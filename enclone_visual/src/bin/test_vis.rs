// Copyright (c) 2022 10X Genomics, Inc. All rights reserved.
//
// Run some tests of enclone visual on a Mac.
//
// As part of the tests, this opens a window.  Ideally we would instead run this without opening
// a window, and that might be eventually possible, given changes on the roadmap for iced.
//
// The current images vary some from laptop to laptop, depending on hardware, OS, and personal
// configuration.  Because of this, perhaps only one person can run this test "officially".
// However other people can create their own regression test results, and test against those, as
// described below.
//
// Create mode.  This mode is invoked by adding the argument CREATE.  This creates PNG images,
// which exist only locally, and are not in git.  When the enclone repo is initialized, this needs
// to be run.  You might also want to run it if you somehow mess up the images.  Otherwise, do
// not run in this mode.  Create mode also creates JPG images.
//
// Create PNG mode.  Similar, invoked by CREATE_PNG, only creates PNG files.
//
// Unofficial mode.  This mode is invoked by adding the argument UNOFFICIAL.  This does not read or
// write files in regression_images, some of which are in git, and instead uses a directory
// unofficial_images.  If you are running "unofficially" you need this.
//
// Local mode.  This mode is invoked by adding the argument LOCAL.  It does not run the remote
// tests, which is currently only possible at 10x.  Local mode is also easier to run because it
// does not require compiling on two machines.
//
// Update mode.  This mode is invoked by adding the argument UPDATE.  This causing failing results
// to be replaced.
//
// If you are outside 10x, then the first time you run you should:
// test_vis CREATE UNOFFICIAL LOCAL
// and thereafter:
// test_vis UNOFFICIAL LOCAL.
// This will allow you to tell if, with fairly high probability, changes you make have an effect
// on enclone visual.  However they cannot tell you with with certainty, and if you deliberately
// make a change, there is no mechanism for you to create the jpg files that need to be updated.
//
// Other arguments:
// QUIET.
// VERBOSE.
// TESTS=... (comma-separated list, subset of 1,2,3,4,5,main).
// PRINTER.
//
//
// This code works by comparing lowest resolution JPEG files.  We use that format to avoid
// having larger files in git.  A better solution would be to use lowest resolution
// JPEG2000 files, which would be even smaller.
//
// The tests use only public data.
//
// Part of these tests assume that you can connect to a server.  Also we test receiving shares,
// and this could fail and mess things up if someone else happened to send a share.
//
// See also show_diffs.

use enclone_core::tilde_expand_me;
use enclone_visual::compare_images::*;
use enclone_visual::messages::*;
use enclone_visual::testsuite::{metatests, TESTS};
use enclone_visual::*;
use image::codecs::jpeg::JpegEncoder;
use image::ColorType::Rgba8;
use io_utils::*;
use jpeg_decoder::Decoder;
use perf_stats::*;
use pretty_trace::*;
use std::env;
use std::fs::{copy, File};
use std::io::{BufReader, BufWriter, Read};
use std::process::Command;
use std::time::Instant;
use string_utils::*;

fn main() {
    PrettyTrace::new().on();
    let tall = Instant::now();
    let args: Vec<String> = env::args().collect();
    let mut update = false;
    let mut quiet = false;
    let mut verbose = false;
    let mut printer = false;
    let mut local = false;
    let mut create = false;
    let mut create_png = false;
    let mut unofficial = false;
    let mut tests = Vec::<String>::new();
    for i in 1..args.len() {
        if args[i] == "UPDATE" {
            update = true;
        } else if args[i] == "QUIET" {
            quiet = true;
        } else if args[i] == "VERBOSE" {
            verbose = true;
        } else if args[i] == "PRINTER" {
            printer = true;
        } else if args[i] == "LOCAL" {
            local = true;
        } else if args[i] == "CREATE" {
            create = true;
        } else if args[i] == "CREATE_PNG" {
            create_png = true;
        } else if args[i] == "UNOFFICIAL" {
            unofficial = true;
        } else if args[i].starts_with("TESTS=") {
            tests = args[i]
                .after("TESTS=")
                .split(',')
                .map(str::to_owned)
                .collect();
        } else {
            eprintln!("\nUnknown argument {}.\n", args[i]);
            std::process::exit(1);
        }
    }
    if !path_exists("enclone_visual") {
        eprintln!("\nYou need to run this from the top level directory of the enclone repo.\n");
        std::process::exit(1);
    }
    let mut all_testnames = Vec::<String>::new();
    let options = fs_extra::dir::CopyOptions::new();
    let source = "enclone_visual/tests/sample_visual";
    let target = "enclone_visual/outputs/sample_visual";
    if path_exists(&target) {
        fs_extra::dir::remove("enclone_visual/outputs/sample_visual").unwrap();
    }
    std::fs::create_dir_all("enclone_visual/outputs").unwrap();
    fs_extra::dir::copy(&source, "enclone_visual/outputs", &options).unwrap();

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
    // TEST THAT REMOTE SHARE DIR HAS WORLD PERMISSIONS
    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // This assumes that you have defined ENCLONE_CONFIG_DEFAULT.
    //
    // The reason for having this test is that if the remote share dir does not have world
    // permissions, then some people may not be able to share.  And the permissions on the remote
    // share dir might get accidentally changed by git.

    if !local && tests.is_empty() {
        let mut found = false;
        for (key, value) in env::vars() {
            if key == "ENCLONE_CONFIG_DEFAULT" {
                found = true;
                let host = value.before(":");
                let rdir = value.after(":").rev_before("/");
                let o = Command::new("ssh")
                    .arg(&host)
                    .arg(&format!("ls -l {}", rdir))
                    .output()
                    .expect("failed to execute ssh for world permissions test");
                let out = strme(&o.stdout);
                let mut found2 = false;
                for line in out.lines() {
                    if line.ends_with(" share") {
                        found2 = true;
                        if line.as_bytes()[7..10].to_vec() != b"rwx".to_vec() {
                            eprintln!("\nThe remote share dir does not have world permissions.\n");
                            std::process::exit(1);
                        }
                    }
                }
                if !found2 {
                    eprintln!("\nFailed to find remote share dir in {}:{}.\n", host, rdir);
                    std::process::exit(1);
                }
            }
        }
        if !found {
            eprintln!("\nYou need to define ENCLONE_CONFIG_DEFAULT.\n");
            std::process::exit(1);
        }
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
    // A UNIT TEST
    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    test_fold();

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
    // RUN A COUPLE TESTS IN LOCAL MODE
    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    for n in [4, 5, 6].iter() {
        if printer {
            println!("running test {}", n);
        }
        let n = *n;
        if tests.is_empty() || tests.contains(&format!("{}", n)) {
            let metas = metatests()[n - 1].clone();
            let mut testnames = Vec::<String>::new();
            for m in metas.iter() {
                match m {
                    Message::SetName(x) => {
                        testnames.push(x.to_string());
                    }
                    _ => {}
                };
            }
            for t in testnames.iter() {
                let png = format!("enclone_visual/outputs/{}.png", t);
                if path_exists(&png) {
                    std::fs::remove_file(&png).unwrap();
                }
            }
            let o = Command::new("enclone")
                .arg(&"VIS")
                .arg(&format!("META={}", n))
                .arg(&"VISUAL_DIR=enclone_visual/outputs/sample_visual")
                .output()
                .expect(&format!("failed to execute enclone visual metatest {}", n));
            if o.status.code() != Some(0) {
                eprintln!("\nnonzero exit code from enclone visual metatest {}\n", n);
                eprintln!("stderr =\n{}", strme(&o.stderr));
                std::process::exit(1);
            }
            all_testnames.append(&mut testnames);
        }
    }

    // Reset sample_visual.

    if path_exists(&target) {
        fs_extra::dir::remove("enclone_visual/outputs/sample_visual").unwrap();
    }
    fs_extra::dir::copy(&source, "enclone_visual/outputs", &options).unwrap();

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
    // RUN ARCHIVE TESTS IN LOCAL MODE
    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    if tests.is_empty() || tests.contains(&"1".to_string()) {
        let metas = metatests()[0].clone();
        let mut testnames = Vec::<String>::new();
        for m in metas.iter() {
            match m {
                Message::SetName(x) => {
                    testnames.push(x.to_string());
                }
                _ => {}
            };
        }
        for t in testnames.iter() {
            let png = format!("enclone_visual/outputs/{}.png", t);
            if path_exists(&png) {
                std::fs::remove_file(&png).unwrap();
            }
        }
        let o = Command::new("enclone")
            .arg(&"VIS")
            .arg(&"META=1")
            .arg(&"VISUAL_DIR=enclone_visual/outputs/sample_visual")
            .output()
            .expect("failed to execute enclone visual metatest 1");
        if o.status.code() != Some(0) {
            eprintln!("\nnonzero exit code from enclone visual metatest 1\n");
            eprintln!("stderr =\n{}", strme(&o.stderr));
            std::process::exit(1);
        }
        all_testnames.append(&mut testnames);
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
    // RUN ARCHIVE TESTS IN REMOTE MODE
    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // This is dangerous, because it tests receiving shares, and someone else could send a share
    // before or during the test.  Should change to restrict to receiving only shares from self.

    if !local && (tests.is_empty() || tests.contains(&"2".to_string())) {
        if path_exists(&target) {
            fs_extra::dir::remove("enclone_visual/outputs/sample_visual").unwrap();
        }
        fs_extra::dir::copy(&source, "enclone_visual/outputs", &options).unwrap();
        let metas = metatests()[1].clone();
        let mut testnames = Vec::<String>::new();
        for m in metas.iter() {
            match m {
                Message::SetName(x) => {
                    testnames.push(x.to_string());
                }
                _ => {}
            };
        }
        for t in testnames.iter() {
            let png = format!("enclone_visual/outputs/{}.png", t);
            if path_exists(&png) {
                std::fs::remove_file(&png).unwrap();
            }
        }
        let o = Command::new("enclone")
            .arg(&"VIS=b")
            .arg(&"META=2")
            .arg(&"REQUIRE_COMPATIBLE")
            .arg(&"VISUAL_DIR=enclone_visual/outputs/sample_visual")
            .output()
            .expect("failed to execute enclone visual metatest 2");
        if o.status.code() != Some(0) {
            eprintln!("\nnonzero exit code from enclone visual metatest 2\n");
            eprintln!("stderr =\n{}", strme(&o.stderr));
            std::process::exit(1);
        }
        all_testnames.append(&mut testnames);
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
    // ANOTHER TEST
    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    if printer {
        println!("starting another test");
    }
    if tests.is_empty() || tests.contains(&"3".to_string()) {
        let metas = metatests()[2].clone();
        let mut testnames = Vec::<String>::new();
        for m in metas.iter() {
            match m {
                Message::SetName(x) => {
                    testnames.push(x.to_string());
                }
                _ => {}
            };
        }
        for t in testnames.iter() {
            let png = format!("enclone_visual/outputs/{}.png", t);
            if path_exists(&png) {
                std::fs::remove_file(&png).unwrap();
            }
        }
        let o = Command::new("enclone")
            .arg(&"VIS")
            .arg("FAIL_ON_ERROR")
            .arg(&"META=3")
            .arg(&"VISUAL_DIR=enclone_visual/outputs/sample_visual")
            .output()
            .expect("failed to execute enclone visual metatest 3");
        if o.status.code() != Some(0) {
            eprintln!("\nnonzero exit code from enclone visual metatest 3\n");
            eprintln!("stderr =\n{}", strme(&o.stderr));
            std::process::exit(1);
        }
        all_testnames.append(&mut testnames);
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
    // PRETEST
    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Run enclone once to get it in cache.  This doesn't totally make sense but seems to improve
    // the reproducibility of timing of the actual work.

    if printer {
        println!("starting pretest");
    }
    if tests.is_empty() {
        let o = Command::new("enclone")
            .arg(&"--version")
            .output()
            .expect("failed to execute enclone visual pretest");
        if o.status.code() != Some(0) {
            eprintln!("\nnonzero exit code from enclone visual pretest\n");
            eprintln!("stderr =\n{}", strme(&o.stderr));
            std::process::exit(1);
        }
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
    // RUN MAIN TESTS
    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    if printer {
        println!("running main tests");
    }
    let t = Instant::now();
    if tests.is_empty() || tests.contains(&"main".to_string()) {
        let o = Command::new("enclone")
            .arg(&"VIS")
            .arg(&"TEST")
            .arg(&"VISUAL_DIR=enclone_visual/tests/sample_visual")
            .output()
            .expect("failed to execute enclone visual test");
        if o.status.code() != Some(0) {
            eprintln!("\nnonzero exit code from enclone visual test\n");
            eprintln!("stderr =\n{}", strme(&o.stderr));
            std::process::exit(1);
        }
        if !quiet {
            print!("{}", strme(&o.stdout));
        }
        for i in 0..TESTS.len() {
            if TESTS[i].2.len() > 0 {
                all_testnames.push(TESTS[i].2.to_string());
            }
        }
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
    // CHECK TESTS VERSUS REGRESSION OUTPUTS
    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    let mut fail = false;
    const MAX_DIFFS: usize = 450;
    let reg = if !unofficial {
        "regression_images"
    } else {
        "unofficial_images"
    };
    if unofficial {
        std::fs::create_dir_all("enclone_visual/unofficial_images").unwrap();
    }
    for i in 0..all_testnames.len() {
        let mut image_new = Vec::<u8>::new();
        let old_png_file = format!("enclone_visual/{reg}/{}.png", all_testnames[i]);
        let new_png_file = format!("enclone_visual/outputs/{}.png", all_testnames[i]);
        let mut f = File::open(&new_png_file).unwrap();
        f.read_to_end(&mut image_new).unwrap();
        let (header, image_data_new0) = png_decoder::decode(&image_new).unwrap();
        let (width, height) = (header.width as usize, header.height as usize);

        // Convert the new png file to a jpg, and read that back in as a bit image.

        let new_jpg_file = format!("{}.jpg", new_png_file.rev_before(".png"));
        {
            let quality = 1 as u8; // lowest quality
            if path_exists(&new_jpg_file) {
                // This file removal is to circumvent a bug in Carbon Black.
                std::fs::remove_file(&new_jpg_file).unwrap();
            }
            let mut f = open_for_write_new![&new_jpg_file];
            let mut buff = BufWriter::new(&mut f);
            let mut encoder = JpegEncoder::new_with_quality(&mut buff, quality);
            encoder
                .encode(&image_data_new0, width as u32, height as u32, Rgba8)
                .unwrap();
        }
        let file = open_for_read![&new_jpg_file];
        let mut decoder_new = Decoder::new(BufReader::new(file));
        let image_data_new = decoder_new.decode().expect("failed to decode image");

        // In create mode, just copy file.

        let old_jpg_file = format!("{}.jpg", old_png_file.rev_before(".png"));
        if create || create_png {
            copy(&new_png_file, &old_png_file).unwrap();
            if create {
                copy(&new_jpg_file, &old_jpg_file).unwrap();
            }
            continue;
        }

        // Check for existence of old jpg file.

        if !path_exists(&old_jpg_file) {
            eprintln!(
                "\nLooks like you've added a test.  Please look at \
                enclone_visual/outputs/{}.png and\n\
                if it's right, copy it and the jpg file to regression_tests and git add the jpg.\n",
                all_testnames[i],
            );
            std::process::exit(1);
        }

        // Read in the old jpg file as a bit image.

        let file = open_for_read![&old_jpg_file];
        let mut decoder_old = Decoder::new(BufReader::new(file));
        let image_data_old = decoder_old.decode().expect("failed to decode image");

        // Test for differences.

        if image_data_old.len() != image_data_new.len() {
            eprintln!(
                "\nimage size for {} = {} changed from {} x {} to {} x {}",
                i,
                all_testnames[i],
                decoder_old.info().unwrap().width,
                decoder_old.info().unwrap().height,
                decoder_new.info().unwrap().width,
                decoder_new.info().unwrap().height,
            );
            fail = true;
        } else {
            let diffs = compare_images(&image_data_old, &image_data_new, width, height, verbose);
            if diffs > MAX_DIFFS {
                eprintln!(
                    "\nThere are {} diffs for {}.  Please open enclone_visual/outputs/\
                    joint.{}.jpg.",
                    diffs, all_testnames[i], all_testnames[i]
                );

                // Create and save concatenated image.  Note that we're depending on the png
                // in regression_images as being current.  We could do the same thing with the jpg
                // versions, which would be guaranteed to be correct, but of lower quality.  We did
                // this in 8949a053552c478af5c952ee407416d0e52ab8a0 of dj/189, if you want to go
                // back to that.

                if !path_exists(&old_png_file) {
                    eprintln!(
                        "\nThe file\n{old_png_file}\nis missing.  If you have never run \
                        test_vis before, then probably you should run either\n\
                        test_vis CREATE [if you're at 10x]\nor\n\
                        test_vis CREATE LOCAL [if not]\n"
                    );
                    std::process::exit(1);
                }
                let mut f = File::open(&old_png_file).unwrap();
                let mut image_old = Vec::<u8>::new();
                f.read_to_end(&mut image_old).unwrap();
                let (_, image_data_old0) = png_decoder::decode(&image_old).unwrap();
                if image_data_old0.len() != width * height * 4 {
                    eprintln!(
                        "\nThe size of {} is wrong, so it is probably corrupted.\n",
                        old_png_file
                    );
                    eprintln!(
                        "You may want to checkout master and run with CREATE_PNG to \
                        regenerate valid PNG files.\n"
                    );
                    std::process::exit(1);
                }
                let mut joint = Vec::<u8>::new();
                for i in 0..height {
                    let start = i * width * 4;
                    let stop = (i + 1) * width * 4;
                    joint.append(&mut image_data_old0[start..stop].to_vec());
                    joint.append(&mut image_data_new0[start..stop].to_vec());
                    for j in start..stop {
                        let diff = image_data_new0[j].wrapping_sub(image_data_old0[j]);
                        joint.push(255 - diff);
                    }
                }
                let new_jpg_file = format!("enclone_visual/outputs/joint.{}.jpg", all_testnames[i]);
                let quality = 80 as u8;
                let mut f = open_for_write_new![&new_jpg_file];
                let mut buff = BufWriter::new(&mut f);
                let mut encoder = JpegEncoder::new_with_quality(&mut buff, quality);
                encoder
                    .encode(&joint, (width * 3) as u32, height as u32, Rgba8)
                    .unwrap();

                // Keep going.

                fail = true;
                if update {
                    copy(&new_png_file, &old_png_file).unwrap();
                    copy(&new_jpg_file, &old_jpg_file).unwrap();
                }
            }
        }
    }
    let state = if fail {
        "UNSUCCESSFULLY"
    } else {
        "successfully"
    };

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
    // TEST COMPUTATIONAL PERFORMANCE
    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    let used = elapsed(&t);
    if tests.is_empty() {
        const EXPECTED_TIME: f64 = 49.4; // this is supposed to be the lowest observed value
        const MAX_PERCENT_OVER: f64 = 4.2;
        let percent_over = 100.0 * (used - EXPECTED_TIME) / EXPECTED_TIME;
        if percent_over > MAX_PERCENT_OVER {
            eprintln!(
                "\nUsed {:.1} seconds, exceeding expected test time of {:.1} seconds by {:.1}%, \
                    versus max allowed = {}%.",
                used, EXPECTED_TIME, percent_over, MAX_PERCENT_OVER,
            );
            eprintln!(
                "You might want to retry a second time.  Note that if you're disconnected from\n\
                the internet, then the Mac Gatekeeper could introduce very long delays.\n"
            );
            std::process::exit(1);
        }

        // Get peak memory.  This is sensitive to changes of a few percent or so.

        #[cfg(not(target_os = "windows"))]
        {
            let maxrss_children;
            unsafe {
                let mut rusage: libc::rusage = std::mem::zeroed();
                let retval = libc::getrusage(libc::RUSAGE_CHILDREN, &mut rusage as *mut _);
                assert_eq!(retval, 0);
                maxrss_children = rusage.ru_maxrss;
            }
            let peak_mem_mb = maxrss_children as f64 / ((1024 * 1024) as f64);
            const MAX_PEAK_MEM: f64 = 168.2; // this is supposed to be the lowest observed value
            const MAX_PERCENT_OVER_MEM: f64 = 17.6;
            let percent_over = 100.0 * (peak_mem_mb - MAX_PEAK_MEM) / MAX_PEAK_MEM;

            eprintln!(
                "\nPeak mem {:.1} MB, exceeded expected peak mem of {:.1} MB by {:.1}%, \
                    versus max allowed = {}%.",
                peak_mem_mb, MAX_PEAK_MEM, percent_over, MAX_PERCENT_OVER_MEM,
            );
            if percent_over > MAX_PERCENT_OVER_MEM {
                eprintln!("That's too high.  This happens occasionally, so please retry.\n");
                eprintln!(
                    "Also the value is probably different if you're using a Mac with an \
                    Apple chip.\n"
                );
                std::process::exit(1);
            }
        }
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
    // TEST FOR STRAY PROCESSES
    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    if tests.is_empty() {
        let o = Command::new("ps")
            .arg(&"ux")
            .output()
            .expect("failed to execute ps");
        let out = strme(&o.stdout);
        for line in out.lines() {
            if line.contains("ssh -p") {
                eprintln!(
                    "\nLooks like you may have a stray process running:\n{}\n\n\
                    Treating this as an error.\n",
                    line
                );
                std::process::exit(1);
            }
        }
        for (key, value) in env::vars() {
            if key == "ENCLONE_CONFIG" {
                let host = value.before(":");
                let o = Command::new("ssh")
                    .arg(&host)
                    .arg(&"ps x")
                    .output()
                    .expect("failed to execute ssh");
                let out = strme(&o.stdout);
                for line in out.lines() {
                    if line.contains("enclone SERVER") {
                        eprintln!(
                            "\nLooks like you may have a stray process running on the \
                            remote server:\n{}\n\n\
                            Treating this as an error.\n",
                            line
                        );
                        std::process::exit(1);
                    }
                }
            }
        }
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
    // MAC TESTS THAT DON'T REALLY BELONG HERE BUT DON'T HAVE A BETTER HOME
    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Test that a particular case of tilde expansion works on a Mac.

    if tests.is_empty() {
        let o = Command::new("enclone")
            .arg("BCR=123085")
            .arg("MIN_CELLS=10")
            .arg("PLOT_BY_ISOTYPE=gui")
            .arg("HONEY_OUT=~/enclone_temp_test_file")
            .arg("NOPRINT")
            .output()
            .expect("failed to execute enclone tilde test");
        if o.status.code() != Some(0) {
            eprintln!("\nnonzero exit code from enclone tilde test\n");
            eprintln!("stderr =\n{}", strme(&o.stderr));
            std::process::exit(1);
        }
        let mut f = "~/enclone_temp_test_file".to_string();
        tilde_expand_me(&mut f);
        std::fs::remove_file(&f).unwrap();
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
    // DONE, REPORT STATUS
    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    println!(
        "\nenclone visual tests completed {} in {:.1} seconds",
        state, used,
    );
    println!("actual total time = {:.1} seconds\n", elapsed(&tall));
    if fail {
        std::process::exit(1);
    }
}
