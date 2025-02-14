// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.

// Read the vars file and export code.  This is a partial implementation.
// Output is {(filename, contents)}.
//
// This writes a temporary file.

use crate::var::*;
use io_utils::*;
use itertools::Itertools;
use std::io::{BufWriter, Write};
use std::process::Command;
use string_utils::*;

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Find upper case strings in var.

fn get_uppers(var: &str) -> Vec<(String, usize)> {
    let mut uppers = Vec::<(String, usize)>::new();
    {
        let mut chars = Vec::<char>::new();
        for c in var.chars() {
            chars.push(c);
        }
        let mut s = String::new();
        let mut start = 0;
        for i in 0..chars.len() {
            if chars[i].is_ascii_uppercase() {
                s.push(chars[i]);
            } else if !s.is_empty() {
                uppers.push((s.clone(), start));
                start = i + 1;
                s.clear();
            }
        }
        if !s.is_empty() {
            uppers.push((s, start));
        }
    }
    if uppers.len() > 1 {
        eprintln!(
            "\nIllegal variable {}, has more than one uppercase string in it.\n",
            var
        );
        std::process::exit(1);
    }
    uppers
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

fn process_var<W: Write>(
    v: &Variable,
    exact: &str,
    cell: &str,
    code: &str,
    f: &mut BufWriter<W>,
    class: &str,
) {
    let var = &v.name;
    let uppers = get_uppers(var);
    let mut rega = false;
    let mut dataset = false;
    let mut name = false;
    let bc = var == "BC";
    let info = var == "INFO";
    for i in 0..uppers.len() {
        if uppers[i].0 == "REGA" {
            rega = true;
        } else if uppers[i].0 == "DATASET" {
            dataset = true;
        } else if uppers[i].0 == "NAME" {
            name = true;
        }
    }
    let upper = !uppers.is_empty();
    if !upper || rega || dataset || name || bc || info {
        let mut passes = 1;
        if v.level == "cell-exact" {
            passes = 2;
        }
        for pass in 1..=passes {
            let mut var = var.clone();
            if pass == 2 {
                var += "_cell";
            }
            emit_code_to_test_for_var(&var, f, class);
            fwriteln!(f, "{}", code);
            if pass == 1 {
                fwriteln!(f, "({}, {}, \"{}\".to_string())", exact, cell, v.level);
            } else {
                fwriteln!(f, "let _exact = {};", exact); // to circumvent warning
                fwriteln!(f, "(String::new(), {}, \"{}\".to_string())", cell, v.level);
            }
        }
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

fn parse_value_return_lines(code: &mut String, level: &str, exact: &mut String, cell: &mut String) {
    *exact = "String::new()".to_string();
    *cell = "Vec::new()".to_string();
    let mut lines = Vec::<String>::new();
    for line in code.lines() {
        lines.push(line.to_string());
    }
    let n = lines.len();
    if n > 0 {
        let mut sub = 0;
        for i in (0..lines.len()).rev() {
            if lines[i].contains("exact: ") {
                *exact = lines[i].after("exact: ").to_string();
                sub += 1;
            } else if lines[i].contains("cell: ") {
                *cell = lines[i].after("cell: ").to_string();
                sub += 1;
            }
        }
        let mut code2 = String::new();
        for i in 0..lines.len() - sub {
            code2 += &mut format!("{}\n", lines[i]);
        }
        *code = code2;
    }
    if level == "cell-exact" {
        assert!(!exact.is_empty());
        assert!(!cell.is_empty());
    }
    if level == "cell" {
        assert!(!cell.is_empty());
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

fn run_rustfmt(f: &str) {
    let new = Command::new("rustfmt")
        .arg(&f)
        .output()
        .unwrap_or_else(|_| panic!("{}", "failed to execute rustfmt".to_string()));
    if new.status.code() != Some(0) {
        eprintln!("\nrustfmt failed\n");
        eprintln!("You can observe the problem by typing rustfmt {}.\n", f);
        std::process::exit(1);
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Emit code that tests for a given variable, allowing for up to three bracket expressions
// in the variable.  See "Test for implemented" before for precisely what is supported.
// also
// ...REGA...
// ...{}...REGA...

fn quote_str_or_char(s: &str) -> String {
    if s.len() == 1 {
        format!("'{}'", s)
    } else {
        format!("\"{}\"", s)
    }
}

fn emit_code_to_test_for_var<W: Write>(var: &str, f: &mut BufWriter<W>, class: &str) {
    let uppers = get_uppers(var);
    let mut rega = None;
    let mut dataset = None;
    let mut name = None;
    let bc = var == "BC";
    let info = var == "INFO";
    for i in 0..uppers.len() {
        if uppers[i].0 == "REGA" {
            rega = Some(uppers[i].1);
        } else if uppers[i].0 == "DATASET" {
            dataset = Some(uppers[i].1);
        } else if uppers[i].0 == "NAME" {
            name = Some(uppers[i].1);
        }
    }
    let nranges = var.matches('{').count();

    // Test for implemented.

    assert_eq!(nranges, var.matches('}').count());
    assert!(uppers.len() <= 1);
    assert!(nranges <= 1 || rega.is_none());
    assert!(nranges == 0 || dataset.is_none());
    assert!(nranges == 0 || name.is_none());
    assert!(nranges <= 3);

    // Proceed.

    if nranges == 0 {
        if rega.is_none() && dataset.is_none() && name.is_none() && !bc && !info {
            fwriteln!(f, r###"}} else if vname == "{}" {{"###, var);
        } else if info && class == "lvar" {
            fwriteln!(
                f,
                r###"}} else if bin_member(&ctl.gen_opt.info_fields, var) {{"###
            );
        } else if bc && class == "lvar" {
            fwriteln!(f, r###"}} else if bin_member(alt_bcs, var) {{"###);
        } else if name.is_some() {
            let (start, stop) = (var.before("NAME"), var.after("NAME"));
            fwriteln!(
                f,
                r###"}} else if vname.starts_with({start_c})
                    && vname.after({start_c}).ends_with({stop_c})
                    && ( bin_member(&ctl.origin_info.dataset_list, 
                        &vname.between2("{start}", "{stop}").to_string())
                        || bin_member(&ctl.origin_info.origin_list, 
                        &vname.between2("{start}", "{stop}").to_string())
                        || bin_member(&ctl.origin_info.donor_list, 
                        &vname.between2("{start}", "{stop}").to_string())
                        || bin_member(&ctl.origin_info.tag_list, 
                        &vname.between2("{start}", "{stop}").to_string()) 
                       ) {{"###,
                start_c = quote_str_or_char(start),
                stop_c = quote_str_or_char(stop),
                start = start,
                stop = stop,
            );
            fwriteln!(
                f,
                r###"let name = vname.between2("{}", "{}");"###,
                start,
                stop
            );
        } else if dataset.is_some() {
            let (start, stop) = (var.before("DATASET"), var.after("DATASET"));
            fwriteln!(
                f,
                r###"}} else if vname.starts_with({start_c})
                    && vname.after("{start}").ends_with({stop_c})
                    && bin_member(&ctl.origin_info.dataset_list, 
                        vname.between2("{start}", "{stop}")) {{"###,
                start_c = quote_str_or_char(start),
                stop_c = quote_str_or_char(stop),
                start = start,
                stop = stop,
            );
            fwriteln!(
                f,
                r###"let dataset = vname.between2("{}", "{}"));"###,
                start,
                stop
            );
        } else {
            let start = var.before("REGA");
            let stop = var.after("REGA");
            fwriteln!(
                f,
                r###"}} else if vname.starts_with({start_c})
                    && vname.after("{start}").ends_with({stop_c})
                    && !vname.between2("{start}", "{stop}").contains('_')
                    && Regex::new(vname.between2("{start}", "{stop}")).is_ok() {{"###,
                start_c = quote_str_or_char(start),
                stop_c = quote_str_or_char(stop),
                start = start,
                stop = stop,
            );
            // Note inefficiency here, as we are instantiating the regular expression for every
            // exact subclonotype, whereas it only needs to be done once (in principle).
            fwriteln!(
                f,
                r###"let reg = Regex::new(vname.between2("{}", "{}")).unwrap();"###,
                start,
                stop
            );
        }
    } else if nranges == 1 {
        if rega.is_none() {
            let begin = var.before("{");
            let end = var.after("}");
            let low = var.after("{").before("..");
            let high = var.after("{").between("..", "}");
            let mut conditions = Vec::<String>::new();
            conditions.push(format!(
                r###"vname.starts_with({})"###,
                quote_str_or_char(begin)
            ));
            conditions.push(format!(
                r###"vname.ends_with({})"###,
                quote_str_or_char(end)
            ));
            conditions.push(format!(
                r###"vname.between2("{}", "{}").parse::<i64>().is_ok()"###,
                begin, end,
            ));
            if !low.is_empty() {
                conditions.push(format!(
                    r###"vname.between2("{}", "{}").force_i64() >= {}"###,
                    begin, end, low,
                ));
            }
            if !high.is_empty() {
                conditions.push(format!(
                    r###"vname.between2("{}", "{}").force_i64() <= {}"###,
                    begin, end, high,
                ));
            }
            fwriteln!(f, "}} else if {} {{ ", conditions.iter().format(" && "));
            fwriteln!(
                f,
                r###"let arg1 = vname.between2("{}", "{}").force_i64();"###,
                begin,
                end,
            );
        } else {
            // <begin>{}<start>REGA<stop>
            let begin = var.before("{");
            let start = var.between("}", "REGA");
            let stop = var.after("REGA");
            let low = var.after("{").before("..");
            let high = var.after("{").between("..", "}");
            let mut conditions = Vec::<String>::new();
            conditions.push(format!(r###"vname.starts_with("{}")"###, begin));
            conditions.push(format!(r###"vname.ends_with("{}")"###, stop));
            conditions.push(format!(
                r###"vname.between2("{}", "{}").contains('{}')"###,
                begin, stop, start,
            ));
            conditions.push(format!(
                r###"vname.between("{}", "{}").parse::<i64>().is_ok()"###,
                begin, start,
            ));
            if !low.is_empty() {
                conditions.push(format!(
                    r###"vname.between("{}", "{}").force_i64() >= {}"###,
                    begin, start, low,
                ));
            }
            if !high.is_empty() {
                conditions.push(format!(
                    r###"vname.between("{}", "{}").force_i64() <= {}"###,
                    begin, start, high,
                ));
            }
            conditions.push(format!(
                r###"!vname.after("{}").between2("{}", "{}").contains('_')"###,
                begin, start, stop,
            ));
            conditions.push(format!(
                r###"Regex::new(vname.between2("{}", "{}")).is_ok()"###,
                start, stop,
            ));
            fwriteln!(f, "}} else if {} {{ ", conditions.iter().format(" && "));
            fwriteln!(
                f,
                r###"let arg1 = vname.between("{}", "{}").force_i64();"###,
                begin,
                start,
            );
            fwriteln!(
                f,
                r###"let reg = Regex::new(vname.after("{}").between2("{}", "{}")).unwrap();"###,
                begin,
                start,
                stop
            );
        }
    } else if nranges == 2 {
        // This code has not been exercised.
        let begin = var.before("{");
        let middle = var.between("}", "{");
        let end = var.rev_after("}").to_string();
        let low1 = var.after("{").before("..");
        let high1 = var.after("{").between("..", "}");
        let low2 = var.rev_after("{").before("..");
        let high2 = var.rev_after("{").between("..", "}");
        let mut conditions = Vec::<String>::new();
        conditions.push(format!(r###"vname.starts_with("{}")"###, begin));
        conditions.push(format!(
            r###"vname.after("{}").contains({})"###,
            begin,
            quote_str_or_char(middle),
        ));
        conditions.push(format!(
            r###"vname.after("{}").after("{}").ends_with("{}")"###,
            begin, middle, end,
        ));
        conditions.push(format!(
            r###"vname.between2("{}", "{}").parse::<i64>().is_ok()"###,
            begin, middle,
        ));
        if !low1.is_empty() {
            conditions.push(format!(
                r###"vname.between2("{}", "{}").force_i64() >= {}"###,
                begin, middle, low1,
            ));
        }
        if !high1.is_empty() {
            conditions.push(format!(
                r###"vname.between2("{}", "{}").force_i64() <= {}"###,
                begin, middle, high1,
            ));
        }
        conditions.push(format!(
            r###"vname.after("{}").between2("{}", "{}").parse::<i64>().is_ok()"###,
            begin, middle, end,
        ));
        if !low2.is_empty() {
            conditions.push(format!(
                r###"vname.after("{}").between2("{}", "{}").force_i64() >= {}"###,
                begin, middle, end, low2,
            ));
        }
        if !high2.is_empty() {
            conditions.push(format!(
                r###"vname.after("{}").between2("{}", "{}").force_i64() <= {}"###,
                begin, middle, end, high2,
            ));
        }
        fwriteln!(f, "}} else if {} {{ ", conditions.iter().format(" && "));
        fwriteln!(
            f,
            r###"let arg1 = vname.between2("{}", "{}").force_i64();"###,
            begin,
            middle,
        );
        fwriteln!(
            f,
            r###"let arg2 = vname.after("{}"),between2("{}", "{}").force_i64();"###,
            begin,
            middle,
            end,
        );
    } else {
        let begin = var.before("{");
        let mid1 = var.between("}", "{");
        let mid2 = var.after("}").between("}", "{");
        let end = var.rev_after("}");
        let low1 = var.after("{").before("..");
        let high1 = var.after("{").between("..", "}");
        let low2 = var.after("{").after("{").before("..");
        let high2 = var.after("{").after("{").between("..", "}");
        let low3 = var.rev_after("{").before("..");
        let high3 = var.rev_after("{").between("..", "}");
        let mut conditions = Vec::<String>::new();
        conditions.push(format!(r###"vname.starts_with("{}")"###, begin));
        conditions.push(format!(
            r###"vname.after("{}").contains({})"###,
            begin,
            quote_str_or_char(mid1),
        ));
        conditions.push(format!(
            r###"vname.after("{}").after("{}").contains({})"###,
            begin,
            mid1,
            quote_str_or_char(mid2),
        ));
        conditions.push(format!(
            r###"vname.after("{}").after("{}").after("{}").ends_with({})"###,
            begin,
            mid1,
            mid2,
            quote_str_or_char(end),
        ));
        conditions.push(format!(
            r###"vname.between("{}", "{}").parse::<i64>().is_ok()"###,
            begin, mid1,
        ));
        if !low1.is_empty() {
            conditions.push(format!(
                r###"vname.between("{}", "{}").force_i64() >= {}"###,
                begin, mid1, low1,
            ));
        }
        if !high1.is_empty() {
            conditions.push(format!(
                r###"vname.between("{}", "{}").force_i64() <= {}"###,
                begin, mid1, high1,
            ));
        }
        if !low2.is_empty() {
            conditions.push(format!(
                r###"vname.after("{}").between("{}", "{}").force_i64() >= {}"###,
                begin, mid1, mid2, low2,
            ));
        }
        if !high2.is_empty() {
            conditions.push(format!(
                r###"vname.after("{}").between("{}", "{}").force_i64() <= {}"###,
                begin, mid1, mid2, high2,
            ));
        }
        conditions.push(format!(
            r###"vname.after("{}").after("{}").between("{}", "{}").parse::<i64>().is_ok()"###,
            begin, mid1, mid2, end,
        ));
        if !low3.is_empty() {
            conditions.push(format!(
                r###"vname.after("{}").after("{}").between("{}", "{}").force_i64() >= {}"###,
                begin, mid1, mid2, end, low3,
            ));
        }
        if !high3.is_empty() {
            conditions.push(format!(
                r###"vname.after("{}").after("{}").between("{}", "{}").force_i64() <= {}"###,
                begin, mid1, mid2, end, high3,
            ));
        }
        fwriteln!(f, "}} else if {} {{ ", conditions.iter().format(" && "));
        fwriteln!(
            f,
            r###"let arg1 = vname.between("{}", "{}").force_i64();"###,
            begin,
            mid1,
        );
        fwriteln!(
            f,
            r###"let arg2 = vname.after("{}").between("{}", "{}").force_i64();"###,
            begin,
            mid1,
            mid2,
        );
        fwriteln!(
            f,
            r###"let arg3 = vname.after("{}").after("{}").between("{}", "{}").force_i64();"###,
            begin,
            mid1,
            mid2,
            end,
        );
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Export code.

pub fn export_code(level: usize) -> Vec<(String, String)> {
    let mut outs = Vec::<(String, String)>::new();

    // Define code start/stop for cvar_vdj.

    let cvar_vdj_start = r###"

        // Copyright (c) 2021 10x Genomics, Inc. All rights reserved.
        // This file is auto-generated by the crate enclone_vars, please do not edit.

        use amino::*;
        use crate::print_utils1::*;
        use crate::print_utils3::*;
        use enclone_core::align_to_vdj_ref::*;
        use enclone_core::defs::*;
        use enclone_core::median::*;
        use enclone_core::opt_d::*;
        use enclone_proto::types::*;
        use itertools::Itertools;
        use stats_utils::*;
        use std::cmp::min;
        use std::collections::HashMap;
        use string_utils::*;
        use vdj_ann::refx::RefData;
        use vector_utils::*;

        pub fn proc_cvar_auto(
            j: usize,
            pass: usize,
            var: &String,
            ex: &ExactClonotype,
            exacts: &Vec<usize>,
            exact_clonotypes: &Vec<ExactClonotype>,
            mid: usize,
            col: usize,
            u: usize,
            rsi: &ColInfo,
            refdata: &RefData,
            dref: &Vec<DonorReferenceItem>,
            ctl: &EncloneControl,
            extra_args: &Vec<String>,
            pcols_sort: &Vec<String>,
            cx: &mut Vec<Vec<String>>,
            varmat: &Vec<Vec<Vec<u8>>>,
            out_data: &mut Vec<HashMap<String, String>>,
            stats: &mut Vec<(String, Vec<String>)>,
            allele_data: &AlleleData,
        ) -> Result<bool, String> {

            let mut vname = var.clone();
            if var.contains(':') {
                vname = var.after(":").to_string();
            }
            let cvars = &ctl.clono_print_opt.cvars;
            let mut abbrc = format!("{}{}", var, col + 1);
            if var.contains(':') {
                abbrc = var.before(":").to_string();
            }
            let val =
            if false {
                (String::new(), Vec::<String>::new(), String::new())

        "###;

    let cvar_vdj_stop = r###"

            } else {
                ("$UNDEFINED".to_string(), Vec::<String>::new(), String::new())
            };
            if val.0 == "$UNDEFINED" {
                Ok(false)
            } else {
                let (exact, cell, _level) = &val;
                let mut varc = format!("{}{}", var, col + 1);
                if !exact.is_empty() {
                    if j < rsi.cvars[col].len() && cvars.contains(var) {
                        cx[col][j] = exact.clone();
                    }
                    if pass == 2
                        && ((!ctl.parseable_opt.pout.is_empty()
                            && (ctl.parseable_opt.pchains == "max"
                                || col < ctl.parseable_opt.pchains.force_usize()))
                            || !extra_args.is_empty())
                    {
                        abbrc = abbrc.replace("_Σ", "_sum");
                        abbrc = abbrc.replace("_μ", "_mean");
                        varc = varc.replace("_Σ", "_sum");
                        varc = varc.replace("_μ", "_mean");
        
                        // Strip escape character sequences from exact.  Can happen in notes, 
                        // maybe other places.
        
                        let mut val_clean = String::new();
                        let mut chars = Vec::<char>::new();
                        let valx = exact.to_string();
                        for c in valx.chars() {
                            chars.push(c);
                        }
                        let mut escaped = false;
                        for l in 0..chars.len() {
                            if chars[l] == '' {
                                escaped = true;
                            }
                            if escaped {
                                if chars[l] == 'm' {
                                    escaped = false;
                                }
                                continue;
                            }
                            val_clean.push(chars[l]);
                        }
        
                        // Proceed.
        
                        // let varc = format!("{}{}", v, col + 1);
                        if pcols_sort.is_empty()
                            || bin_member(pcols_sort, &varc)
                            || bin_member(extra_args, &varc)
                        {
                            out_data[u].insert(abbrc.clone(), val_clean);
                        }
                    }
                    if val.1.is_empty() {
                        stats.push((abbrc, vec![exact.to_string(); ex.ncells()]));
                    } else {
                        stats.push((abbrc, cell.to_vec()));
                    }
                } else if !cell.is_empty()
                    && pass == 2
                    && ((ctl.parseable_opt.pchains == "max"
                        || col < ctl.parseable_opt.pchains.force_usize())
                        || !extra_args.is_empty())
                    && (pcols_sort.is_empty() || bin_member(pcols_sort, &varc))
                {
                    let vals = format!("{}", cell.iter().format(POUT_SEP));
                    out_data[u].insert(abbrc, vals);
                }
                Ok(true)
            }
        }

        "###;

    // Build cvar auto file.

    let actual_out = "enclone_print/src/proc_cvar_auto.rs".to_string();
    let mut temp_out = "enclone_exec/testx/outputs/proc_cvar_auto.rs".to_string();
    let mut vars_loc = "enclone_vars/src/vars".to_string();
    if level == 1 {
        temp_out = format!("../{}", temp_out);
        vars_loc = format!("../{}", vars_loc);
    }
    {
        let mut f = open_for_write_new![&temp_out];
        fwrite!(f, "{}", cvar_vdj_start);
        let vars = std::fs::read_to_string(&vars_loc).unwrap();
        let vars = parse_variables(&vars);
        for v in vars.iter() {
            if v.inputs == "cvar_vdj" {
                let (mut exact, mut cell) = (String::new(), String::new());
                let mut code = v.code.clone();
                parse_value_return_lines(&mut code, &v.level, &mut exact, &mut cell);
                process_var(v, &exact, &cell, &code, &mut f, "cvar");
            }
        }
        fwrite!(f, "{}", cvar_vdj_stop);
    }

    // Rustfmt and save.

    run_rustfmt(&temp_out);
    let f = std::fs::read_to_string(&temp_out).unwrap();
    outs.push((actual_out, f));

    // Define code start/stop for lvar_vdj.

    let lvar_vdj_start = r###"

        // Copyright (c) 2021 10x Genomics, Inc. All rights reserved.
        // This file is auto-generated by the crate enclone_vars, please do not edit.

        use amino::*;
        use enclone_core::defs::*;
        use enclone_core::median::*;
        use enclone_proto::types::*;
        use itertools::Itertools;
        use ndarray::s;
        use regex::Regex;
        use std::cmp::{max, min};
        use std::collections::HashMap;
        use string_utils::*;
        use vdj_ann::refx::RefData;
        use vector_utils::*;
        #[cfg(not(target_os = "windows"))]
        use hdf5x::Reader;
        #[cfg(target_os = "windows")]
        use hdf5::Reader;

        pub fn proc_lvar_auto(
            i: usize,
            pass: usize,
            var: &String,
            exacts: &Vec<usize>,
            exact_clonotypes: &Vec<ExactClonotype>,
            u: usize,
            rsi: &ColInfo,
            refdata: &RefData,
            ctl: &EncloneControl,
            extra_args: &Vec<String>,
            out_data: &mut Vec<HashMap<String, String>>,
            stats: &mut Vec<(String, Vec<String>)>,
            lvars: &Vec<String>,
            row: &mut Vec<String>,
            fate: &Vec<HashMap<String, String>>,
            dref: &Vec<DonorReferenceItem>,
            varmat: &Vec<Vec<Vec<u8>>>,
            fp: &Vec<Vec<usize>>,
            n_vdj_gex: &Vec<usize>,
            vdj_cells: &Vec<Vec<String>>,
            gex_info: &GexInfo,
            groups: &HashMap<usize, Vec<usize>>,
            mults: &Vec<usize>,
            nd_fields: &Vec<String>,
            gex_counts_unsorted: &Vec<usize>,
            gex_fcounts_unsorted: &Vec<f64>,
            n_gexs: &Vec<usize>,
            d_readers: &Vec<Option<Reader>>,
            ind_readers: &Vec<Option<Reader>>,
            h5_data: &Vec<(usize, Vec<u32>, Vec<u32>)>,
            alt_bcs: &Vec<String>,
        ) -> Result<bool, String> {

            let clonotype_id = exacts[u];
            let ex = &exact_clonotypes[clonotype_id];
            let mat = &rsi.mat;
            let cols = varmat[0].len();
            let verbose = ctl.gen_opt.row_fill_verbose;
            let mut vname = var.clone();
            let mut abbr = var.clone();
            if var.contains(':') {
                abbr = var.before(":").to_string();
                vname = var.after(":").to_string();
            }

            macro_rules! speak {
                ($u:expr, $var:expr, $val:expr) => {
                    if pass == 2 && (!ctl.parseable_opt.pout.is_empty() || !extra_args.is_empty()) {
                        let mut v = $var.to_string();
                        v = v.replace("_Σ", "_sum");
                        v = v.replace("_μ", "_mean");
                        if ctl.parseable_opt.pcols.is_empty()
                            || bin_member(&ctl.parseable_opt.pcols_sortx, &v)
                            || bin_member(&extra_args, &v)
                        {
                            out_data[$u].insert(v, $val);
                        }
                    }
                };
            }

            let val =
            if false {
                (String::new(), Vec::<String>::new(), String::new())

        "###;

    let lvar_vdj_stop = r###"

            } else {
                ("$UNDEFINED".to_string(), Vec::<String>::new(), String::new())
            };
            if val.0 == "$UNDEFINED" {
                Ok(false)
            } else {
                let (exact, cell, level) = &val;
                if level == "cell" && !var.ends_with("_cell") {
                    if verbose {
                        eprint!("lvar {} ==> {}; ", var, String::new());
                        eprintln!("i = {}, lvars.len() = {}", i, lvars.len());
                    }
                    if i < lvars.len() {
                        row.push(String::new())
                    }
                    if pass == 2 {
                        speak!(u, abbr, String::new());
                    }
                    stats.push((abbr.to_string(), cell.clone()));
                    if pass == 2 {
                        speak!(u, abbr, format!("{}", cell.iter().format(POUT_SEP)));
                    }
                } else if ( !exact.is_empty() && !var.ends_with("_cell") ) || cell.is_empty() {
                    if verbose {
                        eprint!("lvar {} ==> {}; ", var, exact);
                        eprintln!("i = {}, lvars.len() = {}", i, lvars.len());
                    }
                    if i < lvars.len() {
                        row.push(exact.clone())
                    }
                    if pass == 2 {
                        speak!(u, abbr, exact.to_string());
                    }
                    if cell.is_empty() {
                        stats.push((abbr, vec![exact.to_string(); ex.ncells()]));
                    } else {
                        stats.push((abbr, cell.to_vec()));
                    }
                } else if !cell.is_empty() {
                    if pass == 2 {
                        speak!(u, abbr, format!("{}", cell.iter().format(POUT_SEP)));
                    }
                    stats.push((abbr, cell.to_vec()));
                }
                Ok(true)
            }
        }

        "###;

    // Build lvar auto file.

    let actual_out = "enclone_print/src/proc_lvar_auto.rs".to_string();
    let mut temp_out = "enclone_exec/testx/outputs/proc_lvar_auto.rs".to_string();
    let mut vars_loc = "enclone_vars/src/vars".to_string();
    if level == 1 {
        temp_out = format!("../{}", temp_out);
        vars_loc = format!("../{}", vars_loc);
    }
    {
        let mut f = open_for_write_new![&temp_out];
        fwrite!(f, "{}", lvar_vdj_start);
        let vars = std::fs::read_to_string(&vars_loc).unwrap();
        let vars = parse_variables(&vars);
        for v in vars.iter() {
            if v.inputs == "lvar_vdj" {
                let (mut exact, mut cell) = (String::new(), String::new());
                let mut code = v.code.clone();
                parse_value_return_lines(&mut code, &v.level, &mut exact, &mut cell);
                process_var(v, &exact, &cell, &code, &mut f, "lvar");
            }
        }
        fwrite!(f, "{}", lvar_vdj_stop);
    }

    // Rustfmt and save.

    run_rustfmt(&temp_out);
    let f = std::fs::read_to_string(&temp_out).unwrap();
    outs.push((actual_out, f));

    // Return.

    outs
}
