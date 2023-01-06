use std::fs::File;
use std::io::Write;
use std::str::from_utf8;

use clap::Parser;
use num_cpus;

#[derive(Parser)]
struct Cli {
    #[arg(short, long, value_name = "FILE")]
    bamfile: std::path::PathBuf,

    #[arg(short = 'u', long, value_name = "FILE")]
    unmapped: std::path::PathBuf,

    #[arg(short = 'm', long, value_name = "FILE")]
    mate_mapped: std::path::PathBuf,

    #[arg(short = 'a', long, value_name = "FILE")]
    mate_unmapped: std::path::PathBuf,

    #[arg(short = 't', long, value_name = "INT", default_value_t = 0)]
    extra_threads: usize,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let opts = Cli::parse();

    using_htslib(opts)?;
    Ok(())
}

fn record_to_fastq_string(
    record: &rust_htslib::bam::Record,
) -> Result<String, Box<dyn std::error::Error>> {
    let mut name = from_utf8(&record.qname())?.to_string();
    let qual = record.qual().iter().map(|q| q + 33).collect::<Vec<u8>>();
    let qual_str = from_utf8(&qual)?;
    if record.is_first_in_template() {
        name += "/1";
    } else if record.is_last_in_template() {
        name += "/2";
    }
    return Ok(format!(
        "@{}\n{}\n+\n{}",
        name,
        from_utf8(&record.seq().as_bytes())?,
        qual_str
    ));
}

fn using_htslib(opts: Cli) -> Result<(), Box<dyn std::error::Error>> {
    use rust_htslib::{bam, bam::Read};
    let mut bam = bam::Reader::from_path(opts.bamfile)?;
    if opts.extra_threads > 0 {
        let mut extra_threads = opts.extra_threads;
        if extra_threads > num_cpus::get() - 1 {
            extra_threads = num_cpus::get() - 1;
        }
        println!("Using {} extra threads", extra_threads);
        bam.set_threads(extra_threads)?;
    }

    let mut unmapped_file = File::create(opts.unmapped)?;
    let mut mapped_with_unmapped_mate_file = File::create(opts.mate_unmapped)?;
    let mut unmapped_with_mapped_mate_file = File::create(opts.mate_mapped)?;

    for r in bam.records() {
        let record = r?;
        if record.is_secondary() || record.is_duplicate() || record.is_supplementary() {
            continue;
        }

        if !(record.is_unmapped() || record.is_mate_unmapped()) {
            continue;
        }

        let fastq_string = record_to_fastq_string(&record)?;

        if record.is_unmapped() && record.is_mate_unmapped() {
            writeln!(&mut unmapped_file, "{}", fastq_string)?;
        } else if record.is_mate_unmapped() && !record.is_unmapped() {
            writeln!(&mut mapped_with_unmapped_mate_file, "{}", fastq_string)?;
        } else if record.is_unmapped() && !record.is_mate_unmapped() {
            writeln!(&mut unmapped_with_mapped_mate_file, "{}", fastq_string)?;
        }
    }
    Ok(())
}
