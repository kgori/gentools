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

/// It takes a `rust_htslib::bam::Record` and returns a `Result<String, Box<dyn std::error::Error>>`
///
/// Arguments:
///
/// * `record`: a rust_htslib::bam::Record object
///
/// Returns:
///
/// A string containing the fastq record.
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

#[test]
fn test_record_to_fastq_string() {
    let mut record = rust_htslib::bam::Record::new();
    record.set(b"test", None, b"ACGT", &[33, 34, 35, 36]);
    record.set_first_in_template();
    assert_eq!(
        record_to_fastq_string(&record).unwrap(),
        "@test/1".to_string() + "\n" + "ACGT" + "\n" + "+" + "\n" + "BCDE"
    );
    record.unset_first_in_template();
    record.set_last_in_template();
    assert_eq!(
        record_to_fastq_string(&record).unwrap(),
        "@test/2".to_string() + "\n" + "ACGT" + "\n" + "+" + "\n" + "BCDE"
    );
}

/// It reads a BAM file, and writes out three FASTQ files, one for reads that are unmapped, one for
/// reads that are mapped but have an unmapped mate, and one for reads that are unmapped but have a
/// mapped mate
///
/// Arguments:
///
/// * `opts`: Cli
///
/// Returns:
///
/// A Result<(), Box<dyn std::error::Error>>
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

        match (record.is_unmapped(), record.is_mate_unmapped()) {
            (true, true) => {
                writeln!(&mut unmapped_file, "{}", record_to_fastq_string(&record)?)?;
            }
            (false, true) => {
                writeln!(
                    &mut mapped_with_unmapped_mate_file,
                    "{}",
                    record_to_fastq_string(&record)?
                )?;
            }
            (true, false) => {
                writeln!(
                    &mut unmapped_with_mapped_mate_file,
                    "{}",
                    record_to_fastq_string(&record)?
                )?;
            }
            _ => {
                continue;
            }
        }
    }
    Ok(())
}
