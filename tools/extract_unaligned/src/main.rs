extern crate bio;
extern crate bwa;

use bwa::BwaAligner;
use std::str::from_utf8;

use bio::io::fastq;
use clap::Parser;
use num_cpus;
use rust_htslib::{bam, bam::Read};

#[derive(Parser)]
struct Cli {
    #[arg(short, long, value_name = "FILE")]
    bamfile: std::path::PathBuf,

    #[clap(
        short,
        long,
        value_name = "FILES",
        value_delimiter = ',',
        help = "Reference files to filter reads against"
    )]
    references: Option<Vec<std::path::PathBuf>>,

    #[clap(
        short = 'o',
        long,
        value_name = "FILE",
        help = "Name of output fastq file"
    )]
    outfile: std::path::PathBuf,

    #[clap(
        short,
        long,
        default_value_t = 10.0,
        help = "Filter out reads with average base quality smaller than this"
    )]
    min_quality: f32,

    #[clap(
        short = 's',
        long,
        default_value_t = 5000,
        help = "BWA alignment batch size"
    )]
    batch_size: usize,

    #[clap(
        short = 't',
        long,
        value_name = "INT",
        default_value_t = 1,
        help = "Number of reader threads"
    )]
    reader_threads: usize,

    #[clap(
        short,
        long,
        value_name = "INT",
        default_value_t = 1,
        help = "Number of aligner threads"
    )]
    aligner_threads: usize,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let opts = Cli::parse();

    do_work(opts)?;
    Ok(())
}

fn bam_to_fastq(
    record: &rust_htslib::bam::Record,
) -> Result<fastq::Record, Box<dyn std::error::Error>> {
    let name = from_utf8(&record.qname())?;
    let seq = record.seq().as_bytes();
    let qual = record.qual().iter().map(|q| q + 33).collect::<Vec<u8>>();
    Ok(fastq::Record::with_attrs(name, None, &seq, &qual))
}

fn avg_quality(record: &rust_htslib::bam::Record) -> f32 {
    record.qual().iter().map(|x| *x as f32).sum::<f32>() / record.qual().len() as f32
}

#[test]
fn test_bam_to_fastq() {
    let mut bam_rec = rust_htslib::bam::Record::new();
    bam_rec.set(b"test", None, b"ACGT", &[0, 0, 0, 41]);

    let fastq_rec = fastq::Record::with_attrs("test", None, b"ACGT", b"!!!J");

    assert_eq!(bam_to_fastq(&bam_rec).unwrap(), fastq_rec,);
}

#[test]
fn test_avg_quality() {
    let mut bam_rec = rust_htslib::bam::Record::new();
    bam_rec.set(b"test", None, b"ACGT", &[10, 12, 15, 45]);
    assert_eq!(avg_quality(&bam_rec), 20.5);
}

fn do_work(opts: Cli) -> Result<(), Box<dyn std::error::Error>> {
    println!("Opening bam");
    let mut bam = bam::Reader::from_path(opts.bamfile)?;

    let aligners = if let Some(references) = opts.references {
        println!("Opening BWA alignment reference files");
        references
            .iter()
            .map(|r| BwaAligner::from_path(r).unwrap())
            .collect::<Vec<_>>()
    } else {
        Vec::new()
    };

    if opts.reader_threads > 1 {
        let mut extra_threads = opts.reader_threads - 1;
        if extra_threads > num_cpus::get() - 1 {
            extra_threads = num_cpus::get() - 1;
        }
        println!("Using {} threads for reading", extra_threads + 1);
        bam.set_threads(extra_threads)?;
    } else {
        println!("Using 1 thread for reading");
    }

    let mut aligner_threads = opts.aligner_threads;
    if aligner_threads > num_cpus::get() {
        aligner_threads = num_cpus::get();
    }
    if aligner_threads < 1 {
        aligner_threads = 1;
    }
    if aligner_threads == 1 {
        println!("Using 1 thread for alignment");
    } else {
        println!("Using {} threads for alignment", aligner_threads);
    }

    println!("Opening output file");
    let mut unmapped_file = fastq::Writer::to_file(opts.outfile)?;

    let mut unmapped_reads: Vec<fastq::Record> = Vec::new();
    println!("Looking for unmapped reads");
    for r in bam.records() {
        let record = r?;
        if record.is_secondary() || record.is_duplicate() || record.is_supplementary() {
            continue;
        }

        if record.is_unmapped() {
            if avg_quality(&record) >= opts.min_quality {
                unmapped_reads.push(bam_to_fastq(&record)?);
            }
        }

        if unmapped_reads.len() == opts.batch_size {
            println!("{} unmapped reads collected", unmapped_reads.len());
            let mut mapped = vec![false; unmapped_reads.len()];
            for aligner in &aligners {
                let filters =
                    aligner.get_alignment_status(&unmapped_reads, false, false, aligner_threads)?;
                mapped
                    .iter_mut()
                    .zip(filters)
                    .for_each(|(a, b)| *a = *a || b);
            }
            unmapped_reads
                .iter()
                .zip(mapped)
                .filter(|(_, is_mapped)| !is_mapped)
                .map(|(read, _)| read)
                .for_each(|read| {
                    unmapped_file.write_record(&read).unwrap();
                });
            unmapped_file.flush()?;
            unmapped_reads.clear();
        }
    }

    if !unmapped_reads.is_empty() {
        println!("{} unmapped reads collected", unmapped_reads.len());
        let mut mapped = vec![false; unmapped_reads.len()];
        for aligner in &aligners {
            let filters =
                aligner.get_alignment_status(&unmapped_reads, false, false, aligner_threads)?;
            mapped
                .iter_mut()
                .zip(filters)
                .for_each(|(a, b)| *a = *a || b);
        }
        unmapped_reads
            .iter()
            .zip(mapped)
            .filter(|(_, is_mapped)| !is_mapped)
            .map(|(read, _)| read)
            .for_each(|read| {
                unmapped_file.write_record(&read).unwrap();
            });
        unmapped_reads.clear();
    }

    Ok(())
}
