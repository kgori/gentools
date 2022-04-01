use clap::Parser;
use rust_htslib::bam::{Read, Reader, Record};

#[derive(Parser)]
struct Cli {
    #[clap(parse(from_os_str))]
    bamfile: std::path::PathBuf,
}

struct Region {
    tid: i32,
    start: i64,
    end: i64,
    count: i64,
}

impl Region {
    fn new() -> Self {
        Self {
            tid: 0,
            start: 0,
            end: 0,
            count: 0,
        }
    }

    fn print(&self) {
        println!("{}\t{}\t{}\t{}", self.tid, self.start, self.end, self.count)
    }

    fn print_with_header(&self, header: &rust_htslib::bam::HeaderView) {
        let chr = std::string::String::from_utf8_lossy(header.tid2name(self.tid as u32));
        let max_end = header
            .target_len(self.tid as u32)
            .unwrap_or(self.end as u64);
        let end = std::cmp::min(max_end, self.end as u64);
        println!("{}\t{}\t{}\t{}", chr, self.start, end, self.count)
    }
}

fn main() -> Result<(), anyhow::Error> {
    let args = Cli::parse();
    let mut bam = Reader::from_path(args.bamfile)?;
    let header = bam.header().clone();
    let mut record = Record::new();
    let window_size = 1000;

    let mut region = Region::new();
    region.end = window_size;

    while let Some(result) = bam.read(&mut record) {
        match result {
            Ok(_) => {
                let tid = record.tid();
                let pos = record.pos();

                if tid != region.tid {
                    region.print_with_header(&header);
                    region.tid = tid;
                    region.start = 0;
                    region.end = window_size;
                    region.count = 0;
                }

                if pos > region.end {
                    region.print_with_header(&header);
                    region.start = region.end;
                    region.end += window_size;
                    region.count = 0;
                }

                region.count += 1;
            }
            Err(error) => {
                return Err(error.into());
            }
        }
    }
    region.print_with_header(&header);
    Ok(())
}
