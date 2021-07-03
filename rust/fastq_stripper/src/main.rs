extern crate bio;
extern crate flate2;

use failure::Error;
use docopt::Docopt;
use std::fs::File;
use std::io::{self, prelude::*, BufReader};
use std::collections::HashSet;
//use seq_io::fastq::{Reader,Record};
use bio::io::fastq;
use serde::Deserialize;
use std::path::Path;
use flate2::Compression;
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;


const USAGE: &'static str = "
A tool to extract reads from a fastq given a file containing read names.

Usage:
  extract_reads <in-fastq> <out-fastq> <read-names>
  extract_reads (-h | --help)
Options:
  -h --help                Show this screen.
";

#[derive(Debug, Deserialize, Clone)]
pub struct Args {
    arg_in_fastq: String,
    arg_out_fastq: String,
    arg_read_names: String,
}

/// https://github.com/10XGenomics/vartrix/blob/master/src/main.rs#L708-L723
pub fn read_with_gz<P: AsRef<Path>>(p: P) -> Result<Box<dyn BufRead>, Error> {
    let r = File::open(p.as_ref())?;

    let ext = p.as_ref().extension();

    use std::ffi::OsStr;
    if ext == Some(OsStr::new("gz")) {
        let gz = MultiGzDecoder::new(r);
        let buf_reader = BufReader::new(gz);
        Ok(Box::new(buf_reader))
    } else {
        let buf_reader = BufReader::new(r);
        Ok(Box::new(buf_reader))
    }
}

pub fn write_with_gz<P: AsRef<Path>>(p: P) -> Result<Box<dyn Write>, Error> {
    let r = File::create(p.as_ref())?;

    let ext = p.as_ref().extension();

    use std::ffi::OsStr;
    if ext == Some(OsStr::new("gz")) {
        let gz = GzEncoder::new(r, Compression::default());
        Ok(Box::new(gz))
    } else {
        Ok(Box::new(r))
    }
}

fn main() -> io::Result<()> {
    let args: Args = Docopt::new(USAGE)
        .and_then(|d| d.deserialize())
        .unwrap_or_else(|e| e.exit());
   
    let read_names = Path::new(&args.arg_read_names);
    let in_fastq = Path::new(&args.arg_in_fastq);
    let out_fastq = Path::new(&args.arg_out_fastq);
    println!("Read names: {}", read_names.display());
    println!("Input fastq: {}", in_fastq.display());
    println!("Output fastq: {}", out_fastq.display());
    let mut reads = HashSet::new();
    let reader = read_with_gz(&read_names).unwrap();

    for line in reader.lines() {
        let _fqn1 = line.unwrap();
        let _fqn2 = _fqn1.split_whitespace().next().unwrap().to_string();
        let fqn = _fqn2.strip_prefix("@").unwrap_or(&_fqn2).to_string();
        reads.insert(fqn);
    }

    let reader = fastq::Reader::new(read_with_gz(&in_fastq).unwrap());
    let mut fqout = fastq::Writer::new(write_with_gz(&out_fastq).unwrap());
    
    for result in reader.records() {
        let record = result.unwrap();
        assert!(record.check().is_ok());
        //println!("{}", record.id().unwrap());
        if reads.contains(record.id()) {
            fqout.write_record(&record).expect("writing error");
        }
    }
    Ok(())
}
