use crate::commands::command::Command;
use anyhow::{Result, anyhow};
use clap::Parser;
use log::info;
use split_reads::fastq::FastqReader;
use split_reads::maybe_compressed_io::open_file;
use std::{io::BufReader, num::NonZero, path::PathBuf};

/// Index SAM,BAM, or CRAM. Save to split-index (".si") file for rapid extraction of chunks.
#[derive(Parser, Debug)]
#[command(version, verbatim_doc_comment)]
pub(crate) struct TestFastq {
    /// Input Fastq to read. Use "-" for stdin.
    #[clap(long, short = 'i', required = true)]
    input: PathBuf,

    /// Number of threads to use for reading FASTQ
    #[clap(long, short = 't', required = false, default_value_t = NonZero::new(num_cpus::get()).unwrap_or(NonZero::new(1usize).unwrap()))]
    threads: NonZero<usize>,
}

impl TestFastq {
    /// Build the split index, then downsize to the requested number of bins and write to requested
    /// index path
    pub fn test_count_queries(&self) -> Result<()> {
        // First ensure that the output path is well-specified
        let buf = BufReader::new(open_file(self.input.clone(), false)?);
        let mut reader = FastqReader::new(buf);
        let mut num_records: usize = 0;
        let mut num_queries: usize = 0;
        let mut qname = reader.next().ok_or_else(|| anyhow!("No records"))??.name;
        num_records += 1;
        num_queries += 1;
        for record in reader {
            let record = record?;
            num_records += 1;
            if record.name != qname {
                qname = record.name;
                num_queries += 1;
            }
        }
        info!("Read {num_records} reads and {num_queries} queries.",);
        Ok(())
    }
}

/// Implement the Command trait for `TestFastq` struct.
impl Command for TestFastq {
    /// Execute the test-fastq command to count and display FASTQ records and queries.
    fn execute(&self) -> Result<()> {
        info!("Using {} thread(s)", self.threads);
        self.test_count_queries()?;
        Ok(())
    }
}
