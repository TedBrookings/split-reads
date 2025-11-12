use crate::commands::command::Command;
use anyhow::{Result, anyhow};
use clap::{Parser, builder::PossibleValuesParser, value_parser};
use log::info;
use rust_htslib::bam::Writer as BamWriter;
use split_reads::{
    path_type::PathType,
    sam_writer_spec::SamWriterSpec,
    split_index::{SPLIT_INDEX_EXTENSION, SplitIndex},
    util::{RecordType, get_bam_reader, get_fastq_reader, get_fastq_writer},
};
use std::{num::NonZero, path::PathBuf};

/// Index SAM,BAM, or CRAM. Save to split-index (".si") file for rapid extraction of chunks.
#[derive(Parser, Debug)]
#[command(version, verbatim_doc_comment)]
pub(crate) struct Index {
    /// Input SAM/BAM/CRAM to index. Use "-" for stdin.
    #[clap(long, short = 'i', required = true)]
    input: PathBuf,

    /// Output path for Index file. Use "-" for stdout. Defaults to input path with added ".si"
    /// suffix. Error if unspecified and the input file is stdin.
    #[clap(long, short = 'I', required = false, default_value = None)]
    index: Option<PathBuf>,

    /// Reference FASTA (required for CRAMs)
    #[clap(long, short = 'R', required = false, default_value = None)]
    ref_fasta: Option<PathBuf>,

    /// Output path for pass-through SAM.
    #[clap(long, short = 'o', required = false, default_value = None)]
    output: Option<PathBuf>,

    /// Output format type. When specifying file output file names, the extension (.sam, .bam, or
    /// .cram) determines format, so this setting will only have an effect when writing to stdout
    #[clap(long, short = 'O', required = false, default_value_t = String::from("bam"), value_parser = PossibleValuesParser::new(["sam", "bam", "cram", "fastq"]))]
    output_format: String,

    /// Compression level for output compressed formats. Default to 0 for writing to stdout .
    #[clap(long, short = 'C', required = false, value_parser = value_parser!(u32).range(..=9))]
    compression: Option<u32>,

    /// Number of bins to retain in final index file.
    #[clap(long, short = 'n', required = false, default_value_t = NonZero::new(10000usize).unwrap())]
    num_bins: NonZero<usize>,

    /// Number of threads to use for reading BAM
    #[clap(long, short = 't', required = false, default_value_t = NonZero::new(num_cpus::get()).unwrap_or(NonZero::new(1usize).unwrap()))]
    threads: NonZero<usize>,

    /// Time in seconds between log updates
    #[clap(long, required = false, default_value_t = 30)]
    update_interval: u64,
}

impl Index {
    /// Get the output index path that will be used
    fn get_index_path(&self) -> Result<PathBuf> {
        if let Some(specified_index_path) = self.index.clone() {
            // user specified the index path
            Ok(specified_index_path)
        } else if let Some(actual_output_path) = self.output.clone() {
            PathType::from_path(actual_output_path)?
                .default_index(SPLIT_INDEX_EXTENSION)?
                .ok_or_else(|| {
                    anyhow!("When writing to stdout, must explicitly specify index path.")
                })
        } else {
            PathType::from_path(self.input.clone())?
                .default_index(SPLIT_INDEX_EXTENSION)?
                .ok_or_else(|| {
                    anyhow!("When reading from stdin, must explicitly specify index path.")
                })
        }
    }

    /// Get the type of Record that will be used. Check for consistency if writing pass-through.
    fn get_record_type(&self) -> Result<RecordType> {
        let maybe_input_type = RecordType::from_path(self.input.clone());
        let maybe_output_type = if let Some(ref actual_output_path) = self.output {
            RecordType::from_path(actual_output_path)
        } else {
            None
        };
        match (maybe_input_type, maybe_output_type) {
            (Some(input_type), Some(output_type)) => {
                if input_type == output_type {
                    Ok(input_type)
                } else {
                    Err(anyhow!(
                        "Input type ({input_type}) and output type ({output_type}) do not match."
                    ))
                }
            }
            (Some(input_type), None) => Ok(input_type),
            (None, Some(output_type)) => Ok(output_type),
            (None, None) => {
                if self.output_format == "fastq" {
                    Ok(RecordType::Fastq)
                } else {
                    Ok(RecordType::Bam)
                }
            }
        }
    }

    /// Build the split index, then downsize to the requested number of bins and write to requested
    /// index path
    pub fn index_reads(&self) -> Result<PathBuf> {
        // First ensure that the output path is well-specified
        let index_path = self.get_index_path()?;
        let record_type = self.get_record_type()?;

        // Build and downsample the index
        let split_index = if record_type == RecordType::Bam {
            // read (and possibly write) SAM/BAM/CRAM
            let reader = get_bam_reader(self.input.clone(), self.ref_fasta.clone(), self.threads)?;
            let writer: Option<BamWriter> = if let Some(ref output) = self.output {
                Some(
                    SamWriterSpec::new(output)
                        .header_from_reader(&reader)
                        .format_from_path_or_default(self.output_format.clone())?
                        .threads(self.threads)
                        .reference_fasta(self.ref_fasta.clone().as_ref())
                        .compression(self.compression)
                        .get_bam_writer()?,
                )
            } else {
                None
            };
            SplitIndex::build(reader, writer, self.num_bins, self.update_interval)?
        } else {
            // read (and possibly write) FASTQ
            let reader = get_fastq_reader(self.input.clone(), self.threads)?;
            let writer = if let Some(ref output) = self.output {
                Some(get_fastq_writer(output, self.compression, self.threads)?)
            } else {
                None
            };
            SplitIndex::build(reader, writer, self.num_bins, self.update_interval)?
        };
        info!(
            "Indexed {} reads and {} queries into  {} raw bins.",
            split_index.num_reads(),
            split_index.num_queries(),
            split_index.len()
        );
        let downsized_index = split_index.downsize_reads(self.num_bins)?;
        info!("Downsized index to {} bins", downsized_index.len());

        // Write the downsized index
        downsized_index.write(index_path.clone())?;
        Ok(index_path)
    }
}

/// Implement the Command trait for `Index` struct.
impl Command for Index {
    /// Execute the index command to build and write a split-index file.
    fn execute(&self) -> Result<()> {
        info!("Using {} thread(s)", self.threads);
        self.index_reads()?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::{SplitIndex, get_bam_reader};
    use crate::test_utils::random_bam::QueryType;
    use anyhow::Result;
    use rstest::rstest;
    use rust_htslib::bam::Writer as BamWriter;
    use std::{cmp::min, num::NonZero, path::PathBuf};
    use tempfile::TempDir;

    /// Detailed assertiton of expected SplitIndex structure
    fn assert_valid_split_index(
        split_index: &SplitIndex,
        num_reads: usize,
        num_queries: usize,
        requested_num_bins: usize,
        label: &str,
        is_downsized: bool,
    ) {
        // assert number of reads and queries matches
        assert!(
            split_index.num_reads() == num_reads,
            "{label}: expected {num_reads} reads but got {}",
            split_index.num_reads()
        );
        assert!(
            split_index.num_queries() == num_queries,
            "{label}: expected {num_queries} queries but got {}",
            split_index.num_queries()
        );
        if is_downsized {
            // assert number of bins matches request, or number of queries if they are fewer
            let expected_num_bins = min(requested_num_bins, num_queries);
            assert!(
                split_index.len() == expected_num_bins,
                "{label}: actual bins ({}) != expected {expected_num_bins}",
                split_index.len()
            )
        } else {
            // assert number of bins is in allowable range
            let min_num_bins = min(requested_num_bins, num_queries);
            assert!(
                split_index.len() >= min_num_bins,
                "{label}: num actual bins ({}) < min allowed {min_num_bins}",
                split_index.len()
            );
            assert!(
                split_index.len() <= num_queries,
                "{label}: num actual bins ({}) > num queries {num_queries}",
                split_index.len()
            );
        }

        // Check that change from one bin to the next has allowable range of queries
        let split_record_num_queries: Vec<usize> = split_index.get_split_record_num_queries();
        let expected_mean_delta: usize =
            (num_queries as f64 / requested_num_bins as f64).ceil() as usize;
        let (low_delta, high_delta) = if is_downsized {
            (expected_mean_delta / 2, expected_mean_delta * 2)
        } else {
            (1, expected_mean_delta * 2)
        };
        for idx in 1..split_index.len() {
            let delta = split_record_num_queries.get(idx).unwrap()
                - split_record_num_queries.get(idx - 1).unwrap();
            assert!(
                delta >= low_delta,
                "{label}: actual bin delta-queries {delta} < min allowed {low_delta}"
            );
            assert!(
                delta <= high_delta,
                "{label}: actual bin delta-queries {delta} > max allowed {high_delta}"
            );
        }
    }

    /// Struct for holding a test case twith a specific number of bins and queries, and a label
    /// describing it (for debugging purposes)
    struct TestCase {
        num_queries: usize,
        num_bins: usize,
        label: &'static str,
    }

    impl TestCase {
        pub fn new(num_queries: usize, num_bins: usize, label: &'static str) -> Self {
            TestCase {
                num_queries,
                num_bins,
                label,
            }
        }
    }

    /// Test that we can build a SplitIndex and downsize it, and both the raw and downsized
    /// SplitIndices are valid.
    #[rstest(query_type => [QueryType::Single, QueryType::Paired, QueryType::Grouped],
        test_case => [
            TestCase::new(100, 10, "even_divisions"),
            TestCase::new(101, 13, "uneven divisions"),
            TestCase::new(100, 100, "1 read per bin"),
            TestCase::new(100, 1000, "Too many bins"),
            TestCase::new(0, 10, "No reads"),
        ])]
    fn test_index(query_type: QueryType, test_case: TestCase) -> Result<()> {
        let temp_dir = TempDir::new()?;
        let temp_path: PathBuf = temp_dir.path().to_path_buf();
        let (random_bam, num_reads) = query_type.random_bam(&temp_path, test_case.num_queries)?;
        let reader = get_bam_reader(random_bam, None::<PathBuf>, 1usize.try_into()?)?;
        let raw_split_index = SplitIndex::build(
            reader,
            None::<BamWriter>,
            NonZero::new(test_case.num_bins).unwrap(),
            u64::MAX,
        )?;
        assert_valid_split_index(
            &raw_split_index,
            num_reads,
            test_case.num_queries,
            test_case.num_bins,
            &format!("{}, {}, raw index", query_type.label(), test_case.label),
            false,
        );

        let final_split_index = raw_split_index.downsize_reads(test_case.num_bins.try_into()?)?;
        assert_valid_split_index(
            &final_split_index,
            num_reads,
            test_case.num_queries,
            test_case.num_bins,
            &format!(
                "{}, {}, downsized index",
                query_type.label(),
                test_case.label
            ),
            true,
        );
        Ok(())
    }
}
