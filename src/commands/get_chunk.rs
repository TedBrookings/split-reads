use crate::commands::command::Command;
use anyhow::{Result, anyhow};
use clap::{Parser, builder::PossibleValuesParser, value_parser};
use log::{info, warn};
use split_reads::{
    chunkable::ChunkableRecordReader,
    path_type::PathType,
    sam_writer_spec::SamWriterSpec,
    split_index::{SPLIT_INDEX_EXTENSION, SplitIndex},
    util::{RecordType, get_bam_reader, get_fastq_reader, get_fastq_writer},
};
use std::{
    num::NonZero,
    path::{Path, PathBuf},
};

/// Rapidly extract a chunk from a SAM, BAM, or CRAM that has a split-index (".si") file.
#[derive(Parser, Debug)]
#[command(version, verbatim_doc_comment)]
pub(crate) struct GetChunk {
    /// Input SAM/BAM/CRAM to extract from. Cannot read from stdin, because it is not seekable.
    #[clap(long, short = 'i', required = true)]
    input: PathBuf,

    /// Index for input SAM/BAM/CRAM, built by split-reads index. Use "-" for stdin. Defaults to
    /// input sam path with extra ".si" extension.
    #[clap(long, short = 'I', required = false, default_value = None)]
    index: Option<PathBuf>,

    /// Reference FASTA (required for CRAMs)
    #[clap(long, short = 'R', required = false, default_value = None)]
    ref_fasta: Option<PathBuf>,

    /// Output path for chunk file. Use "-" (or omit) for stdout.
    #[clap(long, short = 'o', required = false, default_value = "-")]
    output: PathBuf,

    /// Compression level for output compressed formats. Default to 0 for writing to stdout .
    #[clap(long, short = 'C', required = false, value_parser = value_parser!(u32).range(..=9))]
    compression: Option<u32>,

    /// Index of chunk to take (0, 1, ..., num_chunks - 1)
    #[clap(long, short = 'c', required = true)]
    chunk_index: usize,

    /// Number of chunks in total input file.
    #[clap(long, short = 'n', required = true)]
    num_chunks: NonZero<usize>,

    /// Output format type. When specifying file output file names, the extension (.sam, .bam, .cram, or .fastq)
    /// determines format, so this setting will only have an effect when writing to stdout. If left unspecified,
    /// use the same format as input.
    #[clap(long, short = 'O', required = false, default_value = None, value_parser = PossibleValuesParser::new(["sam", "bam", "cram", "fastq"]))]
    output_format: Option<String>,

    /// Number of threads to use for reading or writing BAM
    #[clap(long, short = 't', default_value_t = NonZero::new(num_cpus::get()).unwrap_or(NonZero::new(1usize).unwrap()))]
    threads: NonZero<usize>,
}

impl GetChunk {
    /// Load the SplitIndex for the original reads file
    fn load_split_index<P1, P2>(index: Option<P1>, input: P2) -> Result<SplitIndex>
    where
        P1: AsRef<Path>,
        P2: AsRef<Path>,
    {
        if let Some(path_buf) = index {
            SplitIndex::read(path_buf)
        } else {
            let default = PathType::from_path(input)?
                .default_index(SPLIT_INDEX_EXTENSION)?
                .ok_or_else(|| {
                    anyhow!("When reading from stdin, must explicitly specify index path.")
                })?;
            SplitIndex::read(default)
        }
    }

    /// Determine the record type to use for output.
    ///
    /// Uses the output path extension if available, falls back to the output_format
    /// option if specified, otherwise uses the input record type.
    fn get_output_record_type(&self, input_record_type: &RecordType) -> Result<RecordType> {
        if let Some(record_type) = RecordType::from_path(self.input.clone()) {
            Ok(record_type)
        } else if let Some(ref type_string) = self.output_format {
            RecordType::from_extension(Some(type_string))
                .ok_or_else(|| anyhow!("Should be unreachable."))
        } else {
            Ok(input_record_type.clone())
        }
    }

    /// Skip to the beginning of the requested chunk, then write the chunk to the desired output.
    fn write_chunk(&self) -> Result<()> {
        // Load SplitIndex
        let split_index = Self::load_split_index(self.index.clone(), self.input.clone())?;

        // get input record type
        let input_record_type = RecordType::from_path(self.input.clone()).ok_or_else(|| {
            anyhow!("Input type must be FASTQ or SAM/BAM/CRAM. Cannot read from stdin.")
        })?;
        // get output record type
        let output_record_type = self.get_output_record_type(&input_record_type)?;

        if input_record_type == RecordType::Bam {
            // reading from SAM/BAM/CRAM
            let mut reader =
                get_bam_reader(self.input.clone(), self.ref_fasta.clone(), self.threads)?;
            if output_record_type == RecordType::Bam {
                // Reading from SAM/BAM/CRAM and writing to SAM/BAM/CRAM
                let default_format = if let Some(ref output_format) = self.output_format {
                    output_format.clone()
                } else {
                    self.input
                        .clone()
                        .extension()
                        .ok_or_else(|| anyhow!("Input has no extension."))?
                        .to_str()
                        .ok_or_else(|| anyhow!("Input extension cannot convert to str"))?
                        .to_ascii_lowercase()
                };
                let writer_spec = SamWriterSpec::new(self.output.clone())
                    .header_from_reader(&reader)
                    .format_from_path_or_default(default_format)?
                    .threads(self.threads)
                    .reference_fasta(self.ref_fasta.clone())
                    .compression(self.compression)
                    .to_owned();
                let mut writer = writer_spec.get_bam_writer()?;
                // Write the chunk
                let mut fast_forward_info =
                    reader.fast_forward(split_index, self.chunk_index, self.num_chunks)?;
                if let Some(ref mut actual_fast_forward_info) = fast_forward_info {
                    actual_fast_forward_info.write_chunk(&mut writer)?;
                } else {
                    warn!("Chunk {} is empty.", self.chunk_index)
                };
            } else {
                // Reading from SAM/BAM/CRAM and translating to FASTQ
                let mut writer =
                    get_fastq_writer(self.output.clone(), self.compression, self.threads)?;
                // Write the chunk
                let mut fast_forward_info =
                    reader.fast_forward(split_index, self.chunk_index, self.num_chunks)?;
                if let Some(ref mut actual_fast_forward_info) = fast_forward_info {
                    actual_fast_forward_info.translate_and_write_chunk(&mut writer)?;
                } else {
                    warn!("Chunk {} is empty.", self.chunk_index)
                };
            }
        } else {
            // reading from FASTQ
            let mut reader = get_fastq_reader(self.input.clone(), self.threads)?;
            let mut fast_forward_info =
                reader.fast_forward(split_index, self.chunk_index, self.num_chunks)?;

            if output_record_type == RecordType::Fastq {
                // reading from FASTQ and writing to FASTQ
                let mut writer =
                    get_fastq_writer(self.output.clone(), self.compression, self.threads)?;
                // Write the chunk
                if let Some(ref mut actual_fast_forward_info) = fast_forward_info {
                    actual_fast_forward_info.write_chunk(&mut writer)?;
                } else {
                    warn!("Chunk {} is empty.", self.chunk_index)
                };
            } else {
                // Reading from FASTQ and translating to SAM/BAM/CRAM
                // Should only be able to get here if output_format is specified;
                let default_format = self
                    .output_format
                    .clone()
                    .ok_or_else(|| anyhow!("Unspecified output format, should be unreachable."))?;
                // TODO: set minimal header, maybe allow sample ID, set query-group order, or similar?
                let writer_spec = SamWriterSpec::new(self.output.clone())
                    .format_from_path_or_default(default_format)?
                    .threads(self.threads)
                    .reference_fasta(self.ref_fasta.clone())
                    .compression(self.compression)
                    .to_owned();
                let mut writer = writer_spec.get_bam_writer()?;
                // Write the chunk
                if let Some(ref mut actual_fast_forward_info) = fast_forward_info {
                    actual_fast_forward_info.translate_and_write_chunk(&mut writer)?;
                } else {
                    warn!("Chunk {} is empty.", self.chunk_index)
                };
            }
        }
        Ok(())
    }
}

/// Implement the Command trait for `GetChunk` struct.
impl Command for GetChunk {
    /// Execute the get-chunk command to extract a specific chunk from the input file.
    fn execute(&self) -> Result<()> {
        info!("Using {} thread(s)", self.threads);
        self.write_chunk()
    }
}

#[cfg(test)]
mod tests {
    use super::{GetChunk, get_bam_reader};
    use crate::{commands::index::Index, test_utils::random_bam::QueryType};
    use anyhow::Result;
    use clap::Parser;
    use rstest::rstest;
    use rust_htslib::{
        bam::{Header, Read as BamRead, Record as BamRecord},
        errors::Error as HtslibErr,
    };
    use std::{
        collections::HashSet,
        fmt::Debug,
        iter::zip,
        num::NonZero,
        panic,
        path::{Path, PathBuf},
    };
    use tempfile::TempDir;

    /// Load header and records from the original truth BAM
    fn load_truth_bam<P>(bam: P) -> Result<(Header, Vec<BamRecord>)>
    where
        P: AsRef<Path>,
    {
        assert!(bam.as_ref().exists());
        let mut reader = get_bam_reader(bam, None::<PathBuf>, 1usize.try_into()?)?;
        let header = Header::from_template(reader.header());
        let records: Vec<BamRecord> = reader
            .records()
            .collect::<Result<Vec<BamRecord>, HtslibErr>>()?;
        Ok((header, records))
    }

    /// Load header, records, and # records / chunk from chunked BAMs
    fn load_chunk_bams<P>(
        chunk_paths: Vec<P>,
        num_reads: usize,
    ) -> Result<(Vec<Header>, Vec<BamRecord>, Vec<usize>)>
    where
        P: AsRef<Path> + Sized,
    {
        let mut headers: Vec<Header> = Vec::with_capacity(chunk_paths.len());
        let mut records: Vec<BamRecord> = Vec::with_capacity(num_reads);
        let mut chunk_lengths: Vec<usize> = Vec::with_capacity(chunk_paths.len());
        let mut last_qname: Option<String> = None;
        for chunk_path in chunk_paths {
            let mut reader = get_bam_reader(chunk_path, None::<PathBuf>, 1usize.try_into()?)?;
            headers.push(Header::from_template(reader.header()));
            let mut chunk_records: Vec<BamRecord> = reader
                .records()
                .collect::<Result<Vec<BamRecord>, HtslibErr>>()?;
            let chunk_queries = get_chunk_queries(&chunk_records);
            if let Some(previous_qname) = last_qname.clone()
                && let Some(first_chunk_qname) = chunk_queries.first()
            {
                assert!(
                    previous_qname != *first_chunk_qname,
                    "Qname {previous_qname} is split between chunks"
                )
            }
            last_qname = if let Some(qname) = chunk_queries.last() {
                Some(qname.to_owned())
            } else {
                None
            };
            chunk_lengths.push(
                chunk_queries
                    .iter()
                    .to_owned()
                    .collect::<HashSet<_>>()
                    .len(),
            );
            records.append(&mut chunk_records);
        }
        Ok((headers, records, chunk_lengths))
    }

    /// Get query names from chunk
    fn get_chunk_queries(chunk_records: &Vec<BamRecord>) -> Vec<String> {
        chunk_records
            .into_iter()
            .map(|rec| String::from_utf8_lossy(rec.qname()).to_string())
            .collect()
    }

    fn write_test_chunks<P1, P2>(
        test_bam: P1,
        index_path: P2,
        num_chunks: usize,
    ) -> Result<Vec<PathBuf>>
    where
        P1: AsRef<Path>,
        P2: AsRef<Path>,
    {
        let mut chunk_bams: Vec<PathBuf> = Vec::with_capacity(num_chunks);
        for chunk in 0..num_chunks {
            let output = test_bam
                .as_ref()
                .with_extension(format!("chunk_{chunk}_{num_chunks}.bam"));
            let command = GetChunk {
                input: test_bam.as_ref().to_path_buf(),
                index: Some(index_path.as_ref().to_path_buf()),
                ref_fasta: None::<PathBuf>,
                output: output.clone(),
                output_format: Some("bam".to_string()),
                threads: NonZero::<usize>::new(1usize).unwrap(),
                chunk_index: chunk,
                num_chunks: NonZero::<usize>::new(num_chunks).unwrap(),
                compression: Some(0u32),
            };
            command.write_chunk()?;
            chunk_bams.push(output.into_boxed_path().into_path_buf());
        }
        Ok(chunk_bams)
    }

    fn assert_records_equal(test_record: &BamRecord, truth_record: &BamRecord) -> () {
        assert!(
            test_record.qname() == truth_record.qname(),
            "Test qname != truth qname ({:?} != {:?})",
            String::from_utf8_lossy(test_record.qname()),
            String::from_utf8_lossy(truth_record.qname())
        );
        assert!(
            test_record.seq_len() == truth_record.seq_len(),
            "Test sequence length != truth sequence length ({} != {})",
            test_record.seq_len(),
            truth_record.seq_len(),
        );
        let test_sequence = test_record.seq().as_bytes();
        let truth_sequence = truth_record.seq().as_bytes();
        assert!(
            test_sequence == truth_sequence,
            "Test sequence != truth sequence ({:?} != {:?})",
            String::from_utf8_lossy(&test_sequence),
            String::from_utf8_lossy(&truth_sequence),
        );

        let test_quals = test_record.qual();
        let truth_quals = truth_record.qual();
        assert!(
            test_quals == truth_quals,
            "Test quals != truth quals ({:?} != {:?})",
            test_quals,
            truth_quals
        );
    }

    fn assert_vecs_equal<T, F>(test_vec: &Vec<T>, truth_vec: &Vec<T>, check_values_equal: F) -> ()
    where
        T: PartialEq + Debug + std::panic::RefUnwindSafe,
        F: Fn(&T, &T) -> () + std::panic::RefUnwindSafe,
    {
        let idx: usize = 0;
        for (x1, x2) in zip(test_vec, truth_vec) {
            if let Err(err) = panic::catch_unwind(|| check_values_equal(x1, x2)) {
                panic!(
                    "First mismatch at index {idx} ({x1:?} != {x2:?}). Test len: {}, truth len: {}. Error: {err:?}",
                    test_vec.len(),
                    truth_vec.len()
                )
            }
        }
        assert!(
            test_vec.len() == truth_vec.len(),
            "Test len != truth len ({} !=  {})",
            test_vec.len(),
            truth_vec.len()
        );
    }

    struct TestCase {
        num_queries: usize,
        num_bins: usize,
        num_chunks: usize,
        label: &'static str,
    }

    impl TestCase {
        pub fn new(
            num_queries: usize,
            num_bins: usize,
            num_chunks: usize,
            label: &'static str,
        ) -> Self {
            TestCase {
                num_queries,
                num_bins,
                num_chunks,
                label,
            }
        }
    }

    #[rstest(query_type => [QueryType::Single, QueryType::Paired, QueryType::Grouped],
        test_case => [
            TestCase::new(100, 20, 5, "even_divisions"),
            TestCase::new(101, 23, 5, "uneven divisions"),
            TestCase::new(100, 100, 5, "1 read per bin"),
            TestCase::new(100, 1000, 5, "Too few reads"),
            TestCase::new(100, 20, 211, "Too many chunks"),
            TestCase::new(0, 20, 211, "No reads"),
        ],
        output => [None, Some("passthrough.bam")]
        )
    ]
    fn test_chunks_recapitulate_bam(
        query_type: QueryType,
        test_case: TestCase,
        output: Option<&str>,
    ) -> Result<()> {
        let temp_dir = TempDir::new()?;
        let temp_path: PathBuf = temp_dir.path().to_path_buf();
        let output_path: Option<PathBuf> = if let Some(output_filename) = output {
            Some(temp_path.join(output_filename))
        } else {
            None
        };
        let (random_bam, num_reads) = query_type.random_bam(&temp_path, test_case.num_queries)?;

        let num_bins_str = test_case.num_bins.to_string();
        let mut args = vec![
            "index",
            "--input",
            random_bam.to_str().unwrap(),
            "--num-bins",
            num_bins_str.as_str(),
        ];
        if let Some(ref output_path) = output_path {
            args.push("--output");
            args.push(output_path.to_str().unwrap());
        };
        let index_tool = Index::try_parse_from(args)?;
        let index = index_tool.index_reads()?;

        let chunk_source = if let Some(output_path) = output_path {
            output_path
        } else {
            random_bam.clone()
        };
        assert!(chunk_source.is_file());
        let chunk_bams = write_test_chunks(random_bam.clone(), index, test_case.num_chunks)?;

        let (truth_header, truth_records) = load_truth_bam(random_bam)?;
        assert!(
            truth_records.len() == num_reads,
            "{}: Did not read the correct number of test reads (expected {num_reads}, got {})",
            test_case.label,
            truth_records.len()
        );
        let (chunk_headers, chunk_records, chunk_lengths) = load_chunk_bams(chunk_bams, num_reads)?;
        for chunk_header in chunk_headers {
            assert!(chunk_header.to_bytes() == truth_header.to_bytes());
        }
        assert_vecs_equal(&chunk_records, &truth_records, assert_records_equal);
        let expected_low_len = test_case.num_queries / test_case.num_chunks;
        let expected_high_len = expected_low_len + 1;
        for (idx, &chunk_length) in chunk_lengths.iter().enumerate() {
            assert!(
                expected_low_len <= chunk_length && chunk_length <= expected_high_len,
                "Bad chunk length for chunk {idx}: got {chunk_length} but expected [{expected_low_len}, {expected_high_len}]"
            );
        }
        Ok(())
    }
}
