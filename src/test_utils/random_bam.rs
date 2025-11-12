use anyhow::Result;
use bam_builder::{BamBuilder, bam_order::BamSortOrder};
use rand::Rng;
use std::{
    fs,
    ops::Range,
    path::{Path, PathBuf},
};

/// Enum for generating random BAMs with different query types
pub enum QueryType {
    Single,  // single-end, e.g. PacBio, Nanopore
    Paired,  // pair-end, e.g. Illumina
    Grouped, // query-grouped, could be anything
}

impl QueryType {
    /// Generate a random BAM file with the specified number of query groups.
    ///
    /// Creates a BAM file in the given temporary directory with the appropriate query type
    /// (single-end, paired-end, or grouped). The returned tuple contains the path to the
    /// generated BAM file and the total number of reads in it.
    ///
    /// # Arguments
    /// * `temp_path` - Directory where the BAM file will be created
    /// * `num_queries` - Number of query groups to generate
    ///
    /// # Returns
    /// A tuple containing the path to the generated BAM file and the total number of reads
    pub fn random_bam<P>(&self, temp_path: &P, num_queries: usize) -> Result<(PathBuf, usize)>
    where
        P: AsRef<Path>,
    {
        let mut temp_file = temp_path.as_ref().to_path_buf();
        temp_file.push(format!("random-{}-{num_queries}.bam", self.label()));
        if temp_file.exists() {
            let len = fs::metadata(temp_file.clone())?.len();
            assert!(len == 0, "{:?} has size {}", temp_file.clone(), len);
        }
        match self {
            QueryType::Single => Self::random_single_bam(temp_file, num_queries),
            QueryType::Paired => Self::random_paired_bam(temp_file, num_queries),
            QueryType::Grouped => Self::random_grouped_bam(temp_file, num_queries),
        }
    }

    /// Get a descriptive label for this query type.
    ///
    /// # Returns
    /// A string slice containing "single", "paired", or "grouped" depending on the query type.
    pub fn label(&self) -> &'static str {
        match self {
            QueryType::Single => "single",
            QueryType::Paired => "paired",
            QueryType::Grouped => "grouped",
        }
    }

    /// Generate a random BAM file with single-end reads.
    ///
    /// Creates a BAM file with one read per query group.
    ///
    /// # Arguments
    /// * `temp_file` - Path where the BAM file will be written
    /// * `num_queries` - Number of query groups (and reads) to generate
    ///
    /// # Returns
    /// A tuple containing the path to the generated BAM file and the number of reads
    fn random_single_bam<P>(temp_file: P, num_queries: usize) -> Result<(PathBuf, usize)>
    where
        P: AsRef<Path>,
    {
        let mut builder = BamBuilder::new(
            150,                             // default read length
            30,                              // default base quality
            format!("single-{num_queries}"), // name of sample
            None,                            // optional read group id
            BamSortOrder::Unsorted,          // how to sort reads when `.sort` is called
            None,                            // optional sequence dictionary
            None,                            // optional seed used for generating random bases
        );
        for idx in 0..num_queries {
            let read = builder
                .frag_builder()
                .name(format!("Single{idx:06}"))
                .build()?;
            builder.add_frag(read);
        }
        builder.to_path(temp_file.as_ref())?;
        assert!(temp_file.as_ref().exists());
        Ok((temp_file.as_ref().to_owned(), num_queries))
    }

    /// Generate a random BAM file with paired-end reads.
    ///
    /// Creates a BAM file with two reads (forward and reverse) per query group.
    ///
    /// # Arguments
    /// * `temp_file` - Path where the BAM file will be written
    /// * `num_queries` - Number of query groups to generate
    ///
    /// # Returns
    /// A tuple containing the path to the generated BAM file and the total number of reads (2 * num_queries)
    fn random_paired_bam<P>(temp_file: P, num_queries: usize) -> Result<(PathBuf, usize)>
    where
        P: AsRef<Path>,
    {
        let mut builder = BamBuilder::new(
            150,                             // default read length
            30,                              // default base quality
            format!("paired-{num_queries}"), // name of sample
            None,                            // optional read group id
            BamSortOrder::Unsorted,          // how to sort reads when `.sort` is called
            None,                            // optional sequence dictionary
            None,                            // optional seed used for generating random bases
        );
        for idx in 0..num_queries {
            let pair = builder
                .pair_builder()
                .name(format!("Pair{idx:06}"))
                .build()?;
            builder.add_pair(pair);
        }
        builder.to_path(temp_file.as_ref())?;
        assert!(temp_file.as_ref().exists());
        Ok((temp_file.as_ref().to_owned(), 2 * num_queries))
    }

    /// Generate a random BAM file with query-grouped reads.
    ///
    /// Creates a BAM file with a variable number of reads per query group.
    ///
    /// # Arguments
    /// * `temp_file` - Path where the BAM file will be written
    /// * `num_queries` - Number of query groups to generate
    ///
    /// # Returns
    /// A tuple containing the path to the generated BAM file and the total number of reads
    fn random_grouped_bam<P>(temp_file: P, num_queries: usize) -> Result<(PathBuf, usize)>
    where
        P: AsRef<Path>,
    {
        let mut builder = BamBuilder::new(
            150,                              // default read length
            30,                               // default base quality
            format!("grouped-{num_queries}"), // name of sample
            None,                             // optional read group id
            BamSortOrder::Unsorted,           // how to sort reads when `.sort` is called
            None,                             // optional sequence dictionary
            None,                             // optional seed used for generating random bases
        );
        let group_size_range: Range<usize> = Range::<usize> {
            start: 1usize,
            end: 5usize,
        };
        let mut num_reads: usize = 0;
        for group in 0..num_queries {
            let group_size = builder.rng.random_range(group_size_range.clone());
            num_reads += group_size;
            for _ in 0..group_size {
                let name = format!("Group{group:06}");
                let read = builder.frag_builder().name(name.clone()).build()?;
                builder.add_frag(read);
            }
        }
        builder.to_path(temp_file.as_ref())?;
        assert!(temp_file.as_ref().exists());
        Ok((temp_file.as_ref().to_owned(), num_reads))
    }
}
