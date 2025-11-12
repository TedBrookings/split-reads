use crate::{
    chunkable::{
        ChunkableRecord, ChunkableRecordReader, ChunkableRecordWriter, FastForwardIndex, SplitRange,
    },
    path_type::PathType,
};
use anyhow::{Result, anyhow};
use bisection::bisect_left_by;
use log::{debug, info, warn};
use rust_htslib::bgzf::{Reader as BgzfReader, Writer as BgzfWriter};
use serde::{Deserialize, Serialize};
use std::{
    cmp::max,
    io::{Read, Write},
    num::NonZero,
    ops::RangeBounds,
    path::Path,
    time::{Duration, SystemTime},
    vec::Vec,
};

/// Version string for SplitIndex header.
const VERSION: &str = "1.0";

/// Default extension for split index files.
pub const SPLIT_INDEX_EXTENSION: &str = "si";

/// Drain range of bytes from the front of passed Vec, and return it as a new Vec
fn split_off<R>(bytes: &mut Vec<u8>, range: R) -> Result<Vec<u8>>
where
    R: RangeBounds<usize>,
{
    if range.contains(&bytes.len()) {
        Err(anyhow!(
            "Requested range extends past end of bytes. Index record truncated."
        ))
    } else {
        Ok(bytes.drain(range).collect())
    }
}

/// Deserialize a usize from the bytes buffer, and shorten the buffer
fn deserialize_usize(bytes: &mut Vec<u8>) -> Result<usize> {
    let usize_bytes = split_off(bytes, ..size_of::<usize>())?;
    Ok(usize::from_le_bytes(usize_bytes.as_slice().try_into()?))
}

/// Deserialize a u64 from the bytes buffer, and shorten the buffer
fn deserialize_u64(bytes: &mut Vec<u8>) -> Result<u64> {
    let u64_bytes = split_off(bytes, ..size_of::<u64>())?;
    Ok(u64::from_le_bytes(u64_bytes.as_slice().try_into()?))
}

/// Struct for holding records in the SplitIndex. It represents a very small bin in the original
/// reads file.
#[derive(Clone, Copy, Serialize, Deserialize, PartialEq)]
struct SplitRecord {
    /// File offset at the first read in the bin
    pub offset: u64,
    /// Cumulative number of queries in the entire reads file at the end of the bin.
    pub num_queries: usize,
    /// Cumulative number of reads in the entire reads file at the end of the bin.
    pub num_reads: usize,
}

impl SplitRecord {
    /// Serialize by appending to bytes
    pub fn serialize(&self, bytes: &mut Vec<u8>) {
        bytes.extend(self.offset.to_le_bytes());
        bytes.extend(self.num_queries.to_le_bytes());
        bytes.extend(self.num_reads.to_le_bytes());
    }

    /// Deserialize by draining from bytes
    pub fn deserialize(bytes: &mut Vec<u8>) -> Result<Self> {
        Ok(SplitRecord {
            offset: deserialize_u64(bytes)?,
            num_queries: deserialize_usize(bytes)?,
            num_reads: deserialize_usize(bytes)?,
        })
    }
}

/// Struct for holding and manipulating all the SplitRecords for a reads file.
#[derive(Serialize, Deserialize, Clone, PartialEq)]
pub struct SplitIndex {
    split_records: Vec<SplitRecord>,
}

impl SplitIndex {
    /// Create a new empty SplitIndex with specified capacity.
    pub fn with_capacity(num_records: usize) -> Self {
        SplitIndex {
            split_records: Vec::with_capacity(num_records),
        }
    }

    /// Get the length of the index
    pub fn len(&self) -> usize {
        self.split_records.len()
    }

    /// Return true if the index is empty
    pub fn is_empty(&self) -> bool {
        self.split_records.is_empty()
    }

    /// Get the total number of indexed queries. For a complete SplitIndex, this is the number of
    /// queries in the reads file.
    pub fn num_queries(&self) -> usize {
        if let Some(split_record) = self.split_records.last() {
            split_record.num_queries
        } else {
            0
        }
    }

    /// Get the total number of indexed reads. For a complete SplitIndex, this is the number of
    /// reads in the reads file.
    pub fn num_reads(&self) -> usize {
        if let Some(split_record) = self.split_records.last() {
            split_record.num_reads
        } else {
            0
        }
    }

    /// Add a new SplitRecord to the Index
    fn add_record(&mut self, split_record: SplitRecord) {
        self.split_records.push(split_record);
    }

    /// Return a SplitRecord for the next bin
    fn start_next_record(&self, offset: u64) -> SplitRecord {
        SplitRecord {
            offset,
            num_queries: self.num_queries() + 1,
            num_reads: self.num_reads() + 1,
        }
    }

    /// Given the index of a bin, return the corresponding BinRange. Return None if past the end.
    fn index_to_bin_range(&self, index: usize) -> Option<SplitRange> {
        if let Some(split_record) = self.split_records.get(index) {
            if index == 0 {
                Some(SplitRange {
                    offset: split_record.offset,
                    num_previous_queries: 0,
                    num_end_queries: split_record.num_queries,
                    num_previous_reads: 0,
                    num_end_reads: split_record.num_reads,
                })
            } else {
                let previous_record = self.split_records.get(index - 1)?;
                Some(SplitRange {
                    offset: split_record.offset,
                    num_previous_queries: previous_record.num_queries,
                    num_end_queries: split_record.num_queries,
                    num_previous_reads: previous_record.num_reads,
                    num_end_reads: split_record.num_reads,
                })
            }
        } else {
            warn!("Requested index {index} from {} split records.", self.len());
            None
        }
    }

    /// Serialize SplitIndex to bytes.
    pub fn serialize(self) -> Vec<u8> {
        let mut bytes: Vec<u8> = format!("split-index {VERSION}\n").as_bytes().to_vec();
        bytes.extend(&self.len().to_le_bytes());
        for split_record in self.split_records {
            split_record.serialize(&mut bytes);
        }
        bytes
    }

    /// Write SplitIndex to the requested path.
    pub fn write<P>(self, path: P) -> Result<usize>
    where
        P: AsRef<Path>,
    {
        let mut writer = match PathType::from_path(path)? {
            PathType::Pipe => Ok(BgzfWriter::from_stdout()?),
            PathType::FilePath(file_path) => Ok(BgzfWriter::from_path(file_path)?),
            PathType::UrlPath(_) => Err(anyhow!("Cannot write directly to a cloud URL")),
        }?;
        writer
            .write(&self.serialize())
            .map_err(|err| anyhow!("{err}"))
    }

    /// Build the SplitIndex. Never split query groups. Because the total number of records and
    /// query groups is unknown, dynamically space bins as
    /// max(1, running_total_queries / requested_final_number_of_bins)
    /// The number of actual bins grows logarithmically in the limit of large numbers of query
    /// groups. Later on the bins are interpolated down to the requested amount.
    pub fn build<Record, Reader, Writer>(
        mut reader: Reader,
        mut writer: Option<Writer>,
        num_bins: NonZero<usize>,
        update_interval: u64,
    ) -> Result<SplitIndex>
    where
        Record: ChunkableRecord,
        Reader: ChunkableRecordReader<Record>,
        Writer: ChunkableRecordWriter<Record>,
    {
        let mut record = Record::new();
        let mut split_index = SplitIndex::with_capacity(num_bins.into());
        let mut next_query_bin: usize = 1;
        // In this and following calculation of offset, if there is a writer, it we should invoke
        // writer.tell(). However
        // 1. rust_htslib currently does not provide writer.tell
        // 2. if the modality is identical, reader.tell() and writer.tell() yield the same values.
        let mut offset: u64 = reader.tell()?;
        let mut last_update = SystemTime::now();
        let update_duration = Duration::from_secs(update_interval);
        if let Some(result) = reader.read_into(&mut record) {
            result?;
            if let Some(ref mut actual_bam_writer) = writer {
                actual_bam_writer.write(&record)?;
            }
            let mut last_query_name: Vec<u8> = record.qname().to_vec();
            let mut split_record = split_index.start_next_record(offset);
            offset = reader.tell()?;
            while let Some(result) = reader.read_into(&mut record) {
                let now = SystemTime::now();
                if now.duration_since(last_update)? > update_duration {
                    info!(
                        "Indexed {} reads and {} queries.",
                        split_record.num_reads, split_record.num_queries
                    );
                    last_update = now;
                }
                result?;
                if let Some(ref mut actual_bam_writer) = writer {
                    actual_bam_writer.write(&record)?;
                }
                if record.qname() == last_query_name {
                    // inside a query group, do not update bin
                    split_record.num_reads += 1;
                } else if split_record.num_queries < next_query_bin {
                    // new query group, but not time to change the bin yet
                    last_query_name = record.qname().to_vec();
                    split_record.num_reads += 1;
                    split_record.num_queries += 1;
                } else {
                    // time for a new bin and query goal
                    last_query_name = record.qname().to_vec();
                    split_index.add_record(split_record);
                    next_query_bin += max(1usize, split_index.num_queries() / num_bins);
                    split_record = split_index.start_next_record(offset);
                }
                offset = reader.tell()?;
            }
            split_index.add_record(split_record);
        } else {
            warn!("Empty index: no reads");
        }
        Ok(split_index)
    }

    /// Downsize via interpolation to roughly evenly spaced bins of the requested size.
    pub fn downsize_reads(&self, num_bins: NonZero<usize>) -> Result<Self> {
        if usize::from(num_bins) > self.len() {
            // This is a normal thing that can happen when indexing a BAM with very few records,
            // not an error. Just return self
            warn!("Keeping original SplitIndex with fewer bins than requested.");
            return Ok(self.clone());
        }
        let mut downsized = SplitIndex::with_capacity(num_bins.into());
        // the last bin *must* be the same, because it contains the total number of reads and
        // queries. All others are taken as close as possible to evenly-spaced
        let mut last_offset = self
            .split_records
            .first()
            .ok_or_else(|| anyhow!("No bins in original index. Should be unreachable."))?
            .offset;
        let mut last_index: Option<usize> = None;
        for bin in 1..num_bins.into() {
            let target_num_queries: usize = self.get_chunk_query_start(bin, num_bins)?;
            let mut index: usize = bisect_left_by(&self.split_records, |&record| {
                record.num_queries.cmp(&target_num_queries)
            });
            if index > 0 &&
                // index points to the first SplitRecord with num_reads > target_num_reads,
                // but the previous one might be *closer* to the target. Check
                target_num_queries - self.split_records[index - 1].num_queries
                    <= self.split_records[index].num_queries - target_num_queries
            {
                index -= 1;
            }
            if let Some(actual_last_index) = last_index
                && index <= actual_last_index
            {
                warn!("Original SplitIndex has few bins, so down-sizing is sparser than expected.")
            } else {
                let mut new_record = self.split_records[index];
                new_record.offset = last_offset;
                downsized.add_record(new_record);
                if index + 1 < self.len() {
                    last_offset = self.split_records[index + 1].offset;
                } else {
                    // we somehow reached the end of the index early. Warn and return what we have
                    warn!(
                        "Original SplitIndex has few bins, so down-sizing is sparser than expected."
                    );
                    return Ok(downsized);
                }
                last_index = Some(index)
            }
        }
        if let Some(last_split_record) = self.split_records.last() {
            let mut new_record = *last_split_record;
            new_record.offset = last_offset;
            downsized.add_record(new_record);
        }
        Ok(downsized)
    }

    /// Parse the header and extract the version string.
    fn check_header(bytes: &mut Vec<u8>) -> Result<String> {
        let pos = bytes
            .iter()
            .position(|c| *c == b'\n')
            .ok_or_else(|| anyhow!("Unable to parse header. Corrupted index or wrong file."))?;
        let mut header: Vec<u8> = bytes.drain(..=pos).collect();
        let expected_front = b"split-index ";
        if header.len() < expected_front.len() {
            Err(anyhow!(
                "Unable to parse header. Corrupted index or wrong file."
            ))?;
        }
        let front: Vec<u8> = header.drain(..expected_front.len()).collect();
        if front != expected_front {
            Err(anyhow!(
                "Unable to parse header. Corrupted index or wrong file."
            ))
        } else {
            // remainder of header should be version string and newline
            let mut version: String = String::from_utf8(header.to_owned())?;
            version.pop(); // remove newline
            Ok(version)
        }
    }

    /// Deserialize SplitIndex from bytes
    pub fn deserialize(bytes: &mut Vec<u8>) -> Result<Self> {
        let version = Self::check_header(bytes)?;
        // here we could use different loading routines in a hypothetical future with multiple
        // versions of the index. For now we just assert it's equal to the expected
        if version != VERSION {
            return Err(anyhow!("Unknown split-index version: {version}"));
        }
        let len: usize = deserialize_usize(bytes)?;
        debug!("Got {len} records in SplitIndex");
        let mut split_index = SplitIndex::with_capacity(len);
        for _ in 0..len {
            split_index.add_record(SplitRecord::deserialize(bytes)?);
        }
        Ok(split_index)
    }

    /// Read SplitIndex from the requested path or URL.
    pub fn read<P>(path: P) -> Result<Self>
    where
        P: AsRef<Path>,
    {
        let mut reader: BgzfReader = match PathType::from_path(path)? {
            PathType::Pipe => BgzfReader::from_stdin().map_err(|err| anyhow!("{err}")),
            PathType::FilePath(file_path) => Ok(BgzfReader::from_path(file_path)?),
            PathType::UrlPath(url) => Ok(BgzfReader::from_url(&url)?),
        }?;
        let mut buf: Vec<u8> = Vec::new();
        reader.read_to_end(&mut buf)?;
        Self::deserialize(&mut buf)
    }

    /// Only used in tests, but tested in index tool, so can't have cfg(test)
    /// get vec of the num_queries for each record
    pub fn get_split_record_num_queries(&self) -> Vec<usize> {
        self.split_records.iter().map(|sr| sr.num_queries).collect()
    }
}

impl FastForwardIndex for SplitIndex {
    /// Given a number of query groups, return the SplitRange for the bin containing that number.
    fn get_record_for_num_queries(&self, num_queries: usize) -> Option<SplitRange> {
        let index: usize = bisect_left_by(&self.split_records, |&record| {
            record.num_queries.cmp(&num_queries)
        });
        self.index_to_bin_range(index)
    }

    /// Given a chunk index and number of chunks, return the corresponding number of query groups
    /// that should have already been read before that chunk. It could also be viewed as the 0-based
    /// index of the query starting that chunk.
    fn get_chunk_query_start(
        &self,
        chunk_index: usize,
        num_chunks: NonZero<usize>,
    ) -> Result<usize> {
        let num_chunks: usize = num_chunks.into();
        if chunk_index <= num_chunks {
            // do chunk_index * self.num_reads() / num_chunks without rounding error or overflow
            let div_mod: (usize, usize) = (
                self.num_queries() / num_chunks,
                self.num_queries() % num_chunks,
            );
            let start = (chunk_index * div_mod.0) + ((chunk_index * div_mod.1) / num_chunks);
            Ok(start)
        } else {
            Err(anyhow!(
                "Invalid chunk index {chunk_index} for {num_chunks}"
            ))
        }
    }
}

#[cfg(test)]
mod tests {
    use anyhow::Result;
    use std::u64;
    use tempfile::NamedTempFile;

    use crate::split_index::{SplitIndex, SplitRecord};

    /// For testing serialization, etc. Create a random nonsensical SplitRecord.
    fn random_split_record<R>(rng: &mut R) -> SplitRecord
    where
        R: rand::Rng,
    {
        SplitRecord {
            offset: rng.random_range(u64::MIN..u64::MAX),
            num_queries: rng.random_range(0..usize::MAX),
            num_reads: rng.random_range(0..usize::MAX),
        }
    }

    /// For testing serialization, etc. Create a random nonsensical SplitIndex.
    fn random_split_index(num_bins: usize) -> SplitIndex {
        let mut rng = rand::rng();
        let mut split_index = SplitIndex::with_capacity(num_bins);
        for _ in 0..num_bins {
            split_index.add_record(random_split_record(&mut rng));
        }
        split_index
    }

    /// Test that serializing then deserializing recapitulate the original SplitIndex.
    #[test]
    fn test_serialize_round_trip() -> Result<()> {
        let split_index: SplitIndex = random_split_index(10000);
        let deserialized = SplitIndex::deserialize(&mut split_index.clone().serialize())?;
        assert!(deserialized == split_index);
        Ok(())
    }

    /// Test that writing then reading recapitulate the original SplitIndex.
    #[test]
    fn test_write_round_trip() -> Result<()> {
        let index_file = NamedTempFile::new().expect("Could not create temp file");
        let split_index: SplitIndex = random_split_index(10000);
        split_index.clone().write(index_file.path())?;
        let deserialized = SplitIndex::read(index_file.path())?;
        assert!(deserialized == split_index);
        Ok(())
    }
}
