use anyhow::{Result, anyhow};
use log::info;
use rust_htslib::bam::{
    Read as BamRead, Reader as BamReader, Record as BamRecord, Writer as BamWriter,
};
use seq_io::fastq::{
    OwnedRecord as OwnedSeqIoFastqRecord, Position, Reader as SeqIoFastqReader,
    Record as SeqIoFastqRecord,
};
use std::io::{BufRead, Read, Seek, SeekFrom, Write};
use std::num::NonZero;

use crate::fastq::{FastqReader, FastqRecord, FastqWriter};
use crate::maybe_compressed_io::MaybeCompressedWriter;

/// A trait with required functions for records that can be extracte as part of a chunk
pub trait ChunkableRecord {
    fn qname(&self) -> &[u8];
    fn seq(&self) -> &[u8];
    fn qual(&self) -> &[u8];
    fn new() -> Self;
    fn set_fields(&mut self, qname: &[u8], seq: &[u8], qual: &[u8]);

    fn translate<CR: ChunkableRecord>(&mut self, chunkable_record: &CR) {
        self.set_fields(
            chunkable_record.qname(),
            chunkable_record.seq(),
            chunkable_record.qual(),
        );
    }
}

/// Struct that includes all the information in SplitRecord, but includes the counts at the
/// end of the previous bin as well.
#[derive(Debug)]
pub struct SplitRange {
    /// File offset at the first read in this bin
    pub offset: u64,
    /// Cumulative number of queries in the entire reads file at the end of the previous bin
    pub num_previous_queries: usize,
    /// Cumulative number of queries in the entire reads file at the end of this bin.
    pub num_end_queries: usize,
    /// Cumulative number of reads in the entire reads file at the end of the previous bin
    pub num_previous_reads: usize,
    /// Cumulative number of reads in the entire reads file at the end of this bin.
    pub num_end_reads: usize,
}

/// A trait that allows fast-forwarding a chunkable reader. Given a chunk index and number of
/// chunks, get an index struct that yields an offset into the underlying file and reads and queries
/// from index bins.
pub trait FastForwardIndex {
    fn get_chunk_query_start(
        &self,
        chunk_index: usize,
        num_chunks: NonZero<usize>,
    ) -> Result<usize>;
    fn get_record_for_num_queries(&self, num_queries: usize) -> Option<SplitRange>;
}

/// Struct holding information needed to fast-forward a reader to a chunk and write it out
#[derive(Debug)]
pub struct FastForwardInfo<'a, R: ChunkableRecord, Reader: ChunkableRecordReader<R>> {
    num_queries: usize,
    stop_num_queries: usize,
    num_reads: usize,
    hard_stop_num_reads: usize,
    record: R,
    reader: &'a mut Reader,
}

impl<'a, R, Reader> FastForwardInfo<'a, R, Reader>
where
    R: ChunkableRecord + 'a,
    Reader: ChunkableRecordReader<R>,
{
    /// Write a chunk to the writer, reading and writing the same record type
    pub fn write_chunk<Writer>(&mut self, writer: &mut Writer) -> Result<()>
    where
        Writer: ChunkableRecordWriter<R>,
    {
        let mut last_query_name = self.record.qname().to_owned();
        while self.num_queries < self.stop_num_queries {
            // have the 1st record of a new query here
            writer.write(&self.record)?;
            self.reader
                .read_no_missing(&mut self.record, &mut self.num_reads)?;
            while self.record.qname() == last_query_name {
                writer.write(&self.record)?;
                self.reader
                    .read_no_missing(&mut self.record, &mut self.num_reads)?;
            }
            self.num_queries += 1;
            last_query_name = self.record.qname().to_owned();
        }
        // write the last query, being careful to check we don't read past the end of the bin/file
        writer.write(&self.record)?;
        while self.num_reads < self.hard_stop_num_reads {
            self.reader
                .read_no_missing(&mut self.record, &mut self.num_reads)?;
            if self.record.qname() != last_query_name {
                break;
            }
            writer.write(&self.record)?;
        }
        Ok(())
    }

    /// Write a chunk to the writer, translating to a different record type
    pub fn translate_and_write_chunk<WriteRecord, Writer>(
        &mut self,
        writer: &mut Writer,
    ) -> Result<()>
    where
        Writer: ChunkableRecordWriter<WriteRecord>,
        WriteRecord: ChunkableRecord,
    {
        let mut last_query_name = self.record.qname().to_owned();
        let mut write_record = WriteRecord::new();
        while self.num_queries < self.stop_num_queries {
            // have the 1st record of a new query here
            write_record.translate(&self.record);
            writer.write(&write_record)?;
            self.reader
                .read_no_missing(&mut self.record, &mut self.num_reads)?;
            while self.record.qname() == last_query_name {
                write_record.translate(&self.record);
                writer.write(&write_record)?;
                self.reader
                    .read_no_missing(&mut self.record, &mut self.num_reads)?;
            }
            self.num_queries += 1;
            last_query_name = self.record.qname().to_owned();
        }
        // write the last query, being careful to check we don't read past the end of the bin/file
        write_record.translate(&self.record);
        writer.write(&write_record)?;
        while self.num_reads < self.hard_stop_num_reads {
            self.reader
                .read_no_missing(&mut self.record, &mut self.num_reads)?;
            if self.record.qname() != last_query_name {
                break;
            }
            write_record.translate(&self.record);
            writer.write(&write_record)?;
        }
        Ok(())
    }
}

/// Public trait for a reader that can fast-forward to a desired chunk then read only the records
/// from that chunk. Directly tied to the type of ChunkableRecord.
pub trait ChunkableRecordReader<R>
where
    R: ChunkableRecord,
    Self: Sized,
{
    fn tell(&mut self) -> Result<u64>;
    fn seek(&mut self, offset: u64) -> Result<()>;
    // Read into existing record, returning potentially missing record, or Result with anyhow error
    fn read_into(&mut self, record: &mut R) -> Option<Result<()>>;

    /// Read into record that should not be missing, and handle any errors.
    fn read_no_missing(&mut self, record: &mut R, num_reads: &mut usize) -> Result<()> {
        *num_reads += 1;
        self.read_into(record)
            .unwrap_or_else(|| Err(anyhow!("file truncated.")))
            .map_err(|err| anyhow!("Unable to read at record {num_reads}: {err:?}"))
    }

    /// Fast forward the reader to the beginning of the chunk that needs to be read
    /// This may involve reading the first record of that chunk, in which case return it.
    fn fast_forward<'a, SI>(
        &'a mut self,
        split_index: SI,
        chunk_index: usize,
        num_chunks: NonZero<usize>,
    ) -> Result<Option<FastForwardInfo<'a, R, Self>>>
    where
        SI: FastForwardIndex,
    {
        // Number of completed queries that should have been read before this chunk starts
        let mut start_num_queries: usize =
            split_index.get_chunk_query_start(chunk_index, num_chunks)?;
        // Number of completed queries that should have been read by the end of this chunk
        let stop_num_queries: usize =
            split_index.get_chunk_query_start(chunk_index + 1, num_chunks)?;
        if start_num_queries >= stop_num_queries {
            // This will be an empty chunk
            return Ok(None);
        }
        // Get the SplitRange for the bin containing the requested start_num_queries
        let split_range = split_index
            .get_record_for_num_queries(start_num_queries)
            .ok_or_else(|| {
                anyhow!("Requested {start_num_queries} reads is past the end of the index.")
            })?;

        // seek to the file offset
        info!("Seeking to {}", split_range.offset);
        self.seek(split_range.offset)?;
        // if necessary, read until we reach the requested number of queries
        let mut num_reads: usize = split_range.num_previous_reads;
        let mut record = R::new();
        if start_num_queries > split_range.num_previous_queries {
            // Skip records until we *complete* the requested number of query groups.
            // The only way to know this is to *start* the query group AFTER start_num_queries
            let mut num_queries: usize = split_range.num_previous_queries;
            self.read_no_missing(&mut record, &mut num_reads)?;
            let mut last_query_name = record.qname().to_owned();
            num_queries += 1;
            while num_queries <= start_num_queries {
                self.read_no_missing(&mut record, &mut num_reads)?;
                let query_name = record.qname();
                if query_name != last_query_name {
                    num_queries += 1;
                    last_query_name = query_name.to_owned();
                }
            }
            start_num_queries = num_queries;
        } else {
            // Always read the first record for algorithm simplicity
            self.read_no_missing(&mut record, &mut num_reads)?;
            start_num_queries += 1; // this will be the start of a new query, because it's a new bin.
        };
        // We know that bins never split query groups, so set a limit on reads to avoid reading past
        // the end of the bin (or the file!) on the last query group
        let hard_stop_num_reads: usize = split_index
            .get_record_for_num_queries(stop_num_queries)
            .ok_or_else(|| anyhow!("Requested {stop_num_queries} past end of file"))?
            .num_end_reads;

        Ok(Some(FastForwardInfo {
            num_queries: start_num_queries,
            stop_num_queries,
            num_reads,
            hard_stop_num_reads,
            record,
            reader: self,
        }))
    }
}

/// Public trait for a writer that can write records from a chunk. Directly tied to the record type.
pub trait ChunkableRecordWriter<R>
where
    R: ChunkableRecord,
{
    fn write(&mut self, record: &R) -> Result<()>;
}

/// Implement ChunkableRecord trait for BAM/SAM/CRAM records.
impl ChunkableRecord for BamRecord {
    fn qname(&self) -> &[u8] {
        self.qname()
    }

    fn qual(&self) -> &[u8] {
        self.qual()
    }

    fn seq(&self) -> &[u8] {
        self.seq().encoded
    }

    fn new() -> BamRecord {
        BamRecord::new()
    }

    fn set_fields(&mut self, qname: &[u8], seq: &[u8], qual: &[u8]) {
        self.set(qname, None, seq, qual)
    }
}

/// Implement ChunkableRecordReader trait for BAM/SAM/CRAM readers.
impl ChunkableRecordReader<BamRecord> for BamReader {
    fn tell(&mut self) -> Result<u64> {
        Ok(<BamReader as BamRead>::tell(self) as u64)
    }
    fn seek(&mut self, offset: u64) -> Result<()> {
        Ok(<BamReader as BamRead>::seek(self, offset as i64)?)
    }

    fn read_into(&mut self, record: &mut BamRecord) -> Option<Result<()>> {
        match self.read(record) {
            Some(Err(err)) => Some(Err(anyhow!("{err}"))),
            Some(Ok(())) => Some(Ok(())),
            None => None,
        }
    }
}

/// Implement ChunkableRecordWriter trait for BAM/SAM/CRAM writers.
impl ChunkableRecordWriter<BamRecord> for BamWriter {
    fn write(&mut self, record: &BamRecord) -> Result<()> {
        Ok(self.write(record)?)
    }
}

/// Implement ChunkableRecord trait for seq_io FASTQ records.
impl ChunkableRecord for OwnedSeqIoFastqRecord {
    fn new() -> OwnedSeqIoFastqRecord {
        OwnedSeqIoFastqRecord {
            head: Vec::<u8>::new(),
            seq: Vec::<u8>::new(),
            qual: Vec::<u8>::new(),
        }
    }
    fn qname(&self) -> &[u8] {
        self.id_bytes()
    }
    fn qual(&self) -> &[u8] {
        &self.qual
    }
    fn seq(&self) -> &[u8] {
        &self.seq
    }

    fn set_fields(&mut self, qname: &[u8], seq: &[u8], qual: &[u8]) {
        self.head = qname.to_vec();
        self.seq = seq.to_vec();
        self.qual = qual.to_vec();
    }
}

/// Implement ChunkableRecordReader trait for seq_io FASTQ readers.
impl<R> ChunkableRecordReader<OwnedSeqIoFastqRecord> for SeqIoFastqReader<R>
where
    R: Read + Seek,
{
    fn tell(&mut self) -> Result<u64> {
        Ok(self.position().byte())
    }
    fn seek(&mut self, offset: u64) -> Result<()> {
        let pos: Position = Position::new(0, offset);
        //SeqIoFastqReader::<R, P>::seek(self, &pos)
        Ok(self.seek(&pos)?)
    }

    fn read_into(&mut self, record: &mut OwnedSeqIoFastqRecord) -> Option<Result<()>> {
        match self.next() {
            None => None,
            Some(Err(err)) => Some(Err(anyhow!("{err}"))),
            Some(Ok(ref_record)) => {
                ref_record.to_owned_record().clone_into(record);
                Some(Ok(()))
            }
        }
    }
}

/// Implement ChunkableRecordWriter trait for compressed FASTQ writers.
impl ChunkableRecordWriter<OwnedSeqIoFastqRecord> for MaybeCompressedWriter {
    fn write(&mut self, record: &OwnedSeqIoFastqRecord) -> Result<()> {
        Ok(record.write(self)?)
    }
}

/// Implement ChunkableRecord trait for custom FASTQ records.
impl ChunkableRecord for FastqRecord {
    fn new() -> Self {
        FastqRecord::new()
    }

    fn qname(&self) -> &[u8] {
        &self.name
    }

    fn qual(&self) -> &[u8] {
        &self.qualities
    }

    fn seq(&self) -> &[u8] {
        &self.sequence
    }

    fn set_fields(&mut self, qname: &[u8], seq: &[u8], qual: &[u8]) {
        self.name = qname.to_vec();
        self.sequence = seq.to_vec();
        self.qualities = qual.to_vec();
    }
}

/// Implement ChunkableRecordReader trait for custom FASTQ readers.
impl<R: BufRead + Seek> ChunkableRecordReader<FastqRecord> for FastqReader<R> {
    fn tell(&mut self) -> Result<u64> {
        Ok(self.stream_position()?)
        // let offset = self.stream_position()?;
        // info!("offset: {offset}");
        // Ok(offset)
    }
    fn seek(&mut self, offset: u64) -> Result<()> {
        if let Err(err) = <FastqReader<R> as Seek>::seek(self, SeekFrom::Start(offset)) {
            Err(anyhow!("{err}"))
        } else {
            Ok(())
        }
    }

    fn read_into(&mut self, record: &mut FastqRecord) -> Option<Result<()>> {
        match self.next() {
            None => None,
            Some(Err(err)) => Some(Err(anyhow!("{err}"))),
            Some(Ok(fastq_record)) => {
                *record = fastq_record;
                Some(Ok(()))
            }
        }
    }
}

/// Implement ChunkableRecordWriter trait for custom FASTQ writers.
impl<W: Write> ChunkableRecordWriter<FastqRecord> for FastqWriter<W> {
    fn write(&mut self, record: &FastqRecord) -> Result<()> {
        FastqWriter::<W>::write(self, record)
    }
}
