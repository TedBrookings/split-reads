use crate::seekable_split::Split;
use anyhow::{Result, anyhow};
use std::io::{BufRead, Result as IoResult, Seek, Write};

/// Struct for holding fastq records
#[derive(Clone, Debug)]
pub struct FastqRecord {
    pub name: Vec<u8>,
    pub sequence: Vec<u8>,
    pub separator: Vec<u8>,
    pub qualities: Vec<u8>,
}

impl FastqRecord {
    /// Shortcut to get length of the read
    pub fn len(&self) -> usize {
        self.sequence.len()
    }
    /// Unused, should never be true, but keeps clippy happy
    pub fn is_empty(&self) -> bool {
        self.sequence.is_empty()
    }

    pub fn new() -> Self {
        Self {
            name: Vec::<u8>::new(),
            sequence: Vec::<u8>::new(),
            separator: Vec::<u8>::new(),
            qualities: Vec::<u8>::new(),
        }
    }
}

impl Default for FastqRecord {
    fn default() -> Self {
        Self::new()
    }
}

/// Struct for reading individual fastq files, using underlying `BufRead` object
pub struct FastqReader<R: BufRead> {
    split: Split<R>,
}

/// Implement remaining `FastqReader` functions for any `BufRead` underlying reader
impl<R: BufRead> FastqReader<R> {
    /// Create new `FastqReader` from base reader object
    pub fn new(reader: R) -> Self {
        FastqReader {
            split: Split::new(reader, b'\n'),
        }
    }

    /// While reading a record, handle possible missing / incomplete data
    fn unwrap_next(&mut self) -> Result<Vec<u8>> {
        match self.split.next() {
            None => Err(anyhow!("Incomplete fastq record")),
            Some(Ok(vec)) => Ok(vec),
            Some(Err(err)) => Err(anyhow!("{err}")),
        }
    }

    /// Get the next fastq record
    fn next_fastq_record(&mut self, name: Vec<u8>) -> Result<FastqRecord> {
        let sequence = self.unwrap_next()?;
        let separator = self.unwrap_next()?;
        let qualities = self.unwrap_next()?;
        Ok(FastqRecord {
            name,
            sequence,
            separator,
            qualities,
        })
    }
}

/// impl Seek for FastqReader, delegating to underlying Split
impl<R: BufRead + Seek> Seek for FastqReader<R> {
    fn seek(&mut self, pos: std::io::SeekFrom) -> IoResult<u64> {
        self.split.seek(pos)
    }
}

/// impl Iterator for `FastqIterator`: yield Result<FastqRecord>
impl<R: BufRead> Iterator for FastqReader<R> {
    type Item = Result<FastqRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.split.next() {
            None => None,
            Some(Err(err)) => Some(Err(anyhow!("{err}"))),
            Some(Ok(name)) => Some(self.next_fastq_record(name)),
        }
    }
}

/// Public struct for writing fastq records
pub struct FastqWriter<W: Write> {
    inner: W,
}

/// impl FastqWriter, just write out the four lines separated by newlines
impl<W: Write> FastqWriter<W> {
    const NEWLINE: [u8; 1] = [b'\n'];

    pub fn new(writer: W) -> Self {
        FastqWriter { inner: writer }
    }

    pub fn write(&mut self, fastq_record: &FastqRecord) -> Result<()> {
        self.inner.write_all(&fastq_record.name)?;
        self.inner.write_all(&Self::NEWLINE)?;

        self.inner.write_all(&fastq_record.sequence)?;
        self.inner.write_all(&Self::NEWLINE)?;

        self.inner.write_all(&fastq_record.separator)?;
        self.inner.write_all(&Self::NEWLINE)?;

        self.inner.write_all(&fastq_record.qualities)?;
        self.inner.write_all(&Self::NEWLINE)?;
        Ok(())
    }
}
