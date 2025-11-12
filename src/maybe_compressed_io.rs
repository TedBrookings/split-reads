use crate::seekable_chain::Chain;
use anyhow::{Result, anyhow};
use noodles_bgzf::{
    VirtualPosition,
    io::{MultithreadedReader, MultithreadedWriter, Seek as NoodlesSeek},
};
use std::{
    fs::{File, OpenOptions, create_dir_all},
    io::{BufRead, BufReader, BufWriter, Cursor, Read, Seek, SeekFrom, Write},
    num::NonZero,
    path::Path,
};

/// First bytes of gzipped file
const BGZIP_MAGIC_NUMBER: [u8; 2] = [0x1fu8, 0x8bu8];

/// Helper function to get a File object that can be read from or written to, given the supplied
/// path. The path may be "-", in which case we will read from stdin or write to stdout
pub fn open_file<P: AsRef<Path>>(path: P, for_writing: bool) -> Result<File> {
    if path.as_ref().to_str() == Some("-") {
        let default = if for_writing {
            "/dev/stdout"
        } else {
            "/dev/stdin"
        };
        OpenOptions::new()
            .write(for_writing)
            .read(!for_writing)
            .open(default)
            .map_err(|err| anyhow!("Opening {default}: {err}"))
    } else {
        if for_writing && let Some(parent_dir) = path.as_ref().parent() {
            create_dir_all(parent_dir)?
        }
        OpenOptions::new()
            .write(for_writing)
            .read(!for_writing)
            .create(for_writing)
            .open(path.as_ref())
            .map_err(|err| {
                let fq = path.as_ref();
                anyhow!("Opening {fq:?}: {err}")
            })
    }
}

/// Type alias for the ChainReader that is used by Compressed or Uncompressed readers.
type Inner = Chain<Cursor<Vec<u8>>, File>;

/// Enum for a file that may or may not be compressed.
pub enum MaybeCompressedReader {
    Compressed(MultithreadedReader<Inner>),
    Uncompressed(BufReader<Inner>),
}

impl MaybeCompressedReader {
    /// Open a possibly compressed input path. input_path can be set to "-" to read from stdin.
    /// Get a buffered reader that can get read the plaintext.
    pub fn new<P: AsRef<Path>>(
        input_path: P,
        decompression_threads: NonZero<usize>,
    ) -> Result<MaybeCompressedReader> {
        let mut input_file = open_file(input_path, false)?;
        let mut first_bytes = [0u8; 2];
        input_file.read_exact(&mut first_bytes)?;
        let mut first_bytes_cursor = Cursor::new(first_bytes.into());
        first_bytes_cursor.seek(SeekFrom::Start(0))?;
        let chain: Inner = Chain::new(first_bytes_cursor, input_file)?;
        if first_bytes == BGZIP_MAGIC_NUMBER {
            // it's gzipped, unzip with requested number of threads
            Ok(MaybeCompressedReader::Compressed(
                MultithreadedReader::with_worker_count(decompression_threads, chain),
            ))
        } else {
            // it's not gzipped, read plain text single-threaded
            Ok(MaybeCompressedReader::Uncompressed(BufReader::new(chain)))
        }
    }
}

/// impl Seek trait for MaybeCompressedReader
/// - Compressed readers use VirtualPosition for seeking,
/// - Uncompressed readers use normal offset
impl Seek for MaybeCompressedReader {
    fn seek(&mut self, pos: std::io::SeekFrom) -> std::io::Result<u64> {
        match self {
            Self::Compressed(reader) => match pos {
                SeekFrom::Start(start_pos) => {
                    let virtual_pos = VirtualPosition::from(start_pos);
                    if let Err(err) = reader.seek_to_virtual_position(virtual_pos) {
                        Err(err)
                    } else {
                        Ok(start_pos)
                    }
                }
                SeekFrom::Current(0) => Ok(reader.virtual_position().into()),
                _ => Err(std::io::Error::other("Cannot SeekFrom other than Start")),
            },
            Self::Uncompressed(reader) => reader.seek(pos),
        }
    }
}

/// impl Read trait for MaybeCompressedReader
impl Read for MaybeCompressedReader {
    fn read(&mut self, buf: &mut [u8]) -> Result<usize, std::io::Error> {
        match self {
            MaybeCompressedReader::Compressed(inner) => inner.read(buf),
            MaybeCompressedReader::Uncompressed(inner) => inner.read(buf),
        }
    }
}

/// impl BufRead for MaybeCompressedReader
impl BufRead for MaybeCompressedReader {
    fn fill_buf(&mut self) -> std::io::Result<&[u8]> {
        match self {
            MaybeCompressedReader::Compressed(inner) => inner.fill_buf(),
            MaybeCompressedReader::Uncompressed(inner) => inner.fill_buf(),
        }
    }

    fn consume(&mut self, amt: usize) {
        match self {
            MaybeCompressedReader::Compressed(inner) => inner.consume(amt),
            MaybeCompressedReader::Uncompressed(inner) => inner.consume(amt),
        }
    }
}

/// Enum for writing a single fastq (regardless of read organization). Either compressed or not
pub enum MaybeCompressedWriter {
    Compressed(BufWriter<MultithreadedWriter<File>>),
    Uncompressed(BufWriter<File>),
}

impl MaybeCompressedWriter {
    /// Create new writer. Compression is primarily determined by the output path (compressed if it
    /// ends in ".gz" or ".bgz", uncompressed otherwise), but when writing to stdout, the bool
    /// `compressed` will determine if output is compressed or not. Threads will only be used for
    /// compressed output.
    pub fn new<P: AsRef<Path>>(
        input_path: P,
        compressed: bool,
        threads: NonZero<usize>,
    ) -> Result<MaybeCompressedWriter> {
        let fastq_file = open_file(input_path.as_ref(), true)?;
        if MaybeCompressedWriter::is_compressed(input_path, compressed) {
            Ok(MaybeCompressedWriter::Compressed(BufWriter::new(
                MultithreadedWriter::with_worker_count(threads, fastq_file),
            )))
        } else {
            Ok(MaybeCompressedWriter::Uncompressed(BufWriter::new(
                fastq_file,
            )))
        }
    }

    /// Determine if output is compressed. When writing to a real path, make compressed if the path
    /// ends in ".gz" or ".bgz", uncompressed otherwise. When writing to stdout, obey `compressed`
    /// boolean.
    fn is_compressed<P: AsRef<Path>>(input_path: P, compressed: bool) -> bool {
        match input_path.as_ref().extension() {
            Some(os_str) => (os_str == "gz") || (os_str == "bgz"),
            None => compressed,
        }
    }
}

/// impl Write trait for MaybeCompressedWriter
impl Write for MaybeCompressedWriter {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        match self {
            MaybeCompressedWriter::Compressed(inner) => inner.write(buf),
            MaybeCompressedWriter::Uncompressed(inner) => inner.write(buf),
        }
    }

    fn flush(&mut self) -> std::io::Result<()> {
        match self {
            MaybeCompressedWriter::Compressed(inner) => inner.flush(),
            MaybeCompressedWriter::Uncompressed(inner) => inner.flush(),
        }
    }
}
