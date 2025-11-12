use crate::{
    fastq::{FastqReader, FastqWriter},
    maybe_compressed_io::{MaybeCompressedReader, MaybeCompressedWriter},
    path_type::PathType,
};
use anyhow::Result;
use env;
use log::warn;
use rust_htslib::bam::{Read, Reader};
use seq_io::fastq::Reader as SeqIoFastqReader;
use std::{
    fmt::Display,
    num::NonZero,
    path::{Path, PathBuf},
    process::Command,
    str::FromStr,
};

/// Find the path to the system's SSL certificate file.
///
/// This function attempts to locate the CA certificate file needed for HTTPS connections.
/// First checks the standard Linux path, then falls back to using curl to discover the
/// certificate location on macOS.
///
/// # Returns
/// - `Ok(Some(path))` if a certificate file path is found
/// - `Ok(None)` if no certificate file can be found
///
/// # Errors
/// Returns an error if the curl command fails to execute properly.
fn find_cert() -> Result<Option<String>> {
    const STANDARD_LINUX_STR: &str = "/etc/ssl/certs/ca-certificates.crt";
    let Ok(standard_linux_pathbuf) = PathBuf::from_str(STANDARD_LINUX_STR);
    if standard_linux_pathbuf.exists() {
        // If the standard linux path exists, use it
        Ok(Some(STANDARD_LINUX_STR.to_string()))
    } else {
        // Likely we're on a mac. Use 'curl -v' to arbitrary known good URL and extract the location
        // of the certs file it used.
        let output = Command::new("curl")
            .arg("-v")
            .arg("https://www.google.com")
            .output()?;
        for line in String::from_utf8(output.stderr)?.lines() {
            let mut words = line.split_whitespace();
            if Some("*") == words.next()
                && Some("CAfile:") == words.next()
                && let Some(cert_path) = words.next()
            {
                return Ok(Some(cert_path.to_owned()));
            }
        }
        Ok(None)
    }
}

/// Get a BAM reader (also reads SAM and CRAM). Set threads for reading.
pub fn get_bam_reader<P1, P2>(
    input: P1,
    reference_fasta: Option<P2>,
    threads: NonZero<usize>,
) -> Result<Reader>
where
    P1: AsRef<Path>,
    P2: AsRef<Path>,
{
    let mut reader = match PathType::from_path(input)? {
        PathType::Pipe => Reader::from_stdin(),
        PathType::UrlPath(url) => {
            if env::var("CURL_CA_BUNDLE").is_err() {
                // Needed to ensure that certificats are up to date
                if let Some(cert_path) = find_cert()? {
                    if env::set_var("CURL_CA_BUNDLE", cert_path).is_none() {
                        warn!("Unable to update certificates, remote read may fail.");
                    }
                } else {
                    warn!("Unable to find current cert path");
                }
            }
            Reader::from_url(&url)
        }
        PathType::FilePath(file_path) => Reader::from_path(file_path),
    }?;
    reader.set_threads(threads.into())?;
    if let Some(fasta) = reference_fasta {
        reader.set_reference(fasta)?;
    }
    Ok(reader)
}

/// Get a FASTQ reader, set threads for decompression.
pub fn get_seq_io_fastq_reader<P>(
    input: P,
    threads: NonZero<usize>,
) -> Result<SeqIoFastqReader<MaybeCompressedReader>>
where
    P: AsRef<Path>,
{
    let reader = MaybeCompressedReader::new(input, threads)?;
    Ok(SeqIoFastqReader::new(reader))
}

/// Get a FASTQ reader, set threads for decompression.
pub fn get_fastq_reader<P>(
    input: P,
    threads: NonZero<usize>,
) -> Result<FastqReader<MaybeCompressedReader>>
where
    P: AsRef<Path>,
{
    let reader = MaybeCompressedReader::new(input, threads)?;
    Ok(FastqReader::new(reader))
}

/// Get a FASTQ writer, set threads for compression.
pub fn get_seq_io_fastq_writer<P>(
    output: P,
    compression: Option<u32>,
    threads: NonZero<usize>,
) -> Result<MaybeCompressedWriter>
where
    P: AsRef<Path>,
{
    let compressed = if let Some(ref compression_level) = compression {
        *compression_level > 0
    } else {
        false
    };
    MaybeCompressedWriter::new(output, compressed, threads)
}

/// Get a FASTQ writer, set threads for compression.
pub fn get_fastq_writer<P>(
    output: P,
    compression: Option<u32>,
    threads: NonZero<usize>,
) -> Result<FastqWriter<MaybeCompressedWriter>>
where
    P: AsRef<Path>,
{
    let compressed = if let Some(ref compression_level) = compression {
        *compression_level > 0
    } else {
        false
    };
    let inner = MaybeCompressedWriter::new(output, compressed, threads)?;
    Ok(FastqWriter::new(inner))
}

/// Enum for distinguishing between FASTQ and SAM/BAM/CRAM record formats.
#[derive(PartialEq, Debug, Clone)]
pub enum RecordType {
    /// FASTQ format (with extensions .fq, .fastq, .gz, .bgz)
    Fastq,
    /// SAM/BAM/CRAM format (with extensions .bam, .sam, .cram)
    Bam,
}

impl Display for RecordType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            RecordType::Fastq => write!(f, "FASTQ"),
            RecordType::Bam => write!(f, "SAM/BAM/CRAM"),
        }
    }
}

impl RecordType {
    /// Detect the record type from a file path extension.
    ///
    /// # Arguments
    /// * `path` - The file path to analyze
    ///
    /// # Returns
    /// `Some(RecordType)` if the extension is recognized, `None` otherwise.
    pub fn from_path<P>(path: P) -> Option<RecordType>
    where
        P: AsRef<Path>,
    {
        if let Some(extension) = path.as_ref().extension() {
            Self::from_extension(extension.to_str())
        } else {
            None
        }
    }

    /// Detect the record type from a file extension string.
    ///
    /// Recognizes FASTQ extensions (.fq, .fastq, .gz, .bgz) and SAM/BAM/CRAM extensions
    /// (.bam, .sam, .cram).
    ///
    /// # Arguments
    /// * `extension` - The file extension (without leading dot)
    ///
    /// # Returns
    /// `Some(RecordType)` if the extension is recognized, `None` otherwise.
    pub fn from_extension(extension: Option<&str>) -> Option<RecordType> {
        if let Some(extension) = extension {
            match extension.to_ascii_lowercase().as_str() {
                "fq" | "fastq" | "gz" | "bgz" => Some(RecordType::Fastq),
                "bam" | "sam" | "cram" => Some(RecordType::Bam),
                _ => None,
            }
        } else {
            None
        }
    }
}
