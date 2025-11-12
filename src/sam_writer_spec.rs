use crate::path_type::PathType;
use anyhow::{Result, anyhow};
use rust_htslib::bam::{CompressionLevel, Format, Header, Read, Reader, Writer};
use std::{num::NonZero, path::Path};

/// Convert a format string to an htslib Format enum.
///
/// # Arguments
/// * `format_str` - The format string ("bam", "cram", or "sam")
///
/// # Errors
/// Returns an error if the format string is not recognized.
fn str_to_format(format_str: &str) -> Result<Format> {
    match format_str.to_ascii_lowercase().as_str() {
        "bam" => Ok(Format::Bam),
        "cram" => Ok(Format::Cram),
        "sam" => Ok(Format::Sam),
        _ => Err(anyhow!("Unknown Sam format: {format_str}")),
    }
}

/// Get the appropriate output format from the specified output path.
///
/// If the output path has a recognized extension (.bam, .cram, .sam), uses that format.
/// Otherwise falls back to the provided default format string.
///
/// # Arguments
/// * `output` - The output path
/// * `default_format` - The default format string to use if extension is not recognized
///
/// # Errors
/// Returns an error if neither the extension nor default format is recognized.
pub fn get_format<P>(output: P, default_format: String) -> Result<Format>
where
    P: AsRef<Path>,
{
    if let Some(extension) = output.as_ref().extension() {
        if let Some(extension_str) = extension.to_str() {
            str_to_format(extension_str).or_else(|_| str_to_format(&default_format))
        } else {
            str_to_format(&default_format)
        }
    } else {
        str_to_format(&default_format)
    }
}

/// Options for configuring a SAM/BAM/CRAM writer.
///
/// This builder-style struct allows setting optional parameters for writing SAM/BAM/CRAM files.
#[derive(Clone, Debug)]
pub struct SamWriterOptions<P> {
    /// Path to reference FASTA file (required for CRAM format)
    reference_fasta: Option<P>,
    /// Number of threads for compression
    threads: Option<NonZero<usize>>,
    /// Compression level (0-9)
    compression: Option<u32>,
}

/// Builder for creating a SAM/BAM/CRAM writer with custom configuration.
///
/// This struct provides a fluent API for configuring a SAM/BAM/CRAM writer before creation.
#[derive(Clone, Debug)]
pub struct SamWriterSpec<P> {
    /// Output file path
    output: P,
    /// SAM/BAM/CRAM header (required for writer creation)
    header: Option<Header>,
    /// Output format (required for writer creation)
    format: Option<Format>,
    /// Additional writer options that may remain unspecified
    options: SamWriterOptions<P>,
}

impl<P> SamWriterOptions<P>
where
    P: AsRef<Path> + Clone,
{
    /// Create a new SamWriterOptions with all options set to None.
    pub fn new() -> Self {
        Self {
            reference_fasta: None,
            threads: None,
            compression: None,
        }
    }

    /// Set the reference FASTA file path (required for CRAM format).
    pub fn reference_fasta(&mut self, reference_fasta: P) -> &mut Self {
        self.reference_fasta = Some(reference_fasta);
        self
    }

    /// Set the number of threads to use for compression.
    pub fn threads(&mut self, threads: NonZero<usize>) -> &mut Self {
        self.threads = Some(threads);
        self
    }

    /// Set the compression level (0-9).
    pub fn compression(&mut self, compression: u32) -> &mut Self {
        self.compression = Some(compression);
        self
    }
}

impl<P> Default for SamWriterOptions<P>
where
    P: AsRef<Path> + Clone,
{
    fn default() -> Self {
        Self::new()
    }
}

impl<P> SamWriterSpec<P>
where
    P: AsRef<Path> + Clone,
{
    /// Create a new SamWriterSpec for the given output path.
    pub fn new(output: P) -> Self {
        Self {
            output,
            header: None,
            format: None,
            options: SamWriterOptions::new(),
        }
    }

    /// Set the SAM/BAM/CRAM header.
    pub fn header(&mut self, header: Header) -> &mut Self {
        self.header = Some(header);
        self
    }

    /// Set the header by copying from an existing SAM/BAM/CRAM reader.
    pub fn header_from_reader(&mut self, reader: &Reader) -> &mut Self {
        self.header(Header::from_template(reader.header()))
    }

    /// Set the output format (BAM, CRAM, or SAM).
    pub fn format(&mut self, format: Format) -> &mut Self {
        self.format = Some(format);
        self
    }

    /// Set the output format based on the output path, with a fallback default format.
    pub fn format_from_path_or_default(&mut self, default: String) -> Result<&mut Self> {
        Ok(self.format(get_format(self.output.clone(), default)?))
    }

    /// Set the reference FASTA file path (required for CRAM format).
    pub fn reference_fasta(&mut self, reference_fasta: Option<P>) -> &mut Self {
        if let Some(ref fasta) = reference_fasta {
            self.options.reference_fasta(fasta.clone());
        }
        self
    }

    /// Set the number of threads to use for compression.
    pub fn threads(&mut self, threads: NonZero<usize>) -> &mut Self {
        self.options.threads(threads);
        self
    }

    /// Set the compression level (0-9).
    pub fn compression(&mut self, compression: Option<u32>) -> &mut Self {
        if let Some(c) = compression {
            self.options.compression(c);
        }
        self
    }

    /// Create and return a configured SAM/BAM/CRAM writer.
    ///
    /// # Errors
    /// Returns an error if the format or header has not been specified, or if the writer
    /// cannot be created.
    pub fn get_bam_writer(&self) -> Result<Writer> {
        match (self.format, &self.header) {
            (Some(ref format), Some(header)) => {
                let mut compression = self.options.compression;
                let mut writer = match PathType::from_path(self.output.as_ref())? {
                    PathType::Pipe => {
                        if compression.is_none() {
                            compression = Some(0);
                        }
                        Ok(Writer::from_stdout(header, *format)?)
                    }
                    PathType::UrlPath(_) => Err(anyhow!("Cannot write directly to a cloud URL")),
                    PathType::FilePath(file_path) => {
                        Ok(Writer::from_path(file_path, header, *format)?)
                    }
                }?;
                if let Some(threads) = self.options.threads {
                    writer.set_threads(threads.into())?;
                }
                if let Some(ref fasta) = self.options.reference_fasta {
                    writer.set_reference(fasta.as_ref())?;
                }
                if let Some(c) = compression {
                    writer.set_compression_level(CompressionLevel::Level(c))?;
                }
                Ok(writer)
            }
            (None, _) => Err(anyhow!("format was not specified for SamWriterSpec")),
            (_, None) => Err(anyhow!("header was not specified for SamWriterSpec")),
        }
    }
}
