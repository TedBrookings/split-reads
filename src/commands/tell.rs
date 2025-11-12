use crate::commands::command::Command;
use anyhow::Result;
use clap::Parser;
use serde::Serialize;
use split_reads::split_index::SplitIndex;
use std::path::PathBuf;

#[derive(clap::ValueEnum, Clone, Default, Debug, Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum TellWhich {
    NumBins,
    #[default]
    NumQueries,
    NumReads,
}

/// Tell some basic stats as derived from a split-index file.
#[derive(Parser, Debug)]
#[command(version, verbatim_doc_comment)]
pub(crate) struct Tell {
    /// Input path for Index file. Use "-" for stdin.
    #[clap(long, short = 'I', required = true)]
    index: PathBuf,

    /// Number of bins to retain in final index file.
    #[clap(long, short = 't', required = false, default_value_t, value_enum)]
    tell: TellWhich,
}

impl Tell {
    /// Build the split index, then downsize to the requested number of bins and write to requested
    /// output path
    fn tell(&self) -> Result<()> {
        let split_index = SplitIndex::read(self.index.clone())?;
        match self.tell.clone() {
            TellWhich::NumBins => println!("{}", split_index.len()),
            TellWhich::NumQueries => println!("{}", split_index.num_queries()),
            TellWhich::NumReads => println!("{}", split_index.num_reads()),
        }
        Ok(())
    }
}

/// Implement the Command trait for `Tell` struct.
impl Command for Tell {
    /// Execute the tell command to print statistics from a split-index file.
    fn execute(&self) -> Result<()> {
        self.tell()
    }
}
