#![deny(unsafe_code)]

extern crate core;

pub mod commands;

use anyhow::Result;
use clap::Parser;
use commands::command::Command;
use commands::get_chunk::GetChunk;
use commands::index::Index;
use commands::tell::Tell;
use commands::test_fastq::TestFastq;
use commands::test_seq_io::TestSeqIo;
use enum_dispatch::enum_dispatch;
use std::sync::LazyLock;

#[cfg(test)]
mod test_utils;

pub mod built_info {
    // The file has been placed there by the build script.
    include!(concat!(env!("OUT_DIR"), "/built.rs"));
}

/// Get package and git version, including dirty git state.
static VERSION: LazyLock<String> = LazyLock::new(|| {
    let pkg_version = built_info::PKG_VERSION;
    if let Some(git_version) = built_info::GIT_VERSION {
        let suffix = if Some(true) == built_info::GIT_DIRTY {
            "-dirty"
        } else {
            ""
        };
        format!("{pkg_version}-{git_version}{suffix}")
    } else {
        format!("{pkg_version}-no-git")
    }
});

#[derive(Parser, Debug)]
#[clap(version = VERSION.as_str(), term_width=0)]
struct Args {
    #[clap(subcommand)]
    subcommand: Subcommand,
}

#[enum_dispatch(Command)]
#[derive(Parser, Debug)]
#[command(version)]
enum Subcommand {
    Index(Index),
    GetChunk(GetChunk),
    Tell(Tell),
    TestSeqIo(TestSeqIo),
    TestFastq(TestFastq),
}

fn main() -> Result<()> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();
    let args: Args = Args::parse();
    args.subcommand.execute()
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_noop() {}
}
