use anyhow::{Result, anyhow};
use std::{
    path::{Path, PathBuf},
    str::FromStr,
};
use url::Url;

/// Enum to hold kinds of paths: pipe (stdin/stdout), local files (PathBuf) remote URLs:
pub enum PathType {
    Pipe,
    FilePath(PathBuf),
    UrlPath(Url),
}

const URL_PREFIXES: [&str; 5] = ["s3://", "gcs://", "ftp://", "http://", "https://"];

impl PathType {
    /// Form PathType enum from input PathBuf object, e.g. a clap input argument
    pub fn from_path<P>(path: P) -> Result<Self>
    where
        P: AsRef<Path>,
    {
        Ok(if let Some(path_str) = path.as_ref().to_str() {
            if path_str == "-" {
                PathType::Pipe
            } else if URL_PREFIXES
                .into_iter()
                .any(|prefix| path_str.starts_with(prefix))
            {
                PathType::UrlPath(Url::parse(path_str)?)
            } else {
                PathType::FilePath(path.as_ref().into())
            }
        } else {
            PathType::FilePath(path.as_ref().into())
        })
    }

    /// Form default index file location from path to main file
    pub fn default_index(&self, index_extension: &'static str) -> Result<Option<PathBuf>> {
        match self {
            Self::Pipe => Ok(None),
            Self::UrlPath(url) => {
                let local_si = if let Some(last_segment) = &url
                    .path_segments()
                    .and_then(|mut segments| segments.next_back())
                {
                    PathBuf::from_str(last_segment)?.with_added_extension(index_extension)
                } else {
                    Err(anyhow!("Unable to parse url {url}"))?
                };
                if local_si.exists() {
                    Ok(Some(local_si))
                } else {
                    Ok(Some(
                        PathBuf::from_str(url.as_str())?.with_added_extension(index_extension),
                    ))
                }
            }
            Self::FilePath(path_buf) => Ok(Some(path_buf.with_added_extension(index_extension))),
        }
    }
}
