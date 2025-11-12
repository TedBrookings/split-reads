use std::io::{BufRead, Result, Seek};

/// Struct for splitting a buffered reader by a delimiter byte,
#[derive(Debug)]
pub struct Split<B> {
    buf: B,
    delim: u8,
}

impl<B: BufRead> Split<B> {
    /// Create a new Split iterator that yields records delimited by the given byte.
    ///
    /// # Arguments
    /// * `buf` - The buffered reader to split
    /// * `delim` - The delimiter byte (e.g., b'\n' for newline-delimited records)
    pub fn new(buf: B, delim: u8) -> Self {
        Self { buf, delim }
    }
}

/// impl seek
impl<B: BufRead + Seek> Seek for Split<B> {
    fn seek(&mut self, pos: std::io::SeekFrom) -> Result<u64> {
        self.buf.seek(pos)
    }
}

/// impl Iterator yielding Result<Vec<u8>>
impl<B: BufRead> Iterator for Split<B> {
    type Item = Result<Vec<u8>>;

    fn next(&mut self) -> Option<Result<Vec<u8>>> {
        let mut buf = Vec::new();
        match self.buf.read_until(self.delim, &mut buf) {
            Ok(0) => None,
            Ok(_n) => {
                if buf[buf.len() - 1] == self.delim {
                    buf.pop();
                }
                Some(Ok(buf))
            }
            Err(e) => Some(Err(e)),
        }
    }
}
