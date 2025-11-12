use std::io::{Read, Result, Seek, SeekFrom};

/// Struct to make something similar to std::io::Chain that is seekable.
pub struct Chain<R1: Read, R2: Read> {
    front: R1,
    back: R2,
    past_front: bool,
    initial_front_pos: u64,
    front_len: u64,
    initial_back_pos: u64,
}

/// impl new ChainReader
impl<R1: Read + Seek, R2: Read + Seek> Chain<R1, R2> {
    /// Create a new seekable Chain from two readers.
    ///
    /// The Chain will read from the front reader first, then continue with the back reader.
    /// Both readers must support seeking.
    ///
    /// # Arguments
    /// * `front` - The first reader to read from
    /// * `back` - The second reader to read from after the front reader is exhausted
    ///
    /// # Errors
    /// Returns an error if either reader's stream position cannot be determined.
    pub fn new(mut front: R1, mut back: R2) -> Result<Self> {
        let initial_front_pos = front.stream_position()?;
        let front_len = front.seek(SeekFrom::End(0))? - initial_front_pos;
        front.seek(SeekFrom::Start(initial_front_pos))?;
        let initial_back_pos = back.stream_position()?;
        Ok(Chain {
            front,
            back,
            past_front: initial_front_pos >= front_len,
            initial_front_pos,
            front_len,
            initial_back_pos,
        })
    }
}

/// impl Seek trait
impl<R1: Read + Seek, R2: Read + Seek> Seek for Chain<R1, R2> {
    fn stream_position(&mut self) -> Result<u64> {
        if self.past_front {
            Ok(self.front_len + self.back.stream_position()? - self.initial_back_pos)
        } else {
            Ok(self.front.stream_position()? - self.initial_front_pos)
        }
    }

    fn seek(&mut self, pos: std::io::SeekFrom) -> std::io::Result<u64> {
        match pos {
            SeekFrom::Start(pos_from_start) => {
                // we can do this correctly
                if pos_from_start >= self.front_len {
                    self.past_front = true;
                    self.back.seek(SeekFrom::Start(
                        pos_from_start - self.front_len + self.initial_back_pos,
                    ))?;
                } else {
                    self.back.seek(SeekFrom::Start(self.initial_back_pos))?;
                    self.past_front = false;
                    self.front
                        .seek(SeekFrom::Start(pos_from_start + self.initial_front_pos))?;
                }
                Ok(pos_from_start)
            }
            SeekFrom::Current(pos_from_current) => {
                if pos_from_current == 0 {
                    // this just gets stream position, return custom override
                    self.stream_position()
                } else {
                    // Compute position from start and convert to SeekFrom::start
                    let pos_from_start: u64 = if pos_from_current > 0 {
                        self.stream_position()? + (pos_from_current as u64)
                    } else {
                        self.stream_position()? - (-pos_from_current as u64)
                    };
                    self.seek(SeekFrom::Start(pos_from_start))
                }
            }
            SeekFrom::End(pos_from_end) => {
                if pos_from_end >= 0 {
                    // Seek exactly to the end (or past it, which should throw error)
                    self.past_front = true;
                    self.front.seek(SeekFrom::End(0))?;
                    Ok(self.front_len + self.back.seek(pos)? - self.initial_back_pos)
                } else {
                    // Get stream length and convert to SeekFrom::start
                    let stream_length =
                        self.front_len + self.back.seek(SeekFrom::End(0))? - self.initial_back_pos;
                    let pos_from_start = stream_length - (-pos_from_end as u64);
                    self.seek(SeekFrom::Start(pos_from_start))
                }
            }
        }
    }
}

/// impl Read trait
impl<R1: Read, R2: Read> Read for Chain<R1, R2> {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        if self.past_front {
            // Just read from the back
            self.back.read(buf)
        } else {
            // Read from the front
            let num = self.front.read(buf)?;
            if 0 == num && !buf.is_empty() {
                self.past_front = true;
                self.back.read(buf)
            } else {
                Ok(num)
            }
        }
    }
}
