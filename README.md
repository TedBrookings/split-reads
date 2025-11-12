# split-reads

## Summary

Toolkit to construct a split-index (.si) file and split query-grouped SAM/BAM/CRAM or FASTQ into contiguous chunks.

## Rationale

It's often helpful to split unmapped SAM files for alignment or other compute-intensive embarassingly-parallel tasks. The simplest approach is just to split in advance and pass chunked BAMS (e.g. from original `my-reads.bam`, produce `my-reads.1.bam`, `myreads.2.bam`, ..., `my-reads.N.bam`). However this approach uses a lot of IO and as such is slow. In a cloud workflow these files must be delocalized as well. Producing a split index is much faster, uses almost no IO, and allows chunk extraction on the fly at almost no cost.

## Installation

Currently only via `git clone` and `cargo build`. I will publish and add to bioconda when features stabalize.

## Usage

To produce a split index, run

```sh
split-reads index -i my-reads.bam
```

This will produce `my-reads.bam.si`. Tool arguments allow overriding default number of CPUs, output index path, etc. This tool can also index remote files (http, ftp, s3, gcs) similar to `samtools`.

Then to extract e.g. chunk 3/10 from a BAM, run

```sh
split-reads get-chunk -i my-reads.bam -c 3 -n 10 -o my-reads.3.bam
```

Chunk indices run from `0` to `num_chunks - 1`. To stream the chunk to stdout for downstream
processing, just omit the output file or set it to `-`. _This is the usual way you will want to use
this tool._

```sh
split-reads get-chunk -i my-reads.bam -c 3 -n 10 | my-aligner ...
```

## Advanced Usage - Plan chunks by number of reads or queries

If you wish to plan the number of chunks to e.g. be a pre-set number of queries, you can use the
`tell` command to plan the number of chunks.

```sh
num_queries=(split-reads tell -I my-reads.bam.si)
# Plan the least number of chunks with at most 1000 queries per chunk.
num_chunks=$((num_queries / 1000))
```

`tell` can also reveal the number of reads or chunks.

## Advanced Usage - Pass-through indexing

You may want to get a split-indexed bam after some amount of processing. `split-reads index` has
an `--output` argument to pass through the reads to an output BAM file while indexing, to avoid the
needless IO of writing to disk then reading from it. For instance if you want
to split primary alignments from a mapped bam, you could run

```sh
samtools collate -Ouf aligned.bam | split-reads index -i - -o collated.bam
```

Which will produce `collated.bam` and `collated.bam.si`.
At the moment, the output type must be the same as the input (SAM/BAM/CRAM/FASTQ). Compression level can change though.
