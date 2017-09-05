# Assembltrie
Assembltrie is a software tool for compressing collections of (fixed length) Illumina reads, written in C++14 and is availble under an open-source license. Currently, Assembltrie is the only FASTQ compressor that approaches the information theory limit for a given short read collection uniformly sampled from an underlying reference genome. Assembltrie becomes the first FASTQ compressor that achieves both combinatorial optimality and information theoretic optimality under fair assumptions.

## Installation 
Assembltrie is suggested to compile with - GCC version 5.0 or higher (or equivalent GCC version that supports at least C++14)- Intel version 17.0 or higher
to generate reliable compression performance. Note that unlike most existing software tools, Assembltrie does not depend on any down-stream compressors, such as `gzip` or `bzip2` so it is not necessary to install them. The tentative building process is as simple as```cd assembltriemake```which will give an executable program `astrie`.

## Usage
To run Assembltrie from command line, type
```astrie -c -i <input.fastq> -o <result> [options]```
for compression; and
```astrie -d -i <input.out> -o <result> [options]```
for decompression.
 
