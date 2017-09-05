# Assembltrie
Assembltrie is a software tool for compressing collections of (fixed length) Illumina reads, written in C++14 and is availble under an open-source license. Currently, Assembltrie is the only FASTQ compressor that approaches the information theory limit for a given short read collection uniformly sampled from an underlying reference genome. Assembltrie becomes the first FASTQ compressor that achieves both combinatorial optimality and information theoretic optimality under fair assumptions.

## Installation 
Assembltrie is suggested to compile with 
- GCC version 5.0 or higher (or equivalent GCC version that supports at least C++14)
- Intel version 17.0 or higher

to generate reliable compression performance. Note that unlike most existing software tools, Assembltrie does not depend on any down-stream compressors, such as `gzip` or `bzip2` so it is not necessary to install them. The tentative building process is as simple as
```
cd assembltrie
make
```
which will give an executable program `astrie`.

## Usage
To run Assembltrie from command line, type
```
astrie -c -i <input.fastq> -o <result> [options]
```
for compression; and
```
astrie -d -i <input.out> -o <result> [options]
```
for decompression.

 **Compression.** The mode `-c` implies to compress the input FASTQ file, generating two separate binary output (compressed) files: one named `result.out`, containing the encoding of assembled reads; the other named `part.out`, containing the encoding of singletons as well as other meta information necessary for decompression. In addition, `options` specifies the following mandatory and selectable parameters:
 - `-L <integer>` **Mandatory**, specifies the (fixed) read length in one compression run, the maximum value is `L = 250`
 - `-K <integer>` **Mandatory**, specifies the minimum overlap length/hash length, the suggested value is `floor(L / 5)` for `L = 100`
 - `-h 0 | 1 | 2` Optional, `h = 0` ignores any strand correction heuristic; `h = 2` applies our greedy strand correction heuristic
 - `-s <integer>` Optional, accelerates potential children search by ignoring the already processed reads with suffix length less than or equal to `integer`, and the suggested value is `floor(L / cov)`, where `cov` denotes the coverage of the input read collection (FASTQ file)
 - `-e <integer>` Optional, the maximum allowed mismatches for read overlaps. The default value is `3`, but `4` is strongly recommended for `L = 100` and `6` for `L = 150` 
 - `-n <integer>` Optional, the number of working threads, defualt value is `n = 8`
 
 **Decompression.** The mode `-d` implies to decompress the input compressed file `input.out` plus the available `part.out` into `result.fasta`, which contains a permutation (according to their locations in the constructed read forest) of (the sequence content only) of the original uncompressed read collection. To properly decompress `input.out`, Assembltrie expects the following parameters
 - `-L <integer>` **Mandatory**, the (fixed) read length in one compression run, should be the same as what is specified in the compression process.
 - `-h 0 | 1 | 2` **Mandatory**, although in Assembltrie's compression process it's optional. Again, it should follow what is specified in the compression process.

**Sample Usage**
```
(export PATH=.:$PATH)
astrie -c -L100 -K20 -h0 -s4 -iSRR554369_1.fastq -oSRR554369_1.out -e4 -n8
astrie -d -L100 -h0 -iSRR554369_1.out -oSRR554369_1.fasta
```
