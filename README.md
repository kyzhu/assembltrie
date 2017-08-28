# Assembltrie
Assembltrie is a software tool for compressing collections of (fixed length) Illumina reads, written in C++14 and is availble under an open-source license. Currently, Assembltrie is the only FASTQ compressor that approaches the information theory limit for a given short read collection uniformly sampled from an underlying reference genome. Assembltrie becomes the first FASTQ compressor that achieves both combinatorial optimality and information theoretic optimality under fair assumptions.

## Installation 
## Usage
To run Assembltrie from command line, type
```
astrie -c|d -i <input.fastq> -o <result> [options]
```
where the option `-c` implies to compress the input FASTQ file, generating two separate binary output (compressed) files: one named `result.out`, containing the encoding of assembled reads; the other named `part.out`, containing the encoding of singletons as well as other meta informations necessary for decompression. The option `-d` implies to decompress the input compressed file 
