# alamem - Approximate Local Alignment via chained MEMs

## Introduction

**alamem** is a program for finding local alignments between a query and a reference database of genomes above a specified **average nucleotide identity** (ANI) and length. It is designed to work in the ANI > ~90% and length > 40bp regime, to operate as a mostly drop-in replacement for BLAT hits.

alamem uses an approximate mapping method without base-level alignment to get both hits and ANI---notably, it uses MEMs to get base-level resolution at the boundaries, but then uses pseudo-matching statistics to estimate ANI of the region covered by the MEM chains. It is over an order of magnitude faster than BLAT, but seems empirically get get fairly comparable results. alamem offers:

1. **No preprocessing**. We don't require any pre-indexing of the reference database. Instead, for every run, alamem streams over the database in plain FASTA or FASTA.gz format.

2. **Low RAM Usage**. Because alamem streams the database, most of the memory that is used is just keeping track of the query index and any hits that come out. For bacterial genomes as queries, ~4 GB of RAM is more than sufficient; for smaller queries, <1 GB of RAM is required.

3. **Fast computations**. We entirely avoid base-level dynamic programming for extension by using pseudo-matching statistics to compute ANI. This means we are just doing seed+chain, rather than seed+chain+extend. Querying a genome against an (unpreprocessed) database of >85000 prokaryotic genomes takes less than an hour in single-threaded mode with 4GB of RAM. On a 64-core machine, this allows the entire computation to happen in less than a minute.

## Updates

See the [CHANGELOG](https://github.com/yunwilliamyu/alamem/blob/main/CHANGELOG.md) for alamem's full versioning history. [TODO]

## Install

#### Option 1: Build from source

Requirements:
1. [rust](https://www.rust-lang.org/tools/install) programming language and associated tools such as cargo are required and assumed to be in PATH.

Building takes around a minute (depending on # of cores).

```sh
git clone https://github.com/yunwilliamyu/alamem
cd alamem
RUSTFLAGS="-C target-cpu=native" cargo install --path .
alamem -h
```

<!--
#### Option 2: Conda
```sh
conda install -c bioconda skani
```
-->

#### Option 2: Pre-built x86-64 linux statically compiled executable
We offer a pre-built statically compiled executable for x86-64 systems. That is, if you're on an x86-64 Linux system, you can just download the binary and run it without installing anything.

For using the latest version of alamem:
```sh
wget https://github.com/yunwilliamyu/alamem/releases/latest/download/alamem
chmod +x alamem
./alamem -h
```
**Important**: the binary runs slightly slower (~10%) most of the time, but it can sometimes be drastically slower than compiling it from scratch for your machine.

## Quick start
```sh
alamem reference_list.txt query.fna[.gz] out.txt
```

reference_list.txt should be a newline-delimited list of all genomes (gzipped is fine)

out.txt will be where results are written.

Note that alamem is almost symmetric with respect to reference and query in terms of output, but you should have the query be the shorter sequence, so that indexing is fast. alamem is fast enough that we can just stream the much larger reference database.

