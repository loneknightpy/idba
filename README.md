Copyright (C) 2009 Yu Peng (loneknightpy@gmail.com)

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

## Requirements

This software is suitable for all unix-like system with gcc installed.


## Installation Guide

### If you checkout the repo directly.
$ ./build.sh

### If you use the release package.
Extract the package, then use make to compile the source code.
```
$ ./configure
$ make
```

### Install with bioconda
To install with [bioconda](https://bioconda.github.io/), refer to: https://anaconda.org/bioconda/idba.


## Introduction

IDBA is the basic iterative de Bruijn graph assembler for second-generation sequencing reads.
IDBA-UD, an extension of IDBA, is designed to utilize paired-end reads to assemble low-depth
regions and use progressive depth on contigs to reduce errors in high-depth regions. It is a
generic purpose assembler and especially good for single-cell and metagenomic sequencing data.
IDBA-Hybrid is another update version of IDBA-UD, which can make use of a similar reference
genome to improve assembly result. IDBA-Tran is an iterative de Bruijn graph assembler for
RNA-Seq data.

The basic IDBA is included for comparison, you should use more specific assemblers for your data.
If you are assembling genomic data without reference, please use IDBA-UD.
If you are assembling genomic data with a similar reference genome, please use IDBA-Hybrid.
If you are assembling transcriptome data, please use IDBA-Tran.


## Comments

Note that IDBA assemblers are designed for short reads (around 100bp). If you want to assemble
paired-end reads with longer read length, please modify the constant kMaxShortSequence in
src/sequence/short_sequence.h to support longer read length.

Please find the manual by running the assembler without any parameters. For example:
```
$ bin/idba
```

IDBA series assemblers accept FASTA format reads. FASTQ format reads can be converted by
fq2fa program in the package.
```
$ bin/fq2fa read.fq read.fa
```

IDBA-UD, IDBA-Hybrid and IDBA-Tran require paired-end reads stored in the same FASTA file. A pair of
reads should be in two consecutive lines. If not, please use fq2fa to merge two
FASTQ read files to single file.
```
$ bin/fq2fa --merge --filter read_1.fq read_2.fq read.fa
```
or convert a FASTQ read file to FASTA file.
```
$ bin/fq2fa --paired --filter read.fq read.fa
```

This tool assumes that the paired-end reads are in order (->, <-). If your data is in reverse order (<-, ->),
please convert it by yourself.

## IDBA on Docker
A docker image was built for IDBA. Please use the following command to run IDBA-UD on Docker (assuming
the input read file is in the current directory). **If you are using Mac OS and see bus error, please try
this image.**
```
$ docker run -v `pwd`:/data -w /data loneknightpy/idba idba_ud  -r read.fa -o output
```


