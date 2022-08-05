# PGS Catalog utilities

[![CI](https://github.com/PGScatalog/pgscatalog_utils/actions/workflows/main.yml/badge.svg)](https://github.com/PGScatalog/pgscatalog_utils/actions/workflows/main.yml)

This repository is a collection of useful tools for working with data from the
PGS Catalog. This is mostly used internally by the PGS Catalog calculator, but
other users might find some of these tools helpful.

## Overview

* `download_scorefiles`: Download scoring files by PGS ID (accession) in genome builds GRCh37 or GRCh38
* `combine_scorefile`: Combine multiple scoring files into a single scoring file
in 'long' format
* `match_variants`: Match target variants (bim or pvar files) against the output
of `combine_scorefile` to produce scoring files for plink 2

## Installation

```
$ pip install pgscatalog-utils
```

Or clone the repo:

```
$ git clone https://github.com/PGScatalog/pgscatalog_utils.git
```

## Quickstart

```
$ download_scorefiles -i PGS000922 PGS001229 -o . -b GRCh37
$ combine_scorefiles -s PGS*.txt.gz -o combined.txt 
$ match_variants -s combined.txt -t <example.pvar> --min_overlap 0.75 --outdir .
```
