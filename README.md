# PGS Catalog utilities

[![CI](https://github.com/PGScatalog/pgscatalog_utils/actions/workflows/main.yml/badge.svg)](https://github.com/PGScatalog/pgscatalog_utils/actions/workflows/main.yml)

This repository is a collection of useful tools for downloading and working with scoring files from the
PGS Catalog. This is mostly used internally by the PGS Catalog Calculator ([`PGScatalog/pgsc_calc`](https://github.com/PGScatalog/pgsc_calc)); however, other users may find some of these tools helpful.

## Overview

* `download_scorefiles`: Download scoring files by PGS ID (accession) in genome builds GRCh37 or GRCh38
* `combine_scorefile`: Combine multiple scoring files into a single scoring file
in 'long' format
* `match_variants`: Match target variants (bim or pvar files) against the output
of `combine_scorefile` to produce scoring files for plink 2
* `validate_scorefiles`: Check/validate that the scoring files and harmonized scoring files match the PGS Catalog scoring file formats.

## Installation

```
$ pip install pgscatalog-utils
```

## Quickstart

```
$ download_scorefiles -i PGS000922 PGS001229 -o . -b GRCh37
$ combine_scorefiles -s PGS*.txt.gz -o combined.txt 
$ match_variants -s combined.txt -t <example.pvar> --min_overlap 0.75 --outdir .
$ validate_scorefiles -t formatted --dir <scoringfiles_directory> --log_dir <logs_directory>
```

More details are available using the `--help` parameter.

## Install from source

Requirements:

- python 3.10
- [poetry](https://python-poetry.org)

```
$ git clone https://github.com/PGScatalog/pgscatalog_utils.git
$ cd pgscatalog_utils
$ poetry install
$ poetry build
$ pip install --user dist/*.whl 
```

## Credits

The `pgscatalog_utils` package is developed as part of the **Polygenic Score (PGS) Catalog** 
([www.PGSCatalog.org](https://www.PGSCatalog.org)) project, a collaboration between the 
University of Cambridge’s Department of Public Health and Primary Care (Michael Inouye, Samuel Lambert, Laurent Gil) 
and the European Bioinformatics Institute (Helen Parkinson, Aoife McMahon, Ben Wingfield, Laura Harris).

A manuscript describing the tool and larger PGS Catalog Calculator pipeline 
[(`PGSCatalog/pgsc_calc`)](https://github.com/PGScatalog/pgsc_calc) is in preparation. In the meantime 
if you use these tools we ask you to cite the repo(s) and the paper describing the PGS Catalog resource:

- >PGS Catalog utilities _(in development)_. PGS Catalog
  Team. [https://github.com/PGScatalog/pgscatalog_utils](https://github.com/PGScatalog/pgscatalog_utils)
- >PGS Catalog Calculator _(in development)_. PGS Catalog
  Team. [https://github.com/PGScatalog/pgsc_calc](https://github.com/PGScatalog/pgsc_calc)
- >Lambert _et al._ (2021) The Polygenic Score Catalog as an open database for
reproducibility and systematic evaluation.  Nature Genetics. 53:420–425
doi:[10.1038/s41588-021-00783-5](https://doi.org/10.1038/s41588-021-00783-5).

This work has received funding from EMBL-EBI core funds, the Baker Institute, the University of Cambridge, 
Health Data Research UK (HDRUK), and the European Union's Horizon 2020 research and innovation programme 
under grant agreement No 101016775 INTERVENE.
