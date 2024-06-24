# PGS Catalog utilities


> [!IMPORTANT]
> * In June 2024 this repository was archived (made read-only) and the code here deprecated
> * The code was refactored and migrated to a new set of Python packages contained in the [pygscatalog repository](https://github.com/PGScatalog/pygscatalog)
> * [The PGS Catalog Calculator v2-beta](https://github.com/PGScatalog/pgsc_calc/releases/tag/v2.0.0-beta) uses these new Python packages
> * The Python package [pgscatalog-utils](https://pypi.org/project/pgscatalog-utils/) is unaffected by this change, and we will continue to publish updates from the new repository
> * If you experience problems with our Python tools, please [create issues at the new repository](https://github.com/PGScatalog/pygscatalog/issues/new)
 
[![CI](https://github.com/PGScatalog/pgscatalog_utils/actions/workflows/main.yml/badge.svg)](https://github.com/PGScatalog/pgscatalog_utils/actions/workflows/main.yml)
[![DOI](https://zenodo.org/badge/513521373.svg)](https://zenodo.org/badge/latestdoi/513521373)

This repository is a collection of useful tools for downloading and working with scoring files from the
PGS Catalog. This is mostly used internally by the PGS Catalog Calculator ([`PGScatalog/pgsc_calc`](https://github.com/PGScatalog/pgsc_calc)); however, other users may find some of these tools helpful.

## Overview

* `download_scorefiles`: Download scoring files by PGS ID (accession) in genome builds GRCh37 or GRCh38
* `combine_scorefile`: Combine multiple scoring files into a single scoring file
in 'long' format
* `match_variants`: Match target variants (bim or pvar files) against the output
of `combine_scorefile` to produce scoring files for plink 2
* `ancestry_analysis` : use genetic PCA loadings to compare samples to population reference panels, and report PGS adjusted for these axes of genetic ancestry. The PCs will likely have been generated with [FRAPOSA (pgs catalog version)](https://github.com/PGScatalog/fraposa_pgsc)
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
If 
you use the tool we ask you to cite our paper describing software and updated PGS Catalog resource:

- >Lambert, Wingfield _et al._ (2024) The Polygenic Score Catalog: new functionality
  and tools to enable FAIR research.  medRxiv.
  doi:[10.1101/2024.05.29.24307783](https://doi.org/10.1101/2024.05.29.24307783).
  
This work has received funding from EMBL-EBI core funds, the Baker Institute, the University of Cambridge, 
Health Data Research UK (HDRUK), and the European Union's Horizon 2020 research and innovation programme 
under grant agreement No 101016775 INTERVENE.
