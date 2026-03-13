# lkparq: Experimental R package with a TxDb-like resource for Gencode V49 for humans

This package was heavily negotiated with Claude chat Opus 4.5 extended in March 2026.

It takes advantage of works of many people involved with Bioconductor packages such as

- GenomicState
- txdbmaker
- AnnotationDbi

The main purpose is to expose GTF content as directly as possible using methods
of GenomicFeatures.  The back end at present is parquet.
