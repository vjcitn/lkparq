<div id="main" class="col-md-9" role="main">

# Extract and group genomic features from a GTFParquet object

<div class="ref-description section level2">

Generic functions to extract genomic features of a given type grouped
based on another type of genomic feature. These methods extend the
GenomicFeatures generics for GTFParquet objects.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
# S4 method for class 'GTFParquet'
transcriptsBy(x, by="gene", filter=NULL)

# S4 method for class 'GTFParquet'
exonsBy(x, by=c("tx", "gene"), filter=NULL)

# S4 method for class 'GTFParquet'
cdsBy(x, by=c("tx", "gene"), filter=NULL)

# S4 method for class 'GTFParquet'
cdsBy(x, by = c("tx", "gene"), filter = NULL)

# S4 method for class 'GTFParquet'
transcriptsBy(x, by = "gene", filter = NULL)
```

</div>

</div>

<div class="section level2">

## Arguments

-   x:

    A `GTFParquet` object.

-   by:

    One of `"gene"`, `"tx"` (transcript). Determines the grouping. For
    `transcriptsBy`, only `"gene"` is currently supported.

-   filter:

    Optional named list for filtering features before grouping. Names
    should be column names (e.g., `gene_type`, `chrom`), values are
    vectors of acceptable values. Example:
    `filter = list(gene_type = "protein_coding", chrom = "chr1")`

</div>

<div class="section level2">

## Value

A
[GRangesList](https://rdrr.io/pkg/GenomicRanges/man/GRangesList-class.html)
object. The names of the list elements are the IDs of the grouping
features (gene IDs or transcript IDs).

For GTFParquet objects, the names use stripped (unversioned) IDs by
default (e.g., `ENSG00000141510` rather than `ENSG00000141510.18`).

</div>

<div class="section level2">

## Details

These functions return a
[GRangesList](https://rdrr.io/pkg/GenomicRanges/man/GRangesList-class.html)
object where the ranges within each of the elements are ordered
according to the following rule:

When using `exonsBy` or `cdsBy` with `by = "tx"`, the returned exons or
CDS are ordered by ascending exon number for each transcript, that is,
by their position in the transcript. In all other cases, the ranges will
be ordered by chromosome, strand, start, and end values.

Unlike TxDb methods, GTFParquet methods preserve rich metadata columns
including `transcript_name`, `transcript_type`, `exon_number`,
`protein_id`, and `frame`.

The `filter` argument allows efficient server-side filtering before data
is loaded into R, which can dramatically improve performance for large
annotation files.

</div>

<div class="section level2">

## See also

<div class="dont-index">

-   `GTFParquet-class` for the class definition

-   `genes,GTFParquet-method` for extracting ungrouped features

-   `transcriptsBy` for the generic

-   `exonsBy` for the generic

-   `cdsBy` for the generic

</div>

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
if (FALSE) { # \dontrun{
gtf <- GTFParquet(system.file("gc49", package="lkparq"))

# Exons grouped by transcript (sorted by exon_number)
ebt <- exonsBy(gtf, by = "tx")
ebt[[1]]  # Exons for first transcript

# Exons grouped by gene
ebg <- exonsBy(gtf, by = "gene")

# CDS grouped by transcript
cbt <- cdsBy(gtf, by = "tx")

# Transcripts grouped by gene
tbg <- transcriptsBy(gtf, by = "gene")

# Filter to protein-coding only
pc_exons <- exonsBy(gtf, by = "tx", 
                    filter = list(gene_type = "protein_coding"))

# Filter by chromosome
chr1_cds <- cdsBy(gtf, by = "tx", filter = list(chrom = "chr1"))
} # }
```

</div>

</div>

</div>
