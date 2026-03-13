<div id="main" class="col-md-9" role="main">

# Find features overlapping a genomic region

<div class="ref-description section level2">

Efficient region queries that use chromosome-based filtering before
computing overlaps.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
genes_in_region(x, region, ...)

transcripts_in_region(x, region, ...)
```

</div>

</div>

<div class="section level2">

## Arguments

-   x:

    A `GTFParquet` object.

-   region:

    A
    [GRanges](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
    object specifying the query region(s).

-   ...:

    Additional arguments passed to `genes()` or `transcripts()`.

</div>

<div class="section level2">

## Value

A [GRanges](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
object containing features that overlap the query region.

</div>

<div class="section level2">

## Details

These functions first filter by chromosome (using Parquet predicate
pushdown for efficiency), then compute overlaps using `findOverlaps`.

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
if (FALSE) { # \dontrun{
gtf <- GTFParquet(system.file("gc49", package="lkparq"))

# Define a query region
region <- GRanges("chr1", IRanges(1000000, 2000000))

# Find overlapping genes
genes_in_region(gtf, region)

# Find overlapping transcripts (protein-coding only)
transcripts_in_region(gtf, region, 
                      filter = list(transcript_type = "protein_coding"))
} # }
```

</div>

</div>

</div>
