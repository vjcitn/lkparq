<div id="main" class="col-md-9" role="main">

# Convenience functions for common gene types

<div class="ref-description section level2">

Helper functions to quickly extract genes of common biotypes.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
protein_coding_genes(x, ...)

lncRNA_genes(x, ...)
```

</div>

</div>

<div class="section level2">

## Arguments

-   x:

    A `GTFParquet` object.

-   ...:

    Additional arguments passed to `genes,GTFParquet-method`.

</div>

<div class="section level2">

## Value

A [GRanges](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
object.

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
if (FALSE) { # \dontrun{
gtf <- GTFParquet(system.file("gc49", package="lkparq"))
pc <- protein_coding_genes(gtf)
lnc <- lncRNA_genes(gtf)
} # }
```

</div>

</div>

</div>
