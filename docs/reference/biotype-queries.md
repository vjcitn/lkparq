<div id="main" class="col-md-9" role="main">

# List available gene or transcript types

<div class="ref-description section level2">

Query the available biotypes in the annotation and their counts.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
gene_types(x)

transcript_types(x)
```

</div>

</div>

<div class="section level2">

## Arguments

-   x:

    A `GTFParquet` object.

</div>

<div class="section level2">

## Value

A table of biotype counts.

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
if (FALSE) { # \dontrun{
gtf <- GTFParquet(system.file("gc49", package="lkparq"))
gene_types(gtf)
# protein_coding      lncRNA     pseudogene ...
#         19950        16880          15200 ...

transcript_types(gtf)
} # }
```

</div>

</div>

</div>
