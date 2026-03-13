<div id="main" class="col-md-9" role="main">

# Get GTF metadata

<div class="ref-description section level2">

Retrieve metadata from the GTF file header, including provider, version,
date, and genome build.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
gtf_metadata(x)
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

A named character vector of metadata key-value pairs.

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
if (FALSE) { # \dontrun{
gtf <- GTFParquet(system.file("gc49", package="lkparq"))
gtf_metadata(gtf)
#      provider         format           date         genome 
#      "GENCODE"          "gtf"   "2025-07-08"       "GRCh38"
} # }
```

</div>

</div>

</div>
