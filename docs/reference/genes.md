<div id="main" class="col-md-9" role="main">

# Extract genomic features from a GTFParquet object

<div class="ref-description section level2">

Methods to extract genomic features from a GTFParquet object as GRanges.
Unlike TxDb methods, these preserve all GTF attributes as metadata
columns.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
# S4 method for class 'GTFParquet'
genes(x, columns=NULL, filter=NULL, use_versioned_ids=FALSE)

# S4 method for class 'GTFParquet'
transcripts(x, columns=NULL, filter=NULL, use_versioned_ids=FALSE)

# S4 method for class 'GTFParquet'
exons(x, columns=NULL, filter=NULL, use_versioned_ids=FALSE)

# S4 method for class 'GTFParquet'
cds(x, columns=NULL, filter=NULL)

# S4 method for class 'GTFParquet'
transcripts(x, columns = NULL, filter = NULL, use_versioned_ids = FALSE)

# S4 method for class 'GTFParquet'
exons(x, columns = NULL, filter = NULL, use_versioned_ids = FALSE)

# S4 method for class 'GTFParquet'
cds(x, columns = NULL, filter = NULL)
```

</div>

</div>

<div class="section level2">

## Arguments

-   x:

    A `GTFParquet` object.

-   columns:

    Character vector of columns to include in `mcols`. If `NULL`
    (default), includes all available attribute columns. For genes:
    `gene_name`, `gene_type`, `source`, `level`, `tags`, `havana_gene`.
    For transcripts: `transcript_name`, `transcript_type`, `gene_id`,
    `gene_name`, `transcript_support_level`, `ccdsid`, `protein_id`.

-   filter:

    Optional named list for filtering features. Names should be column
    names, values are vectors of acceptable values. Example:
    `filter = list(gene_type = "protein_coding", chrom = "chr1")`

-   use\_versioned\_ids:

    Logical. If `TRUE`, use full versioned IDs (e.g.,
    `ENSG00000141510.18`). If `FALSE` (default), use stripped IDs (e.g.,
    `ENSG00000141510`).

</div>

<div class="section level2">

## Value

A [GRanges](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
object with:

-   Feature IDs as `names`

-   Genomic coordinates (`seqnames`, `ranges`, `strand`)

-   Genome build in `seqinfo` (e.g., "GRCh38")

-   Rich metadata in `mcols`

</div>

<div class="section level2">

## Details

These methods return
[GRanges](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
objects with feature IDs as names and rich metadata columns from the
original GTF file.

The `filter` argument enables efficient server-side filtering through
Arrow/Parquet predicate pushdown, which can dramatically improve
performance compared to subsetting after loading.

Available filter columns include:

-   `chrom`: Chromosome name

-   `gene_type`: Gene biotype (e.g., "protein\_coding", "lncRNA")

-   `transcript_type`: Transcript biotype

-   `level`: Annotation confidence (1=verified, 2=manual, 3=automatic)

-   `source`: Annotation source ("HAVANA", "ENSEMBL")

</div>

<div class="section level2">

## See also

<div class="dont-index">

-   `GTFParquet-class` for the class definition

-   `transcriptsBy,GTFParquet-method` for grouped extraction

-   `genes` for the generic

</div>

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
if (FALSE) { # \dontrun{
gtf <- GTFParquet(system.file("gc49", package="lkparq"))

# Extract all genes with full attributes
gr <- genes(gtf)
mcols(gr)  # gene_name, gene_type, level, tags, source, havana_gene

# Filter by gene type
pc <- genes(gtf, filter = list(gene_type = "protein_coding"))
lnc <- genes(gtf, filter = list(gene_type = "lncRNA"))

# Combine filters
pc_chr1 <- genes(gtf, filter = list(gene_type = "protein_coding", chrom = "chr1"))

# Select specific columns only
gr <- genes(gtf, columns = c("gene_name", "gene_type"))

# Use versioned IDs
gr <- genes(gtf, use_versioned_ids = TRUE)
names(gr)[1]  # "ENSG00000141510.18"

# Transcripts with support level
tx <- transcripts(gtf)
high_conf <- tx[mcols(tx)$transcript_support_level == "1"]

# Exons
ex <- exons(gtf, filter = list(chrom = "chr1"))

# CDS with protein IDs
cds_gr <- cds(gtf)
mcols(cds_gr)$protein_id
} # }
```

</div>

</div>

</div>
