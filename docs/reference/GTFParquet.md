<div id="main" class="col-md-9" role="main">

# Create a GTFParquet object

<div class="ref-description section level2">

Create a GTFParquet object

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
GTFParquet(path)
```

</div>

</div>

<div class="section level2">

## Arguments

-   path:

    Path to directory containing Parquet files from gtf\_to\_parquet.py

</div>

<div class="section level2">

## Value

A GTFParquet S4 object

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
gtf <- GTFParquet(system.file("gc49", package="lkparq"))
genes(gtf)
#> GRanges object with 78691 ranges and 6 metadata columns:
#>                   seqnames            ranges strand |       gene_name
#>                      <Rle>         <IRanges>  <Rle> |     <character>
#>   ENSG00000290825     chr1       11121-24894      + |        DDX11L16
#>   ENSG00000223972     chr1       12010-13670      + |         DDX11L1
#>   ENSG00000310526     chr1       14356-30744      - |          WASH7P
#>   ENSG00000227232     chr1       14696-24886      - |          WASH7P
#>   ENSG00000278267     chr1       17369-17436      - |       MIR6859-1
#>               ...      ...               ...    ... .             ...
#>   ENSG00000292369     chrY 57190738-57208756      + | ENSG00000292369
#>   ENSG00000292370     chrY 57201135-57203405      - |          WASIR1
#>   ENSG00000292372     chrY 57207346-57212230      + |          WASH6P
#>   ENSG00000310542     chrY 57208516-57215634      - | ENSG00000310542
#>   ENSG00000292371     chrY 57212184-57214397      - |        DDX11L16
#>                                gene_type      source     level
#>                              <character> <character> <integer>
#>   ENSG00000290825                 lncRNA      HAVANA      <NA>
#>   ENSG00000223972 transcribed_unproces..      HAVANA      <NA>
#>   ENSG00000310526                 lncRNA      HAVANA      <NA>
#>   ENSG00000227232 transcribed_unproces..      HAVANA      <NA>
#>   ENSG00000278267                  miRNA     ENSEMBL      <NA>
#>               ...                    ...         ...       ...
#>   ENSG00000292369                 lncRNA      HAVANA      <NA>
#>   ENSG00000292370                 lncRNA      HAVANA      <NA>
#>   ENSG00000292372         protein_coding      HAVANA      <NA>
#>   ENSG00000310542                 lncRNA      HAVANA      <NA>
#>   ENSG00000292371 unprocessed_pseudogene      HAVANA      <NA>
#>                                     tags          havana_gene
#>                              <character>          <character>
#>   ENSG00000290825    overlaps_pseudogene                 <NA>
#>   ENSG00000223972                   <NA> OTTHUMG00000000961.2
#>   ENSG00000310526 ncRNA_host;overlappi..                 <NA>
#>   ENSG00000227232      overlapping_locus OTTHUMG00000000958.1
#>   ENSG00000278267                   <NA>                 <NA>
#>               ...                    ...                  ...
#>   ENSG00000292369      overlapping_locus OTTHUMG00000184987.2
#>   ENSG00000292370                   <NA> OTTHUMG00000022676.3
#>   ENSG00000292372      overlapping_locus OTTHUMG00000022677.5
#>   ENSG00000310542 overlapping_locus;ov..                 <NA>
#>   ENSG00000292371                   <NA> OTTHUMG00000022678.1
#>   -------
#>   seqinfo: 25 sequences from GRCh38 genome; no seqlengths
genes(gtf, filter = list(gene_type = "protein_coding"))
#> GRanges object with 20097 ranges and 6 metadata columns:
#>                   seqnames            ranges strand |   gene_name
#>                      <Rle>         <IRanges>  <Rle> | <character>
#>   ENSG00000186092     chr1       65419-71585      + |       OR4F5
#>   ENSG00000284733     chr1     450740-451678      - |      OR4F29
#>   ENSG00000284662     chr1     685716-686654      - |      OR4F16
#>   ENSG00000187634     chr1     923923-944575      + |      SAMD11
#>   ENSG00000188976     chr1     943527-960714      - |       NOC2L
#>               ...      ...               ...    ... .         ...
#>   ENSG00000185894     chrY 25030901-25062548      - |       BPY2C
#>   ENSG00000172288     chrY 25622117-25624902      + |        CDY1
#>   ENSG00000292366     chrY 57067747-57130289      + |       VAMP7
#>   ENSG00000292373     chrY 57184216-57199537      + |        IL9R
#>   ENSG00000292372     chrY 57207346-57212230      + |      WASH6P
#>                        gene_type      source     level              tags
#>                      <character> <character> <integer>       <character>
#>   ENSG00000186092 protein_coding      HAVANA      <NA>              <NA>
#>   ENSG00000284733 protein_coding      HAVANA      <NA>              <NA>
#>   ENSG00000284662 protein_coding      HAVANA      <NA>              <NA>
#>   ENSG00000187634 protein_coding      HAVANA      <NA>              <NA>
#>   ENSG00000188976 protein_coding      HAVANA      <NA>              <NA>
#>               ...            ...         ...       ...               ...
#>   ENSG00000185894 protein_coding      HAVANA      <NA>              <NA>
#>   ENSG00000172288 protein_coding      HAVANA      <NA>         retrogene
#>   ENSG00000292366 protein_coding      HAVANA      <NA>              <NA>
#>   ENSG00000292373 protein_coding      HAVANA      <NA> overlapping_locus
#>   ENSG00000292372 protein_coding      HAVANA      <NA> overlapping_locus
#>                             havana_gene
#>                             <character>
#>   ENSG00000186092  OTTHUMG00000001094.4
#>   ENSG00000284733  OTTHUMG00000002860.3
#>   ENSG00000284662  OTTHUMG00000002581.3
#>   ENSG00000187634 OTTHUMG00000040719.11
#>   ENSG00000188976  OTTHUMG00000040720.2
#>               ...                   ...
#>   ENSG00000185894  OTTHUMG00000045199.3
#>   ENSG00000172288  OTTHUMG00000045274.2
#>   ENSG00000292366  OTTHUMG00000022679.3
#>   ENSG00000292373  OTTHUMG00000022720.1
#>   ENSG00000292372  OTTHUMG00000022677.5
#>   -------
#>   seqinfo: 25 sequences from GRCh38 genome; no seqlengths
```

</div>

</div>

</div>
