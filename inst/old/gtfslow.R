
#' Create a GTFParquet connection object
#' @rawNamespace import(GenomicRanges, except=c(setdiff, intersect, union))
#' @import arrow
#' @import dplyr
#' @param path Path to directory containing Parquet files from gtf_to_parquet.py
#' @return A GTFParquet S3 object
#' @export
#' @examples
#' gtf <- GTFParquet("path/to/parquet_output")
#' genes(gtf)
#' genes(gtf, filter = list(gene_type = "protein_coding"))
GTFParquet <- function(path) {
  path <- normalizePath(path, mustWork = TRUE)
  
  # Check which files exist
  files <- list(
    genes = file.path(path, "genes.parquet"),
    transcripts = file.path(path, "transcripts.parquet"),
    exons = file.path(path, "exons.parquet"),
    cds = file.path(path, "cds.parquet"),
    features = file.path(path, "features.parquet"),
    metadata = file.path(path, "metadata.parquet")
  )
  
  available <- sapply(files, file.exists)
  
  # Also check for partitioned genes directory
  genes_partitioned <- dir.exists(file.path(path, "genes"))
  if (genes_partitioned) {
    files$genes <- file.path(path, "genes")
    available["genes"] <- TRUE
  }
  
  structure(
    list(
      path = path,
      files = files,
      available = available,
      is_partitioned = genes_partitioned
    ),
    class = "GTFParquet"
  )
}

#' @export
print.GTFParquet <- function(x, ...) {
  cat("GTFParquet object\n")
  cat("  Path:", x$path, "\n")
  cat("  Available tables:\n")
  for (name in names(x$available)) {
    status <- if (x$available[[name]]) "✓" else "✗"
    cat(sprintf("    %s %s\n", status, name))
  }
  
  # Show metadata if available
  if (x$available[["metadata"]]) {
    meta <- arrow::read_parquet(x$files$metadata)
    cat("  Metadata:\n")
    for (i in seq_len(nrow(meta))) {
      cat(sprintf("    %s: %s\n", meta$name[i], meta$value[i]))
    }
  }
  invisible(x)
}


# =============================================================================
# Gene Accessor - Full Attributes
# =============================================================================

#' @export
genes <- function(x, ...) {
 UseMethod("genes")
}

#' Extract genes with full GTF attributes
#'
#' Returns GRanges with rich metadata columns including gene_type, gene_name,
#' level, tags, and source - attributes lost in TxDb format.
#'
#' @param x A GTFParquet object
#' @param columns Character vector of columns to include in mcols. 
#'   Default includes all: gene_name, gene_type, level, tags, source, havana_gene
#' @param filter Named list for filtering, e.g., list(gene_type = "protein_coding")
#' @param use_versioned_ids If TRUE, use full versioned IDs (ENSG00000290825.2).
#'   If FALSE (default), use stripped IDs (ENSG00000290825).
#' @return GRanges with gene_id as names and full attribute mcols
#' @export
#' @examples
#' gtf <- GTFParquet("parquet_output")
#' 
#' # All genes with full attributes
#' gr <- genes(gtf)
#' mcols(gr)  # gene_name, gene_type, level, tags, source, havana_gene
#' 
#' # Protein-coding genes only
#' pc_genes <- genes(gtf, filter = list(gene_type = "protein_coding"))
#' 
#' # lncRNAs on chromosome 1
#' lnc_chr1 <- genes(gtf, filter = list(gene_type = "lncRNA", chrom = "chr1"))
#' 
#' # Filter by annotation level (1 = verified, 2 = manual, 3 = automatic)
#' high_conf <- genes(gtf, filter = list(level = 1))
genes.GTFParquet <- function(x, 
                              columns = NULL,
                              filter = NULL,
                              use_versioned_ids = FALSE,
                              ...) {
  if (!x$available[["genes"]]) {
    stop("genes.parquet not found in: ", x$path)
  }
  
  # Default columns - all the rich attributes
 default_cols <- c("gene_name", "gene_type", "source", "level", "tags", "havana_gene")
  if (is.null(columns)) {
    columns <- default_cols
  }
  
  # Always need these for GRanges construction
  id_col <- if (use_versioned_ids) "gene_id" else "gene_id_stripped"
  required_cols <- c(id_col, "chrom", "start", "end", "strand")
  select_cols <- unique(c(required_cols, columns))
  
  # Also include filter columns
 if (!is.null(filter)) {
    select_cols <- unique(c(select_cols, names(filter)))
  }
  
  # Open dataset (handles both partitioned and single file)
  ds <- arrow::open_dataset(x$files$genes)
  
  # Build query
  query <- ds |> dplyr::select(dplyr::any_of(select_cols))
  
  # Apply filters
  if (!is.null(filter)) {
    for (col_name in names(filter)) {
      filter_vals <- filter[[col_name]]
      query <- query |> dplyr::filter(.data[[col_name]] %in% filter_vals)
    }
  }
  
  # Collect data
 df <- query |> 
    dplyr::collect() |>
    dplyr::arrange(chrom, start)
  
  # Build GRanges
  .df_to_granges_rich(df, id_col = id_col, mcol_cols = columns)
}


# =============================================================================
# Transcript Accessor - Full Attributes
# =============================================================================

#' @export
transcripts <- function(x, ...) {
  UseMethod("transcripts")
}

#' Extract transcripts with full GTF attributes
#'
#' @param x A GTFParquet object
#' @param columns Columns to include in mcols
#' @param filter Named list for filtering
#' @param use_versioned_ids Use full versioned IDs
#' @return GRanges with transcript_id as names
#' @export
#' @examples
#' gtf <- GTFParquet("parquet_output")
#' 
#' # All transcripts
#' tx <- transcripts(gtf)
#' 
#' # Protein-coding transcripts with CCDS ID
#' ccds_tx <- transcripts(gtf, filter = list(transcript_type = "protein_coding"))
#' ccds_tx <- ccds_tx[!is.na(mcols(ccds_tx)$ccdsid)]
#' 
#' # Filter by transcript support level
#' high_support <- transcripts(gtf, filter = list(transcript_support_level = "1"))
transcripts.GTFParquet <- function(x,
                                    columns = NULL,
                                    filter = NULL,
                                    use_versioned_ids = FALSE,
                                    ...) {
  if (!x$available[["transcripts"]]) {
    stop("transcripts.parquet not found in: ", x$path)
  }
  
  default_cols <- c("transcript_name", "transcript_type", "gene_id", "gene_name",
                    "source", "level", "tags", "transcript_support_level", 
                    "havana_transcript", "ccdsid", "protein_id")
  if (is.null(columns)) {
    columns <- default_cols
  }
  
  id_col <- if (use_versioned_ids) "transcript_id" else "transcript_id_stripped"
  required_cols <- c(id_col, "chrom", "start", "end", "strand")
  select_cols <- unique(c(required_cols, columns))
  
  if (!is.null(filter)) {
    select_cols <- unique(c(select_cols, names(filter)))
  }
  
  df <- arrow::read_parquet(x$files$transcripts, col_select = dplyr::any_of(select_cols))
  
  # Apply filters
  if (!is.null(filter)) {
    for (col_name in names(filter)) {
      df <- df |> dplyr::filter(.data[[col_name]] %in% filter[[col_name]])
    }
  }
  
  df <- df |> dplyr::arrange(chrom, start)
  
  .df_to_granges_rich(df, id_col = id_col, mcol_cols = columns)
}


# =============================================================================
# Exon Accessor - Full Attributes
# =============================================================================

#' @export
exons <- function(x, ...) {
  UseMethod("exons")
}

#' Extract exons with full GTF attributes
#'
#' @param x A GTFParquet object
#' @param columns Columns to include in mcols
#' @param filter Named list for filtering
#' @param use_versioned_ids Use full versioned IDs
#' @return GRanges with exon_id as names
#' @export
exons.GTFParquet <- function(x,
                              columns = NULL,
                              filter = NULL,
                              use_versioned_ids = FALSE,
                              ...) {
  if (!x$available[["exons"]]) {
    stop("exons.parquet not found in: ", x$path)
  }
  
  default_cols <- c("exon_number", "transcript_id", "gene_id", "source", "level", "tags")
  if (is.null(columns)) {
    columns <- default_cols
  }
  
  id_col <- if (use_versioned_ids) "exon_id" else "exon_id_stripped"
  required_cols <- c(id_col, "chrom", "start", "end", "strand")
  select_cols <- unique(c(required_cols, columns))
  
  if (!is.null(filter)) {
    select_cols <- unique(c(select_cols, names(filter)))
  }
  
  df <- arrow::read_parquet(x$files$exons, col_select = dplyr::any_of(select_cols))
  
  if (!is.null(filter)) {
    for (col_name in names(filter)) {
      df <- df |> dplyr::filter(.data[[col_name]] %in% filter[[col_name]])
    }
  }
  
  df <- df |> dplyr::arrange(chrom, start)
  
  .df_to_granges_rich(df, id_col = id_col, mcol_cols = columns)
}


# =============================================================================
# CDS Accessor
# =============================================================================

#' @export
cds <- function(x, ...) {
  UseMethod("cds")
}

#' Extract CDS regions with full attributes
#'
#' @param x A GTFParquet object
#' @param columns Columns to include in mcols
#' @param filter Named list for filtering
#' @return GRanges
#' @export
cds.GTFParquet <- function(x,
                            columns = NULL,
                            filter = NULL,
                            ...) {
  if (!x$available[["cds"]]) {
    stop("cds.parquet not found in: ", x$path)
  }
  
  default_cols <- c("transcript_id", "gene_id", "protein_id", "exon_number", 
                    "frame", "source", "level", "tags")
  if (is.null(columns)) {
    columns <- default_cols
  }
  
  required_cols <- c("cds_id", "chrom", "start", "end", "strand")
  select_cols <- unique(c(required_cols, columns))
  
  if (!is.null(filter)) {
    select_cols <- unique(c(select_cols, names(filter)))
  }
  
  df <- arrow::read_parquet(x$files$cds, col_select = dplyr::any_of(select_cols))
  
  if (!is.null(filter)) {
    for (col_name in names(filter)) {
      df <- df |> dplyr::filter(.data[[col_name]] %in% filter[[col_name]])
    }
  }
  
  df <- df |> dplyr::arrange(chrom, start)
  
  .df_to_granges_rich(df, id_col = "cds_id", mcol_cols = columns)
}


# =============================================================================
# Grouped Accessors (GRangesList)
# =============================================================================

#' @export
exonsBy <- function(x, ...) {
  UseMethod("exonsBy")
}

#' Extract exons grouped by transcript or gene
#'
#' @param x A GTFParquet object
#' @param by Group by "tx" (transcript) or "gene"
#' @param filter Named list for filtering
#' @return GRangesList
#' @export
exonsBy.GTFParquet <- function(x, by = c("tx", "gene"), filter = NULL, ...) {
  by <- match.arg(by)
  
  # Get all exons with transcript/gene info
  exon_cols <- c("exon_number", "transcript_id", "transcript_id_stripped", 
                 "gene_id", "gene_id_stripped")
  exons_gr <- exons(x, columns = exon_cols, filter = filter)
  
  if (by == "tx") {
    group_col <- "transcript_id_stripped"
  } else {
    group_col <- "gene_id_stripped"
  }
  
  # Split by group
  groups <- mcols(exons_gr)[[group_col]]
  grl <- split(exons_gr, groups)
  
  # Sort exons within each transcript by exon_number
  if (by == "tx") {
    grl <- S4Vectors::endoapply(grl, function(gr) {
      gr[order(mcols(gr)$exon_number)]
    })
  }
  
  grl
}


#' @export
transcriptsBy <- function(x, ...) {
  UseMethod("transcriptsBy")
}

#' Extract transcripts grouped by gene
#'
#' @param x A GTFParquet object
#' @param by Currently only "gene" is supported
#' @param filter Named list for filtering
#' @return GRangesList with gene IDs as names
#' @export
transcriptsBy.GTFParquet <- function(x, by = "gene", filter = NULL, ...) {
  if (by != "gene") {
    stop("transcriptsBy currently only supports by='gene'
")
  }
  
  tx_gr <- transcripts(x, filter = filter)
  
  gene_ids <- mcols(tx_gr)$gene_id
  if ("gene_id_stripped" %in% colnames(mcols(tx_gr))) {
    gene_ids <- mcols(tx_gr)$gene_id_stripped
  }
  
  split(tx_gr, gene_ids)
}


#' @export
cdsBy <- function(x, ...) {
  UseMethod("cdsBy")
}

#' Extract CDS regions grouped by transcript
#'
#' @param x A GTFParquet object
#' @param by Group by "tx" (transcript) or "gene"
#' @param filter Named list for filtering
#' @return GRangesList
#' @export
cdsBy.GTFParquet <- function(x, by = c("tx", "gene"), filter = NULL, ...) {
  by <- match.arg(by)
  
  cds_gr <- cds(x, filter = filter)
  
  if (by == "tx") {
    group_col <- if ("transcript_id_stripped" %in% colnames(mcols(cds_gr))) {
      "transcript_id_stripped"
    } else {
      "transcript_id"
    }
  } else {
    group_col <- if ("gene_id_stripped" %in% colnames(mcols(cds_gr))) {
      "gene_id_stripped"
    } else {
      "gene_id"
    }
  }
  
  groups <- mcols(cds_gr)[[group_col]]
  grl <- split(cds_gr, groups)
  
  # Sort by exon_number within transcript
  if (by == "tx" && "exon_number" %in% colnames(mcols(cds_gr))) {
    grl <- S4Vectors::endoapply(grl, function(gr) {
      gr[order(mcols(gr)$exon_number)]
    })
  }
  
  grl
}


# =============================================================================
# Feature Type Queries
# =============================================================================

#' Extract UTR regions
#'
#' @param x A GTFParquet object
#' @param type "5prime", "3prime", or "both"
#' @param filter Named list for filtering
#' @return GRanges
#' @export
utrs <- function(x, type = c("both", "5prime", "3prime"), filter = NULL) {
  type <- match.arg(type)
  
  if (!x$available[["features"]]) {
    stop("features.parquet not found - no UTR data available")
  }
  
  df <- arrow::read_parquet(x$files$features)
  
  # Filter by UTR type
  utr_types <- switch(type,
    "5prime" = c("five_prime_utr", "UTR"),  # Some GTFs use generic "UTR"
    "3prime" = c("three_prime_utr"),
    "both" = c("five_prime_utr", "three_prime_utr", "UTR")
  )
  
  df <- df |> dplyr::filter(feature_type %in% utr_types)
  
  if (!is.null(filter)) {
    for (col_name in names(filter)) {
      df <- df |> dplyr::filter(.data[[col_name]] %in% filter[[col_name]])
    }
  }
  
  if (nrow(df) == 0) {
    return(GenomicRanges::GRanges())
  }
  
  gr <- GenomicRanges::GRanges(
    seqnames = df$chrom,
    ranges = IRanges::IRanges(start = df$start, end = df$end),
    strand = df$strand
  )
  
  mcols(gr) <- df[, c("feature_type", "transcript_id", "gene_id"), drop = FALSE]
  gr
}


#' Extract start/stop codons
#'
#' @param x A GTFParquet object
#' @param type "start", "stop", or "both"
#' @param filter Named list for filtering
#' @return GRanges
#' @export
codons <- function(x, type = c("both", "start", "stop"), filter = NULL) {
  type <- match.arg(type)
  
  if (!x$available[["features"]]) {
    stop("features.parquet not found")
  }
  
  df <- arrow::read_parquet(x$files$features)
  
  codon_types <- switch(type,
    "start" = "start_codon",
    "stop" = "stop_codon",
    "both" = c("start_codon", "stop_codon")
  )
  
  df <- df |> dplyr::filter(feature_type %in% codon_types)
  
  if (!is.null(filter)) {
    for (col_name in names(filter)) {
      df <- df |> dplyr::filter(.data[[col_name]] %in% filter[[col_name]])
    }
  }
  
  if (nrow(df) == 0) {
    return(GenomicRanges::GRanges())
  }
  
  gr <- GenomicRanges::GRanges(
    seqnames = df$chrom,
    ranges = IRanges::IRanges(start = df$start, end = df$end),
    strand = df$strand
  )
  
  mcols(gr) <- df[, c("feature_type", "transcript_id", "gene_id", "frame"), drop = FALSE]
  gr
}


# =============================================================================
# Convenience Query Functions
# =============================================================================

#' Get protein-coding genes
#'
#' @param x A GTFParquet object
#' @param ... Additional arguments passed to genes()
#' @return GRanges of protein-coding genes
#' @export
protein_coding_genes <- function(x, ...) {
  genes(x, filter = list(gene_type = "protein_coding"), ...)
}


#' Get lncRNA genes
#'
#' @param x A GTFParquet object
#' @param ... Additional arguments passed to genes()
#' @return GRanges of lncRNA genes
#' @export
lncRNA_genes <- function(x, ...) {
  genes(x, filter = list(gene_type = "lncRNA"), ...)
}


#' List available gene types
#'
#' @param x A GTFParquet object
#' @return Character vector of unique gene types with counts
#' @export
gene_types <- function(x) {
  if (!x$available[["genes"]]) {
    stop("genes.parquet not found")
  }
  
  df <- arrow::read_parquet(x$files$genes, col_select = "gene_type")
  table(df$gene_type, useNA = "ifany")
}


#' List available transcript types
#'
#' @param x A GTFParquet object
#' @return Table of transcript types with counts
#' @export
transcript_types <- function(x) {
  if (!x$available[["transcripts"]]) {
    stop("transcripts.parquet not found")
  }
  
  df <- arrow::read_parquet(x$files$transcripts, col_select = "transcript_type")
  table(df$transcript_type, useNA = "ifany")
}


# =============================================================================
# Region Queries
# =============================================================================

#' Find genes overlapping a genomic region
#'
#' @param x A GTFParquet object
#' @param region GRanges specifying query region(s)
#' @param ... Additional arguments passed to genes()
#' @return GRanges of genes overlapping the region
#' @export
genes_in_region <- function(x, region, ...) {
  stopifnot(inherits(region, "GRanges"))
  
  chroms <- as.character(unique(GenomicRanges::seqnames(region)))
  
  # Get genes on relevant chromosomes
  all_genes <- genes(x, filter = list(chrom = chroms), ...)
  
  # Find overlaps
  hits <- GenomicRanges::findOverlaps(all_genes, region)
  all_genes[S4Vectors::queryHits(hits)]
}


#' Find transcripts overlapping a genomic region
#'
#' @param x A GTFParquet object
#' @param region GRanges specifying query region(s)
#' @param ... Additional arguments passed to transcripts()
#' @return GRanges of transcripts overlapping the region
#' @export
transcripts_in_region <- function(x, region, ...) {
  stopifnot(inherits(region, "GRanges"))
  
  chroms <- as.character(unique(GenomicRanges::seqnames(region)))
  
  all_tx <- transcripts(x, filter = list(chrom = chroms), ...)
  
  hits <- GenomicRanges::findOverlaps(all_tx, region)
  all_tx[S4Vectors::queryHits(hits)]
}


# =============================================================================
# Helper Functions
# =============================================================================

#' Convert data frame to GRanges with rich metadata
#' @keywords internal
.df_to_granges_rich <- function(df, id_col, mcol_cols) {
  if (nrow(df) == 0) {
    return(GenomicRanges::GRanges())
  }
  
  gr <- GenomicRanges::GRanges(
    seqnames = df$chrom,
    ranges = IRanges::IRanges(start = df$start, end = df$end),
    strand = df$strand
  )
  
  # Set names from ID column
  names(gr) <- df[[id_col]]
  
  # Add metadata columns (excluding GRanges construction cols)
  exclude_cols <- c(id_col, "chrom", "start", "end", "strand")
  mcol_cols_available <- intersect(mcol_cols, colnames(df))
  mcol_cols_available <- setdiff(mcol_cols_available, exclude_cols)
  
  if (length(mcol_cols_available) > 0) {
    GenomicRanges::mcols(gr) <- df[, mcol_cols_available, drop = FALSE]
  }
  
  gr
}


#' Get GTF metadata
#'
#' @param x A GTFParquet object
#' @return Named character vector of metadata
#' @export
gtf_metadata <- function(x) {
  if (!x$available[["metadata"]]) {
    return(character(0))
  }
  
  meta <- arrow::read_parquet(x$files$metadata)
  setNames(meta$value, meta$name)
}
