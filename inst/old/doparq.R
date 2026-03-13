#' TxDb Parquet Interface for GenomicRanges
#'
#' Functions to work with TxDb data stored in Parquet format,
#' emulating GenomicFeatures accessor methods.
#'

# =============================================================================
# Setup and Connection
# =============================================================================
#' Create a TxDbParquet object pointing to a Parquet directory
#' @import GenomicRanges
#' @import arrow
#' @rawNamespace import(dplyr, except=c(setdiff, intersect, union))
#' @rdname txdb_methods
#' @param path Path to the Parquet output directory (containing normalized/ and/or denormalized/)
#' @return A TxDbParquet S3 object
#' @export
#' @examples
#' txdb <- TxDbParquet("path/to/parquet_output")
#' genes(txdb)
TxDbParquet <- function(path) {
  path <- normalizePath(path, mustWork = TRUE)
  

  structure(
    list(
      path = path,
      normalized_path = file.path(path, "normalized"),
      denormalized_path = file.path(path, "denormalized", "transcript_annotations"),
      has_normalized = dir.exists(file.path(path, "normalized")),
      has_denormalized = dir.exists(file.path(path, "denormalized"))
    ),
    class = "TxDbParquet"
  )
}

#' @export
print.TxDbParquet <- function(x, ...) {
  cat("TxDbParquet object\n")
  cat("  Path:", x$path, "\n")
  cat("  Normalized tables:", if (x$has_normalized) "available" else "not found", "\n")
  cat("  Denormalized view:", if (x$has_denormalized) "available" else "not found", "\n")
  invisible(x)
}


# =============================================================================
# Core Gene Accessor
# =============================================================================

#' Extract gene ranges from TxDbParquet
#'
#' Emulates GenomicFeatures::genes() - returns a GRanges object with one
#' range per gene, spanning from the minimum transcript start to maximum
#' transcript end for each gene.
#'
#' @param x A TxDbParquet object
#' @param columns Additional columns to include as mcols (e.g., "tx_type")
#' @param filter Optional named list for filtering (e.g., list(chrom = "chr1"))
#' @return A GRanges object with gene_id as names
#' @export
#' @examples
#' txdb <- TxDbParquet("path/to/parquet")
#' gr <- genes(txdb)
#' gr <- genes(txdb, filter = list(chrom = "chr1"))
genes.TxDbParquet <- function(x, columns = NULL, filter = NULL) {
  if (x$has_denormalized) {
    .genes_from_denormalized(x, columns, filter)
  } else if (x$has_normalized) {
    .genes_from_normalized(x, columns, filter)
  } else {
    stop("No Parquet data found at: ", x$path)
  }
}

#' @export
genes <- function(x, ...) {
  UseMethod("genes")
}


# Internal: Build genes GRanges from denormalized Parquet
.genes_from_denormalized <- function(x, columns = NULL, filter = NULL) {
  # Open the dataset (handles partitioned data)
  ds <- arrow::open_dataset(x$denormalized_path)
  
  # Select columns needed for gene ranges
  select_cols <- c("gene_id", "chrom", "strand", "tx_start", "tx_end")
  if (!is.null(columns)) {
    select_cols <- unique(c(select_cols, columns))
  }
  
  # Build query
  query <- ds |> dplyr::select(dplyr::any_of(select_cols))
  
  # Apply filters if provided
  if (!is.null(filter)) {
    for (col_name in names(filter)) {
      filter_val <- filter[[col_name]]
      query <- query |> dplyr::filter(.data[[col_name]] %in% filter_val)
    }
  }
  
  # Collect and aggregate by gene
  df <- query |>
    dplyr::collect() |>
    dplyr::filter(!is.na(gene_id)) |>
    dplyr::group_by(gene_id, chrom, strand) |>
    dplyr::summarise(
      start = min(tx_start),
      end = max(tx_end),
      n_transcripts = dplyr::n(),
      .groups = "drop"
    )
  
  .df_to_granges(df, extra_mcols = "n_transcripts")
}


# Internal: Build genes GRanges from normalized Parquet tables
.genes_from_normalized <- function(x, columns = NULL, filter = NULL) {
  # Read gene and transcript tables
  gene_tbl <- arrow::read_parquet(file.path(x$normalized_path, "gene.parquet"))
  tx_tbl <- arrow::read_parquet(file.path(x$normalized_path, "transcript.parquet"))
  
  # Join gene with transcript
  df <- gene_tbl |>
    dplyr::inner_join(tx_tbl, by = "_tx_id")
  
  # Apply filters if provided
  if (!is.null(filter)) {
    if ("chrom" %in% names(filter)) {
      df <- df |> dplyr::filter(tx_chrom %in% filter$chrom)
    }
    if ("strand" %in% names(filter)) {
      df <- df |> dplyr::filter(tx_strand %in% filter$strand)
    }
  }
  
  # Aggregate by gene
  df <- df |>
    dplyr::group_by(gene_id, tx_chrom, tx_strand) |>
    dplyr::summarise(
      start = min(tx_start),
      end = max(tx_end),
      n_transcripts = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::rename(chrom = tx_chrom, strand = tx_strand)
  
  .df_to_granges(df, extra_mcols = "n_transcripts")
}


# =============================================================================
# Helper Functions
# =============================================================================

# Convert a data frame with chrom, start, end, strand, gene_id to GRanges
.df_to_granges <- function(df, extra_mcols = NULL) {
  # Handle empty result

  if (nrow(df) == 0) {
    return(GenomicRanges::GRanges())
  }
  
  # Build GRanges
 gr <- GenomicRanges::GRanges(
    seqnames = df$chrom,
    ranges = IRanges::IRanges(start = df$start, end = df$end),
    strand = df$strand
  )
  
  # Set gene_id as names
  names(gr) <- df$gene_id
  
  # Add metadata columns
  mcols_df <- df[, setdiff(names(df), c("chrom", "start", "end", "strand", "gene_id")), drop = FALSE]
  if (ncol(mcols_df) > 0) {
    GenomicRanges::mcols(gr) <- mcols_df
  }
  
  gr
}


# =============================================================================
# Additional Accessors (Transcripts, Exons, CDS)
# =============================================================================

#' Extract transcript ranges from TxDbParquet
#'
#' @param x A TxDbParquet object
#' @param columns Additional columns to include
#' @param filter Optional named list for filtering
#' @return A GRanges object with tx_id as names
#' @export
transcripts <- function(x, ...) {
  UseMethod("transcripts")
}

#' @export
transcripts.TxDbParquet <- function(x, columns = NULL, filter = NULL) {
  if (x$has_denormalized) {
    ds <- arrow::open_dataset(x$denormalized_path)
    
    select_cols <- c("tx_id", "tx_name", "tx_type", "gene_id", "chrom", "strand", "tx_start", "tx_end")
    if (!is.null(columns)) {
      select_cols <- unique(c(select_cols, columns))
    }
    
    query <- ds |> dplyr::select(dplyr::any_of(select_cols))
    
    if (!is.null(filter)) {
      for (col_name in names(filter)) {
        query <- query |> dplyr::filter(.data[[col_name]] %in% filter[[col_name]])
      }
    }
    
    df <- query |> dplyr::collect()
    
  } else {
    tx_tbl <- arrow::read_parquet(file.path(x$normalized_path, "transcript.parquet"))
    gene_tbl <- arrow::read_parquet(file.path(x$normalized_path, "gene.parquet"))
    
    df <- tx_tbl |>
      dplyr::left_join(gene_tbl, by = "_tx_id") |>
      dplyr::rename(tx_id = `_tx_id`, chrom = tx_chrom, strand = tx_strand)
    
    if (!is.null(filter)) {
      for (col_name in names(filter)) {
        df <- df |> dplyr::filter(.data[[col_name]] %in% filter[[col_name]])
      }
    }
  }
  
  gr <- GenomicRanges::GRanges(
    seqnames = df$chrom,
    ranges = IRanges::IRanges(start = df$tx_start, end = df$tx_end),
    strand = df$strand
  )
  
  names(gr) <- df$tx_id
  
  mcols_df <- df[, intersect(c("tx_name", "tx_type", "gene_id"), names(df)), drop = FALSE]
  if (ncol(mcols_df) > 0) {
    GenomicRanges::mcols(gr) <- mcols_df
  }
  
  gr
}


#' Extract exon ranges from TxDbParquet
#'
#' @param x A TxDbParquet object
#' @param columns Additional columns to include
#' @param filter Optional named list for filtering
#' @return A GRanges object
#' @export
exons <- function(x, ...) {
  UseMethod("exons")
}

#' @export
exons.TxDbParquet <- function(x, columns = NULL, filter = NULL) {
  if (!x$has_normalized) {
    stop("exons() requires normalized Parquet tables")
  }
  
  exon_tbl <- arrow::read_parquet(file.path(x$normalized_path, "exon.parquet"))
  
  df <- exon_tbl |>
    dplyr::rename(chrom = exon_chrom, strand = exon_strand)
  
  if (!is.null(filter)) {
    for (col_name in names(filter)) {
      df <- df |> dplyr::filter(.data[[col_name]] %in% filter[[col_name]])
    }
  }
  
  gr <- GenomicRanges::GRanges(
    seqnames = df$chrom,
    ranges = IRanges::IRanges(start = df$exon_start, end = df$exon_end),
    strand = df$strand
  )
  
  names(gr) <- df$`_exon_id`
  GenomicRanges::mcols(gr)$exon_name <- df$exon_name
  
  gr
}


#' Extract CDS ranges from TxDbParquet
#'
#' @param x A TxDbParquet object
#' @param columns Additional columns to include
#' @param filter Optional named list for filtering
#' @return A GRanges object
#' @export
cds <- function(x, ...) {
  UseMethod("cds")
}

#' @export
cds.TxDbParquet <- function(x, columns = NULL, filter = NULL) {
  if (!x$has_normalized) {
    stop("cds() requires normalized Parquet tables")
  }
  
  cds_tbl <- arrow::read_parquet(file.path(x$normalized_path, "cds.parquet"))
  
  df <- cds_tbl |>
    dplyr::rename(chrom = cds_chrom, strand = cds_strand)
  
  if (!is.null(filter)) {
    for (col_name in names(filter)) {
      df <- df |> dplyr::filter(.data[[col_name]] %in% filter[[col_name]])
    }
  }
  
  gr <- GenomicRanges::GRanges(
    seqnames = df$chrom,
    ranges = IRanges::IRanges(start = df$cds_start, end = df$cds_end),
    strand = df$strand
  )
  
  names(gr) <- df$`_cds_id`
  GenomicRanges::mcols(gr)$cds_name <- df$cds_name
  
  gr
}


# =============================================================================
# Exons/CDS by Transcript (GRangesList)
# =============================================================================

#' Extract exons grouped by transcript
#'
#' @param x A TxDbParquet object
#' @param by Group by "tx" (transcript) or "gene"
#' @param filter Optional named list for filtering
#' @return A GRangesList, one element per transcript or gene
#' @export
exonsBy <- function(x, ...) {
  UseMethod("exonsBy")
}

#' @export
exonsBy.TxDbParquet <- function(x, by = c("tx", "gene"), filter = NULL) {
  by <- match.arg(by)
  
  if (by == "tx") {
    .exons_by_tx(x, filter)
  } else {
    .exons_by_gene(x, filter)
  }
}


.exons_by_tx <- function(x, filter = NULL) {
  if (x$has_denormalized) {
    # Use denormalized view - exons are already nested
    ds <- arrow::open_dataset(x$denormalized_path)
    
    query <- ds |> dplyr::select(tx_id, tx_name, chrom, strand, exons)
    
    if (!is.null(filter)) {
      for (col_name in names(filter)) {
        query <- query |> dplyr::filter(.data[[col_name]] %in% filter[[col_name]])
      }
    }
    
    df <- query |> dplyr::collect()
    
    # Build GRangesList from nested exon data
    grl <- lapply(seq_len(nrow(df)), function(i) {
      row <- df[i, ]
      exon_list <- row$exons[[1]]
      
      if (length(exon_list) == 0) {
        return(GenomicRanges::GRanges())
      }
      
      exon_df <- dplyr::bind_rows(exon_list)
      
      gr <- GenomicRanges::GRanges(
        seqnames = row$chrom,
        ranges = IRanges::IRanges(start = exon_df$exon_start, end = exon_df$exon_end),
        strand = row$strand
      )
      GenomicRanges::mcols(gr)$exon_id <- exon_df$exon_id
      GenomicRanges::mcols(gr)$exon_rank <- exon_df$exon_rank
      GenomicRanges::mcols(gr)$exon_name <- exon_df$exon_name
      
      gr
    })
    
    names(grl) <- df$tx_id
    GenomicRanges::GRangesList(grl)
    
  } else {
    # Use normalized tables
    splicing <- arrow::read_parquet(file.path(x$normalized_path, "splicing.parquet"))
    exon_tbl <- arrow::read_parquet(file.path(x$normalized_path, "exon.parquet"))
    tx_tbl <- arrow::read_parquet(file.path(x$normalized_path, "transcript.parquet"))
    
    # Join and aggregate
    df <- splicing |>
      dplyr::inner_join(exon_tbl, by = "_exon_id") |>
      dplyr::inner_join(
        tx_tbl |> dplyr::select(`_tx_id`, tx_chrom, tx_strand),
        by = "_tx_id"
      ) |>
      dplyr::arrange(`_tx_id`, exon_rank)
    
    if (!is.null(filter)) {
      if ("chrom" %in% names(filter)) {
        df <- df |> dplyr::filter(tx_chrom %in% filter$chrom)
      }
    }
    
    # Split by transcript
    split_df <- split(df, df$`_tx_id`)
    
    grl <- lapply(split_df, function(tx_df) {
      gr <- GenomicRanges::GRanges(
        seqnames = tx_df$exon_chrom,
        ranges = IRanges::IRanges(start = tx_df$exon_start, end = tx_df$exon_end),
        strand = tx_df$exon_strand
      )
      GenomicRanges::mcols(gr)$exon_id <- tx_df$`_exon_id`
      GenomicRanges::mcols(gr)$exon_rank <- tx_df$exon_rank
      GenomicRanges::mcols(gr)$exon_name <- tx_df$exon_name
      gr
    })
    
    GenomicRanges::GRangesList(grl)
  }
}


.exons_by_gene <- function(x, filter = NULL) {
  # Get exons by transcript first, then combine by gene
  exons_by_tx <- .exons_by_tx(x, filter)
  
  # Get tx -> gene mapping
  if (x$has_denormalized) {
    ds <- arrow::open_dataset(x$denormalized_path)
    tx_gene <- ds |>
      dplyr::select(tx_id, gene_id) |>
      dplyr::collect() |>
      dplyr::filter(!is.na(gene_id))
  } else {
    gene_tbl <- arrow::read_parquet(file.path(x$normalized_path, "gene.parquet"))
    tx_gene <- gene_tbl |>
      dplyr::rename(tx_id = `_tx_id`)
  }
  
  # Map tx_ids to gene_ids
  tx_to_gene <- setNames(tx_gene$gene_id, as.character(tx_gene$tx_id))
  
  # Group by gene
  gene_ids <- tx_to_gene[names(exons_by_tx)]
  gene_ids[is.na(gene_ids)] <- "unknown"
  
  # Combine exons by gene and reduce to unique ranges
  by_gene <- split(exons_by_tx, gene_ids)
  
  grl <- lapply(by_gene, function(gene_exons) {
    combined <- unlist(gene_exons)
    GenomicRanges::reduce(combined)
  })
  
  GenomicRanges::GRangesList(grl)
}


# =============================================================================
# Seqinfo
# =============================================================================

#' Get sequence (chromosome) information
#'
#' @param x A TxDbParquet object
#' @return A Seqinfo object
#' @export
#' @note If GenomeInfoDb is loaded, you can use GenomeInfoDb::seqinfo(txdb).
#'       Otherwise use seqinfo.TxDbParquet(txdb) directly.
seqinfo.TxDbParquet <- function(x, ...) {
  if (x$has_normalized) {
    chrominfo <- arrow::read_parquet(file.path(x$normalized_path, "chrominfo.parquet"))
  } else {
    # Extract from denormalized
    ds <- arrow::open_dataset(x$denormalized_path)
    chrominfo <- ds |>
      dplyr::select(chrom, chrom_length, is_circular) |>
      dplyr::distinct() |>
      dplyr::collect()
  }
  
  GenomeInfoDb::Seqinfo(
    seqnames = chrominfo$chrom,
    seqlengths = chrominfo$chrom_length,
    isCircular = as.logical(chrominfo$is_circular)
  )
}


# =============================================================================
# Range Queries
# =============================================================================

#' Find genes overlapping a genomic region
#'
#' Efficient region query using partition pruning when possible.
#'
#' @param x A TxDbParquet object
#' @param region A GRanges object specifying the query region(s)
#' @return A GRanges object of overlapping genes
#' @export
genes_in_region <- function(x, region) {
  stopifnot(inherits(region, "GRanges"))
  
  # Extract region info
  chrom <- as.character(GenomicRanges::seqnames(region))
  region_start <- GenomicRanges::start(region)
  region_end <- GenomicRanges::end(region)
  
  if (x$has_denormalized) {
    ds <- arrow::open_dataset(x$denormalized_path)
    
    # Use partition pruning on chrom, then filter by position
    query <- ds |>
      dplyr::filter(
        chrom %in% !!unique(chrom),
        tx_start <= !!max(region_end),
        tx_end >= !!min(region_start)
      ) |>
      dplyr::select(gene_id, chrom, strand, tx_start, tx_end) |>
      dplyr::collect() |>
      dplyr::filter(!is.na(gene_id))
    
    # Aggregate to gene level
    df <- query |>
      dplyr::group_by(gene_id, chrom, strand) |>
      dplyr::summarise(
        start = min(tx_start),
        end = max(tx_end),
        .groups = "drop"
      )
    
    gene_gr <- .df_to_granges(df)
    
    # Find actual overlaps
    hits <- GenomicRanges::findOverlaps(gene_gr, region)
    gene_gr[S4Vectors::queryHits(hits)]
    
  } else {
    # Fall back to loading all genes and filtering
    all_genes <- genes(x, filter = list(chrom = unique(chrom)))
    hits <- GenomicRanges::findOverlaps(all_genes, region)
    all_genes[S4Vectors::queryHits(hits)]
  }
}

