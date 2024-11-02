#' Convert SnapGene file to Biostrings DNAStringSet
#' @param filepath Path to SnapGene .dna file
#' @return Biostrings DNAStringSet object with features and annotations
snapgene_file_to_dnastringset <- function(filepath) {
  data <- snapgene_file_to_list(filepath)

  # Create DNAStringSet
  seq <- Biostrings::DNAStringSet(data$seq)
  names(seq) <- "sequence"

  # Add features as metadata
  metadata <- list(
    topology = data$dna$topology,
    features = data$features,
    notes = data$notes
  )

  GenomicRanges::mcols(seq)$metadata <- metadata

  return(seq)
}
