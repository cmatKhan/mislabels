#'
#' Create a plink region filter from GRangeList
#'
#' @import GenomicFeatures
#' @importFrom tibble as_tibble
#' @importFrom dplyr filter select
#'
#' @param grl a GRangeList object. must have type with at least factor level
#'   gene and seqnames, start, end and gene_id
#'
#' @return a dataframe which may be written out to act as a plink filter
#'
#' @export
create_plink_region_filter = function(grl){
  grl %>%
    as_tibble() %>%
    filter(type == "gene") %>%
    dplyr::select(seqnames, start, end, gene_id) %>%
    mutate(seqnames = str_remove(seqnames, "chr"))
}

# hamming dist -- exact match to genome. RNAseq reads for heterozygous should
# be ~50/50.
# Question: do we expect there to be 50/50, or maybe a 3rd, to be expressed
# on both allele?
