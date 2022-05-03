#'
#' extract unique SNP locations from a sample (view in the sample_genome_db)
#' as a dataframe
#'
#' @import GenomicRanges iterators RSQLite
#'
#' @inheritParams create_bam_pileup_table
#' @param tablename the table from which to extract the SNPS
#'
#' @export
sample_snp_ranges = function(db_path, tablename){

  con = dbConnect(RSQLite::SQLite(), db_path)

  sql = sprintf("SELECT DISTINCT chr, bp FROM %s", tablename)

  df = dbGetQuery(con, sql) %>%
    as_tibble()

  dbDisconnect(con)

  df
}

#'
#' extract unique SNP locations from RNAseq library
#'
#' @import GenomicRanges iterators RSQLite
#' @importFrom dplyr rename
#'
#' @inheritParams create_bam_pileup_table
#' @param tablename name of the rnaseq table. default is rnaseq_pileup
#'
#' @export
sample_snp_loc_from_pileup = function(db_path, tablename = 'rnaseq_pileup'){

  con = dbConnect(RSQLite::SQLite(), db_path)

  sql = paste0(sprintf("SELECT * FROM %s JOIN
    (SELECT rowid as rid FROM %s
        WHERE random()", tablename, tablename), " % 10 = 0  -- Reduce rowids by 10x\n",
        sprintf("ORDER BY random() LIMIT 1000) AS srid
    ON %s.rowid = srid.rid;", tablename))

  # sql = sprintf("SELECT DISTINCT chr, pos FROM %s", tablename)
  # sql

  df = dbGetQuery(con, sql) %>%
    as_tibble() %>%
    dplyr::select(chr, pos) %>%
    dplyr::rename(bp = pos) %>%
    mutate(chr = as.numeric(chr)) %>%
    mutate(chr = paste0(as.character(chr), ".0"))

  dbDisconnect(con)

  df
}

#'
#' generate Snp vectors from each genome
#' @inheritParams create_bam_pileup_table
#' @param tablename the table from which to extract the SNPS
#'
#' @return a list, length number of individuals in genome database, of snp
#' vectors
#'
#' @export
genome_vectors = function(db_path,tablename){

  con = dbConnect(RSQLite::SQLite(), db_path)

  sql = sprintf("SELECT * FROM %s", tablename)

  df_list = dbGetQuery(con, sql) %>%
    as_tibble()

  dbDisconnect(con)

  list_names = unique(df_list$sample)

  split_geno = str_split(df_list$ref_alt, "/")

  df_list = df_list %>%
    mutate(ref_genotype = unlist(map(seq(1,nrow(df_list)),
                                     ~split_geno[[.]][[1]])),
           alt_genotype = unlist(map(seq(1,nrow(df_list)),
                                     ~split_geno[[.]][[2]]))) %>%
    mutate(nucleotide = ifelse(genotype == "A/A", ref_genotype, "")) %>%
    mutate(nucleotide = ifelse(genotype == "B/B", alt_genotype, nucleotide)) %>%
    mutate(nucleotide = ifelse(genotype == "A/B", ref_alt, nucleotide)) %>%
    mutate(nucleotide = ifelse(genotype == "Uncertain", genotype, nucleotide)) %>%
    group_by(sample) %>%
    arrange(chr,bp, .by_group = TRUE) %>%
    group_split()

  names(df_list) = list_names

  out = map(names(df_list), ~pull(df_list[[.]], nucleotide))

  names(out) = list_names

  out
}

#'
#' get vector of consensus bases over a given set of SNPs in an RNAseq library
#'
#' @param db_path path to an RNAseq library pileup sqlite database
#' @param snp_df path to the snp dataframe (see function)
#' @param tablename the table from which to extract the SNPS
#'
#' @return RNAseq library Snp vector
#'
#'@export
rnaseq_vector = function(db_path, snp_df, tablename = "rnaseq_pileup"){

  snp_df = snp_df %>%
    mutate(chr = str_remove(as.character(chr),".0"))

  snp_where_clause = function(snp_row){

    chromosome = snp_row[['chr']]
    loc = snp_row[['bp']]

    sprintf("chr IS %s AND pos IS %s", chromosome, loc)
  }

  where_clause = paste(apply(snp_df,1,snp_where_clause), collapse = " or ")

  sql = sprintf("SELECT * FROM %s WHERE %s", tablename, where_clause)

  message(sprintf("getting pileup vector from: %s", db_path))
  con = dbConnect(RSQLite::SQLite(), db_path)

  df = dbGetQuery(con, sql)

  dbDisconnect(con)

  snp_df %>%
    left_join(df, by = c('chr' = 'chr', 'bp' = 'pos')) %>%
    group_by(chr, bp) %>%
    summarize(geno = paste(nucleotide, collapse = "/"), .groups = 'keep') %>%
    mutate(geno = str_replace(geno, "NA", "Uncertain")) %>%
    pull(geno)

  # list(db_df = df, join_df = out, vec = pull(out, geno))



}
