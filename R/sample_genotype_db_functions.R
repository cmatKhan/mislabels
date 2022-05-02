#'
#' Create a sqlite database to store sample/genotype data
#'
#' @import RSQLite
#'
#' @param db_path a path, eg /path/to/sample_geno_db.sqlite, to where you'd like
#'   to create the sqlite database file
#'
#' @export
create_sample_geno_db = function(db_path){
  con = dbConnect(RSQLite::SQLite(), db_path)
  sql = paste0("CREATE TABLE sample_genotype(sample TEXT, chr TEXT, ",
               "bp INTEGER, ref_alt TEXT, genotype TEXT);")
  dbExecute(con, sql)
  dbDisconnect(con)
}

#'
#' Create a sqlite database to store sample/genotype data
#'
#' @import RSQLite
#'
#' @param db_path a path, eg /path/to/sample_geno_db.sqlite, to where you'd like
#'   to create the sqlite database file
#'
#' @export
create_bam_pileup_table = function(db_path){
  con = dbConnect(RSQLite::SQLite(), db_path)
  sql = paste0("CREATE TABLE rnaseq_pileup(sample TEXT, visit TEXT, chr TEXT, ",
               "pos INTEGER, nucleotide TEXT, count INTEGER);")
  dbExecute(con, sql)
  dbDisconnect(con)
}

#'
#' Fill the sample_genotype database from VCF file(s)
#'
#' @import RSQLite
#' @importFrom dplyr as_tibble mutate select
#' @importFrom tidyr pivot_longer
#' @importFrom stringr str_extract
#'
#' @inheritParams create_sample_geno_db
#'
#' @param vcf_path path to a vcf file (must already be indexed)
#'
#' @export
add_genotype_table = function(vcf_path, db_path){

  db_conn = dbConnect(RSQLite::SQLite(), db_path)

  vcf_file = VcfFile(file      = vcf_path,
                     index     = paste0(vcf_path, ".tbi"),
                     yieldSize = 100)

  vcf = readVcf(vcf_file)

  geno_mat = genotypeToSnpMatrix(vcf, uncertain=TRUE)

  geno_df = t(as(geno_mat$genotype, "character")) %>%
    as_tibble(rownames = "location") %>%
    pivot_longer(-location,
                 names_to = "sample",
                 values_to = "genotype") %>%
    mutate(chr     = as.numeric(str_extract(location, "(?<=^chr)\\d+")),
           bp      = as.numeric(str_extract(location, "(?<=:)\\d+")),
           ref_alt = str_extract(location, "\\w\\/\\w$")) %>%
    dplyr::select(sample,chr,bp,ref_alt,genotype)

  dbAppendTable(db_conn, "sample_genotype_table", geno_df)

  dbDisconnect(db_conn)

}

#'
#' create virtual table of only highly expressed gene regions in the db
#'
#' @import RSQLite
#'
#' @inheritParams create_sample_geno_db
#' @param highly_expressed_gene_ranges a GRangeList of regions over which to
#'   create a virtual table (a subset of the entire table)
#' @param sample_viewname name of the sample view table in the db
#' @param tablename name of the table from which to make the sample view
#'
#'
#'@export
highly_expressed_gene_virtual_table = function(highly_expressed_gene_ranges,
                                               db_path,
                                               sample_viewname,
                                               tablename = "sample_genotype"){

  db_conn = dbConnect(RSQLite::SQLite(), db_path)

  roi_df = as_tibble(highly_expressed_gene_ranges[highly_expressed_gene_ranges$type=="gene"]) %>%
    dplyr::select(seqnames, start, end, strand) %>%
    sample_n(2000) %>%
    filter(seqnames %in% paste0("chr1", seq(1,22)))
  # %>%
  #   filter(strand == "+", seqnames %in% paste0("chr1", seq(1,22)))

  roi_db_string = function(roi_row){

    chromosome = paste0(str_remove(roi_row[['seqnames']], "chr"), ".0")

    sprintf("chr IS %s AND bp BETWEEN %s AND %s", chromosome,
            roi_row[['start']],
            roi_row[['end']])
  }

  roi_list = apply(roi_df, 1, roi_db_string)

  sql = sprintf('DROP VIEW IF EXISTS "main"."%s"', sample_viewname)
  dbExecute(db_conn, sql)

  sql = sprintf('CREATE VIEW "%s" AS SELECT * FROM %s WHERE (%s)',
                sample_viewname,
                tablename,
                paste(roi_list, collapse = " or "))
  dbExecute(db_conn, sql)

  dbDisconnect(db_conn)
}
