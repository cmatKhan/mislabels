#'
#' Create variants of the distance matrix
#' @description distance between two SNPS vectors are over non NA positions
#' where if lenient is FALSE (default), genotype in genome and rnaseq library
#' pileup must match exactly for a distance score of 0. Note that A/T == T/A.
#' If lenient == TRUE, then if the nucleotide is contained in the genotype,
#' then the distance = 0. For example, A == A, A/T, T/A, A/G, G/A, A/C, C/A.
#'
#' @importFrom stringr str_split str_detect
#' @importFrom purrr map
#'
#' @param lenient Set to TRUE to consider, for example,
#'   A == A, A/T, T/A, A/G, G/A, A/C, C/A for a distance score of 0. Default is
#'   FALSE
#
#' @return a symmetric matrix with rownames corresponding to the 16 possible
#'   genotypes
#'
#' @export
equivalence_class_dist_matrix = function(lenient = FALSE){

  margin_names = c('A','T','G','C', 'A/T', 'T/A', 'A/G', 'G/A', 'A/C','C/A',
                   'T/G', 'G/T', 'T/C', 'C/T', 'G/C', 'C/G')

  if(lenient){
    mat = matrix(
      as.numeric(
        unlist(
          map(
            str_split(margin_names,
                      "/",
                      simplify = TRUE)[,1],
            ~str_detect(margin_names,., negate = TRUE)))),
      nrow = 16, ncol = 16)
  } else{
    mat = matrix(nrow=length(margin_names), ncol =length(margin_names))
    diag(mat) = 0
    mat[is.na(mat)] = 1
  }

  rownames(mat) = margin_names
  colnames(mat) = margin_names

  if(!lenient){
    mat['A/T','T/A'] = 0
    mat['T/A', 'A/T'] = 0

    mat['A/G','G/A'] = 0
    mat['G/A', 'A/G'] = 0

    mat['A/C','C/A'] = 0
    mat['C/A', 'A/C'] = 0

    mat['T/G','G/T'] = 0
    mat['G/T', 'T/G'] = 0

    mat['T/C','C/T'] = 0
    mat['C/T', 'T/C'] = 0

    mat['G/C','C/G'] = 0
    mat['C/G', 'G/C'] = 0

  }

  mat = rbind(mat,rep(NA,16))
  rownames(mat)[17] = "Uncertain"
  mat = cbind(mat, rep(NA,17))
  colnames(mat)[17] = "Uncertain"

  mat

}

#'
#' compare two SNP vectors
#' @param genome_vector a Snp vector of an individual genome
#' @param pileup_vector an Rnaseq Snp vector
#' @param equivalence_mat the scoring matrix which defines what is different
#' from what. Eg, Is A/T == T/A (should be). What about A and A/T?
#'
#' @return fraction of matches over total number of genotyped bases in both
#' vectors
#'
#' @export
compare_vectors = function(genome_vector, pileup_vector, equivalence_mat){

  # dist_vec will look something like this
  # > dist_vec
  # [1] NA NA NA NA  1  1 NA NA  1  1  1  0  1  0 NA NA NA NA NA
  # [20] NA NA NA NA NA  1  1  1  1  1 NA NA NA NA NA NA NA  1  1
  # [39]  1  0  1  1  0 NA NA  1  1  1 NA NA NA NA NA NA NA NA NA
  dist_vec = unlist(map2(genome_vector,
                         pileup_vector,
                         ~equivalence_mat[.x,.y]))

  # ratio of difference position / total positions
  sum(dist_vec[!is.na(dist_vec)]) / length(dist_vec[!is.na(dist_vec)])
}

#'
#' Conduct an experiment to attempt to map RNAseq librarys to genomes
#' @param db_path path to the sample_genotype_db
#' @param pileup_database_dir directory which stores the mislabelled and control
#'   pileup databases
#' @param sample the name of the sample (a view in the sample_geontype database)
#' @param filename_map a dataframe which maps pileup database filenames to
#' id and subject numbers
#' @return the results of the experiment
#'
#' @export
conduct_experiment = function(db_path, pileup_database_dir, sample, filename_map){

  snp_df = sample_snp_ranges(db_path, sample)

  sample_genome_vectors = genome_vectors(db_path,
                                         sample)

  get_diff = function(pileup_db_basename){

    pileup_path = file.path(pileup_database_dir,
                            pileup_db_basename[['filename']])

    pileup_vector = rnaseq_vector(pileup_path, snp_df)

    dist_mat = equivalence_class_dist_matrix(lenient = FALSE)

    out_list = map(sample_genome_vectors,
                   compare_vectors,
                   pileup_vector,
                   dist_mat)

    names(out_list) = names(sample_genome_vectors)

    out_list
  }

  apply(filename_map, 1, get_diff)
}

# run an experiment example:
#   gene_ranges = readRDS(high_expr_genes)
# map(sample_list, ~highly_expressed_gene_virtual_table(gene_ranges,
#                                                       here("data/sample_genotype_db.sqlite"),
#                                                       .))
# sample_list = paste0("sample", seq(1,4))
#
# filename_map = read_csv(here("data/filename_id_subject_mislabel_map.csv")) %>%
#   mutate(filename =
#            paste0(
#              str_remove_all(filename, "_visit|.markdup.sorted.pileup"),
#              '.sqlite')) %>%
#   mutate(filename = str_replace(filename, "_", "."))
#
#
# experiment_res = map(sample_list, conduct_experiment, filename_map)

# sample_genotype_db_path = "/home/oguzkhan/class/bioseq/local_data/sample_genotype_db.sqlite"
#
# sample_list = paste0("sample", seq(1,4))
#
# experiment_res =
#   map(sample_list,
#       ~conduct_experiment(sample_genotype_db_path,
#                           "/home/oguzkhan/class/bioseq/local_data/mapqOver10_no_secondary_databases",
#                           .,
#                           filename_id_subject_mislabel_map))
