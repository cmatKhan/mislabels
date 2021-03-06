% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sample_genotype_db_functions.R
\name{add_genotype_table}
\alias{add_genotype_table}
\title{Fill the sample_genotype database from VCF file(s)}
\usage{
add_genotype_table(vcf_path, db_path)
}
\arguments{
\item{vcf_path}{path to a vcf file (must already be indexed)}

\item{db_path}{a path, eg /path/to/sample_geno_db.sqlite, to where you'd like
to create the sqlite database file}
}
\description{
Fill the sample_genotype database from VCF file(s)
}
