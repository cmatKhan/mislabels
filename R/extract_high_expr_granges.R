library(SNPRelate)

#high_expr_genes_ranges = readRDS("data/high_expr_genes.rds")

plink_files = list(
  bed = "/mnt/scratch/llfs/plink_files/freeze5.allchr.allfam.bed",
  fam = "/mnt/scratch/llfs/plink_files/freeze5.allchr.allfam.fam",
  bim = "/mnt/scratch/llfs/plink_files/freeze5.allchr.allfam.bim"
)



snpgdsBED2GDS(plink_files$bed, plink_files$fam, plink_files$bim, "data/allfam.gds")
