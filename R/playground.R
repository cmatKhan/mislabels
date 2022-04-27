library(VariantAnnotation)

file.gz <- "somaticVcfBeta-HCC1187-H-200-37-ASM-T1-N1.vcf.gz"

stopifnot(file.exists(file.gz))

file.gz.tbi <- paste(file.gz, ".tbi", sep="")

if(!(file.exists(file.gz.tbi))) {
    indexTabix(file.gz, format="vcf")
}


start.loc <- 55000000
end.loc <- 56000000
chr7.gr <- GRanges("7", IRanges(start.loc, end.loc))


params <- ScanVcfParam(which=chr7.gr)


vcf <- readVcf(TabixFile(file.gz), "hg38", params)


writeVcf(vcf, "chr7-sub.vcf")


bgzip("chr7-sub.vcf", overwrite=TRUE)
indexTabix("chr7-sub.vcf.gz", format="vcf")
