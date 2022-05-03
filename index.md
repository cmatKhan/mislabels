---
title: "mislabels"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{mislabels}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = FALSE
)
```

```{r setup, eval = TRUE}
library(finalProject)
library(patchwork)
library(gridExtra)
library(tidyverse)
```

# Introduction

The Long Life Family Study is a genetics/genomics study of families 
who are identified by health records as long-lived. The study tracks family 
members over three generations in a series of three visits each. Each individual's 
genome is fully assembled to generate genotype data,
and the result of each of the three visits are measurements of a large number 
of phenotypes. One such phenotype measurement is whole blood RNA sequencing.

Unfortunately, a mishap has occurred in at least 15 of the RNA sequence
samples. These samples display a gene expression pattern which contrasts 
with the sex of the study subject to which the sample is attributed.

To attempt to resolve this, I hypothesize that there exists a set of SNPs which 
are located in highly and consistently expressed exons. For each individual 
genome, the genotype at each SNP locus will form a vector which, like a finger 
print, will uniquely identify individual genomes. I can then extract the 
consensus base, or bases in the case of heterozygous loci, from RNA sequence 
library alignments and quantify the hamming distance between 
the RNA sequence library genotype vector and all of the individual genome 
genotype vectors. The result will be a ranked list which may be used to 
correct the mislabeled samples in the study metadata.

The purpose of this project is to determine the feasibility of matching RNA 
sequencing libraries to individual genomes.

# Method

The genomes are assembled and genotyped in reference to human genome version hg38. 
A variant call file (VCF) exists for each genome which details all 
variants for each individual. Only high confidence single 
nucleotide polymorphisms over exonic regions are used for this project.

The RNA sequence reads are aligned with STAR, a splice aware aligner,
to the same reference genome. Samtools pileup was used to determine 
the aligned base and depth over the exome for each RNA sequence library. Only 
reads which have a high alignment score and a depth of at least 5 are used in
this project.  

There are 3,688 genes which are expressed at at least 10 log2 counts per 
million in all of the roundly 1000 samples which have been RNA sequenced thus 
far. Using these regions, which are distributed over all autosomal chromosomes, 
I created four genotype vectors (samples of SNP loci, in other words) of length
1,900 each. 

As proof of concept, I selected 7 of the mislabeled samples. Since there are 
two visits to the same individual, I selected both visit 1 and visit 2 for each 
for a total of 14 known mislabeled samples. If the sample is correctly 
identified to a given genome using the method I have outlined above, visit 1 
and visit 2 should agree as to which genome that is. I 
then selected visit 1 and visit 2 libraries for 3 other RNA sequencing samples 
to act as a control. These libraries pass all quality thresholds (all of the
mislabels also pass these thresholds), and are not known mislabels. 
Thus, there is a total of twenty RNA sequence libraries in this project, 
which will be compared against 4,556 individual genomes.

## Nitty Gritty

### Data Preparation

While it may seem trivial to extract this data, for me, at least, it was not. 
The genome data exists in 238 VCF files -- each chromosome is split into 
multiple chunks -- which average 4GB each. Reading these files into memory 
and extracting data is, time-wise, not feasible, nor is it logistically simple 
since each chromosome is split into at least 7 parts. As a result, it was not 
immediately clear how to extract and store the data to make experimentation 
possible.  

I initially tried to use [Plink](https://zzz.bwh.harvard.edu/plink/), 
which seems to be commonly used among the geneticists at Wash U. However, 
for my purposes, this did not provide the quick access to the data that that I 
needed. Simply extracting a SNP matrix (dimensions SNP by samples) ran for more 
than a day, and the output format was not useful.  

The solution I finally arrived at was to extract the data from the VCF files 
into a [SQLite database](https://sqlite.org/index.html). 
Using 
[VariantAnnotation](https://bioconductor.org/packages/release/bioc/html/VariantAnnotation.html), 
an R package, I could fine tune, for speed, how the data streams into my 
parsing program. I could also pre-filter the data so that only high quality 
SNPs are stored.

The table which stores the genome genotype data is called `sample_genotype`, 
and has the following fields: `sample`, `chr` (chromosome), 
`bp` (snp location, 1 indexed from the start to the end of the chromosome), 
`ref_alt` (reference genome genotype, and the alternate allele genotype, 
at that location, eg A/T), and `genotype` (three levels, A/A which 
represents homozygous reference, A/B which represents a heterzygous locus, 
and B/B for homozygous alternate). This table is indexed on `sample`, `chr` and 
`bp` for fast lookup.  

One useful feature of using a database framework is that I can also create 
views, or virtual tables defined by a SQL query of other tables in the database, 
which define the randomly sampled sets of loci to be used to match genomes to 
RNA sequence libraries. For example, the SQL command which creates the first sample looks like this:

```
DROP VIEW IF EXISTS "main"."sample1"; 
CREATE VIEW "sample1" AS 
SELECT * FROM sample_genotype 
WHERE (chr IS 13.0 AND bp BETWEEN 48975912 AND 49209779 OR
       ...
       chr IS 17.0 AND bp BETWEEN 50962174 AND 51120868 OR
       ...
       chr IS 19.0 AND bp BETWEEN 12995475 AND 13098796 OR
       ...
```

For each RNA sequence library, I used 
[Rsamtools](https://bioconductor.org/packages/release/bioc/html/Rsamtools.html), 
an R wrapper for
[Samtools](http://www.htslib.org/doc/samtools-mpileup.html), 
to 'pileup' the alignments from the bam files. Only those reads which met 
the criteria described below were retained in the pileup.

```{r}
# these are the paramters used as input to Rsamtools pileup. scanBamParam 
# and pileupParam are helper functions in Rsamtools which read SAM/BAM flags
# and tags and aid in filtering the alignments

sbp = 
  ScanBamParam(
    flag = 
                  # only count reads for which both R1 and R2 align properly 
                  # (these are reverse strand library, so R1 must align to the 
                  # reverse complement of the reference, while R2 aligns on the 
                  # forward strand with some specified overlap)
      scanBamFlag(isProperPair                = TRUE, 
                  # only accept primary alignments
                  isSecondaryAlignment        = FALSE,
                  isSupplementaryAlignment    = FALSE,
                  # only accept alignments which pass all aligner quality 
                  # metrics
                  isNotPassingQualityControls = FALSE))

p_param = 
              # only count to a depth of 1000 over any given locus (for speed, 
              # and it turned out that few individual bases had anywhere close 
              # to this depth)
  pileupParam(max_depth            = 1000, 
              # only retain in the output bases which have a depth of at least 5
              # reads
              min_nucleotide_depth = 5,
              # count reads from both strands over a given locus
              distinguish_strand   = FALSE,
              # only accept reads with a MAPQ (alignment quality) of at least 10
              min_mapq             = 10,
              # only accept reads with a sequencer quality score of at least 10
              min_base_quality     = 10)
```

### Sample preparation

The genome SNP variants are sampled from the previously identified 4000 highly and 
consistently expressed genes. the result, for one sample, looks like so:

```
> snp_df
# A tibble: 840 × 2
   chr         bp
   <chr>    <int>
 1 12.0  12500010
 2 12.0  12500075
 3 12.0  12500375
 4 12.0  12500383
 5 12.0  12500396
 6 12.0  12500480
 7 12.0  12500556
 8 12.0  12500807
 9 12.0  12500824
10 12.0  12500914
# … with 830 more rows
```
Sampling 1900 high expression genes results, for sample1, in 840 high confidence 
SNPS. A subset of a SNP vector (this one for subject labelled '8') looks like 
this:
```
> sample_genome_vectors$`8`
  [1] "G"         "G"         "G"         "C"         "G"        
  [6] "G"         "A"         "G"         "C"         "G"        
 [11] "G"         "G"         "T"         "G"         "A"        
 [16] "T"         "A"         "G"         "A"         "A"        
 [21] "A"         "T"         "A"         "G"         "G"        
 [26] "A"         "C"         "C"         "T"         "A"        
 [31] "T"         "T"         "C"         "T"         "C"        
 [36] "C"         "T"         "T"         "C"         "G"        
 [41] "C"         "C"         "A"         "G"         "T" 
 ...
```
An RNA sequence library vector looks like this:
```
> pileup_vector
  [1] "G"         "G"         "G"         "C"         "G"        
  [6] "G"         "A"         "G"         "C"         "G"        
 [11] "G"         "Uncertain" "Uncertain" "Uncertain" "Uncertain"
 [16] "Uncertain" "Uncertain" "Uncertain" "Uncertain" "Uncertain"
 [21] "Uncertain" "Uncertain" "Uncertain" "Uncertain" "Uncertain"
 [26] "A"         "Uncertain" "Uncertain" "Uncertain" "Uncertain"
 [31] "Uncertain" "T"         "Uncertain" "Uncertain" "C"        
 [36] "C"         "Uncertain" "Uncertain" "Uncertain" "Uncertain"
 [41] "C"         "Uncertain" "Uncertain" "Uncertain" "Uncertain"
 ...
```
A note here -- I was and am concerned about the number of uncertain bases. So, I
visualized this locus in this sample:

```{r echo=FALSE, eval = TRUE, fig.cap="Figure 1: An RNAseq library SNP locations", out.width = '75%'}
knitr::include_graphics("snp_location.png")
```
The vertical bar in figure 1 is over the first "Uncertain" locus in the 
alignment from which the pileup vector above originates. The detail box 
in the image gives details on one of the reads. It says the following: The read 
is 121 base pairs long. SAM flag 99 means that this read, along with its pair, 
are both mapped (the library is unstranded, so in this case, R1 maps to the 
forward strand while R2 maps to the reverse strand). This location is the 
primary alignment for this read, and it aligns without gaps (Cigar 121M means 
121 matches). There is no soft clipping at the ends of the read and the mate 
is also mapped.  

However, the mapping quality (MAPQ) is 1. That indicates that this 
alignment may occur with ~70% chance in other locations in the genome. And, 
this is the case: the flag NH:3 means that there are 3 other alignments which 
may be as good as this one.  

In thinking about this, I realized that the more stringent MAPQ filter on the 
pileup step, while reasonable for something like RNA sequencing, is not 
necessary for this purpose.  A high quality multi-map does provide evidence of a 
given base's identity. In fact, it increases confidence. Thus, I should only set 
a base_quality score, and accept both primary and secondary mappings which map 
as well as the primary mapping.

So, I generated a second data set with the revised parameters below:

```{r}
sbp = 
  ScanBamParam(
    flag = 
      scanBamFlag(isProperPair                = TRUE, 
                  # accept both primary and secondary alignments which align 
                  # as well as the primary
                  isSecondaryAlignment        = TRUE,
                  # But exclude poor secdonary alignments
                  isSupplementaryAlignment    = FALSE,
                  # only accept alignments which pass all aligner quality 
                  # metrics
                  isNotPassingQualityControls = FALSE))

p_param =
  PileupParam(max_depth            = 250, 
              min_nucleotide_depth = 5,
              distinguish_strand   = FALSE,
              min_base_quality     = 10)

```

The result is significantly fewer 'Uncertain' base locations in the RNA sequence 
library pileups.

```{r, eval = TRUE, fig.height=5, fig.width=10}
filter_comparison_plt
```

### Distance

The hamming distance between an RNA sequence library SNP vector and an 
individual's genome SNP vector is calculated using one of the two 'scoring 
matrices'. Strict only accepts exact matches where all of the homozygous 
bases match each other, and homozygous match in either order, i.e. A/T == T/A.  

The lenient scoring matrix, on the other hand, allows matches between homozygous 
and both homozygous genotypes, and also heterozygous genotypes which include the 
same base pair. For instance, A matches A/A, and any heterozygous genotype 
with A.

```{r, eval = TRUE, echo = FALSE, fig.height=7, fig.width=20, fig.cap="Figure 2: Strict and Lenient Scoring Matricies"}


strict_mat_vis = equivalence_class_dist_matrix(lenient = FALSE) %>% 
  as_tibble(rownames = "genotype") %>%
  tableGrob()

lenient_mat_vis = equivalence_class_dist_matrix(lenient = TRUE) %>%
  as_tibble(rownames = "genotype") %>%
  tableGrob()



patchwork = wrap_elements(panel = strict_mat_vis) + lenient_mat_vis+ 
  plot_annotation(tag_levels = 'A')

patchwork
```

# Results

There are 38 individuals which match another individual over the set of 4,556 
genomes exactly in at least one of the 4 SNP vector samples (Figure 3A). 
Outside of those 38 matches, Figure 3A shows the first percentile, median, 
mean and third percentile of the number of differences between genome SNP 
vectors. This is calculated using the 'strict' difference matrix above (Figure 2).  

Using one sample of _____ SNPs, there are ____ RNA seq libraries which match 
one another exactly (Figure 2B). Figure 2B shows the same statistics for these 
RNA sequence libraries as those described for the genomes in the paragraph above. 

Being unable to differentiate between RNAseq libraries or genomes is obviously 
a problem. However, the remedy is likely to increase the number of regions 
over which I sample SNPs. As this is a proof of concept, the goal is not 
necessarily to unambiguously match libraries to genomes in this experiment, but 
to show that this method is reasonable and yields consistent results. As such, 
what we would expect is for the RNAseq libraries which are the same as one 
another to match the same genomes with the same score, and for those genomes 
which are the same to match RNA seq libraries with the same score.

```{r fig3, eval = TRUE, echo = FALSE, fig.cap="Figure 3: Comparing genome and RNAseq library pileup SNP vectors"}


diff_genomes_grob = difference_between_genomes_df %>% 
  tableGrob()

diff_pileup_grob = equivalence_class_dist_matrix(lenient = TRUE) %>%
  as_tibble(rownames = "genotype") %>%
  tableGrob()



patchwork = wrap_elements(panel = diff_genomes_grob) + diff_pileup_grob+ 
  plot_annotation(tag_levels = 'A')

patchwork

```

Figure 4 shows the results of the experiment in matching pile up SNP vectors 
to genome vectors of the same SNPs using the 'strict' scoring matrix. The top 
scoring (lowest difference ratio) "inferred" subject ID is selected for each 
rnaseq pileup. For instance, if RNAseq library 1 has a difference ratio of 0 
to genome 3, then the 'inferred' label for this the RNAseq library currently 
labelled 1 is 3.

```{r, fig.height=12, fig.width=5}
library(patchwork)
extract_top_n = function(sample_name, pileup_vectors){
  lapply(names(pileup_vectors), function(x) {
                          y = sort(unlist(pileup_vectors[[x]]))[1]
                          tibble(
                            sample = sample_name,
                            subject = x,
                            inferred = names(y),
                            difference_ratio = sort(unlist(y))[1])
  })
}

out = map(names(mapqOver10_no_secondary_results), ~extract_top_n(., mapqOver10_no_secondary_results[[.]]))

z = do.call('rbind', unlist(out, recursive = FALSE)) %>%
  group_by(subject, inferred) %>%
  summarise(mean_diff_ratio = mean(difference_ratio), .groups = 'keep') %>%
  ungroup() %>%
  group_by(subject) %>%
  arrange(mean_diff_ratio, .by_group = TRUE) %>%
  tableGrob()

patchwork = wrap_elements(z)
patchwork
```

What we see is that not a single subject is inferred to be the same as that 
which it is labelled. This is disturbing, since there are ____ samples which 
have been included that are expected to be correctly labelled.  

I am confident in the SNP strings that I am extracting -- I know I am extracting 
genotypes for certain samples with a given label from the sample metadata, and 
matching those genotypes to the corresponding RNAseq base. However, the depth 
over any given base in the RNAseq data is frequently very low. RNA is not as 
stable as DNA, so possibly there are more cases in which a base is called as 
heterozygous due to sequencing error. Another possibility is that, as I mentioned 
above, there simply need to be more SNPs included, though the time to extract 
'genome length' SNPs and do the processing is such that I hesitate to do so 
when this method seems to yield poor results already. Finally, it is also 
possible that there are a lot more mislabels than we expect or hope, though 
I think this is not likely. If it were the case, presumably there are still 
more correctly labelled samples than incorrectly labelled samples. I should, 
therefore, be matching subject to ID more often than chance. I could test this 
by increasing the number of subjects which I examine.

# conclusion

This project tested whether or not it is possible to extract from individual 
genomes a set of SNPs over highly and consistently expressed genes which could 
be used as a 'finger print' with which we may match the individual genome to 
a corresponding RNA sequence library. To test this, I created a database of 
genotypes for 4,556 samples. I next identified abou 3,688 highly and consistently 
expressed genes, and created 4 samples of around ___ base pair locations. With 
these samples, I extracted the consensus base(s) over each of the ___ base 
locations from 40 RNAseq libraries, ____ of which are known to be mislabeled due 
to a mismatch between the gene expression pattern and the sex to which the sample 
is assigned.  

I then calculated the equivalence class based hamming distance (Figure 2) 
between each RNA sequence library SNP vector against all of the 4,556 genomes 
in 4 samples, each with SNP vectors of length _____. I also explored how similar 
the genomes and RNA sequence libraries are within each group (genomes to genomes, 
RNAseq libraries to RNAseq libraries).  

My finding is that my current method does not, with high confidence, match 
genomes to RNAseq libraries. It is possible that increasing the amount 
of data, either by expanding the gene set or by including more RNA seq libraries, 
could yield results. 

```{r}

sample_list = paste0("sample", seq(1,4))

filename_map = read_csv(here("data/filename_id_subject_mislabel_map.csv")) %>%
  mutate(filename = 
           paste0(
             str_remove_all(filename, "_visit|.markdup.sorted.pileup"),
             '.sqlite')) %>%
  mutate(filename = str_replace(filename, "_", "."))


experiment_res = map(sample_list, conduct_experiment, filename_map)
```

# References

https://seqqc.wordpress.com/2015/03/10/calculate-nucelotide-frequency-with-rsamtools-pileup/

# Misc


## mislabels

```
4307 , 1294 , 11221 , 20686 , 20903 , 26715 , 21389 , 25785 , 2333 , 6190 , 6099 , 11185 
```

## supposedly correctly labeled

```
31013 , 1597 , 8161 , 5582 , 15076 , 809 , 15147 , 14639 , 26794 , 13822, 2020,14444 
```
```{r}
sample = 'sample1'

snp_df = sample_snp_ranges(here("data/sample_genotype_db.sqlite"), sample)

sample_genome_vectors = genome_vectors(here("data/sample_genotype_db.sqlite"), 
                                       sample)

pileup_vector = rnaseq_vector(here("data/10367282.1.pileup"), snp_df)

dist_mat = equivalence_class_dist_matrix(lenient = FALSE)

out_list1 = map(sample_genome_vectors, compare_vectors, pileup_vector, dist_mat)
names(out_list1) = names(sample_genome_vectors)
```
