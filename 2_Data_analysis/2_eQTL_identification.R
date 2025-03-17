# using pSTRs on chr21 as an example
# eQLT association analysis
library(tidyverse)
library(data.table)

# pSTRs flanking ±500 kb of genes
pstrFlank500k <- fread("./Demo/chr21_pstrs_flanking_gene_500kb.txt.gz")

# pstr dosage matrix of MAGE samples
mageDosage <- fread("./Demo/chr21_mage_dosage.txt.gz") %>%
  # zscore normalization
  .[, lapply(.SD, scale)]

# gene expression matrix of MAGE samples (inverse normalized TMM)
mageInvNormTmm <- fread("./Demo/chr21_mage_inv_norm.txt.gz")

# covariates of gene expression
mageCovars <- fread("./Demo/mage_covariates.txt.gz")

# exclude pSTRs with missing rate > 0.5
siteFmiss <- map_dbl(mageDosage, ~sum(is.na(.x))/length(.x))
siteMatch <- which(siteFmiss < 0.5)

# association tests

doseGeneLm <- map2(names(mageDosage[, siteMatch, with = F]),
                   mageDosage[, siteMatch, with = F],
                   function(id, dos) {
                     # test for associations across pstr-gene pairs
                     genes <- unique(pstrFlank500k$geneID[pstrFlank500k$id == id])
                     
                     # LM regression
                     pstrGeneCorPairs <- map(genes, function(g) {
                       # combine ancestry, sex, and PEER factors as covariates 
                       dt <- cbind(mageInvNormTmm[, g, with = F], dos, mageCovars[, -1])
                       
                       # perform LM regression: Y=Xβ+Wα+ϵ
                       lmModel <- lm(expr(!!sym(g) ~ .), data = dt) %>%
                         # extract statistics
                         broom::tidy() %>%
                         .[2, 2:5] %>%
                         mutate(site = id, geneID = g)
                       
                       lmModel
                     }) %>%
                       reduce(rbind)
                     # return statistics for each pSTR-Gene pair
                     return(pstrGeneCorPairs)
                   }) %>%
  # combine results across genes
  reduce(rbind)

# Bonferroni correction of p.values
# significant pSTR-Gene association pairs are defined with p.ajd < 0.05
doseGeneLmAdj <- doseGeneLm %>%
  as.data.table() %>%
  .[, p.adj := p.adjust(p.value, method = "bonferroni"), by = "geneID"]

# identification of eSTRs by Benjamini-Hochberg adjustment of p.adjust wit fdr < 0.05
eqtlSummary <- doseGeneLmAdj %>%
  .[, .SD[which.min(p.adj)], by = "geneID"] %>%
  .[, fdr := p.adjust(p.adj, method = "BH")] %>%
  .[fdr < 0.05]

# permuation control
set.seed(123)
mageDosage <- fread("./Demo/chr21_mage_dosage.txt.gz") %>%
  # zscore normalization
  .[, lapply(.SD, scale)] %>%
  set_names(str_remove(names(.), ".V1$")) %>%
  # random shuffling of sample IDs
  .[sample(1:NROW(.), NROW(.)), ]

# repeat the association test
doseGeneLmShuffle <- map2(names(mageDosage[, siteMatch, with = F]),
                          mageDosage[, siteMatch, with = F],
                   function(id, dos) {
                     # test associatins across genes
                     genes <- unique(pstrFlank500k$geneID[pstrFlank500k$id == id])
                     
                     # LM regression
                     pstrGeneCorPairs <- map(genes, function(g) {
                       # combine ancestry, sex, and PEER factors as covariants 
                       dt <- cbind(mageInvNormTmm[, g, with = F], dos, mageCovars[, -1])
                       
                       # perform LM regression: Y=Xβ+Wα+ϵ
                       lmModel <- lm(expr(!!sym(g) ~ .), data = dt) %>%
                         # extract statistics
                         broom::tidy() %>%
                         .[2, 2:5] %>%
                         mutate(site = id, geneID = g)
                       
                       lmModel
                     }) %>%
                       reduce(rbind)
                     # return statistics for each pSTR-Gene pair
                     return(pstrGeneCorPairs)
                   }) %>%
  # combine results across genes
  reduce(rbind)
