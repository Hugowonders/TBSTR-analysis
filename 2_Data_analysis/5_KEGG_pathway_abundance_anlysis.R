library(tidyverse)
library(data.table)
library(readxl)
library(clusterProfiler)

# read in eQTL analysis results
eqtlRes <- fread("./Demo/Data_S3_eQTL.txt.gz",
                 select = c(1, 2, 6, 8, 9)) %>%
  set_names(c("chr", "start", "geneSymbol", "geneType", "beta"))

# read in environment association analysis results
envFctAssocRes <- read_xlsx("./Demo/Data_S5_STR_and_environments_associations.xlsx") %>%
  select(c(1, 2, 6)) %>%
  set_names(c("chr", "start", "cpm"))

# environment and gene expression associated pSTRs (adaptive STRs)
adaptStrRes <- inner_join(envFctAssocRes, eqtlRes,
                          by = c("chr", "start")) %>%
  filter(geneType == "protein_coding")

# convert gene symbol to ENTREZ id
targetGeneEntrez <- adaptStrRes %>%
  .$geneSymbol %>%
  unique() %>%
  clusterProfiler::bitr(fromType = "SYMBOL",
                        toType = "ENTREZID",
                        OrgDb = "org.Hs.eg.db")

# perform KEGG enrichment analysis
targetGeneKegg <- clusterProfiler::enrichKEGG(
  gene = targetGeneEntrez$ENTREZID,
  organism = "hsa",
  pAdjustMethod = "BH") %>%
  clusterProfiler::setReadable(keyType = "ENTREZID",
                               OrgDb = "org.Hs.eg.db")

# significantly enriched pathways (q.value < 0.05)
sigKeggPathway <- targetGeneKegg@result %>%
  filter(qvalue < 0.05)

# background gene list
keggBgGene <- targetGeneKegg@result$geneID %>%
  str_split("/", simplify = T) %>%
  as.character() %>%
  unique()

# putative expression levels of background genes mediated by adaptive STRs
keggBgGeneExpression <- adaptStrRes %>%
  .[, .(effect = sum(beta * cpm)), by = "geneSymbol"] %>%
  .[, trend := ifelse(level > 0, "up", "down")]

# differential abundance scores of significantly enriched pathways
sigKeggPathwayDas <- sigKeggPathway %>%
  plyr::ddply(c("ID"), function(x) {
    # extract genes in the pathway
    genes <- str_split(x$geneID, "/") %>% unlist()
    # count up- and down-regulated genes by adaptive STRs
    upCount <- keggBgGeneExpression %>%
      filter(geneSymbol %in% genes, trend == "up") %>%
      NROW()
    downCount <- keggBgGeneExpression %>%
      filter(geneSymbol %in% genes, trend == "down") %>%
      NROW()
    
    x$upCount <- upCount
    x$downCount <- downCount
    # return updated data
    x
  }) %>%
  # calculate differential abundance score (DAS)
  mutate(das = (upCount - downCount)/sqrt(Count),
         total = upCount + downCount,
         size = Count/mean(Count)) %>%
  
# permutation test
  plyr::ddply("ID", function(pw) {
    ab <- map_dbl(1:2000, function(x) {
      # set random seeds
      set.seed(x)
      
      # randomly sample N genes from the total tested genes
      genes <- sample(targetGeneEntrez$SYMBOL, pw$Count)
      
      # expression levels of sampled genes
      trend <- keggBgGeneExpression$trend[keggBgGeneExpression$geneSymbol %in% genes]
      
      # calculate the null distribution of differential abundance score
      z <- (sum(trend == "up") - sum(trend == "down"))/sqrt(pw$Count)
      return(z)
    })
    
    # expected DAS of the pathway
    exp <- mean(ab)
    
    # calculate empirical P-value
    p <- if(pw$das >= exp) {sum(pw$das < ab)/(2000 - 1)}
    else {sum(pw$das >= ab)/(2000 - 1)}
    
    # update the data
    pw$exp <- exp
    pw$p <- p
    return(pw)
  })