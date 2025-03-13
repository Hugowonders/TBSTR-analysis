library(tidyverse)
library(waddR)

# load pSTR dosage of popA, rows are samples and columns are pSTR loci
dosagePopA <- fread("./Demo/dosage_popA.txt.gz") %>%
  as.data.frame()

# load pSTR dosage of popB, rows are samples and columns are pSTR loci
dosagePopB <- fread("./Demo/dosage_popB.txt.gz") %>%
  as.data.frame()

# calculate Rst
rstPopAB <- map2(dosagePopA, dosagePopB,
                 function(a, b) {
                   nPopA <- NROW(a)
                   nPopB <- NROW(b)
                   
                   # exclude pSTRs with missing rate > 0.5
                   if (sum(is.na(a))/nPopA > 0.5 || sum(is.na(b))/nPopB > 0.5) {
                     rst <- NA
                   } else {
                     # in-population variance
                     sw <- mean(c(var(a, na.rm = T), var(b, na.rm = T)))
                     
                     # total variance
                     st <- var(c(a, b), na.rm = T)
                     
                     rst <- (st - sw)/st
                     
                     # assign negative values as zero
                     
                     rst[rst < 0] <- 0
                     
                     return(rst)
                   }
                 })

write_tsv(trdsPopAB, "./rst.txt")

# calculate Tandem Repeat Disparity Score (TRDS)
trdsPopAB <- map2(dosagePopA, dosagePopB,
                  function(a, b) {
                    # exclude pSTRs with missing rate > 0.5
                    nPopA <- NROW(a)
                    nPopB <- NROW(b)
                    
                    if (sum(is.na(a))/nPopA > 0.5 || sum(is.na(b))/nPopB > 0.5) {
                      trds <- NA
                    } else {
                      # in-population variance
                      trds <- waddR::wasserstein.test(na.omit(a), na.omit(b), method = "ASY")
                    } 
                    return(trds)
                  })

write_tsv(trdsPopAB, "./trds.txt")