# using pSTRs on chr1 of HGDP samples as an example

library(tidyverse)
library(plyr)
library(data.table)
library(mashr)
library(poolr)

# pSTRs on chr21 (missing rate < 50%)
pstrInfoChr21 <- fread("./Demo/chr21_pstr_site_info.txt.gz") %>%
  # only pSTRs with heterozygosity > 0.1 were considered
  filter(het > 0.1)

# STR dosage matrix for HGDP samples
hgdpChr21Dosage <- fread("./Demo/chr21_pstr_dosage_HGDP.txt.gz")

# load environment factors of sampling locations for HGDP samples
hgdpSampleEnvFct <- fread("./Demo/hgdp_env_factors.txt.gz") %>%
  
# apply Z-score normalization for altitude, annual solar radiation,
# and annual temperature data
  mutate(altZ = scale(altitude),
         sradZ = scale(srad),
         tavgZ = scale(tavg))

# perform association analysis on environments and STR dosage
envFctAssoc <- map(c("altZ", "sradZ", "tavgZ"),
                   function(envFct) {
                     envFctDosageAssoc <- map2(hgdpChr21Dosage, pstrInfoChr21$site,
                                               function(x, y) {
                                                 
                                                 # test for associations for environmental factors
                                                 dt <- data.frame(dosage = x) %>%
                                                   cbind(select(hgdpSampleEnvFct, envFct))
                                                 
                                                 # lm regression
                                                 lmModel <- lm(dosage ~ ., data = dt) %>%
                                                   # extract statistics
                                                   broom::tidy() %>%
                                                   .[2, 2:5] %>%
                                                   mutate(site = y)
                                               }) %>%
                       # return association analysis results
                       reduce(rbind) %>%
                       select(6, 2:5) %>%
                       set_names(c("site", paste(c("estimate",
                                                    "std.error",
                                                    "statistic",
                                                    "p.value"), envFct, sep = "_")))
                   }) %>%
  # combine results of tested environmental factors
  reduce(left_join)

# prepare association effects
envFctAssocEstimate <- map(envFctAssoc,
                           ~select(.x, site, startsWith("estimate"))) %>%
  reduce(left_join)

# prepare association effects standard errors
envFctAssocStderror <- map(envFctAssoc,
                           ~select(.x, site, starts_with("std.error"))) %>%
  reduce(left_join)

envFctAssocBhat <- envFctAssocEstimate[, -1]

envFctAssocShat <- envFctAssocStderror[, -1]

# create mash set data
envFctAssocData <- mash_set_data(envFctAssocBhat, envFctAssocShat)

# estimate correlation between environmental factors
V.simple <- estimate_null_correlation_simple(envFctAssocData)

# update mash set data
envFctAssocDataUpdate <- mash_update_data(envFctAssocData, V.simple, ref = NULL)

# set up the covariance matrices
U.c <- cov_canonical(envFctAssocDataUpdate)

# fit the model
envFctAssocDataUpdateMc <-  mash(envFctAssocDataUpdate, U.c)

# extract posterior summaries
envFctAssocDataUpdatePm <- get_pm(envFctAssocDataUpdateMc)

# compute combined posterior means of different environmental factors
envFctAssocCpm <- envFctAssocDataUpdatePm %>%
  as_tibble() %>%
  mutate(site = rownames(envFctAssocDataUpdatePm),
         cpm = rowSums(envFctAssocDataUpdatePm))

# calculate Fisher's P-value for the meta environmental factor
envFctAssocFisherP <- envFctAssoc %>%
  select(site, starts_with("p.value")) %>%
  ddply("site", function(dt) {
    # adjusting P-value by fisher's method
    fisher <- poolr::fisher(unlist(dt[, 2:4]))
    # Fisher's P-values
    dt$fisherp <- fisher$p
    # statistics
    dt$X <- fisher$statistic
    
    return(dt)
  }) %>%
  # adjust Fisher's P
  mutate(fisherPadj = p.adjust(fisherp, method = "BH"))

# combine results together
envFctAssocFinal <- envFctAssoc %>%
  left_join(envFctAssocCpm) %>%
  left_join(envFctAssocFisherP)