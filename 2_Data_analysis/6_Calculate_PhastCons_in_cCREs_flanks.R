library(tidyverse)
library(data.table)
library(phastCons100way.UCSC.hg38)

# read in environment association analysis results
envFctAssocRes <- read_xlsx("./Demo/Data_S5_STR_and_environments_associations.xlsx") %>%
  select(c(1, 2, 3, 12)) %>%
  set_names(c("chr", "start", "end", "ccre"))

# phastCons score for cCRE-overlapping STRs associated with the meta environment (hSTRs)
# phastCons in the 2kb-window upstream hSTRs
hstrCcreUpstreamPhast <- map2(envFctAssocRes$chr[envFctAssocRes$ccre == T],
                              envFctAssocRes$start[envFctAssocRes$ccre == T],
                              function(x, y) {
                                # genomic ranges
                                granges <- GRanges(seqnames = x,
                                                   IRanges(start = (y-2000):y,
                                                           width = 1))
                                # get phastCons score
                                phast <- gscores(phastCons100way.UCSC.hg38,
                                                 granges)$default
                                    })

# phastCons in the 2kb-window downstream hSTRs
hstrCcreDownstreamPhast <- map2(envFctAssocRes$chr[envFctAssocRes$ccre == T],
                                envFctAssocRes$end[envFctAssocRes$ccre == T],
                                function(x, y) {
                                  # genomic ranges
                                  granges <- GRanges(seqnames = x,
                                                     IRanges(start = y:(y+2000),
                                                             width = 1))
                                  # get phastCons score
                                  phast <- gscores(phastCons100way.UCSC.hg38,
                                                   granges)$default
                                    })

# phastCons score for random cCRE-disjoint hSTRs

# set random seed
set.seed(123)

# randomly sample cCRE-disjoint hSTRs
hstrCcreDisjointRandom <- envFctAssocRes %>%
  filter(ccre == F) %>%
# keep the number of sites the same    
  .[sample(1:NROW(.), sum(envFctAssocRes$ccre)), ]

randomUpstreamPhast <- map2(hstrCcreDisjointRandom$chr,
                            hstrCcreDisjointRandom$start,
                            function(x, y) {
                              # genomic ranges
                              granges <- GRanges(seqnames = x,
                                                 IRanges(start = (y-2000):y,
                                                         width = 1))
                              # get phastCons score
                              phast <- gscores(phastCons100way.UCSC.hg38,
                                               granges)$default
                                 })

randomDownstreamPhast <- map2(hstrCcreDisjointRandom$chr,
                              hstrCcreDisjointRandom$start,
                              function(x, y) {
                                # genomic ranges
                                granges <- GRanges(seqnames = x,
                                                   IRanges(start = (y+2000):y,
                                                           width = 1))
                                # get phastCons score
                                phast <- gscores(phastCons100way.UCSC.hg38,
                                                 granges)$default
                              })

# combine results
hstrFlankPhast <- rbind(
  data.frame(phast = c(pmap_dbl(hstrCcreUpstreamPhast,
                                ~mean(c(...), na.rm = T)),
                       pmap_dbl(hstrCcreDownstreamPhast,
                                ~mean(c(...), na.rm = T))),
             pos = c(c(-2001):(-1),1:2001),
             group = "ccr-overlapping"),
  data.frame(phast = c(pmap_dbl(randomUpstreamPhast,
                                ~mean(c(...), na.rm = T)),
                       pmap_dbl(randomDownstreamPhast,
                                ~mean(c(...), na.rm = T))),
             pos = c(c(-2001):(-1),1:2001),
             group = "ccr-disjoint")
)
