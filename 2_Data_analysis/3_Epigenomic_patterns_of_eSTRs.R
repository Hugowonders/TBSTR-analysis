library(tidyverse)
library(data.table)

dir.create("./tmp/")

# load eSTR locations
estrInfo <- fread("./Demo/Suppl/Suppl_Data3_pSTR_gene_association.txt.gz") %>%
  .[eSTR == "Yes", 1:3] %>%
  set_names(c("chr", "start", "end"))

# load pSTR locations
pstrInfo <- fread("./Demo/Suppl/Suppl_Data2_STR_Infomation.txt.gz") %>%
  .[eSTR == "Yes", 1:3] %>%
  set_names(c("chr", "start", "end"))

# randomly shuffle eSTR positions in a ±200 kb window
estrShuffle200k <- map(1:2000, function(x) {
  set.seed(x)
  
  diff <- sample(seq(-200000, 200000, 1), 1)
  estrInfo %>%
    select(chr, start, end) %>%
    unique() %>%
    mutate(start1 = start + diff,
           end1 = end + diff,
           diff = diff,
           op = x) %>%
    select(chr, start1, end1, diff, op) %>%
    return()
}) %>%
  reduce(rbind)

# write results to a bed file
write_tsv(estrShuffle2k[estrShuffle2k$start1 > 0, ],
          "./tmp/estr_shuffle200k.bed.gz",
          col_names = F)

# randomly shuffle pSTR positions in a ±200 kb window
pstrShuffle200k <- map(1:2000, function(x) {
  set.seed(x)
  diff <- sample(seq(-200000, 200000, 1), 1)
  
  set.seed(x)
  pstrInfo %>%
    .[sample(1:NROW(.), length(unique(hipstrGeneEqtl$site))), ] %>%
    mutate(start1 = start + diff,
           end1 = end + diff,
           diff = diff,
           op = x) %>%
    select(chr, start1, end1, diff, op) %>%
    return()
}) %>%
  reduce(rbind)

# write results to a bed file
write_tsv(pstrShuffle200k[pstrShuffle2k$start1 > 0, ], 
          "./tmp/pstr_shuffle200k.bed.gz",
          col_names = F)

# intersect permutated eSTR positions with GM12878 epigenomic marks
system("bedtools intersect -a ./tmp/estr_shuffle200k.bed.gz \
       -b ./Demo/GM12878_epigenomic_peaks.bed.gz -wo | cut -f4,5,9 | \
       sort | uniq -c | awk 'BEGIN{OFS=\t}{print $1,$2,$3,$4}' > ./tmp/estr_shuffle_region_count.txt")

# intersect permutated pSTR positions with GM12878 epigenomic marks
system("bedtools intersect -a ./tmp/pstr_shuffle200k.bed.gz \
       -b ./Demo/GM12878_epigenomic_peaks.bed.gz -wo | cut -f4,5,9 | \
       sort | uniq -c | awk 'BEGIN{OFS=\t}{print $1,$2,$3,$4}' > ./tmp/pstr_shuffle_region_count.txt")


targets <- c("H3K4me1", "H3K4me3", "H3K36me3", "H3K9ac",
             "H3K27ac", "H3K27me3", "CTCF", "DNase")

# calculate empirical enrichment fold change of marks
estrShuffleEmpericalFc <- list(
  fread("./tmp/estr_shuffle_region_count.txt", header = F) %>%
    set_names(c("count", "diff", "op", "target")) %>%
    mutate(type = "pstr"),
  fread("./tmp/pstr_shuffle_region_count.txt", header = F) %>%
    set_names(c("count", "diff", "op", "target")) %>%
    mutate(type = "estr")) %>%
  reduce(left_join, by = c("op", "diff", "target")) %>%
  plyr::ddply("target", function(x) {
    fc <- map_dbl(1:NROW(x), function(i) {
      x$count.y[i] / mean(x$count.x[i])
    })
    
    x$fc <- fc
    return(x)
  })

# plotting
estrPeakEnrich_plot <- estrShuffleEmpericalFC %>%
  filter(target %in% targets) %>%
  mutate(target = factor(target, levels = targets)) %>%
  ggplot() +
  geom_smooth(aes(x = diff, y = log2(fc), color = target, fill = target)) +
  scale_x_continuous(breaks = c(-2e5, -1e5, 0, 1e5, 2e5),
                     labels = c("-200", "-100", "0", "100", "200"),
                     expand = expansion(mult = c(0.01, 0.013))) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.15))) +
  scale_color_manual(NULL, values = c("#96410E", "#FF7F00", "#e498e4", "#7C4B73",
                                      "#E31A1C", "#1F78B4", "#0099B4", "#46bc9b")) +
  scale_fill_manual(NULL, values = c("#96410E", "#FF7F00", "#e498e4", "#7C4B73",
                                     "#E31A1C", "#1F78B4", "#0099B4", "#46bc9b")) +
  labs(x = "Distance from eSTR (kb)",
       y = "LFC") +
  guides(color = guide_legend(nrow = 2)) +
  theme_classic() +
  theme(legend.position = "inside",
        legend.position.inside = c(0.5, 0.94),
        legend.background = element_rect(fill = "transparent"),
        legend.title = element_text(size = unit(12, "pt")),
        legend.key.size = unit(11, "pt"),
        legend.text = element_text(size = unit(11, "pt")),
        axis.title = element_text(size = unit(12, "pt")),
        axis.text = element_text(size = unit(12, "pt")))

ggview(estrPeakEnrich_plot, width = 120, height = 100, units = "mm", dpi = 500)

file.remove("./tmp/", recursive = T)