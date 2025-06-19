# Libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(hrbrthemes)
library(gtools)

#fresh environment
rm(list=ls())

#data import
rer3_data <- read.csv("~/Documents/data/rer/out/epi/rer3/rer3_1_suptr.bedmethyl", sep="", stringsAsFactors = F)

#chromosome list
chrom_ls <- mixedsort(unique(rer3_data$chrom))

#test gene list
name <- c('AKT1', 'ALK', 'APC') 
chr <- c('chr14','chr2','chr5')
start <- c(104769349,29192774,112707498) 
end <- c(104795759,29921586,112846239)
gene_df<- data.frame(name, chr, start, end) 
gene_df 

#subsets
thresh_subset <- rer3_data[rer3_data$pct_mod >= 1, ]
chrcode_subset <- thresh_subset[c("chrom", "code", "pct_mod")]
chr1to5_subset <- rer3_data[rer3_data$chrom %in% c("chr1", "chr2", "chr3", "chr4", "chr5"),]
chr1to10_subset <- rer3_data[rer3_data$chrom %in% chrom_ls[1:10],]
gene_test_subset <- rer3_data[rer3_data$chrom %in% gene_df$chr,]

#alphanumeric sorting
full_data_sorted <- rer3_data[gtools::mixedorder(rer3_data$chrom),]
full_data_sorted2 <- factor(rer3_data$chrom, levels = gtools::mixedsort(rer3_data$chrom))

gene_data_sorted <- gene_test_subset[gtools::mixedorder(gene_test_subset$chrom),]
gene_list_sorted <- gene_df[gtools::mixedorder(gene_df$chr),]

# Plot x as chrom and y as starting base position
p <- ggplot(rer3_data, aes(x=chrom, y=start0)) + geom_point()
chr_labs <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
              "11", "12", "13", "14", "15", "16", "17", "18", "19",
              "20", "21", "22", "X", "Y")

#sub jitter plot
sub_jitterp <- ggplot(chr11_12_13_subset, aes(x = chrom, y = start0, color = pct_mod)) +
  geom_jitter(width = 0.2, height = 0, alpha = 0.6, size = 2) +  
  scale_color_gradient(low =  "yellow", high = "red") +   
  theme_minimal() +
  labs(
    title = "Jitter Plot for Chromosome 11-13 Methylation (REREAD3)",
    x = "Chromosome",
    y = "Base Position",
    color = "Modification"
  )

#sub hex plot (has potential)
sub_hex1 <- ggplot(rer3_data, aes(x=chrom, y=start0, z = pct_mod)) +
  stat_summary_hex(aes(fill = after_stat(value)), bins = 50, fun = mean) +
  scale_fill_gradient(low = "yellow", high = "red") +
  labs(title = "Methylation Map (REREAD3)", subtitle = "50 bins", x = "Chromosome", y = "Base Position") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#sub tile plot (doesn't work, horseshit)
sub_tilep <- ggplot(chr11_12_13_subset, aes(x = chrom, y = start0, fill = pct_mod)) +
  geom_tile(color = "black") +
  scale_fill_gradient(low = "yellow", high = "red") +
  coord_fixed()

#plot with gene markers (test with just 3 genes first)
sub_hex2 <- ggplot(gene_data_sorted, aes(x=chrom, y=start0, z = pct_mod)) +
  stat_summary_hex(aes(fill = after_stat(value)), bins = 25, fun = mean) +
  geom_segment(aes(x = 0.9, y = gene_list_sorted$start[1], xend = 1.1, yend = gene_list_sorted$end[1]), color = "blue", linetype = 2, linewidth = 1) +
  geom_segment(aes(x = 1.9, y = gene_list_sorted$start[2], xend = 2.1, yend = gene_list_sorted$end[2]), color = "blue", linetype = 2, linewidth = 1) +
  geom_segment(aes(x = 2.9, y = gene_list_sorted$start[3], xend = 3.1, yend = gene_list_sorted$end[3]), color = "blue", linetype = 2, linewidth = 1) +
  scale_fill_gradient(low = "green", high = "red") +
  labs(title = "Methylation Map (REREAD3) : AKT1, ALK, APC", subtitle = "25 bins", x = "Chromosome", y = "Base Position") +
  scale_x_discrete(labels = unique(gene_data_sorted$chrom)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))

geom_segment(aes(x = 0.9, y = gene_df$start[1], xend = 1.1, yend = gene_df$end[1]), color = "green", linetype = 2, linewidth = 1) +
  geom_segment(aes(x = 1.9, y = gene_df$start[2], xend = 2.1, yend = gene_df$end[2]), color = "green", linetype = 2, linewidth = 1) +
  geom_segment(aes(x = 2.9, y = gene_df$start[3], xend = 3.1, yend = gene_df$end[3]), color = "green", linetype = 2, linewidth = 1) +
  scale_x_discrete(expand = c(1,0)) +

# Summarize data: count occurrences of 'h' and 'm' for each chromosome
summary_data <- chrcode_subset %>%
  group_by(chrom, code) %>%
  summarise(count = n(), .groups = 'drop') %>%
  pivot_wider(names_from = code, values_from = count, values_fill = 0)
summary_data <- summary_data[mixedorder(summary_data$chrom), ]

# Plot showing h and m modifications per chrom
ggplot(rer3_data, aes(x=x, y=y)) +
  geom_point() + 
  geom_segment( aes(x=x, xend=x, y=0, yend=y))