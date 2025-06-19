library("biomaRt")
library("dplyr")

rm(list=ls())

#list of genes or whatever
gen_pan <- read.csv("mats/rsrcs/cancerGenesLG.csv")

#ensembl connection (archive is more reliable)
human = useMart(biomart="ensembl", dataset = "hsapiens_gene_ensembl",
                verbose = TRUE, host = "https://dec2021.archive.ensembl.org")

#biomaRt pull all information you need at once
gene_info <- getBM(
  attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position"),
  filters = "hgnc_symbol",
  values = gen_pan$Hugo.Symbol,
  mart = human
)

#axe duplicates
gene_info <- gene_info %>%
  distinct(hgnc_symbol, .keep_all = TRUE)

#join with original df from csv
df_updated <- gen_pan %>%
  left_join(gene_info, by = c("Hugo.Symbol" = "hgnc_symbol")) %>%
  mutate(chrom = chromosome_name, start = start_position, end = end_position) %>%
  select(Hugo.Symbol, chrom, start, end, Oncogene, Tumor.Supp)

#gene info filled - now save csv
write_csv(df_updated, "gen_inf_fld.csv")
