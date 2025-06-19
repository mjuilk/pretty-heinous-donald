# Libraries
library("ggplot2")
library("dplyr")

setwd("~/Documents/data/hapmap")
rm(list=ls())
chd7h1 <- read.csv("chd7h1_summ.out")
chd7h2 <- read.csv("chd7h2_summ.out")
ube3a <- read.csv("ube3a_summ.out")
df <- chd7h1 %>% bind_rows(chd7h2) %>% bind_rows(ube3a)
plot <- ggplot(data=df, aes(x=gene, y=pass_count, fill=code)) +
  geom_bar(stat="identity")
plot + scale_fill_manual(values=c("green", "yellow", "red"))
