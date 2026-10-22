#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(dplyr); library(ggplot2)})
args <- commandArgs(trailingOnly=FALSE); script_arg <- sub("^--file=", "", args[grep("^--file=", args)])
pkg_dir <- normalizePath(file.path(dirname(script_arg), ".."))
dat <- read.csv(file.path(pkg_dir,"data","enterosignature_sample_assignments_with_akk.csv"), check.names=FALSE) %>% filter(!is.na(enterosignature))
test <- read.csv(file.path(pkg_dir,"data","enterosignature_akk_abundance_by_es_tests.csv"), check.names=FALSE)
p <- ggplot(dat, aes(enterosignature, log10_akkermansia_muciniphila, fill=enterosignature)) +
  geom_boxplot(width=.65, outlier.alpha=.25, linewidth=.3, show.legend=FALSE) +
  labs(title="Akkermansia muciniphila abundance differs across enterosignatures",
       subtitle=sprintf("Kruskal-Wallis p = %.3g", test$p_value[1]), x=NULL,
       y=expression(log[10](italic("Akkermansia muciniphila")~"relative abundance + 10"^-7))) +
  theme_classic(base_size=7) + theme(plot.title=element_text(face="bold"), axis.text.x=element_text(angle=45,hjust=1))
ggsave(file.path(pkg_dir,"reproduced.pdf"), p, width=5.4, height=4.2)
