#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(dplyr); library(ggplot2)})
args <- commandArgs(trailingOnly=FALSE); script_arg <- sub("^--file=", "", args[grep("^--file=", args)])
pkg_dir <- normalizePath(file.path(dirname(script_arg), ".."))
dat <- read.csv(file.path(pkg_dir,"data","enterosignature_bray_curtis_pcoa.csv"), check.names=FALSE)
perm <- read.csv(file.path(pkg_dir,"data","enterosignature_bray_curtis_permanova.csv"), check.names=FALSE) %>% filter(term=="enterosignature")
p <- ggplot(dat, aes(Bray_PCo1, Bray_PCo2, color=enterosignature)) + geom_point(size=1.15, alpha=.72) +
  labs(title="Bray-Curtis ordination by enterosignature",
       subtitle=sprintf("PERMANOVA ES R2 = %.3f, p = %.3g", perm$R2[1], perm$`Pr(>F)`[1]),
       x="PCo1 (18.2%)", y="PCo2 (9.0%)", color=NULL) +
  theme_classic(base_size=7) + theme(plot.title=element_text(face="bold"), legend.position="right")
ggsave(file.path(pkg_dir,"reproduced.pdf"), p, width=3.79, height=3.29)
