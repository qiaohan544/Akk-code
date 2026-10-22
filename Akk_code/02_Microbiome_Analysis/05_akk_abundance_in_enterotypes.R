#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(dplyr); library(ggplot2)})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grep("^--file=", args)])
pkg_dir <- normalizePath(file.path(dirname(script_arg), ".."))
dat <- read.csv(file.path(pkg_dir, "data", "enterosignature_sample_assignments_with_akk.csv"), check.names=FALSE) %>%
  filter(SPT == "Baseline", Response %in% c("Non_Responder", "Responder"), !is.na(enterosignature)) %>%
  mutate(Response = factor(Response, c("Non_Responder", "Responder")),
         enterosignature = factor(enterosignature, unique(sort(enterosignature))))
ann <- dat %>% group_by(enterosignature) %>%
  summarise(p_value=wilcox.test(log10_akkermansia_muciniphila ~ Response, exact=FALSE)$p.value, .groups="drop") %>%
  mutate(p_fdr=p.adjust(p_value, "BH"), label=sprintf("Wilcoxon p=%.2g\nFDR=%.2g", p_value, p_fdr))
ypos <- dat %>% group_by(enterosignature) %>% summarise(y = max(log10_akkermansia_muciniphila, na.rm=TRUE) + .25, .groups="drop")
ann <- left_join(ann, ypos, by="enterosignature")
p <- ggplot(dat, aes(Response, log10_akkermansia_muciniphila, fill=Response)) +
  geom_boxplot(width=.62, outlier.alpha=.25, linewidth=.3) +
  facet_wrap(~enterosignature, nrow=1) +
  geom_text(data=ann, aes(x=1.5, y=y, label=label), inherit.aes=FALSE, size=2.3) +
  scale_fill_manual(values=c(Non_Responder="#B7C9E2", Responder="#D97A6A")) +
  labs(title="Response-associated A. muciniphila differences within enterosignatures",
       x=NULL, y=expression(log[10](italic("Akkermansia muciniphila")~"relative abundance + 10"^-7))) +
  theme_classic(base_size=7) + theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1), strip.background=element_blank())
ggsave(file.path(pkg_dir, "reproduced.pdf"), p, width=11.4, height=4.6)
