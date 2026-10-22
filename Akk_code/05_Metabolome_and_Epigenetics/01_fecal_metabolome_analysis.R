#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(ggplot2); library(ggrepel); library(patchwork); library(ragg)
})
set.seed(20260818)
args <- commandArgs(trailingOnly=FALSE)
script_dir <- normalizePath(dirname(sub("^--file=", "", args[grep("^--file=", args)][1])))
package_dir <- normalizePath(file.path(script_dir,".."))
out_dir <- file.path(package_dir,"data")
figure_dir <- file.path(package_dir,"figure")
fc_cut <- 2
p_cut <- 0.05

a <- read_csv(file.path(out_dir,"filtered_differential_metabolites_AB_vs_A.csv"),show_col_types=FALSE) %>%
  transmute(feature_id,ion_id,type,Name,log2FC_AB_vs_A=log2FC,p_AB_vs_A=pvalue_directional,
            fdr_AB_vs_A=padj_directional,sig_AB_vs_A=is_sig)
b <- read_csv(file.path(out_dir,"filtered_differential_metabolites_AB_vs_AE.csv"),show_col_types=FALSE) %>%
  transmute(feature_id,ion_id,type,Name,log2FC_AB_vs_AE=log2FC,p_AB_vs_AE=pvalue_directional,
            fdr_AB_vs_AE=padj_directional,sig_AB_vs_AE=is_sig)
d <- inner_join(a,b,by=c("feature_id","ion_id","type","Name")) %>%
  mutate(
    x_class=case_when(log2FC_AB_vs_A > fc_cut ~ "AB>A",log2FC_AB_vs_A < -fc_cut ~ "AB<A",TRUE ~ "No strong change"),
    y_class=case_when(log2FC_AB_vs_AE > fc_cut ~ "AB>AE",log2FC_AB_vs_AE < -fc_cut ~ "AB<AE",TRUE ~ "No strong change"),
    quadrant=case_when(
      x_class=="AB>A" & y_class=="AB>AE" ~ "Q9: AB enriched vs A and AE",
      x_class=="No strong change" & y_class=="AB>AE" ~ "Q8: AB enriched vs AE only",
      x_class=="AB<A" & y_class=="AB>AE" ~ "Q7: A high / AE low",
      x_class=="AB>A" & y_class=="No strong change" ~ "Q6: AB enriched vs A only",
      x_class=="No strong change" & y_class=="No strong change" ~ "Q5: no strong fold change",
      x_class=="AB<A" & y_class=="No strong change" ~ "Q4: A enriched vs AB only",
      x_class=="AB>A" & y_class=="AB<AE" ~ "Q3: A low / AE high",
      x_class=="No strong change" & y_class=="AB<AE" ~ "Q2: AE enriched vs AB only",
      TRUE ~ "Q1: A and AE enriched vs AB"
    ),
    any_nominal_significant=sig_AB_vs_A | sig_AB_vs_AE,
    both_nominal_significant=sig_AB_vs_A & sig_AB_vs_AE,
    label_score=pmin(p_AB_vs_A,p_AB_vs_AE,na.rm=TRUE),
    is_valine=grepl("^(DL-|L-)?valine$",Name,ignore.case=TRUE),
    is_valine_dipeptide=grepl("^(DL-)?(Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|Lys|Met|Phe|Pro|Ser|Thr|Trp|Tyr|Val)-(DL-)?(Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|Lys|Met|Phe|Pro|Ser|Thr|Trp|Tyr|Val)$",Name,ignore.case=TRUE) &
      grepl("(^|-)Val($|-)",Name,ignore.case=TRUE),
    significant_valine_related=(is_valine | is_valine_dipeptide) & any_nominal_significant,
    is_scfa=grepl("^(Acetic acid|Acetate|Propionic acid|Propionate|Butyric acid|Butyrate|Butanoic acid|Isobutyric acid|Isobutyrate|Valeric acid|Valerate|Isovaleric acid|Isovalerate|Caproic acid|Caproate|Hexanoic acid|Hexanoate)$",Name,ignore.case=TRUE),
    significant_scfa=is_scfa & any_nominal_significant
  )

quad_levels <- c("Q1: A and AE enriched vs AB","Q2: AE enriched vs AB only","Q3: A low / AE high",
                 "Q4: A enriched vs AB only","Q5: no strong fold change","Q6: AB enriched vs A only",
                 "Q7: A high / AE low","Q8: AB enriched vs AE only","Q9: AB enriched vs A and AE")
d$quadrant <- factor(d$quadrant,levels=quad_levels)
write_csv(d,file.path(out_dir,"AB_vs_A_AB_vs_AE_joint_nine_quadrant_source_data.csv"))

xlim <- ceiling(max(abs(d$log2FC_AB_vs_A),na.rm=TRUE)*1.06)
ylim <- ceiling(max(abs(d$log2FC_AB_vs_AE),na.rm=TRUE)*1.06)
xb <- c(-xlim,-fc_cut,fc_cut,xlim); yb <- c(-ylim,-fc_cut,fc_cut,ylim)
zones <- expand.grid(xi=1:3,yi=1:3) %>% mutate(
  xmin=xb[xi],xmax=xb[xi+1],ymin=yb[yi],ymax=yb[yi+1],
  q=factor((yi-1)*3+xi,levels=1:9),
  q_label=paste0("Q",(yi-1)*3+xi)
)
zone_cols <- c("#4575B4","#91BFDB","#E0F3F8","#ABD9E9","#F2F2F2","#FEE090","#E0F3A8","#FDAE61","#D73027")
point_cols <- setNames(zone_cols,quad_levels)
valine_label_df <- d %>% filter(significant_valine_related) %>% arrange(label_score)
scfa_label_df <- d %>% filter(is_scfa) %>% arrange(desc(significant_scfa),label_score)
counts <- d %>% count(quadrant,.drop=FALSE)

p9 <- ggplot() +
  geom_rect(data=zones,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=q),alpha=.20,color="white",linewidth=.45) +
  geom_text(data=zones,aes((xmin+xmax)/2,(ymin+ymax)/2,label=q_label),color="grey55",fontface="bold",size=4) +
  geom_point(data=d,aes(log2FC_AB_vs_A,log2FC_AB_vs_AE,color=quadrant,alpha=any_nominal_significant),size=1.45) +
  geom_vline(xintercept=c(-fc_cut,fc_cut),linetype="dashed",linewidth=.35,color="grey30") +
  geom_hline(yintercept=c(-fc_cut,fc_cut),linetype="dashed",linewidth=.35,color="grey30") +
  geom_point(data=valine_label_df,aes(log2FC_AB_vs_A,log2FC_AB_vs_AE),shape=23,size=3.25,fill="#7A0177",color="white",stroke=.55) +
  geom_text_repel(data=valine_label_df,aes(log2FC_AB_vs_A,log2FC_AB_vs_AE,label=Name),size=2.45,fontface="bold",color="#65005F",
                  min.segment.length=0,max.overlaps=Inf,seed=20260818,box.padding=.35,point.padding=.25,segment.color="#9E4A97") +
  geom_point(data=scfa_label_df %>% filter(!significant_scfa),aes(log2FC_AB_vs_A,log2FC_AB_vs_AE),shape=24,size=3.1,fill="white",color="#087E8B",stroke=.85) +
  geom_point(data=scfa_label_df %>% filter(significant_scfa),aes(log2FC_AB_vs_A,log2FC_AB_vs_AE),shape=24,size=3.3,fill="#00A6A6",color="white",stroke=.55) +
  geom_text_repel(data=scfa_label_df,aes(log2FC_AB_vs_A,log2FC_AB_vs_AE,label=Name),size=2.4,fontface="bold",color="#006D77",
                  min.segment.length=0,max.overlaps=Inf,seed=20260819,box.padding=.42,point.padding=.28,segment.color="#4A9DA5") +
  scale_fill_manual(values=setNames(zone_cols,levels(zones$q)),guide="none") +
  scale_color_manual(values=point_cols,drop=FALSE,name="Nine-quadrant class") +
  scale_alpha_manual(values=c(`FALSE`=.30,`TRUE`=.88),guide="none") +
  coord_fixed(xlim=c(-xlim,xlim),ylim=c(-ylim,ylim),clip="off") +
  labs(title="Joint nine-quadrant analysis of AB vs A and AB vs AE",
       subtitle="Purple diamonds: significant valine-related metabolites; teal triangles: SCFAs (filled, significant; open, non-significant); dashed lines: |log2FC|=2",
       x=expression(log[2]~FC~("AB vs A")),y=expression(log[2]~FC~("AB vs AE"))) +
  theme_classic(base_size=8)+theme(plot.title=element_text(face="bold",size=10),plot.subtitle=element_text(size=7,color="grey30"),
    legend.position="right",legend.text=element_text(size=5.8),legend.title=element_text(size=6.5),axis.text=element_text(color="black"),plot.margin=margin(6,8,6,6))

shorten <- function(x,n=38) ifelse(nchar(x)>n,paste0(substr(x,1,n-3),"..."),x)
top_order <- read_csv(file.path(out_dir,"filtered_top20_upregulated_by_log2FC_AB_vs_A.csv"),show_col_types=FALSE) %>%
  arrange(desc(log2FC),pvalue_directional) %>% slice_head(n=20) %>%
  transmute(feature_id,rank=row_number(),lab_txt=sprintf("%02d  %s",rank,shorten(Name)))
top_wide <- top_order %>% inner_join(d,by="feature_id")
top_long <- bind_rows(
  top_wide %>% transmute(feature_id,Name,rank,lab_txt,comparison="AB vs A",log2FC=log2FC_AB_vs_A,
                         p_value=p_AB_vs_A,is_sig=sig_AB_vs_A),
  top_wide %>% transmute(feature_id,Name,rank,lab_txt,comparison="AB vs AE",log2FC=log2FC_AB_vs_AE,
                         p_value=p_AB_vs_AE,is_sig=sig_AB_vs_AE)
) %>% mutate(
  lab=factor(lab_txt,levels=rev(unique(top_order$lab_txt))),
  comparison=factor(comparison,levels=c("AB vs A","AB vs AE")),
  p_label=if_else(p_value < 0.001,sprintf("P=%.1e",p_value),sprintf("P=%.3f",p_value)),
  text_x=log2FC+if_else(log2FC>=0,.08,-.08),text_hjust=if_else(log2FC>=0,0,1)
)
write_csv(top_long,file.path(out_dir,"top20_AB_vs_A_ordered_two_comparisons_source_data.csv"))

dodge <- position_dodge(width=.76)
t_combined <- ggplot(top_long,aes(log2FC,lab,group=comparison))+
  geom_vline(xintercept=0,color="#777777",linewidth=.35)+
  geom_col(aes(fill=comparison,alpha=is_sig),position=dodge,width=.68,color="white",linewidth=.18)+
  geom_text(aes(x=text_x,label=p_label,color=is_sig,hjust=text_hjust),position=dodge,fontface="bold",size=2.25,show.legend=FALSE)+
  scale_fill_manual(values=c("AB vs A"="#B83B3B","AB vs AE"="#E68179"),name=NULL)+
  scale_alpha_manual(values=c(`TRUE`=1,`FALSE`=.48),guide="none")+
  scale_color_manual(values=c(`TRUE`="#7A0177",`FALSE`="#777777"),guide="none")+
  scale_x_continuous(expand=expansion(mult=c(.18,.28)))+
  labs(title="Paired comparison of the top 20 AB-enriched metabolites",
       subtitle="Metabolites and rank follow AB vs A; labels show directional P values",
       x=expression(log[2]~FC),y=NULL)+
  theme_minimal(base_size=7.8)+theme(
    plot.title=element_text(face="bold",size=10),plot.subtitle=element_text(size=6.8,color="grey35"),
    axis.text.y=element_text(size=6.3,color="#202020"),axis.text.x=element_text(color="#4A4A4A"),
    panel.grid.major.y=element_blank(),panel.grid.minor=element_blank(),panel.grid.major.x=element_line(color="#ECECEC",linewidth=.3),
    legend.position="top",legend.key.height=unit(3.2,"mm"),legend.key.width=unit(10,"mm"),
    plot.margin=margin(5,12,5,4))

bottom_inset <- plot_spacer() | t_combined | plot_spacer()
bottom_inset <- bottom_inset + plot_layout(widths=c(.09,.82,.09))
combined <- p9/bottom_inset+plot_layout(heights=c(1.35,1))+plot_annotation(tag_levels=list(c("A","C"))) &
  theme(plot.tag=element_text(face="bold",size=11))

pdf_out <- file.path(figure_dir,"filtered_AB_vs_A_AB_vs_AE_joint_nine_quadrant_top20.pdf")
png_out <- file.path(figure_dir,"filtered_AB_vs_A_AB_vs_AE_joint_nine_quadrant_top20.png")
ggsave(pdf_out,combined,width=13,height=14,device="pdf",useDingbats=FALSE)
ggsave(png_out,combined,width=13,height=14,device=ragg::agg_png,dpi=600,background="white")
message("Wrote: ",pdf_out)
message("Wrote: ",png_out)
