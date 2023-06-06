---
title: "Analysis steps"
author: "JSY"
date: "2023-06-05"
output: html_document
editor_options: 
  markdown: 
---

# 1.Process of raw sequence data

Datasets from different Bioproject were downloaded and denoised independently, then merged for Taxonomic classification

## 1.1 Import raw sequence

### Paired end

```{bash}
	time qiime tools import \
	  --type 'SampleData[PairedEndSequencesWithQuality]' \
	  --input-path manifest \
	  --output-path demux.qza \
	  --input-format PairedEndFastqManifestPhred33V2
```

### Single end

```{bash}
	time qiime tools import \
	  --type 'SampleData[SequencesWithQuality]' \
	  --input-path manifest \
	  --output-path demux.qza \
	  --input-format SingleEndFastqManifestPhred33V2
```

## 1.2 Denoise

### Paired end

```{bash}
	qiime dada2 denoise-paired \
	  --i-demultiplexed-seqs demux.qza \
	  --p-trim-left-f 29 --p-trim-left-r 18 \
	  --p-trunc-len-f 0 --p-trunc-len-r 0 \
	  --o-table dada2-table.qza \
	  --o-representative-sequences dada2-rep-seqs.qza \
	  --o-denoising-stats denoising-stats.qza
```

### Single end

```{bash}
	qiime dada2 denoise-single \
	  --i-demultiplexed-seqs demux.qza \
	  --p-trim-left 29 \
	  --p-trunc-len 0  \
	  --o-table dada2-table.qza \
	  --o-representative-sequences dada2-rep-seqs.qza \
	  --o-denoising-stats denoising-stats.qza
```

the trim- and trunc- length of each dataset were set according to the results from documents "demux.qza" and summarized in the additional file 2

## 1.3 Merge feature-table

```{bash}
  qiime feature-table merge \
  --i-tables table-1.qza \
  --i-tables table-2.qza \
  --o-merged-table table.qza
```

the table name were summarized in the additional file 2

## 1.4 Taxonomic classification

```{bash}
	qiime phylogeny align-to-tree-mafft-fasttree \
	  --i-sequences rep-seqs.qza \
	  --o-alignment aligned-rep-seqs.qza \
	  --o-masked-alignment masked-aligned-rep-seqs.qza \
	  --o-tree unrooted-tree.qza \
	  --o-rooted-tree rooted-tree.qza

	qiime diversity core-metrics-phylogenetic \
	  --i-phylogeny rooted-tree.qza \
	  --i-table table.qza \
	  --p-sampling-depth 214 \
	  --m-metadata-file metadata.txt \
	  --output-dir core-metrics-results

	qiime feature-classifier classify-sklearn \
	  --i-classifier classifier.qza \
	  --i-reads rep-seqs.qza \
	  --o-classification taxonomy.qza

	qiime taxa collapse \
	  --i-table table.qza \
	  --i-taxonomy taxonomy.qza \
	  --p-level 7 \
	  --o-collapsed-table table-l7.qza
```

> classfiers included classifier_V3-V4_785, classifier_V3-V4_805, classifier_V4_515-806, classifier_V4_515-805, classifier_V4_563-926,classifier_V4-V5_926, classifier_V5-V6_1061 trained according to the primer information in each study

## 1.5 Feature table export

```{bash}
	qiime tools export \
  --input-path table-l7.qza \
  --output-path exported-table-l7
  
  biom convert -i table.biom -o table-l7_merge_all.txt --to-tsv

```

# 2.Heterogeneity analysis

table-l7_merge_all.txt was then merged with metadata.txt, added with shannon index, age was stratified as: \<3,1;3-40,2; 41-55,3; 56-70,4; \>70,5; time point was stratified as: -1,1; 0-30,2; 30-100,3; \>100,4 \> Packages used in heterogeneity analysis

```{R}

library(readr)
library(data.table)
library(VIM)
library(car)
library(magrittr)
library(plyr)
library(dplyr)
library (meta)
library(reshape2)
library(ggplot2)
library(vegan)
library(MMUPHin)
library(patchwork)
library(verification)

```

## 2.1 data import and filtering

```{r}
adonis_all_level7 <- read.csv("level-7-all.csv",header=T,stringsAsFactors = FALSE, fileEncoding = 'GBK')

adonis_all <- adonis_all_level7
adonis_all[,2:1513][is.na(adonis_all[,2:1513])] <- 0

sum_abundance <- apply(adonis_all[,2:1513],1,sum)

adonis_all_rel <- array(0,dim=c(2796,1512))
for (j in 1:1512) {
  for (i in 1:2796){
    adonis_all_rel[i,j] <- adonis_all[i,j+1]/sum_abundance[i]*100
  }
}
adonis_all_rel <- as.data.frame(adonis_all_rel)
adonis_all_rel[is.na(adonis_all_rel)] <- 0

sample <- adonis_all$index
species <- colnames(adonis_all)
colnames(adonis_all_rel) <- species[2:1513]
rownames(adonis_all_rel) <- sample

adonis_all_rel_fil_10 <- filtering_taxonomy(adonis_all_rel,0.00000001,10)

adonis_all_summary_taxnomy_filtering_10 <- read.table("summary_taxonomy_filtering.txt",sep="\t",header=T)
adonis_all_filtered_taxonomy_10 <- read.table("filtered_taxonomy.txt",sep="\t",header=T)

adonis_all_filtered_taxonomy_10_1 <- cbind(adonis_all_filtered_taxonomy_10,
                                         index=rownames(adonis_all_filtered_taxonomy_10))

adonis_all_metadata <- as.data.frame(cbind(index=sample,adonis_all[,1514:1584]))
adonis_all_reg <- merge(adonis_all_filtered_taxonomy_10_1,adonis_all_metadata,by="index")
rownames(adonis_all_reg) <- adonis_all_reg$index
adonis_all_reg_name <- adonis_all_reg
colnames(adonis_all_reg_name)[131] <- "scs"
colnames(adonis_all_reg_name)[132] <- "DrH"

```

> Filtering section used the script available at: <https://github.com/WeersmaLabIBD/Microbiome/blob/master/Tools/Filter_taxonomy.R>

## Figure 1 (cohort study β diversity)

```{r}
all_reg<-adonis_all_reg_name[which(adonis_all_reg_name$disease != "N"),]
all_reg$ave<-apply(all_reg[,2:111],1,mean)
all_reg<-all_reg[which(all_reg$ave!="0"),]
all_reg_matrix<-all_reg[,2:111]
all_reg_metadata<-all_reg[,c(1,112:182)]
all_reg_matrix<-all_reg_matrix*0.01

adonis_beta_all_reg_matrix<-adonis2(all_reg_matrix ~ Study,data = all_reg_metadata,permutations = 1000, method="bray")

title_beta_all_cohort_study <- paste0("adonis R2: ",round(adonis_beta_all_reg_matrix$R2,2),"; P-value < ", round(adonis_beta_all_reg_matrix$`Pr(>F)`,4))


beta_all_reg_matrix<-vegdist(all_reg_matrix,method="bray", binary=F) |> 
  cmdscale( k = 5,eig = T) 

beta_all_reg_points <- as.data.frame(beta_all_reg_matrix$points)
beta_sum_eig <- sum(beta_all_reg_matrix$eig)
beta_eig_percent <- round(beta_all_reg_matrix$eig/beta_sum_eig*100,1)

beta_all_reg_plot<-cbind(index=rownames(beta_all_reg_points),beta_all_reg_points)
beta_all_reg_plot<-merge(beta_all_reg_plot,all_reg_metadata,by="index")

beta_all_reg_plot$Cohort<-factor(beta_all_reg_plot$Cohort,levels=c("HSCT","GVHD","Chemo"))
beta_all_reg_plot$Study<-factor(beta_all_reg_plot$Study)
beta_all_cohort_study<-ggplot(beta_all_reg_plot,aes(V1,V2,color=Study))+
  geom_point (alpha=0.5,size=2, aes(shape=Cohort))+
  scale_color_manual(name="Study",
                     labels=levels(beta_all_reg_plot$Study),
                     values = c(hue_pal()(18)))+
  theme(panel.background = element_blank(),
        panel.grid =element_blank(),
        axis.line = element_line(colour = "grey50"),
        axis.title = element_text(size=14,family="serif",face="bold"),
        axis.text = element_text(size=12,family="serif",face="bold"),
        # legend.position = 'none',
        legend.title = element_text(size=16,family="serif",face="bold"),
        legend.text = element_text(size=15,family="serif",face="bold"),
        plot.subtitle = element_text(size=12,family="serif",face="bold"),
        strip.text = element_text(size=12,family="serif",face="bold"),
        title = element_text(size=14,family="serif",face="bold")
  )+
  labs(x=paste("PCoA 1 (", beta_eig_percent[1], "%)", sep=""),
       y=paste("PCoA 2 (", beta_eig_percent[2], "%)", sep=""),
       title=title_beta_all_cohort_study) 
beta_all_cohort_study 

adonis_jaccard_all_reg_matrix<-adonis2(all_reg_matrix ~ Study,data = all_reg_metadata,permutations = 1000,method="jaccard")

title_jaccard_all_cohort_study <- paste0("adonis R2: ",round(adonis_jaccard_all_reg_matrix$R2,2), 
                                      "; P-value < ", round(adonis_jaccard_all_reg_matrix$`Pr(>F)`,4))


jaccard_all_reg_matrix<-vegdist(all_reg_matrix,method="jaccard", binary=F) %>% 
  cmdscale( k = 5,eig = T) 

jaccard_all_reg_points <- as.data.frame(jaccard_all_reg_matrix$points)
sum_eig <- sum(jaccard_all_reg_matrix$eig)
eig_percent <- round(jaccard_all_reg_matrix$eig/sum_eig*100,1)

jaccard_all_reg_plot<-cbind(index=rownames(jaccard_all_reg_points),jaccard_all_reg_points)
jaccard_all_reg_plot<-merge(jaccard_all_reg_plot,all_reg_metadata,by="index")

jaccard_all_reg_plot$Cohort<-factor(jaccard_all_reg_plot$Cohort,levels=c("HSCT","GVHD","Chemo"))
jaccard_all_reg_plot$Study<-factor(jaccard_all_reg_plot$Study)
jaccard_all_cohort_study<-ggplot(jaccard_all_reg_plot,aes(V1,-V2,color=Study))+
  geom_point (alpha=0.5,size=2, aes(shape=Cohort))+
  scale_color_manual(name="Study",
                     labels=levels(jaccard_all_reg_plot$Study),
                     values = c(hue_pal()(18)))+
  theme(panel.background = element_blank(),
        panel.grid =element_blank(),
        axis.line = element_line(colour = "grey50"),
        axis.title = element_text(size=14,family="serif",face="bold"),
        axis.text = element_text(size=12,family="serif",face="bold"),
        # legend.position = 'none',
        legend.title = element_text(size=16,family="serif",face="bold"),
        legend.text = element_text(size=15,family="serif",face="bold"),
        plot.subtitle = element_text(size=12,family="serif",face="bold"),
        strip.text = element_text(size=12,family="serif",face="bold"),
        title = element_text(size=14,family="serif",face="bold")
  )+
  labs(x=paste("PCoA 1 (", eig_percent[1], "%)", sep=""),
       y=paste("PCoA 2 (", eig_percent[2], "%)", sep=""),
       title = title_jaccard_all_cohort_study)
jaccard_all_cohort_study 

all_cohort_study_plot <- beta_all_cohort_study/jaccard_all_cohort_study +plot_layout(guides = 'collect')
ggsave("all_cohort_study.pdf",all_cohort_study_plot,units="in", width=9.0, height=12.0, dpi=600,limitsize = FALSE)

```

## Figure 1 (HSCT GVHD CHEMO cohort β Divresity)

```{r}

all_reg_HSCT<-adonis_all_reg_name[which(adonis_all_reg_name$disease == "C"),]
all_reg_GVHD<-adonis_all_reg_name[which(adonis_all_reg_name$disease == "G"),]
all_reg_chemo<-adonis_all_reg_name[which(adonis_all_reg_name$disease == "D"),]



cohort_plot<-function(a,b,m){
  #a:dataframe
  #b:number fot geom_point_shape 
  #m:vector for scale_color_manual values
  a$ave<-apply(a[,2:111],1,mean)
  a<-a[which(a$ave!="0"),]
  a_matrix<-a[,2:111]
  a_metadata<-a[,c(1,112:182)]

  a_matrix<-a_matrix*0.01
  
  adonis_beta_matrix<-adonis2(a_matrix ~ Study,data = a_metadata,permutations = 1000, method="bray")
  
  title_beta_study <- paste0("adonis R2: ",round(adonis_beta_matrix$R2,2), 
                             "; P-value < ", round(adonis_beta_matrix$`Pr(>F)`,4))
  
beta_a_matrix<-vegdist(a_matrix,method="bray",binary = F) |> 
  cmdscale( k = 5,eig = T) 

beta_a_points <- as.data.frame(beta_a_matrix$points)
beta_sum_eig <- sum(beta_a_matrix$eig)
beta_eig_percent <- round(beta_a_matrix$eig/beta_sum_eig*100,1)

beta_a_plot<-cbind(index=rownames(beta_a_points),beta_a_points)
beta_a_plot<-merge(beta_a_plot,a_metadata,by="index")

beta_a_plot$Study<-factor(beta_a_plot$Study)
beta_a_plot$Cohort<-factor(beta_a_plot$Cohort,levels=c("HSCT","GVHD","Chemo"))
beta_a_study<-ggplot(beta_a_plot,aes(V1,-V2,color=Study))+
  geom_point (alpha=0.5,size=2, aes(group=Study),shape=b)+
  scale_color_manual(name="Study",
                     labels=levels(beta_a_plot$Study),
                     values = m)+
    theme(panel.background = element_blank(),
        panel.grid =element_blank(),
        axis.line = element_line(colour = "grey50"),
        axis.title = element_text(size=14,family="serif",face="bold"),
        axis.text = element_text(size=12,family="serif",face="bold"),
        legend.position = 'none',
        plot.subtitle = element_text(size=12,family="serif",face="bold"),
        strip.text = element_text(size=12,family="serif",face="bold"),
        title = element_text(size=14,family="serif",face="bold")
  )+
  labs(x=paste("PCoA 1 (", beta_eig_percent[1], "%)", sep=""),
       y=paste("PCoA 2 (", beta_eig_percent[2], "%)", sep=""),
       title = title_beta_study) 
beta_a_study

adonis_jaccard_matrix<-adonis2(a_matrix ~ Study,data = a_metadata,permutations = 1000,method="jaccard")

title_jaccard_study <- paste0("adonis R2: ",round(adonis_jaccard_matrix$R2,2),
                              "; P-value < ", round(adonis_jaccard_matrix$`Pr(>F)`,4))

jaccard_a_matrix<-vegdist(a_matrix,method="jaccard",binary = F) %>% 
  cmdscale( k = 5,eig = T) 

jaccard_a_points <- as.data.frame(jaccard_a_matrix$points)
sum_eig <- sum(jaccard_a_matrix$eig)
eig_percent <- round(jaccard_a_matrix$eig/sum_eig*100,1)

jaccard_a_plot<-cbind(index=rownames(jaccard_a_points),jaccard_a_points)
jaccard_a_plot<-merge(jaccard_a_plot,a_metadata,by="index")

jaccard_a_plot$Study<-factor(jaccard_a_plot$Study)
jaccard_a_plot$Cohort<-factor(jaccard_a_plot$Cohort,levels=c("HSCT","GVHD","Chemo"))
jaccard_a_study<-ggplot(jaccard_a_plot,aes(V1,V2,color=Study))+
  geom_point (alpha=0.5,size=2, aes(group=Study),shape=b)+
  scale_color_manual(name="Study",
                     labels=levels(jaccard_a_plot$Study),
                     values = m)+
  theme(panel.background = element_blank(),
        panel.grid =element_blank(),
        axis.line = element_line(colour = "grey50"),
        axis.title = element_text(size=14,family="serif",face="bold"),
        axis.text = element_text(size=12,family="serif",face="bold"),
        legend.position = 'none',
        plot.subtitle = element_text(size=12,family="serif",face="bold"),
        strip.text = element_text(size=12,family="serif",face="bold"),
        title = element_text(size=14,family="serif",face="bold")
  )+
  labs(x=paste("PCoA 1 (", eig_percent[1], "%)", sep=""),
       y=paste("PCoA 2 (", eig_percent[2], "%)", sep=""),
       title = title_jaccard_study) 
jaccard_a_study
 d <- beta_a_study/jaccard_a_study +plot_layout(guides = 'collect')
return(d)
}

HSCT_color<-c(hue_pal()(18)[c(1,2,3,5,6,7,9,10,13,16,18)])
beta_HSCT_plot<-cohort_plot(all_reg_HSCT,18,HSCT_color)
GVHD_color<-c(hue_pal()(18)[c(1,2,3,5,6,7,9,10,13,16,18)])
beta_GVHD_plot<-cohort_plot(all_reg_GVHD,17,GVHD_color)
chemo_color<-c(hue_pal()(18)[c(7,8,11,12,14,15,17)])
beta_chemo_plot<-cohort_plot(all_reg_chemo,15,chemo_color)

```

## Figure 1 (Merged)

```{r}
Figure1<-beta_HSCT_plot|beta_GVHD_plot|beta_chemo_plot|all_cohort_study_plot+
  plot_layout(guides = 'collect')+
  plot_annotation(tag_levels = "A")
ggsave("Figure1.pdf",Figure1,units="in", width=22.5, height=10.0, dpi=600,limitsize = FALSE)


```

## 2.2 Heterogeneity analysis

### 2.2.1 Heterogeneity analysis of common factors in HSCT GVHD CHEMO cohorts

```{r}
all_reg_metadata$Cohort<-factor(all_reg_metadata$Cohort,levels = c("HSCT","GVHD","Chemo"))
all_reg_metadata$Study<-as.factor(all_reg_metadata$Study)
all_reg_metadata$amplicon<-factor(all_reg_metadata$amplicon,levels = c("V3_V4","V4","V4_V5","V5_V6"))
all_reg_metadata$Platform<-factor(all_reg_metadata$Platform,levels = c("Illumina_MiSeq","Illumina_HiSeq_2500",
 "Illumina_HiSeq_3000","Roche_454_Titanium_pyrosequencing"))
all_reg_metadata$Region<-as.factor(all_reg_metadata$Region)
all_reg_metadata$Nationality<-as.factor(all_reg_metadata$Nationality)

all_reg_metadata_permanova1<-all_reg_metadata |> select(index,Cohort,amplicon,Platform,Region,Nationality,Study)
all_reg_matrix001<-all_reg_matrix*0.01

permanova1_heterogeneity_bray<-array(0,dim=c(6,7))
rownames(permanova1_heterogeneity_bray)<-colnames(all_reg_metadata_permanova1)[2:7]
colnames(permanova1_heterogeneity_bray)<-c("Df","SumofSqs","R2","F","Pval","FDR","Sig")
for (i in 2:7){
  tmp<-adonis2(all_reg_matrix001 ~ all_reg_metadata_permanova1[,i], 
               data = all_reg_metadata_permanova1, permutations = 1000, method="bray")
  permanova1_heterogeneity_bray[i-1,1]<-tmp$Df[1]
  permanova1_heterogeneity_bray[i-1,2]<-tmp$SumOfSqs[1]
  permanova1_heterogeneity_bray[i-1,3]<-tmp$R2[1]
  permanova1_heterogeneity_bray[i-1,4]<-tmp$F[1]
  permanova1_heterogeneity_bray[i-1,5]<-tmp$`Pr(>F)`[1]
}

permanova2_heterogeneity_bray<-array(0,dim=c(6,7))
rownames(permanova2_heterogeneity_bray)<-colnames(all_reg_metadata_permanova1)[2:7]
colnames(permanova2_heterogeneity_bray)<-c("Df","SumofSqs","MeanofSqs","F","Pval","FDR","Sig")
for (i in 2:7){
  tmp1<-vegdist(all_reg_matrix001,method="bray", binary=F) |> 
    betadisper(group=all_reg_metadata_permanova1[,i]) |> 
    permutest()
  permanova2_heterogeneity_bray[i-1,1]<-tmp1$tab[1,1]
  permanova2_heterogeneity_bray[i-1,2]<-tmp1$tab[1,2]
  permanova2_heterogeneity_bray[i-1,3]<-tmp1$tab[1,3]
  permanova2_heterogeneity_bray[i-1,4]<-tmp1$tab[1,4]
  permanova2_heterogeneity_bray[i-1,5]<-tmp1$tab[1,6]
}


permanova1_heterogeneity_jaccard<-array(0,dim=c(6,7))
rownames(permanova1_heterogeneity_jaccard)<-colnames(all_reg_metadata_permanova1)[2:7]
colnames(permanova1_heterogeneity_jaccard)<-c("Df","SumofSqs","R2","F","Pval","FDR","Sig")
for (i in 2:7){
  tmp<-adonis2(all_reg_matrix001 ~ all_reg_metadata_permanova1[,i], 
               data = all_reg_metadata_permanova1, permutations = 1000, method="jaccard")
  permanova1_heterogeneity_jaccard[i-1,1]<-tmp$Df[1]
  permanova1_heterogeneity_jaccard[i-1,2]<-tmp$SumOfSqs[1]
  permanova1_heterogeneity_jaccard[i-1,3]<-tmp$R2[1]
  permanova1_heterogeneity_jaccard[i-1,4]<-tmp$F[1]
  permanova1_heterogeneity_jaccard[i-1,5]<-tmp$`Pr(>F)`[1]
}

permanova2_heterogeneity_jaccard<-array(0,dim=c(6,7))
rownames(permanova2_heterogeneity_jaccard)<-colnames(all_reg_metadata_permanova1)[2:7]
colnames(permanova2_heterogeneity_jaccard)<-c("Df","SumofSqs","MeanofSqs","F","Pval","FDR","Sig")
for (i in 2:7){
  tmp1<-vegdist(all_reg_matrix001,method="jaccard", binary=F) |> 
    betadisper(group=all_reg_metadata_permanova1[,i]) |> 
    permutest()
  permanova2_heterogeneity_jaccard[i-1,1]<-tmp1$tab[1,1]
  permanova2_heterogeneity_jaccard[i-1,2]<-tmp1$tab[1,2]
  permanova2_heterogeneity_jaccard[i-1,3]<-tmp1$tab[1,3]
  permanova2_heterogeneity_jaccard[i-1,4]<-tmp1$tab[1,4]
  permanova2_heterogeneity_jaccard[i-1,5]<-tmp1$tab[1,6]
}


```

### Heterogeneity analysis function

```{r}
permanova_heterogeneity<-function(a,f,am){
  #a, dataframe for heterogeneity analysis
  #f, factor for heterogeneity analysis
  #d, analytical methods
  metadata<-a[,c(1,112:182)]
  matrix001<-a[,2:111]*0.01
  permanova1_heterogeneity<-array(0,dim=c(1,7))
  rownames(permanova1_heterogeneity)<-f
  colnames(permanova1_heterogeneity)<-c("Df","SumofSqs","R2","F","Pval","FDR","Sig")
  vari<- metadata[,which(colnames(metadata)== as.character(f))]
  metadata$vari<-vari
  tmp<-adonis2(matrix001 ~ vari,
               data = metadata, permutations = 1000, method=am)
  permanova1_heterogeneity[1,1]<-tmp$Df[1]
  permanova1_heterogeneity[1,2]<-tmp$SumOfSqs[1]
  permanova1_heterogeneity[1,3]<-tmp$R2[1]
  permanova1_heterogeneity[1,4]<-tmp$F[1]
  permanova1_heterogeneity[1,5]<-tmp$`Pr(>F)`[1]
  
  permanova2_heterogeneity<-array(0,dim=c(1,7))
  rownames(permanova2_heterogeneity)<-f
  colnames(permanova2_heterogeneity)<-c("Df","SumofSqs","MeanofSqs","F","Pval","FDR","Sig")
  tmp1<-vegdist(matrix001,method=am, binary=F) |> 
    betadisper(group=vari) |> 
    permutest()
  permanova2_heterogeneity[1,1]<-tmp1$tab[1,1]
  permanova2_heterogeneity[1,2]<-tmp1$tab[1,2]
  permanova2_heterogeneity[1,3]<-tmp1$tab[1,3]
  permanova2_heterogeneity[1,4]<-tmp1$tab[1,4]
  permanova2_heterogeneity[1,5]<-tmp1$tab[1,6]
  
  permanova_heterogeneity<-list(permanova1_heterogeneity,permanova2_heterogeneity)
  return(permanova_heterogeneity)
  
}

```

### 2.2.2 Heterogeneity analysis of common factors in HSCT GVHD CHEMO cohorts (with missing value)：sex age diagnosis

#### sex

```{r}
all_reg_permanova_sex<-all_reg[complete.cases(all_reg$sex),]
permanova_heterogeneity_bray_sex<-permanova_heterogeneity(all_reg_permanova_sex,"sex","bray")
permanova1_heterogeneity_bray_sex<-permanova_heterogeneity_bray_sex[[1]]
permanova2_heterogeneity_bray_sex<-permanova_heterogeneity_bray_sex[[2]]
permanova_heterogeneity_jaccard_sex<-permanova_heterogeneity(all_reg_permanova_sex,"sex","jaccard")
permanova1_heterogeneity_jaccard_sex<-permanova_heterogeneity_jaccard_sex[[1]]
permanova2_heterogeneity_jaccard_sex<-permanova_heterogeneity_jaccard_sex[[2]]

```

#### age

```{r}
all_reg_permanova_age<-all_reg[complete.cases(all_reg$age_level),]
all_reg_permanova_age$Age<-as.factor(all_reg_permanova_age$age_level)
all_reg_permanova_age<-all_reg_permanova_age[,-130]

permanova_heterogeneity_bray_age<-permanova_heterogeneity(all_reg_permanova_age,"Age","bray")
permanova1_heterogeneity_bray_age<-permanova_heterogeneity_bray_age[[1]]
permanova2_heterogeneity_bray_age<-permanova_heterogeneity_bray_age[[2]]
permanova_heterogeneity_jaccard_age<-permanova_heterogeneity(all_reg_permanova_age,"Age","jaccard")
permanova1_heterogeneity_jaccard_age<-permanova_heterogeneity_jaccard_age[[1]]
permanova2_heterogeneity_jaccard_age<-permanova_heterogeneity_jaccard_age[[2]]

```

#### diagnosis

```{r}
all_reg_permanova_diagnosis<-all_reg[complete.cases(all_reg$diagnosis),]
all_reg_permanova_diagnosis$diagnosis<-as.factor(all_reg_permanova_diagnosis$diagnosis)

permanova_heterogeneity_bray_diagnosis<-permanova_heterogeneity(all_reg_permanova_diagnosis,"diagnosis","bray")
permanova1_heterogeneity_bray_diagnosis<-permanova_heterogeneity_bray_diagnosis[[1]]
permanova2_heterogeneity_bray_diagnosis<-permanova_heterogeneity_bray_diagnosis[[2]]
permanova_heterogeneity_jaccard_diagnosis<-permanova_heterogeneity(all_reg_permanova_diagnosis,"diagnosis","jaccard")
permanova1_heterogeneity_jaccard_diagnosis<-permanova_heterogeneity_jaccard_diagnosis[[1]]
permanova2_heterogeneity_jaccard_diagnosis<-permanova_heterogeneity_jaccard_diagnosis[[2]]

```

### 2.2.3 Heterogeneity analysis for GVHD HSCT cohort

```{r}
all_reg_permanova_HSCTGVHD<-all_reg[which(all_reg$disease != "D"),]
```

#### time

```{r}
all_reg_permanova_HSCTGVHD_time<-all_reg_permanova_HSCTGVHD[complete.cases(all_reg_permanova_HSCTGVHD$time_level),]  
all_reg_permanova_HSCTGVHD_time$Time<-as.factor(all_reg_permanova_HSCTGVHD_time$time_level)
all_reg_permanova_HSCTGVHD_time<-all_reg_permanova_HSCTGVHD_time[,-130]

permanova_heterogeneity_bray_HSCTGVHD_time<-permanova_heterogeneity(all_reg_permanova_HSCTGVHD_time,"Time","bray")
permanova1_heterogeneity_bray_HSCTGVHD_time<-permanova_heterogeneity_bray_HSCTGVHD_time[[1]]
permanova2_heterogeneity_bray_HSCTGVHD_time<-permanova_heterogeneity_bray_HSCTGVHD_time[[2]]
permanova_heterogeneity_jaccard_HSCTGVHD_time<-permanova_heterogeneity(all_reg_permanova_HSCTGVHD_time,"Time","jaccard")
permanova1_heterogeneity_jaccard_HSCTGVHD_time<-permanova_heterogeneity_jaccard_HSCTGVHD_time[[1]]
permanova2_heterogeneity_jaccard_HSCTGVHD_time<-permanova_heterogeneity_jaccard_HSCTGVHD_time[[2]]

```

#### scs

```{r}
all_reg_permanova_HSCTGVHD_scs<-all_reg_permanova_HSCTGVHD[complete.cases(all_reg_permanova_HSCTGVHD$scs),]  
all_reg_permanova_HSCTGVHD_scs$scs<-as.factor(all_reg_permanova_HSCTGVHD_scs$scs)

permanova_heterogeneity_bray_HSCTGVHD_scs<-permanova_heterogeneity(all_reg_permanova_HSCTGVHD_scs,"scs","bray")
permanova1_heterogeneity_bray_HSCTGVHD_scs<-permanova_heterogeneity_bray_HSCTGVHD_scs[[1]]
permanova2_heterogeneity_bray_HSCTGVHD_scs<-permanova_heterogeneity_bray_HSCTGVHD_scs[[2]]
permanova_heterogeneity_jaccard_HSCTGVHD_scs<-permanova_heterogeneity(all_reg_permanova_HSCTGVHD_scs,"scs","jaccard")
permanova1_heterogeneity_jaccard_HSCTGVHD_scs<-permanova_heterogeneity_jaccard_HSCTGVHD_scs[[1]]
permanova2_heterogeneity_jaccard_HSCTGVHD_scs<-permanova_heterogeneity_jaccard_HSCTGVHD_scs[[2]]

```

#### DrH

```{r}
all_reg_permanova_HSCTGVHD_DrH<-all_reg_permanova_HSCTGVHD[complete.cases(all_reg_permanova_HSCTGVHD$DrH),]  
all_reg_permanova_HSCTGVHD_DrH$DrH<-as.factor(all_reg_permanova_HSCTGVHD_DrH$DrH)

permanova_heterogeneity_bray_HSCTGVHD_DrH<-permanova_heterogeneity(all_reg_permanova_HSCTGVHD_DrH,"DrH","bray")
permanova1_heterogeneity_bray_HSCTGVHD_DrH<-permanova_heterogeneity_bray_HSCTGVHD_DrH[[1]]
permanova2_heterogeneity_bray_HSCTGVHD_DrH<-permanova_heterogeneity_bray_HSCTGVHD_DrH[[2]]
permanova_heterogeneity_jaccard_HSCTGVHD_DrH<-permanova_heterogeneity(all_reg_permanova_HSCTGVHD_scs,"DrH","jaccard")
permanova1_heterogeneity_jaccard_HSCTGVHD_DrH<-permanova_heterogeneity_jaccard_HSCTGVHD_DrH[[1]]
permanova2_heterogeneity_jaccard_HSCTGVHD_DrH<-permanova_heterogeneity_jaccard_HSCTGVHD_DrH[[2]]

```

#### GVHD class

```{r}
all_reg_permanova_HSCTGVHD_gclass<-all_reg_permanova_HSCTGVHD[complete.cases(all_reg_permanova_HSCTGVHD$GVHD_class.nong_0.1_1.2.4_2.),]  
all_reg_permanova_HSCTGVHD_gclass$GVHDclass<-as.factor(all_reg_permanova_HSCTGVHD_gclass$GVHD_class.nong_0.1_1.2.4_2.)
all_reg_permanova_HSCTGVHD_gclass<-all_reg_permanova_HSCTGVHD_gclass[,-130]

permanova_heterogeneity_bray_HSCTGVHD_gclass<-permanova_heterogeneity(all_reg_permanova_HSCTGVHD_gclass,"GVHDclass","bray")
permanova1_heterogeneity_bray_HSCTGVHD_gclass<-permanova_heterogeneity_bray_HSCTGVHD_gclass[[1]]
permanova2_heterogeneity_bray_HSCTGVHD_gclass<-permanova_heterogeneity_bray_HSCTGVHD_gclass[[2]]
permanova_heterogeneity_jaccard_HSCTGVHD_gclass<-permanova_heterogeneity(all_reg_permanova_HSCTGVHD_gclass,"GVHDclass","jaccard")
permanova1_heterogeneity_jaccard_HSCTGVHD_gclass<-permanova_heterogeneity_jaccard_HSCTGVHD_gclass[[1]]
permanova2_heterogeneity_jaccard_HSCTGVHD_gclass<-permanova_heterogeneity_jaccard_HSCTGVHD_gclass[[2]]

```

### 2.2.4 Heterogeneity analysis for chemo cohort

```{r}
all_reg_permanova_chemo<-all_reg[which(all_reg$disease == "D"),]
```

#### time

```{r}
all_reg_permanova_chemo_time<-all_reg_permanova_chemo[complete.cases(all_reg_permanova_chemo$time_level),]  
all_reg_permanova_chemo_time$Time<-as.factor(all_reg_permanova_chemo_time$time_level)
all_reg_permanova_chemo_time<-all_reg_permanova_chemo_time[,-130]

permanova_heterogeneity_bray_chemo_time<-permanova_heterogeneity(all_reg_permanova_chemo_time,"Time","bray")
permanova1_heterogeneity_bray_chemo_time<-permanova_heterogeneity_bray_chemo_time[[1]]
permanova2_heterogeneity_bray_chemo_time<-permanova_heterogeneity_bray_chemo_time[[2]]
permanova_heterogeneity_jaccard_chemo_time<-permanova_heterogeneity(all_reg_permanova_chemo_time,"Time","jaccard")
permanova1_heterogeneity_jaccard_chemo_time<-permanova_heterogeneity_jaccard_chemo_time[[1]]
permanova2_heterogeneity_jaccard_chemo_time<-permanova_heterogeneity_jaccard_chemo_time[[2]]

```

### 2.2.5 Tables in additional file

#### bray

```{r}
permanova1_heterogeneity_bray_table<-bind_rows(
  as.data.frame(permanova1_heterogeneity_bray),
  as.data.frame(permanova1_heterogeneity_bray_age),
  as.data.frame(permanova1_heterogeneity_bray_sex),
  as.data.frame(permanova1_heterogeneity_bray_diagnosis),
  as.data.frame(permanova1_heterogeneity_bray_HSCTGVHD_DrH),
  as.data.frame(permanova1_heterogeneity_bray_HSCTGVHD_scs),
  as.data.frame(permanova1_heterogeneity_bray_HSCTGVHD_gclass),
  as.data.frame(permanova1_heterogeneity_bray_HSCTGVHD_time),
  as.data.frame(permanova1_heterogeneity_bray_chemo_time)
)
rownames(permanova1_heterogeneity_bray_table)<-c("Cohort","Amplicon","Platform","Region","Nationality","Study","Age","Sex","Diagnosis","DrH","scs","GVHD class","Time(HSCT/GVHD)","Time(Chemo)")
permanova1_heterogeneity_bray_table$FDR=p.adjust(permanova1_heterogeneity_bray_table$Pval,method = "fdr")
permanova1_heterogeneity_bray_table$Sig<-rep("Yes",14)

permanova2_heterogeneity_bray_table<-bind_rows(
  as.data.frame(permanova2_heterogeneity_bray),
  as.data.frame(permanova2_heterogeneity_bray_age),
  as.data.frame(permanova2_heterogeneity_bray_sex),
  as.data.frame(permanova2_heterogeneity_bray_diagnosis),
  as.data.frame(permanova2_heterogeneity_bray_HSCTGVHD_DrH),
  as.data.frame(permanova2_heterogeneity_bray_HSCTGVHD_scs),
  as.data.frame(permanova2_heterogeneity_bray_HSCTGVHD_gclass),
  as.data.frame(permanova2_heterogeneity_bray_HSCTGVHD_time),
  as.data.frame(permanova2_heterogeneity_bray_chemo_time)
)
rownames(permanova2_heterogeneity_bray_table)<-c("Cohort","Amplicon","Platform","Region","Nationality","Study","Age","Sex","Diagnosis","DrH","scs","GVHD class","Time(HSCT/GVHD)","Time(Chemo)")
permanova2_heterogeneity_bray_table$FDR=p.adjust(permanova2_heterogeneity_bray_table$Pval,method = "fdr")
permanova2_heterogeneity_bray_table$Sig<-rep("Yes",14)
permanova2_heterogeneity_bray_table$Sig[8]<-"No"

write.csv(permanova1_heterogeneity_bray_table,file = "heterogeneity analysis for initial table bray.csv")
write.csv(permanova2_heterogeneity_bray_table,file = "heterogeneity analysis for initial table bray (dispersions).csv")

```

#### jaccard

```{r}
permanova1_heterogeneity_jaccard_table<-bind_rows(
  as.data.frame(permanova1_heterogeneity_jaccard),
  as.data.frame(permanova1_heterogeneity_jaccard_age),
  as.data.frame(permanova1_heterogeneity_jaccard_sex),
  as.data.frame(permanova1_heterogeneity_jaccard_diagnosis),
  as.data.frame(permanova1_heterogeneity_jaccard_HSCTGVHD_DrH),
  as.data.frame(permanova1_heterogeneity_jaccard_HSCTGVHD_scs),
  as.data.frame(permanova1_heterogeneity_jaccard_HSCTGVHD_gclass),
  as.data.frame(permanova1_heterogeneity_jaccard_HSCTGVHD_time),
  as.data.frame(permanova1_heterogeneity_jaccard_chemo_time)
)
rownames(permanova1_heterogeneity_jaccard_table)<-c("Cohort","Amplicon","Platform","Region","Nationality","Study","Age","Sex","Diagnosis","DrH","scs","GVHD class","Time(HSCT/GVHD)","Time(Chemo)")
permanova1_heterogeneity_jaccard_table$FDR=p.adjust(permanova1_heterogeneity_jaccard_table$Pval,method = "fdr")
permanova1_heterogeneity_jaccard_table$Sig<-rep("Yes",14)

permanova2_heterogeneity_jaccard_table<-bind_rows(
  as.data.frame(permanova2_heterogeneity_jaccard),
  as.data.frame(permanova2_heterogeneity_jaccard_age),
  as.data.frame(permanova2_heterogeneity_jaccard_sex),
  as.data.frame(permanova2_heterogeneity_jaccard_diagnosis),
  as.data.frame(permanova2_heterogeneity_jaccard_HSCTGVHD_DrH),
  as.data.frame(permanova2_heterogeneity_jaccard_HSCTGVHD_scs),
  as.data.frame(permanova2_heterogeneity_jaccard_HSCTGVHD_gclass),
  as.data.frame(permanova2_heterogeneity_jaccard_HSCTGVHD_time),
  as.data.frame(permanova2_heterogeneity_jaccard_chemo_time)
)
rownames(permanova2_heterogeneity_jaccard_table)<-c("Cohort","Amplicon","Platform","Region","Nationality","Study","Age","Sex","Diagnosis","DrH","scs","GVHD class","Time(HSCT/GVHD)","Time(Chemo)")
permanova2_heterogeneity_jaccard_table$FDR=p.adjust(permanova2_heterogeneity_jaccard_table$Pval,method = "fdr")
permanova2_heterogeneity_jaccard_table$Sig<-rep("Yes",14)
permanova2_heterogeneity_jaccard_table$Sig[8]<-"No"

write.csv(permanova1_heterogeneity_jaccard_table,file = "heterogeneity analysis for initial table jaccard.csv")
write.csv(permanova2_heterogeneity_jaccard_table,file = "heterogeneity analysis for initial table jaccard (dispersions).csv")

```

#### Drug combination

```{r}
library(stringr)

drugcombination<-adonis_all_level7[,c(1,1515,1516,1535:1583)] 

drugcombination1<-array(0,dim = c(nrow(drugcombination),(ncol(drugcombination)-3)))
colnames(drugcombination1)<-colnames(drugcombination)[4:ncol(drugcombination)]

for (j in 1:nrow(drugcombination)){
  temp_a<-NULL
  for (i in 4:(ncol(drugcombination)-1)){
    if (is.na(drugcombination[j,i])){
      drugcombination1[j,i-3]<-""
      
    }else if (drugcombination[j,i]==1){
      drugcombination1[j,i-3]<-colnames(drugcombination)[i]
      
    }else{
      drugcombination1[j,i-3]<-""
    }
    temp_a<-paste(temp_a,drugcombination1[j,i-3],sep = " ")
  }
  drugcombination1[j,ncol(drugcombination1)]<-temp_a
}
drugcombination1<-as.data.frame(drugcombination1)
drugcombination1<-cbind(drugcombination[,1:3],drugcombination1)
drugcombination1<-drugcombination1 |> select(index:Cohort,Drug.Number)
drugcombination1_hsct<-drugcombination1 |> filter(Cohort=="HSCT")
drugcombination2_hsct<-drugcombination1_hsct[!duplicated(drugcombination1_hsct$subject),]

drugcombination1_gvhd<-drugcombination1 |> filter(Cohort=="GVHD")
drugcombination2_gvhd<-drugcombination1_gvhd[!duplicated(drugcombination1_gvhd$subject),]

drugcombination1_chemo<-drugcombination1 |> filter(Cohort=="Chemo")
drugcombination2_chemo<-drugcombination1_chemo[!duplicated(drugcombination1_chemo$subject),]

drugcombination2<-rbind(drugcombination2_hsct,drugcombination2_gvhd,drugcombination2_chemo)
write.csv(drugcombination2,"drugcombination2.csv")#修改联合用药的格式

drugcombination3<-read.csv("drugcombination3.csv")
drugcombination_table<-drugcombination3 |> group_by(Combinations) |> 
  summarise(
    `Combination frequency`=n()
  )
drugcombination_table<-drugcombination_table |> mutate(`Drug number`)
drugcombination_table$`Drug number`<-str_count(drugcombination_table$Combinations,pattern = "\\+")+1

drugcombination_table_hsct<-drugcombination3 |> filter(Cohort=="HSCT") |> group_by(Combinations) |> 
  summarise(
    `Combination frequency HSCT`=n()
  )

drugcombination_table_gvhd<-drugcombination3 |> filter(Cohort=="GVHD") |> group_by(Combinations) |> 
  summarise(
    `Combination frequency GVHD`=n()
  )

drugcombination_table_chemo<-drugcombination3 |> filter(Cohort=="Chemo") |> group_by(Combinations) |> 
  summarise(
    `Combination frequency Chemo`=n()
  )

drugcombination_table_merge<-full_join(drugcombination_table,drugcombination_table_hsct,by="Combinations")
drugcombination_table_merge<-full_join(drugcombination_table_merge,drugcombination_table_gvhd,by="Combinations")
drugcombination_table_merge<-full_join(drugcombination_table_merge,drugcombination_table_chemo,by="Combinations")

drugcombination_table_merge<-drugcombination_table_merge[,c(1,3,2,4,5,6)]

write.csv(drugcombination_table_merge,"drug combinations.csv")

```

# 3.Multidrug\~Diversity

> Packages used in analyzing association between multidrug use and diversity

```{R}

library(readr)
library(data.table)
library(VIM)
library(car)
library(magrittr)
library(plyr)
library(dplyr)
library (meta)
library(reshape2)
library(ggplot2)
library(vegan)
library(stringr)
library(MMUPHin)
library(scales)
library(patchwork)

```

## 3.1 shannon index \~ number of used drug

### HSCT

```{r}
alpha_control_0<-all_reg_HSCT
alpha_control_0[is.na(alpha_control_0)]<-0

sum_abundance<-apply(alpha_control_0[,2:1417],1,sum)

alpha_control_rel<-array(0,dim=c(838,1416))
for (j in 1:1416) {
  for (i in 1:838){
    alpha_control_rel[i,j]<-alpha_control_0[i,j+1]/sum_abundance[i]*100
  }
}
alpha_control_rel <- as.data.frame(alpha_control_rel)

sample<-alpha_control_0$index
species<-colnames(alpha_control_0)
colnames(alpha_control_rel)<-species[2:1417]
rownames(alpha_control_rel)<-sample

#control_rel_fil_10 <- filtering_taxonomy(control_rel,0.00000001,10)
alpha_control_rel_fil_10 <- filtering_taxonomy(alpha_control_rel,0.00000001,10)

alpha_summary_taxnomy_filtering_10<- read.table("summary_taxonomy_filtering.txt",sep="\t",header=T)
alpha_filtered_taxonomy_10<-read.table("filtered_taxonomy.txt",sep="\t",header=T)

alpha_filtered_taxonomy_10_1<-cbind(index=rownames(alpha_filtered_taxonomy_10),alpha_filtered_taxonomy_10)

alpha_metadata<-as.data.frame(cbind(index=sample,alpha_control[,1418:1482]))
alpha_control_reg<-merge(alpha_filtered_taxonomy_10_1,alpha_metadata,by="index")
rownames(alpha_control_reg)<-alpha_control_reg$index
alpha_control_reg_name<-alpha_control_reg
colnames(alpha_control_reg_name)[115]<-"amplicon"
colnames(alpha_control_reg_name)[116]<-"sex"
colnames(alpha_control_reg_name)[119]<-"scs"
colnames(alpha_control_reg_name)[120]<-"DrH"
alpha_control_singledrug<-alpha_control_reg_name
alpha_control_reg_name$Drug.Numberf<-alpha_control_reg_name$Drug.Number
alpha_control_reg_name[alpha_control_reg_name$Drug.Number>5,]$Drug.Numberf=">5"
alpha_control_reg_name$Drug.Numberf <- factor(alpha_control_reg_name$Drug.Numberf,
                                              levels=c("0","1","2","3","4","5",">5"))
alpha_control_reg_name <- alpha_control_reg_name %>% filter(Drug.Number!="0")
alpha_control_reg_name <- alpha_control_reg_name[,c(1,169:171)] 
alpha_control_reg_name <- alpha_control_reg_name[complete.cases(alpha_control_reg_name),]
alpha_control_reg_name$Cohort = "HSCT"

mean(alpha_control_reg_name[complete.cases(alpha_control_reg_name$shannon_index),]$shannon_index)
min(alpha_control_reg_name[complete.cases(alpha_control_reg_name$shannon_index),]$shannon_index)
max(alpha_control_reg_name[complete.cases(alpha_control_reg_name$shannon_index),]$shannon_index)
fit = lm (shannon_index ~ Drug.Number, data= alpha_control_reg_name)
summary(fit)
p=6.71e-13
estimated=-0.28307 
r2=0.1954
mean=4.38(0.09-8.12)


```

### GVHD

```{r}
alpha_gvhd_0<-all_reg_GVHD
alpha_gvhd_0[is.na(alpha_gvhd_0)]<-0

sum_abundance<-apply(alpha_gvhd_0[,2:1417],1,sum)

alpha_gvhd_rel<-array(0,dim=c(715,1416))
for (j in 1:1416) {
  for (i in 1:715){
    alpha_gvhd_rel[i,j]<-alpha_gvhd_0[i,j+1]/sum_abundance[i]*100
  }
}
alpha_gvhd_rel <- as.data.frame(alpha_gvhd_rel)
alpha_gvhd_rel[is.na(alpha_gvhd_rel)]<-0

sample<-alpha_gvhd_0$index
species<-colnames(alpha_gvhd_0)
colnames(alpha_gvhd_rel)<-species[2:1417]
rownames(alpha_gvhd_rel)<-sample

#gvhd_rel_fil_10 <- filtering_taxonomy(gvhd_rel,0.00000001,10)
alpha_gvhd_rel_fil_10 <- filtering_taxonomy(alpha_gvhd_rel,0.00000001,10)

alpha_summary_taxnomy_filtering_10<- read.table("summary_taxonomy_filtering.txt",sep="\t",header=T)
alpha_filtered_taxonomy_10<-read.table("filtered_taxonomy.txt",sep="\t",header=T)

alpha_filtered_taxonomy_10_1<-cbind(alpha_filtered_taxonomy_10,index=rownames(alpha_filtered_taxonomy_10))

alpha_metadata<-as.data.frame(cbind(index=sample,alpha_gvhd[,1418:1482]))
alpha_gvhd_reg<-merge(alpha_filtered_taxonomy_10_1,alpha_metadata,by="index")
rownames(alpha_gvhd_reg)<-alpha_gvhd_reg$index
alpha_gvhd_reg_name<-alpha_gvhd_reg
colnames(alpha_gvhd_reg_name)[122]<-"amplicon"
colnames(alpha_gvhd_reg_name)[123]<-"sex"
colnames(alpha_gvhd_reg_name)[126]<-"scs"
colnames(alpha_gvhd_reg_name)[127]<-"DrH"

alpha_gvhd_singledrug<-alpha_gvhd_reg_name
alpha_gvhd_reg_name$Drug.Numberf<-alpha_gvhd_reg_name$Drug.Number
alpha_gvhd_reg_name[alpha_gvhd_reg_name$Drug.Number>5,]$Drug.Numberf=">5"
alpha_gvhd_reg_name$Drug.Numberf <- factor(alpha_gvhd_reg_name$Drug.Numberf,
                                           levels=c("0","1","2","3","4","5",">5"))
alpha_gvhd_reg_name <- alpha_gvhd_reg_name %>% filter(Drug.Number!="0")
alpha_gvhd_reg_name <- alpha_gvhd_reg_name[,c(1,176:178)] 
alpha_gvhd_reg_name <- alpha_gvhd_reg_name[complete.cases(alpha_gvhd_reg_name),]
alpha_gvhd_reg_name$Cohort = "GVHD"

mean(alpha_gvhd_reg_name[complete.cases(alpha_gvhd_reg_name$shannon_index),]$shannon_index)
min(alpha_gvhd_reg_name[complete.cases(alpha_gvhd_reg_name$shannon_index),]$shannon_index)
max(alpha_gvhd_reg_name[complete.cases(alpha_gvhd_reg_name$shannon_index),]$shannon_index)
fit = lm (shannon_index ~ Drug.Number, data= alpha_gvhd_reg_name)
summary(fit)
p=0.203
estimated=-0.05507 
r2=0.008354
mean=3.46(0.13-7.09)
```

### CHEMO

```{r}
alpha_chemo_0<-all_reg_chemo
alpha_chemo_0[is.na(alpha_chemo_0)]<-0

sum_abundance<-apply(alpha_chemo_0[,2:850],1,sum)

alpha_chemo_rel<-array(0,dim=c(1097,849))
for (j in 1:849) {
  for (i in 1:1097){
    alpha_chemo_rel[i,j]<-alpha_chemo_0[i,j+1]/sum_abundance[i]*100
  }
}
alpha_chemo_rel <- as.data.frame(alpha_chemo_rel)
alpha_chemo_rel[is.na(alpha_chemo_rel)]<-0

sample<-alpha_chemo_0$index
species<-colnames(alpha_chemo_0)
colnames(alpha_chemo_rel)<-species[2:850]
rownames(alpha_chemo_rel)<-sample

alpha_chemo_rel_fil_10 <- filtering_taxonomy(alpha_chemo_rel,0.00000001,10)

alpha_summary_taxnomy_filtering_10<- read.table("summary_taxonomy_filtering.txt",sep="\t",header=T)
alpha_filtered_taxonomy_10<-read.table("filtered_taxonomy.txt",sep="\t",header=T)

alpha_filtered_taxonomy_10_1<-cbind(alpha_filtered_taxonomy_10,index=rownames(alpha_filtered_taxonomy_10))

alpha_metadata<-as.data.frame(cbind(index=sample,alpha_chemo[,851:915]))
alpha_chemo_reg<-merge(alpha_filtered_taxonomy_10_1,alpha_metadata,by="index")
rownames(alpha_chemo_reg)<-alpha_chemo_reg$index
alpha_chemo_reg_name<-alpha_chemo_reg
colnames(alpha_chemo_reg_name)[109]<-"amplicon"
colnames(alpha_chemo_reg_name)[110]<-"sex"
colnames(alpha_chemo_reg_name)[113]<-"scs"
colnames(alpha_chemo_reg_name)[114]<-"DrH"

alpha_chemo_singledrug<-alpha_chemo_reg_name
alpha_chemo_reg_name$Drug.Numberf<-alpha_chemo_reg_name$Drug.Number
alpha_chemo_reg_name[alpha_chemo_reg_name$Drug.Number>5,]$Drug.Numberf=">5"
alpha_chemo_reg_name$Drug.Numberf <- factor(alpha_chemo_reg_name$Drug.Numberf,
                                               levels=c("0","1","2","3","4","5",">5"))
alpha_chemo_reg_name <- alpha_chemo_reg_name %>% filter(Drug.Number!="0")
alpha_chemo_reg_name <- alpha_chemo_reg_name[,c(1,163:165)] 
alpha_chemo_reg_name <- alpha_chemo_reg_name[complete.cases(alpha_chemo_reg_name),]
alpha_chemo_reg_name$Cohort <- "Chemo"

mean(alpha_chemo_reg_name[complete.cases(alpha_chemo_reg_name$shannon_index),]$shannon_index)
min(alpha_chemo_reg_name[complete.cases(alpha_chemo_reg_name$shannon_index),]$shannon_index)
max(alpha_chemo_reg_name[complete.cases(alpha_chemo_reg_name$shannon_index),]$shannon_index)
fit = lm (shannon_index ~ Drug.Number, data= alpha_chemo_reg_name)
summary(fit)
p=0.492
estimated=0.02634 
r2=0.001906
mean=2.99(0.19-5.98)
```

### Figure 2 (shannon index\~number of used drug)

```{r}
alpha_multidrug_tmp<-rbind(alpha_control_reg_name,alpha_gvhd_reg_name)
alpha_multidrug<-rbind(alpha_multidrug_tmp,alpha_chemo_reg_name)
alpha_multidrug$Cohort<-factor(alpha_multidrug$Cohort,levels=c("HSCT","GVHD","Chemo"))
alpha_multidrug_plot <- ggplot(data=alpha_multidrug,
                               aes(x=Drug.Numberf,y=shannon_index,fill=Cohort))+
  geom_boxplot(outlier.shape = NA, alpha=0.7, colour="grey36")+
  scale_fill_manual(values=c(HSCT = "#BC3C294C", GVHD = "#0072B54C", Chemo = "#0081344C"))+
  geom_point(alpha=0.8,position=position_jitter(0.2))+
  geom_smooth(method="lm",color="black",aes(group=1),linewidth=0.4,fill="#7876B14C")+
  theme(panel.background = element_blank(),
        panel.grid =element_blank(),
        axis.line = element_line(colour = "grey50"),
        axis.title = element_text(size=14,family="serif",face="bold"),
        axis.text = element_text(size=12,family="serif",face="bold"),
        #legend.position = 'none',
        legend.title = element_text(size=16,family="serif",face="bold"),
        legend.text = element_text(size=15,family="serif",face="bold"),
        plot.subtitle = element_text(size=12,family="serif",face="bold"),
        strip.background = element_blank(),
        strip.text = element_text(size=12,family="serif",face="bold"),
        title = element_text(size=14,family="serif",face="bold")
  )+
  labs(x="Number of used Medicine",y="Shannon index")+
  facet_grid(~Cohortf)
alpha_multidrug_plot
ggsave("alpha_multidrug_plot.pdf",alpha_multidrug_plot,units="in", width=12.0, height=6.0, dpi=600,limitsize = FALSE)

```

## 3.2 shannon index \~ single drug use analysis

```{r}
alpha_control_singledrug
shannon_singledrug<-function(a,b){
  #a dataframe
  #b col(matrix)number+1
  shannon_singledrug_table<-array(0,dim = c(48,6))
  colnames(shannon_singledrug_table)<-c("Meds", "MeanUsers","SDUsers", "MeanNonUsers", "SDNonUsers", "pvalue")
  rownames(shannon_singledrug_table)<-colnames(a)[b:(b+47)]
  for (i in b:(b+47)){
  a1<-a[complete.cases(a[,c(1,i,ncol(a))]),]
  if (nrow(a1) != 0){
    if(sum(a1$shannon_index[a1[,i]==0])!=0 & 
       sum(a1$shannon_index[a1[,i]==1])!=0){
      #mean_users
      shannon_singledrug_table[(i+1-b),2]=mean(a1$shannon_index[a1[,i]==1])
      #sd_users
      shannon_singledrug_table[(i+1-b),3]=sd(a1$shannon_index[a1[,i]==1])
      #mean_n_users
      shannon_singledrug_table[(i+1-b),4]=mean(a1$shannon_index[a1[,i]==0])
      #sd_n_users
      shannon_singledrug_table[(i+1-b),5]=sd(a1$shannon_index[a1[,i]==0])
      test=wilcox.test(a1$shannon_index[a1[,i]==0],
                       a1$shannon_index[a1[,i]==1] )
      shannon_singledrug_table[(i+1-b),6]=test$p.value

    }else if(sum(a1$shannon_index[a1[,i]==1])!=0){
      #mean_users
      shannon_singledrug_table[(i+1-b),2]=mean(a1$shannon_index[a1[,i]==1])
      #sd_users
      shannon_singledrug_table[(i+1-b),3]=sd(a1$shannon_index[a1[,i]==1])
      #mean_n_users
      shannon_singledrug_table[(i+1-b),4]=NA
      #sd_n_users
      shannon_singledrug_table[(i+1-b),5]=NA
      shannon_singledrug_table[(i+1-b),6]=NA
    }else if(sum(a1$shannon_index[a1[,i]==0])!=0){
      #mean_users
      shannon_singledrug_table[(i+1-b),2]=NA
      #sd_users
      shannon_singledrug_table[(i+1-b),3]=NA
      #mean_n_users
      shannon_singledrug_table[(i+1-b),4]=mean(a1$shannon_index[a1[,i]==0])
      #sd_n_users
      shannon_singledrug_table[(i+1-b),5]=sd(a1$shannon_index[a1[,i]==0])
      shannon_singledrug_table[(i+1-b),6]=NA
    }else{
      #mean_users
      shannon_singledrug_table[(i+1-b),2]=NA
      #sd_users
      shannon_singledrug_table[(i+1-b),3]=NA
      #mean_n_users
      shannon_singledrug_table[(i+1-b),4]=NA
      #sd_n_users
      shannon_singledrug_table[(i+1-b),5]=NA
      shannon_singledrug_table[(i+1-b),6]=NA
    }
    
  }else{
    #mean_users
    shannon_singledrug_table[(i+1-b),2]=NA
    #sd_users
    shannon_singledrug_table[(i+1-b),3]=NA
    #mean_n_users
    shannon_singledrug_table[(i+1-b),4]=NA
    #sd_n_users
    shannon_singledrug_table[(i+1-b),5]=NA
    shannon_singledrug_table[(i+1-b),6]=NA
    
  }
  
  }
  shannon_singledrug_table <- as.data.frame(shannon_singledrug_table)
  shannon_singledrug_table$Meds<-rownames(shannon_singledrug_table)
  shannon_singledrug_table<-shannon_singledrug_table[!is.na(shannon_singledrug_table$MeanNonUsers),]
  
  return(shannon_singledrug_table)

}

alpha_control_singledrug
b_hsct<-ncol(alpha_control_singledrug)-49
print(b_hsct)
shannon_singledrug_hsct_table<-shannon_singledrug(alpha_control_singledrug,b_hsct)
write.csv(shannon_singledrug_hsct_table,file = "shannon_singledrug_hsct_table.csv")

alpha_gvhd_singledrug
b_gvhd<-ncol(alpha_gvhd_singledrug)-49
print(b_gvhd)
shannon_singledrug_gvhd_table<-shannon_singledrug(alpha_gvhd_singledrug,b_gvhd)
write.csv(shannon_singledrug_gvhd_table,file = "shannon_singledrug_gvhd_table.csv")


alpha_chemo_singledrug
b_chemo<-ncol(alpha_chemo_singledrug)-49
print(b_chemo)
shannon_singledrug_chemo_table<-shannon_singledrug(alpha_chemo_singledrug,b_chemo)
write.csv(shannon_singledrug_chemo_table,file = "shannon_singledrug_chemo_table.csv")

```

## 3.3 beta diversity \~ number of used drugs (PERMANOVA)

```{r}
betadrugnumber_plot<-function(a,b,d,m,n){
  #a:dataframe
  #b:distance method 
  #d:legend title name
  #m:vector for scale_color_manual values
  #n:shape of point
  
  a$ave<-apply(a[,2:111],1,mean)
  a<-a[which(a$ave!="0"),]
  a_matrix001<-a[,2:111]*0.01
  a_metadata<-a[,c(1,112:ncol(a))]
  
  adonis<-adonis2(a_matrix001 ~ Drug.Numberf,data = a_metadata,permutations = 1000, method=b)
  
  title <- paste0("adonis R2: ",round(adonis$R2,2), "; P-value < ", round(adonis$`Pr(>F)`,4))
  
  beta_a_matrix<-vegdist(a_matrix001,method=b,binary = F) |> 
    cmdscale( k = 5,eig = T) 
  
  beta_a_points <- as.data.frame(beta_a_matrix$points)
  beta_sum_eig <- sum(beta_a_matrix$eig)
  beta_eig_percent <- round(beta_a_matrix$eig/beta_sum_eig*100,1)
  
  beta_a_plot<-cbind(index=rownames(beta_a_points),beta_a_points)
  beta_a_plot<-merge(beta_a_plot,a_metadata,by="index")
  
  beta_a_plot$Study<-factor(beta_a_plot$Study)
  # beta_a_plot$Cohort<-factor(beta_a_plot$Cohort)
  beta_a_study<-ggplot(beta_a_plot,aes(-V1,V2,color=Drug.Numberf))+
    geom_point (alpha=1.5,size=2, aes(group=Drug.Numberf),shape=n)+
    scale_color_manual(name=d,
                       labels=levels(beta_a_plot$Drug.Numberf),
                       values = m)+
    theme(panel.background = element_blank(),
          panel.grid =element_blank(),
          axis.line = element_line(colour = "grey50"),
          axis.title = element_text(size=14,family="serif",face="bold"),
          axis.text = element_text(size=12,family="serif",face="bold"),
          # legend.position = 'none',
          legend.title = element_text(size=16,family="serif",face="bold"),
          legend.text = element_text(size=15,family="serif",face="bold"),
          plot.subtitle = element_text(size=12,family="serif",face="bold"),
          strip.text = element_text(size=12,family="serif",face="bold"),
          title = element_text(size=14,family="serif",face="bold")
    )+
    labs(x=paste("PCoA 1 (", beta_eig_percent[1], "%)", sep=""),
         y=paste("PCoA 2 (", beta_eig_percent[2], "%)", sep=""),
         title = title) 
  # facet_grid(~ Cohort)
  beta_a_study
  # ggsave(paste("beta_",a,"_study.pdf",sep = ""),beta_a_study,units="in", width=6.0, height=6.0, dpi=600,limitsize = FALSE)
  
  return(beta_a_study)
}

beta_drugnumber<- adonis_all_reg_name|> select(index,k__Bacteria.__.__.__.__.__.__:Unassigned.__.__.__.__.__.__,
                                               Study,Cohort,Drug.Number)
# values=c(HSCT = "#BC3C294C", GVHD = "#0072B54C", Chemo = "#E187274C")

```

### HSCT

```{r}
beta_hsct_drugnumber<-beta_drugnumber |> filter(Cohort == "HSCT") 
beta_hsct_drugnumber <- beta_hsct_drugnumber %>% filter(Drug.Number!="0")
beta_hsct_drugnumber$Drug.Numberf<-beta_hsct_drugnumber$Drug.Number
beta_hsct_drugnumber[beta_hsct_drugnumber$Drug.Number>5,]$Drug.Numberf=">5"
unique(beta_hsct_drugnumber$Drug.Numberf)
beta_hsct_drugnumber$Drug.Numberf <- factor(beta_hsct_drugnumber$Drug.Numberf,
                                       levels=c("1","2","3","4","5",">5"))
# 透明度：100%-ff;85%-D9;70%-b3;55%-8C;40%-66;25%-40
color_hsct<-c("#EECEC9","#E4B1A9","#DA948A","#D07669","#C65949","#BC3C29")
bray_hsct_drugnumber_plot<-betadrugnumber_plot(beta_hsct_drugnumber,"bray","HSCT",color_hsct,18)
jaccard_hsct_drugnumber_plot<-betadrugnumber_plot(beta_hsct_drugnumber,"jaccard","HSCT",color_hsct,18)

```

### gvhd

```{r}
beta_gvhd_drugnumber<-beta_drugnumber |> filter(Cohort == "GVHD") 
beta_gvhd_drugnumber <- beta_gvhd_drugnumber %>% filter(Drug.Number!="0")
beta_gvhd_drugnumber$Drug.Numberf<-beta_gvhd_drugnumber$Drug.Number
beta_gvhd_drugnumber[beta_gvhd_drugnumber$Drug.Number>5,]$Drug.Numberf=">5"
unique(beta_gvhd_drugnumber$Drug.Numberf)
beta_gvhd_drugnumber$Drug.Numberf <- factor(beta_gvhd_drugnumber$Drug.Numberf,
                                            levels=c("1","3","4","5",">5"))
# 透明度：100%-ff;85%-D9;70%-b3;55%-8C;40%-66;25%-40
color_gvhd<-c("#CCE3F0","#98C5DF","#66AAD3","#338EC4","#0072B5")"#0072B566"
bray_gvhd_drugnumber_plot<-betadrugnumber_plot(beta_gvhd_drugnumber,"bray","GVHD",color_gvhd,17)
jaccard_gvhd_drugnumber_plot<-betadrugnumber_plot(beta_gvhd_drugnumber,"jaccard","GVHD",color_gvhd,17)

```

### chemo

```{r}
beta_chemo_drugnumber<-beta_drugnumber |> filter(Cohort == "Chemo") 
beta_chemo_drugnumber <- beta_chemo_drugnumber %>% filter(Drug.Number!="0")
beta_chemo_drugnumber$Drug.Numberf<-beta_chemo_drugnumber$Drug.Number
beta_chemo_drugnumber[beta_chemo_drugnumber$Drug.Number>5,]$Drug.Numberf=">5"
unique(beta_chemo_drugnumber$Drug.Numberf)
beta_chemo_drugnumber$Drug.Numberf <- factor(beta_chemo_drugnumber$Drug.Numberf,
                                            levels=c("1","2","3","4","5",">5"))
# 透明度：100%-ff;85%-D9;70%-b3;55%-8C;40%-66;25%-40
color_chemo<-c("#00813433","#0081345C","#00813485","#008134AD","#008134D6","#008134FF")
show_col(color_chemo)
color_chemo<-c("#CCE6D6","#A3D2B6","#7ABD95","#52AA75","#299555","#008134")
bray_chemo_drugnumber_plot<-betadrugnumber_plot(beta_chemo_drugnumber,"bray","Chemo",color_chemo,15)
jaccard_chemo_drugnumber_plot<-betadrugnumber_plot(beta_chemo_drugnumber,"jaccard","Chemo",color_chemo,15)

```

### Figure 2 (beta diversity \~ number of used drug)

```{r}
bray_all_drugnumber_plot <- bray_hsct_drugnumber_plot + bray_gvhd_drugnumber_plot + bray_chemo_drugnumber_plot +
  plot_layout(guides = 'collect')

jaccard_all_drugnumber_plot <- jaccard_hsct_drugnumber_plot + jaccard_gvhd_drugnumber_plot + jaccard_chemo_drugnumber_plot +
  plot_layout(guides = 'collect')+plot_annotation(tag_levels = "A")



Figure2<-alpha_multidrug_plot/bray_all_drugnumber_plot+plot_annotation(tag_levels = "A")+plot_layout(heights = c(3, 2))
ggsave("Figure2.pdf",Figure2,units="in", width=14.0, height=10.0, dpi=600,limitsize = FALSE)

ggsave("Figure2s.pdf",jaccard_all_drugnumber_plot,units="in", width=14.0, height=4.8, dpi=600,limitsize = FALSE)

```

## 3.4 beta diversity \~ single drug use analysis (PERMANOVA)

### 3.4.1 HSCT

```{r}
beta_hsct_singledrug<-alpha_control_singledrug |> select(index,Prednisolone:Amphotericin.B,Study,k__Bacteria.__.__.__.__.__.__:k__Bacteria.p__Verrucomicrobia.c__Verrucomicrobiae.o__Verrucomicrobiales.f__Verrucomicrobiaceae.g__Akkermansia.s__muciniphila)

```

#### bray

```{r}
adjust_beta_permanova1<-array(0,dim=c(49,8))
rownames(adjust_beta_permanova1)[2:49]<-colnames(beta_hsct_singledrug)[2:49]
rownames(adjust_beta_permanova1)[1]<-"Number of used medicines"
colnames(adjust_beta_permanova1)<-c("Meds","Df","SumofSqs","R2","F","Pval","FDR","Sig")
adjust_beta_permanova2<-array(0,dim=c(49,8))
rownames(adjust_beta_permanova2)[2:49]<-colnames(beta_hsct_singledrug)[2:49]
rownames(adjust_beta_permanova2)[1]<-"Number of used medicines"
colnames(adjust_beta_permanova2)<-c("Meds","Df","SumofSqs","MeanofSqs","F","Pval","FDR","Sig")
  
for (i in 2:49){
    aa<-beta_hsct_singledrug[complete.cases(beta_hsct_singledrug[,i]),]
    aa$Study<-as.factor(aa$Study)
    aa[,i]<-as.factor(aa[,i])
    if(nrow(aa) != 0){
      # # aa$ave<-apply(aa[,51:ncol(aa)],1,mean)
      # # aa<-aa[which(aa$ave!="0"),]
      # aa$Study<-as.factor(aa$Study)
      # aa[,i]<-as.factor(aa[,i])
      
      m<-lapply(levels(aa$Study),function(j){
        length(levels(aa[(aa$Study==j),i]))==2
      }) |> unlist() |> as.vector() |> sum()
      n<-length(levels(aa[,i]))
      
      aa_matrix001<-aa[,51:ncol(aa)]*0.01
      aa_metadata<-aa[,1:50]
      
      if((n==2)&(m > 1)){

        aa_matrix_batch<-aa_matrix001 |> t()
        colnames(aa_matrix_batch)<-NULL
        rownames(aa_matrix_batch)<-NULL
        colnames(aa_matrix_batch)<-rownames(aa_matrix001)
        rownames(aa_matrix_batch)<-colnames(aa_matrix001)
        
        aa_matrix_batch<-as.matrix(aa_matrix_batch)
        adjust_aa_matrix<-adjust_batch(feature_abd = aa_matrix_batch,
                                       batch = "Study",
                                       covariates = colnames(aa_metadata)[i],
                                       data = aa_metadata
        )$feature_abd_ad
        
        adjust_aa_matrix<-adjust_aa_matrix |> t()
        colnames(adjust_aa_matrix)<-colnames(aa_matrix001)
        rownames(adjust_aa_matrix)<-rownames(aa_matrix001)
        adjust_aa_matrix<-as.data.frame(adjust_aa_matrix)
        adjust_aa<-cbind(aa_metadata,adjust_aa_matrix)
        
        adjust_aa$ave<-apply(adjust_aa[51:ncol(adjust_aa)],1,mean)
        adjust_aa<-adjust_aa[which(adjust_aa$ave!="0"),]
        adjust_aa<-adjust_aa[,-ncol(adjust_aa)]
        # adjust_aa_metadata<-adjust_aa[,1:49]
        
        tmp2<-adonis2(adjust_aa[51:ncol(adjust_aa)] ~ adjust_aa[,i],
                      data = adjust_aa, permutations = 1000, method="bray")
        adjust_beta_permanova1[i,2]<-tmp2$Df[1]
        adjust_beta_permanova1[i,3]<-tmp2$SumOfSqs[1]
        adjust_beta_permanova1[i,4]<-tmp2$R2[1]
        adjust_beta_permanova1[i,5]<-tmp2$F[1]
        adjust_beta_permanova1[i,6]<-tmp2$`Pr(>F)`[1]
        
        tmp3<-vegdist(adjust_aa[51:ncol(adjust_aa)],method="bray", binary=F) |> 
          betadisper(group=adjust_aa[,i]) |> 
          permutest()
        adjust_beta_permanova2[i,2]<-tmp3$tab[1,1]
        adjust_beta_permanova2[i,3]<-tmp3$tab[1,2]
        adjust_beta_permanova2[i,4]<-tmp3$tab[1,3]
        adjust_beta_permanova2[i,5]<-tmp3$tab[1,4]
        adjust_beta_permanova2[i,6]<-tmp3$tab[1,6]
        
      }else if((n==2)&(m == 1)){
        tmp2<-adonis2(aa_matrix001 ~ aa_metadata[,i],
                      data = aa_metadata, permutations = 1000, method="bray")
        adjust_beta_permanova1[i,2]<-tmp2$Df[1]
        adjust_beta_permanova1[i,3]<-tmp2$SumOfSqs[1]
        adjust_beta_permanova1[i,4]<-tmp2$R2[1]
        adjust_beta_permanova1[i,5]<-tmp2$F[1]
        adjust_beta_permanova1[i,6]<-tmp2$`Pr(>F)`[1]
        
        tmp3<-vegdist(aa_matrix001,method="bray", binary=F) |> 
          betadisper(group=aa_metadata[,i]) |> 
          permutest()
        adjust_beta_permanova2[i,2]<-tmp3$tab[1,1]
        adjust_beta_permanova2[i,3]<-tmp3$tab[1,2]
        adjust_beta_permanova2[i,4]<-tmp3$tab[1,3]
        adjust_beta_permanova2[i,5]<-tmp3$tab[1,4]
        adjust_beta_permanova2[i,6]<-tmp3$tab[1,6]
        
      }else{
        adjust_beta_permanova1[i,2]<-NA
        adjust_beta_permanova1[i,3]<-NA
        adjust_beta_permanova1[i,4]<-NA
        adjust_beta_permanova1[i,5]<-NA
        adjust_beta_permanova1[i,6]<-NA
        
        adjust_beta_permanova2[i,2]<-NA
        adjust_beta_permanova2[i,3]<-NA
        adjust_beta_permanova2[i,4]<-NA
        adjust_beta_permanova2[i,5]<-NA
        adjust_beta_permanova2[i,6]<-NA
        
      }
    }else{
      adjust_beta_permanova1[i,2]<-NA
      adjust_beta_permanova1[i,3]<-NA
      adjust_beta_permanova1[i,4]<-NA
      adjust_beta_permanova1[i,5]<-NA
      adjust_beta_permanova1[i,6]<-NA
      
      adjust_beta_permanova2[i,2]<-NA
      adjust_beta_permanova2[i,3]<-NA
      adjust_beta_permanova2[i,4]<-NA
      adjust_beta_permanova2[i,5]<-NA
      adjust_beta_permanova2[i,6]<-NA
      
    }
}


beta_hsct_drugnumber$ave<-apply(beta_hsct_drugnumber[,2:111],1,mean)
beta_hsct_drugnumber<-beta_hsct_drugnumber[which(beta_hsct_drugnumber$ave!="0"),]
beta_hsct_drugnumber_matrix001<-beta_hsct_drugnumber[,2:111]*0.01
beta_hsct_drugnumber_metadata<-beta_hsct_drugnumber[,c(1,112:ncol(beta_hsct_drugnumber))]

tmp2<-adonis2(beta_hsct_drugnumber_matrix001 ~ Drug.Numberf,
              data = beta_hsct_drugnumber_metadata, permutations = 1000, method="bray")
adjust_beta_permanova1[1,2]<-tmp2$Df[1]
adjust_beta_permanova1[1,3]<-tmp2$SumOfSqs[1]
adjust_beta_permanova1[1,4]<-tmp2$R2[1]
adjust_beta_permanova1[1,5]<-tmp2$F[1]
adjust_beta_permanova1[1,6]<-tmp2$`Pr(>F)`[1]

tmp3<-vegdist(beta_hsct_drugnumber_matrix001,method="bray", binary=F) |> 
  betadisper(group=beta_hsct_drugnumber_metadata[,"Drug.Numberf"]) |> 
  permutest()
adjust_beta_permanova2[1,2]<-tmp3$tab[1,1]
adjust_beta_permanova2[1,3]<-tmp3$tab[1,2]
adjust_beta_permanova2[1,4]<-tmp3$tab[1,3]
adjust_beta_permanova2[1,5]<-tmp3$tab[1,4]
adjust_beta_permanova2[1,6]<-tmp3$tab[1,6]

adjust_beta_permanova1<-as.data.frame(adjust_beta_permanova1)
adjust_beta_permanova1$Meds<-rownames(adjust_beta_permanova1)
adjust_beta_permanova1<-adjust_beta_permanova1[!is.na(adjust_beta_permanova1$Pval),]
adjust_beta_permanova1$FDR<-p.adjust(adjust_beta_permanova1$Pval,method = "fdr")
adjust_beta_permanova1$Sig<-ifelse(adjust_beta_permanova1$FDR<0.1,"YES","NO")
  
adjust_beta_permanova2<-as.data.frame(adjust_beta_permanova2)
adjust_beta_permanova2$Meds<-rownames(adjust_beta_permanova2)
adjust_beta_permanova2<-adjust_beta_permanova2[!is.na(adjust_beta_permanova2$Pval),]
adjust_beta_permanova2$FDR<-p.adjust(adjust_beta_permanova2$Pval,method = "fdr")
adjust_beta_permanova2$Sig<-ifelse(adjust_beta_permanova2$FDR<0.1,"YES","NO")
  
beta_hsct_singledrug_permanova1<-adjust_beta_permanova1
beta_hsct_singledrug_permanova2<-adjust_beta_permanova2
write.csv(beta_hsct_singledrug_permanova1,file = "beta_hsct_singledrug_permanova1.csv")
write.csv(beta_hsct_singledrug_permanova2,file = "beta_hsct_singledrug_permanova2.csv")

```

#### jaccard

```{r}
adjust_beta_permanova1<-array(0,dim=c(49,8))
rownames(adjust_beta_permanova1)[2:49]<-colnames(beta_hsct_singledrug)[2:49]
rownames(adjust_beta_permanova1)[1]<-"Number of used medicines"
colnames(adjust_beta_permanova1)<-c("Meds","Df","SumofSqs","R2","F","Pval","FDR","Sig")
adjust_beta_permanova2<-array(0,dim=c(49,8))
rownames(adjust_beta_permanova2)[2:49]<-colnames(beta_hsct_singledrug)[2:49]
rownames(adjust_beta_permanova2)[1]<-"Number of used medicines"
colnames(adjust_beta_permanova2)<-c("Meds","Df","SumofSqs","MeanofSqs","F","Pval","FDR","Sig")

for (i in 2:49){
  aa<-beta_hsct_singledrug[complete.cases(beta_hsct_singledrug[,i]),]
  aa$Study<-as.factor(aa$Study)
  aa[,i]<-as.factor(aa[,i])
  if(nrow(aa) != 0){
    # # aa$ave<-apply(aa[,51:ncol(aa)],1,mean)
    # # aa<-aa[which(aa$ave!="0"),]
    # aa$Study<-as.factor(aa$Study)
    # aa[,i]<-as.factor(aa[,i])
    
    m<-lapply(levels(aa$Study),function(j){
      length(levels(aa[(aa$Study==j),i]))==2
    }) |> unlist() |> as.vector() |> sum()
    n<-length(levels(aa[,i]))
    
    aa_matrix001<-aa[,51:ncol(aa)]*0.01
    aa_metadata<-aa[,1:50]
    
    if((n==2)&(m > 1)){
      
      aa_matrix_batch<-aa_matrix001 |> t()
      colnames(aa_matrix_batch)<-NULL
      rownames(aa_matrix_batch)<-NULL
      colnames(aa_matrix_batch)<-rownames(aa_matrix001)
      rownames(aa_matrix_batch)<-colnames(aa_matrix001)
      
      aa_matrix_batch<-as.matrix(aa_matrix_batch)
      adjust_aa_matrix<-adjust_batch(feature_abd = aa_matrix_batch,
                                     batch = "Study",
                                     covariates = colnames(aa_metadata)[i],
                                     data = aa_metadata
      )$feature_abd_ad
      
      adjust_aa_matrix<-adjust_aa_matrix |> t()
      colnames(adjust_aa_matrix)<-colnames(aa_matrix001)
      rownames(adjust_aa_matrix)<-rownames(aa_matrix001)
      adjust_aa_matrix<-as.data.frame(adjust_aa_matrix)
      adjust_aa<-cbind(aa_metadata,adjust_aa_matrix)
      
      adjust_aa$ave<-apply(adjust_aa[51:ncol(adjust_aa)],1,mean)
      adjust_aa<-adjust_aa[which(adjust_aa$ave!="0"),]
      adjust_aa<-adjust_aa[,-ncol(adjust_aa)]
      # adjust_aa_metadata<-adjust_aa[,1:49]
      
      tmp2<-adonis2(adjust_aa[51:ncol(adjust_aa)] ~ adjust_aa[,i],
                    data = adjust_aa, permutations = 1000, method="jaccard")
      adjust_beta_permanova1[(i-1),2]<-tmp2$Df[1]
      adjust_beta_permanova1[(i-1),3]<-tmp2$SumOfSqs[1]
      adjust_beta_permanova1[(i-1),4]<-tmp2$R2[1]
      adjust_beta_permanova1[(i-1),5]<-tmp2$F[1]
      adjust_beta_permanova1[(i-1),6]<-tmp2$`Pr(>F)`[1]
      
      tmp3<-vegdist(adjust_aa[51:ncol(adjust_aa)],method="jaccard", binary=F) |> 
        betadisper(group=adjust_aa[,i]) |> 
        permutest()
      adjust_beta_permanova2[(i-1),2]<-tmp3$tab[1,1]
      adjust_beta_permanova2[(i-1),3]<-tmp3$tab[1,2]
      adjust_beta_permanova2[(i-1),4]<-tmp3$tab[1,3]
      adjust_beta_permanova2[(i-1),5]<-tmp3$tab[1,4]
      adjust_beta_permanova2[(i-1),6]<-tmp3$tab[1,6]
      
    }else if((n==2)&(m == 1)){
      tmp2<-adonis2(aa_matrix001 ~ aa_metadata[,i],
                    data = aa_metadata, permutations = 1000, method="jaccard")
      adjust_beta_permanova1[(i-1),2]<-tmp2$Df[1]
      adjust_beta_permanova1[(i-1),3]<-tmp2$SumOfSqs[1]
      adjust_beta_permanova1[(i-1),4]<-tmp2$R2[1]
      adjust_beta_permanova1[(i-1),5]<-tmp2$F[1]
      adjust_beta_permanova1[(i-1),6]<-tmp2$`Pr(>F)`[1]
      
      tmp3<-vegdist(aa_matrix001,method="jaccard", binary=F) |> 
        betadisper(group=aa_metadata[,i]) |> 
        permutest()
      adjust_beta_permanova2[(i-1),2]<-tmp3$tab[1,1]
      adjust_beta_permanova2[(i-1),3]<-tmp3$tab[1,2]
      adjust_beta_permanova2[(i-1),4]<-tmp3$tab[1,3]
      adjust_beta_permanova2[(i-1),5]<-tmp3$tab[1,4]
      adjust_beta_permanova2[(i-1),6]<-tmp3$tab[1,6]
      
    }else{
      adjust_beta_permanova1[(i-1),2]<-NA
      adjust_beta_permanova1[(i-1),3]<-NA
      adjust_beta_permanova1[(i-1),4]<-NA
      adjust_beta_permanova1[(i-1),5]<-NA
      adjust_beta_permanova1[(i-1),6]<-NA
      
      adjust_beta_permanova2[(i-1),2]<-NA
      adjust_beta_permanova2[(i-1),3]<-NA
      adjust_beta_permanova2[(i-1),4]<-NA
      adjust_beta_permanova2[(i-1),5]<-NA
      adjust_beta_permanova2[(i-1),6]<-NA
      
    }
  }else{
    adjust_beta_permanova1[(i-1),2]<-NA
    adjust_beta_permanova1[(i-1),3]<-NA
    adjust_beta_permanova1[(i-1),4]<-NA
    adjust_beta_permanova1[(i-1),5]<-NA
    adjust_beta_permanova1[(i-1),6]<-NA
    
    adjust_beta_permanova2[(i-1),2]<-NA
    adjust_beta_permanova2[(i-1),3]<-NA
    adjust_beta_permanova2[(i-1),4]<-NA
    adjust_beta_permanova2[(i-1),5]<-NA
    adjust_beta_permanova2[(i-1),6]<-NA
    
  }
}

beta_hsct_drugnumber$ave<-apply(beta_hsct_drugnumber[,2:111],1,mean)
beta_hsct_drugnumber<-beta_hsct_drugnumber[which(beta_hsct_drugnumber$ave!="0"),]
beta_hsct_drugnumber_matrix001<-beta_hsct_drugnumber[,2:111]*0.01
beta_hsct_drugnumber_metadata<-beta_hsct_drugnumber[,c(1,112:ncol(beta_hsct_drugnumber))]

tmp2<-adonis2(beta_hsct_drugnumber_matrix001 ~ Drug.Numberf,
              data = beta_hsct_drugnumber_metadata, permutations = 1000, method="jaccard")
adjust_beta_permanova1[1,2]<-tmp2$Df[1]
adjust_beta_permanova1[1,3]<-tmp2$SumOfSqs[1]
adjust_beta_permanova1[1,4]<-tmp2$R2[1]
adjust_beta_permanova1[1,5]<-tmp2$F[1]
adjust_beta_permanova1[1,6]<-tmp2$`Pr(>F)`[1]

tmp3<-vegdist(beta_hsct_drugnumber_matrix001,method="jaccard", binary=F) |> 
  betadisper(group=beta_hsct_drugnumber_metadata[,"Drug.Numberf"]) |> 
  permutest()
adjust_beta_permanova2[1,2]<-tmp3$tab[1,1]
adjust_beta_permanova2[1,3]<-tmp3$tab[1,2]
adjust_beta_permanova2[1,4]<-tmp3$tab[1,3]
adjust_beta_permanova2[1,5]<-tmp3$tab[1,4]
adjust_beta_permanova2[1,6]<-tmp3$tab[1,6]

adjust_beta_permanova1<-as.data.frame(adjust_beta_permanova1)
adjust_beta_permanova1$Meds<-rownames(adjust_beta_permanova1)
adjust_beta_permanova1<-adjust_beta_permanova1[!is.na(adjust_beta_permanova1$Pval),]
adjust_beta_permanova1$FDR<-p.adjust(adjust_beta_permanova1$Pval,method = "fdr")
adjust_beta_permanova1$Sig<-ifelse(adjust_beta_permanova1$FDR<0.1,"YES","NO")

adjust_beta_permanova2<-as.data.frame(adjust_beta_permanova2)
adjust_beta_permanova2$Meds<-rownames(adjust_beta_permanova2)
adjust_beta_permanova2<-adjust_beta_permanova2[!is.na(adjust_beta_permanova2$Pval),]
adjust_beta_permanova2$FDR<-p.adjust(adjust_beta_permanova2$Pval,method = "fdr")
adjust_beta_permanova2$Sig<-ifelse(adjust_beta_permanova2$FDR<0.1,"YES","NO")

beta_hsct_singledrug_permanova1_jaccard<-adjust_beta_permanova1
beta_hsct_singledrug_permanova2_jaccard<-adjust_beta_permanova2
write.csv(beta_hsct_singledrug_permanova1_jaccard,file = "beta_hsct_singledrug_permanova1_jaccard.csv")
write.csv(beta_hsct_singledrug_permanova2_jaccard,file = "beta_hsct_singledrug_permanova2_jaccard.csv")

```

### 3.4.2 gvhd

```{r}
beta_gvhd_singledrug<-alpha_gvhd_singledrug |> select(index,Prednisolone:Amphotericin.B,Study,k__Bacteria.__.__.__.__.__.__:k__Bacteria.p__Verrucomicrobia.c__Verrucomicrobiae.o__Verrucomicrobiales.f__Verrucomicrobiaceae.g__Akkermansia.s__muciniphila)

```

#### bray

```{r}
adjust_beta_permanova1<-array(0,dim=c(49,8))
rownames(adjust_beta_permanova1)[2:49]<-colnames(beta_hsct_singledrug)[2:49]
rownames(adjust_beta_permanova1)[1]<-"Number of used medicines"
colnames(adjust_beta_permanova1)<-c("Meds","Df","SumofSqs","R2","F","Pval","FDR","Sig")
adjust_beta_permanova2<-array(0,dim=c(49,8))
rownames(adjust_beta_permanova2)[2:49]<-colnames(beta_hsct_singledrug)[2:49]
rownames(adjust_beta_permanova2)[1]<-"Number of used medicines"
colnames(adjust_beta_permanova2)<-c("Meds","Df","SumofSqs","MeanofSqs","F","Pval","FDR","Sig")

for (i in 2:49){
  aa<-beta_gvhd_singledrug[complete.cases(beta_gvhd_singledrug[,i]),]
  aa$Study<-as.factor(aa$Study)
  aa[,i]<-as.factor(aa[,i])
  if(nrow(aa) != 0){
    # # aa$ave<-apply(aa[,51:ncol(aa)],1,mean)
    # # aa<-aa[which(aa$ave!="0"),]
    # aa$Study<-as.factor(aa$Study)
    # aa[,i]<-as.factor(aa[,i])
    
    m<-lapply(levels(aa$Study),function(j){
      length(levels(aa[(aa$Study==j),i]))==2
    }) |> unlist() |> as.vector() |> sum()
    n<-length(levels(aa[,i]))
    
    aa_matrix001<-aa[,51:ncol(aa)]*0.01
    aa_metadata<-aa[,1:50]
    
    if((n==2)&(m > 1)){
      
      aa_matrix_batch<-aa_matrix001 |> t()
      colnames(aa_matrix_batch)<-NULL
      rownames(aa_matrix_batch)<-NULL
      colnames(aa_matrix_batch)<-rownames(aa_matrix001)
      rownames(aa_matrix_batch)<-colnames(aa_matrix001)
      
      aa_matrix_batch<-as.matrix(aa_matrix_batch)
      adjust_aa_matrix<-adjust_batch(feature_abd = aa_matrix_batch,
                                     batch = "Study",
                                     covariates = colnames(aa_metadata)[i],
                                     data = aa_metadata
      )$feature_abd_ad
      
      adjust_aa_matrix<-adjust_aa_matrix |> t()
      colnames(adjust_aa_matrix)<-colnames(aa_matrix001)
      rownames(adjust_aa_matrix)<-rownames(aa_matrix001)
      adjust_aa_matrix<-as.data.frame(adjust_aa_matrix)
      adjust_aa<-cbind(aa_metadata,adjust_aa_matrix)
      
      adjust_aa$ave<-apply(adjust_aa[51:ncol(adjust_aa)],1,mean)
      adjust_aa<-adjust_aa[which(adjust_aa$ave!="0"),]
      adjust_aa<-adjust_aa[,-ncol(adjust_aa)]
      # adjust_aa_metadata<-adjust_aa[,1:49]
      
      tmp2<-adonis2(adjust_aa[51:ncol(adjust_aa)] ~ adjust_aa[,i],
                    data = adjust_aa, permutations = 1000, method="bray")
      adjust_beta_permanova1[(i-1),2]<-tmp2$Df[1]
      adjust_beta_permanova1[(i-1),3]<-tmp2$SumOfSqs[1]
      adjust_beta_permanova1[(i-1),4]<-tmp2$R2[1]
      adjust_beta_permanova1[(i-1),5]<-tmp2$F[1]
      adjust_beta_permanova1[(i-1),6]<-tmp2$`Pr(>F)`[1]
      
      tmp3<-vegdist(adjust_aa[51:ncol(adjust_aa)],method="bray", binary=F) |> 
        betadisper(group=adjust_aa[,i]) |> 
        permutest()
      adjust_beta_permanova2[(i-1),2]<-tmp3$tab[1,1]
      adjust_beta_permanova2[(i-1),3]<-tmp3$tab[1,2]
      adjust_beta_permanova2[(i-1),4]<-tmp3$tab[1,3]
      adjust_beta_permanova2[(i-1),5]<-tmp3$tab[1,4]
      adjust_beta_permanova2[(i-1),6]<-tmp3$tab[1,6]
      
    }else if((n==2)&(m == 1)){
      tmp2<-adonis2(aa_matrix001 ~ aa_metadata[,i],
                    data = aa_metadata, permutations = 1000, method="bray")
      adjust_beta_permanova1[(i-1),2]<-tmp2$Df[1]
      adjust_beta_permanova1[(i-1),3]<-tmp2$SumOfSqs[1]
      adjust_beta_permanova1[(i-1),4]<-tmp2$R2[1]
      adjust_beta_permanova1[(i-1),5]<-tmp2$F[1]
      adjust_beta_permanova1[(i-1),6]<-tmp2$`Pr(>F)`[1]
      
      tmp3<-vegdist(aa_matrix001,method="bray", binary=F) |> 
        betadisper(group=aa_metadata[,i]) |> 
        permutest()
      adjust_beta_permanova2[(i-1),2]<-tmp3$tab[1,1]
      adjust_beta_permanova2[(i-1),3]<-tmp3$tab[1,2]
      adjust_beta_permanova2[(i-1),4]<-tmp3$tab[1,3]
      adjust_beta_permanova2[(i-1),5]<-tmp3$tab[1,4]
      adjust_beta_permanova2[(i-1),6]<-tmp3$tab[1,6]
      
    }else{
      adjust_beta_permanova1[(i-1),2]<-NA
      adjust_beta_permanova1[(i-1),3]<-NA
      adjust_beta_permanova1[(i-1),4]<-NA
      adjust_beta_permanova1[(i-1),5]<-NA
      adjust_beta_permanova1[(i-1),6]<-NA
      
      adjust_beta_permanova2[(i-1),2]<-NA
      adjust_beta_permanova2[(i-1),3]<-NA
      adjust_beta_permanova2[(i-1),4]<-NA
      adjust_beta_permanova2[(i-1),5]<-NA
      adjust_beta_permanova2[(i-1),6]<-NA
      
    }
  }else{
    adjust_beta_permanova1[(i-1),2]<-NA
    adjust_beta_permanova1[(i-1),3]<-NA
    adjust_beta_permanova1[(i-1),4]<-NA
    adjust_beta_permanova1[(i-1),5]<-NA
    adjust_beta_permanova1[(i-1),6]<-NA
    
    adjust_beta_permanova2[(i-1),2]<-NA
    adjust_beta_permanova2[(i-1),3]<-NA
    adjust_beta_permanova2[(i-1),4]<-NA
    adjust_beta_permanova2[(i-1),5]<-NA
    adjust_beta_permanova2[(i-1),6]<-NA
    
  }
}

beta_gvhd_drugnumber$ave<-apply(beta_gvhd_drugnumber[,2:111],1,mean)
beta_gvhd_drugnumber<-beta_gvhd_drugnumber[which(beta_gvhd_drugnumber$ave!="0"),]
beta_gvhd_drugnumber_matrix001<-beta_gvhd_drugnumber[,2:111]*0.01
beta_gvhd_drugnumber_metadata<-beta_gvhd_drugnumber[,c(1,112:ncol(beta_gvhd_drugnumber))]

tmp2<-adonis2(beta_gvhd_drugnumber_matrix001 ~ Drug.Numberf,
              data = beta_gvhd_drugnumber_metadata, permutations = 1000, method="bray")
adjust_beta_permanova1[1,2]<-tmp2$Df[1]
adjust_beta_permanova1[1,3]<-tmp2$SumOfSqs[1]
adjust_beta_permanova1[1,4]<-tmp2$R2[1]
adjust_beta_permanova1[1,5]<-tmp2$F[1]
adjust_beta_permanova1[1,6]<-tmp2$`Pr(>F)`[1]

tmp3<-vegdist(beta_gvhd_drugnumber_matrix001,method="bray", binary=F) |> 
  betadisper(group=beta_gvhd_drugnumber_metadata[,"Drug.Numberf"]) |> 
  permutest()
adjust_beta_permanova2[1,2]<-tmp3$tab[1,1]
adjust_beta_permanova2[1,3]<-tmp3$tab[1,2]
adjust_beta_permanova2[1,4]<-tmp3$tab[1,3]
adjust_beta_permanova2[1,5]<-tmp3$tab[1,4]
adjust_beta_permanova2[1,6]<-tmp3$tab[1,6]

adjust_beta_permanova1<-as.data.frame(adjust_beta_permanova1)
adjust_beta_permanova1$Meds<-rownames(adjust_beta_permanova1)
adjust_beta_permanova1<-adjust_beta_permanova1[!is.na(adjust_beta_permanova1$Pval),]
adjust_beta_permanova1$FDR<-p.adjust(adjust_beta_permanova1$Pval,method = "fdr")
adjust_beta_permanova1$Sig<-ifelse(adjust_beta_permanova1$FDR<0.1,"YES","NO")

adjust_beta_permanova2<-as.data.frame(adjust_beta_permanova2)
adjust_beta_permanova2$Meds<-rownames(adjust_beta_permanova2)
adjust_beta_permanova2<-adjust_beta_permanova2[!is.na(adjust_beta_permanova2$Pval),]
adjust_beta_permanova2$FDR<-p.adjust(adjust_beta_permanova2$Pval,method = "fdr")
adjust_beta_permanova2$Sig<-ifelse(adjust_beta_permanova2$FDR<0.1,"YES","NO")

beta_gvhd_singledrug_permanova1<-adjust_beta_permanova1
beta_gvhd_singledrug_permanova2<-adjust_beta_permanova2
write.csv(beta_gvhd_singledrug_permanova1,file = "beta_gvhd_singledrug_permanova1.csv")
write.csv(beta_gvhd_singledrug_permanova2,file = "beta_gvhd_singledrug_permanova2.csv")

```

#### jaccard

```{r}
adjust_beta_permanova1<-array(0,dim=c(49,8))
rownames(adjust_beta_permanova1)[2:49]<-colnames(beta_hsct_singledrug)[2:49]
rownames(adjust_beta_permanova1)[1]<-"Number of used medicines"
colnames(adjust_beta_permanova1)<-c("Meds","Df","SumofSqs","R2","F","Pval","FDR","Sig")
adjust_beta_permanova2<-array(0,dim=c(49,8))
rownames(adjust_beta_permanova2)[2:49]<-colnames(beta_hsct_singledrug)[2:49]
rownames(adjust_beta_permanova2)[1]<-"Number of used medicines"
colnames(adjust_beta_permanova2)<-c("Meds","Df","SumofSqs","MeanofSqs","F","Pval","FDR","Sig")

for (i in 2:49){
  aa<-beta_gvhd_singledrug[complete.cases(beta_gvhd_singledrug[,i]),]
  aa$Study<-as.factor(aa$Study)
  aa[,i]<-as.factor(aa[,i])
  if(nrow(aa) != 0){
    # # aa$ave<-apply(aa[,51:ncol(aa)],1,mean)
    # # aa<-aa[which(aa$ave!="0"),]
    # aa$Study<-as.factor(aa$Study)
    # aa[,i]<-as.factor(aa[,i])
    
    m<-lapply(levels(aa$Study),function(j){
      length(levels(aa[(aa$Study==j),i]))==2
    }) |> unlist() |> as.vector() |> sum()
    n<-length(levels(aa[,i]))
    
    aa_matrix001<-aa[,51:ncol(aa)]*0.01
    aa_metadata<-aa[,1:50]
    
    if((n==2)&(m > 1)){
      
      aa_matrix_batch<-aa_matrix001 |> t()
      colnames(aa_matrix_batch)<-NULL
      rownames(aa_matrix_batch)<-NULL
      colnames(aa_matrix_batch)<-rownames(aa_matrix001)
      rownames(aa_matrix_batch)<-colnames(aa_matrix001)
      
      aa_matrix_batch<-as.matrix(aa_matrix_batch)
      adjust_aa_matrix<-adjust_batch(feature_abd = aa_matrix_batch,
                                     batch = "Study",
                                     covariates = colnames(aa_metadata)[i],
                                     data = aa_metadata
      )$feature_abd_ad
      
      adjust_aa_matrix<-adjust_aa_matrix |> t()
      colnames(adjust_aa_matrix)<-colnames(aa_matrix001)
      rownames(adjust_aa_matrix)<-rownames(aa_matrix001)
      adjust_aa_matrix<-as.data.frame(adjust_aa_matrix)
      adjust_aa<-cbind(aa_metadata,adjust_aa_matrix)
      
      adjust_aa$ave<-apply(adjust_aa[51:ncol(adjust_aa)],1,mean)
      adjust_aa<-adjust_aa[which(adjust_aa$ave!="0"),]
      adjust_aa<-adjust_aa[,-ncol(adjust_aa)]
      # adjust_aa_metadata<-adjust_aa[,1:49]
      
      tmp2<-adonis2(adjust_aa[51:ncol(adjust_aa)] ~ adjust_aa[,i],
                    data = adjust_aa, permutations = 1000, method="jaccard")
      adjust_beta_permanova1[(i-1),2]<-tmp2$Df[1]
      adjust_beta_permanova1[(i-1),3]<-tmp2$SumOfSqs[1]
      adjust_beta_permanova1[(i-1),4]<-tmp2$R2[1]
      adjust_beta_permanova1[(i-1),5]<-tmp2$F[1]
      adjust_beta_permanova1[(i-1),6]<-tmp2$`Pr(>F)`[1]
      
      tmp3<-vegdist(adjust_aa[51:ncol(adjust_aa)],method="jaccard", binary=F) |> 
        betadisper(group=adjust_aa[,i]) |> 
        permutest()
      adjust_beta_permanova2[(i-1),2]<-tmp3$tab[1,1]
      adjust_beta_permanova2[(i-1),3]<-tmp3$tab[1,2]
      adjust_beta_permanova2[(i-1),4]<-tmp3$tab[1,3]
      adjust_beta_permanova2[(i-1),5]<-tmp3$tab[1,4]
      adjust_beta_permanova2[(i-1),6]<-tmp3$tab[1,6]
      
    }else if((n==2)&(m == 1)){
      tmp2<-adonis2(aa_matrix001 ~ aa_metadata[,i],
                    data = aa_metadata, permutations = 1000, method="jaccard")
      adjust_beta_permanova1[(i-1),2]<-tmp2$Df[1]
      adjust_beta_permanova1[(i-1),3]<-tmp2$SumOfSqs[1]
      adjust_beta_permanova1[(i-1),4]<-tmp2$R2[1]
      adjust_beta_permanova1[(i-1),5]<-tmp2$F[1]
      adjust_beta_permanova1[(i-1),6]<-tmp2$`Pr(>F)`[1]
      
      tmp3<-vegdist(aa_matrix001,method="jaccard", binary=F) |> 
        betadisper(group=aa_metadata[,i]) |> 
        permutest()
      adjust_beta_permanova2[(i-1),2]<-tmp3$tab[1,1]
      adjust_beta_permanova2[(i-1),3]<-tmp3$tab[1,2]
      adjust_beta_permanova2[(i-1),4]<-tmp3$tab[1,3]
      adjust_beta_permanova2[(i-1),5]<-tmp3$tab[1,4]
      adjust_beta_permanova2[(i-1),6]<-tmp3$tab[1,6]
      
    }else{
      adjust_beta_permanova1[(i-1),2]<-NA
      adjust_beta_permanova1[(i-1),3]<-NA
      adjust_beta_permanova1[(i-1),4]<-NA
      adjust_beta_permanova1[(i-1),5]<-NA
      adjust_beta_permanova1[(i-1),6]<-NA
      
      adjust_beta_permanova2[(i-1),2]<-NA
      adjust_beta_permanova2[(i-1),3]<-NA
      adjust_beta_permanova2[(i-1),4]<-NA
      adjust_beta_permanova2[(i-1),5]<-NA
      adjust_beta_permanova2[(i-1),6]<-NA
      
    }
  }else{
    adjust_beta_permanova1[(i-1),2]<-NA
    adjust_beta_permanova1[(i-1),3]<-NA
    adjust_beta_permanova1[(i-1),4]<-NA
    adjust_beta_permanova1[(i-1),5]<-NA
    adjust_beta_permanova1[(i-1),6]<-NA
    
    adjust_beta_permanova2[(i-1),2]<-NA
    adjust_beta_permanova2[(i-1),3]<-NA
    adjust_beta_permanova2[(i-1),4]<-NA
    adjust_beta_permanova2[(i-1),5]<-NA
    adjust_beta_permanova2[(i-1),6]<-NA
    
  }
}

beta_gvhd_drugnumber$ave<-apply(beta_gvhd_drugnumber[,2:111],1,mean)
beta_gvhd_drugnumber<-beta_gvhd_drugnumber[which(beta_gvhd_drugnumber$ave!="0"),]
beta_gvhd_drugnumber_matrix001<-beta_gvhd_drugnumber[,2:111]*0.01
beta_gvhd_drugnumber_metadata<-beta_gvhd_drugnumber[,c(1,112:ncol(beta_gvhd_drugnumber))]

tmp2<-adonis2(beta_gvhd_drugnumber_matrix001 ~ Drug.Numberf,
              data = beta_gvhd_drugnumber_metadata, permutations = 1000, method="jaccard")
adjust_beta_permanova1[1,2]<-tmp2$Df[1]
adjust_beta_permanova1[1,3]<-tmp2$SumOfSqs[1]
adjust_beta_permanova1[1,4]<-tmp2$R2[1]
adjust_beta_permanova1[1,5]<-tmp2$F[1]
adjust_beta_permanova1[1,6]<-tmp2$`Pr(>F)`[1]

tmp3<-vegdist(beta_gvhd_drugnumber_matrix001,method="jaccard", binary=F) |> 
  betadisper(group=beta_gvhd_drugnumber_metadata[,"Drug.Numberf"]) |> 
  permutest()
adjust_beta_permanova2[1,2]<-tmp3$tab[1,1]
adjust_beta_permanova2[1,3]<-tmp3$tab[1,2]
adjust_beta_permanova2[1,4]<-tmp3$tab[1,3]
adjust_beta_permanova2[1,5]<-tmp3$tab[1,4]
adjust_beta_permanova2[1,6]<-tmp3$tab[1,6]


adjust_beta_permanova1<-as.data.frame(adjust_beta_permanova1)
adjust_beta_permanova1$Meds<-rownames(adjust_beta_permanova1)
adjust_beta_permanova1<-adjust_beta_permanova1[!is.na(adjust_beta_permanova1$Pval),]
adjust_beta_permanova1$FDR<-p.adjust(adjust_beta_permanova1$Pval,method = "fdr")
adjust_beta_permanova1$Sig<-ifelse(adjust_beta_permanova1$FDR<0.1,"YES","NO")

adjust_beta_permanova2<-as.data.frame(adjust_beta_permanova2)
adjust_beta_permanova2$Meds<-rownames(adjust_beta_permanova2)
adjust_beta_permanova2<-adjust_beta_permanova2[!is.na(adjust_beta_permanova2$Pval),]
adjust_beta_permanova2$FDR<-p.adjust(adjust_beta_permanova2$Pval,method = "fdr")
adjust_beta_permanova2$Sig<-ifelse(adjust_beta_permanova2$FDR<0.1,"YES","NO")

beta_gvhd_singledrug_permanova1_jaccard<-adjust_beta_permanova1
beta_gvhd_singledrug_permanova2_jaccard<-adjust_beta_permanova2
write.csv(beta_gvhd_singledrug_permanova1_jaccard,file = "beta_gvhd_singledrug_permanova1_jaccard.csv")
write.csv(beta_gvhd_singledrug_permanova2_jaccard,file = "beta_gvhd_singledrug_permanova2_jaccard.csv")

```

### 3.4.3 CHEMO

```{r}
beta_chemo_singledrug<-alpha_chemo_singledrug |> select(index,Prednisolone:Amphotericin.B,Study,k__Bacteria.__.__.__.__.__.__:Unassigned.__.__.__.__.__.__)

```

#### bray

```{r}
adjust_beta_permanova1<-array(0,dim=c(49,8))
rownames(adjust_beta_permanova1)[2:49]<-colnames(beta_hsct_singledrug)[2:49]
rownames(adjust_beta_permanova1)[1]<-"Number of used medicines"
colnames(adjust_beta_permanova1)<-c("Meds","Df","SumofSqs","R2","F","Pval","FDR","Sig")
adjust_beta_permanova2<-array(0,dim=c(49,8))
rownames(adjust_beta_permanova2)[2:49]<-colnames(beta_hsct_singledrug)[2:49]
rownames(adjust_beta_permanova2)[1]<-"Number of used medicines"
colnames(adjust_beta_permanova2)<-c("Meds","Df","SumofSqs","MeanofSqs","F","Pval","FDR","Sig")

for (i in 2:49){
  aa<-beta_chemo_singledrug[complete.cases(beta_chemo_singledrug[,i]),]
  aa$Study<-as.factor(aa$Study)
  aa[,i]<-as.factor(aa[,i])
  if(nrow(aa) != 0){
    # # aa$ave<-apply(aa[,51:ncol(aa)],1,mean)
    # # aa<-aa[which(aa$ave!="0"),]
    # aa$Study<-as.factor(aa$Study)
    # aa[,i]<-as.factor(aa[,i])
    
    m<-lapply(levels(aa$Study),function(j){
      length(levels(aa[(aa$Study==j),i]))==2
    }) |> unlist() |> as.vector() |> sum()
    n<-length(levels(aa[,i]))
    
    aa_matrix001<-aa[,51:ncol(aa)]*0.01
    aa_metadata<-aa[,1:50]
    
    if((n==2)&(m > 1)){
      
      aa_matrix_batch<-aa_matrix001 |> t()
      colnames(aa_matrix_batch)<-NULL
      rownames(aa_matrix_batch)<-NULL
      colnames(aa_matrix_batch)<-rownames(aa_matrix001)
      rownames(aa_matrix_batch)<-colnames(aa_matrix001)
      
      aa_matrix_batch<-as.matrix(aa_matrix_batch)
      adjust_aa_matrix<-adjust_batch(feature_abd = aa_matrix_batch,
                                     batch = "Study",
                                     covariates = colnames(aa_metadata)[i],
                                     data = aa_metadata
      )$feature_abd_ad
      
      adjust_aa_matrix<-adjust_aa_matrix |> t()
      colnames(adjust_aa_matrix)<-colnames(aa_matrix001)
      rownames(adjust_aa_matrix)<-rownames(aa_matrix001)
      adjust_aa_matrix<-as.data.frame(adjust_aa_matrix)
      adjust_aa<-cbind(aa_metadata,adjust_aa_matrix)
      
      adjust_aa$ave<-apply(adjust_aa[51:ncol(adjust_aa)],1,mean)
      adjust_aa<-adjust_aa[which(adjust_aa$ave!="0"),]
      adjust_aa<-adjust_aa[,-ncol(adjust_aa)]
      # adjust_aa_metadata<-adjust_aa[,1:49]
      
      tmp2<-adonis2(adjust_aa[51:ncol(adjust_aa)] ~ adjust_aa[,i],
                    data = adjust_aa, permutations = 1000, method="bray")
      adjust_beta_permanova1[(i-1),2]<-tmp2$Df[1]
      adjust_beta_permanova1[(i-1),3]<-tmp2$SumOfSqs[1]
      adjust_beta_permanova1[(i-1),4]<-tmp2$R2[1]
      adjust_beta_permanova1[(i-1),5]<-tmp2$F[1]
      adjust_beta_permanova1[(i-1),6]<-tmp2$`Pr(>F)`[1]
      
      tmp3<-vegdist(adjust_aa[51:ncol(adjust_aa)],method="bray", binary=F) |> 
        betadisper(group=adjust_aa[,i]) |> 
        permutest()
      adjust_beta_permanova2[(i-1),2]<-tmp3$tab[1,1]
      adjust_beta_permanova2[(i-1),3]<-tmp3$tab[1,2]
      adjust_beta_permanova2[(i-1),4]<-tmp3$tab[1,3]
      adjust_beta_permanova2[(i-1),5]<-tmp3$tab[1,4]
      adjust_beta_permanova2[(i-1),6]<-tmp3$tab[1,6]
      
    }else if((n==2)&(m == 1)){
      tmp2<-adonis2(aa_matrix001 ~ aa_metadata[,i],
                    data = aa_metadata, permutations = 1000, method="bray")
      adjust_beta_permanova1[(i-1),2]<-tmp2$Df[1]
      adjust_beta_permanova1[(i-1),3]<-tmp2$SumOfSqs[1]
      adjust_beta_permanova1[(i-1),4]<-tmp2$R2[1]
      adjust_beta_permanova1[(i-1),5]<-tmp2$F[1]
      adjust_beta_permanova1[(i-1),6]<-tmp2$`Pr(>F)`[1]
      
      tmp3<-vegdist(aa_matrix001,method="bray", binary=F) |> 
        betadisper(group=aa_metadata[,i]) |> 
        permutest()
      adjust_beta_permanova2[(i-1),2]<-tmp3$tab[1,1]
      adjust_beta_permanova2[(i-1),3]<-tmp3$tab[1,2]
      adjust_beta_permanova2[(i-1),4]<-tmp3$tab[1,3]
      adjust_beta_permanova2[(i-1),5]<-tmp3$tab[1,4]
      adjust_beta_permanova2[(i-1),6]<-tmp3$tab[1,6]
      
    }else{
      adjust_beta_permanova1[(i-1),2]<-NA
      adjust_beta_permanova1[(i-1),3]<-NA
      adjust_beta_permanova1[(i-1),4]<-NA
      adjust_beta_permanova1[(i-1),5]<-NA
      adjust_beta_permanova1[(i-1),6]<-NA
      
      adjust_beta_permanova2[(i-1),2]<-NA
      adjust_beta_permanova2[(i-1),3]<-NA
      adjust_beta_permanova2[(i-1),4]<-NA
      adjust_beta_permanova2[(i-1),5]<-NA
      adjust_beta_permanova2[(i-1),6]<-NA
      
    }
  }else{
    adjust_beta_permanova1[(i-1),2]<-NA
    adjust_beta_permanova1[(i-1),3]<-NA
    adjust_beta_permanova1[(i-1),4]<-NA
    adjust_beta_permanova1[(i-1),5]<-NA
    adjust_beta_permanova1[(i-1),6]<-NA
    
    adjust_beta_permanova2[(i-1),2]<-NA
    adjust_beta_permanova2[(i-1),3]<-NA
    adjust_beta_permanova2[(i-1),4]<-NA
    adjust_beta_permanova2[(i-1),5]<-NA
    adjust_beta_permanova2[(i-1),6]<-NA
    
  }
}

beta_chemo_drugnumber$ave<-apply(beta_chemo_drugnumber[,2:111],1,mean)
beta_chemo_drugnumber<-beta_chemo_drugnumber[which(beta_chemo_drugnumber$ave!="0"),]
beta_chemo_drugnumber_matrix001<-beta_chemo_drugnumber[,2:111]*0.01
beta_chemo_drugnumber_metadata<-beta_chemo_drugnumber[,c(1,112:ncol(beta_chemo_drugnumber))]

tmp2<-adonis2(beta_chemo_drugnumber_matrix001 ~ Drug.Numberf,
              data = beta_chemo_drugnumber_metadata, permutations = 1000, method="bray")
adjust_beta_permanova1[1,2]<-tmp2$Df[1]
adjust_beta_permanova1[1,3]<-tmp2$SumOfSqs[1]
adjust_beta_permanova1[1,4]<-tmp2$R2[1]
adjust_beta_permanova1[1,5]<-tmp2$F[1]
adjust_beta_permanova1[1,6]<-tmp2$`Pr(>F)`[1]

tmp3<-vegdist(beta_chemo_drugnumber_matrix001,method="bray", binary=F) |> 
  betadisper(group=beta_chemo_drugnumber_metadata[,"Drug.Numberf"]) |> 
  permutest()
adjust_beta_permanova2[1,2]<-tmp3$tab[1,1]
adjust_beta_permanova2[1,3]<-tmp3$tab[1,2]
adjust_beta_permanova2[1,4]<-tmp3$tab[1,3]
adjust_beta_permanova2[1,5]<-tmp3$tab[1,4]
adjust_beta_permanova2[1,6]<-tmp3$tab[1,6]

adjust_beta_permanova1<-as.data.frame(adjust_beta_permanova1)
adjust_beta_permanova1$Meds<-rownames(adjust_beta_permanova1)
adjust_beta_permanova1<-adjust_beta_permanova1[!is.na(adjust_beta_permanova1$Pval),]
adjust_beta_permanova1$FDR<-p.adjust(adjust_beta_permanova1$Pval,method = "fdr")
adjust_beta_permanova1$Sig<-ifelse(adjust_beta_permanova1$FDR<0.1,"YES","NO")

adjust_beta_permanova2<-as.data.frame(adjust_beta_permanova2)
adjust_beta_permanova2$Meds<-rownames(adjust_beta_permanova2)
adjust_beta_permanova2<-adjust_beta_permanova2[!is.na(adjust_beta_permanova2$Pval),]
adjust_beta_permanova2$FDR<-p.adjust(adjust_beta_permanova2$Pval,method = "fdr")
adjust_beta_permanova2$Sig<-ifelse(adjust_beta_permanova2$FDR<0.1,"YES","NO")

beta_chemo_singledrug_permanova1<-adjust_beta_permanova1
beta_chemo_singledrug_permanova2<-adjust_beta_permanova2

write.csv(beta_chemo_singledrug_permanova1,file = "beta_chemo_singledrug_permanova1.csv")
write.csv(beta_chemo_singledrug_permanova2,file = "beta_chemo_singledrug_permanova2.csv")

```

#### jaccard

```{r}
adjust_beta_permanova1<-array(0,dim=c(49,8))
rownames(adjust_beta_permanova1)[2:49]<-colnames(beta_hsct_singledrug)[2:49]
rownames(adjust_beta_permanova1)[1]<-"Number of used medicines"
colnames(adjust_beta_permanova1)<-c("Meds","Df","SumofSqs","R2","F","Pval","FDR","Sig")
adjust_beta_permanova2<-array(0,dim=c(49,8))
rownames(adjust_beta_permanova2)[2:49]<-colnames(beta_hsct_singledrug)[2:49]
rownames(adjust_beta_permanova2)[1]<-"Number of used medicines"
colnames(adjust_beta_permanova2)<-c("Meds","Df","SumofSqs","MeanofSqs","F","Pval","FDR","Sig")

for (i in 2:49){
  aa<-beta_chemo_singledrug[complete.cases(beta_chemo_singledrug[,i]),]
  aa$Study<-as.factor(aa$Study)
  aa[,i]<-as.factor(aa[,i])
  if(nrow(aa) != 0){
    # # aa$ave<-apply(aa[,51:ncol(aa)],1,mean)
    # # aa<-aa[which(aa$ave!="0"),]
    # aa$Study<-as.factor(aa$Study)
    # aa[,i]<-as.factor(aa[,i])
    
    m<-lapply(levels(aa$Study),function(j){
      length(levels(aa[(aa$Study==j),i]))==2
    }) |> unlist() |> as.vector() |> sum()
    n<-length(levels(aa[,i]))
    
    aa_matrix001<-aa[,51:ncol(aa)]*0.01
    aa_metadata<-aa[,1:50]
    
    if((n==2)&(m > 1)){
      
      aa_matrix_batch<-aa_matrix001 |> t()
      colnames(aa_matrix_batch)<-NULL
      rownames(aa_matrix_batch)<-NULL
      colnames(aa_matrix_batch)<-rownames(aa_matrix001)
      rownames(aa_matrix_batch)<-colnames(aa_matrix001)
      
      aa_matrix_batch<-as.matrix(aa_matrix_batch)
      adjust_aa_matrix<-adjust_batch(feature_abd = aa_matrix_batch,
                                     batch = "Study",
                                     covariates = colnames(aa_metadata)[i],
                                     data = aa_metadata
      )$feature_abd_ad
      
      adjust_aa_matrix<-adjust_aa_matrix |> t()
      colnames(adjust_aa_matrix)<-colnames(aa_matrix001)
      rownames(adjust_aa_matrix)<-rownames(aa_matrix001)
      adjust_aa_matrix<-as.data.frame(adjust_aa_matrix)
      adjust_aa<-cbind(aa_metadata,adjust_aa_matrix)
      
      adjust_aa$ave<-apply(adjust_aa[51:ncol(adjust_aa)],1,mean)
      adjust_aa<-adjust_aa[which(adjust_aa$ave!="0"),]
      adjust_aa<-adjust_aa[,-ncol(adjust_aa)]
      # adjust_aa_metadata<-adjust_aa[,1:49]
      
      tmp2<-adonis2(adjust_aa[51:ncol(adjust_aa)] ~ adjust_aa[,i],
                    data = adjust_aa, permutations = 1000, method="jaccard")
      adjust_beta_permanova1[(i-1),2]<-tmp2$Df[1]
      adjust_beta_permanova1[(i-1),3]<-tmp2$SumOfSqs[1]
      adjust_beta_permanova1[(i-1),4]<-tmp2$R2[1]
      adjust_beta_permanova1[(i-1),5]<-tmp2$F[1]
      adjust_beta_permanova1[(i-1),6]<-tmp2$`Pr(>F)`[1]
      
      tmp3<-vegdist(adjust_aa[51:ncol(adjust_aa)],method="jaccard", binary=F) |> 
        betadisper(group=adjust_aa[,i]) |> 
        permutest()
      adjust_beta_permanova2[(i-1),2]<-tmp3$tab[1,1]
      adjust_beta_permanova2[(i-1),3]<-tmp3$tab[1,2]
      adjust_beta_permanova2[(i-1),4]<-tmp3$tab[1,3]
      adjust_beta_permanova2[(i-1),5]<-tmp3$tab[1,4]
      adjust_beta_permanova2[(i-1),6]<-tmp3$tab[1,6]
      
    }else if((n==2)&(m == 1)){
      tmp2<-adonis2(aa_matrix001 ~ aa_metadata[,i],
                    data = aa_metadata, permutations = 1000, method="jaccard")
      adjust_beta_permanova1[(i-1),2]<-tmp2$Df[1]
      adjust_beta_permanova1[(i-1),3]<-tmp2$SumOfSqs[1]
      adjust_beta_permanova1[(i-1),4]<-tmp2$R2[1]
      adjust_beta_permanova1[(i-1),5]<-tmp2$F[1]
      adjust_beta_permanova1[(i-1),6]<-tmp2$`Pr(>F)`[1]
      
      tmp3<-vegdist(aa_matrix001,method="jaccard", binary=F) |> 
        betadisper(group=aa_metadata[,i]) |> 
        permutest()
      adjust_beta_permanova2[(i-1),2]<-tmp3$tab[1,1]
      adjust_beta_permanova2[(i-1),3]<-tmp3$tab[1,2]
      adjust_beta_permanova2[(i-1),4]<-tmp3$tab[1,3]
      adjust_beta_permanova2[(i-1),5]<-tmp3$tab[1,4]
      adjust_beta_permanova2[(i-1),6]<-tmp3$tab[1,6]
      
    }else{
      adjust_beta_permanova1[(i-1),2]<-NA
      adjust_beta_permanova1[(i-1),3]<-NA
      adjust_beta_permanova1[(i-1),4]<-NA
      adjust_beta_permanova1[(i-1),5]<-NA
      adjust_beta_permanova1[(i-1),6]<-NA
      
      adjust_beta_permanova2[(i-1),2]<-NA
      adjust_beta_permanova2[(i-1),3]<-NA
      adjust_beta_permanova2[(i-1),4]<-NA
      adjust_beta_permanova2[(i-1),5]<-NA
      adjust_beta_permanova2[(i-1),6]<-NA
      
    }
  }else{
    adjust_beta_permanova1[(i-1),2]<-NA
    adjust_beta_permanova1[(i-1),3]<-NA
    adjust_beta_permanova1[(i-1),4]<-NA
    adjust_beta_permanova1[(i-1),5]<-NA
    adjust_beta_permanova1[(i-1),6]<-NA
    
    adjust_beta_permanova2[(i-1),2]<-NA
    adjust_beta_permanova2[(i-1),3]<-NA
    adjust_beta_permanova2[(i-1),4]<-NA
    adjust_beta_permanova2[(i-1),5]<-NA
    adjust_beta_permanova2[(i-1),6]<-NA
    
  }
}

beta_chemo_drugnumber$ave<-apply(beta_chemo_drugnumber[,2:111],1,mean)
beta_chemo_drugnumber<-beta_chemo_drugnumber[which(beta_chemo_drugnumber$ave!="0"),]
beta_chemo_drugnumber_matrix001<-beta_chemo_drugnumber[,2:111]*0.01
beta_chemo_drugnumber_metadata<-beta_chemo_drugnumber[,c(1,112:ncol(beta_chemo_drugnumber))]

tmp2<-adonis2(beta_chemo_drugnumber_matrix001 ~ Drug.Numberf,
              data = beta_chemo_drugnumber_metadata, permutations = 1000, method="jaccard")
adjust_beta_permanova1[1,2]<-tmp2$Df[1]
adjust_beta_permanova1[1,3]<-tmp2$SumOfSqs[1]
adjust_beta_permanova1[1,4]<-tmp2$R2[1]
adjust_beta_permanova1[1,5]<-tmp2$F[1]
adjust_beta_permanova1[1,6]<-tmp2$`Pr(>F)`[1]

tmp3<-vegdist(beta_chemo_drugnumber_matrix001,method="jaccard", binary=F) |> 
  betadisper(group=beta_chemo_drugnumber_metadata[,"Drug.Numberf"]) |> 
  permutest()
adjust_beta_permanova2[1,2]<-tmp3$tab[1,1]
adjust_beta_permanova2[1,3]<-tmp3$tab[1,2]
adjust_beta_permanova2[1,4]<-tmp3$tab[1,3]
adjust_beta_permanova2[1,5]<-tmp3$tab[1,4]
adjust_beta_permanova2[1,6]<-tmp3$tab[1,6]

adjust_beta_permanova1<-as.data.frame(adjust_beta_permanova1)
adjust_beta_permanova1$Meds<-rownames(adjust_beta_permanova1)
adjust_beta_permanova1<-adjust_beta_permanova1[!is.na(adjust_beta_permanova1$Pval),]
adjust_beta_permanova1$FDR<-p.adjust(adjust_beta_permanova1$Pval,method = "fdr")
adjust_beta_permanova1$Sig<-ifelse(adjust_beta_permanova1$FDR<0.1,"YES","NO")

adjust_beta_permanova2<-as.data.frame(adjust_beta_permanova2)
adjust_beta_permanova2$Meds<-rownames(adjust_beta_permanova2)
adjust_beta_permanova2<-adjust_beta_permanova2[!is.na(adjust_beta_permanova2$Pval),]
adjust_beta_permanova2$FDR<-p.adjust(adjust_beta_permanova2$Pval,method = "fdr")
adjust_beta_permanova2$Sig<-ifelse(adjust_beta_permanova2$FDR<0.1,"YES","NO")

beta_chemo_singledrug_permanova1_jaccard<-adjust_beta_permanova1
beta_chemo_singledrug_permanova2_jaccard<-adjust_beta_permanova2

write.csv(beta_chemo_singledrug_permanova1_jaccard,file = "beta_chemo_singledrug_permanova1_jaccard.csv")
write.csv(beta_chemo_singledrug_permanova2_jaccard,file = "beta_chemo_singledrug_permanova2_jaccard.csv")

```

# 4.Linear models

> Packages used in building linear models

```{R}

library(ade4)
library(tidyverse)
library(data.table)
library(VIM)
library(car)
library(magrittr)
library(plyr)
library (meta)
library(reshape2)
library(ggplot2)
library(patchwork)

```

## 4.1 HSCT Cohort

```{r}
control<-all_reg_HSCT
control_0<-control
control_0[is.na(control_0)]<-0

sum_abundance<-apply(control_0[,2:1417],1,sum)

control_rel<-array(0,dim=c(838,1416))
for (j in 1:1416) {
  for (i in 1:838){
    control_rel[i,j]<-control_0[i,j+1]/sum_abundance[i]*100
  }
}
control_rel <- as.data.frame(control_rel)

sample<-control_0$index
species<-colnames(control_0)
colnames(control_rel)<-species[2:1417]
rownames(control_rel)<-sample

control_rel_fil_10 <- filtering_taxonomy(control_rel,0.00000001,10)

summary_taxnomy_filtering_10<- read.table("summary_taxonomy_filtering.txt",sep="\t",header=T)
filtered_taxonomy_10<-read.table("filtered_taxonomy.txt",sep="\t",header=T)

filtered_taxonomy_10_1<-cbind(filtered_taxonomy_10,index=rownames(filtered_taxonomy_10))

metadata<-as.data.frame(cbind(index=sample,control[,1418:1482]))
control_reg<-merge(filtered_taxonomy_10_1,metadata,by="index")
rownames(control_reg)<-sample
control_reg_name<-control_reg
colnames(control_reg_name)[115]<-"amplicon"
colnames(control_reg_name)[116]<-"sex"
colnames(control_reg_name)[119]<-"scs"
colnames(control_reg_name)[120]<-"DrH"
```

### single drug analysis

```{r}
number_hsct<-array(length(unique(control_reg_name$subject)),dim=c(104,48))
rownames(number_hsct)<-colnames(control_reg_name[,2:105])
colnames(number_hsct)<-colnames(control_reg_name[,121:168])
nonzero_hsct<-array(NA,dim=c(104,48))
rownames(nonzero_hsct)<-colnames(control_reg_name[,2:105])
colnames(nonzero_hsct)<-colnames(control_reg_name[,121:168])
users_hsct<-array(NA,dim=c(104,48))
rownames(users_hsct)<-colnames(control_reg_name[,2:105])
colnames(users_hsct)<-colnames(control_reg_name[,121:168])
nonusers_hsct<-array(NA,dim=c(104,48))
rownames(nonusers_hsct)<-colnames(control_reg_name[,2:105])
colnames(nonusers_hsct)<-colnames(control_reg_name[,121:168])

estimate<-array(0,dim=c(104,48))
rownames(estimate)<-colnames(control_reg_name[,2:105])
colnames(estimate)<-colnames(control_reg_name[,121:168])
Std_Error<-array(0,dim=c(104,48))
rownames(Std_Error)<-colnames(control_reg_name[,2:105])
colnames(Std_Error)<-colnames(control_reg_name[,121:168])
t_value<-array(0,dim=c(104,48))
rownames(t_value)<-colnames(control_reg_name[,2:105])
colnames(t_value)<-colnames(control_reg_name[,121:168])
p_value<-array(0,dim=c(104,48))
rownames(p_value)<-colnames(control_reg_name[,2:105])
colnames(p_value)<-colnames(control_reg_name[,121:168])



for (x in 1:104) {
  for (y in 1:48){
    
    genus_name<-colnames(control_reg_name[,2:105])
    drug_name<-colnames(control_reg_name[,121:168])
    
    temp<-select(control_reg_name,genus_name[x],subject,flag,time,amplicon,sex,age,diagnosis,scs,DrH,drug_name[y])
    temp$species<-temp[,1]
    temp_nonzero<- temp |> filter(species!=0) 
    nonzero_hsct[x,y] <-length(unique(temp_nonzero$subject))
    temp$drug<-temp[,(ncol(temp)-1)]
    temp_users<-temp[complete.cases(temp$drug),] |> filter(drug==1) 
    users_hsct[x,y] <- length(unique(temp_users$subject))
    temp_nonusers<-temp[complete.cases(temp$drug),] |> filter(drug==0) 
    nonusers_hsct[x,y] <- length(unique(temp_nonusers$subject))
    
    temp<-temp[complete.cases(temp),]
    drug<-temp[,ncol(temp)]
    if(length(unique(drug))==2){
      if(length(unique(temp$time))!=1 & length(unique(temp$amplicon))!=1 & length(unique(temp$sex))!=1 & length(unique(temp$scs))!=1){
        temp1<-lm(temp[,1] ~ time + amplicon + sex + age + diagnosis + scs + DrH + drug,temp)
        if (rownames(summary(temp1)$coefficients)[nrow(summary(temp1)$coefficients)] == "drug"){
          summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),1]->estimate[x,y]
          summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),2]->Std_Error[x,y]
          summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),3]->t_value[x,y]
          summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),4]->p_value[x,y]
        }else{
          estimate[x,y]<-NA
          Std_Error[x,y]<-NA
          t_value[x,y]<-NA
          p_value[x,y]<-NA
        }
      }else if(length(unique(temp$time))==1 & length(unique(temp$amplicon))!=1 & length(unique(temp$sex))!=1 & length(unique(temp$scs))!=1){
        temp1<-lm(temp[,1] ~ amplicon + sex + age + diagnosis + scs + DrH + drug,temp)
        if (rownames(summary(temp1)$coefficients)[nrow(summary(temp1)$coefficients)] == "drug"){
          estimate[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),1]
          Std_Error[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),2]
          t_value[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),3]
          p_value[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),4]
        }else{
          estimate[x,y]<-NA
          Std_Error[x,y]<-NA
          t_value[x,y]<-NA
          p_value[x,y]<-NA
        }
      }else if(length(unique(temp$time))!=1 & length(unique(temp$amplicon))==1 & length(unique(temp$sex))!=1 & length(unique(temp$scs))!=1){
        temp1<-lm(temp[,1] ~ time + sex + age + diagnosis + scs + DrH + drug,temp)
        if (rownames(summary(temp1)$coefficients)[nrow(summary(temp1)$coefficients)] == "drug"){
          estimate[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),1]
          Std_Error[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),2]
          t_value[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),3]
          p_value[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),4]
        }else{
          estimate[x,y]<-NA
          Std_Error[x,y]<-NA
          t_value[x,y]<-NA
          p_value[x,y]<-NA
        }
      }else if(length(unique(temp$time))!=1 & length(unique(temp$amplicon))!=1 & length(unique(temp$sex))==1 & length(unique(temp$scs))!=1){
        temp1<-lm(temp[,1] ~ time  + age + diagnosis + scs + DrH + drug,temp)
        if (rownames(summary(temp1)$coefficients)[nrow(summary(temp1)$coefficients)] == "drug"){
          estimate[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),1]
          Std_Error[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),2]
          t_value[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),3]
          p_value[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),4]
        }else{
          estimate[x,y]<-NA
          Std_Error[x,y]<-NA
          t_value[x,y]<-NA
          p_value[x,y]<-NA
        }
      }else if(length(unique(temp$time))!=1 & length(unique(temp$amplicon))!=1 & length(unique(temp$sex))!=1 & length(unique(temp$scs))==1){
        temp1<-lm(temp[,1] ~ time + sex + age + diagnosis + DrH + drug,temp)
        if (rownames(summary(temp1)$coefficients)[nrow(summary(temp1)$coefficients)] == "drug"){
          estimate[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),1]
          Std_Error[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),2]
          t_value[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),3]
          p_value[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),4]
        }else{
          estimate[x,y]<-NA
          Std_Error[x,y]<-NA
          t_value[x,y]<-NA
          p_value[x,y]<-NA
        }
      }else if(length(unique(temp$time))==1 & length(unique(temp$amplicon))==1 & length(unique(temp$sex))!=1 & length(unique(temp$scs))!=1){
        temp1<-lm(temp[,1] ~ sex + age + diagnosis + scs + DrH + drug,temp)
        if (rownames(summary(temp1)$coefficients)[nrow(summary(temp1)$coefficients)] == "drug"){
          estimate[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),1]
          Std_Error[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),2]
          t_value[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),3]
          p_value[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),4]
        }else{
          estimate[x,y]<-NA
          Std_Error[x,y]<-NA
          t_value[x,y]<-NA
          p_value[x,y]<-NA
        }
      }else if(length(unique(temp$time))==1 & length(unique(temp$amplicon))!=1 & length(unique(temp$sex))==1 & length(unique(temp$scs))!=1){
        temp1<-lm(temp[,1] ~ amplicon + age + diagnosis + scs + DrH + drug,temp)
        if (rownames(summary(temp1)$coefficients)[nrow(summary(temp1)$coefficients)] == "drug"){
          estimate[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),1]
          Std_Error[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),2]
          t_value[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),3]
          p_value[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),4]
        }else{
          estimate[x,y]<-NA
          Std_Error[x,y]<-NA
          t_value[x,y]<-NA
          p_value[x,y]<-NA
        }
      }else if(length(unique(temp$time))==1 & length(unique(temp$amplicon))!=1 & length(unique(temp$sex))!=1 & length(unique(temp$scs))==1){
        temp1<-lm(temp[,1] ~ amplicon + sex + age + diagnosis + DrH + drug,temp)
        if (rownames(summary(temp1)$coefficients)[nrow(summary(temp1)$coefficients)] == "drug"){
          estimate[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),1]
          Std_Error[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),2]
          t_value[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),3]
          p_value[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),4]
        }else{
          estimate[x,y]<-NA
          Std_Error[x,y]<-NA
          t_value[x,y]<-NA
          p_value[x,y]<-NA
        }
      }else if(length(unique(temp$time))!=1 & length(unique(temp$amplicon))==1 & length(unique(temp$sex))==1 & length(unique(temp$scs))!=1){
        temp1<-lm(temp[,1] ~ time + age + diagnosis + scs + DrH + drug,temp)
        if (rownames(summary(temp1)$coefficients)[nrow(summary(temp1)$coefficients)] == "drug"){
          estimate[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),1]
          Std_Error[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),2]
          t_value[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),3]
          p_value[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),4]
        }else{
          estimate[x,y]<-NA
          Std_Error[x,y]<-NA
          t_value[x,y]<-NA
          p_value[x,y]<-NA
        }
      }else if(length(unique(temp$time))!=1 & length(unique(temp$amplicon))==1 & length(unique(temp$sex))!=1 & length(unique(temp$scs))==1){
        temp1<-lm(temp[,1] ~ time + sex + age + diagnosis + DrH + drug,temp)
        if (rownames(summary(temp1)$coefficients)[nrow(summary(temp1)$coefficients)] == "drug"){
          estimate[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),1]
          Std_Error[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),2]
          t_value[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),3]
          p_value[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),4]
        }else{
          estimate[x,y]<-NA
          Std_Error[x,y]<-NA
          t_value[x,y]<-NA
          p_value[x,y]<-NA
        }
      }else if(length(unique(temp$time))!=1 & length(unique(temp$amplicon))!=1 & length(unique(temp$sex))==1 & length(unique(temp$scs))==1){
        temp1<-lm(temp[,1] ~ amplicon + time + age + diagnosis + DrH + drug,temp)
        if (rownames(summary(temp1)$coefficients)[nrow(summary(temp1)$coefficients)] == "drug"){
          estimate[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),1]
          Std_Error[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),2]
          t_value[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),3]
          p_value[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),4]
        }else{
          estimate[x,y]<-NA
          Std_Error[x,y]<-NA
          t_value[x,y]<-NA
          p_value[x,y]<-NA
        }
      }else if(length(unique(temp$time))==1 & length(unique(temp$amplicon))==1 & length(unique(temp$sex))==1 & length(unique(temp$scs))!=1){
        temp1<-lm(temp[,1] ~ age + diagnosis + scs + DrH + drug,temp)
        if (rownames(summary(temp1)$coefficients)[nrow(summary(temp1)$coefficients)] == "drug"){
          estimate[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),1]
          Std_Error[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),2]
          t_value[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),3]
          p_value[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),4]
        }else{
          estimate[x,y]<-NA
          Std_Error[x,y]<-NA
          t_value[x,y]<-NA
          p_value[x,y]<-NA
        }
      }else if(length(unique(temp$time))==1 & length(unique(temp$amplicon))==1 & length(unique(temp$sex))!=1 & length(unique(temp$scs))==1){
        temp1<-lm(temp[,1] ~ sex + age + diagnosis + DrH + drug,temp)
        if (rownames(summary(temp1)$coefficients)[nrow(summary(temp1)$coefficients)] == "drug"){
          estimate[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),1]
          Std_Error[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),2]
          t_value[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),3]
          p_value[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),4]
        }else{
          estimate[x,y]<-NA
          Std_Error[x,y]<-NA
          t_value[x,y]<-NA
          p_value[x,y]<-NA
        }
      }else if(length(unique(temp$time))==1 & length(unique(temp$amplicon))!=1 & length(unique(temp$sex))==1 & length(unique(temp$scs))==1){
        temp1<-lm(temp[,1] ~ amplicon + age + diagnosis + DrH + drug,temp)
        if (rownames(summary(temp1)$coefficients)[nrow(summary(temp1)$coefficients)] == "drug"){
          estimate[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),1]
          Std_Error[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),2]
          t_value[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),3]
          p_value[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),4]
        }else{
          estimate[x,y]<-NA
          Std_Error[x,y]<-NA
          t_value[x,y]<-NA
          p_value[x,y]<-NA
        }
      }else if(length(unique(temp$time))!=1 & length(unique(temp$amplicon))==1 & length(unique(temp$sex))==1 & length(unique(temp$scs))==1){
        temp1<-lm(temp[,1] ~ time + age + diagnosis + DrH + drug,temp)
        if (rownames(summary(temp1)$coefficients)[nrow(summary(temp1)$coefficients)] == "drug"){
          estimate[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),1]
          Std_Error[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),2]
          t_value[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),3]
          p_value[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),4]
        }else{
          estimate[x,y]<-NA
          Std_Error[x,y]<-NA
          t_value[x,y]<-NA
          p_value[x,y]<-NA
        }
      }else{
        temp1<-lm(temp[,1] ~ age + diagnosis + DrH + drug,temp)
        if (rownames(summary(temp1)$coefficients)[nrow(summary(temp1)$coefficients)] == "drug"){
          estimate[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),1]
          Std_Error[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),2]
          t_value[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),3]
          p_value[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),4]
        }else{
          estimate[x,y]<-NA
          Std_Error[x,y]<-NA
          t_value[x,y]<-NA
          p_value[x,y]<-NA
        }
      }
    }else{
      estimate[x,y]<-NA
      Std_Error[x,y]<-NA
      t_value[x,y]<-NA
      p_value[x,y]<-NA
    }
  }
}

number_hsct<-as.data.frame(number_hsct)
nonzero_hsct<-as.data.frame(nonzero_hsct)
users_hsct<-as.data.frame(users_hsct)
nonusers_hsct<-as.data.frame(nonusers_hsct)

estimate<-as.data.frame(estimate)
p_value<-as.data.frame(p_value)
Std_Error<-as.data.frame(Std_Error)
t_value<-as.data.frame(t_value)

```

### multiple drug analysis

```{r}
number_mul_hsct<-array(length(unique(control_reg_name$subject)),dim=c(104,47))
rownames(number_mul_hsct)<-colnames(control_reg_name[,2:105])
colnames(number_mul_hsct)<-colnames(control_reg_name[,c(121:142,144:168)])
nonzero_mul_hsct<-array(NA,dim=c(104,47))
rownames(nonzero_mul_hsct)<-colnames(control_reg_name[,2:105])
colnames(nonzero_mul_hsct)<-colnames(control_reg_name[,c(121:142,144:168)])
users_mul_hsct<-array(NA,dim=c(104,47))
rownames(users_mul_hsct)<-colnames(control_reg_name[,2:105])
colnames(users_mul_hsct)<-colnames(control_reg_name[,c(121:142,144:168)])
nonusers_mul_hsct<-array(NA,dim=c(104,47))
rownames(nonusers_mul_hsct)<-colnames(control_reg_name[,2:105])
colnames(nonusers_mul_hsct)<-colnames(control_reg_name[,c(121:142,144:168)])


estimate_mul<-array(0,dim=c(104,47))
rownames(estimate_mul)<-colnames(control_reg_name[,2:105])
colnames(estimate_mul)<-colnames(control_reg_name[,c(121:142,144:168)])
Std_Error_mul<-array(0,dim=c(104,47))
rownames(Std_Error_mul)<-colnames(control_reg_name[,2:105])
colnames(Std_Error_mul)<-colnames(control_reg_name[,c(121:142,144:168)])
t_value_mul<-array(0,dim=c(104,47))
rownames(t_value_mul)<-colnames(control_reg_name[,2:105])
colnames(t_value_mul)<-colnames(control_reg_name[,c(121:142,144:168)])
p_value_mul<-array(0,dim=c(104,47))
rownames(p_value_mul)<-colnames(control_reg_name[,2:105])
colnames(p_value_mul)<-colnames(control_reg_name[,c(121:142,144:168)])


for (x in 1:104) {
  for (y in 1:47){
    genus_name<-colnames(control_reg_name[,2:105])
    drug_name<-colnames(control_reg_name[,c(121:142,144:168)])
    
    temp_mul<-select(control_reg_name,genus_name[x],subject,flag,time,amplicon,sex,age,diagnosis,scs,DrH,Antibiotic,drug_name[y])
    temp_mul$species<-temp_mul[,1]
    temp_mul_nonzero<- temp_mul |> filter(species!=0) 
    nonzero_mul_hsct[x,y] <-length(unique(temp_mul_nonzero$subject))
    temp_mul$drug<-temp_mul[,(ncol(temp_mul)-1)]
    temp_mul_users<-temp_mul[complete.cases(temp_mul$drug),] |> filter(drug==1) 
    users_mul_hsct[x,y] <- length(unique(temp_mul_users$subject))
    temp_mul_nonusers<-temp_mul[complete.cases(temp_mul$drug),] |> filter(drug==0) 
    nonusers_mul_hsct[x,y] <- length(unique(temp_mul_nonusers$subject))
    
    
    temp_mul<-temp_mul[complete.cases(temp_mul),]
    drug<-temp_mul[,ncol(temp_mul)]
    Antibiotic<-temp_mul[,"Antibiotic"]
    if(length(unique(drug))==2 & length(unique(Antibiotic))==2){
      if(length(unique(temp_mul$time))!=1 & length(unique(temp_mul$amplicon))!=1 & length(unique(temp_mul$sex))!=1 & length(unique(temp_mul$scs))!=1){
        temp_mul1<-lm(temp_mul[,1] ~ time + amplicon + sex + age + diagnosis + scs + DrH + Antibiotic +drug,temp_mul)
        if (rownames(summary(temp_mul1)$coefficients)[nrow(summary(temp_mul1)$coefficients)] == "drug"){
          summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),1]->estimate_mul[x,y]
          summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),2]->Std_Error_mul[x,y]
          summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),3]->t_value_mul[x,y]
          summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),4]->p_value_mul[x,y]
        }else{
          estimate_mul[x,y]<-NA
          Std_Error_mul[x,y]<-NA
          t_value_mul[x,y]<-NA
          p_value_mul[x,y]<-NA
        }
      }else if(length(unique(temp_mul$time))==1 & length(unique(temp_mul$amplicon))!=1 & length(unique(temp_mul$sex))!=1 & length(unique(temp_mul$scs))!=1){
        temp_mul1<-lm(temp_mul[,1] ~ amplicon + sex + age + diagnosis + scs + DrH + Antibiotic+ drug,temp_mul)
        if (rownames(summary(temp_mul1)$coefficients)[nrow(summary(temp_mul1)$coefficients)] == "drug"){
          estimate_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),1]
          Std_Error_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),2]
          t_value_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),3]
          p_value_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),4]
        }else{
          estimate_mul[x,y]<-NA
          Std_Error_mul[x,y]<-NA
          t_value_mul[x,y]<-NA
          p_value_mul[x,y]<-NA
        }
      }else if(length(unique(temp_mul$time))!=1 & length(unique(temp_mul$amplicon))==1 & length(unique(temp_mul$sex))!=1 & length(unique(temp_mul$scs))!=1){
        temp_mul1<-lm(temp_mul[,1] ~ time + sex + age + diagnosis + scs + DrH + Antibiotic+ drug,temp_mul)
        if (rownames(summary(temp_mul1)$coefficients)[nrow(summary(temp_mul1)$coefficients)] == "drug"){
          estimate_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),1]
          Std_Error_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),2]
          t_value_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),3]
          p_value_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),4]
        }else{
          estimate_mul[x,y]<-NA
          Std_Error_mul[x,y]<-NA
          t_value_mul[x,y]<-NA
          p_value_mul[x,y]<-NA
        }
      }else if(length(unique(temp_mul$time))!=1 & length(unique(temp_mul$amplicon))!=1 & length(unique(temp_mul$sex))==1 & length(unique(temp_mul$scs))!=1){
        temp_mul1<-lm(temp_mul[,1] ~ time  + age + diagnosis + scs + DrH + Antibiotic+ drug,temp_mul)
        if (rownames(summary(temp_mul1)$coefficients)[nrow(summary(temp_mul1)$coefficients)] == "drug"){
          estimate_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),1]
          Std_Error_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),2]
          t_value_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),3]
          p_value_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),4]
        }else{
          estimate_mul[x,y]<-NA
          Std_Error_mul[x,y]<-NA
          t_value_mul[x,y]<-NA
          p_value_mul[x,y]<-NA
        }
      }else if(length(unique(temp_mul$time))!=1 & length(unique(temp_mul$amplicon))!=1 & length(unique(temp_mul$sex))!=1 & length(unique(temp_mul$scs))==1){
        temp_mul1<-lm(temp_mul[,1] ~ time + sex + age + diagnosis + DrH+ Antibiotic + drug,temp_mul)
        if (rownames(summary(temp_mul1)$coefficients)[nrow(summary(temp_mul1)$coefficients)] == "drug"){
          estimate_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),1]
          Std_Error_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),2]
          t_value_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),3]
          p_value_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),4]
        }else{
          estimate_mul[x,y]<-NA
          Std_Error_mul[x,y]<-NA
          t_value_mul[x,y]<-NA
          p_value_mul[x,y]<-NA
        }
      }else if(length(unique(temp_mul$time))==1 & length(unique(temp_mul$amplicon))==1 & length(unique(temp_mul$sex))!=1 & length(unique(temp_mul$scs))!=1){
        temp_mul1<-lm(temp_mul[,1] ~ sex + age + diagnosis + scs + DrH + Antibiotic+ drug,temp_mul)
        if (rownames(summary(temp_mul1)$coefficients)[nrow(summary(temp_mul1)$coefficients)] == "drug"){
          estimate_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),1]
          Std_Error_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),2]
          t_value_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),3]
          p_value_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),4]
        }else{
          estimate_mul[x,y]<-NA
          Std_Error_mul<-NA
          t_value_mul[x,y]<-NA
          p_value_mul[x,y]<-NA
        }
      }else if(length(unique(temp_mul$time))==1 & length(unique(temp_mul$amplicon))!=1 & length(unique(temp_mul$sex))==1 & length(unique(temp_mul$scs))!=1){
        temp_mul1<-lm(temp_mul[,1] ~ amplicon + age + diagnosis + scs + DrH + Antibiotic+ drug,temp_mul)
        if (rownames(summary(temp_mul1)$coefficients)[nrow(summary(temp_mul1)$coefficients)] == "drug"){
          estimate_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),1]
          Std_Error_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),2]
          t_value_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),3]
          p_value_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),4]
        }else{
          estimate_mul[x,y]<-NA
          Std_Error_mul[x,y]<-NA
          t_value_mul[x,y]<-NA
          p_value_mul[x,y]<-NA
        }
      }else if(length(unique(temp_mul$time))==1 & length(unique(temp_mul$amplicon))!=1 & length(unique(temp_mul$sex))!=1 & length(unique(temp_mul$scs))==1){
        temp_mul1<-lm(temp_mul[,1] ~ amplicon + sex + age + diagnosis + DrH + Antibiotic+ drug,temp_mul)
        if (rownames(summary(temp_mul1)$coefficients)[nrow(summary(temp_mul1)$coefficients)] == "drug"){
          estimate_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),1]
          Std_Error_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),2]
          t_value_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),3]
          p_value_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),4]
        }else{
          estimate_mul[x,y]<-NA
          Std_Error_mul[x,y]<-NA
          t_value_mul[x,y]<-NA
          p_value_mul[x,y]<-NA
        }
      }else if(length(unique(temp_mul$time))!=1 & length(unique(temp_mul$amplicon))==1 & length(unique(temp_mul$sex))==1 & length(unique(temp_mul$scs))!=1){
        temp_mul1<-lm(temp_mul[,1] ~ time + age + diagnosis + scs + DrH + Antibiotic+ drug,temp_mul)
        if (rownames(summary(temp_mul1)$coefficients)[nrow(summary(temp_mul1)$coefficients)] == "drug"){
          estimate_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),1]
          Std_Error_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),2]
          t_value_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),3]
          p_value_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),4]
        }else{
          estimate_mul[x,y]<-NA
          Std_Error_mul[x,y]<-NA
          t_value_mul[x,y]<-NA
          p_value_mul[x,y]<-NA
        }
      }else if(length(unique(temp_mul$time))!=1 & length(unique(temp_mul$amplicon))==1 & length(unique(temp_mul$sex))!=1 & length(unique(temp_mul$scs))==1){
        temp_mul1<-lm(temp_mul[,1] ~ time + sex + age + diagnosis + DrH + Antibiotic+ drug,temp_mul)
        if (rownames(summary(temp_mul1)$coefficients)[nrow(summary(temp_mul1)$coefficients)] == "drug"){
          estimate_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),1]
          Std_Error_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),2]
          t_value_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),3]
          p_value_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),4]
        }else{
          estimate_mul[x,y]<-NA
          Std_Error_mul[x,y]<-NA
          t_value_mul[x,y]<-NA
          p_value_mul[x,y]<-NA
        }
      }else if(length(unique(temp_mul$time))!=1 & length(unique(temp_mul$amplicon))!=1 & length(unique(temp_mul$sex))==1 & length(unique(temp_mul$scs))==1){
        temp_mul1<-lm(temp_mul[,1] ~ amplicon + time + age + diagnosis + DrH + Antibiotic+ drug,temp_mul)
        if (rownames(summary(temp_mul1)$coefficients)[nrow(summary(temp_mul1)$coefficients)] == "drug"){
          estimate_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),1]
          Std_Error_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),2]
          t_value_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),3]
          p_value_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),4]
        }else{
          estimate_mul[x,y]<-NA
          Std_Error_mul[x,y]<-NA
          t_value_mul[x,y]<-NA
          p_value_mul[x,y]<-NA
        }
      }else if(length(unique(temp_mul$time))==1 & length(unique(temp_mul$amplicon))==1 & length(unique(temp_mul$sex))==1 & length(unique(temp_mul$scs))!=1){
        temp_mul1<-lm(temp_mul[,1] ~ age + diagnosis + scs + DrH + Antibiotic+ drug,temp_mul)
        if (rownames(summary(temp_mul1)$coefficients)[nrow(summary(temp_mul1)$coefficients)] == "drug"){
          estimate_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),1]
          Std_Error_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),2]
          t_value_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),3]
          p_value_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),4]
        }else{
          estimate_mul[x,y]<-NA
          Std_Error_mul[x,y]<-NA
          t_value_mul[x,y]<-NA
          p_value_mul[x,y]<-NA
        }
      }else if(length(unique(temp_mul$time))==1 & length(unique(temp_mul$amplicon))==1 & length(unique(temp_mul$sex))!=1 & length(unique(temp_mul$scs))==1){
        temp_mul1<-lm(temp_mul[,1] ~ sex + age + diagnosis + DrH + Antibiotic+ drug,temp_mul)
        if (rownames(summary(temp_mul1)$coefficients)[nrow(summary(temp_mul1)$coefficients)] == "drug"){
          estimate_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),1]
          Std_Error_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),2]
          t_value_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),3]
          p_value_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),4]
        }else{
          estimate_mul[x,y]<-NA
          Std_Error_mul[x,y]<-NA
          t_value_mul[x,y]<-NA
          p_value_mul[x,y]<-NA
        }
      }else if(length(unique(temp_mul$time))==1 & length(unique(temp_mul$amplicon))!=1 & length(unique(temp_mul$sex))==1 & length(unique(temp_mul$scs))==1){
        temp_mul1<-lm(temp_mul[,1] ~ amplicon + age + diagnosis + DrH + Antibiotic+ drug,temp_mul)
        if (rownames(summary(temp_mul1)$coefficients)[nrow(summary(temp_mul1)$coefficients)] == "drug"){
          estimate_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),1]
          Std_Error_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),2]
          t_value_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),3]
          p_value_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),4]
        }else{
          estimate_mul[x,y]<-NA
          Std_Error_mul[x,y]<-NA
          t_value_mul[x,y]<-NA
          p_value_mul[x,y]<-NA
        }
      }else if(length(unique(temp_mul$time))!=1 & length(unique(temp_mul$amplicon))==1 & length(unique(temp_mul$sex))==1 & length(unique(temp_mul$scs))==1){
        temp_mul1<-lm(temp_mul[,1] ~ time + age + diagnosis + DrH + Antibiotic+ drug,temp_mul)
        if (rownames(summary(temp_mul1)$coefficients)[nrow(summary(temp_mul1)$coefficients)] == "drug"){
          estimate_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),1]
          Std_Error_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),2]
          t_value_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),3]
          p_value_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),4]
        }else{
          estimate_mul[x,y]<-NA
          Std_Error_mul[x,y]<-NA
          t_value_mul[x,y]<-NA
          p_value_mul[x,y]<-NA
        }
      }else{
        temp_mul1<-lm(temp_mul[,1] ~ age + diagnosis + DrH + Antibiotic+ drug,temp_mul)
        if (rownames(summary(temp_mul1)$coefficients)[nrow(summary(temp_mul1)$coefficients)] == "drug"){
          estimate_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),1]
          Std_Error_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),2]
          t_value_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),3]
          p_value_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),4]
        }else{
          estimate_mul[x,y]<-NA
          Std_Error_mul[x,y]<-NA
          t_value_mul[x,y]<-NA
          p_value_mul[x,y]<-NA
        }
      }
    }else{
      estimate_mul[x,y]<-NA
      Std_Error_mul[x,y]<-NA
      t_value_mul[x,y]<-NA
      p_value_mul[x,y]<-NA
    }
  }
}

number_mul_hsct<-as.data.frame(number_mul_hsct)
nonzero_mul_hsct<-as.data.frame(nonzero_mul_hsct)
users_mul_hsct<-as.data.frame(users_mul_hsct)
nonusers_mul_hsct<-as.data.frame(nonusers_mul_hsct)

estimate_mul<-as.data.frame(estimate_mul)
p_value_mul<-as.data.frame(p_value_mul)
Std_Error_mul<-as.data.frame(Std_Error_mul)
t_value_mul<-as.data.frame(t_value_mul)

```

## 4.2 GVHD Cohort

```{r}

gvhd_0<-all_reg_GVHD
gvhd_0[is.na(gvhd_0)]<-0

sum_abundance_gvhd<-apply(gvhd_0[,2:1417],1,sum)

gvhd_rel<-array(0,dim=c(715,1416))
for (j in 1:1416) {
  for (i in 1:715){
    gvhd_rel[i,j]<-gvhd_0[i,j+1]/sum_abundance_gvhd[i]*100
  }
}

gvhd_rel <- as.data.frame(gvhd_rel)
gvhd_rel[is.na(gvhd_rel)]<-0

sample_gvhd<-gvhd_0$index
species_gvhd<-colnames(gvhd_0)
colnames(gvhd_rel)<-species_gvhd[2:1417]
rownames(gvhd_rel)<-sample_gvhd

gvhd_rel_fil_10 <- filtering_taxonomy(gvhd_rel,0.00000001,10)

summary_taxnomy_filtering_10_gvhd<- read.table("summary_taxonomy_filtering.txt",sep="\t",header=T)
filtered_taxonomy_10_gvhd<-read.table("filtered_taxonomy.txt",sep="\t",header=T)

filtered_taxonomy_10_1_gvhd<-cbind(filtered_taxonomy_10_gvhd,index=rownames(filtered_taxonomy_10_gvhd))

metadata_gvhd<-as.data.frame(cbind(gvhd[,1418:1482],index=sample_gvhd))
gvhd_reg<-merge(filtered_taxonomy_10_1_gvhd,metadata_gvhd,by="index")
rownames(gvhd_reg)<-gvhd_reg$index
gvhd_reg_name<-gvhd_reg
colnames(gvhd_reg_name)[122]<-"amplicon"
colnames(gvhd_reg_name)[123]<-"sex"
colnames(gvhd_reg_name)[126]<-"scs"
colnames(gvhd_reg_name)[127]<-"DrH"
colnames(gvhd_reg_name)[120]<-"gclass"

```

### single drug analysis

```{r}
number_gvhd<-array(length(unique(gvhd_reg_name$subject)),dim=c(111,48))
rownames(number_gvhd)<-colnames(gvhd_reg_name[,2:112])
colnames(number_gvhd)<-colnames(gvhd_reg_name[,128:175])
nonzero_gvhd<-array(NA,dim=c(111,48))
rownames(nonzero_gvhd)<-colnames(gvhd_reg_name[,2:112])
colnames(nonzero_gvhd)<-colnames(gvhd_reg_name[,128:175])
users_gvhd<-array(NA,dim=c(111,48))
rownames(users_gvhd)<-colnames(gvhd_reg_name[,2:112])
colnames(users_gvhd)<-colnames(gvhd_reg_name[,128:175])
nonusers_gvhd<-array(NA,dim=c(111,48))
rownames(nonusers_gvhd)<-colnames(gvhd_reg_name[,2:112])
colnames(nonusers_gvhd)<-colnames(gvhd_reg_name[,128:175])



estimate_gvhd<-array(0,dim=c(111,48))
rownames(estimate_gvhd)<-colnames(gvhd_reg_name[,2:112])
colnames(estimate_gvhd)<-colnames(gvhd_reg_name[,128:175])
Std_Error_gvhd<-array(0,dim=c(111,48))
rownames(Std_Error_gvhd)<-colnames(gvhd_reg_name[,2:112])
colnames(Std_Error_gvhd)<-colnames(gvhd_reg_name[,128:175])
t_value_gvhd<-array(0,dim=c(111,48))
rownames(t_value_gvhd)<-colnames(gvhd_reg_name[,2:112])
colnames(t_value_gvhd)<-colnames(gvhd_reg_name[,128:175])
p_value_gvhd<-array(0,dim=c(111,48))
rownames(p_value_gvhd)<-colnames(gvhd_reg_name[,2:112])
colnames(p_value_gvhd)<-colnames(gvhd_reg_name[,128:175])


for (x in 1:111) {
  for (y in 1:48){
    genus_name<-colnames(gvhd_reg_name[,2:112])
    drug_name<-colnames(gvhd_reg_name[,128:175])
    
    
    temp_gvhd<-select(gvhd_reg_name,genus_name[x],subject,flag,time,amplicon,gclass,sex,age,diagnosis,scs,DrH,drug_name[y])
    temp_gvhd$species<-temp_gvhd[,1]
    temp_gvhd_nonzero<- temp_gvhd |> filter(species!=0) 
    nonzero_gvhd[x,y] <-length(unique(temp_gvhd_nonzero$subject))
    temp_gvhd$drug<-temp_gvhd[,(ncol(temp_gvhd)-1)]
    temp_gvhd_users<-temp_gvhd[complete.cases(temp_gvhd$drug),] |> filter(drug==1) 
    users_gvhd[x,y] <- length(unique(temp_gvhd_users$subject))
    temp_gvhd_nonusers<-temp_gvhd[complete.cases(temp_gvhd$drug),] |> filter(drug==0) 
    nonusers_gvhd[x,y] <- length(unique(temp_gvhd_nonusers$subject))
    
    
    temp_gvhd<-temp_gvhd[complete.cases(temp_gvhd),]
    drug_gvhd<-temp_gvhd[,ncol(temp_gvhd)]
    if(length(unique(drug_gvhd))==2){
      if(length(unique(temp_gvhd$time))!=1 & length(unique(temp_gvhd$amplicon))!=1 & length(unique(temp_gvhd$sex))!=1 & length(unique(temp_gvhd$scs))!=1){
        temp_gvhd1<-lm(temp_gvhd[,1] ~ time + amplicon + sex + age + gclass+diagnosis + scs + DrH + drug_gvhd,temp_gvhd)
        if (rownames(summary(temp_gvhd1)$coefficients)[nrow(summary(temp_gvhd1)$coefficients)] == "drug_gvhd"){
          summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),1]->estimate_gvhd[x,y]
          summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),2]->Std_Error_gvhd[x,y]
          summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),3]->t_value_gvhd[x,y]
          summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),4]->p_value_gvhd[x,y]
        }else{
          estimate_gvhd[x,y]<-NA
          Std_Error_gvhd[x,y]<-NA
          t_value_gvhd[x,y]<-NA
          p_value_gvhd[x,y]<-NA
        }
      }else if(length(unique(temp_gvhd$time))==1 & length(unique(temp_gvhd$amplicon))!=1 & length(unique(temp_gvhd$sex))!=1 & length(unique(temp_gvhd$scs))!=1){
        temp_gvhd1<-lm(temp_gvhd[,1] ~ amplicon + sex + age + gclass+diagnosis + scs + DrH + drug_gvhd,temp_gvhd)
        if (rownames(summary(temp_gvhd1)$coefficients)[nrow(summary(temp_gvhd1)$coefficients)] == "drug_gvhd"){
          estimate_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),1]
          Std_Error_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),2]
          t_value_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),3]
          p_value_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),4]
        }else{
          estimate_gvhd[x,y]<-NA
          Std_Error_gvhd[x,y]<-NA
          t_value_gvhd[x,y]<-NA
          p_value_gvhd[x,y]<-NA
        }
      }else if(length(unique(temp_gvhd$time))!=1 & length(unique(temp_gvhd$amplicon))==1 & length(unique(temp_gvhd$sex))!=1 & length(unique(temp_gvhd$scs))!=1){
        temp_gvhd1<-lm(temp_gvhd[,1] ~ time + sex + age + gclass+diagnosis + scs + DrH + drug_gvhd,temp_gvhd)
        if (rownames(summary(temp_gvhd1)$coefficients)[nrow(summary(temp_gvhd1)$coefficients)] == "drug_gvhd"){
          estimate_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),1]
          Std_Error_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),2]
          t_value_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),3]
          p_value_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),4]
        }else{
          estimate_gvhd[x,y]<-NA
          Std_Error_gvhd[x,y]<-NA
          t_value_gvhd[x,y]<-NA
          p_value_gvhd[x,y]<-NA
        }
      }else if(length(unique(temp_gvhd$time))!=1 & length(unique(temp_gvhd$amplicon))!=1 & length(unique(temp_gvhd$sex))==1 & length(unique(temp_gvhd$scs))!=1){
        temp_gvhd1<-lm(temp_gvhd[,1] ~ time  + age + gclass+diagnosis + scs + DrH + drug_gvhd,temp_gvhd)
        if (rownames(summary(temp_gvhd1)$coefficients)[nrow(summary(temp_gvhd1)$coefficients)] == "drug_gvhd"){
          estimate_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),1]
          Std_Error_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),2]
          t_value_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),3]
          p_value_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),4]
        }else{
          estimate_gvhd[x,y]<-NA
          Std_Error_gvhd[x,y]<-NA
          t_value_gvhd[x,y]<-NA
          p_value_gvhd[x,y]<-NA
        }
      }else if(length(unique(temp_gvhd$time))!=1 & length(unique(temp_gvhd$amplicon))!=1 & length(unique(temp_gvhd$sex))!=1 & length(unique(temp_gvhd$scs))==1){
        temp_gvhd1<-lm(temp_gvhd[,1] ~ time + sex + age + gclass+diagnosis + DrH + drug_gvhd,temp_gvhd)
        if (rownames(summary(temp_gvhd1)$coefficients)[nrow(summary(temp_gvhd1)$coefficients)] == "drug_gvhd"){
          estimate_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),1]
          Std_Error_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),2]
          t_value_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),3]
          p_value_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),4]
        }else{
          estimate_gvhd[x,y]<-NA
          Std_Error_gvhd[x,y]<-NA
          t_value_gvhd[x,y]<-NA
          p_value_gvhd[x,y]<-NA
        }
      }else if(length(unique(temp_gvhd$time))==1 & length(unique(temp_gvhd$amplicon))==1 & length(unique(temp_gvhd$sex))!=1 & length(unique(temp_gvhd$scs))!=1){
        temp_gvhd1<-lm(temp_gvhd[,1] ~ sex + age + gclass+diagnosis + scs + DrH + drug_gvhd,temp_gvhd)
        if (rownames(summary(temp_gvhd1)$coefficients)[nrow(summary(temp_gvhd1)$coefficients)] == "drug_gvhd"){
          estimate_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),1]
          Std_Error_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),2]
          t_value_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),3]
          p_value_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),4]
        }else{
          estimate_gvhd[x,y]<-NA
          Std_Error_gvhd[x,y]<-NA
          t_value_gvhd[x,y]<-NA
          p_value_gvhd[x,y]<-NA
        }
      }else if(length(unique(temp_gvhd$time))==1 & length(unique(temp_gvhd$amplicon))!=1 & length(unique(temp_gvhd$sex))==1 & length(unique(temp_gvhd$scs))!=1){
        temp_gvhd1<-lm(temp_gvhd[,1] ~ amplicon + age + gclass+diagnosis + scs + DrH + drug_gvhd,temp_gvhd)
        if (rownames(summary(temp_gvhd1)$coefficients)[nrow(summary(temp_gvhd1)$coefficients)] == "drug_gvhd"){
          estimate_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),1]
          Std_Error_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),2]
          t_value_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),3]
          p_value_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),4]
        }else{
          estimate_gvhd[x,y]<-NA
          Std_Error_gvhd[x,y]<-NA
          t_value_gvhd[x,y]<-NA
          p_value_gvhd[x,y]<-NA
        }
      }else if(length(unique(temp_gvhd$time))==1 & length(unique(temp_gvhd$amplicon))!=1 & length(unique(temp_gvhd$sex))!=1 & length(unique(temp_gvhd$scs))==1){
        temp_gvhd1<-lm(temp_gvhd[,1] ~ amplicon + sex + age + gclass+diagnosis + DrH + drug_gvhd,temp_gvhd)
        if (rownames(summary(temp_gvhd1)$coefficients)[nrow(summary(temp_gvhd1)$coefficients)] == "drug_gvhd"){
          estimate_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),1]
          Std_Error_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),2]
          t_value_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),3]
          p_value_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),4]
        }else{
          estimate_gvhd[x,y]<-NA
          Std_Error_gvhd[x,y]<-NA
          t_value_gvhd[x,y]<-NA
          p_value_gvhd[x,y]<-NA
        }
      }else if(length(unique(temp_gvhd$time))!=1 & length(unique(temp_gvhd$amplicon))==1 & length(unique(temp_gvhd$sex))==1 & length(unique(temp_gvhd$scs))!=1){
        temp_gvhd1<-lm(temp_gvhd[,1] ~ time + age + gclass+diagnosis + scs + DrH + drug_gvhd,temp_gvhd)
        if (rownames(summary(temp_gvhd1)$coefficients)[nrow(summary(temp_gvhd1)$coefficients)] == "drug_gvhd"){
          estimate_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),1]
          Std_Error_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),2]
          t_value_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),3]
          p_value_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),4]
        }else{
          estimate_gvhd[x,y]<-NA
          Std_Error_gvhd[x,y]<-NA
          t_value_gvhd[x,y]<-NA
          p_value_gvhd[x,y]<-NA
        }
      }else if(length(unique(temp_gvhd$time))!=1 & length(unique(temp_gvhd$amplicon))==1 & length(unique(temp_gvhd$sex))!=1 & length(unique(temp_gvhd$scs))==1){
        temp_gvhd1<-lm(temp_gvhd[,1] ~ time + sex + age + gclass+diagnosis + DrH + drug_gvhd,temp_gvhd)
        if (rownames(summary(temp_gvhd1)$coefficients)[nrow(summary(temp_gvhd1)$coefficients)] == "drug_gvhd"){
          estimate_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),1]
          Std_Error_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),2]
          t_value_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),3]
          p_value_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),4]
        }else{
          estimate_gvhd[x,y]<-NA
          Std_Error_gvhd[x,y]<-NA
          t_value_gvhd[x,y]<-NA
          p_value_gvhd[x,y]<-NA
        }
      }else if(length(unique(temp_gvhd$time))!=1 & length(unique(temp_gvhd$amplicon))!=1 & length(unique(temp_gvhd$sex))==1 & length(unique(temp_gvhd$scs))==1){
        temp_gvhd1<-lm(temp_gvhd[,1] ~ amplicon + time + age + gclass+diagnosis + DrH + drug_gvhd,temp_gvhd)
        if (rownames(summary(temp_gvhd1)$coefficients)[nrow(summary(temp_gvhd1)$coefficients)] == "drug_gvhd"){
          estimate_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),1]
          Std_Error_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),2]
          t_value_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),3]
          p_value_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),4]
        }else{
          estimate_gvhd[x,y]<-NA
          Std_Error_gvhd[x,y]<-NA
          t_value_gvhd[x,y]<-NA
          p_value_gvhd[x,y]<-NA
        }
      }else if(length(unique(temp_gvhd$time))==1 & length(unique(temp_gvhd$amplicon))==1 & length(unique(temp_gvhd$sex))==1 & length(unique(temp_gvhd$scs))!=1){
        temp_gvhd1<-lm(temp_gvhd[,1] ~ age + gclass+diagnosis + scs + DrH + drug_gvhd,temp_gvhd)
        if (rownames(summary(temp_gvhd1)$coefficients)[nrow(summary(temp_gvhd1)$coefficients)] == "drug_gvhd"){
          estimate_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),1]
          Std_Error_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),2]
          t_value_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),3]
          p_value_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),4]
        }else{
          estimate_gvhd[x,y]<-NA
          Std_Error_gvhd[x,y]<-NA
          t_value_gvhd[x,y]<-NA
          p_value_gvhd[x,y]<-NA
        }
      }else if(length(unique(temp_gvhd$time))==1 & length(unique(temp_gvhd$amplicon))==1 & length(unique(temp_gvhd$sex))!=1 & length(unique(temp_gvhd$scs))==1){
        temp_gvhd1<-lm(temp_gvhd[,1] ~ sex + age + gclass+diagnosis + DrH + drug_gvhd,temp_gvhd)
        if (rownames(summary(temp_gvhd1)$coefficients)[nrow(summary(temp_gvhd1)$coefficients)] == "drug_gvhd"){
          estimate_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),1]
          Std_Error_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),2]
          t_value_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),3]
          p_value_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),4]
        }else{
          estimate_gvhd[x,y]<-NA
          Std_Error_gvhd[x,y]<-NA
          t_value_gvhd[x,y]<-NA
          p_value_gvhd[x,y]<-NA
        }
      }else if(length(unique(temp_gvhd$time))==1 & length(unique(temp_gvhd$amplicon))!=1 & length(unique(temp_gvhd$sex))==1 & length(unique(temp_gvhd$scs))==1){
        temp_gvhd1<-lm(temp_gvhd[,1] ~ amplicon + age + gclass+diagnosis + DrH + drug_gvhd,temp_gvhd)
        if (rownames(summary(temp_gvhd1)$coefficients)[nrow(summary(temp_gvhd1)$coefficients)] == "drug_gvhd"){
          estimate_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),1]
          Std_Error_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),2]
          t_value_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),3]
          p_value_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),4]
        }else{
          estimate_gvhd[x,y]<-NA
          Std_Error_gvhd[x,y]<-NA
          t_value_gvhd[x,y]<-NA
          p_value_gvhd[x,y]<-NA
        }
      }else if(length(unique(temp_gvhd$time))!=1 & length(unique(temp_gvhd$amplicon))==1 & length(unique(temp_gvhd$sex))==1 & length(unique(temp_gvhd$scs))==1){
        temp_gvhd1<-lm(temp_gvhd[,1] ~ time + age + gclass+diagnosis + DrH + drug_gvhd,temp_gvhd)
        if (rownames(summary(temp_gvhd1)$coefficients)[nrow(summary(temp_gvhd1)$coefficients)] == "drug_gvhd"){
          estimate_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),1]
          Std_Error_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),2]
          t_value_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),3]
          p_value_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),4]
        }else{
          estimate_gvhd[x,y]<-NA
          Std_Error_gvhd[x,y]<-NA
          t_value_gvhd[x,y]<-NA
          p_value_gvhd[x,y]<-NA
        }
      }else{
        temp_gvhd1<-lm(temp_gvhd[,1] ~ age + gclass+diagnosis + DrH + drug_gvhd,temp_gvhd)
        if (rownames(summary(temp_gvhd1)$coefficients)[nrow(summary(temp_gvhd1)$coefficients)] == "drug_gvhd"){
          estimate_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),1]
          Std_Error_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),2]
          t_value_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),3]
          p_value_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),4]
        }else{
          estimate_gvhd[x,y]<-NA
          Std_Error_gvhd[x,y]<-NA
          t_value_gvhd[x,y]<-NA
          p_value_gvhd[x,y]<-NA
        }
      }
    }else{
      estimate_gvhd[x,y]<-NA
      Std_Error_gvhd[x,y]<-NA
      t_value_gvhd[x,y]<-NA
      p_value_gvhd[x,y]<-NA
    }
  }
}

number_gvhd<-as.data.frame(number_gvhd)
nonzero_gvhd<-as.data.frame(nonzero_gvhd)
users_gvhd<-as.data.frame(users_gvhd)
nonusers_gvhd<-as.data.frame(nonusers_gvhd)


estimate_gvhd<-as.data.frame(estimate_gvhd)
p_value_gvhd<-as.data.frame(p_value_gvhd)
Std_Error_gvhd<-as.data.frame(Std_Error_gvhd)
t_value_gvhd<-as.data.frame(t_value_gvhd)

```

### multiple drug analysis

```{r}

number_mul_gvhd<-array(length(unique(gvhd_reg_name$subject)),dim=c(111,47))
rownames(number_mul_gvhd)<-colnames(gvhd_reg_name[,2:112])
colnames(number_mul_gvhd)<-colnames(gvhd_reg_name[,c(128:149,151:175)])
nonzero_mul_gvhd<-array(NA,dim=c(111,47))
rownames(nonzero_mul_gvhd)<-colnames(gvhd_reg_name[,2:112])
colnames(nonzero_mul_gvhd)<-colnames(gvhd_reg_name[,c(128:149,151:175)])
users_mul_gvhd<-array(NA,dim=c(111,47))
rownames(users_mul_gvhd)<-colnames(gvhd_reg_name[,2:112])
colnames(users_mul_gvhd)<-colnames(gvhd_reg_name[,c(128:149,151:175)])
nonusers_mul_gvhd<-array(NA,dim=c(111,47))
rownames(nonusers_mul_gvhd)<-colnames(gvhd_reg_name[,2:112])
colnames(nonusers_mul_gvhd)<-colnames(gvhd_reg_name[,c(128:149,151:175)])


estimate_mul_gvhd<-array(0,dim=c(111,47))
rownames(estimate_mul_gvhd)<-colnames(gvhd_reg_name[,2:112])
colnames(estimate_mul_gvhd)<-colnames(gvhd_reg_name[,c(128:149,151:175)])
Std_Error_mul_gvhd<-array(0,dim=c(111,47))
rownames(Std_Error_mul_gvhd)<-colnames(gvhd_reg_name[,2:112])
colnames(Std_Error_mul_gvhd)<-colnames(gvhd_reg_name[,c(128:149,151:175)])
t_value_mul_gvhd<-array(0,dim=c(111,47))
rownames(t_value_mul_gvhd)<-colnames(gvhd_reg_name[,2:112])
colnames(t_value_mul_gvhd)<-colnames(gvhd_reg_name[,c(128:149,151:175)])
p_value_mul_gvhd<-array(0,dim=c(111,47))
rownames(p_value_mul_gvhd)<-colnames(gvhd_reg_name[,2:112])
colnames(p_value_mul_gvhd)<-colnames(gvhd_reg_name[,c(128:149,151:175)])


for (x in 1:111) {
  for (y in 1:47){
    genus_name<-colnames(gvhd_reg_name[,2:112])
    drug_name<-colnames(gvhd_reg_name[,c(128:149,151:175)])
    
    temp_mul_gvhd<-select(gvhd_reg_name,genus_name[x],subject,flag,time,amplicon,sex,age,diagnosis,scs,DrH,Antibiotic,drug_name[y])
    temp_mul_gvhd$species<-temp_mul_gvhd[,1]
    temp_mul_gvhd_nonzero<- temp_mul_gvhd |> filter(species!=0) 
    nonzero_mul_gvhd[x,y] <-length(unique(temp_mul_gvhd_nonzero$subject))
    temp_mul_gvhd$drug<-temp_mul_gvhd[,(ncol(temp_mul_gvhd)-1)]
    temp_mul_gvhd_users<-temp_mul_gvhd[complete.cases(temp_mul_gvhd$drug),] |> filter(drug==1) 
    users_mul_gvhd[x,y] <- length(unique(temp_mul_gvhd_users$subject))
    temp_mul_gvhd_nonusers<-temp_mul_gvhd[complete.cases(temp_mul_gvhd$drug),] |> filter(drug==0) 
    nonusers_mul_gvhd[x,y] <- length(unique(temp_mul_gvhd_nonusers$subject))
    
    
    temp_mul_gvhd<-temp_mul_gvhd[complete.cases(temp_mul_gvhd),]
    drug<-temp_mul_gvhd[,ncol(temp_mul_gvhd)]
    Antibiotic<-temp_mul_gvhd[,"Antibiotic"]
    if(length(unique(drug))==2 & length(unique(Antibiotic))==2){
      if(length(unique(temp_mul_gvhd$time))!=1 & length(unique(temp_mul_gvhd$amplicon))!=1 & length(unique(temp_mul_gvhd$sex))!=1 & length(unique(temp_mul_gvhd$scs))!=1){
        temp_mul_gvhd1<-lm(temp_mul_gvhd[,1] ~ time + amplicon + sex + age + diagnosis + scs + DrH + Antibiotic +drug,temp_mul_gvhd)
        if (rownames(summary(temp_mul_gvhd1)$coefficients)[nrow(summary(temp_mul_gvhd1)$coefficients)] == "drug"){
          summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),1]->estimate_mul_gvhd[x,y]
          summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),2]->Std_Error_mul_gvhd[x,y]
          summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),3]->t_value_mul_gvhd[x,y]
          summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),4]->p_value_mul_gvhd[x,y]
        }else{
          estimate_mul_gvhd[x,y]<-NA
          Std_Error_mul_gvhd[x,y]<-NA
          t_value_mul_gvhd[x,y]<-NA
          p_value_mul_gvhd[x,y]<-NA
        }
      }else if(length(unique(temp_mul_gvhd$time))==1 & length(unique(temp_mul_gvhd$amplicon))!=1 & length(unique(temp_mul_gvhd$sex))!=1 & length(unique(temp_mul_gvhd$scs))!=1){
        temp_mul_gvhd1<-lm(temp_mul_gvhd[,1] ~ amplicon + sex + age + diagnosis + scs + DrH + Antibiotic+ drug,temp_mul_gvhd)
        if (rownames(summary(temp_mul_gvhd1)$coefficients)[nrow(summary(temp_mul_gvhd1)$coefficients)] == "drug"){
          estimate_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),1]
          Std_Error_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),2]
          t_value_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),3]
          p_value_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),4]
        }else{
          estimate_mul_gvhd[x,y]<-NA
          Std_Error_mul_gvhd[x,y]<-NA
          t_value_mul_gvhd[x,y]<-NA
          p_value_mul_gvhd[x,y]<-NA
        }
      }else if(length(unique(temp_mul_gvhd$time))!=1 & length(unique(temp_mul_gvhd$amplicon))==1 & length(unique(temp_mul_gvhd$sex))!=1 & length(unique(temp_mul_gvhd$scs))!=1){
        temp_mul_gvhd1<-lm(temp_mul_gvhd[,1] ~ time + sex + age + diagnosis + scs + DrH + Antibiotic+ drug,temp_mul_gvhd)
        if (rownames(summary(temp_mul_gvhd1)$coefficients)[nrow(summary(temp_mul_gvhd1)$coefficients)] == "drug"){
          estimate_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),1]
          Std_Error_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),2]
          t_value_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),3]
          p_value_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),4]
        }else{
          estimate_mul_gvhd[x,y]<-NA
          Std_Error_mul_gvhd[x,y]<-NA
          t_value_mul_gvhd[x,y]<-NA
          p_value_mul_gvhd[x,y]<-NA
        }
      }else if(length(unique(temp_mul_gvhd$time))!=1 & length(unique(temp_mul_gvhd$amplicon))!=1 & length(unique(temp_mul_gvhd$sex))==1 & length(unique(temp_mul_gvhd$scs))!=1){
        temp_mul_gvhd1<-lm(temp_mul_gvhd[,1] ~ time  + age + diagnosis + scs + DrH + Antibiotic+ drug,temp_mul_gvhd)
        if (rownames(summary(temp_mul_gvhd1)$coefficients)[nrow(summary(temp_mul_gvhd1)$coefficients)] == "drug"){
          estimate_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),1]
          Std_Error_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),2]
          t_value_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),3]
          p_value_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),4]
        }else{
          estimate_mul_gvhd[x,y]<-NA
          Std_Error_mul_gvhd[x,y]<-NA
          t_value_mul_gvhd[x,y]<-NA
          p_value_mul_gvhd[x,y]<-NA
        }
      }else if(length(unique(temp_mul_gvhd$time))!=1 & length(unique(temp_mul_gvhd$amplicon))!=1 & length(unique(temp_mul_gvhd$sex))!=1 & length(unique(temp_mul_gvhd$scs))==1){
        temp_mul_gvhd1<-lm(temp_mul_gvhd[,1] ~ time + sex + age + diagnosis + DrH+ Antibiotic + drug,temp_mul_gvhd)
        if (rownames(summary(temp_mul_gvhd1)$coefficients)[nrow(summary(temp_mul_gvhd1)$coefficients)] == "drug"){
          estimate_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),1]
          Std_Error_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),2]
          t_value_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),3]
          p_value_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),4]
        }else{
          estimate_mul_gvhd[x,y]<-NA
          Std_Error_mul_gvhd[x,y]<-NA
          t_value_mul_gvhd[x,y]<-NA
          p_value_mul_gvhd[x,y]<-NA
        }
      }else if(length(unique(temp_mul_gvhd$time))==1 & length(unique(temp_mul_gvhd$amplicon))==1 & length(unique(temp_mul_gvhd$sex))!=1 & length(unique(temp_mul_gvhd$scs))!=1){
        temp_mul_gvhd1<-lm(temp_mul_gvhd[,1] ~ sex + age + diagnosis + scs + DrH + Antibiotic+ drug,temp_mul_gvhd)
        if (rownames(summary(temp_mul_gvhd1)$coefficients)[nrow(summary(temp_mul_gvhd1)$coefficients)] == "drug"){
          estimate_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),1]
          Std_Error_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),2]
          t_value_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),3]
          p_value_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),4]
        }else{
          estimate_mul_gvhd[x,y]<-NA
          Std_Error_mul_gvhd<-NA
          t_value_mul_gvhd[x,y]<-NA
          p_value_mul_gvhd[x,y]<-NA
        }
      }else if(length(unique(temp_mul_gvhd$time))==1 & length(unique(temp_mul_gvhd$amplicon))!=1 & length(unique(temp_mul_gvhd$sex))==1 & length(unique(temp_mul_gvhd$scs))!=1){
        temp_mul_gvhd1<-lm(temp_mul_gvhd[,1] ~ amplicon + age + diagnosis + scs + DrH + Antibiotic+ drug,temp_mul_gvhd)
        if (rownames(summary(temp_mul_gvhd1)$coefficients)[nrow(summary(temp_mul_gvhd1)$coefficients)] == "drug"){
          estimate_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),1]
          Std_Error_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),2]
          t_value_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),3]
          p_value_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),4]
        }else{
          estimate_mul_gvhd[x,y]<-NA
          Std_Error_mul_gvhd[x,y]<-NA
          t_value_mul_gvhd[x,y]<-NA
          p_value_mul_gvhd[x,y]<-NA
        }
      }else if(length(unique(temp_mul_gvhd$time))==1 & length(unique(temp_mul_gvhd$amplicon))!=1 & length(unique(temp_mul_gvhd$sex))!=1 & length(unique(temp_mul_gvhd$scs))==1){
        temp_mul_gvhd1<-lm(temp_mul_gvhd[,1] ~ amplicon + sex + age + diagnosis + DrH + Antibiotic+ drug,temp_mul_gvhd)
        if (rownames(summary(temp_mul_gvhd1)$coefficients)[nrow(summary(temp_mul_gvhd1)$coefficients)] == "drug"){
          estimate_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),1]
          Std_Error_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),2]
          t_value_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),3]
          p_value_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),4]
        }else{
          estimate_mul_gvhd[x,y]<-NA
          Std_Error_mul_gvhd[x,y]<-NA
          t_value_mul_gvhd[x,y]<-NA
          p_value_mul_gvhd[x,y]<-NA
        }
      }else if(length(unique(temp_mul_gvhd$time))!=1 & length(unique(temp_mul_gvhd$amplicon))==1 & length(unique(temp_mul_gvhd$sex))==1 & length(unique(temp_mul_gvhd$scs))!=1){
        temp_mul_gvhd1<-lm(temp_mul_gvhd[,1] ~ time + age + diagnosis + scs + DrH + Antibiotic+ drug,temp_mul_gvhd)
        if (rownames(summary(temp_mul_gvhd1)$coefficients)[nrow(summary(temp_mul_gvhd1)$coefficients)] == "drug"){
          estimate_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),1]
          Std_Error_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),2]
          t_value_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),3]
          p_value_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),4]
        }else{
          estimate_mul_gvhd[x,y]<-NA
          Std_Error_mul_gvhd[x,y]<-NA
          t_value_mul_gvhd[x,y]<-NA
          p_value_mul_gvhd[x,y]<-NA
        }
      }else if(length(unique(temp_mul_gvhd$time))!=1 & length(unique(temp_mul_gvhd$amplicon))==1 & length(unique(temp_mul_gvhd$sex))!=1 & length(unique(temp_mul_gvhd$scs))==1){
        temp_mul_gvhd1<-lm(temp_mul_gvhd[,1] ~ time + sex + age + diagnosis + DrH + Antibiotic+ drug,temp_mul_gvhd)
        if (rownames(summary(temp_mul_gvhd1)$coefficients)[nrow(summary(temp_mul_gvhd1)$coefficients)] == "drug"){
          estimate_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),1]
          Std_Error_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),2]
          t_value_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),3]
          p_value_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),4]
        }else{
          estimate_mul_gvhd[x,y]<-NA
          Std_Error_mul_gvhd[x,y]<-NA
          t_value_mul_gvhd[x,y]<-NA
          p_value_mul_gvhd[x,y]<-NA
        }
      }else if(length(unique(temp_mul_gvhd$time))!=1 & length(unique(temp_mul_gvhd$amplicon))!=1 & length(unique(temp_mul_gvhd$sex))==1 & length(unique(temp_mul_gvhd$scs))==1){
        temp_mul_gvhd1<-lm(temp_mul_gvhd[,1] ~ amplicon + time + age + diagnosis + DrH + Antibiotic+ drug,temp_mul_gvhd)
        if (rownames(summary(temp_mul_gvhd1)$coefficients)[nrow(summary(temp_mul_gvhd1)$coefficients)] == "drug"){
          estimate_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),1]
          Std_Error_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),2]
          t_value_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),3]
          p_value_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),4]
        }else{
          estimate_mul_gvhd[x,y]<-NA
          Std_Error_mul_gvhd[x,y]<-NA
          t_value_mul_gvhd[x,y]<-NA
          p_value_mul_gvhd[x,y]<-NA
        }
      }else if(length(unique(temp_mul_gvhd$time))==1 & length(unique(temp_mul_gvhd$amplicon))==1 & length(unique(temp_mul_gvhd$sex))==1 & length(unique(temp_mul_gvhd$scs))!=1){
        temp_mul_gvhd1<-lm(temp_mul_gvhd[,1] ~ age + diagnosis + scs + DrH + Antibiotic+ drug,temp_mul_gvhd)
        if (rownames(summary(temp_mul_gvhd1)$coefficients)[nrow(summary(temp_mul_gvhd1)$coefficients)] == "drug"){
          estimate_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),1]
          Std_Error_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),2]
          t_value_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),3]
          p_value_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),4]
        }else{
          estimate_mul_gvhd[x,y]<-NA
          Std_Error_mul_gvhd[x,y]<-NA
          t_value_mul_gvhd[x,y]<-NA
          p_value_mul_gvhd[x,y]<-NA
        }
      }else if(length(unique(temp_mul_gvhd$time))==1 & length(unique(temp_mul_gvhd$amplicon))==1 & length(unique(temp_mul_gvhd$sex))!=1 & length(unique(temp_mul_gvhd$scs))==1){
        temp_mul_gvhd1<-lm(temp_mul_gvhd[,1] ~ sex + age + diagnosis + DrH + Antibiotic+ drug,temp_mul_gvhd)
        if (rownames(summary(temp_mul_gvhd1)$coefficients)[nrow(summary(temp_mul_gvhd1)$coefficients)] == "drug"){
          estimate_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),1]
          Std_Error_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),2]
          t_value_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),3]
          p_value_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),4]
        }else{
          estimate_mul_gvhd[x,y]<-NA
          Std_Error_mul_gvhd[x,y]<-NA
          t_value_mul_gvhd[x,y]<-NA
          p_value_mul_gvhd[x,y]<-NA
        }
      }else if(length(unique(temp_mul_gvhd$time))==1 & length(unique(temp_mul_gvhd$amplicon))!=1 & length(unique(temp_mul_gvhd$sex))==1 & length(unique(temp_mul_gvhd$scs))==1){
        temp_mul_gvhd1<-lm(temp_mul_gvhd[,1] ~ amplicon + age + diagnosis + DrH + Antibiotic+ drug,temp_mul_gvhd)
        if (rownames(summary(temp_mul_gvhd1)$coefficients)[nrow(summary(temp_mul_gvhd1)$coefficients)] == "drug"){
          estimate_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),1]
          Std_Error_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),2]
          t_value_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),3]
          p_value_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),4]
        }else{
          estimate_mul_gvhd[x,y]<-NA
          Std_Error_mul_gvhd[x,y]<-NA
          t_value_mul_gvhd[x,y]<-NA
          p_value_mul_gvhd[x,y]<-NA
        }
      }else if(length(unique(temp_mul_gvhd$time))!=1 & length(unique(temp_mul_gvhd$amplicon))==1 & length(unique(temp_mul_gvhd$sex))==1 & length(unique(temp_mul_gvhd$scs))==1){
        temp_mul_gvhd1<-lm(temp_mul_gvhd[,1] ~ time + age + diagnosis + DrH + Antibiotic+ drug,temp_mul_gvhd)
        if (rownames(summary(temp_mul_gvhd1)$coefficients)[nrow(summary(temp_mul_gvhd1)$coefficients)] == "drug"){
          estimate_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),1]
          Std_Error_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),2]
          t_value_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),3]
          p_value_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),4]
        }else{
          estimate_mul_gvhd[x,y]<-NA
          Std_Error_mul_gvhd[x,y]<-NA
          t_value_mul_gvhd[x,y]<-NA
          p_value_mul_gvhd[x,y]<-NA
        }
      }else{
        temp_mul_gvhd1<-lm(temp_mul_gvhd[,1] ~ age + diagnosis + DrH + Antibiotic+ drug,temp_mul_gvhd)
        if (rownames(summary(temp_mul_gvhd1)$coefficients)[nrow(summary(temp_mul_gvhd1)$coefficients)] == "drug"){
          estimate_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),1]
          Std_Error_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),2]
          t_value_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),3]
          p_value_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),4]
        }else{
          estimate_mul_gvhd[x,y]<-NA
          Std_Error_mul_gvhd[x,y]<-NA
          t_value_mul_gvhd[x,y]<-NA
          p_value_mul_gvhd[x,y]<-NA
        }
      }
    }else{
      estimate_mul_gvhd[x,y]<-NA
      Std_Error_mul_gvhd[x,y]<-NA
      t_value_mul_gvhd[x,y]<-NA
      p_value_mul_gvhd[x,y]<-NA
    }
  }
}

number_mul_gvhd<-as.data.frame(number_mul_gvhd)
nonzero_mul_gvhd<-as.data.frame(nonzero_mul_gvhd)
users_mul_gvhd<-as.data.frame(users_mul_gvhd)
nonusers_mul_gvhd<-as.data.frame(nonusers_mul_gvhd)

estimate_mul_gvhd<-as.data.frame(estimate_mul_gvhd)
p_value_mul_gvhd<-as.data.frame(p_value_mul_gvhd)
Std_Error_mul_gvhd<-as.data.frame(Std_Error_mul_gvhd)
t_value_mul_gvhd<-as.data.frame(t_value_mul_gvhd)

```

## 4.3 Chemo Cohort

```{r}

chemo_0<-all_reg_chemo
chemo_0[is.na(chemo_0)]<-0

sum_abundance_chemo<-apply(chemo_0[,2:850],1,sum)

chemo_rel<-array(0,dim=c(1097,849))
for (j in 1:849) {
  for (i in 1:1097){
    chemo_rel[i,j]<-chemo_0[i,j+1]/sum_abundance_chemo[i]*100
  }
}
chemo_rel <- as.data.frame(chemo_rel)
chemo_rel[is.na(chemo_rel)]<-0

sample_chemo<-chemo_0$index
species_chemo<-colnames(chemo_0)
colnames(chemo_rel)<-species_chemo[2:850]
rownames(chemo_rel)<-sample_chemo

chemo_rel_fil_10 <- filtering_taxonomy(chemo_rel,0.00000001,10)

summary_taxnomy_filtering_10_chemo<- read.table("summary_taxonomy_filtering.txt",sep="\t",header=T)
filtered_taxonomy_10_chemo<-read.table("filtered_taxonomy.txt",sep="\t",header=T)

filtered_taxonomy_10_1_chemo<-cbind(filtered_taxonomy_10_chemo,index=rownames(filtered_taxonomy_10_chemo))

metadata_chemo<-as.data.frame(cbind(index=sample_chemo,chemo[,851:915]))
chemo_reg<-merge(filtered_taxonomy_10_1_chemo,metadata_chemo,by="index")
rownames(chemo_reg)<-chemo_reg$index
chemo_reg_name<-chemo_reg
colnames(chemo_reg_name)[109]<-"amplicon"
colnames(chemo_reg_name)[110]<-"sex"
colnames(chemo_reg_name)[113]<-"scs"
colnames(chemo_reg_name)[114]<-"DrH"


```

### single drug analysis

```{r}
number_chemo<-array(length(unique(chemo_reg_name$subject)),dim=c(98,48))
rownames(number_chemo)<-colnames(chemo_reg_name[,2:99])
colnames(number_chemo)<-colnames(chemo_reg_name[,115:162])
nonzero_chemo<-array(NA,dim=c(98,48))
rownames(nonzero_chemo)<-colnames(chemo_reg_name[,2:99])
colnames(nonzero_chemo)<-colnames(chemo_reg_name[,115:162])
users_chemo<-array(NA,dim=c(98,48))
rownames(users_chemo)<-colnames(chemo_reg_name[,2:99])
colnames(users_chemo)<-colnames(chemo_reg_name[,115:162])
nonusers_chemo<-array(NA,dim=c(98,48))
rownames(nonusers_chemo)<-colnames(chemo_reg_name[,2:99])
colnames(nonusers_chemo)<-colnames(chemo_reg_name[,115:162])



estimate_chemo<-array(0,dim=c(98,48))
rownames(estimate_chemo)<-colnames(chemo_reg_name[,2:99])
colnames(estimate_chemo)<-colnames(chemo_reg_name[,115:162])
Std_Error_chemo<-array(0,dim=c(98,48))
rownames(Std_Error_chemo)<-colnames(chemo_reg_name[,2:99])
colnames(Std_Error_chemo)<-colnames(chemo_reg_name[,115:162])
t_value_chemo<-array(0,dim=c(98,48))
rownames(t_value_chemo)<-colnames(chemo_reg_name[,2:99])
colnames(t_value_chemo)<-colnames(chemo_reg_name[,115:162])
p_value_chemo<-array(0,dim=c(98,48))
rownames(p_value_chemo)<-colnames(chemo_reg_name[,2:99])
colnames(p_value_chemo)<-colnames(chemo_reg_name[,115:162])


for (x in 1:98) {
  for (y in 1:48){
    genus_name<-colnames(chemo_reg_name[,2:99])
    drug_name<-colnames(chemo_reg_name[,115:162])
    
    temp_chemo<-select(chemo_reg_name,genus_name[x],subject,flag,time,amplicon,sex,age,diagnosis,drug_name[y])
    
    temp_chemo$species<-temp_chemo[,1]
    temp_chemo_nonzero<- temp_chemo |> filter(species!=0) 
    nonzero_chemo[x,y] <-length(unique(temp_chemo_nonzero$subject))
    temp_chemo$drug<-temp_chemo[,(ncol(temp_chemo)-1)]
    temp_chemo_users<-temp_chemo[complete.cases(temp_chemo$drug),] |> filter(drug==1) 
    users_chemo[x,y] <- length(unique(temp_chemo_users$subject))
    temp_chemo_nonusers<-temp_chemo[complete.cases(temp_chemo$drug),] |> filter(drug==0) 
    nonusers_chemo[x,y] <- length(unique(temp_nonusers$subject))
    
    
    temp_chemo<-temp_chemo[complete.cases(temp_chemo),]
    drug_chemo<-temp_chemo[,ncol(temp_chemo)]
    if(length(unique(drug_chemo))==2){
      if(length(unique(temp_chemo$time))!=1 & length(unique(temp_chemo$amplicon))!=1 & length(unique(temp_chemo$sex))!=1 & length(unique(temp_chemo$diagnosis))!=1){
        temp_chemo1<-lm(temp_chemo[,1] ~ time + amplicon + sex + age + diagnosis + drug_chemo,temp_chemo)
        if (rownames(summary(temp_chemo1)$coefficients)[nrow(summary(temp_chemo1)$coefficients)] == "drug_chemo"){
          summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),1]->estimate_chemo[x,y]
          summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),2]->Std_Error_chemo[x,y]
          summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),3]->t_value_chemo[x,y]
          summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),4]->p_value_chemo[x,y]
        }else{
          estimate_chemo[x,y]<-NA
          Std_Error_chemo[x,y]<-NA
          t_value_chemo[x,y]<-NA
          p_value_chemo[x,y]<-NA
        }
      }else if(length(unique(temp_chemo$time))==1 & length(unique(temp_chemo$amplicon))!=1 & length(unique(temp_chemo$sex))!=1 & length(unique(temp_chemo$diagnosis))!=1){
        temp_chemo1<-lm(temp_chemo[,1] ~ amplicon + sex + age + diagnosis + drug_chemo,temp_chemo)
        if (rownames(summary(temp_chemo1)$coefficients)[nrow(summary(temp_chemo1)$coefficients)] == "drug_chemo"){
          estimate_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),1]
          Std_Error_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),2]
          t_value_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),3]
          p_value_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),4]
        }else{
          estimate_chemo[x,y]<-NA
          Std_Error_chemo[x,y]<-NA
          t_value_chemo[x,y]<-NA
          p_value_chemo[x,y]<-NA
        }
      }else if(length(unique(temp_chemo$time))!=1 & length(unique(temp_chemo$amplicon))==1 & length(unique(temp_chemo$sex))!=1 & length(unique(temp_chemo$diagnosis))!=1){
        temp_chemo1<-lm(temp_chemo[,1] ~ time + sex + age + diagnosis + drug_chemo,temp_chemo)
        if (rownames(summary(temp_chemo1)$coefficients)[nrow(summary(temp_chemo1)$coefficients)] == "drug_chemo"){
          estimate_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),1]
          Std_Error_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),2]
          t_value_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),3]
          p_value_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),4]
        }else{
          estimate_chemo[x,y]<-NA
          Std_Error_chemo[x,y]<-NA
          t_value_chemo[x,y]<-NA
          p_value_chemo[x,y]<-NA
        }
      }else if(length(unique(temp_chemo$time))!=1 & length(unique(temp_chemo$amplicon))!=1 & length(unique(temp_chemo$sex))==1 & length(unique(temp_chemo$diagnosis))!=1){
        temp_chemo1<-lm(temp_chemo[,1] ~ time  + age + diagnosis + drug_chemo,temp_chemo)
        if (rownames(summary(temp_chemo1)$coefficients)[nrow(summary(temp_chemo1)$coefficients)] == "drug_chemo"){
          estimate_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),1]
          Std_Error_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),2]
          t_value_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),3]
          p_value_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),4]
        }else{
          estimate_chemo[x,y]<-NA
          Std_Error_chemo[x,y]<-NA
          t_value_chemo[x,y]<-NA
          p_value_chemo[x,y]<-NA
        }
      }else if(length(unique(temp_chemo$time))!=1 & length(unique(temp_chemo$amplicon))!=1 & length(unique(temp_chemo$sex))!=1 & length(unique(temp_chemo$diagnosis))==1){
        temp_chemo1<-lm(temp_chemo[,1] ~ time + sex + age + drug_chemo,temp_chemo)
        if (rownames(summary(temp_chemo1)$coefficients)[nrow(summary(temp_chemo1)$coefficients)] == "drug_chemo"){
          estimate_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),1]
          Std_Error_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),2]
          t_value_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),3]
          p_value_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),4]
        }else{
          estimate_chemo[x,y]<-NA
          Std_Error_chemo[x,y]<-NA
          t_value_chemo[x,y]<-NA
          p_value_chemo[x,y]<-NA
        }
      }else if(length(unique(temp_chemo$time))==1 & length(unique(temp_chemo$amplicon))==1 & length(unique(temp_chemo$sex))!=1 & length(unique(temp_chemo$diagnosis))!=1){
        temp_chemo1<-lm(temp_chemo[,1] ~ sex + age + diagnosis + drug_chemo,temp_chemo)
        if (rownames(summary(temp_chemo1)$coefficients)[nrow(summary(temp_chemo1)$coefficients)] == "drug_chemo"){
          estimate_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),1]
          Std_Error_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),2]
          t_value_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),3]
          p_value_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),4]
        }else{
          estimate_chemo[x,y]<-NA
          Std_Error_chemo[x,y]<-NA
          t_value_chemo[x,y]<-NA
          p_value_chemo[x,y]<-NA
        }
      }else if(length(unique(temp_chemo$time))==1 & length(unique(temp_chemo$amplicon))!=1 & length(unique(temp_chemo$sex))==1 & length(unique(temp_chemo$diagnosis))!=1){
        temp_chemo1<-lm(temp_chemo[,1] ~ amplicon + age + diagnosis + drug_chemo,temp_chemo)
        if (rownames(summary(temp_chemo1)$coefficients)[nrow(summary(temp_chemo1)$coefficients)] == "drug_chemo"){
          estimate_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),1]
          Std_Error_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),2]
          t_value_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),3]
          p_value_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),4]
        }else{
          estimate_chemo[x,y]<-NA
          Std_Error_chemo[x,y]<-NA
          t_value_chemo[x,y]<-NA
          p_value_chemo[x,y]<-NA
        }
      }else if(length(unique(temp_chemo$time))==1 & length(unique(temp_chemo$amplicon))!=1 & length(unique(temp_chemo$sex))!=1 & length(unique(temp_chemo$diagnosis))==1){
        temp_chemo1<-lm(temp_chemo[,1] ~ amplicon + sex + age + drug_chemo,temp_chemo)
        if (rownames(summary(temp_chemo1)$coefficients)[nrow(summary(temp_chemo1)$coefficients)] == "drug_chemo"){
          estimate_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),1]
          Std_Error_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),2]
          t_value_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),3]
          p_value_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),4]
        }else{
          estimate_chemo[x,y]<-NA
          Std_Error_chemo[x,y]<-NA
          t_value_chemo[x,y]<-NA
          p_value_chemo[x,y]<-NA
        }
      }else if(length(unique(temp_chemo$time))!=1 & length(unique(temp_chemo$amplicon))==1 & length(unique(temp_chemo$sex))==1 & length(unique(temp_chemo$diagnosis))!=1){
        temp_chemo1<-lm(temp_chemo[,1] ~ time + age + diagnosis + drug_chemo,temp_chemo)
        if (rownames(summary(temp_chemo1)$coefficients)[nrow(summary(temp_chemo1)$coefficients)] == "drug_chemo"){
          estimate_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),1]
          Std_Error_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),2]
          t_value_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),3]
          p_value_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),4]
        }else{
          estimate_chemo[x,y]<-NA
          Std_Error_chemo[x,y]<-NA
          t_value_chemo[x,y]<-NA
          p_value_chemo[x,y]<-NA
        }
      }else if(length(unique(temp_chemo$time))!=1 & length(unique(temp_chemo$amplicon))==1 & length(unique(temp_chemo$sex))!=1 & length(unique(temp_chemo$diagnosis))==1){
        temp_chemo1<-lm(temp_chemo[,1] ~ time + sex + age + drug_chemo,temp_chemo)
        if (rownames(summary(temp_chemo1)$coefficients)[nrow(summary(temp_chemo1)$coefficients)] == "drug_chemo"){
          estimate_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),1]
          Std_Error_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),2]
          t_value_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),3]
          p_value_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),4]
        }else{
          estimate_chemo[x,y]<-NA
          Std_Error_chemo[x,y]<-NA
          t_value_chemo[x,y]<-NA
          p_value_chemo[x,y]<-NA
        }
      }else if(length(unique(temp_chemo$time))!=1 & length(unique(temp_chemo$amplicon))!=1 & length(unique(temp_chemo$sex))==1 & length(unique(temp_chemo$diagnosis))==1){
        temp_chemo1<-lm(temp_chemo[,1] ~ amplicon + time + age + drug_chemo,temp_chemo)
        if (rownames(summary(temp_chemo1)$coefficients)[nrow(summary(temp_chemo1)$coefficients)] == "drug_chemo"){
          estimate_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),1]
          Std_Error_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),2]
          t_value_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),3]
          p_value_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),4]
        }else{
          estimate_chemo[x,y]<-NA
          Std_Error_chemo[x,y]<-NA
          t_value_chemo[x,y]<-NA
          p_value_chemo[x,y]<-NA
        }
      }else if(length(unique(temp_chemo$time))==1 & length(unique(temp_chemo$amplicon))==1 & length(unique(temp_chemo$sex))==1 & length(unique(temp_chemo$diagnosis))!=1){
        temp_chemo1<-lm(temp_chemo[,1] ~ age + diagnosis + drug_chemo,temp_chemo)
        if (rownames(summary(temp_chemo1)$coefficients)[nrow(summary(temp_chemo1)$coefficients)] == "drug_chemo"){
          estimate_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),1]
          Std_Error_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),2]
          t_value_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),3]
          p_value_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),4]
        }else{
          estimate_chemo[x,y]<-NA
          Std_Error_chemo[x,y]<-NA
          t_value_chemo[x,y]<-NA
          p_value_chemo[x,y]<-NA
        }
      }else if(length(unique(temp_chemo$time))==1 & length(unique(temp_chemo$amplicon))==1 & length(unique(temp_chemo$sex))!=1 & length(unique(temp_chemo$diagnosis))==1){
        temp_chemo1<-lm(temp_chemo[,1] ~ sex + age + drug_chemo,temp_chemo)
        if (rownames(summary(temp_chemo1)$coefficients)[nrow(summary(temp_chemo1)$coefficients)] == "drug_chemo"){
          estimate_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),1]
          Std_Error_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),2]
          t_value_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),3]
          p_value_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),4]
        }else{
          estimate_chemo[x,y]<-NA
          Std_Error_chemo[x,y]<-NA
          t_value_chemo[x,y]<-NA
          p_value_chemo[x,y]<-NA
        }
      }else if(length(unique(temp_chemo$time))==1 & length(unique(temp_chemo$amplicon))!=1 & length(unique(temp_chemo$sex))==1 & length(unique(temp_chemo$diagnosis))==1){
        temp_chemo1<-lm(temp_chemo[,1] ~ amplicon + age + drug_chemo,temp_chemo)
        if (rownames(summary(temp_chemo1)$coefficients)[nrow(summary(temp_chemo1)$coefficients)] == "drug_chemo"){
          estimate_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),1]
          Std_Error_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),2]
          t_value_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),3]
          p_value_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),4]
        }else{
          estimate_chemo[x,y]<-NA
          Std_Error_chemo[x,y]<-NA
          t_value_chemo[x,y]<-NA
          p_value_chemo[x,y]<-NA
        }
      }else if(length(unique(temp_chemo$time))!=1 & length(unique(temp_chemo$amplicon))==1 & length(unique(temp_chemo$sex))==1 & length(unique(temp_chemo$diagnosis))==1){
        temp_chemo1<-lm(temp_chemo[,1] ~ time + age + drug_chemo,temp_chemo)
        if (rownames(summary(temp_chemo1)$coefficients)[nrow(summary(temp_chemo1)$coefficients)] == "drug_chemo"){
          estimate_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),1]
          Std_Error_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),2]
          t_value_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),3]
          p_value_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),4]
        }else{
          estimate_chemo[x,y]<-NA
          Std_Error_chemo[x,y]<-NA
          t_value_chemo[x,y]<-NA
          p_value_chemo[x,y]<-NA
        }
      }else{
        temp_chemo1<-lm(temp_chemo[,1] ~ age + drug_chemo,temp_chemo)
        if (rownames(summary(temp_chemo1)$coefficients)[nrow(summary(temp_chemo1)$coefficients)] == "drug_chemo"){
          estimate_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),1]
          Std_Error_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),2]
          t_value_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),3]
          p_value_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),4]
        }else{
          estimate_chemo[x,y]<-NA
          Std_Error_chemo[x,y]<-NA
          t_value_chemo[x,y]<-NA
          p_value_chemo[x,y]<-NA
        }
      }
    }else{
      estimate_chemo[x,y]<-NA
      Std_Error_chemo[x,y]<-NA
      t_value_chemo[x,y]<-NA
      p_value_chemo[x,y]<-NA
    }
  }
}

number_chemo<-as.data.frame(number_chemo)
nonzero_chemo<-as.data.frame(nonzero_chemo)
users_chemo<-as.data.frame(users_chemo)
nonusers_chemo<-as.data.frame(nonusers_chemo)

estimate_chemo<-as.data.frame(estimate_chemo)
p_value_chemo<-as.data.frame(p_value_chemo)
Std_Error_chemo<-as.data.frame(Std_Error_chemo)
t_value_chemo<-as.data.frame(t_value_chemo)


```

### multiple drug analysis

```{r}
number_mul_chemo_drug<-array(length(unique(chemo_reg_name$subject)),dim=c(98,22))
rownames(number_mul_chemo_drug)<-colnames(chemo_reg_name[,2:99])
colnames(number_mul_chemo_drug)<-colnames(chemo_reg_name[,115:136])
nonzero_mul_chemo_drug<-array(NA,dim=c(98,22))
rownames(nonzero_mul_chemo_drug)<-colnames(chemo_reg_name[,2:99])
colnames(nonzero_mul_chemo_drug)<-colnames(chemo_reg_name[,115:136])
users_mul_chemo_drug<-array(NA,dim=c(98,22))
rownames(users_mul_chemo_drug)<-colnames(chemo_reg_name[,2:99])
colnames(users_mul_chemo_drug)<-colnames(chemo_reg_name[,115:136])
nonusers_mul_chemo_drug<-array(NA,dim=c(98,22))
rownames(nonusers_mul_chemo_drug)<-colnames(chemo_reg_name[,2:99])
colnames(nonusers_mul_chemo_drug)<-colnames(chemo_reg_name[,115:136])


estimate_mul_chemo_drug<-array(0,dim=c(98,22))
rownames(estimate_mul_chemo_drug)<-colnames(chemo_reg_name[,2:99])
colnames(estimate_mul_chemo_drug)<-colnames(chemo_reg_name[,115:136])
Std_Error_mul_chemo_drug<-array(0,dim=c(98,22))
rownames(Std_Error_mul_chemo_drug)<-colnames(chemo_reg_name[,2:99])
colnames(Std_Error_mul_chemo_drug)<-colnames(chemo_reg_name[,115:136])
t_value_mul_chemo_drug<-array(0,dim=c(98,22))
rownames(t_value_mul_chemo_drug)<-colnames(chemo_reg_name[,2:99])
colnames(t_value_mul_chemo_drug)<-colnames(chemo_reg_name[,115:136])
p_value_mul_chemo_drug<-array(0,dim=c(98,22))
rownames(p_value_mul_chemo_drug)<-colnames(chemo_reg_name[,2:99])
colnames(p_value_mul_chemo_drug)<-colnames(chemo_reg_name[,115:136])


for (x in 1:98) {
  for (y in 1:22){
    genus_name<-colnames(chemo_reg_name[,2:99])
    drug_name<-colnames(chemo_reg_name[,115:136])
    
    temp_mul_chemo_drug<-select(chemo_reg_name,genus_name[x],subject,flag,time,amplicon,sex,diagnosis,Trimethoprim.sulfamethoxazole,drug_name[y])
    temp_mul_chemo_drug$species<-temp_mul_chemo_drug[,1]
    temp_mul_chemo_nonzero_drug<- temp_mul_chemo_drug |> filter(species!=0) 
    nonzero_mul_chemo_drug[x,y] <-length(unique(temp_mul_chemo_nonzero_drug$subject))
    temp_mul_chemo_drug$drug<-temp_mul_chemo_drug[,(ncol(temp_mul_chemo_drug)-1)]
    temp_mul_chemo_users_drug<-temp_mul_chemo_drug[complete.cases(temp_mul_chemo_drug$drug),] |> filter(drug==1) 
    users_mul_chemo_drug[x,y] <- length(unique(temp_mul_chemo_users_drug$subject))
    temp_mul_chemo_nonusers_drug<-temp_mul_chemo_drug[complete.cases(temp_mul_chemo_drug$drug),] |> filter(drug==0) 
    nonusers_mul_chemo_drug[x,y] <- length(unique(temp_mul_chemo_nonusers_drug$subject))
    
    
    temp_mul_chemo_drug<-temp_mul_chemo_drug[complete.cases(temp_mul_chemo_drug),]
    drug<-temp_mul_chemo_drug[,ncol(temp_mul_chemo_drug)]
    Antibiotic<-temp_mul_chemo_drug[,"Trimethoprim.sulfamethoxazole"]
    if(length(unique(drug))==2 & length(unique(Antibiotic))==2){
      if(length(unique(temp_mul_chemo_drug$time))!=1 & length(unique(temp_mul_chemo_drug$amplicon))!=1 & length(unique(temp_mul_chemo_drug$sex))!=1 & length(unique(temp_mul_chemo_drug$diagnosis))!=1){
        temp_mul_chemo_drug1<-lm(temp_mul_chemo_drug[,1] ~ time + amplicon + sex +  diagnosis + Antibiotic +drug,temp_mul_chemo_drug)
        if (rownames(summary(temp_mul_chemo_drug1)$coefficients)[nrow(summary(temp_mul_chemo_drug1)$coefficients)] == "drug"){
          summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),1]->estimate_mul_chemo_drug[x,y]
          summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),2]->Std_Error_mul_chemo_drug[x,y]
          summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),3]->t_value_mul_chemo_drug[x,y]
          summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),4]->p_value_mul_chemo_drug[x,y]
        }else{
          estimate_mul_chemo_drug[x,y]<-NA
          Std_Error_mul_chemo_drug[x,y]<-NA
          t_value_mul_chemo_drug[x,y]<-NA
          p_value_mul_chemo_drug[x,y]<-NA
        }
      }else if(length(unique(temp_mul_chemo_drug$time))==1 & length(unique(temp_mul_chemo_drug$amplicon))!=1 & length(unique(temp_mul_chemo_drug$sex))!=1 & length(unique(temp_mul_chemo_drug$diagnosis))!=1){
        temp_mul_chemo_drug1<-lm(temp_mul_chemo_drug[,1] ~ amplicon + sex +  diagnosis + Antibiotic+ drug,temp_mul_chemo_drug)
        if (rownames(summary(temp_mul_chemo_drug1)$coefficients)[nrow(summary(temp_mul_chemo_drug1)$coefficients)] == "drug"){
          estimate_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),1]
          Std_Error_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),2]
          t_value_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),3]
          p_value_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),4]
        }else{
          estimate_mul_chemo_drug[x,y]<-NA
          Std_Error_mul_chemo_drug[x,y]<-NA
          t_value_mul_chemo_drug[x,y]<-NA
          p_value_mul_chemo_drug[x,y]<-NA
        }
      }else if(length(unique(temp_mul_chemo_drug$time))!=1 & length(unique(temp_mul_chemo_drug$amplicon))==1 & length(unique(temp_mul_chemo_drug$sex))!=1 & length(unique(temp_mul_chemo_drug$diagnosis))!=1){
        temp_mul_chemo_drug1<-lm(temp_mul_chemo_drug[,1] ~ time + sex +  diagnosis + Antibiotic+ drug,temp_mul_chemo_drug)
        if (rownames(summary(temp_mul_chemo_drug1)$coefficients)[nrow(summary(temp_mul_chemo_drug1)$coefficients)] == "drug"){
          estimate_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),1]
          Std_Error_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),2]
          t_value_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),3]
          p_value_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),4]
        }else{
          estimate_mul_chemo_drug[x,y]<-NA
          Std_Error_mul_chemo_drug[x,y]<-NA
          t_value_mul_chemo_drug[x,y]<-NA
          p_value_mul_chemo_drug[x,y]<-NA
        }
      }else if(length(unique(temp_mul_chemo_drug$time))!=1 & length(unique(temp_mul_chemo_drug$amplicon))!=1 & length(unique(temp_mul_chemo_drug$sex))==1 & length(unique(temp_mul_chemo_drug$diagnosis))!=1){
        temp_mul_chemo_drug1<-lm(temp_mul_chemo_drug[,1] ~ time  +  diagnosis + Antibiotic+ drug,temp_mul_chemo_drug)
        if (rownames(summary(temp_mul_chemo_drug1)$coefficients)[nrow(summary(temp_mul_chemo_drug1)$coefficients)] == "drug"){
          estimate_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),1]
          Std_Error_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),2]
          t_value_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),3]
          p_value_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),4]
        }else{
          estimate_mul_chemo_drug[x,y]<-NA
          Std_Error_mul_chemo_drug[x,y]<-NA
          t_value_mul_chemo_drug[x,y]<-NA
          p_value_mul_chemo_drug[x,y]<-NA
        }
      }else if(length(unique(temp_mul_chemo_drug$time))!=1 & length(unique(temp_mul_chemo_drug$amplicon))!=1 & length(unique(temp_mul_chemo_drug$sex))!=1 & length(unique(temp_mul_chemo_drug$diagnosis))==1){
        temp_mul_chemo_drug1<-lm(temp_mul_chemo_drug[,1] ~ time + sex +  Antibiotic + drug,temp_mul_chemo_drug)
        if (rownames(summary(temp_mul_chemo_drug1)$coefficients)[nrow(summary(temp_mul_chemo_drug1)$coefficients)] == "drug"){
          estimate_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),1]
          Std_Error_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),2]
          t_value_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),3]
          p_value_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),4]
        }else{
          estimate_mul_chemo_drug[x,y]<-NA
          Std_Error_mul_chemo_drug[x,y]<-NA
          t_value_mul_chemo_drug[x,y]<-NA
          p_value_mul_chemo_drug[x,y]<-NA
        }
      }else if(length(unique(temp_mul_chemo_drug$time))==1 & length(unique(temp_mul_chemo_drug$amplicon))==1 & length(unique(temp_mul_chemo_drug$sex))!=1 & length(unique(temp_mul_chemo_drug$diagnosis))!=1){
        temp_mul_chemo_drug1<-lm(temp_mul_chemo_drug[,1] ~ sex +  diagnosis + Antibiotic+ drug,temp_mul_chemo_drug)
        if (rownames(summary(temp_mul_chemo_drug1)$coefficients)[nrow(summary(temp_mul_chemo_drug1)$coefficients)] == "drug"){
          estimate_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),1]
          Std_Error_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),2]
          t_value_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),3]
          p_value_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),4]
        }else{
          estimate_mul_chemo_drug[x,y]<-NA
          Std_Error_mul_chemo_drug<-NA
          t_value_mul_chemo_drug[x,y]<-NA
          p_value_mul_chemo_drug[x,y]<-NA
        }
      }else if(length(unique(temp_mul_chemo_drug$time))==1 & length(unique(temp_mul_chemo_drug$amplicon))!=1 & length(unique(temp_mul_chemo_drug$sex))==1 & length(unique(temp_mul_chemo_drug$diagnosis))!=1){
        temp_mul_chemo_drug1<-lm(temp_mul_chemo_drug[,1] ~ amplicon +  diagnosis + Antibiotic+ drug,temp_mul_chemo_drug)
        if (rownames(summary(temp_mul_chemo_drug1)$coefficients)[nrow(summary(temp_mul_chemo_drug1)$coefficients)] == "drug"){
          estimate_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),1]
          Std_Error_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),2]
          t_value_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),3]
          p_value_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),4]
        }else{
          estimate_mul_chemo_drug[x,y]<-NA
          Std_Error_mul_chemo_drug[x,y]<-NA
          t_value_mul_chemo_drug[x,y]<-NA
          p_value_mul_chemo_drug[x,y]<-NA
        }
      }else if(length(unique(temp_mul_chemo_drug$time))==1 & length(unique(temp_mul_chemo_drug$amplicon))!=1 & length(unique(temp_mul_chemo_drug$sex))!=1 & length(unique(temp_mul_chemo_drug$diagnosis))==1){
        temp_mul_chemo_drug1<-lm(temp_mul_chemo_drug[,1] ~ amplicon + sex +  Antibiotic+ drug,temp_mul_chemo_drug)
        if (rownames(summary(temp_mul_chemo_drug1)$coefficients)[nrow(summary(temp_mul_chemo_drug1)$coefficients)] == "drug"){
          estimate_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),1]
          Std_Error_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),2]
          t_value_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),3]
          p_value_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),4]
        }else{
          estimate_mul_chemo_drug[x,y]<-NA
          Std_Error_mul_chemo_drug[x,y]<-NA
          t_value_mul_chemo_drug[x,y]<-NA
          p_value_mul_chemo_drug[x,y]<-NA
        }
      }else if(length(unique(temp_mul_chemo_drug$time))!=1 & length(unique(temp_mul_chemo_drug$amplicon))==1 & length(unique(temp_mul_chemo_drug$sex))==1 & length(unique(temp_mul_chemo_drug$diagnosis))!=1){
        temp_mul_chemo_drug1<-lm(temp_mul_chemo_drug[,1] ~ time +  diagnosis + Antibiotic+ drug,temp_mul_chemo_drug)
        if (rownames(summary(temp_mul_chemo_drug1)$coefficients)[nrow(summary(temp_mul_chemo_drug1)$coefficients)] == "drug"){
          estimate_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),1]
          Std_Error_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),2]
          t_value_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),3]
          p_value_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),4]
        }else{
          estimate_mul_chemo_drug[x,y]<-NA
          Std_Error_mul_chemo_drug[x,y]<-NA
          t_value_mul_chemo_drug[x,y]<-NA
          p_value_mul_chemo_drug[x,y]<-NA
        }
      }else if(length(unique(temp_mul_chemo_drug$time))!=1 & length(unique(temp_mul_chemo_drug$amplicon))==1 & length(unique(temp_mul_chemo_drug$sex))!=1 & length(unique(temp_mul_chemo_drug$diagnosis))==1){
        temp_mul_chemo_drug1<-lm(temp_mul_chemo_drug[,1] ~ time + sex +  Antibiotic+ drug,temp_mul_chemo_drug)
        if (rownames(summary(temp_mul_chemo_drug1)$coefficients)[nrow(summary(temp_mul_chemo_drug1)$coefficients)] == "drug"){
          estimate_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),1]
          Std_Error_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),2]
          t_value_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),3]
          p_value_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),4]
        }else{
          estimate_mul_chemo_drug[x,y]<-NA
          Std_Error_mul_chemo_drug[x,y]<-NA
          t_value_mul_chemo_drug[x,y]<-NA
          p_value_mul_chemo_drug[x,y]<-NA
        }
      }else if(length(unique(temp_mul_chemo_drug$time))!=1 & length(unique(temp_mul_chemo_drug$amplicon))!=1 & length(unique(temp_mul_chemo_drug$sex))==1 & length(unique(temp_mul_chemo_drug$diagnosis))==1){
        temp_mul_chemo_drug1<-lm(temp_mul_chemo_drug[,1] ~ amplicon + time +  Antibiotic+ drug,temp_mul_chemo_drug)
        if (rownames(summary(temp_mul_chemo_drug1)$coefficients)[nrow(summary(temp_mul_chemo_drug1)$coefficients)] == "drug"){
          estimate_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),1]
          Std_Error_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),2]
          t_value_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),3]
          p_value_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),4]
        }else{
          estimate_mul_chemo_drug[x,y]<-NA
          Std_Error_mul_chemo_drug[x,y]<-NA
          t_value_mul_chemo_drug[x,y]<-NA
          p_value_mul_chemo_drug[x,y]<-NA
        }
      }else if(length(unique(temp_mul_chemo_drug$time))==1 & length(unique(temp_mul_chemo_drug$amplicon))==1 & length(unique(temp_mul_chemo_drug$sex))==1 & length(unique(temp_mul_chemo_drug$diagnosis))!=1){
        temp_mul_chemo_drug1<-lm(temp_mul_chemo_drug[,1] ~  diagnosis + Antibiotic+ drug,temp_mul_chemo_drug)
        if (rownames(summary(temp_mul_chemo_drug1)$coefficients)[nrow(summary(temp_mul_chemo_drug1)$coefficients)] == "drug"){
          estimate_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),1]
          Std_Error_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),2]
          t_value_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),3]
          p_value_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),4]
        }else{
          estimate_mul_chemo_drug[x,y]<-NA
          Std_Error_mul_chemo_drug[x,y]<-NA
          t_value_mul_chemo_drug[x,y]<-NA
          p_value_mul_chemo_drug[x,y]<-NA
        }
      }else if(length(unique(temp_mul_chemo_drug$time))==1 & length(unique(temp_mul_chemo_drug$amplicon))==1 & length(unique(temp_mul_chemo_drug$sex))!=1 & length(unique(temp_mul_chemo_drug$diagnosis))==1){
        temp_mul_chemo_drug1<-lm(temp_mul_chemo_drug[,1] ~ sex +  Antibiotic+ drug,temp_mul_chemo_drug)
        if (rownames(summary(temp_mul_chemo_drug1)$coefficients)[nrow(summary(temp_mul_chemo_drug1)$coefficients)] == "drug"){
          estimate_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),1]
          Std_Error_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),2]
          t_value_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),3]
          p_value_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),4]
        }else{
          estimate_mul_chemo_drug[x,y]<-NA
          Std_Error_mul_chemo_drug[x,y]<-NA
          t_value_mul_chemo_drug[x,y]<-NA
          p_value_mul_chemo_drug[x,y]<-NA
        }
      }else if(length(unique(temp_mul_chemo_drug$time))==1 & length(unique(temp_mul_chemo_drug$amplicon))!=1 & length(unique(temp_mul_chemo_drug$sex))==1 & length(unique(temp_mul_chemo_drug$diagnosis))==1){
        temp_mul_chemo_drug1<-lm(temp_mul_chemo_drug[,1] ~ amplicon +  Antibiotic+ drug,temp_mul_chemo_drug)
        if (rownames(summary(temp_mul_chemo_drug1)$coefficients)[nrow(summary(temp_mul_chemo_drug1)$coefficients)] == "drug"){
          estimate_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),1]
          Std_Error_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),2]
          t_value_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),3]
          p_value_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),4]
        }else{
          estimate_mul_chemo_drug[x,y]<-NA
          Std_Error_mul_chemo_drug[x,y]<-NA
          t_value_mul_chemo_drug[x,y]<-NA
          p_value_mul_chemo_drug[x,y]<-NA
        }
      }else if(length(unique(temp_mul_chemo_drug$time))!=1 & length(unique(temp_mul_chemo_drug$amplicon))==1 & length(unique(temp_mul_chemo_drug$sex))==1 & length(unique(temp_mul_chemo_drug$diagnosis))==1){
        temp_mul_chemo_drug1<-lm(temp_mul_chemo_drug[,1] ~ time +  Antibiotic+ drug,temp_mul_chemo_drug)
        if (rownames(summary(temp_mul_chemo_drug1)$coefficients)[nrow(summary(temp_mul_chemo_drug1)$coefficients)] == "drug"){
          estimate_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),1]
          Std_Error_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),2]
          t_value_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),3]
          p_value_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),4]
        }else{
          estimate_mul_chemo_drug[x,y]<-NA
          Std_Error_mul_chemo_drug[x,y]<-NA
          t_value_mul_chemo_drug[x,y]<-NA
          p_value_mul_chemo_drug[x,y]<-NA
        }
      }else{
        temp_mul_chemo_drug1<-lm(temp_mul_chemo_drug[,1] ~  Antibiotic+ drug,temp_mul_chemo_drug)
        if (rownames(summary(temp_mul_chemo_drug1)$coefficients)[nrow(summary(temp_mul_chemo_drug1)$coefficients)] == "drug"){
          estimate_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),1]
          Std_Error_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),2]
          t_value_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),3]
          p_value_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),4]
        }else{
          estimate_mul_chemo_drug[x,y]<-NA
          Std_Error_mul_chemo_drug[x,y]<-NA
          t_value_mul_chemo_drug[x,y]<-NA
          p_value_mul_chemo_drug[x,y]<-NA
        }
      }
    }else{
      estimate_mul_chemo_drug[x,y]<-NA
      Std_Error_mul_chemo_drug[x,y]<-NA
      t_value_mul_chemo_drug[x,y]<-NA
      p_value_mul_chemo_drug[x,y]<-NA
    }
  }
}

number_mul_chemo_anti<-array(length(unique(chemo_reg_name$subject)),dim=c(98,25))
rownames(number_mul_chemo_anti)<-colnames(chemo_reg_name[,2:99])
colnames(number_mul_chemo_anti)<-colnames(chemo_reg_name[,138:162])
nonzero_mul_chemo_anti<-array(NA,dim=c(98,25))
rownames(nonzero_mul_chemo_anti)<-colnames(chemo_reg_name[,2:99])
colnames(nonzero_mul_chemo_anti)<-colnames(chemo_reg_name[,138:162])
users_mul_chemo_anti<-array(NA,dim=c(98,25))
rownames(users_mul_chemo_anti)<-colnames(chemo_reg_name[,2:99])
colnames(users_mul_chemo_anti)<-colnames(chemo_reg_name[,138:162])
nonusers_mul_chemo_anti<-array(NA,dim=c(98,25))
rownames(nonusers_mul_chemo_anti)<-colnames(chemo_reg_name[,2:99])
colnames(nonusers_mul_chemo_anti)<-colnames(chemo_reg_name[,138:162])


estimate_mul_chemo_anti<-array(0,dim=c(98,25))
rownames(estimate_mul_chemo_anti)<-colnames(chemo_reg_name[,2:99])
colnames(estimate_mul_chemo_anti)<-colnames(chemo_reg_name[,138:162])
Std_Error_mul_chemo_anti<-array(0,dim=c(98,25))
rownames(Std_Error_mul_chemo_anti)<-colnames(chemo_reg_name[,2:99])
colnames(Std_Error_mul_chemo_anti)<-colnames(chemo_reg_name[,138:162])
t_value_mul_chemo_anti<-array(0,dim=c(98,25))
rownames(t_value_mul_chemo_anti)<-colnames(chemo_reg_name[,2:99])
colnames(t_value_mul_chemo_anti)<-colnames(chemo_reg_name[,138:162])
p_value_mul_chemo_anti<-array(0,dim=c(98,25))
rownames(p_value_mul_chemo_anti)<-colnames(chemo_reg_name[,2:99])
colnames(p_value_mul_chemo_anti)<-colnames(chemo_reg_name[,138:162])


for (x in 1:98) {
  for (y in 1:25){
    genus_name<-colnames(chemo_reg_name[,2:99])
    drug_name<-colnames(chemo_reg_name[,138:162])
    
    temp_mul_chemo_anti<-select(chemo_reg_name,genus_name[x],subject,flag,time,amplicon,sex,diagnosis,Antibiotic,drug_name[y])
    temp_mul_chemo_anti$species<-temp_mul_chemo_anti[,1]
    temp_mul_chemo_nonzero_anti<- temp_mul_chemo_anti |> filter(species!=0) 
    nonzero_mul_chemo_anti[x,y] <-length(unique(temp_mul_chemo_nonzero_anti$subject))
    temp_mul_chemo_anti$drug<-temp_mul_chemo_anti[,(ncol(temp_mul_chemo_anti)-1)]
    temp_mul_chemo_users_anti<-temp_mul_chemo_anti[complete.cases(temp_mul_chemo_anti$drug),] |> filter(drug==1) 
    users_mul_chemo_anti[x,y] <- length(unique(temp_mul_chemo_users_anti$subject))
    temp_mul_chemo_nonusers_anti<-temp_mul_chemo_anti[complete.cases(temp_mul_chemo_anti$drug),] |> filter(drug==0) 
    nonusers_mul_chemo_anti[x,y] <- length(unique(temp_mul_chemo_nonusers_anti$subject))
    
    temp_mul_chemo_anti<-temp_mul_chemo_anti[complete.cases(temp_mul_chemo_anti),]
    drug<-temp_mul_chemo_anti[,ncol(temp_mul_chemo_anti)]
    Antibiotic<-temp_mul_chemo_anti[,"Antibiotic"]
    if(length(unique(drug))==2 & length(unique(Antibiotic))==2){
      if(length(unique(temp_mul_chemo_anti$time))!=1 & length(unique(temp_mul_chemo_anti$amplicon))!=1 & length(unique(temp_mul_chemo_anti$sex))!=1 & length(unique(temp_mul_chemo_anti$diagnosis))!=1){
        temp_mul_chemo_anti1<-lm(temp_mul_chemo_anti[,1] ~ time + amplicon + sex +  diagnosis + Antibiotic +drug,temp_mul_chemo_anti)
        if (rownames(summary(temp_mul_chemo_anti1)$coefficients)[nrow(summary(temp_mul_chemo_anti1)$coefficients)] == "drug"){
          summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),1]->estimate_mul_chemo_anti[x,y]
          summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),2]->Std_Error_mul_chemo_anti[x,y]
          summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),3]->t_value_mul_chemo_anti[x,y]
          summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),4]->p_value_mul_chemo_anti[x,y]
        }else{
          estimate_mul_chemo_anti[x,y]<-NA
          Std_Error_mul_chemo_anti[x,y]<-NA
          t_value_mul_chemo_anti[x,y]<-NA
          p_value_mul_chemo_anti[x,y]<-NA
        }
      }else if(length(unique(temp_mul_chemo_anti$time))==1 & length(unique(temp_mul_chemo_anti$amplicon))!=1 & length(unique(temp_mul_chemo_anti$sex))!=1 & length(unique(temp_mul_chemo_anti$diagnosis))!=1){
        temp_mul_chemo_anti1<-lm(temp_mul_chemo_anti[,1] ~ amplicon + sex +  diagnosis + Antibiotic+ drug,temp_mul_chemo_anti)
        if (rownames(summary(temp_mul_chemo_anti1)$coefficients)[nrow(summary(temp_mul_chemo_anti1)$coefficients)] == "drug"){
          estimate_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),1]
          Std_Error_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),2]
          t_value_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),3]
          p_value_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),4]
        }else{
          estimate_mul_chemo_anti[x,y]<-NA
          Std_Error_mul_chemo_anti[x,y]<-NA
          t_value_mul_chemo_anti[x,y]<-NA
          p_value_mul_chemo_anti[x,y]<-NA
        }
      }else if(length(unique(temp_mul_chemo_anti$time))!=1 & length(unique(temp_mul_chemo_anti$amplicon))==1 & length(unique(temp_mul_chemo_anti$sex))!=1 & length(unique(temp_mul_chemo_anti$diagnosis))!=1){
        temp_mul_chemo_anti1<-lm(temp_mul_chemo_anti[,1] ~ time + sex +  diagnosis + Antibiotic+ drug,temp_mul_chemo_anti)
        if (rownames(summary(temp_mul_chemo_anti1)$coefficients)[nrow(summary(temp_mul_chemo_anti1)$coefficients)] == "drug"){
          estimate_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),1]
          Std_Error_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),2]
          t_value_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),3]
          p_value_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),4]
        }else{
          estimate_mul_chemo_anti[x,y]<-NA
          Std_Error_mul_chemo_anti[x,y]<-NA
          t_value_mul_chemo_anti[x,y]<-NA
          p_value_mul_chemo_anti[x,y]<-NA
        }
      }else if(length(unique(temp_mul_chemo_anti$time))!=1 & length(unique(temp_mul_chemo_anti$amplicon))!=1 & length(unique(temp_mul_chemo_anti$sex))==1 & length(unique(temp_mul_chemo_anti$diagnosis))!=1){
        temp_mul_chemo_anti1<-lm(temp_mul_chemo_anti[,1] ~ time  +  diagnosis + Antibiotic+ drug,temp_mul_chemo_anti)
        if (rownames(summary(temp_mul_chemo_anti1)$coefficients)[nrow(summary(temp_mul_chemo_anti1)$coefficients)] == "drug"){
          estimate_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),1]
          Std_Error_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),2]
          t_value_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),3]
          p_value_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),4]
        }else{
          estimate_mul_chemo_anti[x,y]<-NA
          Std_Error_mul_chemo_anti[x,y]<-NA
          t_value_mul_chemo_anti[x,y]<-NA
          p_value_mul_chemo_anti[x,y]<-NA
        }
      }else if(length(unique(temp_mul_chemo_anti$time))!=1 & length(unique(temp_mul_chemo_anti$amplicon))!=1 & length(unique(temp_mul_chemo_anti$sex))!=1 & length(unique(temp_mul_chemo_anti$diagnosis))==1){
        temp_mul_chemo_anti1<-lm(temp_mul_chemo_anti[,1] ~ time + sex +  Antibiotic + drug,temp_mul_chemo_anti)
        if (rownames(summary(temp_mul_chemo_anti1)$coefficients)[nrow(summary(temp_mul_chemo_anti1)$coefficients)] == "drug"){
          estimate_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),1]
          Std_Error_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),2]
          t_value_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),3]
          p_value_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),4]
        }else{
          estimate_mul_chemo_anti[x,y]<-NA
          Std_Error_mul_chemo_anti[x,y]<-NA
          t_value_mul_chemo_anti[x,y]<-NA
          p_value_mul_chemo_anti[x,y]<-NA
        }
      }else if(length(unique(temp_mul_chemo_anti$time))==1 & length(unique(temp_mul_chemo_anti$amplicon))==1 & length(unique(temp_mul_chemo_anti$sex))!=1 & length(unique(temp_mul_chemo_anti$diagnosis))!=1){
        temp_mul_chemo_anti1<-lm(temp_mul_chemo_anti[,1] ~ sex +  diagnosis + Antibiotic+ drug,temp_mul_chemo_anti)
        if (rownames(summary(temp_mul_chemo_anti1)$coefficients)[nrow(summary(temp_mul_chemo_anti1)$coefficients)] == "drug"){
          estimate_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),1]
          Std_Error_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),2]
          t_value_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),3]
          p_value_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),4]
        }else{
          estimate_mul_chemo_anti[x,y]<-NA
          Std_Error_mul_chemo_anti<-NA
          t_value_mul_chemo_anti[x,y]<-NA
          p_value_mul_chemo_anti[x,y]<-NA
        }
      }else if(length(unique(temp_mul_chemo_anti$time))==1 & length(unique(temp_mul_chemo_anti$amplicon))!=1 & length(unique(temp_mul_chemo_anti$sex))==1 & length(unique(temp_mul_chemo_anti$diagnosis))!=1){
        temp_mul_chemo_anti1<-lm(temp_mul_chemo_anti[,1] ~ amplicon +  diagnosis + Antibiotic+ drug,temp_mul_chemo_anti)
        if (rownames(summary(temp_mul_chemo_anti1)$coefficients)[nrow(summary(temp_mul_chemo_anti1)$coefficients)] == "drug"){
          estimate_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),1]
          Std_Error_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),2]
          t_value_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),3]
          p_value_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),4]
        }else{
          estimate_mul_chemo_anti[x,y]<-NA
          Std_Error_mul_chemo_anti[x,y]<-NA
          t_value_mul_chemo_anti[x,y]<-NA
          p_value_mul_chemo_anti[x,y]<-NA
        }
      }else if(length(unique(temp_mul_chemo_anti$time))==1 & length(unique(temp_mul_chemo_anti$amplicon))!=1 & length(unique(temp_mul_chemo_anti$sex))!=1 & length(unique(temp_mul_chemo_anti$diagnosis))==1){
        temp_mul_chemo_anti1<-lm(temp_mul_chemo_anti[,1] ~ amplicon + sex +  Antibiotic+ drug,temp_mul_chemo_anti)
        if (rownames(summary(temp_mul_chemo_anti1)$coefficients)[nrow(summary(temp_mul_chemo_anti1)$coefficients)] == "drug"){
          estimate_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),1]
          Std_Error_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),2]
          t_value_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),3]
          p_value_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),4]
        }else{
          estimate_mul_chemo_anti[x,y]<-NA
          Std_Error_mul_chemo_anti[x,y]<-NA
          t_value_mul_chemo_anti[x,y]<-NA
          p_value_mul_chemo_anti[x,y]<-NA
        }
      }else if(length(unique(temp_mul_chemo_anti$time))!=1 & length(unique(temp_mul_chemo_anti$amplicon))==1 & length(unique(temp_mul_chemo_anti$sex))==1 & length(unique(temp_mul_chemo_anti$diagnosis))!=1){
        temp_mul_chemo_anti1<-lm(temp_mul_chemo_anti[,1] ~ time +  diagnosis + Antibiotic+ drug,temp_mul_chemo_anti)
        if (rownames(summary(temp_mul_chemo_anti1)$coefficients)[nrow(summary(temp_mul_chemo_anti1)$coefficients)] == "drug"){
          estimate_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),1]
          Std_Error_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),2]
          t_value_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),3]
          p_value_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),4]
        }else{
          estimate_mul_chemo_anti[x,y]<-NA
          Std_Error_mul_chemo_anti[x,y]<-NA
          t_value_mul_chemo_anti[x,y]<-NA
          p_value_mul_chemo_anti[x,y]<-NA
        }
      }else if(length(unique(temp_mul_chemo_anti$time))!=1 & length(unique(temp_mul_chemo_anti$amplicon))==1 & length(unique(temp_mul_chemo_anti$sex))!=1 & length(unique(temp_mul_chemo_anti$diagnosis))==1){
        temp_mul_chemo_anti1<-lm(temp_mul_chemo_anti[,1] ~ time + sex +  Antibiotic+ drug,temp_mul_chemo_anti)
        if (rownames(summary(temp_mul_chemo_anti1)$coefficients)[nrow(summary(temp_mul_chemo_anti1)$coefficients)] == "drug"){
          estimate_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),1]
          Std_Error_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),2]
          t_value_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),3]
          p_value_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),4]
        }else{
          estimate_mul_chemo_anti[x,y]<-NA
          Std_Error_mul_chemo_anti[x,y]<-NA
          t_value_mul_chemo_anti[x,y]<-NA
          p_value_mul_chemo_anti[x,y]<-NA
        }
      }else if(length(unique(temp_mul_chemo_anti$time))!=1 & length(unique(temp_mul_chemo_anti$amplicon))!=1 & length(unique(temp_mul_chemo_anti$sex))==1 & length(unique(temp_mul_chemo_anti$diagnosis))==1){
        temp_mul_chemo_anti1<-lm(temp_mul_chemo_anti[,1] ~ amplicon + time +  Antibiotic+ drug,temp_mul_chemo_anti)
        if (rownames(summary(temp_mul_chemo_anti1)$coefficients)[nrow(summary(temp_mul_chemo_anti1)$coefficients)] == "drug"){
          estimate_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),1]
          Std_Error_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),2]
          t_value_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),3]
          p_value_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),4]
        }else{
          estimate_mul_chemo_anti[x,y]<-NA
          Std_Error_mul_chemo_anti[x,y]<-NA
          t_value_mul_chemo_anti[x,y]<-NA
          p_value_mul_chemo_anti[x,y]<-NA
        }
      }else if(length(unique(temp_mul_chemo_anti$time))==1 & length(unique(temp_mul_chemo_anti$amplicon))==1 & length(unique(temp_mul_chemo_anti$sex))==1 & length(unique(temp_mul_chemo_anti$diagnosis))!=1){
        temp_mul_chemo_anti1<-lm(temp_mul_chemo_anti[,1] ~  diagnosis + Antibiotic+ drug,temp_mul_chemo_anti)
        if (rownames(summary(temp_mul_chemo_anti1)$coefficients)[nrow(summary(temp_mul_chemo_anti1)$coefficients)] == "drug"){
          estimate_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),1]
          Std_Error_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),2]
          t_value_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),3]
          p_value_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),4]
        }else{
          estimate_mul_chemo_anti[x,y]<-NA
          Std_Error_mul_chemo_anti[x,y]<-NA
          t_value_mul_chemo_anti[x,y]<-NA
          p_value_mul_chemo_anti[x,y]<-NA
        }
      }else if(length(unique(temp_mul_chemo_anti$time))==1 & length(unique(temp_mul_chemo_anti$amplicon))==1 & length(unique(temp_mul_chemo_anti$sex))!=1 & length(unique(temp_mul_chemo_anti$diagnosis))==1){
        temp_mul_chemo_anti1<-lm(temp_mul_chemo_anti[,1] ~ sex +  Antibiotic+ drug,temp_mul_chemo_anti)
        if (rownames(summary(temp_mul_chemo_anti1)$coefficients)[nrow(summary(temp_mul_chemo_anti1)$coefficients)] == "drug"){
          estimate_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),1]
          Std_Error_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),2]
          t_value_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),3]
          p_value_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),4]
        }else{
          estimate_mul_chemo_anti[x,y]<-NA
          Std_Error_mul_chemo_anti[x,y]<-NA
          t_value_mul_chemo_anti[x,y]<-NA
          p_value_mul_chemo_anti[x,y]<-NA
        }
      }else if(length(unique(temp_mul_chemo_anti$time))==1 & length(unique(temp_mul_chemo_anti$amplicon))!=1 & length(unique(temp_mul_chemo_anti$sex))==1 & length(unique(temp_mul_chemo_anti$diagnosis))==1){
        temp_mul_chemo_anti1<-lm(temp_mul_chemo_anti[,1] ~ amplicon +  Antibiotic+ drug,temp_mul_chemo_anti)
        if (rownames(summary(temp_mul_chemo_anti1)$coefficients)[nrow(summary(temp_mul_chemo_anti1)$coefficients)] == "drug"){
          estimate_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),1]
          Std_Error_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),2]
          t_value_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),3]
          p_value_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),4]
        }else{
          estimate_mul_chemo_anti[x,y]<-NA
          Std_Error_mul_chemo_anti[x,y]<-NA
          t_value_mul_chemo_anti[x,y]<-NA
          p_value_mul_chemo_anti[x,y]<-NA
        }
      }else if(length(unique(temp_mul_chemo_anti$time))!=1 & length(unique(temp_mul_chemo_anti$amplicon))==1 & length(unique(temp_mul_chemo_anti$sex))==1 & length(unique(temp_mul_chemo_anti$diagnosis))==1){
        temp_mul_chemo_anti1<-lm(temp_mul_chemo_anti[,1] ~ time +  Antibiotic+ drug,temp_mul_chemo_anti)
        if (rownames(summary(temp_mul_chemo_anti1)$coefficients)[nrow(summary(temp_mul_chemo_anti1)$coefficients)] == "drug"){
          estimate_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),1]
          Std_Error_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),2]
          t_value_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),3]
          p_value_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),4]
        }else{
          estimate_mul_chemo_anti[x,y]<-NA
          Std_Error_mul_chemo_anti[x,y]<-NA
          t_value_mul_chemo_anti[x,y]<-NA
          p_value_mul_chemo_anti[x,y]<-NA
        }
      }else{
        temp_mul_chemo_anti1<-lm(temp_mul_chemo_anti[,1] ~  Antibiotic+ drug,temp_mul_chemo_anti)
        if (rownames(summary(temp_mul_chemo_anti1)$coefficients)[nrow(summary(temp_mul_chemo_anti1)$coefficients)] == "drug"){
          estimate_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),1]
          Std_Error_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),2]
          t_value_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),3]
          p_value_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),4]
        }else{
          estimate_mul_chemo_anti[x,y]<-NA
          Std_Error_mul_chemo_anti[x,y]<-NA
          t_value_mul_chemo_anti[x,y]<-NA
          p_value_mul_chemo_anti[x,y]<-NA
        }
      }
    }else{
      estimate_mul_chemo_anti[x,y]<-NA
      Std_Error_mul_chemo_anti[x,y]<-NA
      t_value_mul_chemo_anti[x,y]<-NA
      p_value_mul_chemo_anti[x,y]<-NA
    }
  }
}


number_mul_chemo<-cbind(as.data.frame(number_mul_chemo_drug),as.data.frame(number_mul_chemo_anti))
nonzero_mul_chemo<-cbind(as.data.frame(nonzero_mul_chemo_drug),as.data.frame(nonzero_mul_chemo_anti))
users_mul_chemo<-cbind(as.data.frame(users_mul_chemo_drug),as.data.frame(users_mul_chemo_anti))
nonusers_mul_chemo<-cbind(as.data.frame(nonusers_mul_chemo_drug),as.data.frame(nonusers_mul_chemo_anti))

estimate_mul_chemo<-cbind(as.data.frame(estimate_mul_chemo_drug),as.data.frame(estimate_mul_chemo_anti))
p_value_mul_chemo<-cbind(as.data.frame(p_value_mul_chemo_drug),as.data.frame(p_value_mul_chemo_anti))
Std_Error_mul_chemo<-cbind(as.data.frame(Std_Error_mul_chemo_drug),as.data.frame(Std_Error_mul_chemo_anti))
t_value_mul_chemo<-cbind(as.data.frame(t_value_mul_chemo_drug),as.data.frame(t_value_mul_chemo_anti))

```

# 5.Meta-analysis

## Data processing

```{r}

taxa_control<-rownames(estimate)
taxa_gvhd<-rownames(estimate_gvhd)
taxa_chemo<-rownames(estimate_chemo)
taxa_mul_control<-rownames(estimate_mul)
taxa_mul_gvhd<-rownames(estimate_mul_gvhd)
taxa_mul_chemo<-rownames(estimate_mul_chemo)

taxa_all<-union(taxa_control,taxa_chemo)
taxa_all<-union(taxa_all,taxa_gvhd)
taxa_mul_all<-union(taxa_mul_control,taxa_mul_chemo)
taxa_mul_all<-union(taxa_mul_all,taxa_mul_gvhd)

drug_all<-union(colnames(estimate),colnames(estimate_gvhd))
drug_all<-union(colnames(estimate_chemo),drug_all)
drug_mul_all<-union(colnames(estimate_mul),colnames(estimate_mul_gvhd))
drug_mul_all<-union(colnames(estimate_mul_chemo),drug_mul_all)


est_control_bu_r<-as.data.frame(array(NA,dim=c(23,48)))
rownames(est_control_bu_r)<-c(setdiff(taxa_all,taxa_control))
colnames(est_control_bu_r)<-colnames(estimate)
estimate_control_full<-rbind(estimate,est_control_bu_r)
std_control_full<-rbind(Std_Error,est_control_bu_r)
pva_control_full<-rbind(p_value,est_control_bu_r)
tva_control_full<-rbind(t_value,est_control_bu_r)
number_hsct_full<-rbind(number_hsct,est_control_bu_r)
nonzero_hsct_full<-rbind(nonzero_hsct,est_control_bu_r)
users_hsct_full<-rbind(users_hsct,est_control_bu_r)
nonusers_hsct_full<-rbind(nonusers_hsct,est_control_bu_r)

est_mul_control_bu_r<-as.data.frame(array(NA,dim=c(23,47)))
rownames(est_mul_control_bu_r)<-c(setdiff(taxa_mul_all,taxa_mul_control))
colnames(est_mul_control_bu_r)<-colnames(estimate_mul)
estimate_mul_control_full<-rbind(estimate_mul,est_mul_control_bu_r)
std_mul_control_full<-rbind(Std_Error_mul,est_mul_control_bu_r)
pva_mul_control_full<-rbind(p_value_mul,est_mul_control_bu_r)
tva_mul_control_full<-rbind(t_value_mul,est_mul_control_bu_r)
number_mul_hsct_full<-rbind(number_mul_hsct,est_mul_control_bu_r)
nonzero_mul_hsct_full<-rbind(nonzero_mul_hsct,est_mul_control_bu_r)
users_mul_hsct_full<-rbind(users_mul_hsct,est_mul_control_bu_r)
nonusers_mul_hsct_full<-rbind(nonusers_mul_hsct,est_mul_control_bu_r)


est_gvhd_bu_r<-as.data.frame(array(NA,dim=c(16,48)))
rownames(est_gvhd_bu_r)<-c(setdiff(taxa_all,taxa_gvhd))
colnames(est_gvhd_bu_r)<-colnames(estimate_gvhd)
estimate_gvhd_full<-rbind(estimate_gvhd,est_gvhd_bu_r)
std_gvhd_full<-rbind(Std_Error_gvhd,est_gvhd_bu_r)
pva_gvhd_full<-rbind(p_value_gvhd,est_gvhd_bu_r)
tva_gvhd_full<-rbind(t_value_gvhd,est_gvhd_bu_r)
number_gvhd_full<-rbind(number_gvhd,est_gvhd_bu_r)
nonzero_gvhd_full<-rbind(nonzero_gvhd,est_gvhd_bu_r)
users_gvhd_full<-rbind(users_gvhd,est_gvhd_bu_r)
nonusers_gvhd_full<-rbind(nonusers_gvhd,est_gvhd_bu_r)

est_mul_gvhd_bu_r<-as.data.frame(array(NA,dim=c(16,47)))
rownames(est_mul_gvhd_bu_r)<-c(setdiff(taxa_mul_all,taxa_mul_gvhd))
colnames(est_mul_gvhd_bu_r)<-colnames(estimate_mul_gvhd)
estimate_mul_gvhd_full<-rbind(estimate_mul_gvhd,est_mul_gvhd_bu_r)
std_mul_gvhd_full<-rbind(Std_Error_mul_gvhd,est_mul_gvhd_bu_r)
pva_mul_gvhd_full<-rbind(p_value_mul_gvhd,est_mul_gvhd_bu_r)
tva_mul_gvhd_full<-rbind(t_value_mul_gvhd,est_mul_gvhd_bu_r)
number_mul_gvhd_full<-rbind(number_mul_gvhd,est_mul_gvhd_bu_r)
nonzero_mul_gvhd_full<-rbind(nonzero_mul_gvhd,est_mul_gvhd_bu_r)
users_mul_gvhd_full<-rbind(users_mul_gvhd,est_mul_gvhd_bu_r)
nonusers_mul_gvhd_full<-rbind(nonusers_mul_gvhd,est_mul_gvhd_bu_r)



est_chemo_bu_r<-as.data.frame(array(NA,dim=c(29,48)))
rownames(est_chemo_bu_r)<-c(setdiff(taxa_all,taxa_chemo))
colnames(est_chemo_bu_r)<-colnames(estimate_chemo)
estimate_chemo_full<-rbind(estimate_chemo,est_chemo_bu_r)
std_chemo_full<-rbind(Std_Error_chemo,est_chemo_bu_r)
pva_chemo_full<-rbind(p_value_chemo,est_chemo_bu_r)
tva_chemo_full<-rbind(t_value_chemo,est_chemo_bu_r)
number_chemo_full<-rbind(number_chemo,est_chemo_bu_r)
nonzero_chemo_full<-rbind(nonzero_chemo,est_chemo_bu_r)
users_chemo_full<-rbind(users_chemo,est_chemo_bu_r)
nonusers_chemo_full<-rbind(nonusers_chemo,est_chemo_bu_r)

est_mul_chemo_bu_r<-as.data.frame(array(NA,dim=c(29,47)))
rownames(est_mul_chemo_bu_r)<-c(setdiff(taxa_mul_all,taxa_mul_chemo))
colnames(est_mul_chemo_bu_r)<-colnames(estimate_mul_chemo)
estimate_mul_chemo_full<-rbind(estimate_mul_chemo,est_mul_chemo_bu_r)
std_mul_chemo_full<-rbind(Std_Error_mul_chemo,est_mul_chemo_bu_r)
pva_mul_chemo_full<-rbind(p_value_mul_chemo,est_mul_chemo_bu_r)
tva_mul_chemo_full<-rbind(t_value_mul_chemo,est_mul_chemo_bu_r)
number_mul_chemo_full<-rbind(number_mul_chemo,est_mul_chemo_bu_r)
nonzero_mul_chemo_full<-rbind(nonzero_mul_chemo,est_mul_chemo_bu_r)
users_mul_chemo_full<-rbind(users_mul_chemo,est_mul_chemo_bu_r)
nonusers_mul_chemo_full<-rbind(nonusers_mul_chemo,est_mul_chemo_bu_r)



estimate_control_full<-arrange(estimate_control_full,rownames(estimate_control_full))
estimate_control_full<-dplyr::select(estimate_control_full,Prednisolone:Idarubicin,everything())
estimate_mul_control_full<-arrange(estimate_mul_control_full,rownames(estimate_mul_control_full))
estimate_mul_control_full<-dplyr::select(estimate_mul_control_full,Prednisolone:Idarubicin,everything())

std_control_full<-arrange(std_control_full,rownames(std_control_full))
std_control_full<-select(std_control_full,Prednisolone:Idarubicin,everything())
std_mul_control_full<-arrange(std_mul_control_full,rownames(std_mul_control_full))
std_mul_control_full<-select(std_mul_control_full,Prednisolone:Idarubicin,everything())

pva_control_full<-arrange(pva_control_full,rownames(pva_control_full))
pva_control_full<-select(pva_control_full,Prednisolone:Idarubicin,everything())
pva_mul_control_full<-arrange(pva_mul_control_full,rownames(pva_mul_control_full))
pva_mul_control_full<-select(pva_mul_control_full,Prednisolone:Idarubicin,everything())

tva_control_full<-arrange(tva_control_full,rownames(tva_control_full))
tva_control_full<-select(tva_control_full,Prednisolone:Idarubicin,everything())
tva_mul_control_full<-arrange(tva_mul_control_full,rownames(tva_mul_control_full))
tva_mul_control_full<-select(tva_mul_control_full,Prednisolone:Idarubicin,everything())

number_hsct_full<-arrange(number_hsct_full,rownames(number_hsct_full))
number_hsct_full<-select(number_hsct_full,Prednisolone:Idarubicin,everything())
number_mul_hsct_full<-arrange(number_mul_hsct_full,rownames(number_mul_hsct_full))
number_mul_hsct_full<-select(number_mul_hsct_full,Prednisolone:Idarubicin,everything())

nonzero_hsct_full<-arrange(nonzero_hsct_full,rownames(nonzero_hsct_full))
nonzero_hsct_full<-select(nonzero_hsct_full,Prednisolone:Idarubicin,everything())
nonzero_mul_hsct_full<-arrange(nonzero_mul_hsct_full,rownames(nonzero_mul_hsct_full))
nonzero_mul_hsct_full<-select(nonzero_mul_hsct_full,Prednisolone:Idarubicin,everything())

users_hsct_full<-arrange(users_hsct_full,rownames(users_hsct_full))
users_hsct_full<-select(users_hsct_full,Prednisolone:Idarubicin,everything())
users_mul_hsct_full<-arrange(users_mul_hsct_full,rownames(users_mul_hsct_full))
users_mul_hsct_full<-select(users_mul_hsct_full,Prednisolone:Idarubicin,everything())

nonusers_hsct_full<-arrange(nonusers_hsct_full,rownames(nonusers_hsct_full))
nonusers_hsct_full<-select(nonusers_hsct_full,Prednisolone:Idarubicin,everything())
nonusers_mul_hsct_full<-arrange(nonusers_mul_hsct_full,rownames(nonusers_mul_hsct_full))
nonusers_mul_hsct_full<-select(nonusers_mul_hsct_full,Prednisolone:Idarubicin,everything())


estimate_gvhd_full<-arrange(estimate_gvhd_full,rownames(estimate_gvhd_full))
estimate_gvhd_full<-select(estimate_gvhd_full,Prednisolone:Idarubicin,everything())
estimate_mul_gvhd_full<-arrange(estimate_mul_gvhd_full,rownames(estimate_mul_gvhd_full))
estimate_mul_gvhd_full<-select(estimate_mul_gvhd_full,Prednisolone:Idarubicin,everything())

std_gvhd_full<-arrange(std_gvhd_full,rownames(std_gvhd_full))
std_gvhd_full<-select(std_gvhd_full,Prednisolone:Idarubicin,everything())
std_mul_gvhd_full<-arrange(std_mul_gvhd_full,rownames(std_mul_gvhd_full))
std_mul_gvhd_full<-select(std_mul_gvhd_full,Prednisolone:Idarubicin,everything())

pva_gvhd_full<-arrange(pva_gvhd_full,rownames(pva_gvhd_full))
pva_gvhd_full<-select(pva_gvhd_full,Prednisolone:Idarubicin,everything())
pva_mul_gvhd_full<-arrange(pva_mul_gvhd_full,rownames(pva_mul_gvhd_full))
pva_mul_gvhd_full<-select(pva_mul_gvhd_full,Prednisolone:Idarubicin,everything())

tva_gvhd_full<-arrange(tva_gvhd_full,rownames(tva_gvhd_full))
tva_gvhd_full<-select(tva_gvhd_full,Prednisolone:Idarubicin,everything())
tva_mul_gvhd_full<-arrange(tva_mul_gvhd_full,rownames(tva_mul_gvhd_full))
tva_mul_gvhd_full<-select(tva_mul_gvhd_full,Prednisolone:Idarubicin,everything())

number_gvhd_full<-arrange(number_gvhd_full,rownames(number_gvhd_full))
number_gvhd_full<-select(number_gvhd_full,Prednisolone:Idarubicin,everything())
number_mul_gvhd_full<-arrange(number_mul_gvhd_full,rownames(number_mul_gvhd_full))
number_mul_gvhd_full<-select(number_mul_gvhd_full,Prednisolone:Idarubicin,everything())

nonzero_gvhd_full<-arrange(nonzero_gvhd_full,rownames(nonzero_gvhd_full))
nonzero_gvhd_full<-select(nonzero_gvhd_full,Prednisolone:Idarubicin,everything())
nonzero_mul_gvhd_full<-arrange(nonzero_mul_gvhd_full,rownames(nonzero_mul_gvhd_full))
nonzero_mul_gvhd_full<-select(nonzero_mul_gvhd_full,Prednisolone:Idarubicin,everything())

users_gvhd_full<-arrange(users_gvhd_full,rownames(users_gvhd_full))
users_gvhd_full<-select(users_gvhd_full,Prednisolone:Idarubicin,everything())
users_mul_gvhd_full<-arrange(users_mul_gvhd_full,rownames(users_mul_gvhd_full))
users_mul_gvhd_full<-select(users_mul_gvhd_full,Prednisolone:Idarubicin,everything())

nonusers_gvhd_full<-arrange(nonusers_gvhd_full,rownames(nonusers_gvhd_full))
nonusers_gvhd_full<-select(nonusers_gvhd_full,Prednisolone:Idarubicin,everything())
nonusers_mul_gvhd_full<-arrange(nonusers_mul_gvhd_full,rownames(nonusers_mul_gvhd_full))
nonusers_mul_gvhd_full<-select(nonusers_mul_gvhd_full,Prednisolone:Idarubicin,everything())


estimate_chemo_full<-arrange(estimate_chemo_full,rownames(estimate_chemo_full))
estimate_mul_chemo_full<-arrange(estimate_mul_chemo_full,rownames(estimate_mul_chemo_full))

std_chemo_full<-arrange(std_chemo_full,rownames(std_chemo_full))
std_mul_chemo_full<-arrange(std_mul_chemo_full,rownames(std_mul_chemo_full))

pva_chemo_full<-arrange(pva_chemo_full,rownames(pva_chemo_full))
pva_mul_chemo_full<-arrange(pva_mul_chemo_full,rownames(pva_mul_chemo_full))

tva_chemo_full<-arrange(tva_chemo_full,rownames(tva_chemo_full))
tva_mul_chemo_full<-arrange(tva_mul_chemo_full,rownames(tva_mul_chemo_full))

number_chemo_full<-arrange(number_chemo_full,rownames(number_chemo_full))
number_mul_chemo_full<-arrange(number_mul_chemo_full,rownames(number_mul_chemo_full))

nonzero_chemo_full<-arrange(nonzero_chemo_full,rownames(nonzero_chemo_full))
nonzero_mul_chemo_full<-arrange(nonzero_mul_chemo_full,rownames(nonzero_mul_chemo_full))

users_chemo_full<-arrange(users_chemo_full,rownames(users_chemo_full))
users_mul_chemo_full<-arrange(users_mul_chemo_full,rownames(users_mul_chemo_full))

nonusers_chemo_full<-arrange(nonusers_chemo_full,rownames(nonusers_chemo_full))
nonusers_mul_chemo_full<-arrange(nonusers_mul_chemo_full,rownames(nonusers_mul_chemo_full))

```

## 5.1 Individual drug-taxonomies association meta-analysis

### Functions

```{r}
taxa_rownames<-rownames(estimate_control_full)
colnames<-c("N_hsct","Nonzeros_hsct","Users_hsct","Nonusers_hsct",
            "est_hsct","std_hsct","pval_hsct","tval_hsct",
            "N_gvhd","Nonzeros_gvhd","Users_gvhd","Nonusers_gvhd",
            "est_gvhd","std_gvhd","pval_gvhd","tval_gvhd",
            "N_chemo","Nonzeros_chemo","Users_chemo","Nonusers_chemo",
            "est_chemo","std_chemo","pval_chemo","tval_chemo")
metaanalysis_data<-array(0,dim=c(127,24,48),dimnames=list(taxa_rownames,colnames,drug_all))
for (i in 1:48) {
  metaanalysis_data[,1,i]<-number_hsct_full[,i]
  metaanalysis_data[,2,i]<-nonzero_hsct_full[,i]
  metaanalysis_data[,3,i]<-users_hsct_full[,i]
  metaanalysis_data[,4,i]<-nonusers_hsct_full[,i]
  metaanalysis_data[,5,i]<-estimate_control_full[,i]
  metaanalysis_data[,6,i]<-std_control_full[,i]
  metaanalysis_data[,7,i]<-pva_control_full[,i]
  metaanalysis_data[,8,i]<-tva_control_full[,i]
  metaanalysis_data[,9,i]<-number_gvhd_full[,i]
  metaanalysis_data[,10,i]<-nonzero_gvhd_full[,i]
  metaanalysis_data[,11,i]<-users_gvhd_full[,i]
  metaanalysis_data[,12,i]<-nonusers_gvhd_full[,i]
  metaanalysis_data[,13,i]<-estimate_gvhd_full[,i]
  metaanalysis_data[,14,i]<-std_gvhd_full[,i]
  metaanalysis_data[,15,i]<-pva_gvhd_full[,i]
  metaanalysis_data[,16,i]<-tva_gvhd_full[,i]
  metaanalysis_data[,17,i]<-number_chemo_full[,i]
  metaanalysis_data[,18,i]<-nonzero_chemo_full[,i]
  metaanalysis_data[,19,i]<-users_chemo_full[,i]
  metaanalysis_data[,20,i]<-nonusers_chemo_full[,i]
  metaanalysis_data[,21,i]<-estimate_chemo_full[,i]
  metaanalysis_data[,22,i]<-std_chemo_full[,i]
  metaanalysis_data[,23,i]<-pva_chemo_full[,i]
  metaanalysis_data[,24,i]<-tva_chemo_full[,i]
  
}

metaanalysis_data_1<-adply(metaanalysis_data,c(1,3))
colnames(metaanalysis_data_1)[2]<-"drug"
colnames(metaanalysis_data_1)[1]<-"taxa"

taxa_rownames_mul<-rownames(estimate_mul_control_full)
colnames_mul<-c("N_hsct","Nonzeros_hsct","Users_hsct","Nonusers_hsct",
                "est_mul_hsct","std_mul_hsct","pval_mul_hsct","tval_mul_hsct",
                "N_gvhd","Nonzeros_gvhd","Users_gvhd","Nonusers_gvhd",
                "est_mul_gvhd","std_mul_gvhd","pval_mul_gvhd","tval_mul_gvhd",
                "N_chemo","Nonzeros_chemo","Users_chemo","Nonusers_chemo",
                "est_mul_chemo","std_mul_chemo","pval_mul_chemo","tval_mul_chemo")
metaanalysis_data_mul<-array(0,dim=c(127,24,47),dimnames=list(taxa_rownames_mul,colnames_mul,drug_mul_all))
for (i in 1:47) {
  metaanalysis_data_mul[,1,i]<-number_mul_hsct_full[,i]
  metaanalysis_data_mul[,2,i]<-nonzero_mul_hsct_full[,i]
  metaanalysis_data_mul[,3,i]<-users_mul_hsct_full[,i]
  metaanalysis_data_mul[,4,i]<-nonusers_mul_hsct_full[,i]
  metaanalysis_data_mul[,5,i]<-estimate_mul_control_full[,i]
  metaanalysis_data_mul[,6,i]<-std_mul_control_full[,i]
  metaanalysis_data_mul[,7,i]<-pva_mul_control_full[,i]
  metaanalysis_data_mul[,8,i]<-tva_mul_control_full[,i]
  metaanalysis_data_mul[,9,i]<-number_mul_gvhd_full[,i]
  metaanalysis_data_mul[,10,i]<-nonzero_mul_gvhd_full[,i]
  metaanalysis_data_mul[,11,i]<-users_mul_gvhd_full[,i]
  metaanalysis_data_mul[,12,i]<-nonusers_mul_gvhd_full[,i]
  metaanalysis_data_mul[,13,i]<-estimate_mul_gvhd_full[,i]
  metaanalysis_data_mul[,14,i]<-std_mul_gvhd_full[,i]
  metaanalysis_data_mul[,15,i]<-pva_mul_gvhd_full[,i]
  metaanalysis_data_mul[,16,i]<-tva_mul_gvhd_full[,i]
  metaanalysis_data_mul[,17,i]<-number_mul_chemo_full[,i]
  metaanalysis_data_mul[,18,i]<-nonzero_mul_chemo_full[,i]
  metaanalysis_data_mul[,19,i]<-users_mul_chemo_full[,i]
  metaanalysis_data_mul[,20,i]<-nonusers_mul_chemo_full[,i]
  metaanalysis_data_mul[,21,i]<-estimate_mul_chemo_full[,i]
  metaanalysis_data_mul[,22,i]<-std_mul_chemo_full[,i]
  metaanalysis_data_mul[,23,i]<-pva_mul_chemo_full[,i]
  metaanalysis_data_mul[,24,i]<-tva_mul_chemo_full[,i]
  
}

metaanalysis_data_mul1<-adply(metaanalysis_data_mul,c(1,3))
colnames(metaanalysis_data_mul1)[2]<-"drug"
colnames(metaanalysis_data_mul1)[1]<-"taxa"


drugmeta_hsctgvhd<-function(a){
  a<-a[,1:20]
  # a<-a[complete.cases(a$N_hsct),]
  a$inverse_var.hsct<-1/a$std_hsct^2
  a$inverse_var.gvhd<-1/a$std_gvhd^2
  #Calculate SE  #53
  a$se=sqrt((a$inverse_var.hsct+a$inverse_var.gvhd))
  
  #Calculate Beta #54
  a$beta=(a$inverse_var.hsct*a$est_hsct+a$inverse_var.gvhd*a$est_gvhd)/(a$inverse_var.hsct+a$inverse_var.gvhd)
  #	selection$beta=(selection$inverse_var.ibd*selection$coeff.correcting.all.IBD+selection$inverse_var.lld*selection$coeff.correcting.all.LLD)/(selection$inverse_var.ibd+selection$inverse_var.lld)
  
  #Calculate Z-score #55
  a$Z=a$beta/a$se
  
  #Calculate meta p-value #56
  a$P=2*pnorm(-abs(a$Z))
  #Remove data only avalable in one cohort
  # a=a[complete.cases(a$P), ]
  
  #Adjust pvalue with FDR #57
  a$FDR=p.adjust(a$P,method = "fdr")
  
  #Create empty columns
  a$Het.Q="NS" 
  a$Het.I2="NS" 
  a$Het.Pval="NS" 
  
  #Heterogeneity using Cochran's Q-test for meta-FDR < 0.1
  for (i in 1:length(rownames(a))){
    if ((a$FDR[i]<0.1)&(!is.na(a$FDR[i]))){
      #Select coefficients
      TE=c( a[i,5], a[i,13])
      #Select Standart error
      SE=c( a[i,6], a[i,14])
      het=metagen(TE,SE)
      #Change the number here depending of the number of columns in your data, should match with column Het.I2
      a[i,29]=het$I2
      #Match Het.Q column number
      a[i,28]=het$Q 
      #Calculate p-value from Q calculation // Het.Pval
      a[i,30]=pchisq(het$Q,df=2,lower.tail=F)
    }else{
      #Change the number here depending of the number of columns in your data, should match with column Het.I2
      a[i,29]=NA
      #Match Het.Q column number
      a[i,28]=NA 
      #Calculate p-value from Q calculation // Het.Pval
      a[i,30]=NA
      
    }
  }
  
  a <- a %>%
    mutate(Taxa=rownames(a),
           Drug="a",
           est_chemo="NS",
           std_chemo="NS",
           pval_chemo="NS",
           tval_chemo="NS",
           inverse_var.chemo="NS") %>%
    select(Taxa,Drug,
           N_hsct,Nonzeros_hsct,Users_hsct,Nonusers_hsct,
           est_hsct,std_hsct,pval_hsct,tval_hsct,
           N_gvhd,Nonzeros_gvhd,Users_gvhd,Nonusers_gvhd,
           est_gvhd,std_gvhd,pval_gvhd,tval_gvhd,
           N_chemo,Nonzeros_chemo,Users_chemo,Nonusers_chemo,
           est_chemo,std_chemo,pval_chemo,tval_chemo,
           inverse_var.hsct,inverse_var.gvhd,inverse_var.chemo,
           se,beta,Z,P,FDR,everything())
  return(a)
}
drugmeta_all<-function(a){
  a$inverse_var.hsct<-1/a$std_hsct^2
  a$inverse_var.gvhd<-1/a$std_gvhd^2
  a$inverse_var.chemo<-1/a$std_chemo^2
  #Calculate SE  #53
  a$se=sqrt(1/(a$inverse_var.hsct+a$inverse_var.gvhd+a$inverse_var.chemo))
  
  #Calculate Beta #54
  a$beta=(a$inverse_var.hsct*a$est_hsct+a$inverse_var.gvhd*a$est_gvhd+a$inverse_var.chemo*a$est_chemo)/(a$inverse_var.hsct+a$inverse_var.gvhd+a$inverse_var.chemo)
  #	selection$beta=(selection$inverse_var.ibd*selection$coeff.correcting.all.IBD+selection$inverse_var.lld*selection$coeff.correcting.all.LLD)/(selection$inverse_var.ibd+selection$inverse_var.lld)
  
  #Calculate Z-score #55
  a$Z=a$beta/a$se
  
  #Calculate meta p-value #56
  a$P=2*pnorm(-abs(a$Z))
  #Remove data only avalable in one cohort
  # a=a[complete.cases(a$P), ]
  
  #Adjust pvalue with FDR #57
  a$FDR=p.adjust(a$P,method = "fdr")
  
  #Create empty columns
  a$Het.Q="NS" #58
  a$Het.I2="NS" #59
  a$Het.Pval="NS" #60
  
  #Heterogeneity using Cochran's Q-test for meta-FDR < 0.1
  for (i in 1:length(rownames(a))){
    if ((a$FDR[i]<0.1)&(!is.na(a$FDR[i]))){
      #Select coefficients
      TE=c( a[i,5], a[i,13], a[i,21])
      #Select Standart error
      SE=c( a[i,6], a[i,14], a[i,22])
      het=metagen(TE,SE)
      #Change the number here depending of the number of columns in your data, should match with column Het.I2
      a[i,34]=het$I2
      #Match Het.Q column number
      a[i,33]=het$Q 
      #Calculate p-value from Q calculation // Het.Pval
      a[i,35]=pchisq(het$Q,df=2,lower.tail=F)
    }else{
      a[i,34]=NA
      #Match Het.Q column number
      a[i,33]=NA 
      #Calculate p-value from Q calculation // Het.Pval
      a[i,35]=NA
      
    }
  }
  
  a <- a %>%
    mutate(Taxa=rownames(a),Drug="a") %>%
    select(Taxa,Drug,
           N_hsct,Nonzeros_hsct,Users_hsct,Nonusers_hsct,
           est_hsct,std_hsct,pval_hsct,tval_hsct,
           N_gvhd,Nonzeros_gvhd,Users_gvhd,Nonusers_gvhd,
           est_gvhd,std_gvhd,pval_gvhd,tval_gvhd,
           N_chemo,Nonzeros_chemo,Users_chemo,Nonusers_chemo,
           est_chemo,std_chemo,pval_chemo,tval_chemo,
           inverse_var.hsct,inverse_var.gvhd,inverse_var.chemo,
           se,beta,Z,P,FDR,everything())
  
}
drugmeta_mul_hsctgvhd<-function(a){
  a<-a[,1:20]
  # a<-a[complete.cases(a$N_hsct),]
  a$inverse_mul_var.hsct<-1/a$std_mul_hsct^2
  a$inverse_mul_var.gvhd<-1/a$std_mul_gvhd^2
  #Calculate SE  #53
  a$se_mul=sqrt((a$inverse_mul_var.hsct+a$inverse_mul_var.gvhd))
  
  #Calculate Beta #54
  a$beta_mul=(a$inverse_mul_var.hsct*a$est_mul_hsct+a$inverse_mul_var.gvhd*a$est_mul_gvhd)/(a$inverse_mul_var.hsct+a$inverse_mul_var.gvhd)
  #	selection$beta=(selection$inverse_var.ibd*selection$coeff.correcting.all.IBD+selection$inverse_var.lld*selection$coeff.correcting.all.LLD)/(selection$inverse_var.ibd+selection$inverse_var.lld)
  
  #Calculate Z-score #55
  a$Z_mul=a$beta_mul/a$se_mul
  
  #Calculate meta p-value #56
  a$P_mul=2*pnorm(-abs(a$Z_mul))
  #Remove data only avalable in one cohort
  # a=a[complete.cases(a$P), ]
  
  #Adjust pvalue with FDR #57
  a$FDR_mul=p.adjust(a$P_mul,method = "fdr")
  
  #Create empty columns
  a$Het.Q_mul="NS" #58
  a$Het.I2_mul="NS" #59
  a$Het.Pval_mul="NS" #60
  
  #Heterogeneity using Cochran's Q-test for meta-FDR < 0.1
  for (i in 1:length(rownames(a))){
    if ((a$FDR_mul[i]<0.1)&(!is.na(a$FDR_mul[i]))){
      #Select coefficients
      TE=c( a[i,5], a[i,13])
      #Select Standart error
      SE=c( a[i,6], a[i,14])
      het=metagen(TE,SE)
      #Change the number here depending of the number of columns in your data, should match with column Het.I2
      a[i,29]=het$I2
      #Match Het.Q column number
      a[i,28]=het$Q 
      #Calculate p-value from Q calculation // Het.Pval
      a[i,30]=pchisq(het$Q,df=2,lower.tail=F)
    }else{
      #Change the number here depending of the number of columns in your data, should match with column Het.I2
      a[i,29]=NA
      #Match Het.Q column number
      a[i,28]=NA 
      #Calculate p-value from Q calculation // Het.Pval
      a[i,30]=NA
      
    }
  }
  
  a <- a %>%
    mutate(Taxa=rownames(a),
           Drug="a",
           est_mul_chemo="NS",
           std_mul_chemo="NS",
           pval_mul_chemo="NS",
           tval_mul_chemo="NS",
           inverse_mul_var.chemo="NS") %>%
    select(Taxa,Drug,
           N_hsct,Nonzeros_hsct,Users_hsct,Nonusers_hsct,
           est_mul_hsct,std_mul_hsct,pval_mul_hsct,tval_mul_hsct,
           N_gvhd,Nonzeros_gvhd,Users_gvhd,Nonusers_gvhd,
           est_mul_gvhd,std_mul_gvhd,pval_mul_gvhd,tval_mul_gvhd,
           N_chemo,Nonzeros_chemo,Users_chemo,Nonusers_chemo,
           est_mul_chemo,std_mul_chemo,pval_mul_chemo,tval_mul_chemo,
           inverse_mul_var.hsct,inverse_mul_var.gvhd,inverse_mul_var.chemo,
           se_mul,beta_mul,Z_mul,P_mul,FDR_mul,everything())
  return(a)
}
drug_merge<-function(a,b){
  mergetable<-merge(a,b[,c(1,7:10,15:18,23:37)],by="Taxa") |> select(
    Taxa,Drug,
    N_hsct,Nonzeros_hsct,Users_hsct,Nonusers_hsct,
    est_hsct,std_hsct,pval_hsct,tval_hsct,
    est_mul_hsct,std_mul_hsct,pval_mul_hsct,tval_mul_hsct,
    N_gvhd,Nonzeros_gvhd,Users_gvhd,Nonusers_gvhd,
    est_gvhd,std_gvhd,pval_gvhd,tval_gvhd,
    est_mul_gvhd,std_mul_gvhd,pval_mul_gvhd,tval_mul_gvhd,
    N_chemo,Nonzeros_chemo,Users_chemo,Nonusers_chemo,
    est_chemo,std_chemo,pval_chemo,tval_chemo,
    est_mul_chemo,std_mul_chemo,pval_mul_chemo,tval_mul_chemo,
    inverse_var.hsct,inverse_var.gvhd,inverse_var.chemo,
    inverse_mul_var.hsct,inverse_mul_var.gvhd,inverse_mul_var.chemo,
    se,beta,Z,P,FDR,Het.Q,Het.I2,Het.Pval,
    se_mul,beta_mul,Z_mul,P_mul,FDR_mul,Het.Q_mul,Het.I2_mul,Het.Pval_mul
  )
  return(mergetable)
}

```

### busulfan

```{r}
busulfan<-as.data.frame(metaanalysis_data[,,"Busulfan"])
busulfan<-drugmeta_hsctgvhd(busulfan)
busulfan$Drug<-"Busulfan"
test_busulfan<-busulfan[complete.cases(busulfan$Het.Pval),] 
test_busulfan_n<-sum((as.numeric(test_busulfan$FDR) < 0.1)&(as.numeric(test_busulfan$Het.Pval)>0.05))
hp_busulfan<- test_busulfan |> filter((as.numeric(test_busulfan$FDR) < 0.1)&(as.numeric(test_busulfan$Het.Pval)>0.05)) |> 
  select(Taxa,Drug,beta,FDR) |> mutate(`Influence`=ifelse(beta<0,"Inhibition","Promotion"),anti="no")

busulfan_mul<-as.data.frame(metaanalysis_data_mul[,,"Busulfan"])
busulfan_mul<-drugmeta_mul_hsctgvhd(busulfan_mul)
busulfan_mul$Drug<-"Busulfan"
test_busulfan_mul<-busulfan_mul[complete.cases(busulfan_mul$Het.Pval_mul),] 
test_busulfan_mul_n<-sum((as.numeric(test_busulfan_mul$FDR_mul) < 0.1)&(as.numeric(test_busulfan_mul$Het.Pval_mul)>0.05))
hp_mul_busulfan<- test_busulfan_mul |> filter((as.numeric(test_busulfan_mul$FDR_mul) < 0.1)&(as.numeric(test_busulfan_mul$Het.Pval_mul)>0.05))|> 
  select(Taxa,Drug,beta_mul,FDR_mul) |> mutate(`Influence`=ifelse(beta_mul<0,"Inhibition","Promotion"),anti="no")

busulfan_merge<-drug_merge(busulfan,busulfan_mul)
write.csv(busulfan_merge,file = "drug/busulfan.csv")

```

### cyclophosphamide

```{r}
cyclophosphamide<-as.data.frame(metaanalysis_data[,,"Cyclophosphamide"])
cyclophosphamide<-drugmeta_hsctgvhd(cyclophosphamide)
cyclophosphamide$Drug<-"Cyclophosphamide"
test_cyclophosphamide<-  cyclophosphamide[complete.cases(cyclophosphamide$Het.Pval),] 
test_cyclophosphamide_n<-sum((as.numeric(test_cyclophosphamide$FDR) < 0.1)&(as.numeric(test_cyclophosphamide$Het.Pval)>0.05))
hp_cyclophosphamide<- test_cyclophosphamide |> filter((as.numeric(test_cyclophosphamide$FDR) < 0.1)&(as.numeric(test_cyclophosphamide$Het.Pval)>0.05))|> 
  select(Taxa,Drug,beta,FDR) |> mutate(`Influence`=ifelse(beta<0,"Inhibition","Promotion"),anti="no")

cyclophosphamide_mul<-as.data.frame(metaanalysis_data_mul[,,"Cyclophosphamide"])
cyclophosphamide_mul<-drugmeta_mul_hsctgvhd(cyclophosphamide_mul)
cyclophosphamide_mul$Drug<-"Cyclophosphamide"
test_cyclophosphamide_mul<-  cyclophosphamide_mul[complete.cases(cyclophosphamide_mul$Het.Pval_mul),] 
test_cyclophosphamide_mul_n<-sum((as.numeric(test_cyclophosphamide_mul$FDR_mul) < 0.1)&(as.numeric(test_cyclophosphamide_mul$Het.Pval_mul)>0.05))
hp_mul_cyclophosphamide<- test_cyclophosphamide_mul |> filter((as.numeric(test_cyclophosphamide_mul$FDR_mul) < 0.1)&(as.numeric(test_cyclophosphamide_mul$Het.Pval_mul)>0.05))|> 
  select(Taxa,Drug,beta_mul,FDR_mul) |> mutate(`Influence`=ifelse(beta_mul<0,"Inhibition","Promotion"),anti="no")

cyclophosphamide_merge<-drug_merge(cyclophosphamide,cyclophosphamide_mul)
write.csv(cyclophosphamide_merge,file = "drug/cyclophosphamide.csv")


```

### fludarabine

```{r}
fludarabine<-as.data.frame(metaanalysis_data[,,"Fludarabine"])
fludarabine<-drugmeta_hsctgvhd(fludarabine)
fludarabine$Drug<-"Fludarabine"
test_fludarabine<-  fludarabine[complete.cases(fludarabine$Het.Pval),] 
test_fludarabine_n<-sum((as.numeric(test_fludarabine$FDR) < 0.1)&(as.numeric(test_fludarabine$Het.Pval)>0.05))
hp_fludarabine<- test_fludarabine |> filter((as.numeric(test_fludarabine$FDR) < 0.1)&(as.numeric(test_fludarabine$Het.Pval)>0.05))|> 
  select(Taxa,Drug,beta,FDR) |> mutate(`Influence`=ifelse(beta<0,"Inhibition","Promotion"),anti="no")

fludarabine_mul<-as.data.frame(metaanalysis_data_mul[,,"Fludarabine"])
fludarabine_mul<-drugmeta_mul_hsctgvhd(fludarabine_mul)
fludarabine_mul$Drug<-"Fludarabine"
test_fludarabine_mul<-  fludarabine_mul[complete.cases(fludarabine_mul$Het.Pval_mul),] 
test_fludarabine_mul_n<-sum((as.numeric(test_fludarabine_mul$FDR_mul) < 0.1)&(as.numeric(test_fludarabine_mul$Het.Pval_mul)>0.05))
hp_mul_fludarabine<- test_fludarabine_mul |> filter((as.numeric(test_fludarabine_mul$FDR_mul) < 0.1)&(as.numeric(test_fludarabine_mul$Het.Pval_mul)>0.05))|> 
  select(Taxa,Drug,beta_mul,FDR_mul) |> mutate(`Influence`=ifelse(beta_mul<0,"Inhibition","Promotion"),anti="no")

fludarabine_merge<-drug_merge(fludarabine,fludarabine_mul)
write.csv(fludarabine_merge,file = "drug/fludarabine.csv")

```

### melphalan

```{r}
melphalan<-as.data.frame(metaanalysis_data[,,"Melphalan"])
melphalan<-drugmeta_hsctgvhd(melphalan)
melphalan$Drug<-"Melphalan"
test_melphalan <-  melphalan[complete.cases(melphalan$Het.Pval),] 
test_melphalan_n<-sum((as.numeric(test_melphalan$FDR) < 0.1)&(as.numeric(test_melphalan$Het.Pval)>0.05))
hp_melphalan<- test_melphalan |> filter((as.numeric(test_melphalan$FDR) < 0.1)&(as.numeric(test_melphalan$Het.Pval)>0.05))|> 
  select(Taxa,Drug,beta,FDR) |> mutate(`Influence`=ifelse(beta<0,"Inhibition","Promotion"),anti="no")

write.csv(melphalan,file = "drug/melphalan.csv")

```

### VP16

```{r}
VP16<-as.data.frame(metaanalysis_data[,,"VP16"])
VP16<-drugmeta_hsctgvhd(VP16)
VP16$Drug<-"VP16"
test_VP16<-  VP16[complete.cases(VP16$Het.Pval),] 
test_VP16_n<-sum((as.numeric(test_VP16$FDR) < 0.1)&(as.numeric(test_VP16$Het.Pval)>0.05))
hp_VP16<- test_VP16 |> filter((as.numeric(test_VP16$FDR) < 0.1)&(as.numeric(test_VP16$Het.Pval)>0.05))|> 
  select(Taxa,Drug,beta,FDR) |> mutate(`Influence`=ifelse(beta<0,"Inhibition","Promotion"),anti="no")

VP16_mul<-as.data.frame(metaanalysis_data_mul[,,"VP16"])
VP16_mul<-drugmeta_mul_hsctgvhd(VP16_mul)
VP16_mul$Drug<-"VP16"
test_VP16_mul<-  VP16_mul[complete.cases(VP16_mul$Het.Pval_mul),] 
test_VP16_mul_n<-sum((as.numeric(test_VP16_mul$FDR_mul) < 0.1)&(as.numeric(test_VP16_mul$Het.Pval_mul)>0.05))
hp_mul_VP16<- test_VP16_mul |> filter((as.numeric(test_VP16_mul$FDR_mul) < 0.1)&(as.numeric(test_VP16_mul$Het.Pval_mul)>0.05))|> 
  select(Taxa,Drug,beta_mul,FDR_mul) |> mutate(`Influence`=ifelse(beta_mul<0,"Inhibition","Promotion"),anti="no")

VP16_merge<-drug_merge(VP16,VP16_mul)
write.csv(VP16_merge,file = "drug/VP16.csv")

```

### TBI(not included in the final version of manuscript)

```{r}
TBI<-as.data.frame(metaanalysis_data[,,"Total.Body.Irradiation"])
TBI<-drugmeta_hsctgvhd(TBI)
TBI$Drug<-"TBI"
test_TBI<-  TBI[complete.cases(TBI$Het.Pval),] 
test_TBI_n<-sum((as.numeric(test_TBI$FDR) < 0.1)&(as.numeric(test_TBI$Het.Pval)>0.05))
hp_TBI<- test_TBI |> filter((as.numeric(test_TBI$FDR) < 0.1)&(as.numeric(test_TBI$Het.Pval)>0.05))|> 
  select(Taxa,Drug,beta,FDR) |> mutate(`Influence`=ifelse(beta<0,"Inhibition","Promotion"),anti="no")

TBI_mul<-as.data.frame(metaanalysis_data_mul[,,"Total.Body.Irradiation"])
TBI_mul<-drugmeta_mul_hsctgvhd(TBI_mul)
TBI_mul$Drug<-"TBI"
test_TBI_mul<-  TBI_mul[complete.cases(TBI_mul$Het.Pval_mul),] 
test_TBI_mul_n<-sum((as.numeric(test_TBI_mul$FDR_mul) < 0.1)&(as.numeric(test_TBI_mul$Het.Pval_mul)>0.05))
hp_mul_TBI<- test_TBI_mul |> filter((as.numeric(test_TBI_mul$FDR_mul) < 0.1)&(as.numeric(test_TBI_mul$Het.Pval_mul)>0.05))|> 
  select(Taxa,Drug,beta_mul,FDR_mul) |> mutate(`Influence`=ifelse(beta_mul<0,"Inhibition","Promotion"),anti="no")

TBI_merge<-drug_merge(TBI,TBI_mul)
write.csv(TBI_merge,file = "drug/TBI.csv")

```

### Thiotepa

```{r}
Thiotepa<-as.data.frame(metaanalysis_data[,,"Thiotepa"])
Thiotepa<-drugmeta_hsctgvhd(Thiotepa)
Thiotepa$Drug<-"Thiotepa"
test_Thiotepa<-  Thiotepa[complete.cases(Thiotepa$Het.Pval),] 
test_Thiotepa_n<-sum((as.numeric(test_Thiotepa$FDR) < 0.1)&(as.numeric(test_Thiotepa$Het.Pval)>0.05))
hp_Thiotepa<- test_Thiotepa |> filter((as.numeric(test_Thiotepa$FDR) < 0.1)&(as.numeric(test_Thiotepa$Het.Pval)>0.05))|> 
  select(Taxa,Drug,beta,FDR) |> mutate(`Influence`=ifelse(beta<0,"Inhibition","Promotion"),anti="no")

Thiotepa_mul<-as.data.frame(metaanalysis_data_mul[,,"Thiotepa"])
Thiotepa_mul<-drugmeta_mul_hsctgvhd(Thiotepa_mul)
Thiotepa_mul$Drug<-"Thiotepa"
test_Thiotepa_mul<-  Thiotepa_mul[complete.cases(Thiotepa_mul$Het.Pval_mul),] 
test_Thiotepa_mul_n<-sum((as.numeric(test_Thiotepa_mul$FDR_mul) < 0.1)&(as.numeric(test_Thiotepa_mul$Het.Pval_mul)>0.05))
hp_mul_Thiotepa<- test_Thiotepa_mul |> filter((as.numeric(test_Thiotepa_mul$FDR_mul) < 0.1)&(as.numeric(test_Thiotepa_mul$Het.Pval_mul)>0.05))|> 
  select(Taxa,Drug,beta_mul,FDR_mul) |> mutate(`Influence`=ifelse(beta_mul<0,"Inhibition","Promotion"),anti="no")

Thiotepa_merge<-drug_merge(Thiotepa,Thiotepa_mul)
write.csv(Thiotepa_merge,file = "drug/Thiotepa.csv")

```

### Methotrexate

```{r}
Methotrexate<-as.data.frame(metaanalysis_data[,,"Methotrexate"])
Methotrexate<-drugmeta_hsctgvhd(Methotrexate)
Methotrexate$Drug<-"Methotrexate"
test_Methotrexate<-  Methotrexate[complete.cases(Methotrexate$Het.Pval),] 
test_Methotrexate_n<-sum((as.numeric(test_Methotrexate$FDR) < 0.1)&(as.numeric(test_Methotrexate$Het.Pval)>0.05))
hp_Methotrexate<- test_Methotrexate |> filter((as.numeric(test_Methotrexate$FDR) < 0.1)&(as.numeric(test_Methotrexate$Het.Pval)>0.05))|> 
  select(Taxa,Drug,beta,FDR) |> mutate(`Influence`=ifelse(beta<0,"Inhibition","Promotion"),anti="no")

Methotrexate_mul<-as.data.frame(metaanalysis_data_mul[,,"Methotrexate"])
Methotrexate_mul<-drugmeta_mul_hsctgvhd(Methotrexate_mul)
Methotrexate_mul$Drug<-"Methotrexate"
test_Methotrexate_mul<-  Methotrexate_mul[complete.cases(Methotrexate_mul$Het.Pval_mul),] 
test_Methotrexate_mul_n<-sum((as.numeric(test_Methotrexate_mul$FDR_mul) < 0.1)&(as.numeric(test_Methotrexate_mul$Het.Pval_mul)>0.05))
hp_mul_Methotrexate<- test_Methotrexate_mul |> filter((as.numeric(test_Methotrexate_mul$FDR_mul) < 0.1)&(as.numeric(test_Methotrexate_mul$Het.Pval_mul)>0.05))|> 
  select(Taxa,Drug,beta_mul,FDR_mul) |> mutate(`Influence`=ifelse(beta_mul<0,"Inhibition","Promotion"),anti="no")

Methotrexate_merge<-drug_merge(Methotrexate,Methotrexate_mul)
write.csv(Methotrexate_merge,file = "drug/Methotrexate.csv")

```

### Trimethoprim.sulfamethoxazole

```{r}
trimethoprim<-as.data.frame(metaanalysis_data[,,"Trimethoprim.sulfamethoxazole"])
trimethoprim<-drugmeta_hsctgvhd(trimethoprim)
trimethoprim$Drug<-"Trimethoprim.sulfamethoxazole"
test_trimethoprim<-  trimethoprim[complete.cases(trimethoprim$Het.Pval),] 
test_trimethoprim_n<-sum((as.numeric(test_trimethoprim$FDR) < 0.1)&(as.numeric(test_trimethoprim$Het.Pval)>0.05))
hp_trimethoprim<- test_trimethoprim |> filter((as.numeric(test_trimethoprim$FDR) < 0.1)&(as.numeric(test_trimethoprim$Het.Pval)>0.05))|> 
  select(Taxa,Drug,beta,FDR) |> mutate(`Influence`=ifelse(beta<0,"Inhibition","Promotion"),anti="yes")

trimethoprim_mul<-as.data.frame(metaanalysis_data_mul[,,"Trimethoprim.sulfamethoxazole"])
trimethoprim_mul<-drugmeta_mul_hsctgvhd(trimethoprim_mul)
trimethoprim_mul$Drug<-"Trimethoprim.sulfamethoxazole"
test_trimethoprim_mul<-  trimethoprim_mul[complete.cases(trimethoprim_mul$Het.Pval_mul),] 
test_trimethoprim_mul_n<-sum((as.numeric(test_trimethoprim_mul$FDR_mul) < 0.1)&(as.numeric(test_trimethoprim_mul$Het.Pval_mul)>0.05))
hp_mul_trimethoprim<- test_trimethoprim_mul |> filter((as.numeric(test_trimethoprim_mul$FDR_mul) < 0.1)&(as.numeric(test_trimethoprim_mul$Het.Pval_mul)>0.05))|> 
  select(Taxa,Drug,beta_mul,FDR_mul) |> mutate(`Influence`=ifelse(beta_mul<0,"Inhibition","Promotion"),anti="yes")

trimethoprim_merge<-drug_merge(trimethoprim,trimethoprim_mul)
write.csv(trimethoprim_merge,file = "drug/trimethoprim.csv")

```

### phenoxymethylpenicillin

```{r}
phenoxymethylpenicillin<-as.data.frame(metaanalysis_data[,,"Phenoxymethylpenicillin"])
phenoxymethylpenicillin<-drugmeta_hsctgvhd(phenoxymethylpenicillin)
phenoxymethylpenicillin$Drug<-"Phenoxymethylpenicillin"
test_phenoxymethylpenicillin<-  phenoxymethylpenicillin[complete.cases(phenoxymethylpenicillin$Het.Pval),] 
test_phenoxymethylpenicillin_n<-sum((as.numeric(test_phenoxymethylpenicillin$FDR) < 0.1)&(as.numeric(test_phenoxymethylpenicillin$Het.Pval)>0.05))
hp_phenoxymethylpenicillin<- test_phenoxymethylpenicillin |> filter((as.numeric(test_phenoxymethylpenicillin$FDR) < 0.1)&(as.numeric(test_phenoxymethylpenicillin$Het.Pval)>0.05))|> 
  select(Taxa,Drug,beta,FDR) |> mutate(`Influence`=ifelse(beta<0,"Inhibition","Promotion"),anti="yes")

phenoxymethylpenicillin_mul<-as.data.frame(metaanalysis_data_mul[,,"Phenoxymethylpenicillin"])
phenoxymethylpenicillin_mul<-drugmeta_mul_hsctgvhd(phenoxymethylpenicillin_mul)
phenoxymethylpenicillin_mul$Drug<-"Phenoxymethylpenicillin"
test_phenoxymethylpenicillin_mul<-  phenoxymethylpenicillin_mul[complete.cases(phenoxymethylpenicillin_mul$Het.Pval_mul),] 
test_phenoxymethylpenicillin_mul_n<-sum((as.numeric(test_phenoxymethylpenicillin_mul$FDR_mul) < 0.1)&(as.numeric(test_phenoxymethylpenicillin_mul$Het.Pval_mul)>0.05))
hp_mul_phenoxymethylpenicillin<- test_phenoxymethylpenicillin_mul |> filter((as.numeric(test_phenoxymethylpenicillin_mul$FDR_mul) < 0.1)&(as.numeric(test_phenoxymethylpenicillin_mul$Het.Pval_mul)>0.05))|> 
  select(Taxa,Drug,beta_mul,FDR_mul) |> mutate(`Influence`=ifelse(beta_mul<0,"Inhibition","Promotion"),anti="yes")

phenoxymethylpenicillin_merge<-drug_merge(phenoxymethylpenicillin,phenoxymethylpenicillin_mul)
write.csv(phenoxymethylpenicillin_merge,file = "drug/phenoxymethylpenicillin.csv")

```

### amoxicillin

```{r}
amoxicillin<-as.data.frame(metaanalysis_data[,,"Amoxicillin"])
amoxicillin<-drugmeta_hsctgvhd(amoxicillin)
amoxicillin$Drug<-"Amoxicillin"
test_amoxicillin<-  amoxicillin[complete.cases(amoxicillin$Het.Pval),] 
test_amoxicillin_n<-sum((as.numeric(test_amoxicillin$FDR) < 0.1)&(as.numeric(test_amoxicillin$Het.Pval)>0.05))
hp_amoxicillin<- test_amoxicillin |> filter((as.numeric(test_amoxicillin$FDR) < 0.1)&(as.numeric(test_amoxicillin$Het.Pval)>0.05))|> 
  select(Taxa,Drug,beta,FDR) |> mutate(`Influence`=ifelse(beta<0,"Inhibition","Promotion"),anti="yes")

amoxicillin_mul<-as.data.frame(metaanalysis_data_mul[,,"Amoxicillin"])
amoxicillin_mul<-drugmeta_mul_hsctgvhd(amoxicillin_mul)
amoxicillin_mul$Drug<-"Amoxicillin"
test_amoxicillin_mul<-  amoxicillin_mul[complete.cases(amoxicillin_mul$Het.Pval_mul),] 
test_amoxicillin_mul_n<-sum((as.numeric(test_amoxicillin_mul$FDR_mul) < 0.1)&(as.numeric(test_amoxicillin_mul$Het.Pval_mul)>0.05))
hp_mul_amoxicillin<- test_amoxicillin_mul |> filter((as.numeric(test_amoxicillin_mul$FDR_mul) < 0.1)&(as.numeric(test_amoxicillin_mul$Het.Pval_mul)>0.05))|> 
  select(Taxa,Drug,beta_mul,FDR_mul) |> mutate(`Influence`=ifelse(beta_mul<0,"Inhibition","Promotion"),anti="yes")

amoxicillin_merge<-drug_merge(amoxicillin,amoxicillin_mul)
write.csv(amoxicillin_merge,file = "drug/amoxicillin.csv")

```

### ceftazidime

```{r}
ceftazidime<-as.data.frame(metaanalysis_data[,,"Ceftazidime"])
ceftazidime<-drugmeta_hsctgvhd(ceftazidime)
ceftazidime$Drug<-"Ceftazidime"
test_ceftazidime<-  ceftazidime[complete.cases(ceftazidime$Het.Pval),] 
test_ceftazidime_n<-sum((as.numeric(test_ceftazidime$FDR) < 0.1)&(as.numeric(test_ceftazidime$Het.Pval)>0.05))
hp_ceftazidime<- test_ceftazidime |> filter((as.numeric(test_ceftazidime$FDR) < 0.1)&(as.numeric(test_ceftazidime$Het.Pval)>0.05))|> 
  select(Taxa,Drug,beta,FDR) |> mutate(`Influence`=ifelse(beta<0,"Inhibition","Promotion"),anti="yes")

ceftazidime_mul<-as.data.frame(metaanalysis_data_mul[,,"Ceftazidime"])
ceftazidime_mul<-drugmeta_mul_hsctgvhd(ceftazidime_mul)
ceftazidime_mul$Drug<-"Ceftazidime"
test_ceftazidime_mul<-  ceftazidime_mul[complete.cases(ceftazidime_mul$Het.Pval_mul),] 
test_ceftazidime_mul_n<-sum((as.numeric(test_ceftazidime_mul$FDR_mul) < 0.1)&(as.numeric(test_ceftazidime_mul$Het.Pval_mul)>0.05))
hp_mul_ceftazidime<- test_ceftazidime_mul |> filter((as.numeric(test_ceftazidime_mul$FDR_mul) < 0.1)&(as.numeric(test_ceftazidime_mul$Het.Pval_mul)>0.05))|> 
  select(Taxa,Drug,beta_mul,FDR_mul) |> mutate(`Influence`=ifelse(beta_mul<0,"Inhibition","Promotion"),anti="yes")

ceftazidime_merge<-drug_merge(ceftazidime,ceftazidime_mul)
write.csv(ceftazidime_merge,file = "drug/ceftazidime.csv")

```

### Ciprofloxacin

```{r}
Ciprofloxacin<-as.data.frame(metaanalysis_data[,,"Ciprofloxacin"])
Ciprofloxacin<-drugmeta_hsctgvhd(Ciprofloxacin)
Ciprofloxacin$Drug<-"Ciprofloxacin"
test_Ciprofloxacin<-  Ciprofloxacin[complete.cases(Ciprofloxacin$Het.Pval),] 
test_Ciprofloxacin_n<-sum((as.numeric(test_Ciprofloxacin$FDR) < 0.1)&(as.numeric(test_Ciprofloxacin$Het.Pval)>0.05))
hp_Ciprofloxacin<- test_Ciprofloxacin |> filter((as.numeric(test_Ciprofloxacin$FDR) < 0.1)&(as.numeric(test_Ciprofloxacin$Het.Pval)>0.05))|> 
  select(Taxa,Drug,beta,FDR) |> mutate(`Influence`=ifelse(beta<0,"Inhibition","Promotion"),anti="yes")

Ciprofloxacin_mul<-as.data.frame(metaanalysis_data_mul[,,"Ciprofloxacin"])
Ciprofloxacin_mul<-drugmeta_mul_hsctgvhd(Ciprofloxacin_mul)
Ciprofloxacin_mul$Drug<-"Ciprofloxacin"
test_Ciprofloxacin_mul<-  Ciprofloxacin_mul[complete.cases(Ciprofloxacin_mul$Het.Pval_mul),] 
test_Ciprofloxacin_mul_n<-sum((as.numeric(test_Ciprofloxacin_mul$FDR_mul) < 0.1)&(as.numeric(test_Ciprofloxacin_mul$Het.Pval_mul)>0.05))
hp_mul_Ciprofloxacin<- test_Ciprofloxacin_mul |> filter((as.numeric(test_Ciprofloxacin_mul$FDR_mul) < 0.1)&(as.numeric(test_Ciprofloxacin_mul$Het.Pval_mul)>0.05))|> 
  select(Taxa,Drug,beta_mul,FDR_mul) |> mutate(`Influence`=ifelse(beta_mul<0,"Inhibition","Promotion"),anti="yes")

Ciprofloxacin_merge<-drug_merge(Ciprofloxacin,Ciprofloxacin_mul)
write.csv(Ciprofloxacin_merge,file = "drug/Ciprofloxacin.csv")

```

### metronidazole

```{r}
metronidazole<-as.data.frame(metaanalysis_data[,,"Metronidazole"])
# metronidazole<-metronidazole[,c(1:4,9:24)]
# metronidazole<-metronidazole[complete.cases(metronidazole),]
metronidazole$inverse_var.chemo<-1/metronidazole$std_chemo^2
metronidazole$inverse_var.gvhd<-1/metronidazole$std_gvhd^2
#Calculate SE  #53
metronidazole$se=sqrt(1/(metronidazole$inverse_var.chemo+metronidazole$inverse_var.gvhd))

#Calculate Beta #54
metronidazole$beta=(metronidazole$inverse_var.chemo*metronidazole$est_chemo+metronidazole$inverse_var.gvhd*metronidazole$est_gvhd)/(metronidazole$inverse_var.chemo+metronidazole$inverse_var.gvhd)
#	selection$beta=(selection$inverse_var.ibd*selection$coeff.correcting.all.IBD+selection$inverse_var.lld*selection$coeff.correcting.all.LLD)/(selection$inverse_var.ibd+selection$inverse_var.lld)

#Calculate Z-score #55
metronidazole$Z=metronidazole$beta/metronidazole$se

#Calculate meta p-value #56
metronidazole$P=2*pnorm(-abs(metronidazole$Z))
#Remove data only avalable in one cohort
# metronidazole=metronidazole[complete.cases(metronidazole$P), ]

#Adjust pvalue with FDR #57
metronidazole$FDR=p.adjust(metronidazole$P,method = "fdr")

#Create empty columns
metronidazole$Het.Q="NS" #58
metronidazole$Het.I2="NS" #59
metronidazole$Het.Pval="NS" #60

#Heterogeneity using Cochran's Q-test for meta-FDR < 0.1
for (i in 1:length(rownames(metronidazole))){
  if ((metronidazole$FDR[i]<0.1)&(!is.na(metronidazole$FDR[i]))){
    #Select coefficients
    TE=c( metronidazole[i,13], metronidazole[i,21])
    #Select Standart error
    SE=c( metronidazole[i,14], metronidazole[i,22])
    het=metagen(TE,SE)
    #Change the number here depending of the number of columns in your data, should match with column Het.I2
    metronidazole[i,33]=het$I2
    #Match Het.Q column number
    metronidazole[i,32]=het$Q 
    #Calculate p-value from Q calculation // Het.Pval
    metronidazole[i,34]=pchisq(het$Q,df=2,lower.tail=F)
  }else{
    metronidazole[i,33]=NA
    #Match Het.Q column number
    metronidazole[i,32]=NA 
    #Calculate p-value from Q calculation // Het.Pval
    metronidazole[i,34]=NA
    
  }
}

metronidazole <- metronidazole %>%
  mutate(Taxa=rownames(metronidazole),Drug="Metronidazole",
         # est_hsct="NS",
         # std_hsct="NS",
         # pval_hsct="NS",
         # tval_hsct="NS",
         inverse_var.hsct="NS") %>%
  select(Taxa,Drug,
         N_hsct,Nonzeros_hsct,Users_hsct,Nonusers_hsct,
         est_hsct,std_hsct,pval_hsct,tval_hsct,
         N_gvhd,Nonzeros_gvhd,Users_gvhd,Nonusers_gvhd,
         est_gvhd,std_gvhd,pval_gvhd,tval_gvhd,
         N_chemo,Nonzeros_chemo,Users_chemo,Nonusers_chemo,
         est_chemo,std_chemo,pval_chemo,tval_chemo,
         inverse_var.hsct,inverse_var.gvhd,inverse_var.chemo,
         se,beta,Z,P,FDR,everything())

test_metronidazole<-  metronidazole[complete.cases(metronidazole$Het.Pval),] 
test_metronidazole_n<-sum((as.numeric(test_metronidazole$FDR) < 0.1)&(as.numeric(test_metronidazole$Het.Pval)>0.05))
hp_metronidazole<- test_metronidazole |> filter((as.numeric(test_metronidazole$FDR) < 0.1)&(as.numeric(test_metronidazole$Het.Pval)>0.05))|> 
  select(Taxa,Drug,beta,FDR) |> mutate(`Influence`=ifelse(beta<0,"Inhibition","Promotion"),anti="yes")


write.csv(metronidazole,file="drug/metronidazole.csv")

```

### piperacillin

```{r}
zpiperacillin<-as.data.frame(metaanalysis_data[,,"Piperacillin.tazobactam"])
# piperacillin<-piperacillin[,c(1:4,9:12)]
# piperacillin<-piperacillin[complete.cases(piperacillin),]
piperacillin$inverse_var.chemo<-1/piperacillin$std_chemo^2
piperacillin$inverse_var.hsct<-1/piperacillin$std_hsct^2
#Calculate SE  #53
piperacillin$se=sqrt(1/(piperacillin$inverse_var.chemo+piperacillin$inverse_var.hsct))

#Calculate Beta #54
piperacillin$beta=(piperacillin$inverse_var.chemo*piperacillin$est_chemo+piperacillin$inverse_var.hsct*piperacillin$est_hsct)/(piperacillin$inverse_var.chemo+piperacillin$inverse_var.hsct)
#	selection$beta=(selection$inverse_var.ibd*selection$coeff.correcting.all.IBD+selection$inverse_var.lld*selection$coeff.correcting.all.LLD)/(selection$inverse_var.ibd+selection$inverse_var.lld)

#Calculate Z-score #55
piperacillin$Z=piperacillin$beta/piperacillin$se

#Calculate meta p-value #56
piperacillin$P=2*pnorm(-abs(piperacillin$Z))
#Remove data only avalable in one cohort
# piperacillin=piperacillin[complete.cases(piperacillin$P), ]

#Adjust pvalue with FDR #57
piperacillin$FDR=p.adjust(piperacillin$P,method = "fdr")

#Create empty columns
piperacillin$Het.Q="NS" #58
piperacillin$Het.I2="NS" #59
piperacillin$Het.Pval="NS" #60

#Heterogeneity using Cochran's Q-test for meta-FDR < 0.1
for (i in 1:length(rownames(piperacillin))){
  if ((piperacillin$FDR[i]<0.1)&(!is.na(piperacillin$FDR[i]))){
    #Select coefficients
    TE=c( piperacillin[i,5], piperacillin[i,21])
    #Select Standart error
    SE=c( piperacillin[i,6], piperacillin[i,22])
    het=metagen(TE,SE)
    #Change the number here depending of the number of columns in your data, should match with column Het.I2
    piperacillin[i,33]=het$I2
    #Match Het.Q column number
    piperacillin[i,32]=het$Q 
    #Calculate p-value from Q calculation // Het.Pval
    piperacillin[i,34]=pchisq(het$Q,df=2,lower.tail=F)
  }else{
    piperacillin[i,33]=NA
    #Match Het.Q column number
    piperacillin[i,32]=NA
    #Calculate p-value from Q calculation // Het.Pval
    piperacillin[i,34]=NA
    
  }
}

piperacillin <- piperacillin %>%
  mutate(Taxa=rownames(piperacillin),Drug="Piperacillin.tazobactam",
         # est_gvhd="NS",
         # std_gvhd="NS",
         # pval_gvhd="NS",
         # tval_gvhd="NS",
         inverse_var.gvhd="NS") %>%
  select(Taxa,Drug,
         N_hsct,Nonzeros_hsct,Users_hsct,Nonusers_hsct,
         est_hsct,std_hsct,pval_hsct,tval_hsct,
         N_gvhd,Nonzeros_gvhd,Users_gvhd,Nonusers_gvhd,
         est_gvhd,std_gvhd,pval_gvhd,tval_gvhd,
         N_chemo,Nonzeros_chemo,Users_chemo,Nonusers_chemo,
         est_chemo,std_chemo,pval_chemo,tval_chemo,
         inverse_var.hsct,inverse_var.gvhd,inverse_var.chemo,
         se,beta,Z,P,FDR,everything())

test_piperacillin<-  piperacillin[complete.cases(piperacillin$Het.Pval),] 
test_piperacillin_n<-sum((as.numeric(test_piperacillin$FDR) < 0.1)&(as.numeric(test_piperacillin$Het.Pval)>0.05))
hp_piperacillin<- test_piperacillin |> filter((as.numeric(test_piperacillin$FDR) < 0.1)&(as.numeric(test_piperacillin$Het.Pval)>0.05))|> 
  select(Taxa,Drug,beta,FDR) |> mutate(`Influence`=ifelse(beta<0,"Inhibition","Promotion"),anti="yes")

write.csv(piperacillin,file="drug/piperacillin.csv")

```

### meropenem

```{r}
meropenem<-as.data.frame(metaanalysis_data[,,"Meropenem"])
meropenem<-drugmeta_all(meropenem)
meropenem$Drug<-"Meropenem"
test_meropenem<-  meropenem[complete.cases(meropenem$Het.Pval),] 
test_meropenem_n<-sum((as.numeric(test_meropenem$FDR) < 0.1)&(as.numeric(test_meropenem$Het.Pval)>0.05))
hp_meropenem<- test_meropenem |> filter((as.numeric(test_meropenem$FDR) < 0.1)&(as.numeric(test_meropenem$Het.Pval)>0.05))|> 
  select(Taxa,Drug,beta,FDR) |> mutate(`Influence`=ifelse(beta<0,"Inhibition","Promotion"),anti="yes")

meropenem_mul<-as.data.frame(metaanalysis_data_mul[,,"Meropenem"])
meropenem_mul<-drugmeta_mul_hsctgvhd(meropenem_mul)
meropenem_mul$Drug<-"Meropenem"
test_meropenem_mul<-  meropenem_mul[complete.cases(meropenem_mul$Het.Pval_mul),] 
test_meropenem_mul_n<-sum((as.numeric(test_meropenem_mul$FDR_mul) < 0.1)&(as.numeric(test_meropenem_mul$Het.Pval_mul)>0.05))
hp_mul_meropenem<- test_meropenem_mul |> filter((as.numeric(test_meropenem_mul$FDR_mul) < 0.1)&(as.numeric(test_meropenem_mul$Het.Pval_mul)>0.05))|> 
  select(Taxa,Drug,beta_mul,FDR_mul) |> mutate(`Influence`=ifelse(beta_mul<0,"Inhibition","Promotion"),anti="yes")

meropenem_merge<-drug_merge(meropenem,meropenem_mul)
write.csv(meropenem_merge,file = "drug/meropenem.csv")

```

### antibiotic

```{r}
antibiotic<-as.data.frame(metaanalysis_data[,,"Antibiotic"])
antibiotic<-drugmeta_all(antibiotic)
antibiotic$Drug<-"Antibiotic"
test_antibiotic<-  antibiotic[complete.cases(antibiotic$Het.Pval),] 
test_antibiotic_n<-sum((as.numeric(test_antibiotic$FDR) < 0.1)&(as.numeric(test_antibiotic$Het.Pval)>0.05))
hp_antibiotic<- test_antibiotic |> filter((as.numeric(test_antibiotic$FDR) < 0.1)&(as.numeric(test_antibiotic$Het.Pval)>0.05))|> 
  select(Taxa,Drug,beta,FDR) |> mutate(`Influence`=ifelse(beta<0,"Inhibition","Promotion"),anti="yes")
write.csv(antibiotic,file = "drug/antibiotic.csv")

```

### vancomycin

```{r}
vancomycin<-as.data.frame(metaanalysis_data[,,"Vancomycin"])
vancomycin<-drugmeta_all(vancomycin)
vancomycin$Drug<-"Vancomycin"
test_vancomycin<-  vancomycin[complete.cases(vancomycin$Het.Pval),] 
test_vancomycin_n<-sum((as.numeric(test_vancomycin$FDR) < 0.1)&(as.numeric(test_vancomycin$Het.Pval)>0.05))
hp_vancomycin<- test_vancomycin |> filter((as.numeric(test_vancomycin$FDR) < 0.1)&(as.numeric(test_vancomycin$Het.Pval)>0.05))|> 
  select(Taxa,Drug,beta,FDR) |> mutate(`Influence`=ifelse(beta<0,"Inhibition","Promotion"),anti="yes")

vancomycin_mul<-as.data.frame(metaanalysis_data_mul[,,"Vancomycin"])
vancomycin_mul<-drugmeta_mul_hsctgvhd(vancomycin_mul)
vancomycin_mul$Drug<-"Vancomycin"
test_vancomycin_mul<-  vancomycin_mul[complete.cases(vancomycin_mul$Het.Pval_mul),] 
test_vancomycin_mul_n<-sum((as.numeric(test_vancomycin_mul$FDR_mul) < 0.1)&(as.numeric(test_vancomycin_mul$Het.Pval_mul)>0.05))
hp_mul_vancomycin<- test_vancomycin_mul |> filter((as.numeric(test_vancomycin_mul$FDR_mul) < 0.1)&(as.numeric(test_vancomycin_mul$Het.Pval_mul)>0.05))|> 
  select(Taxa,Drug,beta_mul,FDR_mul) |> mutate(`Influence`=ifelse(beta_mul<0,"Inhibition","Promotion"),anti="yes")

vancomycin_merge<-drug_merge(vancomycin,vancomycin_mul)
write.csv(vancomycin_merge,file = "drug/vancomycin.csv")

```

## Figure 3 (Summary of the number of associated microbial species)

```{r}
metaanalysis_data_1.1<-metaanalysis_data_1[!rowSums(is.na(metaanalysis_data_1[,c(7:10,15:18,23:26)]))==12,] 
n_hsct<-metaanalysis_data_1.1 %>%
  group_by(drug) %>%
  dplyr::count(pval_hsct<0.05)

n_gvhd<-metaanalysis_data_1.1 %>%
  group_by(drug) %>%
  dplyr::count(pval_gvhd<0.05)

n_chemo<-metaanalysis_data_1.1 %>%
  group_by(drug) %>%
  dplyr::count(pval_chemo<0.05)

metaanalysis_data_mul1.1<-metaanalysis_data_mul1[!rowSums(is.na(metaanalysis_data_mul1[,c(7:10,15:18,23:26)]))==12,] 
n_mul_hsct<-metaanalysis_data_mul1.1 %>%
  group_by(drug) %>%
  dplyr::count(pval_mul_hsct<0.05)

n_mul_gvhd<-metaanalysis_data_mul1.1 %>%
  group_by(drug) %>%
  dplyr::count(pval_mul_gvhd<0.05)

n_mul_chemo<-metaanalysis_data_mul1.1 %>%
  group_by(drug) %>%
  dplyr::count(pval_mul_chemo<0.05)



for (i in 1:68){
  if (is.na(n_hsct$`pval_hsct < 0.05`[i])){
    n_hsct$n[i] = 0
  }
}

test_hsct<-dplyr::filter(n_hsct,`pval_hsct < 0.05`==T |is.na(`pval_hsct < 0.05`)) %>%
  filter(!duplicated(drug) | `pval_hsct < 0.05`==T)
#test_hsct<-rbind(test_hsct,n_hsct[64,])
# test_hsct$n[30]<-0

for (i in 1:67){
  if (is.na(n_gvhd$`pval_gvhd < 0.05`[i])){
    n_gvhd$n[i] = 0
  }
}

test_gvhd<-dplyr::filter(n_gvhd,`pval_gvhd < 0.05`==T |is.na(`pval_gvhd < 0.05`)) %>%
  filter(!duplicated(drug) | `pval_gvhd < 0.05`==T)

for (i in 1:50){
  if (is.na(n_chemo$`pval_chemo < 0.05`[i])){
    n_chemo$n[i] = 0
  }
}

test_chemo<-dplyr::filter(n_chemo,`pval_chemo < 0.05`==T |is.na(`pval_chemo < 0.05`)) %>%
  filter(!duplicated(drug) | `pval_chemo < 0.05`==T)

# setdiff(n_chemo$drug[!duplicated(n_chemo$drug)],test_chemo$drug)
# test_chemo<-rbind(test_chemo,n_chemo[c(1,2,15,16),])
#test_chemo$n[c(27,28,29,30)]<-0

for (i in 1:57){
  if (is.na(n_mul_hsct$`pval_mul_hsct < 0.05`[i])){
    n_mul_hsct$n[i] = 0
  }
}

test_mul_hsct<-dplyr::filter(n_mul_hsct,`pval_mul_hsct < 0.05`==T |is.na(`pval_mul_hsct < 0.05`)) %>%
  filter(!duplicated(drug) | `pval_mul_hsct < 0.05`==T)
# test_mul_hsct<-rbind(test_mul_hsct,n_mul_hsct[56,])
# test_mul_hsct$n[25]<-0


for (i in 1:56){
  if (is.na(n_mul_gvhd$`pval_mul_gvhd < 0.05`[i])){
    n_mul_gvhd$n[i] = 0
  }
}

test_mul_gvhd<-dplyr::filter(n_mul_gvhd,`pval_mul_gvhd < 0.05`==T |is.na(`pval_mul_gvhd < 0.05`)) %>%
  filter(!duplicated(drug) | `pval_mul_gvhd < 0.05`==T)


for (i in 1:25){
  if (is.na(n_mul_chemo$`pval_mul_chemo < 0.05`[i])){
    n_mul_chemo$n[i] = 0
  }
}

test_mul_chemo<-dplyr::filter(n_mul_chemo,`pval_mul_chemo < 0.05`==T |is.na(`pval_mul_chemo < 0.05`)) %>%
  filter(!duplicated(drug) | `pval_mul_chemo < 0.05`==T)



test_allcohort<-merge(test_hsct,test_gvhd,by="drug")
test_allcohort<-merge(test_allcohort,test_chemo,by="drug")
test_allcohort<-test_allcohort[,-c(2,4,6)]
colnames(test_allcohort)<-c("drug","HSCT","GVHD","Chemo")
test_allcohort<-dplyr::mutate(test_allcohort,Meta=c(0,6,0,3,0,10,
                                                    0,4,0,0,0,5,
                                                    7,0,7,0,0,0,
                                                    13,7,20,5,4,1,
                                                    7,0,10,4,15,9))
test_allcohort<-test_allcohort[-c(1,14,16,17,27),]
test_allcohort_1<-melt(test_allcohort,id="drug")
colnames(test_allcohort_1)<-c("Meds","cohort","Associations")
# write.table(test_allcohort_1,file = "test_allcohort_1.csv")

test_allcohort_1$Cohort<-factor(test_allcohort_1$cohort,levels = c("Meta","Chemo","GVHD","HSCT"))
test_allcohort_1$Group<-"Single drug analysis"

test_mul_allcohort<-merge(test_mul_hsct,test_mul_gvhd,by="drug")
test_mul_allcohort<-merge(test_mul_allcohort,test_mul_chemo,by="drug")
test_mul_allcohort<-test_mul_allcohort[,-c(2,4,6)]
colnames(test_mul_allcohort)<-c("drug","HSCT","GVHD","Chemo")
test_mul_allcohort<-dplyr::mutate(test_mul_allcohort,Meta=c(2,0,10,4,0,0,
                                                            2,9,0,8,0,0,
                                                            3,17,0,5,0,9,
                                                            0,9,4,8,8))
test_mul_allcohort<-test_mul_allcohort[-20,]
test_mul_allcohort_1<-melt(test_mul_allcohort,id="drug")
colnames(test_mul_allcohort_1)<-c("Meds","cohort","Associations")
# write.table(test_mul_allcohort_1,file = "test_mul_allcohort_1.csv")

test_mul_allcohort_1$Cohort<-factor(test_mul_allcohort_1$cohort,levels = c("Meta","Chemo","GVHD","HSCT"))
test_mul_allcohort_1$Group<-"Multiple drug analysis"

### bar plot
test_merge_cohort<-rbind(test_allcohort_1,test_mul_allcohort_1)
test_merge_cohort$Group<-factor(test_merge_cohort$Group,levels = c("Single drug analysis","Multiple drug analysis"))
figure3_1<-ggplot(test_merge_cohort,aes(Meds,Associations,fill=Cohort)) +
  geom_bar(stat="identity", position=position_dodge(),alpha=1.0) +
  # scale_fill_manual(values=c( "#241c1f","#CCE6D6","#CCE3F0","#EECEC9"))+
  scale_fill_manual(values=c( "#241c1f","#299555","#66AAD3","#C65949")) +
  theme(panel.background = element_blank(),
        panel.grid =element_line(colour = "grey90"),
        axis.line = element_line(colour = "grey50"),
        axis.title = element_text(size=11,family="serif",face="bold"),
        axis.text = element_text(size=11,family="serif",face="bold"),
        #legend.position = 'none',
        legend.title = element_text(size=11,family="serif",face="bold"),
        legend.text = element_text(size=11,family="serif",face="bold"),
        plot.subtitle = element_text(size=11,family="serif",face="bold"),
        strip.text = element_text(size=11,family="serif",face="bold"),
        strip.background = element_rect(fill = NA),
        title = element_text(size=11,family="serif",face="bold")
  )+
  # labs(title = "Single Drug Analysis")+
  facet_grid(cols = vars(Group), drop=T, scales = "free") +
  coord_flip() #+ scale_x_discrete(limits = rev(levels(test_allcohort$meds)))
ggsave("Figure3.pdf",figure3_1,units="in", width=9, height=5, dpi=600,limitsize = FALSE)

### additional file 19
tmp<-test_merge_cohort |> arrange(Meds) |> mutate(Antibiotic=c(rep("no",68),rep("yes",120)))
tmp1<-tmp |> group_by(Group,cohort,Antibiotic) |> 
  summarise(
    sum=sum(Associations)
  )
tmp2<-dcast(tmp1,Group+Antibiotic~cohort)

```

## Figure 4 (Overview of the microbial species with changed abundance caused by medicines used in independent cohorts)

> Packages used in this section

```{r}
library(ggplot2)
library(dplyr)
library(reshape2)
library(networkD3)
library(ggsci)
library(scales)

```

```{r}
fil_hsct<-metaanalysis_data_1.1 %>%
  group_by(drug) %>%
  dplyr::filter(pval_hsct<0.05) 
fil_hsct_e<-fil_hsct[,c(1,2,7)]
fil_hsct_e<-dcast(fil_hsct_e,taxa~drug)
write.csv(fil_hsct_e,file = "F:/数据/document/学校/硕士/数据库/r数据库/gvhd/sankeyplot/fil_hsct_e_s.csv") #excel修改行名
fil_hsct_s<-fil_hsct[,c(1,2,8)]
fil_hsct_s<-dcast(fil_hsct_s,taxa~drug)
write.csv(fil_hsct_s,file = "F:/数据/document/学校/硕士/数据库/r数据库/gvhd/sankeyplot/fil_hsct_s_s.csv") #excel修改行名

fil_gvhd<-metaanalysis_data_1.1 %>%
  group_by(drug) %>%
  dplyr::filter(pval_gvhd<0.05) 
fil_gvhd_e<-fil_gvhd[,c(1,2,15)]
fil_gvhd_e<-dcast(fil_gvhd_e,taxa~drug) 
write.csv(fil_gvhd_e,file = "F:/数据/document/学校/硕士/数据库/r数据库/gvhd/sankeyplot/fil_gvhd_e_s.csv") #excel修改行名
fil_gvhd_s<-fil_gvhd[,c(1,2,16)]
fil_gvhd_s<-dcast(fil_gvhd_s,taxa~drug) 
write.csv(fil_gvhd_s,file = "F:/数据/document/学校/硕士/数据库/r数据库/gvhd/sankeyplot/fil_gvhd_s_s.csv") #excel修改行名


fil_chemo<-metaanalysis_data_1.1 %>%
  group_by(drug) %>%
  dplyr::filter(pval_chemo<0.05) 
fil_chemo_e<-fil_chemo[,c(1,2,23)]
fil_chemo_e<-dcast(fil_chemo_e,taxa~drug)
write.csv(fil_chemo_e,file = "F:/数据/document/学校/硕士/数据库/r数据库/gvhd/sankeyplot/fil_chemo_e_s.csv") #excel修改行名
fil_chemo_s<-fil_chemo[,c(1,2,24)]
fil_chemo_s<-dcast(fil_chemo_s,taxa~drug)
write.csv(fil_chemo_s,file = "F:/数据/document/学校/硕士/数据库/r数据库/gvhd/sankeyplot/fil_chemo_s_s.csv") #excel修改行名


fil_mul_hsct<-metaanalysis_data_mul1.1 %>%
  group_by(drug) %>%
  dplyr::filter(pval_mul_hsct<0.05) 
fil_mul_hsct_e<-fil_mul_hsct[,c(1,2,7)]
fil_mul_hsct_e<-dcast(fil_mul_hsct_e,taxa~drug)
write.csv(fil_mul_hsct_e,file = "F:/数据/document/学校/硕士/数据库/r数据库/gvhd/sankeyplot/fil_mul_hsct_e_s.csv") #excel修改行名
fil_mul_hsct_s<-fil_mul_hsct[,c(1,2,8)]
fil_mul_hsct_s<-dcast(fil_mul_hsct_s,taxa~drug)
write.csv(fil_mul_hsct_s,file = "F:/数据/document/学校/硕士/数据库/r数据库/gvhd/sankeyplot/fil_mul_hsct_s_s.csv") #excel修改行名


fil_mul_gvhd<-metaanalysis_data_mul1.1 %>%
  group_by(drug) %>%
  dplyr::filter(pval_mul_gvhd<0.05) 
fil_mul_gvhd_e<-fil_mul_gvhd[,c(1,2,15)]
fil_mul_gvhd_e<-dcast(fil_mul_gvhd_e,taxa~drug) 
write.csv(fil_mul_gvhd_e,file = "F:/数据/document/学校/硕士/数据库/r数据库/gvhd/sankeyplot/fil_mul_gvhd_e_s.csv") #excel修改行名
fil_mul_gvhd_s<-fil_mul_gvhd[,c(1,2,16)]
fil_mul_gvhd_s<-dcast(fil_mul_gvhd_s,taxa~drug) 
write.csv(fil_mul_gvhd_s,file = "F:/数据/document/学校/硕士/数据库/r数据库/gvhd/sankeyplot/fil_mul_gvhd_s_s.csv") #excel修改行名

fil_mul_chemo<-metaanalysis_data_mul1.1 %>%
  group_by(drug) %>%
  dplyr::filter(pval_mul_chemo<0.05) 
fil_mul_chemo_e<-fil_mul_chemo[,c(1,2,23)]
fil_mul_chemo_e<-dcast(fil_mul_chemo_e,taxa~drug)
write.csv(fil_mul_chemo_e,file = "F:/数据/document/学校/硕士/数据库/r数据库/gvhd/sankeyplot/fil_mul_chemo_e_s.csv") #excel修改行名
fil_mul_chemo_s<-fil_mul_chemo[,c(1,2,24)]
fil_mul_chemo_s<-dcast(fil_mul_chemo_s,taxa~drug)
write.csv(fil_mul_chemo_s,file = "F:/数据/document/学校/硕士/数据库/r数据库/gvhd/sankeyplot/fil_mul_chemo_s_s.csv") #excel修改行名

```

### Single drug analysis

```{r}
sk_hsct<-read.csv( "F:/数据/document/学校/硕士/数据库/r数据库/gvhd/sankeyplot/sk_hsct.csv",header=T)
sk_hsct<-sk_hsct[,-7]
sk_hsct<-melt(sk_hsct,ID=c("Taxa","NodeColor"))
sk_hsct<-na.omit(sk_hsct)
sk_hsct<-sk_hsct %>% mutate(trend=ifelse(value<0,"neg","pos")) %>% mutate(Cohort="HSCT") %>% mutate(antibioic=c(rep("no",50),rep("yes",64)))

sk_gvhd<-read.csv( "F:/数据/document/学校/硕士/数据库/r数据库/gvhd/sankeyplot/sk_gvhd.csv",header=T)
sk_gvhd<-sk_gvhd[,-c(8,18)]
sk_gvhd<-melt(sk_gvhd,ID=c("Taxa","NodeColor"))
sk_gvhd<-na.omit(sk_gvhd)
sk_gvhd<-sk_gvhd %>% mutate(trend=ifelse(value<0,"neg","pos")) %>% mutate(Cohort="GVHD") %>% mutate(antibioic=c(rep("no",27),rep("yes",66)))

sk_chemo<-read.csv( "F:/数据/document/学校/硕士/数据库/r数据库/gvhd/sankeyplot/sk_chemo.csv",header=T)
sk_chemo<-melt(sk_chemo,ID=c("Taxa","NodeColor"))
sk_chemo<-na.omit(sk_chemo)
sk_chemo<-sk_chemo %>% mutate(trend=ifelse(value<0,"neg","pos")) %>% mutate(Cohort="Chemo") %>% mutate(antibioic=c(rep("yes",27)))


```

#### Sankey plot

```{r}
sk_plt<-rbind(sk_hsct,sk_gvhd,sk_chemo)

sk_plt_link1<-sk_plt[,c(3,1,4,2,5,6,7)]
colnames(sk_plt_link1)<-c("source","target","value","node","trend","cohort","anti")
sk_plt_link2<-sk_plt[,c(1,6,4,2,5,6,7)]
colnames(sk_plt_link2)<-c("source","target","value","node","trend","cohort","anti")

node_micro<-sk_plt_link1[,c("target","node")] %>% filter(!duplicated(sk_plt_link1$target) ==T)
colnames(node_micro)=c("name","group")
node_micro$group<-factor(node_micro$group,levels=c("Bacteroidetes","Firmicutes","Actinobacteria","Proteobacteria","Verrucomicrobia"))
node_micro<-node_micro |> arrange(group)
node_drug<-sk_plt_link1[,c("source","anti")] %>% filter(!duplicated(sk_plt_link1$source) ==T)
colnames(node_drug)=c("name","group")
node_drug$name<-as.character(node_drug$name)
node_drug<-node_drug |> arrange(group,name)
node_cohort<-data.frame(name=c("HSCT","GVHD","Chemo"),group=c("HSCT","GVHD","Chemo"))
sk_plt_nodes<-rbind(node_drug,node_micro,node_cohort)
sk_plt_nodes$name<-as.character(sk_plt_nodes$name)
sk_plt_nodes$group<-as.factor(sk_plt_nodes$group)

sk_plt_link1<-sk_plt_link1 %>% mutate(group=paste(trend,anti, sep = "_", collapse = NULL))
sk_plt_link2<-sk_plt_link2 %>% mutate(group=paste(trend,cohort, sep = "_", collapse = NULL))
sk_plt_link<-rbind(sk_plt_link1,sk_plt_link2)
sk_plt_link<-sk_plt_link[,-c(4:7)]
sk_plt_link$group<-as.factor(sk_plt_link$group)
sk_plt_link$source<-as.character(sk_plt_link$source)
sk_plt_link$value<-abs(log(abs(sk_plt_link$value)))
sk_plt_link$value<-as.character(sk_plt_link$value)

sk_plt_link$IDsource <- match(sk_plt_link$source, sk_plt_nodes$name)-1 
sk_plt_link$IDtarget <- match(sk_plt_link$target, sk_plt_nodes$name)-1

sk_plt_color <- 'd3.scaleOrdinal() .domain(["HSCT","GVHD","Chemo","yes","no",
                                            "Firmicutes","Bacteroidetes","Actinobacteria","Verrucomicrobia","Proteobacteria",
                                            "pos_HSCT","pos_GVHD", "pos_Chemo","neg_HSCT","neg_GVHD", "neg_Chemo",
                                            "pos_no","pos_yes","neg_no","neg_yes"]) 
                                   .range(["#BC3C29FF", "#0072B5FF", "#008134FF", "#FFDC91FF", "#7876B1FF",
                                          "#E64B35FF","#4DBBD5FF","#6CB86FFF","#3C5488FF","#F39B7FFF",
                                          "#BC3C294C", "#0072B54C", "#0081344C","#BC3C294C", "#0072B54C", "#0081344C",
                                          "#F39B7FFF","#F39B7F4C","#8491B4E5","#8491B44C"])'
sk <- sankeyNetwork(Links = sk_plt_link, Nodes = sk_plt_nodes, Source = "IDsource", Target = "IDtarget", 
                    Value = "value", nodePadding = 7,fontFamily = "Times New Roman",
                    NodeID = "name", iterations=0,fontSize = 14,sinksRight = TRUE,
                    colourScale=sk_plt_color, 
                    units = "TWh",height = 900,width = 1200,
                    LinkGroup="group", NodeGroup="group")
sk
```

#### Vennplot

```{r}
library(ggVennDiagram)

vp_hsct<-as.character(fil_hsct_e$taxa)
vp_sk_hsct<-unique(sk_hsct$Taxa)
vp_gvhd<-as.character(fil_gvhd_e$taxa)
vp_sk_gvhd<-unique(sk_gvhd$Taxa)
vp_chemo<-as.character(fil_chemo_e$taxa)
vp_sk_chemo<-unique(sk_chemo$Taxa)

intersect(intersect(vp_hsct,vp_gvhd),vp_chemo)

vp<-list(
  HSCT=vp_hsct,
  GVHD=vp_gvhd,
  Chemo=vp_chemo
)

vp_sk<-list(
  HSCT=vp_sk_hsct,
  GVHD=vp_sk_gvhd,
  Chemo=vp_sk_chemo
)


vp_plt<-ggVennDiagram(vp, 
              category.names = c("HSCT","GVHD","Chemo"),
              set_color = c("#BC3C29FF", "#0072B5FF", "#008134FF"),
              set_size = 6,
              label = "both", 
              label_percent_digit = 1, 
              label_size = 3,
              )+
  scale_x_continuous(expand = expansion(mult = .2))+
  scale_fill_gradient(low="white",high = "#241c1f",limits=c(0,40),name = "Species Number")+
  theme(legend.title = element_text(size=16,family="serif",face="bold"),
        legend.text = element_text(size=15,family="serif",face="bold"),
        plot.subtitle = element_text(size=12,family="serif",face="bold")
        
  )
vp_plt

vp_sk_plt<-ggVennDiagram(vp_sk, 
                      category.names = c("HSCT","GVHD","Chemo"),
                      set_color = c("#BC3C29FF", "#0072B5FF", "#008134FF"),
                      set_size = 6,
                      label = "both", 
                      label_percent_digit = 1, 
                      label_size = 3,
)+
  scale_x_continuous(expand = expansion(mult = .2))+
  scale_fill_gradient(low="white",high = "#241c1f",limits=c(0,40),name = "Species Number")+
  theme(legend.title = element_text(size=16,family="serif",face="bold"),
        legend.text = element_text(size=15,family="serif",face="bold"),
        plot.subtitle = element_text(size=12,family="serif",face="bold")
        
  )
vp_sk_plt
Figure4_venn<-vp_plt|vp_sk_plt+plot_layout(guides = 'collect')
ggsave("Figure4_venn.pdf",Figure4_venn,units="in", width=12, height=6, dpi=600,limitsize = FALSE)

```

### multiple drug analysis

```{r}
sk_mul_hsct<-read.csv( "F:/数据/document/学校/硕士/数据库/r数据库/gvhd/sankeyplot/sk_mul_hsct.csv",header=T)
sk_mul_hsct<-sk_mul_hsct[,-6]
sk_mul_hsct<-melt(sk_mul_hsct,ID=c("Taxa","NodeColor"))
sk_mul_hsct<-na.omit(sk_mul_hsct)
sk_mul_hsct<-sk_mul_hsct %>% mutate(trend=ifelse(value<0,"neg","pos")) %>% mutate(Cohort="HSCT") %>% mutate(antibioic=c(rep("no",23),rep("yes",56)))

sk_mul_gvhd<-read.csv( "F:/数据/document/学校/硕士/数据库/r数据库/gvhd/sankeyplot/sk_mul_gvhd.csv",header=T)
sk_mul_gvhd<-sk_mul_gvhd[,-c(8,17)]
sk_mul_gvhd<-melt(sk_mul_gvhd,ID=c("Taxa","NodeColor"))
sk_mul_gvhd<-na.omit(sk_mul_gvhd)
sk_mul_gvhd<-sk_mul_gvhd %>% mutate(trend=ifelse(value<0,"neg","pos")) %>% mutate(Cohort="GVHD") %>% mutate(antibioic=c(rep("no",32),rep("yes",60)))

sk_mul_chemo<-read.csv( "F:/数据/document/学校/硕士/数据库/r数据库/gvhd/sankeyplot/sk_mul_chemo.csv",header=T)
sk_mul_chemo<-melt(sk_mul_chemo,ID=c("Taxa","NodeColor"))
sk_mul_chemo<-na.omit(sk_mul_chemo)
sk_mul_chemo<-sk_mul_chemo %>% mutate(trend=ifelse(value<0,"neg","pos")) %>% mutate(Cohort="Chemo") %>% mutate(antibioic=c(rep("no",5)))


```

#### Sankey plot

```{r}
sk_mul_plt<-rbind(sk_mul_hsct,sk_mul_gvhd,sk_mul_chemo)

sk_mul_plt_link1<-sk_mul_plt[,c(3,1,4,2,5,6,7)]
colnames(sk_mul_plt_link1)<-c("source","target","value","node","trend","cohort","anti")
sk_mul_plt_link2<-sk_mul_plt[,c(1,6,4,2,5,6,7)]
colnames(sk_mul_plt_link2)<-c("source","target","value","node","trend","cohort","anti")

node_mul_micro<-sk_mul_plt_link1[,c("target","node")] %>% filter(!duplicated(sk_mul_plt_link1$target) ==T)
colnames(node_mul_micro)=c("name","group")
node_mul_micro$group<-factor(node_mul_micro$group,levels=c("Bacteroidetes","Firmicutes","Actinobacteria","Proteobacteria","Verrucomicrobia"))
node_mul_micro<-node_mul_micro |> arrange(group)
node_mul_drug<-sk_mul_plt_link1[,c("source","anti")] %>% filter(!duplicated(sk_mul_plt_link1$source) ==T)
colnames(node_mul_drug)=c("name","group")
node_mul_drug$name<-as.character(node_mul_drug$name)
node_mul_drug<-node_mul_drug |> arrange(group,name)
node_mul_cohort<-data.frame(name=c("HSCT","GVHD","Chemo"),group=c("HSCT","GVHD","Chemo"))
sk_mul_plt_nodes<-rbind(node_mul_drug,node_mul_micro,node_mul_cohort)
sk_mul_plt_nodes$name<-as.character(sk_mul_plt_nodes$name)
sk_mul_plt_nodes$group<-as.factor(sk_mul_plt_nodes$group)


sk_mul_plt_link1<-sk_mul_plt_link1 %>% mutate(group=paste(trend,anti, sep = "_", collapse = NULL))
sk_mul_plt_link2<-sk_mul_plt_link2 %>% mutate(group=paste(trend,cohort, sep = "_", collapse = NULL))
sk_mul_plt_link<-rbind(sk_mul_plt_link1,sk_mul_plt_link2)
sk_mul_plt_link<-sk_mul_plt_link[,-c(4:7)]
sk_mul_plt_link$group<-as.factor(sk_mul_plt_link$group)
sk_mul_plt_link$source<-as.character(sk_mul_plt_link$source)
sk_mul_plt_link$value<-abs(log(abs(sk_mul_plt_link$value)))
sk_mul_plt_link$value<-as.character(sk_mul_plt_link$value)

sk_mul_plt_link$IDsource <- match(sk_mul_plt_link$source, sk_mul_plt_nodes$name)-1 
sk_mul_plt_link$IDtarget <- match(sk_mul_plt_link$target, sk_mul_plt_nodes$name)-1

sk_mul_plt_color <- 'd3.scaleOrdinal() .domain(["HSCT","GVHD","Chemo","yes","no",
                                            "Firmicutes","Bacteroidetes","Actinobacteria","Verrucomicrobia","Proteobacteria",
                                            "pos_HSCT","pos_GVHD", "pos_Chemo","neg_HSCT","neg_GVHD", "neg_Chemo",
                                            "pos_no","pos_yes","neg_no","neg_yes"]) 
                                   .range(["#BC3C29FF", "#0072B5FF", "#008134FF", "#FFDC91FF", "#7876B1FF",
                                          "#E64B35FF","#4DBBD5FF","#6CB86FFF","#3C5488FF","#F39B7FFF",
                                          "#BC3C294C", "#0072B54C", "#0081344C","#BC3C294C", "#0072B54C", "#0081344C",
                                          "#F39B7FFF","#F39B7F4C","#8491B4E5","#8491B44C"])'

sk_mul <- sankeyNetwork(Links = sk_mul_plt_link, Nodes = sk_mul_plt_nodes, Source = "IDsource", Target = "IDtarget", 
                    Value = "value", nodePadding = 7,fontFamily = "Times New Roman",
                    NodeID = "name", iterations=0,fontSize = 14,sinksRight = TRUE,
                    colourScale=sk_mul_plt_color, 
                    units = "TWh",height = 900,width = 1200,
                    LinkGroup="group", NodeGroup="group")

sk_mul
```

#### Vennplot

```{r}
vp_mul_hsct<-as.character(fil_mul_hsct_e$taxa)
vp_mul_sk_hsct<-unique(sk_mul_hsct$Taxa)
vp_mul_gvhd<-as.character(fil_mul_gvhd_e$taxa)
vp_mul_sk_gvhd<-unique(sk_mul_gvhd$Taxa)
vp_mul_chemo<-as.character(fil_mul_chemo_e$taxa)
vp_mul_sk_chemo<-unique(sk_mul_chemo$Taxa)

intersect(intersect(vp_mul_hsct,vp_mul_gvhd),vp_mul_chemo)

vp_mul<-list(
  HSCT=vp_mul_hsct,
  GVHD=vp_mul_gvhd,
  Chemo=vp_mul_chemo
)

vp_mul_sk<-list(
  HSCT=vp_mul_sk_hsct,
  GVHD=vp_mul_sk_gvhd,
  Chemo=vp_mul_sk_chemo
)

vp_mul_plt<-ggVennDiagram(vp_mul, 
                      category.names = c("HSCT","GVHD","Chemo"),
                      set_color = c("#BC3C29FF", "#0072B5FF", "#008134FF"),
                      set_size = 6,
                      label = "both", 
                      label_percent_digit = 1, 
                      label_size = 3,
)+
  scale_x_continuous(expand = expansion(mult = .2))+
  scale_fill_gradient(low="white",high = "#241c1f",limits=c(0,60),name = "Species Number")+
  theme(legend.title = element_text(size=16,family="serif",face="bold"),
        legend.text = element_text(size=15,family="serif",face="bold"),
        plot.subtitle = element_text(size=12,family="serif",face="bold")
        
  )
vp_mul_plt

vp_mul_sk_plt<-ggVennDiagram(vp_mul_sk, 
                         category.names = c("HSCT","GVHD","Chemo"),
                         set_color = c("#BC3C29FF", "#0072B5FF", "#008134FF"),
                         set_size = 6,
                         label = "both", 
                         label_percent_digit = 1, 
                         label_size = 3,
)+
  scale_x_continuous(expand = expansion(mult = .2))+
  scale_fill_gradient(low="white",high = "#241c1f",limits=c(0,60),name = "Species Number")+
  theme(legend.title = element_text(size=16,family="serif",face="bold"),
        legend.text = element_text(size=15,family="serif",face="bold"),
        plot.subtitle = element_text(size=12,family="serif",face="bold")
        
  )
vp_mul_sk_plt
Figure4s_venn<-vp_mul_plt|vp_mul_sk_plt+plot_layout(guides = 'collect')
ggsave("Figure4s_venn.pdf",Figure4s_venn,units="in", width=12, height=6, dpi=600,limitsize = FALSE)

```

## Figure 5

```{r}
library(ade4)
library(ggtree)

tmp <- unique(c(as.character(fil_hsct_e$taxa),as.character(fil_gvhd_e$taxa),as.character(fil_chemo_e$taxa),
                as.character(fil_mul_hsct_e$taxa),as.character(fil_mul_gvhd_e$taxa),as.character(fil_mul_chemo_e$taxa)))
tmp1<-as.data.frame(order="1",
                    tmp)
# write.csv(tmp,"hp_seperate_species.csv")
# 
tree_seperate_species<-read.csv("heatmap/tree_seperate_species.csv")

tree_seperate_species$ID<-factor(tree_seperate_species$ID,levels = tree_seperate_species$ID)


rownames(tree_seperate_species)<-tree_seperate_species$ID
tree_seperate_species<-tree_seperate_species[,-1]
tree_seperate_species$p <-as.factor(tree_seperate_species$p)
tree_seperate_species$c <-as.factor(tree_seperate_species$c)
tree_seperate_species$o <-as.factor(tree_seperate_species$o)
tree_seperate_species$f <-as.factor(tree_seperate_species$f)
tree_seperate_species$g <-as.factor(tree_seperate_species$g)

tmptax1.phy <- taxo2phylog(as.taxo(tree_seperate_species))
print(tmptax1.phy)
tree_sep_plt_name<-ggtree(tmptax1.phy,colour = "grey50")  + 
  geom_tiplab(size=3.0)+
  geom_nodelab(size=3.0)+
  hexpand(.25)#show id out of edge
tree_sep_plt_name

tree_sep_plt_noname<-ggtree(tmptax1.phy,colour = "grey50",size=0.9)  + 
  coord_flip()
tree_sep_plt_noname


hp_sep_species_hsct<-melt(fil_hsct_e,ID=C("taxa","drug")) |> na.omit() |> 
  mutate(Cohort="HSCT",Group="Single drug analysis")
colnames(hp_sep_species_hsct)<-c("Taxa","Meds","Magnitude","Cohort","Group")
pval<-melt(fil_hsct[,c(1,2,9)],ID=C("taxa","drug")) |> na.omit()
hp_sep_species_hsct$`P value`<-pval$value

hp_sep_mul_species_hsct<-melt(fil_mul_hsct_e,ID=C("taxa","drug")) |> na.omit() |> 
  mutate(Cohort="HSCT",Group="Multiple drug analysis")
colnames(hp_sep_mul_species_hsct)<-c("Taxa","Meds","Magnitude","Cohort","Group")
pval<-melt(fil_mul_hsct[,c(1,2,9)],ID=C("taxa","drug")) |> na.omit()
hp_sep_mul_species_hsct$`P value`<-pval$value


hp_sep_species_gvhd<-melt(fil_gvhd_e,ID=C("taxa","drug")) |> na.omit() |> 
  mutate(Cohort="GVHD",Group="Single drug analysis")
colnames(hp_sep_species_gvhd)<-c("Taxa","Meds","Magnitude","Cohort","Group")
pval<-melt(fil_gvhd[,c(1,2,17)],ID=C("taxa","drug")) |> na.omit()
hp_sep_species_gvhd$`P value`<-pval$value

hp_sep_mul_species_gvhd<-melt(fil_mul_gvhd_e,ID=C("taxa","drug")) |> na.omit() |> 
  mutate(Cohort="GVHD",Group="Multiple drug analysis")
colnames(hp_sep_mul_species_gvhd)<-c("Taxa","Meds","Magnitude","Cohort","Group")
pval<-melt(fil_mul_gvhd[,c(1,2,17)],ID=C("taxa","drug")) |> na.omit()
hp_sep_mul_species_gvhd$`P value`<-pval$value


hp_sep_species_chemo<-melt(fil_chemo_e,ID=C("taxa","drug")) |> na.omit() |> 
  mutate(Cohort="Chemo",Group="Single drug analysis")
colnames(hp_sep_species_chemo)<-c("Taxa","Meds","Magnitude","Cohort","Group")
pval<-melt(fil_chemo[,c(1,2,25)],ID=C("taxa","drug")) |> na.omit()
hp_sep_species_chemo$`P value`<-pval$value

hp_sep_mul_species_chemo<-melt(fil_mul_chemo_e,ID=C("taxa","drug")) |> na.omit() |> 
  mutate(Cohort="Chemo",Group="Multiple drug analysis")
colnames(hp_sep_mul_species_chemo)<-c("Taxa","Meds","Magnitude","Cohort","Group")
pval<-melt(fil_mul_chemo[,c(1,2,25)],ID=C("taxa","drug")) |> na.omit()
hp_sep_mul_species_chemo$`P value`<-pval$value

hp_sep_merge_species<-rbind(hp_sep_species_hsct,hp_sep_mul_species_hsct,
                            hp_sep_species_gvhd,hp_sep_mul_species_gvhd,
                            hp_sep_species_chemo,hp_sep_mul_species_chemo) |> 
  filter(Meds!="Total.Body.Irradiation") |> 
  mutate(Influence=ifelse(Magnitude<0,"Inhibition","Promotion"))
write.csv(hp_sep_merge_species,"heatmap/hp_sep_merge_species.csv")


hp_sep_merge_species<-read.csv("heatmap/hp_sep_merge_species_mod.csv")
hp_sep_merge_species$Meds<-factor(hp_sep_merge_species$Meds,levels = levels(test_merge_cohort$Meds))
hp_sep_merge_species$Group<-factor(hp_sep_merge_species$Group,levels = c("Single drug analysis","Multiple drug analysis"))
hp_sep_merge_species$Cohort<-factor(hp_sep_merge_species$Cohort,levels = c("HSCT","GVHD","Chemo"))
hp_sep_merge_species$facet<-paste(hp_sep_merge_species$Cohort,hp_sep_merge_species$Group,sep = " ")
hp_sep_merge_species$facet<-factor(hp_sep_merge_species$facet,levels = c(
  "HSCT Single drug analysis","HSCT Multiple drug analysis",
  "GVHD Single drug analysis","GVHD Multiple drug analysis",
  "Chemo Single drug analysis","Chemo Multiple drug analysis"))
hp_sep_merge_species$ID<-factor(hp_sep_merge_species$ID,levels = c(  "k__Bacteria","Unassigned","o__Streptophyta","g__Methanobrevibacter","g__Fusobacterium","A. muciniphila","g__Sutterella","g__Bilophila","g__Campylobacter","H. parainfluenzae","f__Enterobacteriaceae","g__Shigella","g__Atopobium","C. aerofaciens","E. lenta", "f__Coriobacteriaceae","g__Actinomyces","g__Corynebacterium","R. mucilaginosa","R. dentocariosa","g__Scardovia","g__Bifidobacterium1","g__Bifidobacterium2","B. longum","f__Barnesiellaceae","g__Odoribacter","g__Parabacteroides","P. distasonis","g__Prevotella","P. copri","f__Rikenellaceae","A. finegoldii","A. indistinctus","A. onderdonkii","A. putredinis","g__Bacteroides1","g__Bacteroides2","B. caccae","B. fragilis","B. ovatus","B. uniformis","f__Erysipelotrichaceae","C. spiroforme","g__Coprobacillus","E. dolichum","g__Holdemania","g__Staphylococcus","g__Turicibacter","o__Lactobacillales","g__Granulicatella","f__Enterococcaceae","g__Enterococcus1","g__Enterococcus2","E. casseliflavus","g__Lactococcus","g__Streptococcus1","S. infantis","g__Streptococcus2","f__Lactobacillaceae","g__Pediococcus","g__Lactobacillus1","g__Lactobacillus2","L. delbrueckii","L. zeae","L. salivarius","o__Clostridiales1","o__Clostridiales2","f__Christensenellaceae","g__Anaerofustis","f__Mogibacteriaceae","g__Anaerococcus","g__Finegoldia","f__Peptostreptococcaceae1","f__Peptostreptococcaceae2","C. bartlettii","f__Clostridiaceae","g__SMB53","C. celatum","C. paraputrificum","g__Clostridium1","g__Clostridium2","g__Acidaminococcus","g__Dialister","g__Megasphaera","g__Phascolarctobacterium","V. parvula","V. dispar","f__Ruminococcaceae1","f__Ruminococcaceae2","g__Anaerotruncus","B. pullicaecorum","C. methylpentosum","G. formicilis","g__Oscillospira","S. variabile","g__Faecalibacterium","F. prausnitzii","g__Ruminococcus1","R. bromii","f__Lachnospiraceae1","f__Lachnospiraceae2","g__Coprococcus","g__Lachnospira","g__Blautia","B. obeum","B. producta","g__Dorea","D. longicatena","D. formicigenerans","g__Roseburia1","g__Roseburia2","R. faecis","g__Ruminococcus2","R. gnavus","R. torques","g__Clostridium","C. aldenense","C. citroniae","C. clostridioforme","C. hathewayi","C. lavalense","C. symbiosum"))
hp_sep_merge_species_plt<-hp_sep_merge_species |> filter(facet!="Chemo Multiple drug analysis")
temp<-ggplot(hp_sep_merge_species_plt,aes(x=ID,y=Meds,color=Influence))+
  geom_point(aes(size=abs(Magnitude),alpha=-P.value))+
  scale_size_area(max_size = 7)+
  # scale_size_continuous(name="Magnitude")+
  scale_alpha_continuous(name="P value",
                         limits = c(-1,0), breaks = c(0, -0.01, -0.1,-0.5, -1),
                         guide = guide_colourbar(nbin = 100, draw.ulim = FALSE, draw.llim = FALSE))+
  scale_color_manual(name="Influence",
                     labels=c("Inhibition","Promotion"),
                     values=c("#8491B4FF","#F39B7FFF")) +
  # scale_x_discrete(name="Taxa influenced in three cohorts",
  #                  aes(color=textcol))+
  theme(panel.background = element_blank(),
        panel.grid =element_blank(),
        axis.line = element_blank(),
        # axis.title.y = element_text(size=14,family="serif",face="bold"),
        axis.title = element_blank(),
        axis.text.y =  element_text(size=6,family="serif",face="bold"),
        axis.text.x =  element_text(size=10,family="serif",face="bold",angle=90,hjust = 1,vjust = 0.5),
        #legend.position = 'none',
        legend.title = element_text(size=16,family="serif",face="bold"),
        legend.text = element_text(size=15,family="serif",face="bold"),
        plot.subtitle = element_text(size=12,family="serif",face="bold"),
        strip.text = element_text(size=10,family="serif",face="bold"),
        strip.background = element_rect(fill = NA),
        title = element_text(size=14,family="serif",face="bold")
  )+
  labs(size="Magnitude")+
  facet_grid(rows = vars(facet), drop=T)+
  guides(colour = guide_legend(override.aes = list(size=5)),
         alpha=guide_legend(override.aes = list(size=5)))#修改图例点大小

temp
write.csv(hp_sep_merge_species_plt,"figure4data.csv")
b_sep<-temp/tree_sep_plt_noname+plot_layout(heights = c(10.5, 1.2))
```

## 5.2 Species meta-analysis (Figure 5s)

### Single drug analysis

```{r}
meta_hsct_e<-read.csv( "F:/数据/document/学校/硕士/数据库/r数据库/gvhd/sankeyplot/sk_hsct.csv",header=T)
meta_hsct_e<-meta_hsct_e[,-7]
meta_hsct_e<-melt(meta_hsct_e,ID=c("Taxa","NodeColor"))
meta_hsct_e<-na.omit(meta_hsct_e)
meta_hsct_e<-meta_hsct_e %>% mutate(Cohort="HSCT") %>% mutate(antibioic=c(rep("no",50),rep("yes",64)))

meta_hsct_s<-read.csv( "F:/数据/document/学校/硕士/数据库/r数据库/gvhd/sankeyplot/sk_hsct_s.csv",header=T)
meta_hsct_s<-meta_hsct_s[,-7]
meta_hsct_s<-melt(meta_hsct_s,ID=c("Taxa","NodeColor"))
meta_hsct_s<-na.omit(meta_hsct_s)

meta_hsct<-cbind(meta_hsct_e,meta_hsct_s$value)
colnames(meta_hsct)<-c("Species","Phylum","Medicine","estimate","Cohort","antibiotic","std")

meta_gvhd_e<-read.csv( "F:/数据/document/学校/硕士/数据库/r数据库/gvhd/sankeyplot/sk_gvhd.csv",header=T)
meta_gvhd_e<-meta_gvhd_e[,-c(8,18)]
meta_gvhd_e<-melt(meta_gvhd_e,ID=c("Taxa","NodeColor"))
meta_gvhd_e<-na.omit(meta_gvhd_e)
meta_gvhd_e<-meta_gvhd_e %>% mutate(Cohort="GVHD") %>% mutate(antibioic=c(rep("no",27),rep("yes",66)))

meta_gvhd_s<-read.csv( "F:/数据/document/学校/硕士/数据库/r数据库/gvhd/sankeyplot/sk_gvhd_s.csv",header=T)
meta_gvhd_s<-meta_gvhd_s[,-c(8,18)]
meta_gvhd_s<-melt(meta_gvhd_s,ID=c("Taxa","NodeColor"))
meta_gvhd_s<-na.omit(meta_gvhd_s)

meta_gvhd<-cbind(meta_gvhd_e,meta_gvhd_s$value)
colnames(meta_gvhd)<-c("Species","Phylum","Medicine","estimate","Cohort","antibiotic","std")

meta_chemo_e<-read.csv( "F:/数据/document/学校/硕士/数据库/r数据库/gvhd/sankeyplot/sk_chemo.csv",header=T)
meta_chemo_e<-melt(meta_chemo_e,ID=c("Taxa","NodeColor"))
meta_chemo_e<-na.omit(meta_chemo_e)
meta_chemo_e<-meta_chemo_e %>% mutate(Cohort="Chemo") %>% mutate(antibioic=c(rep("yes",27)))

meta_chemo_s<-read.csv( "F:/数据/document/学校/硕士/数据库/r数据库/gvhd/sankeyplot/sk_chemo_s.csv",header=T)
meta_chemo_s<-melt(meta_chemo_s,ID=c("Taxa","NodeColor"))
meta_chemo_s<-na.omit(meta_chemo_s)

meta_chemo<-cbind(meta_chemo_e,meta_chemo_s$value)
colnames(meta_chemo)<-c("Species","Phylum","Medicine","estimate","Cohort","antibiotic","std")

meta_plt<-rbind(meta_hsct,meta_gvhd,meta_chemo)
meta_species<-c("V. parvula","R. mucilaginosa","B. caccae","B. ovatus","B. uniformis","C. lavalense",
                "C. spiroforme","R. torques","R. faecis","A. onderdonkii","C. citroniae","R. bromii")
meta_plt<-meta_plt |> filter(Species%in%meta_species) |> arrange(Species)

meta_plt$drug<-meta_plt$Medicine
meta_plt<-meta_plt |> group_by(Species,Medicine) |> #metagen(estimate,std) |> 
  summarise(
    size=length(drug),
    Phylum=Phylum,
    Cohort=Cohort,
    antibiotic=antibiotic,
    estimate=estimate,
    std=std
  )

meta_plt1<-meta_plt |> filter(size==1)
meta_plt2<-meta_plt |> filter(size==2)
meta_species2<-unique(meta_plt2$Species)
meta_plt2formerge<-array(0,dim<-c(length(meta_species2),ncol(meta_plt2)))
colnames(meta_plt2formerge)<-colnames(meta_plt2)
x=0
for (i in meta_species2){
  x=x+1
  tmp<-meta_plt2 |> filter(Species==i)
  meta_plt2formerge[x,1]<-i
  meta_plt2formerge[x,2]<-levels(tmp$Medicine)[tmp$Medicine[1]]
  meta_plt2formerge[x,3]<-tmp$size[1]
  meta_plt2formerge[x,4]<-tmp$Phylum[1]
  meta_plt2formerge[x,5]<-paste(tmp$Cohort[1],tmp$Cohort[2],sep="+")
  meta_plt2formerge[x,6]<-tmp$antibiotic[1]
  het<-metagen(tmp$estimate,tmp$std)
  meta_plt2formerge[x,7]<-het$TE.random
  meta_plt2formerge[x,8]<-het$seTE.random
}

meta_plt2formerge<-as.data.frame(meta_plt2formerge)
meta_plt2formerge$size<-as.numeric(meta_plt2formerge$size)
meta_plt2formerge$estimate<-as.numeric(meta_plt2formerge$estimate)
meta_plt2formerge$std<-as.numeric(meta_plt2formerge$std)
meta_pltmerge<-rbind(meta_plt1,meta_plt2formerge)


```

#### non-antibiotic

```{r}
meta_pltmerge_nonanti<-meta_pltmerge |> filter(antibiotic=="no")
meta_species_nonanti<-unique(meta_pltmerge_nonanti$Species)
meta_species_nonantif<-function(i){
  tmp<-meta_pltmerge_nonanti |> filter(Species==i)
  metaspecies<-array(NA,dim=c(nrow(tmp),14))
  colnames(metaspecies)<-c("Medicine","Number of associations","Cohort","Phylum","TE","lci","uci","k",
                           "pval","I2","Q","pval-h","Species","antibiotic")
  het<-metagen(tmp$estimate,tmp$std,tmp$Medicine)
  metaspecies[,1]<-summary(het)$studlab
  metaspecies[,2]<-tmp$size
  metaspecies[,3]<-tmp$Cohort
  metaspecies[,4]<-tmp$Phylum
  metaspecies[,5]<-summary(het)$TE
  metaspecies[,6]<-summary(het)$lower
  metaspecies[,7]<-summary(het)$upper
  metaspecies[,8]<-NA
  metaspecies[,9]<-summary(het)$pval
  metaspecies[,10]<-NA
  metaspecies[,11]<-NA
  metaspecies[,12]<-NA
  metaspecies[,13]<-i
  metaspecies[,14]<-tmp$antibiotic
  metaspecies<-metaspecies |> as.data.frame() |> arrange(Medicine)
  var<-c("IGV Model(non-anti)",summary(het)$k,"Meta",tmp$Phylum[1],het$TE.random,
         het$lower.random,het$upper.random,summary(het)$k,
         het$pval.random,het$I2,het$Q,pchisq(het$Q,df=2,lower.tail=F),i,"Meta")
  metaspecies<-rbind(metaspecies,var)
  metaspecies$Cohort<-factor(metaspecies$Cohort,levels = c("HSCT","GVHD","Chemo","HSCT+GVHD","HSCT+Chemo","Meta"))
  metaspecies$TE<-as.numeric(metaspecies$TE)
  metaspecies$lci<-as.numeric(metaspecies$lci)
  metaspecies$uci<-as.numeric(metaspecies$uci)
  metaspecies$`Number of associations`<-as.numeric(metaspecies$`Number of associations`)
  metaspecies$Medicine<-as.factor(metaspecies$Medicine)
  metaspecies$ynumber<-c(1:(nrow(tmp)+1))
  return(metaspecies)
  
}
species_meta_nonanti<-lapply(meta_species_nonanti,meta_species_nonantif)
species_meta_nonanti1<-adply(species_meta_nonanti,1)[,-1]
```

#### antibitic

```{r}
meta_pltmerge_anti<-meta_pltmerge |> filter(antibiotic=="yes")
meta_species_anti<-unique(meta_pltmerge_anti$Species)
meta_species_antif<-function(i){
  tmp<-meta_pltmerge_anti |> filter(Species==i)
  metaspecies<-array(NA,dim=c(nrow(tmp),14))
  colnames(metaspecies)<-c("Medicine","Number of associations","Cohort","Phylum","TE","lci","uci","k",
                           "pval","I2","Q","pval-h","Species","antibiotic")
  het<-metagen(tmp$estimate,tmp$std,tmp$Medicine)
  metaspecies[,1]<-summary(het)$studlab
  metaspecies[,2]<-tmp$size
  metaspecies[,3]<-tmp$Cohort
  metaspecies[,4]<-tmp$Phylum
  metaspecies[,5]<-summary(het)$TE
  metaspecies[,6]<-summary(het)$lower
  metaspecies[,7]<-summary(het)$upper
  metaspecies[,8]<-NA
  metaspecies[,9]<-summary(het)$pval
  metaspecies[,10]<-NA
  metaspecies[,11]<-NA
  metaspecies[,12]<-NA
  metaspecies[,13]<-i
  metaspecies[,14]<-tmp$antibiotic
  metaspecies<-metaspecies |> as.data.frame() |> arrange(Medicine)
  var<-c("IGV Model(anti)",summary(het)$k,"Meta",tmp$Phylum[1],het$TE.random,
         het$lower.random,het$upper.random,summary(het)$k,
         het$pval.random,het$I2,het$Q,pchisq(het$Q,df=2,lower.tail=F),i,"Meta")
  metaspecies<-rbind(metaspecies,var)
  metaspecies$Cohort<-factor(metaspecies$Cohort,levels = c("HSCT","GVHD","Chemo","HSCT+GVHD","HSCT+Chemo","Meta"))
  metaspecies$TE<-as.numeric(metaspecies$TE)
  metaspecies$lci<-as.numeric(metaspecies$lci)
  metaspecies$uci<-as.numeric(metaspecies$uci)
  metaspecies$`Number of associations`<-as.numeric(metaspecies$`Number of associations`)
  metaspecies$Medicine<-as.factor(metaspecies$Medicine)
  metaspecies$ynumber<-c(10:(nrow(tmp)+10))
  return(metaspecies)
  
}
species_meta_anti<-lapply(meta_species_anti,meta_species_antif)
species_meta_anti1<-adply(species_meta_anti,1)[,-1]
```

#### merge analysis

```{r}
species_meta2<-rbind(species_meta_nonanti1,species_meta_anti1)
species_meta2$ynumber<-as.numeric(species_meta2$ynumber)
species_meta2<-species_meta2 |> arrange(Species,ynumber)
x=0
metaspecies<-array(NA,dim=c(length(meta_species),14))
colnames(metaspecies)<-c("Medicine","Number of associations","Cohort","Phylum","TE","lci","uci","k",
                         "pval","I2","Q","pval-h","Species","antibiotic")
for (i in meta_species){
  x=x+1
  tmp<-meta_pltmerge |> filter(Species==i)
  het<-metagen(tmp$estimate,tmp$std,tmp$Medicine)
  metaspecies[x,1]<-"IGV Model"
  metaspecies[x,2]<-summary(het)$k
  metaspecies[x,3]<-"Meta"
  metaspecies[x,4]<-tmp$Phylum[1]
  metaspecies[x,5]<-het$TE.random
  metaspecies[x,6]<-het$lower.random
  metaspecies[x,7]<-het$upper.random
  metaspecies[x,8]<-summary(het)$k
  metaspecies[x,9]<-het$pval.random
  metaspecies[x,10]<-het$I2
  metaspecies[x,11]<-het$Q
  metaspecies[x,12]<-pchisq(het$Q,df=2,lower.tail=F)
  metaspecies[x,13]<-i
  metaspecies[x,14]<-"Meta"

}
metaspecies<-metaspecies |> as.data.frame() |> mutate(ynumber=rep(20,12))
metaspecies$TE<-as.numeric(metaspecies$TE)
metaspecies$lci<-as.numeric(metaspecies$lci)
metaspecies$uci<-as.numeric(metaspecies$uci)
metaspecies$k<-as.numeric(metaspecies$k)
metaspecies$`Number of associations`<-as.numeric(metaspecies$`Number of associations`)
metaspecies$Medicine<-as.factor(metaspecies$Medicine)

species_meta2<-rbind(species_meta_nonanti1,species_meta_anti1,metaspecies)
species_meta2$ynumber<-as.numeric(species_meta2$ynumber)
species_meta2<-species_meta2 |> arrange(Species,ynumber)
species_meta2<-species_meta2[-c(5,75,82),]
```

#### Figure 5sA

```{r}
species_meta_plt<- ggplot(species_meta2,aes(x=TE,y=reorder(Medicine,-ynumber),color=Cohort))+
  geom_errorbar(aes(xmin=lci,xmax=uci),color = "#666666",width=0.3) +
  geom_point(aes(shape=antibiotic,size=`Number of associations`))+
  scale_shape_manual(name="Antibiotic",
                     labels=levels(as.factor(species_meta2$antibiotic)),
                     values = c(no=17,yes=18,Meta=16))+
  scale_size_area(max_size = 8.5)+
  scale_color_manual(
    name="Cohort",
    labels=levels(species_meta2$Cohort),
    values=c(HSCT="#BC3C29",GVHD="#0072B5",Chemo="#008134",`HSCT+GVHD`="#6f57b4",`HSCT+Chemo`="#9e9f82",Meta="#241c1f"))+
  theme(panel.background = element_blank(),
        panel.grid =element_blank(),
        axis.line = element_line(colour = "grey50"),
        axis.title = element_text(size=15,family="serif",face="bold"),
        axis.text = element_text(size=12,family="serif",face="bold"),
        #legend.position = 'none',
        legend.title = element_text(size=16,family="serif",face="bold"),
        legend.text = element_text(size=15,family="serif",face="bold"),
        plot.subtitle = element_text(size=12,family="serif",face="bold"),
        strip.text = element_text(size=12,family="serif",face="bold"),
        strip.background = element_rect(fill = NA),
        title = element_text(size=15,family="serif",face="bold")
  )+
  labs(x="Regression coefficient",y="Medicine")+
  geom_vline(xintercept = 0,linetype="dashed")+
  facet_wrap(~Species,scales = "free")+
  guides(colour = guide_legend(override.aes = list(size=5)),
         shape=guide_legend(override.aes = list(size=5)))#修改图例点大小

species_meta_plt
```

### Multiple drug analysis

```{r}
meta_mul_hsct_e<-read.csv( "F:/数据/document/学校/硕士/数据库/r数据库/gvhd/sankeyplot/sk_mul_hsct.csv",header=T)
meta_mul_hsct_e<-meta_mul_hsct_e[,-6]
meta_mul_hsct_e<-melt(meta_mul_hsct_e,ID=c("Taxa","NodeColor"))
meta_mul_hsct_e<-na.omit(meta_mul_hsct_e)
meta_mul_hsct_e<-meta_mul_hsct_e %>% mutate(Cohort="HSCT") %>% mutate(antibioic=c(rep("no",23),rep("yes",56)))

meta_mul_hsct_s<-read.csv( "F:/数据/document/学校/硕士/数据库/r数据库/gvhd/sankeyplot/sk_mul_hsct_s.csv",header=T)
meta_mul_hsct_s<-meta_mul_hsct_s[,-6]
meta_mul_hsct_s<-melt(meta_mul_hsct_s,ID=c("Taxa","NodeColor"))
meta_mul_hsct_s<-na.omit(meta_mul_hsct_s)

meta_mul_hsct<-cbind(meta_mul_hsct_e,meta_mul_hsct_s$value)
colnames(meta_mul_hsct)<-c("Species","Phylum","Medicine","estimate","Cohort","antibiotic","std")

meta_mul_gvhd_e<-read.csv( "F:/数据/document/学校/硕士/数据库/r数据库/gvhd/sankeyplot/sk_mul_gvhd.csv",header=T)
meta_mul_gvhd_e<-meta_mul_gvhd_e[,-c(8,17)]
meta_mul_gvhd_e<-melt(meta_mul_gvhd_e,ID=c("Taxa","NodeColor"))
meta_mul_gvhd_e<-na.omit(meta_mul_gvhd_e)
meta_mul_gvhd_e<-meta_mul_gvhd_e %>% mutate(Cohort="GVHD") %>% mutate(antibioic=c(rep("no",32),rep("yes",60)))

meta_mul_gvhd_s<-read.csv( "F:/数据/document/学校/硕士/数据库/r数据库/gvhd/sankeyplot/sk_mul_gvhd_s.csv",header=T)
meta_mul_gvhd_s<-meta_mul_gvhd_s[,-c(8,17)]
meta_mul_gvhd_s<-melt(meta_mul_gvhd_s,ID=c("Taxa","NodeColor"))
meta_mul_gvhd_s<-na.omit(meta_mul_gvhd_s)

meta_mul_gvhd<-cbind(meta_mul_gvhd_e,meta_mul_gvhd_s$value)
colnames(meta_mul_gvhd)<-c("Species","Phylum","Medicine","estimate","Cohort","antibiotic","std")

meta_mul_chemo_e<-read.csv( "F:/数据/document/学校/硕士/数据库/r数据库/gvhd/sankeyplot/sk_mul_chemo.csv",header=T)
meta_mul_chemo_e<-melt(meta_mul_chemo_e,ID=c("Taxa","NodeColor"))
meta_mul_chemo_e<-na.omit(meta_mul_chemo_e)
meta_mul_chemo_e<-meta_mul_chemo_e %>% mutate(Cohort="Chemo") %>% mutate(antibioic=c(rep("no",5)))

meta_mul_chemo_s<-read.csv( "F:/数据/document/学校/硕士/数据库/r数据库/gvhd/sankeyplot/sk_mul_chemo_s.csv",header=T)
meta_mul_chemo_s<-melt(meta_mul_chemo_s,ID=c("Taxa","NodeColor"))
meta_mul_chemo_s<-na.omit(meta_mul_chemo_s)

meta_mul_chemo<-cbind(meta_mul_chemo_e,meta_mul_chemo_s$value)
colnames(meta_mul_chemo)<-c("Species","Phylum","Medicine","estimate","Cohort","antibiotic","std")

meta_mul_plt<-rbind(meta_mul_hsct,meta_mul_gvhd,meta_mul_chemo)
meta_mul_species<-c("V. parvula","E. dolichum","C. methylpentosum","B. producta")
meta_mul_plt<-meta_mul_plt |> filter(Species%in%meta_mul_species) |> arrange(Species)

meta_mul_plt$drug<-meta_mul_plt$Medicine
meta_mul_plt<-meta_mul_plt |> group_by(Species,Medicine) |> 
  summarise(
    size=length(drug),
    Phylum=Phylum,
    Cohort=Cohort,
    antibiotic=antibiotic,
    estimate=estimate,
    std=std
  )

meta_mul_plt1<-meta_mul_plt |> filter(size==1)
meta_mul_plt2<-meta_mul_plt |> filter(size==2)
meta_mul_species2<-unique(meta_mul_plt2$Species)
meta_mul_plt2formerge<-array(0,dim<-c(length(meta_mul_species2),ncol(meta_mul_plt2)))
colnames(meta_mul_plt2formerge)<-colnames(meta_mul_plt2)
x=0
for (i in meta_mul_species2){
  x=x+1
  tmp<-meta_mul_plt2 |> filter(Species==i)
  meta_mul_plt2formerge[x,1]<-i
  meta_mul_plt2formerge[x,2]<-levels(tmp$Medicine)[tmp$Medicine[1]]
  meta_mul_plt2formerge[x,3]<-tmp$size[1]
  meta_mul_plt2formerge[x,4]<-tmp$Phylum[1]
  meta_mul_plt2formerge[x,5]<-paste(tmp$Cohort[1],tmp$Cohort[2],sep="+")
  meta_mul_plt2formerge[x,6]<-tmp$antibiotic[1]
  het<-metagen(tmp$estimate,tmp$std)
  meta_mul_plt2formerge[x,7]<-het$TE.random
  meta_mul_plt2formerge[x,8]<-het$seTE.random
}

meta_mul_plt2formerge<-as.data.frame(meta_mul_plt2formerge)
meta_mul_plt2formerge$size<-as.numeric(meta_mul_plt2formerge$size)
meta_mul_plt2formerge$estimate<-as.numeric(meta_mul_plt2formerge$estimate)
meta_mul_plt2formerge$std<-as.numeric(meta_mul_plt2formerge$std)
meta_mul_pltmerge<-rbind(meta_mul_plt1,meta_mul_plt2formerge)
```

#### non-antibiotic

```{r}
meta_mul_pltmerge_nonanti<-meta_mul_pltmerge |> filter(antibiotic=="no")
meta_mul_species_nonanti<-unique(meta_mul_pltmerge_nonanti$Species)
meta_mul_species_nonantif<-function(i){
  tmp<-meta_mul_pltmerge_nonanti |> filter(Species==i)
  metaspecies<-array(NA,dim=c(nrow(tmp),14))
  colnames(metaspecies)<-c("Medicine","Number of associations","Cohort","Phylum","TE","lci","uci","k",
                           "pval","I2","Q","pval-h","Species","antibiotic")
  het<-metagen(tmp$estimate,tmp$std,tmp$Medicine)
  metaspecies[,1]<-summary(het)$studlab
  metaspecies[,2]<-tmp$size
  metaspecies[,3]<-tmp$Cohort
  metaspecies[,4]<-tmp$Phylum
  metaspecies[,5]<-summary(het)$TE
  metaspecies[,6]<-summary(het)$lower
  metaspecies[,7]<-summary(het)$upper
  metaspecies[,8]<-NA
  metaspecies[,9]<-summary(het)$pval
  metaspecies[,10]<-NA
  metaspecies[,11]<-NA
  metaspecies[,12]<-NA
  metaspecies[,13]<-i
  metaspecies[,14]<-tmp$antibiotic
  metaspecies<-metaspecies |> as.data.frame() |> arrange(Medicine)
  var<-c("IGV Model(non-anti)",summary(het)$k,"Meta",tmp$Phylum[1],het$TE.random,
         het$lower.random,het$upper.random,summary(het)$k,
         het$pval.random,het$I2,het$Q,pchisq(het$Q,df=2,lower.tail=F),i,"Meta")
  metaspecies<-rbind(metaspecies,var)
  metaspecies$Cohort<-factor(metaspecies$Cohort,levels = c("HSCT","GVHD","Chemo","HSCT+GVHD","HSCT+Chemo","Meta"))
  metaspecies$TE<-as.numeric(metaspecies$TE)
  metaspecies$lci<-as.numeric(metaspecies$lci)
  metaspecies$uci<-as.numeric(metaspecies$uci)
  metaspecies$`Number of associations`<-as.numeric(metaspecies$`Number of associations`)
  metaspecies$Medicine<-as.factor(metaspecies$Medicine)
  metaspecies$ynumber<-c(1:(nrow(tmp)+1))
  return(metaspecies)
  
}
species_meta_mul_nonanti<-lapply(meta_mul_species_nonanti,meta_mul_species_nonantif)
species_meta_mul_nonanti1<-adply(species_meta_mul_nonanti,1)[,-1]


```

#### antibitic

```{r}
meta_mul_pltmerge_anti<-meta_mul_pltmerge |> filter(antibiotic=="yes")
meta_mul_species_anti<-unique(meta_mul_pltmerge_anti$Species)
meta_mul_species_antif<-function(i){
  tmp<-meta_mul_pltmerge_anti |> filter(Species==i)
  metaspecies<-array(NA,dim=c(nrow(tmp),14))
  colnames(metaspecies)<-c("Medicine","Number of associations","Cohort","Phylum","TE","lci","uci","k",
                           "pval","I2","Q","pval-h","Species","antibiotic")
  het<-metagen(tmp$estimate,tmp$std,tmp$Medicine)
  metaspecies[,1]<-summary(het)$studlab
  metaspecies[,2]<-tmp$size
  metaspecies[,3]<-tmp$Cohort
  metaspecies[,4]<-tmp$Phylum
  metaspecies[,5]<-summary(het)$TE
  metaspecies[,6]<-summary(het)$lower
  metaspecies[,7]<-summary(het)$upper
  metaspecies[,8]<-NA
  metaspecies[,9]<-summary(het)$pval
  metaspecies[,10]<-NA
  metaspecies[,11]<-NA
  metaspecies[,12]<-NA
  metaspecies[,13]<-i
  metaspecies[,14]<-tmp$antibiotic
  metaspecies<-metaspecies |> as.data.frame() |> arrange(Medicine)
  var<-c("IGV Model(anti)",summary(het)$k,"Meta",tmp$Phylum[1],het$TE.random,
         het$lower.random,het$upper.random,summary(het)$k,
         het$pval.random,het$I2,het$Q,pchisq(het$Q,df=2,lower.tail=F),i,"Meta")
  metaspecies<-rbind(metaspecies,var)
  metaspecies$Cohort<-factor(metaspecies$Cohort,levels = c("HSCT","GVHD","Chemo","HSCT+GVHD","HSCT+Chemo","Meta"))
  metaspecies$TE<-as.numeric(metaspecies$TE)
  metaspecies$lci<-as.numeric(metaspecies$lci)
  metaspecies$uci<-as.numeric(metaspecies$uci)
  metaspecies$`Number of associations`<-as.numeric(metaspecies$`Number of associations`)
  metaspecies$Medicine<-as.factor(metaspecies$Medicine)
  metaspecies$ynumber<-c(10:(nrow(tmp)+10))
  return(metaspecies)
  
}
species_meta_mul_anti<-lapply(meta_mul_species_anti,meta_mul_species_antif)
species_meta_mul_anti1<-adply(species_meta_mul_anti,1)[,-1]

```

#### merge analysis

```{r}
x=0
metaspecies_mul<-array(NA,dim=c(length(meta_mul_species),14))
colnames(metaspecies_mul)<-c("Medicine","Number of associations","Cohort","Phylum","TE","lci","uci","k",
                         "pval","I2","Q","pval-h","Species","antibiotic")
for (i in meta_mul_species){
  x=x+1
  tmp<-meta_mul_pltmerge |> filter(Species==i)
  het<-metagen(tmp$estimate,tmp$std,tmp$Medicine)
  metaspecies_mul[x,1]<-"IGV Model"
  metaspecies_mul[x,2]<-summary(het)$k
  metaspecies_mul[x,3]<-"Meta"
  metaspecies_mul[x,4]<-tmp$Phylum[1]
  metaspecies_mul[x,5]<-het$TE.random
  metaspecies_mul[x,6]<-het$lower.random
  metaspecies_mul[x,7]<-het$upper.random
  metaspecies_mul[x,8]<-summary(het)$k
  metaspecies_mul[x,9]<-het$pval.random
  metaspecies_mul[x,10]<-het$I2
  metaspecies_mul[x,11]<-het$Q
  metaspecies_mul[x,12]<-pchisq(het$Q,df=2,lower.tail=F)
  metaspecies_mul[x,13]<-i
  metaspecies_mul[x,14]<-"Meta"
  
}
metaspecies_mul<-metaspecies_mul |> as.data.frame() |> mutate(ynumber=rep(20,4))
metaspecies_mul$TE<-as.numeric(metaspecies_mul$TE)
metaspecies_mul$lci<-as.numeric(metaspecies_mul$lci)
metaspecies_mul$uci<-as.numeric(metaspecies_mul$uci)
metaspecies_mul$k<-as.numeric(metaspecies_mul$k)
metaspecies_mul$`Number of associations`<-as.numeric(metaspecies_mul$`Number of associations`)
metaspecies_mul$Medicine<-as.factor(metaspecies_mul$Medicine)

species_meta_mul2<-rbind(species_meta_mul_nonanti1,species_meta_mul_anti1,metaspecies_mul)
species_meta_mul2$ynumber<-as.numeric(species_meta_mul2$ynumber)
species_meta_mul2<-species_meta_mul2 |> arrange(Species,ynumber)
```

#### Figure 5sB

```{r}
species_meta_mul_plt<- ggplot(species_meta_mul2,aes(x=TE,y=reorder(Medicine,-ynumber),color=Cohort))+
  geom_errorbar(aes(xmin=lci,xmax=uci),color = "#666666",width=0.3) +
  geom_point(aes(shape=antibiotic,size=`Number of associations`))+
  scale_shape_manual(name="Antibiotic",
                     labels=levels(as.factor(species_meta2$antibiotic)),
                     values = c(no=17,yes=18,Meta=16))+
  scale_size_area(max_size = 7)+
  scale_color_manual(
    name="Cohort",
    labels=levels(species_meta2$Cohort),
    values=c(HSCT="#BC3C29",GVHD="#0072B5",Chemo="#008134",`HSCT+GVHD`="#6f57b4",`HSCT+Chemo`="#9e9f82",Meta="#241c1f"))+
  theme(panel.background = element_blank(),
        panel.grid =element_blank(),
        axis.line = element_line(colour = "grey50"),
        axis.title = element_text(size=15,family="serif",face="bold"),
        axis.text = element_text(size=12,family="serif",face="bold"),
        #legend.position = 'none',
        legend.title = element_text(size=16,family="serif",face="bold"),
        legend.text = element_text(size=15,family="serif",face="bold"),
        plot.subtitle = element_text(size=12,family="serif",face="bold"),
        strip.text = element_text(size=12,family="serif",face="bold"),
        strip.background = element_rect(fill = NA),
        title = element_text(size=15,family="serif",face="bold")
  )+
  labs(x="Regression coefficient",y="Medicine")+
  geom_vline(xintercept = 0,linetype="dashed")+
  facet_wrap(~Species,scales = "free")+
  guides(colour = guide_legend(override.aes = list(size=5)),
         shape=guide_legend(override.aes = list(size=5)))#修改图例点大小

species_meta_mul_plt
```

## Figure 6

```{r}
library(ade4)
library(ggtree)

name_species<-list(hp_busulfan,hp_cyclophosphamide,hp_amoxicillin,hp_antibiotic,hp_ceftazidime,hp_Ciprofloxacin,
                   hp_fludarabine,hp_melphalan,hp_meropenem,hp_Methotrexate,hp_metronidazole,hp_phenoxymethylpenicillin,
                   hp_piperacillin,hp_Thiotepa,hp_trimethoprim,hp_vancomycin,hp_VP16)
hp_species<-Reduce(full_join,name_species)
hp_species$Group<-"Single drug analysis"
length(setdiff(test_allcohort_1$Meds,hp_species$Drug))
hp_species_bu<-data.frame(
  Taxa=rep(NA,8),
  Drug=setdiff(test_allcohort_1$Meds,hp_species$Drug),
  beta=rep(NA,8),
  FDR=rep(NA,8),
  Influence=rep(NA,8),
  anti=c("yes","no",rep("yes",6)),
  Group=rep("Single drug analysis",8)
)
hp_species<-rbind(hp_species,hp_species_bu)
# hp_species$Drug<-factor(hp_species$Drug,levels = levels(test_merge_cohort$Meds))

name_mul_species<-list(hp_mul_busulfan,hp_mul_cyclophosphamide,hp_mul_amoxicillin,hp_mul_ceftazidime,hp_mul_Ciprofloxacin,
                       hp_mul_fludarabine,hp_mul_meropenem,hp_mul_Methotrexate,hp_mul_phenoxymethylpenicillin,
                       hp_mul_Thiotepa,hp_mul_trimethoprim,hp_mul_vancomycin,hp_mul_VP16)

hp_mul_species<-Reduce(full_join,name_mul_species)
hp_mul_species$Group<-"Multiple drug analysis"
length(setdiff(test_mul_allcohort_1$Meds,hp_mul_species$Drug))
hp_mul_species_bu<-data.frame(
  Taxa=rep(NA,9),
  Drug=setdiff(test_mul_allcohort_1$Meds,hp_mul_species$Drug),
  beta_mul=rep(NA,9),
  FDR_mul=rep(NA,9),
  Influence=rep(NA,9),
  anti=c("no","yes","yes","no","yes","no","yes","yes","yes"),
  Group=rep("Multiple drug analysis",9)
)
hp_mul_species<-rbind(hp_mul_species,hp_mul_species_bu)
# hp_mul_species$Drug<-factor(hp_mul_species$Drug,levels = levels(test_merge_cohort$Meds))

hp_mul_species_bu1<-data.frame(
  Taxa=unique(setdiff(hp_species$Taxa,hp_mul_species$Taxa)),
  Drug=rep(NA,17),
  beta_mul=rep(NA,17),
  FDR_mul=rep(NA,17),
  Influence=rep(NA,17),
  anti=rep(NA,17),
  Group=rep("Multiple drug analysis",17)
)
hp_mul_species<-rbind(hp_mul_species,hp_mul_species_bu1)
colnames(hp_mul_species)[3:4]<-c("beta","FDR")

hp_merge_species<-rbind(hp_species,hp_mul_species)
write.csv(hp_merge_species,"hp_merge_species.csv")

tree_merge_species<-read.csv("tree_merge_species.csv")
hp_merge_species<-read.csv("hp_merge_species_mod.csv")
tree_merge_species$ID<-factor(tree_merge_species$ID,levels = tree_merge_species$ID)



rownames(tree_merge_species)<-tree_merge_species$ID
tree_merge_species<-tree_merge_species[,-1]
tree_merge_species$p <-as.factor(tree_merge_species$p)
tree_merge_species$c <-as.factor(tree_merge_species$c)
tree_merge_species$o <-as.factor(tree_merge_species$o)
tree_merge_species$f <-as.factor(tree_merge_species$f)
tree_merge_species$g <-as.factor(tree_merge_species$g)
tmptax.phy <- taxo2phylog(as.taxo(tree_merge_species))
print(tmptax.phy)
tree_plt_name<-ggtree(tmptax.phy,colour = "grey50")  + 
  geom_tiplab()+
  # geom_nodelab()+
  hexpand(.25)#show id out of edge

tree_plt_noname<-ggtree(tmptax.phy,colour = "grey50",size=0.7)  + 
  coord_flip()

hp_merge_species$Drug<-factor(hp_merge_species$Drug,levels = levels(test_merge_cohort$Meds))
hp_merge_species$Group<-factor(hp_merge_species$Group,levels = c("Single drug analysis","Multiple drug analysis"))
hp_merge_species$ID<-factor(hp_merge_species$ID,levels = c("f__Enterobacteriaceae","g__Actinomyces","B. breve","g__Parabacteroides","g__Bacteroides","B. ovatus","B. fragilis","B. caccae","C. ramosum","E. dolichum","o__Lactobacillales","g__Streptococcus","g__Enterococcus1","g__Enterococcus2","g__Lactobacillus","L. zeae","o__Clostridiales1","o__Clostridiales2","g__Clostridium","f__Peptostreptococcaceae","V. parvula","V. dispar","f__Ruminococcaceae","F. prausnitzii","g__Oscillospira","g__Ruminococcus1","S. variabile","f__Lachnospiraceae1","f__Lachnospiraceae2","C. hathewayi","g__Coprococcus","g__Dorea","g__Blautia","B. producta","g__Roseburia","R. faecis","R. gnavus","g__Ruminococcus2"))
hp_merge_species$ID[c(128:135,225:233)]<-"g__Ruminococcus2"
hp_merge_species$Drug[234:250]<-"Tobramycin"
temp<-ggplot(hp_merge_species,aes(x=ID,y=Drug,color=Influence))+
  geom_point(aes(size=abs(beta),alpha=-FDR))+
  scale_size_area(max_size = 8)+
  # scale_size_continuous(name="Magnitude")+
  scale_color_manual(name="Influence",
                     labels=c("Inhibition","Promotion"),
                     values=c("#8491B4FF","#F39B7FFF")) +
  scale_alpha_continuous(name="Adjusted P value (FDR)",
                         limits = c(-1,0), breaks = c(0, -0.01, -0.1,-0.5, -1),
                         guide = guide_colourbar(nbin = 100, draw.ulim = FALSE, draw.llim = FALSE))+
  theme(panel.background = element_blank(),
        panel.grid =element_blank(),
        axis.line = element_blank(),
        # axis.title.y = element_text(size=14,family="serif",face="bold"),
        axis.title = element_blank(),
        axis.text.y =  element_text(size=10,family="serif",face="bold"),
        axis.text.x =  element_text(size=12,family="serif",face="bold",angle=90,hjust = 1,vjust = 0.5),
        #legend.position = 'none',
        legend.title = element_text(size=16,family="serif",face="bold"),
        legend.text = element_text(size=15,family="serif",face="bold"),
        plot.subtitle = element_text(size=12,family="serif",face="bold"),
        strip.text = element_text(size=12,family="serif",face="bold"),
        strip.background = element_rect(fill = NA),
        title = element_text(size=14,family="serif",face="bold")
  )+
  labs(size="Magnitude")+
  facet_grid(rows = vars(Group), drop=T, scales = "free")+
  guides(colour = guide_legend(override.aes = list(size=5)),
         alpha=guide_legend(override.aes = list(size=5)))#修改图例点大小

temp
b<-temp/tree_plt_noname+plot_layout(heights = c(4.2, 0.8))
ggsave("Figure/Figure5_single.pdf",b,units="in", width=16.0, height=10.0, dpi=600,limitsize = FALSE)
```
