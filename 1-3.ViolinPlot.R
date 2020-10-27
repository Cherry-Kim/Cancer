#!/usr/bin/env Rscript
library(ggplot2)
library(ggsignif)
library(tidyverse)
library(ggpubr)
library(stats)
library(gridExtra)

png("ViolinPlot.png", width=3000, height=2000, res=250)

data <- read.csv("BTG2_cms3.txt", header=T,stringsAsFactors=F,sep="\t")
h = ggplot(data, aes(x=CMS, y=TPM, fill=GENE))+geom_violin(trim=FALSE)+geom_jitter(shape=16, position=position_jitter(0.2))+scale_fill_hue(l=50)+coord_cartesian(ylim=c(0,50))
a = h + stat_compare_means(method = "anova",   label.x = 1.5, label.y = 45, label="p.signif")

data2 <- read.csv("IKZF1_cms3.txt", header=T,stringsAsFactors=F,sep="\t")
h = ggplot(data2, aes(x=CMS, y=TPM, fill=GENE))+geom_violin(trim=FALSE)+geom_jitter(shape=16, position=position_jitter(0.2))+scale_fill_hue(l=50)+coord_cartesian(ylim=c(0,50))
b = h + stat_compare_means(method = "anova",   label.x = 1.5, label.y = 45, label="p.signif")

data3 <- read.csv("IRF8_cms3.txt", header=T,stringsAsFactors=F,sep="\t")
h = ggplot(data3, aes(x=CMS, y=TPM, fill=GENE))+geom_violin(trim=FALSE)+geom_jitter(shape=16, position=position_jitter(0.2))+scale_fill_hue(l=50)+coord_cartesian(ylim=c(0,50))
c = h + stat_compare_means(method = "anova",   label.x = 1.5, label.y = 45, label="p.signif")

data4 <- read.csv("NFE2_cms3.txt", header=T,stringsAsFactors=F,sep="\t")
h = ggplot(data4, aes(x=CMS, y=TPM, fill=GENE))+geom_violin(trim=FALSE)+geom_jitter(shape=16, position=position_jitter(0.2))+scale_fill_hue(l=50)+coord_cartesian(ylim=c(0,50))
d = h + stat_compare_means(method = "anova",   label.x = 1.5, label.y = 45, label="p.signif")

data5 <- read.csv("SPIB_cms3.txt", header=T,stringsAsFactors=F,sep="\t")
h = ggplot(data5, aes(x=CMS, y=TPM, fill=GENE))+geom_violin(trim=FALSE)+geom_jitter(shape=16, position=position_jitter(0.2))+scale_fill_hue(l=50)+coord_cartesian(ylim=c(0,50))
e = h + stat_compare_means(method = "anova",   label.x = 1.5, label.y = 45, label="p.signif")

data6 <- read.csv("EAF2_cms3.txt", header=T,stringsAsFactors=F,sep="\t")
h = ggplot(data6, aes(x=CMS, y=TPM, fill=GENE))+geom_violin(trim=FALSE)+geom_jitter(shape=16, position=position_jitter(0.2))+scale_fill_hue(l=50)+coord_cartesian(ylim=c(0,50))
f = h + stat_compare_means(method = "anova",   label.x = 1.5, label.y = 45, label="p.signif")

data7 <- read.csv("TAL1_cms3.txt", header=T,stringsAsFactors=F,sep="\t")
h = ggplot(data7, aes(x=CMS, y=TPM, fill=GENE))+geom_violin(trim=FALSE)+geom_jitter(shape=16, position=position_jitter(0.2))+scale_fill_hue(l=50)+coord_cartesian(ylim=c(0,50))
g = h + stat_compare_means(method = "anova",   label.x = 1.5, label.y = 45, label="p.signif")

data8 <- read.csv("GATA1_cms3.txt", header=T,stringsAsFactors=F,sep="\t")
h = ggplot(data8, aes(x=CMS, y=TPM, fill=GENE))+geom_violin(trim=FALSE)+geom_jitter(shape=16, position=position_jitter(0.2))+scale_fill_hue(l=50)+coord_cartesian(ylim=c(0,50))
i = h + stat_compare_means(method = "anova",   label.x = 1.5, label.y = 45, label="p.signif")

data9 <- read.csv("AKNA_cms3.txt", header=T,stringsAsFactors=F,sep="\t")
h = ggplot(data9, aes(x=CMS, y=TPM, fill=GENE))+geom_violin(trim=FALSE)+geom_jitter(shape=16, position=position_jitter(0.2))+scale_fill_hue(l=50)+coord_cartesian(ylim=c(0,50))
j = h + stat_compare_means(method = "anova",   label.x = 1.5, label.y = 45, label="p.signif")

data10 <- read.csv("POU2F2_cms3.txt", header=T,stringsAsFactors=F,sep="\t")
h = ggplot(data10, aes(x=CMS, y=TPM, fill=GENE))+geom_violin(trim=FALSE)+geom_jitter(shape=16, position=position_jitter(0.2))+scale_fill_hue(l=50)+coord_cartesian(ylim=c(0,50))
k = h + stat_compare_means(method = "anova",   label.x = 1.5, label.y = 45, label="p.signif")

data11 <- read.csv("IKZF3_cms3.txt", header=T,stringsAsFactors=F,sep="\t")
h = ggplot(data11, aes(x=CMS, y=TPM, fill=GENE))+geom_violin(trim=FALSE)+geom_jitter(shape=16, position=position_jitter(0.2))+scale_fill_hue(l=50)+coord_cartesian(ylim=c(0,50))
l = h + stat_compare_means(method = "anova",   label.x = 1.5, label.y = 45, label="p.signif")

data12 <- read.csv("KLF1_cms3.txt", header=T,stringsAsFactors=F,sep="\t")
h = ggplot(data12, aes(x=CMS, y=TPM, fill=GENE))+geom_violin(trim=FALSE)+geom_jitter(shape=16, position=position_jitter(0.2))+scale_fill_hue(l=50)+coord_cartesian(ylim=c(0,50))
m = h + stat_compare_means(method = "anova",   label.x = 1.5, label.y = 45, label="p.signif")

data13 <- read.csv("BCL11A_cms3.txt", header=T,stringsAsFactors=F,sep="\t")
h = ggplot(data13, aes(x=CMS, y=TPM, fill=GENE))+geom_violin(trim=FALSE)+geom_jitter(shape=16, position=position_jitter(0.2))+scale_fill_hue(l=50)+coord_cartesian(ylim=c(0,50))
n = h + stat_compare_means(method = "anova",   label.x = 1.5, label.y = 45, label="p.signif")

data14 <- read.csv("GFI1B_cms3.txt", header=T,stringsAsFactors=F,sep="\t")
h = ggplot(data14, aes(x=CMS, y=TPM, fill=GENE))+geom_violin(trim=FALSE)+geom_jitter(shape=16, position=position_jitter(0.2))+scale_fill_hue(l=50)+coord_cartesian(ylim=c(0,50))
o = h + stat_compare_means(method = "anova",   label.x = 1.5, label.y = 45, label="p.signif")

data15 <- read.csv("IRF5_cms3.txt", header=T,stringsAsFactors=F,sep="\t")
h = ggplot(data15, aes(x=CMS, y=TPM, fill=GENE))+geom_violin(trim=FALSE)+geom_jitter(shape=16, position=position_jitter(0.2))+scale_fill_hue(l=50)+coord_cartesian(ylim=c(0,50))
p = h + stat_compare_means(method = "anova",   label.x = 1.5, label.y = 45, label="p.signif")

data16 <- read.csv("JUNB_cms3.txt", header=T,stringsAsFactors=F,sep="\t")
h = ggplot(data16, aes(x=CMS, y=TPM, fill=GENE))+geom_violin(trim=FALSE)+geom_jitter(shape=16, position=position_jitter(0.2))+scale_fill_hue(l=50)+coord_cartesian(ylim=c(0,50))
q = h + stat_compare_means(method = "anova",   label.x = 1.5, label.y = 45, label="p.signif")

data17 <- read.csv("KLF2_cms3.txt", header=T,stringsAsFactors=F,sep="\t")
h = ggplot(data17, aes(x=CMS, y=TPM, fill=GENE))+geom_violin(trim=FALSE)+geom_jitter(shape=16, position=position_jitter(0.2))+scale_fill_hue(l=50)+coord_cartesian(ylim=c(0,50))
r = h + stat_compare_means(method = "anova",   label.x = 1.5, label.y = 45, label="p.signif")

data18 <- read.csv("MYB_cms3.txt", header=T,stringsAsFactors=F,sep="\t")
h = ggplot(data18, aes(x=CMS, y=TPM, fill=GENE))+geom_violin(trim=FALSE)+geom_jitter(shape=16, position=position_jitter(0.2))+scale_fill_hue(l=50)+coord_cartesian(ylim=c(0,50))
s = h + stat_compare_means(method = "anova",   label.x = 1.5, label.y = 45, label="p.signif")

data19 <- read.csv("TCF7_cms3.txt", header=T,stringsAsFactors=F,sep="\t")
h = ggplot(data19, aes(x=CMS, y=TPM, fill=GENE))+geom_violin(trim=FALSE)+geom_jitter(shape=16, position=position_jitter(0.2))+scale_fill_hue(l=50)+coord_cartesian(ylim=c(0,50))
t = h + stat_compare_means(method = "anova",   label.x = 1.5, label.y = 45, label="p.signif")

data20 <- read.csv("SPI1_cms3.txt", header=T,stringsAsFactors=F,sep="\t")
h = ggplot(data20, aes(x=CMS, y=TPM, fill=GENE))+geom_violin(trim=FALSE)+geom_jitter(shape=16, position=position_jitter(0.2))+scale_fill_hue(l=50)+coord_cartesian(ylim=c(0,50))
u = h + stat_compare_means(method = "anova",   label.x = 1.5, label.y = 45, label="p.signif")

grid.arrange(a,b,c,d,e,f,g,i,j,k,l,m,n,o,p,q,r,s,t,u, nrow=5, ncol=4)
dev.off()
