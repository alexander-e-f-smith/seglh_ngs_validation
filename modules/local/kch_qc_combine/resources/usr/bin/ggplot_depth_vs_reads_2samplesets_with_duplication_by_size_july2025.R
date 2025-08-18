#!/usr/bin/env Rscript

### by Alex Smith. Contact kch-tr.KCHBioinformatics@nhs.net
options(warn = 1)

args <-commandArgs(trailingOnly = TRUE)
dataset1 <- args[1]
name1 <- args[2]
name2 <- args[3]
assay_and_sampleskey <- args[4]
spanin <- args[5]
depth <- args[6]
reads_required <- args[7]
spannumber <-as.numeric(spanin)
reads_required_number <-as.numeric(reads_required)
depthfield <-paste0("COVERED_",depth,"X_PCT")
depthfield1 <-paste0("obs_COVERED_",depth,"X_PCT")
depthfield2 <-paste0("exp_COVERED_",depth,"X_PCT")
print(depthfield)
title<-paste("Depth of Coverage vs Sequence Read Number;\n",assay_and_sampleskey)
print(title)
output<-paste0(assay_and_sampleskey,"depth vs reads",depthfield,"_ggplot.pdf")

file1<-read.delim(dataset1)

library(ggplot2)
library(scales)

df1 = file1[,c('obs_Total_processed_paired_reads',depthfield1,'obs_Duplication_rate')]
df2 = file1[,c('exp_Total_processed_paired_reads',depthfield2,'exp_Duplication_rate')]

gg1<-ggplot(data=df1, aes(obs_Total_processed_paired_reads, get(depthfield1))) + 
    geom_point(aes(size=obs_Duplication_rate), colour="red" )  +
    scale_size(range = c(2, 3)) + 
    geom_smooth(data=df1, aes(obs_Total_processed_paired_reads, get(depthfield1), colour=paste(name1,"\n",depthfield1,"vs Processed sequence read number","\n(","loess smoothing=",spanin,")\n")), method = "loess", formula = 'y ~ x', se=TRUE, linetype="dashed", span=spannumber, fullrange=TRUE, level=0.9, linewidth=0.3) +
    geom_point(data=df2, aes(exp_Total_processed_paired_reads, get(depthfield2), size=exp_Duplication_rate), colour="blue") + guides(size=guide_legend(override.aes=list(colour="grey30"))) +
    geom_smooth(data=df2, aes(exp_Total_processed_paired_reads, get(depthfield2), colour=paste(name2,"\n",depthfield2,"vs Processed sequence read number","\n(","loess smoothing=",spanin,")\n")), method = "loess", formula = 'y ~ x', se=TRUE, linetype="dashed", span=spannumber, fullrange=TRUE, level=0.9, linewidth=0.3)  +
    scale_colour_manual(name="Data key:", values=c("red", "blue")) + 
    labs(x = "Processed sequence reads  (pairs)", y = paste("Proportion of ROI covered >","depth X","\n(filtered ROI)" )) + 
    scale_y_continuous(limits=c(0.7,1), breaks=seq(from = 0.8, to = 1, by = 0.02))  + scale_x_continuous(breaks=pretty_breaks(n=6)) + 
    ggtitle(title) + 
    geom_hline(yintercept = 0.95, linetype="dotted", linewidth=0.3) + 
    geom_vline(xintercept = reads_required_number, linetype="dotted", linewidth=0.3)

pdf(output, width=10, height=6)
par(mar = c(10,4,4,4))
print(gg1)
graphics.off()

