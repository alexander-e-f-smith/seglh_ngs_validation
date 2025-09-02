#!/usr/bin/env Rscript

##ggpubr_0.4.0
##correlation plots
args <-commandArgs(trailingOnly = TRUE)
file_in1 <- args[1]
caller <- args[2] ##variant caller
reading <-args[3] ##eg VAF or VD
expected_reading<-args[4] ##field column name
observed_reading<-args[5] ##field column name
expected_source<-args[6] ##nature of expected readings
observed_source<-args[7] ##nature of observed readings
sample_name<-args[8]
xlimit<-args[9]
ylimit<-args[10]
xnum <-as.numeric(xlimit)
ynum <-as.numeric(ylimit)
output<-paste0(sample_name,"_",expected_source,"_",observed_source,"_",reading,"_correlation.pdf")

library(ggplot2)
library(scales)

file<-read.delim(file_in1)

library(ggpubr)
pdf(output, width=10, height=7)
ggscatter(file, x = observed_reading, y = expected_reading, color = "#5A5A5A", size = 1,
            add = "reg.line",
            add.params = list(color = "blue", fill = "gray"),
            conf.int = TRUE,
            cor.coef = TRUE, 
            cor.coeff.args = list(method = "pearson", r.digits=3, label.x = 0, label.y = ynum, color = "blue"),
            ylim=c(0, ynum),
            xlim=c(0, xnum),
            xlab = paste(reading, "expected","(",expected_source,")"), ylab = paste(reading, "observed","(",observed_source,")"),
            main = paste (reading,"(",caller,")","Trueness (Pearson Correlation 95% CI)\n Observed vs Expected","(",observed_source, 'vs', expected_source,")")) +
            stat_regline_equation(label.x = 0, label.y = ynum * 0.9) 
dev.off()
