## Version R version 3.6.3 (2020-02-29)
###ggxplot2_3.3.3 
### blandr_0.5.1


#####Bland_altmann plots

library(ggplot2)
library(blandr)

args <-commandArgs(trailingOnly = TRUE)
inputDF <- args[1]
data_comparison <-args[2]

input<-read.delim(inputDF)
df_eqa_exp<-(input$exp_AF)
df_eqa_obs<-(input$AF)
statistics_results <- blandr.statistics(df_eqa_exp, df_eqa_obs, sig.level=0.95)
plot_limits <- blandr.plot.limits(statistics_results, lowest_y_axis=FALSE)

##Print plots to a pdf file
output_name1=paste0(data_comparison,"_rplot.pdf")
output_name2=paste0(data_comparison,"_ggplot.pdf")
output_title=paste0("Bland-Altman plot(rplot): ",data_comparison," \nfor  methodology agreement/relative accuracy;\n Expected vs observed VAF")
output_title2=paste0("Bland-Altman plot(ggplot): ",data_comparison," \nfor  methodology agreement/relative accuracy;\n Expected vs observed VAF")
pdf(output_name1, width=10, height=6)

blandr.plot.rplot( statistics_results , plot_limits , plotTitle = output_title, annotate=TRUE)
dev.off()

pdf(output_name2, width=20, height=12)
blandr.plot.ggplot( statistics_results , plotTitle = output_title, overlapping = TRUE)
dev.off()
