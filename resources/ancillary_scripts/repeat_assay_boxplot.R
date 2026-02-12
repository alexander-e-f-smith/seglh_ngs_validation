library(ggplot2)

#input delimited  file and data description - taken as script arguments
args <-commandArgs(trailingOnly = TRUE)
inputDF <-args[1] 
data_set_observed <-args[2]
data_set_expected <-args[3]

##read in delimited txt file for dataframe
df<-read.delim(inputDF)

#construct data labels
output<-paste0(data_set_observed,"_boxplot_ggplot.pdf")
label_title<-paste0("Interassay repeats - VAF spread;",data_set_observed,"(vs ",data_set_expected,")")
legend_label1<-paste0("Observed VAF (", data_set_observed,")\n used in boxplot calc for median/IQ ranges")
legend_label2<-paste0("Expected (Truth) VAF (", data_set_expected,")\n overlaid points not used in boxplot calc")

#boxplot with ggplot
gg1<-ggplot(df, aes(x = reorder(CHROM_POS, exp_AF), y = AF)) +
     geom_boxplot() +
     labs(x = "known variant", y = "VAF (log10)", title = label_title) +
     theme(text = element_text(size = 10),
        axis.text.x = element_text(angle = 90, hjust = 1))+
     geom_point(alpha=0.3, aes(x = reorder(CHROM_POS, exp_AF), y = AF, col="VAF"))+
     geom_point(alpha=0.3, aes(x = reorder(CHROM_POS, exp_AF), y = exp_AF, col="expected_VAF"))+ 
     scale_y_log10() + 
     scale_colour_manual(name="Data key:", breaks = c("VAF", "expected_VAF"), values=c("red", "blue"), labels=c(legend_label1, legend_label2))
#output as pdf
pdf(output, width=10, height=6)
par(mar = c(10,4,4,4))
print(gg1)
graphics.off()

##exmaple running : Rscript repeat_assay_boxplot.R "batch_combined_concordant_variants_A_headed_filtered_formatted.tsv" "DX output for Horizon control variants" "Snappy data comparison"
