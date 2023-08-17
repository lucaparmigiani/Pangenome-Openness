#! /usr/bin/env RScript
#Assuming the R libraries were installed via Conda, the next line can be ignored. If not, use it to manually install the required package
#install.packages("ggplot2", "~/local/R_libs/", "https://cran.us.r-project.org")
library(ggplot2)

# core, accessory and unique -> /vol/data/PanKmer/results/pantools/db/Bacillus_cereus/pangenome_size/gene/gains_losses/average_all.csv
#                            -> /vol/data/PanKmer/results/pantools/db/Bacillus_cereus/pangenome_size/gene/gains_losses/median_all.csv
# core & accessory -> /vol/data/PanKmer/results/pantools/db/Bacillus_cereus/pangenome_size/gene/gains_losses/average_core_accessory.csv
#                  -> /vol/data/PanKmer/results/pantools/db/Bacillus_cereus/pangenome_size/gene/gains_losses/median_core_accessory.csv

df_average = read.csv("/vol/data/PanKmer/results/pantools/db/Bacillus_cereus/pangenome_size/gene/gains_losses/average_all.csv", sep=",", header = TRUE)
df_median = read.csv("/vol/data/PanKmer/results/pantools/db/Bacillus_cereus/pangenome_size/gene/gains_losses/median_all.csv", sep=",", header = TRUE)
plot1 = ggplot(df_average, aes(x=x, y=y, color=Category, linetype=Metric))
if (packageVersion("ggplot2") < "3.4.0") {
   plot1 = plot1 + geom_line(size = 1.0, data = df_median, aes(x =x, y =y, color=Category)) + #increase linewidth for ggplot2 version < 3.4.0
   geom_line(size = 1.7) #increase linewidth for ggplot2 version < 3.4.0
} else {
   plot1 = plot1 + geom_line(linewidth = 1.0, data = df_median, aes(x =x, y =y, color=Category)) + #increase linewidth for ggplot2 version >= 3.4.0
   geom_line(linewidth = 1.7) #increase linewidth for ggplot2 version >= 3.4.0
}
plot1 = plot1 + scale_color_manual(values=c("Accessory" = "#0072B2", "Core" = "#009E73", "Unique" = "#D55E00")) +
   scale_linetype_manual(values=c("Average" = "dotted", "Median" = "solid")) +
   labs(x = "Number of genomes") +
   labs(y = "Homology group gain/loss") +
   geom_hline(yintercept=0, linetype="dashed", color = "darkgrey") +
   #ggtitle("Average number of genes gained/lost per added genome") +
   #ylim(-500, 500) + #set a limit to the y axis if outliers stretch the plot too much

   # Three options for the legend.
   # 1. Uncomment the next two lines for a legend to the right of the plot
   # 2. Only uncomment the first line to have a legend inside the plot
   # 3. Only uncomment the second line to have no legend

   #theme(legend.position = "none") + # line 1
   theme(legend.justification=c(1,1), legend.position=c(1,1)) + # line 2
   guides(colour = guide_legend(override.aes = list(size=5))) +
   theme_bw(base_size = 18) # increase this number of have larger text and numbers

pdf(NULL)
ggsave("/vol/data/PanKmer/results/pantools/db/Bacillus_cereus/pangenome_size/gene/gains_losses_core_accessory_unique.png", plot= plot1, height = 10, width = 10)
cat("\nGain and losses between pangenome sizes written to: /vol/data/PanKmer/results/pantools/db/Bacillus_cereus/pangenome_size/gene/gains_losses_core_accessory_unique.png\n\n")

df_average = read.csv("/vol/data/PanKmer/results/pantools/db/Bacillus_cereus/pangenome_size/gene/gains_losses/average_core_accessory.csv", sep=",", header = TRUE)
df_median = read.csv("/vol/data/PanKmer/results/pantools/db/Bacillus_cereus/pangenome_size/gene/gains_losses/median_core_accessory.csv", sep=",", header = TRUE)
plot1 = ggplot(df_average, aes(x=x, y=y, color=Category, linetype=Metric))
if (packageVersion("ggplot2") < "3.4.0") {
   plot1 = plot1 + geom_line(size = 1.0, data = df_median, aes(x =x, y =y, color=Category)) + #increase linewidth for ggplot2 version < 3.4.0
   geom_line(size = 1.7) #increase linewidth for ggplot2 version < 3.4.0
} else {
   plot1 = plot1 + geom_line(linewidth = 1.0, data = df_median, aes(x =x, y =y, color=Category)) + #increase linewidth for ggplot2 version >= 3.4.0
   geom_line(linewidth = 1.7) #increase linewidth for ggplot2 version >= 3.4.0
}
plot1 = plot1 + scale_color_manual(values=c("Accessory" = "#0072B2", "Core" = "#009E73", "Unique" = "#D55E00")) +
   scale_linetype_manual(values=c("Average" = "dotted", "Median" = "solid")) +
   labs(x = "Number of genomes") +
   labs(y = "Homology group gain/loss") +
   geom_hline(yintercept=0, linetype="dashed", color = "darkgrey") +
   #ggtitle("Average number of genes gained/lost per added genome") +
   #ylim(-500, 500) + #set a limit to the y axis if outliers stretch the plot too much

   # Three options for the legend.
   # 1. Uncomment the next two lines for a legend to the right of the plot
   # 2. Only uncomment the first line to have a legend inside the plot
   # 3. Only uncomment the second line to have no legend

   #theme(legend.position = "none") + # line 1
   theme(legend.justification=c(1,1), legend.position=c(1,1)) + # line 2
   guides(colour = guide_legend(override.aes = list(size=5))) +
   theme_bw(base_size = 18) # increase this number of have larger text and numbers

pdf(NULL)
ggsave("/vol/data/PanKmer/results/pantools/db/Bacillus_cereus/pangenome_size/gene/gains_losses_core_accessory.png", plot= plot1, height = 10, width = 10)
cat("Gain and losses between pangenome sizes written to: /vol/data/PanKmer/results/pantools/db/Bacillus_cereus/pangenome_size/gene/gains_losses_core_accessory.png\n\n")
