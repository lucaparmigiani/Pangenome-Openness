#! /usr/bin/env RScript
#Assuming the R libraries were installed via Conda, the next line can be ignored. If not, use it to manually install the required package
#install.packages("ggplot2", "~/local/R_libs/", "https://cran.us.r-project.org")
library(ggplot2)

# AVERAGE: core, accessory and unique -> /vol/data/Pantools/db/Bacillus_cereus/pangenome_size/gene/gains_losses/average_all.csv
# MEDIAN : core, accessory and unique -> /vol/data/Pantools/db/Bacillus_cereus/pangenome_size/gene/gains_losses/median_all.csv
# AVERAGE: core & accessory -> /vol/data/Pantools/db/Bacillus_cereus/pangenome_size/gene/gains_losses/average_core_accessory.csv
# MEDIAN : core & accessory -> /vol/data/Pantools/db/Bacillus_cereus/pangenome_size/gene/gains_losses/median_core_accessory.csv

df = read.csv("/vol/data/Pantools/db/Bacillus_cereus/pangenome_size/gene/gains_losses/average_all.csv", sep=",",header = TRUE)
plot1 = ggplot(df, aes(x=x, y=y, color=Category)) +
   scale_color_manual(values=c("Accessory" = "#0066CC", "Core" = "#009900", "Unique" = "#990000")) +
   geom_line(size = 1.5) +
   theme(text = element_text(size=12)) +
   labs(x = "Number of genomes") +
   labs(y = "Homology group gain/loss") +
   theme_classic(base_size = 20) +
   theme(legend.position = "none") +
   #theme(legend.justification=c(1,1), legend.position=c(1,1)) +
   geom_hline(yintercept=0, linetype="dashed", color = "darkgrey")#+
   #ggtitle("Average number of genes gained/lost per added genome") +
   #ylim(-500, 500) #set a limit to the y axis if outliers stretch the plot too much

pdf(NULL)
ggsave("/vol/data/Pantools/db/Bacillus_cereus/pangenome_size/gene/gains_losses.png", plot= plot1, height = 10, width = 10)
cat("\nGain and losses between pangenome sizes written to: /vol/data/Pantools/db/Bacillus_cereus/pangenome_size/gene/gains_losses.png\n\n")