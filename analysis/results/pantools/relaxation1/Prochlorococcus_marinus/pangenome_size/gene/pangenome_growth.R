#! /usr/bin/env RScript
#Assuming the R libraries were installed via Conda, the next lines can be ignored. If not, use them to manually install the required packages
#install.packages("ggplot2", "~/local/R_libs/", "https://cran.us.r-project.org")
#install.packages("svglite", "~/local/R_libs/", "https://cran.us.r-project.org")

library(ggplot2)
#library(svglite) #required to save plot as SVG
int_breaks <- function(x, n = 5) { # to make sure line breaks are integers
  l <- pretty(x, n)
  l[abs(l %% 1) < .Machine$double.eps ^ 0.5]
}

df = read.csv("/vol/data/PanKmer/results/pantools/db/Prochlorococcus_marinus/pangenome_size/gene/core_accessory_unique_size.csv", sep=",", header = TRUE)
plot1 = ggplot(df, aes(x=x, y=y, color=Category, size=Category)) +
    scale_color_manual(values=c("Median" = "#000000", "Accessory" = "#0072B2", "Core" = "#009E73", "Unique" = "#D55E00")) +
    scale_size_manual(values=c("Median" = 0.2, "Accessory" = 0.8, "Core" = 0.8, "Unique" = 0.8)) +
    geom_point() + 
    theme(text = element_text(size=12)) + 
    labs(x = "Number of genomes") + 
    labs(y = "Number of homology groups") +
    scale_x_continuous(breaks = int_breaks) +

    # Three options for the legend.
    # 1. Uncomment the next two lines for a legend to the right of the plot
    # 2. Only uncomment the first line to have a legend inside the plot
    # 3. Only uncomment the second line to have no legend

    #theme(legend.position = "none") + # line 1
    theme(legend.justification=c(1,1), legend.position=c(1,1)) + # line 2
    guides(colour = guide_legend(override.aes = list(size=5))) +
    theme_bw(base_size = 18) # increase this number of have larger text and numbers

pdf(NULL)
#ggsave("/vol/data/PanKmer/results/pantools/db/Prochlorococcus_marinus/pangenome_size/gene/core_accessory_unique_growth.svg", plot= plot1, height = 10, width = 10)
ggsave("/vol/data/PanKmer/results/pantools/db/Prochlorococcus_marinus/pangenome_size/gene/core_accessory_unique_growth.png", plot= plot1, height = 10, width = 10)
cat("\nPangenome growth written to: /vol/data/PanKmer/results/pantools/db/Prochlorococcus_marinus/pangenome_size/gene/core_accessory_unique_growth.png\n\n")

df = read.csv("/vol/data/PanKmer/results/pantools/db/Prochlorococcus_marinus/pangenome_size/gene/core_dispensable_size.csv", sep=",", header = TRUE)
plot2 = ggplot(df, aes(x=x, y=y, color=Category, size=Category)) +
    scale_color_manual(values=c("Median" = "#000000", "Core" = "#009E73", "Dispensable" = "#D55E00")) +
    scale_size_manual(values=c("Median" = 0.2, "Core" = 0.8, "Dispensable" = 0.8)) +
    geom_point() + 
    theme(text = element_text(size=12)) + 
    labs(x = "Number of genomes") + 
    labs(y = "Number of homology groups") +
    scale_x_continuous(breaks = int_breaks) +

    # Three options for the legend.
    # 1. Uncomment the next two lines for a legend to the right of the plot
    # 2. Only uncomment the first line to have a legend inside the plot
    # 3. Only uncomment the second line to have no legend

    #theme(legend.position = "none") + # line 1
    theme(legend.justification=c(1,1), legend.position=c(1,1)) + # line 2
    guides(colour = guide_legend(override.aes = list(size=5))) +
    theme_bw(base_size = 18) # increase this number of have larger text and numbers

pdf(NULL)
#ggsave("/vol/data/PanKmer/results/pantools/db/Prochlorococcus_marinus/pangenome_size/gene/core_dispensable_growth.svg", plot= plot2, height = 10, width = 10)
ggsave("/vol/data/PanKmer/results/pantools/db/Prochlorococcus_marinus/pangenome_size/gene/core_dispensable_growth.png", plot= plot2, height = 10, width = 10)
cat("\nPangenome growth written to: /vol/data/PanKmer/results/pantools/db/Prochlorococcus_marinus/pangenome_size/gene/core_dispensable_growth.png\n\n")

df = read.csv("/vol/data/PanKmer/results/pantools/db/Prochlorococcus_marinus/pangenome_size/gene/core_dispensable_total_size.csv", sep=",", header = TRUE)
plot3 = ggplot(df, aes(x=x, y=y, color=Category, size=Category)) +
    scale_color_manual(values=c("Median" = "#000000", "Core" = "#009E73", "Dispensable" = "#D55E00", "Total" = "#0072B2")) +
    scale_size_manual(values=c("Median" = 0.2, "Core" = 0.8, "Dispensable" = 0.8, "Total" = 0.8)) +
    geom_point() + 
    theme(text = element_text(size=12)) + 
    labs(x = "Number of genomes") + 
    labs(y = "Number of homology groups") +
    scale_x_continuous(breaks = int_breaks) +

    # Three options for the legend.
    # 1. Uncomment the next two lines for a legend to the right of the plot
    # 2. Only uncomment the first line to have a legend inside the plot
    # 3. Only uncomment the second line to have no legend

    #theme(legend.position = "none") + # line 1
    theme(legend.justification=c(1,1), legend.position=c(1,1)) + # line 2
    guides(colour = guide_legend(override.aes = list(size=5))) +
    theme_bw(base_size = 18) # increase this number of have larger text and numbers

pdf(NULL)
#ggsave("/vol/data/PanKmer/results/pantools/db/Prochlorococcus_marinus/pangenome_size/gene/core_dispensable_total_growth.svg", plot= plot3, height = 10, width = 10)
ggsave("/vol/data/PanKmer/results/pantools/db/Prochlorococcus_marinus/pangenome_size/gene/core_dispensable_total_growth.png", plot= plot3, height = 10, width = 10)
cat("\nPangenome growth written to: /vol/data/PanKmer/results/pantools/db/Prochlorococcus_marinus/pangenome_size/gene/core_dispensable_total_growth.png\n\n")

