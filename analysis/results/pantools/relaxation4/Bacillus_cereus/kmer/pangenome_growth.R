#! /usr/bin/env RScript
#Assuming the R libraries were installed via Conda, the next lines can be ignored. If not, use them to manually install the required packages
#install.packages("ggplot2", "~/local/R_libs/", "https://cran.us.r-project.org")
#install.packages("dplyr", "~/local/R_libs/", "https://cloud.r-project.org")
#install.packages("scales", "~/local/R_libs/", "https://cran.us.r-project.org")

library(ggplot2)
#library(scales) #uncomment when you want to show the total number of k-mers in millions instead of a scientific number
#library(dplyr) #uncomment when you want to show the total number of k-mers in millions instead of a scientific number

int_breaks <- function(x, n = 5) { # to make sure line breaks are integers
  l <- pretty(x, n)
  l[abs(l %% 1) < .Machine$double.eps ^ 0.5]
}

df = read.csv("/vol/data/Pantools/db/Bacillus_cereus/pangenome_size/kmer/core_accessory_size.csv", sep=",", header = TRUE)
#df <- mutate(df,y = y / 10^6) # convert scientific number to millions (requires dplyr library)
plot1 = ggplot(df, aes(x=x, y=y, color=Category, size=Category)) +
    scale_color_manual(values=c("Median" = "#000000", "Accessory" = "#0066CC", "Core" = "#009900", "Unique" = "#990000")) +
    scale_size_manual(values=c("Median" = 0.2, "Accessory" = 0.8, "Core" = 0.8, "Unique" = 0.8)) +
    geom_point() + 
    theme(text = element_text(size=12)) + 
    labs(x = "Number of genomes") + 
    #scale_y_continuous(labels = unit_format(unit = "M")) + # convert scientific number to millions (requires scales library)
    labs(y = "K-mers") +
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
ggsave("/vol/data/Pantools/db/Bacillus_cereus/pangenome_size/kmer/core_accessory_growth.png", plot= plot1, height = 10, width = 10)
cat("\nPangenome growth written to: /vol/data/Pantools/db/Bacillus_cereus/pangenome_size/kmer/core_accessory_growth.png\n\n")

int_breaks <- function(x, n = 5) { # to make sure line breaks are integers
  l <- pretty(x, n)
  l[abs(l %% 1) < .Machine$double.eps ^ 0.5]
}

df = read.csv("/vol/data/Pantools/db/Bacillus_cereus/pangenome_size/kmer/core_accessory_unique_size.csv", sep=",", header = TRUE)
#df <- mutate(df,y = y / 10^6) # convert scientific number to millions (requires dplyr library)
plot1 = ggplot(df, aes(x=x, y=y, color=Category, size=Category)) +
    scale_color_manual(values=c("Median" = "#000000", "Accessory" = "#0066CC", "Core" = "#009900", "Unique" = "#990000")) +
    scale_size_manual(values=c("Median" = 0.2, "Accessory" = 0.8, "Core" = 0.8, "Unique" = 0.8)) +
    geom_point() + 
    theme(text = element_text(size=12)) + 
    labs(x = "Number of genomes") + 
    #scale_y_continuous(labels = unit_format(unit = "M")) + # convert scientific number to millions (requires scales library)
    labs(y = "K-mers") +
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
ggsave("/vol/data/Pantools/db/Bacillus_cereus/pangenome_size/kmer/core_accessory_unique_growth.png", plot= plot1, height = 10, width = 10)
cat("\nPangenome growth written to: /vol/data/Pantools/db/Bacillus_cereus/pangenome_size/kmer/core_accessory_unique_growth.png\n\n")
