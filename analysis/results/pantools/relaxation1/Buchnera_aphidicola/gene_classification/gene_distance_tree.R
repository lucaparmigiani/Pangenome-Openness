#! /usr/bin/env RScript

# Use one of the following files
# /vol/data/PanKmer/results/pantools/db/Buchnera_aphidicola/gene_classification/distance_distinct_genes.csv
# /vol/data/PanKmer/results/pantools/db/Buchnera_aphidicola/gene_classification/distance_inf_distinct_genes.csv
# /vol/data/PanKmer/results/pantools/db/Buchnera_aphidicola/gene_classification/distance_all_genes.csv

#Assuming the R libraries were installed via Conda, the next line can be ignored. If not, use it to manually install the required package
#install.packages("ape", "~/local/R_libs/", "https://cran.us.r-project.org")
library(ape)

input = read.csv("/vol/data/PanKmer/results/pantools/db/Buchnera_aphidicola/gene_classification/distance_distinct_genes.csv", sep=",", header = TRUE)
input[37] <- NULL
df2 = subset(input, select = -c(Genomes))
df.dist2 =as.matrix(df2, labels=TRUE)
colnames(df.dist2) <- rownames(df.dist2) <- input[['Genomes']]
NJ_tree <- nj(df.dist2)
pdf(NULL)
plot(NJ_tree, main = "Neighbor Joining")
write.tree(NJ_tree, tree.names = TRUE, file="/vol/data/PanKmer/results/pantools/db/Buchnera_aphidicola/gene_classification/gene_distance.tree")
cat("\nGene distance tree written to: /vol/data/PanKmer/results/pantools/db/Buchnera_aphidicola/gene_classification/gene_distance.tree\n\n")