#! /usr/bin/env RScript
df1 <- read.csv("/vol/data/PanKmer/results/pantools/db/Campylobacter_jejuni/pangenome_size/gene/gene_nrs_for_heaps_law.csv", header=FALSE)
df2 <- read.csv(file = "/vol/data/PanKmer/results/pantools/db/Campylobacter_jejuni/pangenome_size/gene/genes_for_heaps_law.csv", header = FALSE)
cat("\nThe Heaps law model is fitted to the number of new genes observed when genomes are ordered in a random way.\nIf alpha is higher than 1.0 the pangenome is closed, if alpha is below 1.0 it is open.\n\n")

x <- as.integer(df1)
y <- as.numeric(df2)
p0 <- c(mean(y[which(x == 2)] ), 1)

objectFun <- function(p, x, y) {
  y.hat <- p[1] * x^(-p[2])
  J <- sqrt(sum((y - y.hat)^2))/length(x)
  return(J)
}

fit <- optim(p0, objectFun, gr = NULL, x, y, method = "L-BFGS-B", lower = c(0, 0), upper = c(10000, 2))

sink("/vol/data/PanKmer/results/pantools/db/Campylobacter_jejuni/pangenome_size/gene/heaps_law_alpha.txt")
cat(fit$par[2])
sink()
cat("alpha =", fit$par[2], "\n\n")