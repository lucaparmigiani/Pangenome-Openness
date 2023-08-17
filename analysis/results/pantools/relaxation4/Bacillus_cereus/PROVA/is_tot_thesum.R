library(dplyr)

core_acc = read.csv("core_accessory_size.csv")
core_acc_uniq = read.csv("core_accessory_unique_size.csv")
head(core_acc)
head(core_acc_uniq)
core_acc$Category = as.factor(core_acc$Category)
core_acc_uniq$Category = as.factor(core_acc_uniq$Category)

core_acc %>% filter(Category == "Accessory")
colSums(core_acc %>% filter((Category == "Core" | Category == "Accessory" ) & x== 2) %>% select(y))
colSums(core_acc %>% filter((Category == "Core" | Category == "Accessory" ) & x== 3) %>% select(y))
colSums(core_acc %>% filter((Category == "Core" | Category == "Accessory" ) & x== 82) %>% select(y))
colSums(core_acc %>% filter((Category == "Core" | Category == "Accessory" ) & x== 83) %>% select(y))
max(core_acc$x)
core_acc[core_acc$x==83,]
#length(core_acc[core_acc$x==82,]) != length(core_acc[core_acc$x==83,])
colSums(core_acc_uniq %>% filter((Category == "Core" | Category == "Accessory" | Category == "Unique") & x== 2) %>% select(y))
colSums(core_acc_uniq %>% filter((Category == "Core" | Category == "Accessory" | Category == "Unique") & x== 3) %>% select(y))
colSums(core_acc_uniq %>% filter((Category == "Core" | Category == "Accessory" | Category == "Unique") & x== 82) %>% select(y))
colSums(core_acc_uniq %>% filter((Category == "Core" | Category == "Accessory" | Category == "Unique") & x== 83) %>% select(y))

colSums(core_acc %>% filter((Category == "Median") & x== 30) %>% select(y))
colSums(core_acc_uniq %>% filter((Category== "Median") & x== 30) %>% select(y))
robustbase::colMedians(as.matrix(core_acc_uniq %>% filter((Category == "Core" | Category == "Accessory" | Category == "Unique") & x== 83) %>% select(y)))
core_acc_uniq %>% filter((Category == "Core" | Category == "Accessory" | Category == "Unique") & x== 82) %>% group_by(Category) %>% count()
colSums(core_acc_uniq %>% filter((Category== "Median") & x== 83) %>% select(y))


pan = read.csv("/home/luca/@focus/code/45_comparison/data/pantools/Bacillus_cereus/gene_classification/classified_groups.csv")
head(pan)
pan[ncol(pan)] = NULL
panmatrix = as.matrix(pan[,3:ncol(pan)])
head(panmatrix)

colSums(panmatrix)
mean(colSums(panmatrix))
