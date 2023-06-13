# Lineage and phenotype correlation pruning 

### read in data
# read in genotype
library(vcfR)
pyseer_snps <- read.vcfR( "~/Desktop/pyseer_tutorial/snps.vcf.gz", verbose = TRUE )
genotype <- extract.gt(pyseer_snps, as.numeric = TRUE, IDtoRowNames = T)

# read in lineage
lineages <- read.csv("~/Desktop/pyseer_tutorial/database/database_clusters.csv")

# read in phenotype
metadata <- read.csv("~/Desktop/Massachusetts Data/sparc_metadata.csv")
resistances <- as.data.frame(metadata$Benzylpenicillin.MIC..ug.mL.)
resistances <- data.frame(MIC=resistances[,1], ID=metadata$Lane.name)

# ensure order/number of samples matches between these 3 data sets
# sort/filter those that dont
index <- match(dimnames(genotype)[[2]], resistances$ID)
resistances <- resistances[index,]

index <- match(dimnames(genotype)[[2]], lineages$Taxon)
lineages <- lineages[index,]

# transpose genotype data (need samples in rows and variants in columns)
genotype <- t(genotype)

# find and replace missing data in genotype file
# will use change all missing SNPs to 0
genotype[is.na(genotype)] <- 0


# Clean up the MIC values in the resistances dataset

# make equal values equal (eg 1.00 -> 1) -- is this losing precision in the data?
resistances$MIC <- replace(resistances$MIC, resistances$MIC=="0.50", 0.5)
resistances$MIC <- replace(resistances$MIC, resistances$MIC=="1.00", 1)
resistances$MIC <- replace(resistances$MIC, resistances$MIC=="2.00", 2)
resistances$MIC <- replace(resistances$MIC, resistances$MIC=="4.00", 4)
resistances$MIC <- replace(resistances$MIC, resistances$MIC==">4.00", ">4")
resistances$MIC <- replace(resistances$MIC, resistances$MIC=="<=0.03", "<0.03")
resistances$MIC <- replace(resistances$MIC, resistances$MIC=="<.016", "<0.016")

# find the order that the MIC categories should be in
MICs <- sort(unique(resistances$MIC))
desired_order <- c(MICs[1:2], MICs[5:6], MICs[3], MICs[7:23], MICs[4], MICs[24])

# convert MIC to sample,s ordinal category (1:24) in resistances
for (i in 1:nrow(resistances)) {
  resistances$MIC[i] <- which(desired_order == resistances$MIC[i])
}

# switch order of columns and rename
resistances <- data.frame(ID=resistances$ID, K=resistances[,1])
resistances$K <- as.numeric(resistances$K)


#can remove constant variants (where there is only the 0 allele, not 1) so no cor errors
gt <- genotype[,-which(colSums(genotype) == 0)]
gt <- gt[,-which(colSums(gt) == 603)]


### simple linear regression

library(tidyverse)
model <- lm(resistances$K ~ gt[,1])
model
ggplot(model, aes(gt[,1], resistances$K)) +
  geom_point() +
  stat_smooth(method = lm)



### simple linear correlation

# SNP v phenotype

#plot
model <- cor(gt, resistances$K, method = "spearman")
model <- as.data.frame(model)

xaxis <- seq(1, ncol(gt), by = 1)
ggplot(model, aes(x=xaxis, y=model$V1)) + 
  geom_point()

# find p values of correlation
p_values <- c()
for (i in 1:ncol(gt)) {
  p_value <- cor.test(gt[,i], resistances$K, method = "spearman")
  p_values <- append(p_values, p_value$p.value)
}

# prune by p values
prune_index <- c()
for (i in 1:ncol(gt)) {
  if (p_values[i] > 0.05) {
    prune_index <- append(prune_index, i)
  }
}
pruned_gt <- gt[,-prune_index]


# prune by correlation (rho) threshold
prune_index <- c()
threshold <- 0.4
for (i in 1:ncol(gt)) {
  if (abs(model$V1[i]) < threshold) {
    prune_index <- append(prune_index, i)
  }
}
pruned_gt <- gt[,-prune_index]

# observe pruned data
model <- cor(pruned_gt, resistances$K, method = "spearman")
model <- as.data.frame(model)

xaxis <- seq(1, ncol(pruned_gt), by = 1)
ggplot(model, aes(x=xaxis, y=model$V1)) + 
  geom_point()


# SNP v lineage

#without SNPvphen pruning
model <- cor(gt, lineages$Cluster, method = "spearman")
model <- as.data.frame(model)

xaxis <- seq(1, ncol(gt), by = 1)
ggplot(model, aes(x=xaxis, y=model$V1)) + 
  geom_point()

#with SNPvphen pruning

# plot
model <- cor(pruned_gt, lineages$Cluster, method = "spearman")
model <- as.data.frame(model)

xaxis <- seq(1, ncol(pruned_gt), by = 1)
ggplot(model, aes(x=xaxis, y=model$V1)) + 
  geom_point()

# find p values of correlation
p_values <- c()
for (i in 1:ncol(pruned_gt)) {
  p_value <- cor.test(pruned_gt[,i], lineages$Cluster, method = "spearman")
  p_values <- append(p_values, p_value$p.value)
}

# prune by p values
prune_index <- c()
for (i in 1:ncol(pruned_gt)) {
  if (p_values[i] > 0.2) {
    prune_index <- append(prune_index, i)
  }
}
pruned_gt <- pruned_gt[,-prune_index]


# prune by correlation (rho)
prune_index <- c()
threshold <- 0.1
for (i in 1:ncol(pruned_gt)) {
  if (abs(model$V1[i]) < threshold) {
    prune_index <- append(prune_index, i)
  }
}
pruned_gt <- pruned_gt[,-prune_index]


# observe pruned data by correlation score (rho)
model <- cor(pruned_gt, resistances$K, method = "spearman")
model <- as.data.frame(model)

xaxis <- seq(1, ncol(pruned_gt), by = 1)
ggplot(model, aes(x=xaxis, y=model$V1)) + 
  geom_point()

# observed pruned data by p value

p_values <- c()
for (i in 1:ncol(pruned_gt)) {
  p_value <- cor.test(pruned_gt[,i], resistances$K, method = "spearman")
  p_values <- append(p_values, p_value$p.value)
}
xaxis <- seq(1, ncol(pruned_gt), by = 1)
ggplot(model, aes(x=xaxis, y=-log10(p_values))) + 
  geom_point()



# Other
k <- kruskal.test(x = linphen$V2, g = linphen$V1)
