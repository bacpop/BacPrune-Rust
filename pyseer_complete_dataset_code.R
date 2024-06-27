# Code for obtaining the genotype data to be used in LD pruning (and then BBGWAS)

# read in genotype
library(vcfR)
pyseer_snps <- read.vcfR( "~/Desktop/pyseer_tutorial/snps.vcf.gz", verbose = TRUE )
gt <- extract.gt(pyseer_snps, as.numeric = TRUE, IDtoRowNames = T)

# read in lineage
lineages <- read.csv("~/Desktop/pyseer_tutorial/database/database_clusters.csv")

# read in phenotype
metadata <- read.csv("~/Desktop/Massachusetts Data/sparc_metadata.csv")
resistances <- as.data.frame(metadata$Benzylpenicillin.MIC..ug.mL.)
resistances <- data.frame(MIC=resistances[,1], ID=metadata$Lane.name)

# ensure order/number of samples matches between these 3 data sets
# sort/filter those that dont
index <- match(dimnames(gt)[[2]], resistances$ID)
resistances <- resistances[index,]

index <- match(dimnames(gt)[[2]], lineages$Taxon)
lineages <- lineages[index,]

# transpose genotype data (need samples in rows and variants in columns)
gt <- t(gt)

# find and replace missing data in genotype file
# will use change all missing SNPs to 0
gt[is.na(gt)] <- 0

# write to folder
write.csv(gt, "~/Documents/GitHub/BacPrune-Rust/pyseer_complete.csv", row.names = FALSE)
