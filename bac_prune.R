# LD pruning homemade

data <- read.csv("~/rust_prune/rust_prune_4/3000_gts.csv")
haplotype_frequencies <- as.data.frame(table(data[,5], data[,2])/603)


### Data wrangle
# read in genotype
library(vcfR)
pyseer_snps <- read.vcfR( "~/Desktop/pyseer_tutorial/snps.vcf.gz", verbose = TRUE )
gt <- extract.gt(pyseer_snps, as.numeric = TRUE, IDtoRowNames = T)

# transpose genotype data (want samples in rows and variants in columns)
gt <- t(gt)

# temporarily - use only 3000 randomly sampled variants to make code run faster
gt <- gt[ , sample(ncol(gt),3000)]
V <- as.numeric(ncol(gt)) #as normal
N <- as.numeric(nrow(gt)) #as normal

# find and replace missing data in genotype file
# will use change all missing SNPs to 0
gt[is.na(gt)] <- as.numeric(0)

write.csv(gt, "3000_gts.csv", row.names = FALSE)

dummydat <- matrix(data=c(0,1), nrow=603, ncol=3000)
write.csv(dummydat, "dummydat.csv", row.names = F, col.names = T)

### Calculate MAF and AF from allele counts

# Compute MAF (and AF if desired)
find_MAF <- function(data, N, exclude_AF = FALSE) {
  MAF <- as.data.frame(colSums(data)/N)
  if (exclude_AF == FALSE) {
    AF <- 1-MAF
    MAF_AF <- cbind(MAF, AF)
    colnames(MAF_AF) <- c("MAF", "AF")
    return(MAF_AF)
  } else {
    return(MAF)
  }
}

# Discard variants with MAF below cutoff
MAF_prune <- function(data, N, cutoff) {
  MAF <- as.data.frame(colSums(data)/N)
  maf_prune <- which(MAF < cutoff)
  return(data[,-maf_prune])
}
#ensure to update V using V <- as.numeric(ncol(data)) after discarding variants


### Define functions

# Calculate the allele frequencies
find_allele_frequencies <- function(data, labels = TRUE) {
  allele_frequencies <- matrix(data=NA, nrow = V, ncol = 2)
  
  allele_frequencies[,2] <- colSums(data)
  allele_frequencies[,1] <- N - allele_frequencies[,2]
  allele_frequencies <- allele_frequencies/N
  
  if (labels == TRUE) {
    row.names(allele_frequencies) <- colnames(data)
    colnames(allele_frequencies) <- c("0", "1")
    return(allele_frequencies)
    } else {
    return(allele_frequencies)
  }
}

# Calculate haplotype frequencies for a given pair of variants
find_haplotype_frequencies <- function(data, variantA, variantB, compact = TRUE) {
  if (compact == FALSE) {
    haplotype_frequencies <- as.data.frame(table(data[,variantA], data[,variantB]))
    haplotype_frequencies[,3] <- haplotype_frequencies[,3]/N
    return(haplotype_frequencies)
  } else {
    haplotype_frequencies <- table(data[,variantA], data[,variantB])/N
    return(haplotype_frequencies)
  }
}

# Calculate D_prime
#standalone function
#note that variantA and variantB are the index number of the variant (1:V)

find_D_prime <- function(data, variantA, variantB, hap_freq_data, allele_freq_data) {
  if (missing(hap_freq_data)) {
    # find haplotype frequencies for the combination
    haplotype_frequencies <- table(data[,variantA], data[,variantB])/N
  } else {
    haplotype_frequencies <- hap_freq_data
  }
  if (missing(allele_freq_data)) {
    # find allele frequency for both variants
    allele_frequencies <- find_allele_frequencies(data)
  } else {
    allele_frequencies <- allele_freq_data
  }
  
  #calculate D
  D <- as.numeric(haplotype_frequencies[1,1]*haplotype_frequencies[2,2] - haplotype_frequencies[2,1]*haplotype_frequencies[1,2])
  #calculate D_max
  if (D < 0) {
    set_elements <- as.numeric(c((-1 * allele_frequencies[variantA,1] * allele_frequencies[variantB,1]) , (-1 * allele_frequencies[variantA,2] * allele_frequencies[variantB,2])))
    D_max <- max(set_elements)
  } else {
    set_elements <- as.numeric(c((allele_frequencies[variantA,1] * allele_frequencies[variantB,2]) , (allele_frequencies[variantA,2] * allele_frequencies[variantB,1])))
    D_max <- min(set_elements)
  }
  #calculate D_prime
  D_prime <- D/D_max
  return(D_prime)
}


### Generate pairwise matrix of D_prime values

D_prime_matrix <- matrix(data=NA, nrow = V, ncol = V) # set up matrix
pb = txtProgressBar(min = 0, max = V, initial = 0, style = 3) # set up progress bar

for (i in 1:V) {
  setTxtProgressBar(pb,i) # for progress bar
  for (j in i:V) {
    D_prime_matrix[i,j] <- find_D_prime(data = gt, variantA = i, variantB = j)
  }
  close(pb)
}

# It works!!


### LD Pruning

LD_prune <- function(data, D_prime_matrix, cutoff) {
  for (i in 1:V) {
    for (j in i:V) {
      if (D_prime_matrix[i,j] >= cutoff) {
        data <- data[,-as.numeric(colnames(gt)[j])]
      }
    }
  }
  return(data)
}

dataaaa <- LD_prune(gt, D_prime_matrix, 1)

### or 


LD_prune <- function(data, D_prime_matrix, cutoff) {
  prune_me_index <- 1
  prune_me <- vector("list", V)
  for (i in 1:V) {
    for (j in i:V) {
      if (D_prime_matrix[i,j] >= cutoff) {
        prune_me[[prune_me_index]] <- as.numeric(j) #add j to list of columns to be pruned
        prune_me_index <- prune_me_index + 1
      }
    }
  }
  print(prune_me)
  data <- data[,-prune_me[[1:length(prune_me)]]]
  return(data)
}

dataaaa <- LD_prune(gt, D_prime_matrix, 1)

# if any value in row is 1 then remove row

for (i in 1:nrow(D_prime_matrix)) {
  if (D_prime_matrix[i,] > 0.99) {
    remove <- i + remove
  }
}

gt <- gt[,-remove]




LD_prune <- function(data, D_prime_matrix, cutoff) {
  prune_me_index <- 1
  prune_me <- vector("list", V)
  for (i in 1:nrow(D_prime_matrix)) {
    if (any(D_prime_matrix[i,] > 0.99)) {
      prune_me[[prune_me_index]] <- as.numeric(i) #add i to list of columns to be pruned
      prune_me_index <- prune_me_index + 1
    }
  }
  data <- data[,-prune_me]
  return(data)
}

dataaaa <- LD_prune(gt, D_prime_matrix, 0.99)


LD_prune <- function(data, D_prime_matrix, cutoff) {
  for (i in 1:nrow(D_prime_matrix)) {
    if (any(D_prime_matrix[i,] > cutoff)) {
      prune_me <- append(prune_me, colnames(data)[i])
    }
  }
  print(length(prune_me))
  data <- data[,-!(names(data) %in% prune_me)]
  return(data)
}

dataaaa <- LD_prune(gt, D_prime_matrix, .99)

colnames(gt)[1]
length(colnames(gt))

# profiling tests

bench::mark(find_haplotype_frequencies(data = gt, variantA = 1, variantB = 2))
bench::mark(find_D_prime(data = gt, variantA = 1, variantB = 2))



