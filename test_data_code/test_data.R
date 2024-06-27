#testing rust LD prune
set.seed(238529)

# test on dataset of two loci in perfect LD

doc <- sample(c(0,1), size = 300, replace = TRUE)
test_data <- as.data.frame(c(doc, doc))


# test on dataset of four loci in perfect LD

test_data <- as.data.frame(doc, doc, doc, doc)

# test on dataset of four loci in perfect LD, four loci in partial LD

doc2a <- sample(c(0,1), size = 300, replace = TRUE)
doc2b <- doc2a
doc2b[5,9] <- 0
doc2c <- doc2a
doc2c[7,1] <- 1
test_data <- as.data.frame(doc, doc, doc2a, doc, doc2b, doc, doc2c, doc2d)