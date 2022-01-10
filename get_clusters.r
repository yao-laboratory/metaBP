library("RBiotools")
library("stringr")

df <- read.csv("sequences.csv")
df <- df[c("accession", "identifier", "protein")]

proteinGrouping <- runLinclust_DF(df)

write.csv(proteinGrouping, "clusters.csv")