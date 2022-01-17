library("RBiotools")
library("stringr")

df <- read.csv("results/sequences.csv")
df <- df[c("accession", "identifier", "protein")]

proteinGrouping <- runLinclust_DF(df)

write.csv(proteinGrouping, "results/clusters.csv")