df <- read.table("../database/vdjdb.txt", header=T, sep="\t")

df$comment <- NULL

cat("Summary:\n")
str(df)
summary(df)

cat("Species:\n")
cat(levels(df$species), "\n")

cat("Ag-species:\n")
cat(levels(df$antigen.species), "\n")
cat("Ag-genes:\n")
cat(levels(df$antigen.gene), "\n")

cat("MHC-A:\n")
cat(levels(df$mhc.a), "\n")
cat("MHC-B:\n")
cat(levels(df$mhc.b), "\n")

cat("References:\n")
cat(levels(df$reference.id), "\n")

cat("Last record id:\n")
cat(as.character(sort(df$record.id, decreasing = T)[1]), "\n")
cat("Last complex id:\n")
cat(as.character(sort(df$complex.id, decreasing = T)[1]), "\n")