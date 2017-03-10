library(dplyr)

#epitopes = c("NLVPMVATV","KRWIILGLNK","GLCTLVAML","ATDALMTGY")

df = read.table("../../database/vdjdb.slim.txt", header=T, sep="\t", stringsAsFactors=F) %>%
  filter(gene == "TRB", species == "HomoSapiens") %>% #, mhc.class == "MHCI", antigen.epitope %in% epitopes) %>%
  distinct() %>%
  group_by(antigen.epitope) %>%
  mutate(cdr3.count = length(unique(cdr3))) %>%
  filter(cdr3.count >= 10) %>%
  select(cdr3, antigen.epitope, antigen.gene, antigen.species, mhc.class)

# remove 'cross-reactive' cdr3s

df.cdr3good = df %>%
  group_by(cdr3) %>%
  summarize(count = n()) %>%
  filter(count == 1)

df = subset(df, cdr3 %in% df.cdr3good$cdr3)

write.table(df, "vdjdb.nodes.txt", quote=F, sep="\t", row.names=F)
