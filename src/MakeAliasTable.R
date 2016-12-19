add_aliases <- function (.vec) {
  res = data.frame(Bad = "", Good = "")
  for (gene in .vec) {
    # TRAV-13-1 -> TRAV-13
    if (length(grep("-", gene, perl = F)) == 0) {
      res = rbind(res, data.frame(Bad = paste0(gene, "-1"), Good = gene))
    }
    
    # TRAV-12-!
    if (substr(gene, nchar(gene), nchar(gene)) == "1") {
      res = rbind(res, data.frame(Bad = paste0(substr(gene, 1, nchar(gene) - 1), "!"), Good = gene))
    }
    
    if (length(grep("/", gene, perl = F)) > 0) {
      pos = gregexpr("/", gene, perl = F)[[1]][1]
      
      # TRAV14 -> TRAV14/DV4
      res = rbind(res, data.frame(Bad = substr(gene, 1, pos - 1), Good = gene))
      
      # TRDV4 -> TRAV14/DV4
      res = rbind(res, data.frame(Bad = paste0(substr(gene, 1, 2), 
                                               substr(gene, pos + 1, nchar(gene))), 
                                  Good = gene))
      
      # TRADV14 -> TRAV14/DV4
      res = rbind(res, data.frame(Bad = paste0(substr(gene, 1, 3), 
                                               substr(gene, pos + 1, pos + 2), 
                                               substr(gene, pos + 3, nchar(gene))), 
                                  Good = gene))
      
      # TRADV4 -> TRAV14/DV4
      res = rbind(res, data.frame(Bad = paste0(substr(gene, 1, 3), 
                                               substr(gene, pos + 1, pos + 2), 
                                               substr(gene, 5, pos - 1)), 
                                  Good = gene))
      
      # TRAV14_DV4   -> TRAV14/DV4
      res = rbind(res, data.frame(Bad = sub("/", "_", gene, perl = F), Good = gene))
    }
  }
  
  res[,1] <- as.character(res[,1])
  res[,2] <- as.character(res[,2])
  res[-1,]
}

segm_df = read.csv("segments.txt", sep = "\t", header = T, stringsAsFactors = F)
segm_df$id = substr(segm_df$id, 1, regexpr("*", segm_df$id, fixed = T) - 1)
res = aggregate(segm_df$id, list(segm_df$X.species, segm_df$gene, segm_df$segment), function (x) {
  add_aliases(x)
}, simplify = F)

logic = sapply(1:length(res$x), function (a) nrow(res$x[a][[1]])) > 0
res = lapply((1:length(res[[1]]))[logic > 0], function (i) { print(i); data.frame(sapply(res, "[[", i), stringsAsFactors = F) })
res = do.call(rbind, res)
res = res[order(res[,1], res[,2], res[,3], res[,5]),]
colnames(res) = c("species", "gene", "segment", "bad", "good")

write.table(res, "segments.alias.txt", sep = "\t", row.names = F, quote = F)