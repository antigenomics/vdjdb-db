# Format the reviewed held-out metric table for the documentation figure; no refit.
data <- read.delim("metrics.tsv", check.names = FALSE)
rows <- list()
for (epitope in c("GLCTLVAML", "YLQPRTFLL")) {
  subset <- data[data$epitope == epitope, ]
  for (metric in c("ap", "accuracy", "balanced_accuracy")) {
    value <- function(name) subset[subset$predictor == name, metric]
    rows[[length(rows) + 1L]] <- data.frame(
      epitope = epitope, metric = metric, iptm = value("ipTM"),
      motif = value("tcremp_TRA"), v3 = value("v3"))
  }
}
write.table(do.call(rbind, rows), "plot.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
