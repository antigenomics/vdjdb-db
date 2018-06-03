library(dplyr)
library(stringr)
library(data.table)

### VDJdb processing

proc_slim_vdjdb <- function(.df) {
  .df %>%
    filter(v.segm != "", j.segm != "") %>%
    mutate(v.segm = str_split_fixed(v.segm, "[*,]", 2)[,1],
           j.segm = str_split_fixed(j.segm, "[*,]", 2)[,1]) %>%
    select(species, gene, antigen.epitope, cdr3, v.segm, j.segm) %>%
    unique
}

### Segment processing
get_segment_parts <- function(.df) {
  .df.v <- .df %>%
    select(species, gene, v.segm, v.end, cdr3) %>%
    filter(v.end > 0) %>%
    mutate(v.segm = str_split_fixed(v.segm, "[*,]", 2)[,1],
           cdr3 = substr(cdr3, 1, v.end)) %>%
    group_by(species, gene, v.segm, cdr3) %>%
    summarise(count = n()) %>%
    group_by(species, gene, cdr3, type = "V") %>%
    summarise(segm = v.segm[which(count == max(count))][1])
  
  .df.j <- .df %>%
    select(species, gene, j.segm, j.start, cdr3) %>%
    filter(j.start > 0) %>%
    mutate(j.segm = str_split_fixed(j.segm, "[*,]", 2)[,1],
           cdr3 = substr(cdr3, j.start, nchar(cdr3))) %>%
    group_by(species, gene, j.segm, cdr3) %>%
    summarise(count = n()) %>%
    group_by(species, gene, cdr3, type = "J") %>%
    summarise(segm = j.segm[which(count == max(count))][1])
  
  rbind(.df.v, .df.j)
}

### VDJtools export

mock_codons <- c('GCT', 'TGT', 'GAT', 'GAA', 'TTT',
                 'GGT', 'ATT', 'CAT', 'AAA', 'TTA',
                 'ATG', 'AAT', 'CCT', 'CAA', 'CGT',
                 'TCT', 'ACT', 'GTT', 'TGG', 'TAT')

names(mock_codons) <- c('A', 'C', 'D', 'E', 'F',
                        'G', 'I', 'H', 'K', 'L',
                        'M', 'N', 'P', 'Q', 'R',
                        'S', 'T', 'V', 'W', 'Y')

mock_back_translate <- function(x) {
  paste0(mock_codons[x], collapse = "")
}

# "CASS" %>% strsplit('') %>% lapply(mock_back_translate)

as.vdjtools.df <- function(.df) {
  .df$cdr3aa <- .df$cdr3
  .df$cdr3nt <- .df$cdr3aa %>% 
    strsplit('') %>% 
    lapply(mock_back_translate)
  
  .df %>%
    mutate(count = 1, v = v.segm, d = "", j = j.segm,
           vend = -1, dstart = -1, dend = -1, jstart = -1) %>%
    group_by(species, gene, antigen.epitope) %>%
    mutate(freq = 1 / n()) %>%
    ungroup %>%
    select(count, freq, cdr3nt, cdr3aa, v, d, j, vend, dstart, dend, jstart,
           species, gene, antigen.epitope)
}