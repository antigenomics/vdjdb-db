---
title: "VDJdb summary statistics"
author: "Mikhail Shugay"
date: "26-02-2022"
output:
  html_document:
    code_folding: hide
params:
  legacy: "../out/legacy"
  reference_years: "reference_years.tsv"
  annotations: "annotations.tsv"
---

<!--
The release dashboard. Its published fragment -- everything between `!summary_embed_start!` and
`!summary_embed_end!` -- is extracted by `MakeEmbedableHtml.py` and injected into vdjdb-web's
`/overview`, so three properties of the output are load-bearing and asserted by
`summary/check_summary.py`: the base64 images are single unbroken lines, there are no `<div>`
wrappers, and the filename does not change.

The paper figures live in `vdjdb_paper_figures.Rmd`, split off at the marker below. That is what
keeps `maps`, `scatterpie` and `ggrepel` out of this render path.

Every input is a parameter, so the render reads the build under `out/` rather than the stale 2024
copy in `database/`, and it makes no network call:

    Rscript -e 'rmarkdown::render("summary/vdjdb_summary.Rmd")'
    Rscript -e 'rmarkdown::render("summary/vdjdb_summary.Rmd", params = list(legacy = "../out/legacy"))'
-->




``` r
library(knitr)
library(ggplot2)
library(RColorBrewer)
library(data.table)
library(forcats)
library(ggalluvial)
library(circlize)
library(tidyverse)
library(stringr)
library(gridExtra)
library(cowplot)
select = dplyr::select

df = fread(file.path(params$legacy, "vdjdb.slim.txt"), header=T, sep="\t")
```

---

!summary_embed_start!


``` r
paste("Last updated on", format(Sys.time(), '%d %B, %Y'))
```

```
## "Last updated on 26 September, 2026"
```

#### Record statistics by species and TCR chain

General statistics. Note that general statistics are computed using the 'slim' database version. Slim version, for example, lists the same TCR sequence found in several donors/studies only once and selects representative V/J for a given CDR3aa clonotype and the best score across all redundant records.


``` r
df.sg = df %>% 
  group_by(species, gene) %>%
  summarize(records = length(complex.id), 
            paired.records = sum(ifelse(complex.id=="0", 0, 1)),
            epitopes = length(unique(antigen.epitope)),
            publications = length(unique(str_split_fixed(reference.id, ",", n = Inf)[,1]))) %>%
  arrange(species, gene)

colnames(df.sg) = c("Species", "Chain", "Records", "Paired records", "Unique epitopes", "Studies")

kable(format = "html", df.sg)
```

<table>
 <thead>
  <tr>
   <th style="text-align:left;"> Species </th>
   <th style="text-align:left;"> Chain </th>
   <th style="text-align:right;"> Records </th>
   <th style="text-align:right;"> Paired records </th>
   <th style="text-align:right;"> Unique epitopes </th>
   <th style="text-align:right;"> Studies </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> HomoSapiens </td>
   <td style="text-align:left;"> TRA </td>
   <td style="text-align:right;"> 60822 </td>
   <td style="text-align:right;"> 41148 </td>
   <td style="text-align:right;"> 1747 </td>
   <td style="text-align:right;"> 348 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HomoSapiens </td>
   <td style="text-align:left;"> TRB </td>
   <td style="text-align:right;"> 120530 </td>
   <td style="text-align:right;"> 71524 </td>
   <td style="text-align:right;"> 1937 </td>
   <td style="text-align:right;"> 431 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MacacaMulatta </td>
   <td style="text-align:left;"> TRB </td>
   <td style="text-align:right;"> 1290 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MusMusculus </td>
   <td style="text-align:left;"> TRA </td>
   <td style="text-align:right;"> 6843 </td>
   <td style="text-align:right;"> 4900 </td>
   <td style="text-align:right;"> 89 </td>
   <td style="text-align:right;"> 71 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MusMusculus </td>
   <td style="text-align:left;"> TRB </td>
   <td style="text-align:right;"> 8185 </td>
   <td style="text-align:right;"> 4742 </td>
   <td style="text-align:right;"> 124 </td>
   <td style="text-align:right;"> 88 </td>
  </tr>
</tbody>
</table>

#### Record statistics by year



Number of unique TCR sequences, epitopes, MHC alleles and studies by **publication** year, cumulative plots.


``` r
# Cumulative counts by first-appearance year.
#
# The previous form cross-joined every (pub_date, pub_date2) pair -- roughly 7M rows over ~35
# distinct years -- and then ran four `length(unique(x[which(pub_date <= pub_date2)]))` per group.
# That is O(n x years) and the only part of this render that could plausibly exhaust a 16 GB
# runner. A key's first year, tabulated and cumulated, gives the same numbers in one pass.
years = sort(unique(dt.pubdate$pub_date))

cumulative = dt.vdjdb.s %>%
  inner_join(dt.pubdate, by = "reference.id", relationship = "many-to-many") %>%
  select(chains, pub_date, tcr = tcr_key, epi = antigen.epitope,
         ref = `reference.id`, mhc = mhc_key) %>%
  pivot_longer(c(tcr, epi, ref, mhc), names_to = "metric", values_to = "key") %>%
  filter(key != "") %>%
  group_by(chains, metric, key) %>%
  summarise(first_year = min(pub_date), .groups = "drop") %>%
  count(chains, metric, first_year, name = "added") %>%
  complete(chains, metric, first_year = years, fill = list(added = 0)) %>%
  arrange(chains, metric, first_year) %>%
  group_by(chains, metric) %>%
  mutate(total = cumsum(added)) %>%
  ungroup()

# Callouts carry no coordinates (#460). The segment ends at the series value in that year and the
# label floats a fixed fraction of the panel maximum above it, so they stay correct as the database
# grows -- the hardcoded y positions they replace were right for a 2022 database and wrong since.
dt.ann = fread(params$annotations, sep = "\t", header = TRUE) %>%
  mutate(label = gsub("\\\\n", "\n", label))
LABEL_FLOAT = 0.05

callouts = function(metric_name) {
  d = cumulative %>% filter(metric == metric_name)
  a = dt.ann %>% filter(panel == metric_name)
  if (nrow(a) == 0 || nrow(d) == 0) return(NULL)
  top = max(d$total)
  a$value = sapply(a$year, function(y) {
    v = d$total[d$first_year == y]
    if (length(v) == 0) 0 else max(v)
  })
  list(annotate("segment", x = a$year, xend = a$year, y = 0, yend = a$value,
                linetype = "solid", colour = "grey25", linewidth = 0.3),
       annotate("text", x = a$year, y = a$value + LABEL_FLOAT * top, label = a$label,
                hjust = a$hjust, vjust = a$vjust))
}

panel = function(metric_name, title, legend_title = NULL) {
  d = cumulative %>% filter(metric == metric_name)
  p = ggplot(d, aes(x = first_year, y = total, color = chains)) +
    callouts(metric_name) +
    geom_line() +
    geom_point() +
    ylab("") +
    scale_x_continuous("", breaks = seq(1995, max(years), by = 2)) +
    theme_classic() + ggtitle(title) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          axis.line = element_line(linewidth = 0.3))
  if (is.null(legend_title))
    p + scale_color_brewer(palette = "Set1")
  else
    p + scale_color_brewer(legend_title, palette = "Set1") + theme(legend.position = "bottom")
}

p1 = panel("tcr", "Number of unique TCRs")
p2 = panel("epi", "Number of unique epitopes")
p3 = panel("ref", "Number of studies")
p4 = panel("mhc", "Number of MHC alleles", "TCR chain(s)")

# `cowplot::get_legend`, not a grep over grob names. The hand-rolled version looked for a grob
# named exactly "guide-box"; that name still matches on ggplot2 4.0.2, but it is an internal detail
# of the gtable and `tmp$grobs[[integer(0)]]` is the error it throws the release it stops matching.
mylegend = cowplot::get_legend(p4)

grid.arrange(arrangeGrob(p1 + theme(legend.position="none"),
                         p2 + theme(legend.position="none"),
                         p3 + theme(legend.position="none"),
                         p4 + theme(legend.position="none"),
                         nrow=2),
             mylegend, nrow=2,heights=c(10, 1)) #-> PX1
```

<img src="vdjdb_summary_files/figure-html/unnamed-chunk-5-1.png" alt="" width="672" />

``` r
#PX1 %>% plot
#pdf("pubyear.pdf")
#PX1 %>% plot
#dev.off()

#fwrite(dt.vdjdb.s2, "vdjdb_stats_pubyear.txt", sep = "\t")
```

#### Summary by antigen and antigen origin

Representative data for Homo Sapiens


``` r
df.a = df %>% 
  filter(species == "HomoSapiens") %>%
  group_by(antigen.species) %>% #, antigen.gene) %>%
  summarize(records = n(), 
            epitopes = length(unique(antigen.epitope)),
            publications = length(unique(str_split_fixed(reference.id, ",", n = Inf)[,1]))) %>%
  arrange(-records)

colnames(df.a) = c("Parent species", #"Parent gene",
                   "Records", "Unique epitopes", "Studies")
kable(format = "html", df.a)
```

<table>
 <thead>
  <tr>
   <th style="text-align:left;"> Parent species </th>
   <th style="text-align:right;"> Records </th>
   <th style="text-align:right;"> Unique epitopes </th>
   <th style="text-align:right;"> Studies </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> CMV </td>
   <td style="text-align:right;"> 54504 </td>
   <td style="text-align:right;"> 154 </td>
   <td style="text-align:right;"> 72 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HomoSapiens </td>
   <td style="text-align:right;"> 41033 </td>
   <td style="text-align:right;"> 718 </td>
   <td style="text-align:right;"> 183 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> EBV </td>
   <td style="text-align:right;"> 36885 </td>
   <td style="text-align:right;"> 56 </td>
   <td style="text-align:right;"> 53 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> InfluenzaA </td>
   <td style="text-align:right;"> 16768 </td>
   <td style="text-align:right;"> 54 </td>
   <td style="text-align:right;"> 41 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SARS-CoV-2 </td>
   <td style="text-align:right;"> 15373 </td>
   <td style="text-align:right;"> 693 </td>
   <td style="text-align:right;"> 44 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HIV-1 </td>
   <td style="text-align:right;"> 3549 </td>
   <td style="text-align:right;"> 89 </td>
   <td style="text-align:right;"> 43 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HBV </td>
   <td style="text-align:right;"> 2806 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 7 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> M.tuberculosis </td>
   <td style="text-align:right;"> 2672 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 7 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> YFV </td>
   <td style="text-align:right;"> 2646 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 8 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HCV </td>
   <td style="text-align:right;"> 1975 </td>
   <td style="text-align:right;"> 25 </td>
   <td style="text-align:right;"> 14 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HTLV-1 </td>
   <td style="text-align:right;"> 571 </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 12 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> DENV </td>
   <td style="text-align:right;"> 419 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> InfluenzaB </td>
   <td style="text-align:right;"> 393 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Wheat </td>
   <td style="text-align:right;"> 329 </td>
   <td style="text-align:right;"> 23 </td>
   <td style="text-align:right;"> 16 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PlasmodiumFalciparum </td>
   <td style="text-align:right;"> 282 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MCPyV </td>
   <td style="text-align:right;"> 274 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SARS-CoV </td>
   <td style="text-align:right;"> 147 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RotavirusA </td>
   <td style="text-align:right;"> 84 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HAdV2 </td>
   <td style="text-align:right;"> 71 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Synthetic </td>
   <td style="text-align:right;"> 70 </td>
   <td style="text-align:right;"> 33 </td>
   <td style="text-align:right;"> 7 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CoxsackievirusB </td>
   <td style="text-align:right;"> 65 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HCoV-OC43 </td>
   <td style="text-align:right;"> 63 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> VZV </td>
   <td style="text-align:right;"> 61 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HCoV-HKU1 </td>
   <td style="text-align:right;"> 48 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> E.Coli </td>
   <td style="text-align:right;"> 41 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 5 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HSV-2 </td>
   <td style="text-align:right;"> 34 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HPV </td>
   <td style="text-align:right;"> 28 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> KlebsiellaOxytoca </td>
   <td style="text-align:right;"> 22 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HPV16 </td>
   <td style="text-align:right;"> 20 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 4 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HPV18 </td>
   <td style="text-align:right;"> 20 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AdV </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HAdV5 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MusMusculus </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 4 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SalmonellaTyphi </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CryptomeriaJaponica </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Unknown </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 4 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HSV-1 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AspergillusOryzae </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> StreptomycesKanamyceticus </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> FusariumOxysporum </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> BacillusSubtilis </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CryptococcusNeoformans </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HHV </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HaemophilusInfluenzae </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> M.avium </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PseudomonasAeruginosa </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PseudomonasFluorescens </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SaccharomycesCerevisiae </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> StaphylococcusAureus </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
</tbody>
</table>

---

#### **COVID-19** data

Summary of antigens and T-cell receptors related to COVID-19 pandemic. Number of records for SARS-CoV-2 epitopes grouped by viral protein and HLA plotted using alluvium plot. Epitopes with less than 30 records in total were not counted.


``` r
df %>% 
  filter(species == "HomoSapiens", 
         startsWith(as.character(antigen.species), "SARS-CoV")) %>%
  mutate(mhc.a = str_split_fixed(mhc.a, "[,:]", 2)[,1],
         mhc.b = str_split_fixed(mhc.b, "[,:]", 2)[,1],
         mhc = ifelse(mhc.class == "MHCI", mhc.a, paste0(mhc.a, '/', substr(mhc.b, 7, 15)))) %>%
  group_by(antigen.gene, mhc, antigen.epitope) %>%
  mutate(publications = length(unique(str_split_fixed(reference.id, ",", n = Inf)[,1]))) %>%
  group_by(antigen.gene, mhc, antigen.epitope, gene, publications) %>%
  summarize(records = n()) -> df.c

colnames(df.c) = c("Gene", "HLA", "Epitope", "TCR chain",
                   "Studies",
                   "Records")

ggplot(df.c %>% 
         ungroup %>% arrange(Records) %>% filter(Records >= 30),
       aes(axis1 = Gene,
           axis2 = gsub("HLA-", "", HLA),
           axis3 = substr(Epitope,1,3),
           axis4 = `TCR chain`,
           y = log2(Records))) +
  geom_alluvium(aes(fill = substr(Epitope,1,3) %>% as.factor %>% as.integer), 
                color = "white", alpha = 0.8, curve_type = "sigmoid") +
  geom_stratum(fill = "grey95", color = "white", linewidth =1.0) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_fill_distiller(palette = "Set3", guide="none", "") +
  #scale_fill_hue(guide="none", "") +
  ylab("") + scale_x_discrete(limits = c("Gene", "HLA", "Epitope", "TCR chain"),
                              expand = c(.1, .1),
                              position = "top") +
  theme_void() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x =  element_text(size = 16, color = "black", vjust = -5),
        axis.ticks.x = element_blank(),
        panel.grid.major.y = element_blank())
```

<img src="vdjdb_summary_files/figure-html/unnamed-chunk-7-1.png" alt="" width="1152" />

Summary of SARS-CoV-2 epitopes and corresponding TCR alpha and beta chain specificity records (cases with 10+ records)


``` r
kable(format = "html", 
      df.c %>% 
        reshape2::dcast(Gene + HLA + Epitope + Studies ~ `TCR chain`, fill = 0) %>%
        mutate(HLA = gsub("*", ".", HLA, fixed = T)) %>%
        filter(TRA+TRB >= 10) %>%
        arrange(-(TRB+TRA)))
```

<table>
 <thead>
  <tr>
   <th style="text-align:left;"> Gene </th>
   <th style="text-align:left;"> HLA </th>
   <th style="text-align:left;"> Epitope </th>
   <th style="text-align:right;"> Studies </th>
   <th style="text-align:right;"> TRA </th>
   <th style="text-align:right;"> TRB </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> N </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> LLLDRLNQL </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 452 </td>
   <td style="text-align:right;"> 1717 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> S </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> YLQPRTFLL </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 708 </td>
   <td style="text-align:right;"> 1118 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> S </td>
   <td style="text-align:left;"> HLA-DPA1.01/B1.04 </td>
   <td style="text-align:left;"> TFEYVSQPFLMDLE </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 249 </td>
   <td style="text-align:right;"> 882 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> N </td>
   <td style="text-align:left;"> HLA-B.07 </td>
   <td style="text-align:left;"> SPRWYFYYL </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 503 </td>
   <td style="text-align:right;"> 626 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> S </td>
   <td style="text-align:left;"> HLA-A.24 </td>
   <td style="text-align:left;"> NYNYLYRLF </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 345 </td>
   <td style="text-align:right;"> 709 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ORF1ab </td>
   <td style="text-align:left;"> HLA-A.01 </td>
   <td style="text-align:left;"> TTDPSFLGRY </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 410 </td>
   <td style="text-align:right;"> 466 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> S </td>
   <td style="text-align:left;"> HLA-A.24 </td>
   <td style="text-align:left;"> QYIKWPWYI </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 302 </td>
   <td style="text-align:right;"> 519 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> S </td>
   <td style="text-align:left;"> HLA-B.15 </td>
   <td style="text-align:left;"> NQKLIANQF </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 218 </td>
   <td style="text-align:right;"> 288 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> N </td>
   <td style="text-align:left;"> HLA-A.03 </td>
   <td style="text-align:left;"> KTFPPTEPK </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 128 </td>
   <td style="text-align:right;"> 278 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> S </td>
   <td style="text-align:left;"> HLA-DPA1.01/B1.04 </td>
   <td style="text-align:left;"> RSFIEDLLFNKVTLA </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 167 </td>
   <td style="text-align:right;"> 169 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> S </td>
   <td style="text-align:left;"> HLA-A.01 </td>
   <td style="text-align:left;"> LTDEMIAQY </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 125 </td>
   <td style="text-align:right;"> 129 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> S </td>
   <td style="text-align:left;"> HLA-A.03 </td>
   <td style="text-align:left;"> KCYGVSPTK </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 214 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> N </td>
   <td style="text-align:left;"> HLA-C.12 </td>
   <td style="text-align:left;"> KAYNVTQAF </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 51 </td>
   <td style="text-align:right;"> 162 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ORF1ab </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> ALWEIQQVV </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 27 </td>
   <td style="text-align:right;"> 153 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> N </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> RTATKQYNV </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 50 </td>
   <td style="text-align:right;"> 93 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> S </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> RLQSLQTYV </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 53 </td>
   <td style="text-align:right;"> 90 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> N </td>
   <td style="text-align:left;"> HLA-B.40 </td>
   <td style="text-align:left;"> MEVTPSGTWL </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 47 </td>
   <td style="text-align:right;"> 69 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ORF3a </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> ALSKGVHFV </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 27 </td>
   <td style="text-align:right;"> 77 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ORF3a </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> LLYDANYFL </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 33 </td>
   <td style="text-align:right;"> 70 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ORF3a </td>
   <td style="text-align:left;"> HLA-A.01 </td>
   <td style="text-align:left;"> FTSDYYQLY </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 46 </td>
   <td style="text-align:right;"> 45 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> M </td>
   <td style="text-align:left;"> HLA-A.24 </td>
   <td style="text-align:left;"> RYRIGNYKL </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 44 </td>
   <td style="text-align:right;"> 34 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ORF3a </td>
   <td style="text-align:left;"> HLA-A.24 </td>
   <td style="text-align:left;"> VYFLQSINF </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 29 </td>
   <td style="text-align:right;"> 29 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ORF1ab </td>
   <td style="text-align:left;"> HLA-A.01 </td>
   <td style="text-align:left;"> PTDNYITTY </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 24 </td>
   <td style="text-align:right;"> 26 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> M </td>
   <td style="text-align:left;"> HLA-B.15 </td>
   <td style="text-align:left;"> RVAGDSGFAAY </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 24 </td>
   <td style="text-align:right;"> 24 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ORF1ab </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> KLWAQCVQL </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 21 </td>
   <td style="text-align:right;"> 25 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> N </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> RTATKAYNV </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 19 </td>
   <td style="text-align:right;"> 23 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> M </td>
   <td style="text-align:left;"> HLA-B.07 </td>
   <td style="text-align:left;"> WPVTLACFV </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 21 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ORF1ab </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> KSVNITFEL </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 16 </td>
   <td style="text-align:right;"> 19 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NSP3 </td>
   <td style="text-align:left;"> HLA-A.01 </td>
   <td style="text-align:left;"> PTDNYITTY </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 19 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> S </td>
   <td style="text-align:left;"> HLA-B.07 </td>
   <td style="text-align:left;"> SIIAYTMSL </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 17 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> S </td>
   <td style="text-align:left;"> HLA-DRA.01/B1.04 </td>
   <td style="text-align:left;"> HWFVTQRNFYEPQII </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 16 </td>
   <td style="text-align:right;"> 16 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ORF1ab </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> VMVELVAEL </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 17 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ORF1ab </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> RQLLFVVEV </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 14 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> M </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> WLLWPVTLA </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 16 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> S </td>
   <td style="text-align:left;"> HLA-B.35 </td>
   <td style="text-align:left;"> FVSNGTHWF </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 19 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ORF1ab </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> TLMNVLTLV </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 11 </td>
   <td style="text-align:right;"> 13 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ORF14 </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> LLLEWLAMA </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 13 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ORF1ab </td>
   <td style="text-align:left;"> HLA-A.24 </td>
   <td style="text-align:left;"> TSQWLTNIF </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 11 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ORF9b </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> KVYPIILRL </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 12 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> S </td>
   <td style="text-align:left;"> HLA-B.07 </td>
   <td style="text-align:left;"> WPWYIWLGF </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 14 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> S </td>
   <td style="text-align:left;"> HLA-DRA.01/B1.04 </td>
   <td style="text-align:left;"> NCTFEYVSQPFLMDL </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 16 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> S </td>
   <td style="text-align:left;"> HLA-DRA.01/B1.11 </td>
   <td style="text-align:left;"> VGGNYNYLYRLFRKS </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 10 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ORF1ab </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> CLAVHECFV </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 10 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ORF1ab </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> FLAHIQWMV </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 11 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ORF1ab </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> FLLNKEMYL </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 10 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> N </td>
   <td style="text-align:left;"> HLA-B.44 </td>
   <td style="text-align:left;"> QELIRQGTDY </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 9 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NSP12 </td>
   <td style="text-align:left;"> HLA-A.24 </td>
   <td style="text-align:left;"> VYIGDPAQL </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 9 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ORF1ab </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> AIFYLITPV </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 11 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ORF1ab </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> FTVLCLTPV </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 11 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> S </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> VFLVLLPLV </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 10 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ORF1ab </td>
   <td style="text-align:left;"> HLA-A.24 </td>
   <td style="text-align:left;"> NYMPYFFTL </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 9 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ORF3a </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> YLYALVYFL </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 10 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> S </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> RLITGRLQSL </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 8 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ORF1ab </td>
   <td style="text-align:left;"> HLA-A.01 </td>
   <td style="text-align:left;"> NSSTCMMCY </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 10 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ORF1ab </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> KLQFTSLEI </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 8 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ORF1ab </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> YVDNSSLTI </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 9 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ORF1ab </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> KLKDCVMYA </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 8 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ORF1ab </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> YVWKSYVHV </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 9 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ORF1ab </td>
   <td style="text-align:left;"> HLA-B.07 </td>
   <td style="text-align:left;"> MPASWVMRI </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 10 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ORF1ab </td>
   <td style="text-align:left;"> HLA-B.07 </td>
   <td style="text-align:left;"> MPYFFTLLL </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 8 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ORF7a </td>
   <td style="text-align:left;"> HLA-B.07 </td>
   <td style="text-align:left;"> SPIFRIVAA </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 11 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> N </td>
   <td style="text-align:left;"> HLA-A.01 </td>
   <td style="text-align:left;"> NTNSSPDDQIGYY </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 7 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ORF10 </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> NVFAFPFTI </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 7 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ORF1ab </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> FLALCADSI </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 7 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ORF1ab </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> FLPRVFSAV </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 8 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ORF1ab </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> KLNEEIAII </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 8 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ORF1ab </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> SVLYYQNNV </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 8 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ORF1ab </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> YMPYFFTLL </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 7 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ORF1ab </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> YTMADLVYA </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 9 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RNP </td>
   <td style="text-align:left;"> HLA-A.01 </td>
   <td style="text-align:left;"> DTDFVNEFY </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 8 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> S </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> IMLCCMTSC </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 8 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> S </td>
   <td style="text-align:left;"> HLA-B.44 </td>
   <td style="text-align:left;"> AEVQIDRLI </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 7 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> S </td>
   <td style="text-align:left;"> HLA-DRA.01/B1.03 </td>
   <td style="text-align:left;"> RISNCVADYSVLYNS </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 7 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> N </td>
   <td style="text-align:left;"> HLA-B.07 </td>
   <td style="text-align:left;"> RIRGGDGKM </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 8 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ORF1ab </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> FLPGVYSVI </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 8 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ORF1ab </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> LLDDFVEII </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 7 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ORF1ab </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> SLLMPILTL </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 7 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> S </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> KQIYKTPPI </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 7 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> S </td>
   <td style="text-align:left;"> HLA-DRA.01/B1.15 </td>
   <td style="text-align:left;"> NLLLQYGSFCTQLNRAL </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 13 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> M </td>
   <td style="text-align:left;"> HLA-A.24 </td>
   <td style="text-align:left;"> LWLLWPVTL </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 6 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> N </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> LLNKHIDAY </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 8 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ORF1ab </td>
   <td style="text-align:left;"> HLA-A.01 </td>
   <td style="text-align:left;"> LMNVLTLVY </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 6 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ORF1ab </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> KLSYGIATV </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 8 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ORF1ab </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> LMIERFVSL </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 7 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ORF1ab </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> RIMTWLDMV </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 8 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ORF1ab </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> SMMILSDDA </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 6 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ORF1ab </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> YLFDESGEF </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 7 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ORF1ab </td>
   <td style="text-align:left;"> HLA-B.07 </td>
   <td style="text-align:left;"> SPYNSQNAV </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 8 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> S </td>
   <td style="text-align:left;"> HLA-DRA.01/B1.03 </td>
   <td style="text-align:left;"> NFSQILPDPSKPSKR </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 6 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> M </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> LACFVLAAV </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 7 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ORF1ab </td>
   <td style="text-align:left;"> HLA-A.01 </td>
   <td style="text-align:left;"> YADVFHLYL </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 6 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ORF1ab </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> FLNRFTTTL </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 6 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ORF1ab </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> TLMNVITLV </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 6 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ORF1ab </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> TTIQTIVEV </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 6 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ORF1ab </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> YADVFHLYL </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 6 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ORF1ab </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> YLGGMSYYC </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 6 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ORF1ab </td>
   <td style="text-align:left;"> HLA-B.07 </td>
   <td style="text-align:left;"> VPYFNMVYM </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 7 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> S </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> FIAGLIAIV </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 7 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ORF1ab </td>
   <td style="text-align:left;"> HLA-A.01 </td>
   <td style="text-align:left;"> FLTENLLLY </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 6 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ORF1ab </td>
   <td style="text-align:left;"> HLA-A.01 </td>
   <td style="text-align:left;"> TSAMQTMLF </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 6 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ORF1ab </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> LMCQPILLL </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 6 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ORF1ab </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> MLDMYSVML </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 5 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ORF1ab </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> SMWALVISV </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 5 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ORF1ab </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> TLKNTVCTV </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 5 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ORF1ab </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> VYIGDPAQL </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 5 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ORF1ab </td>
   <td style="text-align:left;"> HLA-A.24 </td>
   <td style="text-align:left;"> VMHANYIFW </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 5 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ORF7a </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> KLFIRQEEV </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 5 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ORF7b </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> LLFLVLIML </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 5 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> S </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> FLHVTYVPA </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 5 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> S </td>
   <td style="text-align:left;"> HLA-DQA1.01/B1.06 </td>
   <td style="text-align:left;"> NNSYECDIPIGAGIC </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 5 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> S </td>
   <td style="text-align:left;"> HLA-DRA.01/B1.03 </td>
   <td style="text-align:left;"> VGGNYNYLYRLFRKSNLKP </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 10 </td>
  </tr>
</tbody>
</table>

---

#### **Self-antigen** data

Summary of T-cell receptors recognizing self-antigens, including antigens linked to utoimmune diseases and potential neoantigen targets for cancer immunotherapy. Number of records for self-antigens grouped by (mutated) human gene and corresponding HLAs are plotted using alluvium plot. Only self-antigens with at least 10 records are shown.


``` r
df %>% 
  filter(species == "HomoSapiens", 
         startsWith(as.character(antigen.species), "HomoSapiens")) %>%
  mutate(mhc.a = str_split_fixed(mhc.a, "[,:]", 2)[,1],
         mhc.b = str_split_fixed(mhc.b, "[,:]", 2)[,1],
         mhc = ifelse(mhc.class == "MHCI", mhc.a, paste0(mhc.a, '/', substr(mhc.b, 7, 15)))) %>%
  group_by(antigen.gene, mhc, antigen.epitope) %>%
  mutate(publications = length(unique(str_split_fixed(reference.id, ",", n = Inf)[,1]))) %>%
  group_by(antigen.gene, mhc, antigen.epitope, gene, publications) %>%
  summarize(records = n()) -> df.n

colnames(df.n) = c("Gene", "HLA", "Epitope", "TCR chain",
                   "Studies",
                   "Records")

ggplot(df.n %>% 
         ungroup %>% arrange(Records) %>% filter(Records >= 10),
       aes(axis1 = Gene,
           axis2 = gsub("HLA-", "", HLA),
           axis3 = substr(Epitope,1,3),
           axis4 = `TCR chain`,
           y = log2(Records))) +
  geom_alluvium(aes(fill = substr(Epitope,1,3) %>% as.factor %>% as.integer), 
                color = "white", alpha = 0.8, curve_type = "sigmoid") +
  geom_stratum(fill = "grey95", color = "white", linewidth =1.0) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_fill_distiller(palette = "Accent", guide="none", "") +
  #scale_fill_hue(guide="none", "") +
  ylab("") + scale_x_discrete(limits = c("Gene", "HLA", "Epitope", "TCR chain"),
                              expand = c(.1, .1),
                              position = "top") +
  theme_void() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x =  element_text(size = 16, color = "black", vjust = -5),
        axis.ticks.x = element_blank(),
        panel.grid.major.y = element_blank())
```

<img src="vdjdb_summary_files/figure-html/unnamed-chunk-9-1.png" alt="" width="1152" />

Summary of neoantigens and corresponding TCR alpha and beta chain specificity records (cases with 10+ records)


``` r
kable(format = "html", 
      df.n %>% 
        reshape2::dcast(Gene + HLA + Epitope + Studies ~ `TCR chain`, fill = 0) %>%
        mutate(HLA = gsub("*", ".", HLA, fixed = T)) %>%
        filter(TRA+TRB >= 10) %>%
        arrange(-(TRB+TRA)))
```

<table>
 <thead>
  <tr>
   <th style="text-align:left;"> Gene </th>
   <th style="text-align:left;"> HLA </th>
   <th style="text-align:left;"> Epitope </th>
   <th style="text-align:right;"> Studies </th>
   <th style="text-align:right;"> TRA </th>
   <th style="text-align:right;"> TRB </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> NY-ESO-1 </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> SLLMWITQV </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 38 </td>
   <td style="text-align:right;"> 29698 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MLANA </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> ELAGIGILTV </td>
   <td style="text-align:right;"> 19 </td>
   <td style="text-align:right;"> 411 </td>
   <td style="text-align:right;"> 1478 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PGT </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> LLAGIGTVPI </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 479 </td>
   <td style="text-align:right;"> 885 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> BST2 </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> LLLGIGILV </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 384 </td>
   <td style="text-align:right;"> 233 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> KLK3 </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> VISNDVCAQV </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 229 </td>
   <td style="text-align:right;"> 362 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MLANA </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> ALAGIGILTV </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 97 </td>
   <td style="text-align:right;"> 150 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IGF2BP2 </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> NLSALGIFST </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 34 </td>
   <td style="text-align:right;"> 183 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOB </td>
   <td style="text-align:left;"> HLA-DRA.01/B1.07 </td>
   <td style="text-align:left;"> SLFFSAQPFEITAST </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 104 </td>
   <td style="text-align:right;"> 104 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MLANA </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> EAAGIGILTV </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 47 </td>
   <td style="text-align:right;"> 106 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HCRT </td>
   <td style="text-align:left;"> HLA-DQA1.05/B1.06 </td>
   <td style="text-align:left;"> HGAGNHAAGILTL </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 73 </td>
   <td style="text-align:right;"> 70 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> KIF20A </td>
   <td style="text-align:left;"> HLA-DQA1.01/B1.06 </td>
   <td style="text-align:left;"> GTRVIRDMTLHSAPS </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 63 </td>
   <td style="text-align:right;"> 65 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IGRP </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> VLFGLGFAI </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 51 </td>
   <td style="text-align:right;"> 70 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GNL3L </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> NLNCCSVPV </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 48 </td>
   <td style="text-align:right;"> 63 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HCRT </td>
   <td style="text-align:left;"> HLA-DQA1.05/B1.06 </td>
   <td style="text-align:left;"> ASGNHAAGILTM </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 56 </td>
   <td style="text-align:right;"> 51 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TKT </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> AMFWSVPTV </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 95 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PMEL </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> IMDQVPFSV </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 31 </td>
   <td style="text-align:right;"> 62 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MLANA </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> ALGIGILTV </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 40 </td>
   <td style="text-align:right;"> 52 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SEC24A </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> FLYNLLTRV </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 78 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SF3B1 </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> RLPGVLPRA </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 43 </td>
   <td style="text-align:right;"> 48 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SLC30A8 </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> VAANIVLTV </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 43 </td>
   <td style="text-align:right;"> 48 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> BST2 </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> LLLGIGILVL </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 23 </td>
   <td style="text-align:right;"> 59 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> FNDC3B </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> VVLSWAPPV </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 39 </td>
   <td style="text-align:right;"> 42 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NUF2 </td>
   <td style="text-align:left;"> HLA-A.24 </td>
   <td style="text-align:left;"> VYGIRLEHF </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 75 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> FNDC3B </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> VVMSWAPPV </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 37 </td>
   <td style="text-align:right;"> 43 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LY6K </td>
   <td style="text-align:left;"> HLA-A.24 </td>
   <td style="text-align:left;"> RYCNLEGPPI </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 71 </td>
   <td style="text-align:right;"> 8 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NSDHL </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> ILTGLNYEV </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 35 </td>
   <td style="text-align:right;"> 38 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> G6PC2 </td>
   <td style="text-align:left;"> HLA-A.01 </td>
   <td style="text-align:left;"> LTSLTILQLY </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 34 </td>
   <td style="text-align:right;"> 36 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> INS </td>
   <td style="text-align:left;"> HLA-DRA.01/B1.04 </td>
   <td style="text-align:left;"> GIVEQCCTSICSLYQ </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 69 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ZNT8 </td>
   <td style="text-align:left;"> HLA-B.08 </td>
   <td style="text-align:left;"> YAKWKLCSA </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 28 </td>
   <td style="text-align:right;"> 41 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PLA2G6 </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> FLASKIGRLV </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 34 </td>
   <td style="text-align:right;"> 34 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ZNT8 </td>
   <td style="text-align:left;"> HLA-B.08 </td>
   <td style="text-align:left;"> YAYAKWKL </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 31 </td>
   <td style="text-align:right;"> 32 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ZNT8 </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> LLSLFSLWL </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 24 </td>
   <td style="text-align:right;"> 38 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NSDHL </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> ILTGLNYEA </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 27 </td>
   <td style="text-align:right;"> 29 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ZNT8 </td>
   <td style="text-align:left;"> HLA-A.01 </td>
   <td style="text-align:left;"> VTDAAHLLI </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 23 </td>
   <td style="text-align:right;"> 33 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GANAB </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> ALYGFVPVL </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 24 </td>
   <td style="text-align:right;"> 27 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ZNT8 </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> VVTGVLVYL </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 23 </td>
   <td style="text-align:right;"> 27 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> KIF20A </td>
   <td style="text-align:left;"> HLA-DQA1.01/B1.06 </td>
   <td style="text-align:left;"> LHCQERANELMRAMK </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 21 </td>
   <td style="text-align:right;"> 23 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IAPP </td>
   <td style="text-align:left;"> HLA-A.01 </td>
   <td style="text-align:left;"> STNVGSNTY </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 16 </td>
   <td style="text-align:right;"> 27 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> INS </td>
   <td style="text-align:left;"> HLA-DQA1.03/B1.02 </td>
   <td style="text-align:left;"> GQVELGGGNAVEVCK </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 21 </td>
   <td style="text-align:right;"> 21 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PMEL </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> KTWGQYWQV </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 21 </td>
   <td style="text-align:right;"> 21 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GAD65 </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> VMNILLQYV </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 24 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GFAP </td>
   <td style="text-align:left;"> HLA-B.08 </td>
   <td style="text-align:left;"> HLKRNIVV </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 16 </td>
   <td style="text-align:right;"> 23 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> G6PC2 </td>
   <td style="text-align:left;"> HLA-B.08 </td>
   <td style="text-align:left;"> NLIFKWIL </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 19 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> INS </td>
   <td style="text-align:left;"> HLA-B.08 </td>
   <td style="text-align:left;"> MALWMRLL </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 20 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PORCN </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> LLHGFSFYL </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 36 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GAD2 </td>
   <td style="text-align:left;"> HLA-B.08 </td>
   <td style="text-align:left;"> MMIARFKM </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 21 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MLANA </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> ELAGIGLTV </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 22 </td>
   <td style="text-align:right;"> 13 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IGRP </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> FLWSVFWLI </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 22 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PMEL </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> YLEPGPVTV </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 32 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PREINS </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> ALWMRLLPL </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 18 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> KMT2D </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> ALSPVIPHI </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 18 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ZNT8 </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> IQATVMIIV </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 20 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> WT1 </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> RMFPNAPYL </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 17 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> WT1 </td>
   <td style="text-align:left;"> HLA-A.24 </td>
   <td style="text-align:left;"> CYTWNQMNL </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 15 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NY-ESO-1 </td>
   <td style="text-align:left;"> HLA-DQA1.01/B1.07 </td>
   <td style="text-align:left;"> LLEFYLAMPFATP </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 16 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PMEL </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> ITDQVPFSV </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 17 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SEC24A </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> FLYNPLTRV </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 16 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ZNT8 </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> KMYAFTLES </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 20 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AKAP13 </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> KLMNIQQKL </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 25 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HAUS3 </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> ILNAMIAKI </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 14 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> KIF20A </td>
   <td style="text-align:left;"> HLA-DQA1.01/B1.06 </td>
   <td style="text-align:left;"> QASRTVIHSADITFQ </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 10 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PREINS </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> RLLPLLALL </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 15 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> F8 </td>
   <td style="text-align:left;"> HLA-DRA.01/B1.01 </td>
   <td style="text-align:left;"> SYFTNMFATWSPSKARLHLQ </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 26 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MYD88 </td>
   <td style="text-align:left;"> HLA-B.07 </td>
   <td style="text-align:left;"> RPIPIKYKAM </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 13 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PGM5 </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> AVGSYVYSV </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 12 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CTAG1B </td>
   <td style="text-align:left;"> HLA-DRA.01/B3.02 </td>
   <td style="text-align:left;"> LKEFTVSGNILTIRL </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 13 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> G6PC2 </td>
   <td style="text-align:left;"> HLA-B.08 </td>
   <td style="text-align:left;"> YLKTNLFL </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 15 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GAD2 </td>
   <td style="text-align:left;"> HLA-A.01 </td>
   <td style="text-align:left;"> VSATAGTTVY </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 10 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IAPP </td>
   <td style="text-align:left;"> HLA-B.08 </td>
   <td style="text-align:left;"> ILKLQVFL </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 15 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> INS </td>
   <td style="text-align:left;"> HLA-A.01 </td>
   <td style="text-align:left;"> GSHLVEALY </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 14 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> INS </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> RLLPLLALLAL </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 17 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CLK3 </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> RLWGTWVKA </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 11 </td>
   <td style="text-align:right;"> 11 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ZDBF2 </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> YILKYSVFL </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 21 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IAPP </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> KLQVFLIVL </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 14 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> INS </td>
   <td style="text-align:left;"> HLA-A.24 </td>
   <td style="text-align:left;"> LWMRLLPLL </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 12 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PHKA2 </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> LLSIIFFPA </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 10 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SRPX </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> TLWCSPIKV </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 10 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NY-ESO-1 </td>
   <td style="text-align:left;"> HLA-DRA.01/B1.04 </td>
   <td style="text-align:left;"> PGVLLKEFTVSGNIL </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 10 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ZNT8 </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> ILSAHVATA </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 9 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PMEL </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> YLEPGPVTA </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 9 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PTPRN </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> LLPPLLEHL </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 9 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GANAB </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> ALYGSVPVL </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 7 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IGRP </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> LNIDLLWSV </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 8 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IGRP </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> NLFLFLFAV </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 11 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MAGEA1 </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> KVLEYVIKV </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 8 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MLANA </td>
   <td style="text-align:left;"> HLA-DRA.01/B3.02 </td>
   <td style="text-align:left;"> EPVVNAPPAYEKLS </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 8 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> WDR46 </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> FLIYLDVSV </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 8 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ZNT8 </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> LAVDGVLSV </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 7 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ADH </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> RQFGPDWIVA </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 9 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GFAP </td>
   <td style="text-align:left;"> HLA-B.08 </td>
   <td style="text-align:left;"> LRLRLDQL </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 9 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> INS </td>
   <td style="text-align:left;"> HLA-B.08 </td>
   <td style="text-align:left;"> LWMRLLPL </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 8 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MLANA </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> AAGIGILTV </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 10 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NY-ESO-1 </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> SLLMWITQC </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 7 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Plastin-2 </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> NLFNRYPAL </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 8 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ZNT8 </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> LLSILCIWV </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 10 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 5T4 </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> RLARLALVL </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 7 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> BCL2L1 </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> YLNDHLEPWI </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 7 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> INS </td>
   <td style="text-align:left;"> HLA-B.08 </td>
   <td style="text-align:left;"> WMRLLPLL </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 8 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> KMT2D </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> ALSPVIPLI </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 8 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PABPC1 </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> MLGEQLFPL </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 13 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SLC3A2 </td>
   <td style="text-align:left;"> HLA-DRA.01/B1.07 </td>
   <td style="text-align:left;"> DPPALASTNAEVT </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 7 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TERT </td>
   <td style="text-align:left;"> HLA-DRA.01/B3.02 </td>
   <td style="text-align:left;"> GTAFVQMPAHGLFPW </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 7 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TYR </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> CLLWSFQTSA </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 7 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AKAP9 </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> RLSDFSEQL </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 12 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GAD2 </td>
   <td style="text-align:left;"> HLA-B.08 </td>
   <td style="text-align:left;"> EAKQKGFVPF </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 8 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GAD2 </td>
   <td style="text-align:left;"> HLA-B.08 </td>
   <td style="text-align:left;"> TLKKMREI </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 8 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GCN1L1 </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> SLLRSLENV </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 13 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MRM1 </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> LLFGMPPCL </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 7 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PTPRN </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> VIVMLTPLV </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 8 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PTPRN </td>
   <td style="text-align:left;"> HLA-B.08 </td>
   <td style="text-align:left;"> HARIKLKV </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 7 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SMARCD3 </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> KLFEFLVYGV </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 7 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ZNT8 </td>
   <td style="text-align:left;"> HLA-A.01 </td>
   <td style="text-align:left;"> PTEKGANEY </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 7 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> EXOC8 </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> IILVAVPHV </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 12 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GAD2 </td>
   <td style="text-align:left;"> HLA-B.08 </td>
   <td style="text-align:left;"> SRKHKWKL </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 7 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GFAP </td>
   <td style="text-align:left;"> HLA-A.01 </td>
   <td style="text-align:left;"> ALDIEIATY </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 6 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TLR1 </td>
   <td style="text-align:left;"> HLA-A.24 </td>
   <td style="text-align:left;"> SYLDLPWYL </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 6 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GPER </td>
   <td style="text-align:left;"> HLA-B.27 </td>
   <td style="text-align:left;"> GQMWLLAPR </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 4 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HAUS3 </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> ILNAMITKI </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 5 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> INS </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> LALWGPDPAA </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 7 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PTPRN </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> MVWESGCTV </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 8 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TBX3 </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> GMGPLLATV </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 6 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> USP28 </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> LIIPFIHLI </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 6 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CDK4 </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> ALDPHSGHFV </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 5 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GAD2 </td>
   <td style="text-align:left;"> HLA-A.01 </td>
   <td style="text-align:left;"> QQDKHYDLSY </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 6 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GAD2 </td>
   <td style="text-align:left;"> HLA-B.08 </td>
   <td style="text-align:left;"> FQQDKHYDL </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 4 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IAPP </td>
   <td style="text-align:left;"> HLA-B.08 </td>
   <td style="text-align:left;"> LNHLKATPI </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 6 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> INS </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> ALWGPDPAAA </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 6 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> INSDRIP </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> MLYQHLLPL </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 6 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NBPF14 </td>
   <td style="text-align:left;"> HLA-A.24 </td>
   <td style="text-align:left;"> SYKSYSSTF </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 5 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NDC1 </td>
   <td style="text-align:left;"> HLA-A.02 </td>
   <td style="text-align:left;"> CLNEYHLFL </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 10 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NRAS </td>
   <td style="text-align:left;"> HLA-A.11 </td>
   <td style="text-align:left;"> VVVGADGVGK </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 6 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NY-ESO-1 </td>
   <td style="text-align:left;"> HLA-B.07 </td>
   <td style="text-align:left;"> APRGPHGGAASGL </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 5 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PRPF3 </td>
   <td style="text-align:left;"> HLA-B.27 </td>
   <td style="text-align:left;"> TRLALIAPK </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 5 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RNASEH2B </td>
   <td style="text-align:left;"> HLA-B.27 </td>
   <td style="text-align:left;"> GQVMVVAPR </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 4 </td>
  </tr>
</tbody>
</table>

---

#### Distribution VDJdb confidence scores

Legend: 0 - low confidence, 1 - medium confidence, 2 - high confidence, 3 - very high confidence.

Note that this scoring system is currently deprecated due to large amounts of deep multimer+ repertoire sequencing data and 10X scRNAseq with ImmuDEX dextramer multiplex that are almost impossible to validate independently. We suggest applying methods like TCRNET to [infer high-confidence TCR motifs](https://github.com/antigenomics/vdjdb-motifs) from VDJdb.


``` r
df.score <- df[df$species=='HomoSapiens',] %>%
  group_by(mhc.class, gene, vdjdb.score) %>%
  summarize(total = n())

ggplot(df.score, aes(x=paste(mhc.class, gene, sep = " "), y=total, fill = as.factor(vdjdb.score))) + 
  geom_bar(stat = "identity", position = "dodge", color = "black", linewidth = 0.3) +  
  xlab("") + scale_y_log10("Records") +
  scale_fill_brewer("VDJdb score", palette = "PuBuGn") + 
  theme_classic() +
  theme(legend.position="bottom",
        axis.line = element_line(linewidth = 0.3))
```

<img src="vdjdb_summary_files/figure-html/unnamed-chunk-11-1.png" alt="" width="672" />

---

#### Spectratype

Representative data for Homo Sapiens. CDR3 length distribution (spectratype) is colored by cognate epitope. Second plot shows epitope length distribution for MHC class I and II colored by unique CDR3 (alpha or beta) records.


``` r
df.spe = subset(df, species=="HomoSapiens")

ggplot(df.spe %>% mutate(
  epi_len = nchar(antigen.epitope),
  antigen.epitope = as.factor(antigen.epitope)), 
  aes(x=nchar(as.character(cdr3)))) + 
  geom_histogram(aes(fill = antigen.epitope %>% 
                       fct_reorder(epi_len) %>% as.integer(),
                     group = antigen.epitope %>% 
                       fct_reorder(epi_len)), 
                 bins = 21, size = 0,
                 binwidth = 1, alpha = 0.9, color = NA) + 
  geom_density(aes(y = after_stat(count)), adjust = 3.0, linetype = "dotted") +
  scale_x_continuous(limits = c(5,25), breaks = seq(5,25,5)) + 
  facet_wrap(~gene) + 
  scale_fill_distiller(palette = "Spectral", guide="none", "") +
  #scale_fill_viridis_d(guide="none", direction = -1) +
  xlab("CDR3 length") + ylab("Records") +
  theme_classic() +
  theme(axis.line = element_line(linewidth = 0.3),
        strip.background = element_blank())
```

<img src="vdjdb_summary_files/figure-html/unnamed-chunk-12-1.png" alt="" width="672" />

``` r
ggplot(df.spe %>% mutate(
  cdr3_len = nchar(cdr3),
  cdr3 = as.factor(cdr3)), 
  aes(x=nchar(antigen.epitope))) + 
  geom_histogram(aes(fill = cdr3 %>% 
                       fct_reorder(cdr3_len) %>% as.integer(),
                     group = cdr3 %>% 
                       fct_reorder(cdr3_len)), 
                 size = 0,
                 binwidth = 1, alpha = 0.9, color = NA) + 
  scale_x_continuous(breaks = 7:25) + 
  facet_wrap(.~mhc.class, scales = "free") + 
  scale_fill_distiller(palette = "Spectral", guide="none", "") +
  #scale_fill_viridis_d(guide="none", direction = -1) +
  xlab("Epitope length") + ylab("Records") +
  theme_classic() +
  theme(axis.line = element_line(linewidth = 0.3),
        strip.background = element_blank())
```

<img src="vdjdb_summary_files/figure-html/unnamed-chunk-12-2.png" alt="" width="672" />

#### V(D)J usage and MHC alleles

Representative data for Homo Sapiens, Variable gene. Only Variable genes and MHC alleles with at least 10 records are shown.


``` r
df.vhla = df %>% 
  as_tibble %>%
  filter(species == "HomoSapiens" ) %>%
  mutate(id = 1:n()) %>%
  separate_rows(mhc.a, sep = ",") %>%
  separate_rows(mhc.b, sep = ",") %>%
  separate_rows(v.segm, sep = ",") %>%
  mutate(mhc.a.split = str_split_fixed(mhc.a, fixed(":"), n = Inf)[,1],
         mhc.b.split = str_split_fixed(mhc.b, fixed(":"), n = Inf)[,1],
         v.segm.split = str_split_fixed(v.segm, fixed("*"), n = Inf)[,1]) %>%
  group_by(gene, mhc.class, mhc.a.split, mhc.b.split, v.segm.split) %>%
  summarize(records = length(unique(id))) %>%
  group_by(mhc.class, mhc.a.split, mhc.b.split) %>%
  mutate(records.mhc = sum(records)) %>%
  group_by(gene, v.segm.split) %>%
  mutate(records.v = sum(records)) %>%
  filter(records.mhc >= 10, records.v >= 10)
  
#df.vhla$v.segm.split = with(df.vhla, factor(v.segm.split, v.segm.split[order(records.v)]))
#df.vhla$mhc.a.split = with(df.vhla, factor(mhc.a.split, mhc.a.split[order(records.mhc)]))

ggplot(df.vhla, aes(x=gsub("HLA-", "", paste(mhc.a.split, mhc.b.split, sep = " / ")) %>%
                      fct_reorder(records), 
                    y=v.segm.split %>%
                      fct_reorder(records), fill = pmin(records, 1000))) +
  geom_tile() +
  scale_fill_gradientn("Records", colors=colorRampPalette(brewer.pal(9, 'PuBuGn'))(32), 
                       trans="log", breaks = c(1, 10, 100, 1000)) +
  xlab("") + ylab("") +
  facet_grid(gene~mhc.class, scales="free", space="free") +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        axis.text.y = element_text(size = 8),
        panel.grid.major = element_blank(),
        legend.position = "right",
        axis.line = element_line(linewidth = 0.3),
        strip.background = element_blank())
```

<img src="vdjdb_summary_files/figure-html/unnamed-chunk-13-1.png" alt="" width="576" />

Circos plot for correspondence between human TRBV genes and MHC class I alleles for links supported by more than 50 records in terms of the ratio of observed to expected records. Band width and color (from red:highest to light yellow:lowest) are scaled proportional to co-occurrence matrix value divided by row (TRBV) and column (HLA) sums.


``` r
#https://jokergoo.github.io/circlize_book
df.vhla.1b <- df.vhla %>%
  filter(records >= 50) %>%
  filter(gene == "TRB", mhc.class == "MHCI") %>%
  mutate(mhc.a.split = substr(mhc.a.split, 5, nchar(mhc.a.split)),
         v.segm.split = paste0("Vb", substr(v.segm.split, 5, nchar(v.segm.split)))) %>%
  reshape2::dcast(mhc.a.split ~ v.segm.split, value.var = "records", fill = 0)
rownames(df.vhla.1b) <- df.vhla.1b$mhc.a.split
df.vhla.1b$mhc.a.split <- NULL
df.vhla.1b <- df.vhla.1b %>% as.matrix()
df.vhla.1b <- (df.vhla.1b %>% 
  sweep(1, rowSums(df.vhla.1b), `/`) %>%
  sweep(2, colSums(df.vhla.1b), `/`)) * sum(df.vhla.1b)

grid.col <- setNames(c(colorRampPalette(brewer.pal(9, 'YlGn'))(nrow(df.vhla.1b)) %>% rev, 
                       colorRampPalette(brewer.pal(9, 'PuBu'))(ncol(df.vhla.1b))), 
                     union(rownames(df.vhla.1b), colnames(df.vhla.1b)))
band.col <- colorRamp2(range(df.vhla.1b), 
                       colorRampPalette(brewer.pal(9, 'YlOrRd'))(2), 
                       transparency = 0.2)

chordDiagram(df.vhla.1b, annotationTrack = "grid", preAllocateTracks = 1, 
             col = band.col,
             grid.col = grid.col,
             grid.border = "black",
             big.gap = 20)

circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, 
              facing = "clockwise", niceFacing = T, adj = c(0, 0.5), 
              cex = 0.65)
}, bg.border = NA)

title("TRBV ~ HLA-I")
```

<img src="vdjdb_summary_files/figure-html/unnamed-chunk-14-1.png" alt="" width="768" />

``` r
circos.clear()
```

---

#### Detailed summary for HLA

Representative data for Homo Sapiens MHC class I and II


``` r
df.m = df %>% 
  as_tibble %>%
  filter(species == "HomoSapiens") %>%
  mutate(id = 1:n()) %>%
  separate_rows(mhc.a, sep = ",") %>%
  separate_rows(mhc.b, sep = ",") %>%
  mutate(mhc.a.split = str_split_fixed(mhc.a, fixed(":"), n = Inf)[,1],
         mhc.b.split = str_split_fixed(mhc.b, fixed(":"), n = Inf)[,1]) %>%
  group_by(mhc.a.split, mhc.b.split) %>%
  summarize(records = length(unique(id)), 
            antigens = length(unique(antigen.epitope)),
            publications = length(unique(str_split_fixed(reference.id, ",", n = Inf)[,1]))) %>%
  arrange(-records)
  
colnames(df.m) = c("First chain", "Second chain", "Records", "Unique epitopes", "Studies")
kable(format = "html", df.m)
```

<table>
 <thead>
  <tr>
   <th style="text-align:left;"> First chain </th>
   <th style="text-align:left;"> Second chain </th>
   <th style="text-align:right;"> Records </th>
   <th style="text-align:right;"> Unique epitopes </th>
   <th style="text-align:right;"> Studies </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> HLA-A*02 </td>
   <td style="text-align:left;"> B2M </td>
   <td style="text-align:right;"> 113756 </td>
   <td style="text-align:right;"> 1039 </td>
   <td style="text-align:right;"> 210 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HLA-A*03 </td>
   <td style="text-align:left;"> B2M </td>
   <td style="text-align:right;"> 27825 </td>
   <td style="text-align:right;"> 23 </td>
   <td style="text-align:right;"> 19 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HLA-A*11 </td>
   <td style="text-align:left;"> B2M </td>
   <td style="text-align:right;"> 7659 </td>
   <td style="text-align:right;"> 22 </td>
   <td style="text-align:right;"> 19 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HLA-B*07 </td>
   <td style="text-align:left;"> B2M </td>
   <td style="text-align:right;"> 7383 </td>
   <td style="text-align:right;"> 125 </td>
   <td style="text-align:right;"> 35 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HLA-B*08 </td>
   <td style="text-align:left;"> B2M </td>
   <td style="text-align:right;"> 4366 </td>
   <td style="text-align:right;"> 74 </td>
   <td style="text-align:right;"> 27 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HLA-A*01 </td>
   <td style="text-align:left;"> B2M </td>
   <td style="text-align:right;"> 3410 </td>
   <td style="text-align:right;"> 161 </td>
   <td style="text-align:right;"> 28 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HLA-A*24 </td>
   <td style="text-align:left;"> B2M </td>
   <td style="text-align:right;"> 3159 </td>
   <td style="text-align:right;"> 97 </td>
   <td style="text-align:right;"> 37 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HLA-E*01 </td>
   <td style="text-align:left;"> B2M </td>
   <td style="text-align:right;"> 2617 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 8 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HLA-DPA1*01 </td>
   <td style="text-align:left;"> HLA-DPB1*04 </td>
   <td style="text-align:right;"> 1509 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 12 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HLA-DRA*01 </td>
   <td style="text-align:left;"> HLA-DRB1*01 </td>
   <td style="text-align:right;"> 1294 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 15 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HLA-B*27 </td>
   <td style="text-align:left;"> B2M </td>
   <td style="text-align:right;"> 868 </td>
   <td style="text-align:right;"> 74 </td>
   <td style="text-align:right;"> 17 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HLA-DQA1*05 </td>
   <td style="text-align:left;"> HLA-DQB1*06 </td>
   <td style="text-align:right;"> 768 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HLA-B*57 </td>
   <td style="text-align:left;"> B2M </td>
   <td style="text-align:right;"> 752 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 12 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HLA-B*15 </td>
   <td style="text-align:left;"> B2M </td>
   <td style="text-align:right;"> 644 </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 8 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HLA-B*35 </td>
   <td style="text-align:left;"> B2M </td>
   <td style="text-align:right;"> 596 </td>
   <td style="text-align:right;"> 49 </td>
   <td style="text-align:right;"> 28 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HLA-DRA*01 </td>
   <td style="text-align:left;"> HLA-DRB1*15 </td>
   <td style="text-align:right;"> 592 </td>
   <td style="text-align:right;"> 27 </td>
   <td style="text-align:right;"> 15 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HLA-DRA*01 </td>
   <td style="text-align:left;"> HLA-DRB1*04 </td>
   <td style="text-align:right;"> 521 </td>
   <td style="text-align:right;"> 22 </td>
   <td style="text-align:right;"> 16 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HLA-DRA*01 </td>
   <td style="text-align:left;"> HLA-DRB5*01 </td>
   <td style="text-align:right;"> 392 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 7 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HLA-DQA1*01 </td>
   <td style="text-align:left;"> HLA-DQB1*06 </td>
   <td style="text-align:right;"> 358 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 10 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HLA-B*42 </td>
   <td style="text-align:left;"> B2M </td>
   <td style="text-align:right;"> 354 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 6 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HLA-DQA1*05 </td>
   <td style="text-align:left;"> HLA-DQB1*02 </td>
   <td style="text-align:right;"> 332 </td>
   <td style="text-align:right;"> 26 </td>
   <td style="text-align:right;"> 16 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HLA-DRA*01 </td>
   <td style="text-align:left;"> HLA-DRB1*07 </td>
   <td style="text-align:right;"> 301 </td>
   <td style="text-align:right;"> 11 </td>
   <td style="text-align:right;"> 8 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HLA-B*44 </td>
   <td style="text-align:left;"> B2M </td>
   <td style="text-align:right;"> 293 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 7 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HLA-B*37 </td>
   <td style="text-align:left;"> B2M </td>
   <td style="text-align:right;"> 248 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HLA-DRA*01 </td>
   <td style="text-align:left;"> HLA-DRB1*11 </td>
   <td style="text-align:right;"> 224 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 8 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HLA-C*12 </td>
   <td style="text-align:left;"> B2M </td>
   <td style="text-align:right;"> 215 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HLA-A*68 </td>
   <td style="text-align:left;"> B2M </td>
   <td style="text-align:right;"> 212 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 4 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HLA-B*40 </td>
   <td style="text-align:left;"> B2M </td>
   <td style="text-align:right;"> 120 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HLA-A*08 </td>
   <td style="text-align:left;"> B2M </td>
   <td style="text-align:right;"> 117 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HLA-B*18 </td>
   <td style="text-align:left;"> B2M </td>
   <td style="text-align:right;"> 84 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 4 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HLA-A*80 </td>
   <td style="text-align:left;"> B2M </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HLA-DRA*01 </td>
   <td style="text-align:left;"> HLA-DRB1*03 </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 7 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HLA-DQA1*03 </td>
   <td style="text-align:left;"> HLA-DQB1*03 </td>
   <td style="text-align:right;"> 65 </td>
   <td style="text-align:right;"> 20 </td>
   <td style="text-align:right;"> 16 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HLA-DQA1*03 </td>
   <td style="text-align:left;"> HLA-DQB1*02 </td>
   <td style="text-align:right;"> 62 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HLA-DRA*01 </td>
   <td style="text-align:left;"> HLA-DRB3*02 </td>
   <td style="text-align:right;"> 62 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 4 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HLA-B*81 </td>
   <td style="text-align:left;"> B2M </td>
   <td style="text-align:right;"> 48 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HLA-C*07 </td>
   <td style="text-align:left;"> B2M </td>
   <td style="text-align:right;"> 43 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HLA-DQA1*02 </td>
   <td style="text-align:left;"> HLA-DQB1*02 </td>
   <td style="text-align:right;"> 32 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HLA-B*51 </td>
   <td style="text-align:left;"> B2M </td>
   <td style="text-align:right;"> 31 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HLA-DQA1*01 </td>
   <td style="text-align:left;"> HLA-DRB1*07 </td>
   <td style="text-align:right;"> 29 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HLA-DQA1*01 </td>
   <td style="text-align:left;"> HLA-DQB1*02 </td>
   <td style="text-align:right;"> 20 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HLA-DQA1*02 </td>
   <td style="text-align:left;"> HLA-DQB1*06 </td>
   <td style="text-align:right;"> 20 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HLA-DQA1*05 </td>
   <td style="text-align:left;"> HLA-DQB1*03 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 4 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HLA-DPA1*02 </td>
   <td style="text-align:left;"> HLA-DPB1*13 </td>
   <td style="text-align:right;"> 16 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HLA-DRA*01 </td>
   <td style="text-align:left;"> HLA-DRB3*03 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HLA-DQA1*01 </td>
   <td style="text-align:left;"> HLA-DQB1*05 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 5 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HLA-DRA*01 </td>
   <td style="text-align:left;"> HLA-DRB3*01 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HLA-B*53 </td>
   <td style="text-align:left;"> B2M </td>
   <td style="text-align:right;"> 11 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HLA-C*08 </td>
   <td style="text-align:left;"> B2M </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HLA-A*30 </td>
   <td style="text-align:left;"> B2M </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HLA-B*38 </td>
   <td style="text-align:left;"> B2M </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HLA-A*29 </td>
   <td style="text-align:left;"> B2M </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HLA-B*41 </td>
   <td style="text-align:left;"> B2M </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HLA-DPA1*01 </td>
   <td style="text-align:left;"> HLA-DPB1*02 </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 4 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HLA-DPA1*02 </td>
   <td style="text-align:left;"> HLA-DPB1*05 </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HLA-C*05 </td>
   <td style="text-align:left;"> B2M </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HLA-DPA1*02 </td>
   <td style="text-align:left;"> HLA-DPB1*01 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HLA-DPA1*02 </td>
   <td style="text-align:left;"> HLA-DPB1*14 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HLA-DRA*01 </td>
   <td style="text-align:left;"> HLA-DRB1*09 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HLA-C*01 </td>
   <td style="text-align:left;"> B2M </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HLA-C*03 </td>
   <td style="text-align:left;"> B2M </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HLA-A*25 </td>
   <td style="text-align:left;"> B2M </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HLA-B*12 </td>
   <td style="text-align:left;"> B2M </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HLA-B*58 </td>
   <td style="text-align:left;"> B2M </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HLA-C*06 </td>
   <td style="text-align:left;"> B2M </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HLA-C*14 </td>
   <td style="text-align:left;"> B2M </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HLA-DPA1*01 </td>
   <td style="text-align:left;"> HLA-DPB1*03 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HLA-DRA*01 </td>
   <td style="text-align:left;"> HLA-DRB1*08 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HLA-DRA*01 </td>
   <td style="text-align:left;"> HLA-DRB1*13 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HLA-A*32 </td>
   <td style="text-align:left;"> B2M </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HLA-B*14 </td>
   <td style="text-align:left;"> B2M </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HLA-B*52 </td>
   <td style="text-align:left;"> B2M </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HLA-C*04 </td>
   <td style="text-align:left;"> B2M </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HLA-DPA1*01 </td>
   <td style="text-align:left;"> HLA-DRB1*04 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HLA-DRA*01 </td>
   <td style="text-align:left;"> HLA-DPB1*04 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HLA-DRA*01 </td>
   <td style="text-align:left;"> HLA-DRB4*01 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
</tbody>
</table>

!summary_embed_end!
