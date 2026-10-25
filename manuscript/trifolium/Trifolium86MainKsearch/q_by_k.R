# Write the best-run Q matrix for K 3-8 (Trifolium86MainKsearch) to ksearch/Q_K{k}.tsv, rows in .fam order.
suppressPackageStartupMessages({ library(LEA); library(data.table) })
W <- "Trifolium86MainKsearch_results/_work/maf0.1_miss0.3_smiss0.8_rel0.999/ld0.2_win100_step10"
proj <- load.snmfProject(file.path(W, "Trifolium86_chr1-6.snmfProject"))
s <- sub("^0_", "", fread(file.path(W, "Trifolium86_chr1-6.fam"), header = FALSE, colClasses = "character")$V2)
for (k in 3:8) {
  q <- Q(proj, K = k, run = which.min(cross.entropy(proj, K = k)))
  fwrite(data.table(sample = s, q), sprintf("ANALYSIS/Trifolium86MainKsearch/ksearch/Q_K%d.tsv", k), sep = "\t")
}
