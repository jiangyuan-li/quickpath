## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- echo = T----------------------------------------------------------------
library(quickpath)

## ---- message = F, echo = T---------------------------------------------------
data(gene_exp.diff)
print(dim(gene_exp.diff))
head(gene_exp.diff, n = 5)
deg = grab_deg_from_cuffdiff(gene_exp.diff, class = "mmu", criterion = "p_value", cut.off = 0.05)
print(dim(deg))
head(deg[,c(1,4,8,9,10,11,12,13,15,16)],n = 5)

## ---- echo = T----------------------------------------------------------------
sig.genes = deg$external_gene_name
path_res = pathway_analysis(sig.genes, class = "mmu")
head(path_res[[1]])

## ---- echo = T----------------------------------------------------------------
deg_egg = grab_deg_from_cuffdiff(gene_exp.diff_egg, class = "gga")
sig.genes = deg_egg$external_gene_name
path_res_egg = pathway_analysis(sig.genes, class = "gga")

## ---- echo = T----------------------------------------------------------------
pathway <- c("Purine metabolism", "PI3K-Akt signaling pathway", 
             "AMPK signaling pathway", "Choline metabolism in cancer")
deg.list = extract_genes_by_pathway(pathway, deg, class = "mmu")

## ---- echo = T----------------------------------------------------------------
deg_chow = grab_deg_from_cuffdiff(gene_exp.diff_chow, class = "mmu")
sig.genes = deg_chow$external_gene_name
path_res_chow = pathway_analysis(sig.genes, class = "mmu")

deg_e105 = grab_deg_from_cuffdiff(gene_exp.diff_e105, class = "mmu")
sig.genes = deg_e105$external_gene_name
path_res_e105 = pathway_analysis(sig.genes, class = "mmu")

deg_e95 = grab_deg_from_cuffdiff(gene_exp.diff_e95, class = "mmu")
sig.genes = deg_e95$external_gene_name
path_res_e95 = pathway_analysis(sig.genes, class = "mmu")

## ---- fig.height=8, fig.width=8, echo=T---------------------------------------
list.info = list(path_res[[2]], path_res_chow[[2]], path_res_e95[[2]], path_res_e105[[2]])
path.ids = check_pathway_name(pathway, class = "mmu")
group.info = c("MCD","Chow","E9.5","E10.5")
path.names = c("1st path","2nd path","3rd path","4th path")
fig_path(path.ids, list.info, group.info, criterion = "percentage", path.names = path.names)
fig_path(path.ids, list.info, group.info, criterion = "pval", path.names = path.names)

## ---- message = F, echo = T---------------------------------------------------
total.genes = mmu_genes_pathways$gene.symbol
deg = grab_deg_from_cuffdiff(gene_exp.diff, class = "mmu", criterion = "p_value", cut.off = 0.05)
sig.genes = deg$external_gene_name
meth_path_res = pathway_analysis_meth(sig.genes,total.genes, class = "mmu")

