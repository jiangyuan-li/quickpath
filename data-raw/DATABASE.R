## code to prepare `DATABASE` dataset goes here
################### mart #######################
mmu_mart <- readRDS("data-raw/mart_mmu.rds")
usethis::use_data(mmu_mart, overwrite = TRUE)

gga_mart <- readRDS("data-raw/mart_gga.rds")
usethis::use_data(gga_mart, overwrite = TRUE)

hsa_mart <- readRDS("data-raw/mart_hsa.rds")
usethis::use_data(hsa_mart, overwrite = TRUE)

################### mmu #######################
mmu_pathway <- utils::read.csv("data-raw/mmu_pathway.csv", stringsAsFactors = FALSE)
usethis::use_data(mmu_pathway, overwrite = TRUE)

mmu_path_gene_list <- readRDS("data-raw/mmu_path_gene_list.rds")
usethis::use_data(mmu_path_gene_list, overwrite = TRUE)

mmu_genes_pathways <- utils::read.csv("data-raw/genes_pathways_mmu.csv", stringsAsFactors = FALSE)
usethis::use_data(mmu_genes_pathways, overwrite = TRUE)

################### gga #######################
gga_pathway <- utils::read.csv("data-raw/gga_pathway.csv")
usethis::use_data(gga_pathway, overwrite = TRUE)

gga_path_gene_list <- readRDS("data-raw/gga_path_gene_list.rds")
usethis::use_data(gga_path_gene_list, overwrite = TRUE)

gga_genes_pathways <- utils::read.csv("data-raw/genes_pathways_gga.csv")
usethis::use_data(gga_genes_pathways, overwrite = TRUE)

################### hsa #######################
hsa_pathway <- utils::read.csv("data-raw/hsa_pathway.csv")
usethis::use_data(hsa_pathway, overwrite = TRUE)

hsa_path_gene_list <- readRDS("data-raw/hsa_path_gene_list.rds")
usethis::use_data(hsa_path_gene_list, overwrite = TRUE)

hsa_genes_pathways <- utils::read.csv("data-raw/genes_pathways_hsa.csv")
usethis::use_data(hsa_genes_pathways, overwrite = TRUE)

usethis::use_data("DATABASE")
