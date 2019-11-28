## code to prepare `sample_data` dataset goes here
#usethis::use_data("sample_data")

gene_exp.diff <- utils::read.csv("data-raw/gene_exp_MCD.diff", sep="\t",stringsAsFactors = FALSE)
usethis::use_data(gene_exp.diff, overwrite = TRUE)

gene_exp.diff_chow <- utils::read.csv("data-raw/gene_exp_Chow.diff", sep="\t",stringsAsFactors = FALSE)
usethis::use_data(gene_exp.diff_chow, overwrite = TRUE)

gene_exp.diff_egg <- utils::read.csv("data-raw/gene_exp_egg.diff", sep="\t",stringsAsFactors = FALSE)
usethis::use_data(gene_exp.diff_egg, overwrite = TRUE)

gene_exp.diff_e95 <- utils::read.csv("data-raw/E9.5_gene_exp.diff", sep="\t",stringsAsFactors = FALSE)
usethis::use_data(gene_exp.diff_e95, overwrite = TRUE)

gene_exp.diff_e105 <- utils::read.csv("data-raw/E10.5_gene_exp.diff", sep="\t",stringsAsFactors = FALSE)
usethis::use_data(gene_exp.diff_e105, overwrite = TRUE)

