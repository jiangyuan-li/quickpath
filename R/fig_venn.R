#' title
#'
#'@param
#'@return
#'@examples
#'eg
#'@export

fig_venn <- function(list.deg){

  library(VennDiagram)
  mcd <- na.omit(read.csv("../../MCD/MCD_DEG.csv"))
  chow <- na.omit(read.csv("../../Chow/Chow_DEG.csv"))

  num.mcd <- length(unique(mcd$external_gene_name))
  num.chow <- length(unique(chow$external_gene_name))

  over.mcdchow <- length(intersect(mcd$external_gene_name,chow$external_gene_name))

  #grid.newpage()
  tiff("Venn_MCD_Chow.tiff",width = 3600,height = 2400,res=600)
  draw.pairwise.venn(num.mcd,num.chow,over.mcdchow,
                     category = c("Chow","MCD"),
                     lty = rep("blank",2), fill = c("light blue", "pink"), alpha = rep(0.5, 2),
                     cat.pos = c(0,0), cat.dist = rep(0.025, 2),scaled = F)

  dev.off()
}
