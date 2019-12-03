#' Grab differential expressed genes (DEGs) from cuffdiff result
#'
#'@param gene_exp.diff a dataframe from cuffdiff result
#'@param out.name If specified, the output will be written to this file.
#'@param criterion p_value or q_value
#'@param cut.off significance level
#'@return A dataframe contains all DEGs with gene symvbols.
#'@examples
#'head(gene_exp.diff)
#'deg = grab_deg_from_cuffdiff(gene_exp.diff)
#'head(deg)
#'@export
grab_degs_from_cuffdiff <- function(gene_exp.diff, out.name = NULL, class = c("mmu","gga","hsa"), criterion = c("p_value","q_value"), cut.off = 0.05){

  # get DEGs based on criterion and cutoff
  deg <- gene_exp.diff[gene_exp.diff[,criterion] <= cut.off, ]

  # get gene symbol
  mart <- eval(parse(text = paste0(class,"_mart")))
  symbol = biomaRt::getBM(attributes=c("ensembl_gene_id", "external_gene_name","description"),
                 filters="ensembl_gene_id", values=deg$gene_id, mart=mart)
  names(symbol)[1] <- "gene_id"

  # combine symbol info and expression info
  out <- merge(deg,symbol,by="gene_id",all.x = T)

  if(!is.null(out.name)){
    # writting DEGs into a .csv file
    write.csv(out, out.name, row.names = FALSE)
  }

  # return a Dataframe with DEGs' info
  return(out)
}
