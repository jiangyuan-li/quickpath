#' title
#'
#'@param
#'@param
#'@return
#'@examples
#'eg
#'@export

extract_genes_by_pathway <- function(pathway, deg.file, out.name, class=c("mmu","hsa","gga")){

  ref.paths = eval(parse(text = paste0(class,"_pathway")))
  row.names(ref.paths) <- NULL
  path.ids <- c()
  for(i in 1:length(pathway)){
    print(ref.paths[grep(pathway[i],ref.paths[,2],ignore.case = TRUE),])
    print(i)
    path.ids[i] <- ref.paths[grep(pathway[i],ref.paths[,2],ignore.case = TRUE),1]
  }

  #all <- eval(parse(text = paste0("genes_pathways_",class)))

  path_list <- eval(parse(text = paste0(class,"_path_gene_list")))
  wb = xlsx::createWorkbook()

  for (i in 1:length(pathway)){
    sheet = createSheet(wb, pathway[i])

    symbol = as.character(path_list[[path.ids[i]]][,'gene.symbol'])

    out.deg <- deg[deg$gene%in%symbol,]
    print(nrow(out.deg))
    print(length(symbol))


    print(path.ids[i])
    addDataFrame(out.deg, sheet=sheet, startRow =1, row.names=FALSE)


  }

  xslx::saveWorkbook(wb, out.name)

}
