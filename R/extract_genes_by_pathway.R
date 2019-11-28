#' Extract DEGs from pathways
#'
#'@param pathway A vector contains interested pathwat names
#'@param deg A dataframe contains DEGS information
#'@param class Specify which type of gene
#'@param out.name A file name ends with ".xlsx", if specified, the output will be written into this file
#'@returnn A list of DEGS for each pathway
#'@examples
#'pathway <- c("Purine metabolism", "PI3K-Akt signaling pathway", "AMPK signaling pathway", "Choline metabolism in cancer")
#'deg = grab_deg_from_cuffdiff(gene_exp.diff)
#'deg.list = extract_degs_by_pathway(pathway, deg, class = "mmu")
#'@export

extract_degs_by_pathway <- function(pathway, deg, class = c("mmu","hsa","gga"), out.name = NULL){

  n = length(pathway)
  path.ids <- check_pathway_name(pathway, class = class)

  path_list <- eval(parse(text = paste0(class,"_path_gene_list")))

  out.list <- list()
  for(i in 1:n){
    symbol = as.character(path_list[[path.ids[i]]][,'gene.symbol'])
    symbol = symbol[!is.na(symbol)]

    out.deg <- deg[deg$external_gene_name%in%symbol,]
    out.list[[i]] <- out.deg

    print(paste0(pathway[i]," has DEGs ",nrow(out.deg),"/",length(symbol)))
  }
  if(!is.null(out.name)){
    message("Writing the output to a xlsx file...")
    wb = xlsx::createWorkbook()

    for (i in 1:length(pathway)){
      sheet = createSheet(wb, pathway[i])
      print(path.ids[i])
      addDataFrame(out.deg, sheet=sheet, startRow =1, row.names=FALSE)
    }
    xslx::saveWorkbook(wb, out.name)
    message("End writing...")
  }
  return(out.list)
}


####################check pathway name and return path ids ##############
#' Extract DEGs from pathways
#'
#'@param pathway A vector contains interested pathwat names
#'@param class Specify which type of gene
#'@return
#'@examples Path IDs for each pathway
#'pathway <- c("Purine metabolism", "PI3K-Akt signaling pathway", "AMPK signaling pathway", "Choline metabolism in cancer")
#'path.ids = check_pathway_name(pathway, class = "mmu")
#'@export
check_pathway_name <- function(pathway, class = c("mmu","hsa","gga")){

  all.paths = eval(parse(text = paste0(class,"_pathway")))
  ref.paths = all.paths[,2]
  names(ref.paths) = all.paths[,1]

  path.ids <- c()
  n = length(pathway)

  message("Start checking pathway names...")
  for(i in 1:n){
    print(paste0(i,"/",n))
    ind = grep(pathway[i],ref.paths,ignore.case = TRUE)
    tmp.info = ref.paths[ind]
    print(tmp.info)
    if(length(ind) != 1){
      stop("Input contains ambiguous or invalid pathway name!")
    }
    path.ids[i] <- names(tmp.info)
  }
  message("All names are valid!")

  return(path.ids)
}
