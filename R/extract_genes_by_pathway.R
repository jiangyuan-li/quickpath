#' Extract DEGs from pathways
#'
#'@param pathway A vector contains interested pathwat names
#'@param deg A dataframe contains DEGS information
#'@param class Specify which type of gene
#'@param out.name A file name ends with ".xlsx", if specified, the output will be written into this file
#'@return A list of DEGS for each pathway
#'@examples
#'pathway <- c("Purine metabolism", "PI3K-Akt signaling pathway", "AMPK signaling pathway", "Choline metabolism in cancer")
#'deg = grab_deg_from_cuffdiff(gene_exp.diff)
#'deg.list = extract_degs_by_pathway(pathway, deg, class = "mmu")
#'@export

extract_genes_by_pathway <- function(pathway, deg, class = c("mmu","hsa","gga"), out.name = NULL){

  # get number of pathways
  n = length(pathway)

  # check whether pathway names are valid
  path.ids <- check_pathway_name(pathway, class = class)

  # get pathway list
  path_list <- eval(parse(text = paste0(class,"_path_gene_list")))

  # grab DEGs for each pathway
  out.list <- list()

  for(i in 1:n){
    # extract reference genes from pathway list
    symbol = as.character(path_list[[path.ids[i]]][,'gene.symbol'])
    symbol = symbol[!is.na(symbol)]

    # get info of DEGs
    out.deg <- deg[deg$external_gene_name%in%symbol,]
    out.list[[i]] <- out.deg

    print(paste0(pathway[i]," has DEGs ",nrow(out.deg),"/",length(symbol)))
  }

  # writing output into xlsx
  if(!is.null(out.name)){

    message("Writing the output to a xlsx file...")
    wb = xlsx::createWorkbook()

    for (i in 1:n){

      # writing each pathway into one sheet
      sheet = xlsx::createSheet(wb, pathway[i])
      print(path.ids[i])
      xlsx::addDataFrame(out.list[[i]], sheet=sheet, startRow =1, row.names=FALSE)
    }

    # save Workbook
    xlsx::saveWorkbook(wb, out.name)
    message("End writing...")
  }

  # return the output list
  return(out.list)
}

#########################################################################
####################check pathway name and return path ids ##############
#########################################################################

#' Extract DEGs from pathways
#'
#'@param pathway A vector contains interested pathwat names
#'@param class Specify which type of gene
#'@return A vector of path IDs
#'@examples Path IDs for each pathway
#'pathway <- c("Purine metabolism", "PI3K-Akt signaling pathway", "AMPK signaling pathway", "Choline metabolism in cancer")
#'path.ids = check_pathway_name(pathway, class = "mmu")
#'@export
check_pathway_name <- function(pathway, class = c("mmu","hsa","gga")){

  # get all pathways as reference
  all.paths = eval(parse(text = paste0(class,"_pathway")))
  ref.paths = all.paths[,2]
  names(ref.paths) = all.paths[,1]

  # create path ids holder
  path.ids <- c()
  n = length(pathway)

  message("Start checking pathway names...")

  for(i in 1:n){

    print(paste0(i,"/",n))
    # get index after matching pathway with references
    ind = grep(pathway[i],ref.paths,ignore.case = TRUE)
    # get info and print it
    tmp.info = ref.paths[ind]
    print(tmp.info)

    # check whether pathway name is valid
    if(length(ind) != 1){
      stop("Input contains ambiguous or invalid pathway name!")
    }

    # get path id
    path.ids[i] <- names(tmp.info)
  }

  message("All names are valid!")

  # return path ids
  return(path.ids)
}
