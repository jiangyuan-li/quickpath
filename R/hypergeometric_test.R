#' hypergeometric test for RNA-seq
#'
#'@param sig.genes A vector of significantly changed genes
#'@param class Specify which type of gene
#'@param mc.cores Specify how many cores to use
#'@return A list of all significantly changed pathways, and all pathways with p value
#'@examples
#'deg = grab_deg_from_cuffdiff(gene_exp.diff)
#'sig.genes = deg$external_gene_name
#'path_res = pathway_analysis(sig.genes, class = "mmu")
#'@export
pathway_analysis <- function(sig.genes, class=c("mmu","hsa","gga"), mc.cores = 8){

  ############## prepare needed data ####################
  all = eval(parse(text = paste0(class,"_genes_pathways")))
  paths = eval(parse(text = paste0(class,"_pathway")))
  path.gene = eval(parse(text = paste0(class,"_path_gene_list")))

  ############# clean background gene list ##############
  all.genes <- all$gene.symbol
  all.genes <- all.genes[!is.na(all.genes)]
  all.genes <- unique(all.genes)

  ############ clean DEGs list ##########################
  sig.genes <- sig.genes[!is.na(sig.genes)]
  sig.genes <- unique(sig.genes)
  sig.genes <- unique(intersect(sig.genes,all.genes))
  #######################################################

  print(paste0("sig genes in pathway: ",length(sig.genes),"/",length(all.genes)))

  ################### do hypergeometric test ############
  out <- parallel::mclapply(path.gene, each_test,all.genes = all.genes,sig.genes = sig.genes,mc.cores = mc.cores)
  print("Inference part is Done.")

  ############## format results #########################
  tmpname <- names(out)
  out <- t(as.data.frame(out))
  colnames(out) <- c("pval","percentage")
  out <- as.data.frame(out)
  rownames(out) <- tmpname
  qval <- p.adjust(out$pval,"fdr")
  df <- data.frame(name = rownames(out),pval=as.numeric(out$pval),qval=as.numeric(qval),percentage=as.numeric(out$percentage))
  out <- merge(paths,df,by="name")
  out <- na.omit(out)

  return(list(out[out$pval<=0.05,],
              out))

}

######################################################################
######################################################################
######################################################################

#' hypergeometric test for DNA-methylation seq
#'
#'@param sig.genes A vector of significantly changed genes
#'@param total.genes A vector of total genes
#'@param class Specify which type of gene
#'@param mc.cores Specify how many cores to use
#'@return A list of all significantly changed pathways, and all pathways with p value
#'@examples
#'path_res = pathway_analysis_meth(sig.genes, class = "mmu")
#'@export
pathway_analysis_meth <- function(sig.genes, total.genes, class=c("mmu","hsa","gga"), mc.cores = 8){

  ############## prepare needed data ###################
  all = eval(parse(text = paste0(class,"_genes_pathways")))
  paths = eval(parse(text = paste0(class,"_pathway")))
  path.gene = eval(parse(text = paste0(class,"_path_gene_list")))

  ############# clean background gene list #############
  all.genes <- all$gene.symbol
  all.genes <- all.genes[!is.na(all.genes)]
  all.genes <- unique(all.genes)

  all.genes <- intersect(all.genes, total.genes)
  ############ clean DEGs list #########################
  sig.genes <- sig.genes[!is.na(sig.genes)]
  sig.genes <- unique(sig.genes)
  sig.genes <- unique(intersect(sig.genes,all.genes))
  ######################################################

  print(paste0("sig genes in pathway: ",length(sig.genes),"/",length(all.genes)))

  ################### do hypergeometric test ###########
  out <- parallel::mclapply(path.gene, each_test,all.genes = all.genes,sig.genes = sig.genes,mc.cores = mc.cores)
  print("Inference part is Done.")

  ############## format results ########################
  tmpname <- names(out)
  out <- t(as.data.frame(out))
  colnames(out) <- c("pval","percentage")
  out <- as.data.frame(out)
  rownames(out) <- tmpname
  qval <- p.adjust(out$pval,"fdr")
  df <- data.frame(name = rownames(out),pval=as.numeric(out$pval),qval=as.numeric(qval),percentage=as.numeric(out$percentage))
  out <- merge(paths,df,by="name")
  out <- na.omit(out)

  return(list(out[out$pval<=0.05,],
              out))

}

########################################################
################ test for each pathway #################
########################################################
# This function is NOT exported
each_test <- function(x, all.genes, sig.genes){

  ############## check whether x is valid ##############
  if(!is.na(x)){

    ############## clean genes in pathway ##############
    symbol <- na.omit(x[,1])
    symbol = symbol[!is.na(symbol)]
    symbol = intersect(symbol, all.genes)

    ############## perform hypergeometric test #########
    pval.percent = pvalhyper(all.genes,sig.genes,symbol)
    return(pval.percent)
  }
  else{
    return(c(NA,NA))
  }
}

######################################################################
############## get pvalue from hypergeometric distribution ###########
######################################################################
# This function is NOT exported
pvalhyper <- function(all.genes, sig.genes, gs){

  ############## get length of target genes ##############
  N = length(all.genes)
  R = length(sig.genes)

  ############## get length of genes in pathway ##########
  n = length(gs)
  r = length(intersect(gs, sig.genes))

  ############## check whether r is 0 ####################
  if(r==0){
    return(c(1,r/n))
  }
  else{

    ############## get p value ###########################
    pval <- 1 - sum(dhyper(0:r-1,R,N-R,n))
    return(c(pval,r/n))
  }
}
