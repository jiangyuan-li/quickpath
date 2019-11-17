#' hypergeometric test
#'
#'@param sig.genes A vector of significantly changed genes
#'@param class Specify which type of gene
#'@return A list of all significantly changed pathways, and all pathways with p value
#'@examples
#'pathway_analysis(c("Osr1","Tbx5"),"mmu")
#'@export
pathway_analysis <- function(sig.genes,class=c("mmu","hsa","gga")){
  library(readr)
  pvalhyper <- function(all.genes, sig.genes, gs){
    N = length(all.genes)
    n = length(gs)
    R = length(sig.genes)
    r = length(intersect(gs, sig.genes))
    pval <- sum(dhyper(0:r,R,N-R,n))

    #if(r==0){
    #  return(c(pval,r/n))
    #}
    #else{
    #  return(c(min(c(pval,1-pval)),r/n))
    #}

    if(r==0){
      return(c(NA,r/n))
    }
    else{
      return(c(1-pval,r/n))
    }
  }

  if(class == "mmu"){
    all <- read.csv("~/Biodatabase/Pathway/genes_pathways_mmu.csv")
    paths <- read.csv("~/Biodatabase/Pathway/mmu_pathway.csv")
    path.gene <- readRDS("~/Biodatabase/Pathway/mmu_path_gene_list.rds")
    print(class)

  }
  if(class == "hsa"){
    all <- read.csv("~/Biodatabase/Pathway/genes_pathways_hsa.csv")
    paths <- read.csv("~/Biodatabase/Pathway/hsa_pathway.csv")
    path.gene <- readRDS("~/Biodatabase/Pathway/hsa_path_gene_list.rds")
    print(class)

  }
  if(class == "gga"){
    all <- read.csv("~/Biodatabase/Pathway/genes_pathways_gga.csv")
    paths <- read.csv("~/Biodatabase/Pathway/gga_pathway.csv")
    path.gene <- readRDS("~/Biodatabase/Pathway/gga_path_gene_list.rds")
    print(class)

  }
  all.genes <- unique(all$gene.symbol)
  all.genes <- all.genes[!is.na(all.genes)]
  all.genes <- all.genes[!all.genes==""]

  sig.genes <- unique(intersect(sig.genes,all.genes))
  print(paste("sig genes in pathway:",length(sig.genes)))
  out <- parallel::mclapply(path.gene,function(x){
    if(!is.na(x)){
      symbol <- na.omit(x[,1])
      symbol = symbol[!symbol==""]
      pval.percent = pvalhyper(all.genes,sig.genes,symbol)
      return(pval.percent)
    }
    else{
      return(c(NA,NA))
    }
  },mc.cores = 8)

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
