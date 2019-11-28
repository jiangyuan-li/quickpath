#' A line chart for pathways' over-representation level
#'
#'@param path.ids Interested pathway ids
#'@param list.info A list which contains pathway analysis results for all interested groups
#'@param group.info A vector of group names
#'@param criterion Plot is for p value or percentage of DEGs in each pathway
#'@param out.name If specified, the plot will be saved in this file
#'@param path.names It will be used in xlab, if not specified, will use the full name.
#'@return NULL
#'@examples
#' deg = grab_deg_from_cuffdiff(gene_exp.diff, class = "mmu")
#' head(deg)
#' sig.genes = deg$external_gene_name
#' path_res = pathway_analysis(sig.genes, class = "mmu")
#' deg_chow = grab_deg_from_cuffdiff(gene_exp.diff_chow, class = "mmu")
#' sig.genes = deg_chow$external_gene_name
#' path_res_chow = pathway_analysis(sig.genes, class = "mmu")
#'
#' deg_e95 = grab_deg_from_cuffdiff(gene_exp.diff_e95, class = "mmu")
#' sig.genes = deg_e95$external_gene_name
#' path_res_e95 = pathway_analysis(sig.genes, class = "mmu")
#'
#' deg_e105 = grab_deg_from_cuffdiff(gene_exp.diff_e105, class = "mmu")
#' sig.genes = deg_e105$external_gene_name
#' path_res_e105 = pathway_analysis(sig.genes, class = "mmu")
#'
#' list.info = list(path_res[[2]], path_res_chow[[2]], path_res_e95[[2]], path_res_e105[[2]])
#' path.ids = check_pathway_name(pathway, class = "mmu")
#' group.info = c("MCD","Chow","E9.5","E10.5")
#' path.names = c("1st path","2nd path","3rd path","4th path")
#'
#' fig_path(path.ids, list.info, group.info, criterion = "percentage", path.names = path.names)
#' fig_path(path.ids, list.info, group.info, criterion = "pval", path.names = path.names)

#'@export
fig_path <- function(path.ids, list.info, group.info, criterion = c("pval","percentage"),
                     out.name = NULL,
                     path.names = NULL,
                     width=2000, height=2400, res=216){

  n = length(list.info)
  dta.list <- list()
  for(i in 1:n){
    tmp.info = list.info[[i]]
    tmp.dta = tmp.info[tmp.info$name %in% path.ids, c("path", criterion)]

    if(is.null(path.names)){
      tmp.dta$path = unlist(strsplit(tmp.dta$path," - M"))[seq(1,2*length(path.ids),2)]
    }
    else{
      tmp.dta$path = path.names
    }

    if(criterion == "percentage"){
      tmp.dta[,2] = tmp.dta[,2] * 100
    }
    else{
      tmp.dta[,2] = -log(tmp.dta[,2])
    }

    tmp.dta$type = group.info[i]
    dta.list[[i]] = tmp.dta
  }

  final.dta = do.call("rbind", dta.list)
  final.dta$path = factor(final.dta$path)
  final.dta$type = factor(final.dta$type, levels = group.info)
  limit = max(final.dta[,2]+10)

  if(criterion == "percentage"){
    fig <- ggplot(data=final.dta, aes(x=path,y=percentage,group=type,colour=type,linetype=type,shape=type)) +
      geom_line()+
      geom_point(size=1.5)+
      scale_y_continuous(limits=c(0,limit),breaks = c(seq(0,limit,5)),
                         labels = c(paste(seq(0,limit,5),"%",sep="")))+
      labs(x = "Pathway", y = "Percentage") +
      scale_color_manual("Comparison", values=as.character(1:n))+
      scale_linetype_manual("Comparison",
                            values = 1:n)+
      scale_shape_manual("Comparison",
                         values = 15:(14+n))+
      theme(text = element_text(size=20,face="bold"),
            axis.text.x=element_text(angle=75, hjust=1,vjust=1,size = 20,face="bold"),
            axis.text.y = element_text(face="bold", size=20))
  }
  else{
    fig <- ggplot(data=final.dta, aes(x=path,y=pval,group=type,colour=type,linetype=type,shape=type)) +
      geom_line()+
      geom_point(size=1.5)+
      scale_y_continuous(limits=c(0,limit),breaks = c(-log(0.05),seq(0,limit,5),limit+5),
                         labels = c("-log(0.05)",seq(0,limit,5),"Inf"))+
      labs(x = "Pathway", y = "-log pvalue") +
      scale_color_manual("Comparison", values=as.character(1:n))+
      scale_linetype_manual("Comparison",
                            values = 1:n)+
      scale_shape_manual("Comparison",
                         values = 15:(14+n))+
      geom_hline(yintercept = -log(0.05),lty=2)+
      theme(text = element_text(size=20,face="bold"),
            axis.text.x=element_text(angle=75, hjust=1,vjust=1,size = 20,face="bold"),
            axis.text.y = element_text(face="bold", size=20))
  }

  print(fig)

  if(!is.null(out.name)){
    tiff(out.name, width=width, height=height, res=res)
    print(fig)
    dev.off()
  }
}
