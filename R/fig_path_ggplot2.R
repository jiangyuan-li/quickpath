#' title
#'
#'@param
#'@param
#'@return
#'@examples
#'eg
#'@export
fig_path <- function(path.target,out.name,list.info){

  library(ggplot2)

  path.target <- c("path:mmu00562",
                   "path:mmu00564",
                   "path:mmu04070",
                   "path:mmu04072",
                   "path:mmu04151")

  path.info1 <- read.csv("../results/HFvsNF_meth_all.csv",stringsAsFactors = FALSE)
  path.info2 <- read.csv("../results/H2SvsHF_meth_all.csv",stringsAsFactors = FALSE)

  tmp1 <- path.info1[path.info1$name%in% path.target,c(1:2,5)]
  dta1<- tmp1[,2:3]
  dta1$path <- unlist(strsplit(tmp1$path," - M"))[seq(1,2*length(path.target),2)]
  names(dta1)[1] <- "Pathway"
  dta1$percentage = dta1$percentage * 100
  dta1$type = "HF vs. NF"

  tmp2 <- path.info2[path.info2$name%in% path.target,c(1:2,5)]
  dta2<- tmp2[,2:3]
  dta2$path <- unlist(strsplit(tmp2$path," - M"))[seq(1,2*length(path.target),2)]
  names(dta2)[1] <- "Pathway"
  dta2$percentage = dta2$percentage * 100
  dta2$type = "H2S vs. HF"

  path.info3 <- read.csv("../results/H2SvsNF_meth_all.csv",stringsAsFactors = FALSE)
  tmp3 <- path.info3[path.info3$name%in% path.target,c(1:2,5)]
  dta3<- tmp3[,2:3]
  dta3$path <- unlist(strsplit(tmp3$path," - M"))[seq(1,2*length(path.target),2)]
  names(dta3)[1] <- "Pathway"
  dta3$percentage = dta3$percentage * 100
  dta3$type = "H2S vs. NF"

  dat <- rbind(dta1,dta2,dta3)
  #dat <- data.frame(Pathway = rep(sel.kegg.names,2), pval = c(out[,1],out[,2]),type=rep(c("HF vs NF","H9N vs NF"),each= nrow(out)))
  dat$Pathway = factor(dat$Pathway,levels = unique(dat$Pathway)[c(2,5,1,3,4)])
  dat$type <- factor(dat$type,levels = unique(dat$type))

  tiff("../figures/percentage_3comparisons_pathway.tiff", width=2000, height=2400, res=216)
  ggplot(data=dat, aes(x=Pathway,y=percentage,group=type,colour=type,linetype=type,shape=type)) +
    geom_line()+
    geom_point(size=1.5)+
    scale_y_continuous(breaks = c(seq(0,40,5)),
                       labels = c(paste(seq(0,40,5),"%",sep="")))+
    labs(x = "Pathway", y = "Percentage") +
    scale_color_manual("Comparison", values=c("1","2","blue"))+
    scale_linetype_manual("Comparison",
                          values = c(1,2,3))+
    scale_shape_manual("Comparison",
                       values = c(15,16,17))+
    #  scale_x_discrete(labels=c())+
    theme(text = element_text(size=20,face="bold"),
          axis.text.x=element_text(angle=75, hjust=1,vjust=1,size = 20,face="bold"),
          axis.text.y = element_text(face="bold", size=20))
  dev.off()
}
