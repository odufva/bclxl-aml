
FUN_PLOT=function(gene, logicalVectors, namesLV, data=NULL, matrix=NULL, col=NULL, ORDER=F, RANGE=NULL) {
  if(is.null(matrix)&is.null(data))stop("No data to plot, check data/matrix")
  
  GNAME=gsub("N:....:|:::::|DUFVA_", "", gene)
  GNAME=gsub("_", " ", GNAME)
  namesLV=gsub("Cancer_", " ", namesLV)
  
  if(is.null(col)){
    # http://tools.medialab.sciences-po.fr/iwanthue/
    col=c("#d7a85b",
          "#4d56b9",
          "#acb839",
          "#5e2883",
          "#42c87f",
          "#bf4da5",
          "#75b550",
          "#8479e6",
          "#cea632",
          "#5488e3",
          "#d38725",
          "#3e397f",
          "#a4a94e",
          "#be7cde",
          "#4d7122",
          "#8460b5",
          "#62ac6a",
          "#86275d",
          "#43c8ac",
          "#cf483f",
          "#748ed7",
          "#ca6e37",
          "#de91dc",
          "#926a26",
          "#94589d",
          "#822c17",
          "#d76092",
          "#d2745b",
          "#b24258",
          "#d35760")[seq(logicalVectors)]
  }
  
  
  if(!is.null(matrix)){
    gene2=ifelse(grepl("GEXP", gene), gene, paste0("'N:GEXP:", gene, ":::::'"))
    D=as.numeric(read.delim(pipe(paste0("grep -Fw ", gene2, " ", matrix)), row.names = 1, header=F))
  }
  
  if(!is.null(data)){
    D = as.numeric(data[,colnames(data)%in%gene])
  }
  
  bplot_list=lapply(logicalVectors, function(v){
    D[v]
  })
  names(bplot_list)=gsub("_", " ", namesLV)
  
  if(ORDER){
    ord=sapply(bplot_list, median)
    col=col#[order(ord, decreasing = T)]
    bplot_list=bplot_list[c("M0", "M1", "M2", "M3", "M4", "M4eo", "M5", "M6", "M7", "CML", "MDS", "CLL", "pre-B-ALL", "MM", "T-ALL", "HCL", "DLBCL", "FL", "CHL", "NLPHL", "MCL", "MZL", "BL", "MALT", "PTCLNOS", "ALCL", "AITL", "ATL", "CTCL", "HSTCL", "ENKTL")]
  }
  
  plots=FUNCTION_PLOT_LIST(bplot_list, gene, col, ORDER, RANGE)
  return(plots)
}

FUNCTION_PLOT_LIST=function(bplot_list, GNAME, col, ORDER, RANGE){
  
  df=melt(bplot_list)
  
  df$class <- factor(df[,2], levels = unique(as.character(df[,2])),ordered = TRUE)
  
  df$Expression=as.numeric(as.vector(df[,1]))
  p <- ggplot(data=df, aes(x=class, y=Expression)) +  
    geom_boxplot(outlier.shape = NA, fill=col) +
    geom_jitter(width = 0.1, size = 0.1, colour = "grey20") #+
  #stat_compare_means(label = "p.format", label.x.npc = "center", label.y.npc = "top")
  
  p2 <- p +
    
    #theme with white background
    theme_bw() +
    
    # titles
    theme(plot.title = element_text(face="italic", color="black", size=16, hjust=0)) +
    theme(axis.title = element_text(color="black", face=NULL, size=12,angle = 90)) +
    # theme(axis.title.x = element_text(size = 18, angle = 90, color="black")) +
    theme(axis.title.y = element_text(size = 14, angle = 90, color="black", face=NULL)) +
    
    ylab("Expression (log2)") +
    xlab("") +
    labs(title=GNAME) +
    
    
    
    #eliminates background, gridlines, and chart border
    theme(plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    theme(panel.border= element_blank())+
    theme(plot.margin = unit(c(0.1,0.1,0.1,1), "cm"))+
    
    #draws x and y axis line
    theme(axis.line = element_line(),
          axis.line.x = element_line(color="black", size = 0.5),
          axis.line.y = element_line(color="black", size = 0.5)) +
    
    # X - axis text
    theme(axis.text.x = element_text(angle=45, hjust=1, color="black", size = 14, face=NULL),
          axis.text.y = element_text(hjust=1, color="black", size = 12, face=NULL))+ 
    
    # if want to limit to range
    if(!is.null(RANGE))scale_y_continuous(breaks=seq(2,14,2), limits = RANGE)
  
  return(p2)
}

FUN_BINARYFEAT=function(type, annovector){
  logv=annovector%in%type
}

GET_LOGICAL=function(annovector, a=NULL, filterv=NULL, PREFIX=NULL){
  if(sum(filterv)<1&!is.null(filterv))stop("Logical filtering vector is empty")
  if(class(annovector)!="list")stop("Annotation vector should be a list")
  
  binaryfeatures=unlist(lapply(annovector, function(annov, a, filterv, PREFIX){
    
    if(is.null(filterv))filterv=rep(T, length(annov))
    
    if(is.null(a)){  
      a=unique(annov[filterv])
      a=a[!(is.na(a)|a=="na")]
      annov[!filterv]=NA
    }else{
      b=unique(annov[filterv])
      a=a[a%in%b]
    }
    
    # make binary features
    binaryfeats=lapply(a, FUN_BINARYFEAT, annov)
    if(!is.null(PREFIX))names(binaryfeats)=paste0(a ,"_", PREFIX)
    if(is.null(PREFIX))names(binaryfeats)=a
    return(binaryfeats)
  }, a, filterv, PREFIX), recursive = F)
}

