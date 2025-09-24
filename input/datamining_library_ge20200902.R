library(Mfuzz)
ge.na.ratio <- function(x){
  sum(is.na(x))/dim(x)[1]/dim(x)[2]
}

NA_threshold_table = function(matrix) {
  na_ratio_in_each_prot = apply(matrix, 1, function(x) {
    sum(is.na(x))/ncol(matrix)
  })
  
  temp = data.frame(sample = names(na_ratio_in_each_prot),
                    na_ratio = na_ratio_in_each_prot,
                    stringsAsFactors = F)
  
  Table1_na = sapply(10:1, function(x) {
    threshold = x/10
    
    prot_choose = temp$na_ratio <= threshold
    prot_num = sum(prot_choose)
    na = ge.na.ratio(matrix[prot_choose,]) %>% round(4)
    
    return(c(threshold = paste0(10*x,"%"),
             protein_num = prot_num,
             NA_ratio = paste0(100*na,"%")))
  }) %>% t()
  
  return(Table1_na)
}

ge.split <- function(data,split,which=1,start=NULL){
  if(is.null(start)){
    sapply(data,function(v){strsplit(v,split)[[1]][which]})
  }else{
    tmp <- sapply(data,function(v){strsplit(v,split)[[1]][1]})
    sapply(tmp,function(v){strsplit(v,start)[[1]][2]})
  }
}


ge.readtable <- function(data,sep = "\t",header = T){
  read.table(data,sep = sep,header = header,stringsAsFactors = F)
}

ge.writetable <- function(data,filename ,sep = "\t",col.names = T,row.names = T,quote = F){
  write.table(data,filename,sep=sep,col.names = col.names,row.names = row.names,quote = quote)
}

ge.plot.density <- function(data){
  plot(density(na.omit(unlist(data))),main="")
  #axis(1,cex.axis = 3)
}

ge.remove.techrep <- function(data,pattern="_repB",method="mean"){
  repB <- names(data)[grepl(pattern, names(data))]
  for (i in repB) {
    repa <- str_split(i,pattern)[[1]][1]
    df1 <- data[,which(names(data) %in% c(repa,i))]
    data <- data[,-which(names(data) %in% c(repa,i))]
    new_mean <- apply(df1, 1, function(x){ifelse(sum(is.na(x))==ncol(df1),NA, mean(as.numeric(x),na.rm=T))} )
    data <- cbind(data,new_mean)
    names(data)[ncol(data)] <- repa
  }
  return(data)
}

ge.plot.techrep.correlation <- function(cor1,cor2,name="pearson_correlation"){
  pdf(paste0(name,".pdf"))
  r <- cor(cor1, cor2, use = "pairwise.complete.obs")   
  smoothScatter(cor1, cor2, nrpoints = 100,cex = 2,
                colramp = colorRampPalette(c(blues9,"orange", "red")),
                main = name, xlab = "repA", ylab = "repB")
  abline(lm(cor1 ~ cor2), col="red", lwd=2, lty=2)
  text(min(cor1,na.rm = T)*1.3,max(cor2,na.rm = T)*0.8,labels =paste0( "r =", as.character(round(r,4))),cex = 1.2)
  dev.off()
}

ge.plot.pool.correlation <- function(data,name="bio_cor",method="circle",height = 7, width = 7){
  library(corrplot)
  df_cor <- data.frame(data)
  pdf(paste0(name,".pdf"),height = height, width = width)
  mycor=cor(df_cor, use = "pairwise.complete.obs")
  corrplot(mycor, method=method,type = "upper",tl.col = "black",tl.srt = 45, tl.cex = 1.5)
  dev.off()
}




ge.plot.boxplot <- function(data,x,y,type,filename,title="boxplot"){
  a <- ggplot(data=data, aes(x =x, y =y ,color=type,group=type)) +
    geom_jitter(alpha = 0.3,size=3) +
    geom_boxplot(alpha = .5,size=1)+
    labs(x="sample",y="value",fill= "type")+
    ggtitle(title)+
    theme_bw() + 
    theme(panel.border = element_blank())+
    theme(axis.line = element_line(size=1, colour = "black")) +
    theme(panel.grid =element_blank())+  
    theme(axis.text = element_text(size = 15,colour = "black"),text = element_text(size = 15,colour = "black"))+
    theme(axis.text.x = element_text( hjust = 1,angle = 45))
  ggsave(paste0(filename, ".pdf"),plot=a,width=8,height=8)
}

# p <- ggboxplot(df1, x="dose", y="len", color = "dose", 
#                palette = c("#00AFBB", "#E7B800", "#FC4E07"), 
#                add = "jitter", shape="dose")
# my_comparisons <- list(c("0.5", "1"), c("1", "2"), c("0.5", "2"))
# p+stat_compare_means(comparisons = my_comparisons)+
#   stat_compare_means(label.y = 50)


ge.plot.valcano <- function(data, title,fd=1,pvalue=0.05){
  df.t <- data
  cut.fd <- fd
  pdf(paste0(title, "_volcano.pdf"))
  plot(df.t$fd, -log10(df.t$P_value_adjust), col="#00000033", pch=19,
       xlab=paste("log2 (fold change)"),
       ylab="-log10 (P_value_adjust)",
       main= title)
  
  up <- subset(df.t, df.t$P_value_adjust < pvalue & df.t$fd > cut.fd)
  down <- subset(df.t, df.t$P_value_adjust< pvalue & df.t$fd< as.numeric(cut.fd*(-1)))
  write.csv(up,file = paste0(title, "_up.csv"))
  write.csv(down,file = paste0(title, "_down.csv"))
  points(up$fd, -log10(up$P_value_adjust), col=1, bg = brewer.pal(9, "YlOrRd")[6], pch=21, cex=1.5)
  points(down$fd, -log10(down$P_value_adjust), col = 1, bg = brewer.pal(11,"RdBu")[9], pch = 21,cex=1.5)
  abline(h=-log10(pvalue),v=c(-1*fd,fd),lty=2,lwd=1)
  
  dev.off()
}


ge.plot.pca <- function(data,type,title="",ptColors=NULL,label2=NULL,width=12,height=8){ 
  M <- t(data)
  M <- apply(M,2,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)})
  clnames <- row.names(data)
  M[is.na(M)] <- 0
  m1 <- prcomp(M);
  Y  <- scale(M, m1$center, m1$scale) %*% m1$rotation 
  Y  <- Y[,c(1,2)]
  
  Y <- data.frame(Y,type);
  colnames(Y) <- c("PC1","PC2","label")
  eigs <- m1$sdev^2
  percentages <- eigs[1:2] / sum(eigs)
  p <- ggplot(Y, aes(x=PC1, y=PC2, colour=label)) + geom_point(size=4)
  p <- p + theme(  panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   axis.text = element_text(size = 30,color = "black"),
                   panel.border = element_blank(),
                   axis.line.x = element_line(color="black", size = 0.25),
                   axis.line.y = element_line(color="black", size = 0.25),
                   plot.title   = element_text(size=30),
                   axis.title   =element_text(size=30),
                   panel.background = element_blank())
  
  strLabx <- sprintf("PC1(%4.2f%%)",percentages[1]*100)
  p <- p +  labs(x =strLabx,y = sprintf("PC2(%4.2f%%)",percentages[2]*100),
                 title =sprintf("PCA:%d features",length(clnames)))
  if(!is.null(ptColors)){
    p <- p +   scale_colour_manual(values=ptColors)
  }
  if(!is.null(label2)){
    p <- p +   geom_text(aes(label=label2,vjust = -0.8, hjust = 0.5,size=0.5),show.legend = FALSE)
  }
  
  ggsave(paste0(title,"_pca.pdf"),plot =p ,width=width,height=height,device = NULL)
}


ge.plot.tsne <- function(data,type,title="",label=NULL){
  col2=brewer.pal(9,"Set1")[1:length(unique(type))]
  cl <- data.frame(col2,row.names =unique(type),stringsAsFactors = F)
  cl2 <- cl[match(type,row.names(cl)),1]
  
  df10 <- data
  df10[is.na(df10)] <- 0
  df10 <- t(apply(df10, 1, scale))
  colnames(df10) <- type
  df11.tsne <- Rtsne(t(df10), dims = 2, perplexity = (ncol(data)-1)/3-1, verbose = T , check_duplicates = FALSE)
  pdf(paste0(title,"_TNSE.pdf"))
  if(is.null(label)){
  plot(df11.tsne$Y,col=cl2, main = "tsne", pch = 20,cex=2,cex.axis=2,cex.lab=2)
  legend("topright",legend=row.names(cl), fill=cl$col2, lty=1,lwd=1)
  }else{
    plot(df11.tsne$Y,col=cl2, main = "tsne", pch = 20,cex=2,cex.axis=2,cex.lab=2)
    text(df11.tsne$Y,pos = 1, labels = label, col= "DimGrey",cex = 0.8)
    legend("topright",legend=row.names(cl), fill=cl$col2, lty=1,lwd=1)
  }
  dev.off()
}


ge.plot.umap<- function(data,type,title="",label=NULL){
  col2=brewer.pal(9,"Set1")[1:length(unique(type))]
  cl <- data.frame(col2,row.names =unique(type),stringsAsFactors = F)
  cl2 <- cl[match(type,row.names(cl)),1]
  
  df10 <- data
  df10[is.na(df10)] <- 0
  df10 <- t(apply(df10, 1, scale))
  colnames(df10) <- type
  df.umap <- umap(t(df10),n_neighbors=ncol(data)-1)
  pdf(paste0(title,"_UMAP.pdf"))
  if(is.null(label)){
  plot(df.umap$layout,col = cl2, main = "umap", pch = 20,cex=2,cex.axis=2,cex.lab=2)
  legend("topright",legend=row.names(cl), fill=cl$col2, lty=1,lwd=1)
  }else{
    plot(df.umap$layout,col = cl2, main = "umap", pch = 20,cex=2,cex.axis=2,cex.lab=2)
    text(df.umap$layout, pos = 1, labels = label, col= "DimGrey",cex = 0.8)
    legend("topright",legend=row.names(cl), fill=cl$col2, lty=1,lwd=1)
  }
  dev.off()
}



ge.plot.volcano <- function(data, group1, group2, fc= 1, pvalue = 0.05, str1= "grp1",str2= "grp2",pair=F,adjust.bool=T) {
  df8 <- 2^data
  df8[is.na(df8)] <- min(df8,na.rm = T)*0.8
  df8$fd <-
    apply(df8, 1, function(x)
      log2((mean(x[group1], na.rm = T) / mean(x[group2], na.rm = T))))
  x <- c(0.0, 0.0)
  df9 <- data
  df9[is.na(df9)] <- min(df9,na.rm = T)*0.8
  
  df8$P_value <- apply(df9,1,function(y) {p_try = tryCatch(t.test(y[group1],
                                             y[group2],
                                             paired = pair,
                                             var.equal = F)$p.value,
                                      error = function(x) NA)})
  
  df8$P_value_adjust <- p.adjust(df8$P_value, method = "BH")
  if(adjust.bool){
    df8$P <- df8$P_value_adjust
    y.lab <- "-log10 (adjust P)"
  }else{
    df8$P <- df8$P_value
    y.lab <- "-log10 (P_value)"
  }
  pdf(paste0(str1, "_", str2, "_volcano.pdf"),
      width = 4,
      height = 4,
  )
  plot(
    df8$fd,
    -log10(df8$P),
    col = "#00000033",
    pch = 19,
    xlab = paste("log2 (fold change)"),
    ylab = y.lab,
    #xlim = c(-4, 4),
    main = paste0(str1, " vs ", str2)
  )
  
  up <- subset(df8, df8$P < pvalue & df8$fd > fc)
  down <- subset(df8, df8$P < pvalue & df8$fd < -1*fc)
  write.csv(up, file = paste0(str1, "_", str2, "_up_volcano.csv"))
  write.csv(down, file = paste0(str1, "_", str2, "_dw_volcano.csv"))
  points(
    up$fd,
    -log10(up$P),
    col = 1,
    bg = brewer.pal(9, "YlOrRd")[6],
    pch = 21,
    cex = 1.5
  )
  points(
    down$fd,
    -log10(down$P),
    col = 1,
    bg = brewer.pal(11, "RdBu")[9],
    pch = 21,
    cex = 1.5
  )
  abline(
    h = -log10(pvalue),
    v = c(-1*fc, fc),
    lty = 2,
    lwd = 1
  )
  dev.off()
  print(paste0("up:",nrow(up)))
  print(paste0("down:",nrow(down)))
}

ge.plot.bar <- function(data,sample,value,group=group,title="",xlab="sample",ylab="value"){
  a <- ggplot(data,aes(x=sample,y=value,group=group))+ 
    geom_bar(position = "dodge",stat = "identity",width =0.8,alpha=0.8,aes(fill=group))+
    ggtitle(paste0(title,"_barplot"))+
    xlab(xlab)+
    ylab(ylab)+
    theme(legend.text = element_text(size = 15,color = "black"),legend.position = 'top',
          legend.title = element_text(size=15,color="black") ,
          panel.grid.major =element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))+
    theme(panel.grid =element_blank())+
    theme(axis.text = element_text(size = 10,color = "black"))+
    theme(axis.text.x = element_text( hjust = 1,angle = 45))+
    theme(plot.subtitle=element_text(size=30, hjust=0, color="black"))+
    theme(axis.title.x=element_text(size=17, hjust=0.5, color="black"))+
    theme(axis.title.y=element_text(size=17, hjust=0.5, color="black")) + geom_text(aes(x=sample,label=value,vjust = -0.8, hjust = 0.5),position = "dodge",stat = "identity",show.legend = FALSE)
  ggsave(paste0(title,"_barplot.pdf"),plot=a,width=10,height=8)
}


ge.plot.line <- function(data,sample,value,group=group,title="",xlab="sample",ylab="value"){
  a <- ggplot(data,aes(x=sample,y=value,group=group,color=group))+ 
    geom_line()+
    geom_point()+
    ggtitle(paste0(title,"_lineplot"))+
    xlab(xlab)+
    ylab(ylab)+
    theme(legend.text = element_text(size = 15,color = "black"),legend.position = 'top',
          legend.title = element_text(size=15,color="black") ,
          panel.grid.major =element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))+
    theme(panel.grid =element_blank())+
    theme(axis.text = element_text(size = 10,color = "black"))+
    theme(axis.text.x = element_text( hjust = 1,angle = 45))+
    theme(plot.subtitle=element_text(size=30, hjust=0, color="black"))+
    theme(axis.title.x=element_text(size=17, hjust=0.5, color="black"))+
    theme(axis.title.y=element_text(size=17, hjust=0.5, color="black"))+ geom_text(aes(label=group,vjust = -0.8, hjust = 0.5),show.legend = FALSE)
  ggsave(paste0(title,"_lineplot.png"),plot=a,width=10,height=8)
}


# ge.plot.vioplot <- function(sample1,sample2,title="",xlab="sample",ylab="value"){
# pdf(paste0(title, "_violin.pdf"))
# vioplot(sample1 ,sample2  ,
#          areaEqual=FALSE, 
#         # rectCol= color, col= color,
#         lineCol=c("black", "black"),
#         border=c("black","black"),
#         names=c("DIANN","DIANN_quantile"),
#         main="biological replicates", xlab=xlab, ylab=ylab,plotCentre = "point")
# dev.off()
# }


ge.mfuzz.cselection <- function(data,range=seq(5,50,5),repeats = 5){
  df3a<-as.matrix(data)
  df3Ex<- ExpressionSet(assayData = df3a)
  if(interactive()){
    df3F <- filter.NA(df3Ex)
    df3F <- fill.NA(df3F)
    df3F <- standardise(df3F)
  }
  
  df3F <- filter.NA(df3F)
  m<-mestimate(df3F)
  cselection(df3F,m=m,crange = range,repeats = repeats,visu = T)
  return(df3F)
}

ge.mfuzz.getresult <- function(expressionSet, pic,time.label,filename,anova=F,alldata=NULL){
  cl <- mfuzz(expressionSet,c=pic,m=1.25)
  dir.create(path=filename,recursive = TRUE)
  pdf(paste0(filename,".pdf"))
  mfuzz.plot2(expressionSet, cl=cl,time.labels=time.label,mfrow=c(4,4),centre=TRUE,x11=F,centre.lwd=0.2)#min.mem=0.99
  dev.off()
  
  for(i in 1:pic){
    potname<-names(cl$cluster[unname(cl$cluster)==i])
    write.csv(cl[[4]][potname,i],paste0(filename,"/mfuzz_",i,".csv"))
  }
  if(anova){
    for(ii in 1:pic){
      potname<-names(cl$cluster[unname(cl$cluster)==ii])
      tmp <- data.frame(label=time.label,t(alldata[potname,]))
      anova <- c()
      for (n in 1:length(potname)) {
        aov<-(summary(aov(tmp[,n+1] ~ label,tmp))[[1]])$`Pr(>F)`[1]
        anova <- c(anova,aov)
      }
      anova.adjust <-p.adjust(anova, method="BH")
      newdf <- data.frame(prot=names(tmp)[-1],anova,anova.adjust)
      newdf2 <- newdf[newdf$anova.adjust<0.05,]
      write.csv(newdf2,paste0(filename,"/mfuzz_anova_",ii,".csv"),row.names = F)
    }
  }
}


