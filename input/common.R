library("FactoMineR")
library("factoextra")
library("grid")
library("data.table")
library("ggplot2")
library(pheatmap);
library("RColorBrewer");
library("tsne");


###tSNE plot DE: dataframe with rows as sample and columns as features,a column with name 'label' is required, represents the label of samples
# lblColors: a vector containing colors for each label, like  lblColors=c(M = "forestgreen", N = "gray0", P="firebrick", C= "red",A="blue")
drawTSNE <- function(DF,ptColors,rowNormalization=F,colNormalization=F,perplexity=10,strTitle='tSNE'){
    M <- DF[,colnames(DF)!='label']
    if(rowNormalization){M <- data.frame(t(apply(M,1,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)})))}
    if(colNormalization){M <- apply(M,2,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)})}
    M[is.na(M)] <- 0
    indx <- match('label',colnames(M))
    clnames <- colnames(DF)[colnames(DF)!='label']
    tsn = tsne(M,perplexity =perplexity)
    cnames <- colnames(M)
    tsn <- data.frame(tsn,DF$label)
    colnames(tsn) <- c("X","Y","label")
    rownames(tsn) <- rownames(M)
    #tsn <- tsn[-c(which(tsn$X==min(tsn$X)),which(tsn$Y==min(tsn$Y))),]
    #tsn <- tsn[-c(which(tsn$X==max(tsn$X)),which(tsn$Y==max(tsn$Y))),]
   #lblColors <- c(N='#537e35',M='#e17832',A='#f5b819',C='#5992c6',P='#282f89',W='mediumorchid3')

   #lblColors <- c(A='#537e35',M='#e17832',N='#f5b819',B='#5992c6',C='#282f89',W='mediumorchid3')
   p <- ggplot(tsn, aes(x=X, y=Y, colour=label)) + geom_point(size=4)
   p <- p + theme(  panel.grid.major = element_blank(),
	    panel.grid.minor = element_blank(),
	    panel.border = element_blank(),
	     plot.title = element_text(size=15),
	    axis.line.x = element_line(color="black", size = 0.5),
	    axis.line.y = element_line(color="black", size = 0.5),
	    panel.background = element_blank())
   p <-  p +  labs(title=strTitle)
   p <- p +   scale_colour_manual(values=ptColors)
   p
}

drawPCA<- function(DF,ptColors,rowNormalization=F,colNormalization=F,strTitle=NULL){ 
## M is a matrix or dataframe, rows are samples and columns are features, rownames are sample names
   M <- DF[,colnames(DF)!='label']
   if(rowNormalization){
      M <- data.frame(t(apply(M,1,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)})))
      #print('Normalization by row is done!')
   }
   if(colNormalization){
      M <- apply(M,2,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)})
   }
   clnames <- colnames(DF)[colnames(DF)!='label']
   M[is.na(M)] <- 0
   m1 <- prcomp(M,colNormalization);
   Y  <- scale(M, m1$center, m1$scale) %*% m1$rotation 
   Y  <- Y[,c(1,2)]
   
   Y <- data.frame(Y,DF$label);
   colnames(Y) <- c("PC1","PC2","label")
   if(is.null(strTitle)){
      strTitle <- sprintf("PCA:%d features",length(clnames))
	}
   eigs <- m1$sdev^2
   percentages <- eigs[1:2] / sum(eigs)
   # lblColors <- c(N='#537e35',M='#e17832',A='#f5b819',C='#5992c6',P='#282f89',W='mediumorchid3')
   #lblColors <- c(training='#537e35',validation='#e17832',A='#f5b819',C='#5992c6',P='#282f89',W='mediumorchid3')
   p <- ggplot(Y, aes(x=PC1, y=PC2, colour=label)) + geom_point(size=4)
   p <- p + theme(  panel.grid.major = element_blank(),
	    panel.grid.minor = element_blank(),
	    panel.border = element_blank(),
	    axis.line.x = element_line(color="black", size = 0.25),
	    axis.line.y = element_line(color="black", size = 0.25),
	    plot.title   = element_text(size=16),
	    panel.background = element_blank())
            
   strLabx <- sprintf("PC1(%4.2f%%)",percentages[1]*100)
   p <- p +  labs(x =strLabx,y = sprintf("PC2(%4.2f%%)",percentages[2]*100),
	      title =strTitle)
   p <- p +   scale_colour_manual(values=ptColors)
  
   p
}

# parse gene symbol from uniprot website
# acc: protein access Id like ' "O00264" "O00468" "O14773" "O14979" "O15230" "O43175" "O43707"
# return: the gene symbol for the acc
geneName <- function(acc){
   kk <- read.table(sprintf('https://www.uniprot.org/uniprot/%s.fasta',acc),sep="\n",stringsAsFactors=F,header=F)
   a <- strsplit(kk[1,]," ")[[1]]
   b <- a[grepl("GN=",a)]
   strsplit(b,"=")[[1]][2]
}
getSubset <- function(Labels,M0,missRate,imputation=F){
   lbls <- Labels
   lbls <- if(length(Labels)==1){lbls <- unlist(strsplit(lbls,""))}else{lbls <- Labels}
   tmp  <- M0[M0$label %in% lbls,]
   R0   <- apply(tmp,2,function(v){sum(is.na(v))/length(v)*100})
   tmp  <- tmp[,R0<=missRate]

   #if(!imputation){tmp[is.na(tmp)] <- 0}
   if(imputation){ #### kNN imputation 
       K=3
       t1 <- sapply(1:dim(tmp)[1],function(i) {
           lbl <- tmp[i, 'label']
           v0  <- tmp[i, colnames(tmp) != 'label']
           feature0 <- names(v0)[is.na(v0)]
           if(length(feature0)>0){
		     k0  <- tmp[-i,]
		     k0  <- k0[k0$label==lbl,-1]
		     k0  <- rbind(v0,k0)
		     k0 <- apply(k0,2,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)})
		     k0[is.na(k0)] <- 0

		     D0  <- as.matrix(dist(k0))
		     d0  <- order(D0[1,])[2:(K+1)]
		     if(length(feature0)==1){
			v0[feature0] <- mean(tmp[d0,feature0],na.rm=T)
		     }else{
			v0[feature0] <- apply(tmp[d0,feature0],2,function(v){mean(v,na.rm=T)})
		     }
           }    
           v0
       })
      
      t1 <- apply(t(t1),2,unlist)    
      rownames(t1) <- rownames(tmp)
      t1 <- data.frame(t1)
      clnames <- colnames(t1)
      t1$label <- tmp$label
      t1 <- t1[,c('label',clnames)]
      #print(sprintf("get subset for %s with %d rows and %d features",paste0(lbls,collapse=","),dim(t1)[1],dim(t1)[2]))
      tmp <- t1
   }
   return(tmp)
}


 RFScore <- function(feature,M1=t1,nFolds=10,it=10){
   #tmp <- data.frame(t(apply(M1[,feature],1,function(v0){v <- as.numeric(v0); (v-mean(v))/sd(v)})))
   tmp <- M1[,c('label',feature)]
   clsses <- unique(tmp$label)
   #rownames(tmp) <- sapply(1:dim(tmp)[1],function(v){paste0('R',v)})
   sampleNumber <- sapply(clsses,function(v){sum(tmp$label==v)})
   names(sampleNumber) <- clsses
   maxLabel <- names(sampleNumber)[max(sampleNumber)==sampleNumber]
   minLabel <- names(sampleNumber)[min(sampleNumber)==sampleNumber]
   size <- min(sampleNumber)
   
   avgAcc <- 0;
   probs <- rep(0,max(sampleNumber))
   names(probs) <- rownames(tmp)[tmp$label==maxLabel]
   its <- 0
   indx2 <- rownames(tmp)[tmp$label==minLabel]

   prediction <- matrix(0,nrow=dim(M1)[1],ncol=2)
   rownames(prediction) <- rownames(M1)

   while(sum(probs==0)>0){
      indx1 <- names(sort(probs))[1:size]
      probs[indx1] <- probs[indx1]+1
           tmp1 <- tmp[c(indx1,indx2),]
	   for(fd in 1:nFolds){
               folds <- createFolds(tmp1$label,nFolds)
               for(fold in folds){
                   valids <- tmp1[fold,]
		   rownames(valids) <- rownames(tmp1)[fold]

                   trains <- tmp1[setdiff(1:dim(tmp1)[1],fold),]
	           trains$label <- as.factor(trains$label)
		   
                   tmpRF <- randomForest(label ~ . ,data=trains,importance=T,ntree=1000,nodesize=5)
	           predicted <- predict(tmpRF,valids,type='prob')
		   prediction[rownames(predicted),] <- predicted+prediction[rownames(predicted),]
		  
               }
           }
     }#while
      colnames(prediction) <- colnames(predicted)
      prediction <- data.frame(t(apply(prediction,1,function(v){v/sum(v)})))
      prediction$predicted <- as.character(apply(prediction,1,function(v){names(v)[v==max(v)]}))
      prediction$observed  <- M1$label
      rtv <-sprintf("%s  %4.3f",paste0(feature,collapse=","),sum(prediction$observed==prediction$predicted)/dim(M1)[1]*100)
      #prediction
      rtv
 }
#############################################################################################
RFPredict <- function(trainM,testM){
   
   clsses <- unique(trainM$label)
   #rownames(tmp) <- sapply(1:dim(tmp)[1],function(v){paste0('R',v)})
   sampleNumber <- sapply(clsses,function(v){sum(tmp$label==v)})
   names(sampleNumber) <- clsses
   maxLabel <- names(sampleNumber)[max(sampleNumber)==sampleNumber]
   minLabel <- names(sampleNumber)[min(sampleNumber)==sampleNumber]
   size <- min(sampleNumber)
   
   avgAcc <- 0;
   probs <- rep(0,max(sampleNumber))
   names(probs) <- rownames(tmp)[tmp$label==maxLabel]
   its <- 0
   indx2 <- rownames(tmp)[tmp$label==minLabel]

   prediction <- matrix(0,nrow=dim(M1)[1],ncol=2)
   rownames(prediction) <- rownames(M1)

   while(sum(probs==0)>0){
      indx1 <- names(sort(probs))[1:size]
      probs[indx1] <- probs[indx1]+1
           tmp1 <- tmp[c(indx1,indx2),]
	   for(fd in 1:nFolds){
               folds <- createFolds(tmp1$label,nFolds)
               for(fold in folds){
                   valids <- tmp1[fold,]
		   rownames(valids) <- rownames(tmp1)[fold]

                   trains <- tmp1[setdiff(1:dim(tmp1)[1],fold),]
	           trains$label <- as.factor(trains$label)
		   
                   tmpRF <- randomForest(label ~ . ,data=trains,importance=T,ntree=1000,nodesize=5)
	           predicted <- predict(tmpRF,valids,type='prob')
		   prediction[rownames(predicted),] <- predicted+prediction[rownames(predicted),]
		  
               }
           }
     }#while
      colnames(prediction) <- colnames(predicted)
      prediction <- data.frame(t(apply(prediction,1,function(v){v/sum(v)})))
      prediction$predicted <- as.character(apply(prediction,1,function(v){names(v)[v==max(v)]}))
      prediction$observed  <- M1$label
      rtv <-sprintf("%s  %4.3f",paste0(feature,collapse=","),sum(prediction$observed==prediction$predicted)/dim(M1)[1]*100)
      #prediction
      rtv
 }


###############################################################################################################
drawVolcano <- function(df1,pdfPath,strTitle="volcano plot",outFile=NULL){ 
   label <- df1$label
   lbls = unique(df1$label)
   table(df1$label) # A143,C75
   fc <- apply(2^df1[,colnames(df1)!='label'],2,function(x){
       log2( mean(na.omit(x[label==lbls[1]])) /mean(na.omit(x[label==lbls[2]])))
    })
   
   df1[is.na(df1)] <- 0
   pValue <- apply(df1[,colnames(df1)!='label'], 2, function(v) {
       p1 <- t.test(v[label == lbls[1]], v[label == lbls[2]], paired = F, var.equal = F)$p.value
       p1 <-  p.adjust(p1,method="BH")
       p1
    }) 

  pdf(pdfPath)
     plot(fc, -log10(pValue), col = '#00000033', pch = 19,xlab = 'log2(FC)', ylab = '-log10(p-value)', main = strTitle)
     abline(h = 1.3, v = c(-log2(1.5),log2(1.5)), lty = 2, lwd = 1)
  
     up  <- fc >= log2(1.5) & pValue <= 0.05
     points(fc[up], -log10(pValue[up]), col = 1,bg = brewer.pal(9,"YlOrRd")[6], pch = 21, cex = 2)
 
     down <- fc <= -log2(1.5) & pValue <= 0.05
     points(fc[down], -log10(pValue[down]), col = 1,bg = brewer.pal(11,"RdBu")[9], pch = 21, cex = 2)
  dev.off()
  
  name = list(up = data.frame(prot = colnames(df1[,colnames(df1)!='label'])[up],fc = fc[up], p_value = pValue[up],
                            type = rep("Upregulation",sum(up)),stringsAsFactors=F),
            down = data.frame(prot = colnames(df1[,colnames(df1)!='label'])[down],fc = fc[down], p_value = pValue[down],
                              type = rep("Downregulation",sum(down))),stringsAsFactors=F)
   name1 = rbind(name[[1]],name[[2]])
   rownames(name1) <- 1:dim(name1)[1]
   
   if( !is.null(outFile) ){
     write.table(name1,file=outFile,sep="\t",col.names=T,row.names=F,quote=F)
   }
   #write.xlsx(name1, "Table_1_TPD_volcano_prot_AC10_fc1.5_190131.xlsx",showNA = T,row.names = F)
   sum(up)
   sum(down)
   return(name1)
}

drawUMAP <- function(M1,ptColors, strTitle="UMAP",rowNormalization=T,colNormalization=F){
     
     if(!'label' %in% colnames(M1)){
        print('A column with named label must existed in data frame')
	return(NULL)
     }

     tmp <- M1[,colnames(M1)!='label']
     if(rowNormalization){
        tmp <- data.frame(t(apply(tmp,1,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)})),stringsAsFactors=F)
        rownames(tmp) <- rownames(M1)
     }
     if(colNormalization){    tmp <- apply(tmp,2,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)})   }
     	tmp[is.na(tmp)] <- 0

     obj = umap(d=tmp,method='naive')
     clnames <- colnames(tmp)
     df1 <- data.frame(obj$layout)
     df1$label <- M1$label
     colnames(df1) <- c('X','Y','label')
     
     p <- ggplot(df1, aes(x=X, y=Y, colour=label)) + geom_point(size=4)
     p <- p + theme(  panel.grid.major = element_blank(),
	    panel.grid.minor = element_blank(),
	    panel.border = element_blank(),
	    plot.title   = element_text(size=16),
	    axis.line.x = element_line(color="black", size = 0.5),
	    axis.line.y = element_line(color="black", size = 0.5),

	    panel.background = element_blank())
            
   #strLabx <- sprintf("PC1(%4.2f%%)",percentages[1]*100)
   p <- p +  labs(title =strTitle)
   p <- p +   scale_colour_manual(values=ptColors)
   p
 }

Volcano <- function(df1,pdfPath,outFile=NULL,thresholdFC=1.5,thresholdPValue=0.05,strTitle="volcano plot"){ 
   label <- df1$label
   lbls = unique(df1$label)
   table(df1$label) # A143,C75
   fc <- apply(2^df1[,colnames(df1)!='label'],2,function(x){
       log2( mean(na.omit(x[label==lbls[1]])) /mean(na.omit(x[label==lbls[2]])))
    })
   
   df1[is.na(df1)] <- 0
   pValue <- apply(df1[,colnames(df1)!='label'], 2, function(v) {
       p1 <- t.test(v[label == lbls[1]], v[label == lbls[2]], paired = F, var.equal = F)$p.value
       p1 <-  p.adjust(p1,method="BH")
       p1
    }) 
  if(!is.null(pdfPath)){
	   pdf(pdfPath)
	     plot(fc, -log10(pValue), col = '#00000033', pch = 19,xlab = 'log2(FoldChange)', ylab = '-log10(P-value)', main = strTitle)
	     abline(h = -log10(thresholdPValue), v = c(-log2(thresholdFC),log2(thresholdFC)), lty = 2, lwd = 1)
	  
	     up  <- fc >= log2(thresholdFC) & pValue <= thresholdPValue
	     points(fc[up], -log10(pValue[up]), col = 1,bg = brewer.pal(9,"YlOrRd")[6], pch = 21, cex = 2)
	 
	     down <- fc <= -log2(thresholdFC) & pValue <= thresholdPValue
	     points(fc[down], -log10(pValue[down]), col = 1,bg = brewer.pal(11,"RdBu")[9], pch = 21, cex = 2)
	  dev.off()
  }

  name = list(up = data.frame(prot = colnames(df1[,colnames(df1)!='label'])[up],fc = fc[up], p_value = pValue[up],
                            type = rep("Upregulation",sum(up)),stringsAsFactors=F),
            down = data.frame(prot = colnames(df1[,colnames(df1)!='label'])[down],fc = fc[down], p_value = pValue[down],
                              type = rep("Downregulation",sum(down))),stringsAsFactors=F)
   name1 = rbind(name[[1]],name[[2]],stringsAsFactors=F)
   rownames(name1) <- 1:dim(name1)[1]
   name1 <- name1[order(abs(name1$fc),decreasing=T),]
   if( !is.null(outFile) ){
     write.table(name1,file=outFile,sep="\t",col.names=T,row.names=F,quote=F)
   }
   #write.xlsx(name1, "Table_1_TPD_volcano_prot_AC10_fc1.5_190131.xlsx",showNA = T,row.names = F)
   sum(up)
   sum(down)
   return(name1)
 }

RFScore1 <- function(M1,nFolds=10,itors=10){
   if(! 'label' %in% colnames(M1)){
     print('A column with name "label" is required for labeling the class of each sample')
     return(NULL)
   }
   clsses <- unique(M1$label)
   #rownames(tmp) <- sapply(1:dim(tmp)[1],function(v){paste0('R',v)})
   sampleNumber <- sapply(clsses,function(v){sum(M1$label==v)})
   names(sampleNumber) <- clsses
   maxLabel <- names(sampleNumber)[max(sampleNumber)==sampleNumber]
   minLabel <- names(sampleNumber)[min(sampleNumber)==sampleNumber]
   size <- min(sampleNumber)
   
   probs <- rep(0,max(sampleNumber))
   names(probs) <- rownames(M1)[M1$label==maxLabel]
   
   indx2 <- rownames(M1)[M1$label==minLabel]

   prediction <- matrix(0,nrow=dim(M1)[1],ncol=2)
   rownames(prediction) <- rownames(M1)

   while(sum(probs==0)>0){
      indx1 <- names(sort(probs))[1:size]
      probs[indx1] <- probs[indx1]+1
      
      tmp1 <- M1[c(indx1,indx2),]
      for(i in 1:itors){
          folds <- createFolds(tmp1$label,nFolds)
          for(fold in folds){
              valids <- tmp1[fold,]
	      rownames(valids) <- rownames(tmp1)[fold]
              trains <- tmp1[setdiff(1:dim(tmp1)[1],fold),]
	      trains$label <- as.factor(trains$label)
              tmpRF <- randomForest(label ~ . ,data=trains,importance=T,ntree=1000,nodesize=5)
	      predicted <- predict(tmpRF,valids,type='prob')
	      prediction[rownames(predicted),] <- predicted+prediction[rownames(predicted),]
         }
           }
     }#while
      colnames(prediction) <- colnames(predicted)
      prediction <- data.frame(t(apply(prediction,1,function(v){v/sum(v)})))
      prediction$predicted <- as.character(apply(prediction,1,function(v){names(v)[v==max(v)]}))
      prediction$observed  <- M1$label
      feature <- colnames(M1)[colnames(M1) !='label']
      rtv <-sprintf("%s  %4.3f",paste0(feature,collapse=","),sum(prediction$observed==prediction$predicted)/dim(M1)[1]*100)
      print(rtv)
      #prediction
      rtv
}
