setwd("I:/RRA")
list.files(path = "I:/RRA")
files=c("GSE104168 PR8-A549 limma DEG.csv","GSE104168 NY238-A549 limma DEG.csv",
        "GSE106279 PR8-A549 limma DEG.csv","GSE32139 Udorn-hAEC limma DEG.csv",
        "GSE37571 CA04-Calu3 limma DEG.csv","GSE61517 Brisbane-BEAS2B limma DEG.csv",
        "GSE61517 Perth-BEAS2B limma DEG.csv","GSE61517 Udorn-BEAS2B limma DEG.csv",
        "GSE71766 WSN-BEAS2B limma DEG.csv")
pcgenes<-read.csv("HGCN_NCBI.csv")
upList=list()
downList=list()
allFCList=list()
for(i in 1:length(files)){
  inputFile=files[i]
  rt=read.csv(inputFile,header=T)
  rt=rt[,-1]
  rt<-aggregate(x = rt[,2:ncol(rt)],                
                by = list(rt$EntrezID),              
                FUN = mean)
  rt<-merge(rt,pcgenes,by="Group.1")
  rt<-rt[order(rt$logFC),]
  header=unlist(strsplit(inputFile,"_"))
  rt[,1] <- as.character(rt[,1])
  downList[[header[1]]]=as.vector(rt[,1])
  upList[[header[1]]]=rev(as.vector(rt[,1]))
  fcCol=rt[,1:2]
  colnames(fcCol)=c("Gene",header[[1]])
  allFCList[[header[1]]]=fcCol
}

mergeLe=function(x,y){
  merge(x,y,by="Gene",all=T)}
newTab=Reduce(mergeLe,allFCList)
newTab[is.na(newTab)]=0
rownames(newTab)=newTab[,1]
newTab=newTab[,2:ncol(newTab)]
write.csv(newTab,"logFC matrix of 9 instances.csv")

library(RobustRankAggreg)
library(org.Hs.eg.db)
library(annotate)
padj=0.05
logFC=1
upMatrix = rankMatrix(upList)
upAR = aggregateRanks(rmat=upMatrix)
colnames(upAR)=c("EntrezID","Pvalue")
upAdj=p.adjust(upAR$Pvalue,method="bonferroni")
upXls=cbind(upAR,adjPvalue=upAdj)
upFC=newTab[as.vector(upXls[,1]),]
upXls=cbind(upXls,logFC=rowMeans(upFC))
write.table(upXls,file="upProteinCoding-EntrezID.xls",sep="\t",quote=F,row.names=F)
upSig=upXls[(upXls$adjPvalue<padj & upXls$logFC>logFC),]
up<-upSig$EntrezID
upsymbol<-select(org.Hs.eg.db, 
                 keys = up,
                 columns = c("ENTREZID","SYMBOL"),
                 keytype = "ENTREZID")
colnames(upsymbol)=c("EntrezID","GENE.SYMBOL")
upSiganno<-merge(upSig,upsymbol,by="EntrezID",all=T,sort=F)
write.table(upSiganno,file="upSigProteinCoding-EntrezID-GeneSymbol.xls",sep="\t",quote=F,row.names=F)

downMatrix = rankMatrix(downList)
downAR = aggregateRanks(rmat=downMatrix)
colnames(downAR)=c("EntrezID","Pvalue")
downAdj=p.adjust(downAR$Pvalue,method="bonferroni")
downXls=cbind(downAR,adjPvalue=downAdj)
downFC=newTab[as.vector(downXls[,1]),]
downXls=cbind(downXls,logFC=rowMeans(downFC))
write.table(downXls,file="downProteinCoding-EntrezID.xls",sep="\t",quote=F,row.names=F)
downSig=downXls[(downXls$adjPvalue<padj & downXls$logFC< -logFC),]
down<-downSig$EntrezID
downsymbol<-select(org.Hs.eg.db, 
                   keys = down,
                   columns = c("ENTREZID","SYMBOL"),
                   keytype = "ENTREZID")
colnames(downsymbol)=c("EntrezID","GENE.SYMBOL")
downSiganno<-merge(downSig,downsymbol,by="EntrezID",all=T,sort=F)
write.table(downSiganno,file="downSigProteinCoding-EntrezID-GeneSymbol.xls",sep="\t",quote=F,row.names=F)

allSig = rbind(upSiganno,downSiganno)
write.table(allSig,file = 'allSignProteinCoding-EntrezID-GeneSymbol.xls',sep = '\t',quote = F)
