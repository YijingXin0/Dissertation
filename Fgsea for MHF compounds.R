library(fgsea)

drug<-read.csv("I:/logFC_matrix.csv",header=TRUE,colClasses=c("character","character",rep("numeric",159)))
drug<-drug[,-1]

gmt.file <- "I:/IAV_informative_genes.gmt"
IAVlist <- gmtPathways(gmt.file)
IAVlist$cellIAVdown<-IAVlist$cellIAVdown[1:40]

EStable_Statistics <-as.data.frame(c(1:17))
for (i in 2:160) {
  singledrug<-drug[,c(1,i)]
  ranks <- setNames(singledrug[,c(2)], singledrug[,c(1)])
  fgseaRes <- fgsea(IAVlist,ranks,eps=0.0,minSize=0, maxSize=500)
  es<-fgseaRes$ES[2]-fgseaRes$ES[1]
  nes<-fgseaRes$NES[2]-fgseaRes$NES[1]
  EStable<-as.data.frame(c(1:17))
  EStable[1,1]<-colnames(singledrug[2])
  EStable[2,1]<-fgseaRes$ES[2]
  EStable[3,1]<-fgseaRes$ES[1]
  EStable[4,1]<-es
  EStable[5,1]<-fgseaRes$NES[2]
  EStable[6,1]<-fgseaRes$NES[1]
  EStable[7,1]<-nes
  EStable[8,1]<-fgseaRes$pval[2]
  EStable[9,1]<-fgseaRes$pval[1]
  EStable[10,1]<-fgseaRes$padj[2]
  EStable[11,1]<-fgseaRes$padj[1]
  EStable[12,1]<-fgseaRes$log2err[2]
  EStable[13,1]<-fgseaRes$log2err[1]
  EStable[14,1]<-fgseaRes$size[2]
  EStable[15,1]<-fgseaRes$size[1]
  EStable[16,1]<-as.character(fgseaRes$leadingEdge[2])
  EStable[17,1]<-as.character(fgseaRes$leadingEdge[1])
  EStable_Statistics <- cbind(EStable_Statistics,EStable)
  print(i)
  i <- i+1
}
EStable_Statistics <- EStable_Statistics[,-1]
rownames(EStable_Statistics) <-c("drugID","upES","downES","ES","upNES","downNES","NES",
                                 "upPval","downPval","upPadj","downPadj","upLog2err",
                                 "downLog2err","upSize","downSize","upLeadingEdge","downLeadingEdge")
colnames(EStable_Statistics) <- colnames(drug)[2:160]
a<-t(EStable_Statistics)
write.csv(a,"I:/Result_of_fgsea.csv")
