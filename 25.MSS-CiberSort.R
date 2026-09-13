#install.packages("e1071")
#install.packages("parallel")
#install.packages("preprocessCore")
#BiocManager::install("preprocessCore")
#install.packages("devtools")
#install.packages("corrplot")
#install.packages("vioplot")
#devtools::install_github('shenorrlab/bseqsc')
#devtools::install_github('Moonerss/CIBERSORT')

library(devtools)
library(e1071)
library(preprocessCore)
library(parallel)
#library(bseqsc)
library(ggplot2)
library(CIBERSORT)
library(corrplot)
library(vioplot)

setwd("D:/workspace/TCGA/UCEC/input")  

data(LM22)

data = read.table('TCGA_UCEC_TPM.txt',header=T, sep="\t",check.names=F,row.names = 1)
dimnames=list(rownames(data),colnames(data))
data=matrix(as.numeric(as.matrix(data)), nrow=nrow(data), dimnames = dimnames)

setwd("D:/workspace/TCGA/UCEC/output/MSS-2/ML/ML-63_0.0005")
rs = read.table('final_mode_riskscore.txt',header=T, sep="\t",check.names=F,row.names = 1)
rs$riskscore = as.numeric(rs$riskscore)
rs$Type=ifelse(rs[,"riskscore"]>median(rs[,"riskscore"]), "High", "Low")


result <- cibersort(sig_matrix = LM22, mixture_file=data)
results= as.matrix(result[,1:(ncol(result)-3)])
results=rbind(id=colnames(results),results)

setwd("D:/workspace/TCGA/UCEC/output/MSS-2/TME") 

write.table(results,file="IRG_CIBERSORT-Results.txt",sep="\t",quote=F,col.names=F)

# Analysis step; see the surrounding code for details.
data=t(data)

rt=read.table('IRG_CIBERSORT-Results.txt',header=T, sep="\t",check.names=F,row.names = 1)
dimnames=list(rownames(rt),colnames(rt))
rt=matrix(as.numeric(as.matrix(rt)), nrow=nrow(rt), dimnames = dimnames)

rownames(rt)=substr(rownames(rt),1,12)
rownames(rt)=gsub('[.]', '-', rownames(rt))

# Analysis step; see the surrounding code for details.
sameSample=intersect(row.names(rs),row.names(rt))

#
rs=rs[sameSample,"Type",drop=F]
rt=rt[sameSample,,drop=F]
rt=cbind(rt, rs)
#view(rt)
rt$Type=factor(rt$Type, levels=c("Low", "High"))

rt1 = subset(rt,rt$Type=='High')
rt2 = subset(rt,rt$Type=='Low')

highNum=nrow(rt1)
lowNum=nrow(rt2)
Type=c(rep(1,highNum),rep(2,lowNum))


rt=rbind(rt1,rt2)
rt <- subset(rt, select = -c(Type))
#view(rt)
outTab=data.frame()
pdf(file="IRG-CiberSort-Vioplot.pdf",width=13,height=8)
par(las=1,mar=c(10,6,3,3))
x=c(1:ncol(rt))
y=c(1:ncol(rt))

plot(x,y,
     xlim=c(0,63),ylim=c(min(rt),max(rt)+0.05),
     main="",xlab="",ylab="Fraction",
     pch=21,
     col="white",
     xaxt="n")

showcolname <- c()
for(i in 1:ncol(rt)){
  if(sd(rt[1:highNum,i]==0)){
    rt[1,i]=0.00001    
  }
  if(sd(rt[(highNum+1):(highNum+lowNum),i])==0){
    rt[(highNum+1),i]=0.00001
  }
  rt1=rt[1:highNum,i]
  rt2=rt[(highNum+1):(highNum+lowNum),i]
  # Plot the corresponding figure.
  vioplot(rt1,at=3*(i-1),lty=1,add=T,col='blue')
  vioplot(rt2,at=3*(i-1)+1,lty=1,add=T,col='red')
  # Analysis step; see the surrounding code for details.
  wilcoxTest=wilcox.test(rt1,rt2)
  p=wilcoxTest$p.value
  if(p<0.05){
    cellPvalue=cbind(Cell=colnames(rt)[i],pvalue=p)
    outTab=rbind(outTab,cellPvalue)
    #showcolname <- c(showcolname,colnames(rt)[i])
  }
  mx=max(c(rt1,rt2))
  lines(c(x=3*(i-1)+0.2,x=3*(i-1)+0.8),c(mx,mx))
  text(x=3*(i-1)+0.5,y=mx+0.02,labels=ifelse(p<0.001,paste0("p<0.001"),paste0("p=",sprintf("%03f",p))),
       cex=0.8)
}
legend("topright",c("High","Low"),lwd=3,bty="n",cex=1,col=c('blue','red'))
text(seq(1,64,3),-0.04,xpd=NA,labels=colnames(rt),cex=1,srt=45,pos=2)
dev.off()

write.table(outTab,file="IRG-CiberSort-Diff.xls",sep="\t",row.names = F,quote=F)


