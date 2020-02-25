setwd("C:\\Users\\think\\Desktop\\CNV\\mouse2")
chrom <- read.table("chromsize.txt",header=T)

rt <- read.table("M2.txt",header=T)

chrom_length <- c(0)

for (i in c(1:nrow(chrom))){
  
  chrom_length [i] <- sum(chrom[1:i,3])
  
}

new_rt <- matrix(NA,ncol=5, nrow=nrow(rt))

new_rt[which(rt[,1]==1),1] = rt[which(rt[,1]==1),1] 
new_rt[which(rt[,1]==1),2] = rt[which(rt[,1]==1),2]
new_rt[which(rt[,1]==1),3] = rt[which(rt[,1]==1),3]
new_rt[which(rt[,1]==1),4] = rt[which(rt[,1]==1),4]
new_rt[which(rt[,1]==1),5] = rt[which(rt[,1]==1),5]
lines = which(rt[,1]==1)

a <- lines [length(lines)]
b <- a+1

for (i in c(b:nrow(rt))){

 new_rt[i,1] <- rt[i,1]
 new_rt[i,2] <- rt[i,2]+  chrom_length[(rt[i,1]-1)]
 new_rt[i,3] <- rt[i,3]+  chrom_length[(rt[i,1]-1)]
 new_rt[i,4] <- rt[i,4]
 new_rt[i,5] <- rt[i,5]
 
}
colnames(new_rt) <- c("chom","start","end","CopyNumber","Sample")

d <- as.data.frame(new_rt)
library(ggplot2)
p = ggplot(d)+geom_point(aes(x=start,y=CopyNumber),cex=0.001,alpha=0)+scale_y_continuous(limits=c(0, 3), breaks=c(0,1,2,3))+facet_wrap(~Sample,nrow=3)+  theme(panel.grid =element_blank()) 

i =1
x= list(0)
for (i in c(1:nrow(new_rt))){
  print (i)
  a = d[i,2]
  b=d[i,4]
  c=d[i,3]
 
  x[[i]]=d
  x[[i]]$Sample = d[i,5] 


	if	(d[i,1] <=19){
  		if (b==2){
  			p 				=p+geom_segment(data=x[[i]],aes_string(x=a,y=b,xend=c,yend=b),col=	"#B22222",lwd=2)
		}
  		if (b<2){
  			p 		=p+geom_segment(data=x[[i]],aes_string(x=a,y=b,xend=c,yend=b),col=	"#20B2AA",lwd=2)
			}
  		if (b>2){p 		=p+geom_segment(data=x[[i]],aes_string(x=a,y=b,xend=c,yend=b),col=	"red",lwd=2)
			}
	}
	else{
  		if (b==1){
  			p 				=p+geom_segment(data=x[[i]],aes_string(x=a,y=b,xend=c,yend=b),col=	"#B22222",lwd=2)
		}
  		if (b<1){
  			p 		=p+geom_segment(data=x[[i]],aes_string(x=a,y=b,xend=c,yend=b),col=	"#20B2AA",lwd=2)
			}
  		if (b>1){p 		=p+geom_segment(data=x[[i]],aes_string(x=a,y=b,xend=c,yend=b),col=	"red",lwd=2)
			}
	}


}

#aes_string,相比aes是可以将变量转化为数字存入p。

x <- chrom_length[-21]
x <- c(1,x)
y <- (x+chrom_length)/2


p + scale_x_continuous(breaks=y, labels = paste("chr",c(1:19,"X","Y"),sep=""))+theme(axis.text.x = element_text(angle=90,vjust=0.5,size=10,face="bold",colour="black")) +theme(axis.text.y = element_text(angle=90,vjust=0.5,size=10,face="bold",colour="black"))+ theme(axis.text.y = element_text(angle=90,vjust=0.5,size=10,face="bold",colour="black")) 