 mut_mat <- mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)
#²»Òª¸ÄË³Ðò


j=1
for (j in c(1:ncol(mut_mat))){


	filename <- paste(colnames(mut_mat)[j],"pdf",sep=".")
	print(filename)
	freq <- mut_mat[,j]/sum(mut_mat[,j])

	pdf(filename,width=16,height=5)

	bases <- c("A","C","G","T")
  	ctxt16 <- paste(rep(bases,each=4),rep(bases,4),sep=".")
  	mt <- c("CA","CG","CT","TA","TC","TG")
  	types96 <- paste(rep(mt,each=16),rep(ctxt16,6),sep="_")
  	types96 <- sapply(types96,function(z){sub("\\.",substr(z,1,1),z)})
  	context <- substr(types96,4,6)


	col96_original <- c	(rep("#0000EE",16),rep("black",16),rep("#EE0000",16),rep("grey",16),rep("#90EE90",16),rep("#FFE4E1",16))


	col96 <- c	(rep("#0000EE",16),rep("black",16),rep("#EE0000",16),rep("#999999",16),rep("#A1CE63",16),rep("#EBC6C4",16))

	labs <-c(rep("C>A",16),rep("C>G",16),rep("C>T",16),rep("T>	A",16),rep("T>C",16),rep("T>G",16))



	bp <- barplot(freq,col=col96,border=col96,space=	1.5,yaxt="n",xaxt="n",ylim=c(0,0.22),xlim=c(0,50),width=0.2,las=2)

	axis(2,at=round(seq(0,0.2,length.out=3),digits=1),pos=0,las=	1,cex.axis=1.5)

 	axis(1,at=bp,labels=context,pos=0,las=2,cex.axis=1,tick=F,cex.axis=1,font=0.7);

	for(i in seq(1,81,by=16)){
    		rect(bp[i],par()$usr[4],bp[i+15],par()$usr[4]-0.05*diff(par()$usr[3:4]),col=col96[i],border=col96[i]);
    		text((bp[i]+bp[i+15])/2,par()$usr[4]+0.09*diff(par()$usr[3:4]),labels=labs[i], xpd= TRUE, cex = 1.5)
	}

dev.off()
j  = j+1
}



