setwd("E:\\Simon\\Signature\\AAI_mouse_VCF\\MouseTumors")


library(MutationalPatterns)
library("BSgenome.Mmusculus.UCSC.mm10", character.only = TRUE)
ref_genome <- "BSgenome.Mmusculus.UCSC.mm10"

#load data
vcf_files = list.files( pattern = ".vcf", full.names = TRUE)
sample_names= list.files( pattern = ".vcf", full.names = F)

vcfs = read_vcfs_as_granges(vcf_files, sample_names, genome = ref_genome)
auto = extractSeqlevelsByGroup(species="Mus_musculus", 
                               style="UCSC",
                               group="all")
auto <- auto[1:21]
vcfs = lapply(vcfs, function(x) keepSeqlevels(x, auto))
#get 96 matrix
mut_mat <- mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)
colnames(mut_mat) <- c("M1_T1","M1_T2","M1_T3","M10_T1","M11_T1","M2_T1","M2_T2","M2_T3","M2_T4","M3_T1","M9_T1")

load("D:\\0530\\SignatureEstimation\\SignatureEstimation\\data\\signaturesCOSMIC.rda")
#将mut_mat 与 cosmic signature的行名排列一致
mut_mat <- cbind(mut_mat,rownames(mut_mat))
mut_mat <- mut_mat[order(mut_mat[,ncol(mut_mat)],decreasing=F),];
signaturesCOSMIC <- cbind(signaturesCOSMIC,rownames(signaturesCOSMIC))
signaturesCOSMIC <- signaturesCOSMIC [order(signaturesCOSMIC [,ncol(signaturesCOSMIC )],decreasing=F),];
mut_mat <- mut_mat[,-ncol(mut_mat)]
signaturesCOSMIC <- signaturesCOSMIC[,-ncol(signaturesCOSMIC)]
mut_mat <- matrix(as.numeric(mut_mat ),nrow=nrow(mut_mat ), dimnames=dimnames(mut_mat ))
signaturesCOSMIC<- matrix(as.numeric(signaturesCOSMIC),nrow=nrow(signaturesCOSMIC), dimnames=dimnames(signaturesCOSMIC))


sigsLiver = c(1,4,5,6,12,16,17,22,23,24)
# sigsLiver = c(1,2,3,5,6,10,12,13,15,16,17,20,23,26,22)
P = signaturesCOSMIC[, sigsLiver]
liver_signatures <- P

#library(RColorBrewer)
#col = brewer.pal(12,"Paired")
#col2 =brewer.pal(3,"Dark2")
#mycol <- c(col,col2)

liver_signatures <- cbind(liver_signatures[,-8],liver_signatures[,8])
colnames(liver_signatures)[10] <- "Signature 22 (AAI Signature)"

liver_signature_names <- colnames(liver_signatures)
mycol_original <- c("darkgreen","deepskyblue4","grey","orangered1","darkred","goldenrod1","deeppink4","royalblue4","darkolivegreen3","purple4")

mycol <- c("darkgreen","deepskyblue4","goldenrod1","darkred","orangered1","deeppink4","darkolivegreen3","purple4","grey","royalblue4")



# cosine similarity and plot
cos_sim_samples_signatures = cos_sim_matrix(mut_mat, as.matrix(signaturesCOSMIC[,22])) 
colnames(cos_sim_samples_signatures) <- "signature 22"


boxplot(cos_sim_samples_signatures,col="#008B8B",cex=1,pch=16,ylab="Cosine Similarity",xlab="Cosmic Signature 22")
text(1.2,0.66,"M3_T1",col="#EE0000",cex=0.6)
text(1.2,0.6,"M1_T1",col="#EE0000",cex=0.6)



m <- cbind(mut_mat, as.matrix(signaturesCOSMIC[,22]))
colnames(m)[ncol(m)] <- "Signature 22"
cos_sim_samples_signatures = cos_sim_matrix(m, m) 
library(pheatmap)
cols = colorRampPalette(c("#00558C","#F6F023"))(100)
pheatmap(cos_sim_samples_signatures,border_color="white",col=cols,cluster_rows=F,cluster_cols=F)

#strand bias analysis and plot
library("TxDb.Mmusculus.UCSC.mm10.knownGene")
genes_mm10 <- genes(TxDb.Mmusculus.UCSC.mm10.knownGene)

mut_mat_s <- mut_matrix_stranded(vcfs, ref_genome, genes_mm10)
strand_counts <- strand_occurrences(mut_mat_s, by=sample_names)
strand_bias <- strand_bias_test(strand_counts)

#输出，过滤只剩下T>A，命名为strand_bias_TA.txt，再读取
bias <- read.table("strand_bias_TA.txt",header=T)
count <- bias[,3]+bias[,4]
strand <- bias[,3:4]
colnames(strand) <- ("untranscribed","transcribed")


for (j in c(1:nrow(bias))){


	filename <- paste(bias[j,1],"strandbias.pdf",sep="")
	print(filename)


	pdf(filename,width=2,height=3)

	barplot(as.matrix(strand[j,]),ylim=c(0,count[j]),width=	0.5,xlim=c(0,5),space = 0.5,yaxt="n",xaxt="n",col="grey")

	axis(2,at=round(seq(0,count[j],length.out=2),digits=1),pos=0,las=	1,cex.axis=0.5)
	dev.off()
	j = j+1
}
------------------------------------------------------------------------------------------------------------------------------------------------

#plot the signature contriution

#calculate the contribution_counts and contribution_proportions 
fit_res <- fit_to_signatures(mut_mat,liver_signatures)
contribution_counts <- fit_res$contribution
CountsToProportions <- function(x){ 
    z <- matrix(0,nrow=nrow(x),ncol=ncol(x))
    i=1
    for (i in c(1:ncol(x))){
        z[,i] <- x [,i] / colSums(x)[i]
    } 
    dimnames(z) <- dimnames(x) 
    return(z)
}
contribution_proportions <- CountsToProportions(contribution_counts )


mutSign_nums = t(contribution_counts)
mutSign_props=t(contribution_proportions)
sig_cols = mycol

library(ggplot2)
library(reshape2)



  ordering <- order(colSums(t(mutSign_nums)),decreasing=T)
  mutSign_nums <- t(mutSign_nums)
  mutSign_props <- t(mutSign_props)
  mutSign_nums <- mutSign_nums[,ordering]
  mutSign_props <- mutSign_props[,ordering]
  sample.ordering <- colnames(mutSign_nums)
  x1 <- melt(mutSign_nums)
  x2 <- melt(mutSign_props)
  colnames(x1) <- c("Signature","Sample","Activity")
  colnames(x2) <- c("Signature","Sample","Activity")
  x1[,"class0"] <- c("Counts")
  x2[,"class0"] <- c("Proportions")
  df2 <- rbind(x1,x2)
  df2$class0 <- factor(df2$class0,c("Counts","Proportions"))
  df2$Sample <- factor(df2$Sample,sample.ordering)
  
scale<-1          
axis.text.y = element_text(hjust = 0.5,size=10*scale, family="mono")
axis.text = element_text(size = 10*scale, family = "mono")


p = ggplot(df2,aes(x=factor(Sample),y=Activity,fill=Signature))
p = p+geom_bar(stat="identity",position='stack')
p = p + scale_fill_manual(values=sig_cols)
p=p+theme_classic() 


p = p + theme(axis.text.x = element_text(angle=90,vjust=0.5,size=6*scale,face="bold",colour="black"))
p = p + theme(axis.text.y = element_text(size=10*scale,face="bold",colour="black"))


p = p + xlab("Samples") + ylab("Mutational Signature Content")
p = p + theme(axis.title.y = element_text(face="bold",colour="black",size=10*scale))  
p = p + theme(axis.title.x = element_text(face="bold",colour="black",size=10*scale))
 

p = p + ggtitle("Mutational Signature Exposures")
p = p + theme(plot.title=element_text(lineheight=1.0,face="bold",size=10*scale))


