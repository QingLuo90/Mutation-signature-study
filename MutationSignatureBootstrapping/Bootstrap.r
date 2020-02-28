###The scrit is to use bootstrappng method for 

setwd("C:\\Users\\think\\Desktop\\Mouse_otherCancer")

library(pracma)
library(MutationalPatterns)
library("BSgenome.Mmusculus.UCSC.mm10", character.only = TRUE)
ref_genome <- "BSgenome.Mmusculus.UCSC.mm10"

#load data
vcf_files = list.files( pattern = ".vcf", full.names = TRUE)
sample_names= list.files( pattern = ".vcf", full.names = F)

vcfs = read_vcfs_as_granges(vcf_files, sample_names, genome = ref_genome)
auto = extractSeqlevelsByGroup(species="Mus_musculus", 
                               style="UCSC",      group="all")
auto <- auto[1:21]
vcfs = lapply(vcfs, function(x) keepSeqlevels(x, auto))



#get 96 matrix
mut_mat <- mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)
#load the reference cosmic signatures
load("D:\\0530\\SignatureEstimation\\SignatureEstimation\\data\\signaturesCOSMIC.rda")

#make the 96 matrix ready for analysis
mut_mat <- cbind(mut_mat,rownames(mut_mat))
mut_mat <- mut_mat[order(mut_mat[,ncol(mut_mat)],decreasing=F),];
signaturesCOSMIC <- cbind(signaturesCOSMIC,rownames(signaturesCOSMIC))
signaturesCOSMIC <- signaturesCOSMIC [order(signaturesCOSMIC [,ncol(signaturesCOSMIC )],decreasing=F),];
mut_mat <- mut_mat[,-ncol(mut_mat)]
signaturesCOSMIC <- signaturesCOSMIC[,-ncol(signaturesCOSMIC)]
mut_mat <- matrix(as.numeric(mut_mat ),nrow=nrow(mut_mat ), dimnames=dimnames(mut_mat ))
signaturesCOSMIC<- matrix(as.numeric(signaturesCOSMIC),nrow=nrow(signaturesCOSMIC), dimnames=dimnames(signaturesCOSMIC))

#define signatures of liver cancers
sigsLiver = c(1,4,5,6,12,16,17,22,23,24)
P = signaturesCOSMIC[, sigsLiver]

#####define function for bootstrapping###################
bootstrapSigExposures <- function(m, P, R) {
  ## process and check function parameters
  ## m, P
  P = as.matrix(P)
  if(length(m) != nrow(P))
    stop("Length of vector 'm' and number of rows of matrix 'P' must be the same.")
  if(any(names(m) != rownames(P)))
    stop("Elements of vector 'm' and rows of matrix 'P' must have the same names (mutations types).")
  if(ncol(P) == 1)
    stop("Matrices 'P' must have at least 2 columns (signatures).")

  ## normalize m to be a vector of probabilities.
  count = sum(m)
  m = m / sum(m)

  ## find optimal solutions using using provided decomposition method for each bootstrap replicate
  ## matrix of signature exposures per replicate (column)
  K = length(m) #number of mutation types
  exposures = replicate(R, {
      mutations_sampled = sample(seq(K),count, replace = TRUE, prob = m)
      m_sampled = as.numeric(table(factor(mutations_sampled, levels=seq(K))))
      fitsignatures (m_sampled, P)
    })
  rownames(exposures) = colnames(P)
  colnames(exposures) = paste0('Replicate_',seq(R))

  ## compute estimation error for each replicate/trial (Frobenius norm)
  errors = apply(exposures, 2, function(e) FrobeniusNorm(m,P,e))
  names(errors) = colnames(exposures)

  return(list(exposures=exposures, errors=errors))
}
###########################################################################		 
		 

#set empty output matrices
exposure_confidence_interval <- matrix(nrow = ncol(mut_mat),ncol=2)
exposure_boot   <- matrix(nrow = ncol(mut_mat),ncol=1000)
erro_matrix <- matrix(nrow = ncol(mut_mat),ncol=1000)


#calculate bootstrap result for each sample
i=1
for (i in c(1:ncol(mut_mat))) {
    print(i)
    # tumor catelogue
    m = mut_mat[, i]
    ## bootstrap result
    boot = bootstrapSigExposures(m, P, R = 1000)
   
   exposure <- boot$exposures[8,]
   exposure <- sort(exposure) 
   k1 <- 1000*0.05/2
   k2 <- 1000*(1-0.05/2)
   exposure_confidence_interval[i,] <- c(exposure[k1],exposure[k2])
   exposure_boot[i,]  <- boot$exposures[8,]
   erro_matrix[i,] <- boot$errors
   i= i+1
}
rownames(exposure_confidence_interval) <- colnames(mut_mat)
rownames(erro_matrix) <- colnames(mut_mat)
rownames(exposure_boot) <- colnames(mut_mat)


#get output
exposure_confidence_interval <-cbind(exposure_confidence_interval,rowMeans(erro_matrix))

size <- colSums(mut_mat)

exposure_confidence_interval <- cbind(exposure_confidence_interval,size)

colnames(exposure_confidence_interval) <- c("CI_lower","CI_upper","meanErro","size")

write.table(exposure_confidence_interval,"exposure_confidence_interval.txt",row.names = T,quote=F,sep="\t")

write.table(erro_matrix,"erro_matrix.txt",row.names = T,quote=F,sep="\t")


#visualization of some results
plot(exposure_confidence_interval[,4],exposure_confidence_interval[,3],xlab="MutationCounts",ylab="MeanErro",cex=1,pch=16,col="black")

library(ggplot2)

exposure_confidence_interval <- as.data.frame(exposure_confidence_interval)
g <- ggplot(exposure_confidence_interval)
g+geom_point(aes(size,meanErro),pch=16,alpha=0.5,cex=3)

write.table(exposure_boot,"exposure_boot.txt",row.names = T,quote=F,sep="\t")
write.table(erro_matrix,"erro_matrix.txt",row.names = T,quote=F,sep="\t")


e <- t(erro_matrix)
boxplot(e[,1:10], horizontal=T,col="orange",cex=0.1,pch=16)

x <- rowMeans(erro_matrix)
boxplot(x, horizontal=T,col="orange",cex=0.1,pch=16)

