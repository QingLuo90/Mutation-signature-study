#!/usr/bin/env Rscript

library(sciClone)
setwd("E:/BDMi/sciClone ")

v0 = read.table("VAF.txt", sep="\t")
cn0 = read.table("CNV.txt")
clusterParams="empty"
sc = sciClone(vafs=list(v0), sampleNames=c("Demo"), copyNumberCalls=list(cn0), minimumDepth=100, doClustering=TRUE, clusterParams=clusterParams,maximumClusters=10, copyNumberMargins=0.25, useSexChrs=FALSE)


writeClusterTable(sc, "Demo_clusters")
writeClusterSummaryTable(sc, "Demo_cluster.summary")

sc.plot1d(sc,"Demo.pdf", highlightSexChrs=TRUE, highlightsHaveNames=FALSE, overlayClusters=TRUE, showTitle=TRUE, cnToPlot=c(1:3))
