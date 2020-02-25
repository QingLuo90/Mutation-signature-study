setwd("C:\\Users\\think\\Desktop\\CNV")
library(HMMcopy)

sample_names <- list.files( pattern = "_reads.wig", full.names = F)

i=1
sample_names_full <- list()
while (i <= length(sample_names)){
sample_names_full[[i]] <- (list.files( pattern = "_reads.wig", full.names = T)[i])
i <- i+1
}


uncorrected <- lapply (sample_names_full ,wigsToRangedData,gcfile="mm10_gc.wig",mapfile="mm10_map.wig" )

corrected_copy <- lapply(uncorrected, correctReadcount)

somatic_copy1 <- corrected_copy[[1]]
somatic_copy1$copy <- corrected_copy[[1]]$copy - corrected_copy[[4]]$copy


somatic_copy2 <- corrected_copy[[2]]
somatic_copy2$copy <- corrected_copy[[2]]$copy - corrected_copy[[5]]$copy

somatic_copy3 <- corrected_copy[[3]]
somatic_copy3$copy <- corrected_copy[[3]]$copy - corrected_copy[[6]]$copy

somatic_segments1 <- HMMsegment(somatic_copy1)

plotSegments(somatic_copy1, somatic_segments1, pch = ".", chr="chr1", ylab = "Tumour Copy Number",
 xlab = "Chromosome Position")
write.table(somatic_segments1$segs, "M14ALC11_segments.txt",sep="\t")



somatic_segments2 <- HMMsegment(somatic_copy2)

plotSegments(somatic_copy2, somatic_segments2, pch = ".", chr="chr1", ylab = "Tumour Copy Number",
 xlab = "Chromosome Position")
write.table(somatic_segments2$segs, "M14ALC43_segments.txt",sep="\t")


somatic_segments3 <- HMMsegment(somatic_copy3)

plotSegments(somatic_copy3, somatic_segments3, pch = ".", chr="chr1", ylab = "Tumour Copy Number",
 xlab = "Chromosome Position")
write.table(somatic_segments3$segs, "M14ALC7_segments.txt",sep="\t")

---------------------------------------------------------------------
segmented_copy  <- lapply(corrected_copy, HMMsegment)
M14ALC11_segmented <- segmented_copy[[1]]
M14ALC43_segmented <- segmented_copy[[2]]
M14ALC7_segmented <- segmented_copy[[3]]
M14ATA1_segmented <- segmented_copy[[4]]
M14ATA4_segmented <- segmented_copy[[5]]
M14ATA7_segmented <- segmented_copy[[6]]
