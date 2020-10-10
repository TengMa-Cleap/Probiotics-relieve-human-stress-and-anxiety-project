# This script is used to calculate the abundance of sgb
setwd("/userdata1/data_mat/meta/malaysia_bin2_20191228/0.6_all_SGB_bbmap")
library(reshape2)
library(doBy)
library(psych)
library(gplots)
library(RColorBrewer)

dd=read.table("pileup.combined.f", sep = "\t", header = T)
sr=read.table("Sample_reads", header = T, sep = "\t")

dd$Total_SGB_reads=dd$Plus_reads + dd$Minus_reads
dd_sum=summaryBy(list(c("Length", "Covered_bases", "Total_SGB_reads"), c("SGB_ID", "Sample")), data=dd, FUN=sum)
dd_sum$cov=dd_sum$Covered_bases.sum/dd_sum$Length.sum*100
dd_sum=merge(dd_sum, sr, by="Sample")
dd_sum$rpkm=dd_sum$Total_SGB_reads.sum/dd_sum$Length.sum/dd_sum$num_seqs*10^9

reset_reads=''
for(i in 1:nrow(dd_sum)){
  if(dd_sum$cov[i]>=50){
    reset_reads[i]=dd_sum$rpkm[i]
  }else{
    reset_reads[i]=0
  }
}

dd_sum$reset_rpkm=as.numeric(reset_reads)

dd_matrix=dcast(dd_sum, Sample~SGB_ID, value.var = "reset_rpkm")
dd_matrix_t = t(dd_matrix)
write.table(dd_matrix, "rpkm_matrix.txt", quote = F, sep = "\t", row.names = F)

# For abundance 
dd_matrix_ab=melt(dd_matrix)
mp=read.table("mapping", sep = "\t", header = T)
dd_matrix_ab=merge(dd_matrix_ab, mp, by="Sample")

dd_matrix_ab_mean=summaryBy(list("value", c("Time", "variable")), data=dd_matrix_ab, FUN=mean)

write.table(dd_matrix_ab_mean, "dd_matrix_ab_mean.txt", quote = F, sep = "\t", row.names = F)

##### PCoA
library(vegan)
library(ggplot2)
library(ggpubr)

get_plot_data<-function(dist_data){
  plot.list=list()
  dist.pcoa=cmdscale(dist_data, eig=TRUE)
  pc12=dist.pcoa$points[,1:2]
  pc_importance=round(dist.pcoa$eig/sum(dist.pcoa$eig)*100,digits = 2)
  pc12=as.data.frame(pc12)
  pc12[,3]=row.names(pc12)
  x.label=paste("PCoA 1 (", pc_importance[1],digits=4,"%)", sep="")
  y.label=paste("PCoA 2 (", pc_importance[2],digits=4,"%)", sep="")
  plot.list$pc12<-pc12
  plot.list$x.label<-x.label
  plot.list$y.label<-y.label
  return(plot.list)
}

dd.matrix.m=dd_matrix[,2:ncol(dd_matrix)]
row.names(dd.matrix.m)=dd_matrix$Sample

# bray_curtis distance
dist.bray <- vegdist(dd.matrix.m, method = "bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = FALSE) 
bray <- as.matrix(dist.bray)
write.table(bray, "bray.txt", quote = F, sep = "\t", row.names = F)

# jaccard distance
dist.bray <- vegdist(dd.matrix.m, method = "jaccard", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = FALSE) 
jac <- as.matrix(dist.bray)
write.table(jac, "jaccard.txt", quote = F, sep = "\t", row.names = F)


dist.bray.pc=get_plot_data(dist.bray)
dist.bray.pc.plot<-merge(dist.bray.pc$pc12, mp, by.x="V3", by.y="Sample")
pcoa = cmdscale(bray, k=4, eig =T)
points = as.data.frame(pcoa$points)
eig = pcoa$eig
V1 = eig[1]/sum(eig)
V2 = eig[2]/sum(eig)

write.table(V1, "V1.txt", quote = F, sep = "\t", row.names = F)
write.table(V2, "V2.txt", quote = F, sep = "\t", row.names = F)
bray.plot=ggscatter(dist.bray.pc.plot, x = "V1", y = "V2", color = "Time", shape = "Group", ellipse = T, ggtheme = theme_bw())
bray.plot
ggsave(bray.plot, file="bray_pcoa.pdf", width = 7, height = 5)
