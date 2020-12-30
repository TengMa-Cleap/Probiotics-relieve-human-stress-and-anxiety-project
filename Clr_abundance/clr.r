library(microbiome)
library(phyloseq)
library(vegan)

# load data
###### mapping
raw_mp = read.table("mapping.txt", header=T, sep="\t",row.names = 1)
mp = raw_mp 
mp$Sample  =rownames(mp)

###### taxa and abundance table.
meta_ab=read.table("abundance.txt", header=T, row.names = 1, sep = "\t")
meta_tax=read.table("tax.txt", header=T, row.names=1, sep = "\t")

# subset data
meta_ab = subset(meta_ab, select = as.character(mp$Sample))
meta_ab = meta_ab[which(rowSums(meta_ab)!=0),] 
meta_tax = meta_tax[rownames(meta_tax) %in% rownames(meta_ab), ] 


###### Create a phyloseq input
meta_MP = sample_data(raw_mp)
meta_OTU = otu_table(as.matrix(meta_ab), taxa_are_rows = T)
meta_TAX = tax_table(as.matrix(meta_tax))
meta_physeq = phyloseq(meta_OTU, meta_TAX, meta_MP)

# Convert abundance table CLR format
meta_physeq_t = microbiome::transform(meta_physeq, 'clr')

dist.aitch <- vegdist(meta_physeq_t, method = "euclidean", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = FALSE)
Aitch <- as.matrix(dist.aitch)
write.table(Aitch, "Aitch.txt", quote = F, sep = "\t", row.names = F)

