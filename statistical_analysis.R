###################
# This is official R script to post with manuscript


###################
# 1. 16S rRNA sequencing show 1) different in crop system. 2) no different in aggregates
###################




###################
# 2. genomic potential, 1) no different in crop system
###################






###################
# 3. Figure 1: genomic potential of cellulose decomposition enzymes in FP among agg.
###################


###################
# read count table all pro+euk
kegg_all <-  as.data.frame(read.table("~/Box Sync/2018/1January/cobs/merged_kegg_cobs.tsv", header =T, row.names=1))
table_kegg_all <- t(kegg_all)
table <- table_kegg_all 

###################
# Read meta
meta=read.csv("~/Box Sync/2017/10October/cobs/sample_data_cobs.csv")


##############################
# remove poor quality, add R1 and R2
new_meta = data.frame()
fin_table = data.frame()
#uni = unique(meta$filename)
i=1
for ( one_sample in unique(meta$filename) ){
	#one_sample = uni[i]
	one_pair = paste0(gsub("_001","",one_sample),"_R1_001")
	if (meta[ meta$filename_pair == one_pair, ]$Note == "Poor Quality") next
	
	pair_to_remove = paste0(gsub("_001","",one_sample),"_R2_001")
	fin_table = rbind(fin_table, one_sample = rowSums( table[, c(one_pair, pair_to_remove)] ) )
	rownames(fin_table)[i] = toString(one_sample )
	
	new_meta = rbind(new_meta, meta[ meta$filename_pair == one_pair, ])
	i = i + 1
}
colnames(fin_table) <- rownames(table)
dim(new_meta)

dim(fin_table)
full_table <- merge(fin_table, new_meta, by.x = "row.names", by.y = "filename")
dim(full_table)
##############################
# Deseq2
sub_table <- subset(full_table, Crop == "PF")
Y = t(sub_table[, 2:(ncol(fin_table)+1)])
colnames(Y) <- sub_table$FileName.prefix
GRP = sub_table[,(ncol(fin_table)+2):ncol(sub_table)]


library(DESeq2)
library(ggplot2)
## Create a DESeqDataSet object for holding the counts and phenotype data.
deseq2_deseqds <- DESeqDataSetFromMatrix(countData = Y, colData = GRP, design = ~ SoilFrac)

## Differential expression analysis. 
deseq2_de <- DESeq(deseq2_deseqds)

# significantly different between medium and whole sooil
deseq2_de_out <- results(deseq2_de, contrast=c("SoilFrac","MM","WS"))

count <- counts(deseq2_de, normalized=TRUE)


####################
#  cellulose decomposition enzymes
# Figure 1

file_name = "~/Box Sync/2018/1January/cobs/figure1.pdf"
pdf(file_name, width=6, height=9)
one_ec = "ec.3.2.1.4"
d <- plotCounts(deseq2_de, gene = one_ec, intgroup = "SoilFrac", returnData=T)
p1 <- ggplot(d, aes(x=SoilFrac, y=count))+geom_boxplot()+scale_x_discrete(limits=c("WS", "LM","MM","SM","Micro"))+ggtitle("EC 3.2.1.4")+theme(axis.line=element_line(size=1.25), aspect.ratio=1, panel.border=element_blank(),axis.ticks=element_line(size=1.25, colour="black"), panel.background=element_blank(), axis.text=element_text(size=10, face="bold", colour="black"), axis.title=element_text(size=12, face="bold", colour="black"), plot.title = element_text(size=15, face="bold", hjust=0.5))+labs(x="Aggregates Fraction",y="Normalized Count")


for_aov <- data.frame(count = count[one_ec,], frac =GRP$SoilFrac)
fit <- aov(count ~ frac, data=for_aov)
summary(fit)
            # Df  Sum Sq Mean Sq F value   Pr(>F)    
# frac         4 6549154 1637289   21.42 4.64e-06 ***
# Residuals   15 1146379   76425  


one_ec = "ec.3.2.1.91"
d <- plotCounts(deseq2_de, gene = one_ec, intgroup = "SoilFrac", returnData=T)
p2 <- ggplot(d, aes(x=SoilFrac, y=count))+geom_boxplot()+scale_x_discrete(limits=c("WS", "LM","MM","SM","Micro"))+ggtitle("EC 3.2.1.91")+theme(axis.line=element_line(size=1.25), aspect.ratio=1, panel.border=element_blank(),axis.ticks=element_line(size=1.25, colour="black"), panel.background=element_blank(), axis.text=element_text(size=10, face="bold", colour="black"), axis.title=element_text(size=12, face="bold", colour="black"), plot.title = element_text(size=15, face="bold", hjust=0.5))+labs(x="Aggregates Fraction",y="Normalized Count")

for_aov <- data.frame(count = count[one_ec,], frac =GRP$SoilFrac)
fit <- aov(count ~ frac, data=for_aov)
summary(fit)
            # Df Sum Sq Mean Sq F value   Pr(>F)    
# frac         4 309191   77298   34.83 1.98e-07 ***
# Residuals   15  33286    2219 

one_ec = "ec.3.2.1.21"
d <- plotCounts(deseq2_de, gene = one_ec, intgroup = "SoilFrac", returnData=T)
p3 <- ggplot(d, aes(x=SoilFrac, y=count))+geom_boxplot()+scale_x_discrete(limits=c("WS", "LM","MM","SM","Micro"))+ggtitle("EC 3.2.1.21")+theme(axis.line=element_line(size=1.25), aspect.ratio=1, panel.border=element_blank(),axis.ticks=element_line(size=1.25, colour="black"), panel.background=element_blank(), axis.text=element_text(size=10, face="bold", colour="black"), axis.title=element_text(size=12, face="bold", colour="black"), plot.title = element_text(size=15, face="bold", hjust=0.5))+labs(x="Aggregates Fraction",y="Normalized Count")


for_aov <- data.frame(count = count[one_ec,], frac =GRP$SoilFrac)
fit <- aov(count ~ frac, data=for_aov)
summary(fit)
            # Df   Sum Sq  Mean Sq F value   Pr(>F)    
# frac         4 87601920 21900480   9.811 0.000418 ***
# Residuals   15 33483966  2232264 

one_ec = "ec.2.4.1.20"
d <- plotCounts(deseq2_de, gene = one_ec, intgroup = "SoilFrac", returnData=T)
p4 <- ggplot(d, aes(x=SoilFrac, y=count))+geom_boxplot()+scale_x_discrete(limits=c("WS", "LM","MM","SM","Micro"))+ggtitle("EC 2.4.1.20")+theme(axis.line=element_line(size=1.25), aspect.ratio=1, panel.border=element_blank(),axis.ticks=element_line(size=1.25, colour="black"), panel.background=element_blank(), axis.text=element_text(size=10, face="bold", colour="black"), axis.title=element_text(size=12, face="bold", colour="black"), plot.title = element_text(size=15, face="bold", hjust=0.5))+labs(x="Aggregates Fraction",y="Normalized Count")


for_aov <- data.frame(count = count[one_ec,], frac =GRP$SoilFrac)
fit <- aov(count ~ frac, data=for_aov)
summary(fit)
            # Df Sum Sq Mean Sq F value  Pr(>F)    
# frac         4 225011   56253   12.85 9.7e-05 ***
# Residuals   15  65645    4376 


one_ec = "ec.3.2.1.86"
d <- plotCounts(deseq2_de, gene = one_ec, intgroup = "SoilFrac", returnData=T)
p5 <- ggplot(d, aes(x=SoilFrac, y=count))+geom_boxplot()+scale_x_discrete(limits=c("WS", "LM","MM","SM","Micro"))+ggtitle("EC 3.2.1.86")+theme(axis.line=element_line(size=1.25), aspect.ratio=1, panel.border=element_blank(),axis.ticks=element_line(size=1.25, colour="black"), panel.background=element_blank(), axis.text=element_text(size=10, face="bold", colour="black"), axis.title=element_text(size=12, face="bold", colour="black"), plot.title = element_text(size=15, face="bold", hjust=0.5))+labs(x="Aggregates Fraction",y="Normalized Count")

for_aov <- data.frame(count = count[one_ec,], frac =GRP$SoilFrac)
fit <- aov(count ~ frac, data=for_aov)
summary(fit)
            # Df  Sum Sq Mean Sq F value  Pr(>F)   
# frac         4 1199699  299925   5.679 0.00547 **
# Residuals   15  792160   52811 


one_ec = "ec.2.7.1.205"
d <- plotCounts(deseq2_de, gene = one_ec, intgroup = "SoilFrac", returnData=T)
p6 <- ggplot(d, aes(x=SoilFrac, y=count))+geom_boxplot()+scale_x_discrete(limits=c("WS", "LM","MM","SM","Micro"))+ggtitle("EC 2.7.1.205")+theme(axis.line=element_line(size=1.25), aspect.ratio=1, panel.border=element_blank(),axis.ticks=element_line(size=1.25, colour="black"), panel.background=element_blank(), axis.text=element_text(size=10, face="bold", colour="black"), axis.title=element_text(size=12, face="bold", colour="black"), plot.title = element_text(size=15, face="bold", hjust=0.5))+labs(x="Aggregates Fraction",y="Normalized Count")

for_aov <- data.frame(count = count[one_ec,], frac =GRP$SoilFrac)
fit <- aov(count ~ frac, data=for_aov)
summary(fit)
            # Df Sum Sq Mean Sq F value Pr(>F)
# frac         4    3.0   0.750   0.079  0.988
# Residuals   15  142.3   9.484

library("gridExtra")
grid.arrange(p1,p2,p3,p4,p5,p6, nrow=3)
dev.off()


#########################################################
# 4. Heatmap
#########################################################


#heat map


# NMDS??



# start from here
deseq2_de_out <- results(deseq2_de, contrast=c("SoilFrac","MM","WS"))
summary(deseq2_de_out)

# out of 3546 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)     : 375, 11% 
# LFC < 0 (down)   : 342, 9.6% 


deseq2_de_out_sig <- subset(deseq2_de_out, padj < 0.01)
summary(deseq2_de_out_sig)

# out of 346 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)     : 184, 53% 
# LFC < 0 (down)   : 162, 47% 


#install.packages("gplots")
library("gplots")
# heatmap figure

rld <- rlog(deseq2_de)
mat <- assay(rld)
dim(mat)
head(mat)

#only up : High in MM
deseq2_de_out_sig_up <- subset(deseq2_de_out_sig, log2FoldChange > 0)
idx <- rownames(deseq2_de_out_sig_up)[which(deseq2_de_out_sig_up$log2FoldChange > -100)]
range(deseq2_de_out_sig_up$log2FoldChange)
length(idx)

pdf("~/Box Sync/2018/2February/cobs/writing/heat_map_high_in_MM.pdf", width=8, height=6)
heatmap.2(mat[idx,], trace = "none" , col = greenred(10), scale="row", margins=c(12,8))
dev.off()


#only down: High in WS
deseq2_de_out_sig_down <- subset(deseq2_de_out_sig, log2FoldChange < 0)
idx <- rownames(deseq2_de_out_sig_down)[which(deseq2_de_out_sig_down$log2FoldChange <  100)]
pdf("~/Box Sync/2018/2February/cobs/writing/heat_map_high_in_WS.pdf", width=8, height=6)
heatmap.2(mat[idx,], trace = "none" , col = greenred(10), scale="row", margins=c(12,8))
dev.off()


#########################################################
# 5. Co-occurrence
#########################################################

co_result <- read.table("~/Box Sync/2017/11November/cobs/result.tsv")
dim(co_result)
head(co_result)


interest <- c("ec.3.2.1.4", "ec.3.2.1.21","ec.3.2.1.86","ec.3.2.1.91","ec.2.4.1.20","ec.2.7.1.205")
sub1 <- subset(co_result, V2 %in% interest | V3 %in% interest)
dim(sub1)
head(sub1)

sub_all <- subset(sub1, V4 > 0.7 | V4 < -0.7)
length(sub_all)
# high in Medium
deseq2_de_out_sig_up <- subset(deseq2_de_out_sig, log2FoldChange > 0)
high_in_m <- rownames(deseq2_de_out_sig_up)[which(deseq2_de_out_sig_up$log2FoldChange > -100)]
length(high_in_m) #184


sub_up <- subset(sub1, V4 > 0.7 )
dim(sub_up) #high in medium
head(sub_up)
#how many overlap?
posi <- c(as.character(sub_up[,"V2"]), as.character(sub_up[,"V3"]))
length(posi)
uni_posi <- unique(posi)
length(uni_posi) #201

both_posi <- high_in_m[high_in_m %in% uni_posi]
length(both_posi) #59
final_posi <- both_posi[!(both_posi %in% interest) ]
length(final_posi) #54

write(final_posi,"~/Box Sync/2018/2February/cobs/posi.txt")

# high in WS
deseq2_de_out_sig_down <- subset(deseq2_de_out_sig, log2FoldChange < 0)
high_in_ws <- rownames(deseq2_de_out_sig_down)[which(deseq2_de_out_sig_down$log2FoldChange <  100)]
length(high_in_ws)


#how many overlap?
sub_down <- subset(sub1, V4 < -0.7 )
nega <- c(as.character(sub_down[,"V2"]), as.character(sub_down[,"V3"]))
length(nega)
uni_nega <- unique(nega)
length(uni_nega) #126

both_nega <- high_in_ws[high_in_ws %in% uni_nega]
length(both_nega) #42
final_nega <- both_nega[!(both_nega %in% interest) ]
length(final_nega) #42

write(final_nega,"~/Box Sync/2018/2February/cobs/nega.txt")










#########################################################
# 6. organism
#########################################################




interest <- c("3.2.1.4", "3.2.1.21","3.2.1.86","3.2.1.91","2.4.1.20","2.7.1.205")

pdf("~/Box Sync/2018/2February/cobs/writing/high_in_m.pdf")

for(j in 1:length(interest)) {
	ec = interest[5]
ec = interest[j]


path = paste0("~/Box Sync/2017/3March/cobs/data/",ec)
file_abun = paste0(path, ".fasta.summary_count.tsv")
abun=read.table(file_abun,header=T, row.names=1)
#dim(abun)
#abun[1:10,1:10]

#rowSums(abun)
#colSums(abun)

###################
# Read meta
meta=read.csv("~/Box Sync/2017/10October/cobs/sample_data_cobs.csv")

# merge R1 and R2 and remove poor qulaity
table <- abun
new_meta = data.frame()
fin_table = data.frame()
#uni = unique(meta$filename)
i=1
for ( one_sample in unique(meta$filename) ){
	#one_sample = uni[i]
	one_pair = paste0(gsub("_001","",one_sample),"_R1_001")
	if (meta[ meta$filename_pair == one_pair, ]$Note == "Poor Quality") next
	
	pair_to_remove = paste0(gsub("_001","",one_sample),"_R2_001")
	fin_table = rbind(fin_table, one_sample = rowSums( table[, c(one_pair, pair_to_remove)] ) )
	rownames(fin_table)[i] = toString(one_sample )
	
	new_meta = rbind(new_meta, meta[ meta$filename_pair == one_pair, ])
	i = i + 1
}
colnames(fin_table) <- rownames(table)
#dim(new_meta)
#dim(fin_table)
#fin_table[1:10,1:10]
#fin_table[1:10,1200:1214]
#rowSums(fin_table)
#only few bacteria are dominant
#summary(colSums(fin_table))
#hist(colSums(fin_table))

library(ggplot2)
########################
## sum plot
temp <- data.frame(name =  rownames(fin_table), count = rowSums(fin_table))
#ggplot(temp, aes(x=reorder(name, -count), y=count))+geom_bar(stat="identity")+theme(axis.text.x = element_text(angle=90))


temp$Crop = new_meta$Crop
temp$SoilFrac = new_meta$SoilFrac
temp$depth = new_meta$depth

pf <- subset(temp, Crop=="PF")
pf$nor <- pf$count / pf$depth
#ggplot(pf, aes(x=reorder(name, -nor), y=nor))+geom_bar(stat="identity")+theme(axis.text.x = element_text(angle=90))


########################
## plot for each phylum

library(phyloseq)
abun_table <- t(fin_table)
OTU = otu_table(abun_table, taxa_are_rows=TRUE)
row.names(new_meta) <- new_meta$filename
Samp  = sample_data(new_meta)
# read taxonomy
file_tax = paste0(path, ".fasta.summary_count.tsv.org.nb.org.list.tax")
taxmat = as.matrix(read.table(file_tax, sep='\t', row.names=1))
colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
TAX = tax_table(taxmat)

#make phyloseq object
physeq = merge_phyloseq(OTU, Samp, TAX)

#### This step normalize table by depth of sequence
#sample_data(physeq)$depth
#GPr = transform_sample_counts(physeq, function(x) x / sample_data(physeq)$depth )
temp_otu <- otu_table(physeq)

for (i in 1:length(sample_data(physeq)$depth)){
	temp_otu[,i]  <- temp_otu[,i] / sample_data(physeq)$depth[i]
}
otu_table(physeq) <- temp_otu

#sample_data(physeq)$Crop == "PF"
pf <- prune_samples(sample_data(physeq)$Crop == "PF", physeq)

#plot_bar(pf, fill="Phylum")

gl <- tax_glom(pf, "Phylum")
p<-plot_bar(gl, fill = "Phylum")
plot(p)


}

dev.off()
#### this is end of plot




