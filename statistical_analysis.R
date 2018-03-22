###################
# This is official R script to post with manuscript
rm(list=ls())
##########   Load library   ############
library(ape)
library(phyloseq)
library(ggplot2)
library(vegan)

library(vegan)
library(ggplot2)
library(reshape2)
library(corrplot)

library(DESeq2)

##########   Load function   ############


pairwise.adonis <- function(x,factors, sim.function = 'vegdist', sim.method = 'bray', p.adjust.m ='bonferroni')
{

co = combn(unique(as.character(factors)),2)
pairs = c()
F.Model =c()
R2 = c()
p.value = c()

for(elem in 1:ncol(co)){
if(sim.function == 'daisy'){
library(cluster); x1 = daisy(x[factors %in% c(co[1,elem],co[2,elem]),],metric=sim.method)
} else{x1 = vegdist(x[factors %in% c(co[1,elem],co[2,elem]),],method=sim.method)}

ad = adonis(x1 ~ factors[factors %in% c(co[1,elem],co[2,elem])] , permutations=9999);
pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
F.Model =c(F.Model,ad$aov.tab[1,4]);
R2 = c(R2,ad$aov.tab[1,5]);
p.value = c(p.value,ad$aov.tab[1,6])
}
p.adjusted = p.adjust(p.value,method=p.adjust.m)
sig = c(rep('',length(p.adjusted)))
sig[p.adjusted <= 0.05] <-'.'
sig[p.adjusted <= 0.01] <-'*'
sig[p.adjusted <= 0.001] <-'**'
sig[p.adjusted <= 0.0001] <-'***'

pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted,sig)
print("Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1")
return(pairw.res)
} 


###################
# Read 16S rRNA sequencing data make Phyloseq object
###################
biom_file = "~/Box Sync/2017/7July/cobs/otu_table_mc2_w_tax.json.biom"
biom_otu_tax <- import_biom(biom_file)
tree = read.tree("~/Box Sync/2017/7July/cobs/rep_set.tre")
meta  <- read.csv("~/Box Sync/2017/7July/cobs/meta.csv",header=T, row.names=1)
sampledata = sample_data(meta)
physeq = merge_phyloseq(sampledata, biom_otu_tax, tree)
only_july <- prune_samples(as.data.frame(sample_data(physeq))$month_sample == "July", physeq )


remove_corn_m <- prune_samples(!(as.data.frame(sample_data(only_july))$crop_system == "corn" & as.data.frame(sample_data(only_july))$agg_frac == "MM"), only_july)
remove_corn_s <- prune_samples(!(as.data.frame(sample_data(remove_corn_m))$crop_system == "corn" & as.data.frame(sample_data(remove_corn_m))$agg_frac == "SM"),remove_corn_m)

remove_p_m <- prune_samples(!(as.data.frame(sample_data(remove_corn_s ))$crop_system == "prairie" & as.data.frame(sample_data(remove_corn_s ))$agg_frac == "MM"),remove_corn_s )
remove_p_s <- prune_samples(!(as.data.frame(sample_data(remove_p_m))$crop_system == "prairie" & as.data.frame(sample_data(remove_p_m))$agg_frac == "SM"),remove_p_m)

only44sample <- remove_p_s
set.seed(1)
#rare <- rarefy_even_depth(only44sample)
#saveRDS(rare, "~/Box Sync/2018/2February/cobs/cobs_16s_44sample_rarefied.rds")
rare <- readRDS("~/Box Sync/2018/2February/cobs/cobs_16s_44sample_rarefied.rds")
# # this is all 120 files
# #rare <- readRDS("~/Box Sync/2017/7July/cobs/cobs_16s_rarefied.rds")

# glom <- tax_glom(rare, "Rank2")
# plot_bar(glom)

# sample_sums(glom)#13956
# #1081815                        "k__Bacteria" "p__Actinobacteria"   NA


# mean(otu_table(glom)[c("1081815"),]) / 13956
# taxa_sums(glom)




###################
# Supplementary figure 1
###################
rm(list=ls())
rare <- readRDS("~/Box Sync/2018/2February/cobs/cobs_16s_44sample_rarefied.rds")


#NMDS
data.selected <- t(otu_table(rare))
mds.all=metaMDS(data.selected,distance="bray", k=2)
data.sm=as.data.frame(scores(mds.all))
data.sm$crop <- as.factor(sample_data(rare)$crop_system)
data.sm$for_label <- factor(data.sm$crop, labels= c("C", "P", "FP"))

pdf("~/Box Sync/2018/2February/cobs/writing/supplementary_figure1.pdf")
ggplot(data=data.sm, aes(x=NMDS1, y=NMDS2, color=for_label)) + geom_point(size=5) +
 geom_hline(yintercept=0.0, colour="grey", lty=2)+
 geom_vline(xintercept=0.0, colour="grey",lty=2) +
 theme_bw() + scale_color_brewer(palette="Set1") + theme(legend.position = c(0.11, 0.9)) +
 theme(legend.background = element_rect(fill="white", size=0.3, linetype="solid", colour="black"))+stat_ellipse()+labs(color="Crop System")
dev.off()


pairwise.adonis(data.selected, as.factor(sample_data(rare)$crop_system))
# [1] "Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1"
                   # pairs  F.Model         R2 p.value p.adjusted sig
# 1        prairie vs corn 3.976695 0.15308704  0.0004     0.0012   *
# 2 prairie vs prairieFert 1.996555 0.06239905  0.0173     0.0519    
# 3    corn vs prairieFert 5.525448 0.15553492  0.0001     0.0003  **


# #only WS
# only_ws <- prune_samples( as.data.frame(sample_data(rare))$agg_frac == "WS", rare)

# data.selected <- t(otu_table(only_ws))
# mds.all=metaMDS(data.selected,distance="bray", k=2)
# data.sm=as.data.frame(scores(mds.all))
# data.sm$crop <- as.factor(sample_data(only_ws)$crop_system)
# ggplot(data=data.sm, aes(x=NMDS1, y=NMDS2, color=crop)) + geom_point(size=5) +
 # geom_hline(yintercept=0.0, colour="grey", lty=2)+
 # geom_vline(xintercept=0.0, colour="grey",lty=2) +
 # theme_bw() + scale_color_brewer(palette="Set1") + theme(legend.position = c(0.08, 0.82)) +
 # theme(legend.background = element_rect(fill="white", size=0.3, linetype="solid", colour="black"))

# pairwise.adonis(data.selected, as.factor(sample_data(only_ws)$crop_system))
# #only ws
# # [1] "Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1"
                   # # pairs   F.Model        R2 p.value p.adjusted sig
# # 1        prairie vs corn 1.7470238 0.2255090   0.047      0.141    
# # 2 prairie vs prairieFert 0.9227838 0.1332966   0.415      1.000    
# # 3    corn vs prairieFert 1.9232605 0.2427360   0.064      0.192 ''



###################
# Supplementary figure 3A -> Supplementary Figure 5
###################
rm(list=ls())
rare <- readRDS("~/Box Sync/2018/2February/cobs/cobs_16s_44sample_rarefied.rds")
#only PF compare agg
only_pf <- prune_samples( as.data.frame(sample_data(rare))$crop_system == "prairieFert", rare)

data.selected <- t(otu_table(only_pf))
mds.all=metaMDS(data.selected,distance="bray", k=2)
data.sm=as.data.frame(scores(mds.all))
data.sm$agg <- as.factor(sample_data(only_pf)$agg_frac)
data.sm$for_label <- factor(data.sm$agg, labels= c("LM","Micro","MM","SM","WS"))


pdf("~/Box Sync/2018/2February/cobs/writing/supplementary_figure5.pdf")
ggplot(data=data.sm, aes(x=NMDS1, y=NMDS2, color=for_label)) + geom_point(size=5) +
 geom_hline(yintercept=0.0, colour="grey", lty=2)+
 geom_vline(xintercept=0.0, colour="grey",lty=2) +
 theme_bw() + scale_color_brewer(palette="Set1") + theme(legend.position = c(0.9, 0.22)) +
 theme(legend.background = element_rect(fill="white", size=0.3, linetype="solid", colour="black"))+labs(color="Aggregates")
dev.off()



pairwise.adonis(data.selected, as.factor(sample_data(only_pf)$agg))
# [1] "Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1"
         # pairs   F.Model        R2 p.value p.adjusted sig
# 1     WS vs SM 0.7705684 0.1138115  0.7681          1    
# 2     WS vs LM 0.7415927 0.1100026  0.6889          1    
# 3     WS vs MM 0.6975571 0.1041510  0.9168          1    
# 4  WS vs micro 0.9351090 0.1348370  0.5095          1    
# 5     SM vs LM 1.1312410 0.1586317  0.2848          1    
# 6     SM vs MM 0.9426523 0.1357770  0.5141          1    
# 7  SM vs micro 0.7842446 0.1155979  0.7992          1    
# 8     LM vs MM 0.8649800 0.1259989  0.6296          1    
# 9  LM vs micro 1.3383675 0.1823795  0.1400          1    
# 10 MM vs micro 1.1421073 0.1599118  0.3406          1 


# #all compare agg
# data.selected <- t(otu_table(rare))
# mds.all=metaMDS(data.selected,distance="bray", k=2)
# data.sm=as.data.frame(scores(mds.all))
# data.sm$agg <- as.factor(sample_data(rare)$agg_frac)
# ggplot(data=data.sm, aes(x=NMDS1, y=NMDS2, color=agg)) + geom_point(size=5) +
 # geom_hline(yintercept=0.0, colour="grey", lty=2)+
 # geom_vline(xintercept=0.0, colour="grey",lty=2) +
 # theme_bw() + scale_color_brewer(palette="Set1") + theme(legend.position = c(0.08, 0.82)) +
 # theme(legend.background = element_rect(fill="white", size=0.3, linetype="solid", colour="black"))

# pairwise.adonis(data.selected, as.factor(sample_data(rare)$agg))
# # [1] "Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1"
         # # pairs   F.Model         R2 p.value p.adjusted sig
# # 1  micro vs WS 1.1533115 0.04981195   0.208       1.00    
# # 2  micro vs LM 2.2920896 0.09435539   0.011       0.11    
# # 3  micro vs SM 1.0423117 0.06929199   0.309       1.00    
# # 4  micro vs MM 1.1609427 0.07657457   0.230       1.00    
# # 5     WS vs LM 0.9230130 0.04026578   0.434       1.00    
# # 6     WS vs SM 1.1664998 0.07691292   0.215       1.00    
# # 7     WS vs MM 0.9156073 0.06138586   0.420       1.00    
# # 8     LM vs SM 1.8013795 0.11400141   0.042       0.42    
# # 9     LM vs MM 1.4393119 0.09322384   0.104       1.00    
# # 10    SM vs MM 0.9485838 0.13651469   0.525       1.00   

###################
# 2. genomic potential, 1) no different in crop system
###################




### MNDS for all kegg enzyme


data.selected=read.csv("input_all.csv",header=T)
head(data.selected)
dim(data.selected)

sub_data <- subset(data.selected, Sample != "H5")
dim(sub_data)
data.selected <- sub_data

## NMDS

mds.all=metaMDS(data.selected[,2:2726],distance="bray", k=2)

data.sm=as.data.frame(scores(mds.all))

data.sm$Aggregates=as.factor(data.selected[,2727])

data.sm

levels(data.sm$Aggregates)
data.sm$Aggregates <- factor(data.sm$Aggregates, levels=c("Whole Soil", "Large", "Medium", "Small", "Micro"))


ggplot(data=data.sm, aes(x=NMDS1, y=NMDS2, color=Aggregates)) + geom_point(size=5) + xlim(-0.12, 0.12) + ylim(-0.07, 0.07)+
 geom_hline(yintercept=0.0, colour="grey", lty=2)+
 geom_vline(xintercept=0.0, colour="grey",lty=2) +
 theme_bw() + scale_color_brewer(palette="Set1") + theme(legend.position = c(0.08, 0.82)) +
 theme(legend.background = element_rect(fill="white", size=0.3, linetype="solid", colour="black"))





## MRPP & ANOSIM

trt=as.factor(data.selected[,2727])
enzymes=as.matrix(data.selected[,2:2726])

mrpp(enzymes, group=trt, distance="bray")
anosim(enzymes, grouping=trt, permutations=999, distance="bray")


###################
# Supplementary figure 2
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

# sub_table <- subset(full_table, Crop == "PF")
# sub_table <- full_table
# Y = t(sub_table[, 2:(ncol(fin_table)+1)])
# colnames(Y) <- sub_table$FileName.prefix
# GRP = sub_table[,(ncol(fin_table)+2):ncol(sub_table)]


# data.selected <- t(Y)
# dim(data.selected)

# #sub_data <- subset(data.selected, rownames(data.selected) != "H5")
# #sub_meta <- subset(GRP, FileName.prefix != "H5")
# dim(sub_data)
# #data.selected <- sub_data

# mds.all=metaMDS(data.selected,distance="bray", k=2)
# data.sm=as.data.frame(scores(mds.all))
# data.sm$Aggregates = as.factor(sub_meta$SoilFrac)
# data.sm$Crop = as.factor(GRP$Crop)


# ggplot(data=data.sm, aes(x=NMDS1, y=NMDS2, color=Crop)) + geom_point(size=5) + xlim(-0.12, 0.12) + ylim(-0.07, 0.07)+
 # geom_hline(yintercept=0.0, colour="grey", lty=2)+
 # geom_vline(xintercept=0.0, colour="grey",lty=2) +
 # theme_bw() + scale_color_brewer(palette="Set1") + theme(legend.position = c(0.08, 0.82)) +
 # theme(legend.background = element_rect(fill="white", size=0.3, linetype="solid", colour="black"))

# pairwise.adonis(data.selected, as.factor(GRP$Crop))

# interest <- c("ec.3.2.1.4", "ec.3.2.1.21","ec.3.2.1.86","ec.3.2.1.91","ec.2.4.1.20","ec.2.7.1.205")
# carbon <- data.selected[,interest]
# data.selected <- carbon
# mds.all=metaMDS(data.selected,distance="bray", k=2)
# data.sm=as.data.frame(scores(mds.all))
# data.sm$Crop = as.factor(GRP$Crop)

#normalize?
#sub_table <- subset(full_table, Crop == "PF")
#sub_table <- subset(full_table, SoilFrac == "WS")
sub_table <- full_table
Y = t(sub_table[, 2:(ncol(fin_table)+1)])
colnames(Y) <- sub_table$FileName.prefix
GRP = sub_table[,(ncol(fin_table)+2):ncol(sub_table)]


## Create a DESeqDataSet object for holding the counts and phenotype data.
#deseq2_deseqds <- DESeqDataSetFromMatrix(countData = Y, colData = GRP, design = ~ SoilFrac)
deseq2_deseqds <- DESeqDataSetFromMatrix(countData = Y, colData = GRP, design = ~ Crop)
## Differential expression analysis. 
deseq2_de <- DESeq(deseq2_deseqds)

# significantly different between medium and whole sooil
#deseq2_de_out <- results(deseq2_de, contrast=c("SoilFrac","MM","WS"))

count <- counts(deseq2_de, normalized=TRUE)
data.selected <- t(count)



sub_data <- subset(data.selected, rownames(data.selected) != "H5")
sub_meta <- subset(GRP, FileName.prefix != "H5")
data.selected <- sub_data
GRP <- sub_meta


mds.all=metaMDS(data.selected,distance="bray", k=2)
data.sm=as.data.frame(scores(mds.all))
data.sm$Crop = as.factor(GRP$Crop)
data.sm$for_label <- factor(data.sm$Crop, labels= c("C", "P", "FP"))

pdf("~/Box Sync/2018/2February/cobs/writing/supplementary_figure2.pdf")
ggplot(data=data.sm, aes(x=NMDS1, y=NMDS2, color=for_label)) + geom_point(size=5) + 
 geom_hline(yintercept=0.0, colour="grey", lty=2)+
 geom_vline(xintercept=0.0, colour="grey",lty=2) +
 theme_bw() + scale_color_brewer(palette="Set1") + theme(legend.position = c(0.11, 0.9)) +
 theme(legend.background = element_rect(fill="white", size=0.3, linetype="solid", colour="black"))+labs(color="Crop System")
 dev.off()
 
#+stat_ellipse()
pairwise.adonis(data.selected, as.factor(GRP$Crop))
#only whole soil, compare crop
     # pairs  F.Model        R2 p.value p.adjusted sig
# 1  CC vs P 2.148548 0.2636725   0.061      0.183    
# 2 CC vs PF 3.942315 0.3965188   0.031      0.093    
# 3  P vs PF 2.257778 0.2734123   0.066      0.198 

#all
#with H5 <- this is what i use
# [1] "Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1"
     # pairs   F.Model         R2 p.value p.adjusted sig
# 1 PF vs CC 4.5064336 0.13059691  0.0007     0.0021   *
# 2  PF vs P 0.9813639 0.03167594  0.4004     1.0000    
# 3  CC vs P 2.9264497 0.11740339  0.0098     0.0294   .

#without H3
# [1] "Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1"
     # pairs   F.Model        R2 p.value p.adjusted sig
# 1 PF vs CC 4.8254460 0.1426573   0.004      0.012   .
# 2  PF vs P 0.9687456 0.0323252   0.404      1.000    
# 3  CC vs P 2.9264497 0.1174034   0.012      0.036   .


#only PF compare agg
mds.all=metaMDS(data.selected,distance="bray", k=2)
data.sm=as.data.frame(scores(mds.all))
data.sm$agg = as.factor(GRP$SoilFrac)
ggplot(data=data.sm, aes(x=NMDS1, y=NMDS2, color=agg)) + geom_point(size=5) + 
 geom_hline(yintercept=0.0, colour="grey", lty=2)+
 geom_vline(xintercept=0.0, colour="grey",lty=2) +
 theme_bw() + scale_color_brewer(palette="Set1") + theme(legend.position = c(0.08, 0.82)) +
 theme(legend.background = element_rect(fill="white", size=0.3, linetype="solid", colour="black"))
#+stat_ellipse()
pairwise.adonis(data.selected, as.factor(GRP$SoilFrac))
# [1] "Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1"
         # pairs    F.Model        R2 p.value p.adjusted sig
# 1  MM vs Micro  9.8345371 0.6210814   0.028       0.28    
# 2     MM vs SM  6.6725171 0.5716434   0.024       0.24    
# 3     MM vs LM  9.9772689 0.6244665   0.034       0.34    
# 4     MM vs WS 10.7699642 0.6422175   0.028       0.28    
# 5  Micro vs SM  0.6476039 0.1146688   0.745       1.00    
# 6  Micro vs LM  2.6351489 0.3051654   0.039       0.39    
# 7  Micro vs WS  1.0469403 0.1485666   0.439       1.00    
# 8     SM vs LM  1.6638300 0.2496807   0.180       1.00    
# 9     SM vs WS  0.6489030 0.1148724   0.686       1.00    
# 10    LM vs WS  1.1869516 0.1651537   0.361       1.00 

#only carbon

interest <- c("ec.3.2.1.4", "ec.3.2.1.21","ec.3.2.1.86","ec.3.2.1.91","ec.2.4.1.20","ec.2.7.1.205")
carbon <- data.selected[,interest]
data.selected <- carbon

sub_data <- subset(data.selected, rownames(data.selected) != "H5")
sub_meta <- subset(GRP, FileName.prefix != "H5")
data.selected <- sub_data
GRP <- sub_meta


mds.all=metaMDS(data.selected,distance="bray", k=2)
data.sm=as.data.frame(scores(mds.all))
data.sm$frac = as.factor(GRP$SoilFrac)

ggplot(data=data.sm, aes(x=NMDS1, y=NMDS2, color=frac)) + geom_point(size=5) + 
 geom_hline(yintercept=0.0, colour="grey", lty=2)+
 geom_vline(xintercept=0.0, colour="grey",lty=2) +
 theme_bw() + scale_color_brewer(palette="Set1") + theme(legend.position = c(0.08, 0.82)) +
 theme(legend.background = element_rect(fill="white", size=0.3, linetype="solid", colour="black"))

pairwise.adonis(data.selected, as.factor(GRP$SoilFrac))
#PF only carbon only 
# [1] "Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1"
         # pairs     F.Model         R2 p.value p.adjusted sig
# 1  MM vs Micro 76.37241153 0.92716008   0.038       0.38    
# 2     MM vs SM 99.20194262 0.95201625   0.025       0.25    
# 3     MM vs LM 19.19044091 0.76181441   0.018       0.18    
# 4     MM vs WS 52.78186879 0.89792771   0.029       0.29    
# 5  Micro vs SM  0.08884346 0.01745848   0.851       1.00    
# 6  Micro vs LM  0.19981845 0.03222973   0.866       1.00    
# 7  Micro vs WS  0.60484347 0.09157575   0.500       1.00    
# 8     SM vs LM  0.32771802 0.06151189   0.656       1.00    
# 9     SM vs WS  0.63369463 0.11248296   0.478       1.00    
# 10    LM vs WS  0.15414995 0.02504813   0.777       1.00 

#carbon function not diff in crop

for_plot <- data.frame(count = carbon[,"ec.3.2.1.4"], name = rownames(carbon))
ggplot(for_plot, aes(x=name, y=count))+geom_boxplot()


# compare agg only PF
#normalize?
sub_table <- subset(full_table, Crop == "PF")
sub_table <- full_table
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
data.selected <- t(count)


GRP$Crop == "PF"
sub_data <- subset(data.selected, GRP$Crop == "PF")
sub_data <- subset(sub_data, rownames(sub_data) != "H5")


sub_meta <- subset(GRP, Crop == "PF")
sub_meta <- subset(sub_meta , FileName.prefix != "H5")
dim(GRP)

dim(count)

mds.all=metaMDS(sub_data,distance="bray", k=2)
data.sm=as.data.frame(scores(mds.all))

data.sm$Aggregates = as.factor(sub_meta$SoilFrac)

ggplot(data=data.sm, aes(x=NMDS1, y=NMDS2, color=Aggregates)) + geom_point(size=5) + 
 geom_hline(yintercept=0.0, colour="grey", lty=2)+
 geom_vline(xintercept=0.0, colour="grey",lty=2) +
 theme_bw() + scale_color_brewer(palette="Set1") + theme(legend.position = c(0.08, 0.82)) +
 theme(legend.background = element_rect(fill="white", size=0.3, linetype="solid", colour="black"))

pairwise.adonis(sub_data, as.factor(sub_meta$SoilFrac))

#no diff among agg
[# 1] "Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1"
         # pairs    F.Model        R2 p.value p.adjusted sig
# 1  MM vs Micro  9.7708479 0.6195512   0.020       0.20    
# 2     MM vs SM  6.6608214 0.5712137   0.028       0.28    
# 3     MM vs LM  9.8461874 0.6213600   0.034       0.34    
# 4     MM vs WS 10.6630512 0.6399219   0.040       0.40    
# 5  Micro vs SM  0.6406335 0.1135747   0.745       1.00    
# 6  Micro vs LM  2.5611280 0.2991578   0.035       0.35    
# 7  Micro vs WS  1.0515111 0.1491186   0.399       1.00    
# 8     SM vs LM  1.5758623 0.2396435   0.201       1.00    
# 9     SM vs WS  0.6299204 0.1118880   0.722       1.00    
# 10    LM vs WS  1.1370738 0.1593193   0.411       1.00 



#only carbon
interest <- c("ec.3.2.1.4", "ec.3.2.1.21","ec.3.2.1.86","ec.3.2.1.91","ec.2.4.1.20")
carbon <- sub_data[,interest]
mds.all=metaMDS(carbon,distance="bray", k=2)
data.sm=as.data.frame(scores(mds.all))
data.sm$Aggregates = as.factor(sub_meta$SoilFrac)

pairwise.adonis(sub_data, as.factor(sub_meta$SoilFrac))

#no diff among agg in carbon

#########################################################
# Supplementary figure 3B -> Supplementary Figure 6
#########################################################



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

data.selected <- t(Y)
dim(data.selected)

sub_data <- subset(data.selected, rownames(data.selected) != "H5")
sub_meta <- subset(GRP, FileName.prefix != "H5")
dim(sub_data)
data.selected <- sub_data

mds.all=metaMDS(data.selected,distance="bray", k=2)
data.sm=as.data.frame(scores(mds.all))
data.sm$Aggregates = as.factor(sub_meta$SoilFrac)
data.sm$for_label <- factor(data.sm$Aggregates, labels= c("LM", "Micro","MM","SM","WS"))

pdf("~/Box Sync/2018/2February/cobs/writing/supplementary_figure6.pdf")
ggplot(data=data.sm, aes(x=NMDS1, y=NMDS2, color=for_label)) + geom_point(size=5) +
 geom_hline(yintercept=0.0, colour="grey", lty=2)+
 geom_vline(xintercept=0.0, colour="grey",lty=2) +
 theme_bw() + scale_color_brewer(palette="Set1") + theme(legend.position = c(0.1, 0.85)) +
 theme(legend.background = element_rect(fill="white", size=0.3, linetype="solid", colour="black"))+labs(color="Aggregates")
dev.off()



## MRPP & ANOSIM

trt=as.factor(sub_meta$SoilFrac)
enzymes=as.matrix(data.selected)

mrpp(enzymes, group=trt, distance="bray")
anosim(enzymes, grouping=trt, permutations=999, distance="bray")



pairwise.adonis(data.selected, sub_meta$SoilFrac)
# [1] "Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1"
         # pairs    F.Model         R2 p.value p.adjusted sig
# 1  MM vs Micro 2.79395960 0.31771349   0.101       1.00    
# 2     MM vs SM 5.65194747 0.53060227   0.039       0.39    
# 3     MM vs LM 0.94036450 0.13549209   0.451       1.00    
# 4     MM vs WS 0.40820371 0.06370018   0.687       1.00    
# 5  Micro vs SM 0.07160578 0.01411896   0.916       1.00    
# 6  Micro vs LM 0.49284761 0.07590623   0.909       1.00    
# 7  Micro vs WS 1.08128636 0.15269632   0.350       1.00    
# 8     SM vs LM 0.62422001 0.11098784   0.770       1.00    
# 9     SM vs WS 1.47985899 0.22837827   0.310       1.00    
# 10    LM vs WS 0.16316631 0.02647443   0.714       1.00  
#pairwise.adonis(data.selected, GRP$SoilFrac)

###################
# 3. Figure 1: genomic potential of cellulose decomposition enzymes in FP among agg.
###################

#Figure 1C
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
# Figure 1B

file_name = "~/Box Sync/2018/2February/cobs/writing/figure1C.pdf"
pdf(file_name, width=6, height=9)
one_ec = "ec.3.2.1.4"
d <- plotCounts(deseq2_de, gene = one_ec, intgroup = "SoilFrac", returnData=T)
d$for_label <- factor(d$SoilFrac, labels = c("LM", "Micro","MM","SM","WS"))
p1 <- ggplot(d, aes(x=for_label, y=count))+geom_boxplot()+scale_x_discrete(limits=c("WS", "LM","MM","SM","Micro"))+ggtitle("EC 3.2.1.4")+theme(axis.line=element_line(size=1.25), aspect.ratio=1, panel.border=element_blank(),axis.ticks=element_line(size=1.25, colour="black"), panel.background=element_blank(), axis.text=element_text(size=10, face="bold", colour="black"), axis.title=element_text(size=12, face="bold", colour="black"), plot.title = element_text(size=15, face="bold", hjust=0.5))+labs(x="Aggregates Fraction",y="Normalized Count")


for_aov <- data.frame(count = count[one_ec,], frac =GRP$SoilFrac)
fit <- aov(count ~ frac, data=for_aov)
summary(fit)
            # Df  Sum Sq Mean Sq F value   Pr(>F)    
# frac         4 6549154 1637289   21.42 4.64e-06 ***
# Residuals   15 1146379   76425  


one_ec = "ec.3.2.1.91"
d <- plotCounts(deseq2_de, gene = one_ec, intgroup = "SoilFrac", returnData=T)
d$for_label <- factor(d$SoilFrac, labels = c("LM", "Micro","MM","SM","WS"))
p2 <- ggplot(d, aes(x=for_label, y=count))+geom_boxplot()+scale_x_discrete(limits=c("WS", "LM","MM","SM","Micro"))+ggtitle("EC 3.2.1.91")+theme(axis.line=element_line(size=1.25), aspect.ratio=1, panel.border=element_blank(),axis.ticks=element_line(size=1.25, colour="black"), panel.background=element_blank(), axis.text=element_text(size=10, face="bold", colour="black"), axis.title=element_text(size=12, face="bold", colour="black"), plot.title = element_text(size=15, face="bold", hjust=0.5))+labs(x="Aggregates Fraction",y="Normalized Count")

for_aov <- data.frame(count = count[one_ec,], frac =GRP$SoilFrac)
fit <- aov(count ~ frac, data=for_aov)
summary(fit)
            # Df Sum Sq Mean Sq F value   Pr(>F)    
# frac         4 309191   77298   34.83 1.98e-07 ***
# Residuals   15  33286    2219 

one_ec = "ec.3.2.1.21"
d <- plotCounts(deseq2_de, gene = one_ec, intgroup = "SoilFrac", returnData=T)
d$for_label <- factor(d$SoilFrac, labels = c("LM", "Micro","MM","SM","WS"))
p3 <- ggplot(d, aes(x=for_label, y=count))+geom_boxplot()+scale_x_discrete(limits=c("WS", "LM","MM","SM","Micro"))+ggtitle("EC 3.2.1.21")+theme(axis.line=element_line(size=1.25), aspect.ratio=1, panel.border=element_blank(),axis.ticks=element_line(size=1.25, colour="black"), panel.background=element_blank(), axis.text=element_text(size=10, face="bold", colour="black"), axis.title=element_text(size=12, face="bold", colour="black"), plot.title = element_text(size=15, face="bold", hjust=0.5))+labs(x="Aggregates Fraction",y="Normalized Count")


for_aov <- data.frame(count = count[one_ec,], frac =GRP$SoilFrac)
fit <- aov(count ~ frac, data=for_aov)
summary(fit)
            # Df   Sum Sq  Mean Sq F value   Pr(>F)    
# frac         4 87601920 21900480   9.811 0.000418 ***
# Residuals   15 33483966  2232264 

one_ec = "ec.2.4.1.20"
d <- plotCounts(deseq2_de, gene = one_ec, intgroup = "SoilFrac", returnData=T)
d$for_label <- factor(d$SoilFrac, labels = c("LM", "Micro","MM","SM","WS"))
p4 <- ggplot(d, aes(x=for_label, y=count))+geom_boxplot()+scale_x_discrete(limits=c("WS", "LM","MM","SM","Micro"))+ggtitle("EC 2.4.1.20")+theme(axis.line=element_line(size=1.25), aspect.ratio=1, panel.border=element_blank(),axis.ticks=element_line(size=1.25, colour="black"), panel.background=element_blank(), axis.text=element_text(size=10, face="bold", colour="black"), axis.title=element_text(size=12, face="bold", colour="black"), plot.title = element_text(size=15, face="bold", hjust=0.5))+labs(x="Aggregates Fraction",y="Normalized Count")


for_aov <- data.frame(count = count[one_ec,], frac =GRP$SoilFrac)
fit <- aov(count ~ frac, data=for_aov)
summary(fit)
            # Df Sum Sq Mean Sq F value  Pr(>F)    
# frac         4 225011   56253   12.85 9.7e-05 ***
# Residuals   15  65645    4376 


one_ec = "ec.3.2.1.86"
d <- plotCounts(deseq2_de, gene = one_ec, intgroup = "SoilFrac", returnData=T)
d$for_label <- factor(d$SoilFrac, labels = c("LM", "Micro","MM","SM","WS"))
p5 <- ggplot(d, aes(x=for_label, y=count))+geom_boxplot()+scale_x_discrete(limits=c("WS", "LM","MM","SM","Micro"))+ggtitle("EC 3.2.1.86")+theme(axis.line=element_line(size=1.25), aspect.ratio=1, panel.border=element_blank(),axis.ticks=element_line(size=1.25, colour="black"), panel.background=element_blank(), axis.text=element_text(size=10, face="bold", colour="black"), axis.title=element_text(size=12, face="bold", colour="black"), plot.title = element_text(size=15, face="bold", hjust=0.5))+labs(x="Aggregates Fraction",y="Normalized Count")

for_aov <- data.frame(count = count[one_ec,], frac =GRP$SoilFrac)
fit <- aov(count ~ frac, data=for_aov)
summary(fit)
            # Df  Sum Sq Mean Sq F value  Pr(>F)   
# frac         4 1199699  299925   5.679 0.00547 **
# Residuals   15  792160   52811 


# one_ec = "ec.2.7.1.205"
# d <- plotCounts(deseq2_de, gene = one_ec, intgroup = "SoilFrac", returnData=T)
# p6 <- ggplot(d, aes(x=SoilFrac, y=count))+geom_boxplot()+scale_x_discrete(limits=c("WS", "LM","MM","SM","Micro"))+ggtitle("EC 2.7.1.205")+theme(axis.line=element_line(size=1.25), aspect.ratio=1, panel.border=element_blank(),axis.ticks=element_line(size=1.25, colour="black"), panel.background=element_blank(), axis.text=element_text(size=10, face="bold", colour="black"), axis.title=element_text(size=12, face="bold", colour="black"), plot.title = element_text(size=15, face="bold", hjust=0.5))+labs(x="Aggregates Fraction",y="Normalized Count")

# for_aov <- data.frame(count = count[one_ec,], frac =GRP$SoilFrac)
# fit <- aov(count ~ frac, data=for_aov)
# summary(fit)
            # # Df Sum Sq Mean Sq F value Pr(>F)
# # frac         4    3.0   0.750   0.079  0.988
# # Residuals   15  142.3   9.484

library("gridExtra")
grid.arrange(p1,p2,p3,p4,p5, nrow=3)
dev.off()



#########################################################
### figure 1C -> 1B
##############################

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
# Deseq2
sub_table <- full_table
#sub_table <- subset(full_table, SoilFrac == "WS")

Y = t(sub_table[, 2:(ncol(fin_table)+1)])
colnames(Y) <- sub_table$FileName.prefix
GRP = sub_table[,(ncol(fin_table)+2):ncol(sub_table)]


library(DESeq2)
library(ggplot2)
## Create a DESeqDataSet object for holding the counts and phenotype data.
deseq2_deseqds <- DESeqDataSetFromMatrix(countData = Y, colData = GRP, design = ~ Crop)

## Differential expression analysis. 
deseq2_de <- DESeq(deseq2_deseqds)

# significantly different between medium and whole sooil
#deseq2_de_out <- results(deseq2_de, contrast=c("SoilFrac","MM","WS"))

count <- counts(deseq2_de, normalized=TRUE)


####################
#  cellulose decomposition enzymes
# Figure 1C -> 1B

file_name = "~/Box Sync/2018/2February/cobs/writing/figure1B.pdf"
pdf(file_name, width=6, height=9)
one_ec = "ec.3.2.1.4"
d <- plotCounts(deseq2_de, gene = one_ec, intgroup = "Crop", returnData=T)
d$for_label <- factor(d$Crop, labels = c("C", "P","FP"))
p1 <- ggplot(d, aes(x=for_label, y=count))+geom_boxplot()+scale_x_discrete(limits=c("C", "P","FP"))+ggtitle("EC 3.2.1.4")+theme(axis.line=element_line(size=1.25), aspect.ratio=1, panel.border=element_blank(),axis.ticks=element_line(size=1.25, colour="black"), panel.background=element_blank(), axis.text=element_text(size=10, face="bold", colour="black"), axis.title=element_text(size=12, face="bold", colour="black"), plot.title = element_text(size=15, face="bold", hjust=0.5))+labs(x="Crop System",y="Normalized Count")


for_aov <- data.frame(count = count[one_ec,], frac =GRP$Crop)
fit <- aov(count ~ frac, data=for_aov)
summary(fit)

#WS
            # Df Sum Sq Mean Sq F value Pr(>F)
# frac         2  89895   44947   2.272  0.159
# Residuals    9 178012   19779  
#all soil
            # Df   Sum Sq Mean Sq F value Pr(>F)
# frac         2   646703  323351   1.101  0.342
# Residuals   41 12038627  293625 


one_ec = "ec.3.2.1.91"
d <- plotCounts(deseq2_de, gene = one_ec, intgroup = "Crop", returnData=T)
d$for_label <- factor(d$Crop, labels = c("C", "P","FP"))
p2 <-  ggplot(d, aes(x=for_label, y=count))+geom_boxplot()+scale_x_discrete(limits=c("C", "P","FP"))+ggtitle("EC 3.2.1.91")+theme(axis.line=element_line(size=1.25), aspect.ratio=1, panel.border=element_blank(),axis.ticks=element_line(size=1.25, colour="black"), panel.background=element_blank(), axis.text=element_text(size=10, face="bold", colour="black"), axis.title=element_text(size=12, face="bold", colour="black"), plot.title = element_text(size=15, face="bold", hjust=0.5))+labs(x="Crop System",y="Normalized Count")

for_aov <- data.frame(count = count[one_ec,], frac =GRP$Crop)
fit <- aov(count ~ frac, data=for_aov)
summary(fit)

#WS
            # Df Sum Sq Mean Sq F value Pr(>F)
# frac         2    732   366.0   0.388  0.689
# Residuals    9   8495   943.9  
#all soil
            # Df Sum Sq Mean Sq F value Pr(>F)
# frac         2  22186   11093   1.008  0.374
# Residuals   41 451397   11010  

one_ec = "ec.3.2.1.21"
d <- plotCounts(deseq2_de, gene = one_ec, intgroup = "Crop", returnData=T)
d$for_label <- factor(d$Crop, labels = c("C", "P","FP"))
p3 <-  ggplot(d, aes(x=for_label, y=count))+geom_boxplot()+scale_x_discrete(limits=c("C", "P","FP"))+ggtitle("EC 3.2.1.21")+theme(axis.line=element_line(size=1.25), aspect.ratio=1, panel.border=element_blank(),axis.ticks=element_line(size=1.25, colour="black"), panel.background=element_blank(), axis.text=element_text(size=10, face="bold", colour="black"), axis.title=element_text(size=12, face="bold", colour="black"), plot.title = element_text(size=15, face="bold", hjust=0.5))+labs(x="Crop System",y="Normalized Count")


for_aov <- data.frame(count = count[one_ec,], frac =GRP$Crop)
fit <- aov(count ~ frac, data=for_aov)
summary(fit)
            # WS
            # Df  Sum Sq Mean Sq F value Pr(>F)
# frac         2 3832901 1916451   2.131  0.175
# Residuals    9 8092311  899146 

#all soil
            # Df    Sum Sq Mean Sq F value Pr(>F)
# frac         2  10682075 5341037   1.048   0.36
# Residuals   41 208943430 5096181  

one_ec = "ec.2.4.1.20"
d <- plotCounts(deseq2_de, gene = one_ec, intgroup = "Crop", returnData=T)
d$for_label <- factor(d$Crop, labels = c("C", "P","FP"))
p4 <-  ggplot(d, aes(x=for_label, y=count))+geom_boxplot()+scale_x_discrete(limits=c("C", "P","FP"))+ggtitle("EC 2.4.1.20")+theme(axis.line=element_line(size=1.25), aspect.ratio=1, panel.border=element_blank(),axis.ticks=element_line(size=1.25, colour="black"), panel.background=element_blank(), axis.text=element_text(size=10, face="bold", colour="black"), axis.title=element_text(size=12, face="bold", colour="black"), plot.title = element_text(size=15, face="bold", hjust=0.5))+labs(x="Crop System",y="Normalized Count")


for_aov <- data.frame(count = count[one_ec,], frac =GRP$Crop)
fit <- aov(count ~ frac, data=for_aov)
summary(fit)
            # WS
                        # Df Sum Sq Mean Sq F value Pr(>F)
# frac         2   6706    3353   0.658  0.541
# Residuals    9  45885    5098   
#all soil
            # Df Sum Sq Mean Sq F value Pr(>F)  
# frac         2  60364   30182   2.836 0.0702 .
# Residuals   41 436373   10643 


one_ec = "ec.3.2.1.86"
d <- plotCounts(deseq2_de, gene = one_ec, intgroup = "Crop", returnData=T)
d$for_label <- factor(d$Crop, labels = c("C", "P","FP"))
p5 <-  ggplot(d, aes(x=for_label, y=count))+geom_boxplot()+scale_x_discrete(limits=c("C", "P","FP"))+ggtitle("EC 3.2.1.86")+theme(axis.line=element_line(size=1.25), aspect.ratio=1, panel.border=element_blank(),axis.ticks=element_line(size=1.25, colour="black"), panel.background=element_blank(), axis.text=element_text(size=10, face="bold", colour="black"), axis.title=element_text(size=12, face="bold", colour="black"), plot.title = element_text(size=15, face="bold", hjust=0.5))+labs(x="Crop System",y="Normalized Count")

for_aov <- data.frame(count = count[one_ec,], frac =GRP$Crop)
fit <- aov(count ~ frac, data=for_aov)
summary(fit)
            # WS
                       # Df Sum Sq Mean Sq F value Pr(>F)
# frac         2  28711   14356   0.639   0.55
# Residuals    9 202131   22459    
#all soil
            # Df  Sum Sq Mean Sq F value Pr(>F)
# frac         2  136637   68318    0.99   0.38
# Residuals   41 2828435   68986

# one_ec = "ec.2.7.1.205"
# d <- plotCounts(deseq2_de, gene = one_ec, intgroup = "SoilFrac", returnData=T)
# p6 <- ggplot(d, aes(x=SoilFrac, y=count))+geom_boxplot()+scale_x_discrete(limits=c("WS", "LM","MM","SM","Micro"))+ggtitle("EC 2.7.1.205")+theme(axis.line=element_line(size=1.25), aspect.ratio=1, panel.border=element_blank(),axis.ticks=element_line(size=1.25, colour="black"), panel.background=element_blank(), axis.text=element_text(size=10, face="bold", colour="black"), axis.title=element_text(size=12, face="bold", colour="black"), plot.title = element_text(size=15, face="bold", hjust=0.5))+labs(x="Aggregates Fraction",y="Normalized Count")

# for_aov <- data.frame(count = count[one_ec,], frac =GRP$SoilFrac)
# fit <- aov(count ~ frac, data=for_aov)
# summary(fit)
            # # Df Sum Sq Mean Sq F value Pr(>F)
# # frac         4    3.0   0.750   0.079  0.988
# # Residuals   15  142.3   9.484

library("gridExtra")
grid.arrange(p1,p2,p3,p4,p5, nrow=3)
dev.off()







#########################################################
# Supplementary figure 4. Heatmap -> Supplemetary Figure 7
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

#this for class by pathway
kegg_class <- read.table("~/Box Sync/Github/metajinomics/dev/kegg_tools/kegg_full_class.tsv", sep = '\t')
kegg_class$ec = gsub(":",".", kegg_class$V1)


pdf("~/Box Sync/2018/2February/cobs/writing/supplementary_figure_heatmap.pdf", width=8, height=8)
#only up : High in MM
deseq2_de_out_sig_up <- subset(deseq2_de_out_sig, log2FoldChange > 0)
idx <- rownames(deseq2_de_out_sig_up)[which(deseq2_de_out_sig_up$log2FoldChange > 0)]
range(deseq2_de_out_sig_up$log2FoldChange)
length(idx)


p1 <- heatmap.2(mat[idx,], trace = "none" , col = greenred(10), scale="row", margins=c(12,8), labRow=NA, labCol=GRP$SoilFrac)
plot(p1)


subseted_table <- subset(kegg_class, ec %in% idx)
one_category <- subset(subseted_table, V4 == "Carbohydrate metabolism")

p1 <- heatmap.2(mat[one_category$ec,], trace = "none" , col = greenred(10), scale="row", margins=c(10,18), labRow=one_category$V3, labCol=GRP$SoilFrac)
plot(p1)

#only down: High in WS
deseq2_de_out_sig_down <- subset(deseq2_de_out_sig, log2FoldChange < 0)
idx <- rownames(deseq2_de_out_sig_down)[which(deseq2_de_out_sig_down$log2FoldChange <  0)]

p1<-heatmap.2(mat[idx,], trace = "none" , col = greenred(10), scale="row", margins=c(12,8),labRow=NA, labCol=GRP$SoilFrac)
plot(p1)


subseted_table <- subset(kegg_class, ec %in% idx)
one_category <- subset(subseted_table, V4 == "Carbohydrate metabolism")

p1<-heatmap.2(mat[one_category$ec,], trace = "none" , col = greenred(10), scale="row", margins=c(10,18), labRow=one_category$V3, labCol=GRP$SoilFrac)
plot(p1)
dev.off()


#########################################################
# 5. Co-occurrence -> Figure 4
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
# Figure 3. compare agg in organism -> Figure 2
#########################################################
# Read meta
meta=read.csv("~/Box Sync/2017/10October/cobs/sample_data_cobs.csv")


interest <- c("3.2.1.4", "3.2.1.21","3.2.1.86","3.2.1.91","2.4.1.20")

for_plot_table = data.frame()
for(j in 1:length(interest)) {
	ec = interest[1]
ec = interest[j]


path = paste0("~/Box Sync/2017/3March/cobs/data/",ec)
file_abun = paste0(path, ".fasta.summary_count.tsv")
abun=read.table(file_abun,header=T, row.names=1)


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


library(ggplot2)
########################
## sum plot
# temp <- data.frame(name =  rownames(fin_table), count = rowSums(fin_table))

# #use only PF
# pf <- subset(temp, Crop=="PF")
# pf$nor <- pf$count / pf$depth


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
temp_otu <- otu_table(physeq)

for (i in 1:length(sample_data(physeq)$depth)){
	temp_otu[,i]  <- temp_otu[,i] / sample_data(physeq)$depth[i]
}
otu_table(physeq) <- temp_otu

#use only pf
pf <- prune_samples(sample_data(physeq)$Crop == "PF", physeq)
physeq <- pf

temp_table  = data.frame()
for (group in unique(sample_data(physeq)$SoilFrac)){
	temp_physeq = prune_samples(sample_data(physeq)$SoilFrac == group, physeq)

	new_table <- data.frame(taxa_sums(temp_physeq), tax_table(temp_physeq)[,2])
	colnames(new_table) = c("abundance", "phylum")
	sum_table = data.frame()
	for (x in unique(new_table$phylum) ){
		temp = subset(new_table, phylum==x)
		su = sum(temp$abundance)
		sum_table = rbind(sum_table, data.frame(su, temp[1,2]))
	}
	colnames(sum_table) = c("abundance", "phylum")
	perc_table = sum_table
	for (i in 1:nrow(sum_table)){
			perc_table[i,1] = sum_table[i,1] / sum(sum_table[,1])
	}
	
	other_table = data.frame()
	other = c()
	for (i in 1:nrow(perc_table)){
		if(perc_table[i,1] > 0.01) {
			other_table = rbind(other_table, perc_table[i,])
		}else{
			other = c(other, perc_table[i,1])
		}
		
	}
	
	sum(other)
	tep = data.frame(sum(other), "other")
	colnames(tep) = c("abundance", "phylum")
	#tfin = rbind(other_table, tep)
	tfin <- other_table
	ttfin = cbind(tfin,group)
	ttfin = cbind(ttfin,ec)
	temp_table = rbind(temp_table, ttfin)
	phy = unique(temp_table$phylum)
}

temp_table$group = factor(temp_table$group, levels = c("Micro","SM","MM","LM","WS"))
for_plot_table = rbind(for_plot_table, temp_table)


} #end of for

for_plot_table$ec_label = factor(for_plot_table$ec, labels = c("EC 3.2.1.4","EC 3.2.1.21","EC 3.2.1.86","EC 3.2.1.91","EC 2.4.1.20"))

pdf("~/Box Sync/2018/2February/cobs/writing/figure2.pdf", width=11, height=6)
ggplot(for_plot_table, aes(x=group,y=abundance, fill=phylum))+geom_bar(stat="identity")+scale_x_discrete(limits=c("WS", "LM","MM","SM","Micro"))+labs(y="Relative Abundance", x="Aggregates Fraction", fill="Phylum")+facet_grid(~ec_label)+theme_bw() 
dev.off()
#+ scale_fill_brewer(palette="Spectral")

#########################################################
# Supplementary figure 6: compare crop in organism -> Supplementary Figure 4
#########################################################
###################
# Read meta
meta=read.csv("~/Box Sync/2017/10October/cobs/sample_data_cobs.csv")

#interest <- c("3.2.1.4", "3.2.1.21","3.2.1.86","3.2.1.91","2.4.1.20","2.7.1.205")
interest <- c("3.2.1.4", "3.2.1.21","3.2.1.86","3.2.1.91","2.4.1.20")




for_plot_table = data.frame()
for(j in 1:length(interest)) {
	ec = interest[1]
ec = interest[j]


path = paste0("~/Box Sync/2017/3March/cobs/data/",ec)
file_abun = paste0(path, ".fasta.summary_count.tsv")
abun=read.table(file_abun,header=T, row.names=1)


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
temp_otu <- otu_table(physeq)

for (i in 1:length(sample_data(physeq)$depth)){
	temp_otu[,i]  <- temp_otu[,i] / sample_data(physeq)$depth[i]
}
otu_table(physeq) <- temp_otu

#sample_data(physeq)$Crop == "PF"
pf <- prune_samples(sample_data(physeq)$Crop == "PF", physeq)

#plot_bar(pf, fill="Phylum")

gl <- tax_glom(pf, "Phylum")
#here


p<-plot_bar(gl, fill = "Phylum")
plot(p)



temp_table  = data.frame()
for (group in unique(sample_data(physeq)$Crop)){
	temp_physeq = prune_samples(sample_data(physeq)$Crop == group, physeq)

	new_table <- data.frame(taxa_sums(temp_physeq), tax_table(temp_physeq)[,2])
	colnames(new_table) = c("abundance", "phylum")
	sum_table = data.frame()
	for (x in unique(new_table$phylum) ){
		temp = subset(new_table, phylum==x)
		su = sum(temp$abundance)
		sum_table = rbind(sum_table, data.frame(su, temp[1,2]))
	}
	colnames(sum_table) = c("abundance", "phylum")
	perc_table = sum_table
	for (i in 1:nrow(sum_table)){
			perc_table[i,1] = sum_table[i,1] / sum(sum_table[,1])
	}
	
	other_table = data.frame()
	other = c()
	for (i in 1:nrow(perc_table)){
		if(perc_table[i,1] > 0.01) {
			other_table = rbind(other_table, perc_table[i,])
		}else{
			other = c(other, perc_table[i,1])
		}
		
	}
	
	sum(other)
	tep = data.frame(sum(other), "other")
	colnames(tep) = c("abundance", "phylum")
	#tfin = rbind(other_table, tep)
	tfin <- other_table
	ttfin = cbind(tfin,group)
	ttfin = cbind(ttfin,ec)
	temp_table = rbind(temp_table, ttfin)
	phy = unique(temp_table$phylum)
}
temp_table <- subset(temp_table, phylum != "null")
temp_table$group = factor(temp_table$group, levels = c("CC","P","PF"))
for_plot_table = rbind(for_plot_table, temp_table)


} #end of for

for_plot_table$ec_label = factor(for_plot_table$ec, labels = c("EC 3.2.1.4","EC 3.2.1.21","EC 3.2.1.86","EC 3.2.1.91","EC 2.4.1.20"))
for_plot_table$crop_label = factor(for_plot_table$group, labels = c("C","P","FP"))

#### this is end of plot
pdf("~/Box Sync/2018/2February/cobs/writing/supplementary_figure4.pdf")
ggplot(for_plot_table, aes(x=crop_label,y=abundance, fill=phylum))+geom_bar(stat="identity")+labs(y="Relative Abundance", x="Crop System", fill="Phylum")+facet_grid(~ec_label)
dev.off()



interest <- c("3.2.1.4", "3.2.1.21","3.2.1.86","3.2.1.91","2.4.1.20")

tax_name <- tax_table(gl)[,2][2]
one_name <- rownames(tax_table(gl))[2]
tm <- data.frame(count = c(otu_table(gl)[one_name,]), name = colnames(otu_table(gl)), SoilFrac=sample_data(gl)$SoilFrac)
ggplot(tm, aes(x=SoilFrac,y=count))+geom_boxplot()
test <- aov(count ~ SoilFrac, data=tm)
summary(test)


###############################
# figure 4 -> Figure3
###############################
#start here
interest <- c("3.2.1.4", "3.2.1.21","3.2.1.86","3.2.1.91","2.4.1.20")
for_plot <- data.frame()
pdf("~/Box Sync/2018/2February/cobs/writing/figure3.pdf", width=3, heigh=3)
for(j in 1:length(interest)) {
	ec = interest[1]
ec = interest[j]


path = paste0("~/Box Sync/2017/3March/cobs/data/",ec)
file_abun = paste0(path, ".fasta.summary_count.tsv")
abun=read.table(file_abun,header=T, row.names=1)


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
temp_otu <- otu_table(physeq)

for (i in 1:length(sample_data(physeq)$depth)){
	temp_otu[,i]  <- temp_otu[,i] / sample_data(physeq)$depth[i]
}
otu_table(physeq) <- temp_otu

#sample_data(physeq)$Crop == "PF"
pf <- prune_samples(sample_data(physeq)$Crop == "PF", physeq)

#plot_bar(pf, fill="Phylum")

gl <- tax_glom(pf, "Phylum")
#here


for (i in 1:nrow(tax_table(gl))){
otu_name = rownames(tax_table(gl))[i]
phy_name <- tax_table(gl)[otu_name,"Phylum"]
co <- t(otu_table(gl)[otu_name,])
x_name <- sample_data(gl)$SoilFrac
for_aov <- data.frame(coun=co, name=x_name, phy=phy_name, ec=ec)
colnames(for_aov) <- c("count","name","phy","ec")
fit <- aov(count ~ name, data=for_aov)
if( !is.na(summary(fit)[[1]][["Pr(>F)"]][1] )){
	if (summary(fit)[[1]][["Pr(>F)"]][1] < 0.01){
		
		p1 <- ggplot(for_aov, aes(x=name, y=count))+geom_boxplot()+scale_x_discrete(limits=c("WS", "LM","MM","SM","Micro"))+ggtitle(phy_name)+xlab("Aggregate")+ylab("Normalized Count")+theme(axis.line=element_line(size=1.25), aspect.ratio=1, panel.border=element_blank(),axis.ticks=element_line(size=1.25, colour="black"), panel.background=element_blank(), axis.text=element_text(size=10, face="bold", colour="black"), axis.title=element_text(size=12, face="bold", colour="black"), plot.title = element_text(size=15, face="bold", hjust=0.5))
		plot(p1)
		for_plot <- rbind(for_plot, for_aov)
	}
}

} #end of for

}#end of for

dev.off()

#########################################################
# Supplementary figure 5 EC bar -> Supplementary Figure 3
#########################################################
rm(list=ls())
#read kegg
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
#whole soiil?

#full?
sub_table <- full_table

Y = t(sub_table[, 2:(ncol(fin_table)+1)])
colnames(Y) <- sub_table$FileName.prefix
GRP = sub_table[,(ncol(fin_table)+2):ncol(sub_table)]

## Create a DESeqDataSet object for holding the counts and phenotype data.
deseq2_deseqds <- DESeqDataSetFromMatrix(countData = Y, colData = GRP, design = ~ SoilFrac)

## Differential expression analysis. 
deseq2_de <- DESeq(deseq2_deseqds)

#extract normalized count
count <- counts(deseq2_de, normalized=TRUE)


interest <- c("ec.3.2.1.4", "ec.3.2.1.21","ec.3.2.1.86","ec.3.2.1.91","ec.2.4.1.20")
inte_table <- count[interest,]
dim(inte_table)
GRP$Crop == group
inte_table[,GRP$Crop == group]

temp_table  = data.frame()
for (group in unique(GRP$Crop)){
	temp_df = inte_table[,GRP$Crop == group]
	
	sum_table <- data.frame(rowSums(temp_df), rownames(temp_df))
	colnames(sum_table) = c("abundance", "enzyme")

	perc_table = sum_table
	for (i in 1:nrow(sum_table)){
			perc_table[i,1] = sum_table[i,1] / sum(sum_table[,1])
	}
	
	#if anything less than 1% then merge into others
	# other_table = data.frame()
	# other = c()
	# for (i in 1:nrow(perc_table)){
		# if(perc_table[i,1] > 0.001) {
			# other_table = rbind(other_table, perc_table[i,])
		# }else{
			# other = c(other, perc_table[i,1])
		# }
		
	# }
	
	#sum(other)
	# tep = data.frame(sum(other), "other")
	# colnames(tep) = c("abundance", "enzyme")
	
	#tfin = rbind(other_table, tep)
	#ttfin = cbind(tfin,group)
	
	other_table <- perc_table
	ttfin = cbind(other_table, group)

	temp_table = rbind(temp_table, ttfin)
	#phy = unique(temp_table$phylum)
}

temp_table$group = factor(temp_table$group, levels = c("CC","P","PF"))
for_plot_table = temp_table


for_plot_table$ec_label = factor(for_plot_table$enzyme, labels = c("EC 2.4.1.20","EC 3.2.1.21","EC3.2.1.4","EC 3.2.1.86","EC 3.2.1.91"))
for_plot_table$crop_label = factor(for_plot_table$group, labels = c("C", "P","FP"))

dev.off()
#### this is end of plot
pdf("~/Box Sync/2018/2February/cobs/writing/supplementary_figure3.pdf")
ggplot(for_plot_table, aes(x=crop_label,y=abundance, fill=ec_label))+geom_bar(stat="identity")+labs(y="Relative Abundance", x="Crop System", fill = "KEGG Enzyme")+theme_bw()
dev.off()

#########################################################
# Supplementary figure 7
#########################################################
rm(list=ls())
#read kegg
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

## Create a DESeqDataSet object for holding the counts and phenotype data.
deseq2_deseqds <- DESeqDataSetFromMatrix(countData = Y, colData = GRP, design = ~ SoilFrac)

## Differential expression analysis. 
deseq2_de <- DESeq(deseq2_deseqds)

#extract normalized count
count <- counts(deseq2_de, normalized=TRUE)


interest <- c("ec.3.2.1.4", "ec.3.2.1.21","ec.3.2.1.86","ec.3.2.1.91","ec.2.4.1.20")
inte_table <- count[interest,]
dim(inte_table)
GRP$Crop == group
inte_table[,GRP$SoilFrac == group]

temp_table  = data.frame()
for (group in unique(GRP$SoilFrac)){
	temp_df = inte_table[,GRP$SoilFrac == group]
	
	sum_table <- data.frame(rowSums(temp_df), rownames(temp_df))
	colnames(sum_table) = c("abundance", "enzyme")

	perc_table = sum_table
	for (i in 1:nrow(sum_table)){
			perc_table[i,1] = sum_table[i,1] / sum(sum_table[,1])
	}
		
	other_table <- perc_table
	ttfin = cbind(other_table, group)

	temp_table = rbind(temp_table, ttfin)
	#phy = unique(temp_table$phylum)
}

#temp_table$group = factor(temp_table$group, levels = c(""))
for_plot_table = temp_table


for_plot_table$ec_label = factor(for_plot_table$enzyme, labels = c("EC 2.4.1.20","EC 3.2.1.21","EC3.2.1.4","EC 3.2.1.86","EC 3.2.1.91"))
for_plot_table$agg_label = factor(for_plot_table$group, labels = c("MM","Micro","SM","LM","WS"))

dev.off()
#### this is end of plot
pdf("~/Box Sync/2018/2February/cobs/writing/supplementary_figure7.pdf")
ggplot(for_plot_table, aes(x=agg_label,y=abundance, fill=ec_label))+geom_bar(stat="identity")+scale_x_discrete(limits=c("WS", "LM","MM","SM","Micro"))+labs(y="Relative Abundance", x="Aggregates Fraction", fill = "KEGG Enzyme")+theme_bw()
dev.off()
