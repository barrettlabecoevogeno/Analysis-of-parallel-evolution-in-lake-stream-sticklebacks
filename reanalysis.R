library(edgeR)
library(pheatmap)
library(ggbiplot)
library(statmod)
library(VennDiagram)
library(made4)
library(ade4)
library(geomorph)
library(amap)
library(psych)
library(vegan)
library(scales)
setwd("~/Desktop/count")

#Read in sample table
samples=read.csv("samples.csv", header = T, stringsAsFactors = F)
samples$countf=paste(samples$libraryname, "count", sep = ".")
samples

#Read in and merge individual htseq-count files
counts=readDGE(samples$countf)$counts
counts

#First filter weakly expressed and noninformative features in half of samples
noint = rownames(counts) %in% c("__no_feature","__ambiguous","__too_low_aQual", "__not_aligned")
cpms=cpm(counts)
keep=rowSums(cpms>1) >=3 & !noint 
counts=counts[keep,]

#Normalize all libraries withthe trimmed mean of M-value (TMM)
d=DGEList(counts=counts)
d=calcNormFactors(d)

#Check sacling factors for libraries, any libraries with scaling factor not close to 1 should be excluded from further analysis
d$samples

#Inspect the relationships between samples using a multidimensional scaling (MDS) plot 
#Using BCV as a measurment of samples relationship
plotMDS(d,
        method = "bcv",
        labels = samples$libraryname, 
        col = c("grey", "black")[factor(samples$location)]
       )

## PCA analyis ##############################################

#Log-tranform normalized read counts of all expressed genes
logcpm=cpm(d, prior.count = 2, log = T)

#Exchange rows and colums of logcpm for PCA 
logcpm1=data.frame(t(logcpm))

#Define variables
variable=cbind(samples$habitat, samples$location)
variable

#Princple Component Analysis
pcal=prcomp(logcpm1, center = T)

#Summarize PCA results
summary(pcal)

#Define factors
hab=as.factor(variable[,1])
shed=as.factor(variable[,2])

#Visulize PCA results
g1=ggbiplot(pcal, 
            obs.scale = 1, 
            var.scale = 1,
            varname.size = 0,
            var.axes = FALSE
)
g2=g1+geom_point(aes(colour=hab, fill=hab, shape=shed), size=2)+
  scale_shape_manual(values = c(21,22))+
  scale_fill_manual(values = c("grey", "black"))+
  scale_color_manual(values = c("grey", "black"))
g3=g2 + theme_bw()+theme(panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank())
print(g3)

###############
###PERMANOVA###
###############

#Calculate pearson distance between each sample

m<-Dist(logcpm1, method="pearson", diag=T, upper=T)

adonis(m~hab+shed+hab*shed, permutations = 25000)


##########################
##### Method 1 edgeR #####
##########################

design <- model.matrix(~hab+shed+hab*shed, samples) #Make design matrix
rownames(design)=colnames(d) #Give meaningful rownames to design matrix
d2=estimateDisp(d, design) #Estimate common dispersion, trended dispersionds, and tagwise dispersions in one run
f=glmFit(d2, design) #Given the design matrix and dispersion estimates, fit a GLM to each feature

##Likelihood ratio test to see if MS is sig. different from ML
lrt.m <- glmLRT(f, contrast = c(0,1,0,0))
detags.m <- rownames(topTags(lrt.m)) #Get the names of the top 10 de tags (can change 10 to x using n=x)
m.exp<-as.data.frame(t(as.data.frame(cpm(d)[detags.m,]))) #Get the expression levels for the top 10 de tags

#Add habitat and watershed data
m.exp$habitat<-hab
m.exp$shed<-shed 

m.exp$hab.shed<-paste(m.exp$shed, m.exp$habitat) #Combine habitat and watershed for plotting purposes
summary(dt.m <- decideTestsDGE(lrt.m)) #Show how many genes down and up regulated--negative means downregulated in stream

isDE.m <- as.logical(dt.m) #Get index of each de gene
DEnames.m <- rownames(d)[isDE.m] #Get name of each de gene

m.stat<-as.data.frame(lrt.m$table)
m.stat$gene_id<-row.names(m.stat)

##Likelihood ratio test to see if RS is sig. different from RL
lrt.r <- glmLRT(f, contrast = c(0,1,0,1))
detags.r<-rownames(topTags(lrt.r)) #Get the names of the top 10 de tags (can change 10 to x using n=x)
r.exp<-as.data.frame(t(as.data.frame(cpm(d)[detags.r,]))) #Get the expression levels for the top 10 de tags

#Add habitat and watershed data
r.exp$habitat<-hab
r.exp$shed<-shed 

r.exp$hab.shed<-paste(r.exp$shed, r.exp$habitat) #Combine habitat and watershed for plotting purposes
summary(dt.r <- decideTestsDGE(lrt.r)) #Show how many genes down and up regulated
isDE.r <- as.logical(dt.r) #Get index of each de gene
DEnames.r <- rownames(d)[isDE.r] #Get names of each de gene

##Likelihood ratio test to see if the interaction coefficient is sig. different from zero
lrt.int <- glmLRT(f, contrast=c(0,0,0,1))
detags.int<-rownames(topTags(lrt.int)) #Get the names of the top 10 de tags (can change 10 to x using n=x)
int.exp<-as.data.frame(t(as.data.frame(cpm(d)[detags.int,]))) #Get the expression levels for the top 10 de tags
int.exp$habitat<-hab #Add habitat data
int.exp$shed<-shed #Add watershed data
int.exp$hab.shed<-paste(int.exp$shed, int.exp$habitat) #Combine habitat and watershed for plotting purposes
summary(dt.int <- decideTestsDGE(lrt.int)) #Show how many genes down and up regulated
isDE.int <- as.logical(dt.int) #Get index of each de gene
DEnames.int <- rownames(d)[isDE.int] #Get names of each de gene

##Now find genes that are DE in both Misty and Roberts
both<-which(DEnames.r %in% DEnames.m) #Get indices of genes found in both
r.m<-DEnames.r[both] #Make list of gene ids for genes indexed above

##Get list of genes sig. for habitat but not for interaction

getrid<-which(r.m %in% DEnames.int) #Genes that are sig. different for habitat and for interaction
no.int<-r.m[-getrid] #Get only names for genes that have no sig. interaction

#Get list of anti-parallel genes
anti1<-r.m[getrid] 

#Get expression levels for all genes sig. for habitat but not for interaction
exp.no.int<-as.data.frame(cpm(d)[no.int,]) 

#Make new column with gene id
exp.no.int$gene_id<-rownames(exp.no.int) 

#Evaluated the probability of the observed overlap between DE gene in Misty and Robert's
n_r=length(DEnames.r)
n_m=length(DEnames.m)
n_r_m=length(both)
n_all=nrow(counts)
phyper(n_r_m, n_r, n_all-n_r, n_m, lower.tail = F)

##make Venn Diagram for edgeR method#####
venn.diagram(list(DEnames.m, DEnames.r, DEnames.int), filename = "edgeR.png",
             main="GLM Method", main.cex=3.5,
             imagetype = "png", fontface=c(1,2,1,1,1,1,1),
             height=6000, width=6000,
             fill=c("gray60","gray86", "white"),
             cex = c(4.3, 4.6,4.3,4.3,4.3,4.3,4.3), cat.cex=3,
             scaled=F, lwd=c(1,1,1), lty=c(1,1,1), col=c("black","black","black"), 
             cat.pos = c(320,40,180), cat.dist = c(0.08,0.08,0.08),
             category.names=c("Significant in\nMisty",
                              " Significant in\n    Robert's",
                              "Significant\nInteraction"))

#####################################################################################################################################
##Method 2: do 2 individual pairwise tests (1 for M and 1 for R), then get genes found in both with same sign (as in Ghalambor2015)##
#####################################################################################################################################

targets <- samples #import targets file
d2 <- DGEList(counts = counts, group = samples$location) #read in the count files

colnames(d2) <- c("ML1","ML2","ML3","MS1","MS2","MS3","RL1","RL2", "RL3",
                  "RS1", "RS2", "RS3")
d2$samples$group<-as.factor(c("ML", "ML", "ML", "MS", "MS", "MS","RL", "RL", "RL", "RS", "RS", "RS"))

d2 <- calcNormFactors(d2) #calculate the normalization factors

d2 <- estimateCommonDisp(d2) #estimate common dispersion
d2<- estimateTagwiseDisp(d2) #estimate individual gene dispersion

de = exactTest(d2, pair=c("ML","MS")) #test for de genes between Misty L and S
tt = topTags(de, n=nrow(d2)) #get stats for each test including FDR for Misty 
head(tt$table) #see some de genes
detags2 <- rownames(topTags(de)) #get names of top de genes
cpm(d2)[detags2,] #expression data for top de genes

de2<-exactTest(d2, pair=c("RL", "RS")) #test for de genes between Roberts L and S
tt2 = topTags(de2, n=nrow(d2)) #get stats for each test including FDR for Roberts
head(tt2$table) #see some de genes

sum(tt2$table$FDR<0.05) #sig genes in roberts
sum(tt$table$FDR<0.05) #sig genes in misty

sig.mis<-row.names(subset(tt$table, tt$table$FDR < 0.05)) #names of genes sig in misty
sig.rob<-row.names(subset(tt2$table, tt2$table$FDR < 0.05)) #names of genes sig in roberts

all.method2<-merge(tt$table, tt2$table, by="row.names")
all.method2$agreement<-as.factor(ifelse(sign(all.method2$logFC.x)==sign(all.method2$logFC.y), "same","different")) #make new column based on whether log-FC is the same sign in M and R
diff<-subset(all.method2, all.method2$agreement=="same", select=Row.names)

both<-merge(tt$table, tt2$table, by="row.names") #merge test data from Misty and Roberts (all genes are included, even if not DE)

both$sig<-ifelse((both$FDR.x< 0.05) & (both$FDR.y < 0.05), "both",
                 ifelse(both$FDR.x < 0.05, "Misty", 
                        ifelse(both$FDR.y < 0.05, "Roberts", "none"))) #make new column which indicates if gene is DE for M, R, both, or neither
both$sig<-as.factor(both$sig) #make new column a factor
#plot(both$logFC.x, both$logFC.y, xlab="log-fold change in Misty",
#    ylab="log-fold change in Roberts", col=both$sig) # make plot of log-FC in Misty vs. Roberts, with points coloured by significance
con<-subset(both, both$sig=="both") #get subset where genes are DE in both Misty and Roberts
colnames(con)[1] <- "gene_id"
con$agreement<-ifelse(sign(con$logFC.x)==sign(con$logFC.y), "same","different") #make new column based on whether log-FC is the same sign in M and R 
r<-subset(con, con$agreement=="same") #subset where de in same direction for M and R

anti2<-subset(con, con$agreement=="different") #subset where de in opposite direction for M and R
row.names(anti2)<-anti2$gene_id

#Evaluated the probability of the observed overlap between DE gene in Misty and Robert's
n_rob=length(sig.rob)
n_mis=length(sig.mis)
n_rob_mis=nrow(con)
n_all_method2=nrow(counts)
phyper(n_rob_mis, n_rob, n_all_method2-n_rob, n_mis, lower.tail = F)

############################################
###Make Venn Diagram for Two-model method###
############################################

venn.plot <- venn.diagram(list(sig.mis, sig.rob, diff$Row.names), filename = "venn_method2.png",
                          main="Two-model Method", main.cex=3.5,
                          imagetype = "png", fontface=c(1,1,1,1,2,1,1),
                          height=6000, width=6000,
                          fill=c("gray60","gray86", "white"), 
                          cex = c(4.3, 4.3,4.3,4.3,4.6,4.3,4.3), cat.cex=3,
                          scaled=F, lwd=c(1,1,1), lty=c(1,1,1), col=c("black","black","black"), 
                          cat.pos = c(320,40,180), cat.dist = c(0.08,0.08,0.08),
                          category.names=c("Significant in\nMisty",
                                           " Significant in\n    Robert's",
                                           "Directional\nDifference"))


#################################################################
##Find genes that are significant in both method 1 and method 2##
#################################################################

###genes that are parallel

w<-which(r$gene_id %in% exp.no.int$gene_id) #get indices of genes found in both methods
thegenes<-as.vector(r[w,]$gene_id) #names of genes found above
it<-as.data.frame(t(as.data.frame(cpm(d2)[thegenes,]))) #expression data for top de genes
it$site<-substr(row.names(it),1,2) #add column for site data
boxplot(it$ENSGACG00000013480~it$site, ylab=expression("Log"[2]*"CPM"))

rownames(r)=r[,1]
fc.thegenes<-r[thegenes,] #subset of method 2 data for genes sig. in both methods
cor.test(fc.thegenes$logFC.x, fc.thegenes$logFC.y, method = "pearson") #test of correlation between Misty and Roberts over parallel genes

###genes that are anti-parallel

w.anti<-which(anti1 %in% anti2$gene_id) #get indices of anti-parallel genes found in both methods
anti.genes<-as.vector(anti1[w.anti]) #names of genes found above

rownames(con)=con[,1]
fc.antigenes<-con[anti.genes,] #subset of method 2 data for genes sig. in both methods but anti-parallel
cor.test(fc.antigenes$logFC.x, fc.antigenes$logFC.y, method = "pearson") #test of correlation between Misty and Roberts over parallel genes

##Make Venn Diagram!
venn.diagram(list(exp.no.int$gene_id,r$gene_id), filename="both_methods.png",
             imagetype = "png", fontface=c(1,2,1),
             height=6000, width=6000, alpha=c(.2,.2),
             fill=c("gray50", "gray80"), lty=1,
             cex = c(5,5.5,5), cat.cex=3.9,
             scaled=F, cat.pos=c(190,170),cat.dist=c(0.04,0.04),
             category.names=c("GLM Method", "Two-model Method"))

####################################################################
###log2FC between lake and streem in Misty as compared to Robert's##
###across all genes using PPM correlation                         ##
####################################################################

par(mar=c(5,5,2,2))
mycols=c("black","gray90","gray70","gray50" )
mypch=c(1,16,16,16)
names(mycols) <- levels(both$sig) #gives labels to each color
mycols <- mycols[match(both$sig, names(mycols))] #makes a vector giving a color to each individual
names(mypch) <- levels(both$sig) #gives labels to each color
mypch <- mypch[match(both$sig, names(mypch))] #makes a vector giving a color to each individual
plot(both$logFC.x, both$logFC.y, xlab=expression('Log'[2]*'FC in Misty'),
     ylab=expression("Log"[2]*"FC in Robert's"), col=mycols, pch=mypch) # make plot of log-FC in Misty vs. Roberts, with points coloured by significance
points(fc.thegenes$logFC.x, fc.thegenes$logFC.y, col="black", pch=17, cex=1.5) #add DEP genes
points(fc.antigenes$logFC.x, fc.antigenes$logFC.y, col="black", pch=25,bg="white", cex=1.5) #add DEP genes
abline(lm(both$logFC.x~both$logFC.y))
abline(lm(fc.thegenes$logFC.x~fc.thegenes$logFC.y), lty=5)
abline(lm(fc.antigenes$logFC.x~fc.antigenes$logFC.y), lty=3)
legend(5,-5, c("DE in Misty", "DE in Robert's", "DE in neither", "DEP", "ADEP", "\nDEP in 2-model \nmethod only"),
       col=c("gray90", "gray50", "gray70", "black", "black", "black"), pch=c(16,16,16,17,25,1),
       cex=0.9, pt.cex=1.1, y.intersp = -1, bty="n")

cor.test(both$logFC.x, both$logFC.y, method = "pearson") #test of correlation

##############################################################
#Look at properties of thegenes and antigenes
##############################################################


##first for the parallel genes##

exp.thegenes<-as.data.frame(t(as.data.frame(cpm(d)[thegenes,]))) #get cpm of thegenes
exp.thegenes$site<-substr(row.names(exp.thegenes),1,2) #add a site column
exp.thegenes$ws<-substr(exp.thegenes$site,1,1) #add watershed column
exp.thegenes$hab<-substr(exp.thegenes$site,2,2) #add a habitat column

sum(fc.thegenes$logFC.x<0) #how many genes are upregulated in lake
sum(fc.thegenes$logFC.x>0) # how many genes are upregulated in stream
mean(fc.thegenes$logFC.x[fc.thegenes$logFC.x<0]) #calculate mean misty FC of genes upregulated in lake
mean(fc.thegenes$logFC.x[fc.thegenes$logFC.x>0]) #calculate mean misty FC of genes upregulated in stream
mean(fc.thegenes$logFC.y[fc.thegenes$logFC.x<0]) #calculate mean roberts FC of genes upregulated in lake
mean(fc.thegenes$logFC.y[fc.thegenes$logFC.x>0]) #calculate mean roberts FC of genes upregulated in stream

lake.genes<-t(exp.thegenes[,fc.thegenes$gene_id[fc.thegenes$logFC.x<0]])[,c("ML1","ML2","ML3","RL1","RL2","RL3")] #get lake expression for genes upregulated in lake
lake.mean<-mean(rowMeans(lake.genes)) #average lake expression of genes upreg. in lake
mean(rowMeans(lake.genes[,1:3])) #average ML expression of genes upreg. in lake
mean(rowMeans(lake.genes[,4:6])) #average RL expression of genes upreg. in lake

stream.genes<-t(exp.thegenes[,fc.thegenes$gene_id[fc.thegenes$logFC.x>0]])[,c("MS1","MS2","MS3","RS1","RS2","RS3")] #get stream expression for genes upregulated in stream
stream.mean<-mean(rowMeans(stream.genes)) #average stream expression of genes upreg. in stream
mean(rowMeans(stream.genes[,1:3])) #average MS expression of genes upreg. in stream
mean(rowMeans(stream.genes[,4:6])) #average RS expression of genes upreg. in stream

#Across all filtered genes
misty.genes<-as.data.frame(cpm(d))[,c("ML1","ML2","ML3","MS1","MS2","MS3")] #get expression for parallel genes in Misty
mean(rowMeans(misty.genes[,1:3])) #mean ML expression
mean(rowMeans(misty.genes[,4:6])) #mean MS expression

roberts.genes<-as.data.frame(cpm(d))[,c("RL1","RL2","RL3","RS1","RS2","RS3")] #get expression for parallel genes in Roberts
mean(rowMeans(roberts.genes[,1:3])) #mean RL expression
mean(rowMeans(roberts.genes[,4:6])) #mean RS expression

#################
###GO analysis###
#################

library(topGO)
#get GO id for each gene
library(biomaRt)
mart=useDataset("gaculeatus_gene_ensembl", useMart("ensembl"))
#get all go id for universe 
ensembl_genes=rownames(d) 
gene_pool=getBM(filters = "ensembl_gene_id", 
                attributes = c("ensembl_gene_id","go_id"), 
                values = ensembl_genes, 
                mart = mart)
#get all go id for DEP and ADEP genes
ensembl_genes_DEP=fc.thegenes$gene_id
DEP.gene=getBM(filters = "ensembl_gene_id", 
              attributes = c("ensembl_gene_id","go_id"), 
              values = ensembl_genes_DEP, 
              mart = mart)

ensembl_genes_ADEP=fc.antigenes$gene_id
ADEP.gene=getBM(filters = "ensembl_gene_id", 
               attributes = c("ensembl_gene_id","go_id"), 
               values = ensembl_genes_ADEP, 
               mart = mart)

# Remove blank entries
gene_pool <- gene_pool[gene_pool$go_id != '',]
DEP.gene <- DEP.gene[DEP.gene$go_id != '',]
ADEP.gene <- ADEP.gene[ADEP.gene$go_id != '',]

# convert from table format to list format
geneID2GO <- by(gene_pool$go_id,
                gene_pool$ensembl_gene_id,
                function(x) as.character(x))

myInterestingGenes_DEP=by(DEP.gene$go_id,
                      DEP.gene$ensembl_gene_id,
                      function(x) as.character(x))

myInterestingGenes_ADEP=by(ADEP.gene$go_id,
                          ADEP.gene$ensembl_gene_id,
                          function(x) as.character(x))

# examine result
head(geneID2GO)
head(myInterestingGenes_DEP)
head(myInterestingGenes_ADEP)

correction<-"fdr"
geneNames = names(geneID2GO)
myInterestingGenesNames_DEP=names(myInterestingGenes_DEP)
myInterestingGenesNames_ADEP=names(myInterestingGenes_ADEP)

geneList_DEP = factor(as.integer(geneNames %in% myInterestingGenesNames_DEP))
names(geneList_DEP) <- geneNames

geneList_ADEP = factor(as.integer(geneNames %in% myInterestingGenesNames_ADEP))
names(geneList_ADEP) <- geneNames

# GO terms for DEP
ontology=c("MF","BP","CC")
for (i in 1:length(ontology)) {
  tgData = new("topGOdata", 
               ontology = ontology[i], 
               allGenes = geneList_DEP, 
               annot = annFUN.gene2GO, 
               gene2GO = geneID2GO)
  fisherRes = runTest(tgData, algorithm="classic", statistic="fisher")
  fisherResCor = p.adjust(score(fisherRes), method=correction)
  weightRes = runTest(tgData, algorithm="weight01", statistic="fisher")
  weightResCor = p.adjust(score(weightRes), method=correction)
  allRes = GenTable(tgData, 
                    classic=fisherRes, 
                    weight=weightRes, 
                    orderBy="weight", 
                    ranksOf="classic", 
                    topNodes=200)
  allRes$fisher.COR = fisherResCor[allRes$GO.ID]
  allRes$weight.COR = weightResCor[allRes$GO.ID]
  write.csv(allRes, paste("DEPgenes",ontology[i],"csv",sep="."))
} 


# GO terms for ADEP
ontology=c("MF","BP","CC")
for (i in 1:length(ontology)) {
  tgData = new("topGOdata", 
               ontology = ontology[i], 
               allGenes = geneList_ADEP, 
               annot = annFUN.gene2GO, 
               gene2GO = geneID2GO)
  fisherRes = runTest(tgData, algorithm="classic", statistic="fisher")
  fisherResCor = p.adjust(score(fisherRes), method=correction)
  weightRes = runTest(tgData, algorithm="weight01", statistic="fisher")
  weightResCor = p.adjust(score(weightRes), method=correction)
  allRes = GenTable(tgData, 
                    classic=fisherRes, 
                    weight=weightRes, 
                    orderBy="weight", 
                    ranksOf="classic", 
                    topNodes=200)
  allRes$fisher.COR = fisherResCor[allRes$GO.ID]
  allRes$weight.COR = weightResCor[allRes$GO.ID]
  write.csv(allRes, paste("ADEPgenes",ontology[i],"csv",sep="."))
}