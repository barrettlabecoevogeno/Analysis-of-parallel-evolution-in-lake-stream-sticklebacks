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

#Read in and merge individual htseq-count files
counts=readDGE(samples$countf)$counts

#Filter weakly expressed and noninformative features in half of samples
noint = rownames(counts) %in% c("__no_feature","__ambiguous","__too_low_aQual", "__not_aligned")
cpms=cpm(counts)
keep=rowSums(cpms>1) >=3 & !noint 
counts=counts[keep,]

p<-list(1:1000) #make an empty list into which you will put the results of each simulation for dep genes
ap<-list(1:1000) #make an empty list into which you will put the results of each simulation for anti-dep genes
k<-1:1000

for (i in k){

##########################
##### Method 1 edgeR #####
##########################

# Permute habitat for each sample within the same watershed
samples.m=samples[1:6,]
samples.r=samples[7:12,]

samples.m$habitat=sample(samples.m$habitat)
samples.r$habitat=sample(samples.r$habitat)

samples=rbind(samples.m, samples.r)

#Define variables
variable=cbind(samples$habitat, samples$location)

#Define factors
hab=as.factor(variable[,1])
shed=as.factor(variable[,2])

#Normalize all libraries withthe trimmed mean of M-value (TMM)
d=DGEList(counts=counts)
d=calcNormFactors(d)

#Check sacling factors for libraries, any libraries with scaling factor not close to 1 should be excluded from further analysis
d$samples

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
f<-which(DEnames.r %in% DEnames.m) #Get indices of genes found in both
r.m<-DEnames.r[f] #Make list of gene ids for genes indexed above

##Get list of genes sig. for habitat but not for interaction

getrid<-which(r.m %in% DEnames.int) #Genes that are sig. different for habitat and for interaction
no.int<-r.m[-getrid] #Get only names for genes that have no sig. interaction

#Get list of anti-parallel genes
anti1<-r.m[getrid] 

#Get expression levels for all genes sig. for habitat but not for interaction
exp.no.int<-as.data.frame(cpm(d)[no.int,]) 

#Make new column with gene id
exp.no.int$gene_id<-rownames(exp.no.int)

#####################################################################################################################################
##Method 2: do 2 individual pairwise tests (1 for M and 1 for R), then get genes found in both with same sign (as in Ghalambor2015)##
#####################################################################################################################################

d2 <- DGEList(counts = counts, group = samples$location) #read in the count files
colnames(d2) <- c("ML1","ML2","ML3","MS1","MS2","MS3","RL1","RL2", "RL3",
                  "RS1", "RS2", "RS3")

hab2=as.factor(paste(c("M","M","M","M","M","M","R","R","R","R","R","R"), hab, sep = ""))

d2$samples$group<-hab2

d2 <- calcNormFactors(d2) #calculate the normalization factors

d2 <- estimateCommonDisp(d2) #estimate common dispersion
d2<- estimateTagwiseDisp(d2) #estimate individual gene dispersion

de = exactTest(d2, pair=c("MLake","MStream")) #test for de genes between Misty L and S
tt = topTags(de, n=nrow(d2)) #get stats for each test including FDR for Misty 
head(tt$table) #see some de genes
detags2 <- rownames(topTags(de)) #get names of top de genes
cpm(d2)[detags2,] #expression data for top de genes

de2<-exactTest(d2, pair=c("RLake", "RStream")) #test for de genes between Roberts L and S
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

#################################################################
##Find genes that are significant in both method 1 and method 2##
#################################################################

###genes that are parallel

w<-which(r$gene_id %in% exp.no.int$gene_id) #get indices of genes found in both methods
thegenes<-as.vector(r[w,]$gene_id) #names of genes found above
it<-as.data.frame(t(as.data.frame(cpm(d2)[thegenes,]))) #expression data for top de genes
it$site<-substr(row.names(it),1,2) #add column for site data
p[[i]]<-length(thegenes)

###genes that are anti-parallel

w.anti<-which(anti1 %in% anti2$gene_id) #get indices of anti-parallel genes found in both methods
anti.genes<-as.vector(anti1[w.anti]) #names of genes found above
ap[[i]]<-length(anti.genes)
}

v<-unlist(p)
write.csv(v, file="null.p.csv", row.names=F)
hist(v, breaks=90, xlab="Number of Genes", main="", xlim = c(0,30))

x=unlist(ap)
write.csv(x, file="null.ap.csv", row.names=F)
hist(x, breaks=90, xlab="Number of Genes", main="", xlim = c(0,30))
