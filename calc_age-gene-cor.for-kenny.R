library("data.table")
library(parallel)
library(ggplot2)
library(DT)
library(d3heatmap)
library(ggplot2)
library(org.Hs.eg.db)
library(topGO)
numcores <- detectCores() - 1 # Find no. cores
clust <- makeCluster(10) # Initiate cluster

# list of healthy samples
hlt_samples <- read.table("race_healthy1-samples.txt")
hlt_samples <- as.vector(hlt_samples$V1)

# gene ids, symbols, names
gene_sym <- read.delim("human_gene-info_ncbi.txt", header=T, row.names=1, sep="\t")
go_terms <- as.list(org.Hs.egGO2ALLEGS)

# sex
x <- 'f'

# age limits
lwr <- 30; upr <- 40; arng <- lwr:upr

# reading-in the entire pcl (both healthy and disease)
gsm_pcl <- data.frame(fread(paste0('blood.', x, '.pcl'), header=T), row.names=1)

# reading-in sample age annotations
gsm_age <- read.table(paste0('blood.', x, '.txt'), row.names=1); colnames(gsm_age) <- c('age','gse')

ggplot(data=gsm_age, aes(gsm_age[,1])) + 
  geom_histogram(binwidth=2,fill=I("black"),col=I("gray")) + 
  ggtitle("Number of GSM samples by age") +
  labs(x="age",y="# of GSM samples") +
  geom_vline(x=lwr, colour="red",size=0.6, linetype="solid") +
  geom_vline(x=upr, colour="red",size=0.6, linetype="solid") + 
  geom_text(dim(gsm_age[which((gsm_age$age>=lwr) & (gsm_age$age<=upr)),])[1], x=80,y=5)

# intersecting col names with healthy sample ids
hlt_samples <- intersect(hlt_samples,
                         intersect(rownames(gsm_age)[which((gsm_age$age>=lwr) & (gsm_age$age<=upr))],
                                   colnames(gsm_pcl)))

# creating pcl for healthy samples
hlt_gsm_pcl <- gsm_pcl[,hlt_samples]
hlt_gsm_age <- gsm_age[hlt_samples,]

# creating pcl of genes per age with median of gene-expression values per gene
hlt_age_pcl <- array(NaN, c(nrow(hlt_gsm_pcl), length(lwr:upr)))
for(j in 1:length(arng)) {
  age_gsm <- rownames(hlt_gsm_age)[which(hlt_gsm_age$age==arng[j])]
  if(length(age_gsm)==1) {
    hlt_age_pcl[,j] <- hlt_gsm_pcl[,age_gsm]
  } else {
    hlt_age_pcl[,j] <- apply(hlt_gsm_pcl[,age_gsm], 1, median)
  }
}

# scaling genes (rows) to mean 0 and variance 1

idx <- apply(hlt_age_pcl,1,var)
hlt_age_pcl_rowz <- t(apply(hlt_age_pcl, 1, scale))
hlt_age_pcl_rowz[idx==0] <- hlt_age_pcl[idx==0]

colnames(hlt_age_pcl_rowz) <- arng
rownames(hlt_age_pcl_rowz) <- rownames(hlt_gsm_pcl)
rownames(hlt_age_pcl) <- rownames(hlt_gsm_pcl)
# boostrap runs
nboot <- 10
boot_rho <- array(NaN, c(nboot, nrow(hlt_gsm_pcl)))
# boot_rho is the output of the boostrap runs
# It contains <nboot> rows and <#genes> columns

for(n in 1:nboot) {
  cat(n, "...\n")
  bag_gsm <- {}
  for(j in 1:length(arng)) {
    age_gsm <- rownames(hlt_gsm_age)[which(hlt_gsm_age$age==arng[j])]
    if(length(age_gsm)==1) {
      bag_gsm <- c(bag_gsm, age_gsm)
    } else {
      sel_gsm <- sample(age_gsm, replace=T)
      bag_gsm <- c(bag_gsm, unique(sel_gsm))
    }
  }
  
  boot_pcl <- hlt_gsm_pcl[,bag_gsm]
  boot_age <- hlt_gsm_age[bag_gsm,]$age
  
  clusterExport(clust,"boot_pcl"); clusterExport(clust,"boot_age")
  rho <- parApply(clust, boot_pcl, 1, cor, y=boot_age, method="spearman")
  boot_rho[n,] <- rho
}
stopCluster(clust)

# processing rho into comparable correlation scores (fisherz)
boot_fisherz <- t(apply(boot_rho, 1, function(x) { scale(0.5*log((1+x)/(1-x))) }))
bxs_boot_fisherz <- apply(boot_fisherz, 2, quantile, probs=c(1, 0.9, 0.75, 0.5, 0.25, 0.1, 0), na.rm=T)
bxs_boot_fisherz <- data.frame(bxs_boot_fisherz); colnames(bxs_boot_fisherz) <- rownames(hlt_gsm_pcl)

# selecting genes that are highly correlated with age
pos_predg <- colnames(bxs_boot_fisherz)[which((bxs_boot_fisherz[4,]>0) & (bxs_boot_fisherz[6,]>=1))]
neg_predg <- colnames(bxs_boot_fisherz)[which((bxs_boot_fisherz[4,]<0) & (bxs_boot_fisherz[2,]<=-1))]
all_predg <- append(pos_predg, neg_predg)
ann <- annFUN.org("BP",feasibleGenes=rownames(gsm_pcl),mapping = "org.Hs.eg.db",ID=c("entrez"))
geneList <- factor(as.integer (rownames(gsm_pcl) %in% all_predg))
names(geneList) <- rownames(gsm_pcl)
selector <- function(x) {return (x==0)}
sampleGOdata <- new("topGOdata",ontology="BP",allGenes=geneList,geneSel=selector,annot = annFUN.GO2genes,GO2genes=ann)
resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")
resultKS <- runTest(sampleGOdata, algorithm = "classic", statistic = "ks")
# resultKS.elim <- runTest(sampleGOdata, algorithm = "elim", statistic = "ks")
allRes <- GenTable(sampleGOdata, classicFisher = resultFisher,
                   classicKS = resultKS, 
#                    elimKS = resultKS.elim,
                   orderBy = "classicKS", ranksOf = "classicFisher", topNodes = 50)
df <- cbind(data.frame(all_predg), t(abs(bxs_boot_fisherz[4,all_predg])))
ggplot(data=df, aes(x=df[,2])) + 
  geom_histogram(binwidth=.1,fill=I("black"),col=I("gray")) +
  labs(x="abs(score)",y="# of genes")

# subsetting the age_pcl for just the 'predictive' genes from the previous step
pos_predg_pcl_rowz <- cbind(gene_sym[pos_predg,],
                            t(bxs_boot_fisherz[4,pos_predg]), apply(hlt_age_pcl_rowz[pos_predg,], 2, sprintf, fmt="%.6f"))
neg_predg_pcl_rowz <- cbind(gene_sym[neg_predg,],
                            t(bxs_boot_fisherz[4,neg_predg]), apply(hlt_age_pcl_rowz[neg_predg,], 2, sprintf, fmt="%.6f"))
colnames(pos_predg_pcl_rowz) <- c("symbol", "name","score", arng)
colnames(neg_predg_pcl_rowz) <- c("symbol", "name","score", arng)
pos_p <- pos_predg_pcl_rowz[order(-pos_predg_pcl_rowz$score,pos_predg_pcl_rowz$name),]
neg_p <- neg_predg_pcl_rowz[order(neg_predg_pcl_rowz$score,neg_predg_pcl_rowz$name),]

rownames(pos_p) <- pos_p[,1]
pos_p <- pos_p[,-c(1,2,3)]
# convert factor to numeric data frame
indx <- sapply(pos_p, is.factor)
pos_p[indx] <- lapply(pos_p[indx], function(x) as.numeric(as.character(x)))
d3heatmap(pos_p[1:50,],Colv = FALSE,xaxis_font_size=7,yaxis_font_size=8)

