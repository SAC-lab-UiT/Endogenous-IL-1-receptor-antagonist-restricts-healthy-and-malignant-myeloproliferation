pkgs <- c('clusterProfiler','data.table','DESeq2','GEOquery','ggplot2','ggpubr','ggrepel','ggsci','ggsignif','ggthemes','grid','gridExtra','limma','org.Mm.eg.db','plotly','ReactomePA','splitstackshape','stringr','survival','survminer')
for(i in 1:length(pkgs)) {
  library(pkgs[i], character.only = T)
}

library(this.path)

#Change here the path to the analysis folder
analysisPath <- dirname(this.path())
setwd(analysisPath)

########################################################################################
###        Survival AML patients Verhaak dataset (Fig. 1a)                            ##
###                                                                                   ##
########################################################################################

survOrig <- read.table(file = 'survData.txt',sep = '\t',head = T,stringsAsFactors = F)
rownames(survOrig) <-  survOrig[,5]


library(GEOquery)
library(survminer)
library(gridExtra)
library(survival)
gse2 <- getGEO("GSE14468")
expr_mat <- exprs(gse2$GSE14468_series_matrix.txt.gz)
probes <-  c("216243_s_at", "205067_at")
il1rnProbes <- c('212657_s_at', '212659_s_at','216243_s_at','216244_at')
M <- expr_mat[il1rnProbes,]

rownames(survOrig) <- survOrig$geo_accession

modelFrame <- data.frame(OS = survOrig$OS_days,dead = survOrig$OS_bin,
                         age = survOrig$Age,risk = survOrig$risk,row.names = rownames(survOrig))
modelFrame$FAB <- rep('')
modelFrame$FAB[survOrig$FAB %in% c('M0','M1','M2','M3')] <- 'M0M3'
modelFrame$FAB[survOrig$FAB %in% c('M4','M5')] <- 'M4M5'
modelFrame$IL1RN <- colMeans(M[il1rnProbes[1:3],rownames(survOrig)])
keep <- which(modelFrame$FAB %in% c('M0M3', 'M4M5') &  modelFrame$risk != 'unknown')
keepM0M3 <- which(modelFrame$FAB %in% c('M0M3') & modelFrame$risk != 'unknown')
keepM4M5 <- which(modelFrame$FAB %in% c('M4M5') & modelFrame$risk != 'unknown')


####Global Model 

x = as.formula("Surv(OS, dead) ~ IL1RN + risk + age")
xNull <- as.formula("Surv(OS, dead) ~ risk + age")
resNull <- suppressWarnings(coxph(xNull, data = modelFrame[keep,]))
res.cox <- suppressWarnings(coxph(x, data = modelFrame[keep,]))
x <- summary(res.cox)
anovaRes <- anova(resNull,res.cox)

survTest  <- function(modelFrame,cutoffs = 5:95){

    x = as.formula("Surv(OS, dead) ~ il1rnFactor + risk + age")
    xNull <- as.formula("Surv(OS, dead) ~ risk + age")	

    pVals <- matrix(0,length(cutoffs),3)
    colnames(pVals) <- c('pVal','expCoef','')
    il1rnQuants  <- quantile(modelFrame$IL1RN,seq(0.01,1,by = 0.01))
    
    for (ct in 1:length(cutoffs)){
        il1rnFactor <- cut(modelFrame$IL1RN,breaks = c(0,il1rnQuants[cutoffs[ct]],1000),
            labels = c('low','high'))
        modelFrame$il1rnFactor <- il1rnFactor    
        
	resNull <- suppressWarnings(coxph(xNull, data = modelFrame))
	res.cox <- suppressWarnings(coxph(x, data = modelFrame))
	anovaRes <- anova(resNull,res.cox)
	pVals[ct,1] <- anovaRes$'P(>|Chi|)'[2]
	pVals[ct,2] <-  exp(res.cox$coefficients)[1]
    }
    return(list(pVals = pVals,df = modelFrame))
}        

fullMod <- survTest(modelFrame[keep,])$pVals
M3Mod <- survTest(modelFrame[keepM0M3,])$pVals
M5Mod <- survTest(modelFrame[keepM4M5,])$pVals


plot(seq(0,1,by = 0.01)[5:95],fullModel[,1])

pdf(file = 'stability.pdf')
layout(matrix(c(1,2,3),3,1))
plot(seq(0,1,by = 0.01)[5:95],fullMod[,1],type = 'l', col = 'red',
         xlim = c(0,1),ylim = c(0,2),main = 'All samples',xlab = 'quantile cutoff',
         ylab = 'p value [red]   exp(B)[Blue]')
plot.xy(xy.coords(x = seq(0,1,by = 0.01)[5:95],y = fullMod[,2]),
        type = 'l', col = 'blue')
abline(h = 1)
abline(h = 0.05,lty  = 2)

plot(seq(0,1,by = 0.01)[5:95],M3Mod[,1],type = 'l', col = 'red',
         xlim = c(0,1),ylim = c(0,2),main = 'M0M3',xlab = 'quantile cutoff',
        ylab = 'p value [red]   exp(B)[Blue]')
plot.xy(xy.coords(x = seq(0,1,by = 0.01)[5:95],y = M3Mod[,2]),
        type = 'l', col = 'blue')
abline(h = 1)
abline(h = 0.05,lty  = 2)

plot(seq(0,1,by = 0.01)[5:95],M5Mod[,1],type = 'l', col = 'red',
         xlim = c(0,1),ylim = c(0,2),main = 'M4M5',xlab = 'quantile cutoff',
         ylab = 'p value [red]   exp(B)[Blue]')
plot.xy(xy.coords(x = seq(0,1,by = 0.01)[5:95],y = M5Mod[,2]),
        type = 'l', col = 'blue')
abline(h = 1)
abline(h = 0.05,lty  = 2)


##Make plots for paper

##Global model

x = as.formula("Surv(OS, dead) ~ il1rnFactor + risk + age")
xNull <- as.formula("Surv(OS, dead) ~ risk + age")	
il1rnQuants  <- quantile(modelFrame$IL1RN[keep],seq(0.01,1,by = 0.01))
ct <- il1rnQuants[91]   
il1rnFactor <- cut(modelFrame$IL1RN,breaks = c(0,ct,1000),
            labels = c('low','high'))
modelFrame$il1rnFactor <- il1rnFactor
modelFrame$surv <- modelFrame$OS/365.25

resNull <- suppressWarnings(coxph(xNull, data = modelFrame[keep,]))
res.cox <- suppressWarnings(coxph(x, data = modelFrame[keep,]))
anovaRes <- anova(resNull,res.cox)
	
fit1 <- survfit(Surv(surv, dead) ~ il1rnFactor,
               data = modelFrame[keep,])

p1 <- ggsurvplot(fit1, data = modelFrame[keep,],palette = c('red','blue'))
p1 <- p1+ggtitle('Global')

globalPlotData <- summary(fit1)
globalPlotData <- data.frame(time = globalPlotData$time,
                            group = globalPlotData$strata,surv = globalPlotData$surv)


##M4M5 Model
resNull <- suppressWarnings(coxph(xNull, data = modelFrame[keepM4M5,]))
res.cox <- suppressWarnings(coxph(x, data = modelFrame[keepM4M5,]))
anovaRes <- anova(resNull,res.cox)
	
fit2 <- survfit(Surv(surv, dead) ~ il1rnFactor,
               data = modelFrame[keepM4M5,])

p2 <- ggsurvplot(fit2, data = modelFrame[keepM4M5,],palette = c('red','blue'))
p2 <- p2+ggtitle('M4M5 Only')


myplots_out <- arrange_ggsurvplots(
  list(M0M3 = p1,M4M5 = p2),
  print = FALSE,
  ncol = 1,
  nrow = 2,
  title = "")

dir.create("output")

ggsave(
  myplots_out,
  file = "output/SurvivalGlobal2.pdf",
  width = 5,
  height = 7)

M4M5plotData <- summary(fit2)                           
M4M5plotData <- data.frame(time = M4M5plotData$time,group = M4M5plotData$strata,surv = M4M5plotData$surv)


write.table(globalPlotData, file = "output/globalCurve.txt",sep = "\t")
write.table(M4M5plotData, file = "output/M4M5Curve.txt",sep = "\t")

#######################################################################################
############                  19 matched-pair AML patients                 ############
############                      Survival   (Fig. S1B)                    ############
############                                                               ############
#######################################################################################

GSE83533_exprs <- data.table::fread(file = "GSE83533_Dx_Rel_AML_RNAseq_rpkm.txt", sep = '\t')

IL1RN_exprs <- GSE83533_exprs[GSE83533_exprs$gene == "IL1RN", ]
IL1RN_exprs <- IL1RN_exprs[, 3:ncol(IL1RN_exprs)]
IL1RN_exprs <- as.data.frame(t(IL1RN_exprs))
IL1RN_exprs$sample <- row.names(IL1RN_exprs)
IL1RN_exprs$type <- ifelse(grepl("Rel", IL1RN_exprs$sample), "Relapse", "Diagnosis")
IL1RN_exprs$patient <- gsub("AML_|_Dx|_Rel", "", IL1RN_exprs$sample)
IL1RN_exprs <- IL1RN_exprs[order(IL1RN_exprs$patient), ]

write.table(IL1RN_exprs, file = 'IL1RN_expression_Diagnosis-Relapse.csv', sep = '\t')

# SraRunTable from study phs001027 using SRA Run Selector (https://www.ncbi.nlm.nih.gov/Traces/study/?) 
# https://www.ncbi.nlm.nih.gov/Traces/study/?acc=phs001027&f=analyte_type_sam_s%3An%3Arna%3Ac&s=SRR3088324,SRR3088328,SRR3088346,SRR3088348,SRR3088350,SRR3088354,SRR3088357,SRR3088361,SRR3088391,SRR3088395,SRR3088398,SRR3088400,SRR3088423,SRR3088427,SRR3088430,SRR3088432,SRR3088444,SRR3088448,SRR3088456,SRR3088457,SRR3088458,SRR3088459,SRR3088460,SRR3088461,SRR3088462,SRR3088463,SRR3088464,SRR3088465,SRR3088466,SRR3088467,SRR3088468,SRR3088469,SRR3088470,SRR3088471,SRR3088472,SRR3088473,SRR3088474,SRR3088475,SRR3088476,SRR3088477,SRR3088478,SRR3088479,SRR3088480,SRR3088481,SRR3088482,SRR3088483,SRR3088484,SRR3088485,SRR3088486,SRR3088487,SRR3088488,SRR3088489,SRR3088490,SRR3088491,SRR3088492,SRR3088493,SRR3088494,SRR3088495,SRR3088496,SRR3088497,SRR3088498,SRR3088499,SRR3088500,SRR3088501,SRR3088502,SRR3088503,SRR3088504,SRR3088505,SRR3088506,SRR3088507,SRR3088508,SRR3088509,SRR3088510,SRR3088511,SRR3088512,SRR3088513,SRR3088514,SRR3088515,SRR3088516,SRR3088517,SRR3088518,SRR3088519,SRR3088520,SRR3088521,SRR3088522,SRR3088523,SRR3088524,SRR3088525,SRR3088526,SRR3088527,SRR3088528,SRR3088529,SRR3088530,SRR3088531,SRR3088532,SRR3088533,SRR3088534,SRR3088535,SRR3088536,SRR3088537,SRR3088538,SRR3088539,SRR3088540,SRR3088541,SRR3088542,SRR3088543,SRR3088544,SRR3088545,SRR3088546,SRR3088547,SRR3088548,SRR3088549,SRR3088550,SRR3088551,SRR3088564,SRR3088568,SRR3088584,SRR3088588,SRR3088034,SRR3088036,SRR3088058,SRR3088060,SRR3088086,SRR3088088,SRR3088096,SRR3088098,SRR3088130,SRR3088132,SRR3088136,SRR3088138,SRR3088154,SRR3088156,SRR3088166,SRR3088168
# For simplicity, SraRunTable.txt is also included in repository. 
# This file needs to be download to the current working directory before commencing the next command
SraRunTable = fread("SraRunTable.txt", sep = ",")
pheno_subject <- fread("phs001027.v2.pht005216.v1.p1.c1.Epigenetics_AML_Subject_Phenotypes.GRU-PUB.txt", skip = "dbGaP_Subject_ID")

samples_RNAseq = SraRunTable[grep("_RNAseq", SraRunTable$`Sample Name`), ]
samples <- merge(pheno_subject, samples_RNAseq, by.x = "SUBJECT_ID", by.y = "submitted_subject_id", all.y = T, all.x = F)
samples = samples[grep("_scRNAseq", samples$`Sample Name`, invert = T), ]
samples <- samples[grep("N/A", samples$TTR, invert = T), ]
samples$sample <- gsub("_RNAseq", "", samples$`Sample Name`)

samples <- merge(IL1RN_exprs, samples, by = "sample")
samples$TTR <- as.numeric(samples$TTR)
samples$bin <- ifelse(samples$type == "Relapse", 1, 0)

samples$Expression <- factor(ifelse(samples$V1 > mean(samples$V1), "HIGH", "LOW"))
fit <- survfit(Surv(TTR, bin) ~ Expression, data = samples, na.action = na.omit)
res.cox <- suppressWarnings(coxph(Surv(TTR, bin) ~ log2(V1+1) + strata(samples$patient), data = samples, na.action = na.omit))


p <- ggsurvplot(fit, combine = F, censor = FALSE, pval = T, palette = "npg", break.x.by = 365, xscale = "d_y", surv.median.line = "hv")
p <- p + xlab("Years") + ylab("Relapse-free probability") 
p <- p$plot
cols = rev(pal_npg()(2))
#names(cols)
labels = c("High IL1RN expression", "Low IL1B expression")
p <- p + scale_color_manual(values = cols, labels = labels) + theme(legend.background = element_blank())
p <- p + theme(legend.direction = "vertical", legend.position = c(0.7,0.9), legend.title = element_blank())
p <- p + annotate("text", x = (5*365)+100, y = 0.69, label = "69%", size  = 6)
p <- p + annotate("text", x = (5*365)+100, y = 0.1, label = "10%", size  = 6)

pdf(file = 'Relapse-free_survival.pdf')
p
dev.off()

p$data

#value row 32 must be discarded because upper/lower std.error = NA

write.table(p$data,file = 'Relapse-free_probability_IL1RN_expression.csv', sep = '\t')


########################################################################################
###  NFkB activity in hematopoetic and progenitor cells (HSPCs) (Fig. 3b, Fig. 8g,    ##
###  Suppl. Fig. S14G)                                                                ##
########################################################################################
library(this.path)

#Change here the path to the analysis folder
analysisPath <- dirname(this.path())
setwd(analysisPath)
library(DESeq2)
library(org.Mm.eg.db)
library(stringr)
library(ggplot2)
library(ggpubr)


nfkbTargets <- read.table(file = 'List_NFKB_gene_Arranz.txt',stringsAsFactors = F,sep = '\t')
########################
###      Il1rnKO      ##
########################


#read table 
m2   <- read.table(file = 'data_Il1rnko.txt', row.names=1, header=T, sep="\t")
temp <- m2[,grep("_rawCount", colnames(m2))]
keep <- which(rowMeans(temp) > 1) 
raw  <- temp[keep,]
keep <- which(rowMeans(log2(raw+1))>3)
raw  <- raw[keep,]

#annotation 
m2    <- m2[rownames(raw),]
index <- setNames(m2$Description, rownames(m2))
gid   <- setNames(gsub("GeneID:", "", str_extract(m2$Description, "GeneID:[0-9]+")), rownames(m2))
annot <- cbind(as.character(gid[rownames(m2)]), as.character(index[rownames(m2)]))
rownames(annot) <- rownames(m2) 
xx    <- bitr(as.character(gid), "ENTREZID", "SYMBOL",org.Mm.eg.db)
xx    <- setNames(xx[,2], xx[,1])
annot <- data.frame(cbind(annot, as.character(xx[annot[,1]])))
colnames(annot) <- c("ENTREZID", "DESCRIPTION", "SYMBOL")
annot <- annot[, c(1,3,2)]
annot$SYMBOL      <- as.character(annot$SYMBOL)
annot$DESCRIPTION <- as.character(annot$DESCRIPTION)
annot$ENTREZID    <- as.character(annot$ENTREZID)

pheno         <- read.table(file = 'pheno_Il1rn.csv', row.names=1, header=T, sep="\t")
cellLine <- c('LT','LT','LT','LT','LT','LT','ST','ST','ST','ST','ST','ST',
              'MPP','MPP','MPP','MPP','MPP','MPP','CD63','CD63','CD63','CD63','CD63','CD63')
pheno <- cbind(pheno,cellLine)	             
colnames(raw) <- gsub("_rawCount", "", colnames(raw))
raw            <- raw[, rownames(pheno)]

dds <- DESeqDataSetFromMatrix(countData = raw, colData=pheno, design = ~ Comp)

vsd <- vst(dds, blind = TRUE)
M <- vsd@assays@data[[1]]

pheno$treatment <- pheno$Sample
pheno$treatment[pheno$treatment == 'control'] <- 'WT'
pheno$treatment[pheno$treatment == 'IL-1rn KO'] <- 'KO'
pheno$Group <- pheno$cellLine

###Declare function before use


estimateTFactivity <- function(M,targetGenes,pheno,annot,
                               groups = c('LT','ST','MPP'),
                               treatments = c('WT','KO'),n = 3, var.equal = T){
  ##X is a matrix of gene expression with gene IDs on the rowsnames
  ##targetGenes are NFKB target genes 
  ##wt are indexes for the wild type samples
  
  resTable <- matrix(0,nrow = length(groups),(n*2+3))
  rownames(resTable) <- groups
  repPheno <- NULL
  
  for (ct in groups){
    ctInd <- which(pheno$Group == ct)
    trInd <- which(pheno$treatment[ctInd] == treatments[2])
    ctrlInd <- which(pheno$treatment[ctInd] == treatments[1])
    Mtf <- M[annot[,2] %in% targetGenes[,2],ctInd]
    repPheno <- rbind(repPheno,pheno[ctInd,][trInd,],pheno[ctInd,][ctrlInd,])
    
    pcaRes <- prcomp(t(Mtf))
    testRes <- t.test(pcaRes$x[trInd,1],pcaRes$x[ctrlInd,1])
    resTable[ct,] <- c(pcaRes$x[trInd,1],pcaRes$x[ctrlInd,1],testRes$estimate,testRes$p.value)
  }
  
  resTable[,1:6] <- resTable[,1:6]-min(resTable[,1:6])
  exportTable <- cbind(repPheno,as.vector(t(resTable[, 1:6])))
  colnames(exportTable)[ncol(exportTable)] <- 'NFKB_activity'
  
  exportLevels <- as.vector(t(outer(groups,treatments,FUN = "paste",sep = '_')))
  
  exportTable$treatmentGroup <- factor(paste(exportTable$Group,exportTable$treatment,sep = '_'),
                                       levels = exportLevels)
  cp <- ggplot(exportTable, aes(x=treatmentGroup, y=NFKB_activity)) + stat_compare_means(comparisons = list(A = c(1,2),B = c(3,4),C = c(5,6)),
                                                                                         method = 't.test',method.args = list(var.equal = var.equal)) +
    geom_boxplot(fill = rep(c('blue','red'),3))
  return(list(test = resTable,exportTable = exportTable,cp = cp))
}        


out <- estimateTFactivity(M,targetGenes = nfkbTargets,pheno = pheno,annot = annot,
                          groups = c('LT','ST','MPP'),
                          treatments = c('WT','KO'))

pdf(file = 'NFKB_Il1rnKO.pdf')
out$cp
dev.off()


write.table(out$exportTable,file = 'NFKB_activity_Il1rn-KO.csv', sep = '\t')

########################
###      Nras         ##
########################


library(clusterProfiler) 
library(ReactomePA)
library(org.Mm.eg.db)
library(stringr)
library(DESeq2)

nfkbTargets <- read.table(file = 'List_NFKB_gene_Arranz.txt',stringsAsFactors = F,sep = '\t')
#read table 
m2   <- read.table('raw_mat_Nras.csv', row.names=1, header=T, sep="\t")
temp <- m2[,grep("_rawCount", colnames(m2))]
keep <- which(rowMeans(temp) > 1) 
raw  <- temp[keep,]
#how few reads should we consider 
#looks good  boxplot(log2(raw+1))
keep <- which(rowMeans(log2(raw+1))>3)
raw  <- raw[keep,]

#annotation 
m2    <- m2[rownames(raw),]
index <- setNames(m2$Description, rownames(m2))
gid   <- setNames(gsub("GeneID:", "", str_extract(m2$Description, "GeneID:[0-9]+")), rownames(m2))
annot <- cbind(as.character(gid[rownames(m2)]), as.character(index[rownames(m2)]))
rownames(annot) <- rownames(m2) 
xx    <- bitr(as.character(gid), "ENTREZID", "SYMBOL",org.Mm.eg.db)
xx    <- setNames(xx[,2], xx[,1])
annot <- data.frame(cbind(annot, as.character(xx[annot[,1]])))
colnames(annot) <- c("ENTREZID", "DESCRIPTION", "SYMBOL")
annot <- annot[, c(1,3,2)]
annot$SYMBOL      <- as.character(annot$SYMBOL)
annot$DESCRIPTION <- as.character(annot$DESCRIPTION)
annot$ENTREZID    <- as.character(annot$ENTREZID)

#LT -long-term hematopoietic stem cells
#ST - short-term hematopoietic stem cells
#MPP - multipotent progenitors
#CD63 - mesenchymal stromal cells
#For each type of cell e.g. LT; compare control - KO for differential expression
#CD63 - volcano with differentially expressed

#start DESeq2 
pheno <- read.table('pheno_Nras.csv', header=T, sep="\t")
colnames(raw) <- gsub("_rawCount", "", colnames(raw))
rownames(pheno) <- gsub("_rawCount", "", pheno$ID)
raw            <- raw[, rownames(pheno)]
library(DESeq2)
##temp <- raw[,-grep("CD63", colnames(raw))]
dds  <- DESeqDataSetFromMatrix(countData = raw, colData=pheno[colnames(raw),], design = ~ Group+Type)
##dds <- DESeqDataSetFromMatrix(countData = raw, colData=pheno, design = ~ Comp)

vsd <- vst(dds, blind = TRUE)
M <- vsd@assays@data[[1]]

pheno$treatment <- pheno$Type
pheno$treatment[pheno$treatment == 'trans'] <- 'Nras'
pheno$treatment[pheno$treatment == 'control'] <- 'Ctrl'

estimateTFactivity <- function(M,targetGenes,pheno,annot,
                               groups = c('LT','ST','MPP'),
                               treatments = c('Ctrl','Nras'),n = 3, var.equal = T){
  ##X is a matrix of gene expression with gene IDs on the rowsnames
  ##targetGenes are NFKB target genes 
  ##wt are indexes for the wild type samples
  
  resTable <- matrix(0,nrow = length(groups),(n*2+3))
  rownames(resTable) <- groups
  repPheno <- NULL
  
  for (ct in groups){
    ctInd <- which(pheno$Group == ct)
    trInd <- which(pheno$treatment[ctInd] == treatments[2])
    ctrlInd <- which(pheno$treatment[ctInd] == treatments[1])
    Mtf <- M[annot[,2] %in% targetGenes[,2],ctInd]
    repPheno <- rbind(repPheno,pheno[ctInd,][trInd,],pheno[ctInd,][ctrlInd,])
    
    pcaRes <- prcomp(t(Mtf))
    testRes <- t.test(pcaRes$x[trInd,1],pcaRes$x[ctrlInd,1])
    resTable[ct,] <- c(pcaRes$x[trInd,1],pcaRes$x[ctrlInd,1],testRes$estimate,testRes$p.value)
  }
  
  resTable[,1:6] <- resTable[,1:6]-min(resTable[,1:6])
  exportTable <- cbind(repPheno,as.vector(t(resTable[, 1:6])))
  colnames(exportTable)[ncol(exportTable)] <- 'NFKB_activity'
  
  exportLevels <- as.vector(t(outer(groups,treatments,FUN = "paste",sep = '_')))
  
  exportTable$treatmentGroup <- factor(paste(exportTable$Group,exportTable$treatment,sep = '_'),
                                       levels = exportLevels)
  cp <- ggplot(exportTable, aes(x=treatmentGroup, y=NFKB_activity)) + stat_compare_means(comparisons = list(A = c(1,2),B = c(3,4),C = c(5,6)),
                                                                                         method = 't.test',method.args = list(var.equal = var.equal)) +
    geom_boxplot(fill = rep(c('blue','red'),3))
  return(list(test = resTable,exportTable = exportTable,cp = cp))
}        



outRas <- estimateTFactivity(M,targetGenes = nfkbTargets,pheno = pheno,annot = annot,
                             groups = c('LT','ST','MPP'),
                             treatments = c('Ctrl','Nras'))

pdf(file = 'Nfkb_Nras.pdf')
outRas$cp
dev.off()

write.table(outRas$exportTable,file = 'NFKB_activity_Nras.csv', sep = '\t')


########################
###      Vav-Cre      ##
########################

M <- read.table(file = 'm2.txt')
foo <- read.table(file = 'foo.txt',stringsAsFactors = F)
supp <- read.table(file = 'foo.txt',stringsAsFactors = F)
pheno <- read.table(file = 'pheno.txt',stringsAsFactors = F)

supp <- read.table(file = 'List_NFKB_gene_Arranz.txt',stringsAsFactors = F)

use <- foo[foo[,2] %in% supp[,2],1]

pcaRes <- prcomp(t(M[use,]))

x <- pcaRes$x[,1]
#x <- colMeans(M[use,])
x <- -1*x
x <- x-min(x)
#

meanNras <- mean(x[1:3])
seNras <- sd(x[1:3])/sqrt(3)
meanWT <- mean(x[4:6])
seWT <- sd(x[4:6])/sqrt(3)

res <- t.test(x[1:3], x[4:6], var.equal = T)
round(res$p.value, digits = 4)

pdf(file = 'Nfkb_VAV_CRE.pdf',width = 3,height = 4) 
plot(c(0.5,1.5),c(meanWT,meanNras),col = c('blue','red'),xlim = c(0,2),
     pch = 19,xaxt = 'n', ylim =c(0,30 ),xlab = '',ylab = 'NFkB -activity')
plot.xy(xy.coords(x = c(0.5,0.5),y = c((meanWT-seWT),(meanWT+seWT))),
        type = 'l',lwd = 3,col = 'blue')
plot.xy(xy.coords(x = c(1.5,1.5),y = c((meanNras-seNras),(meanNras+seNras))),
        type = 'l',lwd = 3,col = 'red')
axis(1,at = 1.5, 'Vav-Cre')
axis(1,at = 0.5, 'WT')


plot.xy(xy.coords(x = c(0.5,0.5),y = c(22,25)),
        type = 'l')

plot.xy(xy.coords(x = c(1.5,1.5),y = c(22,25)),
        type = 'l')

plot.xy(xy.coords(x = c(0.5,1.5),y = c(25,25)),
        type = 'l')

text(c(1),27,label = round(res$p.value, digits = 4))

dev.off()

resTable <- t(x)
rownames(resTable) <- 'NFkB activity'
write.table(resTable, file = 'Nfkb_plot_Data_VAV_CRE.txt',sep = '\t')

    
##########################################################################################################################################
library(ggplot2)
library(tidyverse)
library(ggvenn)
library(ggrepel)
library(openxlsx)
library(patchwork)
library(RColorBrewer)
library(this.path)
library(edgeR)
library(DESeq2)
library(biomaRt)

color.list <- c(brewer.pal(12, "Paired"),brewer.pal(12, "Set3"),brewer.pal(8, "Pastel2"),colorRampPalette(c("grey20","grey70"))(4))
to <- color.list[22]
ss <- color.list[18]
color.list[22] <- ss
color.list[18] <- to

setwd(dirname(this.path()))

objects <- c("LT_Il1rnKO", "ST_Il1rnKO", "MPP_Il1rnKO")

#Analysis of RNA-sequencing from Il1rn-KO raw counts
GSE126428.rawCounts <- read.delim("GSE126428_eurofins_rawCount.tab",sep = "\t",header = T)
GSE126428.rawCounts$GeneName <- gsub("\\(.*","",GSE126428.rawCounts$Description)
GSE126428.rawCounts$Description <- NULL
head(GSE126428.rawCounts)

GSE126428.rawCounts <- GSE126428.rawCounts %>% dplyr::select(GeneName,ENTREZID,everything()) %>% dplyr::select(-Length,-Effective_Length)

GSE126428.aggCounts <- GSE126428.rawCounts %>% dplyr::select(-ENTREZID)%>% group_by(GeneName) %>% summarise(across(everything(), sum))


#######   #LT Analysis   ############

# Select LT experiment
GSE126428.LT.Counts <- GSE126428.aggCounts %>% dplyr::select(GeneName,contains("_LT"))

# Filter genes: at least 1 cpm in at least 3 samples
idx <- which(rowSums(cpm(GSE126428.LT.Counts[,-c(1,2)]) >= 1) >= 3)

GSE126428.LT.FilterCounts <- GSE126428.LT.Counts[idx,]

GSE126428.LT.FilterCounts <- as.data.frame(GSE126428.LT.FilterCounts)
rownames(GSE126428.LT.FilterCounts) <- GSE126428.LT.FilterCounts$GeneName
GSE126428.LT.FilterCounts$GeneName <- NULL

GSE126428.LT.coldata <- data.frame(Sample=colnames(GSE126428.LT.FilterCounts),
                                   Condition=factor(c("IL1KO","IL1KO","IL1KO","Ctrl","Ctrl","Ctrl")))

rownames(GSE126428.LT.coldata) <- GSE126428.LT.coldata$Sample

GSE126428.LT.dds <- DESeqDataSetFromMatrix(countData = GSE126428.LT.FilterCounts,
                                           colData = GSE126428.LT.coldata,
                                           design= ~ Condition)

GSE126428.LT.dds$Condition <- relevel(GSE126428.LT.dds$Condition, ref = "Ctrl")

GSE126428.LT.dds <- DESeq(GSE126428.LT.dds)

plotDispEsts(GSE126428.LT.dds)

GSE126428.LT.vsd <- vst(GSE126428.LT.dds, blind=FALSE)

plotPCA(GSE126428.LT.vsd,intgroup=c("Condition"))

GSE126428.LT.res <- results(GSE126428.LT.dds,alpha = 0.05)
GSE126428.LT.res <- GSE126428.LT.res[order(GSE126428.LT.res$pvalue),]

summary(GSE126428.LT.res)

GSE126428.LT.res <- as.data.frame(GSE126428.LT.res) %>% mutate(padj2 = p.adjust(pvalue, method = "fdr"))
GSE126428.LT.res$external_gene_name <- rownames(GSE126428.LT.res)

GSE126428.LT.res <- GSE126428.LT.res %>% mutate(score=sqrt(log2FoldChange^2 + (-log10(padj2))^2))

counts_norm <- counts(GSE126428.LT.dds, normalized = TRUE) + 0.5

counts_norm <- counts_norm %>%
  as.data.frame() %>%
  mutate(external_gene_name = rownames(counts_norm)) %>%
  dplyr::select(external_gene_name,everything())

GSE126428.LT.res %>%
  left_join(counts_norm
            , by = c("external_gene_name")) %>%
  write.xlsx(file = "LT_Il1KO_DESeq2_New.xlsx",overwrite = T)

assign(objects[1], GSE126428.LT.res)


############  ST Analysis  ##################

# Select LT experiment
GSE126428.ST.aggCounts <- GSE126428.aggCounts %>% dplyr::select(GeneName,contains("_ST"))

# Filter genes: at least 1 cpm in at least 3 samples
idx <- which(rowSums(cpm(GSE126428.ST.aggCounts[,-c(1,2)]) >= 1) >= 3)
GSE126428.ST.FilterCounts <- GSE126428.ST.aggCounts[idx,]

GSE126428.ST.FilterCounts <- as.data.frame(GSE126428.ST.FilterCounts)
rownames(GSE126428.ST.FilterCounts) <- GSE126428.ST.FilterCounts$GeneName

GSE126428.ST.FilterCounts$GeneName <- NULL

GSE126428.ST.coldata <- data.frame(Sample=colnames(GSE126428.ST.FilterCounts),
                                   Condition=factor(c("IL1KO","IL1KO","IL1KO","Ctrl","Ctrl","Ctrl")))

rownames(GSE126428.ST.coldata) <- GSE126428.ST.coldata$Sample

GSE126428.ST.dds <- DESeqDataSetFromMatrix(countData = GSE126428.ST.FilterCounts,
                                           colData = GSE126428.ST.coldata,
                                           design= ~ Condition)

GSE126428.ST.dds$Condition <- relevel(GSE126428.ST.dds$Condition, ref = "Ctrl")

GSE126428.ST.dds <- DESeq(GSE126428.ST.dds)

plotDispEsts(GSE126428.ST.dds)

GSE126428.ST.vsd <- vst(GSE126428.ST.dds, blind=FALSE)

plotPCA(GSE126428.ST.vsd,intgroup=c("Condition"))

GSE126428.ST.res <- results(GSE126428.ST.dds,alpha = 0.05)
GSE126428.ST.res <- GSE126428.ST.res[order(GSE126428.ST.res$pvalue),]

summary(GSE126428.ST.res)

GSE126428.ST.res <- as.data.frame(GSE126428.ST.res) %>% mutate(padj2 = p.adjust(pvalue, method = "fdr"))
GSE126428.ST.res$external_gene_name <- rownames(GSE126428.ST.res)

GSE126428.ST.res <- GSE126428.ST.res %>% mutate(score=sqrt(log2FoldChange^2 + (-log10(padj2))^2))

counts_norm <- counts(GSE126428.ST.dds, normalized = TRUE) + 0.5

counts_norm <- counts_norm %>%
  as.data.frame() %>%
  mutate(external_gene_name = rownames(counts_norm)) %>%
  dplyr::select(external_gene_name,everything())

GSE126428.ST.res %>%
  left_join(counts_norm
            , by = c("external_gene_name")) %>%
  write.xlsx(file = "ST_Il1KO_DESeq2_New.xlsx",overwrite = T)

assign(objects[2], GSE126428.ST.res)


#########   MPP Analysis   #############

# Select MPP Experiment Samples
GSE126428.MPP.aggCounts <- GSE126428.aggCounts %>% dplyr::select(GeneName,contains("_MPP"))

# Filter genes: at least 1 cpm in at least 3 samples
idx <- which(rowSums(cpm(GSE126428.MPP.aggCounts[,-c(1,2)]) >= 1) >= 3)
GSE126428.MPP.FilterCounts <- GSE126428.MPP.aggCounts[idx,]

GSE126428.MPP.FilterCounts <- as.data.frame(GSE126428.MPP.FilterCounts)
rownames(GSE126428.MPP.FilterCounts) <- GSE126428.MPP.FilterCounts$GeneName

GSE126428.MPP.FilterCounts$GeneName <- NULL

GSE126428.MPP.coldata <- data.frame(Sample=colnames(GSE126428.MPP.FilterCounts),
                                    Condition=factor(c("IL1KO","IL1KO","IL1KO","Ctrl","Ctrl","Ctrl")))

rownames(GSE126428.MPP.coldata) <- GSE126428.MPP.coldata$Sample

GSE126428.MPP.dds <- DESeqDataSetFromMatrix(countData = GSE126428.MPP.FilterCounts,
                                            colData = GSE126428.MPP.coldata,
                                            design= ~ Condition)

GSE126428.MPP.dds$Condition <- relevel(GSE126428.MPP.dds$Condition, ref = "Ctrl")

GSE126428.MPP.dds <- DESeq(GSE126428.MPP.dds)

plotDispEsts(GSE126428.MPP.dds)

GSE126428.MPP.vsd <- vst(GSE126428.MPP.dds, blind=FALSE)

plotPCA(GSE126428.MPP.vsd,intgroup=c("Condition"))

GSE126428.MPP.res <- results(GSE126428.MPP.dds,alpha = 0.05)
GSE126428.MPP.res <- GSE126428.MPP.res[order(GSE126428.MPP.res$pvalue),]

summary(GSE126428.MPP.res)

GSE126428.MPP.res <- as.data.frame(GSE126428.MPP.res) %>% mutate(padj2 = p.adjust(pvalue, method = "fdr"))
GSE126428.MPP.res$external_gene_name <- rownames(GSE126428.MPP.res)

GSE126428.MPP.res <- GSE126428.MPP.res %>% mutate(score=sqrt(log2FoldChange^2 + (-log10(padj2))^2))

counts_norm <- counts(GSE126428.MPP.dds, normalized = TRUE) + 0.5

counts_norm <- counts_norm %>%
  as.data.frame() %>%
  mutate(external_gene_name = rownames(counts_norm)) %>%
  dplyr::select(external_gene_name,everything())

GSE126428.MPP.res %>%
  left_join(counts_norm
            , by = c("external_gene_name")) %>%
  write.xlsx(file = "MPP_Il1KO_DESeq2_New.xlsx", overwrite = T)

assign(objects[3], GSE126428.MPP.res)

#########   CD63 Analysis   #############

objects <- c(objects,"CD63_Il1rnKO")

# Select MPP Experiment Samples
GSE126428.CD63.aggCounts <- GSE126428.aggCounts %>% dplyr::select(GeneName,contains("CD63"))

# Filter genes: at least 1 cpm in at least 3 samples
idx <- which(rowSums(cpm(GSE126428.CD63.aggCounts[,-c(1,2)]) >= 1) >= 3)
GSE126428.CD63.FilterCounts <- GSE126428.CD63.aggCounts[idx,]

GSE126428.CD63.FilterCounts <- as.data.frame(GSE126428.CD63.FilterCounts)
rownames(GSE126428.CD63.FilterCounts) <- GSE126428.CD63.FilterCounts$GeneName

GSE126428.CD63.FilterCounts$GeneName <- NULL

GSE126428.CD63.coldata <- data.frame(Sample=colnames(GSE126428.CD63.FilterCounts),
                                     Condition=factor(c("IL1KO","IL1KO","IL1KO","Ctrl","Ctrl","Ctrl")))

rownames(GSE126428.CD63.coldata) <- GSE126428.CD63.coldata$Sample

GSE126428.CD63.dds <- DESeqDataSetFromMatrix(countData = GSE126428.CD63.FilterCounts,
                                             colData = GSE126428.CD63.coldata,
                                             design= ~ Condition)

GSE126428.CD63.dds$Condition <- relevel(GSE126428.CD63.dds$Condition, ref = "Ctrl")

GSE126428.CD63.dds <- DESeq(GSE126428.CD63.dds)

plotDispEsts(GSE126428.CD63.dds)

GSE126428.CD63.vsd <- vst(GSE126428.CD63.dds, blind=FALSE)

plotPCA(GSE126428.CD63.vsd,intgroup=c("Condition"))


GSE126428.CD63.res <- results(GSE126428.CD63.dds,alpha = 0.05)
GSE126428.CD63.res <- GSE126428.CD63.res[order(GSE126428.CD63.res$pvalue),]

summary(GSE126428.CD63.res)

GSE126428.CD63.res <- as.data.frame(GSE126428.CD63.res) %>% mutate(padj2 = p.adjust(pvalue, method = "fdr"))
GSE126428.CD63.res$external_gene_name <- rownames(GSE126428.CD63.res)

GSE126428.CD63.res <- GSE126428.CD63.res %>% mutate(score=sqrt(log2FoldChange^2 + (-log10(padj2))^2))

counts_norm <- counts(GSE126428.CD63.dds, normalized = TRUE) + 0.5

counts_norm <- counts_norm %>%
  as.data.frame() %>%
  mutate(external_gene_name = rownames(counts_norm)) %>%
  dplyr::select(external_gene_name,everything())

GSE126428.CD63.res %>%
  left_join(counts_norm
            , by = c("external_gene_name")) %>%
  write.xlsx(file = "CD63_Il1KO_DESeq2_New.xlsx", overwrite = T)

assign(objects[4], GSE126428.CD63.res)


#######################################################
##########    Import NfkB Target Genes     ############
#######################################################

all <- rbind(LT_Il1rnKO,ST_Il1rnKO,MPP_Il1rnKO)

nfkbList1 <- read_delim("NfkbTargets1.txt",delim = "\t",col_names = T) %>%
  mutate(Source="Gilmore", MouseGeneName=NA) %>%
  dplyr::select(MouseGeneName,HumanGeneName,Acc,Category,Evidence,Source) %>%
  filter(!is.na(Acc))
nfkbList1 <- nfkbList1[match(unique(sort(nfkbList1$Acc)),nfkbList1$Acc),]

nfkbList2 <- read_delim("NfkbTargets2.txt",delim = "\t",col_names = T) %>%
  mutate(Category = "Gosselin", Source="Gosselin", MouseGeneName=NA) %>%
  rename("HumanGeneName"=Name,"Acc"=RefSeq_Human) %>%
  dplyr::select(MouseGeneName,HumanGeneName,Acc,Category,Evidence,Source) %>%
  filter(!is.na(Acc))
nfkbList2 <- nfkbList2[match(unique(sort(nfkbList2$Acc)),nfkbList2$Acc),]

nfkbList3 <- read_delim("List_NFKB_gene_Arranz.txt",delim = "\t",col_names = T) %>%
  distinct(TF,Target) %>%
  mutate(Evidence="P",Source="Arranz",Acc=NA,HumanGeneName=NA) %>%
  rename("MouseGeneName"=Target,"Category"=TF) %>%
  dplyr::select(MouseGeneName,HumanGeneName,Acc,Category,Evidence,Source)

# Map human to mouse genes
commonName <- "mouse"
wspecies <- "Homo_sapiens"
sp <- "hsapiens"
# spGeneNames <- "mgi_symbol"

genebuildVersion <- "91"
genomeVersion <- "GRCm38"
archive <- "dec2017.archive"

mart <- useMart(host=paste(archive,"ensembl.org",sep = "."), path="/biomart/martservice"
                , biomart = "ENSEMBL_MART_ENSEMBL"
                , port=80,  verbose = T)

mart <- useDataset(paste(sp,"gene","ensembl",sep="_"),mart)

attrList <- listAttributes(mart)

# Source Gilmore
genesMetadataNames <- getBM(attributes = c(
  "external_gene_name"
  ,"ensembl_gene_id"
  ,"refseq_mrna"
)
,uniqueRows=T, mart = mart, filters = "refseq_mrna", values =  nfkbList1$Acc
)
genesMetadata <- getBM(attributes = c("ensembl_gene_id"
                                      ,"mmusculus_homolog_associated_gene_name"
                                      ,"mmusculus_homolog_orthology_type"
                                      ,"mmusculus_homolog_perc_id"
                                      ,"mmusculus_homolog_perc_id_r1"
                                      ,"mmusculus_homolog_orthology_confidence")
                       ,uniqueRows=T, mart = mart, filters = "refseq_mrna", values =  nfkbList1$Acc)
genesMetadata <- genesMetadata %>% filter(mmusculus_homolog_orthology_type %in% c("ortholog_one2one"))
genesMetadata <- genesMetadata %>%
  filter(mmusculus_homolog_orthology_type %in% c("ortholog_one2one")) %>%
  left_join(genesMetadataNames,by = "ensembl_gene_id")

nfkbList1 <- nfkbList1 %>% left_join(genesMetadata,by = c("Acc"="refseq_mrna"))
nfkbList1 <- nfkbList1 %>% mutate(MouseGeneName = nfkbList1$mmusculus_homolog_associated_gene_name)
nfkbList1 <- nfkbList1 %>% filter(!is.na(MouseGeneName))

# Source Gosselin
genesMetadataNames2 <- getBM(attributes = c(
  "external_gene_name"
  ,"ensembl_gene_id"
  ,"refseq_mrna"
)
,uniqueRows=T, mart = mart, filters = "refseq_mrna", values =  nfkbList2$Acc
)
genesMetadata2 <- getBM(attributes = c(
  "ensembl_gene_id"
  ,"mmusculus_homolog_associated_gene_name"
  ,"mmusculus_homolog_orthology_type"
  ,"mmusculus_homolog_perc_id"
  ,"mmusculus_homolog_perc_id_r1"
  ,"mmusculus_homolog_orthology_confidence")
  ,uniqueRows=T, mart = mart, filters = "refseq_mrna", values =  nfkbList2$Acc)
genesMetadata2 <- genesMetadata2 %>%
  filter(mmusculus_homolog_orthology_type %in% c("ortholog_one2one")) %>%
  left_join(genesMetadataNames2,by = "ensembl_gene_id")

nfkbList2 <- nfkbList2 %>% left_join(genesMetadata2,by = c("Acc"="refseq_mrna"))
nfkbList2 <- nfkbList2 %>% mutate(MouseGeneName = nfkbList2$mmusculus_homolog_associated_gene_name)
nfkbList2 <- nfkbList2 %>% filter(!is.na(MouseGeneName))

nfkbGeneNames <- unique(sort(c(nfkbList1$MouseGeneName,nfkbList2$MouseGeneName,nfkbList3$MouseGeneName)))

allTargets <- data.frame(MouseGeneName=nfkbGeneNames)

allTargets <- allTargets %>% left_join(nfkbList1 %>% filter(!is.na(MouseGeneName)) %>% dplyr::select(MouseGeneName,Category), by = c("MouseGeneName"="MouseGeneName")) %>%
  rename("Gilmore"=Category)

allTargets <- allTargets %>% left_join(nfkbList2 %>% filter(!is.na(MouseGeneName)) %>% dplyr::select(MouseGeneName,Category), by = c("MouseGeneName"="MouseGeneName")) %>%
  rename("Gosselin"=Category)

for (arcat in c("NFKB","NFKB1","NFKB2","REL","RELA")) {
  allTargets <- allTargets %>% left_join(nfkbList3 %>% dplyr::select(MouseGeneName,Category) %>% distinct(MouseGeneName,Category) %>% filter(Category == arcat), by = c("MouseGeneName"="MouseGeneName")) %>%
    rename(!!paste0("Arranz_",arcat):=Category)
}

allTargets <- allTargets %>%
  mutate(LT_Il1rnKO=ifelse(!MouseGeneName %in% LT_Il1rnKO$external_gene_name,NA,ifelse(MouseGeneName %in% LT_Il1rnKO[LT_Il1rnKO$padj2 < 0.05,"external_gene_name"],T,F))) %>%
  mutate(ST_Il1rnKO=ifelse(!MouseGeneName %in% ST_Il1rnKO$external_gene_name,NA,ifelse(MouseGeneName %in% ST_Il1rnKO[ST_Il1rnKO$padj2 < 0.05,"external_gene_name"],T,F))) %>%
  mutate(MPP_Il1rnKO=ifelse(!MouseGeneName %in% MPP_Il1rnKO$external_gene_name,NA,ifelse(MouseGeneName %in% MPP_Il1rnKO[MPP_Il1rnKO$padj2 < 0.05,"external_gene_name"],T,F)))

allTargets %>% write.xlsx(file = "NFKB.TargetsLists.xlsx",overwrite = T)

###########################################################################
#######       Annotate DE Tables with Nfkb targets         ################
###########################################################################

#LT
LT_Il1rnKO$isTarget <- F

LT_Il1rnKO <- LT_Il1rnKO %>% left_join(nfkbList1 %>% dplyr::select(MouseGeneName,Category), by = c("external_gene_name"="MouseGeneName")) %>%
  mutate(isTarget=ifelse(!is.na(Category),T,isTarget)) %>%
  rename("Gilmore"=Category)

LT_Il1rnKO <- LT_Il1rnKO %>% left_join(nfkbList2 %>% dplyr::select(MouseGeneName,Category), by = c("external_gene_name"="MouseGeneName")) %>%
  mutate(isTarget=ifelse(!is.na(Category),T,isTarget)) %>%
  rename("Gosselin"=Category)

for (arcat in c("NFKB","NFKB1","NFKB2","REL","RELA")) {
  LT_Il1rnKO <- LT_Il1rnKO %>% left_join(nfkbList3 %>%
                                           dplyr::select(MouseGeneName,Category) %>%
                                           distinct(MouseGeneName,Category) %>%
                                           filter(Category == arcat), by = c("external_gene_name"="MouseGeneName")) %>%
    mutate(isTarget=ifelse(!is.na(Category),T,isTarget)) %>%
    rename(!!paste0("Arranz_",arcat):=Category)
}

LT_Il1rnKO <- LT_Il1rnKO %>% mutate(score=sqrt(log2FoldChange^2 + (-log10(padj2))^2)) %>% arrange(desc(score))
LT_Il1rnKO %>% write.xlsx(file = "LT_Il1rnKO.nfkbAnnot.xlsx",overwrite = T)

#ST
ST_Il1rnKO$isTarget <- F
ST_Il1rnKO <- ST_Il1rnKO %>% left_join(nfkbList1 %>% dplyr::select(MouseGeneName,Category), by = c("external_gene_name"="MouseGeneName")) %>%
  mutate(isTarget=ifelse(!is.na(Category),T,isTarget)) %>%
  rename("Gilmore"=Category)

ST_Il1rnKO <- ST_Il1rnKO %>% left_join(nfkbList2 %>% dplyr::select(MouseGeneName,Category), by = c("external_gene_name"="MouseGeneName")) %>%
  mutate(isTarget=ifelse(!is.na(Category),T,isTarget)) %>%
  rename("Gosselin"=Category)

for (arcat in c("NFKB","NFKB1","NFKB2","REL","RELA")) {
  ST_Il1rnKO <- ST_Il1rnKO %>% left_join(nfkbList3 %>%
                                           dplyr::select(MouseGeneName,Category) %>%
                                           distinct(MouseGeneName,Category) %>%
                                           filter(Category == arcat), by = c("external_gene_name"="MouseGeneName")) %>%
    mutate(isTarget=ifelse(!is.na(Category),T,isTarget)) %>%
    rename(!!paste0("Arranz_",arcat):=Category)
}
ST_Il1rnKO <- ST_Il1rnKO %>% mutate(score=sqrt(log2FoldChange^2 + (-log10(padj2))^2)) %>% arrange(desc(score))
ST_Il1rnKO %>% write.xlsx(file = "ST_Il1rnKO.nfkbAnnot.xlsx",overwrite = T)

#MPP
MPP_Il1rnKO$isTarget <- F
MPP_Il1rnKO <- MPP_Il1rnKO %>% left_join(nfkbList1 %>% dplyr::select(MouseGeneName,Category), by = c("external_gene_name"="MouseGeneName")) %>%
  mutate(isTarget=ifelse(!is.na(Category),T,isTarget)) %>%
  rename("Gilmore"=Category)

MPP_Il1rnKO <- MPP_Il1rnKO %>% left_join(nfkbList2 %>% dplyr::select(MouseGeneName,Category), by = c("external_gene_name"="MouseGeneName")) %>%
  mutate(isTarget=ifelse(!is.na(Category),T,isTarget)) %>%
  rename("Gosselin"=Category)

for (arcat in c("NFKB","NFKB1","NFKB2","REL","RELA")) {
  MPP_Il1rnKO <- MPP_Il1rnKO %>% left_join(nfkbList3 %>%
                                             dplyr::select(MouseGeneName,Category) %>%
                                             distinct(MouseGeneName,Category) %>%
                                             filter(Category == arcat), by = c("external_gene_name"="MouseGeneName")) %>%
    mutate(isTarget=ifelse(!is.na(Category),T,isTarget)) %>%
    rename(!!paste0("Arranz_",arcat):=Category)
}
MPP_Il1rnKO <- MPP_Il1rnKO %>% mutate(score=sqrt(log2FoldChange^2 + (-log10(padj2))^2)) %>% arrange(desc(score))
MPP_Il1rnKO %>% write.xlsx(file = "MPP_Il1rnKO.nfkbAnnot.xlsx",overwrite = T)

#######################################################################################
#######################       Volcanos in HSPCs (Fig. 2N)       #######################
#######################################################################################

# Volcano Plots with Nfkb targets highlighted:
# manually curate genes to label and create a column "toLabel" in .nfkbAnnot.xlsx files created with value =1 for genes to label

LT_Il1rnKO <- as.data.frame(read.xlsx(file.path("LT_Il1rnKO.nfkbAnnot_genes.to.label.xlsx")))
ST_Il1rnKO <- as.data.frame(read.xlsx(file.path("ST_Il1rnKO.nfkbAnnot_genes.to.label.xlsx")))
MPP_Il1rnKO <- as.data.frame(read.xlsx(file.path("MPP_Il1rnKO.nfkbAnnot_genes.to.label.xlsx")))


objects <- c("LT_Il1rnKO","ST_Il1rnKO","MPP_Il1rnKO","CD63_Il1rnKO")


toLabel = list()
for (i in 1:3) {
  toLabel[[objects[i]]] <- read.xlsx(paste0(objects[i],".nfkbAnnot_genes.to.label.xlsx"))
  toLabel[[objects[i]]] <- toLabel[[objects[i]]] %>% filter(toLabel == 1) %>% dplyr::select(external_gene_name)
}
lapply(toLabel, dim)

pCutoff=0.05
FCcutoff=0.5
for(i in 1:3) {
  temp <- get(objects[i])
  # temp <- temp %>% mutate(type = objects[i])
  temp <- temp %>% mutate(toLabel = ifelse(is.na(toLabel),F,T))
  # Get list of Ribosomic related genes
  riboMtOthers = temp[grep("^Rpl|^Rps|^mt-|^Ig[a-z]v", temp$external_gene_name),"external_gene_name"]
  
  # Remove genes with low expression & pvalue over 0.05 | p_adjusted_value over 1
  temp <- temp %>% filter(!is.na(padj)) %>%
    filter(!(baseMean < 500 & pvalue > 0.05)) %>%
    filter(padj < 1)
  
  # Order by distance to 0 in the vulcano plot
  temp <- temp %>% arrange(-score)
  
  # Remove predicted genes (Gm)
  temp <- temp %>% filter(!grepl("Gm[0-9]+", external_gene_name))
  
  # Label selected Nfkb genes except ribo genes
  temp <- temp %>% mutate(anno = ifelse(toLabel == T,external_gene_name,"")) %>%
    mutate(anno = ifelse(external_gene_name %in% riboMtOthers,"",anno))
  
  # Mark as red all genes below 0.05 and with LFC over 0.5
  temp <- temp %>% mutate(fillme = "black", size = 0.1) %>%
    mutate(fillme = ifelse(padj < pCutoff & abs(log2FoldChange) > FCcutoff,"red",fillme))
  
  # Mark as blue all significant Nkfb targets and increase size and border color of labeled genes
  temp <- temp %>% mutate(fillme = ifelse(isTarget & fillme == "red","blue",fillme)) %>%
    # mutate(fillme = ifelse(anno != "" & padj < pCutoff,"blue",fillme)) %>%
    mutate(colorme = ifelse(anno != "","black",fillme)) %>%
    mutate(sizeme = ifelse(anno != "",1,size))
  
  temp <- temp %>% mutate(log2FoldChange = round(log2FoldChange, 2)) %>%
    mutate(padj = signif(temp$padj, 3))
  
  temp <- temp %>% arrange(anno,isTarget,desc(score))
  
  p <- as.data.frame(temp) %>% ggplot(aes(x = log2FoldChange, y = -log10(padj))
                       # ,fill=factor(fillme),color=factor(colorme),size=sizeme
                       ) + 
    # geom_point() +
    geom_point(shape = 21, fill = temp$fillme, color = temp$colorme, size = temp$sizeme) +
    geom_text_repel(size = 2, label = temp$anno,
                    point.padding = unit(1e-10, "lines"),
                    #segment.color = 'transparent',
                    box.padding = 0.1,
                    max.overlaps = 500,
                    seed = 42
    ) + 
    ggtitle(objects[i]) +
    theme(panel.background = element_blank(),
          axis.line = element_line(),
          panel.border = element_rect(fill = "transparent")
    ) +
    ylab("-log10 (Adjusted p value)") + xlab("log2 FC")
  
  assign(paste("p", objects[i], sep = "."), p)
  
}

pdf("Il1rnKO_LT_Vulcanos.pdf")
p.LT_Il1rnKO
dev.off()

pdf("Il1rnKO_ST_Vulcanos.pdf")
p.ST_Il1rnKO
dev.off()

pdf("Il1rnKO_MPP_Vulcanos.pdf")
p.MPP_Il1rnKO
dev.off()

#######################################################################################
########################            IL-1rn-KO                 #########################
########################       Venn in HSPCs (Fig. S4A)       #########################
#######################################################################################

# Overlap between contrasts of significantly upregulated genes
deList <- list(
  LT_Il1rnKO = LT_Il1rnKO %>% filter(padj < 0.05 & log2FoldChange > 0) %>% .$external_gene_name,
  ST_Il1rnKO = ST_Il1rnKO %>% filter(padj < 0.05 & log2FoldChange > 0) %>% .$external_gene_name,
  MPP_Il1rnKO = MPP_Il1rnKO %>% filter(padj < 0.05 & log2FoldChange > 0) %>% .$external_gene_name
)

pdf("Il1rnKO_UpregulatedOverlaps.pdf")
ggvenn(deList, names(deList),
       show_percentage = F,
       fill_color = color.list,
       set_name_size = 4,
       text_size = 5) +
  ggtitle("Upregulated genes")
dev.off()

# Overlap between contrasts of significantly downregulated genes
deList <- list(
  LT_Il1rnKO = LT_Il1rnKO %>% filter(padj < 0.05 & log2FoldChange < 0) %>% .$external_gene_name,
  ST_Il1rnKO = ST_Il1rnKO %>% filter(padj < 0.05 & log2FoldChange < 0) %>% .$external_gene_name,
  MPP_Il1rnKO = MPP_Il1rnKO %>% filter(padj < 0.05 & log2FoldChange < 0) %>% .$external_gene_name
)

pdf("Il1rnKO_DownregulatedOverlaps.pdf")
ggvenn(deList, names(deList),
       show_percentage = F,
       fill_color = color.list,
       set_name_size = 4,
       text_size = 5) +
  ggtitle("Downregulated genes")
dev.off()

#######################################################################################
########################            IL-1rn-KO                 #########################
########################      Volcano in CD63+ (Fig. 4F)      #########################
#######################################################################################

# CD63 Experiment Volcano Plot
CD63_Il1rnKO <- as.data.frame(read.xlsx(file.path("CD63_Il1KO_DESeq2_New.xlsx")))
toLabel <- c("Wif1","Id3","Dab2","Zfp9","Ctgf","Angpt1","Fbln5","Cdh2","Lepr","Adipoq","Thsd4","Cxcl12","Vcam1","Cd68","Ctsg","Mpo","Podxl","Csf2rb","Unc13d","Lrg1","Wfdc17","Mfsd12","Prtn3")

pCutoff=0.05
FCcutoff=0.1

temp <- get(objects[4])

# Get list of Ribosomic related genes
riboMtOthers = temp[grep("^Rpl|^Rps|^mt-|^Ig[a-z]v", temp$external_gene_name),"external_gene_name"]

# Remove genes with low expression & pvalue over 0.05 | p_adjusted_value over 1
temp <- temp %>% filter(!is.na(padj)) %>%
  filter(padj < 1)

# Order by distance to 0 in the vulcano plot
temp <- temp %>% mutate(score=sqrt(log2FoldChange^2 + (-log10(padj2))^2)) %>% arrange(desc(score)) %>% arrange(padj)

# Remove predicted genes (Gm)
temp <- temp %>% filter(!grepl("Gm[0-9]+", external_gene_name))

# Label top 10 up and top 10 down genes plus selected ones
dummy <- temp %>% mutate(up = log2FoldChange > 0) %>% group_split(up)
toLabel <- unique(sort(c(toLabel,unlist(lapply(dummy,function (x) {return(head(x %>% filter(padj < 0.05) %>% .$external_gene_name,10))})))))

# Label selected genes except ribo genes
temp <- temp %>% mutate(anno = ifelse(external_gene_name %in% toLabel,external_gene_name,"")) %>%
  mutate(anno = ifelse(external_gene_name %in% riboMtOthers,"",anno))

temp <- temp %>% mutate(fill = "black", size = 0.1) %>%
  mutate(fill = ifelse(padj < pCutoff & abs(log2FoldChange) > FCcutoff,"red",fill)) %>%
  mutate(color = ifelse(anno != "","black",fill)) %>%
  mutate(size = ifelse(anno != "",1,size)) %>%
  mutate(log2FoldChange = round(log2FoldChange, 2)) %>%
  mutate(padj = signif(temp$padj, 3)) %>%
  arrange(desc(anno),desc(score))

p <- temp %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj))) + 
  geom_point(shape = 21, fill = temp$fill, color = temp$color, size = temp$size) +
  geom_text_repel(size = 2, label = temp$anno,
                  point.padding = unit(1e-10, "lines"),
                  segment.color = 'black',
                  box.padding = 0.1,
                  max.overlaps = 1000,
                  seed = 42
  ) +
  ggtitle(objects[4]) +
  theme(panel.background = element_blank(),
        axis.line = element_line(),
        panel.border = element_rect(fill = "transparent")
  ) +
  ylab("-log10 (Adjusted p value)") + xlab("log2 FC")

assign(paste("p", objects[4], sep = "."), p)

pdf("Il1rnKO_CD63_Vulcano.pdf")
print(p.CD63_Il1rnKO)
dev.off()

#######################################################################################
#########                           Contrasts                 #########################
#########      GSE165810 (Fig. S3I) and GSE166629 (Fig. S3J)  #########################
#######################################################################################

# Overlaps with dataset GSE165810
GSE165810.CountsFile <- "GSE165810.RawCounts.txt"

GSE165810.RawCounts <- read.delim(GSE165810.CountsFile,sep = "\t",header = T)

dim(GSE165810.RawCounts)
head(GSE165810.RawCounts)

# Filter low expressed genes: At least 1 cpm in at least 4 samples
idx <- which(rowSums(cpm(GSE165810.RawCounts[,-1]) >= 1) >= 4)

GSE165810.FilterCounts <- GSE165810.RawCounts[idx,]
rownames(GSE165810.FilterCounts) <- GSE165810.FilterCounts$Geneid
GSE165810.FilterCounts$Geneid <- NULL
GSE165810.coldata <- data.frame(Sample=colnames(GSE165810.FilterCounts),
                                Condition=factor(c("IL1KO","IL1KO","Ctrl","Ctrl","Ctrl","Ctrl","IL1KO","IL1KO")),
                                Batch=factor(rep(c(1,2),each=4))
)
rownames(GSE165810.coldata) <- GSE165810.coldata$Sample

GSE165810.dds <- DESeqDataSetFromMatrix(countData = GSE165810.FilterCounts,
                                        colData = GSE165810.coldata,
                                        design= ~ Batch + Condition)

GSE165810.dds$Condition <- relevel(GSE165810.dds$Condition, ref = "Ctrl")

GSE165810.dds <- DESeq(GSE165810.dds)

GSE165810.vsd <- vst(GSE165810.dds, blind=FALSE)

plotPCA(GSE165810.vsd,intgroup=c("Condition", "Batch"))

GSE165810.res <- results(GSE165810.dds,alpha = 0.05)
GSE165810.res <- GSE165810.res[order(GSE165810.res$pvalue),]


summary(GSE165810.res)

GSE165810.res <- as.data.frame(GSE165810.res) %>% mutate(padj2 = p.adjust(pvalue, method = "fdr"))

GSE165810.res %>% mutate(external_gene_name = rownames(GSE165810.res)) %>% as.data.frame() %>% write.xlsx(file = "GSE165810.res.xlsx",overwrite = T)

objects <- c(objects,"GSE165810")


GSE165810 <- as.data.frame(GSE165810.res)
GSE165810$external_gene_name <- rownames(GSE165810)
GSE165810$padj2 <- GSE165810$padj

allSig <- data.frame(external_gene_name=unique(sort(c(
  GSE165810 %>% filter(padj < 0.05) %>% .$external_gene_name,
  LT_Il1rnKO %>% filter(padj < 0.05) %>% .$external_gene_name,
  ST_Il1rnKO %>% filter(padj < 0.05) %>% .$external_gene_name,
  MPP_Il1rnKO %>% filter(padj < 0.05) %>% .$external_gene_name))))

for (i in 1:length(objects)) {
  allSig <- allSig %>% left_join(get(objects[i]) %>% filter(padj < 0.05) %>% dplyr::select(external_gene_name,log2FoldChange),by = c("external_gene_name")) %>%
    rename(!!objects[i]:=log2FoldChange)
}

allSig$isTarget <- F
allSig[allSig$external_gene_name %in% nfkbGeneNames,"isTarget"] <- T
allSig %>% write.xlsx(file = "GSE165810.AllSig.xlsx",overwrite=T)


#Perform overlap
GSE165810 <- read.xlsx("GSE165810.res.xlsx")
GSE165810 <- as.data.frame(GSE165810)
#GSE165810$external_gene_name <- rownames(GSE165810)
GSE165810$padj2 <- GSE165810$padj

deList <- list(
  LT_Il1rnKO = LT_Il1rnKO %>% filter(padj < 0.05 & log2FoldChange > 0) %>% .$external_gene_name,
  ST_Il1rnKO = ST_Il1rnKO %>% filter(padj < 0.05 & log2FoldChange > 0) %>% .$external_gene_name,
  MPP_Il1rnKO = MPP_Il1rnKO %>% filter(padj < 0.05 & log2FoldChange > 0) %>% .$external_gene_name,
  GSE165810 = GSE165810 %>% filter(padj < 0.05 & log2FoldChange > 0) %>% .$external_gene_name
)

p <- ggvenn(deList, names(deList),
             show_percentage = F,
             fill_color = color.list,
             set_name_size = 4,
             text_size = 5) +
  ggtitle("Upregulated genes")

pdf("Il1rnKO_LTSTMPP_GSE16518_UpregulatedOverlaps.pdf")
print(p)
dev.off()

deList <- list(
  LT_Il1rnKO = LT_Il1rnKO %>% filter(padj < 0.05 & log2FoldChange < 0) %>% .$external_gene_name,
  ST_Il1rnKO = ST_Il1rnKO %>% filter(padj < 0.05 & log2FoldChange < 0) %>% .$external_gene_name,
  MPP_Il1rnKO = MPP_Il1rnKO %>% filter(padj < 0.05 & log2FoldChange < 0) %>% .$external_gene_name,
  GSE165810 = GSE165810 %>% filter(padj < 0.05 & log2FoldChange < 0) %>% .$external_gene_name
)

p <- ggvenn(deList, names(deList),
             show_percentage = F,
             fill_color = color.list,
             set_name_size = 4,
             text_size = 5) +
  ggtitle("Downregulated genes")

pdf("Il1rnKO_LTSTMPP_GSE16518_DownregulatedOverlaps.pdf")
print(p)
dev.off()

# Overlaps with dataset GSE166629

GSE166629.RawCounts <- read.delim(file = "GSE166629_Compiled_counts.txt",header = T,sep = "\t")
dim(GSE166629.RawCounts)

head(GSE166629.RawCounts)[,1:5]

colnames(GSE166629.RawCounts)

rownames(GSE166629.RawCounts) <- GSE166629.RawCounts$Gene
GSE166629.RawCounts$Gene <- NULL
GSE166629.RawCounts <- as.matrix(GSE166629.RawCounts)

sampleMetadataGSE166629 <- read.delim("GSE166629_SampleMetadata.txt",sep = "\t")
sampleMetadataGSE166629 <- as.data.frame(t(sampleMetadataGSE166629))
colnames(sampleMetadataGSE166629) <- sampleMetadataGSE166629[1,]
sampleMetadataGSE166629 <- sampleMetadataGSE166629[-1,]
rownames(sampleMetadataGSE166629) <- sub("\\.","",rownames(sampleMetadataGSE166629))
sampleMetadataGSE166629$SampleName <- rownames(sampleMetadataGSE166629)
rownames(sampleMetadataGSE166629) <- sampleMetadataGSE166629$Sample_description
GSE166629.coldata <- sampleMetadataGSE166629 %>% rename("Condition"=treatment,"Group"=yfp)
rownames(GSE166629.coldata) <- GSE166629.coldata$SampleName

GSE166629.coldata <- GSE166629.coldata %>% filter(mouse_genotype == "WT")

# Extract the values from the "Sample_description" column of GSE166629.coldata and store them in a variable, e.g. "matching_values"
matching_values <- GSE166629.coldata$Sample_description
# Extract the column names of GSE166629.RawCounts
old_colnames <- colnames(GSE166629.RawCounts)
# Use match() function to match the matching_values with old_colnames
match_indices <- match(matching_values, old_colnames)
# Remove non-matching values
match_indices<-match_indices[match_indices>0]
# Subset the GSE166629.RawCounts using the match_indices
GSE166629.RawCounts <- GSE166629.RawCounts[, match_indices]
#Extract the row names of GSE166629.coldata using the match_indices and store them in a variable, e.g. "new_colnames"
new_colnames <- rownames(GSE166629.coldata)[match_indices]
#Use the colnames() function to change the column names of GSE166629.RawCounts to "new_colnames".
colnames(GSE166629.RawCounts) <- new_colnames


## YFP_POS

YFP_POS_Group.coldata <- GSE166629.coldata %>% filter(Group == "YFP_POS" & mouse_genotype == "WT")

# Extract the row names of YFP_POS_Group.coldata
matching_values_YFP_POS <- rownames(YFP_POS_Group.coldata)
# Extract the column names of GSE166629.RawCounts
old_colnames_YFP_POS <- colnames(GSE166629.RawCounts)
# Use match() function to match the matching_values with old_colnames
match_indices_YFP_POS <- match(matching_values_YFP_POS, old_colnames_YFP_POS)
# Remove non-matching values
match_indices_YFP_POS<-match_indices_YFP_POS[match_indices_YFP_POS>0]
# Subset the GSE166629.RawCounts using the match_indices
YFP_POS_Group.RawCounts <- GSE166629.RawCounts[, match_indices_YFP_POS]


#YFP_POS_Group.RawCounts <- GSE166629.RawCounts[,rownames(YFP_POS_Group.coldata)]


# Filter low expressed genes: At least 1 cpm in at least 3 samples
idx <- which(rowSums(cpm(YFP_POS_Group.RawCounts) >= 1) >= 3)

YFP_POS_Group.FilterCounts <- YFP_POS_Group.RawCounts[idx,]

YFP_POS_Group.dds <- DESeqDataSetFromMatrix(countData = YFP_POS_Group.FilterCounts,
                                            colData = YFP_POS_Group.coldata,
                                            design= ~ Condition)

YFP_POS_Group.dds$Condition <- relevel(YFP_POS_Group.dds$Condition, ref = "PBS")

YFP_POS_Group.dds <- DESeq(YFP_POS_Group.dds)

YFP_POS_Group.vsd <- vst(YFP_POS_Group.dds, blind=FALSE)

plotPCA(YFP_POS_Group.vsd,intgroup=c("Condition"))

YFP_POS_Group.res <- results(YFP_POS_Group.dds,alpha = 0.05)
YFP_POS_Group.res <- YFP_POS_Group.res[order(YFP_POS_Group.res$pvalue),]

summary(YFP_POS_Group.res)

YFP_POS_Group.res <- as.data.frame(YFP_POS_Group.res) %>% mutate(padj2 = p.adjust(pvalue, method = "fdr"))

YFP_POS_Group.res %>% mutate(external_gene_name = rownames(YFP_POS_Group.res)) %>%  as.data.frame() %>% write.xlsx(file = "YFP_POS_Group.res.xlsx",overwrite = T)

objects <- c("LT_Il1rnKO","ST_Il1rnKO","MPP_Il1rnKO","YFP_POS_Group")

YFP_POS_Group <- read.xlsx("YFP_POS_Group.res.xlsx")
YFP_POS_Group <- as.data.frame(YFP_POS_Group)
# YFP_POS_Group$external_gene_name <- rownames(YFP_POS_Group)
YFP_POS_Group$padj2 <- YFP_POS_Group$padj

allSig <- data.frame(external_gene_name=unique(sort(c(
  YFP_POS_Group %>% filter(padj < 0.05) %>% .$external_gene_name,
  LT_Il1rnKO %>% filter(padj < 0.05) %>% .$external_gene_name,
  ST_Il1rnKO %>% filter(padj < 0.05) %>% .$external_gene_name,
  MPP_Il1rnKO %>% filter(padj < 0.05) %>% .$external_gene_name))))

for (i in 1:length(objects)) {
  allSig <- allSig %>% left_join(get(objects[i]) %>% filter(padj < 0.05) %>% dplyr::select(external_gene_name,log2FoldChange),by = c("external_gene_name")) %>%
    rename(!!objects[i]:=log2FoldChange)
}

allSig$isTarget <- F
allSig[allSig$external_gene_name %in% nfkbGeneNames,"isTarget"] <- T
allSig %>% write.xlsx(file = "YFP_POS_Group.AllSig.xlsx",overwrite=T)

### Upregulated
deList <- list(
  LT_Il1rnKO = LT_Il1rnKO %>% filter(padj < 0.05 & log2FoldChange > 0) %>% .$external_gene_name,
  ST_Il1rnKO = ST_Il1rnKO %>% filter(padj < 0.05 & log2FoldChange > 0) %>% .$external_gene_name,
  MPP_Il1rnKO = MPP_Il1rnKO %>% filter(padj < 0.05 & log2FoldChange > 0) %>% .$external_gene_name,
  YFP_POS_Group = YFP_POS_Group %>% filter(padj < 0.05 & log2FoldChange > 0) %>% .$external_gene_name
)

p <- ggvenn(deList, names(deList),
             show_percentage = F,
             fill_color = color.list,
             set_name_size = 4,
             text_size = 5) +
  ggtitle("Upregulated genes")

pdf("Il1rnKO_LTSTMPP_GSE166629YPos_UpregulatedOverlaps.pdf")
print(p)
dev.off()

deList <- list(
  LT_Il1rnKO = LT_Il1rnKO %>% filter(padj < 0.05 & log2FoldChange < 0) %>% .$external_gene_name,
  ST_Il1rnKO = ST_Il1rnKO %>% filter(padj < 0.05 & log2FoldChange < 0) %>% .$external_gene_name,
  MPP_Il1rnKO = MPP_Il1rnKO %>% filter(padj < 0.05 & log2FoldChange < 0) %>% .$external_gene_name,
  YFP_POS_Group = YFP_POS_Group %>% filter(padj < 0.05 & log2FoldChange < 0) %>% .$external_gene_name
)

### Downregulated
p <- ggvenn(deList, names(deList),
             show_percentage = F,
             fill_color = color.list,
             set_name_size = 4,
             text_size = 5) +
  ggtitle("Downregulated genes")

pdf("Il1rnKO_LTSTMPP_GSE166629YPos_DownregulatedOverlaps.pdf")
print(p)
dev.off()

## YFP_NEG

YFP_NEG_Group.coldata <- GSE166629.coldata %>% filter(Group == "YFP_NEG" & mouse_genotype == "WT")

# Extract the row names of YFP_POS_Group.coldata
matching_values_YFP_NEG <- rownames(YFP_NEG_Group.coldata)
# Extract the column names of GSE166629.RawCounts
old_colnames_YFP_NEG <- colnames(GSE166629.RawCounts)
# Use match() function to match the matching_values with old_colnames
match_indices_YFP_NEG <- match(matching_values_YFP_NEG, old_colnames_YFP_NEG)
# Remove non-matching values
match_indices_YFP_NEG<-match_indices_YFP_NEG[match_indices_YFP_NEG>0]
# Subset the GSE166629.RawCounts using the match_indices
YFP_NEG_Group.RawCounts <- GSE166629.RawCounts[, match_indices_YFP_NEG]

# Filter low expressed genes: At least 1 cpm in at least 3 samples
idx <- which(rowSums(cpm(YFP_NEG_Group.RawCounts) >= 1) >= 3)

YFP_NEG_Group.FilterCounts <- YFP_NEG_Group.RawCounts[idx,]

YFP_NEG_Group.dds <- DESeqDataSetFromMatrix(countData = YFP_NEG_Group.FilterCounts,
                                            colData = YFP_NEG_Group.coldata,
                                            design= ~ Condition)

YFP_NEG_Group.dds$Condition <- relevel(YFP_NEG_Group.dds$Condition, ref = "PBS")

YFP_NEG_Group.dds <- DESeq(YFP_NEG_Group.dds)

YFP_NEG_Group.vsd <- vst(YFP_NEG_Group.dds, blind=FALSE)

plotPCA(YFP_NEG_Group.vsd,intgroup=c("Condition"))

YFP_NEG_Group.res <- results(YFP_NEG_Group.dds,alpha = 0.05)
YFP_NEG_Group.res <- YFP_NEG_Group.res[order(YFP_NEG_Group.res$pvalue),]


summary(YFP_NEG_Group.res)

YFP_NEG_Group.res <- as.data.frame(YFP_NEG_Group.res) %>% mutate(padj2 = p.adjust(pvalue, method = "fdr"))

YFP_NEG_Group.res %>% mutate(external_gene_name = rownames(YFP_NEG_Group.res)) %>% as.data.frame() %>% write.xlsx(file = "YFP_NEG_Group.res.xlsx",overwrite = T)

objects <- c("LT_Il1rnKO","ST_Il1rnKO","MPP_Il1rnKO","YFP_NEG_Group")

YFP_NEG_Group <- read.xlsx("YFP_NEG_Group.res.xlsx")
YFP_NEG_Group <- as.data.frame(YFP_NEG_Group.res)
YFP_NEG_Group$external_gene_name <- rownames(YFP_NEG_Group)
YFP_NEG_Group$padj2 <- YFP_NEG_Group$padj


allSig <- data.frame(external_gene_name=unique(sort(c(
  YFP_NEG_Group %>% filter(padj < 0.05) %>% .$external_gene_name,
  LT_Il1rnKO %>% filter(padj < 0.05) %>% .$external_gene_name,
  ST_Il1rnKO %>% filter(padj < 0.05) %>% .$external_gene_name,
  MPP_Il1rnKO %>% filter(padj < 0.05) %>% .$external_gene_name))))

for (i in 1:length(objects)) {
  allSig <- allSig %>% left_join(get(objects[i]) %>% filter(padj < 0.05) %>% dplyr::select(external_gene_name,log2FoldChange),by = c("external_gene_name")) %>%
    rename(!!objects[i]:=log2FoldChange)
}

allSig$isTarget <- F
allSig[allSig$external_gene_name %in% nfkbGeneNames,"isTarget"] <- T
allSig %>% write.xlsx(file = "YFP_NEG_Group.AllSig.xlsx",overwrite=T)


### Upregulated
deList <- list(
  LT_Il1rnKO = LT_Il1rnKO %>% filter(padj < 0.05 & log2FoldChange > 0) %>% .$external_gene_name,
  ST_Il1rnKO = ST_Il1rnKO %>% filter(padj < 0.05 & log2FoldChange > 0) %>% .$external_gene_name,
  MPP_Il1rnKO = MPP_Il1rnKO %>% filter(padj < 0.05 & log2FoldChange > 0) %>% .$external_gene_name,
  YFP_NEG_Group = YFP_NEG_Group %>% filter(padj < 0.05 & log2FoldChange > 0) %>% .$external_gene_name
)

p <- ggvenn(deList, names(deList),
             show_percentage = F,
             fill_color = color.list,
             set_name_size = 4,
             text_size = 5) +
  ggtitle("Upregulated genes")

pdf("Il1rnKO_LTSTMPP_GSE166629YNeg_UpregulatedOverlaps.pdf")
print(p)
dev.off()

### Downregulated
deList <- list(
  LT_Il1rnKO = LT_Il1rnKO %>% filter(padj < 0.05 & log2FoldChange < 0) %>% .$external_gene_name,
  ST_Il1rnKO = ST_Il1rnKO %>% filter(padj < 0.05 & log2FoldChange < 0) %>% .$external_gene_name,
  MPP_Il1rnKO = MPP_Il1rnKO %>% filter(padj < 0.05 & log2FoldChange < 0) %>% .$external_gene_name,
  YFP_NEG_Group = YFP_NEG_Group %>% filter(padj < 0.05 & log2FoldChange < 0) %>% .$external_gene_name
)

p <- ggvenn(deList, names(deList),
             show_percentage = F,
             fill_color = color.list,
             set_name_size = 4,
             text_size = 5) +
  ggtitle("Downregulated genes")

pdf("Il1rnKO_LTSTMPP_GSE166629YNeg_DownregulatedOverlaps.pdf")
print(p)
dev.off()




#######################################################################################
#######################                   NRAS                  #######################
#######################          Venn HSPCs (Fig. S13I)          #######################
#######################           PCA HSPCs (Fig. S13J)          #######################
#######################################################################################

#read table 
m2   <- read.table("raw_mat_Nras.csv", row.names=1, header=T, sep="\t")
temp <- m2[,grep("_rawCount", colnames(m2))]
keep <- which(rowMeans(temp) > 1) 
raw  <- temp[keep,]
keep <- which(rowMeans(log2(raw+1))>3)
raw  <- raw[keep,]

#annotation 
m2    <- m2[rownames(raw),]
index <- setNames(m2$Description, rownames(m2))
gid   <- setNames(gsub("GeneID:", "", str_extract(m2$Description, "GeneID:[0-9]+")), rownames(m2))
annot <- cbind(as.character(gid[rownames(m2)]), as.character(index[rownames(m2)]))
rownames(annot) <- rownames(m2) 
xx    <- bitr(as.character(gid), "ENTREZID", "SYMBOL",org.Mm.eg.db)
xx    <- setNames(xx[,2], xx[,1])
annot <- data.frame(cbind(annot, as.character(xx[annot[,1]])))
colnames(annot) <- c("ENTREZID", "DESCRIPTION", "SYMBOL")
annot <- annot[, c(1,3,2)]
annot$SYMBOL      <- as.character(annot$SYMBOL)
annot$DESCRIPTION <- as.character(annot$DESCRIPTION)
annot$ENTREZID    <- as.character(annot$ENTREZID)

pheno <- read.table("pheno_Nras.csv", header=T, sep="\t")
colnames(raw) <- gsub("_rawCount", "", colnames(raw))
rownames(pheno) <- gsub("_rawCount", "", pheno$ID)
raw            <- raw[, rownames(pheno)]

#PCA 
temp <- raw[,-grep("CD63", colnames(raw))]
dds  <- DESeqDataSetFromMatrix(countData = temp, colData=pheno[colnames(temp),], design = ~ Group+Type)
dds  <- DESeq(dds)
cnts <- rlog(counts(dds))
dd <- prcomp(t(cnts))
par(xpd=T)
plot(dd$x, pch=20, col=as.character(colData(dds)$Color), cex=1)
text(dd$x, rownames(dd$x), cex=0.7, pos=2)

getDiffExp <- function(raw, pheno, LT, annot, pval=0.05) 
	{
	dds          <- DESeqDataSetFromMatrix(countData = raw[, LT], colData=pheno[LT,], design = ~ Type)
	dds          <- DESeq(dds)
	perLvl       <- sapply( levels(dds$Type), function(lvl) rowMeans( counts(dds,normalized=TRUE)[,dds$Type == lvl] ) )
	res          <- results(dds) 
	res          <- res[order(abs(res$log2FoldChange), decreasing=T),]
	res          <- cbind(res, log2(perLvl[rownames(res),]))
	res          <- cbind(res, annot[rownames(res),])
	res          <- res[which(res$padj < pval),]
	
	return(res)
	}

MPP  <- getDiffExp(raw, pheno, which(pheno$Group== "MPP"), annot, 0.05)
ST  <- getDiffExp(raw, pheno, which(pheno$Group == "ST"), annot, 0.05)
LT <- getDiffExp(raw, pheno, which(pheno$Group == "LT"), annot, 0.05)
CD63 <-getDiffExp(raw, pheno, which(pheno$Group == "CD63"), annot, 0.05)
#correct for control level
colnames(CD63) <- gsub("trans", "CD63", colnames(CD63))
colnames(LT) <- gsub("trans", "LT", colnames(LT))
colnames(ST) <- gsub("trans", "ST", colnames(ST))
colnames(MPP) <- gsub("trans", "MPP", colnames(MPP))

# LT,ST and MPP venn diagram
all_genes <- unique(c(rownames(LT), rownames(ST), rownames(MPP)))
all_genes <- matrix(nrow=length(all_genes), ncol=3, 0, dimnames=list(all_genes, c("LT", "ST", "MPP")))
all_genes[rownames(ST), "ST"] <- sign(ST$log2FoldChange)
all_genes[rownames(LT), "LT"] <- sign(LT$log2FoldChange)
all_genes[rownames(MPP), "MPP"] <- sign(MPP$log2FoldChange)
all_genes <- as.data.frame(all_genes)

getVennClass <- function(x2){x2[,1]+(x2[,2]*2)+(x2[,3]*4)}
upp <- rownames(all_genes)[unlist(lapply(1:nrow(all_genes), function(x){all(all_genes[x,] >= 0)}))]
upp <- setNames(getVennClass(all_genes[upp,]), upp)
dwn <- rownames(all_genes)[unlist(lapply(1:nrow(all_genes), function(x){all(all_genes[x,] <= 0)}))]
dwn <- setNames(getVennClass(abs(all_genes[dwn,])), dwn)
bth <- setNames(getVennClass(abs(all_genes)), rownames(all_genes))

all_genes$up   <- 0
all_genes$down <- 0
all_genes$both <- 0  
all_genes[names(upp), "up"] = as.integer(upp)
all_genes[names(dwn), "down"] = as.integer(dwn)
all_genes[names(bth), "both"] = as.integer(bth)
all_genes <- cbind(all_genes, "SYMBOL"=as.character(annot[rownames(all_genes),"SYMBOL"]))
all_genes <- cbind(all_genes, "ENTREZID"=as.character(annot[rownames(all_genes),"ENTREZID"]))
all_genes <- cbind(all_genes, "DESCRIPTION"=as.character(annot[rownames(all_genes),"DESCRIPTION"]))

#remove label sectors 
all_genes <- all_genes[, 1:3]
vennDiagram(all_genes, "up", main="Up")
vennDiagram(all_genes, "down", main="Down")
vennDiagram(all_genes, main="All")

#volcano plot CD63  x=logFC,y=-log10(pvalue)
vCD63 <-getDiffExp(raw, pheno, which(pheno$Group == "CD63"), annot, 1)
temp    <- vCD63
vCD63   <- vCD63[order(vCD63$padj),]
topUP   <-vCD63[vCD63$log2FoldChange > 0, "SYMBOL"][1:10]
topBOT  <-vCD63[vCD63$log2FoldChange < 0, "SYMBOL"][1:10]
genes <- vCD63[which(vCD63$padj < 0.05), "SYMBOL"]
found <- which(temp$SYMBOL %in% genes)
vCD63  <- temp 
xx   <- vCD63$log2FoldChange 
yy   <- -log10(vCD63$pvalue)
cols <- rep("black", nrow(vCD63))
sig <- which(vCD63$pvalue < 0.01)
cols[sig] <- "red" 
par(xpd=NA)
plot(xx,yy, pch=20, col=cols, xlab="logFC", ylab="-log10(pvalue)", main="CD63 vs control")
text(xx[found], yy[found], as.character(vCD63[found, "SYMBOL"]), pos=2)
points(xx[found], yy[found], cex=1.5, col="red", pch=20)
points(xx[found], yy[found], cex=1.5, col="black")




#######################################################################################
###########               Single Cell data                                 ############
###########       See scRNA-seq folder in Github                           ############
#######################################################################################


library(Seurat)
library(tidyverse)
library(openxlsx)
library(RColorBrewer)
library(ggExtra)
library(stringr)
library(this.path)

andreaPalete <- c("#ebac23", "#b80058", "#008cf9", "#006e00", "#00bbad", "#d163e6", "#b24502", "#ff9287", "#5954d6", "#00c6f8", "#878500", "#00a76c", "#bdbdbd", "#846b54")
color.list <- c(brewer.pal(12, "Paired"),brewer.pal(12, "Set3"),brewer.pal(8, "Pastel2"),colorRampPalette(c("grey20","grey70"))(4))
to <- color.list[22]
ss <- color.list[18]
color.list[22] <- ss
color.list[18] <- to

color.list <- c(andreaPalete,color.list)

#Change here the path to the analysis folder
analysisPath <- dirname(this.path())
setwd(analysisPath)

#######################################################################################
##############       Myeloid cells (Fig. 5 and S6)                       ##############
#######################################################################################


# This section will generate the plots with the final integrated analysis
bySample.cca <- readRDS("MyeloidFinal.Integrated.seurat.final.rds")

pList <- list()

# Defined cell types
# ImmGen classification by SingleR
pList[[length(pList)+1]] <- DimPlot(bySample.cca,
                                    reduction = "tsne",
                                    group.by = "immGenFTLabels",
                                    label = T, cols = color.list)

# Cell type redefinition based on Neutrophil Markers
DefaultAssay(bySample.cca) <- "RNA"
df <- Embeddings(bySample.cca,reduction = "tsne")
dnames <- colnames(df)
df <- cbind(df,FetchData(bySample.cca,vars = c("Mmp8", "Mpo", "S100a8", "S100a9")))

pList[[length(pList)+1]] <- df %>%
  pivot_longer(cols=-contains("tSNE"),names_to = "Gene", values_to="Expression") %>%
  arrange(Expression) %>%
  ggplot(aes_string(x=dnames[1],y=dnames[2])) +
  geom_point(aes_string(color="Expression")) +
  facet_wrap("Gene") +
  scale_color_gradientn(colors = colorRampPalette(c("grey","orange","red"))(3),
                        name="Log(NormCounts)") +
  theme_classic()

#put here the plot from 2.1.3

# tSNE with final cell type annotation 
pList[[length(pList)+1]] <- DimPlot(bySample.cca,
                                    reduction = "tsne",
                                    group.by = "celltype_cluster_reval",
                                    label = T, cols = color.list)

# The final analysis carries a ManualClustering which split cluster 9 in two
clustering <- "ManualClustering"

# Neutrophils Module Maximum enrichment plot 4
pList[[length(pList)+1]] <- bySample.cca@meta.data %>%
  filter(SampleName == "CD11b_WT" &
           ManualClustering %in% c("C0","C1","C2","C3","C5","C6","C8","C9","C9b","C10")) %>%
  dplyr::select(contains(c(clustering,"G0","G1","G2","G3","G4"))) %>%
  pivot_longer(cols = contains("RNAModule1"),
               names_to = "Module",values_to = "Score") %>%
  rename(Cluster=clustering) %>%
  group_by(Module) %>%
  summarise(Score=scale(Score,center = T,scale = T),
            Cluster=Cluster) %>%
  group_by(Cluster,Module) %>%
  summarise(AverageScore=mean(Score)) %>%
  group_by(Cluster) %>%
  summarise(Module=Module,
            AverageScoreScaled=scale(AverageScore),
            MaxEnrichment=(AverageScore >= (max(AverageScore)-0.02))) %>%
  mutate(Module=sub(".RNAModule1","",Module)) %>%
  mutate(plotBorder=ifelse(MaxEnrichment,1.5,0)) %>%
  ggplot() +
  geom_point(aes_string(x="Module",y="Cluster",
                        fill="AverageScoreScaled",
                        color="MaxEnrichment",
                        stroke = "plotBorder"),size=6,shape=21) +
  scale_fill_gradient(low = "grey90",high = "blue") +
  scale_color_manual(values = c(NA,"red")) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1))

# Monocyte Module Enrichment plot
clustering <- "ManualClustering"
pList[[length(pList)+1]] <- bySample.cca@meta.data %>%
  filter(SampleName == "CD11b_WT" &
           !ManualClustering %in% c("C0","C1","C2","C3","C5","C6","C8","C9","C9b","C10")) %>%
  dplyr::select(contains(c(clustering,"Module1"))) %>%
  dplyr::select(-contains(c("G0","G1","G2","G3","G4"))) %>%
  pivot_longer(cols = contains("RNAModule1"),
               names_to = "Module",values_to = "Score") %>%
  rename(Cluster=clustering) %>%
  group_by(Module) %>%
  summarise(Score=scale(Score,center = T,scale = T),
            Cluster=Cluster) %>%
  group_by(Cluster,Module) %>%
  summarise(AverageScore=mean(Score)) %>%
  group_by(Cluster) %>%
  summarise(Module=Module,
            AverageScoreScaled=scale(AverageScore),
            MaxEnrichment=(AverageScore >= (max(AverageScore)-0.02))) %>%
  mutate(Module=sub(".RNAModule1","",Module)) %>%
  mutate(plotBorder=ifelse(MaxEnrichment,1.5,0)) %>%
  ggplot() +
  geom_point(aes_string(x="Module",y="Cluster",
                        fill="AverageScoreScaled",
                        color="MaxEnrichment",
                        stroke = "plotBorder"),size=6,shape=21) +
  scale_fill_gradient(low = "grey90",high = "blue") +
  scale_color_manual(values = c(NA,"red")) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1))

# Cell Type by Condition
pList[[length(pList)+1]] <- DimPlot(bySample.cca,reduction = "tsne",
                                    group.by = "MyeloidType_ManualClustering_002",
                                    cols = color.list, split.by = "SampleName")

# Cell Type Proportion By Condition
pList[[length(pList)+1]] <- bySample.cca@meta.data %>%
  ggplot(aes(x="",fill=MyeloidType_ManualClustering_002)) +
  geom_bar(stat="count",position="fill") +
  coord_polar("y",start = 0) +
  scale_fill_manual(values = color.list) +
  facet_wrap("SampleName") +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.background = element_blank())

Idents(bySample.cca) <- "MyeloidType_ManualClustering_002"

# IPA Plots
maxzs <- 3
minzs <- -3
ipaList <- list()
for (cluster in levels(Idents(bySample.cca))) {
  if (cluster == "G1") {clusterf <- "G1b"} else {clusterf <- cluster}
  mappedGenes <- read.delim(
    file.path(paste0("ContrastBySampleName"),"MyeloidFinal",
              "KOvsWT","IPA",
              paste0(clusterf,".IPA.Mappings.txt")),
    sep = "\t",header = T)
  mappedGenes$X <- NULL
  mappedGenes$cluster <- cluster
  daf <- read.delim(file = file.path(paste0("ContrastBySampleName"),"MyeloidFinal",
                                     "KOvsWT","IPA",
                                     paste0(clusterf,".IPA.DAF.txt")),
                    sep = "\t",header = T)
  daf$X <- NULL
  daf$cluster <- cluster
  
  dafLong <- daf %>%
    separate_rows(Molecules,sep = ",") %>%
    left_join(mappedGenes,by = c("Molecules"="Symbol","cluster"="cluster"))
  
  mySelectedClusters <- getSheetNames(file.path(paste0("ContrastBySampleName"),"MyeloidFinal",
                                                "KOvsWT","IPA",
                                                paste0("Selected_DAF_Myeloid.xlsx")))
  
  if (paste(cluster,"DAF") %in% mySelectedClusters) {
    dafSelected <- read.xlsx(file.path(paste0("ContrastBySampleName"),"MyeloidFinal",
                                       "KOvsWT","IPA",
                                       paste0("Selected_DAF_Myeloid.xlsx")),
                             sheet = paste(cluster,"DAF"))
    dafSelected <- dafSelected %>% filter(Selected == 1)
    
    myDAFPlot <- dafLong %>%
      filter(Diseases.or.Functions.Annotation %in% dafSelected$Diseases.or.Functions.Annotation) %>%
      group_by(Diseases.or.Functions.Annotation) %>%
      summarise(
        B.H.p.value = -log10(mean(B.H.p.value)),
        zscore = mean(`Activation.z.score`),
        n=n()) %>%
      mutate(cc=ifelse(zscore>0,"r","b"))
    
    maxzs <- ifelse(maxzs < max(myDAFPlot$zscore,na.rm = T),max(myDAFPlot$zscore,na.rm = T),maxzs)
    minzs <- ifelse(minzs > min(myDAFPlot$zscore,na.rm = T),min(myDAFPlot$zscore,na.rm = T),minzs)
    
    myDAFPlot <- myDAFPlot %>%
      arrange(B.H.p.value)
    ipaList[[length(ipaList)+1]] <- myDAFPlot %>%
      arrange(B.H.p.value) %>%
      mutate(Diseases.or.Functions.Annotation=factor(Diseases.or.Functions.Annotation,
                                                     levels = myDAFPlot$Diseases.or.Functions.Annotation)) %>%
      ggplot(aes(y=B.H.p.value,x=Diseases.or.Functions.Annotation,fill=zscore)) +
      geom_bar(stat = "identity") +
      geom_text(color = "black",
                aes(x = Diseases.or.Functions.Annotation, y = B.H.p.value, label = n),
                position = position_stack(0.5)) +
      coord_flip() +
      ylab("-Log10(BH p-value)") +
      ggtitle(paste(cluster,"by BH P-Value"))
  }
}

ipaList <- lapply(ipaList, function (x) {
  return(x  +
           scale_fill_continuous( low=color.list[16],
                                  high = color.list[20],
                                  limits=c(minzs, maxzs),
                                  breaks=seq(round(minzs),round(maxzs),by=1)
           ) +
           theme_classic()
  )
})

pList <- c(pList,ipaList)

DefaultAssay(bySample.cca) <- "RNA"

### Il1B and Il1rn Expression
df <- Embeddings(bySample.cca,reduction = "tsne")
dnames <- colnames(df)
df <- cbind(df,bySample.cca@meta.data)
df <- cbind(df,FetchData(bySample.cca,vars = c("Il1b","Il1rn")))

pList[[length(pList)+1]] <- df %>%
  arrange(Il1b) %>%
  ggplot(aes_string(x=dnames[1],y=dnames[2])) +
  geom_point(aes_string(color="Il1b")) +
  facet_wrap("SampleName") +
  ggtitle("Il1b") +
  scale_color_gradientn(colors = colorRampPalette(c("grey","orange","red"))(3),
                        name="Log(NormCounts)") +
  theme_classic()

pList[[length(pList)+1]] <- df %>%
  filter(SampleName == "CD11b_WT")  %>%
  arrange(Il1rn) %>%
  ggplot(aes_string(x=dnames[1],y=dnames[2])) +
  geom_point(aes_string(color="Il1rn")) +
  facet_wrap("SampleName") +
  ggtitle("Il1rn") +
  scale_color_gradientn(colors = colorRampPalette(c("grey","orange","red"))(3),
                        name="Log(NormCounts)") +
  theme_classic()


pList[[length(pList)+1]] <- df %>%
  filter(SampleName == "CD11b_WT" & !(MyeloidType_ManualClustering_002 %in% c("HSC", "PreDC_II", "CMoP_I"))) %>%
  mutate(Il1Exp = case_when(
    Il1b>0 & Il1rn == 0 ~ "Il1b only",
    Il1b==0 & Il1rn > 0 ~ "Il1rn only",
    Il1b>0 & Il1rn > 0 ~ "Il1b and Il1rn")) %>%
  ggplot(aes(x="",fill=Il1Exp)) +
  geom_bar(stat="count",position="fill") +
  coord_polar("y",start = 0) +
  scale_fill_manual(values = color.list) +
  facet_wrap(c("SampleName","MyeloidType_ManualClustering_002"),ncol = 6) +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.background = element_blank())


# Violin plots on Selected Genes

selectedGenes <- c("Lyn","Hif1a","Lmo4","Csf2rb","Myd88","Cxcr2","Nfkbia","Cebpb")

df <- FetchData(bySample.cca,vars = selectedGenes)
df <- cbind(df,bySample.cca@meta.data)
df <- df %>% mutate(SampleName=factor(SampleName,levels = c("CD11b_WT","CD11b_KO"))) %>%
  filter(!MyeloidType_ManualClustering_002 %in% c("HSC","PreDC_II")) %>%
  mutate(MyeloidType_ManualClustering_002 = factor(MyeloidType_ManualClustering_002,
                                                   levels = c("G0","G1","G2","G3","G4","CMoP_I","BMM_I","BMM_II")))

for (gene in selectedGenes) {
  pList[[length(pList)+1]] <- ggplot(df,aes_string(y=gene,x="SampleName",
                                                   color="SampleName",
                                                   fill="SampleName")) +
    geom_violin(stat = "ydensity",alpha=0.2) +
    scale_color_manual(values = c("grey40","#a80a0b")) +
    scale_fill_manual(values = c("grey40","#a80a0b")) +
    stat_summary(fun = "mean",
                 geom = "crossbar", 
                 width = 0.5,
                 colour = "black") +
    ylab("Expression Level log(NormCounts+1)") +
    ggtitle(gene) +
    facet_wrap("MyeloidType_ManualClustering_002",nrow = 1) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90),
          legend.position = "none")
}

dir.create("Plots")
p <- lapply(seq(length(pList)), function (x) {
  ggsave(file.path("Plots",paste0("AllMyeloidPlots.",x,".pdf")),
         plot = pList[[x]],
         device = pdf())
})



#######################################################################################
##############       LSK cells (Fig. 6 and S8)                           ##############
#######################################################################################

# This section will generate the plots with the final integrated analysis
# The final analysis carries a ManualClustering which split cluster 9 in two
bySample.cca <- readRDS("LSKFinal.Integrated.seurat.rds")
clustering <- "integrated_snn_res.0.5"

pList <- list()

# Module Maximum enrichment plot
clustering <- "integrated_snn_res.0.5"
pList[[length(pList)+1]] <- bySample.cca@meta.data %>%
  filter(SampleName == "LSK_WT") %>%
  dplyr::select(contains(c(clustering,"RNAModule1"))) %>%
  pivot_longer(cols = contains("RNAModule1"),names_to = "Module",values_to = "Score") %>%
  rename(Cluster=clustering) %>%
  group_by(Module) %>%
  summarise(Score=scale(Score,center = T,scale = T),Cluster=Cluster) %>%
  group_by(Cluster,Module) %>%
  summarise(AverageScore=mean(Score)) %>%
  group_by(Cluster) %>%
  summarise(Module=Module,
            AverageScoreScaled=scale(AverageScore),
            MaxEnrichment=(AverageScore >= (max(AverageScore)-0.02))) %>%
  mutate(Module=sub(".RNAModule1","",Module)) %>%
  mutate(plotBorder=ifelse(MaxEnrichment,1.5,0)) %>%
  ggplot() +
  geom_point(aes_string(x="Module",y="Cluster",
                        fill="AverageScoreScaled",
                        color="MaxEnrichment",
                        stroke = "plotBorder"),
             size=6,shape=21) +
  scale_fill_gradient(low = "grey90",high = "blue") +
  scale_color_manual(values = c(NA,"red"))

# Cell Type by Condition
pList[[length(pList)+1]] <- DimPlot(bySample.cca,reduction = "tsne",
                                    group.by = "SommerkampType_05_002",
                                    cols = color.list, split.by = "SampleName")

# Cell Type Proportion By Condition
pList[[length(pList)+1]] <- bySample.cca@meta.data %>%
  ggplot(aes(x="",fill=SommerkampType_05_002)) +
  geom_bar(stat="count",position="fill") +
  coord_polar("y",start = 0) +
  scale_fill_manual(values = color.list) +
  facet_wrap("SampleName") +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.background = element_blank())

Idents(bySample.cca) <- "SommerkampType_05_002"

# IPA Plots
maxzs <- 3
minzs <- -3
ipaList <- list()
for (cluster in levels(Idents(bySample.cca))) {
  clusterf <- cluster
  mappedGenes <- read.delim(
    file.path(paste0("ContrastBySampleName"),"LSKFinal","KOvsWT","IPA",
              paste0(clusterf,".IPA.Mappings.txt")),
    sep = "\t",header = T)
  mappedGenes$X <- NULL
  mappedGenes$cluster <- cluster
  daf <- read.delim(file = file.path(paste0("ContrastBySampleName"),"LSKFinal",
                                     "KOvsWT","IPA",paste0(clusterf,".IPA.DAF.txt")),
                    sep = "\t",header = T)
  daf$X <- NULL
  daf$cluster <- cluster
  
  dafLong <- daf %>%
    separate_rows(Molecules,sep = ",") %>%
    left_join(mappedGenes,by = c("Molecules"="Symbol","cluster"="cluster"))
  
  mySelectedClusters <- getSheetNames(file.path(paste0("ContrastBySampleName"),
                                                "LSKFinal","KOvsWT","IPA",
                                                paste0("Selected_DAF_LSK2.xlsx")))
  
  if (paste(cluster,"DAF") %in% mySelectedClusters) {
    dafSelected <- read.xlsx(file.path(paste0("ContrastBySampleName"),
                                       "LSKFinal","KOvsWT","IPA",
                                       paste0("Selected_DAF_LSK2.xlsx")),
                             sheet = paste(cluster,"DAF"))
    dafSelected <- dafSelected %>% filter(Selected == 1)
    
    myDAFPlot <- dafLong %>%
      filter(Diseases.or.Functions.Annotation %in% dafSelected$Diseases.or.Functions.Annotation) %>%
      group_by(Diseases.or.Functions.Annotation) %>%
      summarise(
        B.H.p.value = -log10(mean(B.H.p.value)),
        zscore = mean(`Activation.z.score`),
        n=n()) %>%
      mutate(cc=ifelse(zscore>0,"r","b"))
    
    maxzs <- ifelse(maxzs < max(myDAFPlot$zscore,na.rm = T),max(myDAFPlot$zscore,na.rm = T),maxzs)
    minzs <- ifelse(minzs > min(myDAFPlot$zscore,na.rm = T),min(myDAFPlot$zscore,na.rm = T),minzs)
    
    myDAFPlot <- myDAFPlot %>%
      arrange(B.H.p.value)
    ipaList[[length(ipaList)+1]] <- myDAFPlot %>%
      arrange(B.H.p.value) %>%
      mutate(Diseases.or.Functions.Annotation=factor(Diseases.or.Functions.Annotation,
                                                     levels = myDAFPlot$Diseases.or.Functions.Annotation)) %>%
      ggplot(aes(y=B.H.p.value,x=Diseases.or.Functions.Annotation,fill=zscore)) +
      geom_bar(stat = "identity") +
      geom_text(color = "black",
                aes(x = Diseases.or.Functions.Annotation, y = B.H.p.value, label = n),
                position = position_stack(0.5)) +
      coord_flip() +
      ylab("-Log10(BH p-value)") +
      ggtitle(paste(cluster,"by BH P-Value"))
  }
}

ipaList <- lapply(ipaList, function (x) {
  return(x  +
           scale_fill_continuous( low=color.list[16],
                                  high = color.list[20],
                                  limits=c(minzs, maxzs),
                                  breaks=seq(round(minzs),round(maxzs),by=1)
           ) +
           theme_classic()
  )
})

pList <- c(pList,ipaList)

DefaultAssay(bySample.cca) <- "RNA"

## Il1B and Il1rn Expression
df <- Embeddings(bySample.cca,reduction = "tsne")
dnames <- colnames(df)
df <- cbind(df,bySample.cca@meta.data)
df <- cbind(df,FetchData(bySample.cca,vars = c("IL1B","IL1RN")))

pList[[length(pList)+1]] <- df %>%
  arrange(IL1B) %>%
  ggplot(aes_string(x=dnames[1],y=dnames[2])) +
  geom_point(aes_string(color="IL1B")) +
  facet_wrap("SampleName") +
  ggtitle("Il1b") +
  scale_color_gradientn(colors = colorRampPalette(c("grey","orange","red"))(3),
                        name="Log(NormCounts)") +
  theme_classic()

pList[[length(pList)+1]] <- df %>%
  arrange(IL1RN) %>%
  filter(SampleName == "LSK_WT")  %>%
  ggplot(aes_string(x=dnames[1],y=dnames[2])) +
  geom_point(aes_string(color="IL1RN")) +
  facet_wrap("SampleName") +
  ggtitle("Il1rn") +
  scale_color_gradientn(colors = colorRampPalette(c("grey","orange","red"))(3),
                        name="Log(NormCounts)") +
  theme_classic()



# Violin plots on Selected Genes
selectedGenes <- c("IFITM1","NMT1","CRIP1","IFITM3","CD52","JUNB","HLF")

df <- FetchData(bySample.cca,vars = selectedGenes)
df <- cbind(df,bySample.cca@meta.data)
df <- df %>% mutate(SampleName=factor(SampleName,levels = c("LSK_WT","LSK_KO")))

for (gene in selectedGenes) {
  pList[[length(pList)+1]] <- ggplot(df,aes_string(y=gene,
                                                   x="SampleName",
                                                   color="SampleName",
                                                   fill="SampleName")) +
    geom_violin(stat = "ydensity",alpha=0.2) +
    scale_color_manual(values = c("grey40","#a80a0b")) +
    scale_fill_manual(values = c("grey40","#a80a0b")) +
    stat_summary(fun = "mean",
                 geom = "crossbar", 
                 width = 0.5,
                 colour = "black") +
    ylab("Expression Level log(NormCounts+1)") +
    ggtitle(str_to_title(gene)) +
    facet_wrap("SommerkampType_05_002",nrow = 1) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90),legend.position = "none")
}

dir.create("Plots")
p <- lapply(seq(length(pList)), function (x) {
  ggsave(file.path("Plots",paste0("AllLSKPlots.",x,".pdf")),
         plot = pList[[x]],
         device = pdf())
})

#######################################################################################
##############       Stromal cells (Fig. 7 and S10)                      ##############
#######################################################################################

# This section will generate the plots with the final integrated analysis
bySample.cca <- readRDS("StromalFinal.IntegratedFilteredDoublets.seurat.rds")
clustering <- "integrated_snn_res.0.1"

pList <- list()

# Module Maximum enrichment plot
clustering <- "integrated_snn_res.0.1"
pList[[length(pList)+1]] <- bySample.cca@meta.data %>%
  filter(SampleName == "CD63_WT") %>%
  dplyr::select(contains(c(clustering,"RNAModule1"))) %>%
  pivot_longer(cols = contains("RNAModule1"),names_to = "Module",values_to = "Score") %>%
  rename(Cluster=clustering) %>%
  group_by(Module) %>%
  summarise(Score=scale(Score,center = T,scale = T),Cluster=Cluster) %>%
  group_by(Cluster,Module) %>%
  summarise(AverageScore=mean(Score)) %>%
  group_by(Cluster) %>%
  summarise(Module=Module,AverageScoreScaled=scale(AverageScore),MaxEnrichment=(AverageScore >= (max(AverageScore)-0.02))) %>%
  mutate(Module=sub(".RNAModule1","",Module)) %>%
  mutate(plotBorder=ifelse(MaxEnrichment,1.5,0)) %>%
  ggplot() +
  geom_point(aes_string(x="Module",y="Cluster",
                        fill="AverageScoreScaled",
                        color="MaxEnrichment", 
                        stroke = "plotBorder"),size=6,shape=21) +
  scale_fill_gradient(low = "grey90",high = "blue") +
  scale_color_manual(values = c(NA,"red"))

# Cell Type by Condition
pList[[length(pList)+1]] <- DimPlot(bySample.cca,reduction = "tsne",
                                    group.by = "ScadenType_01_002",
                                    cols = color.list, split.by = "SampleName")

# Cell Type Proportion By Condition
pList[[length(pList)+1]] <- bySample.cca@meta.data %>%
  ggplot(aes(x="",fill=ScadenType_01_002)) +
  geom_bar(stat="count",position="fill") +
  coord_polar("y",start = 0) +
  scale_fill_manual(values = color.list) +
  facet_wrap("SampleName") +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.background = element_blank())

Idents(bySample.cca) <- "ScadenType_01_002"


# IPA Plots
maxzs <- 3
minzs <- -3
ipaList <- list()
for (cluster in levels(Idents(bySample.cca))) {
  
  mappedGenes <- read.delim(file.path("ContrastBySampleName","StromalFinal",
                                      "KOvsWT","IPA",paste0(cluster,".IPA.Mappings.txt")),sep = "\t",header = T)
  mappedGenes$X <- NULL
  mappedGenes$cluster <- cluster
  daf <- read.delim(file = file.path("ContrastBySampleName","StromalFinal",
                                     "KOvsWT","IPA",paste0(cluster,".IPA.DAF.txt")),sep = "\t",header = T)
  daf$X <- NULL
  daf$cluster <- cluster
  
  dafLong <- daf %>%
    separate_rows(Molecules,sep = ",") %>%
    left_join(mappedGenes,by = c("Molecules"="Symbol","cluster"="cluster"))
  
  mySelectedClusters <- getSheetNames(file.path("ContrastBySampleName","StromalFinal",
                                                "KOvsWT","IPA",
                                                paste0("Selected_DAF_Stromal.xlsx")))
  
  dafSelected <- read.xlsx(file.path("ContrastBySampleName","StromalFinal",
                                     "KOvsWT","IPA",paste0("Selected_DAF_Stromal.xlsx")),
                           sheet = paste(cluster,"DAF"))
  dafSelected <- dafSelected %>% filter(Selected == 1)
  
  myDAFPlot <- dafLong %>%
    filter(Diseases.or.Functions.Annotation %in% dafSelected$Diseases.or.Functions.Annotation) %>%
    group_by(Diseases.or.Functions.Annotation) %>%
    summarise(
      B.H.p.value = -log10(mean(B.H.p.value)),
      zscore = mean(`Activation.z.score`),
      n=n()) %>%
    mutate(cc=ifelse(zscore>0,"r","b"))
  
  maxzs <- ifelse(maxzs < max(myDAFPlot$zscore,na.rm = T),max(myDAFPlot$zscore,na.rm = T),maxzs)
  minzs <- ifelse(minzs > min(myDAFPlot$zscore,na.rm = T),min(myDAFPlot$zscore,na.rm = T),minzs)
  
  myDAFPlot <- myDAFPlot %>%
    arrange(B.H.p.value)
  ipaList[[length(ipaList)+1]] <- myDAFPlot %>%
    arrange(B.H.p.value) %>%
    mutate(Diseases.or.Functions.Annotation=factor(Diseases.or.Functions.Annotation,
                                                   levels = myDAFPlot$Diseases.or.Functions.Annotation)) %>%
    ggplot(aes(y=B.H.p.value,x=Diseases.or.Functions.Annotation,fill=zscore)) +
    geom_bar(stat = "identity") +
    geom_text(color = "black",
              aes(x = Diseases.or.Functions.Annotation, y = B.H.p.value, label = n),
              position = position_stack(0.5)) +
    coord_flip() +
    ylab("-Log10(BH p-value)") +
    ggtitle(paste(cluster,"by BH P-Value"))
  
}

ipaList <- lapply(ipaList, function (x) {
  return(x  +
           scale_fill_continuous( low=color.list[16],
                                  high = color.list[20],
                                  limits=c(minzs, maxzs),
                                  breaks=seq(round(minzs),round(maxzs),by=1)
           ) +
           theme_classic()
  )
})

pList <- c(pList,ipaList)

DefaultAssay(bySample.cca) <- "RNA"

## Il1B and Il1rn Expression
df <- Embeddings(bySample.cca,reduction = "tsne")
dnames <- colnames(df)
df <- cbind(df,bySample.cca@meta.data)
df <- cbind(df,FetchData(bySample.cca,vars = c("Il1b","Il1rn")))

pList[[length(pList)+1]] <- df %>%
  arrange(Il1b) %>%
  ggplot(aes_string(x=dnames[1],y=dnames[2])) +
  geom_point(aes_string(color="Il1b")) +
  facet_wrap("SampleName") +
  ggtitle("IL1B") +
  scale_color_gradientn(colors = colorRampPalette(c("grey","orange","red"))(3),
                        name="Log(NormCounts)") +
  theme_classic()

pList[[length(pList)+1]] <- df %>%
  filter(SampleName == "CD63_WT")  %>%
  arrange(Il1rn) %>%
  ggplot(aes_string(x=dnames[1],y=dnames[2])) +
  geom_point(aes_string(color="Il1rn")) +
  facet_wrap("SampleName") +
  ggtitle("IL1RN") +
  scale_color_gradientn(colors = colorRampPalette(c("grey","orange","red"))(3),
                        name="Log(NormCounts)") +
  theme_classic()

pList[[length(pList)+1]] <- df %>%
  filter(SampleName == "CD63_WT")  %>%
  mutate(Il1Exp = case_when(
    Il1b>0 & Il1rn == 0 ~ "Il1b only",
    Il1b==0 & Il1rn > 0 ~ "Il1rn only",
    Il1b>0 & Il1rn > 0 ~ "Il1b and Il1rn")) %>%
  ggplot(aes(x="",fill=Il1Exp)) +
  geom_bar(stat="count",position="fill") +
  coord_polar("y",start = 0) +
  scale_fill_manual(values = color.list) +
  facet_wrap(c("SampleName","ScadenType_01_002"),ncol = 6) +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.background = element_blank())


# Violin plots on Selected Genes
selectedGenes <- list()
selectedGenes[["Fb_MSC"]] <- c("Nr1d1","Zbtb16","Nfia","Hspg2","Bmp4",
                               "Col3a1","Fbn1","Twist1","Spp1")
selectedGenes[["Sinu"]] <- sort(unique(c("Il6","Il1b","Csf1","Il6",
                                         "Il1b","Osm","Ccl9","Ccl6",
                                         "Nlrp3","Csf1","Ccl4","Il6",
                                         "Il1b","Osm","Ccl9","Ccl6","Nlrp3")))

df <- FetchData(bySample.cca,vars = unlist(selectedGenes))
df <- cbind(df,bySample.cca@meta.data)
df <- df %>% mutate(SampleName=factor(SampleName,levels = c("CD63_WT","CD63_KO")))

## Fibroblasts and MSC
pList[[length(pList)+1]] <- df %>%
  filter(ScadenType_01_002 %in% c("Fibroblast")) %>%
  pivot_longer(cols = selectedGenes$Fb_MSC,names_to = "genes",values_to = "expression") %>%
  ggplot(aes_string(y="expression",x="SampleName",
                    color="SampleName",
                    fill="SampleName")) +
  geom_violin(stat = "ydensity",alpha=0.2) +
  scale_color_manual(values = c("grey40","#a80a0b")) +
  scale_fill_manual(values = c("grey40","#a80a0b")) +
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.5,
               colour = "black") +
  ylab("Expression Level log(NormCounts+1)") +
  ggtitle("Fibroblast") +
  facet_wrap("genes",nrow = 1) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90),legend.position = "none")

pList[[length(pList)+1]] <- df %>%
  filter(ScadenType_01_002 %in% c("MSCs")) %>%
  pivot_longer(cols = selectedGenes$Fb_MSC,names_to = "genes",values_to = "expression") %>%
  ggplot(aes_string(y="expression",x="SampleName",
                    color="SampleName",
                    fill="SampleName")) +
  geom_violin(stat = "ydensity",alpha=0.2) +
  scale_color_manual(values = c("grey40","#a80a0b")) +
  scale_fill_manual(values = c("grey40","#a80a0b")) +
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.5,
               colour = "black") +
  ylab("Expression Level log(NormCounts+1)") +
  ggtitle("MSCs") +
  facet_wrap("genes",nrow = 1) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90),legend.position = "none")

## Sinusoidal
pList[[length(pList)+1]] <- df %>%
  filter(ScadenType_01_002 %in% c("Sinusoidal")) %>%
  pivot_longer(cols = selectedGenes$Sinu,names_to = "genes",values_to = "expression") %>%
  ggplot(aes_string(y="expression",x="SampleName",
                    color="SampleName",
                    fill="SampleName")) +
  geom_violin(stat = "ydensity",alpha=0.2) +
  scale_color_manual(values = c("grey40","#a80a0b")) +
  scale_fill_manual(values = c("grey40","#a80a0b")) +
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.5,
               colour = "black") +
  ylab("Expression Level log(NormCounts+1)") +
  ggtitle("Sinusoidal") +
  facet_wrap("genes",nrow = 1) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90),legend.position = "none")


dir.create("Plots")
p <- lapply(seq(length(pList)), function (x) {
  ggsave(file.path("Plots",paste0("AllStromalPlots.",x,".pdf")),
         plot = pList[[x]],
         device = pdf())
})


