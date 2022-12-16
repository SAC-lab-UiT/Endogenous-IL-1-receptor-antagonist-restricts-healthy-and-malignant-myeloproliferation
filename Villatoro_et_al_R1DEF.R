pkgs <- c('clusterProfiler','data.table','DESeq2','GEOquery','ggplot2','ggpubr','ggrepel','ggsci','ggsignif','ggthemes','grid','gridExtra','limma','org.Mm.eg.db','plotly','ReactomePA','splitstackshape','stringr','survival','survminer')
for(i in 1:length(pkgs)) {
  library(pkgs[i], character.only = T)
}

########################################################################################
###        Survival AML patients Verhaak dataset (Fig. 1A)                            ##
###                                                                                   ##
########################################################################################

setwd( "E:/service/arranz")
survOrig <- read.table(file = 'data/survData.txt',sep = '\t',head = T,stringsAsFactors = F)
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
keepM4M5 <- which(modelFrame$FAB %in% c('M0M3','M4M5') & modelFrame$risk != 'unknown')


####Global Model 

x = as.formula("Surv(OS, dead) ~ IL1RN + risk + age")
xNull <- as.formula("Surv(OS, dead) ~ risk + age")
resNull <- suppressWarnings(coxph(xNull, data = modelFrame[keep,]))
res.cox <- suppressWarnings(coxph(x, data = modelFrame[keep,]))
x <- summary(res.cox)
anovaRes <- anova(resNull,resCox)

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

ggsave(
  myplots_out,
  file = "output/SurvivalGlobal2.pdf",
  width = 5,
  height = 7)

M4M5plotData <- summary(fit2)                           
M4M5plotData <- data.frame(time = M4M5plotData$time,group = M4M5plotData$strata,surv = M4M5plotData$surv)


write.table(globalPlotData,file = 'output/globalCurve.txt;,sep = '\t')
write.table(M4M5plotData,file = 'output/M4M5Curve.txt;,sep = '\t')

#######################################################################################
############                  19 matched-pair AML patients                 ############
############                      Survival   (Fig. S1B)                    ############
############                                                               ############
#######################################################################################

url="ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE83nnn/GSE83533/suppl/GSE83533%5FDx%5FRel%5FAML%5FRNAseq%5Frpkm%2Etxt%2Egz"
GSE83533_exprs <- data.table::fread(paste(paste('URL=',url, sep = ""), ' ; curl "$URL" | gunzip -c', sep = ""))

IL1RN_exprs <- GSE83533_exprs[GSE83533_exprs$gene == "IL1RN", ]
IL1RN_exprs <- IL1RN_exprs[, 3:ncol(IL1RN_exprs)]
IL1RN_exprs <- as.data.frame(t(IL1RN_exprs))
IL1RN_exprs$sample <- row.names(IL1RN_exprs)
IL1RN_exprs$type <- ifelse(grepl("Rel", IL1RN_exprs$sample), "Relapse", "Diagnosis")
IL1RN_exprs$patient <- gsub("AML_|_Dx|_Rel", "", IL1RN_exprs$sample)
IL1RN_exprs <- IL1RN_exprs[order(IL1RN_exprs$patient), ]

# SraRunTable from study phs001027 using SRA Run Selector (https://www.ncbi.nlm.nih.gov/Traces/study/?) 
# https://www.ncbi.nlm.nih.gov/Traces/study/?acc=phs001027&f=analyte_type_sam_s%3An%3Arna%3Ac&s=SRR3088324,SRR3088328,SRR3088346,SRR3088348,SRR3088350,SRR3088354,SRR3088357,SRR3088361,SRR3088391,SRR3088395,SRR3088398,SRR3088400,SRR3088423,SRR3088427,SRR3088430,SRR3088432,SRR3088444,SRR3088448,SRR3088456,SRR3088457,SRR3088458,SRR3088459,SRR3088460,SRR3088461,SRR3088462,SRR3088463,SRR3088464,SRR3088465,SRR3088466,SRR3088467,SRR3088468,SRR3088469,SRR3088470,SRR3088471,SRR3088472,SRR3088473,SRR3088474,SRR3088475,SRR3088476,SRR3088477,SRR3088478,SRR3088479,SRR3088480,SRR3088481,SRR3088482,SRR3088483,SRR3088484,SRR3088485,SRR3088486,SRR3088487,SRR3088488,SRR3088489,SRR3088490,SRR3088491,SRR3088492,SRR3088493,SRR3088494,SRR3088495,SRR3088496,SRR3088497,SRR3088498,SRR3088499,SRR3088500,SRR3088501,SRR3088502,SRR3088503,SRR3088504,SRR3088505,SRR3088506,SRR3088507,SRR3088508,SRR3088509,SRR3088510,SRR3088511,SRR3088512,SRR3088513,SRR3088514,SRR3088515,SRR3088516,SRR3088517,SRR3088518,SRR3088519,SRR3088520,SRR3088521,SRR3088522,SRR3088523,SRR3088524,SRR3088525,SRR3088526,SRR3088527,SRR3088528,SRR3088529,SRR3088530,SRR3088531,SRR3088532,SRR3088533,SRR3088534,SRR3088535,SRR3088536,SRR3088537,SRR3088538,SRR3088539,SRR3088540,SRR3088541,SRR3088542,SRR3088543,SRR3088544,SRR3088545,SRR3088546,SRR3088547,SRR3088548,SRR3088549,SRR3088550,SRR3088551,SRR3088564,SRR3088568,SRR3088584,SRR3088588,SRR3088034,SRR3088036,SRR3088058,SRR3088060,SRR3088086,SRR3088088,SRR3088096,SRR3088098,SRR3088130,SRR3088132,SRR3088136,SRR3088138,SRR3088154,SRR3088156,SRR3088166,SRR3088168
# For simplicity, SraRunTable.txt is also included in repository. 
# This file needs to be download to the current working directory before commencing the next command
SraRunTable = fread("data/SraRunTable.txt", sep = "/")
pheno_subject <- fread("data/phs001027.v2.pht005216.v1.p1.c1.Epigenetics_AML_Subject_Phenotypes.GRU-PUB.txt", skip = "dbGaP_Subject_ID")

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
p <- p + scale_color_npg(labels = c("High IL1RN expression", "Low IL1RN expression")) + theme(legend.background = element_blank())
p <- p + theme(legend.direction = "vertical", legend.position = c(0.7,0.9), legend.title = element_blank())
p <- p + annotate("text", x = (5*365)+100, y = 0.69, label = "69%", size  = 6)
p + annotate("text", x = (5*365)+100, y = 0.1, label = "10%", size  = 6)



########################################################################################
###  NFkB activity in hematopoetic and progenitor cells (HSPCs) (Fig. 2O, Fig. 5G,    ##
###  Suppl. Fig. S7G)                                                                 ##
########################################################################################


#' @title estimateTFactivity 
#' @import ggplot2
#' @import ggpubr
#' @export
#' @param M VST-normalized counts data
#' @param targetGenes NFKB target genes read from upplementaryTableS4.txt
#' @param pheno data frame containing annotation information for each sample should contain columns treatment (WT , KO) and Group (ST LT MPP etc..)
#' @param annot data frame with gene annotation 
#' @param groups sample groups to include in analysis
#' @param treatments tretments to include in analysis
#' @param n number of replicates
#' @return An s3 object with the following slots \cr \cr
#' exportTable:        Estimated activity and pheno input \cr \cr
#' resTable            T-test results \cr \cr
#' cp:                 ggplot object \cr \cr
#' @examples outRas <- estimateTFactivity(M,targetGenes = nfkbTargets,pheno = pheno,annot = annot,groups = c('LT','ST','MPP'),treatments = c('WT','KO'))

#' 
#' See also \href{https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-6-225}{Tomfohr et al., 2005} which was used as an inspiration to generate this function


estimateTFactivity <- function(M,targetGenes,pheno,annot,
         groups = c('LT','ST','MPP'),
         treatments = c('WT','KO'),n = 3, var.equal = F){
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
        geom_boxplot(,fill = rep(c('blue','red'),3))
    return(list(test = resTable,exportTable = exportTable,cp = cp))
}        
    
##########################################################################################################################################
library(ggplot2)
library(tidyverse)
library(ggvenn)
library(ggrepel)
library(openxlsx)
library(patchwork)
library(RColorBrewer)
library(this.path)

color.list <- c(brewer.pal(12, "Paired"),brewer.pal(12, "Set3"),brewer.pal(8, "Pastel2"),colorRampPalette(c("grey20","grey70"))(4))
to <- color.list[22]
ss <- color.list[18]
color.list[22] <- ss
color.list[18] <- to

setwd(dirname(this.path()))

LT_Il1rnKO <- as.data.frame(read.xlsx(file.path("LT_Il1rnKO.nfkbAnnot.xlsx")))
ST_Il1rnKO <- as.data.frame(read.xlsx(file.path("ST_Il1rnKO.nfkbAnnot.xlsx")))
MPP_Il1rnKO <- as.data.frame(read.xlsx(file.path("MPP_Il1rnKO.nfkbAnnot.xlsx")))
CD63_Il1rnKO <- as.data.frame(read.xlsx(file.path("CD63_Il1KO_DESeq2_New.xlsx")))

objects <- c("LT_Il1rnKO","ST_Il1rnKO","MPP_Il1rnKO","CD63_Il1rnKO")


#######################################################################################
#######################       Volcanos in HSPCs (Fig. 2N)       #######################
#######################################################################################

# Volcano Plots with Nfkb targets highlighted
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
                    segment.color = 'transparent',
                    box.padding = 0.1,
                    max.overlaps = 500,
                    seed = 42
    ) +  # min.segment.length = Inf
    ggtitle(objects[i]) +
    theme(panel.background = element_blank(),
          axis.line = element_line(),
          panel.border = element_rect(fill = "transparent")
    ) +
    ylab("-log10 (Adjusted p value)") + xlab("log2 FC")
  
  assign(paste("p", objects[i], sep = "."), p)
  
}

pdf("Il1rnKO_LTSTMPP_Vulcanos.pdf")
p.LT_Il1rnKO + p.ST_Il1rnKO + p.MPP_Il1rnKO
dev.off()


#######################################################################################
########################            IL-1rn-KO                 #########################
########################       Venn in HSPCs (Fig. S3A)       #########################
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
########################      Volcano in CD63+ (Fig. 3F)      #########################
#######################################################################################

# CD63 Experiment Volcano Plot
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
GSE165810 <- read.xlsx("GSE165810.res.xlsx")
GSE165810 <- as.data.frame(GSE165810)
# GSE165810$external_gene_name <- rownames(GSE165810)
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
## YFP_POS
YFP_POS_Group <- read.xlsx("YFP_POS_Group.res.xlsx")
YFP_POS_Group <- as.data.frame(YFP_POS_Group)
# YFP_POS_Group$external_gene_name <- rownames(YFP_POS_Group)
YFP_POS_Group$padj2 <- YFP_POS_Group$padj

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
YFP_NEG_Group <- read.xlsx("YFP_NEG_Group.res.xlsx")
YFP_NEG_Group <- as.data.frame(YFP_NEG_Group)
# YFP_NEG_Group$external_gene_name <- rownames(YFP_NEG_Group)
YFP_NEG_Group$padj2 <- YFP_NEG_Group$padj

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
##############       Enrichment plot IL-1rn-KO LT-HSC (Fig. S3B,C)       ##############
##############       Enrichment plot IL-1rn-KO ST-HSC (Fig. S3D,E)       ##############
#######################################################################################

ridgeplot2 <- function(gsea, n.UP = 10, n.DOWN = 10, 
                       simplify = F, simplifyCutoff = 0.7, 
                       NES = 0.2, pval = 0.05, 
                       selected = NULL, dropTerm = NULL,
                       distanceSort = T, 
                       text.size = 10) {
    
    if(!is.null(selected)) {
        include = gsea@result[gsea@result$Description %in% selected, ]
        include$dist <- sqrt(include$NES^2 + include$p.adjust^2)
    }

    result <- gsea@result
    
    if(simplify) {
        simp <- clusterProfiler::simplify(gsea, simplifyCutoff)
        result <- result[result$Description %in% simp@result$Description, ]
    }
    
    result <- result[result$p.adjust < pval, ]
    
    if(distanceSort) {
        result$dist <- sqrt(result$NES^2 + result$p.adjust^2)
        result <- result[order(-result$dist), ]
    }
    
    if(!is.null(dropTerm)) result <- result[!result$Description %in% dropTerm, ]
    
    up = result[result$NES > NES, ]
    down = result[result$NES < -NES, ]
    
    if(nrow(up) < n.UP) n.UP = nrow(up)
    if(nrow(down) < n.DOWN) n.DOWN = nrow(down)
    
    if(!is.null(selected)) {
        n.UP = ifelse(n.UP < nrow(include[include$NES > NES, ]), 
                n.UP-(n.UP-nrow(include[include$NES > NES, ])), 
                    0)
        n.DOWN = ifelse(n.DOWN < nrow(include[include$NES < -NES, ]),
            n.DOWN-(n.DOWN-nrow(include[include$NES < -NES, ])),
            0)
    }
    
    up <- up[0:n.UP, ]
    down <- down[0:n.DOWN, ]
    
    if(!is.null(selected)) {
        up = rbind(up, include[include$NES > NES, ])
        down = rbind(down, include[include$NES < -NES, ])
    }

    for(i in 1:nrow(up)) {
        dfUP <- data.frame(geneID = unlist(stringr::str_split(up$core_enrichment[i], "/")))
        dfUP$FC <- gseaBP@geneList[match(dfUP$geneID, names(gseaBP@geneList))]
        dfUP$Description <- up$Description[i]
        dfUP$padj <- up$p.adjust[i]
        if(i == 1) allUP <- dfUP
        if(i != 1) allUP <- rbind(allUP, dfUP)
    }
    
    for(i in 1:nrow(down)) {
        dfDOWN <- data.frame(geneID = unlist(stringr::str_split(down$core_enrichment[i], "/")))
        dfDOWN$FC <- gseaBP@geneList[match(dfDOWN$geneID, names(gseaBP@geneList))]
        dfDOWN$Description <- down$Description[i]
        dfDOWN$padj <- down$p.adjust[i]
        if(i == 1) allDOWN <- dfDOWN
        if(i != 1) allDOWN <- rbind(allDOWN, dfDOWN)
    }
    
    all <- rbind(allUP, allDOWN)
    all <- all[!is.na(all$FC), ]
    
    all$Description <- factor(all$Description, levels = unique(all$Description))

    p <- suppressWarnings(ggplot(all, aes(FC, group = Description, fill = padj)) + geom_density() + scale_fill_continuous(low = "red", high = "blue") + 
        theme_few() + 
        facet_wrap(~Description, ncol = 1, switch =  "y", scales = "free_y") +
          theme(text = element_text(size = text.size), strip.text.y.left = element_text(angle = 0), strip.text.y = element_text(hjust = 1, size = text.size),
                strip.background = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(),
                axis.line = element_blank(), panel.border = element_blank(), axis.title.x = element_blank(), axis.ticks.y.left = element_blank(), 
                plot.margin = margin(0,0,0,0), panel.background = element_blank(), plot.background = element_blank(), strip.background.x = element_blank(), 
                strip.background.y = element_blank(), legend.background = element_blank(), axis.title = element_blank(), panel.spacing.y = unit(0, "mm"), panel.ontop = T)
        )
    
    return(p)       
}


LT_Il1rnKO <- fread("unzip -p ../data/LT_Il1rnKO_DESeq2.txt.zip")
ST_Il1rnKO <- fread("unzip -p ../data/ST_Il1rnKO_DESeq2.txt.zip")

selected <- c('cell activation','cytokine production','leukocyte differentiation','regulation of immune system process','oxidative phosphorylation','autophagy','defense response','immune response','cell surface receptor signaling pathway','JAK-STAT cascade','macroautophagy','antigen processing and presentation','regulation of anatomical structure morphogenesis','positive regulation of signaling','regulation of protein modification process','cellular protein-containing complex assembly','regulation of locomotion','peptide biosynthetic process','secretion','positive regulation of response to stimulus','positive regulation of developmental process','positive regulation of multicellular organismal process','positive regulation of protein metabolic process','import into cell','regulation of intracellular signal transduction')

gseaBP_all <- data.frame()
scale=F; nGenesToShrink=10

objects <- c("LT_Il1rnKO", "ST_Il1rnKO")

for (i in 1:length(objects)) {
    message("\n", "Analyzing ", objects[i])
    
    temp <- get(objects[i])[, c(1:9)]
    temp <- temp[!is.na(temp$pvalue) & !is.na(temp$padj) & !is.na(temp$log2FoldChange), ]
    idx = which(temp$baseMean < 2000 & temp$padj > 0.05)
    idx = c(idx, which(temp$padj > 0.05 & abs(temp$log2FoldChange) > 0.5))
    idx2 = which(!1:nrow(temp) %in% idx)
    temp <- temp[idx2, ]
    temp <- temp[temp$baseMean > 100, ]
    temp <- temp[temp$padj < 1, ]
    temp <- temp[temp$gene_biotype == "protein_coding", ]
    
    temp$dist <- sqrt(temp$log2FoldChange^2 + -log10(temp$padj))
    temp$dist <- ifelse(temp$log2FoldChange < 0, -temp$dist, temp$dist)
    geneList <- setNames(temp$dist, temp$ensembl_gene_id)
    geneList <- geneList[order(-geneList)]
    
    counter = nGenesToShrink
    while(counter > 0) { 
        if(counter == nGenesToShrink) message("Shrinking ...")
        
        x <- geneList[counter:length(geneList)]
        outTest <- outliers::grubbs.test(x)
        if (outTest$p.value < 0.05) {
            value = gsub(" is an outlier|highest value |lowest value ", "", outTest$alternative)
            geneList[counter] <- geneList[counter+1]+0.00001
            message("Shrinked gene ", counter, " in list")
        }
        
        x <- geneList[1:(length(geneList)-counter)]
        outTest <- outliers::grubbs.test(x, opposite = T)
        if (outTest$p.value < 0.05) {
            value = gsub(" is an outlier|highest value |lowest value ", "", outTest$alternative)
            idx = (length(geneList)-counter)
            geneList[idx] <- geneList[idx-1]-0.00001
            message("Shrinked gene ", idx, " in list")
        }
        
    counter = counter - 1
    }

    gseaBP <- gseGO(geneList,keyType =  "ENSEMBL", OrgDb = org.Mm.eg.db, ont = "BP", minGSSize = 10, maxGSSize = 500, nPerm = 5000, pvalueCutoff = 0.05)
    
    # gseaplot for leukocyte differention
    text = paste("NES", round(gseaBP@result["leukocyte differentiation" == gseaBP@result$Description, ]$NES,2), 
    "_q", signif(gseaBP@result["leukocyte differentiation" == gseaBP@result$Description, ]$qvalue,3), "_", sep = "")
    gseaplot2(gseaBP, "GO:0002521", pvalue_table = F, subplots = c(1:2))
    
    p <- ridgeplot2(gseaBP, simplify = F, simplifyCutoff = 0.5, selected = unlist(selected, use.names = F), text.size = 10, n.UP = 100, n.DOWN = 100)
    p <- p+scale_fill_continuous(low = "red", high = "blue", limits = c(0,0.04)) 
    file = paste("../res/", objects[i], ".pdf", sep = "")
    pdf(file, height = 4.25, width = 6)
    print(p)
    dev.off()
    
    gseaBP <- clusterProfiler::setReadable(gseaBP, OrgDb = org.Mm.eg.db)
    gseaBP <- clusterProfiler::simplify(gseaBP)
    assign(paste("gseaBP", objects[i], sep = ""), gseaBP)
    gseaBP@result$type <- objects[i]
    gseaBP_all <- rbind(gseaBP_all, gseaBP@result)
}


#######################################################################################
#######################                   NRAS                  #######################
#######################          Venn HSPCs (Fig. S6I)          #######################
#######################           PCA HSPCs (Fig. S6J)          #######################
#######################################################################################

#read table 
unzip("../data/raw_mat_Nras.csv.zip")
m2   <- read.table("../data/raw_mat_Nras.csv", row.names=1, header=T, sep="\t")
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

unzip("../data/pheno_Nras.csv.zip")
pheno <- read.table("../data/pheno_Nras.csv", header=T, sep="\t")
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
##############                          Nes-GFP                          ##############
##############             Gene-set enrichment plot (Fig. S6L)           ##############
#######################################################################################

mart.mouse = biomaRt::useMart("ENSEMBL_MART_ENSEMBL", dataset='mmusculus_gene_ensembl')
msigdb <- msigdbr::msigdbr(species = "Mus musculus")
msigdb <- msigdb[msigdb$gs_cat == "C2" & msigdb$gs_subcat == "CGP", ]
msigdb <- msigdb[, c("gs_name", "gene_symbol")]
length(unique(msigdb$gs_name))
t <- data.table::fread("unzip -p ../data/Nes_GFP_DGE.txt.zip")
df <- getBM(filters = "ensembl_transcript_id_version", attributes = c("ensembl_transcript_id_version", 
        "transcript_biotype"), values = t$ensembl_transcript_id, mart = mart.mouse)
t <- merge(t, df, by.x = "ensembl_transcript_id",by.y = "ensembl_transcript_id_version")
t <- t[t$transcript_biotype == "protein_coding", ]
t <- t[t$mean_obs > 5, ]
t <- t[grep("^mt-|^Rpl|^Rps", t$external_gene_name, invert = T),]
geneList <- setNames(t$b, t$external_gene_name)
geneList <- geneList[which ( !duplicated(names(geneList)) )]
geneList <- geneList[order(-geneList)]
idx1 = which(t$NRAS_TG_rep1 > 10 & t$NRAS_TG_rep2 > 10)
idx2 = which(t$NRAS_WT_rep1 > 10 & t$NRAS_WT_rep2 > 10)
geneList <- geneList[names(geneList) %in% unlist(t[unique(idx1, idx2), "external_gene_name"])]
gsea <- clusterProfiler::GSEA(geneList, TERM2GENE = msigdb)
enrichplot::gseaplot2(gsea, "HAMAI_APOPTOSIS_VIA_TRAIL_UP", title = "Hamai apoptosis via TRAIL up", subplots = c(1,2))



#######################################################################################
###########               Single Cell data                                 ############
###########                                                                ############
#######################################################################################
pkgs <- c('reshape2','cowplot')
for(i in 1:length(pkgs)) {
  library(pkgs[i], character.only = T)
}
blank_background <- theme(panel.grid.major = element_blank()
                          , panel.grid.minor = element_blank()
                          , panel.background = element_blank())
gradient_colors <- colorRampPalette(c("grey"
                                    ,"orange"
                                    ,"red"))(3)

celltype_plot <- function(group){
  ns <- mdata
  ns[which(ns$SampleName!= group),"Celltype"] <- NA
  ns <- ns[order(ns$Celltype, na.last = F),]
  tsne_temp_gene <- ggplot(data=ns) +
    geom_point(aes(x=tSNE_1, y=tSNE_2
                   , colour=Celltype
    )
    , size = 0.7) +
    labs(x="tSNE 1", y="tSNE 2") +
    ggtitle(group) +
    scale_color_manual(values = celltype_cols, na.value = "white") +
    guides(colour = guide_legend(override.aes = list(size=5))) +
    theme(axis.text=element_text(size=8)
          , axis.title=element_text(size=8)
          , legend.text=element_text(size=7)
          , legend.position="right"
          , legend.title = element_blank()
          , axis.line=element_line(color="gray88")
          , panel.border = element_blank()) +
    blank_background
  return(tsne_temp_gene)
}

gene_plot <- function(gene){
  temp_data <- mdata[, c("SampleName", "tSNE_1", "tSNE_2", gene)]
  temp_data$gene_t <- temp_data[, gene]
  gene_min <- min(temp_data$gene)
  gene_max <- max(temp_data$gene)
  gene_plots <- vector('list')
  for(sample in rev(unique(temp_data$SampleName))){
    ns <- temp_data
    ns[which(ns$SampleName != sample), "gene_t"] <- NA
    ns <- ns[order(ns$gene_t, decreasing=F, na.last = F),]
    gene_plots[[sample]] <- ggplot(data=ns) +
      geom_point(aes(x=tSNE_1, y=tSNE_2
                     , colour=gene_t
      ), size = 1) +
      scale_color_gradientn(colours = gradient_colors
                            , na.value = "grey85"
                            , guide = guide_colourbar(barheight=0.5, label.vjust=-1)
                            , breaks = c(gene_min, gene_max)
                            , labels = c(round(gene_min, digits = 2)
                                         , round(gene_max, digits = 2))
                            , limits = c(gene_min, gene_max)) +
      labs(x="t-SNE 1", y="t-SNE 2") +
      ggtitle(sample) +
      blank_background +
      theme(legend.position="bottom"
            , axis.text=element_text(size=8)
            , axis.title=element_text(size=8)
            , axis.line=element_line(color="grey90")
            , panel.border = element_blank()
            , legend.title=element_blank()
            , legend.text=element_text(size=7))
    
  }
  title <- ggdraw() + draw_label(gene, fontface='bold')
  temp_grid <- plot_grid(plotlist = gene_plots, ncol = 3)
  gene_grid <- plot_grid(title, temp_grid, ncol=1, rel_heights=c(0.1, 1))
  return(gene_grid)
}
#######################################################################################
###########                             Figure 4                           ############
#######################################################################################
mdata <- read.csv("./data/cd11b_mdata.csv")
colnames(mdata) <- plyr::revalue(colnames(mdata), c(celltype_cluster_reval="Celltype"))
celltype_cols <- c("#ebac23", "#008cf9", "#006e00", "#00bbad")
names(celltype_cols) <- unique(mdata$Celltype)

# 4A tsnes
plot_list <- lapply(rev(unique(mdata$SampleName)), celltype_plot)
plot_myeloid_separate <- plot_grid(plotlist=plot_list)
plot_myeloid_separate

# 4A pie charts
counts <- table(mdata$Celltype, mdata$SampleName)
pie_chart_data <- setNames(reshape2::melt(counts)
                           , c("predictedID", "Sample", "Cells") 
) 
ggplot(pie_chart_data[which(pie_chart_data$Sample=="CD11b_WT"),], aes(x="", y=Cells, fill = predictedID)) +
  geom_bar(stat="identity", width=1) +
  scale_fill_manual(values = celltype_cols
                    , na.value = "grey80") +
  ggtitle("WT") +
  coord_polar("y", start=0) +
  theme_void() +
  blank_background
ggplot(pie_chart_data[which(pie_chart_data$Sample=="CD11b_KO"),], aes(x="", y=Cells, fill = predictedID)) +
  geom_bar(stat="identity", width=1) +
  scale_fill_manual(values = celltype_cols
                    , na.value = "grey80") +
  ggtitle("KO") +
  coord_polar("y", start=0) +
  theme_void() +
  blank_background

# 4B 
gene_list <- c("Il1b", "Il1rn")
plot_list_il <- lapply(gene_list, gene_plot)
plot_grid(plotlist = plot_list_il, ncol = 1)

# 4E tsnes
tsnedummy_t <- data.frame(Embeddings(object = lorena_lsk_filtered, reduction = "tsne"))
tsnedummy_u <- data.frame(Embeddings(object = lorena_lsk_filtered, reduction = "umap"))
tsnedummy <- merge(tsnedummy_t,tsnedummy_u, by=0)
tsnedummy <- merge(tsnedummy, lorena_lsk_filtered@meta.data, by.x="Row.names", by.y=0)
gene_expression <-  data.frame(lorena_lsk_filtered@assays$RNA@data)[c("Il1b", "Il1rn"),]
colnames(gene_expression) <- gsub("\\.", "-", colnames(gene_expression))
tsnedummy <- merge(tsnedummy, t(gene_expression), by.x = "Row.names", by.y=0)
rownames(tsnedummy) <- tsnedummy$Row.names
tsnedummy <- tsnedummy[,-1]
rownames(tsnedummy) <- 1:nrow(tsnedummy)
colnames(tsnedummy) <- plyr::revalue(colnames(tsnedummy), c(alejo_population="Sommerkamp_annot"))
write.csv(tsnedummy[,c("tSNE_1", "tSNE_2", "SampleName", "Sommerkamp_annot", "Il1b", "Il1rn")]
        , "/home/arubio/NetVolumes/LAB_AH/LAB/Andrea/projects/lorena_sc_il1KO/Nature_code/data/lsk_mdata_1.csv")




mdata <- read.csv("./data/lsk_mdata.csv")
colnames(mdata) <- plyr::revalue(colnames(mdata), c(Sommerkamp_annot="Celltype"))
celltype_cols <- c("#ebac23", "#b80058", "#008cf9", "#006e00")
names(celltype_cols) <- unique(mdata$Celltype)

plot_list <- lapply(rev(unique(mdata$SampleName)), celltype_plot)
plot_grid(plotlist = plot_list)

# 4E pie charts
counts <- table(mdata$Celltype, mdata$SampleName)
pie_chart_data <- setNames(reshape2::melt(counts)
                           , c("predictedID", "Sample", "Cells") 
) 
ggplot(pie_chart_data[which(pie_chart_data$Sample=="LSK_WT"),], aes(x="", y=Cells, fill = predictedID)) +
  geom_bar(stat="identity", width=1) +
  scale_fill_manual(values = celltype_cols
                    , na.value = "grey80") +
  ggtitle("WT") +
  coord_polar("y", start=0) +
  theme_void() +
  blank_background
ggplot(pie_chart_data[which(pie_chart_data$Sample=="LSK_KO"),], aes(x="", y=Cells, fill = predictedID)) +
  geom_bar(stat="identity", width=1) +
  scale_fill_manual(values = celltype_cols
                    , na.value = "grey80") +
  ggtitle("KO") +
  coord_polar("y", start=0) +
  theme_void() +
  blank_background

# 4G
mdata <- read.csv("./data/cd63_mdata.csv")
colnames(mdata) <- plyr::revalue(colnames(mdata), c(predicted.id="Celltype"))
celltype_cols <- c("#ebac23", "#b80058", "#008cf9", "#006e00", "#00bbad", "#d163e6")
names(celltype_cols) <- unique(mdata$Celltype)

# tsnes
plot_list <- lapply(rev(unique(mdata$SampleName)), celltype_plot)
plot_grid(plotlist = plot_list)

# Pie Charts
counts <- table(mdata$Celltype, mdata$SampleName)
pie_chart_data <- setNames(reshape2::melt(counts)
                           , c("predictedID", "Sample", "Cells") 
) 
ggplot(pie_chart_data[which(pie_chart_data$Sample=="CD63_WT"),], aes(x="", y=Cells, fill = predictedID)) +
  geom_bar(stat="identity", width=1) +
  scale_fill_manual(values = celltype_cols
                    , na.value = "grey80") +
  ggtitle("WT") +
  coord_polar("y", start=0) +
  theme_void() +
  blank_background
ggplot(pie_chart_data[which(pie_chart_data$Sample=="CD63_KO"),], aes(x="", y=Cells, fill = predictedID)) +
  geom_bar(stat="identity", width=1) +
  scale_fill_manual(values = celltype_cols
                    , na.value = "grey80") +
  ggtitle("KO") +
  coord_polar("y", start=0) +
  theme_void() +
  blank_background

# 4H
gene_list <- c("Il1b", "Il1rn")
plot_list_il <- lapply(gene_list, gene_plot)
plot_grid(plotlist = plot_list_il, ncol = 1)

#######################################################################################
###########                           Figure S5                            ############
#######################################################################################

# S5 B
# HSC-MPP5
wt_df <- read.csv("./data/hsc_mpp5_df.csv")
wt_df$Categories <- sapply(strsplit(wt_df$Categories, ","), "[[", 1)
wt_df$Activation.z.score <- wt_df$Activation.z.score*(-1)
wt_df$colour <- ifelse(wt_df$Activation.z.score < 0, "negative","positive")
wt_df$hjust <- ifelse(wt_df$Activation.z.score > 0, 1.3, -0.3)
wt_df <- wt_df[order(wt_df$Activation.z.score),]
wt_df$Diseases.or.Functions.Annotation <- factor(wt_df$Diseases.or.Functions.Annotation, levels = wt_df$Diseases.or.Functions.Annotation)
ggplot(wt_df, aes(x=Diseases.or.Functions.Annotation, y=Activation.z.score, label="", hjust=hjust))+
  geom_bar(stat="identity", position = "identity", aes(fill = colour))+
  scale_fill_manual(values=c(positive="firebrick1", negative="steelblue")) +
  ylab("Z-score") +
  xlab("") +
  coord_flip() +
  blank_background +
  theme(axis.text=element_text(size=8)
        , axis.title=element_text(size=8)
        , legend.text=element_text(size=7)
        , legend.position="right"
        , legend.title = element_blank()
        , axis.line=element_line(color="gray88")
        , panel.border = element_blank())

# MPP1
wt_df <- read.csv("./data/mpp1_df.csv", stringsAsFactors = F)
wt_df$Categories <- sapply(strsplit(wt_df$Categories, ","), "[[", 1)
wt_df$Activation.z.score <- wt_df$Activation.z.score*(-1)
wt_df$colour <- ifelse(wt_df$Activation.z.score < 0, "negative","positive")
wt_df$hjust <- ifelse(wt_df$Activation.z.score > 0, 1.3, -0.3)
wt_df <- wt_df[order(wt_df$Activation.z.score),]
wt_df$Diseases.or.Functions.Annotation <- factor(wt_df$Diseases.or.Functions.Annotation, levels = wt_df$Diseases.or.Functions.Annotation)
ggplot(wt_df, aes(x=Diseases.or.Functions.Annotation, y=Activation.z.score, label="", hjust=hjust))+
  geom_bar(stat="identity", position = "identity", aes(fill = colour)) +
  scale_fill_manual(values=c(positive="firebrick1", negative="steelblue")) +
  ylab("Z-score") +
  xlab("") +
  coord_flip() +
  blank_background +
  theme(axis.text=element_text(size=8)
        , axis.title=element_text(size=8)
        , legend.text=element_text(size=7)
        , legend.position="right"
        , legend.title = element_blank()
        , axis.line=element_line(color="gray88")
        , panel.border = element_blank())

# MPP2_MPP3
wt_df <- read.csv("./data/mpp2_mpp3_df.csv", stringsAsFactors = F)
wt_df$Categories <- sapply(strsplit(wt_df$Categories, ","), "[[", 1)
wt_df$Activation.z.score <- wt_df$Activation.z.score*(-1)
wt_df$colour <- ifelse(wt_df$Activation.z.score < 0, "negative","positive")
wt_df$hjust <- ifelse(wt_df$Activation.z.score > 0, 1.3, -0.3)
wt_df <- wt_df[order(wt_df$Activation.z.score),]
wt_df$Diseases.or.Functions.Annotation <- factor(wt_df$Diseases.or.Functions.Annotation, levels = wt_df$Diseases.or.Functions.Annotation)
ggplot(wt_df, aes(x=Diseases.or.Functions.Annotation, y=Activation.z.score, label="", hjust=hjust))+
  geom_bar(stat="identity", position = "identity", aes(fill = colour)) +
  scale_fill_manual(values=c(positive="firebrick1", negative="steelblue")) +
  ylab("Z-score") +
  xlab("") +
  coord_flip() +
  blank_background +
  theme(axis.text=element_text(size=8)
        , axis.title=element_text(size=8)
        , legend.text=element_text(size=7)
        , legend.position="right"
        , legend.title = element_blank()
        , axis.line=element_line(color="gray88")
        , panel.border = element_blank())

# S5 C
gene_list <- c("Il1b", "Il1rn")
mdata <- read.csv("./data/lsk_mdata.csv")
plot_list_il <- lapply(gene_list, gene_plot)
plot_grid(plotlist = plot_list_il, ncol = 1)
