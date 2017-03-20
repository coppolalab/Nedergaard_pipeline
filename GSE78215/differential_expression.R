library(oligo)
library(Biobase)
library(annotate)
library(WGCNA)
library(biomaRt)
library(broom)

library(GEOquery)
library(limma)

library(stringr)
library(readr)
library(openxlsx)

library(rlist)
library(reshape2)
library(dplyr)
library(purrr)
library(magrittr)
library(tidyr)

library(ggplot2)
library(Cairo)

source("../../FRDA project/common_functions.R")
GSE78215.raw <- getGEO("GSE78215", destdir = "./save")[[1]]

supplemental.files <- getGEOSuppFiles('GSE78215', baseDir = '.')
GSE78215.affy <- read.celfiles(list.files("./GSE78215/GSE78215_RAW", full.names = TRUE))
sampleNames(GSE78215.affy) %<>% str_replace("_.*", "")

pheno.data <- pData(GSE78215.raw) %>% select(geo_accession, characteristics_ch1.3)
rownames(pheno.data) <- pheno.data$geo_accession
colnames(pheno.data) <- c("GEO.Accession", "Condition")

pheno.data$Condition %<>% str_replace("sample group: ", "") %>% str_replace(", ", "_") %>% str_replace(" ", "_")
pData(GSE78215.affy) <- pheno.data
GSE78215.norm <- rma(GSE78215.affy)

ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
bm.table <- getBM(attributes = c('affy_mogene_2_1_st_v1', 'ensembl_gene_id', 'entrezgene', 'external_gene_name', 'description'), filters = 'affy_mogene_2_1_st_v1', values = featureNames(GSE78215.norm), mart = ensembl)
bm.table$description %<>% str_replace_all(" \\[.*\\]$", "")
colnames(bm.table) <- c("Affymetrix.ID", "Ensembl.ID", "EntrezGene.ID", "Symbol", "Definition")
bm.filter <- filter(bm.table, !is.na(Symbol)) %>% filter(!duplicated(Ensembl.ID))

GSE78215.annot <- GSE78215.norm[as.character(bm.filter$Affymetrix.ID),]
GSE78215.expr <- exprs(GSE78215.annot)
rownames(GSE78215.expr) <- 1:nrow(GSE78215.annot)
GSE78215.collapse.expr <- collapseRows(GSE78215.expr, bm.filter$Symbol, rownames(GSE78215.expr)) #collapseRows by symbol
GSE78215.collapse <- ExpressionSet(assayData = GSE78215.collapse.expr$datETcollapsed, phenoData = phenoData(GSE78215.annot))

GSE78215.model <- model.matrix( ~ 0 + Condition, pData(GSE78215.collapse))
colnames(GSE78215.model) %<>% str_replace("Condition", "") %>% make.names
contrasts.bound <- makeContrasts("X6hrSD_and.1Hr.Recovery_ZT7 - Circadian_Control_ZT7", "X6HrSD_and.2Hr.Recovery_ZT8 - Circadian_Control_ZT8", "X5HrSD_and.3Hr.Recovery_ZT8 - Circadian_Control_ZT8", "X5HrSD_and.6Hr.Recovery_ZT11 - Circadian_Control_ZT11", levels = GSE78215.model)
GSE78215.lm <- lmFit(GSE78215.collapse, GSE78215.model) %>% eBayes

toptable.bound.zt7 <- topTable(GSE78215.lm, coef = 1, n = Inf) %>% select(AveExpr, logFC, t:B)
toptable.bound.zt7$Symbol <- rownames(toptable.bound.zt7)
colnames(toptable.bound.zt7)[2:6] %<>% str_c(".zt7")

toptable.bound.zt8_1 <- topTable(GSE78215.lm, coef = 2, n = Inf) %>% select(AveExpr, logFC, t:B)
toptable.bound.zt8_1$Symbol <- rownames(toptable.bound.zt8_1)
colnames(toptable.bound.zt8_1)[2:6] %<>% str_c(".zt8_1")

toptable.bound.zt8_2 <- topTable(GSE78215.lm, coef = 3, n = Inf) %>% select(AveExpr, logFC, t:B)
toptable.bound.zt8_2$Symbol <- rownames(toptable.bound.zt8_2)
colnames(toptable.bound.zt8_2)[2:6] %<>% str_c(".zt8_2")

toptable.bound.zt11 <- topTable(GSE78215.lm, coef = 4, n = Inf) %>% select(AveExpr, logFC, t:B)
toptable.bound.zt11$Symbol <- rownames(toptable.bound.zt11)
colnames(toptable.bound.zt11)[2:6] %<>% str_c(".zt11")

toptable.bound.all <- left_join(toptable.bound.zt7, toptable.bound.zt8_1) %>% left_join(toptable.bound.zt8_2) %>% left_join(toptable.bound.zt11)
SaveRDSgz(toptable.bound.all, "./save/toptable.bound.all.rda")
