library(Biobase)
library(annotate)
library(WGCNA)
library(biomaRt)
library(broom)

library(affy)
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
GSE9444.raw <- getGEO("GSE9444", destdir = "./save")[[1]]

supplemental.files <- getGEOSuppFiles('GSE9444', baseDir = '.')
GSE9444.affy <- read.affybatch(list.files("./GSE9444/GSE9444_RAW", full.names = TRUE))
sampleNames(GSE9444.affy) %<>% str_replace("\\.CEL\\.gz", "")

pheno.data <- pData(GSE9444.raw) %>% select(title, geo_accession, source_name_ch1, characteristics_ch1)
GSE9441.pheno <- filter(pheno.data, grepl("^L_|^B_", title)) %>% droplevels
GSE9442.pheno <- filter(pheno.data, grepl("^AK_|^B6_|^D2_", title)) %>% droplevels
GSE9443.pheno <- filter(pheno.data, grepl("pull\\-down|surnSD|RNAtot", title)) %>% droplevels

#C57/BL6, AKR/J, DBA/2J (Tissue x Condition)
GSE9441.source_name_ch1 <- str_split_fixed(GSE9441.pheno$source_name_ch1, ",", 4) %>% data.frame
colnames(GSE9441.source_name_ch1) <- c("Strain", "Condition", "Timepoint", "Tissue")
GSE9441.source_name_ch1$Condition %<>% str_replace("^ ", "")
GSE9441.source_name_ch1$Timepoint %<>% str_split_fixed(" ", 5) %>% magrittr::extract(, 5)
GSE9441.source_name_ch1$Tissue %<>% str_replace("^.*, ", "")
GSE9441.pheno.clean <- mutate(GSE9441.source_name_ch1, GEO.Accession = GSE9441.pheno$geo_accession)
rownames(GSE9441.pheno.clean) <- GSE9441.pheno.clean$GEO.Accession

GSE9441.ESet <- GSE9444.affy[,GSE9441.pheno.clean$GEO.Accession]
pData(GSE9441.ESet) <- GSE9441.pheno.clean

#C57/BL6, AKR/J, DBA/2J (ZT0, ZT6, ZT12, ZT18 x Condition)
GSE9442.source_name_ch1 <- str_split_fixed(GSE9442.pheno$source_name_ch1, ",", 4) %>% data.frame
colnames(GSE9442.source_name_ch1) <- c("Strain", "Condition", "Timepoint", "Tissue")
GSE9442.source_name_ch1$Timepoint %<>% str_split_fixed(" ", 5) %>% magrittr::extract(, 5)
GSE9442.pheno.clean <- select(GSE9442.source_name_ch1, -Tissue) %>% mutate(GEO.Accession = GSE9442.pheno$geo_accession)
rownames(GSE9442.pheno.clean) <- GSE9442.pheno.clean$GEO.Accession

GSE9442.ESet <- GSE9444.affy[,as.character(GSE9442.pheno.clean$GEO.Accession)]
pData(GSE9442.ESet) <- GSE9442.pheno.clean

#Homer1a expression cells, Whole Brain x Condition
GSE9443.source_name_ch1 <- str_split_fixed(GSE9443.pheno$source_name_ch1, ",", 2) %>% data.frame
colnames(GSE9443.source_name_ch1) <- c("Tissue", "Condition")
GSE9443.source_name_ch1$Tissue %<>% str_replace("6 hrs", "6hrs")
GSE9443.pheno.clean <- mutate(GSE9443.source_name_ch1, GEO.Accession = GSE9443.pheno$geo_accession)
rownames(GSE9443.pheno.clean) <- GSE9443.pheno.clean$GEO.Accession

GSE9443.ESet <- GSE9444.affy[,as.character(GSE9443.pheno.clean$GEO.Accession)]
pData(GSE9443.ESet) <- GSE9443.pheno.clean

GSE9441.norm <- rma(GSE9441.ESet)
GSE9442.norm <- rma(GSE9442.ESet)
GSE9443.norm <- rma(GSE9443.ESet)

ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
GSE9441.bm.table <- getBM(attributes = c('affy_mouse430_2', 'ensembl_gene_id', 'entrezgene', 'external_gene_name', 'description'), filters = 'affy_mouse430_2', values = featureNames(GSE9441.norm), mart = ensembl)
GSE9441.bm.table$description %<>% str_replace_all(" \\[.*\\]$", "")
colnames(GSE9441.bm.table) <- c("Affymetrix.ID", "Ensembl.ID", "EntrezGene.ID", "Symbol", "Definition")
GSE9441.bm.filter <- filter(GSE9441.bm.table, !is.na(Symbol))

GSE9441.annot <- GSE9441.norm[GSE9441.bm.filter$Affymetrix.ID,]
GSE9441.expr <- exprs(GSE9441.annot)
rownames(GSE9441.expr) <- 1:nrow(GSE9441.annot)
GSE9441.collapse.expr <- collapseRows(GSE9441.expr, GSE9441.bm.filter$Symbol, rownames(GSE9441.expr)) #collapseRows by symbol
GSE9441.collapse <- ExpressionSet(assayData = GSE9441.collapse.expr$datETcollapsed, phenoData = phenoData(GSE9441.annot))
GSE9441.collapse$Combined <- str_c(GSE9441.collapse$Strain, GSE9441.collapse$Condition, sep = " ") %>% str_replace_all(" ", "_")

GSE9442.bm.table <- getBM(attributes = c('affy_mouse430_2', 'ensembl_gene_id', 'entrezgene', 'external_gene_name', 'description'), filters = 'affy_mouse430_2', values = featureNames(GSE9442.norm), mart = ensembl)
GSE9442.bm.table$description %<>% str_replace_all(" \\[.*\\]$", "")
colnames(GSE9442.bm.table) <- c("Affymetrix.ID", "Ensembl.ID", "EntrezGene.ID", "Symbol", "Definition")
GSE9442.bm.filter <- filter(GSE9442.bm.table, !is.na(Symbol))

GSE9442.annot <- GSE9442.norm[GSE9442.bm.filter$Affymetrix.ID,]
GSE9442.expr <- exprs(GSE9442.annot)
rownames(GSE9442.expr) <- 1:nrow(GSE9442.annot)
GSE9442.collapse.expr <- collapseRows(GSE9442.expr, GSE9442.bm.filter$Symbol, rownames(GSE9442.expr)) #collapseRows by symbol
GSE9442.collapse <- ExpressionSet(assayData = GSE9442.collapse.expr$datETcollapsed, phenoData = phenoData(GSE9442.annot))
GSE9442.collapse$Condition %<>% str_replace(" ", "")
GSE9442.collapse$Combined <- str_c(GSE9442.collapse$Condition, GSE9442.collapse$Timepoint, sep = ".") %>% make.names 

GSE9443.bm.table <- getBM(attributes = c('affy_mouse430_2', 'ensembl_gene_id', 'entrezgene', 'external_gene_name', 'description'), filters = 'affy_mouse430_2', values = featureNames(GSE9443.norm), mart = ensembl)
GSE9443.bm.table$description %<>% str_replace_all(" \\[.*\\]$", "")
colnames(GSE9443.bm.table) <- c("Affymetrix.ID", "Ensembl.ID", "EntrezGene.ID", "Symbol", "Definition")
GSE9443.bm.filter <- filter(GSE9443.bm.table, !is.na(Symbol))

GSE9443.annot <- GSE9443.norm[GSE9443.bm.filter$Affymetrix.ID,]
GSE9443.expr <- exprs(GSE9443.annot)
rownames(GSE9443.expr) <- 1:nrow(GSE9443.annot)
GSE9443.collapse.expr <- collapseRows(GSE9443.expr, GSE9443.bm.filter$Symbol, rownames(GSE9443.expr)) #collapseRows by symbol
GSE9443.collapse <- ExpressionSet(assayData = GSE9443.collapse.expr$datETcollapsed, phenoData = phenoData(GSE9443.annot))

#GSE9441
GSE9441.brain <- GSE9441.collapse[,GSE9441.collapse$Tissue == "Brain"]
GSE9441.liver <- GSE9441.collapse[,GSE9441.collapse$Tissue == "Liver"]

GSE9441.model.brain <- model.matrix( ~ 0 + Combined, pData(GSE9441.brain))
colnames(GSE9441.model.brain) %<>% str_replace("Combined", "") %>% make.names
GSE9441.contrasts.brain <- makeContrasts("AKR.J_sleep_deprived - AKR.J_control", "C57BL.6J_sleep_deprived - C57BL.6J_control", "DBA.2J_sleep_deprived - DBA.2J_control", levels = GSE9441.model.brain)
GSE9441.lm.brain <- lmFit(GSE9441.brain, GSE9441.model.brain) %>% contrasts.fit(GSE9441.contrasts.brain) %>% eBayes

toptable.brain.akr <- topTable(GSE9441.lm.brain, coef = 1, n = Inf) %>% select(AveExpr, logFC, t:B)
colnames(toptable.brain.akr)[2:ncol(toptable.brain.akr)] %<>% str_c('_akr')
toptable.brain.akr$Symbol <- rownames(toptable.brain.akr)

toptable.brain.c57 <- topTable(GSE9441.lm.brain, coef = 2, n = Inf) %>% select(AveExpr, logFC, t:B)
colnames(toptable.brain.c57)[2:ncol(toptable.brain.c57)] %<>% str_c('_c57')
toptable.brain.c57$Symbol <- rownames(toptable.brain.c57)

toptable.brain.dba <- topTable(GSE9441.lm.brain, coef = 3, n = Inf) %>% select(AveExpr, logFC, t:B)
colnames(toptable.brain.dba)[2:ncol(toptable.brain.dba)] %<>% str_c('_dba')
toptable.brain.dba$Symbol <- rownames(toptable.brain.dba)

toptable.brain.GSE9441 <- left_join(toptable.brain.akr, toptable.brain.c57) %>% 
    left_join(toptable.brain.dba) %>% 
    select(Symbol, AveExpr, dplyr::contains("logFC"), dplyr::contains("P.Value"), dplyr::contains("adj.P.Val")) 
SaveRDSgz(toptable.brain.GSE9441, "./save/toptable.brain.GSE9441.rda")

GSE9441.model.liver <- model.matrix( ~ 0 + Combined, pData(GSE9441.liver))
colnames(GSE9441.model.liver) %<>% str_replace("Combined", "") %>% make.names
GSE9441.contrasts.liver <- makeContrasts("AKR.J.sleep.deprived.ZT.0 - AKR.J.control.ZT.0", levels = GSE9441.model.liver)
GSE9441.lm.liver <- lmFit(GSE9441.liver, GSE9441.model.liver) %>% contrasts.fit(GSE9441.contrasts.liver) %>% eBayes

toptable.liver.akr <- topTable(GSE9441.lm.liver, coef = 1, n = Inf) %>% select(AveExpr, logFC, t:B)
colnames(toptable.liver.akr)[2:ncol(toptable.liver.akr)] %<>% str_c('_akr')
toptable.liver.akr$Symbol <- rownames(toptable.liver.akr)

toptable.liver.c57 <- topTable(GSE9441.lm.liver, coef = 2, n = Inf) %>% select(AveExpr, logFC, t:B)
colnames(toptable.liver.c57)[2:ncol(toptable.liver.c57)] %<>% str_c('_c57')
toptable.liver.c57$Symbol <- rownames(toptable.liver.c57)

toptable.liver.dba <- topTable(GSE9441.lm.liver, coef = 3, n = Inf) %>% select(AveExpr, logFC, t:B)
colnames(toptable.liver.dba)[2:ncol(toptable.liver.dba)] %<>% str_c('_dba')
toptable.liver.dba$Symbol <- rownames(toptable.liver.dba)

toptable.liver.GSE9441 <- left_join(toptable.liver.akr, toptable.liver.c57) %>% 
    left_join(toptable.liver.dba) %>% 
    select(Symbol, AveExpr, dplyr::contains("logFC"), dplyr::contains("P.Value"), dplyr::contains("adj.P.Val")) 
SaveRDSgz(toptable.liver.GSE9441, "./save/toptable.liver.GSE9441.rda")

#GSE9442
GSE9442.akr <- GSE9442.collapse[,grepl("AKR", GSE9442.collapse$Strain)]
GSE9442.c57 <- GSE9442.collapse[,grepl("C57", GSE9442.collapse$Strain)]
GSE9442.dba <- GSE9442.collapse[,grepl("DBA", GSE9442.collapse$Strain)]

GSE9442.model.akr <- model.matrix( ~ 0 + Combined, pData(GSE9442.akr))
colnames(GSE9442.model.akr) %<>% str_replace("Combined", "") %>% make.names
GSE9442.contrasts.akr <- makeContrasts("sleep.deprived.ZT.0 - control.ZT.0", "sleep.deprived.ZT.6 - control.ZT.6", "sleep.deprived.ZT.12 - control.ZT.12", "sleep.deprived.ZT.18 - control.ZT.18", levels = GSE9442.model.akr)
GSE9442.lm.akr <- lmFit(GSE9442.akr, GSE9442.model.akr) %>% contrasts.fit(GSE9442.contrasts.akr) %>% eBayes

toptable.akr.zt0 <- topTable(GSE9442.lm.akr, coef = 1, n = Inf) %>% select(AveExpr, logFC, t:B)
colnames(toptable.akr.zt0)[2:ncol(toptable.akr.zt0)] %<>% str_c('.zt0')
toptable.akr.zt0$Symbol <- rownames(toptable.akr.zt0)

toptable.akr.zt6 <- topTable(GSE9442.lm.akr, coef = 2, n = Inf) %>% select(AveExpr, logFC, t:B)
colnames(toptable.akr.zt6)[2:ncol(toptable.akr.zt6)] %<>% str_c('.zt6')
toptable.akr.zt6$Symbol <- rownames(toptable.akr.zt6)

toptable.akr.zt12 <- topTable(GSE9442.lm.akr, coef = 3, n = Inf) %>% select(AveExpr, logFC, t:B)
colnames(toptable.akr.zt12)[2:ncol(toptable.akr.zt12)] %<>% str_c('.zt12')
toptable.akr.zt12$Symbol <- rownames(toptable.akr.zt12)

toptable.akr.zt18 <- topTable(GSE9442.lm.akr, coef = 4, n = Inf) %>% select(AveExpr, logFC, t:B)
colnames(toptable.akr.zt18)[2:ncol(toptable.akr.zt18)] %<>% str_c('.zt18')
toptable.akr.zt18$Symbol <- rownames(toptable.akr.zt18)

toptable.akr.all <- left_join(toptable.akr.zt0, toptable.akr.zt6) %>% 
    left_join(toptable.akr.zt12) %>% 
    left_join(toptable.akr.zt18) %>% 
    select(Symbol, AveExpr, dplyr::contains("logFC"), dplyr::contains("P.Value"), dplyr::contains("adj.P.Val")) 
SaveRDSgz(toptable.akr.all, "./save/toptable.akr.all.rda")

GSE9442.model.c57 <- model.matrix( ~ 0 + Combined, pData(GSE9442.c57))
colnames(GSE9442.model.c57) %<>% str_replace("Combined", "") %>% make.names
GSE9442.contrasts.c57 <- makeContrasts("sleep.deprived.ZT.0 - control.ZT.0", "sleep.deprived.ZT.6 - control.ZT.6", "sleep.deprived.ZT.12 - control.ZT.12", "sleep.deprived.ZT.18 - control.ZT.18", levels = GSE9442.model.c57)
GSE9442.lm.c57 <- lmFit(GSE9442.c57, GSE9442.model.c57) %>% contrasts.fit(GSE9442.contrasts.c57) %>% eBayes

toptable.c57.zt0 <- topTable(GSE9442.lm.c57, coef = 1, n = Inf) %>% select(AveExpr, logFC, t:B)
colnames(toptable.c57.zt0)[2:ncol(toptable.c57.zt0)] %<>% str_c('.zt0')
toptable.c57.zt0$Symbol <- rownames(toptable.c57.zt0)

toptable.c57.zt6 <- topTable(GSE9442.lm.c57, coef = 2, n = Inf) %>% select(AveExpr, logFC, t:B)
colnames(toptable.c57.zt6)[2:ncol(toptable.c57.zt6)] %<>% str_c('.zt6')
toptable.c57.zt6$Symbol <- rownames(toptable.c57.zt6)

toptable.c57.zt12 <- topTable(GSE9442.lm.c57, coef = 3, n = Inf) %>% select(AveExpr, logFC, t:B)
colnames(toptable.c57.zt12)[2:ncol(toptable.c57.zt12)] %<>% str_c('.zt12')
toptable.c57.zt12$Symbol <- rownames(toptable.c57.zt12)

toptable.c57.zt18 <- topTable(GSE9442.lm.c57, coef = 4, n = Inf) %>% select(AveExpr, logFC, t:B)
colnames(toptable.c57.zt18)[2:ncol(toptable.c57.zt18)] %<>% str_c('.zt18')
toptable.c57.zt18$Symbol <- rownames(toptable.c57.zt18)

toptable.c57.all <- left_join(toptable.c57.zt0, toptable.c57.zt6) %>% 
    left_join(toptable.c57.zt12) %>% 
    left_join(toptable.c57.zt18) %>% 
    select(Symbol, AveExpr, dplyr::contains("logFC"), dplyr::contains("P.Value"), dplyr::contains("adj.P.Val")) 
SaveRDSgz(toptable.c57.all, "./save/toptable.c57.all.rda")

GSE9442.model.dba <- model.matrix( ~ 0 + Combined, pData(GSE9442.dba))
colnames(GSE9442.model.dba) %<>% str_replace("Combined", "") %>% make.names
GSE9442.contrasts.dba <- makeContrasts("sleep.deprived.ZT.0 - control.ZT.0", "sleep.deprived.ZT.6 - control.ZT.6", "sleep.deprived.ZT.12 - control.ZT.12", "sleep.deprived.ZT.18 - control.ZT.18", levels = GSE9442.model.dba)
GSE9442.lm.dba <- lmFit(GSE9442.dba, GSE9442.model.dba) %>% contrasts.fit(GSE9442.contrasts.dba) %>% eBayes

toptable.dba.zt0 <- topTable(GSE9442.lm.dba, coef = 1, n = Inf) %>% select(AveExpr, logFC, t:B)
colnames(toptable.dba.zt0)[2:ncol(toptable.dba.zt0)] %<>% str_c('.zt0')
toptable.dba.zt0$Symbol <- rownames(toptable.dba.zt0)

toptable.dba.zt6 <- topTable(GSE9442.lm.dba, coef = 2, n = Inf) %>% select(AveExpr, logFC, t:B)
colnames(toptable.dba.zt6)[2:ncol(toptable.dba.zt6)] %<>% str_c('.zt6')
toptable.dba.zt6$Symbol <- rownames(toptable.dba.zt6)

toptable.dba.zt12 <- topTable(GSE9442.lm.dba, coef = 3, n = Inf) %>% select(AveExpr, logFC, t:B)
colnames(toptable.dba.zt12)[2:ncol(toptable.dba.zt12)] %<>% str_c('.zt12')
toptable.dba.zt12$Symbol <- rownames(toptable.dba.zt12)

toptable.dba.zt18 <- topTable(GSE9442.lm.dba, coef = 4, n = Inf) %>% select(AveExpr, logFC, t:B)
colnames(toptable.dba.zt18)[2:ncol(toptable.dba.zt18)] %<>% str_c('.zt18')
toptable.dba.zt18$Symbol <- rownames(toptable.dba.zt18)

toptable.dba.all <- left_join(toptable.dba.zt0, toptable.dba.zt6) %>% 
    left_join(toptable.dba.zt12) %>% 
    left_join(toptable.dba.zt18) %>% 
    select(Symbol, AveExpr, dplyr::contains("logFC"), dplyr::contains("P.Value"), dplyr::contains("adj.P.Val")) 
SaveRDSgz(toptable.dba.all, "./save/toptable.dba.all.rda")

#GSE9443
GSE9443.homer1 <- GSE9443.collapse[,grepl("Homer1a expressing", GSE9443.collapse$Tissue)]
GSE9443.homer1$Condition %<>% str_replace("6hrs ", "") %>% str_replace("^ ", "") %>% make.names
GSE9443.brain <- GSE9443.collapse[,grepl("Whole brain", GSE9443.collapse$Tissue)]
GSE9443.brain$Condition %<>% str_replace("6 hrs ", "") %>% str_replace("^ ", "") %>% make.names

GSE9443.model.brain <- model.matrix( ~ Condition, pData(GSE9443.brain))
colnames(GSE9443.model.brain) %<>% str_replace("Condition", "") %>% make.names
GSE9443.lm.brain <- lmFit(GSE9443.brain, GSE9443.model.brain) %>% eBayes

toptable.GSE9443.brain <- topTable(GSE9443.lm.brain, coef = 2, n = Inf) %>% select(AveExpr, logFC, t:B)
colnames(toptable.GSE9443.brain)[2:ncol(toptable.GSE9443.brain)] %<>% str_c('.brain')
toptable.GSE9443.brain$Symbol <- rownames(toptable.GSE9443.brain)
SaveRDSgz(toptable.GSE9443.brain, "./save/toptable.GSE9443.brain.rda")

GSE9443.model.homer1 <- model.matrix( ~ Condition, pData(GSE9443.homer1))
colnames(GSE9443.model.homer1) %<>% str_replace("Condition", "") %>% make.names
GSE9443.lm.homer1 <- lmFit(GSE9443.homer1, GSE9443.model.homer1) %>% eBayes

toptable.GSE9443.homer1 <- topTable(GSE9443.lm.homer1, coef = 2, n = Inf) %>% select(AveExpr, logFC, t:B)
colnames(toptable.GSE9443.homer1)[2:ncol(toptable.GSE9443.homer1)] %<>% str_c('.homer1')
toptable.GSE9443.homer1$Symbol <- rownames(toptable.GSE9443.homer1)
SaveRDSgz(toptable.GSE9443.homer1, "./save/toptable.GSE9443.homer1.rda")
