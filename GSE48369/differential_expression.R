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
GSE48369.raw <- getGEO("GSE48369", destdir = "./save")[[1]]

supplemental.files <- getGEOSuppFiles('GSE48369', baseDir = '.')
GSE48369.affy <- read.affybatch(list.files("./GSE48369/GSE48369_RAW", full.names = TRUE))
sampleNames(GSE48369.affy) %<>% str_replace("_.*", "")

pheno.data <- pData(GSE48369.raw) %>% select(title, geo_accession)
rownames(pheno.data) <- pheno.data$geo_accession
pheno.mess <- str_replace(pheno.data$title, "Sleep dep", "Sleep_dep") %>% str_replace(",", "") %>% str_split_fixed(" ", 7) %>% magrittr::extract(, 1:3)
pheno.data$Fraction <- map2_chr(pheno.mess[,1], pheno.mess[,2], str_c, sep = "_")
pheno.data$Condition <- pheno.mess[,3]

pData(GSE48369.affy) <- pheno.data
GSE48369.norm <- rma(GSE48369.affy)

ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
bm.table <- getBM(attributes = c('affy_mouse430_2', 'ensembl_gene_id', 'entrezgene', 'external_gene_name', 'description'), filters = 'affy_mouse430_2', values = featureNames(GSE48369.norm), mart = ensembl)
bm.table$description %<>% str_replace_all(" \\[.*\\]$", "")
colnames(bm.table) <- c("Affymetrix.ID", "Ensembl.ID", "EntrezGene.ID", "Symbol", "Definition")
bm.filter <- filter(bm.table, !is.na(Symbol))

GSE48369.annot <- GSE48369.norm[bm.filter$Affymetrix.ID,]
GSE48369.expr <- exprs(GSE48369.annot)
rownames(GSE48369.expr) <- 1:nrow(GSE48369.annot)
GSE48369.collapse.expr <- collapseRows(GSE48369.expr, bm.filter$Symbol, rownames(GSE48369.expr)) #collapseRows by symbol
GSE48369.collapse <- ExpressionSet(assayData = GSE48369.collapse.expr$datETcollapsed, phenoData = phenoData(GSE48369.annot))

GSE48369.bound <- GSE48369.collapse[,GSE48369.collapse$Fraction == "Bound_fraction"]
GSE48369.unbound <- GSE48369.collapse[,GSE48369.collapse$Fraction == "Unbound_fraction"]

GSE48369.bound.model <- model.matrix( ~ 0 + Condition, pData(GSE48369.bound))
colnames(GSE48369.bound.model) %<>% str_replace("Condition", "") %>% make.names
contrasts.bound <- makeContrasts("Sleep_dep - Sleep", "Sleep_dep - Waking", "Waking - Sleep", levels = GSE48369.bound.model)
GSE48369.bound.lm <- lmFit(GSE48369.bound, GSE48369.bound.model) %>% eBayes

toptable.bound.sds <- topTable(GSE48369.bound.lm, coef = 1, n = Inf) %>% select(AveExpr, logFC, t:B)
toptable.bound.sds$Symbol <- rownames(toptable.bound.sds)
colnames(toptable.bound.sds)[2:6] %<>% str_c(".sds")

toptable.bound.sdw <- topTable(GSE48369.bound.lm, coef = 2, n = Inf) %>% select(AveExpr, logFC, t:B)
toptable.bound.sdw$Symbol <- rownames(toptable.bound.sdw)
colnames(toptable.bound.sdw)[2:6] %<>% str_c(".sdw")

toptable.bound.ws <- topTable(GSE48369.bound.lm, coef = 3, n = Inf) %>% select(AveExpr, logFC, t:B)
toptable.bound.ws$Symbol <- rownames(toptable.bound.ws)
colnames(toptable.bound.ws)[2:6] %<>% str_c(".ws")

toptable.bound.all <- left_join(toptable.bound.sds, toptable.bound.sdw) %>% left_join(toptable.bound.ws)
SaveRDSgz(toptable.bound.all, "./save/toptable.bound.all.rda")

GSE48369.unbound.model <- model.matrix( ~ 0 + Condition, pData(GSE48369.unbound))
colnames(GSE48369.unbound.model) %<>% str_replace("Condition", "") %>% make.names
contrasts.unbound <- makeContrasts("Sleep_dep - Sleep", "Sleep_dep - Waking", "Waking - Sleep", levels = GSE48369.unbound.model)
GSE48369.unbound.lm <- lmFit(GSE48369.unbound, GSE48369.unbound.model) %>% eBayes

toptable.unbound.sds <- topTable(GSE48369.unbound.lm, coef = 1, n = Inf) %>% select(AveExpr, logFC, t:B)
toptable.unbound.sds$Symbol <- rownames(toptable.unbound.sds)
colnames(toptable.unbound.sds)[2:6] %<>% str_c(".sds")

toptable.unbound.sdw <- topTable(GSE48369.unbound.lm, coef = 2, n = Inf) %>% select(AveExpr, logFC, t:B)
toptable.unbound.sdw$Symbol <- rownames(toptable.unbound.sdw)
colnames(toptable.unbound.sdw)[2:6] %<>% str_c(".sdw")

toptable.unbound.ws <- topTable(GSE48369.unbound.lm, coef = 3, n = Inf) %>% select(AveExpr, logFC, t:B)
toptable.unbound.ws$Symbol <- rownames(toptable.unbound.ws)
colnames(toptable.unbound.ws)[2:6] %<>% str_c(".ws")

toptable.unbound.all <- left_join(toptable.unbound.sds, toptable.unbound.sdw) %>% left_join(toptable.unbound.ws)
SaveRDSgz(toptable.unbound.all, "./save/toptable.unbound.all.rda")
