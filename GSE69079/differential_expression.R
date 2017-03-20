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
GSE69079.raw <- getGEO("GSE69079", destdir = "./save")[[1]]

supplemental.files <- getGEOSuppFiles('GSE69079', baseDir = '.')
GSE69079.affy <- read.affybatch(list.files("./GSE69079/GSE69079_RAW", full.names = TRUE))
sampleNames(GSE69079.affy) %<>% str_replace("_.*", "")

pheno.data <- pData(GSE69079.raw) %>% select(title, geo_accession)
rownames(pheno.data) <- pheno.data$geo_accession
pheno.mess <- str_replace(pheno.data$title, "Sleep dep", "Sleep_dep") %>% str_replace(",", "") %>% str_split_fixed(" ", 7) %>% magrittr::extract(, 1:3)
pheno.data$Fraction <- map2_chr(pheno.mess[,1], pheno.mess[,2], str_c, sep = "_")
pheno.data$Condition <- pheno.mess[,3]

pData(GSE69079.affy) <- pheno.data
GSE69079.norm <- rma(GSE69079.affy)

ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
bm.table <- getBM(attributes = c('affy_mouse430_2', 'ensembl_gene_id', 'entrezgene', 'external_gene_name', 'description'), filters = 'affy_mouse430_2', values = featureNames(GSE69079.norm), mart = ensembl)
bm.table$description %<>% str_replace_all(" \\[.*\\]$", "")
colnames(bm.table) <- c("Affymetrix.ID", "Ensembl.ID", "EntrezGene.ID", "Symbol", "Definition")
bm.filter <- filter(bm.table, !is.na(Symbol))

GSE69079.annot <- GSE69079.norm[bm.filter$Affymetrix.ID,]
GSE69079.expr <- exprs(GSE69079.annot)
rownames(GSE69079.expr) <- 1:nrow(GSE69079.annot)
GSE69079.collapse.expr <- collapseRows(GSE69079.expr, bm.filter$Symbol, rownames(GSE69079.expr)) #collapseRows by symbol
GSE69079.collapse <- ExpressionSet(assayData = GSE69079.collapse.expr$datETcollapsed, phenoData = phenoData(GSE69079.annot))

GSE69079.bound <- GSE69079.collapse[,grepl("Bound", GSE69079.collapse$Fraction)]
GSE69079.bound$title %<>% str_replace(" animal.*$", "") %>% str_replace(".*fraction ", "") %>% str_replace(" ", "_")

GSE69079.bound.model <- model.matrix( ~ 0 + title, pData(GSE69079.bound))
colnames(GSE69079.bound.model) %<>% str_replace("title", "") %>% make.names
contrasts.bound <- makeContrasts("Sleep_dep - Sleep", "Sleep_dep - Waking", "Waking - Sleep", levels = GSE69079.bound.model)
GSE69079.bound.lm <- lmFit(GSE69079.bound, GSE69079.bound.model) %>% eBayes

toptable.bound.sds <- topTable(GSE69079.bound.lm, coef = 1, n = Inf) %>% select(AveExpr, logFC, t:B)
toptable.bound.sds$Symbol <- rownames(toptable.bound.sds)
colnames(toptable.bound.sds)[2:6] %<>% str_c(".sds")

toptable.bound.sdw <- topTable(GSE69079.bound.lm, coef = 2, n = Inf) %>% select(AveExpr, logFC, t:B)
toptable.bound.sdw$Symbol <- rownames(toptable.bound.sdw)
colnames(toptable.bound.sdw)[2:6] %<>% str_c(".sdw")

toptable.bound.ws <- topTable(GSE69079.bound.lm, coef = 3, n = Inf) %>% select(AveExpr, logFC, t:B)
toptable.bound.ws$Symbol <- rownames(toptable.bound.ws)
colnames(toptable.bound.ws)[2:6] %<>% str_c(".ws")

toptable.bound.all <- left_join(toptable.bound.sds, toptable.bound.sdw) %>% left_join(toptable.bound.ws)
SaveRDSgz(toptable.bound.all, "./save/toptable.bound.all.rda")
