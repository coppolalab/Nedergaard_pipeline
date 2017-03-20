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
GSE9471.raw <- getGEO("GSE9471", destdir = "./save")[[1]]

supplemental.files <- getGEOSuppFiles('GSE9471', baseDir = '.')
GSE9471.affy <- read.affybatch(list.files("./GSE9471/GSE9471_RAW", full.names = TRUE))
sampleNames(GSE9471.affy) %<>% str_replace_all("\\..*$", "")

pheno.data <- pData(GSE9471.raw) %>% select(title, geo_accession) 
rownames(pheno.data) <- as.character(pheno.data$geo_accession)
pheno.data$Time.Point <- str_replace(pheno.data$title, "bc-", "") %>% str_replace("-.*$", "")
pData(GSE9471.affy) <- pheno.data

GSE9471.norm <- rma(GSE9471.affy)

ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
GSE9471.bm.table <- getBM(attributes = c('affy_mouse430_2', 'ensembl_gene_id', 'entrezgene', 'external_gene_name', 'description'), filters = 'affy_mouse430_2', values = featureNames(GSE9471.norm), mart = ensembl)
GSE9471.bm.table$description %<>% str_replace_all(" \\[.*\\]$", "")
colnames(GSE9471.bm.table) <- c("Affymetrix.ID", "Ensembl.ID", "EntrezGene.ID", "Symbol", "Definition")
GSE9471.bm.filter <- filter(GSE9471.bm.table, !is.na(Symbol))

GSE9471.annot <- GSE9471.norm[GSE9471.bm.filter$Affymetrix.ID,]
GSE9471.expr <- exprs(GSE9471.annot)
rownames(GSE9471.expr) <- 1:nrow(GSE9471.annot)
GSE9471.collapse.expr <- collapseRows(GSE9471.expr, GSE9471.bm.filter$Symbol, rownames(GSE9471.expr)) #collapseRows by symbol
GSE9471.collapse <- ExpressionSet(assayData = GSE9471.collapse.expr$datETcollapsed, phenoData = phenoData(GSE9471.annot))

GSE9471.collapse$Wake <- grepl("ZT15|ZT21", GSE9471.collapse$Time.Point)
GSE9471.model <- model.matrix( ~ Time.Point, pData(GSE9471.collapse))
colnames(GSE9471.model) %<>% str_replace("Condition", "") %>% make.names
GSE9471.lm <- lmFit(GSE9471.collapse, GSE9471.model) %>% eBayes

toptable.GSE9470 <- topTable(GSE9471.lm, coef = 2, n = Inf) %>% select(AveExpr, logFC, t:B)
toptable.GSE9470$Symbol <- rownames(toptable.GSE9471)
SaveRDSgz(toptable.GSE9470, "./save/toptable.GSE9471.rda")

