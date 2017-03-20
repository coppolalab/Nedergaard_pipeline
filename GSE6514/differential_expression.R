library(Biobase)
library(annotate)
library(WGCNA)
library(biomaRt)
library(broom)
library(limma)

library(affy)
library(GEOquery)

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
GSE6514.raw <- getGEO("GSE6514", destdir = "./save")[[1]]

supplemental.files <- getGEOSuppFiles('GSE6514', baseDir = '.')
GSE6514.affy <- read.affybatch(list.files("./GSE6514/GSE6514_RAW", full.names = TRUE))
sampleNames(GSE6514.affy) %<>% str_replace("\\..*$", "")

pheno.data <- pData(GSE6514.raw) %>% select(title, characteristics_ch1.3, geo_accession) 
colnames(pheno.data) <- c("Mess", "Tissue", "GEO.Accession")

pheno.data$Tissue %<>% str_replace("Tissue: ", "") %>% factor
pheno.data$Mess %<>% str_replace("_#.*$", "") %>% str_replace("^+[0-9]_", "")

pheno.cortex <- filter(pheno.data, Tissue == "cerebral cortex")
pheno.cortex$Group <- str_replace(pheno.cortex$Mess, "cerebral_cortex_", "") %>% factor

pheno.hypo <- filter(pheno.data, Tissue == "hypothalamus")
pheno.hypo$Group <- str_replace(pheno.hypo$Mess, "hypothalamus_", "") %>% factor

pheno.combined <- rbind(pheno.cortex, pheno.hypo)
rownames(pheno.combined) <- pheno.combined$GEO.Accession

pData(GSE6514.affy) <- pheno.combined

affy.norm <- rma(GSE6514.affy)

ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
bm.table <- getBM(attributes = c('affy_mouse430_2', 'ensembl_gene_id', 'entrezgene', 'external_gene_name', 'description'), filters = 'affy_mouse430_2', values = featureNames(affy.norm), mart = ensembl)
bm.table$description %<>% str_replace_all(" \\[.*\\]$", "")
colnames(bm.table) <- c("Affymetrix.ID", "Ensembl.ID", "EntrezGene.ID", "Symbol", "Definition")
bm.filter <- filter(bm.table, !is.na(Symbol))

affy.annot <- affy.norm[bm.filter$Affymetrix.ID,]
affy.expr <- exprs(affy.annot)
rownames(affy.expr) <- 1:nrow(affy.annot)
collapse.expr <- collapseRows(affy.expr, bm.filter$Symbol, rownames(affy.expr)) #collapseRows by symbol
collapse.affy <- ExpressionSet(assayData = collapse.expr$datETcollapsed, phenoData = phenoData(affy.annot))

affy.cortex <- collapse.affy[,collapse.affy$Tissue == "cerebral cortex"]
affy.hypo <- collapse.affy[,collapse.affy$Tissue == "hypothalamus"]

model.cortex <- model.matrix( ~ 0 + Group, pData(affy.cortex))
fit.cortex <- makeContrasts("Group3hrs_sleep_deprivation - Group3hrs_sleep", "Group6hrs_sleep_deprivation - Group6hrs_sleep", "Group9hrs_sleep_deprivation - Group9hrs_sleep", "Group12hrs_sleep_deprivation - Group12hrs_sleep", levels = model.cortex)
lm.cortex <- lmFit(affy.cortex, model.cortex) %>% contrasts.fit(fit.cortex) %>% eBayes

toptable.cortex.3hrs <- topTable(lm.cortex, coef = 1, n = Inf) %>% select(AveExpr, logFC, t:B)
colnames(toptable.cortex.3hrs)[2:ncol(toptable.cortex.3hrs)] %<>% str_c('_3hrs')
toptable.cortex.3hrs$Symbol <- rownames(toptable.cortex.3hrs)

toptable.cortex.6hrs <- topTable(lm.cortex, coef = 2, n = Inf) %>% select(AveExpr, logFC, t:B)
colnames(toptable.cortex.6hrs)[2:ncol(toptable.cortex.6hrs)] %<>% str_c('_6hrs')
toptable.cortex.6hrs$Symbol <- rownames(toptable.cortex.6hrs)

toptable.cortex.9hrs <- topTable(lm.cortex, coef = 3, n = Inf) %>% select(AveExpr, logFC, t:B)
colnames(toptable.cortex.9hrs)[2:ncol(toptable.cortex.9hrs)] %<>% str_c('_9hrs')
toptable.cortex.9hrs$Symbol <- rownames(toptable.cortex.9hrs)

toptable.cortex.12hrs <- topTable(lm.cortex, coef = 4, n = Inf) %>% select(AveExpr, logFC, t:B)
colnames(toptable.cortex.12hrs)[2:ncol(toptable.cortex.12hrs)] %<>% str_c('_12hrs')
toptable.cortex.12hrs$Symbol <- rownames(toptable.cortex.12hrs)

toptable.cortex <- left_join(toptable.cortex.3hrs, toptable.cortex.6hrs) %>% 
    left_join(toptable.cortex.9hrs) %>% 
    left_join(toptable.cortex.12hrs) %>%
    select(Symbol, AveExpr, dplyr::contains("logFC"), dplyr::contains("P.Value"), dplyr::contains("adj.P.Val")) 
SaveRDSgz(toptable.cortex, "./save/toptable.cortex.rda")

model.hypo <- model.matrix( ~ 0 + Group, pData(affy.hypo))
fit.hypo <- makeContrasts("Group3hrs_sleep_deprivation - Group3hrs_sleep", "Group6hrs_sleep_deprivation - Group6hrs_sleep", "Group9hrs_sleep_deprivation - Group9hrs_sleep", "Group12hrs_sleep_deprivation - Group12hrs_sleep", levels = model.hypo)
lm.hypo <- lmFit(affy.hypo, model.hypo) %>% contrasts.fit(fit.hypo) %>% eBayes

toptable.hypo.3hrs <- topTable(lm.hypo, coef = 1, n = Inf) %>% select(AveExpr, logFC, t:B)
colnames(toptable.hypo.3hrs)[2:ncol(toptable.hypo.3hrs)] %<>% str_c('_3hrs')
toptable.hypo.3hrs$Symbol <- rownames(toptable.hypo.3hrs)

toptable.hypo.6hrs <- topTable(lm.hypo, coef = 2, n = Inf) %>% select(AveExpr, logFC, t:B)
colnames(toptable.hypo.6hrs)[2:ncol(toptable.hypo.6hrs)] %<>% str_c('_6hrs')
toptable.hypo.6hrs$Symbol <- rownames(toptable.hypo.6hrs)

toptable.hypo.9hrs <- topTable(lm.hypo, coef = 3, n = Inf) %>% select(AveExpr, logFC, t:B)
colnames(toptable.hypo.9hrs)[2:ncol(toptable.hypo.9hrs)] %<>% str_c('_9hrs')
toptable.hypo.9hrs$Symbol <- rownames(toptable.hypo.9hrs)

toptable.hypo.12hrs <- topTable(lm.hypo, coef = 4, n = Inf) %>% select(AveExpr, logFC, t:B)
colnames(toptable.hypo.12hrs)[2:ncol(toptable.hypo.12hrs)] %<>% str_c('_12hrs')
toptable.hypo.12hrs$Symbol <- rownames(toptable.hypo.12hrs)

toptable.hypo <- left_join(toptable.hypo.3hrs, toptable.hypo.6hrs) %>% 
    left_join(toptable.hypo.9hrs) %>% 
    left_join(toptable.hypo.12hrs) %>%
    select(Symbol, AveExpr, dplyr::contains("logFC"), dplyr::contains("P.Value"), dplyr::contains("adj.P.Val")) 
SaveRDSgz(toptable.hypo, "./save/toptable.hypo.rda")

