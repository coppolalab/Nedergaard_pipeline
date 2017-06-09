library(RRHO)

library(openxlsx)
library(Cairo)
library(biomaRt)

library(R.utils)
library(vadr)
#Always load tidyverse last
library(magrittr)
library(stringr)
library(tidyverse)

GetRRHO <- function(colname1, colname2, dataset1, dataset2, symbolname, output.suffix, stepsize = 100) {
    print(colname1)
    print(colname2)
    subset1 <- select_(dataset1, symbolname, str_c("Log.PValue_", colname1))
    subset2 <- select_(dataset2, symbolname, str_c("Log.PValue_", colname2))
    output.dir <- str_c("./rrho", output.suffix, sep = "_")
    dir.create(output.dir, showWarnings = FALSE)
    rrho.out <- RRHO(subset1, subset2, alternative = "enrichment", stepsize = stepsize, labels = c(colname1, colname2), plots = TRUE, outputdir = output.dir, BY = TRUE)
    return(rrho.out)
}

LogPvalue <- function(contrast, top.table) {
    logfc.colname <- str_c("logFC_", contrast)
    pvalue.colname <- str_c("P.Value_", contrast)
    log.pvalue.column <- -log10(top.table[[pvalue.colname]]) * sign(top.table[[logfc.colname]])
    log.pvalue.column
}

source("../../code/common_functions.R")

#Data Preparation

#Nedergaard
nedergaard.top <- ReadRDSgz("../differential_expression/save/top.table.all.rda") %>% select(Symbol, dplyr::contains("logFC"), dplyr::contains("P.Value"))
nedergaard.contrasts <- str_subset(colnames(nedergaard.top), "logFC") %>% str_replace("logFC_", "")
nedergaard.logpval <- map(nedergaard.contrasts, LogPvalue, nedergaard.top) %>% reduce(cbind) %>% set_colnames(str_c("Log.PValue_", nedergaard.contrasts))
nedergaard.reduce <- data.frame(Symbol = nedergaard.top$Symbol, nedergaard.logpval)

#Cirelli 2004 - rejected
cirelli_wake_cortex <- read.delim("../cirelli_2004_neuron/cirelli_2004_wake_cortex", stringsAsFactors = FALSE) %>% 
    as_tibble %>% 
    mutate_all(str_replace_all, pattern = " ", replacement = "")
cirelli_wake_cerebellum <- read.delim("../cirelli_2004_neuron/cirelli_2004_wake_cerebellum", stringsAsFactors = FALSE) %>% 
    as_tibble %>% 
    mutate_all(str_replace_all, pattern = " ", replacement = "")
cirelli_sleep_cortex <- read.delim("../cirelli_2004_neuron/cirelli_2004_sleep_cortex", stringsAsFactors = FALSE) %>% 
    as_tibble %>% 
    mutate_all(str_replace_all, pattern = " ", replacement = "")
cirelli_sleep_cerebellum <- read.delim("../cirelli_2004_neuron/cirelli_2004_sleep_cerebellum", stringsAsFactors = FALSE) %>% 
    as_tibble %>% 
    mutate_all(str_replace_all, pattern = " ", replacement = "")

ensembl <- useMart("ensembl", dataset = "rnorvegicus_gene_ensembl")
ensembl.attributes <- listAttributes(ensembl)
ensembl.features <- listFilters(ensembl)
wake_cortex_ensembl <- getBM(attributes = c('embl', 'external_gene_name'), filters = 'embl', values = cirelli_wake_cortex$Gene, mart = ensembl)

#Gerstner
gerstner_hypo <- read.delim("../Gerstner_2006_Neuroscience/Gerstner_2006_hypo", stringsAsFactors = FALSE) %>%
    as_tibble

GetHyper <- function(group.name, list.genes, top.table) {
  total.genes <- nrow(full.genelist)
  background.genes <- total.genes - length(list.genes) - total.ppi + overlap
  
  hyper.pval <- phyper(overlap, total.ppi, background.genes, length(list.genes), lower.tail = FALSE) %>% signif(3)
  hyper.pval
}

GetSig <- function(top.table, group.name, n.genes) {
    arrange.cond <- str_c("P.Value_", group.name)
    top.sig <- arrange_(top.table, arrange.cond) %>% slice(n.genes)
    top.sig$Symbol
}

#GSE6514
GSE6514.cortex <- ReadRDSgz("../GSE6514/save/toptable.cortex.rda") %>% select(Symbol, dplyr::contains("logFC"), dplyr::contains("P.Value"))
colnames(GSE6514.cortex)[2:ncol(GSE6514.cortex)] %<>% str_c("_cortex")
GSE6514.hypo <- ReadRDSgz("../GSE6514/save/toptable.hypo.rda") %>% select(Symbol, dplyr::contains("logFC"), dplyr::contains("P.Value"))
colnames(GSE6514.hypo)[2:ncol(GSE6514.hypo)] %<>% str_c("_hypo")
GSE6514.combined <- left_join(GSE6514.cortex, GSE6514.hypo)
GSE6514.contrasts <- str_subset(colnames(GSE6514.combined), "logFC") %>% str_replace("logFC_", "")
GSE6514.logpval <- map(GSE6514.contrasts, LogPvalue, GSE6514.combined) %>% reduce(cbind) %>% set_colnames(str_c("Log.PValue_", GSE6514.contrasts))
GSE6514.reduce <- data.frame(Symbol = GSE6514.combined$Symbol, GSE6514.logpval)

#GSE9441
GSE9441.brain <- ReadRDSgz("../GSE9444/save/toptable.GSE9441.brain.rda") %>% select(Symbol, dplyr::contains("logFC"), dplyr::contains("P.Value"))
colnames(GSE9441.brain)[2:ncol(GSE9441.brain)] %<>% str_c("_brain")
GSE9441.liver <- ReadRDSgz("../GSE9444/save/toptable.GSE9441.liver.rda") %>% select(Symbol, dplyr::contains("logFC"), dplyr::contains("P.Value"))
colnames(GSE9441.liver)[2:ncol(GSE9441.liver)] %<>% str_c("_liver")
GSE9441.combined <- left_join(GSE9441.brain, GSE9441.liver)
GSE9441.contrasts <- str_subset(colnames(GSE9441.combined), "logFC") %>% str_replace("logFC_", "")
GSE9441.logpval <- map(GSE9441.contrasts, LogPvalue, GSE9441.combined) %>% reduce(cbind) %>% set_colnames(str_c("Log.PValue_", GSE9441.contrasts))
GSE9441.reduce <- data.frame(Symbol = GSE9441.combined$Symbol, GSE9441.logpval)

#GSE9442
GSE9442.akr <- ReadRDSgz("../GSE9444/save/toptable.GSE9442.akr.rda") %>% select(Symbol, dplyr::contains("logFC"), dplyr::contains("P.Value"))
colnames(GSE9442.akr)[2:ncol(GSE9442.akr)] %<>% str_c("_akr") %>% str_replace("logFC\\.", "logFC_") %>% str_replace("P\\.Value\\.", "P.Value_")
GSE9442.c57 <- ReadRDSgz("../GSE9444/save/toptable.GSE9442.c57.rda") %>% select(Symbol, dplyr::contains("logFC"), dplyr::contains("P.Value"))
colnames(GSE9442.c57)[2:ncol(GSE9442.c57)] %<>% str_c("_c57") %>% str_replace("logFC\\.", "logFC_") %>% str_replace("P\\.Value\\.", "P.Value_")
GSE9442.dba <- ReadRDSgz("../GSE9444/save/toptable.GSE9442.dba.rda") %>% select(Symbol, dplyr::contains("logFC"), dplyr::contains("P.Value"))
colnames(GSE9442.dba)[2:ncol(GSE9442.dba)] %<>% str_c("_dba") %>% str_replace("logFC\\.", "logFC_") %>% str_replace("P\\.Value\\.", "P.Value_")
GSE9442.combined <- left_join(GSE9442.akr, GSE9442.c57) %>% left_join(GSE9442.dba)
GSE9442.contrasts <- str_subset(colnames(GSE9442.combined), "logFC") %>% str_replace("logFC_", "")
GSE9442.logpval <- map(GSE9442.contrasts, LogPvalue, GSE9442.combined) %>% reduce(cbind) %>% set_colnames(str_c("Log.PValue_", GSE9442.contrasts))
GSE9442.reduce <- data.frame(Symbol = GSE9442.combined$Symbol, GSE9442.logpval)

#GSE9443
GSE9443.brain <- ReadRDSgz("../GSE9444/save/toptable.GSE9443.brain.rda") %>% select(Symbol, dplyr::contains("logFC"), dplyr::contains("P.Value"))
colnames(GSE9443.brain)[2:ncol(GSE9443.brain)] %<>% str_replace("logFC\\.", "logFC_") %>% str_replace("P\\.Value\\.", "P.Value_")
GSE9443.homer1 <- ReadRDSgz("../GSE9444/save/toptable.GSE9443.homer1.rda") %>% select(Symbol, dplyr::contains("logFC"), dplyr::contains("P.Value"))
colnames(GSE9443.homer1)[2:ncol(GSE9443.homer1)] %<>% str_replace("logFC\\.", "logFC_") %>% str_replace("P\\.Value\\.", "P.Value_")
GSE9443.combined <- left_join(GSE9443.brain, GSE9443.homer1)
GSE9443.contrasts <- str_subset(colnames(GSE9443.combined), "logFC") %>% str_replace("logFC_", "")
GSE9443.logpval <- map(GSE9443.contrasts, LogPvalue, GSE9443.combined) %>% reduce(cbind) %>% set_colnames(str_c("Log.PValue_", GSE9443.contrasts))
GSE9443.reduce <- data.frame(Symbol = GSE9443.combined$Symbol, GSE9443.logpval)

#GSE9471
GSE9471.top <- ReadRDSgz("../GSE9471/save/toptable.GSE9471.rda") %>% select(Symbol, dplyr::contains("logFC"), dplyr::contains("P.Value"))
GSE9471.top$Log.PValue_sw <- -log10(GSE9471.top$P.Value) * sign(GSE9471.top$logFC)

#GSE48369
GSE48369.bound <- ReadRDSgz("../GSE48369/save/toptable.bound.all.rda") %>% select(Symbol, dplyr::contains("logFC"), dplyr::contains("P.Value"))
colnames(GSE48369.bound)[2:ncol(GSE48369.bound)] %<>% str_c("_bound") %>% str_replace("logFC\\.", "logFC_") %>% str_replace("P\\.Value\\.", "P.Value_")
GSE48369.unbound <- ReadRDSgz("../GSE48369/save/toptable.unbound.all.rda") %>% select(Symbol, dplyr::contains("logFC"), dplyr::contains("P.Value"))
colnames(GSE48369.unbound)[2:ncol(GSE48369.unbound)] %<>% str_c("_unbound") %>% str_replace("logFC\\.", "logFC_") %>% str_replace("P\\.Value\\.", "P.Value_")
GSE48369.combined <- left_join(GSE48369.bound, GSE48369.unbound)
GSE48369.contrasts <- str_subset(colnames(GSE48369.combined), "logFC") %>% str_replace("logFC_", "")
GSE48369.logpval <- map(GSE48369.contrasts, LogPvalue, GSE48369.combined) %>% reduce(cbind) %>% set_colnames(str_c("Log.PValue_", GSE48369.contrasts))
GSE48369.reduce <- data.frame(Symbol = GSE48369.combined$Symbol, GSE48369.logpval)

#GSE69079
GSE69079.top <- ReadRDSgz("../GSE69079/save/toptable.bound.all.rda") %>% select(Symbol, dplyr::contains("logFC"), dplyr::contains("P.Value"))
colnames(GSE69079.top)[2:ncol(GSE69079.top)] %<>% str_replace("logFC\\.", "logFC_") %>% str_replace("P\\.Value\\.", "P.Value_")
GSE69079.contrasts <- str_subset(colnames(GSE69079.top), "logFC") %>% str_replace("logFC_", "")
GSE69079.logpval <- map(GSE69079.contrasts, LogPvalue, GSE69079.top) %>% reduce(cbind) %>% set_colnames(str_c("Log.PValue_", GSE69079.contrasts))
GSE69079.reduce <- data.frame(Symbol = GSE69079.top$Symbol, GSE69079.logpval)

#GSE78215
GSE78215.top <- ReadRDSgz("../GSE78215/save/toptable.bound.all.rda") %>% select(Symbol, dplyr::contains("logFC"), dplyr::contains("P.Value"))
colnames(GSE78215.top)[2:ncol(GSE78215.top)] %<>% str_replace("logFC\\.", "logFC_") %>% str_replace("P\\.Value\\.", "P.Value_")
GSE78215.contrasts <- str_subset(colnames(GSE78215.top), "logFC") %>% str_replace("logFC_", "")
GSE78215.logpval <- map(GSE78215.contrasts, LogPvalue, GSE78215.top) %>% reduce(cbind) %>% set_colnames(str_c("Log.PValue_", GSE78215.contrasts))
GSE78215.reduce <- data.frame(Symbol = GSE78215.top$Symbol, GSE78215.logpval)

#Overlaps
#GSE6514
GSE6514.shared <- intersect(nedergaard.reduce$Symbol, GSE6514.reduce$Symbol)
GSE6514.filtered <- filter(GSE6514.reduce, Symbol %in% GSE6514.shared)
nedergaard.GSE6514 <- filter(nedergaard.reduce, Symbol %in% GSE6514.shared)

GSE6514.rrho <- map(GSE6514.contrasts, mkchain(map(nedergaard.contrasts, GetRRHO, ., nedergaard.GSE6514, GSE6514.filtered, "Symbol", "GSE6514", 100))) %>% flatten
GSE6514.rrho.logpval <- map(GSE6514.rrho, extract2, "hypermat.by") %>% map_dbl(max) 
GSE6514.rrho.pval <- exp(-GSE6514.rrho.logpval) %>% p.adjust(method = "fdr") %>% signif(3)

dim(GSE6514.rrho.pval) <- c(3,8)
rownames(GSE6514.rrho.pval) <- str_replace_all(nedergaard.contrasts, "_", " ")
colnames(GSE6514.rrho.pval) <- str_replace_all(GSE6514.contrasts, "_", " ")
GSE6514.rrho.out <- cbind(Group = rownames(GSE6514.rrho.pval), GSE6514.rrho.pval)

write.csv(GSE6514.rrho.out, "GSE6514.rrho.csv", quote = FALSE, row.names = FALSE)

#GSE9441
GSE9441.shared <- intersect(nedergaard.reduce$Symbol, GSE9441.reduce$Symbol)
GSE9441.filtered <- filter(GSE9441.reduce, Symbol %in% GSE9441.shared)
nedergaard.GSE9441 <- filter(nedergaard.reduce, Symbol %in% GSE9441.shared)

GSE9441.rrho <- map(GSE9441.contrasts, mkchain(map(nedergaard.contrasts, GetRRHO, ., nedergaard.GSE9441, GSE9441.filtered, "Symbol", "GSE9441", 100))) %>% flatten
GSE9441.rrho.logpval <- map(GSE9441.rrho, extract2, "hypermat.by") %>% map_dbl(max) 
GSE9441.rrho.pval <- exp(-GSE9441.rrho.logpval) %>% p.adjust(method = "fdr") %>% signif(3)

dim(GSE9441.rrho.pval) <- c(3,6)
rownames(GSE9441.rrho.pval) <- str_replace_all(nedergaard.contrasts, "_", " ")
GSE9441.contrasts.format <- str_replace_all(GSE9441.contrasts, "akr", "AKR/J") %>% str_replace_all("c57", "C57/BL6") %>% str_replace_all("dba", "DBA/2J")
colnames(GSE9441.rrho.pval) <- str_replace_all(GSE9441.contrasts.format, "_", " ")
GSE9441.rrho.out <- cbind(Group = rownames(GSE9441.rrho.pval), GSE9441.rrho.pval)

write.csv(GSE9441.rrho.out, "GSE9441.rrho.csv", quote = FALSE, row.names = FALSE)

#GSE9442
GSE9442.shared <- intersect(nedergaard.reduce$Symbol, GSE9442.reduce$Symbol)
GSE9442.filtered <- filter(GSE9442.reduce, Symbol %in% GSE9442.shared)
nedergaard.GSE9442 <- filter(nedergaard.reduce, Symbol %in% GSE9442.shared)

GSE9442.rrho <- map(GSE9442.contrasts, mkchain(map(nedergaard.contrasts, GetRRHO, ., nedergaard.GSE9442, GSE9442.filtered, "Symbol", "GSE9442", 100))) %>% flatten
GSE9442.rrho.logpval <- map(GSE9442.rrho, extract2, "hypermat.by") %>% map_dbl(max) 
GSE9442.rrho.pval <- exp(-GSE9442.rrho.logpval) %>% p.adjust(method = "fdr") %>% signif(3)

dim(GSE9442.rrho.pval) <- c(3,12)
rownames(GSE9442.rrho.pval) <- str_replace_all(nedergaard.contrasts, "_", " ")
GSE9442.contrasts.format <- str_replace_all(GSE9442.contrasts, "akr", "AKR/J") %>% str_replace_all("c57", "C57/BL6") %>% str_replace_all("dba", "DBA/2J") %>% str_replace_all("zt", "ZT")
colnames(GSE9442.rrho.pval) <- str_replace_all(GSE9442.contrasts.format, "_", " ")

GSE9442.rrho.akr <- GSE9442.rrho.pval[,1:4]
GSE9442.akr.out <- cbind(Group = rownames(GSE9442.rrho.akr), GSE9442.rrho.akr)

GSE9442.rrho.c57 <- GSE9442.rrho.pval[,5:8]
GSE9442.c57.out <- cbind(Group = rownames(GSE9442.rrho.c57), GSE9442.rrho.c57)

GSE9442.rrho.dba <- GSE9442.rrho.pval[,9:12]
GSE9442.dba.out <- cbind(Group = rownames(GSE9442.rrho.dba), GSE9442.rrho.dba)

write.csv(GSE9442.akr.out, "GSE9442.rrho.akr.csv", quote = FALSE, row.names = FALSE)
write.csv(GSE9442.c57.out, "GSE9442.rrho.c57.csv", quote = FALSE, row.names = FALSE)
write.csv(GSE9442.dba.out, "GSE9442.rrho.dba.csv", quote = FALSE, row.names = FALSE)

#GSE9443
GSE9443.shared <- intersect(nedergaard.reduce$Symbol, GSE9443.reduce$Symbol)
GSE9443.filtered <- filter(GSE9443.reduce, Symbol %in% GSE9443.shared)
nedergaard.GSE9443 <- filter(nedergaard.reduce, Symbol %in% GSE9443.shared)

GSE9443.rrho <- map(GSE9443.contrasts, mkchain(map(nedergaard.contrasts, GetRRHO, ., nedergaard.GSE9443, GSE9443.filtered, "Symbol", "GSE9443", 100))) %>% flatten
GSE9443.rrho.logpval <- map(GSE9443.rrho, extract2, "hypermat.by") %>% map_dbl(max) 
GSE9443.rrho.pval <- exp(-GSE9443.rrho.logpval) %>% p.adjust(method = "fdr") %>% signif(3)

dim(GSE9443.rrho.pval) <- c(3,2)
rownames(GSE9443.rrho.pval) <- str_replace_all(nedergaard.contrasts, "_", " ")
GSE9443.contrasts.format <- capitalize(GSE9443.contrasts) #%>% str_c("GSE9443_", .)
colnames(GSE9443.rrho.pval) <- str_replace_all(GSE9443.contrasts.format, "_", " ")
GSE9443.rrho.out <- cbind(Group = rownames(GSE9443.rrho.pval), GSE9443.rrho.pval)

write.csv(GSE9443.rrho.out, "GSE9443.rrho.csv", quote = FALSE, row.names = FALSE)

#GSE9471
GSE9471.shared <- intersect(nedergaard.reduce$Symbol, GSE9471.reduce$Symbol)
GSE9471.filtered <- filter(GSE9471.top, Symbol %in% GSE9471.shared)
nedergaard.GSE9471 <- filter(nedergaard.reduce, Symbol %in% GSE9471.shared)

GSE9471.rrho <- map(nedergaard.contrasts, GetRRHO, "sw", nedergaard.GSE9471, data.frame(GSE9471.filtered), "Symbol", "GSE9471", 100)
GSE9471.rrho.logpval <- map(GSE9471.rrho, extract2, "hypermat.by") %>% map_dbl(max) 
GSE9471.rrho.pval <- exp(-GSE9471.rrho.logpval) %>% p.adjust(method = "fdr") %>% signif(3)
names(GSE9471.rrho.pval) <- str_replace_all(nedergaard.contrasts, "_", " ")

write.csv(t(GSE9471.rrho.pval), "GSE9471.rrho.csv", quote = FALSE, row.names = FALSE) 

#GSE48369
GSE48369.shared <- intersect(nedergaard.reduce$Symbol, GSE48369.reduce$Symbol)
GSE48369.filtered <- filter(GSE48369.reduce, Symbol %in% GSE48369.shared)
nedergaard.GSE48369 <- filter(nedergaard.reduce, Symbol %in% GSE48369.shared)

GSE48369.rrho <- map(GSE48369.contrasts, mkchain(map(nedergaard.contrasts, GetRRHO, ., nedergaard.GSE48369, GSE48369.filtered, "Symbol", "GSE48369", 100))) %>% flatten
GSE48369.rrho.logpval <- map(GSE48369.rrho, extract2, "hypermat.by") %>% map_dbl(max) 
GSE48369.rrho.pval <- exp(-GSE48369.rrho.logpval) %>% p.adjust(method = "fdr") %>% signif(3) 

dim(GSE48369.rrho.pval) <- c(3,6)
rownames(GSE48369.rrho.pval) <- str_replace_all(nedergaard.contrasts, "_", " ")
GSE48369.contrasts.format <- str_replace_all(GSE48369.contrasts, "sds", "Sleep_dep_vs_Sleep") %>% str_replace_all("sdw", "Sleep_dep_vs_Wake") %>% str_replace_all("ws", "Wake_vs_Sleep")
colnames(GSE48369.rrho.pval) <- str_replace_all(GSE48369.contrasts.format, "_", " ")
GSE48369.rrho.out <- cbind(Group = rownames(GSE48369.rrho.pval), GSE48369.rrho.pval)
GSE48369.rrho.bound <- GSE48369.rrho.out[,1:4]
GSE48369.rrho.unbound <- GSE48369.rrho.out[,c(1, 5:7)]

write.csv(GSE48369.rrho.bound, "GSE48369.rrho.bound.csv", quote = FALSE, row.names = FALSE)
write.csv(GSE48369.rrho.unbound, "GSE48369.rrho.unbound.csv", quote = FALSE, row.names = FALSE)

#GSE69079
GSE69079.shared <- intersect(nedergaard.reduce$Symbol, GSE69079.reduce$Symbol)
GSE69079.filtered <- filter(GSE69079.reduce, Symbol %in% GSE69079.shared)
nedergaard.GSE69079 <- filter(nedergaard.reduce, Symbol %in% GSE69079.shared)

GSE69079.rrho <- map(GSE69079.contrasts, mkchain(map(nedergaard.contrasts, GetRRHO, ., nedergaard.GSE69079, GSE69079.filtered, "Symbol", "GSE69079", 100))) %>% flatten
GSE69079.rrho.logpval <- map(GSE69079.rrho, extract2, "hypermat.by") %>% map_dbl(max) 
GSE69079.rrho.pval <- exp(-GSE69079.rrho.logpval) %>% p.adjust(method = "fdr") %>% signif(3)

dim(GSE69079.rrho.pval) <- c(3,3)
rownames(GSE69079.rrho.pval) <- str_replace_all(nedergaard.contrasts, "_", " ")
GSE69079.contrasts.format <- str_replace_all(GSE69079.contrasts, "sds", "Sleep_dep_vs_Sleep") %>% str_replace_all("sdw", "Sleep_dep_vs_Wake") %>% str_replace_all("ws", "Wake_vs_Sleep")
colnames(GSE69079.rrho.pval) <- str_replace_all(GSE69079.contrasts.format, "_", " ")
GSE69079.rrho.out <- cbind(Group = rownames(GSE69079.rrho.pval), GSE69079.rrho.pval)

write.csv(GSE69079.rrho.out, "GSE69079.rrho.csv", quote = FALSE, row.names = FALSE)

#GSE78215
GSE78215.shared <- intersect(nedergaard.reduce$Symbol, GSE78215.reduce$Symbol)
GSE78215.filtered <- filter(GSE78215.reduce, Symbol %in% GSE78215.shared)
nedergaard.GSE78215 <- filter(nedergaard.reduce, Symbol %in% GSE78215.shared)

GSE78215.rrho <- map(GSE78215.contrasts, mkchain(map(nedergaard.contrasts, GetRRHO, ., nedergaard.GSE78215, GSE78215.filtered, "Symbol", "GSE78215", 100))) %>% flatten
GSE78215.rrho.logpval <- map(GSE78215.rrho, extract2, "hypermat.by") %>% map_dbl(max) 
GSE78215.rrho.pval <- exp(-GSE78215.rrho.logpval) %>% p.adjust(method = "fdr") %>% signif(3)

dim(GSE78215.rrho.pval) <- c(3,4)
rownames(GSE78215.rrho.pval) <- str_replace_all(nedergaard.contrasts, "_", " ")
GSE78215.contrasts.format <- str_replace_all(GSE78215.contrasts, "zt7", "6hrSD_1hrRecov_vs_Cont") %>% str_replace_all("zt8_1", "6hrSD_2hrRecov_vs_Cont") %>% str_replace_all("zt8_2", "5hrSD_3hrRecov_vs_Cont") %>% str_replace_all("zt11", "5hrSD_6hrRecov_vs_Cont")
colnames(GSE78215.rrho.pval) <- str_replace_all(GSE78215.contrasts.format, "_", " ")

GSE78215.rrho.out <- cbind(Group = rownames(GSE78215.rrho.pval), GSE78215.rrho.pval)
write.csv(GSE78215.rrho.out, "GSE78215.rrho.csv", quote = FALSE, row.names = FALSE)
