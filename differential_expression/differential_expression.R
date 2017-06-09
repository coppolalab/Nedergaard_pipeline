library(Biobase)
library(annotate)
library(WGCNA)
library(biomaRt)
library(matrixStats)
library(limma)
library(edgeR)
library(siggenes)

library(openxlsx)
library(Cairo)
library(UpSetR)
library(RColorBrewer)

library(broom)
library(rlist)
library(stringr)
library(magrittr)
library(tidyverse)

FormatTopTable <- function(coef.num, col.suffix, limma.object) {
    top.table <- topTable(limma.object, coef = coef.num, n = Inf) %>% signif(3) %>% data.frame
    colnames(top.table)[c(1,3:ncol(top.table))] %<>% str_c(col.suffix, sep = "_")
    top.table$Symbol <- rownames(top.table)
    top.table
}

GenWorkbook <- function(dataset, filename, pval.name, fdr.name, log.name) {
    pval.cols <- colnames(dataset) %>% str_detect(pval.name) %>% which
    adj.pval.cols <- colnames(dataset) %>% str_detect(fdr.name) %>% which
    logfc.cols <- colnames(dataset) %>% str_detect(log.name) %>% which
    description.cols <- colnames(dataset) %>% str_detect("Description") %>% which

    wb <- createWorkbook()
    addWorksheet(wb = wb, sheetName = "Sheet 1", gridLines = TRUE)
    writeDataTable(wb = wb, sheet = 1, x = dataset)
    sig.pvalues <- createStyle(fontColour = "red")
    conditionalFormatting(wb, 1, cols = pval.cols, rows = 1:nrow(dataset), rule = "<0.005", style = sig.pvalues)
    conditionalFormatting(wb, 1, cols = adj.pval.cols, rows = 1:nrow(dataset), rule = "<0.05", style = sig.pvalues)
    conditionalFormatting(wb, 1, cols = logfc.cols, rows = 1:nrow(dataset), style = c("#63BE7B", "white", "red"), type = "colourScale")
    setColWidths(wb, 1, cols = 1, widths = 15)
    setColWidths(wb, 1, cols = 3, widths = 15)
    setColWidths(wb, 1, cols = 4:ncol(dataset), widths = "auto")
    setColWidths(wb, 1, cols = description.cols, widths = 45)
    pageSetup(wb, 1, orientation = "landscape", fitToWidth = TRUE)
    freezePane(wb, 1, firstRow = TRUE)
    showGridLines(wb, 1, showGridLines = TRUE)

    saveWorkbook(wb, filename, overwrite = TRUE) 
}

#Enrichr
EnrichrSubmit <- function(column.suffix, top.table, pval.prefix, cutoff, enrichr.terms) {
    pval.column <- str_c(pval.prefix, "_", column.suffix)
    log.column <- str_c("logFC_", column.suffix)
    filter.condition <- str_c(pval.column, " < ", cutoff)
    top.filter <- filter_(top.table, filter.condition) %>% select_("Symbol", pval.column, log.column)
    if (nrow(top.filter) > 1000) {
        top.filter <- arrange_(top.filter, str_c("desc(abs(", log.column, "))")) %>% slice(1:1000)
    }
    enrichr.data <- map(enrichr.terms, GetEnrichrData, top.filter)
    names(enrichr.data) <- enrichr.terms
    map(enrichr.terms, EnrichrWorkbook, enrichr.data, column.suffix)
}

EnrichrWorkbook <- function(subindex, full.df, comparison) {
    dataset <- full.df[[subindex]]
    wb <- createWorkbook()
    addWorksheet(wb = wb, sheetName = "sheet 1", gridLines = TRUE)
    writeDataTable(wb = wb, sheet = 1, x = dataset, withFilter = TRUE)
    freezePane(wb, 1, firstRow = TRUE)
    setColWidths(wb, 1, cols = c(1, 3:ncol(dataset)), widths = "auto")
    setColWidths(wb, 1, cols = 2, widths = 45)

    sub.dir <- file.path("./enrichr", comparison)
    dir.create(sub.dir, showWarnings = FALSE, recursive = TRUE)
    filename = str_c(sub.dir, "/", subindex, ".xlsx", sep = "")
    saveWorkbook(wb, filename, overwrite = TRUE) 
}

#MDS
MDSPlot <- function(filename, dataset, targetset, variablename) {
    dataset.plot <- data.frame(rownames(dataset$points), dataset$points)
    target.data <- data.frame(targetset$sample, factor(targetset[[variablename]]))
    colnames(target.data) <- c("Subject.ID", variablename)
    colnames(dataset.plot) <- c("Subject.ID", "Component.1", "Component.2")
    dataset.plot <- left_join(dataset.plot, target.data)

    p <- ggplot(dataset.plot, aes_string(x = "Component.1", y = "Component.2", col = variablename)) + geom_point() 
    p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5))
    p <- p + xlab("Component 1") + ylab("Component 2") + ggtitle(variablename) + theme(plot.background = element_blank(), legend.background = element_blank())
    CairoPDF(filename, height = 6, width = 7, bg = "transparent")
   print(p)
    dev.off()
}

VolcanoPlot <- function(column.suffix, top.table, threshold) {
    log.format <- str_c("logFC_", column.suffix)
    pval.format <- str_c("adj.P.Val_", column.suffix)
    top.table$Log.Pvalue <- -log10(top.table[[pval.format]])
    top.table$Significant <- top.table[[pval.format]] > threshold

    p <- ggplot(top.table, aes_string(x = log.format, y = "Log.Pvalue")) + geom_point(aes(color = Significant))
    p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + theme(legend.position = "none", plot.background = element_blank(), plot.title = element_text(hjust = 0.5))
    p <- p + xlab("Log Fold Change") + ylab("-Log10 FDR P-value") + ggtitle(column.suffix)
    CairoPDF(str_c("volcano_", column.suffix), width = 5, height = 5, bg = "transparent")
    print(p)
    dev.off()
}

DecidePlot <- function(file.name, decide.plot, y.lab, bar.padding = 100, pos.hjust = -0.4, neg.hjust = 1.3) {
    p <- ggplot()
    p <- p + geom_bar(data = filter(decide.plot, Direction == "positive"),  aes(x = Comparison, y = Num.Genes), stat = "identity", colour = "black", fill = "red", position = "dodge")   
    p <- p + geom_text(data = filter(decide.plot, Direction == "positive"), stat = "identity", size = 4, aes(x = Comparison, y = Num.Genes, ymax = max(Num.Genes) + nchar(max(Num.Genes)) * bar.padding, hjust = pos.hjust, label = Num.Genes), position = position_dodge(width = 1))
    p <- p + geom_bar(data = filter(decide.plot, Direction == "negative"),  aes(x = Comparison, y = Num.Genes), stat = "identity", colour = "black", fill = "green", position = "dodge") 
    p <- p + geom_text(data = filter(decide.plot, Direction == "negative"), stat = "identity", size = 4, aes(x = Comparison, y = Num.Genes, ymin = min(Num.Genes) - nchar(abs(min(Num.Genes))) * bar.padding, hjust = neg.hjust, label = abs(Num.Genes)), position = position_dodge(width = 1))
    p <- p + theme_bw() + coord_flip() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_blank()) 
    p <- p + theme(axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_text(hjust = 0)) + ylab(y.lab)
    p <- p + theme(panel.border = element_rect(size = 1, color = "black")) + facet_grid(Threshold + Total.Genes ~ .) 
    CairoPDF(file.name, width = 7, height = 6, bg = "transparent")
    print(p)
    dev.off()
}

GetSizes <- function(p.val, pval.column, log.column, dataset) {
    dataset.sig <- filter_(dataset, str_c(pval.column, " < ", p.val))
    dataset.up <- filter_(dataset.sig, str_c(log.column, " > 0"))
    dataset.down <- filter_(dataset.sig, str_c(log.column, " < 0"))
    return(c(positive = dim(dataset.up)[1], negative = -(dim(dataset.down)[1]), Threshold = p.val))
}

GetThresholds <- function(pval.column, log.column, dataset, thresholds, pval.label) {
    nums.table <- map(thresholds, GetSizes, pval.column, log.column, dataset) %>% reduce(rbind) %>% data.frame
    nums.table$Threshold <- str_c(pval.label, " < ", nums.table$Threshold)
    nums.table$Comparison <- str_replace(pval.column, '^.*?_', '')
    return(nums.table)
}

TopGenesPlot <- function(pval.column, top.table, dge.list, group) {
    format.name <- str_replace(pval.column, "P.Value_", "")
    top5.symbol <- arrange_(top.table, pval.column)$Symbol[1:5]
    top5.expr <- t(dge.list$E[top5.symbol,])
    top5.df <- data.frame(Group = group, top5.expr) %>% gather(Gene, Expression, -Group)
    top5.df$Gene %<>% factor(levels = top5.symbol)

    p <- ggplot(top5.df, aes(x = Group, y = Expression, color = Group)) + geom_boxplot() + geom_jitter() + theme_bw()
    p <- p + facet_wrap(~ Gene, ncol = 5, scales = "free") + theme(legend.position = "none")
    p <- p + theme(plot.title = element_text(hjust = 0.5), axis.title.x = element_blank())
    p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x = element_blank()) 
    p <- p + theme(plot.background = element_blank(), panel.border = element_rect(size = 1, color = "black"))
    p <- p + ggtitle(format.name)
    CairoPDF(str_c("top5.", format.name), height = 3.5, width = 16, bg = "transparent")
    print(p)
    dev.off()
}

EnrichrPlot <- function(enrichr.df, enrichr.expr, filename, plot.title, plot.height = 5, plot.width = 8) {
    enrichr.df$Gene.Count <- map(enrichr.df$Genes, str_split, ",") %>% map(unlist) %>% map_int(length)
    enrichr.df$Adj.P.value <- p.adjust(enrichr.df$P.value, method = "fdr")
    enrichr.df$Log.P.value <- log10(enrichr.df$Adj.P.value)
    log.column <- str_split_fixed(filename, ".", 2)[1]
    enrichr.updown <- map(enrichr.df$Genes, UpDown, enrichr.expr) %>% reduce(rbind)
    colnames(enrichr.updown) <- c("Down", "Up")
    enrichr.df <- cbind(enrichr.df, enrichr.updown)
    enrichr.df$Log.Up <- -enrichr.df$Log.P.value * enrichr.df$Up / enrichr.df$Gene.Count
    enrichr.df$Log.Down <- -enrichr.df$Log.P.value * enrichr.df$Down / enrichr.df$Gene.Count
    enrichr.df$Term %<>% str_replace_all("\\ \\(.*$", "") %>% str_replace_all("\\_Homo.*$", "")  #Remove any thing after the left parenthesis and convert to all lower case
    enrichr.df$Format.Name <- str_c(enrichr.df$Databas, ": ", enrichr.df$Term, " (", enrichr.df$Gene.Count, ")")
    enrichr.df %<>% arrange(-Log.P.value)
    enrichr.df$Format.Name %<>% factor(levels = enrichr.df$Format.Name)
    enrichr.df.plot <- select(enrichr.df, Format.Name, Log.Up, Log.Down) %>% gather(Direction, Length, -Format.Name) 

    p <- ggplot(enrichr.df.plot, aes(Format.Name, Length, fill = Direction)) 
    p <- p + geom_bar(stat = "identity", size = 1) 
    #p <- p + geom_text(label = c(as.character(enrichr.df$Format.Name), rep("", nrow(enrichr.df))), hjust = "left", aes(y = 0.1)) 
    p <- p + scale_fill_discrete(name = "Direction", labels = c("Up", "Down")) 
    p <- p + coord_flip() + theme_bw() + theme(axis.title.y = element_blank(), axis.ticks.y = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ylab(expression(paste(Log[10], ' P-Value')))
    p <- p + theme(plot.background = element_blank(), legend.background = element_blank(), axis.text.y = element_text(size = 12), panel.border = element_rect(color = "black", size = 1))
    CairoPDF(filename, height = plot.height, width = plot.width, bg = "transparent")
    print(p)
    dev.off()
}

UpDown <- function(filter.vector, enrichr.df) {
    split.vector <- str_split(filter.vector, ",")[[1]]
    enrichr.df$Symbol %<>% toupper
    enrichr.filter <- filter(enrichr.df, is.element(Symbol, split.vector))
    enrichr.vector <- c("Up" = length(which(sign(enrichr.filter$logFC) == 1)), "Down" = length(which(sign(enrichr.filter$logFC) == -1)))
    enrichr.vector
}

FilterEnrichr <- function(enrichr.df, size = 200) {
    enrichr.df$Num.Genes <- map(enrichr.df$Genes, str_split, ",") %>% map(extract2, 1) %>% map_int(length)
    enrichr.filter <- filter(enrichr.df, Num.Genes > 4) %>% filter(P.value < 0.05)
    if (nrow(enrichr.df) > size) {
        enrichr.filter %<>% slice(1:size)
    }
    enrichr.filter
}

ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
vega <- useMart("ENSEMBL_MART_VEGA", dataset = "mmusculus_gene_vega")

targets <- read_csv("./target.csv")
targets$group %<>% str_replace("\\(.*\\)", "") %>% factor
SaveRDSgz(targets, "./save/targets.rda")
raw.counts <- read.xlsx("../new_raw_data/2015-9201_metaReadCount_ensembl.xlsx") %>% data.frame
rownames(raw.counts) <- raw.counts$GeneID
counts.only <- select(raw.counts, -(GeneID:Strand)) %>% as.matrix 
counts.only <- counts.only[rowSums(counts.only) > 1, targets$sample]
counts.only <- counts.only[grepl("ENS", rownames(counts.only)),]
counts.mad <- rowMads(counts.only)
counts.only <- counts.only[counts.mad > 0,]
#counts.only <- counts.only[!duplicated(rownames(counts.annot)),]
SaveRDSgz(counts.only, "./save/counts.only.rda")

bm.table <- getBM(attributes = c('ensembl_gene_id', 'mgi_symbol'), filters = 'ensembl_gene_id', values = rownames(counts.only), mart = ensembl) 
colnames(bm.table) <- c("Ensembl.Gene.ID", "Symbol")
bm.filter <- filter(bm.table, nchar(Symbol) > 0 & !duplicated(Ensembl.Gene.ID))

counts.annot <- counts.only[bm.filter$Ensembl.Gene.ID,]
counts.collapse <- collapseRows(counts.annot, bm.filter$Symbol, rownames(counts.annot)) %>% extract2("datETcollapsed")
symbols.real <- !grepl("^Gm[0-9]|Rik$", rownames(counts.collapse))
counts.real <- counts.collapse[symbols.real,]

counts.dge <- DGEList(counts.real)
counts.dge$samples$group <- targets$group
counts.dgenorm <- calcNormFactors(counts.dge)
SaveRDSgz(counts.dgenorm, "./save/counts.dgenorm.rda")

targets$group2 <- str_replace_all(targets$group, "_grp.*$", "")
design.matrix <- model.matrix( ~ 0 + group2, targets)
colnames(design.matrix) %<>% str_replace('group2', "")
#contrasts.all <- makeContrasts(Sleep_dep - Sleep_grp1, Sleep_dep - Sleep_grp2, Sleep_dep - Wake_grp1, Sleep_dep - Wake_grp2, Sleep_grp2 - Sleep_grp1, Wake_grp1 - Sleep_grp1, Wake_grp2 - Sleep_grp1, Wake_grp1 - Sleep_grp2, Wake_grp2 - Sleep_grp2, Wake_grp2 - Wake_grp1, levels = design.matrix)
#format.names <- c("Sleep_dep_vs_Sleep_grp1", "Sleep_dep_vs_Sleep_grp2", "Sleep_dep_vs_Wake_grp1", "Sleep_dep_vs_Wake_grp2", "Sleep_grp2_vs_Sleep_grp1", "Wake_grp1_vs_Sleep_grp1", "Wake_grp2_vs_Sleep_grp1", "Wake_grp1_vs_Sleep_grp2", "Wake_grp2_vs_Sleep_grp2", "Wake_grp2_vs_Wake_grp1")

contrasts.all <- makeContrasts(Sleep_dep - Sleep, Sleep_dep - Wake, Wake - Sleep, levels = design.matrix)
format.names <- c("Sleep_dep_vs_Sleep", "Sleep_dep_vs_Wake", "Wake_vs_Sleep")

#limma voom
voom.raw <- voom(counts.only, design = design.matrix, normalize = "quantile")
voom.tmm <- voom(counts.dgenorm, design = design.matrix, normalize = "quantile")
SaveRDSgz(voom.tmm, "./save/voom.tmm.rda")
SaveRDSgz(voom.raw, "./save/voom.raw.rda")

#Annotate with BioMart (Ensembl and Vega)
#bm.ensembl <- getBM(attributes = c('ensembl_gene_id', 'mgi_symbol', 'description', 'gene_biotype'), filters = 'ensembl_gene_id', values = rownames(voom.tmm), mart = ensembl)
#bm.ensembl$description %<>% str_replace_all(" \\[.*\\]$", "")
#colnames(bm.ensembl) <- c("Ensembl.Gene.ID", "Symbol", "Description", "Gene.Type")
#SaveRDSgz(bm.ensembl, "./save/bm.ensembl.rda")

#bm.vega <- getBM(attributes = c('ens_gene', 'vega_gene_id', 'external_gene_name', 'description', 'gene_biotype'), filters = 'ens_gene', values = rownames(voom.tmm), mart = vega)
#bm.vega$description %<>% str_replace_all(" \\[.*\\]$", "")
#colnames(bm.vega) <- c("Ensembl.Gene.ID", "Vega.ID", "Gene.Name", "Description.Vega", "Gene.Type.Vega")
#SaveRDSgz(bm.vega, "./save/bm.vega.rda")

bm.ensembl <- getBM(attributes = c('mgi_symbol', 'description', 'gene_biotype'), filters = 'mgi_symbol', values = rownames(voom.tmm), mart = ensembl)
bm.ensembl$description %<>% str_replace_all(" \\[.*\\]$", "")
colnames(bm.ensembl) <- c("Symbol", "Description", "Gene.Type")
bm.filter <- filter(bm.ensembl, !duplicated(Symbol))

#combine.annot <- left_join(bm.ensembl, bm.vega) 
#SaveRDSgz(combine.annot, "./save/combine.annot.rda")

fit.voom <- lmFit(voom.tmm, design.matrix) %>% contrasts.fit(contrasts.all) %>% eBayes
top.table.all <- map2(1:3, format.names, FormatTopTable, fit.voom) %>% 
    reduce(left_join) %>%
    as_tibble %>%
    left_join(bm.filter) %>% 
    select(Symbol, Description, Gene.Type, dplyr::contains("logFC"), dplyr::contains("adj.P.Val"), dplyr::contains("P.Value"), AveExpr, dplyr::contains("t_"), dplyr::contains("B_"))
SaveRDSgz(top.table.all, "./save/top.table.all.rda")

GenWorkbook(top.table.all, "table.annot.xlsx", "P.Value", "adj.P.Val", "logFC")

ebam.list <- map(1:10, limma2ebam, fit = fit.voom)


objects.size <- lapply(ls(), function(thing) print(object.size(get(thing)), units = 'auto')) 
names(objects.size) <- ls()
unlist(objects.size) %>% sort

#Figures
mds.all <- t(voom.tmm$E) %>% dist(method = "manhattan") %>% cmdscale(eig = TRUE)
targets$Group <- targets$group2
MDSPlot("MDS.pdf", mds.all, targets, "Group")

thresholds <- c(0.01, 0.005, 0.001)
thresholds.fdr <- c(0.05, 0.01)

#Cutoffs
voom.pval <- str_c("P.Value", format.names, sep = "_")
voom.fdr <- str_c("adj.P.Val", format.names, sep = "_")
voom.logfc <- str_c("logFC", format.names, sep = "_")

thresholds <- c(0.001, 0.01)
thresholds.fdr <- c(0.05, 0.01)

toptable.pval <- map2_df(voom.pval, voom.logfc, GetThresholds, top.table.all, thresholds, "P value") 
toptable.fdr <- map2_df(voom.fdr, voom.logfc, GetThresholds, top.table.all, thresholds.fdr, "FDR")

#melt.voom <- rbind(melt.pval.voom, melt.fdr.voom) 
melt.voom <- gather(toptable.fdr, Direction, Num.Genes, -Threshold, -Comparison)

genetotals.voom <- group_by(melt.voom, Threshold) %>% summarise(Total.Genes = sum(abs(Num.Genes)))

toptable.plot <- left_join(melt.voom, genetotals.voom) #%>% filter(Threshold == "FDR < 0.05")
toptable.plot$Threshold %<>% factor
toptable.plot$Comparison %<>% str_replace_all("_", " ")

DecidePlot("toptable_thresholds", toptable.plot, "Differentially Expressed Genes", bar.padding = 250)

#Top Genes
test <- map(voom.pval, TopGenesPlot, top.table.all, voom.tmm, targets$group2)

#Heatmap
#targets$colors <- labels2colors(targets$group, colorSeq = brewer.pal(n = 5, "Set1"))
#cor.matrix <- bicor(voom.tmm$E, maxPOutliers = 0.05)
#CairoPDF("all.heatmap.pdf", width = 10, height = 10)
#heatmap.plus(cor.matrix, col = heat.colors(40), ColSideColors = as.matrix(targets$color), scale = "none")
#dev.off()

#Volcano
map(format.names, VolcanoPlot, top.table.all, 0.01)

#Upset
top.upset <- select(top.table.all, Symbol, dplyr::contains("adj.P.Val")) 
format.pvals <- str_c("adj.P.Val_", format.names, " < 0.05")
top.list <- map(format.pvals, filter_, .data = top.upset) %>% map(extract2, "Symbol") 
names(top.list) <- format.names

CairoPDF("upset.all", width = 30, height = 10)
upset(fromList(top.list), nsets = 10, nintersects = NA, order.by = "freq")
dev.off()

source("../../code/common_functions.R")

#Enrichr plots
source("../../code/GO/enrichr.R")
enrichr.terms <- c("GO_Biological_Process_2015", "GO_Molecular_Function_2015", "KEGG_2016", "Reactome_2016") 
trap <- map(format.names, EnrichrSubmit, top.table.all, "adj.P.Val", 0.01, enrichr.terms)

sds.symbols <- filter(top.table.all, adj.P.Val_Sleep_dep_vs_Sleep < 0.01) %>%
    arrange(desc(abs(logFC_Sleep_dep_vs_Sleep))) %>% 
    slice(1:1000) %>%
    extract2("Symbol")

sds.gobiol.file <- "./enrichr/Sleep_dep_vs_Sleep/GO_Biological_Process_2015.xlsx"
sds.gobiol <- read.xlsx(sds.gobiol.file) %>% as_tibble
sds.gobiol.filter <- FilterEnrichr(sds.gobiol)
GetKappaCluster(file_path_sans_ext(sds.gobiol.file), sds.gobiol.filter, sds.symbols)
sds.gobiol.final <- slice(sds.gobiol.filter, c(1, 4, 14))
sds.gobiol.final$Database <- "GO BP"

sds.gomole.file <- "./enrichr/Sleep_dep_vs_Sleep/GO_Molecular_Function_2015.xlsx"
sds.gomole <- read.xlsx(sds.gomole.file) %>% as_tibble
sds.gomole.filter <- FilterEnrichr(sds.gomole)
GetKappaCluster(file_path_sans_ext(sds.gomole.file), sds.gomole.filter, sds.symbols)
sds.gomole.final <- slice(sds.gomole.filter, c(9, 17))
sds.gomole.final$Database <- "GO MF"

sds.reactome.file <- "./enrichr/Sleep_dep_vs_Sleep/Reactome_2016.xlsx"
sds.reactome <- read.xlsx(sds.reactome.file) %>% as_tibble()
sds.reactome.filter <- FilterEnrichr(sds.reactome)
GetKappaCluster(file_path_sans_ext(sds.reactome.file), sds.reactome.filter, sds.symbols)
sds.reactome.final <- slice(sds.reactome.filter, c(1))
sds.reactome.final$Database <- "Reactome"

sds.kegg.file <- "./enrichr/Sleep_dep_vs_Sleep/KEGG_2016.xlsx"
sds.kegg <- read.xlsx(sds.kegg.file) %>% as_tibble()
sds.kegg.filter <- FilterEnrichr(sds.kegg)
GetKappaCluster(file_path_sans_ext(sds.kegg.file), sds.kegg.filter, sds.symbols)
sds.kegg.final <- slice(sds.kegg.filter, c(1))
sds.kegg.final$Database <- "KEGG"

sds.enrichr.final <- rbind(sds.gobiol.final, sds.gomole.final, sds.reactome.final, sds.kegg.final)
toptable.sds.enrichr <- select(top.table.all, Symbol, logFC_Sleep_dep_vs_Sleep)
colnames(toptable.sds.enrichr)[2] <- "logFC"
EnrichrPlot(sds.enrichr.final, toptable.sds.enrichr, "sds.enrichr")

sdw.symbols <- filter(top.table.all, adj.P.Val_Sleep_dep_vs_Wake < 0.01) %>%
    arrange(desc(abs(logFC_Sleep_dep_vs_Wake))) %>% 
    slice(1:1000) %>%
    extract2("Symbol")

sdw.gobiol.file <- "./enrichr/Sleep_dep_vs_Wake/GO_Biological_Process_2015.xlsx"
sdw.gobiol <- read.xlsx(sdw.gobiol.file) %>% as_tibble
sdw.gobiol.filter <- FilterEnrichr(sdw.gobiol)
GetKappaCluster(file_path_sans_ext(sdw.gobiol.file), sdw.gobiol.filter, sdw.symbols)
sdw.gobiol.final <- slice(sdw.gobiol.filter, c(1, 3, 6, 9))
sdw.gobiol.final$Database <- "GO BP"

sdw.gomole.file <- "./enrichr/Sleep_dep_vs_Wake/GO_Molecular_Function_2015.xlsx"
sdw.gomole <- read.xlsx(sdw.gomole.file) %>% as_tibble
sdw.gomole.filter <- FilterEnrichr(sdw.gomole)
GetKappaCluster(file_path_sans_ext(sdw.gomole.file), sdw.gomole.filter, sdw.symbols)
sdw.gomole.final <- slice(sdw.gomole.filter, c(1, 2))
sdw.gomole.final$Database <- "GO MF"

sdw.reactome.file <- "./enrichr/Sleep_dep_vs_Wake/Reactome_2016.xlsx"
sdw.reactome <- read.xlsx(sdw.reactome.file) %>% as_tibble
sdw.reactome.filter <- FilterEnrichr(sdw.reactome)
GetKappaCluster(file_path_sans_ext(sdw.reactome.file), sdw.reactome.filter, sdw.symbols)
sdw.reactome.final <- slice(sdw.reactome.filter, c(2, 3))
sdw.reactome.final$Database <- "Reactome"

sdw.kegg.file <- "./enrichr/Sleep_dep_vs_Wake/KEGG_2016.xlsx"
sdw.kegg <- read.xlsx(sdw.kegg.file) %>% as_tibble
sdw.kegg.filter <- FilterEnrichr(sdw.kegg)
GetKappaCluster(file_path_sans_ext(sdw.kegg.file), sdw.kegg.filter, sdw.symbols)
sdw.kegg.final <- slice(sdw.kegg.filter, c(2, 3))
sdw.kegg.final$Database <- "KEGG"

sdw.enrichr.final <- rbind(sdw.gobiol.final, sdw.gomole.final, sdw.reactome.final, sdw.kegg.final)
toptable.sdw.enrichr <- select(top.table.all, Symbol, logFC_Sleep_dep_vs_Sleep)
colnames(toptable.sdw.enrichr)[2] <- "logFC"
EnrichrPlot(sdw.enrichr.final, toptable.sdw.enrichr, "sdw.enrichr")

ws.symbols <- filter(top.table.all, adj.P.Val_Wake_vs_Sleep < 0.01) %>%
    arrange(desc(abs(logFC_Wake_vs_Sleep))) %>% 
    slice(1:1000) %>%
    extract2("Symbol")

ws.gobiol.file <- "./enrichr/Wake_vs_Sleep/GO_Biological_Process_2015.xlsx"
ws.gobiol <- read.xlsx(ws.gobiol.file) %>% as_tibble
ws.gobiol.filter <- FilterEnrichr(ws.gobiol)
GetKappaCluster(file_path_sans_ext(ws.gobiol.file), ws.gobiol.filter, ws.symbols)
ws.gobiol.final <- slice(ws.gobiol.filter, c(1, 4, 6, 13))
ws.gobiol.final$Database <- "GO BP"

ws.gomole.file <- "./enrichr/Wake_vs_Sleep/GO_Molecular_Function_2015.xlsx"
ws.gomole <- read.xlsx(ws.gomole.file) %>% as_tibble
ws.gomole.filter <- FilterEnrichr(ws.gomole)
GetKappaCluster(file_path_sans_ext(ws.gomole.file), ws.gomole.filter, ws.symbols)
ws.gomole.final <- slice(ws.gomole.filter, c(1, 6))
ws.gomole.final$Database <- "GO MF"

ws.reactome.file <- "./enrichr/Wake_vs_Sleep/Reactome_2016.xlsx"
ws.reactome <- read.xlsx(ws.reactome.file) %>% as_tibble
ws.reactome.filter <- FilterEnrichr(ws.reactome)
GetKappaCluster(file_path_sans_ext(ws.reactome.file), ws.reactome.filter, ws.symbols)
ws.reactome.final <- slice(ws.reactome.filter, c(5))
ws.reactome.final$Database <- "Reactome"

ws.kegg.file <- "./enrichr/Wake_vs_Sleep/KEGG_2016.xlsx"
ws.kegg <- read.xlsx(ws.kegg.file) %>% as_tibble
ws.kegg.filter <- FilterEnrichr(ws.kegg)
GetKappaCluster(file_path_sans_ext(ws.kegg.file), ws.kegg.filter, ws.symbols)
ws.kegg.final <- slice(ws.kegg.filter, c(1, 2))
ws.kegg.final$Database <- "KEGG"

ws.enrichr.final <- rbind(ws.gobiol.final, ws.gomole.final, ws.kegg.final, ws.reactome.final)
toptable.ws.enrichr <- select(top.table.all, Symbol, logFC_Sleep_dep_vs_Sleep)
colnames(toptable.ws.enrichr)[2] <- "logFC"
EnrichrPlot(ws.enrichr.final, toptable.ws.enrichr, "ws.enrichr", plot.width = 10)

