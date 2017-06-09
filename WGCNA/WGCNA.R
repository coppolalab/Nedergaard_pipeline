library(WGCNA)
library(flashClust)
enableWGCNAThreads()
library(BayesFactor)
library(PMCMR)
library(biomaRt)

library(Cairo)
library(igraph)
library(openxlsx)

library(R.utils)
library(matrixStats)
library(tools)

library(stringr)
library(broom)
library(magrittr)
library(tidyverse)

EigengeneANOVA <- function(ME.vector, trait.vector) {
    trait.df <- data.frame(Trait = trait.vector, ME = ME.vector)
    trait.anova <- kruskal.test(ME ~ Trait, trait.df) %>% tidy
    trait.anova$p.value
}

EigengeneAOV_old <- function(ME.vector, trait.vector, contrasts.vector) {
    trait.df <- data.frame(Trait = trait.vector, ME = ME.vector)
    trait.contrasts <- posthoc.kruskal.dunn.test(ME ~ Trait, trait.df, p.adjust.method = "none") 
    trait.pvals <- trait.contrasts$p.value

    pval.dep.sleep1 <- trait.pvals["Sleep_grp1","Sleep_dep"]
    pval.dep.sleep2 <- trait.pvals["Sleep_grp2","Sleep_dep"]
    pval.dep.wake1 <- trait.pvals["Wake_grp1","Sleep_dep"]
    pval.dep.wake2 <- trait.pvals["Wake_grp2","Sleep_dep"]

    pval.wake1.sleep1 <- trait.pvals["Wake_grp1","Sleep_grp1"]
    pval.wake2.sleep1 <- trait.pvals["Wake_grp2","Sleep_grp1"]
    pval.wake1.sleep2 <- trait.pvals["Wake_grp1","Sleep_grp2"]
    pval.wake2.sleep2 <- trait.pvals["Wake_grp2","Sleep_grp2"]

    pval.wake2.wake1 <- trait.pvals["Wake_grp2","Wake_grp1"]
    pval.sleep2.sleep1 <- trait.pvals["Sleep_grp2","Sleep_grp1"]

    pvals.subset <- p.adjust(c(pval.dep.sleep1, pval.dep.sleep2, pval.dep.wake1, pval.dep.wake2, pval.wake1.sleep1, pval.wake1.sleep2, pval.wake2.sleep1, pval.wake2.sleep2, pval.wake2.wake1, pval.sleep2.sleep1), method = "fdr")

    trait.medians <- group_by(trait.df, Trait) %>% summarise(median(ME)) %>% data.frame
    colnames(trait.medians)[2] <- "Median"
    rownames(trait.medians) <- trait.medians$Trait

    diff.dep.sleep1 <- trait.medians["Sleep_grp1","Median"] - trait.medians["Sleep_dep","Median"]
    diff.dep.sleep2 <- trait.medians["Sleep_grp2","Median"] - trait.medians["Sleep_dep","Median"]
    diff.dep.wake1 <- trait.medians["Wake_grp1","Median"] - trait.medians["Sleep_dep","Median"]
    diff.dep.wake2 <- trait.medians["Wake_grp2","Median"] - trait.medians["Sleep_dep","Median"]

    diff.wake1.sleep1 <- trait.medians["Wake_grp1","Median"] - trait.medians["Sleep_grp1","Median"]
    diff.wake1.sleep2 <- trait.medians["Wake_grp1","Median"] - trait.medians["Sleep_grp2","Median"]
    diff.wake2.sleep1 <- trait.medians["Wake_grp2","Median"] - trait.medians["Sleep_grp1","Median"]
    diff.wake2.sleep2 <- trait.medians["Wake_grp2","Median"] - trait.medians["Sleep_grp2","Median"]
    
    diff.sleep2.sleep1 <- trait.medians["Sleep_grp2","Median"] - trait.medians["Sleep_grp1","Median"]
    diff.wake2.wake1 <- trait.medians["Wake_grp2","Median"] - trait.medians["Wake_grp1","Median"]
    
    diffs.subset <- c(diff.dep.sleep1, diff.dep.sleep2, diff.dep.wake1, diff.dep.wake2, diff.wake1.sleep1, diff.wake1.sleep2, diff.wake2.sleep1, diff.wake2.sleep2, diff.wake2.wake1, diff.sleep2.sleep1)
    anova.return <- data.frame(Diff = diffs.subset, P.value = pvals.subset)
    rownames(anova.return) <- contrasts.vector

    anova.return
}

EigengeneAOV <- function(ME.vector, trait.vector, contrasts.vector) {
    trait.df <- data.frame(Trait = trait.vector, ME = ME.vector)
    trait.contrasts <- posthoc.kruskal.dunn.test(ME ~ Trait, trait.df, p.adjust.method = "none") 
    trait.pvals <- trait.contrasts$p.value

    pval.dep.sleep <- trait.pvals["Sleep_dep","Sleep"]
    pval.dep.wake <- trait.pvals["Wake","Sleep_dep"]
    pval.wake.sleep <- trait.pvals["Wake","Sleep"]

    pvals.subset <- p.adjust(c(pval.dep.sleep, pval.dep.wake, pval.wake.sleep), method = "fdr")

    trait.medians <- group_by(trait.df, Trait) %>% summarise(median(ME)) %>% data.frame
    colnames(trait.medians)[2] <- "Median"
    rownames(trait.medians) <- trait.medians$Trait

    diff.dep.sleep <- trait.medians["Sleep_dep","Median"] - trait.medians["Sleep","Median"] 
    diff.dep.wake <- trait.medians["Sleep_dep","Median"] - trait.medians["Wake","Median"]
    diff.wake.sleep <- trait.medians["Wake","Median"] - trait.medians["Sleep","Median"]
    
    diffs.subset <- c(diff.dep.sleep, diff.dep.wake, diff.wake.sleep)
    anova.return <- data.frame(Diff = diffs.subset, P.value = pvals.subset)
    rownames(anova.return) <- contrasts.vector

    anova.return
}

ModuleWorkbook <- function(module.table, filename) {
    pval.key <- colnames(module.table) %>% str_detect("pvalue") 
    pval.cols <- which(pval.key)
    MM.key <- colnames(module.table) %>% str_detect("MM.") 
    MM.cols <- which(MM.key & !pval.key)

    wb <- createWorkbook()
    addWorksheet(wb = wb, sheetName = "Sheet 1", gridLines = TRUE)
    writeDataTable(wb = wb, sheet = 1, x = module.table)
    sig.pvalues <- createStyle(fontColour = "red")
    conditionalFormatting(wb, 1, cols = pval.cols, rows = 1:nrow(module.table), rule = "<0.05", style = sig.pvalues)
    conditionalFormatting(wb, 1, cols = MM.cols, rows = 1:nrow(module.table), style = c("#63BE7B", "white", "red"), type = "colourScale")
    setColWidths(wb, 1, cols = 1, widths = 22)
    setColWidths(wb, 1, cols = 2, widths = 45)
    setColWidths(wb, 1, cols = 3:ncol(module.table), widths = "auto")
    pageSetup(wb, 1, orientation = "landscape", fitToWidth = TRUE)
    freezePane(wb, 1, firstRow = TRUE)
    showGridLines(wb, 1, showGridLines = TRUE)
    modifyBaseFont(wb, fontSize = 11)
    saveWorkbook(wb, filename, overwrite = TRUE) 
}

EnrichrSubmit <- function(index, full.df, enrichr.terms, use.weights = FALSE) {
    dataset <- filter(full.df, Module == index)
    dir.create(file.path("./enrichr", index), recursive = TRUE, showWarnings = FALSE)
    enrichr.data <- map(enrichr.terms, GetEnrichrData, dataset, FALSE)
    enrichr.names <- enrichr.terms[!is.na(enrichr.data)]
    enrichr.data <- enrichr.data[!is.na(enrichr.data)]
    names(enrichr.data) <- enrichr.names
    trap1 <- map(names(enrichr.data), EnrichrWorkbook, enrichr.data, index)
}

EnrichrWorkbook <- function(subindex, full.df, index) {
    dataset <- full.df[[subindex]]
    wb <- createWorkbook()
    addWorksheet(wb = wb, sheetName = "sheet 1", gridLines = TRUE)
    writeDataTable(wb = wb, sheet = 1, x = dataset, withFilter = TRUE)
    freezePane(wb, 1, firstRow = TRUE)
    setColWidths(wb, 0, cols = c(1, 3:ncol(dataset)), widths = "auto")
    setColWidths(wb, 1, cols = 2, widths = 45)

    filename = paste("./enrichr/", index, "/", index, "_", subindex, ".xlsx", sep = "")
    saveWorkbook(wb, filename, overwrite = TRUE) 
}

EnrichrPlot <- function(enrichr.df, filename, plot.title, plot.height = 5, plot.width = 8) {
    enrichr.df$Gene.Count <- map(enrichr.df$Genes, str_split, ",") %>% map(unlist) %>% map_int(length)
    enrichr.df$Adj.P.value <- p.adjust(enrichr.df$P.value, method = "fdr")
    enrichr.df$Log.P.value <- -log10(enrichr.df$Adj.P.value)
    enrichr.df$Term %<>% str_replace_all("\\ \\(.*$", "") %>% str_replace_all("\\_Homo.*$", "")  #Remove any thing after the left parenthesis and convert to all lower case
    enrichr.df$Format.Name <- str_c(enrichr.df$Database, ": ", enrichr.df$Term, " (", enrichr.df$Gene.Count, ")")
    enrichr.df %<>% arrange(Log.P.value)
    enrichr.df$Format.Name %<>% factor(levels = enrichr.df$Format.Name)

    p <- ggplot(enrichr.df, aes(Format.Name, Log.P.value)) + geom_bar(stat = "identity", size = 1) 
    p <- p + coord_flip() + theme_bw() + theme(axis.title.y = element_blank(), axis.ticks.y = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ylab(expression(paste(Log[10], ' P-Value')))
    p <- p + theme(plot.background = element_blank(), legend.background = element_blank(), axis.text.y = element_text(size = 12), panel.border = element_rect(color = "black", size = 1))
    CairoPDF(filename, height = plot.height, width = plot.width, bg = "transparent")
    print(p)
    dev.off()
}

FilterEnrichr <- function(enrichr.df, size = 200) {
    enrichr.df$Num.Genes <- map(enrichr.df$Genes, str_split, ",") %>% map(extract2, 1) %>% map_int(length)
    enrichr.filter <- filter(enrichr.df, Num.Genes > 4) %>% filter(P.value < 0.05)
    if (nrow(enrichr.df) > size) {
        enrichr.filter %<>% slice(1:size)
    }
    enrichr.filter
}

GetKscaled <- function(gene.list, module.membership) {
    filter(module.membership, is.element(Symbol, gene.list)) %>% select(Symbol, kscaled) %>% arrange(desc(kscaled))
}

Top5Plot <- function(rank.column, toptable.object, voom.object, pheno.object, pheno.col, levels.vector, maintitle, filename) {
    col.name <- str_c("desc(", rank.column, ")")
    rownames(voom.object) %<>% toupper

    top5.symbol <- arrange_(toptable.object, col.name)$Symbol[1:5] %>% toupper
    top5.expr <- t(voom.object[top5.symbol,])
    colnames(top5.expr) <- top5.symbol
    top5.df <- data.frame(Sample = pheno.object[[pheno.col]], top5.expr) %>%
        as_tibble %>%
        gather_("Gene", "Expression", names(.)[-1])
    top5.df$Gene %<>% factor(levels = make.names(top5.symbol))
    top5.df$Sample %<>% factor(levels = make.names(levels.vector))

    p <- ggplot(top5.df, aes(x = Sample, y = Expression, color = Sample)) + geom_boxplot() + geom_jitter() + theme_bw()
    p <- p + facet_wrap(~ Gene, ncol = 5, scales = "free") + theme(legend.position = "none")
    p <- p + theme(axis.title.x = element_blank())
    p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x = element_blank()) 
    p <- p + theme(plot.background = element_blank(), panel.border = element_rect(size = 1, color = "black"))
    p <- p + theme(plot.title = element_text(hjust = 0.5))
    p <- p + ggtitle(maintitle)
    CairoPDF(filename, height = 3.5, width = 16, bg = "transparent")
    print(p)
    dev.off()
}

GetPPI <- function(gene.list) {
  gene.list %<>% toupper
  biogrid.ppi.reduce <- select(biogrid.ppi, dplyr::contains("Official"))
  biogrid.ppi.reduce$Official.Symbol.Interactor.A %<>% toupper
  biogrid.ppi.reduce$Official.Symbol.Interactor.B %<>% toupper
  
  ppi.biogrid <- filter(biogrid.ppi.reduce, Official.Symbol.Interactor.A %in% gene.list) %>% filter(Official.Symbol.Interactor.B %in% gene.list)
  colnames(ppi.biogrid) <- c("Symbol.A", "Symbol.B")
  
  return(ppi.biogrid)
}

GetHyper <- function(overlap, total.ppi, list.genes, full.genelist) {
  total.genes <- nrow(full.genelist)
  background.genes <- total.genes - length(list.genes) - total.ppi + overlap
  
  hyper.pval <- phyper(overlap, total.ppi, background.genes, length(list.genes), lower.tail = FALSE) %>% signif(3)
  hyper.pval
}

test <- lapply(ls(), function(thing) print(object.size(get(thing)), units = 'auto')) 
names(test) <- ls()
unlist(test) %>% sort

source("../../code/common_functions.R")
voom.tmm <- ReadRDSgz(file = "../differential_expression/save/voom.tmm.rda")
voom.mad <- rowMads(voom.tmm$E)
voom.tmm.expr <- voom.tmm[voom.mad > 0,]

#Calculate scale free topology measures for different values of power adjacency function
powers <- c(c(1:10), seq(from = 12, to = 39, by = 2))

expr.collapse <- voom.tmm.expr$E %>% t
sft <- pickSoftThreshold(expr.collapse, powerVector = powers, verbose = 5, networkType = "signed", corFnc = bicor, corOptions = list(maxPOutliers = 0.05))
sft.df <- as_tibble(sft$fitIndices)
SaveRDSgz(sft, "./save/sft.rda")

#Plot scale indendence and mean connectivity as functions of power
sft.df$multiplied <- sft.df$SFT.R.sq * -sign(sft.df$slope)
p <- ggplot(sft.df, aes(x = Power,  y = multiplied, label = Power)) + geom_point() + geom_text(vjust = -0.6, size = 4, col = "red")
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + geom_hline(aes(yintercept = 0.9))
p <- p + xlab("Soft Threshold") + ylab("Scale Free Topology Model Fit, signed R^2") + ggtitle("Scale Independence")
CairoPDF(file = "./scaleindependence", width = 6, height = 6)
print(p)
dev.off()

p <- ggplot(sft.df, aes(x = Power,  y = mean.k.)) + geom_point() 
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + xlab("Soft Threshold") + ylab("Mean Connectivity") + ggtitle("Mean Connectivity")
CairoPDF(file = "./meanconnectivity", width = 6, height = 6)
print(p)
dev.off()

softPower <- 6
adjacency.expr <- adjacency(expr.collapse, power = softPower, type = "signed", corFnc = "bicor", corOptions = "maxPOutliers = 0.05")

TOM <- TOMsimilarity(adjacency.expr, verbose = 5)
dissimilarity.TOM <- 1 - TOM
rm(TOM)
gc()

geneTree = flashClust(as.dist(dissimilarity.TOM), method = "average")
SaveRDSgz(geneTree, file = "./save/gene.tree.rda")

CairoPDF(file = "./genecluster", height = 10, width = 15)
plot(geneTree, xlab = "", sub = "", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04)
dev.off()

#Identify modules using dynamic tree cutting with hybrid clustering
min.module.size <- 20
dynamic.modules <- cutreeDynamic(dendro = geneTree, cutHeight = 0.995, method = "hybrid", distM = dissimilarity.TOM, pamRespectsDendro = FALSE, minClusterSize = min.module.size, verbose = 2)
dynamic.colors <- labels2colors(dynamic.modules)
SaveRDSgz(dynamic.colors, file = "./save/dynamic.colors.rda")

CairoPDF(file = "./gene_dendrogram_and_module_colors", height = 10, width = 15)
plotDendroAndColors(geneTree, dynamic.colors, "", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "")
dev.off()

#Calculate module eigengenes
ME.list <- moduleEigengenes(expr.collapse, colors = dynamic.colors, softPower = softPower, nPC = 1)
ME.genes <- ME.list$eigengenes
MEDiss <- 1 - bicor(ME.genes, maxPOutliers = 0.05)
METree <- flashClust(as.dist(MEDiss), method = "average")
SaveRDSgz(METree, file = "./save/me.tree.rda")

CairoPDF(file = "./module_eigengene_clustering_tree", height = 10, width = 15)
plot(METree, xlab = "", sub = "", main = "")
dev.off()

#Check if any modules are too similar and merge them. 
ME.dissimilarity.threshold <- 0.3
merge.all <- mergeCloseModules(expr.collapse, dynamic.colors, cutHeight = ME.dissimilarity.threshold, verbose = 3, corFnc = bicor, corOptions = list(maxPOutliers = 0.05)) 
merged.colors <- merge.all$colors
merged.genes <- merge.all$newMEs

CairoPDF("module_eigengene_clustering", height = 8, width = 10, bg = "transparent")
plotDendroAndColors(geneTree, cbind(dynamic.modules, merged.colors), c("Original Modules", "Merged Modules"), dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "", marAll = c(1,8,3,1))
dev.off()

#CairoPDF("module_eigengene_clustering_new", height = 8, width = 10, bg = "transparent")
#plot(geneTree)
#dev.off()

#Use merged eigengenes 
module.colors <- merged.colors
SaveRDSgz(module.colors, file = "./save/module.colors.rda")
color.order <- c("grey", standardColors(50))
modules.labels <- match(module.colors, color.order)
SaveRDSgz(modules.labels, file = "./save/modules.labels.rda")
ME.genes <- as_tibble(merged.genes) %>% select_(., .dots = sort(colnames(.)))
SaveRDSgz(ME.genes, file = "./save/me.genes.rda")

CairoPDF("eigengenes", height = 6, width = 8, bg = "transparent")
plotEigengeneNetworks(ME.genes, "", marDendro = c(0,4,1,2), marHeatmap = c(3,5,1,2), plotPreservation = "standard")
dev.off()

#Create heatmap of eigengene significance
format.names <- c("Sleep_dep_vs_Sleep", "Sleep_dep_vs_Wake", "Wake_vs_Sleep")
color.values <- unique(module.colors)
pdata <- ReadRDSgz("../differential_expression/save/targets.rda")
anova.status <- map_dbl(ME.genes, EigengeneANOVA, pdata$group2) %>% p.adjust("fdr") %>% signif(3)

aov.status <- map(ME.genes, EigengeneAOV, pdata$group2, format.names)
status.diff <- map(aov.status, select, Diff) %>% map(t) %>% reduce(rbind) %>% data.frame
rownames(status.diff) <- names(aov.status)
status.pval <- map(aov.status, select, P.value) %>% map(t) %>% reduce(rbind) %>% data.frame
rownames(status.pval) <- names(aov.status)

text.matrix.traits <- str_c(signif(as.matrix(status.diff), 2), '\n(', signif(as.matrix(status.pval), 1), ')')
dim(text.matrix.traits) = dim(status.diff)

heatmap.range <- c(min(as.matrix(status.diff)) * 1.1, max(as.matrix(status.diff)) * 1.1)
width.dynamic <- 3 + (1 * ncol(text.matrix.traits))

CairoPDF("module_trait_relationships", width = width.dynamic, height = 10, bg = "transparent")
par(mar = c(8, 11, 3, 3))
labeledHeatmap(Matrix = as.matrix(status.diff), xLabels = colnames(status.diff), yLabels = colnames(ME.genes), ySymbols = colnames(ME.genes), yColorLabels = TRUE, colors = greenWhiteRed(50), textMatrix = text.matrix.traits, setStdMargins = F, zlim = heatmap.range, main = "")
dev.off()

group.levels <- c("Sleep", "Wake", "Sleep_dep")
ME.genes.plot <- mutate(ME.genes, Group = factor(pdata$group2, levels = group.levels)) %>% gather(Module.Color, Eigengene, -Group) 
ME.genes.plot$Module.Color %<>% str_replace("ME", "")

p <- ggplot(ME.genes.plot, aes(x = Group, y = Eigengene, color = Group)) + geom_boxplot(width = 0.5) + geom_jitter(width = 0.5)
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
p <- p + scale_fill_manual(values = sort(unique(ME.genes.plot$Module.Color)))
p <- p + facet_wrap(~ Module.Color, ncol = 4, scales = "free") + theme(plot.background = element_blank())
CairoPDF("eigengene_plots", height = 10, width = 15, bg = "transparent")
print(p)
dev.off()

#Generate network statistics
all.degrees <- intramodularConnectivity(adjacency.expr, module.colors)
gene.info <- mutate(all.degrees, Symbol = rownames(all.degrees), Module = module.colors) %>% as_tibble %>% arrange(Module)
kWithin_max <- group_by(gene.info, Module) %>% summarise(kWithin_Max = max(kWithin)) 
gene.info.final <- left_join(gene.info, kWithin_max) %>% 
    mutate(kscaled = kWithin / kWithin_Max) %>% 
    select(Symbol, Module, kTotal:kDiff, kscaled, -kWithin_Max)

gene.module.membership <- bicor(expr.collapse, ME.genes, maxPOutliers = 0.05) %>% 
    as_tibble %>% 
    set_colnames(str_replace(names(ME.genes),"ME", "MM.")) %>%
    mutate(Symbol = colnames(expr.collapse)) 
module.membership.pvalue <- mutate_all(select(gene.module.membership, -Symbol), corPvalueStudent, nSamples = nrow(expr.collapse)) %>% 
    as_tibble %>% 
    set_colnames(str_replace(names(ME.genes),"ME", "MM.pvalue.")) %>%
    mutate(Symbol = colnames(expr.collapse))

ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

#Annotate with BioMart (Ensembl and Vega)
bm.ensembl <- getBM(attributes = c('mgi_symbol', 'description', 'gene_biotype'), filters = 'mgi_symbol', values = gene.info.final$Symbol, mart = ensembl)
bm.ensembl$description %<>% str_replace_all(" \\[.*\\]$", "")
colnames(bm.ensembl) <- c("Symbol", "Description", "Gene.Type")
SaveRDSgz(bm.ensembl, "./save/bm.ensembl.rda")

#bm.vega <- getBM(attributes = c('ens_gene', 'vega_gene_id', 'external_gene_name', 'description', 'gene_biotype'), filters = 'ens_gene', values = gene.info$Ensembl.ID, mart = vega)
#bm.vega$description %<>% str_replace_all(" \\[.*\\]$", "")
#colnames(bm.vega) <- c("Ensembl.ID", "Vega.ID", "Gene.Name", "Description.Vega", "Gene.Type.Vega")
#SaveRDSgz(bm.vega, "./save/bm.vega.rda")

#combine.annot <- left_join(bm.ensembl, bm.vega) 
#SaveRDSgz(combine.annot, "./save/combine.annot.rda")

color.names <- unique(module.colors) %>% sort
mm.select <- str_c("MM.", c(first(color.names), last(color.names))) %>% str_c(collapse = ":")
mm.pvalue.select <- str_c("MM.pvalue.", c(first(color.names), last(color.names))) %>% str_c(collapse = ":")
color.select <- c(mm.select, mm.pvalue.select)

module.membership <- left_join(gene.info.final, bm.ensembl) %>% 
    left_join(gene.module.membership) %>% 
    left_join(module.membership.pvalue) %>%
    select_("Symbol", "Description", "Gene.Type", "Module:kscaled", .dots = color.select) %>%
    arrange(Module, desc(kscaled)) %>%
    mutate_if(is.numeric, signif, digits = 3)
ModuleWorkbook(module.membership, "module_membership.xlsx")

#Enrichr
source("../../code/GO/enrichr.R")

modules.filter <- filter(module.membership, Gene.Type == "protein_coding") %>% select(Symbol, Module)
enrichr.terms <- c("GO_Biological_Process_2015", "GO_Molecular_Function_2015", "KEGG_2016", "Reactome_2016") 
trap1 <- map(color.names, EnrichrSubmit, modules.filter, enrichr.terms)

#Blue
blue.only <- filter(modules.filter, Module == "blue") %>% extract2("Symbol")
blue.gobiol.file <- "./enrichr/blue/blue_GO_Biological_Process_2015.xlsx"
blue.gobiol <- read.xlsx(blue.gobiol.file) %>% as_tibble
blue.gobiol.filter <- FilterEnrichr(blue.gobiol)
GetKappaCluster(file_path_sans_ext(blue.gobiol.file), blue.gobiol.filter, blue.only)
blue.gobiol.final <- slice(blue.gobiol, c(1, 2, 5))
blue.gobiol.final$Database <- "GO BP"

blue.gomole.file <- "./enrichr/blue/blue_GO_Molecular_Function_2015.xlsx"
blue.gomole <- read.xlsx(blue.gomole.file) %>% as_tibble
blue.gomole.filter <- FilterEnrichr(blue.gomole)
GetKappaCluster(file_path_sans_ext(blue.gomole.file), blue.gomole.filter, blue.only)
blue.gomole.final <- slice(blue.gomole, c(16, 18))
blue.gomole.final$Database <- "GO MF"

blue.reactome.file <- "./enrichr/blue/blue_Reactome_2016.xlsx"
blue.reactome <- read.xlsx(blue.reactome.file) %>% as_tibble
blue.reactome.filter <- FilterEnrichr(blue.reactome)
GetKappaCluster(file_path_sans_ext(blue.reactome.file), blue.reactome.filter, blue.only)
blue.reactome.final <- slice(blue.reactome, c(2, 3, 5))
blue.reactome.final$Database <- "Reactome"

blue.kegg.file <- "./enrichr/blue/blue_KEGG_2016.xlsx"
blue.kegg <- read.xlsx(blue.kegg.file) %>% as_tibble
blue.kegg.filter <- FilterEnrichr(blue.kegg)
GetKappaCluster(file_path_sans_ext(blue.kegg.file), blue.kegg.filter, blue.only)
blue.kegg.final <- slice(blue.kegg, c(1, 5))
blue.kegg.final$Database <- "KEGG"

blue.enrichr.final <- rbind(blue.gobiol.final, blue.gomole.final, blue.reactome.final, blue.kegg.final)
EnrichrPlot(blue.enrichr.final, "blue.enrichr")

#Black
black.only <- filter(modules.filter, Module == "black") %>% extract2("Symbol")
black.gobiol.file <- "./enrichr/black/black_GO_Biological_Process_2015.xlsx"
black.gobiol <- read.xlsx(black.gobiol.file) %>% as_tibble
black.gobiol.filter <- FilterEnrichr(black.gobiol)
GetKappaCluster(file_path_sans_ext(black.gobiol.file), black.gobiol.filter, black.only)
black.gobiol.final <- slice(black.gobiol, c(4, 7))
black.gobiol.final$Database <- "GO BP"

black.gomole.file <- "./enrichr/black/black_GO_Molecular_Function_2015.xlsx"
black.gomole <- read.xlsx(black.gomole.file) %>% as_tibble
black.gomole.filter <- FilterEnrichr(black.gomole)
GetKappaCluster(file_path_sans_ext(black.gomole.file), black.gomole.filter, black.only)
black.gomole.final <- slice(black.gomole, c(1))
black.gomole.final$Database <- "GO MF"

black.reactome.file <- "./enrichr/black/black_Reactome_2016.xlsx"
black.reactome <- read.xlsx(black.reactome.file) %>% as_tibble
black.reactome.filter <- FilterEnrichr(black.reactome)
GetKappaCluster(file_path_sans_ext(black.reactome.file), black.reactome.filter, black.only)
black.reactome.final <- slice(black.reactome, c(2, 6))
black.reactome.final$Database <- "Reactome"

black.kegg.file <- "./enrichr/black/black_KEGG_2016.xlsx"
black.kegg <- read.xlsx(black.kegg.file) %>% as_tibble
black.kegg.filter <- FilterEnrichr(black.kegg)
GetKappaCluster(file_path_sans_ext(black.kegg.file), black.kegg.filter, black.only)
black.kegg.final <- slice(black.kegg, c(2))
black.kegg.final$Database <- "KEGG"

black.enrichr.final <- rbind(black.gobiol.final, black.gomole.final, black.reactome.final, black.kegg.final)
EnrichrPlot(black.enrichr.final, "black.enrichr")

#Purple - ENRICHR rejected!
purple.only <- filter(modules.filter, Module == "purple") %>% extract2("Symbol")
purple.gobiol.file <- "./enrichr/purple/purple_GO_Biological_Process_2015.xlsx"
purple.gobiol <- read.xlsx(purple.gobiol.file)  %>% as_tibble
purple.gobiol.filter <- FilterEnrichr(purple.gobiol)
GetKappaCluster(file_path_sans_ext(purple.gobiol.file), purple.gobiol.filter, purple.only)
purple.gobiol.final <- slice(purple.gobiol, c(1, 4, 15, 13))

purple.gomole.file <- "./enrichr/purple/purple_GO_Molecular_Function_2015.xlsx"
purple.gomole <- read.xlsx(purple.gomole.file) %>% as_tibble
purple.gomole.filter <- FilterEnrichr(purple.gomole)
GetKappaCluster(file_path_sans_ext(purple.gomole.file), purple.gomole.filter, purple.only)
purple.gomole.final <- slice(purple.gomole, c(1))

purple.reactome.file <- "./enrichr/purple/purple_Reactome_2016.xlsx"
purple.reactome <- read.xlsx(purple.reactome.file) %>% as_tibble
purple.reactome.filter <- FilterEnrichr(purple.reactome)
GetKappaCluster(file_path_sans_ext(purple.reactome.file), purple.reactome.filter, purple.only)
purple.reactome.final <- slice(purple.reactome, c(35, 54))

purple.kegg.file <- "./enrichr/purple/purple_KEGG_2016.xlsx"
purple.kegg <- read.xlsx(purple.kegg.file) %>% as_tibble
purple.kegg.filter <- FilterEnrichr(purple.kegg)
GetKappaCluster(file_path_sans_ext(purple.kegg.file), purple.kegg.filter, purple.only)
purple.kegg.final <- slice(purple.kegg, c(2, 4))

purple.enrichr.final <- rbind(purple.gobiol.final, purple.gomole.final, purple.reactome.final, purple.kegg.final)
EnrichrPlot(purple.enrichr.final, "purple.enrichr")

#Magenta
red.only <- filter(modules.filter, Module == "red") %>% extract2("Symbol")
red.gobiol.file <- "./enrichr/red/red_GO_Biological_Process_2015.xlsx"
red.gobiol <- read.xlsx(red.gobiol.file) %>% as_tibble
red.gobiol.filter <- FilterEnrichr(red.gobiol)
GetKappaCluster(file_path_sans_ext(red.gobiol.file), red.gobiol.filter, red.only)
red.gobiol.final <- slice(red.gobiol, c(1, 4, 9))
red.gobiol.final$Database <- "GO BP"

red.gomole.file <- "./enrichr/red/red_GO_Molecular_Function_2015.xlsx"
red.gomole <- read.xlsx(red.gomole.file) %>% as_tibble
red.gomole.filter <- FilterEnrichr(red.gomole)
GetKappaCluster(file_path_sans_ext(red.gomole.file), red.gomole.filter, red.only)
red.gomole.final <- slice(red.gomole, c(2))
red.gomole.final$Database <- "GO MF"

red.reactome.file <- "./enrichr/red/red_Reactome_2016.xlsx"
red.reactome <- read.xlsx(red.reactome.file) %>% as_tibble
red.reactome.filter <- FilterEnrichr(red.reactome)
GetKappaCluster(file_path_sans_ext(red.reactome.file), red.reactome.filter, red.only)
red.reactome.final <- slice(red.reactome, c(1))
red.reactome.final$Database <- "Reactome"

red.kegg.file <- "./enrichr/red/red_KEGG_2016.xlsx"
red.kegg <- read.xlsx(red.kegg.file) %>% as_tibble
red.kegg.filter <- FilterEnrichr(red.kegg)
GetKappaCluster(file_path_sans_ext(red.kegg.file), red.kegg.filter, red.only)
red.kegg.final <- slice(red.kegg, c(2))
red.kegg.final$Database <- "KEGG"

red.enrichr.final <- rbind(red.gobiol.final, red.gomole.final, red.reactome.final, red.kegg.final)
EnrichrPlot(red.enrichr.final, "red.enrichr")

#Pink
pink.only <- filter(modules.filter, Module == "pink") %>% extract2("Symbol")
pink.gobiol.file <- "./enrichr/pink/pink_GO_Biological_Process_2015.xlsx"
pink.gobiol <- read.xlsx(pink.gobiol.file) %>% as_tibble
pink.gobiol.filter <- FilterEnrichr(pink.gobiol)
GetKappaCluster(file_path_sans_ext(pink.gobiol.file), pink.gobiol.filter, pink.only)
pink.gobiol.final <- slice(pink.gobiol, c(1, 3, 6))
pink.gobiol.final$Database <- "GO BP"

pink.gomole.file <- "./enrichr/pink/pink_GO_Molecular_Function_2015.xlsx"
pink.gomole <- read.xlsx(pink.gomole.file) %>% as_tibble
pink.gomole.filter <- FilterEnrichr(pink.gomole)
GetKappaCluster(file_path_sans_ext(pink.gomole.file), pink.gomole.filter, pink.only)
pink.gomole.final <- slice(pink.gomole, c(6))
pink.gomole.final$Database <- "GO MF"

pink.reactome.file <- "./enrichr/pink/pink_Reactome_2016.xlsx"
pink.reactome <- read.xlsx(pink.reactome.file) %>% as_tibble
pink.reactome.filter <- FilterEnrichr(pink.reactome)
GetKappaCluster(file_path_sans_ext(pink.reactome.file), pink.reactome.filter, pink.only)
pink.reactome.final <- slice(pink.reactome, c(3, 4))
pink.reactome.final$Database <- "Reactome"

pink.kegg.file <- "./enrichr/pink/pink_KEGG_2016.xlsx"
pink.kegg <- read.xlsx(pink.kegg.file) %>% as_tibble
pink.kegg.filter <- FilterEnrichr(pink.kegg)
GetKappaCluster(file_path_sans_ext(pink.kegg.file), pink.kegg.filter, pink.only)
pink.kegg.final <- slice(pink.kegg, c(1, 3))
pink.kegg.final$Database <- "KEGG"

pink.enrichr.final <- rbind(pink.gobiol.final, pink.gomole.final, pink.reactome.final, pink.kegg.final)
EnrichrPlot(pink.enrichr.final, "pink.enrichr")

#Cyan
cyan.only <- filter(modules.filter, Module == "cyan") %>% extract2("Symbol")
cyan.gobiol.file <- "./enrichr/cyan/cyan_GO_Biological_Process_2015.xlsx"
cyan.gobiol <- read.xlsx(cyan.gobiol.file) %>% as_tibble
cyan.gobiol.filter <- FilterEnrichr(cyan.gobiol)
GetKappaCluster(file_path_sans_ext(cyan.gobiol.file), cyan.gobiol.filter, cyan.only)
cyan.gobiol.final <- slice(cyan.gobiol, c(1, 6))
cyan.gobiol.final$Database <- "GO BP"

cyan.gomole.file <- "./enrichr/cyan/cyan_GO_Molecular_Function_2015.xlsx"
cyan.gomole <- read.xlsx(cyan.gomole.file) %>% as_tibble
cyan.gomole.filter <- FilterEnrichr(cyan.gomole)
GetKappaCluster(file_path_sans_ext(cyan.gomole.file), cyan.gomole.filter, cyan.only)
cyan.gomole.final <- slice(cyan.gomole, c(2, 10))
cyan.gomole.final$Database <- "GO MF"

cyan.reactome.file <- "./enrichr/cyan/cyan_Reactome_2016.xlsx"
cyan.reactome <- read.xlsx(cyan.reactome.file) %>% as_tibble
cyan.reactome.filter <- FilterEnrichr(cyan.reactome)
GetKappaCluster(file_path_sans_ext(cyan.reactome.file), cyan.reactome.filter, cyan.only)
cyan.reactome.final <- slice(cyan.reactome, c(2))
cyan.reactome.final$Database <- "Reactome"

cyan.kegg.file <- "./enrichr/cyan/cyan_KEGG_2016.xlsx"
cyan.kegg <- read.xlsx(cyan.kegg.file) %>% as_tibble
cyan.kegg.filter <- FilterEnrichr(cyan.kegg)
GetKappaCluster(file_path_sans_ext(cyan.kegg.file), cyan.kegg.filter, cyan.only)
cyan.kegg.final <- slice(cyan.kegg, c(8))
cyan.kegg.final$Database <- "KEGG"

cyan.enrichr.final <- rbind(cyan.gobiol.final, cyan.gomole.final, cyan.reactome.final, cyan.kegg.final)
EnrichrPlot(cyan.enrichr.final, "cyan.enrichr")

#Lightcyan
lightcyan.only <- filter(modules.filter, Module == "lightcyan") %>% extract2("Symbol")
lightcyan.gobiol.file <- "./enrichr/lightcyan/lightcyan_GO_Biological_Process_2015.xlsx"
lightcyan.gobiol <- read.xlsx(lightcyan.gobiol.file) %>% as_tibble
lightcyan.gobiol.filter <- FilterEnrichr(lightcyan.gobiol)
GetKappaCluster(file_path_sans_ext(lightcyan.gobiol.file), lightcyan.gobiol.filter, lightcyan.only)
lightcyan.gobiol.final <- slice(lightcyan.gobiol, c(14, 6, 12))
lightcyan.gobiol.final$Database <- "GO BP"

lightcyan.gomole.file <- "./enrichr/lightcyan/lightcyan_GO_Molecular_Function_2015.xlsx"
lightcyan.gomole <- read.xlsx(lightcyan.gomole.file) %>% as_tibble
lightcyan.gomole.filter <- FilterEnrichr(lightcyan.gomole)
GetKappaCluster(file_path_sans_ext(lightcyan.gomole.file), lightcyan.gomole.filter, lightcyan.only)
lightcyan.gomole.final <- slice(lightcyan.gomole, c(2))
lightcyan.gomole.final$Database <- "GO MF"

lightcyan.reactome.file <- "./enrichr/lightcyan/lightcyan_Reactome_2016.xlsx"
lightcyan.reactome <- read.xlsx(lightcyan.reactome.file) %>% as_tibble
lightcyan.reactome.filter <- FilterEnrichr(lightcyan.reactome)
GetKappaCluster(file_path_sans_ext(lightcyan.reactome.file), lightcyan.reactome.filter, lightcyan.only)
lightcyan.reactome.final <- slice(lightcyan.reactome, c(1, 4))
lightcyan.reactome.final$Database <- "Reactome"

lightcyan.kegg.file <- "./enrichr/lightcyan/lightcyan_KEGG_2016.xlsx"
lightcyan.kegg <- read.xlsx(lightcyan.kegg.file) %>% as_tibble
lightcyan.kegg.filter <- FilterEnrichr(lightcyan.kegg)
GetKappaCluster(file_path_sans_ext(lightcyan.kegg.file), lightcyan.kegg.filter, lightcyan.only)
lightcyan.kegg.final <- slice(lightcyan.kegg, c(2))
lightcyan.kegg.final$Database <- "KEGG"

lightcyan.enrichr.final <- rbind(lightcyan.gobiol.final, lightcyan.gomole.final, lightcyan.reactome.final, lightcyan.kegg.final)
EnrichrPlot(lightcyan.enrichr.final, "lightcyan.enrichr")

#Greenyellow
greenyellow.only <- filter(modules.filter, Module == "greenyellow") %>% extract2("Symbol")
greenyellow.gobiol.file <- "./enrichr/greenyellow/greenyellow_GO_Biological_Process_2015.xlsx"
greenyellow.gobiol <- read.xlsx(greenyellow.gobiol.file) %>% as_tibble
greenyellow.gobiol.filter <- FilterEnrichr(greenyellow.gobiol)
GetKappaCluster(file_path_sans_ext(greenyellow.gobiol.file), greenyellow.gobiol.filter, greenyellow.only)
greenyellow.gobiol.final <- slice(greenyellow.gobiol, c(1, 12, 13))
greenyellow.gobiol.final$Database <- "GO BP"

greenyellow.gomole.file <- "./enrichr/greenyellow/greenyellow_GO_Molecular_Function_2015.xlsx"
greenyellow.gomole <- read.xlsx(greenyellow.gomole.file) %>% as_tibble
greenyellow.gomole.filter <- FilterEnrichr(greenyellow.gomole)
GetKappaCluster(file_path_sans_ext(greenyellow.gomole.file), greenyellow.gomole.filter, greenyellow.only)
greenyellow.gomole.final <- slice(greenyellow.gomole, c(1))
greenyellow.gomole.final$Database <- "GO MF"

#greenyellow.reactome.file <- "./enrichr/greenyellow/greenyellow_Reactome_2016.xlsx"
#greenyellow.reactome <- read.xlsx(greenyellow.reactome.file) %>% as_tibble
#greenyellow.reactome.filter <- FilterEnrichr(greenyellow.reactome)
#GetKappaCluster(file_path_sans_ext(greenyellow.reactome.file), greenyellow.reactome.filter, greenyellow.only)
#greenyellow.reactome.final <- slice(greenyellow.reactome, c(8, 3))

#greenyellow.kegg.file <- "./enrichr/greenyellow/greenyellow_KEGG_2016.xlsx"
#greenyellow.kegg <- read.xlsx(greenyellow.kegg.file) %>% as_tibble
#greenyellow.kegg.filter <- FilterEnrichr(greenyellow.kegg)
#GetKappaCluster(file_path_sans_ext(greenyellow.kegg.file), greenyellow.kegg.filter, greenyellow.only)
#greenyellow.kegg.final <- slice(greenyellow.kegg, c(8))

greenyellow.enrichr.final <- rbind(greenyellow.gobiol.final, greenyellow.gomole.final)
EnrichrPlot(greenyellow.enrichr.final, "greenyellow.enrichr")

#Can't use original ppi because it was written for humans
biogrid.ppi <- read_tsv("../../code/BIOGRID-ORGANISM-3.4.136.tab2/BIOGRID-ORGANISM-Mus_musculus-3.4.136.tab2.txt") %>%
    set_colnames(make.names(colnames(.)))
biogrid.ppi$Official.Symbol.Interactor.A %<>% toupper
biogrid.ppi$Official.Symbol.Interactor.B %<>% toupper
biogrid.ppi.reduce <- select(biogrid.ppi, dplyr::contains("Official")) %>%
    set_colnames(c("Symbol.A", "Symbol.B"))

biogrid.graph <- graph_from_edgelist(as.matrix(biogrid.ppi.reduce), directed = FALSE) %>% igraph::simplify()
ppi.incident <- map(V(biogrid.graph), incident, graph = biogrid.graph) %>% map_int(length) %>% sort
ppi.df <- tibble(Symbol = names(ppi.incident), Total = ppi.incident) 

#Top table
top.table <- ReadRDSgz("../differential_expression/save/top.table.all.rda") %>%
    select(Symbol, contains('logFC'), contains('adj.P.Val'))
top.table.upper <- mutate(top.table, Symbol = toupper(Symbol))
module.reduce <- filter(module.membership, Gene.Type == "protein_coding") %>% select(Symbol, Description, Module, kscaled)

red.only <- filter(modules.filter, Module == "red") %>% extract2("Symbol") %>% toupper
pink.only <- filter(modules.filter, Module == "pink") %>% extract2("Symbol") %>% toupper
purple.only <- filter(modules.filter, Module == "purple") %>% extract2("Symbol") %>% toupper
black.only <- filter(modules.filter, Module == "black") %>% extract2("Symbol") %>% toupper
blue.only <- filter(modules.filter, Module == "blue") %>% extract2("Symbol") %>% toupper
cyan.only <- filter(modules.filter, Module == "cyan") %>% extract2("Symbol") %>% toupper
lightcyan.only <- filter(modules.filter, Module == "lightcyan") %>% extract2("Symbol") %>% toupper
greenyellow.only <- filter(modules.filter, Module == "greenyellow") %>% extract2("Symbol") %>% toupper

#PPI and Hub
red.ppi <-  GetPPI(red.only)
pink.ppi <-  GetPPI(pink.only)
purple.ppi <-  GetPPI(purple.only)
black.ppi <-  GetPPI(black.only)
blue.ppi <- GetPPI(blue.only)
cyan.ppi <- GetPPI(cyan.only)
lightcyan.ppi <- GetPPI(lightcyan.only)
greenyellow.ppi <- GetPPI(greenyellow.only)

group.levels <- c("Sleep", "Wake", "Sleep_dep")
#Black
black.top <- filter(module.reduce, Module == "black") %>% 
    left_join(top.table) %>%
    filter(adj.P.Val_Sleep_dep_vs_Sleep < 0.05 & adj.P.Val_Wake_vs_Sleep < 0.05)
Top5Plot("kscaled", black.top, voom.tmm.expr$E, pdata, 'group2', group.levels, "Hub Genes", "black.top5")

black.all.graph <- graph_from_edgelist(as.matrix(black.ppi), directed = FALSE) %>% igraph::simplify()
black.incident <- map(V(black.all.graph), incident, graph = black.all.graph) %>% 
    map_int(length) %>% 
    sort
black.ppi.df <- tibble(Symbol = names(black.incident), PPI = black.incident) %>% 
    left_join(ppi.df) %>%
    arrange(desc(PPI))
black.ppi.df$P.Value = map2_dbl(black.ppi.df$PPI, black.ppi.df$Total, GetHyper, black.only, module.membership) %>% 
    p.adjust(method = "fdr") %>% 
    signif(3)
black.ppi.join <- left_join(black.ppi.df, top.table.upper)
black.ppi.lowsig <- filter(black.ppi.join, adj.P.Val_Sleep_dep_vs_Sleep < 0.05 & adj.P.Val_Wake_vs_Sleep < 0.05 & PPI > 4)
Top5Plot("PPI", black.ppi.lowsig, voom.tmm.expr$E, pdata, 'group2', group.levels, "PPI Hub Genes", "black.ppi.lowsig")

black.ppi.df$color <- factor(black.ppi.df$Symbol %in% black.ppi.lowsig$Symbol) %>% fct_recode(red = 'TRUE', blue = 'FALSE')
black.ppi.edges <- igraph::as_data_frame(black.all.graph, "vertices") %>% as_tibble %>%
    set_colnames("Symbol") %>% 
    left_join(select(black.ppi.df, Symbol, color))
black.ppi.plot <- graph_from_data_frame(igraph::as_data_frame(black.all.graph, "edges"), directed = FALSE, vertices = black.ppi.edges)

CairoPDF("black.ppi", width = 50, height = 50)
    plot.igraph(black.ppi.plot, vertex.size = 0.5, vertex.label.dist = 0.05, vertex.label.degree = -pi/2)
dev.off()

#Purple - no PPI
purple.top <- filter(module.reduce, Module == "purple") %>% 
    left_join(top.table) %>%
    filter(adj.P.Val_Sleep_dep_vs_Sleep < 0.05 & adj.P.Val_Sleep_dep_vs_Wake < 0.05)
Top5Plot("kscaled", purple.top, voom.tmm.expr$E, pdata, 'group2', group.levels, "Hub Genes", "purple.top5")

purple.all.graph <- graph_from_edgelist(as.matrix(purple.ppi), directed = FALSE) %>% igraph::simplify()
purple.incident <- map(V(purple.all.graph), incident, graph = purple.all.graph) %>% 
    map_int(length) %>% 
    sort
purple.ppi.df <- tibble(Symbol = names(purple.incident), PPI = purple.incident) %>% 
    left_join(ppi.df) %>%
    arrange(desc(PPI))
purple.ppi.df$P.Value = map2_dbl(purple.ppi.df$PPI, purple.ppi.df$Total, GetHyper, purple.only, module.membership) %>% 
    p.adjust(method = "fdr") %>% 
    signif(3)
purple.ppi.join <- left_join(purple.ppi.df, top.table.upper)
purple.ppi.lowsig <- filter(purple.ppi.join, PPI > 4 & adj.P.Val_Sleep_dep_vs_Sleep < 0.05 & adj.P.Val_Sleep_dep_vs_Wake < 0.05)
Top5Plot("PPI", purple.ppi.lowsig, voom.tmm.expr$E, pdata, 'group2', group.levels, "PPI Hub Genes", "purple.ppi.lowsig")

purple.ppi.df$color <- factor(purple.ppi.df$Symbol %in% purple.ppi.lowsig$Symbol) %>% fct_recode(red = 'TRUE', blue = 'FALSE')
purple.ppi.edges <- igraph::as_data_frame(purple.all.graph, "vertices") %>% as_tibble %>%
    set_colnames("Symbol") %>% 
    left_join(select(purple.ppi.df, Symbol, color))
purple.ppi.plot <- graph_from_data_frame(igraph::as_data_frame(purple.all.graph, "edges"), directed = FALSE, vertices = purple.ppi.edges)

CairoPDF("purple.ppi", width = 50, height = 50)
    plot.igraph(purple.ppi.plot, vertex.size = 0.5, vertex.label.dist = 0.05, vertex.label.degree = -pi/2)
dev.off()

#Pink
pink.top <- filter(module.reduce, Module == "pink") %>% 
    left_join(top.table) %>%
    filter(adj.P.Val_Sleep_dep_vs_Sleep < 0.05 & adj.P.Val_Wake_vs_Sleep < 0.05)
Top5Plot("kscaled", pink.top, voom.tmm.expr$E, pdata, 'group2', group.levels, "Hub Genes", "pink.top5")

pink.all.graph <- graph_from_edgelist(as.matrix(pink.ppi), directed = FALSE) %>% igraph::simplify()
pink.incident <- map(V(pink.all.graph), incident, graph = pink.all.graph) %>% 
    map_int(length) %>% 
    sort
pink.ppi.df <- tibble(Symbol = names(pink.incident), PPI = pink.incident) %>% 
    left_join(ppi.df) %>%
    arrange(desc(PPI))
pink.ppi.df$P.Value = map2_dbl(pink.ppi.df$PPI, pink.ppi.df$Total, GetHyper, pink.only, module.membership) %>% 
    p.adjust(method = "fdr") %>% 
    signif(3)
pink.ppi.join <- left_join(pink.ppi.df, top.table.upper)
pink.ppi.lowsig <- filter(pink.ppi.join, adj.P.Val_Sleep_dep_vs_Sleep < 0.05 & adj.P.Val_Wake_vs_Sleep < 0.05 & PPI > 4)
Top5Plot("PPI", pink.ppi.lowsig, voom.tmm.expr$E, pdata, 'group2', group.levels, "PPI Hub Genes", "pink.ppi.lowsig")

pink.ppi.df$color <- factor(pink.ppi.df$Symbol %in% pink.ppi.lowsig$Symbol) %>% fct_recode(red = 'TRUE', blue = 'FALSE')
pink.ppi.edges <- igraph::as_data_frame(pink.all.graph, "vertices") %>% as_tibble %>%
    set_colnames("Symbol") %>% 
    left_join(select(pink.ppi.df, Symbol, color))
pink.ppi.plot <- graph_from_data_frame(igraph::as_data_frame(pink.all.graph, "edges"), directed = FALSE, vertices = pink.ppi.edges)

CairoPDF("pink.ppi", width = 50, height = 50)
    plot.igraph(pink.ppi.plot, vertex.size = 0.5, vertex.label.dist = 0.05, vertex.label.degree = -pi/2)
dev.off()

#Red
red.top <- filter(module.reduce, Module == "red") %>% 
    left_join(top.table) %>%
    filter(adj.P.Val_Sleep_dep_vs_Sleep < 0.05 & adj.P.Val_Wake_vs_Sleep < 0.05)
Top5Plot("kscaled", red.top, voom.tmm.expr$E, pdata, 'group2', group.levels, "Hub Genes", "red.top5")

red.all.graph <- graph_from_edgelist(as.matrix(red.ppi), directed = FALSE) %>% igraph::simplify()
red.incident <- map(V(red.all.graph), incident, graph = red.all.graph) %>% 
    map_int(length) %>% 
    sort
red.ppi.df <- tibble(Symbol = names(red.incident), PPI = red.incident) %>% 
    left_join(ppi.df) %>%
    arrange(desc(PPI))
red.ppi.df$P.Value = map2_dbl(red.ppi.df$PPI, red.ppi.df$Total, GetHyper, red.only, module.membership) %>% 
    p.adjust(method = "fdr") %>% 
    signif(3)
red.ppi.join <- left_join(red.ppi.df, top.table.upper)
red.ppi.lowsig <- filter(red.ppi.join, adj.P.Val_Sleep_dep_vs_Sleep < 0.05 & adj.P.Val_Wake_vs_Sleep < 0.05 & PPI > 4)
Top5Plot("PPI", red.ppi.lowsig, voom.tmm.expr$E, pdata, 'group2', group.levels, "PPI Hub Genes", "red.ppi.lowsig")

red.ppi.df$color <- factor(red.ppi.df$Symbol %in% red.ppi.lowsig$Symbol) %>% fct_recode(red = 'TRUE', blue = 'FALSE')
red.ppi.edges <- igraph::as_data_frame(red.all.graph, "vertices") %>% as_tibble %>%
    set_colnames("Symbol") %>% 
    left_join(select(red.ppi.df, Symbol, color))
red.ppi.plot <- graph_from_data_frame(igraph::as_data_frame(red.all.graph, "edges"), directed = FALSE, vertices = red.ppi.edges)

CairoPDF("red.ppi", width = 50, height = 50)
    plot.igraph(red.ppi.plot, vertex.size = 0.5, vertex.label.dist = 0.05, vertex.label.degree = -pi/2)
dev.off()

#Blue
blue.top <- filter(module.reduce, Module == "blue") %>% 
    left_join(top.table) %>%
    filter(adj.P.Val_Sleep_dep_vs_Sleep < 0.05 & adj.P.Val_Sleep_dep_vs_Wake < 0.05)
Top5Plot("kscaled", blue.top, voom.tmm.expr$E, pdata, 'group2', group.levels, "Hub Genes", "blue.top5")

blue.all.graph <- graph_from_edgelist(as.matrix(blue.ppi), directed = FALSE) %>% igraph::simplify()
blue.incident <- map(V(blue.all.graph), incident, graph = blue.all.graph) %>% 
    map_int(length) %>% 
    sort
blue.ppi.df <- tibble(Symbol = names(blue.incident), PPI = blue.incident) %>% 
    left_join(ppi.df) %>%
    arrange(desc(PPI))
blue.ppi.df$P.Value = map2_dbl(blue.ppi.df$PPI, blue.ppi.df$Total, GetHyper, blue.only, module.membership) %>% 
    p.adjust(method = "fdr") %>% 
    signif(3)
blue.ppi.join <- left_join(blue.ppi.df, top.table.upper)
blue.ppi.lowsig <- filter(blue.ppi.join, adj.P.Val_Sleep_dep_vs_Sleep < 0.05 & adj.P.Val_Sleep_dep_vs_Wake < 0.05 & PPI > 4)
Top5Plot("PPI", blue.ppi.lowsig, voom.tmm.expr$E, pdata, 'group2', group.levels, "PPI Hub Genes", "blue.ppi.lowsig")

blue.ppi.df$color <- factor(blue.ppi.df$Symbol %in% blue.ppi.lowsig$Symbol) %>% fct_recode(red = 'TRUE', blue = 'FALSE')
blue.ppi.edges <- igraph::as_data_frame(blue.all.graph, "vertices") %>% as_tibble %>%
    set_colnames("Symbol") %>% 
    left_join(select(blue.ppi.df, Symbol, color))
blue.ppi.plot <- graph_from_data_frame(igraph::as_data_frame(blue.all.graph, "edges"), directed = FALSE, vertices = blue.ppi.edges) %>% 
    delete.vertices(which(clusters(.)$membership != 1))

CairoPDF("blue.ppi", width = 50, height = 50)
    plot.igraph(blue.ppi.plot, vertex.size = 0.5, vertex.label.dist = 0.05, vertex.label.degree = -pi/2)
dev.off()

blue.ppi.hubs <- igraph::as_data_frame(blue.all.graph, "both")
blue.ppi.vertices <- set_colnames(blue.ppi.hubs$vertices, "Symbol") %>% as_tibble
blue.ppi.de <- inner_join(blue.ppi.vertices, mutate(blue.top, Symbol = toupper(Symbol))) %>% left_join(blue.ppi.df)
blue.edge.de <- filter(blue.ppi.hubs$edges, from %in% blue.ppi.de$Symbol) %>% 
    filter(to %in% blue.ppi.de$Symbol) %>% as_tibble
blue.ppi.de.graph <- graph_from_data_frame(blue.edge.de, directed = FALSE, vertices = select(blue.ppi.de, Symbol, color)) %>%
    delete.vertices(which(clusters(.)$membership != 1))

CairoPDF("blue.ppi.de", width = 30, height = 30)
    plot.igraph(blue.ppi.de.graph, vertex.size = 1.0, vertex.label.dist = 0.1, vertex.label.degree = -pi/2)
dev.off()

#light green - PPI REJECTED!
lightcyan.top <- filter(module.reduce, Module == "lightcyan") %>% 
    left_join(top.table) %>%
    filter(adj.P.Val_Sleep_dep_vs_Wake < 0.05)
Top5Plot("kscaled", lightcyan.top, voom.tmm.expr$E, pdata, 'group2', group.levels, "Hub Genes", "lightcyan.top5")

#light green - PPI REJECTED!
greenyellow.top <- filter(module.reduce, Module == "greenyellow") %>% 
    left_join(top.table) %>%
    filter(adj.P.Val_Sleep_dep_vs_Sleep < 0.05)
Top5Plot("kscaled", greenyellow.top, voom.tmm.expr$E, pdata, 'group2', group.levels, "Hub Genes", "greenyellow.top5")

#cyan - PPI REJECTED!
cyan.top <- filter(module.reduce, Module == "cyan") %>% 
    left_join(top.table) %>%
    filter(adj.P.Val_Wake_vs_Sleep < 0.05)
Top5Plot("kscaled", cyan.top, voom.tmm.expr$E, pdata, 'group2', group.levels, "Hub Genes", "cyan.top5")

