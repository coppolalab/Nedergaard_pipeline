#For WGCNA
library(WGCNA)
library(flashClust)
enableWGCNAThreads()
library(BayesFactor)
library(PMCMR)

#For baseline processing
library(limma)
library(R.utils)
library(lumi)
library(biomaRt)
library(matrixStats)

#Plotting
library(ggplot2)
library(extrafont)
library(Cairo)
library(igraph)
library(TeachingDemos)
library(UpSetR)

#Reading and writing tables
library(readr)
library(openxlsx)

#Functional programming
library(magrittr)
library(purrr)

#Data arrangement
library(dplyr)
library(tidyr)
library(broom)

#String operations
library(stringr)
library(tools)

EigengeneANOVA <- function(ME.vector, trait.vector) {
    trait.df <- data.frame(Trait = trait.vector, ME = ME.vector)
    trait.anova <- kruskal.test(ME ~ Trait, trait.df) %>% tidy
    trait.anova$p.value
}

EigengeneBayes <- function(ME.vector, trait.vector) {
    trait.df <- data.frame(Trait = trait.vector, ME = ME.vector)
    trait.anova <- anovaBF(ME ~ Trait, data = trait.df) #%>% anova %>% tidy
    trait.anova
}


EigengeneBoxplot <- function(eigengene.df, status.vector, color) {
    color.column <- str_c("ME", color)
    gene.df <- data.frame(Status = factor(status.vector), Expression = as.vector(eigengene.df[[color.column]]))
    #gene.df$Status %<>% factor(levels = c("Control", "Carrier", "Patient"))

    p <- ggplot(gene.df, aes(x = Status, y = Expression, fill = Status)) + geom_boxplot() + theme_bw()
    p <- p + theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + theme(plot.background = element_blank()) + ggtitle(str_c(capitalize(color), " Module"))
    p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank())
    CairoPDF(color, width = 4, height = 4, bg = "transparent")
    print(p)
    dev.off()
}

EstimateViolinPlot <- function(estimate.df, color, ylimits) {
    estimate.plot <- as.matrix(estimate.df) %>% data.frame %>% select(matches("^Trait")) %>% gather(Group, Estimate) 
    estimate.plot$Group %<>% str_replace("Trait\\.", "") %>% factor()

    p <- ggplot(estimate.plot, aes(x = Group, y = Estimate, fill = Group)) + geom_violin() + geom_boxplot(width = 0.1, outlier.shape = NA)+ theme_bw()
    p <- p + theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(color = "black", size = 1))
    p <- p + theme(plot.background = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank()) + ylim(ylimits) 
    p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    p <- p + ggtitle(str_c(capitalize(color), "Module", sep = " "))
    CairoPDF(str_c(color, "estimate", sep = "."), width = 4, height = 4, bg = "transparent")
    print(p)
    dev.off()
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
    setColWidths(wb, 1, cols = 2, widths = "auto")
    setColWidths(wb, 1, cols = 3, widths = 45)
    setColWidths(wb, 1, cols = 4:ncol(module.table), widths = "auto")
    pageSetup(wb, 1, orientation = "landscape", fitToWidth = TRUE)
    freezePane(wb, 1, firstRow = TRUE)
    showGridLines(wb, 1, showGridLines = TRUE)
    modifyBaseFont(wb, fontSize = 11)
    saveWorkbook(wb, filename, overwrite = TRUE) 
}

PCAPlot <- function(filename, dataset, facet.bool, size.height, size.width) {
    dataset$Module %<>% str_replace("ME", "") 
    p <- ggplot(dataset, aes(x = Gene, y = PCA1, fill = Module, color = Module)) + geom_point()
    p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + theme(axis.ticks.x = element_blank()) + ylab("First Principal Component")
    p <- p + theme(axis.title.x = element_blank(), axis.text.x = element_blank())
    p <- p + scale_color_manual(values = sort(unique(dataset$Module)))
    if (facet.bool == TRUE)
    {
        p <- p + facet_wrap(~ Module)
        p <- p + theme(legend.position = "none")
    } 
    CairoPDF(filename, height = size.height, width = size.width)
    print(p)
    dev.off()
}

MEHeatmap <- function(dataset, ME.genes) {
    color <- as.character(unique(dataset$Module)) %>% str_replace("ME", "")
    colnames(ME.genes) %<>% str_replace("ME", "")
    dataset %<>% select(-Module) %>% scale
    max.dataset <- max(abs(dataset))

    CairoPDF(paste("./modules/", color, sep = ""), width = 21, height = 12)
    par(mar = c(3.5,3,2,3))
    par(oma = c(4,0,2,0))
    plotMat(dataset, zlim = c(-max.dataset, max.dataset))

    ME.genes.plot <- select(ME.genes, Sample.ID, matches(color))
    p <- ggplot(ME.genes.plot, aes_string(x = "Sample.ID", y = color))
    p <- p + geom_bar(stat = "identity") + xlab("Eigengene Expression") 
    p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + theme(axis.text.x = element_blank())  
    p <- p + theme(axis.title.x = element_blank(), axis.ticks.x = element_blank())
    print(p)
    dev.off()
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
    enrichr.df$Log.Bayes.Factor <- log10(enrichr.df$Bayes.Factor)
    enrichr.df$Term %<>% str_replace_all("\\ \\(GO.*$", "") %>% str_replace_all("\\_Homo.*$", "") %>% str_replace_all(",.*$", "")  #Remove any thing after the left parenthesis and convert to all lower case
    enrichr.df$Format.Name <- str_c(enrichr.df$Database, ": ", enrichr.df$Term, " (", enrichr.df$Gene.Count, ")")
    enrichr.df %<>% arrange(Log.Bayes.Factor)
    enrichr.df$Format.Name %<>% factor(levels = enrichr.df$Format.Name)
    enrichr.df.plot <- select(enrichr.df, Format.Name, Log.Bayes.Factor)  

    p <- ggplot(enrichr.df.plot, aes(Format.Name, Log.Bayes.Factor)) + geom_bar(stat = "identity") + coord_flip() + theme_bw() + theme(legend.position = "none") 
    p <- p + theme(axis.title.y = element_blank(), axis.ticks.y = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ylab(expression(paste(Log[10], ' Bayes Factor'))) + theme(plot.background = element_blank())
    p <- p + theme(panel.border = element_rect(color = "black", size = 1))
    p <- p + ggtitle(plot.title)
    CairoPDF(filename, height = plot.height, width = plot.width, bg = "transparent")
    print(p)
    dev.off()
}

GetKscaled <- function(gene.list, module.membership) {
    filter(module.membership, is.element(Symbol, gene.list)) %>% select(Symbol, kscaled) %>% arrange(desc(kscaled))
}

test <- lapply(ls(), function(thing) print(object.size(get(thing)), units = 'auto')) 
names(test) <- ls()
unlist(test) %>% sort

source("../../FRDA project/common_functions.R")
voom.tmm <- ReadRDSgz(file = "../differential_expression/save/voom.tmm.rda")
voom.mad <- rowMads(voom.tmm$E)
voom.tmm.expr <- voom.tmm[voom.mad > 0,]

#Calculate scale free topology measures for different values of power adjacency function
powers <- c(c(1:10), seq(from = 12, to = 39, by = 2))

expr.collapse <- voom.tmm.expr$E %>% t
sft <- pickSoftThreshold(expr.collapse, powerVector = powers, verbose = 5, networkType = "signed", corFnc = bicor, corOptions = list(maxPOutliers = 0.05))
sft.df <- sft$fitIndices
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

#Check if any modules are too similar and merge them.  Possibly not working.
ME.dissimilarity.threshold <- 0.2
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
ME.genes <- merged.genes
SaveRDSgz(ME.genes, file = "./save/me.genes.rda")

CairoPDF("eigengenes", height = 6, width = 8, bg = "transparent")
plotEigengeneNetworks(ME.genes, "", marDendro = c(0,4,1,2), marHeatmap = c(3,5,1,2), plotPreservation = "standard")
dev.off()

EigengeneAOV <- function(ME.vector, trait.vector, contrasts.vector) {
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

#Create heatmap of eigengene significance
format.names <- c("Sleep_grp1_vs_Sleep_dep", "Sleep_grp2_vs_Sleep_dep", "Wake_grp1_vs_Sleep_dep", "Wake_grp2_vs_Sleep_dep", "Sleep_grp2_vs_Sleep_grp1", "Wake_grp1_vs_Sleep_grp1", "Wake_grp2_vs_Sleep_grp1", "Wake_grp1_vs_Sleep_grp2", "Wake_grp2_vs_Sleep_grp2", "Wake_grp2_vs_Wake_grp1")
color.values <- unique(module.colors)
pdata <- ReadRDSgz("../differential_expression/save/targets.rda")
anova.status <- map_dbl(ME.genes, EigengeneANOVA, pdata$group) %>% p.adjust("fdr") %>% signif(3)
bayes.status <- map(ME.genes, EigengeneBayes, pdata$group) 
bf.status <- map(bayes.status, extractBF) %>% map_dbl(extract2, "bf")
posterior.status <- map(bayes.status, posterior, iterations = 100000) 

aov.status <- map(ME.genes, EigengeneAOV, pdata$group, format.names)
status.diff <- map(aov.status, select, Diff) %>% map(t) %>% reduce(rbind) %>% data.frame
rownames(status.diff) <- names(aov.status)
status.pval <- map(aov.status, select, P.value) %>% map(t) %>% reduce(rbind) %>% data.frame
rownames(status.pval) <- names(aov.status)

pval.adjust <- map(status.pval, p.adjust, method = "fdr", n = length(color.values)) %>% reduce(cbind) %>% data.frame
rownames(pval.adjust) <- rownames(status.pval)
colnames(pval.adjust) <- paste(colnames(status.pval), ".pval")

text.matrix.traits <- str_c(signif(as.matrix(status.diff), 2), '\n(', signif(as.matrix(pval.adjust), 1), ')')
dim(text.matrix.traits) = dim(status.diff)

heatmap.range <- c(min(as.matrix(status.diff)) * 1.1, max(as.matrix(status.diff)) * 1.1)
width.dynamic <- 3 + (1 * ncol(text.matrix.traits))

CairoPDF("module_trait_relationships", width = width.dynamic, height = 10, bg = "transparent")
par(mar = c(8, 11, 3, 3))
labeledHeatmap(Matrix = as.matrix(status.diff), xLabels = colnames(status.diff), yLabels = colnames(ME.genes), ySymbols = colnames(ME.genes), yColorLabels = TRUE, colors = greenWhiteRed(50), textMatrix = text.matrix.traits, setStdMargins = F, zlim = heatmap.range, main = "")
dev.off()

#Eigengene plots
EigengeneBoxplot(ME.genes, pdata$group, "purple")
EigengeneBoxplot(ME.genes, pdata$group, "magenta")

EstimateViolinPlot(posterior.status$MEpurple, "purple", c(-0.25, 0.25))
EstimateViolinPlot(posterior.status$MEmagenta, "magenta", c(-0.30, 0.30))

ME.genes.plot <- mutate(ME.genes, Group = pdata$group) %>% gather(Module.Color, Eigengene, -Group) 
ME.genes.plot$Module.Color %<>% str_replace("ME", "")

p <- ggplot(ME.genes.plot, aes(x = Group, y = Eigengene, color = Group)) + geom_boxplot(width = 0.5) + geom_jitter(width = 0.5)
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
p <- p + scale_fill_manual(values = sort(unique(ME.genes.plot$Module.Color)))
p <- p + facet_wrap(~ Module.Color, ncol = 4, scales = "free") + theme(plot.background = element_blank())
CairoPDF("eigengene_plots", height = 12, width = 18, bg = "transparent")
print(p)
dev.off()

#Generate network statistics
all.degrees <- intramodularConnectivity(adjacency.expr, module.colors)
gene.info <- data.frame(Ensembl.ID = rownames(all.degrees), Module = module.colors, all.degrees, stringsAsFactors = FALSE) %>% arrange(Module)
gene.info$kscaled <- by(gene.info, gene.info$Module, select, kWithin) %>% map(function(kWithin) kWithin / max(kWithin)) %>% reduce(c)

gene.module.membership <- data.frame(bicor(expr.collapse, ME.genes, maxPOutliers = 0.05))
module.membership.pvalue <- data.frame(corPvalueStudent(as.matrix(gene.module.membership), nrow(expr.collapse)))
colnames(gene.module.membership) <- str_replace(names(ME.genes),"ME", "MM.")
colnames(module.membership.pvalue) <- str_replace(names(ME.genes),"ME", "MM.pvalue.")
gene.module.membership$Ensembl.ID <- rownames(gene.module.membership)
module.membership.pvalue$Ensembl.ID <- rownames(module.membership.pvalue)

ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
vega <- useMart("ENSEMBL_MART_VEGA", dataset = "mmusculus_gene_vega")

#Annotate with BioMart (Ensembl and Vega)
bm.ensembl <- getBM(attributes = c('ensembl_gene_id', 'mgi_symbol', 'description', 'gene_biotype'), filters = 'ensembl_gene_id', values = gene.info$Ensembl.ID, mart = ensembl)
bm.ensembl$description %<>% str_replace_all(" \\[.*\\]$", "")
colnames(bm.ensembl) <- c("Ensembl.ID", "Symbol", "Description", "Gene.Type")
SaveRDSgz(bm.ensembl, "./save/bm.ensembl.rda")

bm.vega <- getBM(attributes = c('ens_gene', 'vega_gene_id', 'external_gene_name', 'description', 'gene_biotype'), filters = 'ens_gene', values = gene.info$Ensembl.ID, mart = vega)
bm.vega$description %<>% str_replace_all(" \\[.*\\]$", "")
colnames(bm.vega) <- c("Ensembl.ID", "Vega.ID", "Gene.Name", "Description.Vega", "Gene.Type.Vega")
SaveRDSgz(bm.vega, "./save/bm.vega.rda")

combine.annot <- left_join(bm.ensembl, bm.vega) 
SaveRDSgz(combine.annot, "./save/combine.annot.rda")

top.table.all <- map2(1:10, format.names, FormatTopTable, fit.voom) %>% 
    reduce(left_join) %>%
    left_join(combine.annot) %>% 
    select(Ensembl.Gene.ID, Symbol, Description, Gene.Type, dplyr::contains("logFC"), dplyr::contains("adj.P.Val"), dplyr::contains("P.Value"), AveExpr, Vega.ID:Gene.Type.Vega, dplyr::contains("t_"), dplyr::contains("B_"))
SaveRDSgz(top.table.all, "./save/top.table.all.rda")

module.membership <- left_join(gene.info, combine.annot) %>% 
    left_join(gene.module.membership) %>% 
    left_join(module.membership.pvalue) %>%
    select(Ensembl.ID, Symbol:Gene.Type, Module:kscaled, MM.blue:MM.salmon, MM.pvalue.blue:MM.pvalue.salmon, Vega.ID:Gene.Type.Vega) %>%
    arrange(Module, desc(kscaled))
ModuleWorkbook(module.membership, "module_membership.xlsx")

#PCA plots
all.smooth <- apply(ME.genes, 2, smooth.spline, spar = 0.4) %>% map(extract2, "y")
smooth.df <- data.frame(all.smooth)
colnames(smooth.df) <- names(all.smooth)
smooth.df$Gene <- as.factor(1:nrow(smooth.df))
smooth.plot <- gather(smooth.df, Module, PCA1, -Gene)

PCAPlot("all_principal_components", smooth.plot, FALSE, 10, 15)
PCAPlot("facet_principal_components", smooth.plot, TRUE, 13, 20)

#Heatmaps
sample.ids <- factor(rownames(expr.collapse), levels = rownames(expr.collapse))
ME.genes.plot <- mutate(data.frame(ME.genes), Sample.ID = sample.ids)
expr.data.plot <- data.frame(t(expr.collapse), Module = module.colors)
split(expr.data.plot, expr.data.plot$Module) %>% map(MEHeatmap, ME.genes.plot)

#Enrichr
source("../../code/GO/enrichr.R")

modules.filter <- filter(module.membership, Gene.Type == "protein_coding") %>% select(Ensembl.ID, Symbol, Module)
enrichr.terms <- c("GO_Biological_Process_2015", "GO_Molecular_Function_2015", "KEGG_2016", "Reactome_2016", "Human_Gene_Atlas", "GTEx_Tissue_Sample_Gene_Expression_Profiles_up") 
color.names <- unique(module.colors) %>% sort
trap1 <- map(color.names, EnrichrSubmit, modules.filter, enrichr.terms, FALSE)

EnrichrPlot <- function(enrichr.df, filename, plot.title, plot.height = 5, plot.width = 8) {
    enrichr.df$Gene.Count <- map(enrichr.df$Genes, str_split, ",") %>% map(unlist) %>% map_int(length)
    enrichr.df$Adj.P.value <- p.adjust(enrichr.df$P.value, method = "fdr")
    enrichr.df$Log.P.value <- -log10(enrichr.df$Adj.P.value)
    enrichr.df$Term %<>% str_replace_all("\\ \\(.*$", "") %>% str_replace_all("\\_Homo.*$", "")  #Remove any thing after the left parenthesis and convert to all lower case
    enrichr.df$Format.Name <- str_c(enrichr.df$Term, " (", enrichr.df$Gene.Count, ")")
    enrichr.df %<>% arrange(Log.P.value)
    enrichr.df$Format.Name %<>% factor(levels = enrichr.df$Format.Name)

    p <- ggplot(enrichr.df, aes(Format.Name, Log.P.value)) + geom_bar(stat = "identity", size = 1) 
    p <- p + coord_flip() + theme_bw() + theme(axis.title.y = element_blank(), axis.ticks.y = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ylab(expression(paste(Log[10], ' P-Value')))
    p <- p + theme(plot.background = element_blank(), legend.background = element_blank(), axis.text.y = element_text(size = 12), panel.border = element_rect(color = "black", size = 1))
    CairoPDF(filename, height = plot.height, width = plot.width, bg = "transparent")
    print(p)
    dev.off()
}

#Enrichr plots
FilterEnrichr <- function(enrichr.df, size = 100) {
    enrichr.df$Num.Genes <- map(enrichr.df$Genes, str_split, ",") %>% map(extract2, 1) %>% map_int(length)
    enrichr.filter <- filter(enrichr.df, Num.Genes > 4) %>% filter(P.value < 0.05)
    if (nrow(enrichr.df) > size) {
        enrichr.filter %<>% slice(1:size)
    }
    enrichr.filter
}

#Purple
purple.only <- filter(module.filter, Module == "purple")$Symbol
purple.gobiol.file <- "./enrichr/purple/purple_GO_Biological_Process_2015.xlsx"
purple.gobiol <- read.xlsx(purple.gobiol.file) 
purple.gobiol.filter <- FilterEnrichr(purple.gobiol)
GetKappaCluster(file_path_sans_ext(purple.gobiol.file), purple.gobiol.filter, purple.only)
purple.gobiol.final <- slice(purple.gobiol, c(1, 4, 15, 13))

purple.gomole.file <- "./enrichr/purple/purple_GO_Molecular_Function_2015.xlsx"
purple.gomole <- read.xlsx(purple.gomole.file) 
purple.gomole.filter <- FilterEnrichr(purple.gomole)
GetKappaCluster(file_path_sans_ext(purple.gomole.file), purple.gomole.filter, purple.only)
purple.gomole.final <- slice(purple.gomole, c(1))

purple.reactome.file <- "./enrichr/purple/purple_Reactome_2016.xlsx"
purple.reactome <- read.xlsx(purple.reactome.file) 
purple.reactome.filter <- FilterEnrichr(purple.reactome)
GetKappaCluster(file_path_sans_ext(purple.reactome.file), purple.reactome.filter, purple.only)
purple.reactome.final <- slice(purple.reactome, c(35, 54))

purple.kegg.file <- "./enrichr/purple/purple_KEGG_2016.xlsx"
purple.kegg <- read.xlsx(purple.kegg.file) 
purple.kegg.filter <- FilterEnrichr(purple.kegg)
GetKappaCluster(file_path_sans_ext(purple.kegg.file), purple.kegg.filter, purple.only)
purple.kegg.final <- slice(purple.kegg, c(2, 4))

purple.enrichr.final <- rbind(purple.gobiol.final, purple.gomole.final, purple.reactome.final, purple.kegg.final)
EnrichrPlot(purple.enrichr.final, "purple.enrichr")

#Magenta
magenta.only <- filter(module.filter, Module == "magenta")$Symbol
magenta.gobiol.file <- "./enrichr/magenta/magenta_GO_Biological_Process_2015.xlsx"
magenta.gobiol <- read.xlsx(magenta.gobiol.file) 
magenta.gobiol.filter <- FilterEnrichr(magenta.gobiol)
GetKappaCluster(file_path_sans_ext(magenta.gobiol.file), magenta.gobiol.filter, magenta.only)
magenta.gobiol.final <- slice(magenta.gobiol, c(1, 32))

magenta.gomole.file <- "./enrichr/magenta/magenta_GO_Molecular_Function_2015.xlsx"
magenta.gomole <- read.xlsx(magenta.gomole.file) 
magenta.gomole.filter <- FilterEnrichr(magenta.gomole)
GetKappaCluster(file_path_sans_ext(magenta.gomole.file), magenta.gomole.filter, magenta.only)
magenta.gomole.final <- slice(magenta.gomole, c(1))

magenta.reactome.file <- "./enrichr/magenta/magenta_Reactome_2016.xlsx"
magenta.reactome <- read.xlsx(magenta.reactome.file) 
magenta.reactome.filter <- FilterEnrichr(magenta.reactome)
GetKappaCluster(file_path_sans_ext(magenta.reactome.file), magenta.reactome.filter, magenta.only)
magenta.reactome.final <- slice(magenta.reactome, c(8, 3))

magenta.kegg.file <- "./enrichr/magenta/magenta_KEGG_2016.xlsx"
magenta.kegg <- read.xlsx(magenta.kegg.file) 
magenta.kegg.filter <- FilterEnrichr(magenta.kegg)
GetKappaCluster(file_path_sans_ext(magenta.kegg.file), magenta.kegg.filter, magenta.only)
magenta.kegg.final <- slice(magenta.kegg, c(8))

magenta.enrichr.final <- rbind(magenta.gobiol.final, magenta.reactome.final, magenta.kegg.final)
EnrichrPlot(magenta.enrichr.final, "magenta.enrichr")
