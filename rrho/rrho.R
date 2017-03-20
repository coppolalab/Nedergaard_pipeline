library(Biobase)
library(annotate)
library(WGCNA)
library(biomaRt)
library(broom)
library(RRHO)

library(stringr)
library(readr)
library(openxlsx)

library(rlist)
library(reshape2)
library(dplyr)
library(purrr)
library(magrittr)
library(tidyr)
library(vadr)

library(ggplot2)
library(Cairo)

GetRRHO <- function(colname1, colname2, dataset1, dataset2, symbolname, output.suffix, stepsize = 100) {
    print(colname1)
    print(colname2)
    subset1 <- select_(dataset1, symbolname, str_c("Log.Pvalue.", colname1))
    subset2 <- select_(dataset2, symbolname, str_c("Log.Pvalue.", colname2))
    output.dir <- str_c("./rrho", output.suffix, sep = "_")
    dir.create(output.dir, showWarnings = FALSE)
    rrho.out <- RRHO(subset1, subset2, alternative = "enrichment", stepsize = stepsize, labels = c(colname1, colname2), plots = TRUE, outputdir = output.dir, BY = TRUE)
    return(rrho.out)
}


