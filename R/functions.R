library(fs)
library(vroom)
library(Biobase)
library(here)
library(magrittr)
library(data.table)
library(clusterProfiler)
library(edgeR)
library(openxlsx)
library(pheatmap)
library(gridExtra)
library(PoiClaClu)
library(RColorBrewer)
library(limma)
library(GO.db)
library(stringr)
library(dplyr)
library(ggplot2)
library(org.Mm.eg.db)


#' count_matrix_loader
#'
#' @param file_type your raw data file type (presently optimized for txt)
#'
#' @return a count matrix containing raw count values for the experiment with group annotations.

count_matrix_assembly <- function(file_type){
    count_file <- fs::dir_ls(here::here("data-raw/"),
                             regexp = file_type,
                             recurse = TRUE)
    count_matrix <- openxlsx::read.xlsx(xlsxFile = count_file)
    rownames(count_matrix) <- count_matrix$ACCESSION
    count_matrix <- count_matrix %>%
        dplyr::select(-ACCESSION, -SYMBOL)

    return(count_matrix)
}

#' Load metadata and sort them according to count matrix input
#'
#' @param file_name the name of teh metadata file (default "metadata.csv")
#'
#' @return metadata file sorted according to count matrix order

load_metadata <- function(file_name) {
    data_file <- fs::dir_ls(here::here("data-raw/"),
                            regexp = file_name,
                            recurse = T)
    metadata <- openxlsx::read.xlsx(xlsxFile = data_file)
    return(metadata)
}



#Run the GO analysis
GO_analysis <- function(data_sheet, ontology){
    bg <- clusterProfiler::bitr(data_sheet[[1]]$SYMBOL,
                                fromType = "SYMBOL",
                                toType = "ENTREZID",
                                OrgDb = "org.Mm.eg.db",
                                drop = T)
    for (i in 1:length(data_sheet)){
        data_sig <- data_sheet[[i]] %>%
            dplyr::filter(adj.P.Val < 0.05)
        entrez_id <- clusterProfiler::bitr(data_sig$SYMBOL,
                                           fromType = "SYMBOL",
                                           toType = "ENTREZID",
                                           OrgDb = "org.Mm.eg.db",
                                           drop = T)
        GO_data[[i]]<- clusterProfiler::enrichGO(gene = entrez_id$ENTREZID,
                                                 universe = bg$ENTREZID,
                                                 OrgDb = org.Mm.eg.db,
                                                 ont = ontology)
    }
    for (i in 1:5){
        GO_data[[i]]<- clusterProfiler::setReadable(GO_data[[i]], OrgDb = "org.Mm.eg.db", keyType = "ENTREZID")
    }
    return(GO_data)
}

