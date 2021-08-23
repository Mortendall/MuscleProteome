counts <- count_matrix_assembly("CPM.xlsx")
metadata <- load_metadata("setup.xlsx")
metadata <- metadata %>%
    tidyr::unite("Group", Genotype:Training, remove = F)
res <- normalizeBetweenArrays(counts, method = "quantile")


# Missing samples are not present so not really necessary
# missingSamples <- data.table(is.na(res), keep.rownames = TRUE) %>%
#     melt(measure.vars = colnames(res), variable.name = "ID")
# setnames(missingSamples, "rn", "Accession")
# metadata$ID <- as.character(metadata$ID)
#
# missingSamples <- merge(metadata, missingSamples, by = "ID")
#
# missingSamples_DF <- as.data.table(missingSamples)
# missingSamples_DF <- missingSamples_DF[, .(nMissing = sum(value)), by = c("rn", "Group")]
#
#     missingSamples_DF <- reshape2::dcast(missingSamples_DF, rn ~ Group, value.var = "nMissing")


#Check MDS plot

mdsData <- limma::plotMDS(res, plot = F)

all(rownames(mdsData) == metadata$sample)
mdsData <- cbind(metadata, mdsData)
ggplot(mdsData, aes(x = eigen.vectors.1, y = eigen.vectors.2, colour = Group)) +
    geom_point() +
    scale_color_brewer(name = "Significant", type = "qual", palette = "Set1")

ggplot(mdsData, aes(x = eigen.vectors.1, y = eigen.vectors.3, colour = Group)) +
    geom_point() +
    scale_color_brewer(name = "Significant", type = "qual", palette = "Set1")

ggplot(mdsData, aes(x = eigen.vectors.2, y = eigen.vectors.3, colour = Group)) +
    geom_point() +
    scale_color_brewer(name = "Significant", type = "qual", palette = "Set1")

#Generate design matrix for linear model
design <- model.matrix(~ 0 + Group, metadata)
colnames(design) <- str_remove_all(colnames(design), "Group")
cont.matrix <- makeContrasts(
    Exercise_in_WT = wt_tr - wt_sed,
    Exercise_in_KO = ko_tr - ko_sed,
    KO_in_Control = ko_sed - wt_sed,
    KO_in_Exercise = ko_tr - wt_tr,
    Interaction = (ko_tr - ko_sed) - (wt_tr - wt_sed),
    levels = design
)

#Run linear modeling
fit <- lmFit(res, design = design, method = "robust")
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)

#Generate
resultTables <- list(
    Exercise_in_WT = topTable(fit2, coef = "Exercise_in_WT", number = Inf, p.value = 1) %>% data.table(keep.rownames = TRUE),
    Exercise_in_KO = topTable(fit2, coef = "Exercise_in_KO", number = Inf, p.value = 1) %>% data.table(keep.rownames = TRUE),
    KO_in_Control = topTable(fit2, coef = "KO_in_Control", number = Inf, p.value = 1) %>% data.table(keep.rownames = TRUE),
    KO_in_Exercise = topTable(fit2, coef = "KO_in_Exercise", number = Inf, p.value = 1) %>% data.table(keep.rownames = TRUE),
    Interaction = topTable(fit2, coef = "Interaction", number = Inf, p.value = 1) %>% data.table(keep.rownames = TRUE)
)

lapply(resultTables, setnames, "rn", "Accession")
conv <- clusterProfiler::bitr(resultTables[[1]]$Accession,
                              fromType = "ACCNUM",
                              toType = "SYMBOL",
                              OrgDb = "org.Mm.eg.db")

conv <- as.data.table(conv)

for (i in 1:5){
    resultTables[[i]]<- dplyr::left_join(resultTables[[i]], conv, by = c("Accession" = "ACCNUM"))
}
resultTables_sig <- resultTables

for (i in 1:5){
    resultTables_sig[[i]]<- resultTables_sig[[i]] %>%
        dplyr::filter(adj.P.Val < 0.05)
}
#write.xlsx(resultTables, file = here::here("data/limma_results.xlsx"))
