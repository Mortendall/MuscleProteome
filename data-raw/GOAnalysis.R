counts <- count_matrix_assembly("CPM.xlsx")
metadata <- load_metadata("setup.xlsx")
metadata <- metadata %>%
    tidyr::unite("Group", Genotype:Training, remove = F)
res <- normalizeBetweenArrays(counts, method = "quantile")

#####Go analysis####

limma_data <- list ("Exercise_in_WT" = NA,
                 "Exercise_in_KO" = NA,
                 "KO_in_Control"= NA,
                 "KO_in_Exercise" = NA,
                 "Interaction" = NA
)
for (i in 1:5){
  limma_data[[i]]<-  openxlsx::read.xlsx(here::here("data/limma_results.xlsx"), sheet = i)
}
#generate universe for GO analysis
GO_data <- list ("Exercise_in_WT" = NA,
                 "Exercise_in_KO" = NA,
                 "KO_in_Control"= NA,
                 "KO_in_Exercise" = NA,
                 "Interaction" = NA
)
GO_data <- GO_analysis(limma_data, "MF")

dotplot_MF <- list ("Exercise_in_WT" = NA,
                    "Exercise_in_KO" = NA,
                    "KO_in_Control"= NA,
                    "KO_in_Exercise" = NA,
                    "Interaction" = NA
)


for (i in 1:5){
    dotplot_MF[[i]] <- enrichplot::dotplot(GO_data[[i]])+ ggplot2::ggtitle(str_replace_all(names(GO_data[i]), "_", " "))
}

layout <- (dotplot_MF[[1]] + dotplot_MF[[2]] + dotplot_MF[[3]] + dotplot_MF[[4]]) +
    patchwork::plot_annotation("GO analysis: Molecular Function", theme = theme(plot.title = element_text(hjust = 0.5)))

tiff(here::here("data/GO_MF.tif"), unit = "cm", height = 20, width = 50, res = 300)
layout
dev.off()

#repeat for CC analysis####
GO_data_CC <- list ("Exercise_in_WT" = NA,
                 "Exercise_in_KO" = NA,
                 "KO_in_Control"= NA,
                 "KO_in_Exercise" = NA,
                 "Interaction" = NA
)
GO_data_CC <- GO_analysis(limma_data, "CC")

dotplot_CC <- list ("Exercise_in_WT" = NA,
                    "Exercise_in_KO" = NA,
                    "KO_in_Control"= NA,
                    "KO_in_Exercise" = NA,
                    "Interaction" = NA
)


for (i in 1:5){
    dotplot_CC[[i]] <- enrichplot::dotplot(GO_data_CC[[i]])+ ggplot2::ggtitle(str_replace_all(names(GO_data_CC[i]), "_", " "))
}

layout_CC <- (dotplot_CC[[1]] + dotplot_CC[[2]])  +
    patchwork::plot_annotation("GO analysis: Cellular Compartment", theme = theme(plot.title = element_text(hjust = 0.5)))

tiff(here::here("data/GO_CC.tif"), unit = "cm", height = 20, width = 50, res = 300)
layout_CC
dev.off()

#save data

#openxlsx::write.xlsx(GO_data, here::here("data/GO_data_MF.xlsx"))

#####Generate upset plot for exercise figures#####
GO_data <- list ("Exercise_in_WT" = NA,
                 "Exercise_in_KO" = NA,
                 "KO_in_Control"= NA,
                 "KO_in_Exercise" = NA,
                 "Interaction" = NA
)
for (i in 1:5){
    GO_data[[i]]<-  openxlsx::read.xlsx(here::here("data/GO_data_MF.xlsx"), sheet = i)
}


