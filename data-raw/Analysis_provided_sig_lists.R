counts <- count_matrix_assembly("CPM.xlsx")

data_WT <- openxlsx::read.xlsx(here::here("data-raw/Exercise_effect.xlsx"), sheet = 1)
data_KO <- openxlsx::read.xlsx(here::here("data-raw/Exercise_effect.xlsx"), sheet = 2)

bg <-clusterProfiler::bitr(rownames(counts),
                           fromType = "ACCNUM",
                           toType = "ENTREZID",
                           OrgDb = "org.Mm.eg.db",
                           drop = T)
data_WT_entrez <- clusterProfiler::bitr(data_WT$ACCES,
                                        fromType = "ACCNUM",
                                        toType = "ENTREZID",
                                        OrgDb = "org.Mm.eg.db",
                                        drop = T)
data_KO_entrez <- clusterProfiler::bitr(data_KO$ACCES,
                                        fromType = "ACCNUM",
                                        toType = "ENTREZID",
                                        OrgDb = "org.Mm.eg.db",
                                        drop = T)

goResults_WT <- clusterProfiler::enrichGO(gene = data_WT_entrez$ENTREZID,
                                          universe = bg$ENTREZID,
                                          OrgDb = org.Mm.eg.db,
                                          ont = "MF")
goResults_KO <- clusterProfiler::enrichGO(gene = data_KO_entrez$ENTREZID,
                                          universe = bg$ENTREZID,
                                          OrgDb = org.Mm.eg.db,
                                          ont = "MF")
goResults_WT_CC <- clusterProfiler::enrichGO(gene = data_WT_entrez$ENTREZID,
                                             universe = bg$ENTREZID,
                                             OrgDb = org.Mm.eg.db,
                                             ont = "CC")
goResults_KO_CC <- clusterProfiler::enrichGO(gene = data_KO_entrez$ENTREZID,
                                             universe = bg$ENTREZID,
                                             OrgDb = org.Mm.eg.db,
                                             ont = "CC")
enrichplot::dotplot(goResults_WT)
enrichplot::dotplot(goResults_KO)
enrichplot::dotplot(goResults_WT_CC)
enrichplot::dotplot(goResults_KO_CC)
