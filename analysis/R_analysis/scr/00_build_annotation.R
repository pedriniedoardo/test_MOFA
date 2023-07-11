# python version ----------------------------------------------------------
Sys.setenv(RETICULATE_PYTHON = "/home/edo/miniconda3/envs/env_MOFA/bin/python")
library(reticulate)
reticulate::use_condaenv(condaenv = "env_MOFA")

py_config()

# 2Load libraries and data ------------------------------------------------
library(tidyverse)
library(msigdbr)
library(multiGSEA)
library(metaboliteIDmapping)

# read in the annotations -------------------------------------------------
all_gene_sets <- msigdbr(species = "Homo sapiens")
msigdbr_collections() %>% 
  print(n=30)

# REACTOME ----------------------------------------------------------------
# filter only the reactome annotation
df_reactome <- all_gene_sets %>%
  dplyr::filter(gs_subcat == "CP:REACTOME") %>% 
  dplyr::select(gs_name,gene_symbol)

# define the unique genes in the dataset
df_genes_reactome <- data.frame(gene=unique(df_reactome$gene_symbol))

# make a list by terms
list_reactome_terms <- df_reactome %>% 
  split(f=.$gs_name)

# iteratively build a matrix of genes x terms
list_reactome_df <- pmap(list(list_reactome_terms,names(list_reactome_terms)),function(x,name){
  df_genes_reactome %>%
    mutate(test = as.numeric(gene %in% x$gene_symbol)) %>% 
    dplyr::rename(!!name := "test")
})

# reduce the list to a dataframe
reactome_t2g <- list_reactome_df %>% 
  purrr::reduce(left_join,by=c("gene")) %>% 
  pivot_longer(names_to = "terms",values_to = "presence",-gene) %>% 
  pivot_wider(names_from = "gene",values_from = "presence") %>% 
  column_to_rownames("terms") %>% 
  as.matrix()

# save the matrix as an object
saveRDS(reactome_t2g,"../../out/object/reactome_t2g.rds")

# check validity of the matrix --------------------------------------------
# check that the matrix is correct
test_matrix <- reactome_t2g %>% 
  data.frame() %>% 
  rownames_to_column("term") %>% 
  pivot_longer(names_to = "gene",values_to = "presence",-term) %>% 
  # pull(term) %>% unique()
  filter(term=="REACTOME_ABC_TRANSPORTER_DISORDERS") %>% 
  filter(presence == 1) %>% 
  pull(gene) %>% 
  unique()

# confrime they are the same genes
test_orig <- df_reactome %>% 
  filter(gs_name=="REACTOME_ABC_TRANSPORTER_DISORDERS") %>% pull(gene_symbol) %>% unique()
# they are different because they have some duplicated gene namest
df_reactome %>% 
  filter(gs_name=="REACTOME_ABC_TRANSPORTER_DISORDERS") %>% pull(gene_symbol) %>% table() %>% sort()

sum(!test_matrix %in% test_orig)

# KEGG --------------------------------------------------------------------
# filter only the reactome annotation
df_kegg <- all_gene_sets %>%
  dplyr::filter(gs_subcat == "CP:KEGG") %>% 
  dplyr::select(gs_name,gene_symbol)

# define the unique genes in the dataset
df_genes_kegg <- data.frame(gene=unique(df_kegg$gene_symbol))

# make a list by terms
list_kegg_terms <- df_kegg %>% 
  split(f=.$gs_name)

# iteratively build a matrix of genes x terms
list_kegg_df <- pmap(list(list_kegg_terms,names(list_kegg_terms)),function(x,name){
  df_genes_kegg %>%
    mutate(test = as.numeric(gene %in% x$gene_symbol)) %>% 
    dplyr::rename(!!name := "test")
})

# reduce the list to a dataframe
kegg_t2g <- list_kegg_df %>% 
  purrr::reduce(left_join,by=c("gene")) %>% 
  pivot_longer(names_to = "terms",values_to = "presence",-gene) %>% 
  pivot_wider(names_from = "gene",values_from = "presence") %>% 
  column_to_rownames("terms") %>% 
  as.matrix()

# save the matrix as an object
saveRDS(kegg_t2g,"../../out/object/kegg_t2g.rds")

# check validity of the matrix --------------------------------------------
# check that the matrix is correct
test_matrix <- kegg_t2g %>% 
  data.frame() %>% 
  rownames_to_column("term") %>% 
  pivot_longer(names_to = "gene",values_to = "presence",-term) %>% 
  # pull(term) %>% unique()
  filter(term=="KEGG_ABC_TRANSPORTERS") %>% 
  filter(presence == 1) %>% 
  pull(gene) %>% 
  unique()

# confrime they are the same genes
test_orig <- df_kegg %>% 
  filter(gs_name=="KEGG_ABC_TRANSPORTERS") %>% pull(gene_symbol) %>% unique()
# they are different because they have some duplicated gene namest
df_kegg %>% 
  filter(gs_name=="KEGG_ABC_TRANSPORTERS") %>% pull(gene_symbol) %>% table() %>% sort()

sum(!test_matrix %in% test_orig)

# BP ----------------------------------------------------------------------
# filter only the reactome annotation
df_bp <- all_gene_sets %>%
  dplyr::filter(gs_subcat == "GO:BP") %>% 
  dplyr::select(gs_name,gene_symbol)

# define the unique genes in the dataset
df_genes_bp <- data.frame(gene=unique(df_bp$gene_symbol))

# make a list by terms
list_bp_terms <- df_bp %>% 
  split(f=.$gs_name)

# iteratively build a matrix of genes x terms
list_bp_df <- pmap(list(list_bp_terms,names(list_bp_terms)),function(x,name){
  df_genes_bp %>%
    mutate(test = as.numeric(gene %in% x$gene_symbol)) %>% 
    dplyr::rename(!!name := "test")
})

# reduce the list to a dataframe
bp_t2g <- list_bp_df %>% 
  purrr::reduce(left_join,by=c("gene")) %>% 
  column_to_rownames("gene") %>% 
  as.matrix() %>% 
  t()

bp_t2g[1:10,1:10]

# save the matrix as an object
saveRDS(bp_t2g,"../../out/object/bp_t2g.rds")

# check validity of the matrix --------------------------------------------
# check that the matrix is correct
test_matrix <- bp_t2g %>% 
  data.frame() %>% 
  rownames_to_column("term") %>% 
  pivot_longer(names_to = "gene",values_to = "presence",-term) %>% 
  # pull(term) %>% unique()
  filter(term=="GOBP_2FE_2S_CLUSTER_ASSEMBLY") %>% 
  filter(presence == 1) %>% 
  pull(gene) %>% 
  unique()

# confrime they are the same genes
test_orig <- df_bp %>% 
  filter(gs_name=="GOBP_2FE_2S_CLUSTER_ASSEMBLY") %>% pull(gene_symbol) %>% unique()
# they are different because they have some duplicated gene namest
df_bp %>% 
  filter(gs_name=="GOBP_2FE_2S_CLUSTER_ASSEMBLY") %>% pull(gene_symbol) %>% table() %>% sort()

sum(!test_matrix %in% test_orig)

# metabolites annotation --------------------------------------------------
# filter only the reactome annotation
# MOFAobject@data$metabolic
df_metabolites <- read_tsv("../../data/path_metadata.tsv") %>%
  dplyr::select(sub_pathway,metabolite)

# define the unique genes in the dataset
df_feature_metabolites <- data.frame(metabolite=unique(df_metabolites$metabolite))

# make a list by terms
list_metabolites_terms <- df_metabolites %>% 
  split(f=.$sub_pathway)

# iteratively build a matrix of genes x terms
list_metabolites_df <- pmap(list(list_metabolites_terms,names(list_metabolites_terms)),function(x,name){
  df_feature_metabolites %>%
    mutate(test = as.numeric(metabolite %in% x$metabolite)) %>% 
    dplyr::rename(!!name := "test")
})

# reduce the list to a dataframe
metabolites_t2g <- list_metabolites_df %>% 
  purrr::reduce(left_join,by=c("metabolite")) %>% 
  pivot_longer(names_to = "terms",values_to = "presence",-metabolite) %>% 
  pivot_wider(names_from = "metabolite",values_from = "presence") %>% 
  column_to_rownames("terms") %>% 
  as.matrix()

# save the matrix as an object
saveRDS(metabolites_t2g,"../../out/object/metabolites_t2g.rds")

# check validity of the matrix --------------------------------------------
# check that the matrix is correct
test_matrix <- metabolites_t2g %>% 
  data.frame() %>% 
  rownames_to_column("term") %>% 
  pivot_longer(names_to = "gene",values_to = "presence",-term) %>% 
  # pull(term) %>% unique()
  filter(term=="Acetylated Peptides") %>% 
  filter(presence == 1) %>% 
  pull(gene) %>% 
  unique()

# confrime they are the same genes
test_orig <- df_metabolites %>% 
  filter(sub_pathway=="Acetylated Peptides") %>% pull(metabolite) %>% unique()

sum(!test_matrix %in% test_orig)

# build annotation for KEGG multiGSEA -------------------------------------
# read in the compounds IDs
LUT_metabol <- read_tsv("../../data/path_metadata.tsv")

# modify the table to link HMDB and metabolites
LUT_metabolitc_fix <- LUT_metabol %>% 
  split(f = .$metabolite) %>% 
  lapply(function(x){
    # split the kegg value
    hmdb_id <- str_split(x$group_hmdb,pattern = ",") %>% unlist()
    data.frame(hmdb = hmdb_id) %>%
      mutate(metabolite = x$metabolite)
  }) %>% 
  bind_rows()

# use the annotation fo the kegg from multiGSEA
# pathways_kegg_multiGSEA <- getMultiOmicsFeatures(dbs = "kegg",
#                                                  layer = c("transcriptome","metabolome"),
#                                                  returnTranscriptome = "SYMBOL",
#                                                  # returnProteome = "SYMBOL",
#                                                  returnMetabolome = "HMDB",
#                                                  organism = "hsapiens",
#                                                  useLocal = FALSE)
# # save the object
# saveRDS(pathways_kegg_multiGSEA,"../../out/object/pathways_hs_METABOLhmdb_RNAsymbol_kegg.rds")
pathways_kegg_multiGSEA <- readRDS("../../out/object/pathways_hs_METABOLhmdb_RNAsymbol_kegg.rds")

# define the template for the metabolites and transcript
df_metabol_kegg_multiGSEA <- data.frame(hmdb=unique(unlist(pathways_kegg_multiGSEA$metabolome)))
#
df_gene_kegg_multiGSEA <- data.frame(gene=unique(unlist(pathways_kegg_multiGSEA$transcriptome)))

# build the matrix metabolite and transcript
list_kegg_df_metabol <- pmap(list(pathways_kegg_multiGSEA$metabolome,
                                  names(pathways_kegg_multiGSEA$metabolome)),
                             function(x,name){
                               df_metabol_kegg_multiGSEA %>%
                                 mutate(test = as.numeric(hmdb %in% x)) %>% 
                                 dplyr::rename(!!name := "test")
                             })
lapply(list_kegg_df_metabol,dim)

# # for each element convert the HMDB to metabolites
# list_kegg_df_metabol_fix <- lapply(list_kegg_df_metabol,function(x){
#   x %>%
#     left_join(LUT_metabolitc_fix %>% filter(!is.na(hmdb)),
#               by = c("hmdb")) %>%
#     dplyr::filter(!is.na(metabolite)) %>% 
#     # if I use the metabolite the reduce funciton will take forever to run
#     dplyr::select(-metabolite)
# })
# lapply(list_kegg_df_metabol_fix,dim)

#
list_kegg_df_gene <- pmap(list(pathways_kegg_multiGSEA$transcriptome,
                               names(pathways_kegg_multiGSEA$transcriptome)),
                          function(x,name){
                            df_gene_kegg_multiGSEA %>%
                              mutate(test = as.numeric(gene %in% x)) %>%
                              dplyr::rename(!!name := "test")
                          })

# length(list_kegg_df_metabol)
# reduce the list to a dataframe for both gene and metabolites
kegg_t2m_multiGSEA_long <- list_kegg_df_metabol %>%
  purrr::reduce(left_join,by=c("hmdb")) %>%
  pivot_longer(names_to = "terms",values_to = "presence",-hmdb)

kegg_t2m_multiGSEA_wide <- list_kegg_df_metabol %>%
  purrr::reduce(left_join,by=c("hmdb")) %>%
  pivot_longer(names_to = "terms",values_to = "presence",-hmdb) %>% 
  pivot_wider(names_from = "hmdb",values_from = "presence") %>%
  column_to_rownames("terms") %>%
  as.matrix()

dim(kegg_t2m_multiGSEA_wide)
kegg_t2m_multiGSEA_wide[1:10,1:10]
#
# list_test <- LUT_metabolitc_fix %>% filter(!is.na(hmdb)) %>% 
#   split(f = .$hmdb)

list_test2 <- LUT_metabolitc_fix %>% filter(!is.na(hmdb)) %>% 
  split(f = .$metabolite)

# keep only the first one of the hmdb
LUT_metabolitc_fix_shortlist <- lapply(list_test2, function(x){
  x %>% 
    slice(1)
}) %>% 
  bind_rows()

# # remnove those hmdb to avoid issues
# hmdb_test <- lapply(list_test, nrow) %>% unlist() %>% 
#   data.frame(row = .) %>% 
#   rownames_to_column("hmdb") %>% 
#   filter(row > 1) %>% 
#   pull(hmdb)
# 
# metabolite_test <- lapply(list_test2, nrow) %>% unlist() %>% 
#   data.frame(row = .) %>% 
#   rownames_to_column("metabolite") %>% 
#   filter(row > 1) %>% 
#   pull(hmdb)
# 
# LUT_metabolitc_fix %>% 
#   filter(metabolite == "isobar_hexose_diphosphates")

kegg_t2m_multiGSEA <- kegg_t2m_multiGSEA_long %>% 
  # dplyr::filter(!hmdb %in% hmdb_test) %>% 
  left_join(LUT_metabolitc_fix_shortlist %>% filter(!is.na(hmdb)),by = c("hmdb")) %>% 
  # filter out some of the metabolites
  # filter the instances where the metabolite is not present
  dplyr::filter(!is.na(metabolite)) %>%
  # filter(hmdb == "HMDB0000568",metabolite == "arabitol_xylitol")
  dplyr::select(-hmdb) %>% 
  pivot_wider(names_from = "metabolite",values_from = "presence") %>%
  column_to_rownames("terms") %>%
  as.matrix()

dim(kegg_t2m_multiGSEA)
kegg_t2m_multiGSEA[1:10,1:10]

#
kegg_t2g_multiGSEA <- list_kegg_df_gene %>%
  purrr::reduce(left_join,by=c("gene")) %>%
  pivot_longer(names_to = "terms",values_to = "presence",-gene) %>%
  pivot_wider(names_from = "gene",values_from = "presence") %>%
  column_to_rownames("terms") %>%
  as.matrix()

#
kegg_t2m_multiGSEA[1:10,1:10]
kegg_t2g_multiGSEA[1:10,1:10]

# save the matrix as an object
saveRDS(kegg_t2m_multiGSEA,"../../out/object/kegg_t2m_multiGSEA.rds")
saveRDS(kegg_t2g_multiGSEA,"../../out/object/kegg_t2g_multiGSEA.rds")

# check validity of the matrix --------------------------------------------
# check that the matrix is correct
test_matrix_metabolite <- kegg_t2m_multiGSEA %>% 
  data.frame() %>% 
  rownames_to_column("term") %>% 
  pivot_longer(names_to = "metabolite",values_to = "presence",-term) %>% 
  # pull(term) %>% unique()
  filter(term=="(KEGG) Citrate cycle (TCA cycle)") %>% 
  filter(presence == 1) %>% 
  pull(metabolite) %>% 
  unique()
#
test_matrix_gene <- kegg_t2g_multiGSEA %>% 
  data.frame() %>% 
  rownames_to_column("term") %>% 
  pivot_longer(names_to = "gene",values_to = "presence",-term) %>% 
  # pull(term) %>% unique()
  filter(term=="(KEGG) Citrate cycle (TCA cycle)") %>% 
  filter(presence == 1) %>% 
  pull(gene) %>% 
  unique()

# confrime they are the same genes
test_orig_metabol <- LUT_metabolitc_fix_shortlist[LUT_metabolitc_fix_shortlist$hmdb %in% (pathways_kegg_multiGSEA$metabolome[["(KEGG) Citrate cycle (TCA cycle)"]] %>% unique()),] %>% pull(metabolite)
test_orig_gene <- pathways_kegg_multiGSEA$transcriptome[["(KEGG) Citrate cycle (TCA cycle)"]] %>% unique()

#
sum(!test_matrix_metabolite %in% test_orig_metabol)
sum(!test_matrix_gene %in% test_orig_gene)
