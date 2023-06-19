# python version ----------------------------------------------------------
Sys.setenv(RETICULATE_PYTHON = "/home/edo/miniconda3/envs/env_MOFA/bin/python")
library(reticulate)
reticulate::use_condaenv(condaenv = "env_MOFA")

py_config()

# 1Introduction -----------------------------------------------------------
# This vignette show how to use MOFA+ on the bulk multi-omics data set that was used in the first publication of MOFA and the original vignette of the MOFA package.

# Briefly, the data consists of four omics including DNA methylation, RNA-seq, somatic mutations and drug response data from blood for N=200 patients with Chronic Lymphocytic Leukemia (CLL). The data set is explained in detail in this article and is publicly available here

# 2Load libraries and data ------------------------------------------------
library(MOFA2)
library(MOFAdata)
library(data.table)
library(ggplot2)
library(tidyverse)
# Data is stored as a list of matrices. Features are stored in the rows and samples in the columns
# utils::data("CLL_data")

# save the object in local
# saveRDS(CLL_data,file = "../../out/object/CLL_data.rds")
CLL_data <- readRDS("../../out/object/CLL_data.rds")

lapply(CLL_data,dim)

str(CLL_data)

CLL_data$mRNA[1:10,1:10]

# run the test ------------------------------------------------------------
# in this test I want to see the effect of adding same sample (all NA to all the layers) to see the effect of having NA samples across all the layers
CLL_test <- lapply(CLL_data, function(x){
  test <- x %>% 
    data.frame() %>% 
    mutate(sample_test01 = NA,
           sample_test02 = NA,
           sample_test03 = NA,
           sample_test04 = NA) %>% 
    as.matrix()
  return(test)
})
# confirm the addition of the samples
lapply(CLL_test,dim)

# create a face metadata for the new test samples
test_meta <- data.frame(sample = c("sample_test01",
                      "sample_test02",
                      "sample_test03",
                      "sample_test04"),
           Gender = c("f","f","m","m"),
           age = NA,
           TTT = NA,
           TTD = NA,
           treatedAfter = NA,
           died = NA,
           IGHV = NA,
           trisomy12 = NA)

# read in the metadata
CLL_metadata <- read_tsv("../../out/table/CLL_metadata.tsv") %>% 
  # add the three samples with a fake Gender but with NA for all the other metadata
  bind_rows(test_meta)

# 3Create the MOFA obejct and train the model -----------------------------
# Create the MOFA object
MOFAobject <- create_mofa(CLL_test)
MOFAobject

# 3.1Plot data overview ---------------------------------------------------
# Visualise the number of views (rows) and the number of groups (columns) exist, what are their corresponding dimensionalities and how many missing information they have (grey bars).
plot_data_overview(MOFAobject)
ggsave("../../out/image/plot_data_overview_addNA.pdf",width = 5,height = 5)

# 3.2Define MOFA options --------------------------------------------------
data_opts <- get_default_data_options(MOFAobject)
data_opts

# 3.2.2Model options ------------------------------------------------------
model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 15

model_opts

# 3.2.3Training options ---------------------------------------------------
train_opts <- get_default_training_options(MOFAobject)
train_opts$convergence_mode <- "slow"
train_opts$seed <- 42

train_opts

# 3.3Train the MOFA model -------------------------------------------------
# Prepare the MOFA object
MOFAobject <- prepare_mofa(MOFAobject,
                           data_options = data_opts,
                           model_options = model_opts,
                           training_options = train_opts
)

MOFAobject <- run_mofa(MOFAobject, outfile="../../out/object/MOFA2_CLL_test_addNA.hdf5")
saveRDS(MOFAobject,"../../out/object/MOFA2_CLL_test_addNA.rds")
MOFAobject <- readRDS("../../out/object/MOFA2_CLL_test_addNA.rds")

# 4Overview of the trained MOFA model -------------------------------------

# 4.1Slots ----------------------------------------------------------------
# confirm the size
dim(MOFAobject@data$Drugs$group1)

# 4.2Add sample metadata to the model -------------------------------------
# Add sample metadata to the model
samples_metadata(MOFAobject) <- CLL_metadata
samples_metadata(MOFAobject)

# 4.3Correlation between factors ------------------------------------------
plot_factor_cor(MOFAobject)

# 4.4Plot variance decomposition ------------------------------------------

# 4.4.1Variance decomposition by Factor -----------------------------------
# The most important insight that MOFA generates is the variance decomposition analysis. This plot shows the percentage of variance explained by each factor across each data modality (and group, if provided). It summarises the sources of variation from a complex heterogeneous data set in a single figure.
plot_variance_explained(MOFAobject, max_r2=15)
test <- plot_variance_explained(MOFAobject, max_r2=15)
test$data

# 4.4.2Total variance explained per view ----------------------------------
plot_variance_explained(MOFAobject, plot_total = T)[[2]]
ggsave("../../out/image/plot_variance_explained_addNA.pdf",width = 5,height = 5)

# 5Characterisation of Factor 1 -------------------------------------------

# 5.1Association analysis -------------------------------------------------
# Let’s test the association between MOFA Factors and Gender, survival outcome (dead vs alive) and age:
pdf("../../out/image/correlate_factors_with_covariates_addNA.pdf",width = 5,height = 5)
correlate_factors_with_covariates(MOFAobject,
                                  covariates = c("Gender","died","age"),
                                  plot="log_pval"
)
dev.off()

# 5.2Plot factor values ---------------------------------------------------
plot_factor(MOFAobject, 
            factors = 1, 
            color_by = "Factor1"
)

# 5.3Plot feature weights -------------------------------------------------

# 5.3.1Plot feature weights for somatic mutations -------------------------
plot_weights(MOFAobject,
             view = "Mutations",
             factor = 1,
             nfeatures = 10,     # Top number of features to highlight
             scale = T           # Scale weights from -1 to 1
)

# An alternative visualistion to the full distribution of weights is to do a line plot that displays only the top features with the corresponding weight sign on the right:
plot_top_weights(MOFAobject,
                 view = "Mutations",
                 factor = 1,
                 nfeatures = 10,     # Top number of features to highlight
                 scale = T           # Scale weights from -1 to 1
)
ggsave("../../out/image/plot_weightsF1_Mut_addNA.pdf",width = 5,height = 5)

# IGHV has a positve weight. This means that samples with positive Factor 1 values have IGHV mutation whereas samples with negative Factor 1 values do not have the IGHV mutation. To confirm this, let’s plot the Factor values and colour the IGHV mutation status.
plot_factor(MOFAobject, 
            factors = 1, 
            color_by = "IGHV",
            add_violin = TRUE,
            dodge = TRUE
)

# 5.3.2Plot gene weights for mRNA expression ------------------------------
# From the variance explained plot we know that Factor 1 drives variation across all data modalities. Let’s visualise the mRNA expression changes that are associated with Factor 1:
plot_weights(MOFAobject,
             view = "mRNA",
             factor = 1,
             nfeatures = 10
)

plot_top_weights(MOFAobject,
             view = "mRNA",
             factor = 1,
             nfeatures = 10
)
ggsave("../../out/image/plot_weightsF1_RNA_addNA.pdf",width = 5,height = 5)

# 5.3.3Plot molecular signatures in the input data ------------------------
# In this case we have a large amount of genes that have large positive and negative weights. Genes with large positive values will be more expressed in the samples with IGHV mutation, whereas genes with large negative values will be more expressed in the samples without the IGHV mutation. Let’s verify this. The function plot_data_scatter generates a scatterplot of Factor 1 values (x-axis) versus expression values (y-axis) for the top 4 genes with largest positive weight. Samples are coloured by IGHV status:

# This function generates a scatterplot of Factor 1 values (x-axis) versus expression values (y-axis) for the top 4 genes with largest negative weight. Samples are coloured by IGHV status:
plot_data_scatter(MOFAobject,
                  view = "mRNA",
                  factor = 1,
                  features = 4,
                  sign = "negative",
                  color_by = "IGHV"
) +
  labs(y="RNA expression")

# An alternative visualisation is to use a heatmap
plot_data_heatmap(MOFAobject,
                  view = "mRNA",
                  factor = 1,  
                  features = 25,
                  cluster_rows = FALSE, cluster_cols = FALSE,
                  show_rownames = TRUE, show_colnames = FALSE,
                  scale = "row"
)

# plot_data_heatmap has an interesting argument to “beautify” the heatmap: denoise = TRUE. Instead of plotting the (noisy) input data, we can plot the data reconstructed by the model, where noise has been removed:
plot_data_heatmap(MOFAobject,
                  view = "mRNA",
                  factor = 1,
                  features = 25,
                  denoise = TRUE,
                  cluster_rows = FALSE, cluster_cols = FALSE,
                  show_rownames = TRUE, show_colnames = FALSE,
                  scale = "row"
)

# 9Gene set enrichment analysis (GSEA) ------------------------------------
# In addition to exploring the individual weights for each factor, we can use enrichment analysis to look for signiificant associations of factors to genesets. Here, we use the Reactome genesets for illustrations, which is contained in the MOFAdata package. For more details on how the GSEA works we encourage the users to read the GSEA vignette

# 9.1Load Reactome gene set annotations. ----------------------------------
# Gene set annotations are provided as a binary membership matrix. Genes are stored in the rows, pathways are stored in the columns. A value of 1 indicates that gene j belongs to the pathway i.

utils::data(reactomeGS)

head(colnames(reactomeGS))
head(rownames(reactomeGS))

# 9.2Run enrichment analysis ----------------------------------------------
# GSEA on positive weights, with default options
res.positive <- run_enrichment(MOFAobject, 
                               feature.sets = reactomeGS, 
                               view = "mRNA",
                               sign = "positive"
)

# GSEA on negative weights, with default options
res.negative <- run_enrichment(MOFAobject, 
                               feature.sets = reactomeGS, 
                               view = "mRNA",
                               sign = "negative"
)

names(res.positive)

# 9.2.1Plot enrichment analysis results -----------------------------------
plot_enrichment_heatmap(res.positive)
plot_enrichment_heatmap(res.negative)

# Let’s plot the GSEA results for Factor 5. It seems that this Factor is capturing differences in the stress response of the blood cells.
pdf("../../out/image/plot_enrichment_res.pos_addNA.pdf",width = 10,height = 3)
plot_enrichment(res.positive, factor = 5, max.pathways = 15)
dev.off()
