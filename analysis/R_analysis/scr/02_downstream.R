# python version ----------------------------------------------------------
Sys.setenv(RETICULATE_PYTHON = "/home/edo/miniconda3/envs/env_MOFA/bin/python")
library(reticulate)
reticulate::use_condaenv(condaenv = "env_MOFA")

py_config()

# Introduction ------------------------------------------------------------
# In the MOFA2 R package we provide a wide range of downstream analysis to visualise and interpret the model output. Here we provide a brief description of the main functionalities. This vignette is made of simulated data and we do not highlight biologically relevant results.

# 2Load libraries ---------------------------------------------------------
library(ggplot2)
library(MOFA2)

# 3Load trained model -----------------------------------------------------
# filepath <- system.file("extdata", "model.hdf5", package = "MOFA2")
filepath <- "out/object/model.hdf5"
print(filepath)

model <- load_model(filepath)

# 3.1Overview of data -----------------------------------------------------
# The function plot_data_overview can be used to obtain an overview of the input data. It shows how many views (rows) and how many groups (columns) exist, what are their corresponding dimensionalities and how many missing information they have (grey bars).
plot_data_overview(model)

# 4Add metadata to the model ----------------------------------------------
# The metadata is stored as a data.frame object in model@samples_metadata, and it requires at least the column sample. The column group is required only if you are doing multi-group inference.
# The number of rows must match the total number of samples in the model (sum(model@dimensions$N)).

# Let’s add some artifical metadata…

Nsamples <- sum(model@dimensions$N)

sample_metadata <- data.frame(
  sample = samples_names(model)[[1]],
  condition = sample(c("A","B"), size = Nsamples, replace = T),
  age = sample(1:100, size = Nsamples, replace = T)
)

samples_metadata(model) <- sample_metadata
head(model@samples_metadata, n=3)

# 5Variance decomposition -------------------------------------------------
# The first step in the MOFA analysis is to quantify the amount of variance explained (R2) by each factor in each data modality.
# The variance explained estimates are stored in the hdf5 file and loaded in model@cache[["variance_explained"]]:
  
# Total variance explained per view and group
head(model@cache$variance_explained$r2_total[[1]]) # group 1

# Variance explained for every factor in per view and group
head(model@cache$variance_explained$r2_per_factor[[1]]) # group 1

# Variance explained estimates can be plotted using plot_variance_explained(model, ...). Options:
# factors: character vector with a factor name(s), or numeric vector with the index(es) of the factor(s). Default is “all”.
# x: character specifying the dimension for the x-axis (“view”, “factor”, or “group”).
# y: character specifying the dimension for the y-axis (“view”, “factor”, or “group”).
# split_by: character specifying the dimension to be faceted (“view”, “factor”, or “group”).
# plot_total: logical value to indicate if to plot the total variance explained (for the variable in the x-axis)
plot_variance_explained(model, x="view", y="factor")

plot_variance_explained(model, x="group", y="factor", plot_total = T)[[2]]

# 5.0.1Visualisation of samples in the latent space -----------------------
# Each MOFA factor captures a different dimension of heterogeneity in the data. Mathematically, each factor ordinates cells along a one-dimensional axis centered at zero. Samples with different signs manifest opposite phenotypes along the inferred axis of variation, with higher absolute value indicating a stronger effect. Note that the interpretation of factors is analogous to the interpretation of the principal components in PCA.

# 5.1Visualisation of single factors --------------------------------------
# Factors can be plotted using plot_factor (for beeswarm plots of individual factors) or plot_factors (for scatter plots of factor combinations)
plot_factor(model, 
            factor = 1:3,
            color_by = "age",
            shape_by = "condition"
)


# Adding more options
p <- plot_factor(model, 
                 factors = c(1,2,3),
                 color_by = "condition",
                 dot_size = 3,        # change dot size
                 dodge = T,           # dodge points with different colors
                 legend = F,          # remove legend
                 add_violin = T,      # add violin plots,
                 violin_alpha = 0.25  # transparency of violin plots
)

# The output of plot_factor is a ggplot2 object that we can edit
p <- p + 
  scale_color_manual(values=c("A"="black", "B"="red")) +
  scale_fill_manual(values=c("A"="black", "B"="red"))

print(p)

# 5.2Visualisation of combinations of factors -----------------------------
# Scatter plots
plot_factors(model, 
             factors = 1:3,
             color_by = "condition"
)

# 5.3Visualisation of feature weights -------------------------------------
# The weights provide a score for how strong each feature relates to each factor. Features with no association with the factor have values close to zero, while features with strong association with the factor have large absolute values. The sign of the weight indicates the direction of the effect: a positive weight indicates that the feature has higher levels in the cells with positive factor values, and vice versa.

# Weights can be plotted using plot_weights (beeswarm plots) or plot_top_weights (scatter plots)
plot_weights(model,
             view = "view_0",
             factor = 1,
             nfeatures = 10,     # Number of features to highlight
             scale = T,          # Scale weights from -1 to 1
             abs = F             # Take the absolute value?
)


plot_top_weights(model,
                 view = "view_0",
                 factor = 1,
                 nfeatures = 10
)

# 6Visualisation of patterns in the input data ----------------------------
# Instead of looking at weights, it is useful to observe the coordinated heterogeneity that MOFA captures in the original data. This can be done using the plot_data_heatmap and plot_data_scatter function.

# 6.1Heatmaps -------------------------------------------------------------
# Heatmap of observations. Top features are selected by its weight in the selected factor. By default, samples are ordered according to their corresponding factor value.
plot_data_heatmap(model,
                  view = "view_1",         # view of interest
                  factor = 1,             # factor of interest
                  features = 20,          # number of features to plot (they are selected by weight)
                  
                  # extra arguments that are passed to the `pheatmap` function
                  cluster_rows = TRUE, cluster_cols = FALSE,
                  show_rownames = TRUE, show_colnames = FALSE
)

# 6.2Scatter plots --------------------------------------------------------
# Scatter plots of observations vs factor values. It is useful to add a linear regression estimate to visualise if the relationship between (top) features and factor values is linear.
plot_data_scatter(model,
                  view = "view_1",         # view of interest
                  factor = 1,             # factor of interest
                  features = 5,           # number of features to plot (they are selected by weight)
                  add_lm = TRUE,          # add linear regression
                  color_by = "condition"
)

# 6.3Non-linear dimensionality reduction ----------------------------------
# The MOFA factors are linear (as in Principal Component analysis). Nevertheless, the MOFA factors can be used as input to other methods that learn compact nonlinear manifolds (t-SNE or UMAP).

# Run UMAP and t-SNE

set.seed(42)
# model <- run_umap(model)
model <- run_tsne(model)
# Plot non-linear dimensionality reduction

plot_dimred(model,
            method = "TSNE",  # method can be either "TSNE" or "UMAP"
            color_by = "condition"
)

# 7Other functionalities --------------------------------------------------

# 7.1Renaming dimensions --------------------------------------------------
# The user can rename the dimensions of the model
views_names(model) <- c("Transcriptomics", "Proteomics")
factors_names(model) <- paste("Factor", 1:model@dimensions[["K"]], sep=" ")
views_names(model)

# 7.2Extracting data for downstream analysis ------------------------------
# The user can extract the feature weights, the data and the factors to generate their own plots.
# Extract factors

# "factors" is a list of matrices, one matrix per group with dimensions (nsamples, nfactors)
factors <- get_factors(model, factors = "all")
lapply(factors,dim)

# Extract weights
# "weights" is a list of matrices, one matrix per view with dimensions (nfeatures, nfactors)
weights <- get_weights(model, views = "all", factors = "all")
lapply(weights,dim)

# Extract data

# "data" is a nested list of matrices, one matrix per view and group with dimensions (nfeatures, nsamples)
data <- get_data(model)
lapply(data, function(x) lapply(x, dim))[[1]]

# For convenience, the user can extract the data in long data.frame format:
factors <- get_factors(model, as.data.frame = T)
head(factors, n=3)

weights <- get_weights(model, as.data.frame = T)
head(weights, n=3)

data <- get_data(model, as.data.frame = T)
head(data, n=3)
