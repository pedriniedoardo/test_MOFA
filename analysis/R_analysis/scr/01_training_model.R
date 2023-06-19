# python version ----------------------------------------------------------
Sys.setenv(RETICULATE_PYTHON = "/home/edo/miniconda3/envs/env_MOFA/bin/python")
library(reticulate)
reticulate::use_condaenv(condaenv = "env_MOFA")

py_config()

# 1Load libraries ---------------------------------------------------------
library(data.table)
library(MOFA2)

# This vignette contains a detailed tutorial on how to train a MOFA model using R. A concise template script can be found here


# 2What is MOFA? ----------------------------------------------------------
# MOFA is an unsupervised statistical framework for the integration of multi-omic data sets.
# Intuitively, MOFA can be viewed as a versatile and statistically rigorous generalization of principal component analysis (PCA) to multi-omics data. Given several data matrices with measurements of multiple -omics data types on the same or on overlapping sets of samples, MOFA infers an interpretable low-dimensional representation in terms of a few latent factors that (hopefully) captures the relevant signal in the input data.
# Effectively MOFA disentangles the sources of variation in the data, identifying the factors that are shared across multiple data modalities from the factors that drive variability in a single data modality.

# 2.1Is MOFA appropriate for my data? -------------------------------------
# MOFA (and factor analysis models in general) are useful to uncover variation in complex data sets that are expected to contain multiple sources of heterogeneity. This requires a relatively large sample size (at least ~15 samples). In addition, MOFA needs the multi-modal measurements to be derived from the same samples. It is fine if you have samples that are missing some data modality, but there has to be a significant degree of matched measurements.

# 3Load data --------------------------------------------------------------
# To create a MOFA object you need to specify four dimensions: samples, features, view(s) and group(s). MOFA objects can be created from a wide range of input formats, including:

# a list of matrices: this is recommended for relatively simple data.
# a long data.frame: this is recommended for complex data sets with multiple views and/or groups.
# MultiAssayExperiment: to connect with Bioconductor objects.
# Seurat: only for single-cell genomics users.

# 3.1List of matrices -----------------------------------------------------
# A list of matrices, where each entry corresponds to one view. Samples are stored in columns and features in rows.

# Let’s simulate some data to start with

data <- make_example_data(
  n_views = 2, 
  n_samples = 200, 
  n_features = 1000, 
  n_factors = 10
)[[1]]

lapply(data,dim)

# Create the MOFA object:
MOFAobject <- create_mofa(data)

# In case you are using the multi-group functionality, the groups can be specified using a vector with the group ID for each sample.

# Please keep in mind that this is a rather advanced option that we disencourage for MOFA beginners. For more details on how the multi-group inference works, read the FAQ section.

N <- ncol(data[[1]])
groups <- c(rep("A",N/2), rep("B",N/2))

MOFAobject <- create_mofa(data, groups=groups)

# 3.2Long data.frame ------------------------------------------------------
# A long data.frame with columns sample, feature, view, group (optional), value.
# In my opinion this is the best format for complex data sets with multiple omics and potentially multiple groups of data. Also, there is no need to add rows that correspond to missing data:
dt <- fread("ftp://ftp.ebi.ac.uk/pub/databases/mofa/getting_started/data.txt.gz")
head(dt)

# Let's ignore the grouping information to start
dt[,group:=NULL]

# Create the MOFA object
MOFAobject <- create_mofa(dt)
print(MOFAobject)

# 3.3Visualise the structure of the data ----------------------------------
plot_data_overview(MOFAobject)

# 4Define options ---------------------------------------------------------

# 4.1Define data options --------------------------------------------------
# scale_groups: if groups have different ranges/variances, it is good practice to scale each group to unit variance. Default is FALSE
# scale_views: if views have different ranges/variances, it is good practice to scale each view to unit variance. Default is FALSE
data_opts <- get_default_data_options(MOFAobject)
head(data_opts)

# 4.2Define model options -------------------------------------------------
# num_factors: number of factors
# likelihoods: likelihood per view (options are “gaussian”, “poisson”, “bernoulli”). By default they are learnt automatically. We advise users to use “gaussian” whenever possible!
# spikeslab_factors: use spike-slab sparsity prior in the factors? default is FALSE.
# spikeslab_weights: use spike-slab sparsity prior in the weights? default is TRUE.
# ard_factors: use ARD prior in the factors? Default is TRUE if using multiple groups.
# ard_weights: use ARD prior in the weights? Default is TRUE if using multiple views.
# Only change the default model options if you are familiar with the underlying mathematical model!
model_opts <- get_default_model_options(MOFAobject)
head(model_opts)

# 4.3Define train options -------------------------------------------------
# maxiter: number of iterations. Default is 1000.
# convergence_mode: “fast”, “medium”, “slow”. For exploration, the fast mode is good enough.
# startELBO: initial iteration to compute the ELBO (the objective function used to assess convergence).
# freqELBO: frequency of computations of the ELBO.
# gpu_mode: use GPU mode? (needs cupy installed and a functional GPU).
# stochastic: use stochastic inference? (default is FALSE).
# verbose: verbose mode?
# seed: random seed
train_opts <- get_default_training_options(MOFAobject)
head(train_opts)

# 5Build and train the MOFA object ----------------------------------------
# Prepare the MOFA object
MOFAobject <- prepare_mofa(
  object = MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

# Train the MOFA model
outfile <- file.path(getwd(),"../../out/object/model.hdf5")
MOFAobject.trained <- run_mofa(MOFAobject, outfile)

# 6Downstream analysis ----------------------------------------------------
# This finishes the tutorial on how to train a MOFA object from R. To continue with the downstream analysis, follow this tutorial

# edo recommended to input the z score do the scaled data