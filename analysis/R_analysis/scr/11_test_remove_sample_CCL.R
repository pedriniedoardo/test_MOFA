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
library(patchwork)

# Data is stored as a list of matrices. Features are stored in the rows and samples in the columns
# utils::data("CLL_data")

# save the object in local
# saveRDS(CLL_data,file = "../../out/object/CLL_data.rds")
CLL_data <- readRDS("../../out/object/CLL_data.rds")

lapply(CLL_data,dim)

str(CLL_data)

CLL_data$mRNA[1:10,1:10]

# run the test ------------------------------------------------------------
# in this test I want to generate different datasets with a progressive reduction of samples in the RNA layer to see the effect on the whole analysis. I am planning to track the relative proportion of variance explained by each layer as comparative metrics.

# # identify the samples with a RNA profiling
# id_samplesRNA <- apply(is.na(CLL_data$mRNA),MARGIN = 2,FUN = function(x){
#   sum(x)
# })
# 
# # how many are with data
# sum(id_samplesRNA==0)
# 
# # define the number of sample to keep per condition. round to integer
# id_sample_keep <- c(
#   round(1*sum(id_samplesRNA==0),digits = 0),
#   round(3/4*sum(id_samplesRNA==0),digits = 0),
#   round(1/2*sum(id_samplesRNA==0),digits = 0),
#   round(1/4*sum(id_samplesRNA==0),digits = 0),
#   round(1/10*sum(id_samplesRNA==0),digits = 0)) %>% 
#   rep(c(1,3,3,3,3)) %>% 
#   setNames(c("all",
#              "test75_01",
#              "test75_02",
#              "test75_03",
#              "test50_01",
#              "test50_02",
#              "test50_03",
#              "test25_01",
#              "test25_02",
#              "test25_03",
#              "test10_01",
#              "test10_02",
#              "test10_03"))
# 
# # among the sample with an RNA profile randomly select sets of sample with a decreased number of samples
# # set a seed to reporduce the random sampling
# set.seed(21)
# 
# list_CLL_test <- lapply(id_sample_keep, function(x){
#   test <- CLL_data
#   # the one below is working because all the sample are the first 136. if this is not the case we need to implement differently
#   id_keep <- sample(x = 1:sum(id_samplesRNA==0),size = x,replace = F) %>% 
#     sort()
#   id_remove <- (1:sum(id_samplesRNA==0))[!(1:sum(id_samplesRNA==0) %in% id_keep)]
#   
#   # make the id_remove all NA in the RNA layer only
#   test$mRNA[,id_remove] <- NA
#   
#   return(test)
# })
# # save the list of sample input datasets
# saveRDS(list_CLL_test,file = "../../out/object/list_CLL_test_remove_sample.rds")

# read in the test dataset generated
list_CLL_test <- readRDS(file = "../../out/object/list_CLL_test_remove_sample.rds")

# confirm the addition modification of the samples
lapply(list_CLL_test,function(x){
  data.frame(NA_sample = apply(is.na(x$mRNA),MARGIN = 2,any) %>% 
    sum()) %>% 
    mutate(non_NA_sample = dim(x$mRNA)[2] - NA_sample)
    
})

# read in the metadata
CLL_metadata <- read_tsv("../../out/table/CLL_metadata.tsv")

# 3Create the MOFA obejct and train the model -----------------------------
# do the regular processing in an iterative loop
pmap(list(list_CLL_test,names(list_CLL_test)),function(CLL_test,id_name){
  # Create the MOFA object
  MOFAobject <- create_mofa(CLL_test)
  MOFAobject
  
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
  
  MOFAobject <- run_mofa(MOFAobject, outfile=paste0("../../out/object/test_remove_sample/MOFA2_CLL_test_remove_sample_",id_name,".hdf5"))
  saveRDS(MOFAobject,paste0("../../out/object/test_remove_sample/MOFA2_CLL_test_remove_sample_",id_name,".rds"))
})


# read in the results -----------------------------------------------------
# loop read the results
sample_id <- c("all",
               "test75_01",
               "test75_02",
               "test75_03",
               "test50_01",
               "test50_02",
               "test50_03",
               "test25_01",
               "test25_02",
               "test25_03",
               "test10_01",
               "test10_02",
               "test10_03")

# define the sample
id_MOFA <- paste0("../../out/object/test_remove_sample/MOFA2_CLL_test_remove_sample_",sample_id,".rds")

# read in the results
list_MOFA_result <- lapply(id_MOFA, function(x){
  test <- readRDS(x)
  return(test)
}) %>% 
  setNames(sample_id)

# plot the sample organization
list_plot_MOFA <- pmap(list(list_MOFA_result,names(list_MOFA_result)), function(x,name){
  plot_data_overview(x)+ggtitle(name)
})

wrap_plots(list_plot_MOFA[1:10])
ggsave("../../out/image/plot_data_overview_test_remove_sample.pdf",width = 12,height = 9)

# extract the variance explained
list_data_variance <- lapply(list_MOFA_result, function(x){
  test <- plot_variance_explained(x, max_r2=15)
  return(test$data)
})

df_var <- list_data_variance %>% 
  bind_rows(.id = "dataset") %>% 
  mutate(condition = str_extract(dataset,pattern = "75|50|25|10")) %>%
  mutate(condition = str_pad(condition,width = 3,pad = "0")) %>% 
  mutate(condition = case_when(is.na(condition)~"100",
                               T~condition)) %>% 
  mutate(condition_fix = paste0("t_",condition)) %>% 
  mutate(replicate = str_extract(dataset,pattern = "01|02|03")) %>% 
  mutate(replicate = case_when(is.na(replicate)~"01",
                               T~replicate))

# plot the variance across dataset
# df_var %>% 
#   group_by(condition_fix,factor,view) %>% 
#   summarise(avg_value = mean(value)) %>% 
#   ggplot() + 
#   # geom_point() +
#   geom_line(aes(x = condition_fix, y = avg_value, group = view,col=view))+facet_wrap(~factor)

df_var %>% 
  filter(factor %in% c("Factor1","Factor2","Factor3","Factor4","Factor5","Factor6","Factor7","Factor8","Factor9")) %>% 
  filter(! condition_fix %in% c("t_010")) %>% 
  group_by(condition_fix,factor,view) %>% 
  mutate(avg_value = mean(value)) %>% 
  ggplot(aes(x = condition_fix, y = value, color = view)) + 
  geom_point(alpha = 0.5) +
  facet_wrap(~factor)+
  theme_bw() +
  theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 90)) +
  geom_line(aes(x = condition_fix, y = avg_value, group = view,col=view))
ggsave("../../out/image/plotVar_test_remove_sample.pdf",width = 10,height = 9)

# -------------------------------------------------------------------------
# extract the variance explained
list_data_tot_variance <- lapply(list_MOFA_result, function(x){
  test <- plot_variance_explained(x, plot_total = T)[[2]]
  return(test$data)
})

df_tot_var <- list_data_tot_variance %>% 
  bind_rows(.id = "dataset") %>% 
  mutate(condition = str_extract(dataset,pattern = "75|50|25|10")) %>%
  mutate(condition = str_pad(condition,width = 3,pad = "0")) %>% 
  mutate(condition = case_when(is.na(condition)~"100",
                               T~condition)) %>% 
  mutate(condition_fix = paste0("t_",condition)) %>% 
  mutate(replicate = str_extract(dataset,pattern = "01|02|03")) %>% 
  mutate(replicate = case_when(is.na(replicate)~"01",
                               T~replicate))

df_tot_var %>% 
  filter(! condition_fix %in% c("t_010")) %>% 
  group_by(condition_fix,view) %>% 
  mutate(avg_value = mean(R2)) %>% 
  ggplot(aes(x = condition_fix, y = R2, color = view)) + 
  geom_point(alpha = 0.5) +
  theme_bw() +
  theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 90)) +
  geom_line(aes(x = condition_fix, y = avg_value, group = view,col=view))
ggsave("../../out/image/plotTotVar_test_remove_sample.pdf",width = 5,height = 4)
