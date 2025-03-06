rm(list = ls()) # erase previous Global Environment, import manually the file

# libraries
if (!require(foreach)) install.packages("foreach"); library(foreach)
if (!require(doParallel)) install.packages("doParallel"); library(doParallel)

# import data and functions
db <- read.csv("data/data_modelling_v2.csv")
source("functions.R")

# in order to not refit the same again
model_fitted <- list.files("data/",pattern = ".RData")
# remove recovery files
model_fitted <- model_fitted[!grepl("recovery",model_fitted)]
# extract only ID
model_fitted <- substr(model_fitted,6,nchar(model_fitted)-6)
# get the outersect to identify etra IDs
extraIDs <- outersect(model_fitted,unique(db$src_subject_id))

# vector of subjects
if (length(extraIDs) == 0) {
  src_subject_id <- unique(db$src_subject_id) #"M042322"
} else {
  src_subject_id <- extraIDs
}

# start time
start_time <- Sys.time()

# Set up parallel processing
num_cores <- detectCores() - 6 # Adjust as needed
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Parallelize the loop using foreach
foreach(i = 1:length(src_subject_id)) %dopar% {
# for (i in 1:length(src_subject_id)) {message(i)
  
  # # # # # get one subject # # # # #
  oneSubj <- db[db$src_subject_id == src_subject_id[i],]
  
  # # # # # training protocol # # # # #
  x <- list()
  # inputs
  x$inp <- oneSubj[,grepl('in.',colnames(oneSubj)) & nchar(colnames(oneSubj)) < 6]
  # outputs
  x$out <- ifelse(matrix(oneSubj$outcome)==0,-1,1)
  # trial structure (trial types, blocks, phases)
  x$trial_structure <- oneSubj[,c("test_part","begin_timestamp","trial_type",
                                  "change_phase","trial_block","trial_number",
                                  "prediction")]
  x$trial_structure$change_phase[1] <- 0
  
  
  # # # # # Rosa's Model # # # # #
  bins <- c(2,2,2,2) #paper=c(23,23,23,23)
  
  # prepare list for fitting function
  toFit <- list(src_subject_id = src_subject_id[i], 
                phenotype = unique(oneSubj$phenotype),
                site = unique(oneSubj$site), x = x, bins = bins)
  
  # # # # # fit # # # # #
  fit <- fit_mod0(toFit)
  
  # # # # # save # # # # #
  save(fit, file = paste0("data/mod0_", src_subject_id[i], ".RData"))
}

# Clean up parallel processing
stopCluster(cl)
registerDoSEQ() # Revert to single-threaded processing

# end time
end_time <- Sys.time()

# how long was it?
end_time - start_time

# c(23,23,23,23); num_cores = 8; Time difference of ~7.515847 hours; 0.999