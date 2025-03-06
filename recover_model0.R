rm(list = ls()) # erase previous Global Environment, import manually the file

# libraries
if (!require(foreach)) install.packages("foreach"); library(foreach)
if (!require(doParallel)) install.packages("doParallel"); library(doParallel)

# import data and functions
source("functions.R")

# fitted files names
fitFiles <- paste0("data/",list.files("data/",pattern = "mod0_"))

# participants with recovery file
recovery <- fitFiles[grepl("recovery",fitFiles)]
# remove "recovery" and add .RData to get the outersect with fitFiles
if (length(recovery)!=0) {
  recovery <- paste0(substr(recovery,1,nchar(recovery)-15),".RData")
}

# remove "recovery" from fit files
fitFiles <- fitFiles[!grepl("recovery",fitFiles)]
# is there any extra IDs with no "recovery file?
extraIDs <- outersect(fitFiles,recovery)
if (length(extraIDs)!=0) {
  fitFiles <- extraIDs
}

# number of simulations per participants' parameters
nSimEachSubj <- 10

# start time
start_time <- Sys.time()

# Set up parallel processing
num_cores <- detectCores() - 4 # Adjust as needed
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Parallelize the loop using foreach
foreach(i = 1:length(fitFiles)) %dopar% {
# for (i in 1:length(fitFiles)) {
  
  # load files from one participant
  load(fitFiles[i])
  
  # add response variance to bestMod 
  fit$bestMod$param$var <- fit$parameters$resp_var
  
  for (j in 1:nSimEachSubj) {
    # simulate 
    sim <- mod0(fit$bestMod$param, fit$bestMod$x)
    
    # # # # # Rosa's Model # # # # #
    bins <- c(2,2,2,2) #paper=c(13,13,13,13)
    
    # prepare list for fitting function
    toFit <- list(src_subject_id = fit$modCompar$src_subject_id, 
                  phenotype = fit$modCompar$phenotype,
                  site = fit$modCompar$site, x = fit$bestMod$x, bins = bins)
    toFit$x$trial_structure$prediction <- sim$sim_resp
    
    # # # # # fit # # # # #
    refit <- fit_mod0(toFit)
    
    # combine 
    if (j == 1) {
      recovery <- data.frame(nSim=j,refit$parameters,refit$modCompar)
    } else {
      recovery <- rbind(recovery,data.frame(nSim=j,refit$parameters,refit$modCompar))
    }
  } # end j
  save(recovery, file = paste0("data/mod0_",refit$modCompar$src_subject_id,"_recovery.RData"))
}

# Clean up parallel processing
stopCluster(cl)
registerDoSEQ() # Revert to single-threaded processing

# end time
end_time <- Sys.time()

# how long was it?
end_time - start_time

# c(13,13,13,13); nSim = 10; num_cores = 8; Time difference of 9.581508 hours; 0.999