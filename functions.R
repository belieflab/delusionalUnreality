# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # Models# # # # # # # # # # # # # # # # # # # # # ####
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# learning model 1 (Rosa's model)
mod1 <- function (param,x) {
  # parameters
  lr <- param$lr
  invT <- param$invT
  tau <- param$tau
  eta <- param$eta
  startV <- param$startV
  # cues
  input <- as.matrix(x$inp)
  nCues <- ncol(input)
  # outcomes or feedback
  output <- x$out
  nOut <- ncol(output)
  nTrials <- nrow(x$trial_structure)
  # phases
  phases <- x$trial_structure
  # weights (w) and change (delta) of w
  w <- delta <- array(0,dim=c(nTrials,nCues)); w[1,] <- startV
  # responses, prob, and softmax vectors
  pResp <- netInput <- resp <- xi <- PE <- rep(NA,nTrials)
  for (t in 1:nTrials) {
    if (phases$change_phase[t]==1) {
      w[t,] <- eta * w[t,]
    }
    # softmax input
    netInput[t] <- w[t,] %*% input[t,]
    # p(response)
    pResp[t] <- f_logistic(invT, netInput[t])
    # predictive response
    resp[t] <- rbinom(1, 1, pResp[t])
    # overall expectation with noisy-MAX function
    # xi[t] <- f_noisyMax(w[t,],input[t,],tau)
    if (sum(input[t,])==1) {
      xi[t] <- netInput[t]
    } else {
      xi[t] <- f_rosaAdditivity(w[t,],input[t,],tau)
    }
    # prediction error
    PE[t] <- as.vector(output[t] - xi[t])
    # learning
    delta[t,] <- lr * input[t,] * PE[t]
    # update weights (learning)
    if (t < nTrials) {
      w[t+1,] <- w[t,] + delta[t,]
    }
  } # end t
 
  # function output
  return(list(netInput = netInput, pResp = pResp, resp = resp, xi = xi, 
              PE = PE, delta = delta, w = w, x = x, param = param)) 
} # end

mod0 <- function (param,x) {
  # parameters
  lr <- param$lr
  gamma <- param$gamma
  eta <- param$eta
  lambda <- param$lambda
  var <- param$var
  # cues
  input <- as.matrix(x$inp)
  nCues <- 16#ncol(input)
  # outcomes or feedback
  output <- x$out
  nTrials <- nrow(x$trial_structure)
  # phases
  change_phase <- x$trial_structure$change_phase
  pred <- x$trial_structure$prediction
  # weights (w) and change (delta) of w
  w <- delta <- array(0,dim=c(nTrials,nCues)); w[1,] <- 0
  # responses, prob, and softmax vectors
  r_hat <- PE_p <- PE_s <- sim_resp <- rep(NA,nTrials)
  for (t in 1:nTrials) {
    if (change_phase[t]==1) {
      w[t,] <- eta * w[t,]
    }
    if (sum(input[t,])==1) {
      # prediction
      r_hat[t] <- w[t,] %*% input[t,]
      # prediction error
      PE_p[t] <- as.vector(output[t] - r_hat[t])
      # learning
      delta[t,] <- lr * input[t,] * PE_p[t]
    } else {
      weight <- w[t,]*input[t,]
      weight <- weight[input[t,]==1]
      w_p <- max(weight)
      w_s <- min(weight)
      if (w_p >= 0) lambda_t <- lambda else lambda_t <- 1
      # prediction
      r_hat[t] <- w_p + gamma * w_s
      # prediction error
      PE_p[t] <- as.vector(output[t] - r_hat[t])
      PE_s[t] <- as.vector(lambda_t * output[t] - w_s)
      # locate position from max (p) and min (s)
      max_min <- f_getMaxMinWeights(input[t,],weight)
      # learning
      delta[t,max_min] <- lr * c(PE_p[t],PE_s[t])
    }
    
    # produce simulated responses
    if (!is.na(pred[t]) && !is.null(var)) {
      sim_resp[t] <- rnorm(1,r_hat[t],var)
    }
    
    # update weights (learning)
    if (t < nTrials) {
      w[t+1,] <- w[t,] + delta[t,]
    }
  } # end t
  
  # function output
  return(list(r_hat = r_hat, PE_p = PE_p, PE_s = PE_s, delta = delta, w = w, 
              param = param, x = x, sim_resp = sim_resp)) 
} # end



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # Model Fit functions # # # # # # # # # # # # # # ####
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# learning model 1 (Rosa's model)
fit_mod1 <- function(toFit) {
  # start time
  start_time <- Sys.time()
  # number of parameters
  nPar <- 4
  # dichotomous actions {-1,1}
  resp <- ifelse(toFit$resp==1,1,0)
  x <- toFit$x
  nTrials <- nrow(x$out)
  bin1 <- toFit$bins[1]
  bin2 <- toFit$bins[2]
  bin3 <- toFit$bins[3]
  bin4 <- toFit$bins[4]
  # parametric space
  vLR <- seq(f_logit(0.01,0),f_logit(0.99,0), (f_logit(0.99,0)-f_logit(0.01,0))/(bin1-1))
  # vTau <- seq(log(0.1),log(100), (log(100)-log(0.1))/(bin2-1))
  vTau <- seq(f_logit(0.01,0),f_logit(0.99,0), (f_logit(0.99,0)-f_logit(0.01,0))/(bin2-1))
  vEta <- seq(f_logit(0.01,0),f_logit(0.99,0), (f_logit(0.99,0)-f_logit(0.01,0))/(bin3-1))
  vInvT <-  seq(log(0.1),log(100), (log(100)-log(0.1))/(bin4-1))
  # run parameters loop and extract softmax input
  V_inSm <- array(NA,dim = c(nTrials,bin1,bin2,bin3))
  for (i in 1:bin1) {
    for (j in 1:bin2) {
      for (k in 1:bin3) {
        param <- list(lr=f_logit(vLR[i],1), tau=f_logit(vTau[j],1),#exp(vTau[j]), 
                      eta=f_logit(vEta[k],1), invT=1, startV=0)
        temp <- mod1(param,x)
        V_inSm[,i,j,k] <- temp$netInput
      } # end k
    } # end j 
  } # end i
  # increase invT dimension
  V_inSm <- array(V_inSm, dim = c(nTrials,bin1,bin2,bin3,bin4))
  # create 3d invT space and get it in the correct order
  invTMat <- array(exp(vInvT), dim = c(bin4,nTrials,bin1,bin2,bin3))
  invTMat <- aperm(invTMat, c(2,3,4,5,1))
  # probability of response 
  pAs <- 1/(1+exp(-invTMat*V_inSm))
  # array with fitted responses
  respArr <- array(resp, dim = c(nTrials,bin1,bin2,bin3,bin4))
  
  # remove NAs and adjust due to missing responses
  pAs <- pAs[!is.na(resp),,,,]
  respArr <- respArr[!is.na(resp),,,,]
  
  # likelihoods
  LLd <- apply((pAs * respArr) + ((1-pAs) * (1-respArr)),c(2,3,4,5),prod)*10^10
  # obtain posterior probability with uniform prior
  post_prob <- (LLd*1)/sum(LLd)
  
  # learning rate marginal
  lr_marg <- apply(post_prob,1,sum)
  # learning rate ml
  lr_ml <- as.vector(f_logit(vLR[which.max(lr_marg)],1))
  # learning rate weighted mean
  lr_wm <- as.vector(f_logit(vLR %*% lr_marg, 1))
  # learning rate variance
  lr_var <- as.vector(f_logit( ((vLR - f_logit(lr_wm,0))^2) %*% lr_marg,1 ))
 
  # tau marginal
  tau_marg <- apply(post_prob,2,sum)
  # tau ml
  tau_ml <- as.vector(f_logit(vTau[which.max(tau_marg)],1))#as.vector(exp(vTau[which.max(tau_marg)]))
  # tau weighted mean
  tau_wm <- as.vector(f_logit(vTau %*% tau_marg, 1))#as.vector(exp(vTau %*% tau_marg))
  # tau variance
  tau_var <- as.vector(f_logit( ((vTau - f_logit(tau_wm,0))^2) %*% tau_marg,1 ))#as.vector(exp( ((vTau-log(tau_wm))^2) %*% tau_marg))
  
  # eta marginal
  eta_marg <- apply(post_prob,3,sum)
  # learning rate ml
  eta_ml <- as.vector(f_logit(vEta[which.max(eta_marg)],1))
  # learning rate weighted mean
  eta_wm <- as.vector(f_logit(vEta %*% eta_marg, 1))
  # learning rate variance
  eta_var <- as.vector(f_logit( ((vEta - f_logit(eta_wm,0))^2) %*% eta_marg,1 ))
  
  # inverse temperature marginal
  invT_marg <- apply(post_prob,4,sum)
  # inverse temperature ml
  invT_ml <- as.vector(exp(vInvT[which.max(invT_marg)]))
  # inverse temperature weighted mean
  invT_wm <- as.vector(exp(vInvT %*% invT_marg))
  # inverse temperature variance
  invT_var <- as.vector(exp( ((vInvT-log(invT_wm))^2) %*% invT_marg))
  
  # use best parameters to estimate maxLL, BIC, and hitRate
  bestMod <- mod1(list(lr=lr_wm, tau=tau_wm, eta=eta_wm, invT=invT_wm, 
                       startV = 0),x)
  
  # remove NAs and adjust due to missing responses
  nTrials <- sum(!is.na(resp))
  bestModPResp <- bestMod$pResp[!is.na(resp)]
  resp <- resp[!is.na(resp)]
  
  # maximum LLd
  maxLL <- prod((bestModPResp * resp) + ((1-bestModPResp) * (1-resp)))*10^10
  # Bayesian Information Criterion
  BIC <- (nPar * log(nTrials)) - (2 * log(maxLL)) 
  # hit rate
  hitRate <- mean((bestModPResp * resp) + (1-bestModPResp) * (1-resp))
  
  # set output list
  parameters <- data.frame(lr_ml=lr_ml,lr_wm=lr_wm,lr_var=lr_var,
                           tau_ml=tau_ml,tau_wm=tau_wm,tau_var=tau_var,
                           eta_ml=eta_ml,eta_wm=eta_wm,eta_var=eta_var,
                           invT_ml=invT_ml,invT_wm=invT_wm,invT_var=invT_var)
  # end time
  end_time <- Sys.time()
  fitDuration <- end_time - start_time
  
  # model specifications
  modCompar <- data.frame(src_subject_id = toFit$src_subject_id, 
                          phenotype = toFit$phenotype,
                          site = toFit$site,
                          maxLL,hitRate, BIC, fitDuration)
  # create output list
  output <- list(parameters = parameters, modCompar = modCompar, 
                 post_prob = post_prob, bestMod = bestMod)
  
  return(output)
} # end fitting function

fit_mod0 <- function(toFit) {
  # start time
  start_time <- Sys.time()
  # number of parameters
  nPar <- 4
  # responses are predictions bounded between -1 and 1
  pred <- toFit$x$trial_structure$prediction
  valid_trials <- !is.na(pred); valid_trials[1:7] <- F
  predClean <- pred[valid_trials]
  x <- toFit$x
  nTrials <- sum(valid_trials)
  bin1 <- toFit$bins[1]
  bin2 <- toFit$bins[2]
  bin3 <- toFit$bins[3]
  bin4 <- toFit$bins[4]
  # parametric space
  min <- 0.001
  max <- 0.999
  # vLR <- seq(f_logit(0.001,0),f_logit(0.999,0), (f_logit(0.999,0)-f_logit(0.001,0))/(bin1-1))
  # vGamma <- seq(f_logit(0.001,0),f_logit(0.999,0), (f_logit(0.999,0)-f_logit(0.001,0))/(bin2-1))
  # vEta <- seq(f_logit(0.001,0),f_logit(0.999,0), (f_logit(0.999,0)-f_logit(0.001,0))/(bin3-1))
  # vLambda <- seq(atanh(-0.999),atanh(0.999), (atanh(0.999)-atanh(-0.999))/(bin4-1))
  vLR <- seq(f_logit(min,0),f_logit(max,0), (f_logit(max,0)-f_logit(min,0))/(bin1-1))
  vGamma <- seq(f_logit(min,0),f_logit(max,0), (f_logit(max,0)-f_logit(min,0))/(bin2-1))
  vEta <- seq(f_logit(min,0),f_logit(max,0), (f_logit(max,0)-f_logit(min,0))/(bin3-1))
  vLambda <- seq(atanh(-max),atanh(max), (atanh(max)-atanh(-max))/(bin4-1))
 
  # run parameters loop and extract softmax input
  modelPredictions <- likelihood <- array(NA,dim = c(nTrials,bin1,bin2,bin3,bin4))
  for (i in 1:bin1) {
    for (j in 1:bin2) {
      for (k in 1:bin3) {
        for (l in 1:bin4) {
          param <- list(lr=f_logit(vLR[i],1), gamma=f_logit(vGamma[j],1), 
                        eta=f_logit(vEta[k],1), lambda=tanh(vLambda[l]),
                        var=NULL)
          temp <- mod0(param,x)
          modPred <- temp$r_hat[valid_trials]
          var <- sum((predClean-modPred)^2)/nTrials
          # p(D|M)
          likelihood[,i,j,k,l] <- dnorm(predClean,modPred,sqrt(var),log = F)
        } # end l
      } # end k
    } # end j 
  } # end i

  # sum log-likelihood density: p(D|M)
  prodlikelihood <- apply(likelihood,c(2,3,4,5),prod)*10^10
  # obtain posterior probability with uniform prior
  # p(M|D) = (p(D|M) * p(M)) / p(D)
  post_prob <- (prodlikelihood*1)/sum(prodlikelihood)
  
  # image(apply(post_prob,c(1,2),sum))
  
  # learning rate marginal
  lr_marg <- apply(post_prob,1,sum)
  # learning rate maximum likelihood (ml)
  lr_ml <- as.vector(f_logit(vLR[which.max(lr_marg)], 1))
  # learning rate weighted mean
  lr_wm <- as.vector(f_logit(vLR %*% lr_marg, 1))
  # learning rate variance
  lr_var <- as.vector(f_logit( ((vLR - f_logit(lr_wm,0))^2) %*% lr_marg, 1))
  
  # gamma marginal
  gamma_marg <- apply(post_prob,2,sum)
  # gamma maximum likelihood (ml)
  gamma_ml <- as.vector(f_logit(vGamma[which.max(gamma_marg)], 1))
  # gamma weighted mean
  gamma_wm <- as.vector(f_logit(vGamma %*% gamma_marg, 1))#as.vector(exp(vGamma %*% tau_marg))
  # gamma variance
  gamma_var <- as.vector(f_logit( ((vGamma - f_logit(gamma_wm,0))^2) %*% gamma_marg, 1))#as.vector(exp( ((vGamma-log(tau_wm))^2) %*% tau_marg))
  
  # eta marginal
  eta_marg <- apply(post_prob,3,sum)
  # eta maximum likelihood (ml)
  eta_ml <- as.vector(f_logit(vEta[which.max(eta_marg)], 1))
  # eta weighted mean
  eta_wm <- as.vector(f_logit(vEta %*% eta_marg, 1))
  # eta variance
  eta_var <- as.vector(f_logit( ((vEta - f_logit(eta_wm,0))^2) %*% eta_marg, 1))
  
  # lambda marginal
  lambda_marg <- apply(post_prob,4,sum)
  # lambda maximum likelihood (ml)
  lambda_ml <- as.vector(tanh(vLambda[which.max(lambda_marg)]))
  # lambda weighted mean
  lambda_wm <- as.vector(tanh(vLambda %*% lambda_marg))
  # lambda variance
  lambda_var <- as.vector(tanh( ((vLambda - atanh(lambda_wm))^2) %*% lambda_marg))
  
  # use best parameters to estimate maxLL, BIC, and hitRate
  bestMod <- mod0(list(lr=lr_wm, gamma=gamma_wm, eta=eta_wm, lambda=lambda_wm,
                       var=NULL),x)
  
  # remove NAs and adjust due to missing responses
  bestModPred <- bestMod$r_hat[valid_trials]
  
  # response model variance
  resp_var <- sqrt(sum((predClean-bestModPred)^2)/nTrials)
  # sum log-likelihood (model evidence)
  sumLL <- sum(dnorm(predClean, bestModPred, resp_var, log = T))
  # Bayesian Information Criterion
  BIC <- (nPar * log(nTrials)) - (2 * sumLL) 
  
  # set output list
  parameters <- data.frame(lr_ml=lr_ml,lr_wm=lr_wm,lr_var=lr_var,
                           gamma_ml=gamma_ml,gamma_wm=gamma_wm,gamma_var=gamma_var,
                           eta_ml=eta_ml,eta_wm=eta_wm,eta_var=eta_var,
                           lambda_ml=lambda_ml,lambda_wm=lambda_wm,lambda_var=lambda_var,
                           resp_var=resp_var)
  # end time
  end_time <- Sys.time()
  fitDuration <- end_time - start_time
  
  # model specifications
  modCompar <- data.frame(src_subject_id = toFit$src_subject_id, 
                          phenotype = toFit$phenotype,
                          site = toFit$site,
                          sumLL, BIC, fitDuration)
  # create output list
  output <- list(parameters = parameters, modCompar = modCompar, 
                 post_prob = post_prob, bestMod = bestMod)
  
  return(output)
} # end fitting function



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # Other Functions # # # # # # # # # # # # # # # # ####
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# noisy max used for weights the weights (V)
f_noisyMax <- function (v,i,tau=1) {
  weight <- v*i
  weight[i==1] <- exp(v[i==1]/tau)/sum(exp(v[i==1]/tau))
  return(v%*%weight)
}

# softmax function for two actions
f_softmax <- function(x,invT) {return(exp(invT*x)/sum(exp(invT*x)))}

# logistic function (or activation function, similar to backpropagation
f_logistic <- function(invT,netInput) {1/(1+exp(-invT*netInput))}

# logit and inverse logit function
f_logit <- function (innum, inverse) {
  #logit returns the logit or inverse logit of input. If inverse = 0 then perform logit
  if (inverse == 0) { 
    if (sum(innum < 0 | innum > 1) > 0) {
      message("Logit only defined for numbers between 0 and 1")
    }
    out <- -log((1/innum)-1)
  } else if (inverse == 1) {
    out <- 1/(1+exp(-innum))
  } else {
    warning("Inverse paramter must be 0 or 1")
    out <- NULL
  }
  return(out)
}

Phi_approx <- function(x) {
  # Bowling, Shannon R., Mohammad T. Khasawneh, Sittichai Kaewkuekool, and Byung Rae Cho. 2009.
  #"A Logistic Approximation to the Cumulative Normal Distribution."
  # Journal of Industrial Engineering and Management 2 (1): 114-27.
  return(1/(1+exp(-((0.07056*(x^3))+(1.5976*x)))))
}


outersect <- function(x, y) {
  sort(c(setdiff(x, y),
         setdiff(y, x)))
}

f_getMaxMinWeights <- function(inp,weight) {
  return(c(which(names(inp)==names(which.max(weight))),
           which(names(inp)==names(which.min(weight)))))
}


# combine participants
mergeParticipants <- function(estimateFiles,lf) {
  # read data
  for (i in 1:length(estimateFiles)) {
    # message(i)
    # load fitted data
    load(estimateFiles[i])
    # number of trials
    nTrials <- nrow(fit$bestMod$x$trial_structure)
    
    # get one long format participant
    temp0 <- lf[lf$src_subject_id==fit$modCompar$src_subject_id,] 
    
    # combine trial structure with model predictions (r_hat)
    temp <- cbind(temp0, r_hat=fit$bestMod$r_hat)
    
    if (i == 1) {
      longFormat <- temp
      param <- cbind(fit$modCompar,fit$parameters)
      pred_error <- data.frame(src_subject_id=fit$modCompar$src_subject_id,
                               phenotype=fit$modCompar$phenotype,
                               site=fit$modCompar$site,
                               fit$bestMod$x$trial_structure,
                               PE_p=fit$bestMod$PE_p,
                               PE_s=fit$bestMod$PE_s)
    } else {
      longFormat <- rbind(longFormat,temp)
      param <- rbind(param,cbind(fit$modCompar,fit$parameters))
      pred_error <- rbind(pred_error,data.frame(src_subject_id=fit$modCompar$src_subject_id,
                                                phenotype=fit$modCompar$phenotype,
                                                site=fit$modCompar$site,
                                                fit$bestMod$x$trial_structure,
                                                PE_p=fit$bestMod$PE_p,
                                                PE_s=fit$bestMod$PE_s))
    }
  }
  # add group (controls, psychics, adn psychosis)
  param$group <- recode(param$phenotype, "ncv"="psychics","hc"="controls","nc"="controls",
                        "szc"="psychosis","sze"="psychosis","szl"="psychosis")
  pred_error$group <- recode(pred_error$phenotype, "ncv"="psychics","hc"="controls","nc"="controls",
                             "szc"="psychosis","sze"="psychosis","szl"="psychosis")
  # output
  return(list(param=param,pred_error=pred_error,longFormat=longFormat))
}

mergeParticipants_recovery <- function(recoveryFiles) {
  # read data
  for (i in 1:length(recoveryFiles)) {
    # message(i)
    # load fitted data
    load(recoveryFiles[i])
  
    if (i == 1) {
      param <- recovery
    } else {
      param <- rbind(param,recovery)
    }
  }
  return(param)
}

cleanStanModel <- function(mFit,src_subject_id) {
  parMod <- summary(mFit)$summary
  if (sum(grepl("diff",rownames(parMod))) == 0) {
    other <- parMod[c(1:8,nrow(parMod)),]
    parMod <- parMod[9:(nrow(parMod)-1),]
  } else {
    other <- parMod[c(1:16,(nrow(parMod)-24):nrow(parMod)),]
    parMod <- parMod[17:(nrow(parMod)-25),]
  }
  parMod <- data.frame(parMod,src_subject_id=rep(src_subject_id,8))
  # numbers <- gsub("\\D", "", rownames(parMod))
  parMod$params <- gsub("[[:digit:][:punct:]]", "", rownames(parMod))
  parMod <- parMod[!grepl("tilde",parMod$params),]
  alpha <- parMod[parMod$params == "alpha",]
  gamma <- parMod[parMod$params == "gamma",]
  eta <- parMod[parMod$params == "eta",]
  lambda <- parMod[parMod$params == "lambda",]
  stanParams <- data.frame(src_subject_id,alpha_stan=Phi_approx(alpha$mean),
                           gamma_stan=Phi_approx(gamma$mean),eta_stan=Phi_approx(eta$mean),
                           lambda_stan=tanh(lambda$mean))
  return(stanParams)
}


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # Visualization functions # # # # # # # # # # # # ####
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
f_figure2 <- function(wf, deltas, pdi_type = "conv", psychics = T) {
  # wf <- wf[wf$group!="psychics",]
  # deltas <- deltas[deltas$group!="psychics",]
  
  # # # Behaviour # # #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  # melt prediction (for B and D)
  wf2 <- reshape2::melt(wf, measure.vars = c("B1+","B2-","D1+","D2-"))
  wf2 <- wf2[order(wf2$src_subject_id),]
  colnames(wf2)[(ncol(wf2)-1):ncol(wf2)] <- c("cue","pred")
  wf2$cue2 <- substr(wf2$cue,1,1)
  if (pdi_type == "conv") {
    wf$PDI_high <- factor(wf$pdi40_con_high, levels = c("low","high"))
    wf2$PDI_high <- factor(wf2$pdi40_con_high, levels = c("low","high"))
    wf2$PDI_cont <- wf2$pdi40_con
    pdi_label <- "PDI-C"
    deltas$PDI_high <- factor(deltas$pdi40_con_high, levels = c("low","high"))
  } else if (pdi_type == "dist")  {
    wf$PDI_high <- factor(wf$pdi40_dis_high, levels = c("low","high"))
    wf2$PDI_high <- factor(wf2$pdi40_con_high, levels = c("low","high"))
    wf2$PDI_cont <- wf2$pdi40_dis
    pdi_label <- "PDI-D"
    deltas$PDI_high <- factor(deltas$pdi40_dis_high, levels = c("low","high"))
  }

  range(wf2$pred,na.rm = T)
  # PDI conviction
  m<-lmer(pred~cue2*PDI_high+(1|site/src_subject_id),REML=F,wf2)
  anova(m); summary(m)
  figA <- ggplot(wf2[!is.na(wf2$PDI_high),],aes(x=PDI_high,y=pred,
                                                 col=cue2,fill=cue2)) + 
    labs(x=pdi_label, 
         y= "Prediction",
         col = "Trial Type:", fill = "Trial Type:") +
    geom_hline(yintercept = 0, col="black") + 
    geom_violin(alpha=0.1,position=position_dodge(0.8)) +
    stat_summary(fun=mean,position=position_dodge(0.8),geom="bar",alpha=0.7) +
    stat_summary(fun.data=mean_cl_normal,position=position_dodge(0.8),
                 geom="errorbar", col="black", width=0.1) +
    scale_y_continuous(limits = c(-1,1.14), breaks = seq(-1,1,by=1)) +
    # scale_color_manual(values = c("skyblue", "blueviolet")) +
    # scale_fill_manual(values = c("skyblue", "blueviolet")) +
    # geom_signif(comparisons = list(c("low", "high")),
    #             annotations = "ns interaction",
    #             col="black",  y_position = 1,
    #             tip_length = 0, vjust = 0.2) +
    theme_classic()
  figA
  # group
  m<-lmer(pred~cue2*group+(1|site/src_subject_id),REML=F,wf2)
  anova(m); summary(m)
  figB <- ggplot(wf2,aes(x=group,
                         y=pred,col=cue2,fill=cue2)) + 
    labs(x="Group", 
         y= "Prediction",
         col = "Trial Type:", fill = "Trial Type:") +
    geom_hline(yintercept = 0, col="black") + 
    geom_violin(alpha=0.1,position=position_dodge(0.8)) +
    stat_summary(fun=mean,position=position_dodge(0.8),geom="bar",alpha=0.7) +
    stat_summary(fun.data=mean_cl_normal,position=position_dodge(0.8),
                 geom="errorbar", col="black", width=0.1) +
    scale_y_continuous(limits = c(-1,1.14), breaks = seq(-1,1,by=1)) +
    geom_signif(comparisons=list(c("controls","psychics"), c("psychics","psychosis"),
                                 c("controls","psychosis")), 
                annotations = c("ns","ns","*"), 
                col="black", y_position = c(0.9,0.97,1.04), 
                tip_length = 0.01, vjust = 0.2) +
    theme_classic()
  figB
  
  
  # # # Prediction Error # # #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  # get deltas (data.frame) from above
  range(deltas$pe,na.rm = T)
  # PDI conviction
  # m<-lmer(pe~pe_type*trial_type*PDIc_high+(1|site/src_subject_id),REML=F,deltas)
  # anova(m); #summary(m)
  m<-lmer(pe~trial_type*PDI_high+(1|site/src_subject_id),REML=F,deltas[deltas$pe_type=="primary",])
  anova(m); summary(m)
  figC.1 <- ggplot(deltas[!is.na(deltas$PDI_high)&deltas$pe_type=="primary",],
                   aes(x=PDI_high,y=pe,col=trial_type,fill=trial_type)) + 
    labs(x=pdi_label, 
         y= expression(delta[primary]),
         col = "Trial Type:", fill = "Trial Type:") +
    geom_hline(yintercept = 0, col="black") + 
    geom_violin(alpha=0.1,position=position_dodge(0.8)) +
    stat_summary(fun=mean,position=position_dodge(0.8),geom="bar",alpha=0.7) +
    stat_summary(fun.data=mean_cl_normal,position=position_dodge(0.8),
                 geom="errorbar", col="black", width=0.1) +
    scale_y_continuous(breaks = seq(-1,2,by=1)) +
    coord_cartesian(ylim = c(-1,1.8)) +
    geom_signif(comparisons = list(c("low", "high")),
                annotations = ".",
                col="black",  y_position = 1.5,
                tip_length = 0, vjust = 0.2) +
    theme_classic() + theme(axis.title.x = element_blank())
  figC.1
  m<-lmer(pe~trial_type*PDI_high+(1|site/src_subject_id),REML=F,deltas[deltas$pe_type=="secondary",])
  anova(m); summary(m)
  figC.2 <- ggplot(deltas[!is.na(deltas$PDI_high)&deltas$pe_type=="secondary",],
                   aes(x=PDI_high,y=pe,col=trial_type,fill=trial_type)) + 
    labs(x=pdi_label, 
         y= expression(delta[secondary]),
         col = "Trial Type:", fill = "Trial Type:") +
    geom_hline(yintercept = 0, col="black") + 
    geom_violin(alpha=0.1,position=position_dodge(0.8)) +
    stat_summary(fun=mean,position=position_dodge(0.8),geom="bar",alpha=0.7) +
    stat_summary(fun.data=mean_cl_normal,position=position_dodge(0.8),
                 geom="errorbar", col="black", width=0.1) +
    scale_y_continuous(breaks = seq(-1,2,by=1)) +
    coord_cartesian(ylim = c(-1,1.8)) +
    # geom_signif(comparisons = list(c("low", "high")),
    #             annotations = "ns", 
    #             col="black",  y_position = 1.5,
    #             tip_length = 0, vjust = 0.2) +
    theme_classic() + theme(axis.title.x = element_blank())
  figC.2
  # group
  # m<-lmer(pe~pe_type*trial_type*group+(1|site/src_subject_id),REML=F,deltas)
  # anova(m); #summary(m)
  m<-lmer(pe~trial_type*group+(1|site/src_subject_id),REML=F,deltas[deltas$pe_type=="primary",])
  anova(m); summary(m)
  figD.1 <- ggplot(deltas[deltas$pe_type=="primary",],
                   aes(x=group,y=pe,col=trial_type,fill=trial_type)) + 
    labs(x="Group", 
         y= expression(delta[primary]),
         col = "Trial Type:", fill = "Trial Type:") +
    geom_hline(yintercept = 0, col="black") + 
    geom_violin(alpha=0.1,position=position_dodge(0.8)) +
    stat_summary(fun=mean,position=position_dodge(0.8),geom="bar",alpha=0.7) +
    stat_summary(fun.data=mean_cl_normal,position=position_dodge(0.8),
                 geom="errorbar", col="black", width=0.1) +
    scale_y_continuous(breaks = seq(-1,2,by=1)) +
    coord_cartesian(ylim = c(-1,1.8)) +
    geom_signif(comparisons=list(c("controls","psychics"), c("psychics","psychosis"),
                                 c("controls","psychosis")), 
                annotations = c("ns","ns","*"), 
                col="black", y_position = c(1.5,1.6,1.7), 
                tip_length = 0.01, vjust = 0.2) +
    theme_classic() + theme(axis.title.x = element_blank())
  figD.1
  m<-lmer(pe~trial_type*group+(1|site/src_subject_id),REML=F,deltas[deltas$pe_type=="secondary",])
  anova(m); summary(m)
  figD.2 <- ggplot(deltas[deltas$pe_type=="secondary",],
                   aes(x=group,y=pe,col=trial_type,fill=trial_type)) + 
    labs(x="Group", 
         y= expression(delta[secondary]),
         col = "Trial Type:", fill = "Trial Type:") +
    geom_hline(yintercept = 0, col="black") + 
    geom_violin(alpha=0.1,position=position_dodge(0.8)) +
    stat_summary(fun=mean,position=position_dodge(0.8),geom="bar",alpha=0.7) +
    stat_summary(fun.data=mean_cl_normal,position=position_dodge(0.8),
                 geom="errorbar", col="black", width=0.1) +
    scale_y_continuous(breaks = seq(-1,2,by=1)) +
    coord_cartesian(ylim = c(-1,1.8)) +
    # geom_signif(comparisons=list(c("controls","psychics"), c("psychics","psychosis"),
    #                              c("controls","psychosis")), 
    #             annotations = c("ns","ns","*"), 
    #             col="black", y_position = c(0.9,0.97,1.04), 
    #             tip_length = 0.01, vjust = 0.2) +
    theme_classic() + theme(axis.title.x = element_blank())
  figD.2
  
  
  # # # Alpha # # #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  range(f_logit(wf$lr_wm,0),na.rm = T)
  # PDI
  m<-lmer(f_logit(lr_wm,0)~PDI_high+(1|site),REML=F,wf)
  anova(m); summary(m)
  figE <- ggplot(wf[!is.na(wf$PDI_high),],aes(x=PDI_high,
                                              y=f_logit(lr_wm,0),
                                              col=PDI_high,fill=PDI_high)) + 
    labs(x=pdi_label, 
         y= expression(alpha), col=pdi_label,fill=pdi_label) +
    geom_hline(yintercept = 0, col="black") + 
    geom_violin(alpha=0.1,position=position_dodge(0.8)) +
    stat_summary(fun=mean,position=position_dodge(0.8),geom="bar",alpha=0.7) +
    stat_summary(fun.data=mean_cl_normal,position=position_dodge(0.8),
                 geom="errorbar", col="black", width=0.1) +
    scale_y_continuous(limits = c(-7,1), breaks = seq(-6,6,by=2)) +
    scale_color_manual(values = c("blueviolet","skyblue")) +
    scale_fill_manual(values = c("blueviolet","skyblue")) +
    geom_signif(comparisons = list(c("low", "high")), 
                annotations = "ns",
                col="black", y_position = 0.5,
                tip_length = 0, vjust = 0.2) +
    theme_classic() + theme(axis.title.x = element_blank())
  figE
  # group
  m<-lmer(f_logit(lr_wm,0)~group+(1|site),REML=F,wf)
  anova(m); summary(m)
  # anova(lmer(f_logit(lr_wm,0)~group+(1|site),REML=F,wf[wf$group!="controls",]))
  figF <- ggplot(wf[,],aes(x=group,
                           y=f_logit(lr_wm,0),
                           col=group,fill=group)) + 
    labs(x="Group", 
         y= expression(alpha), col="Group",fill="Group") +
    geom_hline(yintercept = 0, col="black") + 
    geom_violin(alpha=0.1,position=position_dodge(0.8)) +
    stat_summary(fun=mean,position=position_dodge(0.8),geom="bar",alpha=0.7) +
    stat_summary(fun.data=mean_cl_normal,position=position_dodge(0.8),
                 geom="errorbar", col="black", width=0.1) +
    scale_y_continuous(limits = c(-7,1.15), breaks = seq(-6,6,by=2)) +
    scale_color_manual(values = c("orange","yellow","red")) +
    scale_fill_manual(values = c("orange","yellow","red")) +
    geom_signif(comparisons=list(c("controls","psychics"), c("controls","psychosis"),
                                 c("psychics","psychosis")), 
                annotations = c("*","***","ns"),
                col="black", y_position = c(-0,0.7,0.35), 
                tip_length = 0.01, vjust = 0.2) +
    theme_classic() + theme(axis.title.x = element_blank())
  figF
  
  
  # # # Gamma # # #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  range(f_logit(wf$gamma_wm,0),na.rm = T)
  # PDI
  m<-lmer(f_logit(gamma_wm,0)~PDI_high+(1|site),REML=F,wf)
  anova(m); summary(m)
  figG <- ggplot(wf[!is.na(wf$PDI_high),],aes(x=PDI_high,
                                              y=f_logit(gamma_wm,0),
                                              col=PDI_high,fill=PDI_high)) + 
    labs(x=pdi_label, 
         y= expression(gamma), col=pdi_label,fill=pdi_label) +
    geom_hline(yintercept = 0, col="black") + 
    geom_violin(alpha=0.1,position=position_dodge(0.8)) +
    stat_summary(fun=mean,position=position_dodge(0.8),geom="bar",alpha=0.7) +
    stat_summary(fun.data=mean_cl_normal,position=position_dodge(0.8),
                 geom="errorbar", col="black", width=0.1) +
    scale_y_continuous(limits = c(-6,4.5), breaks = seq(-6,6,by=2)) +
    scale_color_manual(values = c("blueviolet","skyblue")) +
    scale_fill_manual(values = c("blueviolet","skyblue")) +
    geom_signif(comparisons = list(c("low", "high")), 
                annotations = "ns",
                col="black", y_position = 4,
                tip_length = 0, vjust = 0.2) +
    theme_classic() + theme(axis.title.x = element_blank())
  figG
  # group
  m<-lmer(f_logit(gamma_wm,0)~group+(1|site),REML=F,wf)
  anova(m); summary(m)
  figH <- ggplot(wf[,],aes(x=group,
                           y=f_logit(gamma_wm,0),
                           col=group,fill=group)) + 
    labs(x="Group", 
         y= expression(gamma), col="Group",fill="Group") +
    geom_hline(yintercept = 0, col="black") + 
    geom_violin(alpha=0.1,position=position_dodge(0.8)) +
    stat_summary(fun=mean,position=position_dodge(0.8),geom="bar",alpha=0.7) +
    stat_summary(fun.data=mean_cl_normal,position=position_dodge(0.8),
                 geom="errorbar", col="black", width=0.1) +
    scale_y_continuous(limits = c(-6,4.5), breaks = seq(-6,4.5,by=2)) +
    scale_color_manual(values = c("orange","yellow","red")) +
    scale_fill_manual(values = c("orange","yellow","red")) +
    geom_signif(comparisons = list(c("controls", "psychosis")), 
                annotations = "ns",
                col="black", y_position = 4,
                tip_length = 0, vjust = 0.2) +
    theme_classic() + theme(axis.title.x = element_blank())
  figH
  
  
  # # # Eta # # #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  range(f_logit(wf$gamma_wm,0),na.rm = T)
  # PDI
  m<-lmer(f_logit(eta_wm,0)~PDI_high+(1|site),REML=F,wf)
  anova(m); summary(m)
  figI <- ggplot(wf[!is.na(wf$PDI_high),],aes(x=PDI_high,
                                              y=f_logit(eta_wm,0),
                                              col=PDI_high,fill=PDI_high)) + 
    labs(x=pdi_label, 
         y= expression(eta), col=pdi_label,fill=pdi_label) +
    geom_hline(yintercept = 0, col="black") + 
    geom_violin(alpha=0.1,position=position_dodge(0.8)) +
    stat_summary(fun=mean,position=position_dodge(0.8),geom="bar",alpha=0.7) +
    stat_summary(fun.data=mean_cl_normal,position=position_dodge(0.8),
                 geom="errorbar", col="black", width=0.1) +
    scale_y_continuous(limits = c(-5.5,4.5), breaks = seq(-6,6,by=2)) +
    scale_color_manual(values = c("blueviolet","skyblue")) +
    scale_fill_manual(values = c("blueviolet","skyblue")) +
    geom_signif(comparisons = list(c("low", "high")), 
                annotations = "ns",
                col="black", y_position = 4,
                tip_length = 0, vjust = 0.2) +
    theme_classic() + theme(axis.title.x = element_blank())
  figI
  # group
  m<-lmer(f_logit(eta_wm,0)~group+(1|site),REML=F,wf)
  anova(m); summary(m)
  figJ <- ggplot(wf[,],aes(x=group,
                           y=f_logit(eta_wm,0),
                           col=group,fill=group)) + 
    labs(x="Group", 
         y= expression(eta), col="Group",fill="Group") +
    geom_hline(yintercept = 0, col="black") + 
    geom_violin(alpha=0.1,position=position_dodge(0.8)) +
    stat_summary(fun=mean,position=position_dodge(0.8),geom="bar",alpha=0.7) +
    stat_summary(fun.data=mean_cl_normal,position=position_dodge(0.8),
                 geom="errorbar", col="black", width=0.1) +
    scale_y_continuous(limits = c(-5.5,4.5), breaks = seq(-6,6,by=2)) +
    scale_color_manual(values = c("orange","yellow","red")) +
    scale_fill_manual(values = c("orange","yellow","red")) +
    geom_signif(comparisons = list(c("controls", "psychosis")), 
                annotations = "ns",
                col="black", y_position = 4,
                tip_length = 0, vjust = 0.2) +
    theme_classic() + theme(axis.title.x = element_blank())
  figJ
  
  
  # # # Lambda # # #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  range(atanh(wf$lambda_wm),na.rm = T)
  # PDI
  m<-lmer(atanh(lambda_wm)~PDI_high+(1|site),REML=F,wf)
  anova(m); summary(m)
  figK <- ggplot(wf[!is.na(wf$PDI_high),],aes(x=PDI_high,
                                              y=atanh(lambda_wm),
                                              col=PDI_high,fill=PDI_high)) + 
    labs(x=pdi_label, 
         y= expression(lambda), col=pdi_label,fill=pdi_label) +
    geom_violin(alpha=0.1,position=position_dodge(0.8)) +
    stat_summary(fun=mean,position=position_dodge(0.8),geom="bar",alpha=0.7) +
    stat_summary(fun.data=mean_cl_normal,position=position_dodge(0.8),
                 geom="errorbar", col="black", width=0.1) +
    scale_y_continuous(limits = c(-3,3), breaks = seq(-6,6,by=1.5)) +
    scale_color_manual(values = c("blueviolet","skyblue")) +
    scale_fill_manual(values = c("blueviolet","skyblue")) +
    geom_signif(comparisons = list(c("low", "high")), 
                annotations = ".",
                col="black", y_position = 2.5,
                tip_length = 0, vjust = 0.2) +
    theme_classic() + theme(axis.title.x = element_blank())
  figK
  
  # group
  m<-lmer(atanh(lambda_wm)~group+(1|site),REML=F,wf)
  anova(m); summary(m)
  figL <- ggplot(wf[,],aes(x=group,
                           y=atanh(lambda_wm),
                           col=group,fill=group)) + 
    labs(x="Group", 
         y= expression(lambda), col="Group",fill="Group") +
    geom_violin(alpha=0.1,position=position_dodge(0.8)) +
    stat_summary(fun=mean,position=position_dodge(0.8),geom="bar",alpha=0.7) +
    stat_summary(fun.data=mean_cl_normal,position=position_dodge(0.8),
                 geom="errorbar", col="black", width=0.1) +
    scale_y_continuous(limits = c(-3,3), breaks = seq(-6,6,by=1.5)) +
    scale_color_manual(values = c("orange","yellow","red")) +
    scale_fill_manual(values = c("orange","yellow","red")) +
    geom_signif(comparisons = list(c("controls", "psychosis")), 
                annotations = "ns",
                col="black", y_position = 2.5,
                tip_length = 0, vjust = 0.2) +
    theme_classic() + theme(axis.title.x = element_blank())
  figL
  
  
  
  # # # Combine plots # # #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  
  
  # create left panel
  upperLeft <- annotate_figure(ggarrange(figA, figB, ncol=2, 
                                         align = "hv", widths = c(2,3), 
                                         common.legend = T, labels = c("A","B")),
                               top = text_grob("Behaviour (Phase 3)", color = "black", 
                                               face = "bold", size = 14))
  bottomLeft <- annotate_figure(ggarrange(figC.1, figD.1, figC.2, figD.2, nrow=2, 
                                          ncol=2, align = "hv", widths = c(2,3), 
                                          common.legend = T, labels = c("C","D","E","F")),
                                top = text_grob("Prediction Errors (Phase 2)", color = "black", 
                                                face = "bold", size = 14))
  left <- ggarrange(upperLeft,bottomLeft, nrow=2, heights = c(1,2))
  # create right panel
  right <- annotate_figure(ggarrange(
    ggarrange(figE,figG,figI,figK,nrow=4,ncol=1,align="hv",common.legend = T,labels = c("G","I","K","M")),
    ggarrange(figF,figH,figJ,figL,nrow=4,ncol=1,align="hv",common.legend = T,labels = c("H","J","L","N")),
    ncol = 2, widths = c(2,3)),
    top = text_grob("Rosa's Kamin Model", color = "black",face = "bold", size = 14))
  # combine left and right panels
  figure2 <- ggarrange(left,NULL,right,ncol=3,widths = c(10,1,10))
  
  return(figure2)
}

f_figure2_noPsychics <- function(wf, deltas, pdi_type = "conv") {
  wf <- wf[wf$group!="psychics",]
  deltas <- deltas[deltas$group!="psychics",]
  
  # # # Behaviour # # #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  # melt prediction (for B and D)
  wf2 <- reshape2::melt(wf, measure.vars = c("B1+","B2-","D1+","D2-"))
  wf2 <- wf2[order(wf2$src_subject_id),]
  colnames(wf2)[(ncol(wf2)-1):ncol(wf2)] <- c("cue","pred")
  wf2$cue2 <- substr(wf2$cue,1,1)
  if (pdi_type == "conv") {
    wf$PDI_high <- factor(wf$PDIc_high, levels = c("low","high"))
    wf2$PDI_high <- factor(wf2$PDIc_high, levels = c("low","high"))
    wf2$PDI_cont <- wf2$PDI_conv
    pdi_label <- "PDI-C"
    deltas$PDI_high <- factor(deltas$PDIc_high, levels = c("low","high"))
  } else if (pdi_type == "dist")  {
    wf$PDI_high <- factor(wf$PDId_high, levels = c("low","high"))
    wf2$PDI_high <- factor(wf2$PDId_high, levels = c("low","high"))
    wf2$PDI_cont <- wf2$PDI_dist
    pdi_label <- "PDI-D"
    deltas$PDI_high <- factor(deltas$PDId_high, levels = c("low","high"))
  }
  
  range(wf2$pred,na.rm = T)
  # PDI conviction
  m<-lmer(pred~cue2*PDI_high+(1|site/src_subject_id),REML=F,wf2)
  anova(m); summary(m)
  figA <- ggplot(wf2[!is.na(wf2$PDI_high),],aes(x=PDI_high,y=pred,
                                                col=cue2,fill=cue2)) + 
    labs(x=pdi_label, 
         y= "Prediction",
         col = "Trial Type:", fill = "Trial Type:") +
    geom_hline(yintercept = 0, col="black") + 
    geom_violin(alpha=0.1,position=position_dodge(0.8)) +
    stat_summary(fun=mean,position=position_dodge(0.8),geom="bar",alpha=0.7) +
    stat_summary(fun.data=mean_cl_normal,position=position_dodge(0.8),
                 geom="errorbar", col="black", width=0.1) +
    scale_y_continuous(limits = c(-1,1.14), breaks = seq(-1,1,by=1)) +
    # scale_color_manual(values = c("skyblue", "blueviolet")) +
    # scale_fill_manual(values = c("skyblue", "blueviolet")) +
    # geom_signif(comparisons = list(c("low", "high")),
    #             annotations = "ns interaction",
    #             col="black",  y_position = 1,
    #             tip_length = 0, vjust = 0.2) +
    theme_classic()
  figA
  # group
  m<-lmer(pred~cue2*group+(1|site/src_subject_id),REML=F,wf2)
  anova(m); summary(m)
  figB <- ggplot(wf2,aes(x=group,
                         y=pred,col=cue2,fill=cue2)) + 
    labs(x="Group", 
         y= "Prediction",
         col = "Trial Type:", fill = "Trial Type:") +
    geom_hline(yintercept = 0, col="black") + 
    geom_violin(alpha=0.1,position=position_dodge(0.8)) +
    stat_summary(fun=mean,position=position_dodge(0.8),geom="bar",alpha=0.7) +
    stat_summary(fun.data=mean_cl_normal,position=position_dodge(0.8),
                 geom="errorbar", col="black", width=0.1) +
    scale_y_continuous(limits = c(-1,1.14), breaks = seq(-1,1,by=1)) +
    geom_signif(comparisons=list(c("controls","psychosis")), 
                annotations = "*", 
                col="black", y_position = c(1), 
                tip_length = 0.01, vjust = 0.2) +
    theme_classic()
  figB
  
  
  # # # Prediction Error # # #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  # get deltas (data.frame) from above
  range(deltas$pe,na.rm = T)
  # PDI conviction
  # m<-lmer(pe~pe_type*trial_type*PDI_high+(1|site/src_subject_id),REML=F,deltas)
  # anova(m); #summary(m)
  m<-lmer(pe~trial_type*PDI_high+(1|site/src_subject_id),REML=F,deltas[deltas$pe_type=="primary",])
  anova(m); summary(m)
  figC.1 <- ggplot(deltas[!is.na(deltas$PDI_high)&deltas$pe_type=="primary",],
                   aes(x=PDI_high,y=pe,col=trial_type,fill=trial_type)) + 
    labs(x=pdi_label, 
         y= expression(delta[primary]),
         col = "Trial Type:", fill = "Trial Type:") +
    geom_hline(yintercept = 0, col="black") + 
    geom_violin(alpha=0.1,position=position_dodge(0.8)) +
    stat_summary(fun=mean,position=position_dodge(0.8),geom="bar",alpha=0.7) +
    stat_summary(fun.data=mean_cl_normal,position=position_dodge(0.8),
                 geom="errorbar", col="black", width=0.1) +
    scale_y_continuous(breaks = seq(-1,2,by=1)) +
    coord_cartesian(ylim = c(-1,1.8)) +
    # geom_signif(comparisons = list(c("low", "high")),
    #             annotations = "*",
    #             col="black",  y_position = 1.5,
    #             tip_length = 0.01, vjust = 0.2) +
    theme_classic() + theme(axis.title.x = element_blank())
  figC.1
  m<-lmer(pe~trial_type*PDI_high+(1|site/src_subject_id),REML=F,deltas[deltas$pe_type=="secondary",])
  anova(m); summary(m)
  figC.2 <- ggplot(deltas[!is.na(deltas$PDI_high)&deltas$pe_type=="secondary",],
                   aes(x=PDI_high,y=pe,col=trial_type,fill=trial_type)) + 
    labs(x=pdi_label, 
         y= expression(delta[secondary]),
         col = "Trial Type:", fill = "Trial Type:") +
    geom_hline(yintercept = 0, col="black") + 
    geom_violin(alpha=0.1,position=position_dodge(0.8)) +
    stat_summary(fun=mean,position=position_dodge(0.8),geom="bar",alpha=0.7) +
    stat_summary(fun.data=mean_cl_normal,position=position_dodge(0.8),
                 geom="errorbar", col="black", width=0.1) +
    scale_y_continuous(breaks = seq(-1,2,by=1)) +
    coord_cartesian(ylim = c(-1,1.8)) +
    # geom_signif(comparisons = list(c("low", "high")),
    #             annotations = "ns", 
    #             col="black",  y_position = 1.5,
    #             tip_length = 0, vjust = 0.2) +
    theme_classic() + theme(axis.title.x = element_blank())
  figC.2
  # group
  # m<-lmer(pe~pe_type*trial_type*group+(1|site/src_subject_id),REML=F,deltas)
  # anova(m); #summary(m)
  m<-lmer(pe~trial_type*group+(1|site/src_subject_id),REML=F,deltas[deltas$pe_type=="primary",])
  anova(m); summary(m)
  figD.1 <- ggplot(deltas[deltas$pe_type=="primary",],
                   aes(x=group,y=pe,col=trial_type,fill=trial_type)) + 
    labs(x="Group", 
         y= expression(delta[primary]),
         col = "Trial Type:", fill = "Trial Type:") +
    geom_hline(yintercept = 0, col="black") + 
    geom_violin(alpha=0.1,position=position_dodge(0.8)) +
    stat_summary(fun=mean,position=position_dodge(0.8),geom="bar",alpha=0.7) +
    stat_summary(fun.data=mean_cl_normal,position=position_dodge(0.8),
                 geom="errorbar", col="black", width=0.1) +
    scale_y_continuous(breaks = seq(-1,2,by=1)) +
    coord_cartesian(ylim = c(-1,1.8)) +
    geom_signif(comparisons=list(c("controls","psychosis")), 
                annotations = c("*"), 
                col="black", y_position = c(1.5), 
                tip_length = 0.01, vjust = 0.2) +
    theme_classic() + theme(axis.title.x = element_blank())
  figD.1
  m<-lmer(pe~trial_type*group+(1|site/src_subject_id),REML=F,deltas[deltas$pe_type=="secondary",])
  anova(m); summary(m)
  figD.2 <- ggplot(deltas[deltas$pe_type=="secondary",],
                   aes(x=group,y=pe,col=trial_type,fill=trial_type)) + 
    labs(x="Group", 
         y= expression(delta[secondary]),
         col = "Trial Type:", fill = "Trial Type:") +
    geom_hline(yintercept = 0, col="black") + 
    geom_violin(alpha=0.1,position=position_dodge(0.8)) +
    stat_summary(fun=mean,position=position_dodge(0.8),geom="bar",alpha=0.7) +
    stat_summary(fun.data=mean_cl_normal,position=position_dodge(0.8),
                 geom="errorbar", col="black", width=0.1) +
    scale_y_continuous(breaks = seq(-1,2,by=1)) +
    coord_cartesian(ylim = c(-1,1.8)) +
    # geom_signif(comparisons=list(c("controls","psychics"), c("psychics","psychosis"),
    #                              c("controls","psychosis")), 
    #             annotations = c("ns","ns","*"), 
    #             col="black", y_position = c(0.9,0.97,1.04), 
    #             tip_length = 0.01, vjust = 0.2) +
    theme_classic() + theme(axis.title.x = element_blank())
  figD.2
  
  
  # # # Alpha # # #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  range(f_logit(wf$lr_wm,0),na.rm = T)
  # PDI
  m<-lmer(f_logit(lr_wm,0)~PDI_high+(1|site),REML=F,wf)
  anova(m); summary(m)
  figE <- ggplot(wf[!is.na(wf$PDI_high),],aes(x=PDI_high,
                                              y=f_logit(lr_wm,0),
                                              col=PDI_high,fill=PDI_high)) + 
    labs(x=pdi_label, 
         y= expression(alpha), col=pdi_label,fill=pdi_label) +
    geom_hline(yintercept = 0, col="black") + 
    geom_violin(alpha=0.1,position=position_dodge(0.8)) +
    stat_summary(fun=mean,position=position_dodge(0.8),geom="bar",alpha=0.7) +
    stat_summary(fun.data=mean_cl_normal,position=position_dodge(0.8),
                 geom="errorbar", col="black", width=0.1) +
    scale_y_continuous(limits = c(-7,1), breaks = seq(-6,6,by=2)) +
    scale_color_manual(values = c("blueviolet","skyblue")) +
    scale_fill_manual(values = c("blueviolet","skyblue")) +
    geom_signif(comparisons = list(c("low", "high")), 
                annotations = "ns",
                col="black", y_position = 0.5,
                tip_length = 0, vjust = 0.2) +
    theme_classic() + theme(axis.title.x = element_blank())
  figE
  # group
  m<-lmer(f_logit(lr_wm,0)~group+(1|site),REML=F,wf)
  anova(m); summary(m)
  # anova(lmer(f_logit(lr_wm,0)~group+(1|site),REML=F,wf[wf$group!="controls",]))
  figF <- ggplot(wf[,],aes(x=group,
                           y=f_logit(lr_wm,0),
                           col=group,fill=group)) + 
    labs(x="Group", 
         y= expression(alpha), col="Group",fill="Group") +
    geom_hline(yintercept = 0, col="black") + 
    geom_violin(alpha=0.1,position=position_dodge(0.8)) +
    stat_summary(fun=mean,position=position_dodge(0.8),geom="bar",alpha=0.7) +
    stat_summary(fun.data=mean_cl_normal,position=position_dodge(0.8),
                 geom="errorbar", col="black", width=0.1) +
    scale_y_continuous(limits = c(-7,1.15), breaks = seq(-6,6,by=2)) +
    scale_color_manual(values = c("orange","red")) +
    scale_fill_manual(values = c("orange","red")) +
    geom_signif(comparisons=list(c("controls","psychosis")), 
                annotations = c("***"),
                col="black", y_position = c(0.5), 
                tip_length = 0.01, vjust = 0.2) +
    theme_classic() + theme(axis.title.x = element_blank())
  figF
  
  
  # # # Gamma # # #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  range(f_logit(wf$gamma_wm,0),na.rm = T)
  # PDI
  m<-lmer(f_logit(gamma_wm,0)~PDI_high+(1|site),REML=F,wf)
  anova(m); summary(m)
  figG <- ggplot(wf[!is.na(wf$PDI_high),],aes(x=PDI_high,
                                              y=f_logit(gamma_wm,0),
                                              col=PDI_high,fill=PDI_high)) + 
    labs(x=pdi_label, 
         y= expression(gamma), col=pdi_label,fill=pdi_label) +
    geom_hline(yintercept = 0, col="black") + 
    geom_violin(alpha=0.1,position=position_dodge(0.8)) +
    stat_summary(fun=mean,position=position_dodge(0.8),geom="bar",alpha=0.7) +
    stat_summary(fun.data=mean_cl_normal,position=position_dodge(0.8),
                 geom="errorbar", col="black", width=0.1) +
    scale_y_continuous(limits = c(-6,4.5), breaks = seq(-6,6,by=2)) +
    scale_color_manual(values = c("blueviolet","skyblue")) +
    scale_fill_manual(values = c("blueviolet","skyblue")) +
    geom_signif(comparisons = list(c("low", "high")), 
                annotations = ".",
                col="black", y_position = 4,
                tip_length = 0, vjust = 0.2) +
    theme_classic() + theme(axis.title.x = element_blank())
  figG
  # group
  m<-lmer(f_logit(gamma_wm,0)~group+(1|site),REML=F,wf)
  anova(m); summary(m)
  figH <- ggplot(wf[,],aes(x=group,
                           y=f_logit(gamma_wm,0),
                           col=group,fill=group)) + 
    labs(x="Group", 
         y= expression(gamma), col="Group",fill="Group") +
    geom_hline(yintercept = 0, col="black") + 
    geom_violin(alpha=0.1,position=position_dodge(0.8)) +
    stat_summary(fun=mean,position=position_dodge(0.8),geom="bar",alpha=0.7) +
    stat_summary(fun.data=mean_cl_normal,position=position_dodge(0.8),
                 geom="errorbar", col="black", width=0.1) +
    scale_y_continuous(limits = c(-6,4.5), breaks = seq(-6,4.5,by=2)) +
    scale_color_manual(values = c("orange","red")) +
    scale_fill_manual(values = c("orange","red")) +
    geom_signif(comparisons = list(c("controls", "psychosis")), 
                annotations = "ns",
                col="black", y_position = 4,
                tip_length = 0, vjust = 0.2) +
    theme_classic() + theme(axis.title.x = element_blank())
  figH
  
  
  # # # Eta # # #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  range(f_logit(wf$gamma_wm,0),na.rm = T)
  # PDI
  m<-lmer(f_logit(eta_wm,0)~PDI_high+(1|site),REML=F,wf)
  anova(m); summary(m)
  figI <- ggplot(wf[!is.na(wf$PDI_high),],aes(x=PDI_high,
                                              y=f_logit(eta_wm,0),
                                              col=PDI_high,fill=PDI_high)) + 
    labs(x=pdi_label, 
         y= expression(eta), col=pdi_label,fill=pdi_label) +
    geom_hline(yintercept = 0, col="black") + 
    geom_violin(alpha=0.1,position=position_dodge(0.8)) +
    stat_summary(fun=mean,position=position_dodge(0.8),geom="bar",alpha=0.7) +
    stat_summary(fun.data=mean_cl_normal,position=position_dodge(0.8),
                 geom="errorbar", col="black", width=0.1) +
    scale_y_continuous(limits = c(-5.5,4.5), breaks = seq(-6,6,by=2)) +
    scale_color_manual(values = c("blueviolet","skyblue")) +
    scale_fill_manual(values = c("blueviolet","skyblue")) +
    geom_signif(comparisons = list(c("low", "high")), 
                annotations = "ns",
                col="black", y_position = 4,
                tip_length = 0, vjust = 0.2) +
    theme_classic() + theme(axis.title.x = element_blank())
  figI
  # group
  m<-lmer(f_logit(eta_wm,0)~group+(1|site),REML=F,wf)
  anova(m); summary(m)
  figJ <- ggplot(wf[,],aes(x=group,
                           y=f_logit(eta_wm,0),
                           col=group,fill=group)) + 
    labs(x="Group", 
         y= expression(eta), col="Group",fill="Group") +
    geom_hline(yintercept = 0, col="black") + 
    geom_violin(alpha=0.1,position=position_dodge(0.8)) +
    stat_summary(fun=mean,position=position_dodge(0.8),geom="bar",alpha=0.7) +
    stat_summary(fun.data=mean_cl_normal,position=position_dodge(0.8),
                 geom="errorbar", col="black", width=0.1) +
    scale_y_continuous(limits = c(-5.5,4.5), breaks = seq(-6,6,by=2)) +
    scale_color_manual(values = c("orange","red")) +
    scale_fill_manual(values = c("orange","red")) +
    geom_signif(comparisons = list(c("controls", "psychosis")), 
                annotations = "ns",
                col="black", y_position = 4,
                tip_length = 0, vjust = 0.2) +
    theme_classic() + theme(axis.title.x = element_blank())
  figJ
  
  
  # # # Lambda # # #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  range(atanh(wf$lambda_wm),na.rm = T)
  # PDI
  m<-lmer(atanh(lambda_wm)~PDI_high+(1|site),REML=F,wf)
  anova(m); summary(m)
  figK <- ggplot(wf[!is.na(wf$PDI_high),],aes(x=PDI_high,
                                              y=atanh(lambda_wm),
                                              col=PDI_high,fill=PDI_high)) + 
    labs(x=pdi_label, 
         y= expression(lambda), col=pdi_label,fill=pdi_label) +
    geom_violin(alpha=0.1,position=position_dodge(0.8)) +
    stat_summary(fun=mean,position=position_dodge(0.8),geom="bar",alpha=0.7) +
    stat_summary(fun.data=mean_cl_normal,position=position_dodge(0.8),
                 geom="errorbar", col="black", width=0.1) +
    scale_y_continuous(limits = c(-3,3), breaks = seq(-6,6,by=1.5)) +
    scale_color_manual(values = c("blueviolet","skyblue")) +
    scale_fill_manual(values = c("blueviolet","skyblue")) +
    geom_signif(comparisons = list(c("low", "high")), 
                annotations = "ns",
                col="black", y_position = 2.5,
                tip_length = 0, vjust = 0.2) +
    theme_classic() + theme(axis.title.x = element_blank())
  figK
  
  # group
  m<-lmer(atanh(lambda_wm)~group+(1|site),REML=F,wf)
  anova(m); summary(m)
  figL <- ggplot(wf[,],aes(x=group,
                           y=atanh(lambda_wm),
                           col=group,fill=group)) + 
    labs(x="Group", 
         y= expression(lambda), col="Group",fill="Group") +
    geom_violin(alpha=0.1,position=position_dodge(0.8)) +
    stat_summary(fun=mean,position=position_dodge(0.8),geom="bar",alpha=0.7) +
    stat_summary(fun.data=mean_cl_normal,position=position_dodge(0.8),
                 geom="errorbar", col="black", width=0.1) +
    scale_y_continuous(limits = c(-3,3), breaks = seq(-6,6,by=1.5)) +
    scale_color_manual(values = c("orange","red")) +
    scale_fill_manual(values = c("orange","red")) +
    geom_signif(comparisons = list(c("controls", "psychosis")), 
                annotations = "ns",
                col="black", y_position = 2.5,
                tip_length = 0, vjust = 0.2) +
    theme_classic() + theme(axis.title.x = element_blank())
  figL
  
  
  
  # # # Combine plots # # #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  
  
  # create left panel
  upperLeft <- annotate_figure(ggarrange(figA, figB, ncol=2, 
                                         align = "hv", widths = c(1,1), 
                                         common.legend = T, labels = c("A","B")),
                               top = text_grob("Behaviour (Phase 3)", color = "black", 
                                               face = "bold", size = 14))
  bottomLeft <- annotate_figure(ggarrange(figC.1, figD.1, figC.2, figD.2, nrow=2, 
                                          ncol=2, align = "hv", widths = c(1,1), 
                                          common.legend = T, labels = c("C","D","E","F")),
                                top = text_grob("Prediction Errors (Phase 2)", color = "black", 
                                                face = "bold", size = 14))
  left <- ggarrange(upperLeft,bottomLeft, nrow=2, heights = c(1,2))
  # create right panel
  right <- annotate_figure(ggarrange(
    ggarrange(figE,figG,figI,figK,nrow=4,ncol=1,align="hv",common.legend = T,labels = c("G","I","K","M")),
    ggarrange(figF,figH,figJ,figL,nrow=4,ncol=1,align="hv",common.legend = T,labels = c("H","J","L","N")),
    ncol = 2, widths = c(1,1)),
    top = text_grob("Rosa's Kamin Model", color = "black",face = "bold", size = 14))
  # combine left and right panels
  figure2 <- ggarrange(left,NULL,right,ncol=3,widths = c(10,1,10))
  
  return(figure2)
}

seeBGGM <- function (param,vec,vec_labs) {
  dims <- param[,vec]
  for (i in 1:ncol(dims)) {
    dims[,i] <- scale(dims[,i])[1:nrow(dims)]
  }
  # dims$pdi40_i10_unreality <- param$pdi40_i10_unreality
  # remove nas
  dims <- dims[rowSums(is.na(dims))==0,]
  # estimate network
  comp.net <- estimate(dims, iter = 10000, type = "continuous")
  # comp.net <- estimate(dims, iter = 10000, type = "mixed", mixed_type = c(0,1,0,0,0))
  # select cool connections
  comp.sel.net <- BGGM::select(comp.net, cred = 0.95)
  # BGGM
  qgraph(comp.sel.net$pcor_adj, layout = "spring",
         labels = vec_labs)
  # output
  return(list(summary(comp.net),comp.sel.net))
}
