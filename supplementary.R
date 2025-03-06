rm(list = ls()) # erase previous Global Environment, import manually the file

# libraries
if (!require(reshape2)) install.packages("reshape2"); library(reshape2)
if (!require(ggplot2)) install.packages("ggplot2"); library(ggplot2)
if (!require(ggpubr)) install.packages("ggpubr"); library(ggpubr)
if (!require(dplyr)) install.packages("dplyr"); library(dplyr)
if (!require(lmerTest)) install.packages("lmerTest"); library(lmerTest)

# print figures?
print_figures <- 0

# import data and functions
source("functions.R")

# merge all participants and create csv
if (!file.exists("data/mod0_summary.csv")) {
  # all, recover and estimates files names
  all <- paste0("data/",list.files("data/", pattern = "mod0"))
  recoveryFiles <- all[grepl("_recovery",all)]
  recoveryFiles <- recoveryFiles[!grepl(".csv",recoveryFiles)]
  estimateFiles <- outersect(all,recoveryFiles); remove(all)
  estimateFiles <- estimateFiles[!grepl(".csv",estimateFiles)]
  
  # add r_hat to it
  lf <- read.csv("data/data_modelling_v2.csv")
  
  # merge participants
  estimate <- mergeParticipants(estimateFiles,lf)
  recovery <- mergeParticipants_recovery(recoveryFiles)
  
  
  # read data
  write.csv(recovery,"data/mod0_recovery.csv",row.names = F)
  write.csv(estimate$param,"data/mod0_summary.csv",row.names = F)
  write.csv(estimate$pred_error,"data/mod0_pred_error.csv",row.names = F)
  write.csv(estimate$longFormat, "data/data_modelling_v2.csv",row.names = F)
  r_param <- read.csv("data/mod0_recovery.csv")
  e_param <- read.csv("data/mod0_summary.csv")
  e_pred_error <- read.csv("data/mod0_pred_error.csv")
  lf <- read.csv("data/data_modelling_v2.csv")
} else {
  r_param <- read.csv("data/mod0_recovery.csv")
  e_param <- read.csv("data/mod0_summary.csv")
  e_pred_error <- read.csv("data/mod0_pred_error.csv")
  lf <- read.csv("data/data_modelling_v2.csv")
}


# subject ID vectors
src_subject_id <- unique(e_param$src_subject_id)
src_subject_id <- src_subject_id[order(src_subject_id)]

# stan model clean
# stanParams <- cleanStanModel(m0Fit,src_subject_id)

# average iterations
temp <- data.frame(r_param %>% group_by(src_subject_id) %>%
                     summarise(lr_rec = max(lr_wm), gamma_rec = max(gamma_wm), 
                               eta_rec = max(eta_wm), lambda_rec = max(lambda_wm)))
# get the best recovery
for (i in 1:length(src_subject_id)) {
  temp0 <- r_param[r_param$src_subject_id == src_subject_id[i],]
  temp0 <- temp0[temp0$sumLL == max(temp0$sumLL),]
  if (i == 1) {
    temp <- data.frame(src_subject_id=src_subject_id[i], lr_rec=temp0$lr_wm, 
                       gamma_rec=temp0$gamma_wm,eta_rec=temp0$eta_wm, 
                       lambda_rec=temp0$lambda_wm)
  } else {
    temp <- rbind(temp,data.frame(src_subject_id=src_subject_id[i], lr_rec=temp0$lr_wm, 
                                  gamma_rec=temp0$gamma_wm,eta_rec=temp0$eta_wm, 
                                  lambda_rec=temp0$lambda_wm))
  }
}

# order src_subject_id
temp <- temp[order(temp$src_subject_id),]
e_param <- e_param[order(e_param$src_subject_id),]

# sanity check, it must be 142
sum(temp$src_subject_id == e_param$src_subject_id)

# combine estimate and recovery
param <- cbind(e_param,temp[,-1])#,stanParams[,-1])
param$group <- recode(param$phenotype, "ncv"="psychics","hc"="controls",
                      "nc"="controls","szc"="psychosis","sze"="psychosis",
                      "szl"="psychosis")
e_pred_error$group <- recode(e_pred_error$phenotype, "ncv"="psychics","hc"="controls",
                      "nc"="controls","szc"="psychosis","sze"="psychosis",
                      "szl"="psychosis")

# how is the sample per site
table(param$phenotype,param$site)



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
if (!require(MASS)) install.packages("MASS"); library(MASS)
fitdistr(f_logit(param$lr_wm,0), "normal")
fitdistr(f_logit(param$gamma_wm,0), "normal")
fitdistr(f_logit(param$eta_wm,0), "normal")
fitdistr(atanh(param$lambda_wm), "normal")
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #



# range(f_logit(param$lr_wm,0));range(f_logit(param$lr_rec,0))
rango <- range(f_logit(c(param$lr_wm,param$lr_rec),0))
rango[1] <- floor(rango[1]); rango[2] <- ceiling(rango[2])
p_alpha <- ggplot(param, aes(x=f_logit(lr_wm,0),y=f_logit(lr_rec,0),
                  col=site)) +
  labs(title = "learning rate",
       subtitle = "more = more learning",
       x = expression(alpha), y= expression(hat(alpha))) +
  coord_cartesian(xlim = rango, ylim = rango) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point(alpha=0.5) +
  geom_smooth(method="lm", se = F) +
  stat_cor(method = "pearson", 
           label.x.npc = 0,
           label.y.npc = 1,
           size=2) +
  theme_classic()

# range(log(param$tau_wm));range(log(param$tau_rec))
# range(f_logit(param$tau_wm,0));range(f_logit(param$tau_rec,0))
rango <- range(f_logit(c(param$gamma_wm,param$gamma_rec),0))
rango[1] <- floor(rango[1]); rango[2] <- ceiling(rango[2])
p_gamma <- ggplot(param, aes(x=f_logit(gamma_wm,0),y=f_logit(gamma_rec,0),
                             col=site)) +
  labs(title = "additivity",
       subtitle = "more = more additivity",
       x = expression(gamma), y= expression(hat(gamma))) +
  coord_cartesian(xlim = rango, ylim = rango) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point(alpha=0.5) +
  geom_smooth(method="lm", se = F) +
  stat_cor(method = "pearson", 
           label.x.npc = 0,
           label.y.npc = 1,
           size=2) +
  theme_classic()
  
# range(f_logit(param$eta_wm,0));range(f_logit(param$eta_rec,0))
rango <- range(f_logit(c(param$eta_wm,param$eta_rec),0))
rango[1] <- floor(rango[1]); rango[2] <- ceiling(rango[2])
p_eta <- ggplot(param, aes(x=f_logit(eta_wm,0),y=f_logit(eta_rec,0),
                col=site)) + 
  labs(title = "new phase, w decrement",
       subtitle = "more = more forgetting",
       x = expression(eta), y= expression(hat(eta))) +
  geom_abline(slope = 1, intercept = 0) +
  coord_cartesian(xlim = rango, ylim = rango) +
  geom_point(alpha=0.5) +
  geom_smooth(method="lm", se = F) +
  stat_cor(method = "pearson",
           label.x.npc = 0,
           label.y.npc = 1,
           size=2) +
  theme_classic()

# range(log(param$invT_wm));range(log(param$invT_rec))
rango <- range(atanh(c(param$lambda_wm,param$lambda_rec)))
rango[1] <- floor(rango[1]); rango[2] <- ceiling(rango[2])
p_lambda <- ggplot(param, aes(x=atanh(lambda_wm),y=atanh(lambda_rec),
                            col=site)) + 
  labs(title = "counterfactual",
       subtitle = "less = downward counterfactual update",
       x = expression(lambda), y= expression(hat(lambda))) +
  geom_abline(slope = 1, intercept = 0) +
  coord_cartesian(xlim = rango, ylim = rango) +
  geom_point(alpha=0.5) +
  geom_smooth(method="lm", se = F) +
  stat_cor(method = "pearson", 
           label.x.npc = 0,
           label.y.npc = 1,
           size=2) +
  theme_classic()

figS1 <- ggarrange(p_alpha,p_gamma,p_eta,p_lambda, labels=LETTERS[1:4],
                   common.legend = T)
figS1
if (print_figures == 1) {
  ggsave(filename = "figures/figS1.pdf",plot = figS1,dpi = 600,units = "in",
         width = 6,height = 6)
}



# # # Validation, r_hat against predictions # # #
validation <- lf[substr(lf$trial_type,1,1)=="B"|substr(lf$trial_type,1,1)=="D",]
validation$condition <- ifelse(substr(validation$trial_type,1,1)=="B","Blocking","Control")
figS2 <- ggplot(validation, aes(x=prediction,y=r_hat)) + 
  labs(title = "Model Blocking Validation", 
       y = expression(Model~Prediction~`(`*hat(r)*`)`), 
       x = "Participants' Prediction") +
  geom_point(alpha = 0.1) + 
  geom_smooth(method = "lm", se = F, col = "black") +
  facet_grid(.~condition) + 
  stat_cor() +
  theme_classic()
figS2

if (print_figures == 1) {
  ggsave(filename = "figures/figS2.pdf",plot = figS4, dpi = 600, units = "in",
         width = 6, height = 4)
}