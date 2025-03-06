rm(list = ls())

# import libraries
if (!require(reshape2)) {install.packages("reshape2")}; library(reshape2)
if (!require(ggplot2)) {install.packages("ggplot2")}; library(ggplot2)
if (!require(ggpubr)) {install.packages("ggpubr")}; library(ggpubr)
if (!require(lmerTest)) {install.packages("lmerTest")}; library(lmerTest)

# read long format
lf <- read.csv("data/longFormat.csv")
wf <- read.csv("data/wideFormat.csv")
deltas <- read.csv("data/lfDeltas.csv")

# general parameters
avg_size <- 3
ind_size <- 1.5
avg_stroke <- 0.9
sem_width <- 0.35
y_text_size <- 12
sig_line_size <- 0.5
sig_text_size <- 6

# general theme
theme_mine <- theme_classic()

# how many have unreality?
table(!is.na(wf$pdi40_i10_unreality),wf$group)

# remove participants with no reality item
wf <- wf[!is.na(wf$pdi40_i10_unreality),]
lf <- lf[!is.na(lf$pdi40_i10_unreality),]
deltas <- deltas[!is.na(deltas$pdi40_i10_unreality),]

# how many have unreality? (after removing ids with no item 10)
table(!is.na(wf$pdi40_i10_unreality), wf$group)






# # # # # # # # # # # Unreality (PDI-40)# # # # # # # # # # # # # # # # # # #### 
# melt unreality points
for_plot <- reshape2::melt(wf, measure.vars = c("pdi40_i10_dis_unreality",
                                                "pdi40_i10_pre_unreality",
                                                "pdi40_i10_con_unreality"))
for_plot$variable <- factor(for_plot$variable, levels = c("pdi40_i10_dis_unreality",
                                                          "pdi40_i10_pre_unreality",
                                                          "pdi40_i10_con_unreality"))

# unreality item 10 from PDI-10 is
# Do things around you ever feel unreal, as though it was all part of an experiment?

# plot labels
anno <- data.frame(x1 = c(2,2), x2 = c(3,3), y1 = c(4.3,4.8), y2 = c(4.4,4.9), 
                   xstar = c(2.5,2.5), ystar = c(4.4,4.9), lab = c("*","*"), 
                   variable = factor(c("pdi40_i10_dis_unreality",
                                       "pdi40_i10_pre_unreality"),
                                     levels = c("pdi40_i10_dis_unreality",
                                                "pdi40_i10_pre_unreality",
                                                "pdi40_i10_con_unreality")))

fig_unreality <- ggplot(for_plot[,], 
                        aes(x=group, y=value, col=group, shape=group)) + 
  labs(subtitle = "'Do things around you ever feel unreal,\nas though it was all part of an experiment?'\nIf Yes...",#Unreality item (from PDI-40)
       col="Group",shape="Group",fill="Group",
       x="Group",y="Intensity") +
  stat_summary(geom="bar", fill="white") +
  geom_jitter(alpha = 0.2, size = ind_size, fill="white", width = 0.1, height = 0.1) +
  stat_summary(fun.data = mean_cl_boot, geom="errorbar", width = sem_width) +
  stat_summary(geom="point", size = avg_size, stroke = avg_stroke, fill="white") +
  scale_colour_manual(values = c("orange","purple","red"),
                      labels = c("HC (n=42)",
                                 "EEB (n=31)",
                                 "SZ (n=62)")) +
  scale_shape_manual(values = c(21,22,23),
                     labels = c("HC (n=42)",
                                "EEB (n=31)",
                                "SZ (n=62)")) +
  # scale_colour_manual(values = c("orange","purple","red"),
  #                     labels = c("Healthy Controls (HC; n=42)",
  #                                "Elevated Esoteric Beliefs (EEB; n=31)",
  #                                "Schizophrenia (SZ; n=62)")) +
  # scale_shape_manual(values = c(21,22,23),
  #                    labels = c("Healthy Controls (HC; n=42)",
  #                               "Elevated Esoteric Beliefs (EEB; n=31)",
  #                               "Schizophrenia (SZ; n=62)")) +
  scale_x_discrete(labels = c("HC","EEB","SZ")) +
  geom_segment(data=anno, aes(x=x1,xend=x2,y=y2,yend=y2),inherit.aes=F, size=sig_line_size) +
  geom_segment(data=anno, aes(x=x1,xend=x1,y=y1,yend=y2),inherit.aes=F, size=sig_line_size) +
  geom_segment(data=anno, aes(x=x2,xend=x2,y=y1,yend=y2),inherit.aes=F, size=sig_line_size) +
  geom_text(data = anno, aes(x = xstar,  y = ystar, label = lab), inherit.aes=F, size = sig_text_size) +
  scale_y_continuous(breaks = c(1,3,5)) +
  coord_cartesian(ylim=c(0.9,5.1)) +
  facet_grid(.~variable,
             labeller = labeller(variable= c(`pdi40_i10_dis_unreality` = "Distress",
                                             `pdi40_i10_pre_unreality` = "Preoccupation",
                                             `pdi40_i10_con_unreality` = "Conviction"))) + 
  theme_mine + 
  theme(plot.margin = margin(t=5.5, l=13, r=25, b=5.5),
        legend.position = "none",
        legend.background = element_blank(),
        legend.title = element_blank(),
        axis.title.x = element_blank())
# fig_unreality


# Wilcoxon Rank-Sum test for unpaired samples 
# statistical analysis
temp <- for_plot[for_plot$variable == "pdi40_i10_dis_unreality" & 
                   for_plot$group != "controls" & !is.na(for_plot$value),]
wilcox.test(temp$value[temp$group=="psychosis"],
            temp$value[temp$group=="psychics"], exact = F, paired = F, conf.int = F)
report::report_table(t.test(temp$value[temp$group=="psychosis"],temp$value[temp$group=="psychics"]))

temp <- for_plot[for_plot$variable == "pdi40_i10_pre_unreality" & 
                   for_plot$group != "controls" & !is.na(for_plot$value),]
wilcox.test(temp$value[temp$group=="psychosis"],
            temp$value[temp$group=="psychics"], exact = F, paired = F, conf.int = F)
report::report_table(t.test(temp$value[temp$group=="psychosis"],temp$value[temp$group=="psychics"]))

temp <- for_plot[for_plot$variable == "pdi40_i10_con_unreality" & 
                   for_plot$group != "controls" & !is.na(for_plot$value),]
wilcox.test(temp$value[temp$group=="psychosis"],
            temp$value[temp$group=="psychics"], exact = F, paired = F, conf.int = F)
report::report_table(t.test(temp$value[temp$group=="psychosis"],
                            temp$value[temp$group=="psychics"]))





# # # # # # # # # # # Brief Psychiatric Rating Scale# # # # # # # # # # # # #### 
# melt BRPS scores
for_plot <- reshape2::melt(wf, measure.vars = c("bprs_dis","bprs_neg","bprs_pos"))
# remove absent scores
for_plot <- for_plot[!is.na(for_plot$value),]

# plot labels
anno <- data.frame(x1 = c(1,1,1), x2 = c(3,3,3), y1 = c(4.1,4.3,4.5), y2 = c(4.2,4.4,4.6), 
                   xstar = c(2,2,2), ystar = c(4.3,4.5,4.7), lab = c("***","***","***"), 
                   variable = factor(c("bprs_dis","bprs_neg","bprs_pos"),
                                     levels = c("bprs_dis","bprs_neg","bprs_pos")),
                   pdi40_i10_unreality = factor(c(0,0,0),levels = c(0,1)))

fig_bprs <- ggplot(for_plot, aes(x=group, y=value, col=group, shape=group)) + 
  labs(subtitle = "Symptoms by Clinician (BPRS)",
       #subtitle = "Symptoms assessed by a clinician:\n Brief Psychiatric Rating Scale (BPRS)",
       col="Group",shape="Group",fill="Group",
       x="Group",y="Score") +
  stat_summary(geom="bar", fill="white") +
  geom_jitter(alpha = 0.2, size = ind_size, fill="white", width = 0.1, height = 0.1) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = sem_width) +
  stat_summary(geom = "point", size = avg_size, stroke = avg_stroke, fill = "white") +
  # scale_color_manual(
  #   labels = c("Healthy Controls (HC; n=42)",
  #              "Elevated Esoteric Beliefs (EEB; n=31)",
  #              "Schizophrenia (SZ; n=62)"),
  #   values = c("orange","purple","red")) +
  # scale_shape_manual(
  #   labels = c("Healthy Controls (HC; n=42)",
  #              "Elevated Esoteric Beliefs (EEB; n=31)",
  #              "Schizophrenia (SZ; n=62)"),
  #   values = c(21,22,23)) +
  scale_colour_manual(values = c("orange","purple","red"),
                      labels = c("HC (n=42)",
                                 "EEB (n=31)",
                                 "SZ (n=62)")) +
  scale_shape_manual(values = c(21,22,23),
                     labels = c("HC (n=42)",
                                "EEB (n=31)",
                                "SZ (n=62)")) +
  geom_segment(data=anno, aes(x=x1,xend=x2,y=y2,yend=y2),inherit.aes=F, size=sig_line_size) +
  geom_segment(data=anno, aes(x=x1,xend=x1,y=y1,yend=y2),inherit.aes=F, size=sig_line_size) +
  geom_segment(data=anno, aes(x=x2,xend=x2,y=y1,yend=y2),inherit.aes=F, size=sig_line_size) +
  geom_text(data = anno, aes(x = xstar,  y = ystar, label = lab), inherit.aes=F, size = sig_text_size) +
  scale_x_discrete(labels = c("HC","EEB","SZ")) + 
  scale_y_continuous(breaks = c(1,3,5)) +
  coord_cartesian(ylim=c(1,5)) +
  facet_grid(pdi40_i10_unreality~variable,
             labeller = labeller(variable= c(`bprs_dis` = "Disorganized",
                                             `bprs_neg` = "Negative",
                                             `bprs_pos` = "Positive"),
                                 pdi40_i10_unreality = c(`0` = "Unreality: No",
                                                         `1` = "Unreality: Yes"))) + 
  theme_mine +
  theme(legend.position = c(0.3,0.3),#"none",
        legend.background = element_blank(),
        legend.title = element_blank(),
        axis.title.x = element_blank())
# fig_bprs


# Non-parameter tests: Wilcoxon Rank-Sum test for unpaired samples 
shapiro.test(wf$bprs_dis); shapiro.test(wf$bprs_neg); shapiro.test(wf$bprs_pos)
# statistical analysis
# none unreality in bprs distress
kruskal.test(bprs_dis~group,wf[wf$pdi40_i10_unreality==0,])
wilcox.test(bprs_dis~group,wf[wf$pdi40_i10_unreality==0 & wf$group!="psychosis",],exact = F, paired = F)
# none unreality in bprs negative
kruskal.test(bprs_neg~group,wf[wf$pdi40_i10_unreality==0,])
wilcox.test(bprs_neg~group,wf[wf$pdi40_i10_unreality==0 & wf$group!="psychosis",],exact = F, paired = F)
# none unreality in bprs positive
kruskal.test(bprs_pos~group,wf[wf$pdi40_i10_unreality==0,])
wilcox.test(bprs_pos~group,wf[wf$pdi40_i10_unreality==0 & wf$group!="controls",],exact = F, paired = F)

# none unreality
# unreality in bprs distress
wilcox.test(bprs_dis~group,wf[wf$pdi40_i10_unreality==1,],exact = F, paired = F, conf.int = F)
report::report_table(t.test(bprs_dis~group,wf[wf$pdi40_i10_unreality==1,]))
# unreality in bprs negative
wilcox.test(bprs_neg~group,wf[wf$pdi40_i10_unreality==1,],exact = F, paired = F, conf.int = F)
report::report_table(t.test(bprs_neg~group,wf[wf$pdi40_i10_unreality==1,]))
# unreality in bprs positive
wilcox.test(bprs_pos~group,wf[wf$pdi40_i10_unreality==1,],exact = F, paired = F, conf.int = F)
report::report_table(t.test(bprs_pos~group,wf[wf$pdi40_i10_unreality==1,]))





# # # # # # # # # # # Schneider's First Rank Symptoms # # # # # # # # # # # #### 
# melt unreality points
for_plot <- wf[!is.na(wf$pdi40_i10_unreality),]

# plot labels
anno <- data.frame(x1 = c(2), x2 = c(3), y1 = c(84), y2 = c(88), 
                   xstar = c(2.5), ystar = c(92), lab = c("***"), 
                   pdi40_i10_unreality = factor(c("0"),levels = c("0","1")))

fig_frs <- ggplot(for_plot[,], 
                  aes(x=group, y=frs_perc, col=group, shape=group)) + 
  labs(subtitle = "Schneider's First Rank Symptoms (FRS)", #subtitle = "Schneider's First Rank Symptoms,\n Group, and Unreality",
       col="Group",shape="Group",fill="Group",
       x="Group",y="Percentage (%)") +
  geom_hline(yintercept = 0, col = "black") +
  stat_summary(geom="bar", fill="white") +
  # geom_point(alpha = 0.2, size = ind_size, fill="white") +
  geom_jitter(alpha = 0.2, size = ind_size, fill="white", width = 0.1, height = 0.1) +
  stat_summary(fun.data = mean_cl_boot, geom="errorbar", width = sem_width) +
  stat_summary(geom="point", size = avg_size, stroke = avg_stroke, fill="white") +
  scale_color_manual(
    labels = c("Healthy Controls (HC; n=42)",
               "Elevated Esoteric Beliefs (EEB; n=31)",
               "Clinical Psychosis (SZ; n=62)"),
    values = c("orange","purple","red")) +
  scale_shape_manual(
    labels = c("Healthy Controls (HC; n=42)",
               "Elevated Esoteric Beliefs (EEB; n=31)",
               "Clinical Psychosis (SZ; n=62)"),
    values = c(21,22,23)) +
  geom_segment(data=anno, aes(x=x1,xend=x2,y=y2,yend=y2),inherit.aes=F, size=sig_line_size) +
  geom_segment(data=anno, aes(x=x1,xend=x1,y=y1,yend=y2),inherit.aes=F, size=sig_line_size) +
  geom_segment(data=anno, aes(x=x2,xend=x2,y=y1,yend=y2),inherit.aes=F, size=sig_line_size) +
  geom_text(data = anno, aes(x = xstar,  y = ystar, label = lab), inherit.aes=F, size = sig_text_size) +
  scale_x_discrete(labels = c("HC","EEB","SZ")) +
  scale_y_continuous(breaks = c(0,50,100)) +
  coord_cartesian(ylim=c(0,100)) +
  facet_grid(.~pdi40_i10_unreality,
             labeller = labeller(
               pdi40_i10_unreality = c(`0` = "Unreality: No",
                                       `1` = "Unreality: Yes"))) + 
  theme_mine +
  theme(legend.background = element_blank(),
        legend.title = element_blank(),
        axis.title.x = element_blank())
fig_legend <- get_legend(fig_frs)
fig_frs <- fig_frs + theme(legend.position = "none")
# fig_frs


# Non-parameter tests: Wilcoxon Rank-Sum test for unpaired samples 
shapiro.test(wf$frs_perc)
# statistical analysis
# psychics vs psychosis in unreality
wilcox.test(wf$frs_perc[wf$group=="psychics" & wf$pdi40_i10_unreality==1],
            wf$frs_perc[wf$group=="psychosis" & wf$pdi40_i10_unreality==1],exact = F, paired = F, conf.int = F)
report::report_table(t.test(wf$frs_perc[wf$group=="psychics" & wf$pdi40_i10_unreality==1],
                            wf$frs_perc[wf$group=="psychosis" & wf$pdi40_i10_unreality==1]))
# psychics vs psychosis in none unreality
wilcox.test(wf$frs_perc[wf$group=="psychics" & wf$pdi40_i10_unreality==0],
            wf$frs_perc[wf$group=="psychosis" & wf$pdi40_i10_unreality==0],exact = F, paired = F, conf.int = F)
report::report_table(t.test(wf$frs_perc[wf$group=="psychics" & wf$pdi40_i10_unreality==0],
                            wf$frs_perc[wf$group=="psychosis" & wf$pdi40_i10_unreality==0]))





# # # # # # # # # # # Blocking Behaviour# # # # # # # # # # # # # # # # # # #### 
# filter first trials from last phase
for_plot <- lf[lf$trial_block == 1 & (substr(lf$trial_type,1,1) == "B" |
                                        substr(lf$trial_type,1,1) == "D"),]
# create new variable
for_plot$comparison <- ifelse(substr(for_plot$trial_type,1,1)=="B","B","D")
# add qualitatively Response
for_plot$Response <- ifelse(for_plot$response == 1, "Allergy",
                            ifelse(for_plot$response == -1, "No Allergy",NA))
# add condition name
for_plot$Condition <- ifelse(for_plot$comparison=="B","Blocking","Control")
# remove rows without response
for_plot <- for_plot[!is.na(for_plot$Response),]
# remove absent unreality items
for_plot <- for_plot[!is.na(for_plot$pdi40_i10_unreality),]

# plot labels
anno <- data.frame(x1 = c(1), x2 = c(2), y1 = c(0.45), y2 = c(0.47), 
                   xstar = c(1.5), ystar = c(.485), lab = c("***"), 
                   group = "controls", pdi40_i10_unreality = 0)

fig_blocking <- ggplot(for_plot, 
               aes(x=comparison,y=prediction,col=group,shape=group)) + 
  labs(subtitle = "Kamin Blocking Behaviour",
       x = "Food (B: Blocking, D: Control)",
       y = "Allergy Prediction",
       col = "Group:",shape = "Group:") +
  geom_hline(yintercept = 0, col = "black") +
  stat_summary(geom="bar", fill="white") +
  stat_summary(fun.data = mean_cl_boot, geom="errorbar", width = sem_width) +
  stat_summary(geom="point", size = avg_size, stroke = avg_stroke, fill="white") +
  scale_color_manual(labels = c("HC (n=42)","EEB (n=31)","SZ (n=62)"),
                     values = c("orange","purple","red")) +
  scale_shape_manual(labels = c("HC (n=42)","EEB (n=31)","SZ (n=62)"),
                     values = c(21,22,23)) +
  scale_y_continuous(breaks = c(-0.2,0,0.5)) +
  scale_x_discrete(labels = c("Blocking","Control")) + 
  coord_cartesian(ylim=c(-0.25,0.55)) +
  geom_segment(data=anno, aes(x=x1,xend=x2,y=y2,yend=y2),inherit.aes=F, size=sig_line_size) +
  geom_segment(data=anno, aes(x=x1,xend=x1,y=y1,yend=y2),inherit.aes=F, size=sig_line_size) +
  geom_segment(data=anno, aes(x=x2,xend=x2,y=y1,yend=y2),inherit.aes=F, size=sig_line_size) +
  geom_text(data = anno, aes(x = xstar, y = ystar, label = lab), inherit.aes=F, size = sig_text_size) +
  facet_grid(pdi40_i10_unreality~group,
             labeller = labeller(group= c(`controls` = "HC",
                                          `psychics` = "EEB",
                                          `psychosis` = "SZ"),
                                 pdi40_i10_unreality = c(`0` = "Unreality: No",
                                                         `1` = "Unreality: Yes"))) + 
  theme_mine +
  theme(legend.position = c(0.22,0.41),
        legend.background = element_blank(),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1))
# fig_blocking

# linear mixed model predicting 'prediction' trials nested in participants
# statistical analysis
temp <- for_plot[for_plot$group == "controls" & for_plot$pdi40_i10_unreality == 0,]
report::report_table(lmer(prediction ~ comparison + (1|src_subject_id), REML = T, temp))

temp <- for_plot[for_plot$group == "psychics" & for_plot$pdi40_i10_unreality == 0,]
report::report_table(lmer(prediction ~ comparison + (1|src_subject_id), REML = T, temp))

temp <- for_plot[for_plot$group == "psychics" & for_plot$pdi40_i10_unreality == 1,]
report::report_table(lmer(prediction ~ comparison + (1|src_subject_id), REML = T, temp))

temp <- for_plot[for_plot$group == "psychosis" & for_plot$pdi40_i10_unreality == 0,]
report::report_table(lmer(prediction ~ comparison + (1|src_subject_id), REML = T, temp))

temp <- for_plot[for_plot$group == "psychosis" & for_plot$pdi40_i10_unreality == 1,]
report::report_table(lmer(prediction ~ comparison + (1|src_subject_id), REML = T, temp))






# # # # # # # # # # # Prediction Errors' types# # # # # # # # # # # # # # # #### 
# remove absent unreality items
for_plot <- deltas[!is.na(deltas$pdi40_i10_unreality),]
# remove CDs
for_plot <- for_plot[for_plot$trial_type == "AB",]

# plot labels
anno <- data.frame(x1 = c(2.2,2.2), x2 = c(2.8,2.8), 
                   y1 = c(mean(for_plot$pe[for_plot$group == "psychics" & for_plot$pe_type == "primary" & for_plot$pdi40_i10_unreality == 1]),
                          mean(for_plot$pe[for_plot$group == "psychics" & for_plot$pe_type == "secondary" & for_plot$pdi40_i10_unreality == 1])), 
                   y2 = c(mean(for_plot$pe[for_plot$group == "psychosis" & for_plot$pe_type == "primary" & for_plot$pdi40_i10_unreality == 1]),
                          mean(for_plot$pe[for_plot$group == "psychosis" & for_plot$pe_type == "secondary" & for_plot$pdi40_i10_unreality == 1])), 
                   xstar = c(2.7), ystar = c(.39), 
                   lab = c("***",""), pdi40_i10_unreality = 1)

rowSums(table(for_plot$group,for_plot$src_subject_id)!=0)
fig_pe <- ggplot(for_plot,
                aes(x=group, y=pe, col=pe_type, shape=group, fill=pe_type)) + 
  labs(subtitle = "Prediction Errors",
       col = expression(delta~type), fill= expression(delta~type),
       x = "Group", y = expression(delta)) +
  geom_hline(yintercept = 0, col = "black") +
  stat_summary(fun.data = mean_cl_boot, geom="errorbar", width = sem_width) +
  stat_summary(geom="point", size = avg_size, stroke = avg_stroke) +
  scale_shape_manual(name = NULL,labels=c(bquote(primary~delta),bquote(secondary~delta)),
                     values = c(21,22,23)) +
  scale_color_discrete(name = NULL,labels=c("primary","secondary")) +
  scale_fill_discrete(name = NULL,labels=c(bquote(primary~delta),bquote(secondary~delta))) +
  scale_y_continuous(breaks = c(0,0.3,0.6)) +
  scale_x_discrete(labels = c("HC","EEB","SZ")) +
  coord_cartesian(ylim=c(-.05,0.7)) +
  geom_segment(data=anno, aes(x=x1, xend=x2, y=y1, yend=y2),inherit.aes=F, size=sig_line_size) +
  geom_text(data = anno, aes(x = xstar,  y = ystar, label = lab), inherit.aes=F, size = sig_text_size) +
  facet_grid(pdi40_i10_unreality~.,
             labeller = labeller(pdi40_i10_unreality = c(`0` = "Unreality: No",
                                                         `1` = "Unreality: Yes"))) + 
  guides(shape="none", fill="none") +
  theme_mine +
  theme(legend.position = c(0.3,0.15),
        legend.background = element_blank(),
        axis.title.x = element_blank())
# fig_pe

# statistics
shapiro.test(for_plot$pe)

# full models
summary(lmer(pe ~ group * pe_type * pdi40_i10_unreality + (1|src_subject_id),
             REML = T, for_plot[,]))
# with anova and no controls
anova(lmer(pe ~ group * pe_type * pdi40_i10_unreality + (1|src_subject_id), 
           REML = T, for_plot[for_plot$group != "controls",]))
# linear model effect sizes no controls
report::report_table(lmer(pe ~ group * pe_type * pdi40_i10_unreality + (1|src_subject_id), 
                          REML = T, for_plot[for_plot$group != "controls",]))

# linear model effect sizes only unreality
report::report_table(lmer(pe ~ group * pe_type + (1|src_subject_id),
                           REML = T, for_plot[for_plot$pdi40_i10_unreality == 1,]))



# # # # # # # # # # # Task made in Power Point# # # # # # # # # # # # # # # #### 
# get it in power point, hand made

# # # custom layer for ggplot2 # # #
# https://stackoverflow.com/questions/44688623/adding-custom-images-to-ggplot-facets
library(ggpubr); library(png)
# get the task png
task_img <- readPNG("figures/task_figure1B.png")
# create an object
fig_task <- ggplot() + labs(subtitle = "Kamin Blocking Task") + background_image(task_img) + 
  # This ensures that the image leaves some space at the edges
  theme(plot.margin = margin(t=0, l=.5, r=0, b=0, unit = "cm")) # margin(t=0.5, l=1.5, r=1.5, b=0.5, unit = "cm")






# # # # # # # # # # # Extinction Behaviour# # # # # # # # # # # # # # # # # #### 
# filter last phase
for_plot <- lf[(lf$trial_type == "B2-" | lf$trial_type == "D2-"),]
# create new variable
for_plot$comparison <- ifelse(substr(for_plot$trial_type,1,1)=="B","B","D")
# add qualitatively Response
for_plot$Response <- ifelse(for_plot$response == 1, "Allergy",
                            ifelse(for_plot$response == -1, "No Allergy",NA))
# add condition name
for_plot$Condition <- ifelse(for_plot$comparison=="B","Blocking","Control")
# remove rows without response
for_plot <- for_plot[!is.na(for_plot$Response),]
# remove absent unreality items
for_plot <- for_plot[!is.na(for_plot$pdi40_i10_unreality),]

# plot labels
anno <- data.frame(x1 = 6.2, x2 = 6.35, y1 = -.75, y2 = -.2, xstar = 6.5, ystar = -.2+(-.75-(-.15))/2, 
                   lab = "*", pdi40_i10_unreality = 1)

fig_ext <- ggplot(for_plot,aes(x=trial_block, y=prediction, col=group, shape=group)) + 
  labs(subtitle = "Extinction", 
       col = "Group", shape = "Group",
       x = "Blocks", y = "Allergy Prediction") +
  geom_hline(yintercept = 0, col = "black") +
  # stat_summary(geom="line", size = 1.2, position = position_dodge(0.2)) +
  geom_smooth(method = "lm", se = F, size = 1) +
  stat_summary(fun.data = mean_se, geom="errorbar", width=sem_width, position = position_dodge(0.2)) +
  stat_summary(geom="point", size = avg_size, stroke = avg_stroke, fill="white", position = position_dodge(0.2)) +
  scale_color_manual(labels = c("HC (n=42)","EEB (n=31)","SZ (n=62)"),
                     values = c("orange","purple","red")) +
  scale_shape_manual(labels = c("HC (n=42)","EEB (n=31)","SZ (n=62)"),
                     values = c(21,22,23)) +
  scale_x_continuous(breaks = 1:6) +
  scale_y_continuous(breaks = c(-0.7,0,0.3)) +
  coord_cartesian(ylim=c(-0.75,0.35), xlim=c(1,6.5)) +
  geom_segment(data=anno, aes(x=x2, xend=x2, y=y1, yend=y2),inherit.aes=F, size=sig_line_size) +
  geom_segment(data=anno, aes(x=x1, xend=x2, y=y2, yend=y2),inherit.aes=F, size=sig_line_size) +
  geom_segment(data=anno, aes(x=x1, xend=x2, y=y1, yend=y1),inherit.aes=F, size=sig_line_size) +
  geom_text(data = anno, aes(x = xstar,  y = ystar, label = lab), inherit.aes=F, size = sig_text_size) +
  facet_grid(pdi40_i10_unreality ~ .,
             labeller = labeller(pdi40_i10_unreality = c(`0` = "Unreality: No",
                                                         `1` = "Unreality: Yes"))) +
  theme_mine +
  theme(legend.position = c(.6,.5), # c(.2,.1) "none"
        legend.background = element_blank(),
        legend.title = element_blank())
# fig_ext

# statistical analysis
# full stats, did not work, small sample size
m <- lmer(prediction ~ trial_block * group * pdi40_i10_unreality + (1|src_subject_id),
          REML = F, for_plot[for_plot$group != "controls",])
summary(m); anova(m)

# what about the last two blocks (when Extinction is stronger)?
temp <- for_plot[for_plot$pdi40_i10_unreality == 1 & for_plot$trial_block > 4,]
mean(temp$prediction[temp$group =="psychosis"]); sd(temp$prediction[temp$group =="psychosis"])
mean(temp$prediction[temp$group =="psychics"]); sd(temp$prediction[temp$group =="psychics"])
# group to binary
ggplot(temp, aes(x=group,y=prediction)) + stat_summary()

# t.test with the last 2 trials
report::report_table(t.test(prediction~group,temp))

# t.test with the slopes
ggplot(wf[wf$pdi40_i10_unreality==1,], aes(x=group,y=exS)) + stat_summary()
t.test(wf$exS[wf$pdi40_i10_unreality==1 & wf$group=="psychics"],
       wf$exS[wf$pdi40_i10_unreality==1 & wf$group=="psychosis"])

# "cherry on the cake" 
# extinction and control
cor.test(wf$exLast2, wf$chat_control, method = "spearman", exact = F)
summary(lm(exLast2 ~ chat_control * group, wf[wf$pdi40_i10_unreality==1,]))
summary(lm(chat_control ~ group, wf[wf$pdi40_i10_unreality==1,]))

# what about extinction and control. It is significant only for unreality 
summary(lm(exLast2 ~ chat_control, wf[wf$pdi40_i10_unreality==1,]))
cor.test(wf$exLast2[wf$pdi40_i10_unreality==1], 
         wf$chat_control[wf$pdi40_i10_unreality==1], method = "spearman", exact = F)
cor.test(wf$exLast2[wf$pdi40_i10_unreality==1], 
         wf$chat_control[wf$pdi40_i10_unreality==1], method = "pearson", exact = F)
# but nor for not unreality 
summary(lm(exLast2 ~ chat_control, wf[wf$pdi40_i10_unreality==0,]))
cor.test(wf$exLast2[wf$pdi40_i10_unreality==0], 
         wf$chat_control[wf$pdi40_i10_unreality==0], method = "spearman", exact = F)
cor.test(wf$exLast2[wf$pdi40_i10_unreality==0], 
         wf$chat_control[wf$pdi40_i10_unreality==0], method = "pearson", exact = F)
# even with the interaction
summary(lm(exLast2 ~ chat_control * pdi40_i10_unreality, wf[,]))


# differences in correlations
spearman_diff_test <- function(r1, n1, r2, n2) {
  # Fisher r-to-z transformation
  z1 <- 0.5 * log((1 + r1) / (1 - r1))
  z2 <- 0.5 * log((1 + r2) / (1 - r2))
  
  # Standard errors
  SE1 <- 1 / sqrt(n1 - 3)
  SE2 <- 1 / sqrt(n2 - 3)
  
  # Difference and its standard error
  z_diff <- z1 - z2
  SE_diff <- sqrt(SE1^2 + SE2^2)
  
  # z-test statistic
  z <- z_diff / SE_diff
  
  # Two-tailed p-value
  p_value <- 2 * (1 - pnorm(abs(z)))
  
  ci95 <- paste0("dif: ",z_diff, ", from ", z_diff + SE_diff*1.96, " to ",z_diff - SE_diff*1.96)
  
  return(list(z_score = z, p_value = p_value, ci95 = ci95))
}
# differences in correlations unrreality vs none unreality
spearman_diff_test(r1=-0.7487018,n1=12,r2=0.1493659,n2=123)

fig_ext2 <- ggplot(wf, aes(x=chat_control,y=exLast2, col=as.factor(pdi40_i10_unreality))) + 
  labs(x="Control",y="Allergy Prediction\n(last 2 trials)",col="Unreality") +
  geom_hline(yintercept = 0, col = "black") +
  geom_point(alpha=0.2) + 
  geom_smooth(method = "lm",se = F) +
  scale_color_manual(values=scales::hue_pal()(4)[c(2,4)], labels=c("No","Yes")) +
  scale_y_continuous(breaks = c(-1,0,1)) +
  stat_cor(method = "spearman") +
  theme_mine + theme(legend.position = c(.8,.7),
                     legend.background = element_blank())
# fig_ext2

# lambda and extinction
shapiro.test(wf$lambda_wm) # is it normal? no
cor.test(wf$exLast2, wf$lambda_wm, method = "spearman")
cor.test(wf$exLast2, wf$lambda_wm, method = "pearson")
fig_ext3 <- ggplot(wf, aes(x=lambda_wm, y=exLast2)) + 
  labs(x=expression(lambda*`:`~Counterfactual),y="Allergy Prediction\n(last 2 trials)") +
  geom_hline(yintercept = 0, col = "black") +
  geom_point(alpha=0.2) + geom_smooth(method = "lm",col="black",se = F) +
  scale_y_continuous(breaks = c(-1,0,1)) +
  stat_cor(method = "spearman") +
  theme_mine
# fig_ext3






# # # # # # # # # # # model parameters# # # # # # # # # # # # # # # # # # # ####
# melt transform parameters
for_plot <- wf

# normal? no
shapiro.test(wf$chat_control)

# melt
for_plot <- reshape2::melt(for_plot, measure.vars = c("lr_wm","gamma_wm","eta_wm","lambda_wm"))
levels(for_plot$variable) <- c("alpha*`:`~Learning~Rate","gamma*`:`~Additivity","eta*`:`~Forgetting","lambda*`:`~Counterfactual")

# unreality item 10 from PDI-10 is
# Do things around you ever feel unreal, as though it was all part of an experiment?

# plot labels
anno <- data.frame(x1 = c(1,2,2), x2 = c(3,3,3), y1 = c(.35,.65,1.05), y2 = c(.4,.7,1.1),
                   xstar = c(2,2.5,2.5), ystar = c(.45,.75,1.15), lab = c("**","*","*"),
                   variable = factor(c("alpha*`:`~Learning~Rate","gamma*`:`~Additivity","lambda*`:`~Counterfactual"),
                                     levels = c("alpha*`:`~Learning~Rate","gamma*`:`~Additivity","eta*`:`~Forgetting","lambda*`:`~Counterfactual")),
                   unreality = c("Unreality: No","Unreality: Yes","Unreality: Yes"))
# change labels
for_plot$unreality <- ifelse(for_plot$pdi40_i10_unreality==1,"Unreality: Yes",
                             ifelse(for_plot$pdi40_i10_unreality==0,"Unreality: No",NA))

fig_param <- ggplot(for_plot[!is.na(for_plot$pdi40_i10_unreality),], 
                    aes(x=group, y=value, col=group, shape=group)) + 
  labs(subtitle = "Parameters",x="Group",y="Estimates") +
  geom_hline(yintercept = c(-1,0,1), col = "grey") +
  # geom_point(alpha = 0.2, size = ind_size, fill="white") +
  stat_summary(fun.data = mean_cl_boot, geom="errorbar", width = sem_width) +
  stat_summary(geom="point", size = avg_size, stroke = avg_stroke, fill="white") +
  geom_smooth(method = "lm", se = F) +
  scale_color_manual(values = c("orange","purple","red"),
                     labels = c("HC","EEB","SZ")) +
  scale_shape_manual(values = c(21,22,23)) +
  scale_x_discrete(labels = c("HC","EEB","SZ")) +
  scale_y_continuous(breaks = seq(-1,1,by=.5)) +
  coord_cartesian(ylim = c(-0.32,1.2)) +
  facet_grid(unreality~variable, labeller = label_parsed, scales = "free_y") +
  geom_segment(data=anno, aes(x=x1,xend=x2,y=y2,yend=y2),inherit.aes=F, size=sig_line_size) +
  geom_segment(data=anno, aes(x=x1,xend=x1,y=y1,yend=y2),inherit.aes=F, size=sig_line_size) +
  geom_segment(data=anno, aes(x=x2,xend=x2,y=y1,yend=y2),inherit.aes=F, size=sig_line_size) +
  geom_text(data = anno, aes(x = xstar, y = ystar, label = lab), inherit.aes=F, size = sig_text_size) +
  theme_mine + 
  theme(legend.position = "none",
        legend.background = element_blank(),
        legend.title = element_blank(),
        axis.title.x = element_blank())
# fig_param

# statistical analysis

# learning rate, anova, none unreality
report::report_table(anova(lmer(lr_wm~group+(1|site),wf[wf$pdi40_i10_unreality==1,])))
# learning rate, anova, unreality
report::report_table(anova(lmer(lr_wm~group+(1|site),wf[wf$pdi40_i10_unreality==0,])))
# learning rate, linear model, none unreality
report::report_table(lmer(lr_wm~group+(1|site),wf[wf$pdi40_i10_unreality==1,]))
# learning rate, linear model, unreality
report::report_table(lmer(lr_wm~group+(1|site),wf[wf$pdi40_i10_unreality==0,]))

# gamma, anova, none unreality
report::report_table(anova(lmer(gamma_wm~group+(1|site),wf[wf$pdi40_i10_unreality==1,])))
# gamma, anova, unreality
report::report_table(anova(lmer(gamma_wm~group+(1|site),wf[wf$pdi40_i10_unreality==0,])))
# gamma, linear model, none unreality
report::report_table(lmer(gamma_wm~group+(1|site),wf[wf$pdi40_i10_unreality==1,]))
# gamma, linear model, unreality
report::report_table(lmer(gamma_wm~group+(1|site),wf[wf$pdi40_i10_unreality==0,]))

# eta, anova, none unreality
report::report_table(anova(lmer(eta_wm~group+(1|site),wf[wf$pdi40_i10_unreality==1,])))
# eta, anova, unreality
report::report_table(anova(lmer(eta_wm~group+(1|site),wf[wf$pdi40_i10_unreality==0,])))
# eta, linear model, none unreality
report::report_table(lmer(eta_wm~group+(1|site),wf[wf$pdi40_i10_unreality==1,]))
# eta, linear model, unreality
report::report_table(lmer(eta_wm~group+(1|site),wf[wf$pdi40_i10_unreality==0,]))

# lambda, anova, none unreality
report::report_table(anova(lmer(lambda_wm~group+(1|site),wf[wf$pdi40_i10_unreality==1,])))
# lambda, anova, unreality
report::report_table(anova(lmer(lambda_wm~group+(1|site),wf[wf$pdi40_i10_unreality==0,])))
# lambda, linear model, none unreality
report::report_table(lmer(lambda_wm~group+(1|site),wf[wf$pdi40_i10_unreality==1,]))
# lambda, linear model, unreality
report::report_table(lmer(lambda_wm~group+(1|site),wf[wf$pdi40_i10_unreality==0,]))



# PDI paranoia items
anno <- data.frame(x1 = c(1), x2 = c(3), y1 = c(.5,.4), y2 = c(.55,.45),
                   xstar = c(2), ystar = c(.6,.5), lab = c("**","*"),
                   variable = factor(c("alpha*`:`~Learning~Rate","alpha*`:`~Learning~Rate"),
                                     levels = c("alpha*`:`~Learning~Rate","gamma*`:`~Additivity","eta*`:`~Forgetting","lambda*`:`~Counterfactual")),
                   pdiParanoia = factor(c("Paranoia: High","Paranoia: Low"),
                                  levels = c("Paranoia: Low","Paranoia: High")))
# change labels
wf$pdiParanoia <- factor(ifelse(wf$pdi40P > median(wf$pdi40P), "Paranoia: High" ,
                          ifelse(wf$pdi40P <= median(wf$pdi40P), "Paranoia: Low",NA)),
                   levels = c("Paranoia: Low","Paranoia: High"))
for_plot$pdiParanoia <- factor(ifelse(for_plot$pdi40P > median(for_plot$pdi40P), "Paranoia: High" ,
                                ifelse(for_plot$pdi40P <= median(for_plot$pdi40P), "Paranoia: Low",NA)),
                         levels = c("Paranoia: Low","Paranoia: High"))

rowSums(table(for_plot$group[!is.na(for_plot$value)],
              for_plot$src_subject_id[!is.na(for_plot$value)])!=0)
fig_param2 <- ggplot(for_plot[!is.na(for_plot$pdiParanoia),], 
                     aes(x=group, y=value, col=group, shape=group)) + 
  labs(subtitle = "Paranoia from PDI-40",x="Group",y="Parameters Values") +
  geom_hline(yintercept = c(-1,0,1), col = "grey") +
  stat_summary(fun.data = mean_se, geom="errorbar", width = sem_width) +
  stat_summary(geom="point", size = avg_size, stroke = avg_stroke, fill="white") +
  geom_smooth(method = "lm", se = F) +
  scale_color_manual(values = c("orange","purple","red"),
                     labels = c("HC","EEB","SZ")) +
  scale_shape_manual(values = c(21,22,23)) +
  scale_x_discrete(labels = c("HC","EEB","SZ")) +
  scale_y_continuous(breaks = seq(-1,1,by=.5)) +
  coord_cartesian(ylim = c(0,1)) +
  facet_grid(pdiParanoia~variable, labeller = label_parsed, scales = "free_y") +
  geom_segment(data=anno, aes(x=x1,xend=x2,y=y2,yend=y2),inherit.aes=F, size=sig_line_size) +
  geom_segment(data=anno, aes(x=x1,xend=x1,y=y1,yend=y2),inherit.aes=F, size=sig_line_size) +
  geom_segment(data=anno, aes(x=x2,xend=x2,y=y1,yend=y2),inherit.aes=F, size=sig_line_size) +
  geom_text(data = anno, aes(x = xstar, y = ystar, label = lab), inherit.aes=F, size = sig_text_size) +
  theme_mine + 
  theme(plot.subtitle = element_text(face="bold"),
        legend.position = "none",
        legend.background = element_blank(),
        legend.title = element_blank(),
        axis.title.x = element_blank())
# fig_param2

# statistical analysis
report::report_table(anova(lmer(lr_wm~group+(1|site),wf[wf$pdiParanoia=="Paranoia: High",])))
report::report_table(anova(lmer(lr_wm~group+(1|site),wf[wf$pdiParanoia=="Paranoia: Low",])))
report::report_table(lmer(lr_wm~group+(1|site),wf[wf$pdiParanoia=="Paranoia: High",]))
report::report_table(lmer(lr_wm~group+(1|site),wf[wf$pdiParanoia=="Paranoia: Low",]))

report::report_table(anova(lmer(gamma_wm~group+(1|site),wf[wf$pdiParanoia=="Paranoia: High",])))
report::report_table(anova(lmer(gamma_wm~group+(1|site),wf[wf$pdiParanoia=="Paranoia: Low",])))
report::report_table(lmer(gamma_wm~group+(1|site),wf[wf$pdiParanoia=="Paranoia: High",]))
report::report_table(lmer(gamma_wm~group+(1|site),wf[wf$pdiParanoia=="Paranoia: Low",]))

report::report_table(anova(lmer(eta_wm~group+(1|site),wf[wf$pdiParanoia=="Paranoia: High",])))
report::report_table(anova(lmer(eta_wm~group+(1|site),wf[wf$pdiParanoia=="Paranoia: Low",])))
report::report_table(lmer(eta_wm~group+(1|site),wf[wf$pdiParanoia=="Paranoia: High",]))
report::report_table(lmer(eta_wm~group+(1|site),wf[wf$pdiParanoia=="Paranoia: Low",]))

report::report_table(anova(lmer(lambda_wm~group+(1|site),wf[wf$pdiParanoia=="Paranoia: High",])))
report::report_table(anova(lmer(lambda_wm~group+(1|site),wf[wf$pdiParanoia=="Paranoia: Low",])))
report::report_table(lmer(lambda_wm~group+(1|site),wf[wf$pdiParanoia=="Paranoia: High",]))
report::report_table(lmer(lambda_wm~group+(1|site),wf[wf$pdiParanoia=="Paranoia: Low",]))






# # # # # # # # # # # Figure 1# # # # # # # # # # # # # # # # # # # # # # # #### 
fig1_bp <- annotate_figure(ggarrange(ggarrange(fig_unreality,fig_frs,nrow=2,labels = c("A","C")),
                                     fig_bprs,ncol=2,labels = c("","B")),
                     top = text_grob("Questionnaires and Symptoms", color = "black", face = "bold", size = 14,
                                     hjust = 0.5, x = 0.5))
fig1_bp

fig2_bp <- annotate_figure(ggarrange(ggarrange(NULL,fig_task,NULL,nrow=3,heights = c(1,6,1)),
                                     fig_blocking,ncol=2,labels = c("A","B")),
                     top = text_grob("Task and Behaviour", color = "black", face = "bold", size = 14,
                                     hjust = 0.5, x = 0.5))
fig2_bp

fig3_bp <- annotate_figure(ggarrange(fig_param,fig_pe,ncol=2,labels = c("A","B"),widths = c(2,1)),
                           top = text_grob("Computational Model", color = "black", face = "bold", size = 14,
                                           hjust = 0.5, x = 0.5))
fig3_bp

fig4_bp <- annotate_figure(ggarrange(fig_ext,ggarrange(fig_ext2,fig_ext3,nrow=2,labels = c("B","C")),ncol=2,labels = c("A","")),
                           top = text_grob("Exploratory Analysis", color = "black", face = "bold", size = 14,
                                           hjust = 0.5, x = 0.5))
fig4_bp

# print
print_figures <- 1
if (print_figures == 1) {
  ggsave(filename = "figures/fig1_bp.pdf", plot = fig1_bp, dpi = 2400, units = "in",
         width = 6, height = 4, scale = 1.25)
  ggsave(filename = "figures/fig2_bp.pdf", plot = fig2_bp, dpi = 2400, units = "in",
         width = 6, height = 4, scale = 1.25)
  ggsave(filename = "figures/fig3_bp.pdf", plot = fig3_bp, dpi = 2400, units = "in",
         width = 6, height = 4, scale = 1.25)
  ggsave(filename = "figures/fig4_bp.pdf", plot = fig4_bp, dpi = 2400, units = "in",
         width = 6, height = 4, scale = 1.25)
  ggsave(filename = "figures/figS3_bp.pdf", plot = fig_param2, dpi = 2400, units = "in",
         width = 4, height = 4, scale = 1.25)
}







# # # # # # # # # # # Supplementary Figures # # # # # # # # # # # # # # # # #### 
# # # # # # # # # # # Figure# # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # Figure# # # # # # # # # # # # # # # # # # # # # # # # # # 


# new descriptive and statistics for a continuous variable
f_continiousDescriptive <- function(txt,group,conVar,stats) {
  # parameters: txt (common features), group (which variable will be used for grouping),
  # conVar (continious variable to be analized), stats (conduct statistics?)
  
  # to factors
  group <- as.factor(group)
  
  # groups in ascending order
  grp <- levels(group)
  nGrp <- length(grp)
  
  # output matrix
  out <- matrix(NA,nrow=1,ncol=10+nGrp)
  
  # to data frame
  out <- as.data.frame(out)
  
  # combine out with txt
  out[1,1] <- txt[1]
  out[1,2] <- txt[2]
  out[1,3] <- txt[3]
  out[1,4] <- txt[4]
  out[1,5] <- "Continuous; M (SD)"
  
  # add total
  out[1,7] <- paste0(round(mean(conVar,na.rm = T),2), " (",
                     round(sd(conVar,na.rm = T),2),")")
  for (i in 1:nGrp) {
    out[1,7+i] <- paste0(round(mean(conVar[group == grp[i]],na.rm = T),2), " (",
                         round(sd(conVar[group == grp[i]],na.rm = T),2),")")
  }
  
  # # # # only works when length(grp) == 2 # # # #
  # filter 
  filtro <- !is.na(conVar) & !is.na(group)
  # do stats if: stats == 1 AND if all conVar are not NAs AND 
  # if non-NAs are at least 2 groups
  if (stats == 1 & (sum(is.na(conVar)) != length(conVar)) & 
      length(unique(group[filtro])) != 1) {
    # statistical test
    test <- anova(lm(conVar[filtro] ~ group[filtro]))
    
    # combine out with test
    out[1,ncol(out)-2] <- round(test$`F value`[1],6)
    if (!is.null(test$parameter)) {
      out[1,ncol(out)-1] <- test$parameter
    } else {
      out[1,ncol(out)-1] <- paste(test$Df, collapse = ", ")
    }
    out[1,ncol(out)] <- round(test$`Pr(>F)`[1],6)
  }
  
  # column names
  names(out) <- c("Table_name","Char.","Sub. Char.","item","Var. Type","Lvl. Names",
                  "Total",grp,"Statistic","Parameter","p value")
  
  return(out)
}
# new descriptive and statistics for a categorical variable
f_categoricalDescriptive <- function(txt,group,catVar,stats,per_type) {
  # to factors
  catVar <- as.factor(catVar)
  group <- as.factor(group)
  
  # levels in catVar and in group
  names <- levels(catVar)
  nLev <- length(names)
  grp <- levels(group)
  nGrp <- length(grp)
  
  
  
  # # # # frequency and percentages
  # percentages per row
  if (per_type == 1) {
    out <- matrix(paste0(table(catVar,group)," (",
                         round(table(catVar,group)/rowSums(table(catVar,group))*100,1)
                         ,")"),nrow = nLev, ncol = nGrp)
    # total
    out <- cbind(paste0(table(catVar)," (100)"),out)
    
  } else if (per_type == 2) {
    # percentages per column
    out <- matrix(paste0(table(catVar,group)," (",
                         t(round(t(table(catVar,group))/colSums(table(catVar,group))*100,1))
                         ,")"),nrow = nLev, ncol = nGrp)
    # total
    out <- cbind(paste0(table(catVar)," (",round(table(catVar)/sum(table(catVar))*100,1),")"),out)
  } else {
    # frequencies
    out <- matrix(table(catVar,group), nrow = nLev, ncol = nGrp)
    # total
    out <- cbind(matrix(table(catVar)),out)
  }
  
  
  # create space in output matrix
  out <- cbind(matrix(NA,nLev,6),out,matrix(NA,nLev,3))
  
  # combine out with txt
  for (j in 1:4) {
    out[,j] <- txt[j]
  }
  out[,5] <- "Categorical; # (%)"
  out[,6] <- names
  
  if (stats == 1) {
    # statistical test
    test <- chisq.test(table(catVar,group))
    
    # combine out with test
    out[1,ncol(out)-2] <- round(test$statistic,6)
    out[1,ncol(out)-1] <- test$parameter
    out[1,ncol(out)] <- round(test$p.value,6)
  }
  
  # to data frame
  out <- as.data.frame(out)
  
  # column names
  names(out) <- c("Table_name","Char.","Sub. Char.","item","Var. Type","Lvl. Names",
                  "Total",grp,"Statistic","Parameter","p value")
  # out <- out[,c(1:7,9:11,8,12:14)]
  
  # if there are missing data (NAs), create "Lost Data" output
  # if (sum(is.na(catVar)) > 0) {
  out <- rbind(out,NA)
  out[nLev+1,"Lvl. Names"] <- "Missing Data"
  # out[nLev+1,7:(7+nGrp)] <- c(sum((table(group)-colSums(table(catVar,group)))),
  #                             table(group)-colSums(table(catVar,group)))
  for (i in 1:nGrp) {
    if (i == 1) {
      temp <- sum(is.na(group[group == grp[i]|is.na(group)]))
    } else {
      temp <- c(temp,sum(is.na(group[group == grp[i]|is.na(group)])))
    }
  }
  temp <- c(sum(is.na(catVar)),temp)
  out[nLev+1,7:(7+nGrp)] <- temp
  
  for (j in 1:4) {
    out[nLev+1,j] <- txt[j]
  }
  out[nLev+1,5] <- "Categorical; # (%)"
  # }
  return(out)
}
# function to build table given a specific dataset (useful for stratification)
f_buildTable <- function (db, table_variables, group_by, do_stats, per_type) {
  # run the whole table
  for (i in 1:nrow(table_variables)) {
    message(paste("writing table, var:",table_variables$Sub.[i]))
    
    # if group and catVar/conVar are empty then do not run
    if ((sum(table(db[,group_by])) + 
         sum(table(db[,table_variables$item[i]]))) == 0) {
      message(paste("ERROR",table_variables$Sub.[i]))
    # if  catVar/conVar are not empty then run stats
    } else {
      # categorical or continuous
      if (table_variables$cont[i] == 1) {
        # continuous
        temp <- f_continiousDescriptive(txt = c(table_variables$nTable[i],table_variables$Char.[i],
                                                table_variables$Sub.[i],table_variables$item[i]),
                                        group = db[,group_by], conVar = db[,table_variables$item[i]],
                                        stats = do_stats)
      } else {
        # only if not all of them are NAs
        if (sum(is.na(db[,table_variables$item[i]])) != nrow(db)) {
          # categorical
          temp <- f_categoricalDescriptive(txt = c(table_variables$nTable[i],table_variables$Char.[i],
                                                   table_variables$Sub.[i],table_variables$item[i]),
                                           group = db[,group_by], catVar = db[,table_variables$item[i]],
                                           stats = do_stats, per_type = per_type)
        }
      }
      # combine
      if (i == 1) {
        Temp <- temp
      } else {
        Temp <- rbind(Temp,temp)
      }
    }; remove(temp)
  } # end nrow table_variables
  # output 
  return(Temp)
}


feedTheCode <- data.frame(
  t(matrix(c(1, "Gen. Char.", "Age (years)",             "age",              1,
             1, "Gen. Char.", "Sex",                     "curr_sex",         0,
             1, "Gen. Char.", "Duration of Illness",     "duration_illness", 1,
             1, "Gen. Char.", "Tot. Chlorpromazine",     "TotChlorpromazine", 1,
             1, "Gen. Char.", "Years of Education (YoE)","educat_level",     1,
             1, "Gen. Char.", "Race",                    "race",             0,
             1, "Gen. Char.", "Handedness",              "handedness",       0,
             1, "Gen. Char.", "Age at First Episde",     "age_first_ah",     1,
             1, "Gen. Char.", "Mother YoE",              "F_educat_level",   1,
             1, "Gen. Char.", "Father YoE",              "M_educat_level",   1,
             2, "Gen. Char.", "KB score",                "kb_score",         1,
             2, "Gen. Char.", "PDI-40 total",            "pdi40",            1,
             2, "Gen. Char.", "Unreality: PDI-item-10","pdi40_i10_unreality",0,
             2, "Gen. Char.", "BPRS Positive",           "bprs_pos",         1,
             2, "Gen. Char.", "BPRS Negative",           "bprs_neg",         1,
             2, "Gen. Char.", "BPRS Distress",           "bprs_dis",         1,
             2, "Gen. Char.", "FRS Percentaje",          "frs_perc",         1,
             2, "Gen. Char.", "Learning Rate (alpha)" ,  "lr_wm",            1,
             2, "Gen. Char.", "Additivity (gamma)",      "gamma_wm",         1,
             2, "Gen. Char.", "Forgetting (eta)",        "eta_wm",           1,
             2, "Gen. Char.", "Counterfactual (lambda)", "lambda_wm",        1,
             2, "Gen. Char.", "Chat control",            "chat_control",     1),
           nrow=5)))
colnames(feedTheCode) <- c("nTable","Char.","Sub.","item","cont")

# two tables
tabS1feed <- feedTheCode[feedTheCode$nTable == 1,]
tabS2feed <- feedTheCode[feedTheCode$nTable == 2,]

tableS1 <- f_buildTable(db = wf, table_variables = tabS1feed,
                        group_by = "group", do_stats = 1, per_type = 2)
tableS2 <- f_buildTable(db = wf, table_variables = tabS2feed,
                        group_by = "group", do_stats = 1, per_type = 2)


print_csv <- 1
if (print_csv == 1) {
  write.csv(tableS1, "figures/tableS1.csv", row.names = F, na = "")
  write.csv(tableS2, "figures/tableS2.csv", row.names = F, na = "")
}

table(wf$group[wf$handedness!=""],wf$handedness[wf$handedness!=""])
chisq.test(wf$group[wf$handedness!=""],wf$handedness[wf$handedness!=""])


# sex analysis
# unreality and sex
table(wf$pdi40_i10_unreality,wf$curr_sex)
chisq.test(table(wf$pdi40_i10_unreality,wf$curr_sex))
# sex and paranoia
wilcox.test(pdi40~curr_sex,wf)
report::report_table(t.test(pdi40~curr_sex,wf))
# is kb score normally distributed?
shapiro.test(wf$kb_score)
# kb_score as a function of sex
report::report_table(t.test(kb_score~curr_sex,wf))

# group and sex
table(wf$group,wf$curr_sex)
chisq.test(table(wf$group,wf$curr_sex))

# Medication normally distributed
shapiro.test(wf$TotChlorpromazine)
# and counterfactual does not correlates with medication
cor.test(wf$TotChlorpromazine,wf$lambda_wm,method = "spearman")
cor.test(wf$TotChlorpromazine,wf$gamma_wm,method = "spearman")
cor.test(wf$TotChlorpromazine[wf$pdi40_i10_unreality==1],wf$lambda_wm[wf$pdi40_i10_unreality==1],method = "spearman")
cor.test(wf$TotChlorpromazine[wf$pdi40_i10_unreality==1],wf$gamma_wm[wf$pdi40_i10_unreality==1],method = "spearman")

# Intelligence normally distributed
shapiro.test(wf$wtarss)
# counterfactuals does not correlate with intelligence
cor.test(wf$wtarss, wf$lambda_wm,method = "spearman")
cor.test(wf$wtarss, wf$gamma_wm,method = "spearman")
cor.test(wf$wtarss[wf$pdi40_i10_unreality==1],wf$lambda_wm[wf$pdi40_i10_unreality==1],method = "spearman")
cor.test(wf$wtarss[wf$pdi40_i10_unreality==1],wf$gamma_wm[wf$pdi40_i10_unreality==1],method = "spearman")

