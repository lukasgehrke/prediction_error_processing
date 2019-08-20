# bei mir funktioniert so alles korrekt bis zu den ### (Zeile 75)! Die eingelesenen Werte in R sind korrekt auch wenn Sie im csv anders dargestellt werden.
# Vielleicht darf man im csv keine Bearbeitungen vornehmen da dann die bearbeiteten Felder einen falschen Typ haben, vielleicht...

# read in csv data
erp1 <- read.csv(file="P:\\Project_Sezen\\data_processing\\ERP\\lg_confirmed\\peaks_locs_ERP_normal_FCz_100.csv", header=TRUE, sep=";")
erp1$Congruency <- "normal"

erp2 <- read.csv(file="P:\\Project_Sezen\\data_processing\\ERP\\lg_confirmed\\peaks_locs_ERP_conflict_FCz_100.csv", header=TRUE, sep=";")
erp2$Congruency <- "conflict"

# merge two data frames
erp <- rbind(erp1, erp2)
head(erp)

## make factors
erp$Participant <- as.factor(erp$Participant)
# sensory congruence condition: conflict or normal
erp$Condition <- as.factor(erp$Condition)
# sensory feedback condition: vis or vibro
erp$Congruency <- as.factor(erp$Congruency) 

# sum up means of study design for overview
aggregate(. ~ Condition + Congruency, erp, mean)

# linear model of regression on Min Peak Amplitude
lm_model_FCz_amp <- lm(Min_Peak_Amplitude ~ Condition*Congruency, data = erp)
summary(aov(lm_model_FCz_amp))

library(lsr)
etaSquared(lm_model_FCz_amp, type = 2, anova = TRUE)
# es interessiert nur die Spalte mit eta.sq.part (partielles eta quadrat)

# Boxplot
# specify labs titles, colors
xlab_title <- 'feedback condition'
ylab_title <- 'min peak amplitude (?V)'
colors <- c("#C2D7FF", "#CC7F72")
require(ggpubr)
erp_min_peak <- ggplot(erp, aes(x = Condition,
                                y = Min_Peak_Amplitude,
                                fill = Congruency,
                                colour = Congruency,
                                shape = Congruency)) +
  stat_summary(fun.y=median,
               geom="line", 
               position=position_nudge(x = c(-.25, -.25, .25, .25), y = c(0,0)),
               aes(group=Congruency),
               size = 1.5) + 
  ggbeeswarm::geom_quasirandom(aes(fill=Congruency),
                               alpha = .6,
                               shape = 22,
                               colour = "grey10",
                               #fill = colors,
                               size = 2.3,
                               width = .1) +
  geom_boxplot(aes(colour=Congruency, fill = Congruency),
               width = .3,
               position=position_dodge(1),
               alpha=.2,
               outlier.shape = NA) +
  theme_bw() +
  theme(axis.title.y = element_text(colour="grey20",size=12),
        axis.title.x = element_text(colour="grey20",size=12),
        axis.text.x = element_text(colour="grey20",size=10),
        axis.text.y = element_text(colour="grey20",size=10)) +
  ylab(ylab_title) +
  xlab(xlab_title) +
  scale_fill_manual(values=colors, name="congruency") +
  scale_shape_manual(values=colors, name="congruency") +
  scale_colour_manual(values=colors, name="congruency") +
  scale_x_discrete(limits=c("visual", "vibro")) +
  theme(legend.background = element_rect(size=0.5, linetype="solid"))
erp_min_peak


# 
# ### 
# # linear model of regression on Min Peak Latency
# lm_model_FCz_lat <- lm(erp$Min_Peak_Latency~erp$Congruency*erp$Feedback)
# require(lsr)
# summary(aov(lm_model_FCz_lat))
# 
# # FCz                         Df Sum Sq Mean Sq F value Pr(>F)  
# # erp$Congruency               1   7050    7050   2.937 0.0909 .
# # erp$Feedback                 1   8337    8337   3.473 0.0664 .
# # erp$Congruency:erp$Feedback  1   1122    1122   0.467 0.4964  
# # Residuals                   72 172822    2400                 
# # ---
# #   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# # Pz                          Df Sum Sq Mean Sq F value Pr(>F)  
# # erp$Congruency               1   7050    7050   2.937 0.0909 .
# # erp$Feedback                 1   8337    8337   3.473 0.0664 .
# # erp$Congruency:erp$Feedback  1   1122    1122   0.467 0.4964  
# # Residuals                   72 172822    2400                 
# # ---
# #   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# # Cz
# 
# 
# etaSquared(lm_model_FCz_lat, type = 2, anova = TRUE)
# 
# # FCz                              eta.sq eta.sq.part         SS df       MS         F          p
# # erp$Congruency              0.037238011 0.039196248   7050.316  1 7050.316 2.9372594 0.09085994
# # erp$Feedback                0.044034235 0.046020651   8337.053  1 8337.053 3.4733318 0.06644324
# # erp$Congruency:erp$Feedback 0.005925568 0.006449754   1121.895  1 1121.895 0.4673969 0.49638153
# # Residuals                   0.912802186          NA 172821.895 72 2400.304        NA         NA
# 
# # Pz                               eta.sq eta.sq.part         SS df       MS         F          p
# # erp$Congruency              0.037238011 0.039196248   7050.316  1 7050.316 2.9372594 0.09085994
# # erp$Feedback                0.044034235 0.046020651   8337.053  1 8337.053 3.4733318 0.06644324
# # erp$Congruency:erp$Feedback 0.005925568 0.006449754   1121.895  1 1121.895 0.4673969 0.49638153
# # Residuals                   0.912802186          NA 172821.895 72 2400.304        NA         NA
# 
# # Cz
# 
# # TO DO lUKAS ON 30.04.19 --> done, but ERP values have huuuuge variance... why? check matlab script!
# 
# ### plot taken over from erp_results_CHI.R 
# 
# # post hoc significance testing
# # leave out here
# # require(ggpubr)
# # my_comparisons <- list(c("visual", "vibro"))
# 
# # Boxplot
# # specify labs titles, colors
# xlab_title <- 'feedback condition'
# ylab_title <- 'negative peak amplitude (?V)'
# colors <- c("#C2D7FF", "#CC7F72")
# 
# require(ggpubr)
# 
# 
# erp_min_peak <- ggplot(erp, aes(x = Feedback,
#                                 y = Min_Peak_Amplitude,
#                                 fill = Congruency,
#                                 colour = Congruency,
#                                 shape = Congruency)) +
#   stat_summary(fun.y=median,
#                geom="line", 
#                position=position_nudge(x = c(-.25, -.25, .25, .25), y = c(0,0)),
#                aes(group=Congruency),
#                size = 1.5) + 
#   ggbeeswarm::geom_quasirandom(aes(fill=Congruency),
#                                alpha = .6,
#                                shape = 22,
#                                colour = "grey10",
#                                #fill = colors,
#                                size = 2.3,
#                                width = .1) +
#   geom_boxplot(aes(colour=Congruency, fill = Congruency),
#                width = .3,
#                position=position_dodge(1),
#                alpha=.2,
#                outlier.shape = NA) +
#   theme_bw() +
#   theme(axis.title.y = element_text(colour="grey20",size=12),
#         axis.title.x = element_text(colour="grey20",size=12),
#         axis.text.x = element_text(colour="grey20",size=10),
#         axis.text.y = element_text(colour="grey20",size=10)) +
#   ylab(ylab_title) +
#   xlab(xlab_title) +
#   scale_fill_manual(values=colors, name="congruency") +
#   scale_shape_manual(values=colors, name="congruency") +
#   scale_colour_manual(values=colors, name="congruency") +
#   scale_x_discrete(limits=c("visual", "vibro")) +
#   theme(legend.background = element_rect(size=0.5, linetype="solid"))
#   #theme(legend.position = c(0.5, 0.85)) + 
#   #stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") +
# # erp_min_peak
# 
# ylab_title <- 'negative peak latency (ms)'
# erp_min_lat <- ggplot(erp, aes(x = Feedback,
#                                 y = Min_Peak_Latency,
#                                 fill = Congruency,
#                                 colour = Congruency,
#                                 shape = Congruency)) +
#   stat_summary(fun.y=median,
#                geom="line", 
#                position=position_nudge(x = c(-.25, -.25, .25, .25), y = c(0,0)),
#                aes(group=Congruency),
#                size = 1.5) + 
#   ggbeeswarm::geom_quasirandom(aes(fill=Congruency),
#                                alpha = .6,
#                                shape = 22,
#                                colour = "grey10",
#                                #fill = colors,
#                                size = 2.3,
#                                width = .1) +
#   geom_boxplot(aes(colour=Congruency, fill = Congruency),
#                width = .3,
#                position=position_dodge(1),
#                alpha=.2,
#                outlier.shape = NA) +
#   theme_bw() +
#   theme(axis.title.y = element_text(colour="grey20",size=12),
#         axis.title.x = element_text(colour="grey20",size=12),
#         axis.text.x = element_text(colour="grey20",size=10),
#         axis.text.y = element_text(colour="grey20",size=10)) +
#   ylab(ylab_title) +
#   xlab(xlab_title) +
#   scale_fill_manual(values=colors, name="congruency") +
#   scale_shape_manual(values=colors, name="congruency") +
#   scale_colour_manual(values=colors, name="congruency") +
#   scale_x_discrete(limits=c("visual", "vibro")) +
#   theme(legend.background = element_rect(size=0.5, linetype="solid"))
#   #theme(legend.position = c(0.5, 0.85)) + 
#   #stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") +
# 
# require(cowplot)
# two_plot <- plot_grid(erp_min_peak, erp_min_lat, labels=c("AUTO"), nrow=1, ncol = 2, align = "hv")
# two_plot
