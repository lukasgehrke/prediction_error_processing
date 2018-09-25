
# read in csv data
#erp <- read.csv(file="P:\\Project_Sezen\\data_processing\\ERP\\peaks_locs_ERP_incl_EMS_FCz_95.csv", header=TRUE, sep=";")
erp <- read.csv(file="P:\\Project_Sezen\\data_processing\\ERP\\peaks_locs_ERP_incl_EMS_FCz_100.csv", header=TRUE, sep=";")
#erp <- read.csv(file="P:\\Project_Sezen\\data_processing\\ERP\\peaks_locs_ERP_incl_EMS_FCz_105.csv", header=TRUE, sep=";")

nr_of_subjects <- 11

# make factors
erp$Participant <- as.factor(erp$Participant)
erp$Condition <- as.factor(erp$Condition)

# inspect
f1 <- function(x) c(Mean = mean(x), Max = max(x), SD = sd(x))

do.call(data.frame, aggregate(Min_Peak_Amplitude~Condition, erp, f1))
# Condition Min_Peak_Amplitude.Mean Min_Peak_Amplitude.Max Min_Peak_Amplitude.SD
# 1       ems               -6.216436               -2.80250              2.093379
# 2     vibro               -4.670238               -0.34552              2.377441
# 3    visual               -3.983173               -1.04120              1.734263

do.call(data.frame, aggregate(Min_Peak_Latency~Condition, erp, f1))
# Condition Min_Peak_Latency.Mean Min_Peak_Latency.Max Min_Peak_Latency.SD
# 1       ems              176.0000                  228            32.39753
# 2     vibro              183.6364                  276            37.11407
# 3    visual              173.0909                  236            52.46228

# ANOVA on aggregated data Min Peak Amplitude
require(lsr)
fit_amp <- aov(Min_Peak_Amplitude ~ Condition, data=erp)
#summary(fit_amp)
etaSquared(fit_amp, type = 2, anova = TRUE)
# eta.sq eta.sq.part        SS df        MS        F          p
# Condition 0.1807994   0.1807994  28.78427  2 14.392133 3.310533 0.05021739
# Residuals 0.8192006          NA 130.42130 30  4.347377       NA         NA

# ANOVA on aggregated data Min Peak Latency
fit_lat <- aov(Min_Peak_Latency ~ Condition, data=erp)
#summary(fit_lat)
etaSquared(fit_lat, type = 2, anova = TRUE)
# eta.sq eta.sq.part         SS df       MS         F         p
# Condition 0.01244338  0.01244338   652.6061  2  326.303 0.1890025 0.8287625
# Residuals 0.98755662          NA 51793.4545 30 1726.448        NA        NA

### plot
# Boxplot Duration
# specify title, labs titles, legend title
xlab_title <- 'feedback condition'

# post hoc significance testing
require(ggpubr)
my_comparisons <- list(c("visual", "ems"), c("vibro", "ems"))

ylab_title <- 'negative peak amplitude (µV)'
erp_min_peak <- ggplot(erp, aes(x = Condition,
                                        y = Min_Peak_Amplitude, 
                                        fill = Condition,
                                        colour = Condition,
                                        shape = Condition)) +
  #geom_line(aes(group=Participant), colour="grey10", size=.5, alpha=.2) +
  ggbeeswarm::geom_quasirandom(alpha = .4,
                               shape = 22,
                               colour = "grey10",
                               fill = "#185994",
                               size = 2.3,
                               width = .1) +
  geom_boxplot(colour="grey10", fill = "#185994", alpha=.2, outlier.shape = NA) +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.title.y = element_text(colour="grey20",size=12),
        axis.title.x = element_text(colour="grey20",size=12),
        axis.text.x = element_text(colour="grey20",size=11),
        axis.text.y = element_text(colour="grey20",size=11)) +
  ylab(ylab_title) +
  xlab(xlab_title) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") +
  scale_x_discrete(name="feedback condition", limits = rev(levels(erp$Condition)))
#erp_min_peak

ylab_title <- 'negative peak latency (ms)'
erp_min_lat <- ggplot(erp, aes(x = Condition,
                                y = Min_Peak_Latency, 
                                fill = Condition,
                                colour = Condition,
                                shape = Condition)) +
  #geom_line(aes(group=Participant), colour="grey10", size=.5, alpha=.2) +
  ggbeeswarm::geom_quasirandom(alpha = .4,
                               shape = 22,
                               colour = "grey10",
                               fill = "#185994",
                               size = 2.3,
                               width = .1) +
  geom_boxplot(colour="grey10", fill = "#185994", alpha=.2, outlier.shape = NA) +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.title.y = element_text(colour="grey20",size=12),
        axis.title.x = element_text(colour="grey20",size=12),
        axis.text.x = element_text(colour="grey20",size=11),
        axis.text.y = element_text(colour="grey20",size=11)) +
  ylab(ylab_title) +
  xlab(xlab_title) +
  #stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") +
  scale_x_discrete(name="feedback condition", limits = rev(levels(erp$Condition)))
#erp_min_lat

require(cowplot)
two_plot <- plot_grid(erp_min_peak, erp_min_lat, labels=c("AUTO"), nrow=1, ncol = 2, align = "hv")
two_plot

