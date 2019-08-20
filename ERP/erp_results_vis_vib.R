
# read in csv data
#erp <- read.csv(file="P:\\Project_Sezen\\data_processing\\ERP\\peaks_locs_ERP_incl_EMS_FCz_95.csv", header=TRUE, sep=";")
#erp <- read.csv(file="P:\\Project_Sezen\\data_processing\\ERP\\peak_locs_baseline_250\\peaks_locs_ERP_FCz_100.csv", header=TRUE, sep=";")
#erp <- read.csv(file="P:\\Project_Sezen\\data_processing\\ERP\\peak_locs_baseline_250\\peaks_locs_ERP_Pz_100.csv", header=TRUE, sep=";")
#erp <- read.csv(file="P:\\Project_Sezen\\data_processing\\ERP\\peak_locs_baseline_250\\peaks_locs_ERP_Oz_100.csv", header=TRUE, sep=";")
erp <- read.csv(file="P:\\Project_Sezen\\data_processing\\ERP\\peak_locs_baseline_100\\peaks_locs_ERP_Fz_100.csv", header=TRUE, sep=";")
#erp <- read.csv(file="P:\\Project_Sezen\\data_processing\\ERP\\peak_locs_baseline_250\\peaks_locs_ERP_Cz_100.csv", header=TRUE, sep=";")
#erp <- read.csv(file="P:\\Project_Sezen\\data_processing\\ERP\\peaks_locs_ERP_incl_EMS_FCz_105.csv", header=TRUE, sep=";")

nr_of_subjects <- 19

# make factors
erp$Participant <- as.factor(erp$Participant)
erp$Condition <- as.factor(erp$Condition)

# inspect
f1 <- function(x) c(Mean = mean(x), Max = max(x), SD = sd(x))

do.call(data.frame, aggregate(Min_Peak_Amplitude~Condition, erp, f1))
# Condition Min_Peak_Amplitude.Mean Min_Peak_Amplitude.Max Min_Peak_Amplitude.SD
# 1     vibro               -3.372541                1.76100              4.062581
# 2    visual               -2.809810                0.18644              2.516000

do.call(data.frame, aggregate(Min_Peak_Latency~Condition, erp, f1))
# Condition Min_Peak_Latency.Mean Min_Peak_Latency.Max Min_Peak_Latency.SD
# 1     vibro              193.0526                  300            49.43848
# 2    visual              198.9474                  300            63.27144

# ANOVA on aggregated data Min Peak Amplitude
require(lsr)
fit_amp <- aov(Min_Peak_Amplitude ~ Condition, data=erp)
#summary(fit_amp)
etaSquared(fit_amp, type = 2, anova = TRUE)
# eta.sq eta.sq.part        SS df        MS        F          p
# Condition 0.007265869 0.007265869   3.008325  1  3.008325 0.2634858 0.6108714
# Residuals 0.992734131          NA 411.026822 36 11.417412        NA        NA

# ANOVA on aggregated data Min Peak Latency
fit_lat <- aov(Min_Peak_Latency ~ Condition, data=erp)
#summary(fit_lat)
etaSquared(fit_lat, type = 2, anova = TRUE)
# eta.sq eta.sq.part         SS df       MS         F         p
# Condition 0.002836346 0.002836346    330.1053  1  330.1053 0.1023989 0.7508195
# Residuals 0.997163654          NA 116053.8947 36 3223.7193        NA        NA

### plot
# Boxplot Duration
# specify title, labs titles, legend title
xlab_title <- 'feedback condition'

# post hoc significance testing
require(ggpubr)
my_comparisons <- list(c("visual", "vibro"))

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

