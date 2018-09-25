
# read in csv data
rt <- read.csv(file="P:\\Project_Sezen\\data_processing\\RT\\RT_results_ALL.csv", header=TRUE, sep=";")
nr_of_subjects <- 11

# make factors
rt$RT <- as.numeric(rt$RT)
rt$subject <- as.factor(rt$subject)
rt$condition <- as.factor(rt$condition)
rt$congruency <- as.factor(rt$congruency)

# remove outliers
outlier_thresh <- min(boxplot.stats(rt$RT)$out)
rt_no_outliers <- subset(rt, rt[ , 4] < outlier_thresh)

# inspect
aggregate(RT ~ condition + congruency, rt_no_outliers, mean)
# condition congruency         RT
# 1       ems   conflict 0.22473725
# 2     vibro   conflict 0.24741854
# 3    visual   conflict 0.26353916
# 4       ems     normal 0.04646336
# 5     vibro     normal 0.06927854
# 6    visual     normal 0.07673378

# aggreate data by condition
rt_no_outliers_normal <- rt_no_outliers[rt_no_outliers$congruency %in% c("normal"),]
rt_no_outliers_conflict <- rt_no_outliers[rt_no_outliers$congruency %in% c("conflict"),]
rt_agg_norm <- aggregate(RT ~ subject + condition, rt_no_outliers_normal, mean)
rt_agg_conf <- aggregate(RT ~ subject + condition, rt_no_outliers_conflict, mean)

# normal
f1 <- function(x) c(Mean = mean(x), Max = max(x), SD = sd(x))
do.call(data.frame, aggregate(RT ~ condition, rt_no_outliers_normal, f1))
# condition    RT.Mean RT.Max      RT.SD
# 1       ems 0.04646336  0.412 0.05389574
# 2     vibro 0.06927854  0.414 0.06331435
# 3    visual 0.07673378  0.394 0.06283353

# conflict
do.call(data.frame, aggregate(RT ~ condition, rt_no_outliers_conflict, f1))
# condition   RT.Mean RT.Max      RT.SD
# 1       ems 0.2247373  0.430 0.10689765
# 2     vibro 0.2474185  0.428 0.07814822
# 3    visual 0.2635392  0.430 0.06921642

# stats aggregate data
library(lsr)
rt_no_outliers_conflict_agg <- aggregate(RT ~ subject + condition, rt_no_outliers_conflict, mean)
fit_conf_no_outliers <- aov(RT ~ condition, data=rt_no_outliers_conflict_agg)#[-c(1,12),])
summary(fit_conf_no_outliers)
etaSquared(fit_conf_no_outliers, type = 2, anova = TRUE)
# eta.sq eta.sq.part         SS df          MS        F         p
# condition 0.04754661  0.04754661 0.00828848  2 0.004144240 0.723842 0.4934398
# Residuals 0.95245339          NA 0.16603479 29 0.005725338       NA        NA

rt_with_outliers_normal <- rt[rt$congruency %in% c("normal"),]
# normal
aggregate(RT ~ condition, rt_with_outliers_normal, mean)
# condition         RT
# 1       ems 0.06095139
# 2     vibro 0.08926677
# 3    visual 0.08277067

rt_with_outliers_conflict <- rt[rt$congruency %in% c("conflict"),]
aggregate(RT ~ condition, rt_with_outliers_conflict, mean)
# condition        RT
# 1       ems 0.2471927
# 2     vibro 0.2942874
# 3    visual 0.2920187

### plot
# Boxplot Duration
# specify title, labs titles, legend title
ylab_title <- 'reaction time in s. (mismatch trials)'
xlab_title <- 'feedback condition'
# post hoc significance testing

require(ggpubr)
my_comparisons <- list(c("visual", "vibro"), c("visual", "ems"), c("vibro", "ems"))

rt_plot_conf <- ggplot(rt_no_outliers_conflict_agg, aes(x = condition,
                      y = RT, 
                      fill = condition,
                      colour = condition,
                      shape = condition)) +
  geom_line(aes(group=subject), colour="grey10", size=.5, alpha=.2) +
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
  #stat_compare_means(comparisons = my_comparisons) +
  ylim(0,0.4) + 
  scale_x_discrete(name="feedback condition", limits = rev(levels(rt_no_outliers_conflict_agg$condition)))
rt_plot_conf

###
rt_plot_norm <- ggplot(rt_agg_norm[-c(1,12),], aes(x = condition,
                                        y = RT, 
                                        fill = condition,
                                        colour = condition,
                                        shape = condition)) +
  geom_line(aes(group=subject), colour="grey10", size=.5, alpha=.2) +
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
  #stat_compare_means(comparisons = my_comparisons) +
  ylim(0,0.4) +
  scale_x_discrete(name="feedback condition", limits = rev(levels(erp$Condition)))
rt_plot_norm

require(cowplot)
all_plot <- plot_grid(rt_plot_norm, rt_plot_conf, labels=c("A", "B"), nrow=1, ncol = 2, align = "hv")
all_plot
