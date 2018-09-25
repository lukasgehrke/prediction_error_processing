
# read in csv data
ipq <- read.csv(file="P:\\Project_Sezen\\data_processing\\IPQ\\ipq_long.csv", header=TRUE, sep=";")
nr_of_subjects <- 11

# add factor 
questions <- rep(c("G1", "SP4", "INV1", "REAL2"), nr_of_subjects)
ipq$questions <- questions

library(tidyr)
ipq_long <- gather(ipq, condition, ipq_score, ems:visual, factor_key=TRUE)

# summary and stats
aggregate(ipq_score ~ condition + questions, ipq_long, mean)
# condition questions ipq_score
# 1     visual        G1  4.545455
# 2      vibro        G1  4.818182
# 3        ems        G1  4.636364
# 4     visual      INV1  5.090909
# 5      vibro      INV1  4.545455
# 6        ems      INV1  4.727273
# 7     visual     REAL2  3.545455
# 8      vibro     REAL2  3.818182
# 9        ems     REAL2  3.363636
# 10    visual       SP4  5.363636
# 11     vibro       SP4  5.454545
# 12       ems       SP4  5.181818

# ANOVA
library(lsr)
fit_ipq_long <- aov(ipq_score ~ condition, data=ipq_long)
summary(fit_ipq_long)
etaSquared(fit_ipq_long, type = 2, anova = TRUE)
# eta.sq eta.sq.part          SS  df        MS         F         p
# condition 0.00419426  0.00419426   0.8636364   2 0.4318182 0.2716693 0.7625412
# Residuals 0.99580574          NA 205.0454545 129 1.5894996        NA        NA

# significance testing
library(ggpubr)
my_comparisons <- list(c("visual", "vibro"), c("visual", "ems"), c("vibro", "ems"))

### plot ALL questions
# Boxplot Duration
# specify title, labs titles, legend title
ylab_title <- 'ipq score (G1)'
xlab_title <- 'feedback condition'

ipq_long$condition <- factor(ipq_long$condition, levels = rev(levels(ipq_long$condition)))
ipq_long$subject <- as.factor(ipq_long$subject)

ipq_all <- ggplot(ipq_long, aes(x = condition,
                                     y = ipq_score, fill = condition, colour = condition, shape = condition)) +
  geom_line(aes(group=subject), colour="grey10", size=.5, alpha=.2) +
  ggbeeswarm::geom_quasirandom(alpha = .4,
                               shape = 22,
                               colour = "grey10",
                               fill = "#185994",
                               size = 2.3,
                               width = .1) +
  geom_boxplot(colour="grey10", fill = "#185994", alpha=.2, outlier.shape = NA) +
  facet_grid(.~questions) + 
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.title.y = element_text(colour="grey20",size=12),
        axis.title.x = element_text(colour="grey20",size=12),
        axis.text.x = element_text(colour="grey20",size=11),
        axis.text.y = element_text(colour="grey20",size=11)) +
  ylab(ylab_title) +
  xlab(xlab_title) +
  stat_compare_means(comparisons = my_comparisons) +
  ylim(1, 7)
ipq_all

# select only 1 question and plot it
this_question <- c("G1")
ipq_single_question <- ipq_long[ipq_long$questions %in% this_question,]
ipq_1 <- ggplot(ipq_single_question, aes(x = condition,
                                y = ipq_score, fill = condition, colour = condition, shape = condition)) +
  geom_line(aes(group=subject), colour="grey10", size=.5, alpha=.2) +
  ggbeeswarm::geom_quasirandom(alpha = .4,
                               shape = 22,
                               colour = "grey10",
                               fill = "#185994",
                               size = 2.3,
                               width = .1) +
  geom_boxplot(colour="grey10", fill = "#185994", alpha=.2, outlier.shape = NA) +
  #facet_grid(.~questions) + 
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.title.y = element_text(colour="grey20",size=12),
        axis.title.x = element_text(colour="grey20",size=12),
        axis.text.x = element_text(colour="grey20",size=11),
        axis.text.y = element_text(colour="grey20",size=11)) +
  ylab(ylab_title) +
  xlab(xlab_title) +
  #stat_compare_means(comparisons = my_comparisons) +
  ylim(1, 7)
ipq_1

# inspect
aggregate(ipq_score ~ condition, ipq_single_question, mean)
# condition ipq_score
# 1    visual  4.545455
# 2     vibro  4.818182
# 3       ems  4.636364

# ANOVA
library(lsr)
fit_ipq_single_question <- aov(ipq_score ~ condition, data=ipq_single_question)
summary(fit_ipq_single_question)
etaSquared(fit_ipq_single_question, type = 2, anova = TRUE)
# eta.sq eta.sq.part         SS df        MS         F         p
# condition 0.01674641  0.01674641  0.4242424  2 0.2121212 0.2554745 0.7762169
# Residuals 0.98325359          NA 24.9090909 30 0.8303030        NA        NA
