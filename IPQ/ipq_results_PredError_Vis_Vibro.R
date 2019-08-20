
# read in csv data
ipq <- read.csv(file="P:\\Project_Sezen\\data_processing\\IPQ\\ipq_long.csv", header=TRUE, sep=";")
nr_of_subjects <- 19

# add factor 
questions <- rep(c("G1", "SP4", "INV1", "REAL2"), nr_of_subjects)
ipq$questions <- questions

require(tidyr)
library(tidyr)
#ipq_long <- gather(ipq, condition, ipq_score, ems:visual, factor_key=TRUE)
ipq_long <- gather(ipq, condition, ipq_score, visual:vibro, factor_key=TRUE)

# summary and stats
aggregate(ipq_score ~ condition + questions, ipq_long, mean)
# condition questions ipq_score
# 1    visual        G1  4.578947
# 2     vibro        G1  4.684211
# 3    visual      INV1  4.736842
# 4     vibro      INV1  4.736842
# 5    visual     REAL2  3.684211
# 6     vibro     REAL2  4.157895
# 7    visual       SP4  4.894737
# 8     vibro       SP4  5.263158

# ANOVA
require(lsr)
library(lsr)
fit_ipq_long <- aov(ipq_score ~ condition, data=ipq_long)
summary(fit_ipq_long)
# Df Sum Sq Mean Sq F value Pr(>F)
# condition     1   2.13   2.132   1.375  0.243
# Residuals   150 232.58   1.550  
etaSquared(fit_ipq_long, type = 2, anova = TRUE)
#             eta.sq    eta.sq.part         SS  df       MS        F        p
# condition 0.009081736 0.009081736   2.131579   1 2.131579 1.374745 0.242857
# Residuals 0.990918264          NA 232.578947 150 1.550526       NA       NA

# significance testing
require(ggpubr)
library(ggpubr)
my_comparisons <- list(c("visual", "vibro"))

### plot ALL questions
# Boxplot Duration
# specify title, labs titles, legend title
ylab_title <- 'ipq score'
xlab_title <- 'feedback condition'

ipq_long$condition <- factor(ipq_long$condition, levels = rev(levels(ipq_long$condition)))
ipq_long$subject <- as.factor(ipq_long$subject)

require(ggplot2)
#ipq_all <- ggplot(ipq_long, aes(x = condition, y = ipq_score, fill = condition, colour = condition, shape = condition)) +
#  geom_line(aes(group=subject), colour="grey10", size=.5, alpha=.2) +
#  ggbeeswarm::geom_quasirandom(alpha = .4,
#                               shape = 22,
#                               colour = "grey10",
#                               fill = "#185994",
#                               size = 2.3,
#                               width = .1) +
#  geom_boxplot(colour="grey10", fill = "#185994", alpha=.2, outlier.shape = NA) +
#  facet_grid(.~questions) + 
#  theme_bw() +
#  theme(legend.position = "none") +
#  theme(axis.title.y = element_text(colour="grey20",size=12),
#        axis.title.x = element_text(colour="grey20",size=12),
#        axis.text.x = element_text(colour="grey20",size=11),
#        axis.text.y = element_text(colour="grey20",size=11)) +
#  ylab(ylab_title) +
#  xlab(xlab_title) +
#  stat_compare_means(comparisons = my_comparisons) +
#  ylim(1, 7)
#ipq_all

#barplot
ipq_all <- ggplot(ipq_long, aes(x = condition, y = ipq_score, fill = condition, colour = condition, shape = condition)) +
  geom_line(aes(group=subject), colour="grey10", size=.5, alpha=.2) +
  ggbeeswarm::geom_quasirandom(alpha = .4,
                               shape = 22,
                               colour = "grey10",
                               fill = "#185994",
                               size = 2.3,
                               width = .1) +
  geom_col(colour="grey10", fill = "#185994", alpha=.2) +
  #geom_bar(colour="grey10", fill = "#185994", alpha=.2, outlier.shape = NA) +
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
ylab_title <- 'ipq score (G1)'
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
  stat_compare_means(comparisons = my_comparisons) +
  ylim(1, 7)
ipq_1

# select only 1 question and plot it
this_question <- c("REAL2")
ylab_title <- 'ipq score (REAL2)'
ipq_single_question <- ipq_long[ipq_long$questions %in% this_question,]
ipq_2 <- ggplot(ipq_single_question, aes(x = condition,
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
ipq_2

# inspect
aggregate(ipq_score ~ condition, ipq_single_question, mean)
# condition ipq_score
# 1    visual  4.545455
# 2     vibro  4.818182
# 3       ems  4.636364

# ANOVA
fit_ipq_single_question <- aov(ipq_score ~ condition, data=ipq_single_question)
summary(fit_ipq_single_question)
etaSquared(fit_ipq_single_question, type = 2, anova = TRUE)
# eta.sq eta.sq.part         SS df        MS         F         p
# condition 0.01674641  0.01674641  0.4242424  2 0.2121212 0.2554745 0.7762169
# Residuals 0.98325359          NA 24.9090909 30 0.8303030        NA        NA
