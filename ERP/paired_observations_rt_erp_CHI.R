
nr_of_subjects <- 11

# read in csv data ERP
erp <- read.csv(file="P:\\Project_Sezen\\data_processing\\ERP\\peaks_locs_ERP_incl_EMS_FCz.csv", header=TRUE, sep=";")
# make factors
erp$Participant <- as.factor(erp$Participant)
erp$Condition <- as.factor(erp$Condition)
erp <- erp[order(erp$Condition),]

# read in csv data RT
rt <- read.csv(file="P:\\Project_Sezen\\data_processing\\RT\\RT_results_ALL.csv", header=TRUE, sep=";")
# make factors
rt$RT <- as.numeric(rt$RT)
rt$subject <- as.factor(rt$subject)
rt$condition <- as.factor(rt$condition)
rt$congruency <- as.factor(rt$congruency)
# remove outliers
outlier_thresh <- min(boxplot.stats(rt$RT)$out)
rt_no_outliers <- subset(rt, rt[ , 4] < outlier_thresh)
rt_no_outliers_conflict <- rt_no_outliers[rt_no_outliers$congruency %in% c("conflict"),]
rt_agg_conf <- aggregate(RT ~ subject + condition, rt_no_outliers_conflict, mean)
rt_agg_conf <- rt_agg_conf[order(rt_agg_conf$condition),]

# correlation
library(ggpubr)
cor <- cor.test(erp$Min_Peak_Amplitude[-23], rt_agg_conf$RT, method="spearman")
# Spearman's rank correlation rho
# 
# data:  erp$Min_Peak_Amplitude[-23] and rt_agg_conf$RT
# S = 3622, p-value = 0.06053
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.3361437

df <- cbind(erp[-c(23),c(1:3)], rt_agg_conf$RT)

# make paired observations plot
ylab_title <- 'reaction time in s. (mismatch trials)'
xlab_title <- 'min. peak amplitude (µV)'

names(df) <- c('participant', 'group', 'condition1', 'condition2')

library(plyr)
df$group <- revalue(df$group, c("ems"="EMS+Vibro+Visual", "vibro"="Vibro+Visual", "visual"="Visual"))

# abline coeffs
coef(lm(df$condition2 ~ df$condition1))
# (Intercept) df$condition1 
# 0.3134883     0.0105075

# scatterplot of paired observations -----------------
paired_obs <- ggplot(df, aes(x=condition1,y=condition2,group=group,fill=group,colour=group,shape=group)) + 
  geom_point(size=4.5,stroke=.5) +
  theme_bw() +
  scale_shape_manual(values=c(22,21,24)) +
  scale_fill_manual(values = c("#185994", "#D4140A", "#09B567")) +
  scale_colour_manual(values = c("grey5","grey5", "grey5")) +
  theme(axis.text.x = element_text(colour="grey20",size=11),
        axis.text.y = element_text(colour="grey20",size=11),  
        axis.title.x = element_text(colour="grey20",size=14),
        axis.title.y = element_text(colour="grey20",size=14),
        legend.title = element_blank(),
        legend.position = c(.80, .20),
        plot.margin = unit(c(5.5,5.5,5.5,5.5), "pt"),
        legend.text = element_text(colour="grey20",size=12),
        legend.background = element_blank())+ 
  #labs(title="Paired observations") +
  xlab(xlab_title)+
  ylab(ylab_title)+
  ylim(c(0,0.5))+
  xlim(c(-11.5,0))+
  geom_abline(slope = 0.0105075, intercept = 0.3134883)
paired_obs


