
nr_of_subjects <- 11

# read in csv data ERP
erp <- read.csv(file="P:\\Project_Sezen\\data_processing\\ERP\\peaks_locs_ERP_incl_EMS_FCz.csv", header=TRUE, sep=";")
# make factors
erp$Participant <- as.factor(erp$Participant)
erp$Condition <- as.factor(erp$Condition)
erp <- erp[order(erp$Condition),]
erp

# read in csv data
ipq <- read.csv(file="P:\\Project_Sezen\\data_processing\\IPQ\\ipq_long.csv", header=TRUE, sep=";")

# add factor 
questions <- rep(c("G1", "SP4", "INV1", "REAL2"), nr_of_subjects)
ipq$questions <- questions

library(tidyr)
ipq_long <- gather(ipq, condition, ipq_score, ems:visual, factor_key=TRUE)

# summary and stats
# remove/revers inv1 item?
select <- c("G1")
ipq_sel <- ipq_long[ipq_long$question %in% select,]
ipq_dat <- aggregate(ipq_score ~ subject + condition, ipq_long, mean)
#ipq_dat <- ipq_dat[order(ipq_dat$condition),]
ipq_dat

# correlation
library(ggpubr)
cor <- cor.test(erp$Min_Peak_Amplitude, ipq_sel$ipq_score, method="pearson")
cor
# mean of all IPQ items
# Pearson's product-moment correlation
# 
# data:  erp$Min_Peak_Amplitude and ipq_dat$ipq_score
# t = -0.42544, df = 31, p-value = 0.6735
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
# -0.4088047  0.2742950
# sample estimates:
# cor 
# -0.07618848 

### only IPQ item G1
# Pearson's product-moment correlation
# 
# data:  erp$Min_Peak_Amplitude and ipq_sel$ipq_score
# t = 0.48169, df = 31, p-value = 0.6334
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.2649573  0.4171570
# sample estimates:
#        cor 
# 0.08619192 

df <- cbind(erp[,c(1:3)], ipq_dat$ipq_score)
df

# make paired observations plot
ylab_title <- 'ipq score (mean all items)'
xlab_title <- 'min. peak amplitude (µV)'

names(df) <- c('participant', 'group', 'condition1', 'condition2')

library(plyr)
df$group <- revalue(df$group, c("ems"="EMS+Vibro+Visual", "vibro"="Vibro+Visual", "visual"="Visual"))

# abline coeffs
coef(lm(df$condition2 ~ df$condition1))
# (Intercept) df$condition1 
# 4.4562118    -0.0248958 

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
  ylim(c(1,7))+
  xlim(c(-11.5,0))+
  geom_abline(slope = -0.0248958 , intercept = 4.4562118)
paired_obs


