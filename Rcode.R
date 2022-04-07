######read nhance demo and oral examination
library(readr)
library(dplyr)
library(haven)
setwd("/")
#read 2013 data 
demo2013<-read_xpt('DEMO_H_2014.XPT')
peri2013<-read_xpt('OHXPER_H_2014.XPT')
#read 2011 data
demo2011<-read_xpt('DEMO_G_2012.XPT')
peri2011<-read_xpt('OHXPER_G_2012.XPT')
age18_2011<-subset(demo2011, demo2011$RIDAGEYR>=18)
#read 2009 data
demo2009<-read_xpt('DEMO_F_2010.XPT')
peri2009<-read_xpt('OHXPER_F_2010.XPT')
#merge 6 year data
perio_demo_6<-rbind(perio_complete_2013_2,perio_complete_2011_2,perio_complete_2009_2)
# read clustering results
# read 2012 definition and 2018 classification 
library(pROC)
library(readxl)
cluster3<-read_xlsx("3clsuters.xlsx", ) 
#calculate multiclass AUC
auc124 <- multiclass.roc(data$Def_2012_4group,data$CLUSTER)
print(auc124)  
# plot confusion matrix
library(ggplot2)
conf_matrix <- function(df.true, df.pred, title = "", true.lab ="", pred.lab ="",
                        high.col = 'red', low.col = 'white') {
  df.true <- as.factor(df.true)
  df.pred <- factor(df.pred, levels = levels(df.true))
  df.cm <- table(True = df.true, Pred = df.pred)
  df.cm.col <- df.cm / rowSums(df.cm)
  df.table <- reshape2::melt(df.cm)
  df.table.col <- reshape2::melt(df.cm.col)
  df.table <- left_join(df.table, df.table.col, by =c("True", "Pred"))
  acc.vector <- c(diag(df.cm)) / c(rowSums(df.cm))
  class.acc <- data.frame(Pred = "Class Acc.", True = names(acc.vector), value = acc.vector)
  acc <- sum(diag(df.cm)) / sum(df.cm)
  ggplot() +
    geom_tile(aes(x=Pred, y=True, fill=value.y),
              data=df.table, size=0.2, color=grey(0.5)) +
    geom_tile(aes(x=Pred, y=True),
              data=df.table[df.table$True==df.table$Pred, ], size=1, color="black", fill = 'transparent') +
    scale_x_discrete(position = "top",  limits = c(levels(df.table$Pred), "")) +
    scale_y_discrete(limits = rev(unique(levels(df.table$Pred)))) +
    labs(x=pred.lab, y=true.lab, fill=NULL
    ) +
    geom_text(aes(x=Pred, y=True, label=value.x),
              data=df.table, size=4, colour="black") +
    geom_text(data = class.acc, aes(Pred, True, label = paste0(round(100*value), "%"))) +
    scale_fill_gradient(low=low.col, high=high.col, labels = scales::percent,
                        limits = c(0,1), breaks = c(0,0.5,1)) +
    guides(size=F) +
    theme_bw() +
    theme(panel.border = element_blank(), legend.position = "bottom",
          axis.text = element_text(color='black'), axis.ticks = element_blank(),
          panel.grid = element_blank(), axis.text.x.top = element_text(angle = 30, vjust = 0, hjust = 0)) +
    coord_fixed()} 
conf_matrix(data$Def_2012_4group, data$Clas_2018_3group, title='')

# sankey plot between NEWDEF and CLUSTER
library(networkD3)
library(dplyr)
links0<-data.frame(
  source=c(paste0(dataCLUSTER3,'_3clusters'),
           paste0(data$Clas_2018_3group,'_2018'),
           paste0(data$Def_2012_4group,'_2012')),
  target=c( paste0(data$Clas_2018_3group,'_2018'),
            paste0(data$Def_2012_4group,'_2012'),
            paste0(data$CLUSTER,'_4clusters')),
  value=c(1))
links0$source<-as.character(links0$source)
links0$target<-as.character(links0$target)
nodes0<- data.frame(
  name = c(as.character(links0$source), 
           as.character(links0$target)) %>% 
    unique())
links0$IDsource <- match(links0$source, nodes0$name)-1
links0$IDtarget <- match(links0$target, nodes0$name)-1
my_color0<-'d3.scaleOrdinal().domain([ ]).
                                      range([ ])'
sankey0<-sankeyNetwork(Links = links0, Nodes = nodes0,
                       Source = "IDsource", 
                       Target = "IDtarget",
                       Value = "value", 
                       NodeID = "name", 
                       fontSize =16,
                       colourScale = my_color0,
                       LinkGroup = "group",
                       units = 'votes',          
                       nodeWidth = 15, 
                       iterations = 0,
                       sinksRight=FALSE)
print(sankey0)
##read medical data 
#read 2009 medical data
med2009<-read_xpt('MCQ_F_2010.XPT')
dia2009<-read_xpt("DIQ_F_2010.XPT")
med2009_2<-select(med2009,c("SEQN","MCQ010","MCQ070","MCQ080","MCQ082",
                            "MCQ160A","MCQ160B","MCQ160C","MCQ160D","MCQ160E",
                            "MCQ160F","MCQ160G","MCQ160K","MCQ160L","MCQ160M",
                            "MCQ160N",
                            "MCQ220"))
#read 2012 medical data
med2011<-read_xpt('MCQ_G_2012.XPT')
dia2011<-read_xpt("DIQ_G_2012.XPT")
colnames(med2011)
med2011_2<-select(med2011,c("SEQN","MCQ010","MCQ070","MCQ080","MCQ082",
                            "MCQ160A","MCQ160B","MCQ160C","MCQ160D","MCQ160E",
                            "MCQ160F","MCQ160G","MCQ160K","MCQ160L","MCQ160M",
                            "MCQ160N",
                            "MCQ220"))
# read 2014 data 
med2013<-read_xpt('MCQ_H_2014.XPT')
dia2013<-read_xpt("DIQ_H_2014.XPT")
colnames(med2013)
med2013_2<-select(med2013,c("SEQN","MCQ010","MCQ070","MCQ080","MCQ082",
                            "MCQ160A","MCQ160B","MCQ160C","MCQ160D","MCQ160E",
                            "MCQ160F","MCQ160G","MCQ160K","MCQ160L","MCQ160M",
                            "MCQ160N",
                            "MCQ220"))
# merge medical data 6 years 
med_data_6<-rbind(dia_med2009,dia_med2011,dia_med2013)
#merge med data and perio-demo data for regression analysis
data_reg<-merge(perio_data_regression,med_data_regression,by="SEQN",all = F)
# regression
glm1<-glm(diabetes ~ RIAGENDR + RIDAGEYR + INDFMIN2 + RIDRETH1 + DMDEDUC2 + Def_2012_4group, 
          family = binomial(),
          data=data_reg_r)
summary(glm1)
