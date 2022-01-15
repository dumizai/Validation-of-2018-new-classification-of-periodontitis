library(dplyr)
library(matrixStats)
library(readr)
library(dplyr)
library(haven)
getwd()
setwd("")
#read 2013 data 
demo2013<-read_xpt('DEMO_H_2014.XPT')
peri2013<-read_xpt('OHXPER_H_2014.XPT')
age18_2013<-subset(demo2013, demo2013$RIDAGEYR>=18)
# keep age， gender，  "SEQN"
age18_2013_2<-age18_2013[,c("RIAGENDR","RIDAGEYR","SEQN")]
age18_2013_3<-age18_2013[,c("DMDEDUC2","INDFMIN2","RIDRETH1","RIAGENDR","RIDAGEYR","SEQN")]  # DMDEDUC2 (education level aged 20+)
demo_peri_2013<-merge(age18_2013_3, peri2013, by.x = 'SEQN', by.y = 'SEQN', all=F)
#those who complete perio examination
# OHASCST4: 1 Complete,  2 Partial, 3 Not done, .Missing
perio_complete_2013<-subset(demo_peri_2013,demo_peri_2013$OHDPDSTS=="1")
# remove useless variables # "OHAEXSTS"  ”OHASCST4“
variable.names(perio_complete_2013)
perio_complete_2013_2<-perio_complete_2013[,-c(7:9)]

#read 2011 data
demo2011<-read_xpt('DEMO_G_2012.XPT')
peri2011<-read_xpt('OHXPER_G_2012.XPT')
age18_2011<-subset(demo2011, demo2011$RIDAGEYR>=18)
# keep age， gender，  "SEQN"
age18_2011_2<-age18_2011[,c("RIAGENDR","RIDAGEYR","SEQN")]
age18_2011_3<-age18_2011[,c("DMDEDUC2","INDFMIN2","RIDRETH1","RIAGENDR","RIDAGEYR","SEQN")]
demo_peri_2011<-merge(age18_2011_3, peri2011, by.x = 'SEQN', by.y = 'SEQN', all=F)
#those who complete perio examination
# OHASCST4: 1 Complete,  2 Partial, 3 Not done, .Missing
perio_complete_2011<-subset(demo_peri_2011,demo_peri_2011$OHDPDSTS=="1")
# remove useless variables # "OHAEXSTS"  ”OHASCST4“
variable.names(perio_complete_2011)
perio_complete_2011_2<-perio_complete_2011[,-c(7:9)]

#read 2009 data
demo2009<-read_xpt('DEMO_F_2010.XPT')
peri2009<-read_xpt('OHXPER_F_2010.XPT')
age18_2009<-subset(demo2009, demo2009$RIDAGEYR>=18)
# keep age， gender，  "SEQN"
age18_2009_2<-age18_2009[,c("RIAGENDR","RIDAGEYR","SEQN")]
age18_2009_3<-age18_2009[,c("DMDEDUC2","INDFMIN2","RIDRETH1","RIAGENDR","RIDAGEYR","SEQN")]
demo_peri_2009<-merge(age18_2009_3, peri2009, by.x = 'SEQN', by.y = 'SEQN', all=F)
#those who complete perio examination
# OHASCST4: 1 Complete,  2 Partial, 3 Not done, .Missing
perio_complete_2009<-subset(demo_peri_2009,demo_peri_2009$OHDPDSTS=="1")
# remove useless variables # "OHAEXSTS"  ”OHASCST4“
variable.names(perio_complete_2009)
perio_complete_2009_2<-perio_complete_2009[,-c(7:9)]

#merge 6 year data
perio_demo_6<-rbind(perio_complete_2013_2,perio_complete_2011_2,perio_complete_2009_2)

# using matlab file to compute the mean.
#read 6teeth data
mo<-read_csv("")
mo_contains <- mo[ , grepl('SEQN|LA', names(mo))]
mo_containsPD <- mo[ , grepl('SEQN|PC', names(mo))]
mo_contains01<-mo[ , grepl('SEQN|LAM', names(mo))]
mo_contains02<-mo[ , grepl('SEQN|LAL', names(mo))]
mo_contains0<-merge(mo_contains01,mo_contains02,by="SEQN",all=F)
mo_contains$RS0<-rowSums(mo_contains0[,c(2:57)]>=3)
mo_contains1<-mo[ , grepl('SEQN|LAD', names(mo))]
mo_contains2<-mo[ , grepl('SEQN|LAS', names(mo))]
mo_contains3<-mo[ , grepl('SEQN|LAP', names(mo))]
mo_contains4<-mo[ , grepl('SEQN|LAA', names(mo))]
mo_contains00<-cbind(mo_contains1,mo_contains2[,c(2:29)],
                     mo_contains3[,c(2:29)],mo_contains4[,c(2:29)])
mo_contains$maxcal= rowMaxs(as.matrix(mo_contains00[,c(2:113)]))
mo_contains$maxpd= rowMaxs(as.matrix(mo_containsPD[,c(2:169)]))
mo_contains$newdef = ifelse((mo_contains$RS0 <=2 | mo_contains$maxcal <1),
                            "Periodontal_healthy", "StageI") # perio_demo_6_contains$maxpd <=3 
table(mo_contains$newdef)
d0<-subset(mo_contains,mo_contains$newdef=="Periodontal_healthy")
d1<-subset(mo_contains,mo_contains$newdef=="StageI")
d1$newdef= ifelse((d1$maxcal <=2 | d1$maxpd <4),"StageI", "StageII") # 
table(d1$newdef)
d2<-subset(d1,d1$newdef=="StageII")
d2$newdef= ifelse((d2$maxcal <=4 | d2$maxpd <6),"StageII", "StageIII-IV") #
table(d2$newdef)
perio_demo_60<-rbind(d0,d1[which(d1$newdef=="StageI"),],d2)
mo1<-merge(perio_demo_6[,c(1:6)],perio_demo_60[,c(1,173)],by="SEQN",all = F)
table(mo1$newdef)

#read cluster
library(readxl)
clusters<-read_xlsx("")
table(clusters$CLUSTER)
#merge cluster number and perio- and newdef
perio_demo_6_cluster<-merge(mo1,clusters,by='SEQN',all = F)
perio_data_regression<-perio_demo_6_cluster[,c('SEQN',"RIAGENDR","RIDAGEYR","INDFMIN2","RIDRETH1","DMDEDUC2",
                                               "newdef", "CLUSTER" )]
colnames(perio_data_regression)[which(names(perio_data_regression) == "newdef")] <- "NEWDEF"

# multiclass ROC curve
library(pROC)
dataROC<-perio_data_regression[,c(7,8)]
levels(dataROC$NEWDEF)<-c(1,2,3,4)
dataROC$NEWDEF[which(dataROC$NEWDEF=="Periodontal_healthy")]<-1
dataROC$NEWDEF[which(dataROC$NEWDEF=="StageI")]<-2
dataROC$NEWDEF[which(dataROC$NEWDEF=="StageII")]<-3
dataROC$NEWDEF[which(dataROC$NEWDEF=="StageIII-IV")]<-4
dataROC$NEWDEF<-as.ordered(dataROC$NEWDEF)
table(dataROC$NEWDEF)
levels(dataROC$CLUSTER)<-c(1,2,3,4)
dataROC$CLUSTER[which(dataROC$CLUSTER==1)]<-1
dataROC$CLUSTER[which(dataROC$CLUSTER==2)]<-2
dataROC$CLUSTER[which(dataROC$CLUSTER==3)]<-3
dataROC$CLUSTER[which(dataROC$CLUSTER==4)]<-4
dataROC$CLUSTER<-as.ordered(dataROC$CLUSTER)
table(dataROC$CLUSTER)
auc <- multiclass.roc(dataROC$NEWDEF,dataROC$CLUSTER)
print(auc) 

library(ggplot2)
conf_matrix <- function(df.true, df.pred, title = "", true.lab ="2018 case definition", pred.lab ="Periodontitis clusters",
                        high.col = 'red', low.col = 'white') {
  #convert input vector to factors, and ensure they have the same levels
  df.true <- as.factor(df.true)
  df.pred <- factor(df.pred, levels = levels(df.true))
  #generate confusion matrix, and confusion matrix as a pecentage of each true class (to be used for color) 
  df.cm <- table(True = df.true, Pred = df.pred)
  df.cm.col <- df.cm / rowSums(df.cm)
  #convert confusion matrices to tables, and binding them together
  df.table <- reshape2::melt(df.cm)
  df.table.col <- reshape2::melt(df.cm.col)
  df.table <- left_join(df.table, df.table.col, by =c("True", "Pred"))
  #calculate accuracy and class accuracy
  acc.vector <- c(diag(df.cm)) / c(rowSums(df.cm))
  class.acc <- data.frame(Pred = "Class Acc.", True = names(acc.vector), value = acc.vector)
  acc <- sum(diag(df.cm)) / sum(df.cm)
  #plot
  ggplot() +
    geom_tile(aes(x=Pred, y=True, fill=value.y),
              data=df.table, size=0.2, color=grey(0.5)) +
    geom_tile(aes(x=Pred, y=True),
              data=df.table[df.table$True==df.table$Pred, ], size=1, color="black", fill = 'transparent') +
    scale_x_discrete(position = "top",  limits = c(levels(df.table$Pred), "")) +#Class Acc.
    scale_y_discrete(limits = rev(unique(levels(df.table$Pred)))) +
    labs(x=pred.lab, y=true.lab, fill=NULL) +
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
conf_matrix(dataROC$NEWDEF, dataROC$CLUSTER, title='')
table(dataROC$NEWDEF,dataROC$CLUSTER)

# sankey plot between NEWDEF and CLUSTER
library(networkD3)
links0<-data.frame(
  source=perio_data_regression$NEWDEF,
  target=perio_data_regression$CLUSTER,
  value=c(1))
nodes0<- data.frame(
  name = c(as.character(links0$source), 
           as.character(links0$target)) %>% 
    unique())
links0$IDsource <- match(links0$source, nodes0$name)-1 
links0$IDtarget <- match(links0$target, nodes0$name)-1
links0$group[links0$source=="Periodontal_healthy"]<-"Periodontal_healthy"
links0$group[links0$source=="StageI"]<-"StageI"
links0$group[links0$source=="StageII"]<-"StageII"
links0$group[links0$source=="StageIII-IV"]<-"StageIII-IV"
my_color0<-'d3.scaleOrdinal().domain(["Periodontal_healthy","StageI","StageII","StageIII-IV",
                                      "1","2","3","4"]).
                                      range([ "#FDDBC7" ,"#F4A582","#D6604D","#B2182B", 
                                            "#D1E5F0", "#92C5DE" ,"#4393C3","#2166AC" ])'
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
                       sinksRight=FALSE)
print(sankey0)

# read medical data
#read 2009 medical data
med2009<-read_xpt('MCQ_F_2010.XPT')
dia2009<-read_xpt("DIQ_F_2010.XPT")
med2009_2<-select(med2009,c("SEQN","MCQ010","MCQ070","MCQ080","MCQ082",
                            "MCQ160A","MCQ160B","MCQ160C","MCQ160D","MCQ160E",
                            "MCQ160F","MCQ160G","MCQ160K","MCQ160L","MCQ160M",
                            "MCQ160N","MCQ220"))
dia2009_2<-select(dia2009,c("SEQN","DIQ010"))
dia_med2009<-merge(dia2009_2, med2009_2, by.x = 'SEQN', by.y = 'SEQN', all=F)

# example codes for regression analysis
glm1<-glm(diabetes ~ RIAGENDR + RIDAGEYR + INDFMIN2 + RIDRETH1 + DMDEDUC2 + NEWDEF, 
          family = binomial(),
          data=data)
summary(glm1)
#odds ratios and 95% CI
diab_newdef<- data.frame(exp(cbind(OR = coef(glm1), confint(glm1))))
round(diab_newdef[c(28:30),],2)