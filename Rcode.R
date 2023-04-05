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
# GAP  
library(cluster)
library(factoextra)
library(phyloseq)
# Compute GAP statistics for km clustering
gap_statistic = function(data, min_num_clusters = 1, max_num_clusters = 6, num_reference_bootstraps = 10) {
  num_clusters = min_num_clusters:max_num_clusters
  actual_dispersions = maply(num_clusters, function(n) dispersion(data, n))
  ref_dispersions = maply(num_clusters, function(n) reference_dispersion(data, n, num_reference_bootstraps))
  mean_ref_dispersions = ref_dispersions[ , 1]
  stddev_ref_dispersions = ref_dispersions[ , 2]
  gaps = mean_ref_dispersions - actual_dispersions
  print(plot_gap_statistic(gaps, stddev_ref_dispersions, num_clusters))
  print(paste("The estimated number of clusters is ", num_clusters[which.max(gaps)], ".", sep = ""))
  list(gaps = gaps, gap_stddevs = stddev_ref_dispersions)
}
gap_statistic(mo[,-1])
# Compute Silhourre score for k-medoids analysis
pam.res<-pam(mo01,3)
pam.sil<-silhouette(pam.res,dist(mo01))
fviz_silhouette(pam.res) 
# FCM analysis
library(ppclust)
fcm.res<-fcm(mo01,centers=3)
#Compute Sihhouette score for FCM  
si <- function(x, u, v, m, t=NULL, eta, av=1, tidx="f"){
  if(missing(x))
    stop("Missing input argument. A ppclust object or a numeric data set is required")
  tidx <- match.arg(tidx, c("e","f","g"))
  if(inherits(x, "ppclust")){
    X <- as.matrix(x$x)
    if(!is.null(x$u)){
      U <- as.matrix(x$u)
      m <- x$m
    }
    else if(!is.null(x$t)){
      U <- as.matrix(x$t)
      m <- x$eta
    }
    else{
      stop("Argument 'x' does not have the fuzzy membership or typicality matrix")
    }
    V <- as.matrix(x$v)
    if(tidx == "e" || tidx == "g"){
      if(!is.null(x$t)){
        T <- x$t
        eta <- x$eta
      }
      else
        stop("Argument 'x' does not have the typicality matrix")
    }
  }
  else{
    if(!missing(x))
      if(is.matrix(x) || is.data.frame(x) || is.vector(x))
        X <- as.matrix(x)
      else
        stop("Argument 'x' must be a valid instance of the 'ppclust', a numeric vector, data frame or matrix")
      else
        stop("Missing argument 'x'")
      if(!missing(u))
        if(is.matrix(u) || is.data.frame(u))
          U <- as.matrix(u)
        else
          stop("Argument 'u' must be a numeric data frame or matrix")
        else
          stop("Missing argument 'u'")
        if(!missing(v))
          if(is.matrix(v) || is.data.frame(v))
            V <- as.matrix(v)
          else
            stop("Argument 'v' must be a numeric data frame or matrix")
          else
            stop("Missing argument 'v'")
          if(tidx != "f")
            if(!is.null(t))
              if(is.matrix(t) || is.data.frame(t))
                T <- as.matrix(t)
              else
                stop("Argument 't' must be a numeric data frame or matrix")
              else
                stop("Argument 't' is null")
              if(tidx != "f"){
                if(missing(eta))
                  eta <- 2
                if(!is.numeric(eta))
                  stop("Argument 'eta' must be number")
                if(eta < 1)
                  stop("Argument 'eta' should be a positive number equals to or greater than 1")
              }
              if(missing(m))
                m <- 2
              if(!is.numeric(m))
                stop("Argument 'm' must be number")
              if(m < 1)
                stop("Argument 'm' should be a positive number equals to or greater than 1")
  }
  if(nrow(X) != nrow(U))
    stop("The number of rows of data set is not equal to the number of rows of the membership matrix")
  if(ncol(X) != ncol(V))
    stop("The number of columns of the data set matrix is not equal to the number of columns of prototypes matrix")
  if(ncol(U) != nrow(V))
    stop("The number of columns of the membership matrix is not equal to the number of rows of prototypes matrix")
  if(tidx == "g"){
    if(!is.null(T))
      U <- T/rowSums(T)
    else
      stop("Typicality matrix is required to compute the generalized SI index")
  }
  n <- nrow(U)
  k <- ncol(U)
  vm <- vector(length(n), mode = "numeric")
  for(i in 1:n)
    vm[i] <- which.max(U[i,])
  counts <- c()
  for(j in 1:k)
    counts[j] <- length(which(vm == j))
  D <- matrix(nrow = n, ncol = n, 0)
  for(i1 in 1:(n-1))
    for(i2 in (i1+1):n)
      D[i2,i1] <- D[i1,i2] <- sum((X[i1,]-X[i2,])^2)
  a <- b <- si.obj <- rep(0, n)
  Z <- matrix(nrow = n, ncol = k, 0)
  for(i in 1:n)
    for(j in 1:k)
      for(i2 in 1:n)
        if(vm[i2] == j)
          Z[i,j] = Z[i,j] + D[i,i2]
  for(i in 1:n){
    for(j in 1:k){
      if(vm[i] == j){
        if(counts[j] != 1){
          Z[i,j] = Z[i,j]/(counts[j] - 1)
          a[i] = Z[i,j]
          Z[i,j] = max(Z[i,]) + 1
        }
      }
      else{
        Z[i,j] = Z[i,j]/counts[j]
      }
    }
    if(counts[vm[i]] != 1){
      b[i] = min(Z[i,])
      si.obj[i] = (b[i]-a[i])/max(a[i], b[i])
    }
    else
      si.obj[i] = 0
  }
  idx1 <- mean(si.obj)
  idx2 <- vector(length = n, mode = "numeric")
  Y <- rep(0, n)
  for(i in 1:n)
    Y[i] <- (max(U[i,]) - max(U[i,][-(which.max(U[i,]))]))^av
  idx2 <- sum(Y * si.obj)/sum(Y)
  result = list()
  result$si.obj <- si.obj
  result$sih <- idx1
  if(tidx == "f")
    result$sif <- idx2
  else if(tidx == "e")
    result$sif.e <- idx2
  else if(tidx == "g")
    result$sif.g <- idx2
  return(result)
}
idx<-si(fcm.res)
print(idx)
mean(idx$si.obj) 
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
