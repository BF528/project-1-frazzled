#!/usr/bin/env Rscript
​
​
#installing R-package
#install.packages("ggplot2")
#install.packages("ggpubr")
library(affy)
library(affyPLM)
library(sva)
library(AnnotationDbi)
library(hgu133plus2.db)
#library(ggplot2)
#library(ggarrange)
library(ggpubr)
library(dplyr)
​
​
# reading the csv file/data into final_3
#Set row names
#rownames(final_data3)=final_data3[,1];
final_3 <- read.csv("/projectnb/bf528/users/frazzled/project_1/final_expression.csv")
rownames(final_3)<-final_3[,1]
final_3=final_3[,-1] 
​
final_data3 <- final_3
#4.1 implementing first filter
# at least 20% expressed gene to the value of log2(15) 
final_data3=final_data3 #for statistical computation
​
final_data_1 <- function(fdata){
filterdata<-fdata
N<-nrow(fdata)
column_number<-ncol(fdata)
logval<-log2(15)
logpercentofgenes<- (column_number)*0.2
filter1result<-c(); #Empty list to store the genes that pass the log(15) test
counter=0;
index=1;
​
for (x in 1:N)
{
  counter=0
  for (y in 1:column_number){
    if(fdata[x,y]>logval){
      counter=counter+1
    }
  }
  if(counter>(logpercentofgenes)){
    filter1result[index]<-x;
    index<-index+1;
  }
}
return(filter1result)
}
finalfilter_1 <-final_data_1(final_data3)
​
#chi-squared-test/t-test
#filter 2
​
T_function_Chisq<- function(Cdata) {
  N<-nrow(Cdata)
  column_number<-ncol(Cdata)
  filter2result<-c();
  counter=0;
  index=1;
  degfreedom<-column_number-1
  alpha<-0.01
  chisq_low<-qchisq(alpha/2, degfreedom);
  chisq_high<-qchisq((1-alpha)/2, degfreedom);
  
  # T=(N−1)(s/σ0)2
  stdsample<-sd(as.matrix(Cdata))
  for (x in 1:N){
    T_test<-(column_number-1)*((stdsample/sd(Cdata[x,])**2))
    if ((T_test<chisq_low) || (T_test>chisq_high)){
      filter2result[index]<-x
      index=index+1
    }
    
    
  }
  return(filter2result)
  
}
finalfilter_2<-T_function_Chisq(final_data3)
​
​
#covariant test
#filter 3 
​
convariantfunc<-function(Cdata){
  N<-nrow(Cdata)
  column_number<-ncol(Cdata)
  logval<-log2(15)
  logpercentofgenes<- (column_number)*0.2
  filter3result<-c();
  counter=0;
  index=1;
  for (x in 1:N){
    std<-sd(Cdata[x,])
    meandata<-mean(as.numeric(Cdata[x,]))
    cvvalue<-std/meandata
    if (cvvalue>0.186){
      filter3result[index]<-x
      index=index+1
    }
      
      
  }
  return(filter3result)
  }
  finalfilter_3<-convariantfunc(final_data3)
  
  
#intersect genes passes all the three filters
passall3<-intersect(intersect(finalfilter_1,finalfilter_2), finalfilter_3)
passall3<-as.data.frame(passall3)
​
final_all3<- final_data3[passall3[[1]],]
datapassfilter2<- final_data3[finalfilter_2,]
​
write.csv(final_all3,"/projectnb/bf528/users/frazzled/project_1/final_Noise_filter.csv")
print("number of final that pass all the filters")
print(count(final_all3))
​
#Section_5
#Cluster the dataset of patients 
set.seed(112)
Clusteringdata<-as.data.frame(t(final_all3)) 
#Clusteringdata
Ndata<-as.data.frame(scale(Clusteringdata))
Mdistance<-dist(Ndata,method = "euclidean") #distance matrix 
Cluster<-hclust(Mdistance,method = "average") # method to cluster
plot(Cluster)
​
#cuting the cluster in two groups 
Newclustercut<- cutree(Cluster, k=2)
plot(Cluster)
rect.hclust(Cluster, k=2,border = 2:4)
table(Newclustercut)
  
#heatmap of the gene expression 
#using data from metadata read above cluster
#read meta data file, /project/bf528/project_1/doc/proj_metadata.csv
matadata <- read.csv("/project/bf528/project_1/doc/proj_metadata.csv")
​
#heatmap(as.matrix(t(datawith_cluster[,-1])), #ColSideColors=red, label=names(Clusteringdata))
#function to set the subtype colors
setcolsidecolors<-function(typedata){
  colside_colors=c();#nrow gave errors due to data type of form vector
  u_limit=length(typedata);
  for (i in 1:u_limit)
  {
    if (typedata[i]=="C3"){
      colside_colors[i]<-"red"
    }else{
      colside_colors[i]<-"blue"
    }
  }
  return(colside_colors)
}
​
cancer_mol_subtype<-matadata[ ,"cit.coloncancermolecularsubtype"]
heatmapcolors<-setcolsidecolors(cancer_mol_subtype)
datawith_cluster<-mutate(as.data.frame(Clusteringdata), cluster=Newclustercut)
nameheatmap<-names(Clusteringdata)
heatmap(as.matrix(t(datawith_cluster[,-1])), ColSideColors=heatmapcolors, label=nameheatmap)
​
#5.4
ttestdata<-as.matrix(final_all3)
Newclustercut<-data.frame(cluster = Newclustercut)
​
p_values_run<-function (ttestdata)t.test(ttestdata[Newclustercut['cluster']==1],ttestdata[Newclustercut['cluster']==2], var.equal = T)$p.value
T_Stat_run<- function (ttestdata)t.test(ttestdata[Newclustercut['cluster']==1],ttestdata [Newclustercut['cluster']==2], var.equal = T)$statistic
p_value<-apply(ttestdata,1,p_values_run)
T_stat<-apply(ttestdata,1,T_Stat_run)
​
padjust.value<-p.adjust(p_value,method = "fdr")
  
result_section5<-data.frame(t.stat=T_stat, p_val=p_value, p.adjust=padjust.value)
​
result_section5<-arrange(result_section5,desc(result_section5["t.stat"]),by_group = FALSE)
​
pvalue0.05<-p_value<0.05 
pvalue0.05<-pvalue0.05[pvalue0.05==TRUE]
length(pvalue0.05)
​
​
#5.6 - number of expressed gene at p<0.05
​
datapassfilter2<- final_data3[finalfilter_2,]
p_values_filter2<-function (tdata)t.test(tdata[Newclustercut['cluster']==1],tdata[Newclustercut['cluster']==2], var.equal = T)$p.value
T_Stat_filter2<- function (tdata)t.test(tdata[Newclustercut['cluster']==1],tdata [Newclustercut['cluster']==2], var.equal = T)$statistic
p_value_f2<-apply(datapassfilter2,1,p_values_filter2)
T_stat_f2<-apply(datapassfilter2,1,T_Stat_filter2)
​
padjust.value_f2<-p.adjust(p_value_f2,method = "fdr")
​
result_section5.6<-data.frame(t.stat=T_stat_f2, p_val=p_value_f2, p.adjust=padjust.value_f2)
​
result_section5.6<-arrange(result_section5.6,desc(result_section5.6["t.stat"]),by_group = FALSE)
​
write.csv(result_section5.6,"/projectnb/bf528/users/frazzled/project_1/final_section5.6_result.csv")
​
write.csv(result_section5,"/projectnb/bf528/users/frazzled/project_1/final_section5_result.csv")
