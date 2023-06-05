library(ade4)
library(tidyverse)
library(data.table)
library(VIM)
library(car)
library(magrittr)
library(plyr)
# library(dplyr)
library (meta)
library(reshape2)
library(ggplot2)
library(patchwork)
install.packages("plyr")
#修改药物名称首字母大写，antibiotic列，修改药物赋值2变为0
group_hsct_level7<-read.csv("F:/数据/document/学校/硕士/MMF/meta/HCT gvhd/srr/level-7-hsct.csv",header=T)



# hsct --------------------------------------------------------------------

control<-group_hsct_level7[which(group_hsct_level7$disease == "C"),]
control_0<-control
control_0[is.na(control_0)]<-0

sum_abundance<-apply(control_0[,2:1417],1,sum)

control_rel<-array(0,dim=c(838,1416))
for (j in 1:1416) {
  for (i in 1:838){
    control_rel[i,j]<-control_0[i,j+1]/sum_abundance[i]*100
  }
}
control_rel <- as.data.frame(control_rel)

sample<-control_0$index
species<-colnames(control_0)
colnames(control_rel)<-species[2:1417]
rownames(control_rel)<-sample

control_rel_fil_10 <- filtering_taxonomy(control_rel,0.00000001,10)
# control_rel_fil_20 <- filtering_taxonomy(control_rel,0.00000001,20)

summary_taxnomy_filtering_10<- read.table("summary_taxonomy_filtering.txt",sep="\t",header=T)
filtered_taxonomy_10<-read.table("filtered_taxonomy.txt",sep="\t",header=T)

filtered_taxonomy_10_1<-cbind(filtered_taxonomy_10,index=rownames(filtered_taxonomy_10))

metadata<-as.data.frame(cbind(index=sample,control[,1418:1482]))
control_reg<-merge(filtered_taxonomy_10_1,metadata,by="index")
rownames(control_reg)<-sample
control_reg_name<-control_reg
colnames(control_reg_name)[115]<-"amplicon"
colnames(control_reg_name)[116]<-"sex"
colnames(control_reg_name)[119]<-"scs"
colnames(control_reg_name)[120]<-"DrH"





## single drug analysis ----------------------------------------------------



number_hsct<-array(length(unique(control_reg_name$subject)),dim=c(104,48))
rownames(number_hsct)<-colnames(control_reg_name[,2:105])
colnames(number_hsct)<-colnames(control_reg_name[,121:168])
nonzero_hsct<-array(NA,dim=c(104,48))
rownames(nonzero_hsct)<-colnames(control_reg_name[,2:105])
colnames(nonzero_hsct)<-colnames(control_reg_name[,121:168])
users_hsct<-array(NA,dim=c(104,48))
rownames(users_hsct)<-colnames(control_reg_name[,2:105])
colnames(users_hsct)<-colnames(control_reg_name[,121:168])
nonusers_hsct<-array(NA,dim=c(104,48))
rownames(nonusers_hsct)<-colnames(control_reg_name[,2:105])
colnames(nonusers_hsct)<-colnames(control_reg_name[,121:168])

estimate<-array(0,dim=c(104,48))
rownames(estimate)<-colnames(control_reg_name[,2:105])
colnames(estimate)<-colnames(control_reg_name[,121:168])
Std_Error<-array(0,dim=c(104,48))
rownames(Std_Error)<-colnames(control_reg_name[,2:105])
colnames(Std_Error)<-colnames(control_reg_name[,121:168])
t_value<-array(0,dim=c(104,48))
rownames(t_value)<-colnames(control_reg_name[,2:105])
colnames(t_value)<-colnames(control_reg_name[,121:168])
p_value<-array(0,dim=c(104,48))
rownames(p_value)<-colnames(control_reg_name[,2:105])
colnames(p_value)<-colnames(control_reg_name[,121:168])



for (x in 1:104) {
  for (y in 1:48){
    
    genus_name<-colnames(control_reg_name[,2:105])
    drug_name<-colnames(control_reg_name[,121:168])
    
    temp<-select(control_reg_name,genus_name[x],subject,flag,time,amplicon,sex,age,diagnosis,scs,DrH,drug_name[y])
    temp$species<-temp[,1]
    temp_nonzero<- temp |> filter(species!=0) 
    nonzero_hsct[x,y] <-length(unique(temp_nonzero$subject))
    temp$drug<-temp[,(ncol(temp)-1)]
    temp_users<-temp[complete.cases(temp$drug),] |> filter(drug==1) 
    users_hsct[x,y] <- length(unique(temp_users$subject))
    temp_nonusers<-temp[complete.cases(temp$drug),] |> filter(drug==0) 
    nonusers_hsct[x,y] <- length(unique(temp_nonusers$subject))
    
    temp<-temp[complete.cases(temp),]
    drug<-temp[,ncol(temp)]
    if(length(unique(drug))==2){
      if(length(unique(temp$time))!=1 & length(unique(temp$amplicon))!=1 & length(unique(temp$sex))!=1 & length(unique(temp$scs))!=1){
        temp1<-lm(temp[,1] ~ time + amplicon + sex + age + diagnosis + scs + DrH + drug,temp)
        if (rownames(summary(temp1)$coefficients)[nrow(summary(temp1)$coefficients)] == "drug"){
          summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),1]->estimate[x,y]
          summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),2]->Std_Error[x,y]
          summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),3]->t_value[x,y]
          summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),4]->p_value[x,y]
        }else{
          estimate[x,y]<-NA
          Std_Error[x,y]<-NA
          t_value[x,y]<-NA
          p_value[x,y]<-NA
        }
      }else if(length(unique(temp$time))==1 & length(unique(temp$amplicon))!=1 & length(unique(temp$sex))!=1 & length(unique(temp$scs))!=1){
        temp1<-lm(temp[,1] ~ amplicon + sex + age + diagnosis + scs + DrH + drug,temp)
        if (rownames(summary(temp1)$coefficients)[nrow(summary(temp1)$coefficients)] == "drug"){
          estimate[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),1]
          Std_Error[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),2]
          t_value[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),3]
          p_value[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),4]
        }else{
          estimate[x,y]<-NA
          Std_Error[x,y]<-NA
          t_value[x,y]<-NA
          p_value[x,y]<-NA
        }
      }else if(length(unique(temp$time))!=1 & length(unique(temp$amplicon))==1 & length(unique(temp$sex))!=1 & length(unique(temp$scs))!=1){
        temp1<-lm(temp[,1] ~ time + sex + age + diagnosis + scs + DrH + drug,temp)
        if (rownames(summary(temp1)$coefficients)[nrow(summary(temp1)$coefficients)] == "drug"){
          estimate[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),1]
          Std_Error[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),2]
          t_value[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),3]
          p_value[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),4]
        }else{
          estimate[x,y]<-NA
          Std_Error[x,y]<-NA
          t_value[x,y]<-NA
          p_value[x,y]<-NA
        }
      }else if(length(unique(temp$time))!=1 & length(unique(temp$amplicon))!=1 & length(unique(temp$sex))==1 & length(unique(temp$scs))!=1){
        temp1<-lm(temp[,1] ~ time  + age + diagnosis + scs + DrH + drug,temp)
        if (rownames(summary(temp1)$coefficients)[nrow(summary(temp1)$coefficients)] == "drug"){
          estimate[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),1]
          Std_Error[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),2]
          t_value[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),3]
          p_value[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),4]
        }else{
          estimate[x,y]<-NA
          Std_Error[x,y]<-NA
          t_value[x,y]<-NA
          p_value[x,y]<-NA
        }
      }else if(length(unique(temp$time))!=1 & length(unique(temp$amplicon))!=1 & length(unique(temp$sex))!=1 & length(unique(temp$scs))==1){
        temp1<-lm(temp[,1] ~ time + sex + age + diagnosis + DrH + drug,temp)
        if (rownames(summary(temp1)$coefficients)[nrow(summary(temp1)$coefficients)] == "drug"){
          estimate[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),1]
          Std_Error[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),2]
          t_value[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),3]
          p_value[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),4]
        }else{
          estimate[x,y]<-NA
          Std_Error[x,y]<-NA
          t_value[x,y]<-NA
          p_value[x,y]<-NA
        }
      }else if(length(unique(temp$time))==1 & length(unique(temp$amplicon))==1 & length(unique(temp$sex))!=1 & length(unique(temp$scs))!=1){
        temp1<-lm(temp[,1] ~ sex + age + diagnosis + scs + DrH + drug,temp)
        if (rownames(summary(temp1)$coefficients)[nrow(summary(temp1)$coefficients)] == "drug"){
          estimate[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),1]
          Std_Error[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),2]
          t_value[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),3]
          p_value[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),4]
        }else{
          estimate[x,y]<-NA
          Std_Error[x,y]<-NA
          t_value[x,y]<-NA
          p_value[x,y]<-NA
        }
      }else if(length(unique(temp$time))==1 & length(unique(temp$amplicon))!=1 & length(unique(temp$sex))==1 & length(unique(temp$scs))!=1){
        temp1<-lm(temp[,1] ~ amplicon + age + diagnosis + scs + DrH + drug,temp)
        if (rownames(summary(temp1)$coefficients)[nrow(summary(temp1)$coefficients)] == "drug"){
          estimate[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),1]
          Std_Error[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),2]
          t_value[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),3]
          p_value[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),4]
        }else{
          estimate[x,y]<-NA
          Std_Error[x,y]<-NA
          t_value[x,y]<-NA
          p_value[x,y]<-NA
        }
      }else if(length(unique(temp$time))==1 & length(unique(temp$amplicon))!=1 & length(unique(temp$sex))!=1 & length(unique(temp$scs))==1){
        temp1<-lm(temp[,1] ~ amplicon + sex + age + diagnosis + DrH + drug,temp)
        if (rownames(summary(temp1)$coefficients)[nrow(summary(temp1)$coefficients)] == "drug"){
          estimate[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),1]
          Std_Error[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),2]
          t_value[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),3]
          p_value[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),4]
        }else{
          estimate[x,y]<-NA
          Std_Error[x,y]<-NA
          t_value[x,y]<-NA
          p_value[x,y]<-NA
        }
      }else if(length(unique(temp$time))!=1 & length(unique(temp$amplicon))==1 & length(unique(temp$sex))==1 & length(unique(temp$scs))!=1){
        temp1<-lm(temp[,1] ~ time + age + diagnosis + scs + DrH + drug,temp)
        if (rownames(summary(temp1)$coefficients)[nrow(summary(temp1)$coefficients)] == "drug"){
          estimate[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),1]
          Std_Error[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),2]
          t_value[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),3]
          p_value[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),4]
        }else{
          estimate[x,y]<-NA
          Std_Error[x,y]<-NA
          t_value[x,y]<-NA
          p_value[x,y]<-NA
        }
      }else if(length(unique(temp$time))!=1 & length(unique(temp$amplicon))==1 & length(unique(temp$sex))!=1 & length(unique(temp$scs))==1){
        temp1<-lm(temp[,1] ~ time + sex + age + diagnosis + DrH + drug,temp)
        if (rownames(summary(temp1)$coefficients)[nrow(summary(temp1)$coefficients)] == "drug"){
          estimate[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),1]
          Std_Error[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),2]
          t_value[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),3]
          p_value[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),4]
        }else{
          estimate[x,y]<-NA
          Std_Error[x,y]<-NA
          t_value[x,y]<-NA
          p_value[x,y]<-NA
        }
      }else if(length(unique(temp$time))!=1 & length(unique(temp$amplicon))!=1 & length(unique(temp$sex))==1 & length(unique(temp$scs))==1){
        temp1<-lm(temp[,1] ~ amplicon + time + age + diagnosis + DrH + drug,temp)
        if (rownames(summary(temp1)$coefficients)[nrow(summary(temp1)$coefficients)] == "drug"){
          estimate[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),1]
          Std_Error[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),2]
          t_value[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),3]
          p_value[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),4]
        }else{
          estimate[x,y]<-NA
          Std_Error[x,y]<-NA
          t_value[x,y]<-NA
          p_value[x,y]<-NA
        }
      }else if(length(unique(temp$time))==1 & length(unique(temp$amplicon))==1 & length(unique(temp$sex))==1 & length(unique(temp$scs))!=1){
        temp1<-lm(temp[,1] ~ age + diagnosis + scs + DrH + drug,temp)
        if (rownames(summary(temp1)$coefficients)[nrow(summary(temp1)$coefficients)] == "drug"){
          estimate[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),1]
          Std_Error[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),2]
          t_value[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),3]
          p_value[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),4]
        }else{
          estimate[x,y]<-NA
          Std_Error[x,y]<-NA
          t_value[x,y]<-NA
          p_value[x,y]<-NA
        }
      }else if(length(unique(temp$time))==1 & length(unique(temp$amplicon))==1 & length(unique(temp$sex))!=1 & length(unique(temp$scs))==1){
        temp1<-lm(temp[,1] ~ sex + age + diagnosis + DrH + drug,temp)
        if (rownames(summary(temp1)$coefficients)[nrow(summary(temp1)$coefficients)] == "drug"){
          estimate[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),1]
          Std_Error[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),2]
          t_value[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),3]
          p_value[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),4]
        }else{
          estimate[x,y]<-NA
          Std_Error[x,y]<-NA
          t_value[x,y]<-NA
          p_value[x,y]<-NA
        }
      }else if(length(unique(temp$time))==1 & length(unique(temp$amplicon))!=1 & length(unique(temp$sex))==1 & length(unique(temp$scs))==1){
        temp1<-lm(temp[,1] ~ amplicon + age + diagnosis + DrH + drug,temp)
        if (rownames(summary(temp1)$coefficients)[nrow(summary(temp1)$coefficients)] == "drug"){
          estimate[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),1]
          Std_Error[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),2]
          t_value[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),3]
          p_value[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),4]
        }else{
          estimate[x,y]<-NA
          Std_Error[x,y]<-NA
          t_value[x,y]<-NA
          p_value[x,y]<-NA
        }
      }else if(length(unique(temp$time))!=1 & length(unique(temp$amplicon))==1 & length(unique(temp$sex))==1 & length(unique(temp$scs))==1){
        temp1<-lm(temp[,1] ~ time + age + diagnosis + DrH + drug,temp)
        if (rownames(summary(temp1)$coefficients)[nrow(summary(temp1)$coefficients)] == "drug"){
          estimate[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),1]
          Std_Error[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),2]
          t_value[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),3]
          p_value[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),4]
        }else{
          estimate[x,y]<-NA
          Std_Error[x,y]<-NA
          t_value[x,y]<-NA
          p_value[x,y]<-NA
        }
      }else{
        temp1<-lm(temp[,1] ~ age + diagnosis + DrH + drug,temp)
        if (rownames(summary(temp1)$coefficients)[nrow(summary(temp1)$coefficients)] == "drug"){
          estimate[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),1]
          Std_Error[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),2]
          t_value[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),3]
          p_value[x,y]<-summary(temp1)$coefficients[nrow(summary(temp1)$coefficients),4]
        }else{
          estimate[x,y]<-NA
          Std_Error[x,y]<-NA
          t_value[x,y]<-NA
          p_value[x,y]<-NA
        }
      }
    }else{
      estimate[x,y]<-NA
      Std_Error[x,y]<-NA
      t_value[x,y]<-NA
      p_value[x,y]<-NA
    }
  }
}

number_hsct<-as.data.frame(number_hsct)
nonzero_hsct<-as.data.frame(nonzero_hsct)
users_hsct<-as.data.frame(users_hsct)
nonusers_hsct<-as.data.frame(nonusers_hsct)

estimate<-as.data.frame(estimate)
p_value<-as.data.frame(p_value)
Std_Error<-as.data.frame(Std_Error)
t_value<-as.data.frame(t_value)


## multiple drug analysis --------------------------------------------------

number_mul_hsct<-array(length(unique(control_reg_name$subject)),dim=c(104,47))
rownames(number_mul_hsct)<-colnames(control_reg_name[,2:105])
colnames(number_mul_hsct)<-colnames(control_reg_name[,c(121:142,144:168)])
nonzero_mul_hsct<-array(NA,dim=c(104,47))
rownames(nonzero_mul_hsct)<-colnames(control_reg_name[,2:105])
colnames(nonzero_mul_hsct)<-colnames(control_reg_name[,c(121:142,144:168)])
users_mul_hsct<-array(NA,dim=c(104,47))
rownames(users_mul_hsct)<-colnames(control_reg_name[,2:105])
colnames(users_mul_hsct)<-colnames(control_reg_name[,c(121:142,144:168)])
nonusers_mul_hsct<-array(NA,dim=c(104,47))
rownames(nonusers_mul_hsct)<-colnames(control_reg_name[,2:105])
colnames(nonusers_mul_hsct)<-colnames(control_reg_name[,c(121:142,144:168)])


estimate_mul<-array(0,dim=c(104,47))
rownames(estimate_mul)<-colnames(control_reg_name[,2:105])
colnames(estimate_mul)<-colnames(control_reg_name[,c(121:142,144:168)])
Std_Error_mul<-array(0,dim=c(104,47))
rownames(Std_Error_mul)<-colnames(control_reg_name[,2:105])
colnames(Std_Error_mul)<-colnames(control_reg_name[,c(121:142,144:168)])
t_value_mul<-array(0,dim=c(104,47))
rownames(t_value_mul)<-colnames(control_reg_name[,2:105])
colnames(t_value_mul)<-colnames(control_reg_name[,c(121:142,144:168)])
p_value_mul<-array(0,dim=c(104,47))
rownames(p_value_mul)<-colnames(control_reg_name[,2:105])
colnames(p_value_mul)<-colnames(control_reg_name[,c(121:142,144:168)])


for (x in 1:104) {
  for (y in 1:47){
    genus_name<-colnames(control_reg_name[,2:105])
    drug_name<-colnames(control_reg_name[,c(121:142,144:168)])
    
    temp_mul<-select(control_reg_name,genus_name[x],subject,flag,time,amplicon,sex,age,diagnosis,scs,DrH,Antibiotic,drug_name[y])
    temp_mul$species<-temp_mul[,1]
    temp_mul_nonzero<- temp_mul |> filter(species!=0) 
    nonzero_mul_hsct[x,y] <-length(unique(temp_mul_nonzero$subject))
    temp_mul$drug<-temp_mul[,(ncol(temp_mul)-1)]
    temp_mul_users<-temp_mul[complete.cases(temp_mul$drug),] |> filter(drug==1) 
    users_mul_hsct[x,y] <- length(unique(temp_mul_users$subject))
    temp_mul_nonusers<-temp_mul[complete.cases(temp_mul$drug),] |> filter(drug==0) 
    nonusers_mul_hsct[x,y] <- length(unique(temp_mul_nonusers$subject))
    
    
    temp_mul<-temp_mul[complete.cases(temp_mul),]
    drug<-temp_mul[,ncol(temp_mul)]
    Antibiotic<-temp_mul[,"Antibiotic"]
    if(length(unique(drug))==2 & length(unique(Antibiotic))==2){
      if(length(unique(temp_mul$time))!=1 & length(unique(temp_mul$amplicon))!=1 & length(unique(temp_mul$sex))!=1 & length(unique(temp_mul$scs))!=1){
        temp_mul1<-lm(temp_mul[,1] ~ time + amplicon + sex + age + diagnosis + scs + DrH + Antibiotic +drug,temp_mul)
        if (rownames(summary(temp_mul1)$coefficients)[nrow(summary(temp_mul1)$coefficients)] == "drug"){
          summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),1]->estimate_mul[x,y]
          summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),2]->Std_Error_mul[x,y]
          summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),3]->t_value_mul[x,y]
          summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),4]->p_value_mul[x,y]
        }else{
          estimate_mul[x,y]<-NA
          Std_Error_mul[x,y]<-NA
          t_value_mul[x,y]<-NA
          p_value_mul[x,y]<-NA
        }
      }else if(length(unique(temp_mul$time))==1 & length(unique(temp_mul$amplicon))!=1 & length(unique(temp_mul$sex))!=1 & length(unique(temp_mul$scs))!=1){
        temp_mul1<-lm(temp_mul[,1] ~ amplicon + sex + age + diagnosis + scs + DrH + Antibiotic+ drug,temp_mul)
        if (rownames(summary(temp_mul1)$coefficients)[nrow(summary(temp_mul1)$coefficients)] == "drug"){
          estimate_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),1]
          Std_Error_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),2]
          t_value_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),3]
          p_value_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),4]
        }else{
          estimate_mul[x,y]<-NA
          Std_Error_mul[x,y]<-NA
          t_value_mul[x,y]<-NA
          p_value_mul[x,y]<-NA
        }
      }else if(length(unique(temp_mul$time))!=1 & length(unique(temp_mul$amplicon))==1 & length(unique(temp_mul$sex))!=1 & length(unique(temp_mul$scs))!=1){
        temp_mul1<-lm(temp_mul[,1] ~ time + sex + age + diagnosis + scs + DrH + Antibiotic+ drug,temp_mul)
        if (rownames(summary(temp_mul1)$coefficients)[nrow(summary(temp_mul1)$coefficients)] == "drug"){
          estimate_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),1]
          Std_Error_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),2]
          t_value_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),3]
          p_value_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),4]
        }else{
          estimate_mul[x,y]<-NA
          Std_Error_mul[x,y]<-NA
          t_value_mul[x,y]<-NA
          p_value_mul[x,y]<-NA
        }
      }else if(length(unique(temp_mul$time))!=1 & length(unique(temp_mul$amplicon))!=1 & length(unique(temp_mul$sex))==1 & length(unique(temp_mul$scs))!=1){
        temp_mul1<-lm(temp_mul[,1] ~ time  + age + diagnosis + scs + DrH + Antibiotic+ drug,temp_mul)
        if (rownames(summary(temp_mul1)$coefficients)[nrow(summary(temp_mul1)$coefficients)] == "drug"){
          estimate_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),1]
          Std_Error_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),2]
          t_value_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),3]
          p_value_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),4]
        }else{
          estimate_mul[x,y]<-NA
          Std_Error_mul[x,y]<-NA
          t_value_mul[x,y]<-NA
          p_value_mul[x,y]<-NA
        }
      }else if(length(unique(temp_mul$time))!=1 & length(unique(temp_mul$amplicon))!=1 & length(unique(temp_mul$sex))!=1 & length(unique(temp_mul$scs))==1){
        temp_mul1<-lm(temp_mul[,1] ~ time + sex + age + diagnosis + DrH+ Antibiotic + drug,temp_mul)
        if (rownames(summary(temp_mul1)$coefficients)[nrow(summary(temp_mul1)$coefficients)] == "drug"){
          estimate_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),1]
          Std_Error_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),2]
          t_value_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),3]
          p_value_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),4]
        }else{
          estimate_mul[x,y]<-NA
          Std_Error_mul[x,y]<-NA
          t_value_mul[x,y]<-NA
          p_value_mul[x,y]<-NA
        }
      }else if(length(unique(temp_mul$time))==1 & length(unique(temp_mul$amplicon))==1 & length(unique(temp_mul$sex))!=1 & length(unique(temp_mul$scs))!=1){
        temp_mul1<-lm(temp_mul[,1] ~ sex + age + diagnosis + scs + DrH + Antibiotic+ drug,temp_mul)
        if (rownames(summary(temp_mul1)$coefficients)[nrow(summary(temp_mul1)$coefficients)] == "drug"){
          estimate_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),1]
          Std_Error_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),2]
          t_value_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),3]
          p_value_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),4]
        }else{
          estimate_mul[x,y]<-NA
          Std_Error_mul<-NA
          t_value_mul[x,y]<-NA
          p_value_mul[x,y]<-NA
        }
      }else if(length(unique(temp_mul$time))==1 & length(unique(temp_mul$amplicon))!=1 & length(unique(temp_mul$sex))==1 & length(unique(temp_mul$scs))!=1){
        temp_mul1<-lm(temp_mul[,1] ~ amplicon + age + diagnosis + scs + DrH + Antibiotic+ drug,temp_mul)
        if (rownames(summary(temp_mul1)$coefficients)[nrow(summary(temp_mul1)$coefficients)] == "drug"){
          estimate_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),1]
          Std_Error_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),2]
          t_value_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),3]
          p_value_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),4]
        }else{
          estimate_mul[x,y]<-NA
          Std_Error_mul[x,y]<-NA
          t_value_mul[x,y]<-NA
          p_value_mul[x,y]<-NA
        }
      }else if(length(unique(temp_mul$time))==1 & length(unique(temp_mul$amplicon))!=1 & length(unique(temp_mul$sex))!=1 & length(unique(temp_mul$scs))==1){
        temp_mul1<-lm(temp_mul[,1] ~ amplicon + sex + age + diagnosis + DrH + Antibiotic+ drug,temp_mul)
        if (rownames(summary(temp_mul1)$coefficients)[nrow(summary(temp_mul1)$coefficients)] == "drug"){
          estimate_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),1]
          Std_Error_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),2]
          t_value_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),3]
          p_value_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),4]
        }else{
          estimate_mul[x,y]<-NA
          Std_Error_mul[x,y]<-NA
          t_value_mul[x,y]<-NA
          p_value_mul[x,y]<-NA
        }
      }else if(length(unique(temp_mul$time))!=1 & length(unique(temp_mul$amplicon))==1 & length(unique(temp_mul$sex))==1 & length(unique(temp_mul$scs))!=1){
        temp_mul1<-lm(temp_mul[,1] ~ time + age + diagnosis + scs + DrH + Antibiotic+ drug,temp_mul)
        if (rownames(summary(temp_mul1)$coefficients)[nrow(summary(temp_mul1)$coefficients)] == "drug"){
          estimate_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),1]
          Std_Error_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),2]
          t_value_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),3]
          p_value_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),4]
        }else{
          estimate_mul[x,y]<-NA
          Std_Error_mul[x,y]<-NA
          t_value_mul[x,y]<-NA
          p_value_mul[x,y]<-NA
        }
      }else if(length(unique(temp_mul$time))!=1 & length(unique(temp_mul$amplicon))==1 & length(unique(temp_mul$sex))!=1 & length(unique(temp_mul$scs))==1){
        temp_mul1<-lm(temp_mul[,1] ~ time + sex + age + diagnosis + DrH + Antibiotic+ drug,temp_mul)
        if (rownames(summary(temp_mul1)$coefficients)[nrow(summary(temp_mul1)$coefficients)] == "drug"){
          estimate_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),1]
          Std_Error_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),2]
          t_value_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),3]
          p_value_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),4]
        }else{
          estimate_mul[x,y]<-NA
          Std_Error_mul[x,y]<-NA
          t_value_mul[x,y]<-NA
          p_value_mul[x,y]<-NA
        }
      }else if(length(unique(temp_mul$time))!=1 & length(unique(temp_mul$amplicon))!=1 & length(unique(temp_mul$sex))==1 & length(unique(temp_mul$scs))==1){
        temp_mul1<-lm(temp_mul[,1] ~ amplicon + time + age + diagnosis + DrH + Antibiotic+ drug,temp_mul)
        if (rownames(summary(temp_mul1)$coefficients)[nrow(summary(temp_mul1)$coefficients)] == "drug"){
          estimate_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),1]
          Std_Error_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),2]
          t_value_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),3]
          p_value_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),4]
        }else{
          estimate_mul[x,y]<-NA
          Std_Error_mul[x,y]<-NA
          t_value_mul[x,y]<-NA
          p_value_mul[x,y]<-NA
        }
      }else if(length(unique(temp_mul$time))==1 & length(unique(temp_mul$amplicon))==1 & length(unique(temp_mul$sex))==1 & length(unique(temp_mul$scs))!=1){
        temp_mul1<-lm(temp_mul[,1] ~ age + diagnosis + scs + DrH + Antibiotic+ drug,temp_mul)
        if (rownames(summary(temp_mul1)$coefficients)[nrow(summary(temp_mul1)$coefficients)] == "drug"){
          estimate_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),1]
          Std_Error_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),2]
          t_value_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),3]
          p_value_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),4]
        }else{
          estimate_mul[x,y]<-NA
          Std_Error_mul[x,y]<-NA
          t_value_mul[x,y]<-NA
          p_value_mul[x,y]<-NA
        }
      }else if(length(unique(temp_mul$time))==1 & length(unique(temp_mul$amplicon))==1 & length(unique(temp_mul$sex))!=1 & length(unique(temp_mul$scs))==1){
        temp_mul1<-lm(temp_mul[,1] ~ sex + age + diagnosis + DrH + Antibiotic+ drug,temp_mul)
        if (rownames(summary(temp_mul1)$coefficients)[nrow(summary(temp_mul1)$coefficients)] == "drug"){
          estimate_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),1]
          Std_Error_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),2]
          t_value_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),3]
          p_value_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),4]
        }else{
          estimate_mul[x,y]<-NA
          Std_Error_mul[x,y]<-NA
          t_value_mul[x,y]<-NA
          p_value_mul[x,y]<-NA
        }
      }else if(length(unique(temp_mul$time))==1 & length(unique(temp_mul$amplicon))!=1 & length(unique(temp_mul$sex))==1 & length(unique(temp_mul$scs))==1){
        temp_mul1<-lm(temp_mul[,1] ~ amplicon + age + diagnosis + DrH + Antibiotic+ drug,temp_mul)
        if (rownames(summary(temp_mul1)$coefficients)[nrow(summary(temp_mul1)$coefficients)] == "drug"){
          estimate_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),1]
          Std_Error_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),2]
          t_value_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),3]
          p_value_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),4]
        }else{
          estimate_mul[x,y]<-NA
          Std_Error_mul[x,y]<-NA
          t_value_mul[x,y]<-NA
          p_value_mul[x,y]<-NA
        }
      }else if(length(unique(temp_mul$time))!=1 & length(unique(temp_mul$amplicon))==1 & length(unique(temp_mul$sex))==1 & length(unique(temp_mul$scs))==1){
        temp_mul1<-lm(temp_mul[,1] ~ time + age + diagnosis + DrH + Antibiotic+ drug,temp_mul)
        if (rownames(summary(temp_mul1)$coefficients)[nrow(summary(temp_mul1)$coefficients)] == "drug"){
          estimate_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),1]
          Std_Error_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),2]
          t_value_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),3]
          p_value_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),4]
        }else{
          estimate_mul[x,y]<-NA
          Std_Error_mul[x,y]<-NA
          t_value_mul[x,y]<-NA
          p_value_mul[x,y]<-NA
        }
      }else{
        temp_mul1<-lm(temp_mul[,1] ~ age + diagnosis + DrH + Antibiotic+ drug,temp_mul)
        if (rownames(summary(temp_mul1)$coefficients)[nrow(summary(temp_mul1)$coefficients)] == "drug"){
          estimate_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),1]
          Std_Error_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),2]
          t_value_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),3]
          p_value_mul[x,y]<-summary(temp_mul1)$coefficients[nrow(summary(temp_mul1)$coefficients),4]
        }else{
          estimate_mul[x,y]<-NA
          Std_Error_mul[x,y]<-NA
          t_value_mul[x,y]<-NA
          p_value_mul[x,y]<-NA
        }
      }
    }else{
      estimate_mul[x,y]<-NA
      Std_Error_mul[x,y]<-NA
      t_value_mul[x,y]<-NA
      p_value_mul[x,y]<-NA
    }
  }
}

number_mul_hsct<-as.data.frame(number_mul_hsct)
nonzero_mul_hsct<-as.data.frame(nonzero_mul_hsct)
users_mul_hsct<-as.data.frame(users_mul_hsct)
nonusers_mul_hsct<-as.data.frame(nonusers_mul_hsct)

estimate_mul<-as.data.frame(estimate_mul)
p_value_mul<-as.data.frame(p_value_mul)
Std_Error_mul<-as.data.frame(Std_Error_mul)
t_value_mul<-as.data.frame(t_value_mul)



# GVHD -------------------------------------------------------------------


gvhd<-group_hsct_level7[which(group_hsct_level7$disease == "G"),]
gvhd_0<-gvhd
gvhd_0[is.na(gvhd_0)]<-0

sum_abundance_gvhd<-apply(gvhd_0[,2:1417],1,sum)

gvhd_rel<-array(0,dim=c(715,1416))
for (j in 1:1416) {
  for (i in 1:715){
    gvhd_rel[i,j]<-gvhd_0[i,j+1]/sum_abundance_gvhd[i]*100
  }
}

gvhd_rel <- as.data.frame(gvhd_rel)
gvhd_rel[is.na(gvhd_rel)]<-0

sample_gvhd<-gvhd_0$index
species_gvhd<-colnames(gvhd_0)
colnames(gvhd_rel)<-species_gvhd[2:1417]
rownames(gvhd_rel)<-sample_gvhd

gvhd_rel_fil_10 <- filtering_taxonomy(gvhd_rel,0.00000001,10)

summary_taxnomy_filtering_10_gvhd<- read.table("summary_taxonomy_filtering.txt",sep="\t",header=T)
filtered_taxonomy_10_gvhd<-read.table("filtered_taxonomy.txt",sep="\t",header=T)

filtered_taxonomy_10_1_gvhd<-cbind(filtered_taxonomy_10_gvhd,index=rownames(filtered_taxonomy_10_gvhd))

metadata_gvhd<-as.data.frame(cbind(gvhd[,1418:1482],index=sample_gvhd))
gvhd_reg<-merge(filtered_taxonomy_10_1_gvhd,metadata_gvhd,by="index")
rownames(gvhd_reg)<-gvhd_reg$index
gvhd_reg_name<-gvhd_reg
colnames(gvhd_reg_name)[122]<-"amplicon"
colnames(gvhd_reg_name)[123]<-"sex"
colnames(gvhd_reg_name)[126]<-"scs"
colnames(gvhd_reg_name)[127]<-"DrH"
colnames(gvhd_reg_name)[120]<-"gclass"


## single drug analysis -------------------------------------------------

number_gvhd<-array(length(unique(gvhd_reg_name$subject)),dim=c(111,48))
rownames(number_gvhd)<-colnames(gvhd_reg_name[,2:112])
colnames(number_gvhd)<-colnames(gvhd_reg_name[,128:175])
nonzero_gvhd<-array(NA,dim=c(111,48))
rownames(nonzero_gvhd)<-colnames(gvhd_reg_name[,2:112])
colnames(nonzero_gvhd)<-colnames(gvhd_reg_name[,128:175])
users_gvhd<-array(NA,dim=c(111,48))
rownames(users_gvhd)<-colnames(gvhd_reg_name[,2:112])
colnames(users_gvhd)<-colnames(gvhd_reg_name[,128:175])
nonusers_gvhd<-array(NA,dim=c(111,48))
rownames(nonusers_gvhd)<-colnames(gvhd_reg_name[,2:112])
colnames(nonusers_gvhd)<-colnames(gvhd_reg_name[,128:175])



estimate_gvhd<-array(0,dim=c(111,48))
rownames(estimate_gvhd)<-colnames(gvhd_reg_name[,2:112])
colnames(estimate_gvhd)<-colnames(gvhd_reg_name[,128:175])
Std_Error_gvhd<-array(0,dim=c(111,48))
rownames(Std_Error_gvhd)<-colnames(gvhd_reg_name[,2:112])
colnames(Std_Error_gvhd)<-colnames(gvhd_reg_name[,128:175])
t_value_gvhd<-array(0,dim=c(111,48))
rownames(t_value_gvhd)<-colnames(gvhd_reg_name[,2:112])
colnames(t_value_gvhd)<-colnames(gvhd_reg_name[,128:175])
p_value_gvhd<-array(0,dim=c(111,48))
rownames(p_value_gvhd)<-colnames(gvhd_reg_name[,2:112])
colnames(p_value_gvhd)<-colnames(gvhd_reg_name[,128:175])


for (x in 1:111) {
  for (y in 1:48){
    genus_name<-colnames(gvhd_reg_name[,2:112])
    drug_name<-colnames(gvhd_reg_name[,128:175])
    
    
    temp_gvhd<-select(gvhd_reg_name,genus_name[x],subject,flag,time,amplicon,gclass,sex,age,diagnosis,scs,DrH,drug_name[y])
    temp_gvhd$species<-temp_gvhd[,1]
    temp_gvhd_nonzero<- temp_gvhd |> filter(species!=0) 
    nonzero_gvhd[x,y] <-length(unique(temp_gvhd_nonzero$subject))
    temp_gvhd$drug<-temp_gvhd[,(ncol(temp_gvhd)-1)]
    temp_gvhd_users<-temp_gvhd[complete.cases(temp_gvhd$drug),] |> filter(drug==1) 
    users_gvhd[x,y] <- length(unique(temp_gvhd_users$subject))
    temp_gvhd_nonusers<-temp_gvhd[complete.cases(temp_gvhd$drug),] |> filter(drug==0) 
    nonusers_gvhd[x,y] <- length(unique(temp_gvhd_nonusers$subject))
    
    
    temp_gvhd<-temp_gvhd[complete.cases(temp_gvhd),]
    drug_gvhd<-temp_gvhd[,ncol(temp_gvhd)]
    if(length(unique(drug_gvhd))==2){
      if(length(unique(temp_gvhd$time))!=1 & length(unique(temp_gvhd$amplicon))!=1 & length(unique(temp_gvhd$sex))!=1 & length(unique(temp_gvhd$scs))!=1){
        temp_gvhd1<-lm(temp_gvhd[,1] ~ time + amplicon + sex + age + gclass+diagnosis + scs + DrH + drug_gvhd,temp_gvhd)
        if (rownames(summary(temp_gvhd1)$coefficients)[nrow(summary(temp_gvhd1)$coefficients)] == "drug_gvhd"){
          summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),1]->estimate_gvhd[x,y]
          summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),2]->Std_Error_gvhd[x,y]
          summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),3]->t_value_gvhd[x,y]
          summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),4]->p_value_gvhd[x,y]
        }else{
          estimate_gvhd[x,y]<-NA
          Std_Error_gvhd[x,y]<-NA
          t_value_gvhd[x,y]<-NA
          p_value_gvhd[x,y]<-NA
        }
      }else if(length(unique(temp_gvhd$time))==1 & length(unique(temp_gvhd$amplicon))!=1 & length(unique(temp_gvhd$sex))!=1 & length(unique(temp_gvhd$scs))!=1){
        temp_gvhd1<-lm(temp_gvhd[,1] ~ amplicon + sex + age + gclass+diagnosis + scs + DrH + drug_gvhd,temp_gvhd)
        if (rownames(summary(temp_gvhd1)$coefficients)[nrow(summary(temp_gvhd1)$coefficients)] == "drug_gvhd"){
          estimate_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),1]
          Std_Error_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),2]
          t_value_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),3]
          p_value_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),4]
        }else{
          estimate_gvhd[x,y]<-NA
          Std_Error_gvhd[x,y]<-NA
          t_value_gvhd[x,y]<-NA
          p_value_gvhd[x,y]<-NA
        }
      }else if(length(unique(temp_gvhd$time))!=1 & length(unique(temp_gvhd$amplicon))==1 & length(unique(temp_gvhd$sex))!=1 & length(unique(temp_gvhd$scs))!=1){
        temp_gvhd1<-lm(temp_gvhd[,1] ~ time + sex + age + gclass+diagnosis + scs + DrH + drug_gvhd,temp_gvhd)
        if (rownames(summary(temp_gvhd1)$coefficients)[nrow(summary(temp_gvhd1)$coefficients)] == "drug_gvhd"){
          estimate_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),1]
          Std_Error_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),2]
          t_value_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),3]
          p_value_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),4]
        }else{
          estimate_gvhd[x,y]<-NA
          Std_Error_gvhd[x,y]<-NA
          t_value_gvhd[x,y]<-NA
          p_value_gvhd[x,y]<-NA
        }
      }else if(length(unique(temp_gvhd$time))!=1 & length(unique(temp_gvhd$amplicon))!=1 & length(unique(temp_gvhd$sex))==1 & length(unique(temp_gvhd$scs))!=1){
        temp_gvhd1<-lm(temp_gvhd[,1] ~ time  + age + gclass+diagnosis + scs + DrH + drug_gvhd,temp_gvhd)
        if (rownames(summary(temp_gvhd1)$coefficients)[nrow(summary(temp_gvhd1)$coefficients)] == "drug_gvhd"){
          estimate_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),1]
          Std_Error_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),2]
          t_value_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),3]
          p_value_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),4]
        }else{
          estimate_gvhd[x,y]<-NA
          Std_Error_gvhd[x,y]<-NA
          t_value_gvhd[x,y]<-NA
          p_value_gvhd[x,y]<-NA
        }
      }else if(length(unique(temp_gvhd$time))!=1 & length(unique(temp_gvhd$amplicon))!=1 & length(unique(temp_gvhd$sex))!=1 & length(unique(temp_gvhd$scs))==1){
        temp_gvhd1<-lm(temp_gvhd[,1] ~ time + sex + age + gclass+diagnosis + DrH + drug_gvhd,temp_gvhd)
        if (rownames(summary(temp_gvhd1)$coefficients)[nrow(summary(temp_gvhd1)$coefficients)] == "drug_gvhd"){
          estimate_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),1]
          Std_Error_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),2]
          t_value_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),3]
          p_value_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),4]
        }else{
          estimate_gvhd[x,y]<-NA
          Std_Error_gvhd[x,y]<-NA
          t_value_gvhd[x,y]<-NA
          p_value_gvhd[x,y]<-NA
        }
      }else if(length(unique(temp_gvhd$time))==1 & length(unique(temp_gvhd$amplicon))==1 & length(unique(temp_gvhd$sex))!=1 & length(unique(temp_gvhd$scs))!=1){
        temp_gvhd1<-lm(temp_gvhd[,1] ~ sex + age + gclass+diagnosis + scs + DrH + drug_gvhd,temp_gvhd)
        if (rownames(summary(temp_gvhd1)$coefficients)[nrow(summary(temp_gvhd1)$coefficients)] == "drug_gvhd"){
          estimate_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),1]
          Std_Error_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),2]
          t_value_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),3]
          p_value_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),4]
        }else{
          estimate_gvhd[x,y]<-NA
          Std_Error_gvhd[x,y]<-NA
          t_value_gvhd[x,y]<-NA
          p_value_gvhd[x,y]<-NA
        }
      }else if(length(unique(temp_gvhd$time))==1 & length(unique(temp_gvhd$amplicon))!=1 & length(unique(temp_gvhd$sex))==1 & length(unique(temp_gvhd$scs))!=1){
        temp_gvhd1<-lm(temp_gvhd[,1] ~ amplicon + age + gclass+diagnosis + scs + DrH + drug_gvhd,temp_gvhd)
        if (rownames(summary(temp_gvhd1)$coefficients)[nrow(summary(temp_gvhd1)$coefficients)] == "drug_gvhd"){
          estimate_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),1]
          Std_Error_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),2]
          t_value_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),3]
          p_value_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),4]
        }else{
          estimate_gvhd[x,y]<-NA
          Std_Error_gvhd[x,y]<-NA
          t_value_gvhd[x,y]<-NA
          p_value_gvhd[x,y]<-NA
        }
      }else if(length(unique(temp_gvhd$time))==1 & length(unique(temp_gvhd$amplicon))!=1 & length(unique(temp_gvhd$sex))!=1 & length(unique(temp_gvhd$scs))==1){
        temp_gvhd1<-lm(temp_gvhd[,1] ~ amplicon + sex + age + gclass+diagnosis + DrH + drug_gvhd,temp_gvhd)
        if (rownames(summary(temp_gvhd1)$coefficients)[nrow(summary(temp_gvhd1)$coefficients)] == "drug_gvhd"){
          estimate_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),1]
          Std_Error_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),2]
          t_value_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),3]
          p_value_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),4]
        }else{
          estimate_gvhd[x,y]<-NA
          Std_Error_gvhd[x,y]<-NA
          t_value_gvhd[x,y]<-NA
          p_value_gvhd[x,y]<-NA
        }
      }else if(length(unique(temp_gvhd$time))!=1 & length(unique(temp_gvhd$amplicon))==1 & length(unique(temp_gvhd$sex))==1 & length(unique(temp_gvhd$scs))!=1){
        temp_gvhd1<-lm(temp_gvhd[,1] ~ time + age + gclass+diagnosis + scs + DrH + drug_gvhd,temp_gvhd)
        if (rownames(summary(temp_gvhd1)$coefficients)[nrow(summary(temp_gvhd1)$coefficients)] == "drug_gvhd"){
          estimate_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),1]
          Std_Error_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),2]
          t_value_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),3]
          p_value_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),4]
        }else{
          estimate_gvhd[x,y]<-NA
          Std_Error_gvhd[x,y]<-NA
          t_value_gvhd[x,y]<-NA
          p_value_gvhd[x,y]<-NA
        }
      }else if(length(unique(temp_gvhd$time))!=1 & length(unique(temp_gvhd$amplicon))==1 & length(unique(temp_gvhd$sex))!=1 & length(unique(temp_gvhd$scs))==1){
        temp_gvhd1<-lm(temp_gvhd[,1] ~ time + sex + age + gclass+diagnosis + DrH + drug_gvhd,temp_gvhd)
        if (rownames(summary(temp_gvhd1)$coefficients)[nrow(summary(temp_gvhd1)$coefficients)] == "drug_gvhd"){
          estimate_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),1]
          Std_Error_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),2]
          t_value_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),3]
          p_value_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),4]
        }else{
          estimate_gvhd[x,y]<-NA
          Std_Error_gvhd[x,y]<-NA
          t_value_gvhd[x,y]<-NA
          p_value_gvhd[x,y]<-NA
        }
      }else if(length(unique(temp_gvhd$time))!=1 & length(unique(temp_gvhd$amplicon))!=1 & length(unique(temp_gvhd$sex))==1 & length(unique(temp_gvhd$scs))==1){
        temp_gvhd1<-lm(temp_gvhd[,1] ~ amplicon + time + age + gclass+diagnosis + DrH + drug_gvhd,temp_gvhd)
        if (rownames(summary(temp_gvhd1)$coefficients)[nrow(summary(temp_gvhd1)$coefficients)] == "drug_gvhd"){
          estimate_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),1]
          Std_Error_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),2]
          t_value_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),3]
          p_value_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),4]
        }else{
          estimate_gvhd[x,y]<-NA
          Std_Error_gvhd[x,y]<-NA
          t_value_gvhd[x,y]<-NA
          p_value_gvhd[x,y]<-NA
        }
      }else if(length(unique(temp_gvhd$time))==1 & length(unique(temp_gvhd$amplicon))==1 & length(unique(temp_gvhd$sex))==1 & length(unique(temp_gvhd$scs))!=1){
        temp_gvhd1<-lm(temp_gvhd[,1] ~ age + gclass+diagnosis + scs + DrH + drug_gvhd,temp_gvhd)
        if (rownames(summary(temp_gvhd1)$coefficients)[nrow(summary(temp_gvhd1)$coefficients)] == "drug_gvhd"){
          estimate_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),1]
          Std_Error_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),2]
          t_value_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),3]
          p_value_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),4]
        }else{
          estimate_gvhd[x,y]<-NA
          Std_Error_gvhd[x,y]<-NA
          t_value_gvhd[x,y]<-NA
          p_value_gvhd[x,y]<-NA
        }
      }else if(length(unique(temp_gvhd$time))==1 & length(unique(temp_gvhd$amplicon))==1 & length(unique(temp_gvhd$sex))!=1 & length(unique(temp_gvhd$scs))==1){
        temp_gvhd1<-lm(temp_gvhd[,1] ~ sex + age + gclass+diagnosis + DrH + drug_gvhd,temp_gvhd)
        if (rownames(summary(temp_gvhd1)$coefficients)[nrow(summary(temp_gvhd1)$coefficients)] == "drug_gvhd"){
          estimate_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),1]
          Std_Error_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),2]
          t_value_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),3]
          p_value_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),4]
        }else{
          estimate_gvhd[x,y]<-NA
          Std_Error_gvhd[x,y]<-NA
          t_value_gvhd[x,y]<-NA
          p_value_gvhd[x,y]<-NA
        }
      }else if(length(unique(temp_gvhd$time))==1 & length(unique(temp_gvhd$amplicon))!=1 & length(unique(temp_gvhd$sex))==1 & length(unique(temp_gvhd$scs))==1){
        temp_gvhd1<-lm(temp_gvhd[,1] ~ amplicon + age + gclass+diagnosis + DrH + drug_gvhd,temp_gvhd)
        if (rownames(summary(temp_gvhd1)$coefficients)[nrow(summary(temp_gvhd1)$coefficients)] == "drug_gvhd"){
          estimate_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),1]
          Std_Error_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),2]
          t_value_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),3]
          p_value_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),4]
        }else{
          estimate_gvhd[x,y]<-NA
          Std_Error_gvhd[x,y]<-NA
          t_value_gvhd[x,y]<-NA
          p_value_gvhd[x,y]<-NA
        }
      }else if(length(unique(temp_gvhd$time))!=1 & length(unique(temp_gvhd$amplicon))==1 & length(unique(temp_gvhd$sex))==1 & length(unique(temp_gvhd$scs))==1){
        temp_gvhd1<-lm(temp_gvhd[,1] ~ time + age + gclass+diagnosis + DrH + drug_gvhd,temp_gvhd)
        if (rownames(summary(temp_gvhd1)$coefficients)[nrow(summary(temp_gvhd1)$coefficients)] == "drug_gvhd"){
          estimate_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),1]
          Std_Error_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),2]
          t_value_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),3]
          p_value_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),4]
        }else{
          estimate_gvhd[x,y]<-NA
          Std_Error_gvhd[x,y]<-NA
          t_value_gvhd[x,y]<-NA
          p_value_gvhd[x,y]<-NA
        }
      }else{
        temp_gvhd1<-lm(temp_gvhd[,1] ~ age + gclass+diagnosis + DrH + drug_gvhd,temp_gvhd)
        if (rownames(summary(temp_gvhd1)$coefficients)[nrow(summary(temp_gvhd1)$coefficients)] == "drug_gvhd"){
          estimate_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),1]
          Std_Error_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),2]
          t_value_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),3]
          p_value_gvhd[x,y]<-summary(temp_gvhd1)$coefficients[nrow(summary(temp_gvhd1)$coefficients),4]
        }else{
          estimate_gvhd[x,y]<-NA
          Std_Error_gvhd[x,y]<-NA
          t_value_gvhd[x,y]<-NA
          p_value_gvhd[x,y]<-NA
        }
      }
    }else{
      estimate_gvhd[x,y]<-NA
      Std_Error_gvhd[x,y]<-NA
      t_value_gvhd[x,y]<-NA
      p_value_gvhd[x,y]<-NA
    }
  }
}

number_gvhd<-as.data.frame(number_gvhd)
nonzero_gvhd<-as.data.frame(nonzero_gvhd)
users_gvhd<-as.data.frame(users_gvhd)
nonusers_gvhd<-as.data.frame(nonusers_gvhd)


estimate_gvhd<-as.data.frame(estimate_gvhd)
p_value_gvhd<-as.data.frame(p_value_gvhd)
Std_Error_gvhd<-as.data.frame(Std_Error_gvhd)
t_value_gvhd<-as.data.frame(t_value_gvhd)



## multiple drug analysis --------------------------------------------------

number_mul_gvhd<-array(length(unique(gvhd_reg_name$subject)),dim=c(111,47))
rownames(number_mul_gvhd)<-colnames(gvhd_reg_name[,2:112])
colnames(number_mul_gvhd)<-colnames(gvhd_reg_name[,c(128:149,151:175)])
nonzero_mul_gvhd<-array(NA,dim=c(111,47))
rownames(nonzero_mul_gvhd)<-colnames(gvhd_reg_name[,2:112])
colnames(nonzero_mul_gvhd)<-colnames(gvhd_reg_name[,c(128:149,151:175)])
users_mul_gvhd<-array(NA,dim=c(111,47))
rownames(users_mul_gvhd)<-colnames(gvhd_reg_name[,2:112])
colnames(users_mul_gvhd)<-colnames(gvhd_reg_name[,c(128:149,151:175)])
nonusers_mul_gvhd<-array(NA,dim=c(111,47))
rownames(nonusers_mul_gvhd)<-colnames(gvhd_reg_name[,2:112])
colnames(nonusers_mul_gvhd)<-colnames(gvhd_reg_name[,c(128:149,151:175)])


estimate_mul_gvhd<-array(0,dim=c(111,47))
rownames(estimate_mul_gvhd)<-colnames(gvhd_reg_name[,2:112])
colnames(estimate_mul_gvhd)<-colnames(gvhd_reg_name[,c(128:149,151:175)])
Std_Error_mul_gvhd<-array(0,dim=c(111,47))
rownames(Std_Error_mul_gvhd)<-colnames(gvhd_reg_name[,2:112])
colnames(Std_Error_mul_gvhd)<-colnames(gvhd_reg_name[,c(128:149,151:175)])
t_value_mul_gvhd<-array(0,dim=c(111,47))
rownames(t_value_mul_gvhd)<-colnames(gvhd_reg_name[,2:112])
colnames(t_value_mul_gvhd)<-colnames(gvhd_reg_name[,c(128:149,151:175)])
p_value_mul_gvhd<-array(0,dim=c(111,47))
rownames(p_value_mul_gvhd)<-colnames(gvhd_reg_name[,2:112])
colnames(p_value_mul_gvhd)<-colnames(gvhd_reg_name[,c(128:149,151:175)])


for (x in 1:111) {
  for (y in 1:47){
    genus_name<-colnames(gvhd_reg_name[,2:112])
    drug_name<-colnames(gvhd_reg_name[,c(128:149,151:175)])
    
    temp_mul_gvhd<-select(gvhd_reg_name,genus_name[x],subject,flag,time,amplicon,sex,age,diagnosis,scs,DrH,Antibiotic,drug_name[y])
    temp_mul_gvhd$species<-temp_mul_gvhd[,1]
    temp_mul_gvhd_nonzero<- temp_mul_gvhd |> filter(species!=0) 
    nonzero_mul_gvhd[x,y] <-length(unique(temp_mul_gvhd_nonzero$subject))
    temp_mul_gvhd$drug<-temp_mul_gvhd[,(ncol(temp_mul_gvhd)-1)]
    temp_mul_gvhd_users<-temp_mul_gvhd[complete.cases(temp_mul_gvhd$drug),] |> filter(drug==1) 
    users_mul_gvhd[x,y] <- length(unique(temp_mul_gvhd_users$subject))
    temp_mul_gvhd_nonusers<-temp_mul_gvhd[complete.cases(temp_mul_gvhd$drug),] |> filter(drug==0) 
    nonusers_mul_gvhd[x,y] <- length(unique(temp_mul_gvhd_nonusers$subject))
    
    
    temp_mul_gvhd<-temp_mul_gvhd[complete.cases(temp_mul_gvhd),]
    drug<-temp_mul_gvhd[,ncol(temp_mul_gvhd)]
    Antibiotic<-temp_mul_gvhd[,"Antibiotic"]
    if(length(unique(drug))==2 & length(unique(Antibiotic))==2){
      if(length(unique(temp_mul_gvhd$time))!=1 & length(unique(temp_mul_gvhd$amplicon))!=1 & length(unique(temp_mul_gvhd$sex))!=1 & length(unique(temp_mul_gvhd$scs))!=1){
        temp_mul_gvhd1<-lm(temp_mul_gvhd[,1] ~ time + amplicon + sex + age + diagnosis + scs + DrH + Antibiotic +drug,temp_mul_gvhd)
        if (rownames(summary(temp_mul_gvhd1)$coefficients)[nrow(summary(temp_mul_gvhd1)$coefficients)] == "drug"){
          summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),1]->estimate_mul_gvhd[x,y]
          summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),2]->Std_Error_mul_gvhd[x,y]
          summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),3]->t_value_mul_gvhd[x,y]
          summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),4]->p_value_mul_gvhd[x,y]
        }else{
          estimate_mul_gvhd[x,y]<-NA
          Std_Error_mul_gvhd[x,y]<-NA
          t_value_mul_gvhd[x,y]<-NA
          p_value_mul_gvhd[x,y]<-NA
        }
      }else if(length(unique(temp_mul_gvhd$time))==1 & length(unique(temp_mul_gvhd$amplicon))!=1 & length(unique(temp_mul_gvhd$sex))!=1 & length(unique(temp_mul_gvhd$scs))!=1){
        temp_mul_gvhd1<-lm(temp_mul_gvhd[,1] ~ amplicon + sex + age + diagnosis + scs + DrH + Antibiotic+ drug,temp_mul_gvhd)
        if (rownames(summary(temp_mul_gvhd1)$coefficients)[nrow(summary(temp_mul_gvhd1)$coefficients)] == "drug"){
          estimate_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),1]
          Std_Error_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),2]
          t_value_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),3]
          p_value_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),4]
        }else{
          estimate_mul_gvhd[x,y]<-NA
          Std_Error_mul_gvhd[x,y]<-NA
          t_value_mul_gvhd[x,y]<-NA
          p_value_mul_gvhd[x,y]<-NA
        }
      }else if(length(unique(temp_mul_gvhd$time))!=1 & length(unique(temp_mul_gvhd$amplicon))==1 & length(unique(temp_mul_gvhd$sex))!=1 & length(unique(temp_mul_gvhd$scs))!=1){
        temp_mul_gvhd1<-lm(temp_mul_gvhd[,1] ~ time + sex + age + diagnosis + scs + DrH + Antibiotic+ drug,temp_mul_gvhd)
        if (rownames(summary(temp_mul_gvhd1)$coefficients)[nrow(summary(temp_mul_gvhd1)$coefficients)] == "drug"){
          estimate_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),1]
          Std_Error_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),2]
          t_value_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),3]
          p_value_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),4]
        }else{
          estimate_mul_gvhd[x,y]<-NA
          Std_Error_mul_gvhd[x,y]<-NA
          t_value_mul_gvhd[x,y]<-NA
          p_value_mul_gvhd[x,y]<-NA
        }
      }else if(length(unique(temp_mul_gvhd$time))!=1 & length(unique(temp_mul_gvhd$amplicon))!=1 & length(unique(temp_mul_gvhd$sex))==1 & length(unique(temp_mul_gvhd$scs))!=1){
        temp_mul_gvhd1<-lm(temp_mul_gvhd[,1] ~ time  + age + diagnosis + scs + DrH + Antibiotic+ drug,temp_mul_gvhd)
        if (rownames(summary(temp_mul_gvhd1)$coefficients)[nrow(summary(temp_mul_gvhd1)$coefficients)] == "drug"){
          estimate_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),1]
          Std_Error_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),2]
          t_value_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),3]
          p_value_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),4]
        }else{
          estimate_mul_gvhd[x,y]<-NA
          Std_Error_mul_gvhd[x,y]<-NA
          t_value_mul_gvhd[x,y]<-NA
          p_value_mul_gvhd[x,y]<-NA
        }
      }else if(length(unique(temp_mul_gvhd$time))!=1 & length(unique(temp_mul_gvhd$amplicon))!=1 & length(unique(temp_mul_gvhd$sex))!=1 & length(unique(temp_mul_gvhd$scs))==1){
        temp_mul_gvhd1<-lm(temp_mul_gvhd[,1] ~ time + sex + age + diagnosis + DrH+ Antibiotic + drug,temp_mul_gvhd)
        if (rownames(summary(temp_mul_gvhd1)$coefficients)[nrow(summary(temp_mul_gvhd1)$coefficients)] == "drug"){
          estimate_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),1]
          Std_Error_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),2]
          t_value_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),3]
          p_value_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),4]
        }else{
          estimate_mul_gvhd[x,y]<-NA
          Std_Error_mul_gvhd[x,y]<-NA
          t_value_mul_gvhd[x,y]<-NA
          p_value_mul_gvhd[x,y]<-NA
        }
      }else if(length(unique(temp_mul_gvhd$time))==1 & length(unique(temp_mul_gvhd$amplicon))==1 & length(unique(temp_mul_gvhd$sex))!=1 & length(unique(temp_mul_gvhd$scs))!=1){
        temp_mul_gvhd1<-lm(temp_mul_gvhd[,1] ~ sex + age + diagnosis + scs + DrH + Antibiotic+ drug,temp_mul_gvhd)
        if (rownames(summary(temp_mul_gvhd1)$coefficients)[nrow(summary(temp_mul_gvhd1)$coefficients)] == "drug"){
          estimate_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),1]
          Std_Error_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),2]
          t_value_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),3]
          p_value_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),4]
        }else{
          estimate_mul_gvhd[x,y]<-NA
          Std_Error_mul_gvhd<-NA
          t_value_mul_gvhd[x,y]<-NA
          p_value_mul_gvhd[x,y]<-NA
        }
      }else if(length(unique(temp_mul_gvhd$time))==1 & length(unique(temp_mul_gvhd$amplicon))!=1 & length(unique(temp_mul_gvhd$sex))==1 & length(unique(temp_mul_gvhd$scs))!=1){
        temp_mul_gvhd1<-lm(temp_mul_gvhd[,1] ~ amplicon + age + diagnosis + scs + DrH + Antibiotic+ drug,temp_mul_gvhd)
        if (rownames(summary(temp_mul_gvhd1)$coefficients)[nrow(summary(temp_mul_gvhd1)$coefficients)] == "drug"){
          estimate_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),1]
          Std_Error_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),2]
          t_value_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),3]
          p_value_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),4]
        }else{
          estimate_mul_gvhd[x,y]<-NA
          Std_Error_mul_gvhd[x,y]<-NA
          t_value_mul_gvhd[x,y]<-NA
          p_value_mul_gvhd[x,y]<-NA
        }
      }else if(length(unique(temp_mul_gvhd$time))==1 & length(unique(temp_mul_gvhd$amplicon))!=1 & length(unique(temp_mul_gvhd$sex))!=1 & length(unique(temp_mul_gvhd$scs))==1){
        temp_mul_gvhd1<-lm(temp_mul_gvhd[,1] ~ amplicon + sex + age + diagnosis + DrH + Antibiotic+ drug,temp_mul_gvhd)
        if (rownames(summary(temp_mul_gvhd1)$coefficients)[nrow(summary(temp_mul_gvhd1)$coefficients)] == "drug"){
          estimate_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),1]
          Std_Error_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),2]
          t_value_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),3]
          p_value_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),4]
        }else{
          estimate_mul_gvhd[x,y]<-NA
          Std_Error_mul_gvhd[x,y]<-NA
          t_value_mul_gvhd[x,y]<-NA
          p_value_mul_gvhd[x,y]<-NA
        }
      }else if(length(unique(temp_mul_gvhd$time))!=1 & length(unique(temp_mul_gvhd$amplicon))==1 & length(unique(temp_mul_gvhd$sex))==1 & length(unique(temp_mul_gvhd$scs))!=1){
        temp_mul_gvhd1<-lm(temp_mul_gvhd[,1] ~ time + age + diagnosis + scs + DrH + Antibiotic+ drug,temp_mul_gvhd)
        if (rownames(summary(temp_mul_gvhd1)$coefficients)[nrow(summary(temp_mul_gvhd1)$coefficients)] == "drug"){
          estimate_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),1]
          Std_Error_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),2]
          t_value_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),3]
          p_value_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),4]
        }else{
          estimate_mul_gvhd[x,y]<-NA
          Std_Error_mul_gvhd[x,y]<-NA
          t_value_mul_gvhd[x,y]<-NA
          p_value_mul_gvhd[x,y]<-NA
        }
      }else if(length(unique(temp_mul_gvhd$time))!=1 & length(unique(temp_mul_gvhd$amplicon))==1 & length(unique(temp_mul_gvhd$sex))!=1 & length(unique(temp_mul_gvhd$scs))==1){
        temp_mul_gvhd1<-lm(temp_mul_gvhd[,1] ~ time + sex + age + diagnosis + DrH + Antibiotic+ drug,temp_mul_gvhd)
        if (rownames(summary(temp_mul_gvhd1)$coefficients)[nrow(summary(temp_mul_gvhd1)$coefficients)] == "drug"){
          estimate_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),1]
          Std_Error_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),2]
          t_value_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),3]
          p_value_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),4]
        }else{
          estimate_mul_gvhd[x,y]<-NA
          Std_Error_mul_gvhd[x,y]<-NA
          t_value_mul_gvhd[x,y]<-NA
          p_value_mul_gvhd[x,y]<-NA
        }
      }else if(length(unique(temp_mul_gvhd$time))!=1 & length(unique(temp_mul_gvhd$amplicon))!=1 & length(unique(temp_mul_gvhd$sex))==1 & length(unique(temp_mul_gvhd$scs))==1){
        temp_mul_gvhd1<-lm(temp_mul_gvhd[,1] ~ amplicon + time + age + diagnosis + DrH + Antibiotic+ drug,temp_mul_gvhd)
        if (rownames(summary(temp_mul_gvhd1)$coefficients)[nrow(summary(temp_mul_gvhd1)$coefficients)] == "drug"){
          estimate_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),1]
          Std_Error_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),2]
          t_value_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),3]
          p_value_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),4]
        }else{
          estimate_mul_gvhd[x,y]<-NA
          Std_Error_mul_gvhd[x,y]<-NA
          t_value_mul_gvhd[x,y]<-NA
          p_value_mul_gvhd[x,y]<-NA
        }
      }else if(length(unique(temp_mul_gvhd$time))==1 & length(unique(temp_mul_gvhd$amplicon))==1 & length(unique(temp_mul_gvhd$sex))==1 & length(unique(temp_mul_gvhd$scs))!=1){
        temp_mul_gvhd1<-lm(temp_mul_gvhd[,1] ~ age + diagnosis + scs + DrH + Antibiotic+ drug,temp_mul_gvhd)
        if (rownames(summary(temp_mul_gvhd1)$coefficients)[nrow(summary(temp_mul_gvhd1)$coefficients)] == "drug"){
          estimate_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),1]
          Std_Error_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),2]
          t_value_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),3]
          p_value_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),4]
        }else{
          estimate_mul_gvhd[x,y]<-NA
          Std_Error_mul_gvhd[x,y]<-NA
          t_value_mul_gvhd[x,y]<-NA
          p_value_mul_gvhd[x,y]<-NA
        }
      }else if(length(unique(temp_mul_gvhd$time))==1 & length(unique(temp_mul_gvhd$amplicon))==1 & length(unique(temp_mul_gvhd$sex))!=1 & length(unique(temp_mul_gvhd$scs))==1){
        temp_mul_gvhd1<-lm(temp_mul_gvhd[,1] ~ sex + age + diagnosis + DrH + Antibiotic+ drug,temp_mul_gvhd)
        if (rownames(summary(temp_mul_gvhd1)$coefficients)[nrow(summary(temp_mul_gvhd1)$coefficients)] == "drug"){
          estimate_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),1]
          Std_Error_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),2]
          t_value_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),3]
          p_value_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),4]
        }else{
          estimate_mul_gvhd[x,y]<-NA
          Std_Error_mul_gvhd[x,y]<-NA
          t_value_mul_gvhd[x,y]<-NA
          p_value_mul_gvhd[x,y]<-NA
        }
      }else if(length(unique(temp_mul_gvhd$time))==1 & length(unique(temp_mul_gvhd$amplicon))!=1 & length(unique(temp_mul_gvhd$sex))==1 & length(unique(temp_mul_gvhd$scs))==1){
        temp_mul_gvhd1<-lm(temp_mul_gvhd[,1] ~ amplicon + age + diagnosis + DrH + Antibiotic+ drug,temp_mul_gvhd)
        if (rownames(summary(temp_mul_gvhd1)$coefficients)[nrow(summary(temp_mul_gvhd1)$coefficients)] == "drug"){
          estimate_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),1]
          Std_Error_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),2]
          t_value_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),3]
          p_value_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),4]
        }else{
          estimate_mul_gvhd[x,y]<-NA
          Std_Error_mul_gvhd[x,y]<-NA
          t_value_mul_gvhd[x,y]<-NA
          p_value_mul_gvhd[x,y]<-NA
        }
      }else if(length(unique(temp_mul_gvhd$time))!=1 & length(unique(temp_mul_gvhd$amplicon))==1 & length(unique(temp_mul_gvhd$sex))==1 & length(unique(temp_mul_gvhd$scs))==1){
        temp_mul_gvhd1<-lm(temp_mul_gvhd[,1] ~ time + age + diagnosis + DrH + Antibiotic+ drug,temp_mul_gvhd)
        if (rownames(summary(temp_mul_gvhd1)$coefficients)[nrow(summary(temp_mul_gvhd1)$coefficients)] == "drug"){
          estimate_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),1]
          Std_Error_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),2]
          t_value_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),3]
          p_value_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),4]
        }else{
          estimate_mul_gvhd[x,y]<-NA
          Std_Error_mul_gvhd[x,y]<-NA
          t_value_mul_gvhd[x,y]<-NA
          p_value_mul_gvhd[x,y]<-NA
        }
      }else{
        temp_mul_gvhd1<-lm(temp_mul_gvhd[,1] ~ age + diagnosis + DrH + Antibiotic+ drug,temp_mul_gvhd)
        if (rownames(summary(temp_mul_gvhd1)$coefficients)[nrow(summary(temp_mul_gvhd1)$coefficients)] == "drug"){
          estimate_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),1]
          Std_Error_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),2]
          t_value_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),3]
          p_value_mul_gvhd[x,y]<-summary(temp_mul_gvhd1)$coefficients[nrow(summary(temp_mul_gvhd1)$coefficients),4]
        }else{
          estimate_mul_gvhd[x,y]<-NA
          Std_Error_mul_gvhd[x,y]<-NA
          t_value_mul_gvhd[x,y]<-NA
          p_value_mul_gvhd[x,y]<-NA
        }
      }
    }else{
      estimate_mul_gvhd[x,y]<-NA
      Std_Error_mul_gvhd[x,y]<-NA
      t_value_mul_gvhd[x,y]<-NA
      p_value_mul_gvhd[x,y]<-NA
    }
  }
}

number_mul_gvhd<-as.data.frame(number_mul_gvhd)
nonzero_mul_gvhd<-as.data.frame(nonzero_mul_gvhd)
users_mul_gvhd<-as.data.frame(users_mul_gvhd)
nonusers_mul_gvhd<-as.data.frame(nonusers_mul_gvhd)

estimate_mul_gvhd<-as.data.frame(estimate_mul_gvhd)
p_value_mul_gvhd<-as.data.frame(p_value_mul_gvhd)
Std_Error_mul_gvhd<-as.data.frame(Std_Error_mul_gvhd)
t_value_mul_gvhd<-as.data.frame(t_value_mul_gvhd)


# chemo ------------------------------------------------------------------


group_chemo_level7<-read.csv("F:/数据/document/学校/硕士/MMF/meta/HCT gvhd/srr/level-7-chemo.csv",header=T)

chemo<-group_chemo_level7[which(group_chemo_level7$disease == "D"),]

chemo_0<-chemo
chemo_0[is.na(chemo_0)]<-0

sum_abundance_chemo<-apply(chemo_0[,2:850],1,sum)

chemo_rel<-array(0,dim=c(1097,849))
for (j in 1:849) {
  for (i in 1:1097){
    chemo_rel[i,j]<-chemo_0[i,j+1]/sum_abundance_chemo[i]*100
  }
}
chemo_rel <- as.data.frame(chemo_rel)
chemo_rel[is.na(chemo_rel)]<-0

sample_chemo<-chemo_0$index
species_chemo<-colnames(chemo_0)
colnames(chemo_rel)<-species_chemo[2:850]
rownames(chemo_rel)<-sample_chemo

chemo_rel_fil_10 <- filtering_taxonomy(chemo_rel,0.00000001,10)

summary_taxnomy_filtering_10_chemo<- read.table("summary_taxonomy_filtering.txt",sep="\t",header=T)
filtered_taxonomy_10_chemo<-read.table("filtered_taxonomy.txt",sep="\t",header=T)

filtered_taxonomy_10_1_chemo<-cbind(filtered_taxonomy_10_chemo,index=rownames(filtered_taxonomy_10_chemo))

metadata_chemo<-as.data.frame(cbind(index=sample_chemo,chemo[,851:915]))
chemo_reg<-merge(filtered_taxonomy_10_1_chemo,metadata_chemo,by="index")
rownames(chemo_reg)<-chemo_reg$index
chemo_reg_name<-chemo_reg
colnames(chemo_reg_name)[109]<-"amplicon"
colnames(chemo_reg_name)[110]<-"sex"
colnames(chemo_reg_name)[113]<-"scs"
colnames(chemo_reg_name)[114]<-"DrH"



## single drug analysis ----------------------------------------------------


number_chemo<-array(length(unique(chemo_reg_name$subject)),dim=c(98,48))
rownames(number_chemo)<-colnames(chemo_reg_name[,2:99])
colnames(number_chemo)<-colnames(chemo_reg_name[,115:162])
nonzero_chemo<-array(NA,dim=c(98,48))
rownames(nonzero_chemo)<-colnames(chemo_reg_name[,2:99])
colnames(nonzero_chemo)<-colnames(chemo_reg_name[,115:162])
users_chemo<-array(NA,dim=c(98,48))
rownames(users_chemo)<-colnames(chemo_reg_name[,2:99])
colnames(users_chemo)<-colnames(chemo_reg_name[,115:162])
nonusers_chemo<-array(NA,dim=c(98,48))
rownames(nonusers_chemo)<-colnames(chemo_reg_name[,2:99])
colnames(nonusers_chemo)<-colnames(chemo_reg_name[,115:162])



estimate_chemo<-array(0,dim=c(98,48))
rownames(estimate_chemo)<-colnames(chemo_reg_name[,2:99])
colnames(estimate_chemo)<-colnames(chemo_reg_name[,115:162])
Std_Error_chemo<-array(0,dim=c(98,48))
rownames(Std_Error_chemo)<-colnames(chemo_reg_name[,2:99])
colnames(Std_Error_chemo)<-colnames(chemo_reg_name[,115:162])
t_value_chemo<-array(0,dim=c(98,48))
rownames(t_value_chemo)<-colnames(chemo_reg_name[,2:99])
colnames(t_value_chemo)<-colnames(chemo_reg_name[,115:162])
p_value_chemo<-array(0,dim=c(98,48))
rownames(p_value_chemo)<-colnames(chemo_reg_name[,2:99])
colnames(p_value_chemo)<-colnames(chemo_reg_name[,115:162])


for (x in 1:98) {
  for (y in 1:48){
    genus_name<-colnames(chemo_reg_name[,2:99])
    drug_name<-colnames(chemo_reg_name[,115:162])
    
    temp_chemo<-select(chemo_reg_name,genus_name[x],subject,flag,time,amplicon,sex,age,diagnosis,drug_name[y])
    
    temp_chemo$species<-temp_chemo[,1]
    temp_chemo_nonzero<- temp_chemo |> filter(species!=0) 
    nonzero_chemo[x,y] <-length(unique(temp_chemo_nonzero$subject))
    temp_chemo$drug<-temp_chemo[,(ncol(temp_chemo)-1)]
    temp_chemo_users<-temp_chemo[complete.cases(temp_chemo$drug),] |> filter(drug==1) 
    users_chemo[x,y] <- length(unique(temp_chemo_users$subject))
    temp_chemo_nonusers<-temp_chemo[complete.cases(temp_chemo$drug),] |> filter(drug==0) 
    nonusers_chemo[x,y] <- length(unique(temp_nonusers$subject))
    
    
    temp_chemo<-temp_chemo[complete.cases(temp_chemo),]
    drug_chemo<-temp_chemo[,ncol(temp_chemo)]
    if(length(unique(drug_chemo))==2){
      if(length(unique(temp_chemo$time))!=1 & length(unique(temp_chemo$amplicon))!=1 & length(unique(temp_chemo$sex))!=1 & length(unique(temp_chemo$diagnosis))!=1){
        temp_chemo1<-lm(temp_chemo[,1] ~ time + amplicon + sex + age + diagnosis + drug_chemo,temp_chemo)
        if (rownames(summary(temp_chemo1)$coefficients)[nrow(summary(temp_chemo1)$coefficients)] == "drug_chemo"){
          summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),1]->estimate_chemo[x,y]
          summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),2]->Std_Error_chemo[x,y]
          summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),3]->t_value_chemo[x,y]
          summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),4]->p_value_chemo[x,y]
        }else{
          estimate_chemo[x,y]<-NA
          Std_Error_chemo[x,y]<-NA
          t_value_chemo[x,y]<-NA
          p_value_chemo[x,y]<-NA
        }
      }else if(length(unique(temp_chemo$time))==1 & length(unique(temp_chemo$amplicon))!=1 & length(unique(temp_chemo$sex))!=1 & length(unique(temp_chemo$diagnosis))!=1){
        temp_chemo1<-lm(temp_chemo[,1] ~ amplicon + sex + age + diagnosis + drug_chemo,temp_chemo)
        if (rownames(summary(temp_chemo1)$coefficients)[nrow(summary(temp_chemo1)$coefficients)] == "drug_chemo"){
          estimate_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),1]
          Std_Error_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),2]
          t_value_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),3]
          p_value_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),4]
        }else{
          estimate_chemo[x,y]<-NA
          Std_Error_chemo[x,y]<-NA
          t_value_chemo[x,y]<-NA
          p_value_chemo[x,y]<-NA
        }
      }else if(length(unique(temp_chemo$time))!=1 & length(unique(temp_chemo$amplicon))==1 & length(unique(temp_chemo$sex))!=1 & length(unique(temp_chemo$diagnosis))!=1){
        temp_chemo1<-lm(temp_chemo[,1] ~ time + sex + age + diagnosis + drug_chemo,temp_chemo)
        if (rownames(summary(temp_chemo1)$coefficients)[nrow(summary(temp_chemo1)$coefficients)] == "drug_chemo"){
          estimate_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),1]
          Std_Error_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),2]
          t_value_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),3]
          p_value_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),4]
        }else{
          estimate_chemo[x,y]<-NA
          Std_Error_chemo[x,y]<-NA
          t_value_chemo[x,y]<-NA
          p_value_chemo[x,y]<-NA
        }
      }else if(length(unique(temp_chemo$time))!=1 & length(unique(temp_chemo$amplicon))!=1 & length(unique(temp_chemo$sex))==1 & length(unique(temp_chemo$diagnosis))!=1){
        temp_chemo1<-lm(temp_chemo[,1] ~ time  + age + diagnosis + drug_chemo,temp_chemo)
        if (rownames(summary(temp_chemo1)$coefficients)[nrow(summary(temp_chemo1)$coefficients)] == "drug_chemo"){
          estimate_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),1]
          Std_Error_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),2]
          t_value_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),3]
          p_value_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),4]
        }else{
          estimate_chemo[x,y]<-NA
          Std_Error_chemo[x,y]<-NA
          t_value_chemo[x,y]<-NA
          p_value_chemo[x,y]<-NA
        }
      }else if(length(unique(temp_chemo$time))!=1 & length(unique(temp_chemo$amplicon))!=1 & length(unique(temp_chemo$sex))!=1 & length(unique(temp_chemo$diagnosis))==1){
        temp_chemo1<-lm(temp_chemo[,1] ~ time + sex + age + drug_chemo,temp_chemo)
        if (rownames(summary(temp_chemo1)$coefficients)[nrow(summary(temp_chemo1)$coefficients)] == "drug_chemo"){
          estimate_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),1]
          Std_Error_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),2]
          t_value_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),3]
          p_value_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),4]
        }else{
          estimate_chemo[x,y]<-NA
          Std_Error_chemo[x,y]<-NA
          t_value_chemo[x,y]<-NA
          p_value_chemo[x,y]<-NA
        }
      }else if(length(unique(temp_chemo$time))==1 & length(unique(temp_chemo$amplicon))==1 & length(unique(temp_chemo$sex))!=1 & length(unique(temp_chemo$diagnosis))!=1){
        temp_chemo1<-lm(temp_chemo[,1] ~ sex + age + diagnosis + drug_chemo,temp_chemo)
        if (rownames(summary(temp_chemo1)$coefficients)[nrow(summary(temp_chemo1)$coefficients)] == "drug_chemo"){
          estimate_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),1]
          Std_Error_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),2]
          t_value_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),3]
          p_value_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),4]
        }else{
          estimate_chemo[x,y]<-NA
          Std_Error_chemo[x,y]<-NA
          t_value_chemo[x,y]<-NA
          p_value_chemo[x,y]<-NA
        }
      }else if(length(unique(temp_chemo$time))==1 & length(unique(temp_chemo$amplicon))!=1 & length(unique(temp_chemo$sex))==1 & length(unique(temp_chemo$diagnosis))!=1){
        temp_chemo1<-lm(temp_chemo[,1] ~ amplicon + age + diagnosis + drug_chemo,temp_chemo)
        if (rownames(summary(temp_chemo1)$coefficients)[nrow(summary(temp_chemo1)$coefficients)] == "drug_chemo"){
          estimate_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),1]
          Std_Error_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),2]
          t_value_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),3]
          p_value_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),4]
        }else{
          estimate_chemo[x,y]<-NA
          Std_Error_chemo[x,y]<-NA
          t_value_chemo[x,y]<-NA
          p_value_chemo[x,y]<-NA
        }
      }else if(length(unique(temp_chemo$time))==1 & length(unique(temp_chemo$amplicon))!=1 & length(unique(temp_chemo$sex))!=1 & length(unique(temp_chemo$diagnosis))==1){
        temp_chemo1<-lm(temp_chemo[,1] ~ amplicon + sex + age + drug_chemo,temp_chemo)
        if (rownames(summary(temp_chemo1)$coefficients)[nrow(summary(temp_chemo1)$coefficients)] == "drug_chemo"){
          estimate_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),1]
          Std_Error_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),2]
          t_value_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),3]
          p_value_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),4]
        }else{
          estimate_chemo[x,y]<-NA
          Std_Error_chemo[x,y]<-NA
          t_value_chemo[x,y]<-NA
          p_value_chemo[x,y]<-NA
        }
      }else if(length(unique(temp_chemo$time))!=1 & length(unique(temp_chemo$amplicon))==1 & length(unique(temp_chemo$sex))==1 & length(unique(temp_chemo$diagnosis))!=1){
        temp_chemo1<-lm(temp_chemo[,1] ~ time + age + diagnosis + drug_chemo,temp_chemo)
        if (rownames(summary(temp_chemo1)$coefficients)[nrow(summary(temp_chemo1)$coefficients)] == "drug_chemo"){
          estimate_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),1]
          Std_Error_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),2]
          t_value_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),3]
          p_value_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),4]
        }else{
          estimate_chemo[x,y]<-NA
          Std_Error_chemo[x,y]<-NA
          t_value_chemo[x,y]<-NA
          p_value_chemo[x,y]<-NA
        }
      }else if(length(unique(temp_chemo$time))!=1 & length(unique(temp_chemo$amplicon))==1 & length(unique(temp_chemo$sex))!=1 & length(unique(temp_chemo$diagnosis))==1){
        temp_chemo1<-lm(temp_chemo[,1] ~ time + sex + age + drug_chemo,temp_chemo)
        if (rownames(summary(temp_chemo1)$coefficients)[nrow(summary(temp_chemo1)$coefficients)] == "drug_chemo"){
          estimate_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),1]
          Std_Error_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),2]
          t_value_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),3]
          p_value_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),4]
        }else{
          estimate_chemo[x,y]<-NA
          Std_Error_chemo[x,y]<-NA
          t_value_chemo[x,y]<-NA
          p_value_chemo[x,y]<-NA
        }
      }else if(length(unique(temp_chemo$time))!=1 & length(unique(temp_chemo$amplicon))!=1 & length(unique(temp_chemo$sex))==1 & length(unique(temp_chemo$diagnosis))==1){
        temp_chemo1<-lm(temp_chemo[,1] ~ amplicon + time + age + drug_chemo,temp_chemo)
        if (rownames(summary(temp_chemo1)$coefficients)[nrow(summary(temp_chemo1)$coefficients)] == "drug_chemo"){
          estimate_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),1]
          Std_Error_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),2]
          t_value_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),3]
          p_value_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),4]
        }else{
          estimate_chemo[x,y]<-NA
          Std_Error_chemo[x,y]<-NA
          t_value_chemo[x,y]<-NA
          p_value_chemo[x,y]<-NA
        }
      }else if(length(unique(temp_chemo$time))==1 & length(unique(temp_chemo$amplicon))==1 & length(unique(temp_chemo$sex))==1 & length(unique(temp_chemo$diagnosis))!=1){
        temp_chemo1<-lm(temp_chemo[,1] ~ age + diagnosis + drug_chemo,temp_chemo)
        if (rownames(summary(temp_chemo1)$coefficients)[nrow(summary(temp_chemo1)$coefficients)] == "drug_chemo"){
          estimate_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),1]
          Std_Error_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),2]
          t_value_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),3]
          p_value_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),4]
        }else{
          estimate_chemo[x,y]<-NA
          Std_Error_chemo[x,y]<-NA
          t_value_chemo[x,y]<-NA
          p_value_chemo[x,y]<-NA
        }
      }else if(length(unique(temp_chemo$time))==1 & length(unique(temp_chemo$amplicon))==1 & length(unique(temp_chemo$sex))!=1 & length(unique(temp_chemo$diagnosis))==1){
        temp_chemo1<-lm(temp_chemo[,1] ~ sex + age + drug_chemo,temp_chemo)
        if (rownames(summary(temp_chemo1)$coefficients)[nrow(summary(temp_chemo1)$coefficients)] == "drug_chemo"){
          estimate_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),1]
          Std_Error_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),2]
          t_value_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),3]
          p_value_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),4]
        }else{
          estimate_chemo[x,y]<-NA
          Std_Error_chemo[x,y]<-NA
          t_value_chemo[x,y]<-NA
          p_value_chemo[x,y]<-NA
        }
      }else if(length(unique(temp_chemo$time))==1 & length(unique(temp_chemo$amplicon))!=1 & length(unique(temp_chemo$sex))==1 & length(unique(temp_chemo$diagnosis))==1){
        temp_chemo1<-lm(temp_chemo[,1] ~ amplicon + age + drug_chemo,temp_chemo)
        if (rownames(summary(temp_chemo1)$coefficients)[nrow(summary(temp_chemo1)$coefficients)] == "drug_chemo"){
          estimate_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),1]
          Std_Error_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),2]
          t_value_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),3]
          p_value_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),4]
        }else{
          estimate_chemo[x,y]<-NA
          Std_Error_chemo[x,y]<-NA
          t_value_chemo[x,y]<-NA
          p_value_chemo[x,y]<-NA
        }
      }else if(length(unique(temp_chemo$time))!=1 & length(unique(temp_chemo$amplicon))==1 & length(unique(temp_chemo$sex))==1 & length(unique(temp_chemo$diagnosis))==1){
        temp_chemo1<-lm(temp_chemo[,1] ~ time + age + drug_chemo,temp_chemo)
        if (rownames(summary(temp_chemo1)$coefficients)[nrow(summary(temp_chemo1)$coefficients)] == "drug_chemo"){
          estimate_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),1]
          Std_Error_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),2]
          t_value_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),3]
          p_value_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),4]
        }else{
          estimate_chemo[x,y]<-NA
          Std_Error_chemo[x,y]<-NA
          t_value_chemo[x,y]<-NA
          p_value_chemo[x,y]<-NA
        }
      }else{
        temp_chemo1<-lm(temp_chemo[,1] ~ age + drug_chemo,temp_chemo)
        if (rownames(summary(temp_chemo1)$coefficients)[nrow(summary(temp_chemo1)$coefficients)] == "drug_chemo"){
          estimate_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),1]
          Std_Error_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),2]
          t_value_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),3]
          p_value_chemo[x,y]<-summary(temp_chemo1)$coefficients[nrow(summary(temp_chemo1)$coefficients),4]
        }else{
          estimate_chemo[x,y]<-NA
          Std_Error_chemo[x,y]<-NA
          t_value_chemo[x,y]<-NA
          p_value_chemo[x,y]<-NA
        }
      }
    }else{
      estimate_chemo[x,y]<-NA
      Std_Error_chemo[x,y]<-NA
      t_value_chemo[x,y]<-NA
      p_value_chemo[x,y]<-NA
    }
  }
}

number_chemo<-as.data.frame(number_chemo)
nonzero_chemo<-as.data.frame(nonzero_chemo)
users_chemo<-as.data.frame(users_chemo)
nonusers_chemo<-as.data.frame(nonusers_chemo)

estimate_chemo<-as.data.frame(estimate_chemo)
p_value_chemo<-as.data.frame(p_value_chemo)
Std_Error_chemo<-as.data.frame(Std_Error_chemo)
t_value_chemo<-as.data.frame(t_value_chemo)



## multiple drug analysis --------------------------------------------------

number_mul_chemo_drug<-array(length(unique(chemo_reg_name$subject)),dim=c(98,22))
rownames(number_mul_chemo_drug)<-colnames(chemo_reg_name[,2:99])
colnames(number_mul_chemo_drug)<-colnames(chemo_reg_name[,115:136])
nonzero_mul_chemo_drug<-array(NA,dim=c(98,22))
rownames(nonzero_mul_chemo_drug)<-colnames(chemo_reg_name[,2:99])
colnames(nonzero_mul_chemo_drug)<-colnames(chemo_reg_name[,115:136])
users_mul_chemo_drug<-array(NA,dim=c(98,22))
rownames(users_mul_chemo_drug)<-colnames(chemo_reg_name[,2:99])
colnames(users_mul_chemo_drug)<-colnames(chemo_reg_name[,115:136])
nonusers_mul_chemo_drug<-array(NA,dim=c(98,22))
rownames(nonusers_mul_chemo_drug)<-colnames(chemo_reg_name[,2:99])
colnames(nonusers_mul_chemo_drug)<-colnames(chemo_reg_name[,115:136])


estimate_mul_chemo_drug<-array(0,dim=c(98,22))
rownames(estimate_mul_chemo_drug)<-colnames(chemo_reg_name[,2:99])
colnames(estimate_mul_chemo_drug)<-colnames(chemo_reg_name[,115:136])
Std_Error_mul_chemo_drug<-array(0,dim=c(98,22))
rownames(Std_Error_mul_chemo_drug)<-colnames(chemo_reg_name[,2:99])
colnames(Std_Error_mul_chemo_drug)<-colnames(chemo_reg_name[,115:136])
t_value_mul_chemo_drug<-array(0,dim=c(98,22))
rownames(t_value_mul_chemo_drug)<-colnames(chemo_reg_name[,2:99])
colnames(t_value_mul_chemo_drug)<-colnames(chemo_reg_name[,115:136])
p_value_mul_chemo_drug<-array(0,dim=c(98,22))
rownames(p_value_mul_chemo_drug)<-colnames(chemo_reg_name[,2:99])
colnames(p_value_mul_chemo_drug)<-colnames(chemo_reg_name[,115:136])


for (x in 1:98) {
  for (y in 1:22){
    genus_name<-colnames(chemo_reg_name[,2:99])
    drug_name<-colnames(chemo_reg_name[,115:136])
    
    temp_mul_chemo_drug<-select(chemo_reg_name,genus_name[x],subject,flag,time,amplicon,sex,diagnosis,Trimethoprim.sulfamethoxazole,drug_name[y])
    temp_mul_chemo_drug$species<-temp_mul_chemo_drug[,1]
    temp_mul_chemo_nonzero_drug<- temp_mul_chemo_drug |> filter(species!=0) 
    nonzero_mul_chemo_drug[x,y] <-length(unique(temp_mul_chemo_nonzero_drug$subject))
    temp_mul_chemo_drug$drug<-temp_mul_chemo_drug[,(ncol(temp_mul_chemo_drug)-1)]
    temp_mul_chemo_users_drug<-temp_mul_chemo_drug[complete.cases(temp_mul_chemo_drug$drug),] |> filter(drug==1) 
    users_mul_chemo_drug[x,y] <- length(unique(temp_mul_chemo_users_drug$subject))
    temp_mul_chemo_nonusers_drug<-temp_mul_chemo_drug[complete.cases(temp_mul_chemo_drug$drug),] |> filter(drug==0) 
    nonusers_mul_chemo_drug[x,y] <- length(unique(temp_mul_chemo_nonusers_drug$subject))
    
    
    temp_mul_chemo_drug<-temp_mul_chemo_drug[complete.cases(temp_mul_chemo_drug),]
    drug<-temp_mul_chemo_drug[,ncol(temp_mul_chemo_drug)]
    Antibiotic<-temp_mul_chemo_drug[,"Trimethoprim.sulfamethoxazole"]
    if(length(unique(drug))==2 & length(unique(Antibiotic))==2){
      if(length(unique(temp_mul_chemo_drug$time))!=1 & length(unique(temp_mul_chemo_drug$amplicon))!=1 & length(unique(temp_mul_chemo_drug$sex))!=1 & length(unique(temp_mul_chemo_drug$diagnosis))!=1){
        temp_mul_chemo_drug1<-lm(temp_mul_chemo_drug[,1] ~ time + amplicon + sex +  diagnosis + Antibiotic +drug,temp_mul_chemo_drug)
        if (rownames(summary(temp_mul_chemo_drug1)$coefficients)[nrow(summary(temp_mul_chemo_drug1)$coefficients)] == "drug"){
          summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),1]->estimate_mul_chemo_drug[x,y]
          summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),2]->Std_Error_mul_chemo_drug[x,y]
          summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),3]->t_value_mul_chemo_drug[x,y]
          summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),4]->p_value_mul_chemo_drug[x,y]
        }else{
          estimate_mul_chemo_drug[x,y]<-NA
          Std_Error_mul_chemo_drug[x,y]<-NA
          t_value_mul_chemo_drug[x,y]<-NA
          p_value_mul_chemo_drug[x,y]<-NA
        }
      }else if(length(unique(temp_mul_chemo_drug$time))==1 & length(unique(temp_mul_chemo_drug$amplicon))!=1 & length(unique(temp_mul_chemo_drug$sex))!=1 & length(unique(temp_mul_chemo_drug$diagnosis))!=1){
        temp_mul_chemo_drug1<-lm(temp_mul_chemo_drug[,1] ~ amplicon + sex +  diagnosis + Antibiotic+ drug,temp_mul_chemo_drug)
        if (rownames(summary(temp_mul_chemo_drug1)$coefficients)[nrow(summary(temp_mul_chemo_drug1)$coefficients)] == "drug"){
          estimate_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),1]
          Std_Error_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),2]
          t_value_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),3]
          p_value_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),4]
        }else{
          estimate_mul_chemo_drug[x,y]<-NA
          Std_Error_mul_chemo_drug[x,y]<-NA
          t_value_mul_chemo_drug[x,y]<-NA
          p_value_mul_chemo_drug[x,y]<-NA
        }
      }else if(length(unique(temp_mul_chemo_drug$time))!=1 & length(unique(temp_mul_chemo_drug$amplicon))==1 & length(unique(temp_mul_chemo_drug$sex))!=1 & length(unique(temp_mul_chemo_drug$diagnosis))!=1){
        temp_mul_chemo_drug1<-lm(temp_mul_chemo_drug[,1] ~ time + sex +  diagnosis + Antibiotic+ drug,temp_mul_chemo_drug)
        if (rownames(summary(temp_mul_chemo_drug1)$coefficients)[nrow(summary(temp_mul_chemo_drug1)$coefficients)] == "drug"){
          estimate_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),1]
          Std_Error_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),2]
          t_value_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),3]
          p_value_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),4]
        }else{
          estimate_mul_chemo_drug[x,y]<-NA
          Std_Error_mul_chemo_drug[x,y]<-NA
          t_value_mul_chemo_drug[x,y]<-NA
          p_value_mul_chemo_drug[x,y]<-NA
        }
      }else if(length(unique(temp_mul_chemo_drug$time))!=1 & length(unique(temp_mul_chemo_drug$amplicon))!=1 & length(unique(temp_mul_chemo_drug$sex))==1 & length(unique(temp_mul_chemo_drug$diagnosis))!=1){
        temp_mul_chemo_drug1<-lm(temp_mul_chemo_drug[,1] ~ time  +  diagnosis + Antibiotic+ drug,temp_mul_chemo_drug)
        if (rownames(summary(temp_mul_chemo_drug1)$coefficients)[nrow(summary(temp_mul_chemo_drug1)$coefficients)] == "drug"){
          estimate_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),1]
          Std_Error_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),2]
          t_value_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),3]
          p_value_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),4]
        }else{
          estimate_mul_chemo_drug[x,y]<-NA
          Std_Error_mul_chemo_drug[x,y]<-NA
          t_value_mul_chemo_drug[x,y]<-NA
          p_value_mul_chemo_drug[x,y]<-NA
        }
      }else if(length(unique(temp_mul_chemo_drug$time))!=1 & length(unique(temp_mul_chemo_drug$amplicon))!=1 & length(unique(temp_mul_chemo_drug$sex))!=1 & length(unique(temp_mul_chemo_drug$diagnosis))==1){
        temp_mul_chemo_drug1<-lm(temp_mul_chemo_drug[,1] ~ time + sex +  Antibiotic + drug,temp_mul_chemo_drug)
        if (rownames(summary(temp_mul_chemo_drug1)$coefficients)[nrow(summary(temp_mul_chemo_drug1)$coefficients)] == "drug"){
          estimate_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),1]
          Std_Error_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),2]
          t_value_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),3]
          p_value_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),4]
        }else{
          estimate_mul_chemo_drug[x,y]<-NA
          Std_Error_mul_chemo_drug[x,y]<-NA
          t_value_mul_chemo_drug[x,y]<-NA
          p_value_mul_chemo_drug[x,y]<-NA
        }
      }else if(length(unique(temp_mul_chemo_drug$time))==1 & length(unique(temp_mul_chemo_drug$amplicon))==1 & length(unique(temp_mul_chemo_drug$sex))!=1 & length(unique(temp_mul_chemo_drug$diagnosis))!=1){
        temp_mul_chemo_drug1<-lm(temp_mul_chemo_drug[,1] ~ sex +  diagnosis + Antibiotic+ drug,temp_mul_chemo_drug)
        if (rownames(summary(temp_mul_chemo_drug1)$coefficients)[nrow(summary(temp_mul_chemo_drug1)$coefficients)] == "drug"){
          estimate_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),1]
          Std_Error_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),2]
          t_value_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),3]
          p_value_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),4]
        }else{
          estimate_mul_chemo_drug[x,y]<-NA
          Std_Error_mul_chemo_drug<-NA
          t_value_mul_chemo_drug[x,y]<-NA
          p_value_mul_chemo_drug[x,y]<-NA
        }
      }else if(length(unique(temp_mul_chemo_drug$time))==1 & length(unique(temp_mul_chemo_drug$amplicon))!=1 & length(unique(temp_mul_chemo_drug$sex))==1 & length(unique(temp_mul_chemo_drug$diagnosis))!=1){
        temp_mul_chemo_drug1<-lm(temp_mul_chemo_drug[,1] ~ amplicon +  diagnosis + Antibiotic+ drug,temp_mul_chemo_drug)
        if (rownames(summary(temp_mul_chemo_drug1)$coefficients)[nrow(summary(temp_mul_chemo_drug1)$coefficients)] == "drug"){
          estimate_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),1]
          Std_Error_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),2]
          t_value_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),3]
          p_value_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),4]
        }else{
          estimate_mul_chemo_drug[x,y]<-NA
          Std_Error_mul_chemo_drug[x,y]<-NA
          t_value_mul_chemo_drug[x,y]<-NA
          p_value_mul_chemo_drug[x,y]<-NA
        }
      }else if(length(unique(temp_mul_chemo_drug$time))==1 & length(unique(temp_mul_chemo_drug$amplicon))!=1 & length(unique(temp_mul_chemo_drug$sex))!=1 & length(unique(temp_mul_chemo_drug$diagnosis))==1){
        temp_mul_chemo_drug1<-lm(temp_mul_chemo_drug[,1] ~ amplicon + sex +  Antibiotic+ drug,temp_mul_chemo_drug)
        if (rownames(summary(temp_mul_chemo_drug1)$coefficients)[nrow(summary(temp_mul_chemo_drug1)$coefficients)] == "drug"){
          estimate_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),1]
          Std_Error_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),2]
          t_value_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),3]
          p_value_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),4]
        }else{
          estimate_mul_chemo_drug[x,y]<-NA
          Std_Error_mul_chemo_drug[x,y]<-NA
          t_value_mul_chemo_drug[x,y]<-NA
          p_value_mul_chemo_drug[x,y]<-NA
        }
      }else if(length(unique(temp_mul_chemo_drug$time))!=1 & length(unique(temp_mul_chemo_drug$amplicon))==1 & length(unique(temp_mul_chemo_drug$sex))==1 & length(unique(temp_mul_chemo_drug$diagnosis))!=1){
        temp_mul_chemo_drug1<-lm(temp_mul_chemo_drug[,1] ~ time +  diagnosis + Antibiotic+ drug,temp_mul_chemo_drug)
        if (rownames(summary(temp_mul_chemo_drug1)$coefficients)[nrow(summary(temp_mul_chemo_drug1)$coefficients)] == "drug"){
          estimate_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),1]
          Std_Error_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),2]
          t_value_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),3]
          p_value_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),4]
        }else{
          estimate_mul_chemo_drug[x,y]<-NA
          Std_Error_mul_chemo_drug[x,y]<-NA
          t_value_mul_chemo_drug[x,y]<-NA
          p_value_mul_chemo_drug[x,y]<-NA
        }
      }else if(length(unique(temp_mul_chemo_drug$time))!=1 & length(unique(temp_mul_chemo_drug$amplicon))==1 & length(unique(temp_mul_chemo_drug$sex))!=1 & length(unique(temp_mul_chemo_drug$diagnosis))==1){
        temp_mul_chemo_drug1<-lm(temp_mul_chemo_drug[,1] ~ time + sex +  Antibiotic+ drug,temp_mul_chemo_drug)
        if (rownames(summary(temp_mul_chemo_drug1)$coefficients)[nrow(summary(temp_mul_chemo_drug1)$coefficients)] == "drug"){
          estimate_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),1]
          Std_Error_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),2]
          t_value_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),3]
          p_value_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),4]
        }else{
          estimate_mul_chemo_drug[x,y]<-NA
          Std_Error_mul_chemo_drug[x,y]<-NA
          t_value_mul_chemo_drug[x,y]<-NA
          p_value_mul_chemo_drug[x,y]<-NA
        }
      }else if(length(unique(temp_mul_chemo_drug$time))!=1 & length(unique(temp_mul_chemo_drug$amplicon))!=1 & length(unique(temp_mul_chemo_drug$sex))==1 & length(unique(temp_mul_chemo_drug$diagnosis))==1){
        temp_mul_chemo_drug1<-lm(temp_mul_chemo_drug[,1] ~ amplicon + time +  Antibiotic+ drug,temp_mul_chemo_drug)
        if (rownames(summary(temp_mul_chemo_drug1)$coefficients)[nrow(summary(temp_mul_chemo_drug1)$coefficients)] == "drug"){
          estimate_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),1]
          Std_Error_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),2]
          t_value_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),3]
          p_value_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),4]
        }else{
          estimate_mul_chemo_drug[x,y]<-NA
          Std_Error_mul_chemo_drug[x,y]<-NA
          t_value_mul_chemo_drug[x,y]<-NA
          p_value_mul_chemo_drug[x,y]<-NA
        }
      }else if(length(unique(temp_mul_chemo_drug$time))==1 & length(unique(temp_mul_chemo_drug$amplicon))==1 & length(unique(temp_mul_chemo_drug$sex))==1 & length(unique(temp_mul_chemo_drug$diagnosis))!=1){
        temp_mul_chemo_drug1<-lm(temp_mul_chemo_drug[,1] ~  diagnosis + Antibiotic+ drug,temp_mul_chemo_drug)
        if (rownames(summary(temp_mul_chemo_drug1)$coefficients)[nrow(summary(temp_mul_chemo_drug1)$coefficients)] == "drug"){
          estimate_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),1]
          Std_Error_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),2]
          t_value_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),3]
          p_value_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),4]
        }else{
          estimate_mul_chemo_drug[x,y]<-NA
          Std_Error_mul_chemo_drug[x,y]<-NA
          t_value_mul_chemo_drug[x,y]<-NA
          p_value_mul_chemo_drug[x,y]<-NA
        }
      }else if(length(unique(temp_mul_chemo_drug$time))==1 & length(unique(temp_mul_chemo_drug$amplicon))==1 & length(unique(temp_mul_chemo_drug$sex))!=1 & length(unique(temp_mul_chemo_drug$diagnosis))==1){
        temp_mul_chemo_drug1<-lm(temp_mul_chemo_drug[,1] ~ sex +  Antibiotic+ drug,temp_mul_chemo_drug)
        if (rownames(summary(temp_mul_chemo_drug1)$coefficients)[nrow(summary(temp_mul_chemo_drug1)$coefficients)] == "drug"){
          estimate_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),1]
          Std_Error_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),2]
          t_value_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),3]
          p_value_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),4]
        }else{
          estimate_mul_chemo_drug[x,y]<-NA
          Std_Error_mul_chemo_drug[x,y]<-NA
          t_value_mul_chemo_drug[x,y]<-NA
          p_value_mul_chemo_drug[x,y]<-NA
        }
      }else if(length(unique(temp_mul_chemo_drug$time))==1 & length(unique(temp_mul_chemo_drug$amplicon))!=1 & length(unique(temp_mul_chemo_drug$sex))==1 & length(unique(temp_mul_chemo_drug$diagnosis))==1){
        temp_mul_chemo_drug1<-lm(temp_mul_chemo_drug[,1] ~ amplicon +  Antibiotic+ drug,temp_mul_chemo_drug)
        if (rownames(summary(temp_mul_chemo_drug1)$coefficients)[nrow(summary(temp_mul_chemo_drug1)$coefficients)] == "drug"){
          estimate_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),1]
          Std_Error_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),2]
          t_value_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),3]
          p_value_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),4]
        }else{
          estimate_mul_chemo_drug[x,y]<-NA
          Std_Error_mul_chemo_drug[x,y]<-NA
          t_value_mul_chemo_drug[x,y]<-NA
          p_value_mul_chemo_drug[x,y]<-NA
        }
      }else if(length(unique(temp_mul_chemo_drug$time))!=1 & length(unique(temp_mul_chemo_drug$amplicon))==1 & length(unique(temp_mul_chemo_drug$sex))==1 & length(unique(temp_mul_chemo_drug$diagnosis))==1){
        temp_mul_chemo_drug1<-lm(temp_mul_chemo_drug[,1] ~ time +  Antibiotic+ drug,temp_mul_chemo_drug)
        if (rownames(summary(temp_mul_chemo_drug1)$coefficients)[nrow(summary(temp_mul_chemo_drug1)$coefficients)] == "drug"){
          estimate_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),1]
          Std_Error_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),2]
          t_value_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),3]
          p_value_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),4]
        }else{
          estimate_mul_chemo_drug[x,y]<-NA
          Std_Error_mul_chemo_drug[x,y]<-NA
          t_value_mul_chemo_drug[x,y]<-NA
          p_value_mul_chemo_drug[x,y]<-NA
        }
      }else{
        temp_mul_chemo_drug1<-lm(temp_mul_chemo_drug[,1] ~  Antibiotic+ drug,temp_mul_chemo_drug)
        if (rownames(summary(temp_mul_chemo_drug1)$coefficients)[nrow(summary(temp_mul_chemo_drug1)$coefficients)] == "drug"){
          estimate_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),1]
          Std_Error_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),2]
          t_value_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),3]
          p_value_mul_chemo_drug[x,y]<-summary(temp_mul_chemo_drug1)$coefficients[nrow(summary(temp_mul_chemo_drug1)$coefficients),4]
        }else{
          estimate_mul_chemo_drug[x,y]<-NA
          Std_Error_mul_chemo_drug[x,y]<-NA
          t_value_mul_chemo_drug[x,y]<-NA
          p_value_mul_chemo_drug[x,y]<-NA
        }
      }
    }else{
      estimate_mul_chemo_drug[x,y]<-NA
      Std_Error_mul_chemo_drug[x,y]<-NA
      t_value_mul_chemo_drug[x,y]<-NA
      p_value_mul_chemo_drug[x,y]<-NA
    }
  }
}

number_mul_chemo_anti<-array(length(unique(chemo_reg_name$subject)),dim=c(98,25))
rownames(number_mul_chemo_anti)<-colnames(chemo_reg_name[,2:99])
colnames(number_mul_chemo_anti)<-colnames(chemo_reg_name[,138:162])
nonzero_mul_chemo_anti<-array(NA,dim=c(98,25))
rownames(nonzero_mul_chemo_anti)<-colnames(chemo_reg_name[,2:99])
colnames(nonzero_mul_chemo_anti)<-colnames(chemo_reg_name[,138:162])
users_mul_chemo_anti<-array(NA,dim=c(98,25))
rownames(users_mul_chemo_anti)<-colnames(chemo_reg_name[,2:99])
colnames(users_mul_chemo_anti)<-colnames(chemo_reg_name[,138:162])
nonusers_mul_chemo_anti<-array(NA,dim=c(98,25))
rownames(nonusers_mul_chemo_anti)<-colnames(chemo_reg_name[,2:99])
colnames(nonusers_mul_chemo_anti)<-colnames(chemo_reg_name[,138:162])


estimate_mul_chemo_anti<-array(0,dim=c(98,25))
rownames(estimate_mul_chemo_anti)<-colnames(chemo_reg_name[,2:99])
colnames(estimate_mul_chemo_anti)<-colnames(chemo_reg_name[,138:162])
Std_Error_mul_chemo_anti<-array(0,dim=c(98,25))
rownames(Std_Error_mul_chemo_anti)<-colnames(chemo_reg_name[,2:99])
colnames(Std_Error_mul_chemo_anti)<-colnames(chemo_reg_name[,138:162])
t_value_mul_chemo_anti<-array(0,dim=c(98,25))
rownames(t_value_mul_chemo_anti)<-colnames(chemo_reg_name[,2:99])
colnames(t_value_mul_chemo_anti)<-colnames(chemo_reg_name[,138:162])
p_value_mul_chemo_anti<-array(0,dim=c(98,25))
rownames(p_value_mul_chemo_anti)<-colnames(chemo_reg_name[,2:99])
colnames(p_value_mul_chemo_anti)<-colnames(chemo_reg_name[,138:162])


for (x in 1:98) {
  for (y in 1:25){
    genus_name<-colnames(chemo_reg_name[,2:99])
    drug_name<-colnames(chemo_reg_name[,138:162])
    
    temp_mul_chemo_anti<-select(chemo_reg_name,genus_name[x],subject,flag,time,amplicon,sex,diagnosis,Antibiotic,drug_name[y])
    temp_mul_chemo_anti$species<-temp_mul_chemo_anti[,1]
    temp_mul_chemo_nonzero_anti<- temp_mul_chemo_anti |> filter(species!=0) 
    nonzero_mul_chemo_anti[x,y] <-length(unique(temp_mul_chemo_nonzero_anti$subject))
    temp_mul_chemo_anti$drug<-temp_mul_chemo_anti[,(ncol(temp_mul_chemo_anti)-1)]
    temp_mul_chemo_users_anti<-temp_mul_chemo_anti[complete.cases(temp_mul_chemo_anti$drug),] |> filter(drug==1) 
    users_mul_chemo_anti[x,y] <- length(unique(temp_mul_chemo_users_anti$subject))
    temp_mul_chemo_nonusers_anti<-temp_mul_chemo_anti[complete.cases(temp_mul_chemo_anti$drug),] |> filter(drug==0) 
    nonusers_mul_chemo_anti[x,y] <- length(unique(temp_mul_chemo_nonusers_anti$subject))
    
    temp_mul_chemo_anti<-temp_mul_chemo_anti[complete.cases(temp_mul_chemo_anti),]
    drug<-temp_mul_chemo_anti[,ncol(temp_mul_chemo_anti)]
    Antibiotic<-temp_mul_chemo_anti[,"Antibiotic"]
    if(length(unique(drug))==2 & length(unique(Antibiotic))==2){
      if(length(unique(temp_mul_chemo_anti$time))!=1 & length(unique(temp_mul_chemo_anti$amplicon))!=1 & length(unique(temp_mul_chemo_anti$sex))!=1 & length(unique(temp_mul_chemo_anti$diagnosis))!=1){
        temp_mul_chemo_anti1<-lm(temp_mul_chemo_anti[,1] ~ time + amplicon + sex +  diagnosis + Antibiotic +drug,temp_mul_chemo_anti)
        if (rownames(summary(temp_mul_chemo_anti1)$coefficients)[nrow(summary(temp_mul_chemo_anti1)$coefficients)] == "drug"){
          summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),1]->estimate_mul_chemo_anti[x,y]
          summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),2]->Std_Error_mul_chemo_anti[x,y]
          summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),3]->t_value_mul_chemo_anti[x,y]
          summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),4]->p_value_mul_chemo_anti[x,y]
        }else{
          estimate_mul_chemo_anti[x,y]<-NA
          Std_Error_mul_chemo_anti[x,y]<-NA
          t_value_mul_chemo_anti[x,y]<-NA
          p_value_mul_chemo_anti[x,y]<-NA
        }
      }else if(length(unique(temp_mul_chemo_anti$time))==1 & length(unique(temp_mul_chemo_anti$amplicon))!=1 & length(unique(temp_mul_chemo_anti$sex))!=1 & length(unique(temp_mul_chemo_anti$diagnosis))!=1){
        temp_mul_chemo_anti1<-lm(temp_mul_chemo_anti[,1] ~ amplicon + sex +  diagnosis + Antibiotic+ drug,temp_mul_chemo_anti)
        if (rownames(summary(temp_mul_chemo_anti1)$coefficients)[nrow(summary(temp_mul_chemo_anti1)$coefficients)] == "drug"){
          estimate_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),1]
          Std_Error_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),2]
          t_value_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),3]
          p_value_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),4]
        }else{
          estimate_mul_chemo_anti[x,y]<-NA
          Std_Error_mul_chemo_anti[x,y]<-NA
          t_value_mul_chemo_anti[x,y]<-NA
          p_value_mul_chemo_anti[x,y]<-NA
        }
      }else if(length(unique(temp_mul_chemo_anti$time))!=1 & length(unique(temp_mul_chemo_anti$amplicon))==1 & length(unique(temp_mul_chemo_anti$sex))!=1 & length(unique(temp_mul_chemo_anti$diagnosis))!=1){
        temp_mul_chemo_anti1<-lm(temp_mul_chemo_anti[,1] ~ time + sex +  diagnosis + Antibiotic+ drug,temp_mul_chemo_anti)
        if (rownames(summary(temp_mul_chemo_anti1)$coefficients)[nrow(summary(temp_mul_chemo_anti1)$coefficients)] == "drug"){
          estimate_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),1]
          Std_Error_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),2]
          t_value_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),3]
          p_value_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),4]
        }else{
          estimate_mul_chemo_anti[x,y]<-NA
          Std_Error_mul_chemo_anti[x,y]<-NA
          t_value_mul_chemo_anti[x,y]<-NA
          p_value_mul_chemo_anti[x,y]<-NA
        }
      }else if(length(unique(temp_mul_chemo_anti$time))!=1 & length(unique(temp_mul_chemo_anti$amplicon))!=1 & length(unique(temp_mul_chemo_anti$sex))==1 & length(unique(temp_mul_chemo_anti$diagnosis))!=1){
        temp_mul_chemo_anti1<-lm(temp_mul_chemo_anti[,1] ~ time  +  diagnosis + Antibiotic+ drug,temp_mul_chemo_anti)
        if (rownames(summary(temp_mul_chemo_anti1)$coefficients)[nrow(summary(temp_mul_chemo_anti1)$coefficients)] == "drug"){
          estimate_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),1]
          Std_Error_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),2]
          t_value_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),3]
          p_value_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),4]
        }else{
          estimate_mul_chemo_anti[x,y]<-NA
          Std_Error_mul_chemo_anti[x,y]<-NA
          t_value_mul_chemo_anti[x,y]<-NA
          p_value_mul_chemo_anti[x,y]<-NA
        }
      }else if(length(unique(temp_mul_chemo_anti$time))!=1 & length(unique(temp_mul_chemo_anti$amplicon))!=1 & length(unique(temp_mul_chemo_anti$sex))!=1 & length(unique(temp_mul_chemo_anti$diagnosis))==1){
        temp_mul_chemo_anti1<-lm(temp_mul_chemo_anti[,1] ~ time + sex +  Antibiotic + drug,temp_mul_chemo_anti)
        if (rownames(summary(temp_mul_chemo_anti1)$coefficients)[nrow(summary(temp_mul_chemo_anti1)$coefficients)] == "drug"){
          estimate_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),1]
          Std_Error_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),2]
          t_value_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),3]
          p_value_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),4]
        }else{
          estimate_mul_chemo_anti[x,y]<-NA
          Std_Error_mul_chemo_anti[x,y]<-NA
          t_value_mul_chemo_anti[x,y]<-NA
          p_value_mul_chemo_anti[x,y]<-NA
        }
      }else if(length(unique(temp_mul_chemo_anti$time))==1 & length(unique(temp_mul_chemo_anti$amplicon))==1 & length(unique(temp_mul_chemo_anti$sex))!=1 & length(unique(temp_mul_chemo_anti$diagnosis))!=1){
        temp_mul_chemo_anti1<-lm(temp_mul_chemo_anti[,1] ~ sex +  diagnosis + Antibiotic+ drug,temp_mul_chemo_anti)
        if (rownames(summary(temp_mul_chemo_anti1)$coefficients)[nrow(summary(temp_mul_chemo_anti1)$coefficients)] == "drug"){
          estimate_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),1]
          Std_Error_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),2]
          t_value_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),3]
          p_value_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),4]
        }else{
          estimate_mul_chemo_anti[x,y]<-NA
          Std_Error_mul_chemo_anti<-NA
          t_value_mul_chemo_anti[x,y]<-NA
          p_value_mul_chemo_anti[x,y]<-NA
        }
      }else if(length(unique(temp_mul_chemo_anti$time))==1 & length(unique(temp_mul_chemo_anti$amplicon))!=1 & length(unique(temp_mul_chemo_anti$sex))==1 & length(unique(temp_mul_chemo_anti$diagnosis))!=1){
        temp_mul_chemo_anti1<-lm(temp_mul_chemo_anti[,1] ~ amplicon +  diagnosis + Antibiotic+ drug,temp_mul_chemo_anti)
        if (rownames(summary(temp_mul_chemo_anti1)$coefficients)[nrow(summary(temp_mul_chemo_anti1)$coefficients)] == "drug"){
          estimate_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),1]
          Std_Error_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),2]
          t_value_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),3]
          p_value_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),4]
        }else{
          estimate_mul_chemo_anti[x,y]<-NA
          Std_Error_mul_chemo_anti[x,y]<-NA
          t_value_mul_chemo_anti[x,y]<-NA
          p_value_mul_chemo_anti[x,y]<-NA
        }
      }else if(length(unique(temp_mul_chemo_anti$time))==1 & length(unique(temp_mul_chemo_anti$amplicon))!=1 & length(unique(temp_mul_chemo_anti$sex))!=1 & length(unique(temp_mul_chemo_anti$diagnosis))==1){
        temp_mul_chemo_anti1<-lm(temp_mul_chemo_anti[,1] ~ amplicon + sex +  Antibiotic+ drug,temp_mul_chemo_anti)
        if (rownames(summary(temp_mul_chemo_anti1)$coefficients)[nrow(summary(temp_mul_chemo_anti1)$coefficients)] == "drug"){
          estimate_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),1]
          Std_Error_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),2]
          t_value_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),3]
          p_value_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),4]
        }else{
          estimate_mul_chemo_anti[x,y]<-NA
          Std_Error_mul_chemo_anti[x,y]<-NA
          t_value_mul_chemo_anti[x,y]<-NA
          p_value_mul_chemo_anti[x,y]<-NA
        }
      }else if(length(unique(temp_mul_chemo_anti$time))!=1 & length(unique(temp_mul_chemo_anti$amplicon))==1 & length(unique(temp_mul_chemo_anti$sex))==1 & length(unique(temp_mul_chemo_anti$diagnosis))!=1){
        temp_mul_chemo_anti1<-lm(temp_mul_chemo_anti[,1] ~ time +  diagnosis + Antibiotic+ drug,temp_mul_chemo_anti)
        if (rownames(summary(temp_mul_chemo_anti1)$coefficients)[nrow(summary(temp_mul_chemo_anti1)$coefficients)] == "drug"){
          estimate_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),1]
          Std_Error_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),2]
          t_value_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),3]
          p_value_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),4]
        }else{
          estimate_mul_chemo_anti[x,y]<-NA
          Std_Error_mul_chemo_anti[x,y]<-NA
          t_value_mul_chemo_anti[x,y]<-NA
          p_value_mul_chemo_anti[x,y]<-NA
        }
      }else if(length(unique(temp_mul_chemo_anti$time))!=1 & length(unique(temp_mul_chemo_anti$amplicon))==1 & length(unique(temp_mul_chemo_anti$sex))!=1 & length(unique(temp_mul_chemo_anti$diagnosis))==1){
        temp_mul_chemo_anti1<-lm(temp_mul_chemo_anti[,1] ~ time + sex +  Antibiotic+ drug,temp_mul_chemo_anti)
        if (rownames(summary(temp_mul_chemo_anti1)$coefficients)[nrow(summary(temp_mul_chemo_anti1)$coefficients)] == "drug"){
          estimate_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),1]
          Std_Error_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),2]
          t_value_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),3]
          p_value_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),4]
        }else{
          estimate_mul_chemo_anti[x,y]<-NA
          Std_Error_mul_chemo_anti[x,y]<-NA
          t_value_mul_chemo_anti[x,y]<-NA
          p_value_mul_chemo_anti[x,y]<-NA
        }
      }else if(length(unique(temp_mul_chemo_anti$time))!=1 & length(unique(temp_mul_chemo_anti$amplicon))!=1 & length(unique(temp_mul_chemo_anti$sex))==1 & length(unique(temp_mul_chemo_anti$diagnosis))==1){
        temp_mul_chemo_anti1<-lm(temp_mul_chemo_anti[,1] ~ amplicon + time +  Antibiotic+ drug,temp_mul_chemo_anti)
        if (rownames(summary(temp_mul_chemo_anti1)$coefficients)[nrow(summary(temp_mul_chemo_anti1)$coefficients)] == "drug"){
          estimate_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),1]
          Std_Error_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),2]
          t_value_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),3]
          p_value_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),4]
        }else{
          estimate_mul_chemo_anti[x,y]<-NA
          Std_Error_mul_chemo_anti[x,y]<-NA
          t_value_mul_chemo_anti[x,y]<-NA
          p_value_mul_chemo_anti[x,y]<-NA
        }
      }else if(length(unique(temp_mul_chemo_anti$time))==1 & length(unique(temp_mul_chemo_anti$amplicon))==1 & length(unique(temp_mul_chemo_anti$sex))==1 & length(unique(temp_mul_chemo_anti$diagnosis))!=1){
        temp_mul_chemo_anti1<-lm(temp_mul_chemo_anti[,1] ~  diagnosis + Antibiotic+ drug,temp_mul_chemo_anti)
        if (rownames(summary(temp_mul_chemo_anti1)$coefficients)[nrow(summary(temp_mul_chemo_anti1)$coefficients)] == "drug"){
          estimate_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),1]
          Std_Error_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),2]
          t_value_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),3]
          p_value_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),4]
        }else{
          estimate_mul_chemo_anti[x,y]<-NA
          Std_Error_mul_chemo_anti[x,y]<-NA
          t_value_mul_chemo_anti[x,y]<-NA
          p_value_mul_chemo_anti[x,y]<-NA
        }
      }else if(length(unique(temp_mul_chemo_anti$time))==1 & length(unique(temp_mul_chemo_anti$amplicon))==1 & length(unique(temp_mul_chemo_anti$sex))!=1 & length(unique(temp_mul_chemo_anti$diagnosis))==1){
        temp_mul_chemo_anti1<-lm(temp_mul_chemo_anti[,1] ~ sex +  Antibiotic+ drug,temp_mul_chemo_anti)
        if (rownames(summary(temp_mul_chemo_anti1)$coefficients)[nrow(summary(temp_mul_chemo_anti1)$coefficients)] == "drug"){
          estimate_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),1]
          Std_Error_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),2]
          t_value_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),3]
          p_value_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),4]
        }else{
          estimate_mul_chemo_anti[x,y]<-NA
          Std_Error_mul_chemo_anti[x,y]<-NA
          t_value_mul_chemo_anti[x,y]<-NA
          p_value_mul_chemo_anti[x,y]<-NA
        }
      }else if(length(unique(temp_mul_chemo_anti$time))==1 & length(unique(temp_mul_chemo_anti$amplicon))!=1 & length(unique(temp_mul_chemo_anti$sex))==1 & length(unique(temp_mul_chemo_anti$diagnosis))==1){
        temp_mul_chemo_anti1<-lm(temp_mul_chemo_anti[,1] ~ amplicon +  Antibiotic+ drug,temp_mul_chemo_anti)
        if (rownames(summary(temp_mul_chemo_anti1)$coefficients)[nrow(summary(temp_mul_chemo_anti1)$coefficients)] == "drug"){
          estimate_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),1]
          Std_Error_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),2]
          t_value_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),3]
          p_value_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),4]
        }else{
          estimate_mul_chemo_anti[x,y]<-NA
          Std_Error_mul_chemo_anti[x,y]<-NA
          t_value_mul_chemo_anti[x,y]<-NA
          p_value_mul_chemo_anti[x,y]<-NA
        }
      }else if(length(unique(temp_mul_chemo_anti$time))!=1 & length(unique(temp_mul_chemo_anti$amplicon))==1 & length(unique(temp_mul_chemo_anti$sex))==1 & length(unique(temp_mul_chemo_anti$diagnosis))==1){
        temp_mul_chemo_anti1<-lm(temp_mul_chemo_anti[,1] ~ time +  Antibiotic+ drug,temp_mul_chemo_anti)
        if (rownames(summary(temp_mul_chemo_anti1)$coefficients)[nrow(summary(temp_mul_chemo_anti1)$coefficients)] == "drug"){
          estimate_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),1]
          Std_Error_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),2]
          t_value_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),3]
          p_value_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),4]
        }else{
          estimate_mul_chemo_anti[x,y]<-NA
          Std_Error_mul_chemo_anti[x,y]<-NA
          t_value_mul_chemo_anti[x,y]<-NA
          p_value_mul_chemo_anti[x,y]<-NA
        }
      }else{
        temp_mul_chemo_anti1<-lm(temp_mul_chemo_anti[,1] ~  Antibiotic+ drug,temp_mul_chemo_anti)
        if (rownames(summary(temp_mul_chemo_anti1)$coefficients)[nrow(summary(temp_mul_chemo_anti1)$coefficients)] == "drug"){
          estimate_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),1]
          Std_Error_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),2]
          t_value_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),3]
          p_value_mul_chemo_anti[x,y]<-summary(temp_mul_chemo_anti1)$coefficients[nrow(summary(temp_mul_chemo_anti1)$coefficients),4]
        }else{
          estimate_mul_chemo_anti[x,y]<-NA
          Std_Error_mul_chemo_anti[x,y]<-NA
          t_value_mul_chemo_anti[x,y]<-NA
          p_value_mul_chemo_anti[x,y]<-NA
        }
      }
    }else{
      estimate_mul_chemo_anti[x,y]<-NA
      Std_Error_mul_chemo_anti[x,y]<-NA
      t_value_mul_chemo_anti[x,y]<-NA
      p_value_mul_chemo_anti[x,y]<-NA
    }
  }
}


number_mul_chemo<-cbind(as.data.frame(number_mul_chemo_drug),as.data.frame(number_mul_chemo_anti))
nonzero_mul_chemo<-cbind(as.data.frame(nonzero_mul_chemo_drug),as.data.frame(nonzero_mul_chemo_anti))
users_mul_chemo<-cbind(as.data.frame(users_mul_chemo_drug),as.data.frame(users_mul_chemo_anti))
nonusers_mul_chemo<-cbind(as.data.frame(nonusers_mul_chemo_drug),as.data.frame(nonusers_mul_chemo_anti))

estimate_mul_chemo<-cbind(as.data.frame(estimate_mul_chemo_drug),as.data.frame(estimate_mul_chemo_anti))
p_value_mul_chemo<-cbind(as.data.frame(p_value_mul_chemo_drug),as.data.frame(p_value_mul_chemo_anti))
Std_Error_mul_chemo<-cbind(as.data.frame(Std_Error_mul_chemo_drug),as.data.frame(Std_Error_mul_chemo_anti))
t_value_mul_chemo<-cbind(as.data.frame(t_value_mul_chemo_drug),as.data.frame(t_value_mul_chemo_anti))


# meta analysis--------------------------------------------------------------------

taxa_control<-rownames(estimate)
taxa_gvhd<-rownames(estimate_gvhd)
taxa_chemo<-rownames(estimate_chemo)
taxa_mul_control<-rownames(estimate_mul)
taxa_mul_gvhd<-rownames(estimate_mul_gvhd)
taxa_mul_chemo<-rownames(estimate_mul_chemo)

taxa_all<-union(taxa_control,taxa_chemo)
taxa_all<-union(taxa_all,taxa_gvhd)
taxa_mul_all<-union(taxa_mul_control,taxa_mul_chemo)
taxa_mul_all<-union(taxa_mul_all,taxa_mul_gvhd)

drug_all<-union(colnames(estimate),colnames(estimate_gvhd))
drug_all<-union(colnames(estimate_chemo),drug_all)
drug_mul_all<-union(colnames(estimate_mul),colnames(estimate_mul_gvhd))
drug_mul_all<-union(colnames(estimate_mul_chemo),drug_mul_all)


est_control_bu_r<-as.data.frame(array(NA,dim=c(23,48)))
rownames(est_control_bu_r)<-c(setdiff(taxa_all,taxa_control))
colnames(est_control_bu_r)<-colnames(estimate)
estimate_control_full<-rbind(estimate,est_control_bu_r)
std_control_full<-rbind(Std_Error,est_control_bu_r)
pva_control_full<-rbind(p_value,est_control_bu_r)
tva_control_full<-rbind(t_value,est_control_bu_r)
number_hsct_full<-rbind(number_hsct,est_control_bu_r)
nonzero_hsct_full<-rbind(nonzero_hsct,est_control_bu_r)
users_hsct_full<-rbind(users_hsct,est_control_bu_r)
nonusers_hsct_full<-rbind(nonusers_hsct,est_control_bu_r)

est_mul_control_bu_r<-as.data.frame(array(NA,dim=c(23,47)))
rownames(est_mul_control_bu_r)<-c(setdiff(taxa_mul_all,taxa_mul_control))
colnames(est_mul_control_bu_r)<-colnames(estimate_mul)
estimate_mul_control_full<-rbind(estimate_mul,est_mul_control_bu_r)
std_mul_control_full<-rbind(Std_Error_mul,est_mul_control_bu_r)
pva_mul_control_full<-rbind(p_value_mul,est_mul_control_bu_r)
tva_mul_control_full<-rbind(t_value_mul,est_mul_control_bu_r)
number_mul_hsct_full<-rbind(number_mul_hsct,est_mul_control_bu_r)
nonzero_mul_hsct_full<-rbind(nonzero_mul_hsct,est_mul_control_bu_r)
users_mul_hsct_full<-rbind(users_mul_hsct,est_mul_control_bu_r)
nonusers_mul_hsct_full<-rbind(nonusers_mul_hsct,est_mul_control_bu_r)



est_gvhd_bu_r<-as.data.frame(array(NA,dim=c(16,48)))
rownames(est_gvhd_bu_r)<-c(setdiff(taxa_all,taxa_gvhd))
colnames(est_gvhd_bu_r)<-colnames(estimate_gvhd)
estimate_gvhd_full<-rbind(estimate_gvhd,est_gvhd_bu_r)
std_gvhd_full<-rbind(Std_Error_gvhd,est_gvhd_bu_r)
pva_gvhd_full<-rbind(p_value_gvhd,est_gvhd_bu_r)
tva_gvhd_full<-rbind(t_value_gvhd,est_gvhd_bu_r)
number_gvhd_full<-rbind(number_gvhd,est_gvhd_bu_r)
nonzero_gvhd_full<-rbind(nonzero_gvhd,est_gvhd_bu_r)
users_gvhd_full<-rbind(users_gvhd,est_gvhd_bu_r)
nonusers_gvhd_full<-rbind(nonusers_gvhd,est_gvhd_bu_r)

est_mul_gvhd_bu_r<-as.data.frame(array(NA,dim=c(16,47)))
rownames(est_mul_gvhd_bu_r)<-c(setdiff(taxa_mul_all,taxa_mul_gvhd))
colnames(est_mul_gvhd_bu_r)<-colnames(estimate_mul_gvhd)
estimate_mul_gvhd_full<-rbind(estimate_mul_gvhd,est_mul_gvhd_bu_r)
std_mul_gvhd_full<-rbind(Std_Error_mul_gvhd,est_mul_gvhd_bu_r)
pva_mul_gvhd_full<-rbind(p_value_mul_gvhd,est_mul_gvhd_bu_r)
tva_mul_gvhd_full<-rbind(t_value_mul_gvhd,est_mul_gvhd_bu_r)
number_mul_gvhd_full<-rbind(number_mul_gvhd,est_mul_gvhd_bu_r)
nonzero_mul_gvhd_full<-rbind(nonzero_mul_gvhd,est_mul_gvhd_bu_r)
users_mul_gvhd_full<-rbind(users_mul_gvhd,est_mul_gvhd_bu_r)
nonusers_mul_gvhd_full<-rbind(nonusers_mul_gvhd,est_mul_gvhd_bu_r)



est_chemo_bu_r<-as.data.frame(array(NA,dim=c(29,48)))
rownames(est_chemo_bu_r)<-c(setdiff(taxa_all,taxa_chemo))
colnames(est_chemo_bu_r)<-colnames(estimate_chemo)
estimate_chemo_full<-rbind(estimate_chemo,est_chemo_bu_r)
std_chemo_full<-rbind(Std_Error_chemo,est_chemo_bu_r)
pva_chemo_full<-rbind(p_value_chemo,est_chemo_bu_r)
tva_chemo_full<-rbind(t_value_chemo,est_chemo_bu_r)
number_chemo_full<-rbind(number_chemo,est_chemo_bu_r)
nonzero_chemo_full<-rbind(nonzero_chemo,est_chemo_bu_r)
users_chemo_full<-rbind(users_chemo,est_chemo_bu_r)
nonusers_chemo_full<-rbind(nonusers_chemo,est_chemo_bu_r)

est_mul_chemo_bu_r<-as.data.frame(array(NA,dim=c(29,47)))
rownames(est_mul_chemo_bu_r)<-c(setdiff(taxa_mul_all,taxa_mul_chemo))
colnames(est_mul_chemo_bu_r)<-colnames(estimate_mul_chemo)
estimate_mul_chemo_full<-rbind(estimate_mul_chemo,est_mul_chemo_bu_r)
std_mul_chemo_full<-rbind(Std_Error_mul_chemo,est_mul_chemo_bu_r)
pva_mul_chemo_full<-rbind(p_value_mul_chemo,est_mul_chemo_bu_r)
tva_mul_chemo_full<-rbind(t_value_mul_chemo,est_mul_chemo_bu_r)
number_mul_chemo_full<-rbind(number_mul_chemo,est_mul_chemo_bu_r)
nonzero_mul_chemo_full<-rbind(nonzero_mul_chemo,est_mul_chemo_bu_r)
users_mul_chemo_full<-rbind(users_mul_chemo,est_mul_chemo_bu_r)
nonusers_mul_chemo_full<-rbind(nonusers_mul_chemo,est_mul_chemo_bu_r)



estimate_control_full<-arrange(estimate_control_full,rownames(estimate_control_full))
estimate_control_full<-dplyr::select(estimate_control_full,Prednisolone:Idarubicin,everything())
estimate_mul_control_full<-arrange(estimate_mul_control_full,rownames(estimate_mul_control_full))
estimate_mul_control_full<-dplyr::select(estimate_mul_control_full,Prednisolone:Idarubicin,everything())

std_control_full<-arrange(std_control_full,rownames(std_control_full))
std_control_full<-select(std_control_full,Prednisolone:Idarubicin,everything())
std_mul_control_full<-arrange(std_mul_control_full,rownames(std_mul_control_full))
std_mul_control_full<-select(std_mul_control_full,Prednisolone:Idarubicin,everything())

pva_control_full<-arrange(pva_control_full,rownames(pva_control_full))
pva_control_full<-select(pva_control_full,Prednisolone:Idarubicin,everything())
pva_mul_control_full<-arrange(pva_mul_control_full,rownames(pva_mul_control_full))
pva_mul_control_full<-select(pva_mul_control_full,Prednisolone:Idarubicin,everything())

tva_control_full<-arrange(tva_control_full,rownames(tva_control_full))
tva_control_full<-select(tva_control_full,Prednisolone:Idarubicin,everything())
tva_mul_control_full<-arrange(tva_mul_control_full,rownames(tva_mul_control_full))
tva_mul_control_full<-select(tva_mul_control_full,Prednisolone:Idarubicin,everything())

number_hsct_full<-arrange(number_hsct_full,rownames(number_hsct_full))
number_hsct_full<-select(number_hsct_full,Prednisolone:Idarubicin,everything())
number_mul_hsct_full<-arrange(number_mul_hsct_full,rownames(number_mul_hsct_full))
number_mul_hsct_full<-select(number_mul_hsct_full,Prednisolone:Idarubicin,everything())

nonzero_hsct_full<-arrange(nonzero_hsct_full,rownames(nonzero_hsct_full))
nonzero_hsct_full<-select(nonzero_hsct_full,Prednisolone:Idarubicin,everything())
nonzero_mul_hsct_full<-arrange(nonzero_mul_hsct_full,rownames(nonzero_mul_hsct_full))
nonzero_mul_hsct_full<-select(nonzero_mul_hsct_full,Prednisolone:Idarubicin,everything())

users_hsct_full<-arrange(users_hsct_full,rownames(users_hsct_full))
users_hsct_full<-select(users_hsct_full,Prednisolone:Idarubicin,everything())
users_mul_hsct_full<-arrange(users_mul_hsct_full,rownames(users_mul_hsct_full))
users_mul_hsct_full<-select(users_mul_hsct_full,Prednisolone:Idarubicin,everything())

nonusers_hsct_full<-arrange(nonusers_hsct_full,rownames(nonusers_hsct_full))
nonusers_hsct_full<-select(nonusers_hsct_full,Prednisolone:Idarubicin,everything())
nonusers_mul_hsct_full<-arrange(nonusers_mul_hsct_full,rownames(nonusers_mul_hsct_full))
nonusers_mul_hsct_full<-select(nonusers_mul_hsct_full,Prednisolone:Idarubicin,everything())


estimate_gvhd_full<-arrange(estimate_gvhd_full,rownames(estimate_gvhd_full))
estimate_gvhd_full<-select(estimate_gvhd_full,Prednisolone:Idarubicin,everything())
estimate_mul_gvhd_full<-arrange(estimate_mul_gvhd_full,rownames(estimate_mul_gvhd_full))
estimate_mul_gvhd_full<-select(estimate_mul_gvhd_full,Prednisolone:Idarubicin,everything())

std_gvhd_full<-arrange(std_gvhd_full,rownames(std_gvhd_full))
std_gvhd_full<-select(std_gvhd_full,Prednisolone:Idarubicin,everything())
std_mul_gvhd_full<-arrange(std_mul_gvhd_full,rownames(std_mul_gvhd_full))
std_mul_gvhd_full<-select(std_mul_gvhd_full,Prednisolone:Idarubicin,everything())

pva_gvhd_full<-arrange(pva_gvhd_full,rownames(pva_gvhd_full))
pva_gvhd_full<-select(pva_gvhd_full,Prednisolone:Idarubicin,everything())
pva_mul_gvhd_full<-arrange(pva_mul_gvhd_full,rownames(pva_mul_gvhd_full))
pva_mul_gvhd_full<-select(pva_mul_gvhd_full,Prednisolone:Idarubicin,everything())

tva_gvhd_full<-arrange(tva_gvhd_full,rownames(tva_gvhd_full))
tva_gvhd_full<-select(tva_gvhd_full,Prednisolone:Idarubicin,everything())
tva_mul_gvhd_full<-arrange(tva_mul_gvhd_full,rownames(tva_mul_gvhd_full))
tva_mul_gvhd_full<-select(tva_mul_gvhd_full,Prednisolone:Idarubicin,everything())

number_gvhd_full<-arrange(number_gvhd_full,rownames(number_gvhd_full))
number_gvhd_full<-select(number_gvhd_full,Prednisolone:Idarubicin,everything())
number_mul_gvhd_full<-arrange(number_mul_gvhd_full,rownames(number_mul_gvhd_full))
number_mul_gvhd_full<-select(number_mul_gvhd_full,Prednisolone:Idarubicin,everything())

nonzero_gvhd_full<-arrange(nonzero_gvhd_full,rownames(nonzero_gvhd_full))
nonzero_gvhd_full<-select(nonzero_gvhd_full,Prednisolone:Idarubicin,everything())
nonzero_mul_gvhd_full<-arrange(nonzero_mul_gvhd_full,rownames(nonzero_mul_gvhd_full))
nonzero_mul_gvhd_full<-select(nonzero_mul_gvhd_full,Prednisolone:Idarubicin,everything())

users_gvhd_full<-arrange(users_gvhd_full,rownames(users_gvhd_full))
users_gvhd_full<-select(users_gvhd_full,Prednisolone:Idarubicin,everything())
users_mul_gvhd_full<-arrange(users_mul_gvhd_full,rownames(users_mul_gvhd_full))
users_mul_gvhd_full<-select(users_mul_gvhd_full,Prednisolone:Idarubicin,everything())

nonusers_gvhd_full<-arrange(nonusers_gvhd_full,rownames(nonusers_gvhd_full))
nonusers_gvhd_full<-select(nonusers_gvhd_full,Prednisolone:Idarubicin,everything())
nonusers_mul_gvhd_full<-arrange(nonusers_mul_gvhd_full,rownames(nonusers_mul_gvhd_full))
nonusers_mul_gvhd_full<-select(nonusers_mul_gvhd_full,Prednisolone:Idarubicin,everything())


estimate_chemo_full<-arrange(estimate_chemo_full,rownames(estimate_chemo_full))
#estimate_chemo_full<-select(estimate_chemo_full,prednisolone:idarubicin,everything())
estimate_mul_chemo_full<-arrange(estimate_mul_chemo_full,rownames(estimate_mul_chemo_full))
#estimate_chemo_full<-select(estimate_chemo_full,prednisolone:idarubicin,everything())

std_chemo_full<-arrange(std_chemo_full,rownames(std_chemo_full))
#std_chemo_full<-select(std_chemo_full,prednisolone:idarubicin,everything())
std_mul_chemo_full<-arrange(std_mul_chemo_full,rownames(std_mul_chemo_full))
#std_chemo_full<-select(std_chemo_full,prednisolone:idarubicin,everything())

pva_chemo_full<-arrange(pva_chemo_full,rownames(pva_chemo_full))
#pva_chemo_full<-select(pva_chemo_full,prednisolone:idarubicin,everything())
pva_mul_chemo_full<-arrange(pva_mul_chemo_full,rownames(pva_mul_chemo_full))
#pva_chemo_full<-select(pva_chemo_full,prednisolone:idarubicin,everything())

tva_chemo_full<-arrange(tva_chemo_full,rownames(tva_chemo_full))
tva_mul_chemo_full<-arrange(tva_mul_chemo_full,rownames(tva_mul_chemo_full))

number_chemo_full<-arrange(number_chemo_full,rownames(number_chemo_full))
# number_chemo_full<-select(number_chemo_full,Prednisolone:Idarubicin,everything())
number_mul_chemo_full<-arrange(number_mul_chemo_full,rownames(number_mul_chemo_full))
# number_mul_chemo_full<-select(number_mul_chemo_full,Prednisolone:Idarubicin,everything())

nonzero_chemo_full<-arrange(nonzero_chemo_full,rownames(nonzero_chemo_full))
# nonzero_chemo_full<-select(nonzero_chemo_full,Prednisolone:Idarubicin,everything())
nonzero_mul_chemo_full<-arrange(nonzero_mul_chemo_full,rownames(nonzero_mul_chemo_full))
# nonzero_mul_chemo_full<-select(nonzero_mul_chemo_full,Prednisolone:Idarubicin,everything())

users_chemo_full<-arrange(users_chemo_full,rownames(users_chemo_full))
# users_chemo_full<-select(users_chemo_full,Prednisolone:Idarubicin,everything())
users_mul_chemo_full<-arrange(users_mul_chemo_full,rownames(users_mul_chemo_full))
# users_mul_chemo_full<-select(users_mul_chemo_full,Prednisolone:Idarubicin,everything())

nonusers_chemo_full<-arrange(nonusers_chemo_full,rownames(nonusers_chemo_full))
# nonusers_chemo_full<-select(nonusers_chemo_full,Prednisolone:Idarubicin,everything())
nonusers_mul_chemo_full<-arrange(nonusers_mul_chemo_full,rownames(nonusers_mul_chemo_full))
# nonusers_mul_chemo_full<-select(nonusers_mul_chemo_full,Prednisolone:Idarubicin,everything())



## supplementary table -----------------------------------------------------


taxa_rownames<-rownames(estimate_control_full)
colnames<-c("N_hsct","Nonzeros_hsct","Users_hsct","Nonusers_hsct",
            "est_hsct","std_hsct","pval_hsct","tval_hsct",
            "N_gvhd","Nonzeros_gvhd","Users_gvhd","Nonusers_gvhd",
            "est_gvhd","std_gvhd","pval_gvhd","tval_gvhd",
            "N_chemo","Nonzeros_chemo","Users_chemo","Nonusers_chemo",
            "est_chemo","std_chemo","pval_chemo","tval_chemo")
metaanalysis_data<-array(0,dim=c(127,24,48),dimnames=list(taxa_rownames,colnames,drug_all))
for (i in 1:48) {
  metaanalysis_data[,1,i]<-number_hsct_full[,i]
  metaanalysis_data[,2,i]<-nonzero_hsct_full[,i]
  metaanalysis_data[,3,i]<-users_hsct_full[,i]
  metaanalysis_data[,4,i]<-nonusers_hsct_full[,i]
  metaanalysis_data[,5,i]<-estimate_control_full[,i]
  metaanalysis_data[,6,i]<-std_control_full[,i]
  metaanalysis_data[,7,i]<-pva_control_full[,i]
  metaanalysis_data[,8,i]<-tva_control_full[,i]
  metaanalysis_data[,9,i]<-number_gvhd_full[,i]
  metaanalysis_data[,10,i]<-nonzero_gvhd_full[,i]
  metaanalysis_data[,11,i]<-users_gvhd_full[,i]
  metaanalysis_data[,12,i]<-nonusers_gvhd_full[,i]
  metaanalysis_data[,13,i]<-estimate_gvhd_full[,i]
  metaanalysis_data[,14,i]<-std_gvhd_full[,i]
  metaanalysis_data[,15,i]<-pva_gvhd_full[,i]
  metaanalysis_data[,16,i]<-tva_gvhd_full[,i]
  metaanalysis_data[,17,i]<-number_chemo_full[,i]
  metaanalysis_data[,18,i]<-nonzero_chemo_full[,i]
  metaanalysis_data[,19,i]<-users_chemo_full[,i]
  metaanalysis_data[,20,i]<-nonusers_chemo_full[,i]
  metaanalysis_data[,21,i]<-estimate_chemo_full[,i]
  metaanalysis_data[,22,i]<-std_chemo_full[,i]
  metaanalysis_data[,23,i]<-pva_chemo_full[,i]
  metaanalysis_data[,24,i]<-tva_chemo_full[,i]
  
}

metaanalysis_data_1<-adply(metaanalysis_data,c(1,3))
colnames(metaanalysis_data_1)[2]<-"drug"
colnames(metaanalysis_data_1)[1]<-"taxa"

taxa_rownames_mul<-rownames(estimate_mul_control_full)
colnames_mul<-c("N_hsct","Nonzeros_hsct","Users_hsct","Nonusers_hsct",
                "est_mul_hsct","std_mul_hsct","pval_mul_hsct","tval_mul_hsct",
                "N_gvhd","Nonzeros_gvhd","Users_gvhd","Nonusers_gvhd",
                "est_mul_gvhd","std_mul_gvhd","pval_mul_gvhd","tval_mul_gvhd",
                "N_chemo","Nonzeros_chemo","Users_chemo","Nonusers_chemo",
                "est_mul_chemo","std_mul_chemo","pval_mul_chemo","tval_mul_chemo")
metaanalysis_data_mul<-array(0,dim=c(127,24,47),dimnames=list(taxa_rownames_mul,colnames_mul,drug_mul_all))
for (i in 1:47) {
  metaanalysis_data_mul[,1,i]<-number_mul_hsct_full[,i]
  metaanalysis_data_mul[,2,i]<-nonzero_mul_hsct_full[,i]
  metaanalysis_data_mul[,3,i]<-users_mul_hsct_full[,i]
  metaanalysis_data_mul[,4,i]<-nonusers_mul_hsct_full[,i]
  metaanalysis_data_mul[,5,i]<-estimate_mul_control_full[,i]
  metaanalysis_data_mul[,6,i]<-std_mul_control_full[,i]
  metaanalysis_data_mul[,7,i]<-pva_mul_control_full[,i]
  metaanalysis_data_mul[,8,i]<-tva_mul_control_full[,i]
  metaanalysis_data_mul[,9,i]<-number_mul_gvhd_full[,i]
  metaanalysis_data_mul[,10,i]<-nonzero_mul_gvhd_full[,i]
  metaanalysis_data_mul[,11,i]<-users_mul_gvhd_full[,i]
  metaanalysis_data_mul[,12,i]<-nonusers_mul_gvhd_full[,i]
  metaanalysis_data_mul[,13,i]<-estimate_mul_gvhd_full[,i]
  metaanalysis_data_mul[,14,i]<-std_mul_gvhd_full[,i]
  metaanalysis_data_mul[,15,i]<-pva_mul_gvhd_full[,i]
  metaanalysis_data_mul[,16,i]<-tva_mul_gvhd_full[,i]
  metaanalysis_data_mul[,17,i]<-number_mul_chemo_full[,i]
  metaanalysis_data_mul[,18,i]<-nonzero_mul_chemo_full[,i]
  metaanalysis_data_mul[,19,i]<-users_mul_chemo_full[,i]
  metaanalysis_data_mul[,20,i]<-nonusers_mul_chemo_full[,i]
  metaanalysis_data_mul[,21,i]<-estimate_mul_chemo_full[,i]
  metaanalysis_data_mul[,22,i]<-std_mul_chemo_full[,i]
  metaanalysis_data_mul[,23,i]<-pva_mul_chemo_full[,i]
  metaanalysis_data_mul[,24,i]<-tva_mul_chemo_full[,i]
  
}

metaanalysis_data_mul1<-adply(metaanalysis_data_mul,c(1,3))
colnames(metaanalysis_data_mul1)[2]<-"drug"
colnames(metaanalysis_data_mul1)[1]<-"taxa"


drugmeta_hsctgvhd<-function(a){
  a<-a[,1:20]
  # a<-a[complete.cases(a$N_hsct),]
  a$inverse_var.hsct<-1/a$std_hsct^2
  a$inverse_var.gvhd<-1/a$std_gvhd^2
  #Calculate SE  #53
  a$se=sqrt((a$inverse_var.hsct+a$inverse_var.gvhd))
  
  #Calculate Beta #54
  a$beta=(a$inverse_var.hsct*a$est_hsct+a$inverse_var.gvhd*a$est_gvhd)/(a$inverse_var.hsct+a$inverse_var.gvhd)
  #	selection$beta=(selection$inverse_var.ibd*selection$coeff.correcting.all.IBD+selection$inverse_var.lld*selection$coeff.correcting.all.LLD)/(selection$inverse_var.ibd+selection$inverse_var.lld)
  
  #Calculate Z-score #55
  a$Z=a$beta/a$se
  
  #Calculate meta p-value #56
  a$P=2*pnorm(-abs(a$Z))
  #Remove data only avalable in one cohort
  # a=a[complete.cases(a$P), ]
  
  #Adjust pvalue with FDR #57
  a$FDR=p.adjust(a$P,method = "fdr")
  
  #Create empty columns
  a$Het.Q="NS" 
  a$Het.I2="NS" 
  a$Het.Pval="NS" 
  
  #Heterogeneity using Cochran's Q-test for meta-FDR < 0.1
  for (i in 1:length(rownames(a))){
    if ((a$FDR[i]<0.1)&(!is.na(a$FDR[i]))){
      #Select coefficients
      TE=c( a[i,5], a[i,13])
      #Select Standart error
      SE=c( a[i,6], a[i,14])
      het=metagen(TE,SE)
      #Change the number here depending of the number of columns in your data, should match with column Het.I2
      a[i,29]=het$I2
      #Match Het.Q column number
      a[i,28]=het$Q 
      #Calculate p-value from Q calculation // Het.Pval
      a[i,30]=pchisq(het$Q,df=2,lower.tail=F)
    }else{
      #Change the number here depending of the number of columns in your data, should match with column Het.I2
      a[i,29]=NA
      #Match Het.Q column number
      a[i,28]=NA 
      #Calculate p-value from Q calculation // Het.Pval
      a[i,30]=NA
      
    }
  }
  
  a <- a %>%
    mutate(Taxa=rownames(a),
           Drug="a",
           est_chemo="NS",
           std_chemo="NS",
           pval_chemo="NS",
           tval_chemo="NS",
           inverse_var.chemo="NS") %>%
    select(Taxa,Drug,
           N_hsct,Nonzeros_hsct,Users_hsct,Nonusers_hsct,
           est_hsct,std_hsct,pval_hsct,tval_hsct,
           N_gvhd,Nonzeros_gvhd,Users_gvhd,Nonusers_gvhd,
           est_gvhd,std_gvhd,pval_gvhd,tval_gvhd,
           N_chemo,Nonzeros_chemo,Users_chemo,Nonusers_chemo,
           est_chemo,std_chemo,pval_chemo,tval_chemo,
           inverse_var.hsct,inverse_var.gvhd,inverse_var.chemo,
           se,beta,Z,P,FDR,everything())
  return(a)
}
drugmeta_all<-function(a){
  a$inverse_var.hsct<-1/a$std_hsct^2
  a$inverse_var.gvhd<-1/a$std_gvhd^2
  a$inverse_var.chemo<-1/a$std_chemo^2
  #Calculate SE  #53
  a$se=sqrt(1/(a$inverse_var.hsct+a$inverse_var.gvhd+a$inverse_var.chemo))
  
  #Calculate Beta #54
  a$beta=(a$inverse_var.hsct*a$est_hsct+a$inverse_var.gvhd*a$est_gvhd+a$inverse_var.chemo*a$est_chemo)/(a$inverse_var.hsct+a$inverse_var.gvhd+a$inverse_var.chemo)
  #	selection$beta=(selection$inverse_var.ibd*selection$coeff.correcting.all.IBD+selection$inverse_var.lld*selection$coeff.correcting.all.LLD)/(selection$inverse_var.ibd+selection$inverse_var.lld)
  
  #Calculate Z-score #55
  a$Z=a$beta/a$se
  
  #Calculate meta p-value #56
  a$P=2*pnorm(-abs(a$Z))
  #Remove data only avalable in one cohort
  # a=a[complete.cases(a$P), ]
  
  #Adjust pvalue with FDR #57
  a$FDR=p.adjust(a$P,method = "fdr")
  
  #Create empty columns
  a$Het.Q="NS" #58
  a$Het.I2="NS" #59
  a$Het.Pval="NS" #60
  
  #Heterogeneity using Cochran's Q-test for meta-FDR < 0.1
  for (i in 1:length(rownames(a))){
    if ((a$FDR[i]<0.1)&(!is.na(a$FDR[i]))){
      #Select coefficients
      TE=c( a[i,5], a[i,13], a[i,21])
      #Select Standart error
      SE=c( a[i,6], a[i,14], a[i,22])
      het=metagen(TE,SE)
      #Change the number here depending of the number of columns in your data, should match with column Het.I2
      a[i,34]=het$I2
      #Match Het.Q column number
      a[i,33]=het$Q 
      #Calculate p-value from Q calculation // Het.Pval
      a[i,35]=pchisq(het$Q,df=2,lower.tail=F)
    }else{
      a[i,34]=NA
      #Match Het.Q column number
      a[i,33]=NA 
      #Calculate p-value from Q calculation // Het.Pval
      a[i,35]=NA
      
    }
  }
  
  a <- a %>%
    mutate(Taxa=rownames(a),Drug="a") %>%
    select(Taxa,Drug,
           N_hsct,Nonzeros_hsct,Users_hsct,Nonusers_hsct,
           est_hsct,std_hsct,pval_hsct,tval_hsct,
           N_gvhd,Nonzeros_gvhd,Users_gvhd,Nonusers_gvhd,
           est_gvhd,std_gvhd,pval_gvhd,tval_gvhd,
           N_chemo,Nonzeros_chemo,Users_chemo,Nonusers_chemo,
           est_chemo,std_chemo,pval_chemo,tval_chemo,
           inverse_var.hsct,inverse_var.gvhd,inverse_var.chemo,
           se,beta,Z,P,FDR,everything())
  
}
drugmeta_mul_hsctgvhd<-function(a){
  a<-a[,1:20]
  # a<-a[complete.cases(a$N_hsct),]
  a$inverse_mul_var.hsct<-1/a$std_mul_hsct^2
  a$inverse_mul_var.gvhd<-1/a$std_mul_gvhd^2
  #Calculate SE  #53
  a$se_mul=sqrt((a$inverse_mul_var.hsct+a$inverse_mul_var.gvhd))
  
  #Calculate Beta #54
  a$beta_mul=(a$inverse_mul_var.hsct*a$est_mul_hsct+a$inverse_mul_var.gvhd*a$est_mul_gvhd)/(a$inverse_mul_var.hsct+a$inverse_mul_var.gvhd)
  #	selection$beta=(selection$inverse_var.ibd*selection$coeff.correcting.all.IBD+selection$inverse_var.lld*selection$coeff.correcting.all.LLD)/(selection$inverse_var.ibd+selection$inverse_var.lld)
  
  #Calculate Z-score #55
  a$Z_mul=a$beta_mul/a$se_mul
  
  #Calculate meta p-value #56
  a$P_mul=2*pnorm(-abs(a$Z_mul))
  #Remove data only avalable in one cohort
  # a=a[complete.cases(a$P), ]
  
  #Adjust pvalue with FDR #57
  a$FDR_mul=p.adjust(a$P_mul,method = "fdr")
  
  #Create empty columns
  a$Het.Q_mul="NS" #58
  a$Het.I2_mul="NS" #59
  a$Het.Pval_mul="NS" #60
  
  #Heterogeneity using Cochran's Q-test for meta-FDR < 0.1
  for (i in 1:length(rownames(a))){
    if ((a$FDR_mul[i]<0.1)&(!is.na(a$FDR_mul[i]))){
      #Select coefficients
      TE=c( a[i,5], a[i,13])
      #Select Standart error
      SE=c( a[i,6], a[i,14])
      het=metagen(TE,SE)
      #Change the number here depending of the number of columns in your data, should match with column Het.I2
      a[i,29]=het$I2
      #Match Het.Q column number
      a[i,28]=het$Q 
      #Calculate p-value from Q calculation // Het.Pval
      a[i,30]=pchisq(het$Q,df=2,lower.tail=F)
    }else{
      #Change the number here depending of the number of columns in your data, should match with column Het.I2
      a[i,29]=NA
      #Match Het.Q column number
      a[i,28]=NA 
      #Calculate p-value from Q calculation // Het.Pval
      a[i,30]=NA
      
    }
  }
  
  a <- a %>%
    mutate(Taxa=rownames(a),
           Drug="a",
           est_mul_chemo="NS",
           std_mul_chemo="NS",
           pval_mul_chemo="NS",
           tval_mul_chemo="NS",
           inverse_mul_var.chemo="NS") %>%
    select(Taxa,Drug,
           N_hsct,Nonzeros_hsct,Users_hsct,Nonusers_hsct,
           est_mul_hsct,std_mul_hsct,pval_mul_hsct,tval_mul_hsct,
           N_gvhd,Nonzeros_gvhd,Users_gvhd,Nonusers_gvhd,
           est_mul_gvhd,std_mul_gvhd,pval_mul_gvhd,tval_mul_gvhd,
           N_chemo,Nonzeros_chemo,Users_chemo,Nonusers_chemo,
           est_mul_chemo,std_mul_chemo,pval_mul_chemo,tval_mul_chemo,
           inverse_mul_var.hsct,inverse_mul_var.gvhd,inverse_mul_var.chemo,
           se_mul,beta_mul,Z_mul,P_mul,FDR_mul,everything())
  return(a)
}
drug_merge<-function(a,b){
  mergetable<-merge(a,b[,c(1,7:10,15:18,23:37)],by="Taxa") |> select(
    Taxa,Drug,
    N_hsct,Nonzeros_hsct,Users_hsct,Nonusers_hsct,
    est_hsct,std_hsct,pval_hsct,tval_hsct,
    est_mul_hsct,std_mul_hsct,pval_mul_hsct,tval_mul_hsct,
    N_gvhd,Nonzeros_gvhd,Users_gvhd,Nonusers_gvhd,
    est_gvhd,std_gvhd,pval_gvhd,tval_gvhd,
    est_mul_gvhd,std_mul_gvhd,pval_mul_gvhd,tval_mul_gvhd,
    N_chemo,Nonzeros_chemo,Users_chemo,Nonusers_chemo,
    est_chemo,std_chemo,pval_chemo,tval_chemo,
    est_mul_chemo,std_mul_chemo,pval_mul_chemo,tval_mul_chemo,
    inverse_var.hsct,inverse_var.gvhd,inverse_var.chemo,
    inverse_mul_var.hsct,inverse_mul_var.gvhd,inverse_mul_var.chemo,
    se,beta,Z,P,FDR,Het.Q,Het.I2,Het.Pval,
    se_mul,beta_mul,Z_mul,P_mul,FDR_mul,Het.Q_mul,Het.I2_mul,Het.Pval_mul
  )
  return(mergetable)
}


### busulfan ----------------------------------------------------------------


busulfan<-as.data.frame(metaanalysis_data[,,"Busulfan"])
busulfan<-drugmeta_hsctgvhd(busulfan)
busulfan$Drug<-"Busulfan"
test_busulfan<-busulfan[complete.cases(busulfan$Het.Pval),] 
test_busulfan_n<-sum((as.numeric(test_busulfan$FDR) < 0.1)&(as.numeric(test_busulfan$Het.Pval)>0.05))
hp_busulfan<- test_busulfan |> filter((as.numeric(test_busulfan$FDR) < 0.1)&(as.numeric(test_busulfan$Het.Pval)>0.05)) |> 
  select(Taxa,Drug,beta,FDR) |> mutate(`Influence`=ifelse(beta<0,"Inhibition","Promotion"),anti="no")

busulfan_mul<-as.data.frame(metaanalysis_data_mul[,,"Busulfan"])
busulfan_mul<-drugmeta_mul_hsctgvhd(busulfan_mul)
busulfan_mul$Drug<-"Busulfan"
test_busulfan_mul<-busulfan_mul[complete.cases(busulfan_mul$Het.Pval_mul),] 
test_busulfan_mul_n<-sum((as.numeric(test_busulfan_mul$FDR_mul) < 0.1)&(as.numeric(test_busulfan_mul$Het.Pval_mul)>0.05))
hp_mul_busulfan<- test_busulfan_mul |> filter((as.numeric(test_busulfan_mul$FDR_mul) < 0.1)&(as.numeric(test_busulfan_mul$Het.Pval_mul)>0.05))|> 
  select(Taxa,Drug,beta_mul,FDR_mul) |> mutate(`Influence`=ifelse(beta_mul<0,"Inhibition","Promotion"),anti="no")

busulfan_merge<-drug_merge(busulfan,busulfan_mul)
write.csv(busulfan_merge,file = "drug/busulfan.csv")


### cyclophosphamide --------------------------------------------------------


cyclophosphamide<-as.data.frame(metaanalysis_data[,,"Cyclophosphamide"])
cyclophosphamide<-drugmeta_hsctgvhd(cyclophosphamide)
cyclophosphamide$Drug<-"Cyclophosphamide"
test_cyclophosphamide<-  cyclophosphamide[complete.cases(cyclophosphamide$Het.Pval),] 
test_cyclophosphamide_n<-sum((as.numeric(test_cyclophosphamide$FDR) < 0.1)&(as.numeric(test_cyclophosphamide$Het.Pval)>0.05))
hp_cyclophosphamide<- test_cyclophosphamide |> filter((as.numeric(test_cyclophosphamide$FDR) < 0.1)&(as.numeric(test_cyclophosphamide$Het.Pval)>0.05))|> 
  select(Taxa,Drug,beta,FDR) |> mutate(`Influence`=ifelse(beta<0,"Inhibition","Promotion"),anti="no")

cyclophosphamide_mul<-as.data.frame(metaanalysis_data_mul[,,"Cyclophosphamide"])
cyclophosphamide_mul<-drugmeta_mul_hsctgvhd(cyclophosphamide_mul)
cyclophosphamide_mul$Drug<-"Cyclophosphamide"
test_cyclophosphamide_mul<-  cyclophosphamide_mul[complete.cases(cyclophosphamide_mul$Het.Pval_mul),] 
test_cyclophosphamide_mul_n<-sum((as.numeric(test_cyclophosphamide_mul$FDR_mul) < 0.1)&(as.numeric(test_cyclophosphamide_mul$Het.Pval_mul)>0.05))
hp_mul_cyclophosphamide<- test_cyclophosphamide_mul |> filter((as.numeric(test_cyclophosphamide_mul$FDR_mul) < 0.1)&(as.numeric(test_cyclophosphamide_mul$Het.Pval_mul)>0.05))|> 
  select(Taxa,Drug,beta_mul,FDR_mul) |> mutate(`Influence`=ifelse(beta_mul<0,"Inhibition","Promotion"),anti="no")

cyclophosphamide_merge<-drug_merge(cyclophosphamide,cyclophosphamide_mul)
write.csv(cyclophosphamide_merge,file = "drug/cyclophosphamide.csv")


### fludarabine -------------------------------------------------------------


fludarabine<-as.data.frame(metaanalysis_data[,,"Fludarabine"])
fludarabine<-drugmeta_hsctgvhd(fludarabine)
fludarabine$Drug<-"Fludarabine"
test_fludarabine<-  fludarabine[complete.cases(fludarabine$Het.Pval),] 
test_fludarabine_n<-sum((as.numeric(test_fludarabine$FDR) < 0.1)&(as.numeric(test_fludarabine$Het.Pval)>0.05))
hp_fludarabine<- test_fludarabine |> filter((as.numeric(test_fludarabine$FDR) < 0.1)&(as.numeric(test_fludarabine$Het.Pval)>0.05))|> 
  select(Taxa,Drug,beta,FDR) |> mutate(`Influence`=ifelse(beta<0,"Inhibition","Promotion"),anti="no")

fludarabine_mul<-as.data.frame(metaanalysis_data_mul[,,"Fludarabine"])
fludarabine_mul<-drugmeta_mul_hsctgvhd(fludarabine_mul)
fludarabine_mul$Drug<-"Fludarabine"
test_fludarabine_mul<-  fludarabine_mul[complete.cases(fludarabine_mul$Het.Pval_mul),] 
test_fludarabine_mul_n<-sum((as.numeric(test_fludarabine_mul$FDR_mul) < 0.1)&(as.numeric(test_fludarabine_mul$Het.Pval_mul)>0.05))
hp_mul_fludarabine<- test_fludarabine_mul |> filter((as.numeric(test_fludarabine_mul$FDR_mul) < 0.1)&(as.numeric(test_fludarabine_mul$Het.Pval_mul)>0.05))|> 
  select(Taxa,Drug,beta_mul,FDR_mul) |> mutate(`Influence`=ifelse(beta_mul<0,"Inhibition","Promotion"),anti="no")

fludarabine_merge<-drug_merge(fludarabine,fludarabine_mul)
write.csv(fludarabine_merge,file = "drug/fludarabine.csv")


### melphalan ---------------------------------------------------------------


melphalan<-as.data.frame(metaanalysis_data[,,"Melphalan"])
melphalan<-drugmeta_hsctgvhd(melphalan)
melphalan$Drug<-"Melphalan"
test_melphalan <-  melphalan[complete.cases(melphalan$Het.Pval),] 
test_melphalan_n<-sum((as.numeric(test_melphalan$FDR) < 0.1)&(as.numeric(test_melphalan$Het.Pval)>0.05))
hp_melphalan<- test_melphalan |> filter((as.numeric(test_melphalan$FDR) < 0.1)&(as.numeric(test_melphalan$Het.Pval)>0.05))|> 
  select(Taxa,Drug,beta,FDR) |> mutate(`Influence`=ifelse(beta<0,"Inhibition","Promotion"),anti="no")

write.csv(melphalan,file = "drug/melphalan.csv")


### VP16 --------------------------------------------------------------------


VP16<-as.data.frame(metaanalysis_data[,,"VP16"])
VP16<-drugmeta_hsctgvhd(VP16)
VP16$Drug<-"VP16"
test_VP16<-  VP16[complete.cases(VP16$Het.Pval),] 
test_VP16_n<-sum((as.numeric(test_VP16$FDR) < 0.1)&(as.numeric(test_VP16$Het.Pval)>0.05))
hp_VP16<- test_VP16 |> filter((as.numeric(test_VP16$FDR) < 0.1)&(as.numeric(test_VP16$Het.Pval)>0.05))|> 
  select(Taxa,Drug,beta,FDR) |> mutate(`Influence`=ifelse(beta<0,"Inhibition","Promotion"),anti="no")

VP16_mul<-as.data.frame(metaanalysis_data_mul[,,"VP16"])
VP16_mul<-drugmeta_mul_hsctgvhd(VP16_mul)
VP16_mul$Drug<-"VP16"
test_VP16_mul<-  VP16_mul[complete.cases(VP16_mul$Het.Pval_mul),] 
test_VP16_mul_n<-sum((as.numeric(test_VP16_mul$FDR_mul) < 0.1)&(as.numeric(test_VP16_mul$Het.Pval_mul)>0.05))
hp_mul_VP16<- test_VP16_mul |> filter((as.numeric(test_VP16_mul$FDR_mul) < 0.1)&(as.numeric(test_VP16_mul$Het.Pval_mul)>0.05))|> 
  select(Taxa,Drug,beta_mul,FDR_mul) |> mutate(`Influence`=ifelse(beta_mul<0,"Inhibition","Promotion"),anti="no")

VP16_merge<-drug_merge(VP16,VP16_mul)
write.csv(VP16_merge,file = "drug/VP16.csv")


### TBI ---------------------------------------------------------------------


TBI<-as.data.frame(metaanalysis_data[,,"Total.Body.Irradiation"])
TBI<-drugmeta_hsctgvhd(TBI)
TBI$Drug<-"TBI"
test_TBI<-  TBI[complete.cases(TBI$Het.Pval),] 
test_TBI_n<-sum((as.numeric(test_TBI$FDR) < 0.1)&(as.numeric(test_TBI$Het.Pval)>0.05))
hp_TBI<- test_TBI |> filter((as.numeric(test_TBI$FDR) < 0.1)&(as.numeric(test_TBI$Het.Pval)>0.05))|> 
  select(Taxa,Drug,beta,FDR) |> mutate(`Influence`=ifelse(beta<0,"Inhibition","Promotion"),anti="no")

TBI_mul<-as.data.frame(metaanalysis_data_mul[,,"Total.Body.Irradiation"])
TBI_mul<-drugmeta_mul_hsctgvhd(TBI_mul)
TBI_mul$Drug<-"TBI"
test_TBI_mul<-  TBI_mul[complete.cases(TBI_mul$Het.Pval_mul),] 
test_TBI_mul_n<-sum((as.numeric(test_TBI_mul$FDR_mul) < 0.1)&(as.numeric(test_TBI_mul$Het.Pval_mul)>0.05))
hp_mul_TBI<- test_TBI_mul |> filter((as.numeric(test_TBI_mul$FDR_mul) < 0.1)&(as.numeric(test_TBI_mul$Het.Pval_mul)>0.05))|> 
  select(Taxa,Drug,beta_mul,FDR_mul) |> mutate(`Influence`=ifelse(beta_mul<0,"Inhibition","Promotion"),anti="no")

TBI_merge<-drug_merge(TBI,TBI_mul)
write.csv(TBI_merge,file = "drug/TBI.csv")


### Thiotepa ----------------------------------------------------------------


Thiotepa<-as.data.frame(metaanalysis_data[,,"Thiotepa"])
Thiotepa<-drugmeta_hsctgvhd(Thiotepa)
Thiotepa$Drug<-"Thiotepa"
test_Thiotepa<-  Thiotepa[complete.cases(Thiotepa$Het.Pval),] 
test_Thiotepa_n<-sum((as.numeric(test_Thiotepa$FDR) < 0.1)&(as.numeric(test_Thiotepa$Het.Pval)>0.05))
hp_Thiotepa<- test_Thiotepa |> filter((as.numeric(test_Thiotepa$FDR) < 0.1)&(as.numeric(test_Thiotepa$Het.Pval)>0.05))|> 
  select(Taxa,Drug,beta,FDR) |> mutate(`Influence`=ifelse(beta<0,"Inhibition","Promotion"),anti="no")

Thiotepa_mul<-as.data.frame(metaanalysis_data_mul[,,"Thiotepa"])
Thiotepa_mul<-drugmeta_mul_hsctgvhd(Thiotepa_mul)
Thiotepa_mul$Drug<-"Thiotepa"
test_Thiotepa_mul<-  Thiotepa_mul[complete.cases(Thiotepa_mul$Het.Pval_mul),] 
test_Thiotepa_mul_n<-sum((as.numeric(test_Thiotepa_mul$FDR_mul) < 0.1)&(as.numeric(test_Thiotepa_mul$Het.Pval_mul)>0.05))
hp_mul_Thiotepa<- test_Thiotepa_mul |> filter((as.numeric(test_Thiotepa_mul$FDR_mul) < 0.1)&(as.numeric(test_Thiotepa_mul$Het.Pval_mul)>0.05))|> 
  select(Taxa,Drug,beta_mul,FDR_mul) |> mutate(`Influence`=ifelse(beta_mul<0,"Inhibition","Promotion"),anti="no")

Thiotepa_merge<-drug_merge(Thiotepa,Thiotepa_mul)
write.csv(Thiotepa_merge,file = "drug/Thiotepa.csv")


### Methotrexate ------------------------------------------------------------


Methotrexate<-as.data.frame(metaanalysis_data[,,"Methotrexate"])
Methotrexate<-drugmeta_hsctgvhd(Methotrexate)
Methotrexate$Drug<-"Methotrexate"
test_Methotrexate<-  Methotrexate[complete.cases(Methotrexate$Het.Pval),] 
test_Methotrexate_n<-sum((as.numeric(test_Methotrexate$FDR) < 0.1)&(as.numeric(test_Methotrexate$Het.Pval)>0.05))
hp_Methotrexate<- test_Methotrexate |> filter((as.numeric(test_Methotrexate$FDR) < 0.1)&(as.numeric(test_Methotrexate$Het.Pval)>0.05))|> 
  select(Taxa,Drug,beta,FDR) |> mutate(`Influence`=ifelse(beta<0,"Inhibition","Promotion"),anti="no")

Methotrexate_mul<-as.data.frame(metaanalysis_data_mul[,,"Methotrexate"])
Methotrexate_mul<-drugmeta_mul_hsctgvhd(Methotrexate_mul)
Methotrexate_mul$Drug<-"Methotrexate"
test_Methotrexate_mul<-  Methotrexate_mul[complete.cases(Methotrexate_mul$Het.Pval_mul),] 
test_Methotrexate_mul_n<-sum((as.numeric(test_Methotrexate_mul$FDR_mul) < 0.1)&(as.numeric(test_Methotrexate_mul$Het.Pval_mul)>0.05))
hp_mul_Methotrexate<- test_Methotrexate_mul |> filter((as.numeric(test_Methotrexate_mul$FDR_mul) < 0.1)&(as.numeric(test_Methotrexate_mul$Het.Pval_mul)>0.05))|> 
  select(Taxa,Drug,beta_mul,FDR_mul) |> mutate(`Influence`=ifelse(beta_mul<0,"Inhibition","Promotion"),anti="no")

Methotrexate_merge<-drug_merge(Methotrexate,Methotrexate_mul)
write.csv(Methotrexate_merge,file = "drug/Methotrexate.csv")


### Trimethoprim.sulfamethoxazole -------------------------------------------


trimethoprim<-as.data.frame(metaanalysis_data[,,"Trimethoprim.sulfamethoxazole"])
trimethoprim<-drugmeta_hsctgvhd(trimethoprim)
trimethoprim$Drug<-"Trimethoprim.sulfamethoxazole"
test_trimethoprim<-  trimethoprim[complete.cases(trimethoprim$Het.Pval),] 
test_trimethoprim_n<-sum((as.numeric(test_trimethoprim$FDR) < 0.1)&(as.numeric(test_trimethoprim$Het.Pval)>0.05))
hp_trimethoprim<- test_trimethoprim |> filter((as.numeric(test_trimethoprim$FDR) < 0.1)&(as.numeric(test_trimethoprim$Het.Pval)>0.05))|> 
  select(Taxa,Drug,beta,FDR) |> mutate(`Influence`=ifelse(beta<0,"Inhibition","Promotion"),anti="yes")

trimethoprim_mul<-as.data.frame(metaanalysis_data_mul[,,"Trimethoprim.sulfamethoxazole"])
trimethoprim_mul<-drugmeta_mul_hsctgvhd(trimethoprim_mul)
trimethoprim_mul$Drug<-"Trimethoprim.sulfamethoxazole"
test_trimethoprim_mul<-  trimethoprim_mul[complete.cases(trimethoprim_mul$Het.Pval_mul),] 
test_trimethoprim_mul_n<-sum((as.numeric(test_trimethoprim_mul$FDR_mul) < 0.1)&(as.numeric(test_trimethoprim_mul$Het.Pval_mul)>0.05))
hp_mul_trimethoprim<- test_trimethoprim_mul |> filter((as.numeric(test_trimethoprim_mul$FDR_mul) < 0.1)&(as.numeric(test_trimethoprim_mul$Het.Pval_mul)>0.05))|> 
  select(Taxa,Drug,beta_mul,FDR_mul) |> mutate(`Influence`=ifelse(beta_mul<0,"Inhibition","Promotion"),anti="yes")

trimethoprim_merge<-drug_merge(trimethoprim,trimethoprim_mul)
write.csv(trimethoprim_merge,file = "drug/trimethoprim.csv")


### phenoxymethylpenicillin -------------------------------------------------


phenoxymethylpenicillin<-as.data.frame(metaanalysis_data[,,"Phenoxymethylpenicillin"])
phenoxymethylpenicillin<-drugmeta_hsctgvhd(phenoxymethylpenicillin)
phenoxymethylpenicillin$Drug<-"Phenoxymethylpenicillin"
test_phenoxymethylpenicillin<-  phenoxymethylpenicillin[complete.cases(phenoxymethylpenicillin$Het.Pval),] 
test_phenoxymethylpenicillin_n<-sum((as.numeric(test_phenoxymethylpenicillin$FDR) < 0.1)&(as.numeric(test_phenoxymethylpenicillin$Het.Pval)>0.05))
hp_phenoxymethylpenicillin<- test_phenoxymethylpenicillin |> filter((as.numeric(test_phenoxymethylpenicillin$FDR) < 0.1)&(as.numeric(test_phenoxymethylpenicillin$Het.Pval)>0.05))|> 
  select(Taxa,Drug,beta,FDR) |> mutate(`Influence`=ifelse(beta<0,"Inhibition","Promotion"),anti="yes")

phenoxymethylpenicillin_mul<-as.data.frame(metaanalysis_data_mul[,,"Phenoxymethylpenicillin"])
phenoxymethylpenicillin_mul<-drugmeta_mul_hsctgvhd(phenoxymethylpenicillin_mul)
phenoxymethylpenicillin_mul$Drug<-"Phenoxymethylpenicillin"
test_phenoxymethylpenicillin_mul<-  phenoxymethylpenicillin_mul[complete.cases(phenoxymethylpenicillin_mul$Het.Pval_mul),] 
test_phenoxymethylpenicillin_mul_n<-sum((as.numeric(test_phenoxymethylpenicillin_mul$FDR_mul) < 0.1)&(as.numeric(test_phenoxymethylpenicillin_mul$Het.Pval_mul)>0.05))
hp_mul_phenoxymethylpenicillin<- test_phenoxymethylpenicillin_mul |> filter((as.numeric(test_phenoxymethylpenicillin_mul$FDR_mul) < 0.1)&(as.numeric(test_phenoxymethylpenicillin_mul$Het.Pval_mul)>0.05))|> 
  select(Taxa,Drug,beta_mul,FDR_mul) |> mutate(`Influence`=ifelse(beta_mul<0,"Inhibition","Promotion"),anti="yes")

phenoxymethylpenicillin_merge<-drug_merge(phenoxymethylpenicillin,phenoxymethylpenicillin_mul)
write.csv(phenoxymethylpenicillin_merge,file = "drug/phenoxymethylpenicillin.csv")


### amoxicillin -------------------------------------------------------------


amoxicillin<-as.data.frame(metaanalysis_data[,,"Amoxicillin"])
amoxicillin<-drugmeta_hsctgvhd(amoxicillin)
amoxicillin$Drug<-"Amoxicillin"
test_amoxicillin<-  amoxicillin[complete.cases(amoxicillin$Het.Pval),] 
test_amoxicillin_n<-sum((as.numeric(test_amoxicillin$FDR) < 0.1)&(as.numeric(test_amoxicillin$Het.Pval)>0.05))
hp_amoxicillin<- test_amoxicillin |> filter((as.numeric(test_amoxicillin$FDR) < 0.1)&(as.numeric(test_amoxicillin$Het.Pval)>0.05))|> 
  select(Taxa,Drug,beta,FDR) |> mutate(`Influence`=ifelse(beta<0,"Inhibition","Promotion"),anti="yes")

amoxicillin_mul<-as.data.frame(metaanalysis_data_mul[,,"Amoxicillin"])
amoxicillin_mul<-drugmeta_mul_hsctgvhd(amoxicillin_mul)
amoxicillin_mul$Drug<-"Amoxicillin"
test_amoxicillin_mul<-  amoxicillin_mul[complete.cases(amoxicillin_mul$Het.Pval_mul),] 
test_amoxicillin_mul_n<-sum((as.numeric(test_amoxicillin_mul$FDR_mul) < 0.1)&(as.numeric(test_amoxicillin_mul$Het.Pval_mul)>0.05))
hp_mul_amoxicillin<- test_amoxicillin_mul |> filter((as.numeric(test_amoxicillin_mul$FDR_mul) < 0.1)&(as.numeric(test_amoxicillin_mul$Het.Pval_mul)>0.05))|> 
  select(Taxa,Drug,beta_mul,FDR_mul) |> mutate(`Influence`=ifelse(beta_mul<0,"Inhibition","Promotion"),anti="yes")

amoxicillin_merge<-drug_merge(amoxicillin,amoxicillin_mul)
write.csv(amoxicillin_merge,file = "drug/amoxicillin.csv")


### ceftazidime -------------------------------------------------------------


ceftazidime<-as.data.frame(metaanalysis_data[,,"Ceftazidime"])
ceftazidime<-drugmeta_hsctgvhd(ceftazidime)
ceftazidime$Drug<-"Ceftazidime"
test_ceftazidime<-  ceftazidime[complete.cases(ceftazidime$Het.Pval),] 
test_ceftazidime_n<-sum((as.numeric(test_ceftazidime$FDR) < 0.1)&(as.numeric(test_ceftazidime$Het.Pval)>0.05))
hp_ceftazidime<- test_ceftazidime |> filter((as.numeric(test_ceftazidime$FDR) < 0.1)&(as.numeric(test_ceftazidime$Het.Pval)>0.05))|> 
  select(Taxa,Drug,beta,FDR) |> mutate(`Influence`=ifelse(beta<0,"Inhibition","Promotion"),anti="yes")

ceftazidime_mul<-as.data.frame(metaanalysis_data_mul[,,"Ceftazidime"])
ceftazidime_mul<-drugmeta_mul_hsctgvhd(ceftazidime_mul)
ceftazidime_mul$Drug<-"Ceftazidime"
test_ceftazidime_mul<-  ceftazidime_mul[complete.cases(ceftazidime_mul$Het.Pval_mul),] 
test_ceftazidime_mul_n<-sum((as.numeric(test_ceftazidime_mul$FDR_mul) < 0.1)&(as.numeric(test_ceftazidime_mul$Het.Pval_mul)>0.05))
hp_mul_ceftazidime<- test_ceftazidime_mul |> filter((as.numeric(test_ceftazidime_mul$FDR_mul) < 0.1)&(as.numeric(test_ceftazidime_mul$Het.Pval_mul)>0.05))|> 
  select(Taxa,Drug,beta_mul,FDR_mul) |> mutate(`Influence`=ifelse(beta_mul<0,"Inhibition","Promotion"),anti="yes")

ceftazidime_merge<-drug_merge(ceftazidime,ceftazidime_mul)
write.csv(ceftazidime_merge,file = "drug/ceftazidime.csv")


### Ciprofloxacin -----------------------------------------------------------


Ciprofloxacin<-as.data.frame(metaanalysis_data[,,"Ciprofloxacin"])
Ciprofloxacin<-drugmeta_hsctgvhd(Ciprofloxacin)
Ciprofloxacin$Drug<-"Ciprofloxacin"
test_Ciprofloxacin<-  Ciprofloxacin[complete.cases(Ciprofloxacin$Het.Pval),] 
test_Ciprofloxacin_n<-sum((as.numeric(test_Ciprofloxacin$FDR) < 0.1)&(as.numeric(test_Ciprofloxacin$Het.Pval)>0.05))
hp_Ciprofloxacin<- test_Ciprofloxacin |> filter((as.numeric(test_Ciprofloxacin$FDR) < 0.1)&(as.numeric(test_Ciprofloxacin$Het.Pval)>0.05))|> 
  select(Taxa,Drug,beta,FDR) |> mutate(`Influence`=ifelse(beta<0,"Inhibition","Promotion"),anti="yes")

Ciprofloxacin_mul<-as.data.frame(metaanalysis_data_mul[,,"Ciprofloxacin"])
Ciprofloxacin_mul<-drugmeta_mul_hsctgvhd(Ciprofloxacin_mul)
Ciprofloxacin_mul$Drug<-"Ciprofloxacin"
test_Ciprofloxacin_mul<-  Ciprofloxacin_mul[complete.cases(Ciprofloxacin_mul$Het.Pval_mul),] 
test_Ciprofloxacin_mul_n<-sum((as.numeric(test_Ciprofloxacin_mul$FDR_mul) < 0.1)&(as.numeric(test_Ciprofloxacin_mul$Het.Pval_mul)>0.05))
hp_mul_Ciprofloxacin<- test_Ciprofloxacin_mul |> filter((as.numeric(test_Ciprofloxacin_mul$FDR_mul) < 0.1)&(as.numeric(test_Ciprofloxacin_mul$Het.Pval_mul)>0.05))|> 
  select(Taxa,Drug,beta_mul,FDR_mul) |> mutate(`Influence`=ifelse(beta_mul<0,"Inhibition","Promotion"),anti="yes")

Ciprofloxacin_merge<-drug_merge(Ciprofloxacin,Ciprofloxacin_mul)
write.csv(Ciprofloxacin_merge,file = "drug/Ciprofloxacin.csv")


### metronidazole -----------------------------------------------------------


metronidazole<-as.data.frame(metaanalysis_data[,,"Metronidazole"])
# metronidazole<-metronidazole[,c(1:4,9:24)]
# metronidazole<-metronidazole[complete.cases(metronidazole),]
metronidazole$inverse_var.chemo<-1/metronidazole$std_chemo^2
metronidazole$inverse_var.gvhd<-1/metronidazole$std_gvhd^2
#Calculate SE  #53
metronidazole$se=sqrt(1/(metronidazole$inverse_var.chemo+metronidazole$inverse_var.gvhd))

#Calculate Beta #54
metronidazole$beta=(metronidazole$inverse_var.chemo*metronidazole$est_chemo+metronidazole$inverse_var.gvhd*metronidazole$est_gvhd)/(metronidazole$inverse_var.chemo+metronidazole$inverse_var.gvhd)
#	selection$beta=(selection$inverse_var.ibd*selection$coeff.correcting.all.IBD+selection$inverse_var.lld*selection$coeff.correcting.all.LLD)/(selection$inverse_var.ibd+selection$inverse_var.lld)

#Calculate Z-score #55
metronidazole$Z=metronidazole$beta/metronidazole$se

#Calculate meta p-value #56
metronidazole$P=2*pnorm(-abs(metronidazole$Z))
#Remove data only avalable in one cohort
# metronidazole=metronidazole[complete.cases(metronidazole$P), ]

#Adjust pvalue with FDR #57
metronidazole$FDR=p.adjust(metronidazole$P,method = "fdr")

#Create empty columns
metronidazole$Het.Q="NS" #58
metronidazole$Het.I2="NS" #59
metronidazole$Het.Pval="NS" #60

#Heterogeneity using Cochran's Q-test for meta-FDR < 0.1
for (i in 1:length(rownames(metronidazole))){
  if ((metronidazole$FDR[i]<0.1)&(!is.na(metronidazole$FDR[i]))){
    #Select coefficients
    TE=c( metronidazole[i,13], metronidazole[i,21])
    #Select Standart error
    SE=c( metronidazole[i,14], metronidazole[i,22])
    het=metagen(TE,SE)
    #Change the number here depending of the number of columns in your data, should match with column Het.I2
    metronidazole[i,33]=het$I2
    #Match Het.Q column number
    metronidazole[i,32]=het$Q 
    #Calculate p-value from Q calculation // Het.Pval
    metronidazole[i,34]=pchisq(het$Q,df=2,lower.tail=F)
  }else{
    metronidazole[i,33]=NA
    #Match Het.Q column number
    metronidazole[i,32]=NA 
    #Calculate p-value from Q calculation // Het.Pval
    metronidazole[i,34]=NA
    
  }
}

metronidazole <- metronidazole %>%
  mutate(Taxa=rownames(metronidazole),Drug="Metronidazole",
         # est_hsct="NS",
         # std_hsct="NS",
         # pval_hsct="NS",
         # tval_hsct="NS",
         inverse_var.hsct="NS") %>%
  select(Taxa,Drug,
         N_hsct,Nonzeros_hsct,Users_hsct,Nonusers_hsct,
         est_hsct,std_hsct,pval_hsct,tval_hsct,
         N_gvhd,Nonzeros_gvhd,Users_gvhd,Nonusers_gvhd,
         est_gvhd,std_gvhd,pval_gvhd,tval_gvhd,
         N_chemo,Nonzeros_chemo,Users_chemo,Nonusers_chemo,
         est_chemo,std_chemo,pval_chemo,tval_chemo,
         inverse_var.hsct,inverse_var.gvhd,inverse_var.chemo,
         se,beta,Z,P,FDR,everything())

test_metronidazole<-  metronidazole[complete.cases(metronidazole$Het.Pval),] 
test_metronidazole_n<-sum((as.numeric(test_metronidazole$FDR) < 0.1)&(as.numeric(test_metronidazole$Het.Pval)>0.05))
hp_metronidazole<- test_metronidazole |> filter((as.numeric(test_metronidazole$FDR) < 0.1)&(as.numeric(test_metronidazole$Het.Pval)>0.05))|> 
  select(Taxa,Drug,beta,FDR) |> mutate(`Influence`=ifelse(beta<0,"Inhibition","Promotion"),anti="yes")


write.csv(metronidazole,file="drug/metronidazole.csv")


### piperacillin ------------------------------------------------------------


piperacillin<-as.data.frame(metaanalysis_data[,,"Piperacillin.tazobactam"])
# piperacillin<-piperacillin[,c(1:4,9:12)]
# piperacillin<-piperacillin[complete.cases(piperacillin),]
piperacillin$inverse_var.chemo<-1/piperacillin$std_chemo^2
piperacillin$inverse_var.hsct<-1/piperacillin$std_hsct^2
#Calculate SE  #53
piperacillin$se=sqrt(1/(piperacillin$inverse_var.chemo+piperacillin$inverse_var.hsct))

#Calculate Beta #54
piperacillin$beta=(piperacillin$inverse_var.chemo*piperacillin$est_chemo+piperacillin$inverse_var.hsct*piperacillin$est_hsct)/(piperacillin$inverse_var.chemo+piperacillin$inverse_var.hsct)
#	selection$beta=(selection$inverse_var.ibd*selection$coeff.correcting.all.IBD+selection$inverse_var.lld*selection$coeff.correcting.all.LLD)/(selection$inverse_var.ibd+selection$inverse_var.lld)

#Calculate Z-score #55
piperacillin$Z=piperacillin$beta/piperacillin$se

#Calculate meta p-value #56
piperacillin$P=2*pnorm(-abs(piperacillin$Z))
#Remove data only avalable in one cohort
# piperacillin=piperacillin[complete.cases(piperacillin$P), ]

#Adjust pvalue with FDR #57
piperacillin$FDR=p.adjust(piperacillin$P,method = "fdr")

#Create empty columns
piperacillin$Het.Q="NS" #58
piperacillin$Het.I2="NS" #59
piperacillin$Het.Pval="NS" #60

#Heterogeneity using Cochran's Q-test for meta-FDR < 0.1
for (i in 1:length(rownames(piperacillin))){
  if ((piperacillin$FDR[i]<0.1)&(!is.na(piperacillin$FDR[i]))){
    #Select coefficients
    TE=c( piperacillin[i,5], piperacillin[i,21])
    #Select Standart error
    SE=c( piperacillin[i,6], piperacillin[i,22])
    het=metagen(TE,SE)
    #Change the number here depending of the number of columns in your data, should match with column Het.I2
    piperacillin[i,33]=het$I2
    #Match Het.Q column number
    piperacillin[i,32]=het$Q 
    #Calculate p-value from Q calculation // Het.Pval
    piperacillin[i,34]=pchisq(het$Q,df=2,lower.tail=F)
  }else{
    piperacillin[i,33]=NA
    #Match Het.Q column number
    piperacillin[i,32]=NA
    #Calculate p-value from Q calculation // Het.Pval
    piperacillin[i,34]=NA
    
  }
}

piperacillin <- piperacillin %>%
  mutate(Taxa=rownames(piperacillin),Drug="Piperacillin.tazobactam",
         # est_gvhd="NS",
         # std_gvhd="NS",
         # pval_gvhd="NS",
         # tval_gvhd="NS",
         inverse_var.gvhd="NS") %>%
  select(Taxa,Drug,
         N_hsct,Nonzeros_hsct,Users_hsct,Nonusers_hsct,
         est_hsct,std_hsct,pval_hsct,tval_hsct,
         N_gvhd,Nonzeros_gvhd,Users_gvhd,Nonusers_gvhd,
         est_gvhd,std_gvhd,pval_gvhd,tval_gvhd,
         N_chemo,Nonzeros_chemo,Users_chemo,Nonusers_chemo,
         est_chemo,std_chemo,pval_chemo,tval_chemo,
         inverse_var.hsct,inverse_var.gvhd,inverse_var.chemo,
         se,beta,Z,P,FDR,everything())

test_piperacillin<-  piperacillin[complete.cases(piperacillin$Het.Pval),] 
test_piperacillin_n<-sum((as.numeric(test_piperacillin$FDR) < 0.1)&(as.numeric(test_piperacillin$Het.Pval)>0.05))
hp_piperacillin<- test_piperacillin |> filter((as.numeric(test_piperacillin$FDR) < 0.1)&(as.numeric(test_piperacillin$Het.Pval)>0.05))|> 
  select(Taxa,Drug,beta,FDR) |> mutate(`Influence`=ifelse(beta<0,"Inhibition","Promotion"),anti="yes")

write.csv(piperacillin,file="drug/piperacillin.csv")




### meropenem ---------------------------------------------------------------


meropenem<-as.data.frame(metaanalysis_data[,,"Meropenem"])
meropenem<-drugmeta_all(meropenem)
meropenem$Drug<-"Meropenem"
test_meropenem<-  meropenem[complete.cases(meropenem$Het.Pval),] 
test_meropenem_n<-sum((as.numeric(test_meropenem$FDR) < 0.1)&(as.numeric(test_meropenem$Het.Pval)>0.05))
hp_meropenem<- test_meropenem |> filter((as.numeric(test_meropenem$FDR) < 0.1)&(as.numeric(test_meropenem$Het.Pval)>0.05))|> 
  select(Taxa,Drug,beta,FDR) |> mutate(`Influence`=ifelse(beta<0,"Inhibition","Promotion"),anti="yes")

meropenem_mul<-as.data.frame(metaanalysis_data_mul[,,"Meropenem"])
meropenem_mul<-drugmeta_mul_hsctgvhd(meropenem_mul)
meropenem_mul$Drug<-"Meropenem"
test_meropenem_mul<-  meropenem_mul[complete.cases(meropenem_mul$Het.Pval_mul),] 
test_meropenem_mul_n<-sum((as.numeric(test_meropenem_mul$FDR_mul) < 0.1)&(as.numeric(test_meropenem_mul$Het.Pval_mul)>0.05))
hp_mul_meropenem<- test_meropenem_mul |> filter((as.numeric(test_meropenem_mul$FDR_mul) < 0.1)&(as.numeric(test_meropenem_mul$Het.Pval_mul)>0.05))|> 
  select(Taxa,Drug,beta_mul,FDR_mul) |> mutate(`Influence`=ifelse(beta_mul<0,"Inhibition","Promotion"),anti="yes")

meropenem_merge<-drug_merge(meropenem,meropenem_mul)
write.csv(meropenem_merge,file = "drug/meropenem.csv")


### antibiotic --------------------------------------------------------------


antibiotic<-as.data.frame(metaanalysis_data[,,"Antibiotic"])
antibiotic<-drugmeta_all(antibiotic)
antibiotic$Drug<-"Antibiotic"
test_antibiotic<-  antibiotic[complete.cases(antibiotic$Het.Pval),] 
test_antibiotic_n<-sum((as.numeric(test_antibiotic$FDR) < 0.1)&(as.numeric(test_antibiotic$Het.Pval)>0.05))
hp_antibiotic<- test_antibiotic |> filter((as.numeric(test_antibiotic$FDR) < 0.1)&(as.numeric(test_antibiotic$Het.Pval)>0.05))|> 
  select(Taxa,Drug,beta,FDR) |> mutate(`Influence`=ifelse(beta<0,"Inhibition","Promotion"),anti="yes")
write.csv(antibiotic,file = "drug/antibiotic.csv")


### vancomycin --------------------------------------------------------------


vancomycin<-as.data.frame(metaanalysis_data[,,"Vancomycin"])
vancomycin<-drugmeta_all(vancomycin)
vancomycin$Drug<-"Vancomycin"
test_vancomycin<-  vancomycin[complete.cases(vancomycin$Het.Pval),] 
test_vancomycin_n<-sum((as.numeric(test_vancomycin$FDR) < 0.1)&(as.numeric(test_vancomycin$Het.Pval)>0.05))
hp_vancomycin<- test_vancomycin |> filter((as.numeric(test_vancomycin$FDR) < 0.1)&(as.numeric(test_vancomycin$Het.Pval)>0.05))|> 
  select(Taxa,Drug,beta,FDR) |> mutate(`Influence`=ifelse(beta<0,"Inhibition","Promotion"),anti="yes")

vancomycin_mul<-as.data.frame(metaanalysis_data_mul[,,"Vancomycin"])
vancomycin_mul<-drugmeta_mul_hsctgvhd(vancomycin_mul)
vancomycin_mul$Drug<-"Vancomycin"
test_vancomycin_mul<-  vancomycin_mul[complete.cases(vancomycin_mul$Het.Pval_mul),] 
test_vancomycin_mul_n<-sum((as.numeric(test_vancomycin_mul$FDR_mul) < 0.1)&(as.numeric(test_vancomycin_mul$Het.Pval_mul)>0.05))
hp_mul_vancomycin<- test_vancomycin_mul |> filter((as.numeric(test_vancomycin_mul$FDR_mul) < 0.1)&(as.numeric(test_vancomycin_mul$Het.Pval_mul)>0.05))|> 
  select(Taxa,Drug,beta_mul,FDR_mul) |> mutate(`Influence`=ifelse(beta_mul<0,"Inhibition","Promotion"),anti="yes")

vancomycin_merge<-drug_merge(vancomycin,vancomycin_mul)
write.csv(vancomycin_merge,file = "drug/vancomycin.csv")


## Figure3-------------------------------------------------------------------------


metaanalysis_data_1.1<-metaanalysis_data_1[!rowSums(is.na(metaanalysis_data_1[,c(7:10,15:18,23:26)]))==12,] 
n_hsct<-metaanalysis_data_1.1 %>%
  group_by(drug) %>%
  dplyr::count(pval_hsct<0.05)

n_gvhd<-metaanalysis_data_1.1 %>%
  group_by(drug) %>%
  dplyr::count(pval_gvhd<0.05)

n_chemo<-metaanalysis_data_1.1 %>%
  group_by(drug) %>%
  dplyr::count(pval_chemo<0.05)

metaanalysis_data_mul1.1<-metaanalysis_data_mul1[!rowSums(is.na(metaanalysis_data_mul1[,c(7:10,15:18,23:26)]))==12,] 
n_mul_hsct<-metaanalysis_data_mul1.1 %>%
  group_by(drug) %>%
  dplyr::count(pval_mul_hsct<0.05)

n_mul_gvhd<-metaanalysis_data_mul1.1 %>%
  group_by(drug) %>%
  dplyr::count(pval_mul_gvhd<0.05)

n_mul_chemo<-metaanalysis_data_mul1.1 %>%
  group_by(drug) %>%
  dplyr::count(pval_mul_chemo<0.05)



for (i in 1:68){
  if (is.na(n_hsct$`pval_hsct < 0.05`[i])){
    n_hsct$n[i] = 0
  }
}

test_hsct<-dplyr::filter(n_hsct,`pval_hsct < 0.05`==T |is.na(`pval_hsct < 0.05`)) %>%
  filter(!duplicated(drug) | `pval_hsct < 0.05`==T)
#test_hsct<-rbind(test_hsct,n_hsct[64,])
# test_hsct$n[30]<-0

for (i in 1:67){
  if (is.na(n_gvhd$`pval_gvhd < 0.05`[i])){
    n_gvhd$n[i] = 0
  }
}

test_gvhd<-dplyr::filter(n_gvhd,`pval_gvhd < 0.05`==T |is.na(`pval_gvhd < 0.05`)) %>%
  filter(!duplicated(drug) | `pval_gvhd < 0.05`==T)

for (i in 1:50){
  if (is.na(n_chemo$`pval_chemo < 0.05`[i])){
    n_chemo$n[i] = 0
  }
}

test_chemo<-dplyr::filter(n_chemo,`pval_chemo < 0.05`==T |is.na(`pval_chemo < 0.05`)) %>%
  filter(!duplicated(drug) | `pval_chemo < 0.05`==T)

# setdiff(n_chemo$drug[!duplicated(n_chemo$drug)],test_chemo$drug)
# test_chemo<-rbind(test_chemo,n_chemo[c(1,2,15,16),])
#test_chemo$n[c(27,28,29,30)]<-0

for (i in 1:57){
  if (is.na(n_mul_hsct$`pval_mul_hsct < 0.05`[i])){
    n_mul_hsct$n[i] = 0
  }
}

test_mul_hsct<-dplyr::filter(n_mul_hsct,`pval_mul_hsct < 0.05`==T |is.na(`pval_mul_hsct < 0.05`)) %>%
  filter(!duplicated(drug) | `pval_mul_hsct < 0.05`==T)
# test_mul_hsct<-rbind(test_mul_hsct,n_mul_hsct[56,])
# test_mul_hsct$n[25]<-0


for (i in 1:56){
  if (is.na(n_mul_gvhd$`pval_mul_gvhd < 0.05`[i])){
    n_mul_gvhd$n[i] = 0
  }
}

test_mul_gvhd<-dplyr::filter(n_mul_gvhd,`pval_mul_gvhd < 0.05`==T |is.na(`pval_mul_gvhd < 0.05`)) %>%
  filter(!duplicated(drug) | `pval_mul_gvhd < 0.05`==T)


for (i in 1:25){
  if (is.na(n_mul_chemo$`pval_mul_chemo < 0.05`[i])){
    n_mul_chemo$n[i] = 0
  }
}

test_mul_chemo<-dplyr::filter(n_mul_chemo,`pval_mul_chemo < 0.05`==T |is.na(`pval_mul_chemo < 0.05`)) %>%
  filter(!duplicated(drug) | `pval_mul_chemo < 0.05`==T)



test_allcohort<-merge(test_hsct,test_gvhd,by="drug")
test_allcohort<-merge(test_allcohort,test_chemo,by="drug")
test_allcohort<-test_allcohort[,-c(2,4,6)]
colnames(test_allcohort)<-c("drug","HSCT","GVHD","Chemo")
test_allcohort<-dplyr::mutate(test_allcohort,Meta=c(0,6,0,3,0,10,
                                                    0,4,0,0,0,5,
                                                    7,0,7,0,0,0,
                                                    13,7,20,5,4,1,
                                                    7,0,10,4,15,9))
test_allcohort<-test_allcohort[-c(1,14,16,17,27),]
test_allcohort_1<-melt(test_allcohort,id="drug")
colnames(test_allcohort_1)<-c("Meds","cohort","Associations")
# write.table(test_allcohort_1,file = "test_allcohort_1.csv")

test_allcohort_1$Cohort<-factor(test_allcohort_1$cohort,levels = c("Meta","Chemo","GVHD","HSCT"))
test_allcohort_1$Group<-"Single drug analysis"

test_mul_allcohort<-merge(test_mul_hsct,test_mul_gvhd,by="drug")
test_mul_allcohort<-merge(test_mul_allcohort,test_mul_chemo,by="drug")
test_mul_allcohort<-test_mul_allcohort[,-c(2,4,6)]
colnames(test_mul_allcohort)<-c("drug","HSCT","GVHD","Chemo")
test_mul_allcohort<-dplyr::mutate(test_mul_allcohort,Meta=c(2,0,10,4,0,0,
                                                            2,9,0,8,0,0,
                                                            3,17,0,5,0,9,
                                                            0,9,4,8,8))
test_mul_allcohort<-test_mul_allcohort[-20,]
test_mul_allcohort_1<-melt(test_mul_allcohort,id="drug")
colnames(test_mul_allcohort_1)<-c("Meds","cohort","Associations")
# write.table(test_mul_allcohort_1,file = "test_mul_allcohort_1.csv")

test_mul_allcohort_1$Cohort<-factor(test_mul_allcohort_1$cohort,levels = c("Meta","Chemo","GVHD","HSCT"))
test_mul_allcohort_1$Group<-"Multiple drug analysis"

### bar plot 3_1 -----------------------------------------------------------------
test_merge_cohort<-rbind(test_allcohort_1,test_mul_allcohort_1)
test_merge_cohort$Group<-factor(test_merge_cohort$Group,levels = c("Single drug analysis","Multiple drug analysis"))
figure3_1<-ggplot(test_merge_cohort,aes(Meds,Associations,fill=Cohort)) +
  geom_bar(stat="identity", position=position_dodge(),alpha=1.0) +
  # scale_fill_manual(values=c( "#241c1f","#CCE6D6","#CCE3F0","#EECEC9"))+
  scale_fill_manual(values=c( "#241c1f","#299555","#66AAD3","#C65949")) +
  theme(panel.background = element_blank(),
        panel.grid =element_line(colour = "grey90"),
        axis.line = element_line(colour = "grey50"),
        axis.title = element_text(size=11,family="serif",face="bold"),
        axis.text = element_text(size=11,family="serif",face="bold"),
        #legend.position = 'none',
        legend.title = element_text(size=11,family="serif",face="bold"),
        legend.text = element_text(size=11,family="serif",face="bold"),
        plot.subtitle = element_text(size=11,family="serif",face="bold"),
        strip.text = element_text(size=11,family="serif",face="bold"),
        strip.background = element_rect(fill = NA),
        title = element_text(size=11,family="serif",face="bold")
  )+
  # labs(title = "Single Drug Analysis")+
  facet_grid(cols = vars(Group), drop=T, scales = "free") +
  coord_flip() #+ scale_x_discrete(limits = rev(levels(test_allcohort$meds)))
figure3_1
ggsave("Figure0111/Figure3.pdf",figure3_1,units="in", width=9, height=5, dpi=600,limitsize = FALSE)

tmp<-test_merge_cohort |> arrange(Meds) |> mutate(Antibiotic=c(rep("no",68),rep("yes",120)))
tmp1<-tmp |> group_by(Group,cohort,Antibiotic) |> 
  summarise(
    sum=sum(Associations)
  )
tmp2<-dcast(tmp1,Group+Antibiotic~cohort)
write.csv(tmp2,"Supplementary data 17.csv")

### sankey plot 3_2 -------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(reshape2)
library(networkD3)
library(ggsci)
library(scales)

fil_hsct<-metaanalysis_data_1.1 %>%
  group_by(drug) %>%
  dplyr::filter(pval_hsct<0.05) 
fil_hsct_e<-fil_hsct[,c(1,2,7)]
fil_hsct_e<-dcast(fil_hsct_e,taxa~drug)
write.csv(fil_hsct_e,file = "F:/数据/document/学校/硕士/数据库/r数据库/gvhd/sankeyplot/fil_hsct_e_s.csv") #excel修改行名
fil_hsct_s<-fil_hsct[,c(1,2,8)]
fil_hsct_s<-dcast(fil_hsct_s,taxa~drug)
write.csv(fil_hsct_s,file = "F:/数据/document/学校/硕士/数据库/r数据库/gvhd/sankeyplot/fil_hsct_s_s.csv") #excel修改行名

fil_gvhd<-metaanalysis_data_1.1 %>%
  group_by(drug) %>%
  dplyr::filter(pval_gvhd<0.05) 
fil_gvhd_e<-fil_gvhd[,c(1,2,15)]
fil_gvhd_e<-dcast(fil_gvhd_e,taxa~drug) 
write.csv(fil_gvhd_e,file = "F:/数据/document/学校/硕士/数据库/r数据库/gvhd/sankeyplot/fil_gvhd_e_s.csv") #excel修改行名
fil_gvhd_s<-fil_gvhd[,c(1,2,16)]
fil_gvhd_s<-dcast(fil_gvhd_s,taxa~drug) 
write.csv(fil_gvhd_s,file = "F:/数据/document/学校/硕士/数据库/r数据库/gvhd/sankeyplot/fil_gvhd_s_s.csv") #excel修改行名


fil_chemo<-metaanalysis_data_1.1 %>%
  group_by(drug) %>%
  dplyr::filter(pval_chemo<0.05) 
fil_chemo_e<-fil_chemo[,c(1,2,23)]
fil_chemo_e<-dcast(fil_chemo_e,taxa~drug)
write.csv(fil_chemo_e,file = "F:/数据/document/学校/硕士/数据库/r数据库/gvhd/sankeyplot/fil_chemo_e_s.csv") #excel修改行名
fil_chemo_s<-fil_chemo[,c(1,2,24)]
fil_chemo_s<-dcast(fil_chemo_s,taxa~drug)
write.csv(fil_chemo_s,file = "F:/数据/document/学校/硕士/数据库/r数据库/gvhd/sankeyplot/fil_chemo_s_s.csv") #excel修改行名


fil_mul_hsct<-metaanalysis_data_mul1.1 %>%
  group_by(drug) %>%
  dplyr::filter(pval_mul_hsct<0.05) 
fil_mul_hsct_e<-fil_mul_hsct[,c(1,2,7)]
fil_mul_hsct_e<-dcast(fil_mul_hsct_e,taxa~drug)
write.csv(fil_mul_hsct_e,file = "F:/数据/document/学校/硕士/数据库/r数据库/gvhd/sankeyplot/fil_mul_hsct_e_s.csv") #excel修改行名
fil_mul_hsct_s<-fil_mul_hsct[,c(1,2,8)]
fil_mul_hsct_s<-dcast(fil_mul_hsct_s,taxa~drug)
write.csv(fil_mul_hsct_s,file = "F:/数据/document/学校/硕士/数据库/r数据库/gvhd/sankeyplot/fil_mul_hsct_s_s.csv") #excel修改行名


fil_mul_gvhd<-metaanalysis_data_mul1.1 %>%
  group_by(drug) %>%
  dplyr::filter(pval_mul_gvhd<0.05) 
fil_mul_gvhd_e<-fil_mul_gvhd[,c(1,2,15)]
fil_mul_gvhd_e<-dcast(fil_mul_gvhd_e,taxa~drug) 
write.csv(fil_mul_gvhd_e,file = "F:/数据/document/学校/硕士/数据库/r数据库/gvhd/sankeyplot/fil_mul_gvhd_e_s.csv") #excel修改行名
fil_mul_gvhd_s<-fil_mul_gvhd[,c(1,2,16)]
fil_mul_gvhd_s<-dcast(fil_mul_gvhd_s,taxa~drug) 
write.csv(fil_mul_gvhd_s,file = "F:/数据/document/学校/硕士/数据库/r数据库/gvhd/sankeyplot/fil_mul_gvhd_s_s.csv") #excel修改行名

fil_mul_chemo<-metaanalysis_data_mul1.1 %>%
  group_by(drug) %>%
  dplyr::filter(pval_mul_chemo<0.05) 
fil_mul_chemo_e<-fil_mul_chemo[,c(1,2,23)]
fil_mul_chemo_e<-dcast(fil_mul_chemo_e,taxa~drug)
write.csv(fil_mul_chemo_e,file = "F:/数据/document/学校/硕士/数据库/r数据库/gvhd/sankeyplot/fil_mul_chemo_e_s.csv") #excel修改行名
fil_mul_chemo_s<-fil_mul_chemo[,c(1,2,24)]
fil_mul_chemo_s<-dcast(fil_mul_chemo_s,taxa~drug)
write.csv(fil_mul_chemo_s,file = "F:/数据/document/学校/硕士/数据库/r数据库/gvhd/sankeyplot/fil_mul_chemo_s_s.csv") #excel修改行名


#### single drug analysis -------------------------------------------------------------


sk_hsct<-read.csv( "F:/数据/document/学校/硕士/数据库/r数据库/gvhd/sankeyplot/sk_hsct.csv",header=T)
sk_hsct<-sk_hsct[,-7]
sk_hsct<-melt(sk_hsct,ID=c("Taxa","NodeColor"))
sk_hsct<-na.omit(sk_hsct)
sk_hsct<-sk_hsct %>% mutate(trend=ifelse(value<0,"neg","pos")) %>% mutate(Cohort="HSCT") %>% mutate(antibioic=c(rep("no",50),rep("yes",64)))

sk_gvhd<-read.csv( "F:/数据/document/学校/硕士/数据库/r数据库/gvhd/sankeyplot/sk_gvhd.csv",header=T)
sk_gvhd<-sk_gvhd[,-c(8,18)]
sk_gvhd<-melt(sk_gvhd,ID=c("Taxa","NodeColor"))
sk_gvhd<-na.omit(sk_gvhd)
sk_gvhd<-sk_gvhd %>% mutate(trend=ifelse(value<0,"neg","pos")) %>% mutate(Cohort="GVHD") %>% mutate(antibioic=c(rep("no",27),rep("yes",66)))

sk_chemo<-read.csv( "F:/数据/document/学校/硕士/数据库/r数据库/gvhd/sankeyplot/sk_chemo.csv",header=T)
sk_chemo<-melt(sk_chemo,ID=c("Taxa","NodeColor"))
sk_chemo<-na.omit(sk_chemo)
sk_chemo<-sk_chemo %>% mutate(trend=ifelse(value<0,"neg","pos")) %>% mutate(Cohort="Chemo") %>% mutate(antibioic=c(rep("yes",27)))

sk_plt<-rbind(sk_hsct,sk_gvhd,sk_chemo)

sk_plt_link1<-sk_plt[,c(3,1,4,2,5,6,7)]
colnames(sk_plt_link1)<-c("source","target","value","node","trend","cohort","anti")
sk_plt_link2<-sk_plt[,c(1,6,4,2,5,6,7)]
colnames(sk_plt_link2)<-c("source","target","value","node","trend","cohort","anti")

node_micro<-sk_plt_link1[,c("target","node")] %>% filter(!duplicated(sk_plt_link1$target) ==T)
colnames(node_micro)=c("name","group")
node_micro$group<-factor(node_micro$group,levels=c("Bacteroidetes","Firmicutes","Actinobacteria","Proteobacteria","Verrucomicrobia"))
node_micro<-node_micro |> arrange(group)
node_drug<-sk_plt_link1[,c("source","anti")] %>% filter(!duplicated(sk_plt_link1$source) ==T)
colnames(node_drug)=c("name","group")
node_drug$name<-as.character(node_drug$name)
node_drug<-node_drug |> arrange(group,name)
node_cohort<-data.frame(name=c("HSCT","GVHD","Chemo"),group=c("HSCT","GVHD","Chemo"))
sk_plt_nodes<-rbind(node_drug,node_micro,node_cohort)
sk_plt_nodes$name<-as.character(sk_plt_nodes$name)
sk_plt_nodes$group<-as.factor(sk_plt_nodes$group)

sk_plt_link1<-sk_plt_link1 %>% mutate(group=paste(trend,anti, sep = "_", collapse = NULL))
sk_plt_link2<-sk_plt_link2 %>% mutate(group=paste(trend,cohort, sep = "_", collapse = NULL))
sk_plt_link<-rbind(sk_plt_link1,sk_plt_link2)
sk_plt_link<-sk_plt_link[,-c(4:7)]
sk_plt_link$group<-as.factor(sk_plt_link$group)
sk_plt_link$source<-as.character(sk_plt_link$source)
sk_plt_link$value<-abs(log(abs(sk_plt_link$value)))
sk_plt_link$value<-as.character(sk_plt_link$value)

sk_plt_link$IDsource <- match(sk_plt_link$source, sk_plt_nodes$name)-1 
sk_plt_link$IDtarget <- match(sk_plt_link$target, sk_plt_nodes$name)-1

library(scales)
library(ggsci)
mypal_drugcohort = pal_nejm("default", alpha = 1)(8)
mypal = pal_npg("nrc", alpha =1)(3)
show_col(pal_aaas("default", alpha = 1)(5))
show_col(pal_lancet("lanonc", alpha =1)(10))
show_col(mypal)


#### Figure3_2 -----------------------------------------------------------------


sk_plt_color <- 'd3.scaleOrdinal() .domain(["HSCT","GVHD","Chemo","yes","no",
                                            "Firmicutes","Bacteroidetes","Actinobacteria","Verrucomicrobia","Proteobacteria",
                                            "pos_HSCT","pos_GVHD", "pos_Chemo","neg_HSCT","neg_GVHD", "neg_Chemo",
                                            "pos_no","pos_yes","neg_no","neg_yes"]) 
                                   .range(["#BC3C29FF", "#0072B5FF", "#008134FF", "#FFDC91FF", "#7876B1FF",
                                          "#E64B35FF","#4DBBD5FF","#6CB86FFF","#3C5488FF","#F39B7FFF",
                                          "#BC3C294C", "#0072B54C", "#0081344C","#BC3C294C", "#0072B54C", "#0081344C",
                                          "#F39B7FFF","#F39B7F4C","#8491B4E5","#8491B44C"])'
sk <- sankeyNetwork(Links = sk_plt_link, Nodes = sk_plt_nodes, Source = "IDsource", Target = "IDtarget", 
                    Value = "value", nodePadding = 7,fontFamily = "Times New Roman",
                    NodeID = "name", iterations=0,fontSize = 14,sinksRight = TRUE,
                    colourScale=sk_plt_color, 
                    units = "TWh",height = 900,width = 1200,
                    LinkGroup="group", NodeGroup="group")
sk

write.csv(sk_plt_link,"sk_plt_link.csv")
write.csv(sk_plt_nodes,"sk_plt_nodes.csv")

#### multiple drug analysis --------------------------------------------------

sk_mul_hsct<-read.csv( "F:/数据/document/学校/硕士/数据库/r数据库/gvhd/sankeyplot/sk_mul_hsct.csv",header=T)
sk_mul_hsct<-sk_mul_hsct[,-6]
sk_mul_hsct<-melt(sk_mul_hsct,ID=c("Taxa","NodeColor"))
sk_mul_hsct<-na.omit(sk_mul_hsct)
sk_mul_hsct<-sk_mul_hsct %>% mutate(trend=ifelse(value<0,"neg","pos")) %>% mutate(Cohort="HSCT") %>% mutate(antibioic=c(rep("no",23),rep("yes",56)))

sk_mul_gvhd<-read.csv( "F:/数据/document/学校/硕士/数据库/r数据库/gvhd/sankeyplot/sk_mul_gvhd.csv",header=T)
sk_mul_gvhd<-sk_mul_gvhd[,-c(8,17)]
sk_mul_gvhd<-melt(sk_mul_gvhd,ID=c("Taxa","NodeColor"))
sk_mul_gvhd<-na.omit(sk_mul_gvhd)
sk_mul_gvhd<-sk_mul_gvhd %>% mutate(trend=ifelse(value<0,"neg","pos")) %>% mutate(Cohort="GVHD") %>% mutate(antibioic=c(rep("no",32),rep("yes",60)))

sk_mul_chemo<-read.csv( "F:/数据/document/学校/硕士/数据库/r数据库/gvhd/sankeyplot/sk_mul_chemo.csv",header=T)
sk_mul_chemo<-melt(sk_mul_chemo,ID=c("Taxa","NodeColor"))
sk_mul_chemo<-na.omit(sk_mul_chemo)
sk_mul_chemo<-sk_mul_chemo %>% mutate(trend=ifelse(value<0,"neg","pos")) %>% mutate(Cohort="Chemo") %>% mutate(antibioic=c(rep("no",5)))

sk_mul_plt<-rbind(sk_mul_hsct,sk_mul_gvhd,sk_mul_chemo)

sk_mul_plt_link1<-sk_mul_plt[,c(3,1,4,2,5,6,7)]
colnames(sk_mul_plt_link1)<-c("source","target","value","node","trend","cohort","anti")
sk_mul_plt_link2<-sk_mul_plt[,c(1,6,4,2,5,6,7)]
colnames(sk_mul_plt_link2)<-c("source","target","value","node","trend","cohort","anti")

node_mul_micro<-sk_mul_plt_link1[,c("target","node")] %>% filter(!duplicated(sk_mul_plt_link1$target) ==T)
colnames(node_mul_micro)=c("name","group")
node_mul_micro$group<-factor(node_mul_micro$group,levels=c("Bacteroidetes","Firmicutes","Actinobacteria","Proteobacteria","Verrucomicrobia"))
node_mul_micro<-node_mul_micro |> arrange(group)
node_mul_drug<-sk_mul_plt_link1[,c("source","anti")] %>% filter(!duplicated(sk_mul_plt_link1$source) ==T)
colnames(node_mul_drug)=c("name","group")
node_mul_drug$name<-as.character(node_mul_drug$name)
node_mul_drug<-node_mul_drug |> arrange(group,name)
node_mul_cohort<-data.frame(name=c("HSCT","GVHD","Chemo"),group=c("HSCT","GVHD","Chemo"))
sk_mul_plt_nodes<-rbind(node_mul_drug,node_mul_micro,node_mul_cohort)
sk_mul_plt_nodes$name<-as.character(sk_mul_plt_nodes$name)
sk_mul_plt_nodes$group<-as.factor(sk_mul_plt_nodes$group)


sk_mul_plt_link1<-sk_mul_plt_link1 %>% mutate(group=paste(trend,anti, sep = "_", collapse = NULL))
sk_mul_plt_link2<-sk_mul_plt_link2 %>% mutate(group=paste(trend,cohort, sep = "_", collapse = NULL))
sk_mul_plt_link<-rbind(sk_mul_plt_link1,sk_mul_plt_link2)
sk_mul_plt_link<-sk_mul_plt_link[,-c(4:7)]
sk_mul_plt_link$group<-as.factor(sk_mul_plt_link$group)
sk_mul_plt_link$source<-as.character(sk_mul_plt_link$source)
sk_mul_plt_link$value<-abs(log(abs(sk_mul_plt_link$value)))
sk_mul_plt_link$value<-as.character(sk_mul_plt_link$value)

sk_mul_plt_link$IDsource <- match(sk_mul_plt_link$source, sk_mul_plt_nodes$name)-1 
sk_mul_plt_link$IDtarget <- match(sk_mul_plt_link$target, sk_mul_plt_nodes$name)-1



#### Figure3_2 s ---------------------------------------------------------------


sk_mul_plt_color <- 'd3.scaleOrdinal() .domain(["HSCT","GVHD","Chemo","yes","no",
                                            "Firmicutes","Bacteroidetes","Actinobacteria","Verrucomicrobia","Proteobacteria",
                                            "pos_HSCT","pos_GVHD", "pos_Chemo","neg_HSCT","neg_GVHD", "neg_Chemo",
                                            "pos_no","pos_yes","neg_no","neg_yes"]) 
                                   .range(["#BC3C29FF", "#0072B5FF", "#008134FF", "#FFDC91FF", "#7876B1FF",
                                          "#E64B35FF","#4DBBD5FF","#6CB86FFF","#3C5488FF","#F39B7FFF",
                                          "#BC3C294C", "#0072B54C", "#0081344C","#BC3C294C", "#0072B54C", "#0081344C",
                                          "#F39B7FFF","#F39B7F4C","#8491B4E5","#8491B44C"])'

sk_mul <- sankeyNetwork(Links = sk_mul_plt_link, Nodes = sk_mul_plt_nodes, Source = "IDsource", Target = "IDtarget", 
                    Value = "value", nodePadding = 7,fontFamily = "Times New Roman",
                    NodeID = "name", iterations=0,fontSize = 14,sinksRight = TRUE,
                    colourScale=sk_mul_plt_color, 
                    units = "TWh",height = 900,width = 1200,
                    LinkGroup="group", NodeGroup="group")

sk_mul

write.csv(sk_mul_plt_link,"sk_mul_plt_link.csv")



### vennplot 3_3 ----------------------------------------------------------------

library(ggVennDiagram)

vp_hsct<-as.character(fil_hsct_e$taxa)
vp_sk_hsct<-unique(sk_hsct$Taxa)
vp_gvhd<-as.character(fil_gvhd_e$taxa)
vp_sk_gvhd<-unique(sk_gvhd$Taxa)
vp_chemo<-as.character(fil_chemo_e$taxa)
vp_sk_chemo<-unique(sk_chemo$Taxa)

intersect(intersect(vp_hsct,vp_gvhd),vp_chemo)

vp<-list(
  HSCT=vp_hsct,
  GVHD=vp_gvhd,
  Chemo=vp_chemo
)

vp_sk<-list(
  HSCT=vp_sk_hsct,
  GVHD=vp_sk_gvhd,
  Chemo=vp_sk_chemo
)

# vp_plt<-ggvenn(vp,c("HSCT","GVHD","Chemo"),show_percentage = T,
#                stroke_color = "white",
#                fill_color = c("#BC3C29","#0072B5","#008134"),
#                set_name_color =c("#BC3C29","#0072B5","#008134"))+
#   theme(text = element_text(size=16,family="serif",face="bold"))
# vp_plt
# ggsave("Figure/Figure3_3l.pdf",vp_plt,units="in", width=3, height=3, dpi=600,limitsize = FALSE)
# vp_sk_plt<-ggvenn(vp_sk,c("HSCT","GVHD","Chemo"),show_percentage = T,
# stroke_color = "white",
# fill_color = c("#BC3C29","#0072B5","#008134"),
# set_name_color =c("#BC3C29","#0072B5","#008134"))+
#   theme(text = element_text(size=16,family="serif",face="bold"))
# vp_sk_plt
# ggsave("Figure/Figure3_3r.pdf",vp_sk_plt,units="in", width=3, height=3, dpi=600,limitsize = FALSE)


vp_plt<-ggVennDiagram(vp, 
              category.names = c("HSCT","GVHD","Chemo"),
              set_color = c("#BC3C29FF", "#0072B5FF", "#008134FF"),
              set_size = 6,
              label = "both", 
              label_percent_digit = 1, 
              label_size = 3,
              )+
  scale_x_continuous(expand = expansion(mult = .2))+
  scale_fill_gradient(low="white",high = "#241c1f",limits=c(0,40),name = "Species Number")+
  theme(legend.title = element_text(size=16,family="serif",face="bold"),
        legend.text = element_text(size=15,family="serif",face="bold"),
        plot.subtitle = element_text(size=12,family="serif",face="bold")
        
  )
vp_plt

vp_sk_plt<-ggVennDiagram(vp_sk, 
                      category.names = c("HSCT","GVHD","Chemo"),
                      set_color = c("#BC3C29FF", "#0072B5FF", "#008134FF"),
                      set_size = 6,
                      label = "both", 
                      label_percent_digit = 1, 
                      label_size = 3,
)+
  scale_x_continuous(expand = expansion(mult = .2))+
  scale_fill_gradient(low="white",high = "#241c1f",limits=c(0,40),name = "Species Number")+
  theme(legend.title = element_text(size=16,family="serif",face="bold"),
        legend.text = element_text(size=15,family="serif",face="bold"),
        plot.subtitle = element_text(size=12,family="serif",face="bold")
        
  )
vp_sk_plt
Figure4_venn<-vp_plt|vp_sk_plt+plot_layout(guides = 'collect')
ggsave("Figure0111/Figure4_venn.pdf",Figure4_venn,units="in", width=12, height=6, dpi=600,limitsize = FALSE)



vp_mul_hsct<-as.character(fil_mul_hsct_e$taxa)
vp_mul_sk_hsct<-unique(sk_mul_hsct$Taxa)
vp_mul_gvhd<-as.character(fil_mul_gvhd_e$taxa)
vp_mul_sk_gvhd<-unique(sk_mul_gvhd$Taxa)
vp_mul_chemo<-as.character(fil_mul_chemo_e$taxa)
vp_mul_sk_chemo<-unique(sk_mul_chemo$Taxa)

intersect(intersect(vp_mul_hsct,vp_mul_gvhd),vp_mul_chemo)

vp_mul<-list(
  HSCT=vp_mul_hsct,
  GVHD=vp_mul_gvhd,
  Chemo=vp_mul_chemo
)

vp_mul_sk<-list(
  HSCT=vp_mul_sk_hsct,
  GVHD=vp_mul_sk_gvhd,
  Chemo=vp_mul_sk_chemo
)

# vp_mul_plt<-ggvenn(vp_mul,c("HSCT","GVHD","Chemo"),show_percentage = T,
#                    stroke_color = "white",
#                    fill_color = c("#BC3C29","#0072B5","#008134"),
#                    set_name_color =c("#BC3C29","#0072B5","#008134"))+
#   theme(text = element_text(size=16,family="serif",face="bold"))
# vp_mul_plt
# ggsave("Figure/Figure3_3ls.pdf",vp_mul_plt,units="in", width=3, height=3, dpi=600,limitsize = FALSE)
# 
# vp_mul_sk_plt<-ggvenn(vp_mul_sk,c("HSCT","GVHD","Chemo"),show_percentage = T,
#                       stroke_color = "white",
#                       fill_color = c("#BC3C29","#0072B5","#008134"),
#                       set_name_color =c("#BC3C29","#0072B5","#008134"))+
#   theme(text = element_text(size=16,family="serif",face="bold"))
# vp_mul_sk_plt
# ggsave("Figure/Figure3_3rs.pdf",vp_mul_sk_plt,units="in", width=3, height=3, dpi=600,limitsize = FALSE)

vp_mul_plt<-ggVennDiagram(vp_mul, 
                      category.names = c("HSCT","GVHD","Chemo"),
                      set_color = c("#BC3C29FF", "#0072B5FF", "#008134FF"),
                      set_size = 6,
                      label = "both", 
                      label_percent_digit = 1, 
                      label_size = 3,
)+
  scale_x_continuous(expand = expansion(mult = .2))+
  scale_fill_gradient(low="white",high = "#241c1f",limits=c(0,60),name = "Species Number")+
  theme(legend.title = element_text(size=16,family="serif",face="bold"),
        legend.text = element_text(size=15,family="serif",face="bold"),
        plot.subtitle = element_text(size=12,family="serif",face="bold")
        
  )
vp_mul_plt

vp_mul_sk_plt<-ggVennDiagram(vp_mul_sk, 
                         category.names = c("HSCT","GVHD","Chemo"),
                         set_color = c("#BC3C29FF", "#0072B5FF", "#008134FF"),
                         set_size = 6,
                         label = "both", 
                         label_percent_digit = 1, 
                         label_size = 3,
)+
  scale_x_continuous(expand = expansion(mult = .2))+
  scale_fill_gradient(low="white",high = "#241c1f",limits=c(0,60),name = "Species Number")+
  theme(legend.title = element_text(size=16,family="serif",face="bold"),
        legend.text = element_text(size=15,family="serif",face="bold"),
        plot.subtitle = element_text(size=12,family="serif",face="bold")
        
  )
vp_mul_sk_plt
Figure4s_venn<-vp_mul_plt|vp_mul_sk_plt+plot_layout(guides = 'collect')
ggsave("Figure0111/Figure4s_venn.pdf",Figure4s_venn,units="in", width=12, height=6, dpi=600,limitsize = FALSE)



### Figure 3s ---------------------------------------------------------------
library(ade4)
library(ggtree)

tmp <- unique(c(as.character(fil_hsct_e$taxa),as.character(fil_gvhd_e$taxa),as.character(fil_chemo_e$taxa),
                as.character(fil_mul_hsct_e$taxa),as.character(fil_mul_gvhd_e$taxa),as.character(fil_mul_chemo_e$taxa)))
tmp1<-as.data.frame(order="1",
                    tmp)
write.csv(tmp,"hp_seperate_species.csv")

tree_seperate_species<-read.csv("heatmap/tree_seperate_species.csv")
# hp_merge_species<-read.csv("hp_merge_species_mod.csv")
tree_seperate_species$ID<-factor(tree_seperate_species$ID,levels = tree_seperate_species$ID)


rownames(tree_seperate_species)<-tree_seperate_species$ID
tree_seperate_species<-tree_seperate_species[,-1]
tree_seperate_species$p <-as.factor(tree_seperate_species$p)
tree_seperate_species$c <-as.factor(tree_seperate_species$c)
tree_seperate_species$o <-as.factor(tree_seperate_species$o)
tree_seperate_species$f <-as.factor(tree_seperate_species$f)
tree_seperate_species$g <-as.factor(tree_seperate_species$g)

tmptax1.phy <- taxo2phylog(as.taxo(tree_seperate_species))
print(tmptax1.phy)
tree_sep_plt_name<-ggtree(tmptax1.phy,colour = "grey50")  + 
  geom_tiplab(size=3.0)+
  geom_nodelab(size=3.0)+
  hexpand(.25)#show id out of edge
tree_sep_plt_name

tree_sep_plt_noname<-ggtree(tmptax1.phy,colour = "grey50",size=0.9)  + 
  coord_flip()
tree_sep_plt_noname


hp_sep_species_hsct<-melt(fil_hsct_e,ID=C("taxa","drug")) |> na.omit() |> 
  mutate(Cohort="HSCT",Group="Single drug analysis")
colnames(hp_sep_species_hsct)<-c("Taxa","Meds","Magnitude","Cohort","Group")
pval<-melt(fil_hsct[,c(1,2,9)],ID=C("taxa","drug")) |> na.omit()
hp_sep_species_hsct$`P value`<-pval$value

hp_sep_mul_species_hsct<-melt(fil_mul_hsct_e,ID=C("taxa","drug")) |> na.omit() |> 
  mutate(Cohort="HSCT",Group="Multiple drug analysis")
colnames(hp_sep_mul_species_hsct)<-c("Taxa","Meds","Magnitude","Cohort","Group")
pval<-melt(fil_mul_hsct[,c(1,2,9)],ID=C("taxa","drug")) |> na.omit()
hp_sep_mul_species_hsct$`P value`<-pval$value


hp_sep_species_gvhd<-melt(fil_gvhd_e,ID=C("taxa","drug")) |> na.omit() |> 
  mutate(Cohort="GVHD",Group="Single drug analysis")
colnames(hp_sep_species_gvhd)<-c("Taxa","Meds","Magnitude","Cohort","Group")
pval<-melt(fil_gvhd[,c(1,2,17)],ID=C("taxa","drug")) |> na.omit()
hp_sep_species_gvhd$`P value`<-pval$value

hp_sep_mul_species_gvhd<-melt(fil_mul_gvhd_e,ID=C("taxa","drug")) |> na.omit() |> 
  mutate(Cohort="GVHD",Group="Multiple drug analysis")
colnames(hp_sep_mul_species_gvhd)<-c("Taxa","Meds","Magnitude","Cohort","Group")
pval<-melt(fil_mul_gvhd[,c(1,2,17)],ID=C("taxa","drug")) |> na.omit()
hp_sep_mul_species_gvhd$`P value`<-pval$value


hp_sep_species_chemo<-melt(fil_chemo_e,ID=C("taxa","drug")) |> na.omit() |> 
  mutate(Cohort="Chemo",Group="Single drug analysis")
colnames(hp_sep_species_chemo)<-c("Taxa","Meds","Magnitude","Cohort","Group")
pval<-melt(fil_chemo[,c(1,2,25)],ID=C("taxa","drug")) |> na.omit()
hp_sep_species_chemo$`P value`<-pval$value

hp_sep_mul_species_chemo<-melt(fil_mul_chemo_e,ID=C("taxa","drug")) |> na.omit() |> 
  mutate(Cohort="Chemo",Group="Multiple drug analysis")
colnames(hp_sep_mul_species_chemo)<-c("Taxa","Meds","Magnitude","Cohort","Group")
pval<-melt(fil_mul_chemo[,c(1,2,25)],ID=C("taxa","drug")) |> na.omit()
hp_sep_mul_species_chemo$`P value`<-pval$value

hp_sep_merge_species<-rbind(hp_sep_species_hsct,hp_sep_mul_species_hsct,
                            hp_sep_species_gvhd,hp_sep_mul_species_gvhd,
                            hp_sep_species_chemo,hp_sep_mul_species_chemo) |> 
  filter(Meds!="Total.Body.Irradiation") |> 
  mutate(Influence=ifelse(Magnitude<0,"Inhibition","Promotion"))
write.csv(hp_sep_merge_species,"heatmap/hp_sep_merge_species.csv")


hp_sep_merge_species<-read.csv("heatmap/hp_sep_merge_species_mod.csv")
hp_sep_merge_species$Meds<-factor(hp_sep_merge_species$Meds,levels = levels(test_merge_cohort$Meds))
hp_sep_merge_species$Group<-factor(hp_sep_merge_species$Group,levels = c("Single drug analysis","Multiple drug analysis"))
hp_sep_merge_species$Cohort<-factor(hp_sep_merge_species$Cohort,levels = c("HSCT","GVHD","Chemo"))
hp_sep_merge_species$facet<-paste(hp_sep_merge_species$Cohort,hp_sep_merge_species$Group,sep = " ")
hp_sep_merge_species$facet<-factor(hp_sep_merge_species$facet,levels = c("HSCT Single drug analysis","HSCT Multiple drug analysis",
                                                                         "GVHD Single drug analysis","GVHD Multiple drug analysis",
                                                                         "Chemo Single drug analysis","Chemo Multiple drug analysis"))
hp_sep_merge_species$ID<-factor(hp_sep_merge_species$ID,levels = c("k__Bacteria","Unassigned","o__Streptophyta","g__Methanobrevibacter",
                                                                   "g__Fusobacterium","A. muciniphila","g__Sutterella","g__Bilophila",
                                                                   "g__Campylobacter","H. parainfluenzae","f__Enterobacteriaceae",
                                                                   "g__Shigella","g__Atopobium","C. aerofaciens","E. lenta",
                                                                   "f__Coriobacteriaceae","g__Actinomyces","g__Corynebacterium","R. mucilaginosa",
                                                                   "R. dentocariosa","g__Scardovia","g__Bifidobacterium1","g__Bifidobacterium2",
                                                                   "B. longum","f__Barnesiellaceae","g__Odoribacter","g__Parabacteroides","P. distasonis",
                                                                   "g__Prevotella","P. copri","f__Rikenellaceae","A. finegoldii","A. indistinctus",
                                                                   "A. onderdonkii","A. putredinis","g__Bacteroides1","g__Bacteroides2","B. caccae",
                                                                   "B. fragilis","B. ovatus","B. uniformis","f__Erysipelotrichaceae","C. spiroforme",
                                                                   "g__Coprobacillus","E. dolichum","g__Holdemania","g__Staphylococcus",
                                                                   "g__Turicibacter","o__Lactobacillales","g__Granulicatella","f__Enterococcaceae",
                                                                   "g__Enterococcus1","g__Enterococcus2","E. casseliflavus","g__Lactococcus",
                                                                   "g__Streptococcus1","S. infantis","g__Streptococcus2","f__Lactobacillaceae",
                                                                   "g__Pediococcus","g__Lactobacillus1","g__Lactobacillus2","L. delbrueckii","L. zeae",
                                                                   "L. salivarius","o__Clostridiales1","o__Clostridiales2","f__Christensenellaceae",
                                                                   "g__Anaerofustis","f__Mogibacteriaceae","g__Anaerococcus","g__Finegoldia",
                                                                   "f__Peptostreptococcaceae1","f__Peptostreptococcaceae2","C. bartlettii","f__Clostridiaceae",
                                                                   "g__SMB53","C. celatum","C. paraputrificum","g__Clostridium1","g__Clostridium2",
                                                                   "g__Acidaminococcus","g__Dialister","g__Megasphaera","g__Phascolarctobacterium","V. parvula",
                                                                   "V. dispar","f__Ruminococcaceae1","f__Ruminococcaceae2","g__Anaerotruncus","B. pullicaecorum",
                                                                   "C. methylpentosum","G. formicilis","g__Oscillospira","S. variabile","g__Faecalibacterium",
                                                                   "F. prausnitzii","g__Ruminococcus1","R. bromii","f__Lachnospiraceae1","f__Lachnospiraceae2",
                                                                   "g__Coprococcus","g__Lachnospira","g__Blautia","B. obeum","B. producta","g__Dorea",
                                                                   "D. longicatena","D. formicigenerans","g__Roseburia1","g__Roseburia2","R. faecis",
                                                                   "g__Ruminococcus2","R. gnavus","R. torques","g__Clostridium","C. aldenense","C. citroniae",
                                                                   "C. clostridioforme","C. hathewayi","C. lavalense","C. symbiosum"))
hp_sep_merge_species_plt<-hp_sep_merge_species |> filter(facet!="Chemo Multiple drug analysis")
# hp_sep_merge_species_plt$textcol<-ifelse(hp_sep_merge_species_plt$Taxa %in% intersect(intersect(vp_hsct,vp_gvhd),vp_chemo),"1",
#                                          ifelse(hp_sep_merge_species_plt$Taxa %in% intersect(intersect(vp_mul_hsct,vp_mul_gvhd),vp_mul_chemo),"2","3"))
temp<-ggplot(hp_sep_merge_species_plt,aes(x=ID,y=Meds,color=Influence))+
  geom_point(aes(size=abs(Magnitude),alpha=-P.value))+
  scale_size_area(max_size = 7)+
  # scale_size_continuous(name="Magnitude")+
  scale_alpha_continuous(name="P value",
                         limits = c(-1,0), breaks = c(0, -0.01, -0.1,-0.5, -1),
                         guide = guide_colourbar(nbin = 100, draw.ulim = FALSE, draw.llim = FALSE))+
  scale_color_manual(name="Influence",
                     labels=c("Inhibition","Promotion"),
                     values=c("#8491B4FF","#F39B7FFF")) +
  # scale_x_discrete(name="Taxa influenced in three cohorts",
  #                  aes(color=textcol))+
  theme(panel.background = element_blank(),
        panel.grid =element_blank(),
        axis.line = element_blank(),
        # axis.title.y = element_text(size=14,family="serif",face="bold"),
        axis.title = element_blank(),
        axis.text.y =  element_text(size=6,family="serif",face="bold"),
        axis.text.x =  element_text(size=10,family="serif",face="bold",angle=90,hjust = 1,vjust = 0.5),
        #legend.position = 'none',
        legend.title = element_text(size=16,family="serif",face="bold"),
        legend.text = element_text(size=15,family="serif",face="bold"),
        plot.subtitle = element_text(size=12,family="serif",face="bold"),
        strip.text = element_text(size=10,family="serif",face="bold"),
        strip.background = element_rect(fill = NA),
        title = element_text(size=14,family="serif",face="bold")
  )+
  labs(size="Magnitude")+
  facet_grid(rows = vars(facet), drop=T)+
  guides(colour = guide_legend(override.aes = list(size=5)),
         alpha=guide_legend(override.aes = list(size=5)))#修改图例点大小

temp
write.csv(hp_sep_merge_species_plt,"figure4data.csv")
b_sep<-temp/tree_sep_plt_noname+plot_layout(heights = c(10.5, 1.2))
ggsave("Figure/figure4.tiff",b_sep,units="in", width=20, height=17.6, dpi=600,limitsize = FALSE)                                     

 
# Species meta-analysis ---------------------------------------------------

## single drug analysis ----------------------------------------------------

meta_hsct_e<-read.csv( "F:/数据/document/学校/硕士/数据库/r数据库/gvhd/sankeyplot/sk_hsct.csv",header=T)
meta_hsct_e<-meta_hsct_e[,-7]
meta_hsct_e<-melt(meta_hsct_e,ID=c("Taxa","NodeColor"))
meta_hsct_e<-na.omit(meta_hsct_e)
meta_hsct_e<-meta_hsct_e %>% mutate(Cohort="HSCT") %>% mutate(antibioic=c(rep("no",50),rep("yes",64)))

meta_hsct_s<-read.csv( "F:/数据/document/学校/硕士/数据库/r数据库/gvhd/sankeyplot/sk_hsct_s.csv",header=T)
meta_hsct_s<-meta_hsct_s[,-7]
meta_hsct_s<-melt(meta_hsct_s,ID=c("Taxa","NodeColor"))
meta_hsct_s<-na.omit(meta_hsct_s)

meta_hsct<-cbind(meta_hsct_e,meta_hsct_s$value)
colnames(meta_hsct)<-c("Species","Phylum","Medicine","estimate","Cohort","antibiotic","std")

meta_gvhd_e<-read.csv( "F:/数据/document/学校/硕士/数据库/r数据库/gvhd/sankeyplot/sk_gvhd.csv",header=T)
meta_gvhd_e<-meta_gvhd_e[,-c(8,18)]
meta_gvhd_e<-melt(meta_gvhd_e,ID=c("Taxa","NodeColor"))
meta_gvhd_e<-na.omit(meta_gvhd_e)
meta_gvhd_e<-meta_gvhd_e %>% mutate(Cohort="GVHD") %>% mutate(antibioic=c(rep("no",27),rep("yes",66)))

meta_gvhd_s<-read.csv( "F:/数据/document/学校/硕士/数据库/r数据库/gvhd/sankeyplot/sk_gvhd_s.csv",header=T)
meta_gvhd_s<-meta_gvhd_s[,-c(8,18)]
meta_gvhd_s<-melt(meta_gvhd_s,ID=c("Taxa","NodeColor"))
meta_gvhd_s<-na.omit(meta_gvhd_s)

meta_gvhd<-cbind(meta_gvhd_e,meta_gvhd_s$value)
colnames(meta_gvhd)<-c("Species","Phylum","Medicine","estimate","Cohort","antibiotic","std")

meta_chemo_e<-read.csv( "F:/数据/document/学校/硕士/数据库/r数据库/gvhd/sankeyplot/sk_chemo.csv",header=T)
meta_chemo_e<-melt(meta_chemo_e,ID=c("Taxa","NodeColor"))
meta_chemo_e<-na.omit(meta_chemo_e)
meta_chemo_e<-meta_chemo_e %>% mutate(Cohort="Chemo") %>% mutate(antibioic=c(rep("yes",27)))

meta_chemo_s<-read.csv( "F:/数据/document/学校/硕士/数据库/r数据库/gvhd/sankeyplot/sk_chemo_s.csv",header=T)
meta_chemo_s<-melt(meta_chemo_s,ID=c("Taxa","NodeColor"))
meta_chemo_s<-na.omit(meta_chemo_s)

meta_chemo<-cbind(meta_chemo_e,meta_chemo_s$value)
colnames(meta_chemo)<-c("Species","Phylum","Medicine","estimate","Cohort","antibiotic","std")

meta_plt<-rbind(meta_hsct,meta_gvhd,meta_chemo)
meta_species<-c("V. parvula","R. mucilaginosa","B. caccae","B. ovatus","B. uniformis","C. lavalense",
                "C. spiroforme","R. torques","R. faecis","A. onderdonkii","C. citroniae","R. bromii")
meta_plt<-meta_plt |> filter(Species%in%meta_species) |> arrange(Species)

meta_plt$drug<-meta_plt$Medicine
meta_plt<-meta_plt |> group_by(Species,Medicine) |> #metagen(estimate,std) |> 
  summarise(
    size=length(drug),
    Phylum=Phylum,
    Cohort=Cohort,
    antibiotic=antibiotic,
    estimate=estimate,
    std=std
  )

meta_plt1<-meta_plt |> filter(size==1)
meta_plt2<-meta_plt |> filter(size==2)
meta_species2<-unique(meta_plt2$Species)
meta_plt2formerge<-array(0,dim<-c(length(meta_species2),ncol(meta_plt2)))
colnames(meta_plt2formerge)<-colnames(meta_plt2)
x=0
for (i in meta_species2){
  x=x+1
  tmp<-meta_plt2 |> filter(Species==i)
  meta_plt2formerge[x,1]<-i
  meta_plt2formerge[x,2]<-levels(tmp$Medicine)[tmp$Medicine[1]]
  meta_plt2formerge[x,3]<-tmp$size[1]
  meta_plt2formerge[x,4]<-tmp$Phylum[1]
  meta_plt2formerge[x,5]<-paste(tmp$Cohort[1],tmp$Cohort[2],sep="+")
  meta_plt2formerge[x,6]<-tmp$antibiotic[1]
  het<-metagen(tmp$estimate,tmp$std)
  meta_plt2formerge[x,7]<-het$TE.random
  meta_plt2formerge[x,8]<-het$seTE.random
}

meta_plt2formerge<-as.data.frame(meta_plt2formerge)
meta_plt2formerge$size<-as.numeric(meta_plt2formerge$size)
meta_plt2formerge$estimate<-as.numeric(meta_plt2formerge$estimate)
meta_plt2formerge$std<-as.numeric(meta_plt2formerge$std)
meta_pltmerge<-rbind(meta_plt1,meta_plt2formerge)


### non-antibiotic ----------------------------------------------------------



meta_pltmerge_nonanti<-meta_pltmerge |> filter(antibiotic=="no")
meta_species_nonanti<-unique(meta_pltmerge_nonanti$Species)
meta_species_nonantif<-function(i){
  tmp<-meta_pltmerge_nonanti |> filter(Species==i)
  metaspecies<-array(NA,dim=c(nrow(tmp),14))
  colnames(metaspecies)<-c("Medicine","Number of associations","Cohort","Phylum","TE","lci","uci","k",
                           "pval","I2","Q","pval-h","Species","antibiotic")
  het<-metagen(tmp$estimate,tmp$std,tmp$Medicine)
  metaspecies[,1]<-summary(het)$studlab
  metaspecies[,2]<-tmp$size
  metaspecies[,3]<-tmp$Cohort
  metaspecies[,4]<-tmp$Phylum
  metaspecies[,5]<-summary(het)$TE
  metaspecies[,6]<-summary(het)$lower
  metaspecies[,7]<-summary(het)$upper
  metaspecies[,8]<-NA
  metaspecies[,9]<-summary(het)$pval
  metaspecies[,10]<-NA
  metaspecies[,11]<-NA
  metaspecies[,12]<-NA
  metaspecies[,13]<-i
  metaspecies[,14]<-tmp$antibiotic
  metaspecies<-metaspecies |> as.data.frame() |> arrange(Medicine)
  var<-c("IGV Model(non-anti)",summary(het)$k,"Meta",tmp$Phylum[1],het$TE.random,
         het$lower.random,het$upper.random,summary(het)$k,
         het$pval.random,het$I2,het$Q,pchisq(het$Q,df=2,lower.tail=F),i,"Meta")
  metaspecies<-rbind(metaspecies,var)
  metaspecies$Cohort<-factor(metaspecies$Cohort,levels = c("HSCT","GVHD","Chemo","HSCT+GVHD","HSCT+Chemo","Meta"))
  metaspecies$TE<-as.numeric(metaspecies$TE)
  metaspecies$lci<-as.numeric(metaspecies$lci)
  metaspecies$uci<-as.numeric(metaspecies$uci)
  metaspecies$`Number of associations`<-as.numeric(metaspecies$`Number of associations`)
  metaspecies$Medicine<-as.factor(metaspecies$Medicine)
  metaspecies$ynumber<-c(1:(nrow(tmp)+1))
  return(metaspecies)
  
}
species_meta_nonanti<-lapply(meta_species_nonanti,meta_species_nonantif)
species_meta_nonanti1<-adply(species_meta_nonanti,1)[,-1]

# species_meta_nonanti_plt<- ggplot(species_meta_nonanti1,aes(x=TE,y=reorder(Medicine,-ynumber),color=Cohort))+
#   geom_point(aes(shape=antibiotic,size=`Number of associations`))+
#   scale_shape_manual(name="Antibiotic",
#                      labels=levels(as.factor(species_meta_nonanti1$antibiotic)),
#                      values = c(no=17,Meta=16))+
#   scale_size_area(max_size = 8)+
#   scale_color_manual(
#     name="Cohort",
#     labels=c("HSCT","GVHD","HSCT+GVHD","Meta"),
#     values=c(HSCT="#BC3C29",GVHD="#0072B5",`HSCT+GVHD`="#5E576F",Meta="#8A6757"))+
#   geom_errorbar(aes(xmin=lci,xmax=uci),color = "#666666",width=0.3) +
#   theme(panel.background = element_blank(),
#         panel.grid =element_blank(),
#         axis.line = element_line(colour = "grey50"),
#         axis.title = element_text(size=15,family="serif",face="bold"),
#         axis.text = element_text(size=12,family="serif",face="bold"),
#         #legend.position = 'none',
#         legend.title = element_text(size=16,family="serif",face="bold"),
#         legend.text = element_text(size=15,family="serif",face="bold"),
#         plot.subtitle = element_text(size=12,family="serif",face="bold"),
#         strip.text = element_text(size=12,family="serif",face="bold"),
#         strip.background = element_rect(fill = NA),
#         title = element_text(size=15,family="serif",face="bold")
#   )+
#   labs(x="Regression coefficient",y="Medicine")+
#   geom_vline(xintercept = 0,linetype="dashed")+
#   facet_wrap(~Species,scales = "free")+
#   guides(colour = guide_legend(override.aes = list(size=5)),
#          shape=guide_legend(override.aes = list(size=5)))#修改图例点大小
# species_meta_nonanti_plt
# ggsave("Figure5s1.pdf",species_meta_nonanti_plt,units="in", width=12.75, height=8.0, dpi=600,limitsize = FALSE)



### antibitic ---------------------------------------------------------------


meta_pltmerge_anti<-meta_pltmerge |> filter(antibiotic=="yes")
meta_species_anti<-unique(meta_pltmerge_anti$Species)
meta_species_antif<-function(i){
  tmp<-meta_pltmerge_anti |> filter(Species==i)
  metaspecies<-array(NA,dim=c(nrow(tmp),14))
  colnames(metaspecies)<-c("Medicine","Number of associations","Cohort","Phylum","TE","lci","uci","k",
                           "pval","I2","Q","pval-h","Species","antibiotic")
  het<-metagen(tmp$estimate,tmp$std,tmp$Medicine)
  metaspecies[,1]<-summary(het)$studlab
  metaspecies[,2]<-tmp$size
  metaspecies[,3]<-tmp$Cohort
  metaspecies[,4]<-tmp$Phylum
  metaspecies[,5]<-summary(het)$TE
  metaspecies[,6]<-summary(het)$lower
  metaspecies[,7]<-summary(het)$upper
  metaspecies[,8]<-NA
  metaspecies[,9]<-summary(het)$pval
  metaspecies[,10]<-NA
  metaspecies[,11]<-NA
  metaspecies[,12]<-NA
  metaspecies[,13]<-i
  metaspecies[,14]<-tmp$antibiotic
  metaspecies<-metaspecies |> as.data.frame() |> arrange(Medicine)
  var<-c("IGV Model(anti)",summary(het)$k,"Meta",tmp$Phylum[1],het$TE.random,
         het$lower.random,het$upper.random,summary(het)$k,
         het$pval.random,het$I2,het$Q,pchisq(het$Q,df=2,lower.tail=F),i,"Meta")
  metaspecies<-rbind(metaspecies,var)
  metaspecies$Cohort<-factor(metaspecies$Cohort,levels = c("HSCT","GVHD","Chemo","HSCT+GVHD","HSCT+Chemo","Meta"))
  metaspecies$TE<-as.numeric(metaspecies$TE)
  metaspecies$lci<-as.numeric(metaspecies$lci)
  metaspecies$uci<-as.numeric(metaspecies$uci)
  metaspecies$`Number of associations`<-as.numeric(metaspecies$`Number of associations`)
  metaspecies$Medicine<-as.factor(metaspecies$Medicine)
  metaspecies$ynumber<-c(10:(nrow(tmp)+10))
  return(metaspecies)
  
}
species_meta_anti<-lapply(meta_species_anti,meta_species_antif)
species_meta_anti1<-adply(species_meta_anti,1)[,-1]

# species_meta_anti_plt<- ggplot(species_meta_anti1,aes(x=TE,y=reorder(Medicine,-ynumber),color=Cohort))+
#   geom_point(aes(shape=antibiotic,size=`Number of associations`))+
#   scale_shape_manual(name="Antibiotic",
#                      labels=levels(as.factor(species_meta_anti1$antibiotic)),
#                      values = c(yes=18,Meta=16))+
#   scale_size_area(max_size = 8)+
#   scale_color_manual(
#     name="Cohort",
#     labels=levels(species_meta_anti1$Cohort),
#     values=c(HSCT="#BC3C29",GVHD="#0072B5",Chemo="#E18727",`HSCT+GVHD`="#5E576F",`HSCT+Chemo`="#CF6228",Meta="#8A6757"))+
#   geom_errorbar(aes(xmin=lci,xmax=uci),color = "#666666",width=0.3) +
#   theme(panel.background = element_blank(),
#         panel.grid =element_blank(),
#         axis.line = element_line(colour = "grey50"),
#         axis.title = element_text(size=15,family="serif",face="bold"),
#         axis.text = element_text(size=12,family="serif",face="bold"),
#         #legend.position = 'none',
#         legend.title = element_text(size=16,family="serif",face="bold"),
#         legend.text = element_text(size=15,family="serif",face="bold"),
#         plot.subtitle = element_text(size=12,family="serif",face="bold"),
#         strip.text = element_text(size=12,family="serif",face="bold"),
#         strip.background = element_rect(fill = NA),
#         title = element_text(size=15,family="serif",face="bold")
#   )+
#   labs(x="Regression coefficient",y="Medicine")+
#   geom_vline(xintercept = 0,linetype="dashed")+
#   facet_wrap(~Species,scales = "free")+
#   guides(colour = guide_legend(override.aes = list(size=5)),
#          shape=guide_legend(override.aes = list(size=5)))#修改图例点大小
# species_meta_anti_plt
# ggsave("Figure5s2.pdf",species_meta_anti_plt,units="in", width=19.0, height=8.0, dpi=600,limitsize = FALSE)




### merge analysis ----------------------------------------------------------

species_meta2<-rbind(species_meta_nonanti1,species_meta_anti1)
species_meta2$ynumber<-as.numeric(species_meta2$ynumber)
species_meta2<-species_meta2 |> arrange(Species,ynumber)
x=0
metaspecies<-array(NA,dim=c(length(meta_species),14))
colnames(metaspecies)<-c("Medicine","Number of associations","Cohort","Phylum","TE","lci","uci","k",
                         "pval","I2","Q","pval-h","Species","antibiotic")
for (i in meta_species){
  x=x+1
  tmp<-meta_pltmerge |> filter(Species==i)
  het<-metagen(tmp$estimate,tmp$std,tmp$Medicine)
  metaspecies[x,1]<-"IGV Model"
  metaspecies[x,2]<-summary(het)$k
  metaspecies[x,3]<-"Meta"
  metaspecies[x,4]<-tmp$Phylum[1]
  metaspecies[x,5]<-het$TE.random
  metaspecies[x,6]<-het$lower.random
  metaspecies[x,7]<-het$upper.random
  metaspecies[x,8]<-summary(het)$k
  metaspecies[x,9]<-het$pval.random
  metaspecies[x,10]<-het$I2
  metaspecies[x,11]<-het$Q
  metaspecies[x,12]<-pchisq(het$Q,df=2,lower.tail=F)
  metaspecies[x,13]<-i
  metaspecies[x,14]<-"Meta"

}
metaspecies<-metaspecies |> as.data.frame() |> mutate(ynumber=rep(20,12))
metaspecies$TE<-as.numeric(metaspecies$TE)
metaspecies$lci<-as.numeric(metaspecies$lci)
metaspecies$uci<-as.numeric(metaspecies$uci)
metaspecies$k<-as.numeric(metaspecies$k)
metaspecies$`Number of associations`<-as.numeric(metaspecies$`Number of associations`)
metaspecies$Medicine<-as.factor(metaspecies$Medicine)

species_meta2<-rbind(species_meta_nonanti1,species_meta_anti1,metaspecies)
species_meta2$ynumber<-as.numeric(species_meta2$ynumber)
species_meta2<-species_meta2 |> arrange(Species,ynumber)
species_meta2<-species_meta2[-c(5,75,82),]


### Figure4 -----------------------------------------------------------------


species_meta_plt<- ggplot(species_meta2,aes(x=TE,y=reorder(Medicine,-ynumber),color=Cohort))+
  geom_errorbar(aes(xmin=lci,xmax=uci),color = "#666666",width=0.3) +
  geom_point(aes(shape=antibiotic,size=`Number of associations`))+
  scale_shape_manual(name="Antibiotic",
                     labels=levels(as.factor(species_meta2$antibiotic)),
                     values = c(no=17,yes=18,Meta=16))+
  scale_size_area(max_size = 8.5)+
  scale_color_manual(
    name="Cohort",
    labels=levels(species_meta2$Cohort),
    values=c(HSCT="#BC3C29",GVHD="#0072B5",Chemo="#008134",`HSCT+GVHD`="#6f57b4",`HSCT+Chemo`="#9e9f82",Meta="#241c1f"))+
  theme(panel.background = element_blank(),
        panel.grid =element_blank(),
        axis.line = element_line(colour = "grey50"),
        axis.title = element_text(size=15,family="serif",face="bold"),
        axis.text = element_text(size=12,family="serif",face="bold"),
        #legend.position = 'none',
        legend.title = element_text(size=16,family="serif",face="bold"),
        legend.text = element_text(size=15,family="serif",face="bold"),
        plot.subtitle = element_text(size=12,family="serif",face="bold"),
        strip.text = element_text(size=12,family="serif",face="bold"),
        strip.background = element_rect(fill = NA),
        title = element_text(size=15,family="serif",face="bold")
  )+
  labs(x="Regression coefficient",y="Medicine")+
  geom_vline(xintercept = 0,linetype="dashed")+
  facet_wrap(~Species,scales = "free")+
  guides(colour = guide_legend(override.aes = list(size=5)),
         shape=guide_legend(override.aes = list(size=5)))#修改图例点大小

species_meta_plt

ggsave("Figure/Figure4.pdf",species_meta_plt,units="in", width=20.0, height=10.0, dpi=600,limitsize = FALSE)



## multiple drug analysis --------------------------------------------------

meta_mul_hsct_e<-read.csv( "F:/数据/document/学校/硕士/数据库/r数据库/gvhd/sankeyplot/sk_mul_hsct.csv",header=T)
meta_mul_hsct_e<-meta_mul_hsct_e[,-6]
meta_mul_hsct_e<-melt(meta_mul_hsct_e,ID=c("Taxa","NodeColor"))
meta_mul_hsct_e<-na.omit(meta_mul_hsct_e)
meta_mul_hsct_e<-meta_mul_hsct_e %>% mutate(Cohort="HSCT") %>% mutate(antibioic=c(rep("no",23),rep("yes",56)))

meta_mul_hsct_s<-read.csv( "F:/数据/document/学校/硕士/数据库/r数据库/gvhd/sankeyplot/sk_mul_hsct_s.csv",header=T)
meta_mul_hsct_s<-meta_mul_hsct_s[,-6]
meta_mul_hsct_s<-melt(meta_mul_hsct_s,ID=c("Taxa","NodeColor"))
meta_mul_hsct_s<-na.omit(meta_mul_hsct_s)

meta_mul_hsct<-cbind(meta_mul_hsct_e,meta_mul_hsct_s$value)
colnames(meta_mul_hsct)<-c("Species","Phylum","Medicine","estimate","Cohort","antibiotic","std")

meta_mul_gvhd_e<-read.csv( "F:/数据/document/学校/硕士/数据库/r数据库/gvhd/sankeyplot/sk_mul_gvhd.csv",header=T)
meta_mul_gvhd_e<-meta_mul_gvhd_e[,-c(8,17)]
meta_mul_gvhd_e<-melt(meta_mul_gvhd_e,ID=c("Taxa","NodeColor"))
meta_mul_gvhd_e<-na.omit(meta_mul_gvhd_e)
meta_mul_gvhd_e<-meta_mul_gvhd_e %>% mutate(Cohort="GVHD") %>% mutate(antibioic=c(rep("no",32),rep("yes",60)))

meta_mul_gvhd_s<-read.csv( "F:/数据/document/学校/硕士/数据库/r数据库/gvhd/sankeyplot/sk_mul_gvhd_s.csv",header=T)
meta_mul_gvhd_s<-meta_mul_gvhd_s[,-c(8,17)]
meta_mul_gvhd_s<-melt(meta_mul_gvhd_s,ID=c("Taxa","NodeColor"))
meta_mul_gvhd_s<-na.omit(meta_mul_gvhd_s)

meta_mul_gvhd<-cbind(meta_mul_gvhd_e,meta_mul_gvhd_s$value)
colnames(meta_mul_gvhd)<-c("Species","Phylum","Medicine","estimate","Cohort","antibiotic","std")

meta_mul_chemo_e<-read.csv( "F:/数据/document/学校/硕士/数据库/r数据库/gvhd/sankeyplot/sk_mul_chemo.csv",header=T)
meta_mul_chemo_e<-melt(meta_mul_chemo_e,ID=c("Taxa","NodeColor"))
meta_mul_chemo_e<-na.omit(meta_mul_chemo_e)
meta_mul_chemo_e<-meta_mul_chemo_e %>% mutate(Cohort="Chemo") %>% mutate(antibioic=c(rep("no",5)))

meta_mul_chemo_s<-read.csv( "F:/数据/document/学校/硕士/数据库/r数据库/gvhd/sankeyplot/sk_mul_chemo_s.csv",header=T)
meta_mul_chemo_s<-melt(meta_mul_chemo_s,ID=c("Taxa","NodeColor"))
meta_mul_chemo_s<-na.omit(meta_mul_chemo_s)

meta_mul_chemo<-cbind(meta_mul_chemo_e,meta_mul_chemo_s$value)
colnames(meta_mul_chemo)<-c("Species","Phylum","Medicine","estimate","Cohort","antibiotic","std")

meta_mul_plt<-rbind(meta_mul_hsct,meta_mul_gvhd,meta_mul_chemo)
meta_mul_species<-c("V. parvula","E. dolichum","C. methylpentosum","B. producta")
meta_mul_plt<-meta_mul_plt |> filter(Species%in%meta_mul_species) |> arrange(Species)

meta_mul_plt$drug<-meta_mul_plt$Medicine
meta_mul_plt<-meta_mul_plt |> group_by(Species,Medicine) |> 
  summarise(
    size=length(drug),
    Phylum=Phylum,
    Cohort=Cohort,
    antibiotic=antibiotic,
    estimate=estimate,
    std=std
  )

meta_mul_plt1<-meta_mul_plt |> filter(size==1)
meta_mul_plt2<-meta_mul_plt |> filter(size==2)
meta_mul_species2<-unique(meta_mul_plt2$Species)
meta_mul_plt2formerge<-array(0,dim<-c(length(meta_mul_species2),ncol(meta_mul_plt2)))
colnames(meta_mul_plt2formerge)<-colnames(meta_mul_plt2)
x=0
for (i in meta_mul_species2){
  x=x+1
  tmp<-meta_mul_plt2 |> filter(Species==i)
  meta_mul_plt2formerge[x,1]<-i
  meta_mul_plt2formerge[x,2]<-levels(tmp$Medicine)[tmp$Medicine[1]]
  meta_mul_plt2formerge[x,3]<-tmp$size[1]
  meta_mul_plt2formerge[x,4]<-tmp$Phylum[1]
  meta_mul_plt2formerge[x,5]<-paste(tmp$Cohort[1],tmp$Cohort[2],sep="+")
  meta_mul_plt2formerge[x,6]<-tmp$antibiotic[1]
  het<-metagen(tmp$estimate,tmp$std)
  meta_mul_plt2formerge[x,7]<-het$TE.random
  meta_mul_plt2formerge[x,8]<-het$seTE.random
}

meta_mul_plt2formerge<-as.data.frame(meta_mul_plt2formerge)
meta_mul_plt2formerge$size<-as.numeric(meta_mul_plt2formerge$size)
meta_mul_plt2formerge$estimate<-as.numeric(meta_mul_plt2formerge$estimate)
meta_mul_plt2formerge$std<-as.numeric(meta_mul_plt2formerge$std)
meta_mul_pltmerge<-rbind(meta_mul_plt1,meta_mul_plt2formerge)

### non-antibiotic ----------------------------------------------------------



meta_mul_pltmerge_nonanti<-meta_mul_pltmerge |> filter(antibiotic=="no")
meta_mul_species_nonanti<-unique(meta_mul_pltmerge_nonanti$Species)
meta_mul_species_nonantif<-function(i){
  tmp<-meta_mul_pltmerge_nonanti |> filter(Species==i)
  metaspecies<-array(NA,dim=c(nrow(tmp),14))
  colnames(metaspecies)<-c("Medicine","Number of associations","Cohort","Phylum","TE","lci","uci","k",
                           "pval","I2","Q","pval-h","Species","antibiotic")
  het<-metagen(tmp$estimate,tmp$std,tmp$Medicine)
  metaspecies[,1]<-summary(het)$studlab
  metaspecies[,2]<-tmp$size
  metaspecies[,3]<-tmp$Cohort
  metaspecies[,4]<-tmp$Phylum
  metaspecies[,5]<-summary(het)$TE
  metaspecies[,6]<-summary(het)$lower
  metaspecies[,7]<-summary(het)$upper
  metaspecies[,8]<-NA
  metaspecies[,9]<-summary(het)$pval
  metaspecies[,10]<-NA
  metaspecies[,11]<-NA
  metaspecies[,12]<-NA
  metaspecies[,13]<-i
  metaspecies[,14]<-tmp$antibiotic
  metaspecies<-metaspecies |> as.data.frame() |> arrange(Medicine)
  var<-c("IGV Model(non-anti)",summary(het)$k,"Meta",tmp$Phylum[1],het$TE.random,
         het$lower.random,het$upper.random,summary(het)$k,
         het$pval.random,het$I2,het$Q,pchisq(het$Q,df=2,lower.tail=F),i,"Meta")
  metaspecies<-rbind(metaspecies,var)
  metaspecies$Cohort<-factor(metaspecies$Cohort,levels = c("HSCT","GVHD","Chemo","HSCT+GVHD","HSCT+Chemo","Meta"))
  metaspecies$TE<-as.numeric(metaspecies$TE)
  metaspecies$lci<-as.numeric(metaspecies$lci)
  metaspecies$uci<-as.numeric(metaspecies$uci)
  metaspecies$`Number of associations`<-as.numeric(metaspecies$`Number of associations`)
  metaspecies$Medicine<-as.factor(metaspecies$Medicine)
  metaspecies$ynumber<-c(1:(nrow(tmp)+1))
  return(metaspecies)
  
}
species_meta_mul_nonanti<-lapply(meta_mul_species_nonanti,meta_mul_species_nonantif)
species_meta_mul_nonanti1<-adply(species_meta_mul_nonanti,1)[,-1]




### antibitic ---------------------------------------------------------------


meta_mul_pltmerge_anti<-meta_mul_pltmerge |> filter(antibiotic=="yes")
meta_mul_species_anti<-unique(meta_mul_pltmerge_anti$Species)
meta_mul_species_antif<-function(i){
  tmp<-meta_mul_pltmerge_anti |> filter(Species==i)
  metaspecies<-array(NA,dim=c(nrow(tmp),14))
  colnames(metaspecies)<-c("Medicine","Number of associations","Cohort","Phylum","TE","lci","uci","k",
                           "pval","I2","Q","pval-h","Species","antibiotic")
  het<-metagen(tmp$estimate,tmp$std,tmp$Medicine)
  metaspecies[,1]<-summary(het)$studlab
  metaspecies[,2]<-tmp$size
  metaspecies[,3]<-tmp$Cohort
  metaspecies[,4]<-tmp$Phylum
  metaspecies[,5]<-summary(het)$TE
  metaspecies[,6]<-summary(het)$lower
  metaspecies[,7]<-summary(het)$upper
  metaspecies[,8]<-NA
  metaspecies[,9]<-summary(het)$pval
  metaspecies[,10]<-NA
  metaspecies[,11]<-NA
  metaspecies[,12]<-NA
  metaspecies[,13]<-i
  metaspecies[,14]<-tmp$antibiotic
  metaspecies<-metaspecies |> as.data.frame() |> arrange(Medicine)
  var<-c("IGV Model(anti)",summary(het)$k,"Meta",tmp$Phylum[1],het$TE.random,
         het$lower.random,het$upper.random,summary(het)$k,
         het$pval.random,het$I2,het$Q,pchisq(het$Q,df=2,lower.tail=F),i,"Meta")
  metaspecies<-rbind(metaspecies,var)
  metaspecies$Cohort<-factor(metaspecies$Cohort,levels = c("HSCT","GVHD","Chemo","HSCT+GVHD","HSCT+Chemo","Meta"))
  metaspecies$TE<-as.numeric(metaspecies$TE)
  metaspecies$lci<-as.numeric(metaspecies$lci)
  metaspecies$uci<-as.numeric(metaspecies$uci)
  metaspecies$`Number of associations`<-as.numeric(metaspecies$`Number of associations`)
  metaspecies$Medicine<-as.factor(metaspecies$Medicine)
  metaspecies$ynumber<-c(10:(nrow(tmp)+10))
  return(metaspecies)
  
}
species_meta_mul_anti<-lapply(meta_mul_species_anti,meta_mul_species_antif)
species_meta_mul_anti1<-adply(species_meta_mul_anti,1)[,-1]


### merge analysis ----------------------------------------------------------


# species_meta_mul2<-rbind(species_meta_mul_nonanti1,species_meta_mul_anti1)
# species_meta_mul2$ynumber<-as.numeric(species_meta_mul2$ynumber)
# species_meta_mul2<-species_meta_mul2 |> arrange(Species,ynumber)
x=0
metaspecies_mul<-array(NA,dim=c(length(meta_mul_species),14))
colnames(metaspecies_mul)<-c("Medicine","Number of associations","Cohort","Phylum","TE","lci","uci","k",
                         "pval","I2","Q","pval-h","Species","antibiotic")
for (i in meta_mul_species){
  x=x+1
  tmp<-meta_mul_pltmerge |> filter(Species==i)
  het<-metagen(tmp$estimate,tmp$std,tmp$Medicine)
  metaspecies_mul[x,1]<-"IGV Model"
  metaspecies_mul[x,2]<-summary(het)$k
  metaspecies_mul[x,3]<-"Meta"
  metaspecies_mul[x,4]<-tmp$Phylum[1]
  metaspecies_mul[x,5]<-het$TE.random
  metaspecies_mul[x,6]<-het$lower.random
  metaspecies_mul[x,7]<-het$upper.random
  metaspecies_mul[x,8]<-summary(het)$k
  metaspecies_mul[x,9]<-het$pval.random
  metaspecies_mul[x,10]<-het$I2
  metaspecies_mul[x,11]<-het$Q
  metaspecies_mul[x,12]<-pchisq(het$Q,df=2,lower.tail=F)
  metaspecies_mul[x,13]<-i
  metaspecies_mul[x,14]<-"Meta"
  
}
metaspecies_mul<-metaspecies_mul |> as.data.frame() |> mutate(ynumber=rep(20,4))
metaspecies_mul$TE<-as.numeric(metaspecies_mul$TE)
metaspecies_mul$lci<-as.numeric(metaspecies_mul$lci)
metaspecies_mul$uci<-as.numeric(metaspecies_mul$uci)
metaspecies_mul$k<-as.numeric(metaspecies_mul$k)
metaspecies_mul$`Number of associations`<-as.numeric(metaspecies_mul$`Number of associations`)
metaspecies_mul$Medicine<-as.factor(metaspecies_mul$Medicine)

species_meta_mul2<-rbind(species_meta_mul_nonanti1,species_meta_mul_anti1,metaspecies_mul)
species_meta_mul2$ynumber<-as.numeric(species_meta_mul2$ynumber)
species_meta_mul2<-species_meta_mul2 |> arrange(Species,ynumber)


### Figure4s ----------------------------------------------------------------


species_meta_mul_plt<- ggplot(species_meta_mul2,aes(x=TE,y=reorder(Medicine,-ynumber),color=Cohort))+
  geom_errorbar(aes(xmin=lci,xmax=uci),color = "#666666",width=0.3) +
  geom_point(aes(shape=antibiotic,size=`Number of associations`))+
  scale_shape_manual(name="Antibiotic",
                     labels=levels(as.factor(species_meta2$antibiotic)),
                     values = c(no=17,yes=18,Meta=16))+
  scale_size_area(max_size = 7)+
  scale_color_manual(
    name="Cohort",
    labels=levels(species_meta2$Cohort),
    values=c(HSCT="#BC3C29",GVHD="#0072B5",Chemo="#008134",`HSCT+GVHD`="#6f57b4",`HSCT+Chemo`="#9e9f82",Meta="#241c1f"))+
  theme(panel.background = element_blank(),
        panel.grid =element_blank(),
        axis.line = element_line(colour = "grey50"),
        axis.title = element_text(size=15,family="serif",face="bold"),
        axis.text = element_text(size=12,family="serif",face="bold"),
        #legend.position = 'none',
        legend.title = element_text(size=16,family="serif",face="bold"),
        legend.text = element_text(size=15,family="serif",face="bold"),
        plot.subtitle = element_text(size=12,family="serif",face="bold"),
        strip.text = element_text(size=12,family="serif",face="bold"),
        strip.background = element_rect(fill = NA),
        title = element_text(size=15,family="serif",face="bold")
  )+
  labs(x="Regression coefficient",y="Medicine")+
  geom_vline(xintercept = 0,linetype="dashed")+
  facet_wrap(~Species,scales = "free")+
  guides(colour = guide_legend(override.aes = list(size=5)),
         shape=guide_legend(override.aes = list(size=5)))#修改图例点大小

species_meta_mul_plt

ggsave("Figure/Figure4s.pdf",species_meta_mul_plt,units="in", width=11.5, height=6.6, dpi=600,limitsize = FALSE)

## Figure5 -----------------------------------------------------------------


### barplot 5_1 -------------------------------------------------------------


figure5_1<-ggplot(test_merge_cohort,aes(Meds,Associations,fill=Cohort)) +
  geom_bar(stat="identity", position=position_dodge(),alpha=1.0) +
  scale_fill_manual(values=c( "#241c1f","#CCE6D6","#CCE3F0","#EECEC9"))+
  # scale_fill_manual(values=c( "#D3D2D2","#299555","#66AAD3","#C65949")) +
  theme(panel.background = element_blank(),
        panel.grid =element_blank(),
        axis.line = element_line(colour = "grey50"),
        axis.title = element_text(size=11,family="serif",face="bold"),
        axis.text = element_text(size=11,family="serif",face="bold"),
        #legend.position = 'none',
        legend.title = element_text(size=11,family="serif",face="bold"),
        legend.text = element_text(size=11,family="serif",face="bold"),
        plot.subtitle = element_text(size=11,family="serif",face="bold"),
        strip.text = element_text(size=11,family="serif",face="bold"),
        strip.background = element_rect(fill = NA),
        title = element_text(size=11,family="serif",face="bold")
  )+
  # labs(title = "Single Drug Analysis")+
  facet_grid(rows = vars(Group), drop=T, scales = "free") +
  coord_flip() #+ scale_x_discrete(limits = rev(levels(test_allcohort$meds)))

### heatmap 5_2-----------------------------------------------------------------


library(ade4)
library(ggtree)

name_species<-list(hp_busulfan,hp_cyclophosphamide,hp_amoxicillin,hp_antibiotic,hp_ceftazidime,hp_Ciprofloxacin,
                   hp_fludarabine,hp_melphalan,hp_meropenem,hp_Methotrexate,hp_metronidazole,hp_phenoxymethylpenicillin,
                   hp_piperacillin,hp_Thiotepa,hp_trimethoprim,hp_vancomycin,hp_VP16)
hp_species<-Reduce(full_join,name_species)
hp_species$Group<-"Single drug analysis"
length(setdiff(test_allcohort_1$Meds,hp_species$Drug))
hp_species_bu<-data.frame(
  Taxa=rep(NA,8),
  Drug=setdiff(test_allcohort_1$Meds,hp_species$Drug),
  beta=rep(NA,8),
  FDR=rep(NA,8),
  Influence=rep(NA,8),
  anti=c("yes","no",rep("yes",6)),
  Group=rep("Single drug analysis",8)
)
hp_species<-rbind(hp_species,hp_species_bu)
# hp_species$Drug<-factor(hp_species$Drug,levels = levels(test_merge_cohort$Meds))

name_mul_species<-list(hp_mul_busulfan,hp_mul_cyclophosphamide,hp_mul_amoxicillin,hp_mul_ceftazidime,hp_mul_Ciprofloxacin,
                       hp_mul_fludarabine,hp_mul_meropenem,hp_mul_Methotrexate,hp_mul_phenoxymethylpenicillin,
                       hp_mul_Thiotepa,hp_mul_trimethoprim,hp_mul_vancomycin,hp_mul_VP16)

hp_mul_species<-Reduce(full_join,name_mul_species)
hp_mul_species$Group<-"Multiple drug analysis"
length(setdiff(test_mul_allcohort_1$Meds,hp_mul_species$Drug))
hp_mul_species_bu<-data.frame(
  Taxa=rep(NA,9),
  Drug=setdiff(test_mul_allcohort_1$Meds,hp_mul_species$Drug),
  beta_mul=rep(NA,9),
  FDR_mul=rep(NA,9),
  Influence=rep(NA,9),
  anti=c("no","yes","yes","no","yes","no","yes","yes","yes"),
  Group=rep("Multiple drug analysis",9)
)
hp_mul_species<-rbind(hp_mul_species,hp_mul_species_bu)
# hp_mul_species$Drug<-factor(hp_mul_species$Drug,levels = levels(test_merge_cohort$Meds))

hp_mul_species_bu1<-data.frame(
  Taxa=unique(setdiff(hp_species$Taxa,hp_mul_species$Taxa)),
  Drug=rep(NA,17),
  beta_mul=rep(NA,17),
  FDR_mul=rep(NA,17),
  Influence=rep(NA,17),
  anti=rep(NA,17),
  Group=rep("Multiple drug analysis",17)
)
hp_mul_species<-rbind(hp_mul_species,hp_mul_species_bu1)
colnames(hp_mul_species)[3:4]<-c("beta","FDR")

hp_merge_species<-rbind(hp_species,hp_mul_species)
write.csv(hp_merge_species,"hp_merge_species.csv")

tree_merge_species<-read.csv("tree_merge_species.csv")
hp_merge_species<-read.csv("hp_merge_species_mod.csv")
tree_merge_species$ID<-factor(tree_merge_species$ID,levels = tree_merge_species$ID)



rownames(tree_merge_species)<-tree_merge_species$ID
tree_merge_species<-tree_merge_species[,-1]
tree_merge_species$p <-as.factor(tree_merge_species$p)
tree_merge_species$c <-as.factor(tree_merge_species$c)
tree_merge_species$o <-as.factor(tree_merge_species$o)
tree_merge_species$f <-as.factor(tree_merge_species$f)
tree_merge_species$g <-as.factor(tree_merge_species$g)
tmptax.phy <- taxo2phylog(as.taxo(tree_merge_species))
print(tmptax.phy)
tree_plt_name<-ggtree(tmptax.phy,colour = "grey50")  + 
  geom_tiplab()+
  # geom_nodelab()+
  hexpand(.25)#show id out of edge

tree_plt_noname<-ggtree(tmptax.phy,colour = "grey50",size=0.7)  + 
  coord_flip()

hp_merge_species$Drug<-factor(hp_merge_species$Drug,levels = levels(test_merge_cohort$Meds))
hp_merge_species$Group<-factor(hp_merge_species$Group,,levels = c("Single drug analysis","Multiple drug analysis"))
hp_merge_species$ID<-factor(hp_merge_species$ID,levels = c("f__Enterobacteriaceae","g__Actinomyces","B. breve","g__Parabacteroides",
                                                           "g__Bacteroides","B. ovatus","B. fragilis","B. caccae","C. ramosum",
                                                           "E. dolichum","o__Lactobacillales","g__Streptococcus","g__Enterococcus1",
                                                           "g__Enterococcus2","g__Lactobacillus","L. zeae","o__Clostridiales1",
                                                           "o__Clostridiales2","g__Clostridium","f__Peptostreptococcaceae","V. parvula",
                                                           "V. dispar","f__Ruminococcaceae","F. prausnitzii","g__Oscillospira","g__Ruminococcus1",
                                                           "S. variabile","f__Lachnospiraceae1","f__Lachnospiraceae2","C. hathewayi","g__Coprococcus",
                                                           "g__Dorea","g__Blautia","B. producta","g__Roseburia","R. faecis","R. gnavus","g__Ruminococcus2"))
hp_merge_species$ID[c(128:135,225:233)]<-"g__Ruminococcus2"
hp_merge_species$Drug[234:250]<-"Tobramycin"
temp<-ggplot(hp_merge_species,aes(x=ID,y=Drug,color=Influence))+
  geom_point(aes(size=abs(beta),alpha=-FDR))+
  scale_size_area(max_size = 8)+
  # scale_size_continuous(name="Magnitude")+
  scale_color_manual(name="Influence",
                     labels=c("Inhibition","Promotion"),
                     values=c("#8491B4FF","#F39B7FFF")) +
  scale_alpha_continuous(name="Adjusted P value (FDR)",
                         limits = c(-1,0), breaks = c(0, -0.01, -0.1,-0.5, -1),
                         guide = guide_colourbar(nbin = 100, draw.ulim = FALSE, draw.llim = FALSE))+
  theme(panel.background = element_blank(),
        panel.grid =element_blank(),
        axis.line = element_blank(),
        # axis.title.y = element_text(size=14,family="serif",face="bold"),
        axis.title = element_blank(),
        axis.text.y =  element_text(size=10,family="serif",face="bold"),
        axis.text.x =  element_text(size=12,family="serif",face="bold",angle=90,hjust = 1,vjust = 0.5),
        #legend.position = 'none',
        legend.title = element_text(size=16,family="serif",face="bold"),
        legend.text = element_text(size=15,family="serif",face="bold"),
        plot.subtitle = element_text(size=12,family="serif",face="bold"),
        strip.text = element_text(size=12,family="serif",face="bold"),
        strip.background = element_rect(fill = NA),
        title = element_text(size=14,family="serif",face="bold")
  )+
  labs(size="Magnitude")+
  facet_grid(rows = vars(Group), drop=T, scales = "free")+
  guides(colour = guide_legend(override.aes = list(size=5)),
         alpha=guide_legend(override.aes = list(size=5)))#修改图例点大小

temp
# a<-figure3/plot_spacer()+plot_layout(heights = c(3.45, 1.55))
b<-temp/tree_plt_noname+plot_layout(heights = c(4.2, 0.8))
# m<-(figure3|b)+plot_layout(widths = c(1, 4),guides = 'collect')+ plot_annotation(tag_levels = "A")
ggsave("Figure/Figure5_single.pdf",b,units="in", width=16.0, height=10.0, dpi=600,limitsize = FALSE)


