rm(list=ls())
gc()

####################
## Load libraries ##
install.packages('pheatmap')
install.packages("lmerTest")
install.packages("robustHD")
library(data.table)
library(dplyr)
library(readxl)
library(pheatmap)
library(lmerTest)
library(robustHD)
setwd('/share/inspurStorage/home1/zixing/ABCDR/Project_1')

dataset_win<-readRDS('/share/inspurStorage/home1/zixing/ABCDR/Project_1/dataset_win.RDS')
dataset_baseline = dataset_win[dataset_win$eventname == "baseline_year_1_arm_1",]

########exclude missing data based on baseline model######
exclusion <- read_excel("/share/inspurStorage/home1/zixing/ABCDR/Project_1/xlsx/preprocessing.xlsx", col_names = TRUE);exclusion=t(data.frame(exclusion));ind_exclu<-c()
for(i in 1:length(exclusion['exclusion',])){
  ind_exclu = c(ind_exclu,c(which(names(dataset_baseline)==exclusion['exclusion',][i])))
}
na_data <- exclusion['exclusion',]

dataset_baseline_rmNA=dataset_baseline[complete.cases(dataset_baseline[,ind_exclu]),]


# qc-smri
# 
dataset_baseline_rmNA_smriqc<-dataset_baseline_rmNA[!(dataset_baseline_rmNA$imgincl_t1w_include==0),]

saveRDS(dataset_baseline_rmNA_smriqc, file = "dataset_baseline_rmNA_smriqc.RDS")

dim(dataset_baseline_rmNA_smriqc)
#9838 subjects

subjects <- subset(dataset_baseline_rmNA_smriqc, select=subjectkey)
#cbcl data has data at baseline and 1/2/3 year followup
#mri data only has data at baseline and 2year followup
dataset_baseline_rmNA_smriqc[dataset_baseline_rmNA_smriqc$subjectkey=='NDAR_INV46P89PTE',]$bmi=NA
dataset_baseline_rmNA_smriqc[dataset_baseline_rmNA_smriqc$subjectkey=='NDAR_INV0P34UPZ9',]$bmi=NA
bmilower10<-as.character(dataset_baseline_rmNA_smriqc[dataset_baseline_rmNA_smriqc$bmi<10 & is.na(dataset_baseline_rmNA_smriqc$bmi)==FALSE,]$subjectkey)
dataset_baseline_rmNA_smriqc[dataset_baseline_rmNA_smriqc$subjectkey %in% bmilower10,]$bmi=NA

dataset_baseline_rmNA_smriqc$bmi<-winsor(dataset_baseline_rmNA_smriqc$bmi, trim = 0.02,na.rm=TRUE)

dataset_baseline_followup<-merge(subjects,dataset_win,by="subjectkey")
saveRDS(dataset_baseline_followup, file = "dataset_baseline_followup.RDS")
