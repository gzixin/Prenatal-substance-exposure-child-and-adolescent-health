library(dplyr)
library(optimx)
library(performance)


setwd('/share/inspurStorage/home1/zixing/ABCDR/Project_1')
project1_dataset_baseline_rmNA_win<-readRDS("dataset_baseline_rmNA_smriqc.RDS")

########### add PCs ##################
setwd('/share/inspurStorage/home1/zixing/ABCDR/Project_1/gwas/scz2022')

pcs <- read.table("ABCD.eigenvec", header=F)

colnames(pcs) <- c("FID", "src_subject_id", paste0("PC",1:15))
project1_dataset_baseline_rmNA_phillip_win_prs <- left_join(project1_dataset_baseline_rmNA_win, pcs, by="src_subject_id")


# setwd('/home1/zixing/ABCDR/Project_1/gwas/adhd2019')
# 
# pcs1 <- read.table("ABCD.eigenvec", header=F)
# 
# # colnames(pcs1) <- c("FID", "src_subject_id", paste0("PCs",1:6))
# project1_dataset_baseline_rm_NA_phillip_win_prs2norm_beforemeancentered <- left_join(project1_dataset_baseline_rm_NA_phillip_win_prs2norm_beforemeancentered, pcs1[1:8], by="src_subject_id")
############ add prs ###########
setwd_list <- c('/home1/zixing/ABCDR/Project_1/gwas/scz2022','/home1/zixing/ABCDR/Project_1/gwas/adhd2019','/home1/zixing/ABCDR/Project_1/gwas/cud','/home1/zixing/ABCDR/Project_1/gwas/mdd2018','/home1/zixing/ABCDR/Project_1/gwas/alcdep_aug2018','/home1/zixing/ABCDR/Project_1/gwas/edu2022')
p.threshold <- c(0.001,0.05,0.1,0.2,0.3,0.4,0.5)


for(i in setwd_list){
  setwd(i)
  # Now merge the files

  for(j in p.threshold){
    prs <- read.table(paste0("ABCD.",j,".profile"), header=T)

    if(i==setwd_list[1]) {
      colnames(prs) <- c("FID","src_subject_id","PHENO","CNT","CNT2",paste0("SCZ_SCORE_",j))
      project1_dataset_baseline_rmNA_phillip_win_prs <- left_join(project1_dataset_baseline_rmNA_phillip_win_prs, prs[,c("src_subject_id",paste0("SCZ_SCORE_",j))], by="src_subject_id")
    }
    if(i==setwd_list[2]) {
      colnames(prs) <- c("FID","src_subject_id","PHENO","CNT","CNT2",paste0("ADHD_SCORE_",j))
      project1_dataset_baseline_rmNA_phillip_win_prs <- left_join(project1_dataset_baseline_rmNA_phillip_win_prs, prs[,c("src_subject_id",paste0("ADHD_SCORE_",j))], by="src_subject_id")
    }
    if(i==setwd_list[3]) {
      colnames(prs) <- c("FID","src_subject_id","PHENO","CNT","CNT2",paste0("CUD_SCORE_",j))
      project1_dataset_baseline_rmNA_phillip_win_prs <- left_join(project1_dataset_baseline_rmNA_phillip_win_prs, prs[,c("src_subject_id",paste0("CUD_SCORE_",j))], by="src_subject_id")
    }
    if(i==setwd_list[4]) {
      colnames(prs) <- c("FID","src_subject_id","PHENO","CNT","CNT2",paste0("MDD_SCORE_",j))
      project1_dataset_baseline_rmNA_phillip_win_prs <- left_join(project1_dataset_baseline_rmNA_phillip_win_prs, prs[,c("src_subject_id",paste0("MDD_SCORE_",j))], by="src_subject_id")
    }
    if(i==setwd_list[5]) {
      colnames(prs) <- c("FID","src_subject_id","PHENO","CNT","CNT2",paste0("ALCDEP_SCORE_",j))
      project1_dataset_baseline_rmNA_phillip_win_prs <- left_join(project1_dataset_baseline_rmNA_phillip_win_prs, prs[,c("src_subject_id",paste0("ALCDEP_SCORE_",j))], by="src_subject_id")
    }
    if(i==setwd_list[6]) {
      colnames(prs) <- c("FID","src_subject_id","PHENO","CNT","CNT2",paste0("EDU_SCORE_",j))
      project1_dataset_baseline_rmNA_phillip_win_prs <- left_join(project1_dataset_baseline_rmNA_phillip_win_prs, prs[,c("src_subject_id",paste0("EDU_SCORE_",j))], by="src_subject_id")
    }
  }
}

SCZ_name<-paste0("SCZ_SCORE_",p.threshold)
ADHD_name<-paste0("ADHD_SCORE_",p.threshold)
CUD_name<-paste0("CUD_SCORE_",p.threshold)
MDD_name<-paste0("MDD_SCORE_",p.threshold)
ALCDEP_name<-paste0("ALCDEP_SCORE_",p.threshold)
EDU_name<-paste0("EDU_SCORE_",p.threshold)

project1_dataset_baseline_rmNA_phillip_win_prs$adhd_avg<-rowMeans(project1_dataset_baseline_rmNA_phillip_win_prs[ADHD_name],na.rm=TRUE)
project1_dataset_baseline_rmNA_phillip_win_prs$scz_avg<-rowMeans(project1_dataset_baseline_rmNA_phillip_win_prs[SCZ_name],na.rm=TRUE)
project1_dataset_baseline_rmNA_phillip_win_prs$mdd_avg<-rowMeans(project1_dataset_baseline_rmNA_phillip_win_prs[MDD_name],na.rm=TRUE)

project1_dataset_baseline_rmNA_phillip_win_prs$EA_avg<-rowMeans(project1_dataset_baseline_rmNA_phillip_win_prs[EDU_name],na.rm=TRUE)
project1_dataset_baseline_rmNA_phillip_win_prs$cud_avg<-rowMeans(project1_dataset_baseline_rmNA_phillip_win_prs[CUD_name],na.rm=TRUE)
project1_dataset_baseline_rmNA_phillip_win_prs$alcdep_avg<-rowMeans(project1_dataset_baseline_rmNA_phillip_win_prs[ALCDEP_name],na.rm=TRUE)



name_cbcl_6scale<-c('CBCL.social.problems.score','CBCL.thought.problems','CBCL.attention.problems.score','CBCL.internalizing.factors','CBCL.externalizing.factors','CBCL.total.problems.score')


setwd('/share/inspurStorage/home1/zixing/ABCDR/Project_1')
saveRDS(project1_dataset_baseline_rmNA_phillip_win_prs, file = "project1_dataset_baseline_rmNA_phillip_win_prs.RDS")
