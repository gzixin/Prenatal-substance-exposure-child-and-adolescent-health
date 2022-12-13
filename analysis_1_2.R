#find the assoiations between exposures and confounding variables
rm(list=ls())
gc()

####################
## Load libraries ##
install.packages('pheatmap')
install.packages("lmerTest")
install.packages("robustHD")
install.packages("fmsb")
library(data.table)
library(dplyr)
library(readxl)
library(pheatmap)
library(lmerTest)
library(robustHD)
library(RColorBrewer)
library(MuMIn)
library(lme4)
library(fmsb)


setwd('/share/inspurStorage/home1/zixing/ABCDR/Project_1')
project1_dataset_baseline_rm_NA_phillip_win_prs<-readRDS('project1_dataset_baseline_rmNA_phillip_win_prs.RDS')

project1_dataset_baseline_rm_NA_phillip_win_prs2<-project1_dataset_baseline_rm_NA_phillip_win_prs
#define Min-Max normalization function
min_max_norm <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

# center with 'apply()'
center_apply <- function(x) {
  apply(x, 2, function(y) y - mean(y))
}


#危险
name_pes_risky<-c(name_pes_f[c(3,4,5,9)],name_pes_p[c(4:6,8:20)],name_pes_s[c(1:4,10,13)])
name_pes_protective<-name_pes[!name_pes %in% c(name_pes_risky,'obsteric_complication')]

pesp_risky<-name_pes_p[name_pes_p %in% name_pes_risky]
pesp_protective<-name_pes_p[name_pes_p %in% name_pes_protective]

pesf_risky<-name_pes_f[name_pes_f %in% name_pes_risky]
pesf_protective<-name_pes_f[name_pes_f %in% name_pes_protective]

pess_risky<-name_pes_s[name_pes_s %in% name_pes_risky]
pess_protective<-name_pes_s[name_pes_s %in% name_pes_protective]


name_pes_continuous<- name_pes[name_pes %in% continuous_v]
for (i in name_pes_continuous ){
  project1_dataset_baseline_rm_NA_phillip_win_prs2[,i]<-min_max_norm(project1_dataset_baseline_rm_NA_phillip_win_prs2[,i])
}
otr_substance<-c("coc_crack","her_morph","oxycont")
namepesdmmy<-name_pes[!name_pes %in% name_pes_continuous] 
for (i in namepesdmmy[!namepesdmmy %in% otr_substance] ){
  project1_dataset_baseline_rm_NA_phillip_win_prs2[,i]<-as.numeric(as.character(project1_dataset_baseline_rm_NA_phillip_win_prs2[,i]))
}
for (i in otr_substance ){
  project1_dataset_baseline_rm_NA_phillip_win_prs2[,i]<-as.numeric(project1_dataset_baseline_rm_NA_phillip_win_prs2[,i])
  project1_dataset_baseline_rm_NA_phillip_win_prs2[,i]<-ifelse(project1_dataset_baseline_rm_NA_phillip_win_prs2[,i]==1,0,1)
}

#special for income and education and urbanicity
for (i in c('demo_comb_income_v2','demo_prnt_prtnr_ed_v2','reshist_addr1_urban_area') ){
  project1_dataset_baseline_rm_NA_phillip_win_prs2[,i]<-min_max_norm(project1_dataset_baseline_rm_NA_phillip_win_prs2[,i])
}

for (i in name_genetic[31:36] ){
  project1_dataset_baseline_rm_NA_phillip_win_prs2[,i]<-min_max_norm(project1_dataset_baseline_rm_NA_phillip_win_prs2[,i])
}
for (i in name_genetic[1:30] ){
  project1_dataset_baseline_rm_NA_phillip_win_prs2[,i]<-as.numeric(as.character(project1_dataset_baseline_rm_NA_phillip_win_prs2[,i]))
}


project1_dataset_baseline_rm_NA_phillip_win_prs2$score_pes_risky<-rowSums(project1_dataset_baseline_rm_NA_phillip_win_prs2[name_pes_risky],na.rm=TRUE)
project1_dataset_baseline_rm_NA_phillip_win_prs2$score_pes_protective<-rowSums(project1_dataset_baseline_rm_NA_phillip_win_prs2[name_pes_protective],na.rm=TRUE)
project1_dataset_baseline_rm_NA_phillip_win_prs2$sumscore_pes<-project1_dataset_baseline_rm_NA_phillip_win_prs2$score_pes_risky-project1_dataset_baseline_rm_NA_phillip_win_prs2$score_pes_protective

project1_dataset_baseline_rm_NA_phillip_win_prs2$score_pesp_risky<-rowSums(project1_dataset_baseline_rm_NA_phillip_win_prs2[pesp_risky],na.rm=TRUE)
project1_dataset_baseline_rm_NA_phillip_win_prs2$score_pesp_protective<-rowSums(project1_dataset_baseline_rm_NA_phillip_win_prs2[pesp_protective],na.rm=TRUE)
project1_dataset_baseline_rm_NA_phillip_win_prs2$sumscore_pesp<-project1_dataset_baseline_rm_NA_phillip_win_prs2$score_pesp_risky-project1_dataset_baseline_rm_NA_phillip_win_prs2$score_pesp_protective


project1_dataset_baseline_rm_NA_phillip_win_prs2$score_pesf_risky<-rowSums(project1_dataset_baseline_rm_NA_phillip_win_prs2[pesf_risky],na.rm=TRUE)
project1_dataset_baseline_rm_NA_phillip_win_prs2$score_pesf_protective<-rowSums(project1_dataset_baseline_rm_NA_phillip_win_prs2[pesf_protective],na.rm=TRUE)
project1_dataset_baseline_rm_NA_phillip_win_prs2$sumscore_pesf<-project1_dataset_baseline_rm_NA_phillip_win_prs2$score_pesf_risky-project1_dataset_baseline_rm_NA_phillip_win_prs2$score_pesf_protective


project1_dataset_baseline_rm_NA_phillip_win_prs2$score_pess_risky<-rowSums(project1_dataset_baseline_rm_NA_phillip_win_prs2[pess_risky],na.rm=TRUE)
project1_dataset_baseline_rm_NA_phillip_win_prs2$score_pess_protective<-rowSums(project1_dataset_baseline_rm_NA_phillip_win_prs2[pess_protective],na.rm=TRUE)
project1_dataset_baseline_rm_NA_phillip_win_prs2$sumscore_pess<-project1_dataset_baseline_rm_NA_phillip_win_prs2$score_pess_risky-project1_dataset_baseline_rm_NA_phillip_win_prs2$score_pess_protective

project1_dataset_baseline_rm_NA_phillip_win_prs2$newscore_pesf_risky<-rowSums(project1_dataset_baseline_rm_NA_phillip_win_prs2[pesf_risky[3:4]],na.rm=TRUE)
project1_dataset_baseline_rm_NA_phillip_win_prs2$newsumscore_pesf<-project1_dataset_baseline_rm_NA_phillip_win_prs2$newscore_pesf_risky-project1_dataset_baseline_rm_NA_phillip_win_prs2$score_pesf_protective

pessum<-c('sumscore_pes','sumscore_pesp','sumscore_pesf','sumscore_pess')

project1_dataset_baseline_rm_NA_phillip_win_prs2$EA_avg<- -project1_dataset_baseline_rm_NA_phillip_win_prs2$EA_avg
project1_dataset_baseline_rm_NA_phillip_win_prs2$score_genetic<-rowSums(project1_dataset_baseline_rm_NA_phillip_win_prs2[name_genetic],na.rm=TRUE) 

project1_dataset_baseline_rm_NA_phillip_win_prs2$score_fam<-rowSums(project1_dataset_baseline_rm_NA_phillip_win_prs2[name_genetic[c(1:30)]],na.rm=TRUE) 

project1_dataset_baseline_rm_NA_phillip_win_prs2$famhx_ss_alc_p<-ifelse(rowSums(project1_dataset_baseline_rm_NA_phillip_win_prs2[name_genetic[c(1:3)]],na.rm=TRUE)>1,1,0)
project1_dataset_baseline_rm_NA_phillip_win_prs2$famhx_ss_dg_p<-ifelse(rowSums(project1_dataset_baseline_rm_NA_phillip_win_prs2[name_genetic[c(4:6)]],na.rm=TRUE)>1,1,0)
project1_dataset_baseline_rm_NA_phillip_win_prs2$famhx_ss_dprs_p<-ifelse(rowSums(project1_dataset_baseline_rm_NA_phillip_win_prs2[name_genetic[c(7:9)]],na.rm=TRUE)>1,1,0)
project1_dataset_baseline_rm_NA_phillip_win_prs2$famhx_ss_ma_p<-ifelse(rowSums(project1_dataset_baseline_rm_NA_phillip_win_prs2[name_genetic[c(10:12)]],na.rm=TRUE)>1,1,0)
project1_dataset_baseline_rm_NA_phillip_win_prs2$famhx_ss_vs_p<-ifelse(rowSums(project1_dataset_baseline_rm_NA_phillip_win_prs2[name_genetic[c(13:15)]],na.rm=TRUE)>1,1,0)
project1_dataset_baseline_rm_NA_phillip_win_prs2$famhx_ss_trb_p<-ifelse(rowSums(project1_dataset_baseline_rm_NA_phillip_win_prs2[name_genetic[c(16:18)]],na.rm=TRUE)>1,1,0)
project1_dataset_baseline_rm_NA_phillip_win_prs2$famhx_ss_nrv_p<-ifelse(rowSums(project1_dataset_baseline_rm_NA_phillip_win_prs2[name_genetic[c(19:21)]],na.rm=TRUE)>1,1,0)
project1_dataset_baseline_rm_NA_phillip_win_prs2$famhx_ss_prf_p<-ifelse(rowSums(project1_dataset_baseline_rm_NA_phillip_win_prs2[name_genetic[c(22:24)]],na.rm=TRUE)>1,1,0)
project1_dataset_baseline_rm_NA_phillip_win_prs2$famhx_ss_hspd_p<-ifelse(rowSums(project1_dataset_baseline_rm_NA_phillip_win_prs2[name_genetic[c(25:27)]],na.rm=TRUE)>1,1,0)
project1_dataset_baseline_rm_NA_phillip_win_prs2$famhx_ss_scd_p<-ifelse(rowSums(project1_dataset_baseline_rm_NA_phillip_win_prs2[name_genetic[c(28:30)]],na.rm=TRUE)>1,1,0)


name_newfam<-c('famhx_ss_alc_p','famhx_ss_dg_p','famhx_ss_dprs_p','famhx_ss_ma_p','famhx_ss_vs_p','famhx_ss_trb_p','famhx_ss_nrv_p','famhx_ss_prf_p','famhx_ss_hspd_p','famhx_ss_scd_p')
project1_dataset_baseline_rm_NA_phillip_win_prs2$newscore_fam<-rowSums(project1_dataset_baseline_rm_NA_phillip_win_prs2[name_newfam],na.rm=TRUE) 

project1_dataset_baseline_rm_NA_phillip_win_prs2$newscore_genetic<-rowSums(project1_dataset_baseline_rm_NA_phillip_win_prs2[c(name_newfam,name_genetic[31:36])],na.rm=TRUE) 

EGsum<-c(pessum,'score_genetic','score_fam','newscore_genetic','newsumscore_pesf')

project1_dataset_baseline_rm_NA_phillip_win_prs2[,EGsum] <-center_apply(project1_dataset_baseline_rm_NA_phillip_win_prs2[,EGsum])

# project1_dataset_baseline_rm_NA_phillip_win_prs2[,EGsum] <-center_apply(project1_dataset_baseline_rm_NA_phillip_win_prs2[,EGsum])
# 
# project1_dataset_baseline_rm_NA_phillip_win_prs2$newsumscore_pesf<-center_apply(project1_dataset_baseline_rm_NA_phillip_win_prs2[,EGsum])
levels(project1_dataset_baseline_rm_NA_phillip_win_prs2[,'caffeine'])<-c("FALSE", "TRUE")

############ scale ###################

nonmri_outcomes_v <- nonmri_outcomes_v[c(-32:-40,-44,-46:-48,-50:-53,-56,-84:-86)]

project1_dataset_baseline_rm_NA_phillip_win_prs2[,c(name_area_cdk,name_vol_cdk,name_thickness_cdk,'smri_vol_scs_intracranialv')]<-lapply(project1_dataset_baseline_rm_NA_phillip_win_prs2[,c(name_area_cdk,name_vol_cdk,name_thickness_cdk,'smri_vol_scs_intracranialv')],scale, center = FALSE)

saveRDS(project1_dataset_baseline_rm_NA_phillip_win_prs2, file = "project1_dataset_baseline_rm_NA_phillip_win_prs2norm_beforemeancentered.RDS")

####################### load ################

project1_dataset_baseline_rm_NA_phillip_win_prs2norm_beforemeancentered<-readRDS('project1_dataset_baseline_rm_NA_phillip_win_prs2norm_beforemeancentered.RDS')

################# statistical analysis #####################
library('optimx')

analysis1_2_outcome<-rep()
analysis1_2_exposure<-rep()
analysis1_2_out_beta<-rep()
analysis1_2_out_tvalue <- rep()
analysis1_2_out_pvalue <- rep()
analysis1_2_out_df <- rep()
analysis1_2_stdCoef<-rep()
analysis1_2_CI_low<-rep()
analysis1_2_CI_high<-rep()

analysis1_2_alcohol_exposure<-rep()
analysis1_2_alcohol_out_beta<-rep()
analysis1_2_alcohol_out_tvalue <- rep()
analysis1_2_alcohol_out_pvalue <- rep()
analysis1_2_alcohol_out_df <- rep()
analysis1_2_alcohol_stdCoef<-rep()
analysis1_2_alcohol_CI_low<-rep()
analysis1_2_alcohol_CI_high<-rep()

analysis1_2_tobacco_exposure<-rep()
analysis1_2_tobacco_out_beta<-rep()
analysis1_2_tobacco_out_tvalue <- rep()
analysis1_2_tobacco_out_pvalue <- rep()
analysis1_2_tobacco_out_df <- rep()
analysis1_2_tobacco_stdCoef<-rep()
analysis1_2_tobacco_CI_low<-rep()
analysis1_2_tobacco_CI_high<-rep()

analysis1_2_marijuana_exposure<-rep()
analysis1_2_marijuana_out_beta<-rep()
analysis1_2_marijuana_out_tvalue <- rep()
analysis1_2_marijuana_out_pvalue <- rep()
analysis1_2_marijuana_out_df <- rep()
analysis1_2_marijuana_stdCoef<-rep()
analysis1_2_marijuana_CI_low<-rep()
analysis1_2_marijuana_CI_high<-rep()

analysis1_2_caffeine_exposure<-rep()
analysis1_2_caffeine_out_beta<-rep()
analysis1_2_caffeine_out_tvalue <- rep()
analysis1_2_caffeine_out_pvalue <- rep()
analysis1_2_caffeine_out_df <- rep()
analysis1_2_caffeine_stdCoef<-rep()
analysis1_2_caffeine_CI_low<-rep()
analysis1_2_caffeine_CI_high<-rep()

analysis1_2_r2m<-rep()
analysis1_2_r2c<-rep()
analysis1_2_alcohol_etasq<-rep()
analysis1_2_caffeine_etasq<-rep()
analysis1_2_tobacco_etasq<-rep()
analysis1_2_marijuana_etasq<-rep()
analysis1_2_number=1

for(i in EGsum[c(2,8,4,7)]){ 
  project1_dataset_baseline_rm_NA_phillip_win_prs2norm_beforemeancentered$y = project1_dataset_baseline_rm_NA_phillip_win_prs2norm_beforemeancentered[,i]
  out <- formula(paste( 'y', "~",paste(c(exposure_name,name_fixed,'(1 | site_id_l/rel_family_id)'),collapse="+")))
  lmer.fit <- lmer(out,data = project1_dataset_baseline_rm_NA_phillip_win_prs2norm_beforemeancentered)
  analysis1_2=summary(lmer.fit)$coefficients
  ci<-confint(lmer.fit,method='Wald')
  
  analysis1_2_alcohol_exposure[analysis1_2_number]<-'alcohol'
  analysis1_2_alcohol_out_beta[analysis1_2_number]<-analysis1_2['alcoholTRUE',1]
  analysis1_2_alcohol_out_tvalue[analysis1_2_number] <- analysis1_2['alcoholTRUE',4]
  analysis1_2_alcohol_out_pvalue[analysis1_2_number] <- analysis1_2['alcoholTRUE',5]
  analysis1_2_alcohol_out_df[analysis1_2_number] <- analysis1_2['alcoholTRUE',3]
  analysis1_2_alcohol_stdCoef[analysis1_2_number]<-stdCoef.merMod(lmer.fit)['alcoholTRUE',1]
  analysis1_2_alcohol_CI_low[analysis1_2_number]<-as.numeric(ci['alcoholTRUE',1])
  analysis1_2_alcohol_CI_high[analysis1_2_number]<-as.numeric(ci['alcoholTRUE',2])
  
  analysis1_2_caffeine_exposure[analysis1_2_number]<-'caffeine'
  analysis1_2_caffeine_out_beta[analysis1_2_number]<-analysis1_2['caffeineTRUE',1]
  analysis1_2_caffeine_out_tvalue[analysis1_2_number] <- analysis1_2['caffeineTRUE',4]
  analysis1_2_caffeine_out_pvalue[analysis1_2_number] <- analysis1_2['caffeineTRUE',5]
  analysis1_2_caffeine_out_df[analysis1_2_number] <- analysis1_2['caffeineTRUE',3]
  analysis1_2_caffeine_stdCoef[analysis1_2_number]<-stdCoef.merMod(lmer.fit)['caffeineTRUE',1]
  analysis1_2_caffeine_CI_low[analysis1_2_number]<-as.numeric(ci['caffeineTRUE',1])
  analysis1_2_caffeine_CI_high[analysis1_2_number]<-as.numeric(ci['caffeineTRUE',2])
  
  analysis1_2_tobacco_exposure[analysis1_2_number]<-'tobacco'
  analysis1_2_tobacco_out_beta[analysis1_2_number]<-analysis1_2['tobaccoTRUE',1]
  analysis1_2_tobacco_out_tvalue[analysis1_2_number] <- analysis1_2['tobaccoTRUE',4]
  analysis1_2_tobacco_out_pvalue[analysis1_2_number] <- analysis1_2['tobaccoTRUE',5]
  analysis1_2_tobacco_out_df[analysis1_2_number] <- analysis1_2['tobaccoTRUE',3]
  analysis1_2_tobacco_stdCoef[analysis1_2_number]<-stdCoef.merMod(lmer.fit)['tobaccoTRUE',1]
  analysis1_2_tobacco_CI_low[analysis1_2_number]<-as.numeric(ci['tobaccoTRUE',1])
  analysis1_2_tobacco_CI_high[analysis1_2_number]<-as.numeric(ci['tobaccoTRUE',2])
  
  analysis1_2_marijuana_exposure[analysis1_2_number]<-'marijuana'
  analysis1_2_marijuana_out_beta[analysis1_2_number]<-analysis1_2['marijuanaTRUE',1]
  analysis1_2_marijuana_out_tvalue[analysis1_2_number] <- analysis1_2['marijuanaTRUE',4]
  analysis1_2_marijuana_out_pvalue[analysis1_2_number] <- analysis1_2['marijuanaTRUE',5]
  analysis1_2_marijuana_out_df[analysis1_2_number] <- analysis1_2['marijuanaTRUE',3]
  analysis1_2_marijuana_stdCoef[analysis1_2_number]<-stdCoef.merMod(lmer.fit)['marijuanaTRUE',1]
  analysis1_2_marijuana_CI_low[analysis1_2_number]<-as.numeric(ci['marijuanaTRUE',1])
  analysis1_2_marijuana_CI_high[analysis1_2_number]<-as.numeric(ci['marijuanaTRUE',2])
  
  analysis1_2_r2c[analysis1_2_number] = as.numeric(performance::r2(lmer.fit)[1])
  analysis1_2_r2m[analysis1_2_number] = as.numeric(performance::r2(lmer.fit)[2])
  analysis1_2_alcohol_etasq[analysis1_2_number] = as.numeric(effectsize::eta_squared(lmer.fit)['1',2])
  analysis1_2_tobacco_etasq[analysis1_2_number] = as.numeric(effectsize::eta_squared(lmer.fit)['2',2])
  analysis1_2_marijuana_etasq[analysis1_2_number] = as.numeric(effectsize::eta_squared(lmer.fit)['3',2])
  analysis1_2_caffeine_etasq[analysis1_2_number] = as.numeric(effectsize::eta_squared(lmer.fit)['4',2])
  
  analysis1_2_outcome[analysis1_2_number] = i
  analysis1_2_number = analysis1_2_number + 1

}
results_analysis1_2_alcohol<-data.frame(analysis1_2_outcome,analysis1_2_alcohol_exposure,analysis1_2_alcohol_out_beta,analysis1_2_alcohol_out_tvalue,analysis1_2_alcohol_out_pvalue,analysis1_2_alcohol_out_df,analysis1_2_alcohol_stdCoef,analysis1_2_alcohol_CI_low,analysis1_2_alcohol_etasq,analysis1_2_alcohol_CI_high,analysis1_2_r2c,analysis1_2_r2m)
results_analysis1_2_tobacco<-data.frame(analysis1_2_outcome,analysis1_2_tobacco_exposure,analysis1_2_tobacco_out_beta,analysis1_2_tobacco_out_tvalue,analysis1_2_tobacco_out_pvalue,analysis1_2_tobacco_out_df,analysis1_2_tobacco_stdCoef,analysis1_2_tobacco_CI_low,analysis1_2_tobacco_etasq,analysis1_2_tobacco_CI_high,analysis1_2_r2c,analysis1_2_r2m)
results_analysis1_2_marijuana<-data.frame(analysis1_2_outcome,analysis1_2_marijuana_exposure,analysis1_2_marijuana_out_beta,analysis1_2_marijuana_out_tvalue,analysis1_2_marijuana_out_pvalue,analysis1_2_marijuana_out_df,analysis1_2_marijuana_stdCoef,analysis1_2_marijuana_CI_low,analysis1_2_marijuana_etasq,analysis1_2_marijuana_CI_high,analysis1_2_r2c,analysis1_2_r2m)
results_analysis1_2_caffeine<-data.frame(analysis1_2_outcome,analysis1_2_caffeine_exposure,analysis1_2_caffeine_out_beta,analysis1_2_caffeine_out_tvalue,analysis1_2_caffeine_out_pvalue,analysis1_2_caffeine_out_df,analysis1_2_caffeine_stdCoef,analysis1_2_caffeine_CI_low,analysis1_2_caffeine_etasq,analysis1_2_caffeine_CI_high,analysis1_2_r2c,analysis1_2_r2m)


colnames(results_analysis1_2_alcohol)<-c('analysis1_2_outcome','analysis1_2_exposure','analysis1_2_out_beta','analysis1_2_out_tvalue','analysis1_2_out_pvalue','analysis1_2_out_df','analysis1_2_stdCoef','analysis1_2_CI_low','analysis1_2_alcohol_etasq','analysis1_2_CI_high','analysis1_2_r2c','analysis1_2_r2m')
colnames(results_analysis1_2_marijuana)<-colnames(results_analysis1_2_alcohol)
colnames(results_analysis1_2_caffeine)<-colnames(results_analysis1_2_alcohol)
colnames(results_analysis1_2_tobacco)<-colnames(results_analysis1_2_alcohol)

results0804_analysis1_2<-rbind(results_analysis1_2_alcohol,results_analysis1_2_caffeine,results_analysis1_2_tobacco,results_analysis1_2_marijuana)
results0804_analysis1_2_new <- results0804_analysis1_2 %>% mutate(value=ifelse(results0804_analysis1_2$analysis1_2_out_pvalue<0.05,analysis1_2_stdCoef,NA))

write.csv(results0804_analysis1_2_new,'/share/inspurStorage/home1/zixing/ABCDR/Project_1/new_results/analysis1_2/results0804_analysis1_2_new.csv')








##################
analysis1_2_outcome<-rep()
analysis1_2_exposure<-rep()
analysis1_2_out_beta<-rep()
analysis1_2_out_tvalue <- rep()
analysis1_2_out_pvalue <- rep()
analysis1_2_out_df <- rep()
analysis1_2_stdCoef<-rep()
analysis1_2_CI_low<-rep()
analysis1_2_CI_high<-rep()

analysis1_2_tobacco_exposure<-rep()
analysis1_2_tobacco_out_beta<-rep()
analysis1_2_tobacco_out_tvalue <- rep()
analysis1_2_tobacco_out_pvalue <- rep()
analysis1_2_tobacco_out_df <- rep()
analysis1_2_tobacco_stdCoef<-rep()
analysis1_2_tobacco_CI_low<-rep()
analysis1_2_tobacco_CI_high<-rep()

analysis1_2_marijuana_exposure<-rep()
analysis1_2_marijuana_out_beta<-rep()
analysis1_2_marijuana_out_tvalue <- rep()
analysis1_2_marijuana_out_pvalue <- rep()
analysis1_2_marijuana_out_df <- rep()
analysis1_2_marijuana_stdCoef<-rep()
analysis1_2_marijuana_CI_low<-rep()
analysis1_2_marijuana_CI_high<-rep()

analysis1_2_r2m<-rep()
analysis1_2_r2c<-rep()
analysis1_2_alcohol_etasq<-rep()
analysis1_2_caffeine_etasq<-rep()
analysis1_2_tobacco_etasq<-rep()
analysis1_2_marijuana_etasq<-rep()
analysis1_2_number=1

for(i in name_pes_f){ 
  project1_dataset_baseline_rm_NA_phillip_win_prs2norm_beforemeancentered$y = project1_dataset_baseline_rm_NA_phillip_win_prs2norm_beforemeancentered[,i]
  out <- formula(paste( 'y', "~",paste(c(exposure_name,name_fixed,'(1 | site_id_l/rel_family_id)'),collapse="+")))
  lmer.fit <- lmer(out,data = project1_dataset_baseline_rm_NA_phillip_win_prs2norm_beforemeancentered)
  analysis1_2=summary(lmer.fit)$coefficients
  ci<-confint(lmer.fit,method='Wald')
  
  analysis1_2_alcohol_exposure[analysis1_2_number]<-'alcohol'
  analysis1_2_alcohol_out_beta[analysis1_2_number]<-analysis1_2['alcoholTRUE',1]
  analysis1_2_alcohol_out_tvalue[analysis1_2_number] <- analysis1_2['alcoholTRUE',4]
  analysis1_2_alcohol_out_pvalue[analysis1_2_number] <- analysis1_2['alcoholTRUE',5]
  analysis1_2_alcohol_out_df[analysis1_2_number] <- analysis1_2['alcoholTRUE',3]
  analysis1_2_alcohol_stdCoef[analysis1_2_number]<-stdCoef.merMod(lmer.fit)['alcoholTRUE',1]
  analysis1_2_alcohol_CI_low[analysis1_2_number]<-as.numeric(ci['alcoholTRUE',1])
  analysis1_2_alcohol_CI_high[analysis1_2_number]<-as.numeric(ci['alcoholTRUE',2])
  
  analysis1_2_caffeine_exposure[analysis1_2_number]<-'caffeine'
  analysis1_2_caffeine_out_beta[analysis1_2_number]<-analysis1_2['caffeineTRUE',1]
  analysis1_2_caffeine_out_tvalue[analysis1_2_number] <- analysis1_2['caffeineTRUE',4]
  analysis1_2_caffeine_out_pvalue[analysis1_2_number] <- analysis1_2['caffeineTRUE',5]
  analysis1_2_caffeine_out_df[analysis1_2_number] <- analysis1_2['caffeineTRUE',3]
  analysis1_2_caffeine_stdCoef[analysis1_2_number]<-stdCoef.merMod(lmer.fit)['caffeineTRUE',1]
  analysis1_2_caffeine_CI_low[analysis1_2_number]<-as.numeric(ci['caffeineTRUE',1])
  analysis1_2_caffeine_CI_high[analysis1_2_number]<-as.numeric(ci['caffeineTRUE',2])
  
  analysis1_2_tobacco_exposure[analysis1_2_number]<-'tobacco'
  analysis1_2_tobacco_out_beta[analysis1_2_number]<-analysis1_2['tobaccoTRUE',1]
  analysis1_2_tobacco_out_tvalue[analysis1_2_number] <- analysis1_2['tobaccoTRUE',4]
  analysis1_2_tobacco_out_pvalue[analysis1_2_number] <- analysis1_2['tobaccoTRUE',5]
  analysis1_2_tobacco_out_df[analysis1_2_number] <- analysis1_2['tobaccoTRUE',3]
  analysis1_2_tobacco_stdCoef[analysis1_2_number]<-stdCoef.merMod(lmer.fit)['tobaccoTRUE',1]
  analysis1_2_tobacco_CI_low[analysis1_2_number]<-as.numeric(ci['tobaccoTRUE',1])
  analysis1_2_tobacco_CI_high[analysis1_2_number]<-as.numeric(ci['tobaccoTRUE',2])
  
  analysis1_2_marijuana_exposure[analysis1_2_number]<-'marijuana'
  analysis1_2_marijuana_out_beta[analysis1_2_number]<-analysis1_2['marijuanaTRUE',1]
  analysis1_2_marijuana_out_tvalue[analysis1_2_number] <- analysis1_2['marijuanaTRUE',4]
  analysis1_2_marijuana_out_pvalue[analysis1_2_number] <- analysis1_2['marijuanaTRUE',5]
  analysis1_2_marijuana_out_df[analysis1_2_number] <- analysis1_2['marijuanaTRUE',3]
  analysis1_2_marijuana_stdCoef[analysis1_2_number]<-stdCoef.merMod(lmer.fit)['marijuanaTRUE',1]
  analysis1_2_marijuana_CI_low[analysis1_2_number]<-as.numeric(ci['marijuanaTRUE',1])
  analysis1_2_marijuana_CI_high[analysis1_2_number]<-as.numeric(ci['marijuanaTRUE',2])
  
  analysis1_2_r2c[analysis1_2_number] = as.numeric(performance::r2(lmer.fit)[1])
  analysis1_2_r2m[analysis1_2_number] = as.numeric(performance::r2(lmer.fit)[2])
  analysis1_2_alcohol_etasq[analysis1_2_number] = as.numeric(effectsize::eta_squared(lmer.fit)['1',2])
  analysis1_2_tobacco_etasq[analysis1_2_number] = as.numeric(effectsize::eta_squared(lmer.fit)['2',2])
  analysis1_2_marijuana_etasq[analysis1_2_number] = as.numeric(effectsize::eta_squared(lmer.fit)['3',2])
  analysis1_2_caffeine_etasq[analysis1_2_number] = as.numeric(effectsize::eta_squared(lmer.fit)['4',2])
  
  analysis1_2_outcome[analysis1_2_number] = i
  analysis1_2_number = analysis1_2_number + 1
  
}
results_analysis1_2_alcohol<-data.frame(analysis1_2_outcome,analysis1_2_alcohol_exposure,analysis1_2_alcohol_out_beta,analysis1_2_alcohol_out_tvalue,analysis1_2_alcohol_out_pvalue,analysis1_2_alcohol_out_df,analysis1_2_alcohol_stdCoef,analysis1_2_alcohol_CI_low,analysis1_2_alcohol_etasq,analysis1_2_alcohol_CI_high,analysis1_2_r2c,analysis1_2_r2m)
results_analysis1_2_tobacco<-data.frame(analysis1_2_outcome,analysis1_2_tobacco_exposure,analysis1_2_tobacco_out_beta,analysis1_2_tobacco_out_tvalue,analysis1_2_tobacco_out_pvalue,analysis1_2_tobacco_out_df,analysis1_2_tobacco_stdCoef,analysis1_2_tobacco_CI_low,analysis1_2_tobacco_etasq,analysis1_2_tobacco_CI_high,analysis1_2_r2c,analysis1_2_r2m)
results_analysis1_2_marijuana<-data.frame(analysis1_2_outcome,analysis1_2_marijuana_exposure,analysis1_2_marijuana_out_beta,analysis1_2_marijuana_out_tvalue,analysis1_2_marijuana_out_pvalue,analysis1_2_marijuana_out_df,analysis1_2_marijuana_stdCoef,analysis1_2_marijuana_CI_low,analysis1_2_marijuana_etasq,analysis1_2_marijuana_CI_high,analysis1_2_r2c,analysis1_2_r2m)
results_analysis1_2_caffeine<-data.frame(analysis1_2_outcome,analysis1_2_caffeine_exposure,analysis1_2_caffeine_out_beta,analysis1_2_caffeine_out_tvalue,analysis1_2_caffeine_out_pvalue,analysis1_2_caffeine_out_df,analysis1_2_caffeine_stdCoef,analysis1_2_caffeine_CI_low,analysis1_2_caffeine_etasq,analysis1_2_caffeine_CI_high,analysis1_2_r2c,analysis1_2_r2m)


colnames(results_analysis1_2_alcohol)<-c('analysis1_2_outcome','analysis1_2_exposure','analysis1_2_out_beta','analysis1_2_out_tvalue','analysis1_2_out_pvalue','analysis1_2_out_df','analysis1_2_stdCoef','analysis1_2_CI_low','analysis1_2_alcohol_etasq','analysis1_2_CI_high','analysis1_2_r2c','analysis1_2_r2m')
colnames(results_analysis1_2_marijuana)<-colnames(results_analysis1_2_alcohol)
colnames(results_analysis1_2_caffeine)<-colnames(results_analysis1_2_alcohol)
colnames(results_analysis1_2_tobacco)<-colnames(results_analysis1_2_alcohol)

results0804_analysis1_2<-rbind(results_analysis1_2_alcohol,results_analysis1_2_caffeine,results_analysis1_2_tobacco,results_analysis1_2_marijuana)
results0804_analysis1_2_new <- results0804_analysis1_2 %>% mutate(value=ifelse(results0804_analysis1_2$analysis1_2_out_pvalue<0.05,analysis1_2_stdCoef,NA))

write.csv(results0804_analysis1_2_new,'/share/inspurStorage/home1/zixing/ABCDR/Project_1/new_results/analysis1_2/results0804_analysis1_2_new.csv')


analysis1_2_outcome<-rep()
analysis1_2_exposure<-rep()
analysis1_2_out_beta<-rep()
analysis1_2_out_tvalue <- rep()
analysis1_2_out_pvalue <- rep()
analysis1_2_out_df <- rep()
analysis1_2_stdCoef<-rep()
analysis1_2_CI_low<-rep()
analysis1_2_CI_high<-rep()

analysis1_2_r2m<-rep()
analysis1_2_r2c<-rep()
analysis1_2_etasq<-rep()
analysis1_2_number=1

for(i in EGsum[c(2,8,4,7)]){ 
  project1_dataset_baseline_rm_NA_phillip_win_prs2norm_beforemeancentered$y = project1_dataset_baseline_rm_NA_phillip_win_prs2norm_beforemeancentered[,i]
  for (j in exposure_name){
    project1_dataset_baseline_rm_NA_phillip_win_prs2norm_beforemeancentered$x = project1_dataset_baseline_rm_NA_phillip_win_prs2norm_beforemeancentered[,j]
    out <- formula(paste( 'y', "~",paste(c('x',name_fixed,'(1 | site_id_l/rel_family_id)'),collapse="+")))
    lmer.fit <- lmer(out,data = project1_dataset_baseline_rm_NA_phillip_win_prs2norm_beforemeancentered)
    analysis1_2=summary(lmer.fit)$coefficients
    ci<-confint(lmer.fit,method='Wald')
    
    analysis1_2_exposure[analysis1_2_number]<-j
    analysis1_2_out_beta[analysis1_2_number]<-analysis1_2['xTRUE',1]
    analysis1_2_out_tvalue[analysis1_2_number] <- analysis1_2['xTRUE',4]
    analysis1_2_out_pvalue[analysis1_2_number] <- analysis1_2['xTRUE',5]
    analysis1_2_out_df [analysis1_2_number]<- analysis1_2['xTRUE',3]
    analysis1_2_stdCoef[analysis1_2_number]<-stdCoef.merMod(lmer.fit)['xTRUE',1]
    analysis1_2_CI_low[analysis1_2_number]<-as.numeric(ci['xTRUE',1])
    analysis1_2_CI_high[analysis1_2_number]<-as.numeric(ci['xTRUE',2])
    
    analysis1_2_r2c[analysis1_2_number] = as.numeric(performance::r2(lmer.fit)[1])
    analysis1_2_r2m[analysis1_2_number] = as.numeric(performance::r2(lmer.fit)[2])
    analysis1_2_etasq[analysis1_2_number] = as.numeric(effectsize::eta_squared(lmer.fit)['1',2])
    
    analysis1_2_outcome[analysis1_2_number] = i
    analysis1_2_number = analysis1_2_number + 1
  }
  
}

results0602_analysis1_2 = data.frame(analysis1_2_outcome, analysis1_2_exposure, analysis1_2_out_tvalue, analysis1_2_out_beta, analysis1_2_stdCoef, analysis1_2_r2c,analysis1_2_r2m,analysis1_2_etasq,analysis1_2_CI_low,analysis1_2_CI_high,analysis1_2_out_pvalue,analysis1_2_out_df)

results0602_analysis1_2new <- results0602_analysis1_2 %>% mutate(value=ifelse(results0602_analysis1_2$analysis1_2_out_pvalue<0.05,analysis1_2_stdCoef,NA))

write.csv(results0602_analysis1_2new,'/share/inspurStorage/home1/zixing/ABCDR/Project_1/new_results/analysis1_2/results0602_analysis1_2new.csv')
# load packages
library(ggplot2) # ggplot() for plotting
library(dplyr) # data reformatting
library(tidyr) # data reformatting
library(stringr) # string manipulation


# #modified ggplot
# p <- ggplot(results_analysis1_2new, aes(x=analysis1_2_exposure, y=analysis1_2_outcome, fill=value))+
#   #add border white colour of line thickness 0.2
#   geom_tile(colour="white", size=0.2)+scale_colour_gradient(palette = "Set1",na.value = "grey90")+
#   #remove x and y axis labels
#   labs(x="", y="")+
#   #remove extra space
#   scale_y_discrete(expand=c(0, 0))+
#   #set a base size for all fonts
#   theme_grey(base_size=12)+
#   guides(fill=guide_legend(title="Effect size"))+
#   labs(x="", y="", title="Association between prenatal substance exposure and PES and genetic factors")+
#   #theme options
#   theme(
#     #bold font for legend text
#     legend.text=element_text(face="bold"),
#     #set thickness of axis ticks
#     axis.ticks=element_line(size=0.4),
#     #remove plot background
#     plot.background=element_blank(),
#     #remove plot border
#     panel.border=element_blank()
#   )
# 
# ggsave(p, filename="x-z.png",height=5.5, width=9.8, units="in", dpi=200)
# 
# 






# 
# # Define line colors
# colors_line1_2 <-  c(scales::alpha("darkorange", 0.9),
#                   scales::alpha("darkgreen", 0.9),
#                   scales::alpha("darkblue", 0.9),
#                   scales::alpha("purple", 0.9))
# legend_colors_line <-  c(scales::alpha("darkorange", 0.9),
#                          scales::alpha("darkgreen", 0.9),
#                          scales::alpha("darkblue", 0.9),
#                          scales::alpha("purple", 0.9))
# 
# create_beautiful_radarchart1_2 <- function(data, color = "#00AFBB", 
#                                             vlabels = colnames(data),vlcex = NULL,
#                                             caxislabels = NULL, title = NULL, ...){
#   radarchart(
#     data, axistype = 1,
#     # Customize the polygon
#     pcol = colors_line1_2, pfcol = colors_fill, plwd = 2,
#     # Customize the grid
#     cglcol = "grey", cglty = 1, cglwd = 0.8,
#     # Customize the axis
#     axislabcol = "grey", 
#     # Variable labels
#     vlcex = vlcex, vlabels = vlabels,
#     caxislabels = caxislabels, title = title, ...
#   )
# }
# 
# #x<-exposure_name
# #z<-name_pes_f,name_pes_p,name_pes_s
# name_pes<-c(name_pes_f,name_pes_p,name_pes_s)
# analysis1_2_outcome<-rep()
# analysis1_2_exposure<-rep()
# analysis1_2_out_beta<-rep()
# analysis1_2_stdCoef<-rep()
# analysis1_2_r.squared<-rep()
# analysis1_2_out_tvalue <- rep()
# analysis1_2_number=1
# 
project1_dataset_baseline_rm_NA_phillip_win_prs3<-project1_dataset_baseline_rm_NA_phillip_win_prs 

for (i in exposure_name){
  project1_dataset_baseline_rm_NA_phillip_win_prs3[,i]<-ifelse(project1_dataset_baseline_rm_NA_phillip_win_prs[,i]==FALSE,0,1)
}

# 
analysis1_2_outcome<-rep()
analysis1_2_exposure<-rep()
analysis1_2_out_beta<-rep()
analysis1_2_stdCoef<-rep()
analysis1_2_number=1

for(i in exposure_name){ 
  for(j in name_genetic[31:36]){
    project1_dataset_baseline_rm_NA_phillip_win_prs3$y = project1_dataset_baseline_rm_NA_phillip_win_prs3[,i]
    project1_dataset_baseline_rm_NA_phillip_win_prs3$x = project1_dataset_baseline_rm_NA_phillip_win_prs3[,j]
    out <- formula(paste( 'y', "~",paste(c('x',name_fixed,name_PC,'(1 | site_id_l/rel_family_id)'),collapse="+")))
    lmer.fit <- lmer(out,data = project1_dataset_baseline_rm_NA_phillip_win_prs3)
    analysis1_2=summary(lmer.fit)$coefficients
    analysis1_2_out_beta[analysis1_2_number] = analysis1_2['x',1]
    analysis1_2_stdCoef[analysis1_2_number] = stdCoef.merMod(lmer.fit)['x',1]
    analysis1_2_outcome[analysis1_2_number] = i
    analysis1_2_exposure[analysis1_2_number] = j
    analysis1_2_number = analysis1_2_number + 1
    
  }
}
results_exposure_prs_stdcoef<-data.frame(analysis1_2_outcome,analysis1_2_exposure,analysis1_2_out_beta,analysis1_2_stdCoef)

results_exposure_prs_stdcoef_alcohol<-results_exposure_prs_stdcoef[results_exposure_prs_stdcoef$analysis1_2_outcome=='alcohol',]
results_exposure_prs_stdcoef_tobacco<-results_exposure_prs_stdcoef[results_exposure_prs_stdcoef$analysis1_2_outcome=='tobacco',]
results_exposure_prs_stdcoef_marijuana<-results_exposure_prs_stdcoef[results_exposure_prs_stdcoef$analysis1_2_outcome=='marijuana',]
results_exposure_prs_stdcoef_caffeine<-results_exposure_prs_stdcoef[results_exposure_prs_stdcoef$analysis1_2_outcome=='caffeine',]


name_alc_prs<-paste0("alc_",name_genetic[31:36])
project1_dataset_baseline_rm_NA_phillip_win_prs2[,name_alc_prs]<-project1_dataset_baseline_rm_NA_phillip_win_prs2[,name_genetic[31:36]]

for ( i in 1:length(name_alc_prs)){
  project1_dataset_baseline_rm_NA_phillip_win_prs2[,name_alc_prs[i]]<-min_max_norm(results_exposure_prs_stdcoef_alcohol$analysis1_2_stdCoef[i]*project1_dataset_baseline_rm_NA_phillip_win_prs2[,name_alc_prs[i]])
}
project1_dataset_baseline_rm_NA_phillip_win_prs2[,name_alc_prs[6]]<- -project1_dataset_baseline_rm_NA_phillip_win_prs2[,name_alc_prs[6]]

project1_dataset_baseline_rm_NA_phillip_win_prs2$alc_genetic<-rowSums(project1_dataset_baseline_rm_NA_phillip_win_prs2[c(name_newfam,name_alc_prs)],na.rm=TRUE) 


name_tob_prs<-paste0("tob_",name_genetic[31:36])
project1_dataset_baseline_rm_NA_phillip_win_prs2[,name_tob_prs]<-project1_dataset_baseline_rm_NA_phillip_win_prs2[,name_genetic[31:36]]

for ( i in 1:length(name_tob_prs)){
  project1_dataset_baseline_rm_NA_phillip_win_prs2[,name_tob_prs[i]]<-min_max_norm(results_exposure_prs_stdcoef_tobacco$analysis1_2_stdCoef[i]*project1_dataset_baseline_rm_NA_phillip_win_prs2[,name_tob_prs[i]])
}
project1_dataset_baseline_rm_NA_phillip_win_prs2[,name_tob_prs[6]]<- -project1_dataset_baseline_rm_NA_phillip_win_prs2[,name_tob_prs[6]]

project1_dataset_baseline_rm_NA_phillip_win_prs2$tob_genetic<-rowSums(project1_dataset_baseline_rm_NA_phillip_win_prs2[c(name_newfam,name_tob_prs)],na.rm=TRUE) 


name_caf_prs<-paste0("caf_",name_genetic[31:36])
project1_dataset_baseline_rm_NA_phillip_win_prs2[,name_caf_prs]<-project1_dataset_baseline_rm_NA_phillip_win_prs2[,name_genetic[31:36]]

for ( i in 1:length(name_caf_prs)){
  project1_dataset_baseline_rm_NA_phillip_win_prs2[,name_caf_prs[i]]<-min_max_norm(results_exposure_prs_stdcoef_caffeine$analysis1_2_stdCoef[i]*project1_dataset_baseline_rm_NA_phillip_win_prs2[,name_caf_prs[i]])
}
project1_dataset_baseline_rm_NA_phillip_win_prs2[,name_caf_prs[6]]<- -project1_dataset_baseline_rm_NA_phillip_win_prs2[,name_caf_prs[6]]

project1_dataset_baseline_rm_NA_phillip_win_prs2$caf_genetic<-rowSums(project1_dataset_baseline_rm_NA_phillip_win_prs2[c(name_newfam,name_caf_prs)],na.rm=TRUE) 




name_maj_prs<-paste0("maj_",name_genetic[31:36])
project1_dataset_baseline_rm_NA_phillip_win_prs2[,name_maj_prs]<-project1_dataset_baseline_rm_NA_phillip_win_prs2[,name_genetic[31:36]]

for ( i in 1:length(name_maj_prs)){
  project1_dataset_baseline_rm_NA_phillip_win_prs2[,name_maj_prs[i]]<-min_max_norm(results_exposure_prs_stdcoef_marijuana$analysis1_2_stdCoef[i]*project1_dataset_baseline_rm_NA_phillip_win_prs2[,name_maj_prs[i]])
}
project1_dataset_baseline_rm_NA_phillip_win_prs2[,name_maj_prs[6]]<- -project1_dataset_baseline_rm_NA_phillip_win_prs2[,name_maj_prs[6]]

project1_dataset_baseline_rm_NA_phillip_win_prs2$maj_genetic<-rowSums(project1_dataset_baseline_rm_NA_phillip_win_prs2[c(name_newfam,name_maj_prs)],na.rm=TRUE) 


project1_dataset_baseline_rm_NA_phillip_win_prs4<-project1_dataset_baseline_rm_NA_phillip_win_prs2 

for (i in exposure_name){
  project1_dataset_baseline_rm_NA_phillip_win_prs4[,i]<-ifelse(project1_dataset_baseline_rm_NA_phillip_win_prs2[,i]==FALSE,0,1)
}
