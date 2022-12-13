rm(list=ls())
gc()

####################
## Load libraries ##
install.packages('pheatmap')
install.packages("lmerTest")
install.packages("robustHD")
library(pheatmap)
library(lmerTest)
library(gamm4)
library(nlme)
library(psych)
library(mice)
library(ggplot2)
library(data.table)
library(dplyr)
library(readxl)
library(purrr)
library(broom)
library(ggpubr)
library(mgcv)
library(RColorBrewer)
library(robustHD)
library(foreign)
library(lme4)

setwd('/share/inspurStorage/home1/zixing/ABCDR/Project_1')
# dataset_baseline<-readRDS("dataset_baseline.Rds")
dataset<-readRDS("dataset4.0.Rds")
# datasetraw<-readRDS("dataset4.0raw.Rds")

dataset$devhx_caffeine_11.x[dataset$devhx_caffeine_11.x==-1] = NA
dataset$devhx_caffeine_11.x[dataset$devhx_caffeine_11.x>0] = TRUE
dataset$devhx_caffeine_11.x[dataset$devhx_caffeine_11.x==0] = FALSE
#integration
#substnace use
dataset$tobacco<-as.factor(dataset$devhx_8_tobacco.x|dataset$devhx_9_tobacco.x)
dataset$alcohol<-as.factor(dataset$devhx_8_alcohol.x|dataset$devhx_9_alcohol.x)
dataset$marijuana<-as.factor(dataset$devhx_8_marijuana.x|dataset$devhx_9_marijuana.x)
dataset$caffeine<-as.factor(dataset$devhx_caffeine_11.x)

dataset$coc_crack<-as.factor(dataset$devhx_8_coc_crack.x|dataset$devhx_9_coc_crack.x)
dataset$her_morph<-as.factor(dataset$devhx_8_her_morph.x|dataset$devhx_9_her_morph.x)
dataset$oxycont<-as.factor(dataset$devhx_8_oxycont.x|dataset$devhx_9_oxycont.x)

dataset$rel_family_id<-as.factor(dataset$rel_family_id)
#selecting nesting variables
ind_nest = c(which(names(dataset)=="rel_family_id"), which(names(dataset)=="site_id_l"));summary(dataset[,ind_nest])

#########################################
## Select & process data for analyses  ##
#########################################
setwd('/share/inspurStorage/home1/zixing/ABCDR/Project_1/xlsx')
total_number=dim(dataset)[1]

## Select exposures
exposures <- read_excel("/share/inspurStorage/home1/zixing/ABCDR/Project_1/xlsx/exposures_x2.xlsx",1, col_names = TRUE)
exposures=t(data.frame(exposures));ind_expo <- c()
for(i in 1:length(exposures['Name',])){
  ind_expo = c(ind_expo,c(which(names(dataset)==exposures['Name',][i])))
}
exposure_name<-names(dataset)[ind_expo];summary(dataset[exposure_name])


#separated by sex
ind_male=which(dataset$sex=='M')
ind_female=which(dataset$sex=='F')


#######ZZZZZZZZZZZ#######################
#########################################
#########Physical health#################
#########################################

#interview age
dataset$interview_age<-(dataset$interview_age)/12

#body mass index

# Average Measured Weight (lbs):If three measurements were obtained, 
# the two closest measurements will be averaged. 
# Should the third measurement fall equally between the first two measurements, 
# all three will be averaged.
anthroweight_3m<-c('anthroweight1lb','anthroweight2lb','anthroweight3lb')
anthroweight_2m<-c('anthroweight1lb','anthroweight2lb')

dataset$anthroweight12lb<-abs(dataset$anthroweight1lb-dataset$anthroweight2lb)
dataset$anthroweight13lb<-abs(dataset$anthroweight1lb-dataset$anthroweight3lb)
dataset$anthroweight23lb<-abs(dataset$anthroweight2lb-dataset$anthroweight3lb)

dataset$anthroweight_avg<-ifelse(is.na(dataset$anthroweight3lb)==TRUE,rowMeans(dataset[c('anthroweight1lb','anthroweight2lb')],na.rm=TRUE),
                                 ifelse(is.na(dataset$anthroweight3lb)==TRUE
                                      & dataset$anthroweight1lb<dataset$anthroweight13lb & dataset$anthroweight12lb<dataset$anthroweight23lb, rowMeans(dataset[c('anthroweight1lb','anthroweight2lb')],na.rm=TRUE)
                                        ,ifelse(is.na(dataset$anthroweight3lb)==TRUE&
                                                  dataset$anthroweight13lb<dataset$anthroweight12lb & dataset$anthroweight13lb<dataset$anthroweight23lb, rowMeans(dataset[c('anthroweight1lb','anthroweight3lb')],na.rm=TRUE)
                                                ,ifelse(is.na(dataset$anthroweight3lb)==TRUE&
                                                          dataset$anthroweight23lb<dataset$anthroweight12lb & dataset$anthroweight23lb<dataset$anthroweight13lb, rowMeans(dataset[c('anthroweight2lb','anthroweight3lb')],na.rm=TRUE)
                                                        ,rowMeans(dataset[c('anthroweight1lb','anthroweight2lb','anthroweight3lb')],na.rm=TRUE)))))


dataset$bmi<-703*(dataset$anthroweight_avg)/((dataset$anthroheightcalc)^2)

#screentime
sqt_weekday<-c('screen1_wkdy_y','screen2_wkdy_y','screen3_wkdy_y','screen4_wkdy_y','screen5_wkdy_y','screen_wkdy_y')
sqt_weekend<-c('screen7_wknd_y','screen8_wknd_y','screen9_wknd_y','screen10_wknd_y','screen11_wknd_y','screen12_wknd_y')
dataset$screen_weekday<-rowSums(dataset[sqt_weekday],na.rm=TRUE)
dataset$screen_weekend<-rowSums(dataset[sqt_weekend],na.rm=TRUE)

#ABCD Youth Pubertal Development Scale and Menstrual Cycle Survey History (PDMS)
for (indx in 1:total_number) {
  #male
  if (indx %in% ind_male){
    dataset$ypdms_male_sum[indx]=sum(as.numeric(unlist(dataset$pds_ht2_y))[indx],as.numeric(unlist(dataset$pds_bdyhair_y))[indx],as.numeric(unlist(dataset$pds_skin2_y))[indx],as.numeric(unlist(dataset$pds_m4_y))[indx],as.numeric(unlist(dataset$pds_m5_y))[indx],na.rm = TRUE)
    dataset$ypdms_male_mean[indx]=mean(dataset$pds_male_sum[indx],na.rm = TRUE)
    dataset$ypdms_sum_y_m[indx]=sum(as.numeric(unlist(dataset$pds_bdyhair_y[indx])), as.numeric(unlist(dataset$pds_m4_y[indx])),as.numeric(unlist(dataset$pds_m5_y[indx])),na.rm = TRUE)
  }
  #female
  else{
    dataset$ypdms_female_sum[indx]=sum(as.numeric(unlist(dataset$pds_ht2_y))[indx],as.numeric(unlist(dataset$pds_bdyhair_y))[indx],as.numeric(unlist(dataset$pds_skin2_y))[indx],as.numeric(unlist(dataset$pds_f4_2_y))[indx],as.numeric(unlist(dataset$pds_f5_y))[indx],na.rm = TRUE)
    dataset$ypdms_female_mean[indx]=mean(dataset$pds_female_sum[indx],na.rm = TRUE)
    dataset$ypdms_sum_y_m[indx]=sum(as.numeric(unlist(dataset$pds_bdyhair_y[indx])), as.numeric(unlist(dataset$pds_f4_2_y[indx])),as.numeric(unlist(dataset$pds_f5_y[indx])),na.rm = TRUE)
  }
}

###################
#1:prepuberty        3
#2:early puberty  >= 4 and <= 5
#3:mid puberty    >= 6 and <= 8
#4:late puberty   >= 9 and <=11
#5:post puberty      12
puberty.category <- function(s){
  ifelse(s==3,return(1),ifelse(s>= 4 & s<= 5,return(2),ifelse(s>= 6 & s<= 8,return(3),ifelse(s>= 9 & s<= 11,return(4),ifelse(s==12,return(5),return(NA))))))
}

ypdms_cat<-c()
for (i in dataset$ypdms_sum_y_m){
  ypdms_cat<-c(ypdms_cat,puberty.category(as.numeric(unlist(i))))
}
dataset$ypdms_y_ss_cat_sum<-as.factor(ypdms_cat)

#take vitamins
#-1 = Not applicable No aplica
dataset$devhx_10.x[dataset$devhx_10.x==-1] = NA
dataset$devhx_vitamin<-as.factor(dataset$devhx_10.x)

#planned pregnancy
#-1 = Not applicable No aplica
dataset$planned_pregnancy<-as.factor(dataset$devhx_6_p.x)

dataset$maternal_age<-dataset$devhx_3_p.x#maternal age
dataset$paternal_age<-dataset$devhx_4_p.x#paternal age

dataset$born_premature<-as.factor(dataset$devhx_12a_p.x)#born_premature 

dataset$incubator_days<-dataset$devhx_15.x
dataset$firstyear_fever_104degrees_days<-dataset$devhx_16_p.x
dataset$firstyear_infections_ser_ill_days<-dataset$devhx_17_p.x
dataset$breast_fed_months<-dataset$devhx_18_p.x

#obstetric complications (at least 1 condition = present)
obs_comps<-c('devhx_10a3_p.x','devhx_10b3_p.x','devhx_10c3_p.x','devhx_10d3_p.x','devhx_10e3_p.x','devhx_10f3_p.x','devhx_10g3_p.x','devhx_10h3_p.x','devhx_10i3_p.x','devhx_10j3_p.x','devhx_10k3_p.x','devhx_10l3_p.x','devhx_10m3_p.x')
dataset$obsteric_complication<-as.factor(ifelse(rowSums(dataset[obs_comps],na.rm=TRUE)>0,1,0))

names(dataset)[which(names(dataset)=='devhx_10a3_p.x'):which(names(dataset)=='devhx_10m3_p.x')]<-c('Severe.nausea.vomiting','Heavy.bleeding','Pre_eclampsia.eclampsia.toxemia','Severe.gall.bladder.attack','Persistent.proteinuria','Rubella','Severe.anemia','UTI','Preg_related.diabetes','Preg_related.HBP','problems.with.the.placenta','accident.or.injury','other.conditions')

#birth complications (at least 1 condition = present)
birth_comps<-c('devhx_14a3_p.x','devhx_14b3_p.x','devhx_14c3_p.x','devhx_14d3_p.x','devhx_14e3_p.x','devhx_14f3_p.x','devhx_14g3_p.x','devhx_14h3_p.x')
dataset$birth_complication<-as.factor(ifelse(rowSums(dataset[birth_comps],na.rm=TRUE)>0,1,0))
names(dataset)[which(names(dataset)=='devhx_14a3_p.x'):which(names(dataset)=='devhx_14h3_p.x')]<-c('Blue.at.birth','slow.heart.beat','not.breathe.at.first','convulsions','Jaundice.needing.treatment','Required.oxygen','Required.blood.transfusion','Rh.incompatibility')


#########################################
###########Mental health#################
#########################################

#family history problem :substance use and mental health problem
famhx_ss_relative <- read_excel("/share/inspurStorage/home1/zixing/ABCDR/Project_1/xlsx/famhx_relative.xlsx", col_names = TRUE);famhx_ss_relative=t(data.frame(famhx_ss_relative));rownames(famhx_ss_relative)

name_relative_alcohol <- famhx_ss_relative[1,];name_relative_dg <- famhx_ss_relative[2,];name_relative_dprs <- famhx_ss_relative[3,];name_relative_ma <- famhx_ss_relative[4,];name_relative_vs<- famhx_ss_relative[5,];name_relative_trb <- famhx_ss_relative[6,];name_relative_nrv <- famhx_ss_relative[7,];name_relative_prf <- famhx_ss_relative[8,];name_relative_hspd <- famhx_ss_relative[9,];name_relative_scd <- famhx_ss_relative[10,]

dataset$famhx_ss_relative_prob_alc_p<-ifelse(rowSums(dataset[name_relative_alcohol],na.rm=TRUE)>0,1,0)
dataset$famhx_ss_relative_prob_dg_p<-ifelse(rowSums(dataset[name_relative_dg],na.rm=TRUE)>0,1,0)
dataset$famhx_ss_relative_prob_dprs_p<-ifelse(rowSums(dataset[name_relative_dprs],na.rm=TRUE)>0,1,0)
dataset$famhx_ss_relative_prob_ma_p<-ifelse(rowSums(dataset[name_relative_ma],na.rm=TRUE)>0,1,0)
dataset$famhx_ss_relative_prob_vs_p<-ifelse(rowSums(dataset[name_relative_vs],na.rm=TRUE)>0,1,0)
dataset$famhx_ss_relative_prob_trb_p<-ifelse(rowSums(dataset[name_relative_trb],na.rm=TRUE)>0,1,0)
dataset$famhx_ss_relative_prob_nrv_p<-ifelse(rowSums(dataset[name_relative_nrv],na.rm=TRUE)>0,1,0)
dataset$famhx_ss_relative_prob_prf_p<-ifelse(rowSums(dataset[name_relative_prf],na.rm=TRUE)>0,1,0)
dataset$famhx_ss_relative_prob_hspd_p<-ifelse(rowSums(dataset[name_relative_hspd],na.rm=TRUE)>0,1,0)
dataset$famhx_ss_relative_prob_scd_p<-ifelse(rowSums(dataset[name_relative_scd],na.rm=TRUE)>0,1,0)


####################Substance use###################
#ABCD Youth Substance Use Interview
tlfb_alc_v<-c('tlfb_alc_use','tlfb_alc_last_calc')
tlfb_tob_v<-c('tlfb_cig_use','tlfb_ecig_use','tlfb_hookah_use','tlfb_chew_use')
tlfb_mj_v<-c('tlfb_mj_use','tlfb_blunt_use','tlfb_mj_use', 'tlfb_edible_use','tlfb_mj_conc_use')

dataset$tlfb_alc<-ifelse(rowSums(dataset[tlfb_alc_v],na.rm=TRUE)>0,1,0)
dataset$tlfb_tob<-ifelse(rowSums(dataset[tlfb_tob_v],na.rm=TRUE)>0,1,0)
dataset$tlfb_mj<-ifelse(rowSums(dataset[tlfb_mj_v],na.rm=TRUE)>0,1,0)
dataset$tlfb_caf<-dataset$caff_max_type
dataset$tlfb_caf[dataset$tlfb_caf>0]=1

#ABCD Parent Community Risk and Protective Factors (CRPF)
su_risk_p<-c('su_risk_p_1','su_risk_p_2','su_risk_p_3','su_risk_p_4','su_risk_p_5','su_risk_p_7','su_risk_p_8','su_risk_p_10','su_risk_p_11','su_risk_p_12','su_risk_p_13')

#how hard to get those substances
#the higher the score, more risks it has
#0 = Very hard ; 1 = Sort of hard ; 2 = Sort of easy ; 3 = Very easy; 4 = Don't know
dataset[su_risk_p][dataset[su_risk_p] ==4 ] <- NA

#Is "medical marijuana" (marijuana prescribed by a doctor) legal in your state? 
#0 = Yes; 1 = No; 2 = I do not know 
#change to 0 = No; 1 = Yes; NA = I do not know 
dataset$su_risk_p_6[dataset$su_risk_p_6 ==2 ] <- NA;dataset$su_risk_p_6[dataset$su_risk_p_6 ==0 ] <- 2;dataset$su_risk_p_6[dataset$su_risk_p_6 ==1 ] <- 0;dataset$su_risk_p_6[dataset$su_risk_p_6 ==2 ] <- 1

#su_risk_p_9
#How many of your friends or family members have a medical marijuana prescription?
#0 = 0; 8 = 1-3; 9 = 4-6; 10 = 7-9; 11 = 10 or more 
su_risks<-c(su_risk_p,'su_risk_p_6','su_risk_p_9')
dataset$su_risks_Perceived_availability<-rowSums(dataset[su_risks],na.rm=TRUE)

#ABCD Parent Demographics Survey
dataset['demo_prnt_prtnr_v2'][dataset['demo_prnt_prtnr_v2'] ==2 ] <- 0

dataset$demo_prnt_prtnr_ed_v2<-pmax(dataset$demo_prnt_ed_v2,dataset$demo_prtnr_ed_v2,na.rm=TRUE)

######PHYSICAL HEALTH##############
dataset$Asthma<-as.factor(dataset$medhx_2a)
dataset$Allergies<-as.factor(dataset$medhx_2b)
dataset$Brain_Injury<-as.factor(dataset$medhx_2c)
dataset$Bronchitis<-as.factor(dataset$medhx_2d)
dataset$Cancer_or_Leukemia<-as.factor(dataset$medhx_2e)
dataset$Cerebral_Palsy<-as.factor(dataset$medhx_2f)
dataset$Diabetes<-as.factor(dataset$medhx_2g)

dataset$Epilepsy_or_Seizures<-as.factor(dataset$medhx_2h)
dataset$Hearing_Problem<-as.factor(dataset$medhx_2i)
dataset$Kidney_Disease<-as.factor(dataset$medhx_2j)
dataset$Lead_Poisoning<-as.factor(dataset$medhx_2k)
dataset$Muscular_Dystrophy_Since_LAST_meeting<-as.factor(dataset$medhx_2l)
dataset$Multiple_Sclerosis<-as.factor(dataset$medhx_2m)
dataset$Problems_with_Vision<-as.factor(dataset$medhx_2n)
dataset$Problems_with_Heart<-as.factor(dataset$medhx_2o)
dataset$Sickle_Cell_Anemia<-as.factor(dataset$medhx_2p)
dataset$Very_Bad_Headaches<-as.factor(dataset$medhx_2q)
dataset$Operation<-as.factor(dataset$medhx_2r)
dataset$Other_Illness<-as.factor(dataset$medhx_2s)
#######################################################################
dataset$BrokenBones_times<-as.factor(dataset$medhx_ss_6a_times_p)
dataset$Sprains_times<-as.factor(dataset$medhx_ss_6b_times_p)
dataset$CutsScrapes_times<-as.factor(dataset$medhx_ss_6c_times_p)
dataset$Stitches_times<-as.factor(dataset$medhx_ss_6d_times_p)
dataset$OTRSeriousWounds_times<-as.factor(dataset$medhx_ss_6e_times_p)
dataset$Falls_times<-as.factor(dataset$medhx_ss_6f_times_p)
dataset$Burns_times<-as.factor(dataset$medhx_ss_6g_times_p)
dataset$HighFever_times<-as.factor(dataset$medhx_ss_6h_times_p)
dataset$HeadInjury_times<-as.factor(dataset$medhx_ss_6i_times_p)
dataset$KnockedUnconscious_times<-as.factor(dataset$medhx_ss_6j_times_p)
dataset$Bruises_times<-as.factor(dataset$medhx_ss_6k_times_p)
dataset$AsthmaAttack_times<-as.factor(dataset$medhx_ss_6l_times_p)
dataset$BrokenTeeth_times<-as.factor(dataset$medhx_ss_6m_times_p)
dataset$AnimalBite_times<-as.factor(dataset$medhx_ss_6n_times_p)
dataset$Overdose_times<-as.factor(dataset$medhx_ss_6o_times_p)
dataset$Seizure_times<-as.factor(dataset$medhx_ss_6p_times_p)
dataset$AccidentalPoisoning_times<-as.factor(dataset$medhx_ss_6q_times_p)
dataset$GunShotWound_times<-as.factor(dataset$medhx_ss_6r_times_p)
dataset$WoundfromOTRWeapon_times<-as.factor(dataset$medhx_ss_6s_times_p)
dataset$OTRHospitalizations_times<-as.factor(dataset$medhx_ss_6t_times_p)

dataset$ehi_y_ss_scoreb<-as.factor(dataset$ehi_y_ss_scoreb)


#########YYYYYYYYYYYYYYYYYY############

#ABCD Parent Diagnostic Interview for DSM-5 Full (KSADS-5)
ksad_cat_fle <- read_excel("/share/inspurStorage/home1/zixing/ABCDR/Project_1/xlsx/ksad_cat.xlsx", col_names = TRUE)
ksad_cat_fle=t(data.frame(ksad_cat_fle))
rownames(ksad_cat_fle)

name_ksad_all <- ksad_cat_fle[1,]
summary(dataset[name_ksad_all])

#remove na values larger than 1 such as 555,888
dataset[name_ksad_all][dataset[name_ksad_all] > 1 ] <- NA

dataset$KSADS.Depression<-ifelse(rowSums(dataset[ksad_cat_fle[3,][complete.cases(ksad_cat_fle[3,])]],na.rm=TRUE)>0,1,0)

dataset$KSADS.Bipolar.and.Related.Disorder<-ifelse(rowSums(dataset[ksad_cat_fle[5,][complete.cases(ksad_cat_fle[5,])]],na.rm=TRUE)>0,1,0)

dataset$KSADS.Hallucinations<-ifelse(rowSums(dataset[ksad_cat_fle[13,][complete.cases(ksad_cat_fle[13,])]],na.rm=TRUE)>0,1,0)
dataset$KSADS.Auditory.Hallucinations<-ifelse(rowSums(dataset[ksad_cat_fle[11,][complete.cases(ksad_cat_fle[11,])]],na.rm=TRUE)>0,1,0)
dataset$KSADS.Delusions<-ifelse(rowSums(dataset[ksad_cat_fle[12,][complete.cases(ksad_cat_fle[12,])]],na.rm=TRUE)>0,1,0)
dataset$KSADS.Associated.Psychotic.Symptoms<-ifelse(rowSums(dataset[ksad_cat_fle[10,][complete.cases(ksad_cat_fle[10,])]],na.rm=TRUE)>0,1,0)

dataset$KSADS.Panic.and.Other.Specified.Anxiety.Disorder<-ifelse(rowSums(dataset[ksad_cat_fle[15,][complete.cases(ksad_cat_fle[15,])]],na.rm=TRUE)>0,1,0)

dataset$KSADS.Agoraphobia.and.Other.Specified.Anxiety.Disorder<-ifelse(rowSums(dataset[ksad_cat_fle[17,][complete.cases(ksad_cat_fle[17,])]],na.rm=TRUE)>0,1,0)

dataset$KSADS.Separation.and.Other.Specified.Anxiety.Disorder<-ifelse(rowSums(dataset[ksad_cat_fle[19,][complete.cases(ksad_cat_fle[19,])]],na.rm=TRUE)>0,1,0)

dataset$KSADS.Social.and.Other.Specified.Anxiety.Disorder<-ifelse(rowSums(dataset[ksad_cat_fle[21,][complete.cases(ksad_cat_fle[21,])]],na.rm=TRUE)>0,1,0)

dataset$KSADS.Specific.Phobia<-ifelse(rowSums(dataset[ksad_cat_fle[23,][complete.cases(ksad_cat_fle[23,])]],na.rm=TRUE)>0,1,0)

dataset$KSADS.Generalized.Anxiety.Disorder.and.other.Specified.Anxiety.Disorder<-ifelse(rowSums(dataset[ksad_cat_fle[25,][complete.cases(ksad_cat_fle[25,])]],na.rm=TRUE)>0,1,0)

dataset$KSADS.Obsessive.Compulsive.and.Related.Disorder<-ifelse(rowSums(dataset[ksad_cat_fle[27,][complete.cases(ksad_cat_fle[27,])]],na.rm=TRUE)>0,1,0)

dataset$KSADS.Encopresis<-ifelse(rowSums(dataset[ksad_cat_fle[29,][complete.cases(ksad_cat_fle[29,])]],na.rm=TRUE)>0,1,0)

dataset$KSADS.Feeding.or.Eating.Disorder<-ifelse(rowSums(dataset[ksad_cat_fle[31,][complete.cases(ksad_cat_fle[31,])]],na.rm=TRUE)>0,1,0)

dataset$KSADS.Attention.Deficit.Hyperactivity.Disorder<-ifelse(rowSums(dataset[ksad_cat_fle[33,][complete.cases(ksad_cat_fle[33,])]],na.rm=TRUE)>0,1,0)

dataset$KSADS.Oppositional.Defiant.Disorder<-ifelse(rowSums(dataset[ksad_cat_fle[35,][complete.cases(ksad_cat_fle[35,])]],na.rm=TRUE)>0,1,0)

dataset$KSADS.Conduct.Disorder<-ifelse(rowSums(dataset[ksad_cat_fle[37,][complete.cases(ksad_cat_fle[37,])]],na.rm=TRUE)>0,1,0)

dataset$KSADS.Alcohol.Use.and.Alcohol.Related.Disorder<-ifelse(rowSums(dataset[ksad_cat_fle[39,][complete.cases(ksad_cat_fle[39,])]],na.rm=TRUE)>0,1,0)

dataset$KSADS.Substance.Related.Disorder<-ifelse(rowSums(dataset[ksad_cat_fle[41,][complete.cases(ksad_cat_fle[41,])]],na.rm=TRUE)>0,1,0)

dataset$KSADS.PTSD.and.Other.Specified.Trauma.and.Stressor.Related.Disorder<-ifelse(rowSums(dataset[ksad_cat_fle[43,][complete.cases(ksad_cat_fle[43,])]],na.rm=TRUE)>0,1,0)

dataset$KSADS.Insomnia<-ifelse(rowSums(dataset[ksad_cat_fle[45,][complete.cases(ksad_cat_fle[45,])]],na.rm=TRUE)>0,1,0)

dataset$KSADS.homicidal.ideation.and.behavior<-ifelse(rowSums(dataset[ksad_cat_fle[47,][complete.cases(ksad_cat_fle[47,])]],na.rm=TRUE)>0,1,0)

dataset$KSADS.Selective.Mutism.and.Other.Specified.Anxiety.Disorder<-ifelse(rowSums(dataset[ksad_cat_fle[49,][complete.cases(ksad_cat_fle[49,])]],na.rm=TRUE)>0,1,0)

#sleep disturbance scale
add_sds_name<-names(dataset)[which(colnames(dataset)=='sleepdisturb1_p'):which(colnames(dataset)=='sleepdisturb26_p')]
sds_dims_name<-add_sds_name[c(1:5,11)];sds_sbd_name<-add_sds_name[13:15];sds_da_name<-add_sds_name[c(17,20,21)];sds_swtd_name<-add_sds_name[c(6:8,12,18,19)];sds_does_name<-add_sds_name[22:26];sds_shy_name<-add_sds_name[c(9,16)];

dataset$sds_dims <- rowSums(dataset[sds_dims_name],na.rm=TRUE)
dataset$sds_sbd <- rowSums(dataset[sds_sbd_name],na.rm=TRUE)
dataset$sds_da <- rowSums(dataset[sds_da_name],na.rm=TRUE)
dataset$sds_swtd <- rowSums(dataset[sds_swtd_name],na.rm=TRUE)
dataset$sds_does <- rowSums(dataset[sds_does_name],na.rm=TRUE)
dataset$sds_shy <- rowSums(dataset[sds_shy_name],na.rm=TRUE)
sds_subscale<-names(dataset)[which(names(dataset)=='sds_dims'):which(names(dataset)=='sds_shy')]
dataset$total.sleep.problems<-rowSums(dataset[sds_subscale],na.rm=TRUE)

##ABCD Pearson Scores
#Rey Auditory Verbal Learning Test, Matrix Reasoning Test, and Rey Delayed Recall Test
ravlt_verbal_learn_v<-c('pea_ravlt_sd_trial_i_tc','pea_ravlt_sd_trial_ii_tc','pea_ravlt_sd_trial_iii_tc','pea_ravlt_sd_trial_iv_tc','pea_ravlt_sd_trial_v_tc')
dataset$ravlt_verbal_learning<-rowSums(dataset[ravlt_verbal_learn_v],na.rm=TRUE)


#ABCD Little Man Task Summary Scores
#LMT Efficiency	=	Percentage correct divided by average reaction for correct trials
dataset$lmt_scr_efficiency<-dataset$lmt_scr_perc_correct / dataset$lmt_scr_rt_correct



all_variables <- read_excel("/share/inspurStorage/home1/zixing/ABCDR/Project_1/xlsx/covariates_x2.xlsx",1,col_names = TRUE)
all_variables=t(data.frame(all_variables))
dummy_v<-all_variables['dummy',][complete.cases(all_variables['dummy',])]
name_all<-all_variables['name_variable_in_abcdV4',][complete.cases(all_variables['name_variable_in_abcdV4',])]

ind_dummy<-c()
for(i in 1:length(dummy_v)){
  ind_dummy<-c(ind_dummy,which(names(dataset)==dummy_v[i]))
}
dataset[,ind_dummy]<-lapply(dataset[,ind_dummy],as.factor)
continuous_v<-name_all[!name_all %in% dummy_v]
# ind_continuous<-c()
# for(i in 1:length(continuous_v)){
#   ind_continuous<-c(ind_continuous,which(names(dataset)==continuous_v[i]))
# }

name_continuous<-continuous_v

setwd('/share/inspurStorage/home1/zixing/ABCDR/Project_1')
saveRDS(dataset, file = "dataset_beforewin.RDS")

dataset_beforewin<-readRDS("dataset_beforewin.RDS")
dataset_win<-dataset_beforewin
# ##winsorizaion- 0.6% of data  to move from the top and bottom of the distributions, no scores were ± 3 SDs
for (i in name_continuous){
  dataset_win[,i]=winsor(dataset_beforewin[,i], trim = 0.006,na.rm=TRUE)
}

ind_cbcl<-c()
for (i in name_continuous[which(name_continuous=="cbcl_scr_syn_anxdep_r"):which(name_continuous=="cbcl_scr_07_stress_r")]){
  ind_cbcl<-c(ind_cbcl,which(names(dataset_win)==i))
}

name_cbcl<-c('CBCL.anxious.or.depressed.score','CBCL.withdrawn.or.depressed.score','CBCL.somatic.complaints','CBCL.social.problems.score','CBCL.thought.problems','CBCL.attention.problems.score','CBCL.rule.breaking.behavior.score','CBCL.aggressive.behavior.score','CBCL.internalizing.factors','CBCL.externalizing.factors','CBCL.total.problems.score','CBCL.Depress','CBCL.AnxDisord','CBCL.SomaticPro','CBCL.ADHD','CBCL.Opposite','CBCL.Conduct','CBCL.Sluggish.Cognitive.Tempo(SCT)','CBCL.Obsessive.Compulsive.Problems(OCD)','CBCL.Stress')

for (i in 1:length(name_cbcl)){
  names(dataset_win)[ind_cbcl[i]]=name_cbcl[i]
}
name_continuous[which(name_continuous=="cbcl_scr_syn_anxdep_r"):which(name_continuous=="cbcl_scr_07_stress_r")]=name_cbcl

library('readxl')
confounding_v<-read_excel("/share/inspurStorage/home1/zixing/ABCDR/Project_1/xlsx/covariates_x2.xlsx",2,col_names = TRUE)
confounding_v=t(data.frame(confounding_v))
name_pes_f<-confounding_v['pes_f',][complete.cases(confounding_v['pes_f',])]
name_pes_p<-confounding_v['pes_p',][complete.cases(confounding_v['pes_p',])]
name_pes_s<-confounding_v['pes_s',][complete.cases(confounding_v['pes_s',])]
name_pes <- c(name_pes_f,name_pes_p,name_pes_s)
name_genetic<-confounding_v['family.history.problem',][complete.cases(confounding_v['family.history.problem',])]
name_additional<-confounding_v['additional',][complete.cases(confounding_v['additional',])]

outcomes_v<-read_excel("/share/inspurStorage/home1/zixing/ABCDR/Project_1/xlsx/covariates_x2.xlsx",3,col_names = TRUE)
outcomes_v=t(data.frame(outcomes_v))
outcomes_v<-outcomes_v['name_variable_in_abcdV4',][complete.cases(outcomes_v['name_variable_in_abcdV4',])]
nonmri_outcomes_v<-outcomes_v[1:147]

brain_outcomes_v<-outcomes_v[!outcomes_v %in% nonmri_outcomes_v]
outcomes_cortical<-brain_outcomes_v[!(brain_outcomes_v %in% brain_outcomes_v[c(72:101,103:115)])]


## Select exposures
exposures <- read_excel("/share/inspurStorage/home1/zixing/ABCDR/Project_1/xlsx/exposures_x2.xlsx",1, col_names = TRUE)
exposures=t(data.frame(exposures));ind_expo <- c()
for(i in 1:length(exposures['Name',])){
  ind_expo = c(ind_expo,c(which(names(dataset_win)==exposures['Name',][i])))
}
exposure_name<-names(dataset_win)[ind_expo];summary(dataset_win[exposure_name])

levels(dataset_win$caffeine)<-c('FALSE','TRUE');summary(dataset_win[exposure_name])
dataset_win$sex<-as.factor(dataset_win$sex)

dataset_win$race_ethnicity[dataset_win$race_ethnicity==1]<-'White';dataset_win$race_ethnicity[dataset_win$race_ethnicity==2]<-'Black';dataset_win$race_ethnicity[dataset_win$race_ethnicity==3]<-'Hispanic';dataset_win$race_ethnicity[dataset_win$race_ethnicity==4]<-'Asian';dataset_win$race_ethnicity[dataset_win$race_ethnicity==5]<-'Other'
dataset_win$race_ethnicity<-as.factor(dataset_win$race_ethnicity)
dataset_win$demo_prnt_marital_v2<-as.factor(dataset_win$demo_prnt_marital_v2)
dataset_win$demo_prnt_marital_living<-ifelse(dataset_win$demo_prnt_marital_v2 %in% c(1,6),1,0)
dataset_win$demo_prnt_marital_living<-as.factor(dataset_win$demo_prnt_marital_living)

dataset_win$rel_family_id<-as.factor(dataset_win$rel_family_id)
saveRDS(dataset_win, file = "dataset_win.RDS")

name_fixed<-c('sex','race_ethnicity','interview_age')
