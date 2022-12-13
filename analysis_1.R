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
library(RColorBrewer)
library(MuMIn)
library(lme4)
library(sjstats)
library(effectsize)
setwd('/share/inspurStorage/home1/zixing/ABCDR/Project_1')
project1_dataset_baseline_rmNA_win<-readRDS("dataset_baseline_rmNA_smriqc.RDS")

#selecting nesting variables
ind_nest = c(which(names(project1_dataset_baseline_rmNA_win)=="rel_family_id"), which(names(project1_dataset_baseline_rmNA_win)=="site_id_l"));summary(project1_dataset_baseline_rmNA_win[,ind_nest])

ind_fixed = c(which(names(project1_dataset_baseline_rmNA_win)=="sex"), which(names(project1_dataset_baseline_rmNA_win)=="race_ethnicity"));summary(project1_dataset_baseline_rmNA_win[,ind_fixed])


levels (project1_dataset_baseline_rmNA_win$caffeine)<-c('FALSE','TRUE')


project1_dataset_baseline_rmNA_win[,exposure_name]<-lapply(project1_dataset_baseline_rmNA_win[,exposure_name], as.logical)
project1_dataset_baseline_rmNA_win[,exposure_name]<-lapply(project1_dataset_baseline_rmNA_win[,exposure_name], as.numeric)
saveRDS(project1_dataset_baseline_rmNA_win, file = "project1_dataset_baseline_rmNA_win.RDS")


stdCoef.merMod <- function(object) {
  sdy <- sd(getME(object,"y"))
  sdx <- apply(getME(object,"X"), 2, sd)
  sc <- fixef(object)*sdx/sdy
  se.fixef <- coef(summary(object))[,"Std. Error"]
  se <- se.fixef*sdx/sdy
  return(data.frame(stdcoef=sc, stdse=se))
}

#########################################
## Select & process data for analyses  ##
#########################################

analysis1_outcome<-rep()
analysis1_exposure<-rep()
analysis1_out_beta<-rep()
analysis1_out_tvalue <- rep()
analysis1_out_pvalue <- rep()
analysis1_out_df <- rep()
analysis1_stdCoef<-rep()
analysis1_r2m<-rep()
analysis1_r2c<-rep()
analysis1_CI_low<-rep()
analysis1_CI_high<-rep()
analysis1_etasq<-rep()
analysis1_number=1


for (i in exposure_name){
  for(j in exposure_name){
    if(i!=j){
      project1_dataset_baseline_rmNA_win$y = as.numeric(unlist(project1_dataset_baseline_rmNA_win[,i]))
      project1_dataset_baseline_rmNA_win$x = as.numeric(unlist(project1_dataset_baseline_rmNA_win[,j]))
      lmer.fit <- lmer(y~x+(1 | site_id_l/rel_family_id),data = project1_dataset_baseline_rmNA_win)
      ci<-confint(lmer.fit,method='Wald')
      analysis1=summary(lmer.fit)$coefficients
      analysis1_out_tvalue[analysis1_number] = analysis1['x',4]
      analysis1_out_pvalue[analysis1_number] = analysis1['x',5]
      analysis1_out_df[analysis1_number] = analysis1['x',3]
      analysis1_out_beta[analysis1_number] = analysis1['x',1]
      analysis1_stdCoef[analysis1_number] = stdCoef.merMod(lmer.fit)['x',1]
      analysis1_r2c[analysis1_number] = as.numeric(performance::r2(lmer.fit)[1])
      analysis1_r2m[analysis1_number] = as.numeric(performance::r2(lmer.fit)[2])
      analysis1_etasq[analysis1_number] = as.numeric(effectsize::eta_squared(lmer.fit)[2])
      analysis1_CI_low[analysis1_number] = as.numeric(ci['x',1])
      analysis1_CI_high[analysis1_number] = as.numeric(ci['x',2])
      analysis1_outcome[analysis1_number] = i
      analysis1_exposure[analysis1_number] = j
      analysis1_number = analysis1_number + 1
    }
    else{
      analysis1=summary(lmer.fit)$coefficients
      analysis1_out_tvalue[analysis1_number] = 0
      analysis1_out_pvalue[analysis1_number] = 0
      analysis1_out_df[analysis1_number] = 0
      analysis1_out_beta[analysis1_number] = 0
      analysis1_stdCoef[analysis1_number] = 0
      analysis1_r2c[analysis1_number] = 0
      analysis1_r2m[analysis1_number] = 0
      analysis1_etasq[analysis1_number] = 0
      analysis1_CI_low[analysis1_number] = 0
      analysis1_CI_high[analysis1_number] = 0
      analysis1_outcome[analysis1_number] = i
      analysis1_exposure[analysis1_number] = j
      analysis1_number = analysis1_number + 1
    }
  }
}
results_analysis1 = data.frame(analysis1_outcome, analysis1_exposure, analysis1_out_tvalue, analysis1_out_beta, analysis1_stdCoef, analysis1_r2c,analysis1_r2m,analysis1_etasq,analysis1_CI_low,analysis1_CI_high,analysis1_out_pvalue,analysis1_out_df)

write.csv(results_analysis1,'/share/inspurStorage/home1/zixing/ABCDR/Project_1/new_results/analysis1/results_analysis1.csv');cat("results_analysis1.csv is exported can ready to download\n")

#sort descending and simplify the association matrix
results_pairwise_substance_use<-results_analysis1[,c(1,2,5,8,11)]


##### prepost data ######
pre_beforemeancentered<-project1_dataset_baseline_rm_NA_phillip_win_prs2norm_beforemeancentered_drop
pre_beforemeancentered$alcohol<-ifelse(project1_dataset_baseline_rm_NA_phillip_win_prs2norm_beforemeancentered_drop$devhx_9_alcohol.x==1,NA,ifelse(project1_dataset_baseline_rm_NA_phillip_win_prs2norm_beforemeancentered_drop$devhx_8_alcohol.x==1,'Pre','Nonexposed'))
pre_beforemeancentered$tobacco<-ifelse(project1_dataset_baseline_rm_NA_phillip_win_prs2norm_beforemeancentered_drop$devhx_9_tobacco.x==1,NA,ifelse(project1_dataset_baseline_rm_NA_phillip_win_prs2norm_beforemeancentered_drop$devhx_8_tobacco.x==1,'Pre','Nonexposed'))
pre_beforemeancentered$marijuana<-ifelse(project1_dataset_baseline_rm_NA_phillip_win_prs2norm_beforemeancentered_drop$devhx_9_marijuana.x==1,NA,ifelse(project1_dataset_baseline_rm_NA_phillip_win_prs2norm_beforemeancentered_drop$devhx_8_marijuana.x==1,'Pre','Nonexposed'))


pre_beforemeancentered[,exposure_name]<-lapply(pre_beforemeancentered[,exposure_name],as.factor)
pre_beforemeancentered$alcohol <- factor(pre_beforemeancentered$alcohol, levels = c('Nonexposed',"Pre"))
pre_beforemeancentered$tobacco <- factor(pre_beforemeancentered$tobacco, levels = c('Nonexposed',"Pre"))
pre_beforemeancentered$marijuana <- factor(pre_beforemeancentered$marijuana, levels = c('Nonexposed',"Pre"))

post_beforemeancentered<-project1_dataset_baseline_rm_NA_phillip_win_prs2norm_beforemeancentered_drop
post_beforemeancentered$alcohol<-ifelse(project1_dataset_baseline_rm_NA_phillip_win_prs2norm_beforemeancentered_drop$devhx_9_alcohol.x==1,'Post',ifelse(project1_dataset_baseline_rm_NA_phillip_win_prs2norm_beforemeancentered_drop$devhx_8_alcohol.x==1,NA,'Nonexposed'))
post_beforemeancentered$tobacco<-ifelse(project1_dataset_baseline_rm_NA_phillip_win_prs2norm_beforemeancentered_drop$devhx_9_tobacco.x==1,'Post',ifelse(project1_dataset_baseline_rm_NA_phillip_win_prs2norm_beforemeancentered_drop$devhx_8_tobacco.x==1,NA,'Nonexposed'))
post_beforemeancentered$marijuana<-ifelse(project1_dataset_baseline_rm_NA_phillip_win_prs2norm_beforemeancentered_drop$devhx_9_marijuana.x==1,'Post',ifelse(project1_dataset_baseline_rm_NA_phillip_win_prs2norm_beforemeancentered_drop$devhx_8_marijuana.x==1,NA,'Nonexposed'))


post_beforemeancentered[,exposure_name]<-lapply(post_beforemeancentered[,exposure_name],as.factor)
post_beforemeancentered$alcohol <- factor(post_beforemeancentered$alcohol, levels = c('Nonexposed',"Post"))
post_beforemeancentered$tobacco <- factor(post_beforemeancentered$tobacco, levels = c('Nonexposed',"Post"))
post_beforemeancentered$marijuana <- factor(post_beforemeancentered$marijuana, levels = c('Nonexposed',"Post"))


######### pre #########

analysis1_outcome<-rep()
analysis1_exposure<-rep()
analysis1_out_beta<-rep()
analysis1_out_tvalue <- rep()
analysis1_out_pvalue <- rep()
analysis1_out_df <- rep()
analysis1_stdCoef<-rep()
analysis1_r2m<-rep()
analysis1_r2c<-rep()
analysis1_CI_low<-rep()
analysis1_CI_high<-rep()
analysis1_etasq<-rep()
analysis1_number=1


for (i in exposure_name){
  for(j in exposure_name){
    if(i!=j){
      pre_beforemeancentered$y = as.numeric(unlist(pre_beforemeancentered[,i]))
      pre_beforemeancentered$x = as.numeric(unlist(pre_beforemeancentered[,j]))
      lmer.fit <- lmer(y~x+(1 | site_id_l/rel_family_id),data = pre_beforemeancentered)
      ci<-confint(lmer.fit,method='Wald')
      analysis1=summary(lmer.fit)$coefficients
      analysis1_out_tvalue[analysis1_number] = analysis1['x',4]
      analysis1_out_pvalue[analysis1_number] = analysis1['x',5]
      analysis1_out_df[analysis1_number] = analysis1['x',3]
      analysis1_out_beta[analysis1_number] = analysis1['x',1]
      analysis1_stdCoef[analysis1_number] = stdCoef.merMod(lmer.fit)['x',1]
      analysis1_r2c[analysis1_number] = as.numeric(performance::r2(lmer.fit)[1])
      analysis1_r2m[analysis1_number] = as.numeric(performance::r2(lmer.fit)[2])
      analysis1_etasq[analysis1_number] = as.numeric(effectsize::eta_squared(lmer.fit)[2])
      analysis1_CI_low[analysis1_number] = as.numeric(ci['x',1])
      analysis1_CI_high[analysis1_number] = as.numeric(ci['x',2])
      analysis1_outcome[analysis1_number] = i
      analysis1_exposure[analysis1_number] = j
      analysis1_number = analysis1_number + 1
    }
    else{
      analysis1=summary(lmer.fit)$coefficients
      analysis1_out_tvalue[analysis1_number] = 0
      analysis1_out_pvalue[analysis1_number] = 0
      analysis1_out_df[analysis1_number] = 0
      analysis1_out_beta[analysis1_number] = 0
      analysis1_stdCoef[analysis1_number] = 0
      analysis1_r2c[analysis1_number] = 0
      analysis1_r2m[analysis1_number] = 0
      analysis1_etasq[analysis1_number] = 0
      analysis1_CI_low[analysis1_number] = 0
      analysis1_CI_high[analysis1_number] = 0
      analysis1_outcome[analysis1_number] = i
      analysis1_exposure[analysis1_number] = j
      analysis1_number = analysis1_number + 1
    }
  }
}
results_analysis1_pre = data.frame(analysis1_outcome, analysis1_exposure, analysis1_out_tvalue, analysis1_out_beta, analysis1_stdCoef, analysis1_r2c,analysis1_r2m,analysis1_etasq,analysis1_CI_low,analysis1_CI_high,analysis1_out_pvalue,analysis1_out_df)

write.csv(results_analysis1_pre,'/share/inspurStorage/home1/zixing/ABCDR/Project_1/new_results/analysis1/results_analysis1_pre.csv');cat("results_analysis1_pre.csv is exported can ready to download\n")

#sort descending and simplify the association matrix
results_pairwise_presubstance_use<-results_analysis1_pre[,c(1,2,5,8,11)]
######### post #########

analysis1_outcome<-rep()
analysis1_exposure<-rep()
analysis1_out_beta<-rep()
analysis1_out_tvalue <- rep()
analysis1_out_pvalue <- rep()
analysis1_out_df <- rep()
analysis1_stdCoef<-rep()
analysis1_r2m<-rep()
analysis1_r2c<-rep()
analysis1_CI_low<-rep()
analysis1_CI_high<-rep()
analysis1_etasq<-rep()
analysis1_number=1


for (i in exposure_name){
  for(j in exposure_name){
    if(i!=j){
      post_beforemeancentered$y = as.numeric(unlist(post_beforemeancentered[,i]))
      post_beforemeancentered$x = as.numeric(unlist(post_beforemeancentered[,j]))
      lmer.fit <- lmer(y~x+(1 | site_id_l/rel_family_id),data = post_beforemeancentered)
      ci<-confint(lmer.fit,method='Wald')
      analysis1=summary(lmer.fit)$coefficients
      analysis1_out_tvalue[analysis1_number] = analysis1['x',4]
      analysis1_out_pvalue[analysis1_number] = analysis1['x',5]
      analysis1_out_df[analysis1_number] = analysis1['x',3]
      analysis1_out_beta[analysis1_number] = analysis1['x',1]
      analysis1_stdCoef[analysis1_number] = stdCoef.merMod(lmer.fit)['x',1]
      analysis1_r2c[analysis1_number] = as.numeric(performance::r2(lmer.fit)[1])
      analysis1_r2m[analysis1_number] = as.numeric(performance::r2(lmer.fit)[2])
      analysis1_etasq[analysis1_number] = as.numeric(effectsize::eta_squared(lmer.fit)[2])
      analysis1_CI_low[analysis1_number] = as.numeric(ci['x',1])
      analysis1_CI_high[analysis1_number] = as.numeric(ci['x',2])
      analysis1_outcome[analysis1_number] = i
      analysis1_exposure[analysis1_number] = j
      analysis1_number = analysis1_number + 1
    }
    else{
      analysis1=summary(lmer.fit)$coefficients
      analysis1_out_tvalue[analysis1_number] = 0
      analysis1_out_pvalue[analysis1_number] = 0
      analysis1_out_df[analysis1_number] = 0
      analysis1_out_beta[analysis1_number] = 0
      analysis1_stdCoef[analysis1_number] = 0
      analysis1_r2c[analysis1_number] = 0
      analysis1_r2m[analysis1_number] = 0
      analysis1_etasq[analysis1_number] = 0
      analysis1_CI_low[analysis1_number] = 0
      analysis1_CI_high[analysis1_number] = 0
      analysis1_outcome[analysis1_number] = i
      analysis1_exposure[analysis1_number] = j
      analysis1_number = analysis1_number + 1
    }
  }
}
results_analysis1_post= data.frame(analysis1_outcome, analysis1_exposure, analysis1_out_tvalue, analysis1_out_beta, analysis1_stdCoef, analysis1_r2c,analysis1_r2m,analysis1_etasq,analysis1_CI_low,analysis1_CI_high,analysis1_out_pvalue,analysis1_out_df)

write.csv(results_analysis1_post,'/share/inspurStorage/home1/zixing/ABCDR/Project_1/new_results/analysis1/results_analysis1_post.csv');cat("results_analysis1_post.csv is exported can ready to download\n")

#sort descending and simplify the association matrix
results_pairwise_postsubstance_use<-results_analysis1_post[,c(1,2,5,8,11)]
##############################################################
############draw the association from the strongest###########
##############################################################
library(ggplot2)
library(viridis)
library(psych)
library(plotly)

pairwise_substance_use.data <- results_pairwise_substance_use$analysis1_stdCoef
pairwise_substance_use=matrix(data = pairwise_substance_use.data, nrow =4, ncol = 4)
pairwise_substance_use<-t(pairwise_substance_use)
colnames(pairwise_substance_use) =exposure_name;rownames(pairwise_substance_use) =exposure_name

#do this before the transformation!
pairwise_substance_use[lower.tri(pairwise_substance_use, diag = TRUE)] <- NA
pairwise_substance_use <- pairwise_substance_use[-1, -ncol(pairwise_substance_use)]

#Store our variable names for later use
x_labels <- colnames(pairwise_substance_use)
y_labels <- rownames(pairwise_substance_use)
#Change the variable names to numeric for the grid
colnames(pairwise_substance_use) <- 1:ncol(pairwise_substance_use)
rownames(pairwise_substance_use) <- nrow(pairwise_substance_use):1


library(reshape)
#Melt the data into the desired format
plotdata <- melt(pairwise_substance_use)
colnames(plotdata)[3]<-'β'
# fig <- plot_ly(data = plotdata, width = 500, height = 500)
# fig <- fig %>% add_trace(x = ~Var2, y = ~Var1, type = "scatter",   mode = "scatter", color = ~etasq, symbol = I("square"))

#Adding the size variable & scaling it
plotdata$size <-(abs(plotdata$β))
scaling <- 1000 / ncol(pairwise_substance_use) 
plotdata$size <- plotdata$size * scaling

fig <- plot_ly(data = plotdata, width = 500, height = 500)
fig <- fig %>% add_trace(x = ~X2, y = ~X1, type = "scatter", mode = "markers", color = ~value, symbol = I("square"))


xAx1 <- list(showgrid = FALSE,
             showline = FALSE,
             zeroline = FALSE,
             tickvals = colnames(pairwise_substance_use),
             ticktext = x_labels,
             title = FALSE)
yAx1 <- list(autoaxis = FALSE,
             showgrid = FALSE,
             showline = FALSE,
             zeroline = FALSE,
             tickvals = rownames(pairwise_substance_use),
             ticktext = y_labels,
             title = FALSE)

fig <- fig %>% layout(xaxis = xAx1,
                      yaxis = yAx1,
                      plot_bgcolor = "rgba(0,0,0,0)",
                      paper_bgcolor = "rgba(0, 0, 0, 0.01)")

# fig <- plot_ly(data = plotdata, width = 500, height = 500)
# fig <- fig %>% add_trace(x = ~Var2, y = ~Var1, type = "scatter", mode = "markers", color = ~etasq, marker = list(size = ~size, opacity = 1), symbol = I("square"))
# # fig <- fig %>% layout(xaxis = xAx1, yaxis = yAx1)
# fig <- fig %>% layout(xaxis = xAx1,
#                       yaxis = yAx1,
#                       plot_bgcolor = "rgba(0,0,0,0)",
#                       paper_bgcolor = "rgba(0, 0, 0, 0.01)")
# fig <- fig %>% colorbar(limits = c(0,0.4), x = 0.95, y = 1)

Dataexposurealcohol<-project1_dataset_baseline_rm_NA_phillip_win_prs2norm_beforemeancentered[,c("src_subject_id",'alcohol')]
Dataexposurealcohol<-Dataexposurealcohol[Dataexposurealcohol$alcohol==TRUE,]
Dataexposurealcohol$Knowledge<-'Exposure'
Dataexposurealcohol$PSE<-'alcohol'
Dataexposurecoffee<-project1_dataset_baseline_rm_NA_phillip_win_prs2norm_beforemeancentered[,c("src_subject_id",'caffeine')]
Dataexposurecoffee<-Dataexposurecoffee[Dataexposurecoffee$caffeine==TRUE,]
Dataexposurecoffee$Knowledge<-'Exposure'
Dataexposurecoffee$PSE<-'coffee'
Dataexposuremarijuana<-project1_dataset_baseline_rm_NA_phillip_win_prs2norm_beforemeancentered[,c("src_subject_id",'marijuana')]
Dataexposuremarijuana<-Dataexposuremarijuana[Dataexposuremarijuana$marijuana==TRUE,]
Dataexposuremarijuana$Knowledge<-'Exposure'
Dataexposuremarijuana$PSE<-'marijuana'
Dataexposuretobacco<-project1_dataset_baseline_rm_NA_phillip_win_prs2norm_beforemeancentered[,c("src_subject_id",'tobacco')]
Dataexposuretobacco<-Dataexposuretobacco[Dataexposuretobacco$tobacco==TRUE,]
Dataexposuretobacco$Knowledge<-'Exposure'
Dataexposuretobacco$PSE<-'tobacco'
Dataexposurealcohol<-Dataexposurealcohol[complete.cases(Dataexposurealcohol), ] 
Dataexposurecoffee<-Dataexposurecoffee[complete.cases(Dataexposurecoffee), ] 
Dataexposuremarijuana<-Dataexposuremarijuana[complete.cases(Dataexposuremarijuana), ] 
Dataexposuretobacco<-Dataexposuretobacco[complete.cases(Dataexposuretobacco), ] 


Dataexposurealcohol<-Dataexposurealcohol[,-2]
Dataexposurecoffee<-Dataexposurecoffee[,-2]
Dataexposuremarijuana<-Dataexposuremarijuana[,-2]
Dataexposuretobacco<-Dataexposuretobacco[,-2]


Dataexposurealcoholpre<-project1_dataset_baseline_rm_NA_phillip_win_prs2norm_beforemeancentered[,c("src_subject_id",'devhx_8_alcohol.x')]
Dataexposurealcoholpre<-Dataexposurealcoholpre[Dataexposurealcoholpre$devhx_8_alcohol.x==1,]
Dataexposurealcoholpre$Knowledge<-'Before knowledge'
Dataexposurealcoholpre$PSE<-'alcohol'
Dataexposuremarijuanapre<-project1_dataset_baseline_rm_NA_phillip_win_prs2norm_beforemeancentered[,c("src_subject_id",'devhx_8_marijuana.x')]
Dataexposuremarijuanapre<-Dataexposuremarijuanapre[Dataexposuremarijuanapre$devhx_8_marijuana.x==1,]
Dataexposuremarijuanapre$Knowledge<-'Before knowledge'
Dataexposuremarijuanapre$PSE<-'marijuana'
Dataexposuretobaccopre<-project1_dataset_baseline_rm_NA_phillip_win_prs2norm_beforemeancentered[,c("src_subject_id",'devhx_8_tobacco.x')]
Dataexposuretobaccopre<-Dataexposuretobaccopre[Dataexposuretobaccopre$devhx_8_tobacco.x==1,]
Dataexposuretobaccopre$Knowledge<-'Before knowledge'
Dataexposuretobaccopre$PSE<-'tobacco'
Dataexposurealcoholpre<-Dataexposurealcoholpre[complete.cases(Dataexposurealcoholpre), ] 
Dataexposuremarijuanapre<-Dataexposuremarijuanapre[complete.cases(Dataexposuremarijuanapre), ] 
Dataexposuretobaccopre<-Dataexposuretobaccopre[complete.cases(Dataexposuretobaccopre), ] 

Dataexposurealcoholpre<-Dataexposurealcoholpre[,-2]
Dataexposuremarijuanapre<-Dataexposuremarijuanapre[,-2]
Dataexposuretobaccopre<-Dataexposuretobaccopre[,-2]



Dataexposurealcoholpost<-project1_dataset_baseline_rm_NA_phillip_win_prs2norm_beforemeancentered[,c("src_subject_id",'devhx_9_alcohol.x')]
Dataexposurealcoholpost<-Dataexposurealcoholpost[Dataexposurealcoholpost$devhx_9_alcohol.x==1,]
Dataexposurealcoholpost$Knowledge<-'After knowledge'
Dataexposurealcoholpost$PSE<-'alcohol'
Dataexposurealcoholpost<-Dataexposurealcoholpost[complete.cases(Dataexposurealcoholpost), ] 
Dataexposuremarijuanapost<-project1_dataset_baseline_rm_NA_phillip_win_prs2norm_beforemeancentered[,c("src_subject_id",'devhx_9_marijuana.x')]
Dataexposuremarijuanapost<-Dataexposuremarijuanapost[Dataexposuremarijuanapost$devhx_9_marijuana.x==1,]
Dataexposuremarijuanapost$Knowledge<-'After knowledge'
Dataexposuremarijuanapost$PSE<-'marijuana'
Dataexposuremarijuanapost <- Dataexposuremarijuanapost[complete.cases(Dataexposuremarijuanapost), ]
Dataexposuretobaccopost<-project1_dataset_baseline_rm_NA_phillip_win_prs2norm_beforemeancentered[,c("src_subject_id",'devhx_9_tobacco.x')]
Dataexposuretobaccopost<-Dataexposuretobaccopost[Dataexposuretobaccopost$devhx_9_tobacco.x==1,]
Dataexposuretobaccopost$Knowledge<-'After knowledge'
Dataexposuretobaccopost$PSE<-'tobacco'
Dataexposuretobaccopost<-Dataexposuretobaccopost[complete.cases(Dataexposuretobaccopost), ] 

Dataexposurealcoholpost<-Dataexposurealcoholpost[,-2]
Dataexposuremarijuanapost<-Dataexposuremarijuanapost[,-2]
Dataexposuretobaccopost<-Dataexposuretobaccopost[,-2]



df_Dataprepare<-rbind(Dataexposurecoffee,Dataexposurealcohol,Dataexposurealcoholpre,Dataexposurealcoholpost,Dataexposuretobacco,Dataexposuretobaccopre,Dataexposuretobaccopost,Dataexposuremarijuana,Dataexposuremarijuanapre,Dataexposuremarijuanapost)

Dataexposurecaffeinepre<-Dataexposurealcoholpre
Dataexposurecaffeinepost<-Dataexposurealcoholpost
Dataexposurecaffeinepre$PSE<-'coffee'
Dataexposurecaffeinepost$PSE<-'coffee'


df_Dataprepare1<-rbind(Dataexposurecoffee,Dataexposurecaffeinepre,Dataexposurecaffeinepost,Dataexposurealcohol,Dataexposurealcoholpre,Dataexposurealcoholpost,Dataexposuretobacco,Dataexposuretobaccopre,Dataexposuretobaccopost,Dataexposuremarijuana,Dataexposuremarijuanapre,Dataexposuremarijuanapost)

# df_alcohol<-rbind(Dataexposurealcohol,Dataexposurealcoholpre,Dataexposurealcoholpost)
# df_coffee<-rbind(Dataexposurecoffee)
# df_marijuana<-rbind(Dataexposuremarijuana,Dataexposuremarijuanapre,Dataexposuremarijuanapost)
# df_tobacco<-rbind(Dataexposuretobacco,Dataexposuretobaccopre,Dataexposuretobaccopost)



  
df_Dataprepare$PSE<-factor(df_Dataprepare$PSE,levels = c('coffee','alcohol','tobacco','marijuana'))
df_Dataprepare$Knowledge<-factor(df_Dataprepare$Knowledge,levels = c('Exposure','Before knowledge','After knowledge'))

df_Dataprepare1$PSE<-factor(df_Dataprepare1$PSE,levels = c('coffee','alcohol','tobacco','marijuana'))
df_Dataprepare1$Knowledge<-factor(df_Dataprepare1$Knowledge,levels = c('Exposure','Before knowledge','After knowledge'))

# df_alcohol$PSE<-factor(df_alcohol$PSE,levels = c('alcohol','coffee','marijuana','tobacco'))
# df_alcohol$Knowledge<-factor(df_alcohol$Knowledge,levels = c('Exposure','Before knowledge','After knowledge'))
# df_coffee$PSE<-factor(df_coffee$PSE,levels = c('alcohol','coffee','marijuana','tobacco'))
# df_coffee$Knowledge<-factor(df_coffee$Knowledge,levels = c('Exposure','Before knowledge','After knowledge'))
# df_marijuana$PSE<-factor(df_marijuana$PSE,levels = c('alcohol','coffee','marijuana','tobacco'))
# df_marijuana$Knowledge<-factor(df_marijuana$Knowledge,levels = c('Exposure','Before knowledge','After knowledge'))
# df_tobacco$PSE<-factor(df_tobacco$PSE,levels = c('alcohol','coffee','marijuana','tobacco'))
# df_tobacco$Knowledge<-factor(df_tobacco$Knowledge,levels = c('Exposure','Before knowledge','After knowledge'))


install.packages('ggridges')
install.packages("tidyverse")

library("ggridges")
library("tidyverse")
ggplot(df_Dataprepare, aes(x = Knowledge)) +
  geom_histogram(fill = "white", colour = "black") +
  facet_grid(PSE ~ ., scales = "free")

ggplot(df_Dataprepare, aes(x=Knowledge)) +
  geom_bar() +
  facet_grid(PSE ~ ., scales = "free")

# Bar chart side by side
ggplot(df_Dataprepare, aes(x = Knowledge, fill = PSE)) +
  geom_bar(position = position_dodge(0.7), width = 0.6) +
  theme_classic()+scale_fill_manual(values = c("#7E6148BF", "#E64B35BF","#00A087BF", "#3C5488BF"))+theme(legend.position="bottom")

ggplot(df_Dataprepare, aes(x = Knowledge, fill = PSE)) +
  geom_bar(position = position_dodge(0.7), width = 0.6) +
  theme_classic()+scale_fill_manual(values = c(rgb(0.3,0.1,0.4,0.6) , rgb(0.3,0.5,0.4,0.6) , rgb(0.3,0.9,0.4,0.6) ,  rgb(0.3,0.9,0.4,0.6)))+theme(legend.position="bottom")

install.packages('viridis')
library(viridis)
ggplot(df_Dataprepare, aes(x = Knowledge, fill = PSE)) +
  geom_bar(position = position_dodge(0.7), width = 0.6) +
  theme_classic()+scale_fill_viridis(discrete=TRUE)+theme(legend.position="bottom")

ggplot(df_Dataprepare, aes(x = Knowledge, fill = PSE)) +
  geom_bar(position = position_dodge(0.7), width = 0.6) +
  theme_classic()+scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9", "#009E73"))+theme(legend.position="bottom")

install.packages('wesanderson')
library(wesanderson)
ggplot(df_Dataprepare, aes(x = Knowledge, fill = PSE)) +
  geom_bar(position = position_dodge(0.7), width = 0.6) +
  theme_classic()+scale_fill_manual(values = wes_palette("Moonrise2", n = 4))+theme(legend.position="bottom")

ggplot(df_Dataprepare, aes(x = Knowledge, fill = PSE)) +
  geom_bar( width = 0.6) +
  theme_classic()+scale_fill_manual(values = wes_palette("Moonrise2", n = 4))+theme(legend.position="bottom")

ggplot(df_Dataprepare1, aes(x = Knowledge, fill = PSE)) +
  geom_bar(position = position_dodge(0.7), width = 0.6) +
  theme_classic()+scale_fill_manual(values = wes_palette("Moonrise2", n = 4))+theme(legend.position="bottom")




# newdf_Dataprepare<-as.data.frame(rep(c('Alcohol','Coffee','Marijuana','Tobacco'),times=3))
# colnames(newdf_Dataprepare)[1]<-'PSE'
# newdf_Dataprepare$Knowledge<-c(rep('Exposure',4),rep('Before knowledge',4),rep('After knowledge',4))
# 
# newdf_Dataprepare$Prevalence<-c(2524/9838*100,5880/9838*100,547/9838*100,1300/9838*100,2488/9838*100,NA,541/9838*100,1289/9838*100,254/9838*100,NA,187/9838*100,473/9838*100)
# newdf_Dataprepare<-newdf_Dataprepare[c(-6,-10),]
# 
# newdf_Dataprepare$PSE<-factor(newdf_Dataprepare$PSE,levels = c('Alcohol','Coffee','Marijuana','Tobacco'))
# newdf_Dataprepare$Knowledge<-factor(newdf_Dataprepare$Knowledge,levels = c('Exposure','Before knowledge','After knowledge'))





