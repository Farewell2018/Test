getwd()
setwd("E:/ybw")
############Load data##############
library(data.table)
library(foreign)
data_all <- as.data.frame(fread('0329data_machine_final.csv')) #data_all pre all data
data_all <- data_all[,-1]
cov <- as.character(colnames(data_all))
fixcov <- cov[c(1:12,29)]
olicov <- cov[13:26]
data_uu <- data_all #data to calculate

data_all <- as.data.frame(read.spss("lipid_1.sav"))
data_uu <- data_all #data to calculate

data_all1 <- as.data.frame(read.spss("lipid_11.sav"))
data_uu$frmscore <- data_all1$frmrisknew
data_uu$s2d <- data_all1$S2Drisk_1
data_uu$ascvd <- data_all1$ASCVDrisk
#############建立模型##############
library(caret)
data_uu <- as.data.frame(data_uu)
data_uu <- filter(data_uu,followuptime>=1)
set.seed(79782)#2223
intrain <- createDataPartition(y = data_uu$incidentCVD,
                             p = 0.7,#改变比例
                             list = FALSE)
data_train <- data_uu[intrain,]  #划定训练集
data_test  <- data_uu[-intrain,] #划定测试集

colMeans(data_test)
apply(data_test, 2, sd)
summary(data_uu)
table(data_test$lipid)
write.csv(data_test,"data_test.csv")
#########With olink model#########
library(timeROC)
library(survival)
train_set  <- data_train[,c(olicov,covtest,"followuptime","incidentCVD")] #train set to calculate
test_set   <- data_test[,c(olicov,covtest,"followuptime","incidentCVD")]  #test set  to calculate

#cox
coxph_fit <- coxph(Surv(followuptime, incidentCVD) ~., data = train_set, model=T, x=T,y=T)
test_set$pred <- predict(coxph_fit, newdata = test_set,type ="lp")
tROC_cox<-timeROC(T=test_set$followuptime,delta=test_set$incidentCVD,marker=test_set$pred,weighting="marginal",
                  times=c(1:15), cause=1, ROC = TRUE, iid=T)
tROC_cox$AUC

#random forest
library(randomForestSRC)
tune.nodesize(Surv(followuptime,incidentCVD) ~ ., train_set)
n <- 5000
RF_fit<-tune.rfsrc(Surv(followuptime, incidentCVD==1) ~., data=train_set, ntreeTry = n)
{RF.final <- rfsrc(Surv(followuptime, incidentCVD==1) ~., data = train_set,
                  nodedepth = 3,
                  ntree = 4500,
                  mtry = 1,#RF_fit[["optimal"]][["mtry"]],#500.6.1.3
                  nodesize =1,#RF_fit[["optimal"]][["nodesize"]],#1
                  importance = T,
                  seed= 1234
                  )
RFS.test<- predict(RF.final, newdata = test_set,type="lp")
tmp.test <- test_set
tmp.test$pred <- RFS.test$predicted
tROC_RF<-timeROC(T=tmp.test$followuptime,delta=tmp.test$incidentCVD,marker=tmp.test$pred,weighting="marginal",
                 times=c(1:15), cause=1, ROC = TRUE, iid=TRUE)
plot(RF.final)
tROC_RF$AUC}
tROC_cox$AUC

tROC_FRM <- timeROC(T=data_test$followuptime,delta=data_test$incidentCVD,marker=data_test$frmscore,weighting="marginal",
                    times=c(1:15), cause=1, ROC = TRUE, iid=TRUE)
tROC_FRM$AUC

tROC_s2d <- timeROC(T=data_test$followuptime,delta=data_test$incidentCVD,marker=data_test$s2d,weighting="marginal",
                    times=c(1:15), cause=1, ROC = TRUE, iid=TRUE)
tROC_s2d$AUC

tROC_ascvd <- timeROC(T=data_test$followuptime,delta=data_test$incidentCVD,marker=data_test$ascvd,weighting="marginal",
                    times=c(1:15), cause=1, ROC = TRUE, iid=TRUE)
tROC_ascvd$AUC

saveRDS(RF.final,"oliRF0408.rds")#save model
#saveRDS(RF.final,"oliRF0401.rds")#save model 4500_1_25,seed=2223_233,p=0.7

library(survex)  
plot(RF.final)
explainer <- explain(RF.final, 
                     data =train_set[, -c(23, 24)],
                     y = Surv(train_set$followuptime, train_set$incidentCVD))
model_performance(explainer)
oli_imp_p <- model_parts(explainer)
plot(oli_imp_p)

########Without olink model#######
# "age"                "sex"                "bmi_1"              "smoke_1"            "drinking"          
# "sbp_1"              "hba1c_1"            "tc_1"               "diabete_medicine"   "chole_medicine"    
# "bloodpre_medicine"  "heartdiseasefamily" "lipid"
covtest <- c("age","sex","bmi_1","hba1c_1","sbp","lipid"
             ,"smoke_1","diabete_medicine"
             )
train_set_1<-data_train[,c(covtest,
                           "followuptime","incidentCVD")]#train_set_1 without olink
test_set_1<-data_test[,c(covtest,
                         "followuptime","incidentCVD")]#test_set_1  without olink
#1)cox
coxph_fit_1 <- coxph(Surv(followuptime, incidentCVD==1) ~., data = train_set_1, model=T, x=T,y=T)
test_set_1$pred <- predict(coxph_fit_1, newdata = test_set_1,type ="lp")#根据cox回归拟合测试集的结局
#计算time-dependent ROC
tROC_cox_1<-timeROC(T=test_set_1$followuptime,delta=test_set_1$incidentCVD,marker=test_set_1$pred,
                    weighting="marginal",times=c(1:15), cause=1, ROC = TRUE, iid=T)
tROC_cox_1$AUC

#2)RF
#randomForestSRC
n <- 8000
RF_fit_1<-tune.rfsrc(Surv(followuptime, incidentCVD==1) ~., data=train_set_1, ntreeTry = n)
{RF.final_1 <- rfsrc(Surv(followuptime, incidentCVD==1) ~., data = train_set_1,
                    nodedepth = 2,
                    ntree =6000,
                    mtry = 7,#RF_fit_1[["optimal"]][["mtry"]],#
                    nodesize = 30,#RF_fit_1[["optimal"]][["nodesize"]],#
                    importance = T,
                    seed=1234
                    )
RFS.test_1<- predict(RF.final_1, newdata = test_set_1)
tmp.test_1<- test_set_1
tmp.test_1$pred <- RFS.test_1$predicted
tROC_RF_1<-timeROC(T=tmp.test_1$followuptime,delta=tmp.test_1$incidentCVD,marker=tmp.test_1$pred,
                   weighting="marginal",times=c(1:15), cause=1, ROC = TRUE, iid=TRUE)
tROC_RF_1$AUC}
tROC_cox_1$AUC
plot(RF.final_1)

###确定模型参数进行保存
saveRDS(RF.final_1,"0409basicRF.rds")#save全部基础变量的模型6000 1 60 depth=2

library("magrittr")
library("dplyr")
library("tibble")
RFmodel_oli <- readRDS("oliRF0408.rds") #load RFmodel with olink
RFmodel_out <- readRDS("0409basicRF.rds")  #load without olink

plot(RFmodel_oli)
importance_olink <- data.frame(RFmodel_oli$importance) %>%
  rownames_to_column("cov") %>%
  arrange(- RFmodel_oli.importance)
importance_olink

library(stringr)
library(ggplot2)
library(viridis)
library(ggsci)
reorder(importance_olink,cov)
write.csv(as.data.frame(importance_olink),"0409ipm.csv")
test_imp <- read.csv("0409ipm.csv")

ggplot(data=test_imp, aes(x=reorder(cov,-order),
                          y=RFmodel_oli.importance, fill= RFmodel_oli.importance))+
  geom_bar(stat="identity")+
  theme_classic()+
  theme(legend.position ='none')+
  coord_flip()+ 
  scale_fill_viridis(begin = 1, end = 0.5, option = "G")
  # scale_fill_hue(h = c(0, 360), c = 80)#变量重要性排序作图

# importance_olink1 <- importance_olink[c(4:10,13:15,17:19),]
# olink_importance <- ggplot(data=importance_olink1, aes(x=reorder(cov,RFmodel_oli.importance),
#                                   y=RFmodel_oli.importance, fill=cov))+
#   geom_bar(stat="identity")+
#   theme_classic()+
#   theme(legend.position ='none')+
#   coord_flip()#变量重要性排序作图

##########Time dependent ROC######
library(timeROC)
library(survival)
library(ggsci)
library(cowplot)
library(irtoys)
library(randomForest)
dev.off()
RFS.test <- predict(RFmodel_oli, newdata = test_set)
tmp.test <- test_set
tmp.test$pred <- RFS.test$predicted
tROC_RF  <- timeROC(T=tmp.test$followuptime,delta=tmp.test$incidentCVD,marker=tmp.test$pred,weighting="marginal",
                    times=c(1:15), cause=1, ROC = TRUE, iid=TRUE)
tROC_RF$AUC#RF with olink time dependent ROC

RFS.test_1<- predict(RFmodel_out, newdata = test_set_1)
tmp.test_1<- test_set_1
tmp.test_1$pred <- RFS.test_1$predicted
tROC_RF_1<-timeROC(T=tmp.test_1$followuptime,delta=tmp.test_1$incidentCVD,marker=tmp.test_1$pred,weighting="marginal",
                   times=c(1:15), cause=1, ROC = TRUE, iid=TRUE)
tROC_RF_1$AUC#RF without olink time dependent ROC

sRocFuction = function(td = null){
  par(mar= c(5,5,1,1),cex.lab=1.2,cex.axis= 1.2)
  plot(td,time=3,lwd = 2, cex.main=1.3, cex.lab=1.5, cex.axis=1.2, font=1.2)
  plot(td,time=5,lwd = 2,add=TRUE,col="blue")
  plot(td,time=10,lwd = 2,add=TRUE,col="grey50")
  legend("bottomright",c(paste0("3   years  ", "(AUC=", sprintf("%0.2f",td$AUC[3]) ,")"),
                         paste0("5   years  ", "(AUC=", sprintf("%0.2f",td$AUC[5]) ,")"),
                         paste0("10 years "  ," (AUC=", sprintf("%0.2f",td$AUC[10]),")")
                         ),
  col=c("red","blue","grey50"),lty=1,lwd=2)
}
sRocFuction(td = tROC_RF)
Polink_all <- recordPlot()#保存为with_all

sRocFuction(td = tROC_RF_1)
Pnooli_all <- recordPlot()#保存为without_all

##########另一种做timeROC图###############
result <- with(tmp.test, timeROC(
  T = followuptime,
  delta = incidentCVD,
  marker = pred,
  cause = 1,
  weighting="marginal",
  times = c(3, 5, 10),
  ROC = TRUE,
  iid = TRUE
))

tROC_FRM$AUC
tROC_s2d$AUC
tROC_ascvd$AUC

result_frm <- with(data_test, timeROC(
  T = followuptime,
  delta = incidentCVD,
  marker = frmscore,
  cause = 1,
  weighting = "marginal",
  times = c(3,5,10),
  ROC = TRUE,
  iid = TRUE
))
result_frm$AUC
dat_frm = data.frame(fpr = as.numeric(result_frm$FP),
                   tpr = as.numeric(result_frm$TP),
                   time = rep(as.factor(c(3,5,10)),each = nrow(result_frm$TP)))

result_s2d <- with(data_test, timeROC(
  T = followuptime,
  delta = incidentCVD,
  marker = s2d,
  cause = 1,
  weighting = "marginal",
  times = c(3,5,10),
  ROC = TRUE,
  iid = TRUE
))
dat_s2d = data.frame(fpr = as.numeric(result_s2d$FP),
                     tpr = as.numeric(result_s2d$TP),
                     time = rep(as.factor(c(3,5,10)),each = nrow(result_s2d$TP)))

result_ascvd <- with(data_test, timeROC(
  T = followuptime,
  delta = incidentCVD,
  marker = ascvd,
  cause = 1,
  weighting = "marginal",
  times = c(3,5,10),
  ROC = TRUE,
  iid = TRUE
))
dat_ascvd = data.frame(fpr = as.numeric(result_ascvd$FP),
                     tpr = as.numeric(result_ascvd$TP),
                     time = rep(as.factor(c(3,5,10)),each = nrow(result_ascvd$TP)))

dat_1 = data.frame(fpr = as.numeric(result$FP),
                   tpr = as.numeric(result$TP),
                   time = rep(as.factor(c(3,5,10)),each = nrow(result$TP)))
p1<-ggplot() +
  geom_line(data = dat_1,aes(x = fpr, y = tpr,color = time),size = 1) +
  scale_color_manual(name = NULL,values = c("#92C5DE", "#F4A582", "#66C2A5"),
                     labels = paste0("AUC of ",c(3,5,10),"-y CVD incidence: ",
                                     format(round(result$AUC,3),nsmall = 2)))+
  geom_line(aes(x=c(0,1),y=c(0,1)),color = "grey")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.background = element_rect(linetype = 1, size = 0.2, colour = "black"),
        legend.position = c(0.765,0.125))+
  scale_x_continuous(expand = c(0.005,0.005))+
  scale_y_continuous(expand = c(0.005,0.005))+
  labs(x = "1 - Specificity",
       y = "Sensitivity")+
  coord_fixed()
p1

result_1 <- with(tmp.test_1, timeROC(
  T = followuptime,
  delta = incidentCVD,
  marker = pred,
  cause = 1,
  weighting="marginal",
  times = c(3,5, 10),
  ROC = TRUE,
  iid = TRUE
))

dat_2 = data.frame(fpr = as.numeric(result_1$FP),
                   tpr = as.numeric(result_1$TP),
                   time = rep(as.factor(c(3,5,10)),each = nrow(result_1$TP)))
p2<-ggplot() +
  geom_line(data = dat_2,aes(x = fpr, y = tpr,color = time),size = 1) +
  scale_color_manual(name = NULL,values = c("#92C5DE", "#F4A582", "#66C2A5"),
                     labels = paste0("AUC of ",c(3,5,10),"-y DM incidence: ",
                                     format(round(result_1$AUC,3),nsmall = 2)))+
  geom_line(aes(x=c(0,1),y=c(0,1)),color = "grey")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.background = element_rect(linetype = 1, size = 0.2, colour = "black"),
        legend.position = c(0.765,0.125))+
  scale_x_continuous(expand = c(0.005,0.005))+
  scale_y_continuous(expand = c(0.005,0.005))+
  labs(x = "1 - Specificity",
       y = "Sensitivity")+
  coord_fixed()
p2

#######TimeROC_with_vs_without_3_5_10##########
library(patchwork)
dat_1$model<-'Full  model'
dat_2$model<-'Basic model'
dat_3$model<-'FRM   model'
datt<-rbind(dat_1,dat_2)
datt$model<-factor(datt$model)
datt2<-rbind(dat_1,dat_3)
datt2$model<-factor(datt2$model)
dat_frm$model <- 'Framingham'
dat_s2d$model <- 'Score2-Diabetes'
dat_ascvd$model <-'ASCVD'
datt3 <- rbind(dat_1, dat_frm, dat_s2d, dat_ascvd)

vs_3<-ggplot() + 
  geom_line(data = subset(datt,time==3),aes(x = fpr, y = tpr,color = model),size = 1) + 
  scale_color_manual(name = NULL,values = c("blue", "red"),
                     labels = paste0(c('Basic model','Full    model'),": ",c(sprintf("%0.2f",tROC_RF_1$AUC[3]),sprintf("%0.2f",tROC_RF$AUC[3]))))+
  geom_line(aes(x=c(0,1),y=c(0,1)),color = "grey",linetype ="dashed")+
  theme_classic()+
  theme(panel.grid = element_blank(),
        #legend.background = element_rect(linetype = 1, size = 0.2, colour = "black"),
        legend.position = c(0.7,0.125), legend.text=element_text(size =16), axis.title = element_text(size = 16),
        axis.line = element_line(linewidth = 1.0, color = "black"),
        axis.ticks = element_line(size = 1.0),
        axis.ticks.length = unit(0.2, "cm"),
        axis.text = element_text(size = 15))+
  scale_x_continuous(expand = c(0.005,0.005))+
  scale_y_continuous(expand = c(0.005,0.005))+
  labs(x = "1 - Specificity",
       y = "Sensitivity")
vs_3
vs_5<-ggplot() + 
  geom_line(data = subset(datt,time==5),aes(x = fpr, y = tpr,color = model),size = 1) + 
  scale_color_manual(name = NULL,values = c("blue", "red"),
                     labels = paste0(c('Basic model','Full    model'),": ",c(sprintf("%0.2f",tROC_RF_1$AUC[5]),sprintf("%0.2f",tROC_RF$AUC[5]))))+
  geom_line(aes(x=c(0,1),y=c(0,1)),color = "grey",linetype ="dashed")+
  theme_classic()+
  theme(panel.grid = element_blank(),
        #legend.background = element_rect(linetype = 1, size = 0.2, colour = "black"),
        legend.position = c(0.7,0.125),legend.text=element_text(size =16), axis.title = element_text(size = 16),
        axis.line = element_line(linewidth = 1.0, color = "black"),
        axis.ticks = element_line(size = 1.0),
        axis.ticks.length = unit(0.2, "cm"),
        axis.text = element_text(size = 15))+
  scale_x_continuous(expand = c(0.005,0.005))+
  scale_y_continuous(expand = c(0.005,0.005))+
  labs(x = "1 - Specificity",
       y = "Sensitivity")
vs_5
vs_10<-ggplot() + 
  geom_line(data = subset(datt,time==10),aes(x = fpr, y = tpr,color = model),size = 1) + 
  scale_color_manual(name = NULL,values = c("blue", "red"),
                     labels = paste0(c('Basic model','Full    model'),": ",c(sprintf("%0.2f",tROC_RF_1$AUC[10]),sprintf("%0.2f",tROC_RF$AUC[10]))))+
  geom_line(aes(x=c(0,1),y=c(0,1)),color = "grey",linetype ="dashed")+
  theme_classic()+
  theme(panel.grid = element_blank(),
        #legend.background = element_rect(linetype = 1, size = 0.2, colour = "black"),
        legend.position = c(0.7,0.125),legend.text=element_text(size =16), axis.title = element_text(size = 16),
        axis.line = element_line(linewidth = 1.0, color = "black"),
        axis.ticks = element_line(size = 1.0),
        axis.ticks.length = unit(0.2, "cm"),
        axis.text = element_text(size = 15))+
  scale_x_continuous(expand = c(0.005,0.005))+
  scale_y_continuous(expand = c(0.005,0.005))+
  labs(x = "1 - Specificity",
       y = "Sensitivity")
vs_10

vs_3+vs_5+vs_10+plot_layout(ncol=3,nrow=1)+plot_annotation(tag_levels="A")&theme(plot.tag = element_text(size = 14,face = 'bold'))
dev.off()### RF_with_vs_without_3_5_10

dat_frm$model <- 'Framingham'
dat_s2d$model <- 'Score2-Diabetes'
dat_ascvd$model <-'ASCVD'

vs_3_1<-ggplot() + 
  geom_line(data = subset(datt3,time==3),aes(x = fpr, y = tpr,color = model),size = 1)+ 
  scale_color_manual(name = NULL,values = c("#14517C", "#D8383A","#63E398","#F3D266"),
                     labels = paste0(c('ASCVD','Framingham','Full model','Score2-Diabetes'),": ",
                                     c(sprintf("%0.2f",tROC_ascvd$AUC[3]),
                                       sprintf("%0.2f",tROC_FRM$AUC[3]),
                                       sprintf("%0.2f",tROC_RF$AUC[3]),
                                       sprintf("%0.2f",tROC_s2d$AUC[3]))))+
  geom_line(aes(x=c(0,1),y=c(0,1)),color = "grey",linetype ="dashed")+
  theme_classic()+
  theme(panel.grid = element_blank(),
        #legend.background = element_rect(linetype = 1, size = 0.2, colour = "black"),
        legend.position = c(0.7,0.125), legend.text=element_text(size =16), axis.title = element_text(size = 16),
        axis.line = element_line(linewidth = 1.0, color = "black"),
        axis.ticks = element_line(size = 1.0),
        axis.ticks.length = unit(0.2, "cm"),
        axis.text = element_text(size = 15))+
  scale_x_continuous(expand = c(0.005,0.005))+
  scale_y_continuous(expand = c(0.005,0.005))+
  labs(x = "1 - Specificity",
       y = "Sensitivity")
vs_3_1

vs_5_1<-ggplot() + 
  geom_line(data = subset(datt3,time==5),aes(x = fpr, y = tpr,color = model),size = 1)+ 
  scale_color_manual(name = NULL,values = c("#14517C", "#D8383A","#63E398","#F3D266"),
                     labels = paste0(c('ASCVD','Framingham','Full model','Score2-Diabetes'),": ",
                                     c(sprintf("%0.2f",tROC_ascvd$AUC[5]),
                                       sprintf("%0.2f",tROC_FRM$AUC[5]),
                                       sprintf("%0.2f",tROC_RF$AUC[5]),
                                       sprintf("%0.2f",tROC_s2d$AUC[5]))))+
  geom_line(aes(x=c(0,1),y=c(0,1)),color = "grey",linetype ="dashed")+
  theme_classic()+
  theme(panel.grid = element_blank(),
        #legend.background = element_rect(linetype = 1, size = 0.2, colour = "black"),
        legend.position = c(0.7,0.125),legend.text=element_text(size =16), axis.title = element_text(size = 16),
        axis.line = element_line(linewidth = 1.0, color = "black"),
        axis.ticks = element_line(size = 1.0),
        axis.ticks.length = unit(0.2, "cm"),
        axis.text = element_text(size = 15))+
  scale_x_continuous(expand = c(0.005,0.005))+
  scale_y_continuous(expand = c(0.005,0.005))+
  labs(x = "1 - Specificity",
       y = "Sensitivity")
vs_5_1
vs_10_1<-ggplot() + 
  geom_line(data = subset(datt3,time==10),aes(x = fpr, y = tpr,color = model),size = 1)+ 
  scale_color_manual(name = NULL,values = c("#14517C", "#D8383A","#63E398","#F3D266"),
                     labels = paste0(c('ASCVD','Framingham','Full model','Score2-Diabetes'),": ",
                                     c(sprintf("%0.2f",tROC_ascvd$AUC[10]),
                                       sprintf("%0.2f",tROC_FRM$AUC[10]),
                                       sprintf("%0.2f",tROC_RF$AUC[10]),
                                       sprintf("%0.2f",tROC_s2d$AUC[10]))))+
  geom_line(aes(x=c(0,1),y=c(0,1)),color = "grey",linetype ="dashed")+
  theme_classic()+
  theme(panel.grid = element_blank(),
        #legend.background = element_rect(linetype = 1, size = 0.2, colour = "black"),
        legend.position = c(0.7,0.125),legend.text=element_text(size =16), axis.title = element_text(size = 16),
        axis.line = element_line(linewidth = 1.0, color = "black"),
        axis.ticks = element_line(size = 1.0),
        axis.ticks.length = unit(0.2, "cm"),
        axis.text = element_text(size = 15))+
  scale_x_continuous(expand = c(0.005,0.005))+
  scale_y_continuous(expand = c(0.005,0.005))+
  labs(x = "1 - Specificity",
       y = "Sensitivity")
vs_10_1

vs_3_1+vs_5_1+vs_10_1+plot_layout(ncol=3,nrow=1)+plot_annotation(tag_levels="A")&theme(plot.tag = element_text(size = 14,face = 'bold'))
dev.off()### RF_with_vs_without_3_5_10

#########Roc of 15 years#############
AUC   <- as.data.frame(tROC_RF$AUC)
AUC_1 <- as.data.frame(tROC_RF_1$AUC)
AUC$time   <- c(1:15)
AUC_1$time <- c(1:15)
AUC <- AUC[2:15,]
AUC_1 <- AUC_1[2:15,]

ROC_15years <- ggplot()+
  geom_point(mapping = aes(x= time, y= tROC_RF$AUC[2:15]),
             data = AUC,
             size = 3,
             shape = 16,
             color = "red")+
  geom_line(data = AUC, 
            aes(x= time, y= tROC_RF$AUC[2:15]),
            color="red",
            size= 1.1)+
  # geom_smooth(data = AUC, 
  #             aes(x= time, y= tROC_RF$AUC[3:15]),
  #             method="gam",se= FALSE,color="red")+
  geom_point(mapping = aes(x= time, y= tROC_RF_1$AUC[2:15]),
             data = AUC_1,
             size = 3,
             shape = 15,
             color = "blue")+
  geom_line(data = AUC_1, 
            aes(x= time, y= tROC_RF_1$AUC[2:15]),
            color="blue",
            size= 1.1)+
  # geom_smooth(data = AUC_1, 
  #             aes(x= time, y= tROC_RF_1$AUC[3:15]),
  #             method="gam",se= FALSE,color="blue")+
  scale_x_continuous(limits = c(1, 15), breaks= c(2:15))+
  theme_classic()
ROC_15years

###############Birer Score##############
library(randomForestSRC)

bs_km <- get.brier.survival(RFmodel_oli, 
                            cens.model = "km")$brier.score
head(bs_km)
bs_rsf <- get.brier.survival(RFmodel_oli, 
                             cens.model = "rfsrc")$brier.score
head(bs_rsf)
bs_rsf$model <- "Full Model"

bs_km_1 <- get.brier.survival(RFmodel_out, 
                              cens.model = "km")$brier.score
head(bs_km_1)
bs_rsf_1 <- get.brier.survival(RFmodel_out, 
                               cens.model = "rfsrc")$brier.score
head(bs_rsf_1)
bs_rsf_1$model <- "Basic Model"

bs <- rbind(bs_rsf, bs_rsf_1)
birer_socre <- ggplot(data = bs)+
  geom_line(mapping = aes(x= time, y= brier.score,color = model),
             size = 0.7
            )+
  scale_color_manual(values = c("blue","red"),
                     labels = c("Basic model","Full    model"))+
  # geom_smooth(aes(x= time, y= brier.score,color = model),
  #             method="gam",se= TRUE,size=1.1)+
  theme_classic()
birer_socre
t.test(bs_km$brier.score, bs_km_1$brier.score, paired = TRUE, conf.level = 0.95)
dev.off()

########寻找最佳截点画生存曲线#########
install.packages("survminer")
library(survminer)
library(ggpubr)

# library(survival)
# library(survminer)
# data(myeloma)
# head(myeloma)
# res.cut1 <- surv_cutpoint(myeloma, time = "time", event = "event",
#                          variables = c("CCND1", "CRIM1", "DEPDC1") # 找这3个变量的最佳切点
# )
# 
# summary(res.cut1)
# plot(res.cut1, "DEPDC1", palette = "npg")

RSF_train   <- predict(RFmodel_oli, newdata = train_set)#训练集预测
RSF_train_1 <- predict(RFmodel_out, newdata = train_set_1)
train_set$pred    <- RSF_train$predicted  
train_set_1$pred  <- RSF_train_1$predicted
                                      
ggplot(train_set, aes(x=pred))+
  geom_histogram(fill="steelblue", color="white")+
  theme_classic()
res.cut <- surv_cutpoint(train_set, time="followuptime", event="incidentCVD",
                         variables=c("pred"))
test_cut <- tmp.test[,c("followuptime","incidentCVD","pred")]
test_cut$pred <- ifelse(test_cut$pred>=26.66554,"high","low")
fit <- survfit(Surv(followuptime, incidentCVD)~ pred,data=test_cut)
# res.cat <- surv_categorize(res.cut)
km_cutpoint <-ggsurvplot(fit,
                         data=test_cut, 
                         surv.median.line = "hv",
                         legend.title ="Risk Group",
                         legend.labs = c("High Risk","Low risk"),
                         pval = TRUE,
                         risk.table = "abs_pct",
                         conf.int = TRUE,
                         palette ="npg",
                         ncensor.plot = TRUE,
                         ggtheme = theme_classic2(),
                         xlim = c(0,15),
                         break.time.by =1
) #>cutpoint(26.66554) is High
km_cutpoint

ggplot(train_set_1, aes(x=pred))+
  geom_histogram(fill="steelblue", color="white")+
  theme_classic()
res.cut_1 <- surv_cutpoint(train_set_1, time="followuptime", event="incidentCVD",
                         variables=c("pred")
                        )
test_cut_1 <- tmp.test_1[,c("followuptime","incidentCVD","pred")]
test_cut_1$pred <- ifelse(test_cut_1$pred>=23.31984,"high","low")
fit_1 <- survfit(Surv(followuptime, incidentCVD)~ pred,data=test_cut_1)
summary(fit_1)
km_cutpoint_1 <-ggsurvplot(fit_1,
                          data=test_cut_1,
                          surv.median.line = "hv",
                          legend.title ="Risk Group",
                          legend.labs = c("High Risk","Low risk"),
                          pval = TRUE,
                          risk.table = "abs_pct",
                          conf.int = TRUE,
                          palette ="npg",
                          ncensor.plot = TRUE,
                          ggtheme = theme_classic2(),
                          xlim = c(0,15),
                          break.time.by =1
                          ) #>cutpoint(23.31984) is High
km_cutpoint_1

table(test_cut$incidentCVD,test_cut$pred)
table(test_cut_1$incidentCVD,test_cut_1$pred)

########NRI/IDI 净重分类指数NRI；综合判别改善指数IDI######
library(nricens)
library(gt)
uu_RF <- predict(RFmodel_oli, newdata = data_uu)#全人群中进行模型预测
uu_RF_1 <- predict(RFmodel_out, newdata = data_uu)
uu_with <- subset(data_uu, select = c(covtest,olicov,"followuptime","incidentCVD"))
uu_out  <- subset(data_uu, select = c(covtest,"followuptime","incidentCVD"))
uu_with$pred <- uu_RF$predicted  #全人群 With Olink
uu_out$pred  <- uu_RF_1$predicted#全人群 Without Olink

#建立全人群cox回归模型
coxfit.rsf.test <- coxph(Surv(followuptime, incidentCVD==1) ~ pred, data=uu_with, x=T)
#rsf with olink
coxfit.std.test <- coxph(Surv(followuptime, incidentCVD==1) ~ pred, data=uu_out,  x=T)
#std without olink

NRIvalue <- matrix(nrow = 16, ncol = 10, byrow = FALSE, 
            dimnames = list(1:16,c("NRI","NRI1","NRI3","NRI+","NRI+1","NRI+3","NRI-","NRI1-","NRI3-","P")))
NRIvalue <- as.data.frame(NRIvalue)
par(mfrow = c(3,5))
dev.off()
for (i in 2:16) {
set.seed(1234)
NRI <- nricens(mdl.std = coxfit.std.test,
          mdl.new = coxfit.rsf.test,
          t0 = 16,
          cut = c(0.3,0.7),
          niter = 1000,
          updown = "category")
NRIvalue[i-1,1:3] <- NRI$nri[1,1:3]
NRIvalue[i-1,4:6] <- NRI$nri[2,1:3]
NRIvalue[i-1,7:9] <- NRI$nri[3,1:3]
}#NRI绘图4*4
write.csv(NRIvalue, "0409NRIvalue.csv", row.names = FALSE)#save NRIvalue
par(mfrow = c(1,1))

# tROC_RF$AUC
# nricens(time = data_uu$followuptime,  
#         event = data_uu$incidentCVD,
#         p.std = uu_out$pred/100, 
#         p.new = uu_with$pred/100,
#         t0 = 1, 
#         cut = 34.14652/100, 
#         niter = 100)

library(survIDINRI)
z.new.test <- as.matrix(subset(uu_with, select = pred))
z.std.test <- as.matrix(subset(uu_out, select = pred))
IDIvalue <- matrix(nrow = 16, ncol = 8, byrow = FALSE, 
                   dimnames = list(1:16,c("IDI","IDILower","IDIUpper","IDIP","NRI","NRILower","NRIUpper","NRIp")))
IDIvalue <- as.data.frame(IDIvalue)

for (i in 2:16) {
  res.test <- IDI.INF(indata = data_uu[,c(28,29)],
                     covs0 = z.std.test,
                     covs1 = z.new.test,
                     t0= i,
                     npert = 10000,
                     seed1 = 1234
                     )
IDI <- IDI.INF.OUT(res.test)
IDI.INF.GRAPH(res.test,
              size = 1)
IDIvalue[i,1:4] <- IDI[1,]
IDIvalue[i,5:8] <- IDI[2,]
}
# M1 IDI;M2 NRI
write.csv(IDIvalue, "0409IDI_NRI_value.csv", row.names = FALSE)#save IDIvalue

##############校准曲线#################
install.packages("Hmisc")
install.packages("Matrix")
library(rms)
library(survival)
library(riskRegression)
library(magrittr)
library(dplyr)
library(ggplot2)

tmp.train <- train_set
tmp.train$pred <- RFmodel_oli$predicted
tmp.train_1 <- train_set_1
tmp.train_1$pred <- RFmodel_out$predicted
set.seed(1234)

cox_fit_train   <- coxph(Surv(followuptime, incidentCVD) ~ pred, data=tmp.train  , x=T, y=T)
cox_fit_train_1 <- coxph(Surv(followuptime, incidentCVD) ~ pred, data=tmp.train_1, x=T, y=T)

cox_fit_1 <- Score(list("fit1" = cox_fit_train,
                        "fit2" = cox_fit_train_1
                        ),
                   formula = Surv(followuptime, incidentCVD) ~1,
                   data = tmp.train,
                   plots = "calibration",
                   B= 100,
                   M= 50,
                   time= c(3),
                   seed = 1234
                   )
plotCalibration(cox_fit_1,
                rug=T,
                cens.method="local",
                xlab="Predicted Risk",
                #xlim = c(0.0,0.25),
                ylab="Observed  Risk",
                #ylim = c(0.0,0.8)
                #method = "quantile"
                )

data_plot <- plotCalibration(cox_fit_1, plot = T) 
plot_df <- bind_rows(data_plot$plotFrames) %>%
  mutate(fits = rep(c("fit1","fit2"),c(23,28)))
plot_df <- plot_df[c(1,4,10,14,16,18,24,28,32,36,39,43),]
ggplot(plot_df)+
  #geom_smooth(aes(x= Pred, y= Obs,color=fits),
  #            method="lm", se=F)+
  geom_line(aes(x= Pred, y= Obs,color=fits),size=1)+
  geom_point(aes(x= Pred, y= Obs,color=fits),size = 1.2) +
  scale_color_manual(values = c("red","blue"),name=NULL)+
  scale_x_continuous(limits = c(0,0.45),name = "Predicted Risk")+
  scale_y_continuous(limits = c(0,0.45),name = "Observerd Risk")+
  geom_abline(slope = 1,intercept = 0,lty=2,lwd=1)+
  geom_rug(aes(color=fits))+
  theme_bw()#3 year train

cox_fit_2 <- Score(list("fit1" = cox_fit_train,
                        "fit2" = cox_fit_train_1
                        ),
                   formula = Surv(followuptime, incidentCVD) ~1,
                   data = tmp.test,
                   plots = "calibration",
                   B= 100,
                   M= 50,
                   time= c(3),
                   seed = 1234
                   )
plotCalibration(cox_fit_2,
                rug=T,
                cens.method="local",
                xlab="Predicted Risk",
                #xlim = c(0.0,0.25),
                ylab="Observed  Risk",
                #ylim = c(0.0,0.8)
                #method = "quantile"
                )

data_plot_1 <- plotCalibration(cox_fit_2, plot = T) 
plot_df_1 <- bind_rows(data_plot_1$plotFrames) %>%
  mutate(fits = rep(c("fit1","fit2"),c(24,30)))
plot_df_1 <- plot_df_1[c(1,4,11,16,18,25,30,36,37,47),]
ggplot(plot_df_1)+
  #geom_smooth(aes(x= Pred, y= Obs,color=fits),
  #            method="lm", se=F)+
  geom_line(aes(x= Pred, y= Obs,color=fits),size=1)+
  geom_point(aes(x= Pred, y= Obs,color=fits),size = 1.2) +
  scale_color_manual(values = c("red","blue"),name=NULL)+
  scale_x_continuous(limits = c(0,0.3),name = "Predicted Risk")+
  scale_y_continuous(limits = c(0,0.3),name = "Observerd Risk")+
  geom_abline(slope = 1,intercept = 0,lty=2,lwd=1)+
  geom_rug(aes(color=fits))+
  theme_bw()#3 year test

cox_fit_1 <- Score(list("fit1" = cox_fit_train,
                        "fit2" = cox_fit_train_1
                        ),
                   formula = Surv(followuptime, incidentCVD) ~1,
                   data = tmp.train,
                   plots = "calibration",
                   B= 100,
                   M= 50,
                   time= c(5),
                   seed = 1234
                   )
plotCalibration(cox_fit_1,
                rug=T,
                cens.method="local",
                xlab="Predicted Risk",
                #xlim = c(0.0,0.25),
                ylab="Observed  Risk",
                #ylim = c(0.0,0.8)
                #method = "quantile"
)

data_plot <- plotCalibration(cox_fit_1, plot = T) 
plot_df <- bind_rows(data_plot$plotFrames) %>%
  mutate(fits = rep(c("fit1","fit2"),c(32,35)))
plot_df <- plot_df[c(1,5,12,17,23,33,38,45,50,56),]
ggplot(plot_df)+
  #geom_smooth(aes(x= Pred, y= Obs,color=fits),
  #            method="lm", se=F)+
  geom_line(aes(x= Pred, y= Obs,color=fits),size=1)+
  geom_point(aes(x= Pred, y= Obs,color=fits),size = 1.2) +
  scale_color_manual(values = c("red","blue"),name=NULL)+
  scale_x_continuous(limits = c(0,0.5),name = "Predicted Risk")+
  scale_y_continuous(limits = c(0,0.5),name = "Observerd Risk")+
  geom_abline(slope = 1,intercept = 0,lty=2,lwd=1)+
  geom_rug(aes(color=fits))+
  theme_bw()#5 year train

cox_fit_2 <- Score(list("fit1" = cox_fit_train,
                        "fit2" = cox_fit_train_1
                        ),
                   formula = Surv(followuptime, incidentCVD) ~1,
                   data = tmp.test,
                   plots = "calibration",
                   B= 100,
                   M= 50,
                   time= c(5),
                   seed = 1234
                   )
plotCalibration(cox_fit_2,
                rug=T,
                cens.method="local",
                xlab="Predicted Risk",
                #xlim = c(0.0,0.25),
                ylab="Observed  Risk",
                #ylim = c(0.0,0.8)
                #method = "quantile"
)

data_plot_1 <- plotCalibration(cox_fit_2, plot = T) 
plot_df_1 <- bind_rows(data_plot_1$plotFrames) %>%
  mutate(fits = rep(c("fit1","fit2"),c(30,36)))
plot_df_1 <- filter(plot_df_1,Pred<=0.18)
plot_df_1 <- plot_df_1[c(1,4,6,10,14,15,18,21,25,29),]
ggplot(plot_df_1)+
  #geom_smooth(aes(x= Pred, y= Obs,color=fits),
  #            method="lm", se=F)+
  geom_line(aes(x= Pred, y= Obs,color=fits),size=1)+
  geom_point(aes(x= Pred, y= Obs,color=fits),size = 1.2) +
  scale_color_manual(values = c("red","blue"),name=NULL)+
  scale_x_continuous(limits = c(0,0.2),name = "Predicted Risk")+
  scale_y_continuous(limits = c(0,0.2),name = "Observerd Risk")+
  geom_abline(slope = 1,intercept = 0,lty=2,lwd=1)+
  geom_rug(aes(color=fits))+
  theme_bw()#5 year test

cox_fit_1 <- Score(list("fit1" = cox_fit_train,
                        "fit2" = cox_fit_train_1
                        ),
                   formula = Surv(followuptime, incidentCVD) ~1,
                   data = tmp.train,
                   plots = "calibration",
                   B= 100,
                   M= 50,
                   time= c(10),
                   seed = 1234
                   )
plotCalibration(cox_fit_1,
                rug=T,
                cens.method="local",
                xlab="Predicted Risk",
                #xlim = c(0.0,0.25),
                ylab="Observed  Risk",
                #ylim = c(0.0,0.8)
                #method = "quantile"
)

data_plot <- plotCalibration(cox_fit_1, plot = T) 
plot_df <- bind_rows(data_plot$plotFrames) %>%
  mutate(fits = rep(c("fit1","fit2"),c(52,53)))
plot_df <- plot_df[c(1,4,9,13,26,43,54,57,62,65,83,98),]
ggplot(plot_df)+
  #geom_smooth(aes(x= Pred, y= Obs,color=fits),
  #            method="gam", se=F)+
  geom_line(aes(x= Pred, y= Obs,color=fits),size=1)+
  geom_point(aes(x= Pred, y= Obs,color=fits),size = 1.2) +
  scale_color_manual(values = c("red","blue"),name=NULL)+
  scale_x_continuous(limits = c(0,0.8),name = "Predicted Risk")+
  scale_y_continuous(limits = c(0,0.8),name = "Observerd Risk")+
  geom_abline(slope = 1,intercept = 0,lty=2,lwd=1)+
  geom_rug(aes(color=fits))+
  theme_bw()#10 year train

cox_fit_2 <- Score(list("fit1" = cox_fit_train,
                        "fit2" = cox_fit_train_1
                        ),
                   formula = Surv(followuptime, incidentCVD) ~1,
                   data = tmp.test,
                   plots = "calibration",
                   B= 100,
                   M= 50,
                   time= c(10),
                   seed = 1234
                   )
plotCalibration(cox_fit_2,
                rug=T,
                cens.method="local",
                xlab="Predicted Risk",
                #xlim = c(0.0,0.25),
                ylab="Observed  Risk",
                #ylim = c(0.0,0.8)
                method = "quantile"
)

data_plot_1 <- plotCalibration(cox_fit_2, plot = T) 
plot_df_1 <- bind_rows(data_plot_1$plotFrames) %>%
  mutate(fits = rep(c("fit1","fit2"),c(49,50)))
plot_df_1 <- plot_df_1[c(1,10,18,20,31,51,56,58,68,83),]
ggplot(plot_df_1)+
  #geom_smooth(aes(x= Pred, y= Obs,color=fits),
  #            method="lm", se=F)+
  geom_line(aes(x= Pred, y= Obs,color=fits),size=1)+
  geom_point(aes(x= Pred, y= Obs,color=fits),size = 1.2) +
  scale_color_manual(values = c("red","blue"),name=NULL)+
  scale_x_continuous(limits = c(0,0.5),name = "Predicted Risk")+
  scale_y_continuous(limits = c(0,0.5),name = "Observerd Risk")+
  geom_abline(slope = 1,intercept = 0,lty=2,lwd=1)+
  geom_rug(aes(color=fits))+
  theme_bw()#10 year test

# coxfit1 <- cph(Surv(followuptime, incidentCVD) ~ pred,
#                data = tmp.train, x=T,y=T,surv = T,
#                time.inc = 10
#                )
# cal1 <- calibrate(coxfit1, cmethod="KM", method="boot",u=10, m=200,B=500) 
# plot(cal1,
#      #lwd = 1.5, # 误差线粗细
#      #lty = 1, # 误差线类型，可选0-6
#      #errbar.col = c("#2166AC"), # 误差线颜色
#      xlim = c(0.0,1),ylim= c(0.0,1),
#      xlab = "Nomogram-prediced OS (%)",ylab = "Observed OS (%)",
#      cex.lab=1.2, cex.axis=1, cex.main=1.2, cex.sub=0.6) # 字体大小
# lines(cal1[,c('mean.predicted',"KM")], 
#       type = 'b', # 连线的类型，可以是"p","b","o"
#       lwd = 3, # 连线的粗细
#       pch = 16, # 点的形状，可以是0-20
#       col = "tomato") # 连线的颜色
# box(lwd = 2) # 边框粗细
# abline(0,1,lty = 3, # 对角线为虚线
#        lwd = 2, # 对角线的粗细
#        col = "grey70" # 对角线的颜色
# ) 
# 
# coxfit2 <- cph(Surv(followuptime, incidentCVD) ~ pred,
#                data = tmp.train_1, x=T,y=T,surv = T,
#                time.inc = 5
# )
# 
# cal2 <- calibrate(coxfit2, cmethod="KM", method="boot",u=5,m=200,B=500) 
# plot(cal1,lwd = 2,lty = 0,errbar.col = c("#2166AC"),
#      xlim = c(0.5,1),ylim= c(0.5,1),
#      xlab = "Nomogram-prediced OS (%)",ylab = "Observed OS (%)",
#      col = c("#2166AC"),
#      cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6)
# lines(cal1[,c('mean.predicted',"KM")],
#       type = 'b', lwd = 3, col = c("#2166AC"), pch = 16)
# plot(cal2,lwd = 2,lty = 0,errbar.col = c("#B2182B"),
#      xlim = c(0,1),ylim= c(0,1),col = c("#B2182B"),add = T)
# lines(cal2[,c('mean.predicted',"KM")],
#       type = 'b', lwd = 3, col = c("#B2182B"), pch = 16)
# abline(0,1, lwd = 2, lty = 3, col = c("#224444"))
# legend("bottomright", #图例的位置
#        legend = c("5-year","10-year"), #图例文字
#        col =c("#2166AC","#B2182B"), #图例线的颜色，与文字对应
#        lwd = 2,#图例中线的粗细
#        cex = 1.2,#图例字体大小
#        bty = "n")#不显示图例边框
# 
# cox_fit_2 <- Score(list("fit" = cox_fit_train),
#                    formula = Surv(followuptime, incidentCVD) ~1,
#                    data = tmp.test,
#                    plots = "calibration",
#                    B= 100,
#                    M= 50,
#                    time= c(10)
#                    )
# plotCalibration(cox_fit_2,
#                 rug=T,
#                 cens.method="local",
#                 xlab="Predicted Risk",
#                 ylab="Observed  Risk",
#                 method = "quantile")
# 
# train_prob <- c(predictRisk(RFmodel_oli, newdata = tmp.test, times = 3))
# cuts <- unique(quantile(c(0,1,train_prob),seq(0,1,length=20),na.rm=TRUE))
# cuts
# 
# km_surv <- groupkm(train_prob,
#                    Srv = Surv(tmp.test$followuptime,tmp.test$incidentCVD),
#                    u = 3,
#                    cuts = cuts)
# km_surv[,4] <- 1-km_surv[,4]
# km_surv
# #等渗回归
# plot(km_surv[,1], km_surv[,4],
#      xlim=c(0,1),ylim=c(0,1),
#      xlab = 'Predicted probability of 3-year',
#      ylab = 'Observed probability of 3-year',
#      cex.axis = 1.4, # 调整坐标轴标签的大小
#      cex.lab = 1.5 # 调整坐标轴标题的大小
# )
# lines(km_surv[,1], km_surv[,4])
# 
# # 计算误差线范围
# errl <- ifelse(km_surv[,"KM"] == 0, 0,  
#                km_surv[,"KM"] * exp(1.959964 * (-km_surv[,"std.err"])))
# errh <- ifelse(km_surv[,"KM"] == 0, 0, 
#                pmin(1, km_surv[,"KM"] * exp(1.959964 * km_surv[,"std.err"])))
# # 添加误差线
# errbar(x = km_surv[,"x"],
#        y = km_surv[,"KM"],
#        yminus = errl,yplus = errh,
#        add = T,
#        pch=16,cex=1,
#        asp=1,xaxs='i',yaxs='i'
# )
# # 添加对角线
# abline(a = 0,b = 1,col='grey')
# #plot(smooth(ROC),col="red",print.auc=T,legacy.axes=T)

##########决策曲线#################
RFmodel

if (!require(rmda)) {
  install.packages("rmda")
}
## Warning: 程辑包'rmda'是用R版本4.1.3 来建造的

if (!require(ggDCA)) {
  devtools::install_github("yikeshu0611/ggDCA")
}
if (!require(rms)) {
  install.packages("rms")
}
if (!require(caret)) {
  install.packages("caret")
}
install.packages("ggsci")
library(rmda)
library(ggDCA)
library(ggplot2)
library(rms)
library(caret)
library(riskRegression)
library(tidyverse)
library(ggsci)
library(utf8)
library(tidyr)
library(tidyverse)
tmp.train$pred2 <- tmp.train_1$pred
rf_dca.test <- stdca(data=tmp.train,
                     outcome = "incidentCVD",
                     ttoutcome = "followuptime",
                     timepoint = 3,
                     predictors = c("pred","pred2"),
                     probability = c(F,F),
                     smooth = T,
                     graph = F
)

rf_dca.test$net.benefit %>%
  pivot_longer(cols = c(all, none, contains("sm")),
               names_to ="models",
               values_to = "net_benefit") %>%
  dplyr::mutate(models = factor(models,
                                levels =c("all","none","pred_sm","pred2_sm")))%>%
  ggplot(aes(threshold, net_benefit))+
  geom_line(aes(color=models),size=1.2)+
  scale_color_jama(name="Models Types",
                   labels=c("All","None","Full Model","Basic Model"))+
  scale_x_continuous(labels= scales::label_percent(accuracy=1),
                     name="Threshold Probility",
                     limits = c(0,0.25))+
  scale_y_continuous(limits = c(-0.01,0.03),name="Net Benefit")+
  theme_bw(base_size = 14)+
  theme(legend.background = element_blank(),
        legend.position = c(0.85,0.75)
  )#3year train

rf_dca.test2 <- stdca(data=tmp.train,
                     outcome = "incidentCVD",
                     ttoutcome = "followuptime",
                     timepoint = 5,
                     predictors = c("pred","pred2"),
                     probability = c(F,F),
                     smooth = T,
                     graph = F
)

rf_dca.test2$net.benefit %>%
  pivot_longer(cols = c(all, none, contains("sm")),
               names_to ="models",
               values_to = "net_benefit") %>%
  dplyr::mutate(models = factor(models,
                                levels =c("all","none","pred_sm","pred2_sm")))%>%
  ggplot(aes(threshold, net_benefit))+
  geom_line(aes(color=models),size=1.2)+
  scale_color_jama(name="Models Types",
                   labels=c("All","None","Full Model","Basic Model"))+
  scale_x_continuous(labels= scales::label_percent(accuracy=1),
                     name="Threshold Probility",
                     limits = c(0,0.2))+
  scale_y_continuous(limits = c(-0.01,0.075),name="Net Benefit")+
  theme_bw(base_size = 14)+
  theme(legend.background = element_blank(),
        legend.position = c(0.85,0.75)
  )#5year train

rf_dca.test3 <- stdca(data=tmp.train,
                     outcome = "incidentCVD",
                     ttoutcome = "followuptime",
                     timepoint = 10,
                     predictors = c("pred","pred2"),
                     probability = c(F,F),
                     smooth = T,
                     graph = F
)

rf_dca.test3$net.benefit %>%
  pivot_longer(cols = c(all, none, contains("sm")),
               names_to ="models",
               values_to = "net_benefit") %>%
  dplyr::mutate(models = factor(models,
                                levels =c("all","none","pred_sm","pred2_sm")))%>%
  ggplot(aes(threshold, net_benefit))+
  geom_line(aes(color=models),size=1.2)+
  scale_color_jama(name="Models Types",
                   labels=c("All","None","Full Model","Basic Model"))+
  scale_x_continuous(labels= scales::label_percent(accuracy=1),
                     name="Threshold Probility",
                     limits = c(0,0.5))+
  scale_y_continuous(limits = c(-0.01,0.2),name="Net Benefit")+
  theme_bw(base_size = 14)+
  theme(legend.background = element_blank(),
        legend.position = c(0.85,0.75)
  )#10year train

tmp.test$pred2 <- tmp.test_1$pred
cox.fit.dca1 <- coxph(Surv(followuptime, incidentCVD) ~ pred, data=tmp.test)
tmp.test$pred3 <- c((1-(summary(survfit(cox_fit_train, newdata= tmp.test),time=3)$surv)))
cox.fit.dca2 <- coxph(Surv(followuptime, incidentCVD) ~ pred2, data=tmp.test)
tmp.test$pred4 <- c((1-(summary(survfit(cox_fit_train_1, newdata= tmp.test),time=3)$surv)))
rf_dca.test41 <- stdca(data=tmp.test,
                     outcome = "incidentCVD",
                     ttoutcome = "followuptime",
                     timepoint = 3,
                     predictors = c("pred"),
                     probability = c(F),
                     smooth = T,
                     graph = T
)

rf_dca.test4$net.benefit %>%
  pivot_longer(cols = c(all, none, contains("sm")),
               names_to ="models",
               values_to = "net_benefit") %>%
  dplyr::mutate(models = factor(models,
                                levels =c("all","none","pred3_sm","pred4_sm")))%>%
  ggplot(aes(threshold, net_benefit))+
  geom_line(aes(color=models),size=1.2)+
  scale_color_jama(name="Models Types",
                   labels=c("All","None","Full Model","Basic Model"))+
  scale_x_continuous(labels= scales::label_percent(accuracy=1),
                     name="Threshold Probility",
                     limits = c(0,0.2))+
  scale_y_continuous(limits = c(-0.001,0.01),name="Net Benefit")+
  theme_bw(base_size = 14)+
  theme(legend.background = element_blank(),
        legend.position = c(0.85,0.75)
  )#3year train

library(rmda)
plot_clinical_impact(rf_dca.test41,population.size= 1000,
                     cost.benefit.axis = T,
                     n.cost.benefits= 8,
                     col =c('red','blue'),
                     confidence.intervals= T,
                     ylim=c(0,1000),
                     legend.position="topright")

rf_dca.test5 <- stdca(data=tmp.test,
                      outcome = "incidentCVD",
                      ttoutcome = "followuptime",
                      timepoint = 5,
                      predictors = c("pred","pred2"),
                      probability = c(F,F),
                      smooth = T,
                      graph = T
)

rf_dca.test5$net.benefit %>%
  pivot_longer(cols = c(all, none, contains("sm")),
               names_to ="models",
               values_to = "net_benefit") %>%
  dplyr::mutate(models = factor(models,
                                levels =c("all","none","pred_sm","pred2_sm")))%>%
  ggplot(aes(threshold, net_benefit))+
  geom_line(aes(color=models),size=1.2)+
  scale_color_jama(name="Models Types",
                   labels=c("All","None","Full Model","Basic Model"))+
  scale_x_continuous(labels= scales::label_percent(accuracy=1),
                     name="Threshold Probility",
                     limits = c(0,0.2))+
  scale_y_continuous(limits = c(-0.01,0.08),name="Net Benefit")+
  theme_bw(base_size = 14)+
  theme(legend.background = element_blank(),
        legend.position = c(0.85,0.75)
  )#5year train

rf_dca.test6 <- stdca(data=tmp.test,
                      outcome = "incidentCVD",
                      ttoutcome = "followuptime",
                      timepoint = 10,
                      predictors = c("pred","pred2"),
                      probability = c(F,F),
                      smooth = T,
                      graph = T
)

rf_dca.test6$net.benefit %>%
  pivot_longer(cols = c(all, none, contains("sm")),
               names_to ="models",
               values_to = "net_benefit") %>%
  dplyr::mutate(models = factor(models,
                                levels =c("all","none","pred_sm","pred2_sm")))%>%
  ggplot(aes(threshold, net_benefit))+
  geom_line(aes(color=models),size=1.2)+
  scale_color_jama(name="Models Types",
                   labels=c("All","None","Full Model","Basic Model"))+
  scale_x_continuous(labels= scales::label_percent(accuracy=1),
                     name="Threshold Probility",
                     limits = c(0,0.3))+
  scale_y_continuous(limits = c(-0.01,0.18),name="Net Benefit")+
  theme_bw(base_size = 14)+
  theme(legend.background = element_blank(),
        legend.position = c(0.85,0.75)
  )#10year train

# train_set$incidentCVD <- as.numeric(train_set$incidentCVD)
# df <- ggDCA::dca(cox_fit_train,
#                  times=3
# )
# ggplot(df)+
#   scale_x_continuous(labels= scales::label_percent(accuracy=1),
#                      name="Threshold Probility",
#                      limits = c(0,0.15),)+
#   scale_y_continuous(limits = c(-0.025,0.025),name="Net Benefit")+
#   theme_bw(base_size = 14)+
#   theme(legend.background = element_blank(),
#         legend.position = c(0.85,0.75))+
#   geom_smooth(method = "loess", se = FALSE)
# 
# df_rf .<- dca(data= cox_fit_train
#              # outcome = "incidentCVD",
#              # predictors ="pred",
#              # probability =T,
#              # graph =T
# )
# ggplot(df_rf)
# install.packages("dcurves")
# library(dcurves)
# library(survival)
# library(ggplot2)
# library(dplyr)
# train_set$pred1 <- c(1-summary(survfit(cox_fit_train,newdata = train_set),times=10)$surv)
# dca(Surv(followuptime,incidentCVD)~pred1,
#     data=train_set,time =10) %>% plot()

################相关性################
library(corrplot)
library(ggplot2)
library(ggcorrplot)
library(vcd)
library(psych)
library(ggrepel)
library(nortest)
library(pheatmap)
library(Hmisc)
library(grDevices)

data_spearman <- data_uu[,c(covtest,olicov)]
head(data_spearman)
data_spearman$age <- scale(data_spearman$age)
data_spearman$bmi_1 <- scale(data_spearman$bmi_1)
data_spearman$hba1c_1 <- scale(data_spearman$hba1c_1)

col3 <- colorRampPalette(c("blue", "white","red"))
spearman <- cor(data_spearman, method = "spearman")
res <- rcorr(as.matrix(data_spearman))
corrplot(res$r, method="color", #square方形，ellipse, 椭圆形，number数值，shade阴影，color颜色，pie饼图
         title = "pearson",   #指定标题
         type="lower",  #full完全(默认)，lower下三角，upper上三角
         col=col3(100), #指定图形展示的颜色，默认以均匀的颜色展示。支持grDevices包中的调色板，也支持RColorBrewer包中调色板。
         outline = T,  #是否添加圆形、方形或椭圆形的外边框，默认为FALSE。
         diag = FALSE,  #是否展示对角线上的结果，默认为TRUE
         mar = c(0,0,0,0), #设置图形的四边间距。数字分别对应(bottom, left, top, right)。
         bg="white", #指定背景颜色
         add = FALSE, #表示是否添加到已经存在的plot中。默认FALSE生成新plot。
         is.corr = TRUE, #是否为相关系数绘图，默认为TRUE,FALSE则可将其它数字矩阵进行可视化。
         addgrid.col = "darkgray", #设置网格线颜色，当指定method参数为color或shade时， 默认的网格线颜色为白色，其它method则默认为灰色，也可以自定义颜色。
         addCoef.col = NULL, #设置相关系数值的颜色，只有当method不是number时才有效
         addCoefasPercent = TRUE, #是否将相关系数转化为百分比形式，以节省空间，默认为FALSE。
         order = "original", #指定相关系数排序的方法, 可以是original原始顺序，AOE特征向量角序，FPC第一主成分顺序，hclust层次聚类顺序，alphabet字母顺序。
         hclust.method = "complete", # 指定hclust中细分的方法，只有当指定order参数为hclust时有效。有7种可选：complete,ward,single,average,mcquitty,median,centroid。
         addrect = NULL, #是否添加矩形框，只有当指定order参数为hclust时有效， 默认不添加， 用整数指定即可添加。
         rect.col = "black", #指定矩形框的颜色。
         rect.lwd = 2, #指定矩形框的线宽。
         tl.pos = "ld",  #指定文本标签(变量名称)相对绘图区域的位置，为"lt"(左侧和顶部),"ld"(左侧和对角线),"td"(顶部和对角线),"d"(对角线),"n"(无);当type="full"时默认"lt"。当type="lower"时默认"ld"。当type="upper"时默认"td"。
         tl.cex = 0.8,  #设置文本标签的大小
         tl.col = "black", #设置文本标签的颜色。
         cl.pos = NULL, #设置图例位置，为"r"(右边),"b"(底部),"n"(无)之一。当type="full"/"upper"时，默认"r"; 当type="lower"时，默认"b"。
         # addshade = c("negative", "positive", "all"), # 表示给增加阴影，只有当method="shade"时有效。#为"negative"(对负相关系数增加阴影135度)；"positive"(对正相关系数增加阴影45度)；"all"(对所有相关系数增加阴影)。
         # shade.lwd = 0,  #指定阴影线宽。
         # shade.col = "white",#指定阴影线的颜色。   
         p.mat= res$P,sig.level= 0.01/22,insig= "label_sig",pch.cex = 1.5#只有指定矩阵的P值，sig.level，pch等参数才有效。只有当insig = "pch"时，pch.col和pch.cex参数才有效。
         )
corrplot(res$r,type="lower",add=TRUE,method="number", 
         tl.pos = "n",cl.pos = "n",diag=FALSE,
         col=col3(100), 
         number.cex = 0.7)

dev.off()

install.packages("Hmisc")
library(Hmisc)
library(lava)
library(stats)
ir <- isoreg(tmp.test$pred)
plot(ir, plot.type = "row")
ir$yf
