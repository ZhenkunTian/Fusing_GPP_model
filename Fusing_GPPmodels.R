# ===================
# Author: Zhenkun Tian
# Purpose: Calculating GPP using five LUE models: VPM, EC, GLO_PEM, CHJ, and C-Fix
#          Fusing these five individual models using BMA, SVM and RF
# Created: 2022-4-1
# Modified: 2023-1-31
# ===================

setwd("D:\\GPPmodel")


#n: sample size
#return the index set of n-fold validation
CV<-function(n,Z=5,seed=888){
  z<-rep(1:Z,each=ceiling(n/Z))[1:n]
  set.seed(seed)
  z=sample(z,n)
  mm=list()
  for(i in 1:Z){
    mm[[i]]<-(1:n)[z==i]
  }
  return (mm)
}

##C-Fix temperature stress function
f_t_CFix<-function(Ta){
  C1=21.77
  delt_HaP=52750
  Rg=8.31
  delt_S=704.98
  delt_HdP=211000
  Ta<-Ta+273.15
  X<-exp(C1-delt_HaP/(Rg*Ta))
  Y<-1+exp((delt_S*Ta-delt_HdP)/(Rg*Ta))
  return (X/Y)
}
x<-c(-20:50)
y<-f_t_CFix(x)
plot(x,y)

#C-Fix CO2 fertilization effect
CO2fert<-function(CO2){
  O2=20.9
  tao=2550
  CO2ref=281
  Km=948
  K0=30
  X<-(CO2-O2/(tao*2))/(CO2ref-O2/(tao*2))
  Y<-(Km*(1+O2/K0)+CO2ref)/(Km*(1+O2/K0)+CO2)
  return(X*Y)
}
x<-c(281:500)
y<-CO2fert(x)
plot(x,y)

#temperature stress fuction
f_t <- function(Ta, T_opt,T_min,T_max) {
  f=vector(length=length(Ta))
  for(i in 1:length(Ta)){
    if(Ta[i]<T_min | Ta[i]>T_max){
      f[i]<-0
    }
    else{
      f[i]<-(Ta[i]-T_min)*(Ta[i]-T_max)/((Ta[i]-T_min)*(Ta[i]-T_max)-(Ta[i]-T_opt)*(Ta[i]-T_opt))
    }
  }
  return(f)
}
x<-c(-20:50)
y<-f_t(x,20,0,40)
plot(x,y)

#vpm water stress function
f_w_vpm <- function(LSWI,LSWI_max){
  f=vector(length=length(LSWI))
  for(i in 1:length(LSWI)){
    if(LSWI[i]>LSWI_max){
      f[i]<-1
    }
    else{
      f[i] <- (1+LSWI[i])/(1+LSWI_max)
    }
  }
  P_scalar <- 1
  return(f*P_scalar)
}
x<-seq(-1,1,0.1)
y<-f_w_vpm(x,0.8)
plot(x,y)

#vpm P scalar function
#before leaf expansion: Ps=(1+lswi)/2; after leaf expansion: Ps=1
#DOY:day of year
f_vpm_Ps<-function(PTF,DOY,lswi){
  if(PTF=="EBF"|PTF=="ENF") ps=rep(1,times=length(lswi))
  else {
    ps=vector(length=length(DOY))
    for(i in 1:length(DOY)){
      if(DOY[i]>=105){
        ps[i]<-1
      }
      else{
        ps[i] <- (1+lswi[i])/2
      }
    }
  }
  return(ps)
}

#EC-LUE moisture stress function
f_w_EF<-function(LE,H){
  EF<-abs(LE)/(abs(LE)+abs(H))
  return(EF)
}


#CASA water stress function
f_w_CASA<-function(EET,PET){
  f=vector(length=length(EET))
  for(i in 1:length(EET)){
    if(is.na(EET[i])|is.na(PET[i])){
      f[i]=NA
    }
    else if(EET[i]>=PET[i]|EET[i]<0){
      f[i]=1
    }
    else if(EET[i]<PET[i]){
      f[i]<-EET[i]/PET[i]*0.5+0.5
    }
  }
  return(f)
}


#Priestley-Taylor model PET(W m-2)
PT_PET<-function(Ta,Rn,G=0){
  #Psychometric constant gamma=0.667
  gamma=0.667
  #Slope of vapour pressure curve (delta) for different temperatures (T)
  delta_Numerator<-4098*(0.6108*exp(17.27*Ta/(Ta+237.3)))
  delta_Denominator<-(Ta+237.3)*(Ta+237.3)
  delta<-delta_Numerator/delta_Denominator
  PET<-1.26*(delta/(delta+gamma))*(Rn-G)
  PET<-ifelse(PET<0,0,PET)
  return(PET)
}

LUE_parameters<-read.csv("sites_LUE_parameters.csv",header = T)
Temperature_parameters<-read.csv("IGBP_LUE_parameters.csv",header = T)

sitenames<-read.csv("sites-selected-56.csv",header=T)
df_metrics<-as.data.frame(matrix(nrow=0,ncol=8))
names(df_metrics)<-c("site","PFT","Lat","Lon","model","R2","RMSE","RPE")

#compute daily gpp for each site
for(i in 1:length(sitenames$SITE_ID)){
  flux_path<-paste("flux2015-sites-process-shared\\",sitenames$SITE_ID[i],sep="")
  flux_name<-paste(flux_path,"\\",sitenames$SITE_ID[i],"_flux2015_modis_daily_fullGPP.csv",sep="")
  print(flux_name)
  df <- read.csv(flux_name, header = T)
  df$TA_F[df$TA_F==-9999]<-NA
  df$SW_IN_F[df$SW_IN_F==-9999]<-NA
  df$VPD_F[df$VPD_F==-9999]<-NA
  df$PA_F[df$PA_F==-9999]<-NA
  df$NETRAD[df$NETRAD==-9999]<-NA
  df$CO2_F_MDS[df$CO2_F_MDS==-9999]<-NA
  df$LE_F_MDS[df$LE_F_MDS==-9999]<-NA
  df$LE_CORR[df$LE_CORR==-9999]<-NA
  df$H_F_MDS[df$H_F_MDS==-9999]<-NA
  df$H_CORR[df$H_CORR==-9999]<-NA
  df$G_F_MDS[df$G_F_MDS==-9999]<-NA
  df$NEE_VUT_REF[df$NEE_VUT_REF==-9999]<-NA
  df$GPP_DT_VUT_REF[df$GPP_DT_VUT_REF==-9999]<-NA
  df$GPP_NT_VUT_REF[df$GPP_NT_VUT_REF==-9999]<-NA
  df$GPP_DTNT_VUT_REF<-(df$GPP_DT_VUT_REF+df$GPP_NT_VUT_REF)/2
  
  df$PAR<-df$SW_IN_F*0.0864*0.45  # PAR
  #modis GPP unit to fluxtnet2015 gpp unit
  df$GPP<-df$GPP_modis_8days*1000/8  
  df2<-data.frame(calendar_date=df$calendar_date,gpp_obs=df$GPP_DTNT_VUT_REF,gpp_mod=df$GPP)
  
  # #retriev LUEmax(gC/m2/d/MJ)
  LUE_max<-LUE_parameters[LUE_parameters$SITE_ID==sitenames$SITE_ID[i],c("LUE")]
  ENF_Tmin<-Temperature_parameters[Temperature_parameters$IGBP==sitenames$IGBP[i],6]
  ENF_Tmax<-Temperature_parameters[Temperature_parameters$IGBP==sitenames$IGBP[i],7]
  ENF_Topt<-Temperature_parameters[Temperature_parameters$IGBP==sitenames$IGBP[i],8]
  
  #vpm
  lswi_max<-max(df$LSWI)
  #lswi_max
  ft<-f_t(Ta=df$TA_F,T_opt =ENF_Topt,T_min = ENF_Tmin,
          T_max =ENF_Tmax )
  fw_vpm<-f_w_vpm(LSWI = df$LSWI,LSWI_max =lswi_max )
  fw_vpm_Ps<-f_vpm_Ps(sitenames$IGBP[i],df$DOY,lswi=df$LSWI)
  LUE_vpm<-LUE_max*ft*fw_vpm*fw_vpm_Ps
  gpp_vpm<-LUE_vpm*df$PAR*df$Fpar
  df2$VPM<-gpp_vpm
  
  #GLO_PEM
  ft<-f_t(Ta=df$TA_F,T_opt =ENF_Topt,T_min = ENF_Tmin,
          T_max =ENF_Tmax )
  fw_GLO_PEM<-f_w_CASA(EET = df$LE_F_MDS,PET =df$NETRAD )
  LUE_GLO_PEM<-LUE_max*ft*fw_GLO_PEM
  gpp_GLO_PEM<-LUE_GLO_PEM*df$PAR*df$Fpar
  df2$GLO_PEM<-gpp_GLO_PEM
  
  #GPP_EC
  ft_EM<-f_t(Ta=df$TA_F,T_opt =ENF_Topt,T_min = ENF_Tmin,
             T_max =ENF_Tmax )
  fw_EM<-f_w_EF(df$LE_F_MDS,df$H_F_MDS)
  gpp_ec<-2.14*pmin(ft_EM,fw_EM)*df$PAR*df$Fpar
  df2$EC<-gpp_ec
  
  
  #GPP_max
  GPP_max<-max(df$GPP_DT_VUT_REF,na.rm = TRUE)
  ft<-f_t(Ta=df$TA_F,T_opt =ENF_Topt,T_min = ENF_Tmin,
          T_max =ENF_Tmax )
  lswi_max<-max(df$LSWI)
  fw_vpm<-f_w_vpm(LSWI = df$LSWI,LSWI_max =lswi_max )
  gpp_chj<-GPP_max*ft*fw_vpm
  df2$CHJ<-gpp_chj
  
  #C-Fix
  ft_C_Fix<-f_t_CFix(Ta=df$TA_F)
  fw_CO2<-CO2fert(df$CO2_F_MDS)
  gpp_C_Fix<-LUE_max*ft_C_Fix*fw_CO2*df$PAR*df$Fpar
  df2$C_Fix<-gpp_C_Fix
  
  #remove NA
  sum(is.na(df2$VPM))
  sum(is.na(df2$GLO_PEM))
  sum(is.na(df2$EC))
  sum(is.na(df2$CHJ))
  df2<-df2[(!is.na(df2$VPM))&(!is.na(df2$GLO_PEM))&(!is.na(df2$EC))&(!is.na(df2$CHJ))&(!is.na(df2$C_Fix)),]
  
  #vpm metrics
  SSE<-sum((df2$gpp_obs-df2$VPM)^2,na.rm=TRUE)
  SST<-sum((df2$gpp_obs-mean(df2$gpp_obs,na.rm = TRUE))^2,na.rm=TRUE)
  R2_vpm<-cor(df2$gpp_obs,df2$VPM,method = "pearson")^2
  RMSE_vpm<-sqrt(SSE/length(df2$gpp_obs))
  RPE_vpm<-(mean(df2$VPM,na.rm = TRUE)-mean(df2$gpp_obs,na.rm = TRUE))/mean(df2$gpp_obs,na.rm = TRUE)
  df_metrics<- rbind(df_metrics[,],c(sitenames$SITE_ID[i],sitenames$IGBP[i],sitenames$LOCATION_LAT[i],sitenames$LOCATION_LONG[i],"VPM",R2_vpm,RMSE_vpm,RPE_vpm))
  
  #GOL_PEM metrics
  SSE<-sum((df2$gpp_obs-df2$GLO_PEM)^2,na.rm=TRUE)
  SST<-sum((df2$gpp_obs-mean(df2$gpp_obs,na.rm = TRUE))^2,na.rm=TRUE)
  R2_GLO_PEM<-cor(df2$gpp_obs,df2$GLO_PEM,method = "pearson")^2
  RMSE_GLO_PEM<-sqrt(SSE/length(df2$gpp_obs))
  RPE_GLO_PEM<-(mean(df2$GLO_PEM,na.rm = TRUE)-mean(df2$gpp_obs,na.rm = TRUE))/mean(df2$gpp_obs,na.rm = TRUE)
  df_metrics<- rbind(df_metrics[,],c(sitenames$SITE_ID[i],sitenames$IGBP[i],sitenames$LOCATION_LAT[i],sitenames$LOCATION_LONG[i],"GLO_PEM",R2_GLO_PEM,RMSE_GLO_PEM,RPE_GLO_PEM))
  
  #EC metrics
  SSE<-sum((df2$gpp_obs-df2$EC)^2,na.rm=TRUE)
  SST<-sum((df2$gpp_obs-mean(df2$gpp_obs,na.rm = TRUE))^2,na.rm=TRUE)
  R2_ec<-cor(df2$gpp_obs,df2$EC,method = "pearson")^2
  RMSE_ec<-sqrt(SSE/length(df2$gpp_obs))
  RPE_ec<-(mean(df2$EC,na.rm = TRUE)-mean(df2$gpp_obs,na.rm = TRUE))/mean(df2$gpp_obs,na.rm = TRUE)
  df_metrics<-rbind(df_metrics[,],c(sitenames$SITE_ID[i],sitenames$IGBP[i],sitenames$LOCATION_LAT[i],sitenames$LOCATION_LONG[i],"EC_LUE",R2_ec,RMSE_ec,RPE_ec))
  
  #chj metrics
  SSE<-sum((df2$gpp_obs-df2$CHJ)^2,na.rm=TRUE)
  SST<-sum((df2$gpp_obs-mean(df2$gpp_obs,na.rm = TRUE))^2,na.rm=TRUE)
  R2_chj<-cor(df2$gpp_obs,df2$CHJ,method = "pearson")^2
  RMSE_chj<-sqrt(SSE/length(df2$gpp_obs))
  RPE_chj<-(mean(df2$CHJ,na.rm = TRUE)-mean(df2$gpp_obs,na.rm = TRUE))/mean(df2$gpp_obs,na.rm = TRUE)
  df_metrics<-rbind(df_metrics[,],c(sitenames$SITE_ID[i],sitenames$IGBP[i],sitenames$LOCATION_LAT[i],sitenames$LOCATION_LONG[i],"CHJ",R2_chj,RMSE_chj,RPE_chj))
  
  #C_Fix metrics
  SSE<-sum((df2$gpp_obs-df2$C_Fix)^2,na.rm=TRUE)
  SST<-sum((df2$gpp_obs-mean(df2$gpp_obs,na.rm = TRUE))^2,na.rm=TRUE)
  R2_C_Fix<-cor(df2$gpp_obs,df2$C_Fix,method = "pearson")^2
  RMSE_C_Fix<-sqrt(SSE/length(df2$gpp_obs))
  RPE_C_Fix<-(mean(df2$C_Fix,na.rm = TRUE)-mean(df2$gpp_obs,na.rm = TRUE))/mean(df2$gpp_obs,na.rm = TRUE)
  df_metrics<-rbind(df_metrics[,],c(sitenames$SITE_ID[i],sitenames$IGBP[i],sitenames$LOCATION_LAT[i],sitenames$LOCATION_LONG[i],"C_Fix",R2_C_Fix,RMSE_C_Fix,RPE_C_Fix))
  
  
  #modis metrics
  SSE<-sum((df2$gpp_obs-df2$gpp_mod)^2,na.rm=TRUE)
  SST<-sum((df2$gpp_obs-mean(df2$gpp_obs,na.rm = TRUE))^2,na.rm=TRUE)
  R2_mod<-cor(df2$gpp_obs,df2$gpp_mod,method = "pearson")^2
  RMSE_mod<-sqrt(SSE/length(df2$gpp_obs))
  RPE_mod<-(mean(df2$gpp_mod,na.rm = TRUE)-mean(df2$gpp_obs,na.rm = TRUE))/mean(df2$gpp_obs,na.rm = TRUE)
  df_metrics<-rbind(df_metrics[,],c(sitenames$SITE_ID[i],sitenames$IGBP[i],sitenames$LOCATION_LAT[i],sitenames$LOCATION_LONG[i],"MODIS",R2_mod,RMSE_mod,RPE_mod))
  
  ## fusing with RF, 5 fold cross validation
  library(randomForest)
  rf_mod<-randomForest(gpp_obs ~ VPM+GLO_PEM+EC+CHJ+C_Fix,data =df2)
  rf_mod_pre<-predict(rf_mod,df2)
  df2$RF<-rf_mod_pre
  n<-length(df2$gpp_obs)
  Z=5
  mm<-CV(n,Z)
  R2_rf<-rep(0,Z)
  RMSE_rf<-rep(0,Z)
  RPE_rf<-rep(0,Z)
  for(j in 1:Z){
    m=mm[[j]]
    validate_data<-df2[m,]
    train_data<-df2[-m,]
    rf_mod<-randomForest(gpp_obs ~ VPM+GLO_PEM+EC+CHJ+C_Fix,data =train_data)
    rf_mod_pre<-predict(rf_mod,validate_data)
    validate_data$rf<-rf_mod_pre
    SSE<-sum((validate_data$gpp_obs-validate_data$rf)^2,na.rm=TRUE)
    SST<-sum((validate_data$gpp_obs-mean(validate_data$gpp_obs))^2,na.rm=TRUE)
    R2_rf[j]<-cor(validate_data$gpp_obs,validate_data$rf,method="pearson")^2
    RMSE_rf[j]<-sqrt(SSE/length(validate_data$gpp_obs))
    #Relative predictive error (RPE)
    RPE_rf[j]<-(mean(validate_data$rf)-mean(validate_data$gpp_obs))/mean(validate_data$gpp_obs)
  }
  c(mean(R2_rf),mean(RMSE_rf),mean(RPE_rf))
  df_metrics<-rbind(df_metrics[,],c(sitenames$SITE_ID[i],sitenames$IGBP[i],sitenames$LOCATION_LAT[i],sitenames$LOCATION_LONG[i],"RF",c(mean(R2_rf),mean(RMSE_rf),mean(RPE_rf))))
  
  # fusing with SVM, 5 fold cross validation
  library(e1071)
  svm_mod<-svm(gpp_obs ~ VPM+GLO_PEM+EC+CHJ+C_Fix, data=df2)
  svm_mod_pre<-predict(svm_mod,df2)
  df2$SVM<-svm_mod_pre
  n<-length(df2$gpp_obs)
  Z=5
  mm<-CV(n,Z)
  R2_svm<-rep(0,Z)
  RMSE_svm<-rep(0,Z)
  RPE_svm<-rep(0,Z)
  for(j in 1:Z){
    m=mm[[j]]
    validate_data<-df2[m,]
    train_data<-df2[-m,]
    svm_mod<-svm(gpp_obs ~ VPM+GLO_PEM+EC+CHJ+C_Fix, data=train_data)
    svm_mod_pre<-predict(svm_mod,validate_data)
    validate_data$svm<-svm_mod_pre
    SSE<-sum((validate_data$gpp_obs-validate_data$svm)^2,na.rm=TRUE)
    SST<-sum((validate_data$gpp_obs-mean(validate_data$gpp_obs))^2,na.rm=TRUE)
    R2_svm[j]<-cor(validate_data$gpp_obs,validate_data$svm,method="pearson")^2
    RMSE_svm[j]<-sqrt(SSE/length(validate_data$gpp_obs))
    #Relative predictive error (RPE)
    RPE_svm[j]<-(mean(validate_data$svm)-mean(validate_data$gpp_obs))/mean(validate_data$gpp_obs)
  }
  c(mean(R2_svm),mean(RMSE_svm),mean(RPE_svm))
  df_metrics<-rbind(df_metrics[,],c(sitenames$SITE_ID[i],sitenames$IGBP[i],sitenames$LOCATION_LAT[i],sitenames$LOCATION_LONG[i],"SVM",c(mean(R2_svm),mean(RMSE_svm),mean(RPE_svm))))

  # fusing with BMA, 5 fold cross validation
  library("BMA")
  bma_mod<-bic.glm(gpp_obs ~ VPM+GLO_PEM+EC+CHJ+C_Fix, data=df2,
                   glm.family=gaussian(link="identity"))
  bma_mod_pre<-predict(bma_mod,df2)
  df2$BMA<-bma_mod_pre
  n<-length(df2$gpp_obs)
  Z=5
  mm<-CV(n,Z)
  R2_bma<-rep(0,Z)
  RMSE_bma<-rep(0,Z)
  RPE_bma<-rep(0,Z)
  for(j in 1:Z){
    m=mm[[j]]
    validate_data<-df2[m,]
    train_data<-df2[-m,]
    bma_mod<-bic.glm(gpp_obs ~ VPM+GLO_PEM+EC+CHJ+C_Fix, data=train_data,
                     glm.family=gaussian(link="identity"))
    bma_mod_pre<-predict(bma_mod,validate_data)
    validate_data$bma<-bma_mod_pre
    SSE<-sum((validate_data$gpp_obs-validate_data$bma)^2,na.rm=TRUE)
    SST<-sum((validate_data$gpp_obs-mean(validate_data$gpp_obs))^2,na.rm=TRUE)
    R2_bma[j]<-cor(validate_data$gpp_obs,validate_data$bma,method="pearson")^2
    RMSE_bma[j]<-sqrt(SSE/length(validate_data$gpp_obs))
    #Relative predictive error (RPE)
    RPE_bma[j]<-(mean(validate_data$bma)-mean(validate_data$gpp_obs))/mean(validate_data$gpp_obs)
  }
  c(mean(R2_bma),mean(RMSE_bma),mean(RPE_bma))
  df_metrics<-rbind(df_metrics[,],c(sitenames$SITE_ID[i],sitenames$IGBP[i],sitenames$LOCATION_LAT[i],sitenames$LOCATION_LONG[i],"BMA",c(mean(R2_bma),mean(RMSE_bma),mean(RPE_bma))))
  
  #write the results for each site
  write_name<-paste(flux_path,"\\",sitenames$SITE_ID[i],"_model_gpp_DTNT.csv",sep="")
  write.csv(df2,write_name)
}

# write the metrics
names(df_metrics)<-c("site","PFT","Lat","Lon","model","R2","RMSE","RPE")
View(df_metrics)
write.csv(df_metrics,"sites_sta_DTNT.csv")



