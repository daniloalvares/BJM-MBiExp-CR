format_change <- function(data){
 
  # Data by biomarker
  data_M <- data[which(!is.na(data$m_spike_serum)),]
  data_F <- data[which(!is.na(data$flc_serum)),]
  
  # Unique ID
  u_ID <- unique(data$patientid)
  
  # Start and stop positions of longitudinal measurements per patient
  start_M <- stop_M <- rep(NA,length(u_ID))
  start_F <- stop_F <- rep(NA,length(u_ID))
  aux_M <- aux_F <- 0
  
  # Variables for short/survival format
  os <- surv <- status_os <- status <- lastLoT <- rep(NA,length(u_ID))
  sex <- race <- ecog <- iss <- rep(NA,length(u_ID))
  age <- albumin <- beta2_microglobulin <- creatinine <- hemoglobin <- rep(NA,length(u_ID))
  ldh <- lymphocyte <- neutrophil <- platelet <- rep(NA,length(u_ID))
  immunoglobulin_a <- immunoglobulin_g <- immunoglobulin_m <- rep(NA,length(u_ID))
  
  for(i in 1:length(u_ID)){
    pos <- which(data$patientid==u_ID[i])
    pos_M <- which(data_M$patientid==u_ID[i])
    pos_F <- which(data_F$patientid==u_ID[i])
    
    # M-spike
    if(length(pos_M)>0){
      start_M[i] <- aux_M + 1
      stop_M[i] <- aux_M + length(pos_M)
      aux_M <- stop_M[i]
    }
    
    # FLC
    if(length(pos_F)>0){
      start_F[i] <- aux_F + 1
      stop_F[i] <- aux_F + length(pos_F)
      aux_F <- stop_F[i]
    }
    
    # Survival information
    os[i] <- data$os[pos][1]
    surv[i] <- data$surv[pos][1]
    status_os[i] <- data$status_os[pos][1]
    status[i] <- data$status_lot[pos][1]
    lastLoT[i] <- data$lastLoT[pos][1]
    
    # Categorical baseline variables
    sex[i] <- data$gender[pos][1]
    race[i] <- data$race_ethnicity[pos][1]
    ecog[i] <- data$ecoggrp[pos][1]
    iss[i] <- data$iss_stage[pos][1]
    
    # Continuous baseline variables
    age[i] <- data$age[pos][1]
    albumin[i] <- data$albumin_bl[pos][1]
    beta2_microglobulin[i] <- data$beta2_microglobulin_bl[pos][1]
    creatinine[i] <- data$creatinine_bl[pos][1]
    hemoglobin[i] <- data$hemoglobin_bl[pos][1]
    ldh[i] <- data$ldh_bl[pos][1]
    lymphocyte[i] <- data$lymphocyte_bl[pos][1]
    neutrophil[i] <- data$neutrophil_count_bl[pos][1]
    platelet[i] <- data$platelet_count_bl[pos][1]
    immunoglobulin_a[i] <- data$immunoglobulin_a_bl[pos][1]
    immunoglobulin_g[i] <- data$immunoglobulin_g_bl[pos][1]
    immunoglobulin_m[i] <- data$immunoglobulin_m_bl[pos][1]
    
  }

  
  # Reference: Male
  sex <- ifelse(sex=="Female",1,0)
  
  # Reference: Non-Hispanic White
  race_0 <- ifelse(race=="Non-Hispanic White",1,0)
  race_1 <- ifelse(race=="Non-Hispanic Black",1,0)
  race_2 <- ifelse(race=="Hispanic or Latino (any)" | race=="Non-Hispanic Asian" | race=="Other Non-Hispanic",1,0)
  race_3 <- ifelse(race=="Not reported",1,0)
  
  # Reference: 0
  ecog_0 <- ifelse(ecog=="0",1,0)
  ecog_1 <- ifelse(ecog=="1",1,0)
  ecog_2 <- ifelse(ecog=="2" | ecog=="3+",1,0)
  ecog_3 <- ifelse(ecog=="Not reported",1,0)
  
  # Reference: Stage I
  iss_0 <- ifelse(iss=="Stage I",1,0)
  iss_1 <- ifelse(iss=="Stage II",1,0)
  iss_2 <- ifelse(iss=="Stage III",1,0)
  iss_3 <- ifelse(iss=="Unknown/not documented",1,0)


  Short <- data.frame(patientid=data[!duplicated(data$patientid),"patientid"],
                      ID1=data[!duplicated(data$patientid),"ID1"],
                      ID2=data[!duplicated(data$patientid),"ID2"],
                      os=os/365.25, surv=surv/365.25, status_os=status_os, status=status, lastLoT=lastLoT,
                      sex=sex, race_0=race_0, race_1=race_1, 
                      race_2=race_2, race_3=race_3, ecog_0=ecog_0, ecog_1=ecog_1, ecog_2=ecog_2, 
                      ecog_3=ecog_3, iss_0=iss_0, iss_1=iss_1, iss_2=iss_2, iss_3=iss_3,
                      age=age, albumin=albumin, beta2_microglobulin=beta2_microglobulin,
                      creatinine=creatinine, hemoglobin=hemoglobin, ldh=ldh, lymphocyte=lymphocyte,
                      neutrophil=neutrophil, platelet=platelet, immunoglobulin_a=immunoglobulin_a, 
                      immunoglobulin_g=immunoglobulin_g, immunoglobulin_m=immunoglobulin_m,
                      start_M=start_M, stop_M=stop_M, start_F=start_F, stop_F=stop_F)


  # Longitudinal information
  M_Spike <- data.frame(patientid=data_M$patientid, ID1=data_M$ID1, ID2=data_M$ID2,
                        y=data_M$m_spike_serum, time=as.numeric(data_M$timefromlot)/365.25)
  FLC <- data.frame(patientid=data_F$patientid, ID1=data_F$ID1, ID2=data_F$ID2,
                    y=data_F$flc_serum, time=as.numeric(data_F$timefromlot)/365.25)

  Long <- list(M_Spike=M_Spike, FLC=FLC)
  
  return( list(Short=Short, Long=Long) )
  
}


# Test sets
testlot1 <- format_change(testlot1)
testlot2 <- format_change(testlot2)
testlot3 <- format_change(testlot3)
testlot4 <- format_change(testlot4)

# Training sets
trainlot1 <- format_change(trainlot1)
trainlot2 <- format_change(trainlot2)
trainlot3 <- format_change(trainlot3)
trainlot4 <- format_change(trainlot4)


# Previous LoT duration
trainlot2$Short$duration <- trainlot1$Short$surv[which(!is.na(trainlot1$Short$ID2))]
trainlot3$Short$duration <- trainlot2$Short$surv[which(!is.na(trainlot2$Short$ID2))]
trainlot4$Short$duration <- trainlot3$Short$surv[which(!is.na(trainlot3$Short$ID2))]
testlot2$Short$duration <- testlot1$Short$surv[which(!is.na(testlot1$Short$ID2))]
testlot3$Short$duration <- testlot2$Short$surv[which(!is.na(testlot2$Short$ID2))]
testlot4$Short$duration <- testlot3$Short$surv[which(!is.na(testlot3$Short$ID2))]


rm(list=setdiff(ls(), c("testlot1","testlot2","testlot3","testlot4","trainlot1","trainlot2","trainlot3","trainlot4")))
