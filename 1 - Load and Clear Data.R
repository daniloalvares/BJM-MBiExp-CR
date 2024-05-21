# Establish connection with Redshift
con <- DBI::dbConnect(odbc::odbc(), Driver="redshift",              
                      Server="redshift-01-eu.dap.apollo.roche.com",
                      Port="5439",
                      Database="flatiron_edm",
                      UID=Sys.getenv("REDSHIFT_USERNAME"),
                      PWD=Sys.getenv("REDSHIFT_PASSWORD"))

modeldata <- dbGetQuery(con, "SELECT * FROM scr_scr_361_joint_models.modeldata_05012023")
modeldata <- modeldata %>% filter(ismaintenancetherapy == 0 & is_transplant == 0)
modeldata <- modeldata[with(modeldata,order(patientid,testdate)),]


# Selecting variables that will be used
modeldata <- modeldata[,c("patientid","linename","linenumber","LineOfTherapy","testdate","start_lot_date","end_lot_date",
                          "first_trt_date","last_trt_date","dateofdeath","dateofdeath_imp","death_indicator",
                          "m_spike_serum","kappa_flc_serum","lambda_flc_serum",
                          "gender","race_ethnicity","ecoggrp","iss_stage",
                          "age","albumin_bl","beta2_microglobulin_bl","creatinine_bl","hemoglobin_bl",
                          "ldh_bl","lymphocyte_bl","neutrophil_count_bl","platelet_count_bl",
                          "immunoglobulin_a_bl","immunoglobulin_g_bl","immunoglobulin_m_bl")]


# g/dL -> g/L
modeldata$kappa_flc_serum <- modeldata$kappa_flc_serum/10
modeldata$lambda_flc_serum <- modeldata$lambda_flc_serum/10
modeldata$flc_serum <- NA


# Combine Kappa-FLC and Lambda-FLC into a single biomarker
# Whatever the higher (first observation) value in LoT 1 to follow through
uniqueID <- unique(modeldata$patientid)
for(i in 1:length(uniqueID)){
  pos1 <- which(modeldata$patientid==uniqueID[i] & modeldata$linenumber==1)
  pos2 <- which(modeldata$patientid==uniqueID[i] & modeldata$linenumber==2)
  pos3 <- which(modeldata$patientid==uniqueID[i] & modeldata$linenumber==3)
  pos4 <- which(modeldata$patientid==uniqueID[i] & modeldata$linenumber==4)
  pick1 <- NULL
  
  # LoT 1
  if(length(pos1)>0){
    bio_K <- sum(is.na(modeldata$kappa_flc_serum[pos1])) != length(pos1)
    bio_L <- sum(is.na(modeldata$lambda_flc_serum[pos1])) != length(pos1)
    
    if(bio_K & bio_L){
      first_K <- modeldata$kappa_flc_serum[pos1][!is.na(modeldata$kappa_flc_serum[pos1])][1]
      first_L <- modeldata$lambda_flc_serum[pos1][!is.na(modeldata$lambda_flc_serum[pos1])][1]
      pick1 <- ifelse(first_K >= first_L, "kappa_flc_serum", "lambda_flc_serum")
    }else{
      if(bio_K){ pick1 <- "kappa_flc_serum" }
      else{ 
        if(bio_L){ pick1 <- "lambda_flc_serum" } 
      }
    }
    if(!is.null(pick1)){ modeldata$flc_serum[pos1] <- modeldata[pos1,pick1] }
  }
  
  if(!is.null(pick1)){
    modeldata$flc_serum[pos2] <- modeldata[pos2,pick1]
    modeldata$flc_serum[pos3] <- modeldata[pos3,pick1]
    modeldata$flc_serum[pos4] <- modeldata[pos4,pick1]
  }else{
    pick2 <- NULL
    
    # LoT 2
    if(length(pos2)>0){
      bio_K <- sum(is.na(modeldata$kappa_flc_serum[pos2])) != length(pos2)
      bio_L <- sum(is.na(modeldata$lambda_flc_serum[pos2])) != length(pos2)
      
      if(bio_K & bio_L){
        first_K <- modeldata$kappa_flc_serum[pos2][!is.na(modeldata$kappa_flc_serum[pos2])][1]
        first_L <- modeldata$lambda_flc_serum[pos2][!is.na(modeldata$lambda_flc_serum[pos2])][1]
        pick2 <- ifelse(first_K >= first_L, "kappa_flc_serum", "lambda_flc_serum")
      }else{
        if(bio_K){ pick2 <- "kappa_flc_serum" }
        else{ 
          if(bio_L){ pick2 <- "lambda_flc_serum" } 
        }
      }
      if(!is.null(pick2)){ modeldata$flc_serum[pos2] <- modeldata[pos2,pick2] }
    }
    
    if(!is.null(pick2)){
      modeldata$flc_serum[pos3] <- modeldata[pos3,pick2]
      modeldata$flc_serum[pos4] <- modeldata[pos4,pick2]
    }else{
      pick3 <- NULL
      
      # LoT 3
      if(length(pos3)>0){
        bio_K <- sum(is.na(modeldata$kappa_flc_serum[pos3])) != length(pos3)
        bio_L <- sum(is.na(modeldata$lambda_flc_serum[pos3])) != length(pos3)
        
        if(bio_K & bio_L){
          first_K <- modeldata$kappa_flc_serum[pos3][!is.na(modeldata$kappa_flc_serum[pos3])][1]
          first_L <- modeldata$lambda_flc_serum[pos3][!is.na(modeldata$lambda_flc_serum[pos3])][1]
          pick3 <- ifelse(first_K >= first_L, "kappa_flc_serum", "lambda_flc_serum")
        }else{
          if(bio_K){ pick3 <- "kappa_flc_serum" }
          else{ 
            if(bio_L){ pick3 <- "lambda_flc_serum" } 
          }
        }
        if(!is.null(pick3)){ modeldata$flc_serum[pos3] <- modeldata[pos3,pick3] }
      }
      
      if(!is.null(pick3)){
        modeldata$flc_serum[pos4] <- modeldata[pos4,pick3]
      }else{
        pick4 <- NULL
        
        # LoT 4
        if(length(pos4)>0){
          bio_K <- sum(is.na(modeldata$kappa_flc_serum[pos4])) != length(pos4)
          bio_L <- sum(is.na(modeldata$lambda_flc_serum[pos4])) != length(pos4)
          
          if(bio_K & bio_L){
            first_K <- modeldata$kappa_flc_serum[pos4][!is.na(modeldata$kappa_flc_serum[pos4])][1]
            first_L <- modeldata$lambda_flc_serum[pos4][!is.na(modeldata$lambda_flc_serum[pos4])][1]
            pick4 <- ifelse(first_K >= first_L, "kappa_flc_serum", "lambda_flc_serum")
          }else{
            if(bio_K){ pick4 <- "kappa_flc_serum" }
            else{ 
              if(bio_L){ pick4 <- "lambda_flc_serum" } 
            }
          }
          if(!is.null(pick4)){ modeldata$flc_serum[pos4] <- modeldata[pos4,pick4] }
        }
      }
    }
  }
}


# Cut-off: 99th percentile
cutoff_m_spike_serum <- quantile(modeldata$m_spike_serum, probs = 0.99, na.rm = TRUE)
# M-spikes that exceed the threshold
modeldata <- modeldata[-which(!is.na(modeldata$m_spike_serum) & modeldata$m_spike_serum > cutoff_m_spike_serum),]

cutoff_flc_serum <- quantile(modeldata$flc_serum, probs = 0.99, na.rm = TRUE)
# FLC that exceed the threshold
modeldata <- modeldata[-which(!is.na(modeldata$flc_serum) & modeldata$flc_serum > cutoff_flc_serum),]


# Check LoT jumps
uniqueID <- unique(modeldata$patientid)
for(i in 1:length(uniqueID)){
    pos <- which(modeldata$patientid==uniqueID[i])
    if(sum(modeldata$linenumber[pos]==1)==0){
        modeldata <- modeldata[-pos,]   
    }else{
        pos1 <- which(modeldata$linenumber[pos]==1)
        if(sum(modeldata$linenumber[pos]==2)==0){
            if(length(pos1) < length(pos)){ modeldata <- modeldata[-pos[-pos1],] }
        }else{
            pos2 <- which(modeldata$linenumber[pos]==2)
            if(sum(modeldata$linenumber[pos]==3)==0){
                if(length(c(pos1,pos2)) < length(pos)){ modeldata <- modeldata[-pos[-c(pos1,pos2)],] }
            }else{
                pos3 <- which(modeldata$linenumber[pos]==3)
                if(sum(modeldata$linenumber[pos]==4)==0){
                    if(length(c(pos1,pos2,pos3)) < length(pos)){ modeldata <- modeldata[-pos[-c(pos1,pos2,pos3)],] }
                }
            }
        }
    }
}


# Last LoT identifier and date of death
uniqueID <- unique(modeldata$patientid)
modeldata$lastLoT <- 0
for(i in 1:length(uniqueID)){
    pos <- which(modeldata$patientid==uniqueID[i])
    pos1 <- which(modeldata$linenumber[pos]==max(modeldata$linenumber[pos]))
    modeldata$lastLoT[pos][pos1] <- 1
    if(!is.na(modeldata$dateofdeath[pos][pos1][1]) & (modeldata$dateofdeath[pos][pos1][1]==format(modeldata$end_lot_date[pos][pos1][1],format="%Y-%m"))){
        modeldata$dateofdeath_imp[pos] <- modeldata$end_lot_date[pos][pos1][1]
    }
}


# Redefining the end of LoT
uniqueID <- unique(modeldata$patientid)
for(i in 1:length(uniqueID)){
    pos <- which(modeldata$patientid==uniqueID[i])
    if(max(modeldata$linenumber[pos]) > 1){
        modeldata$end_lot_date[pos][modeldata$linenumber[pos] == 1] <- min(modeldata$start_lot_date[pos][modeldata$linenumber[pos] == 2]) - 1
    }
    if(max(modeldata$linenumber[pos]) > 2){
        modeldata$end_lot_date[pos][modeldata$linenumber[pos] == 2] <- min(modeldata$start_lot_date[pos][modeldata$linenumber[pos] == 3]) - 1
    }
    if(max(modeldata$linenumber[pos]) > 3){
        modeldata$end_lot_date[pos][modeldata$linenumber[pos] == 3] <- min(modeldata$start_lot_date[pos][modeldata$linenumber[pos] == 4]) - 1
    }
    if(max(modeldata$linenumber[pos]) > 4){
        modeldata$end_lot_date[pos][modeldata$linenumber[pos] == 4] <- min(modeldata$start_lot_date[pos][modeldata$linenumber[pos] >= 5]) - 1
    }
}


# LoTs 1, 2, and 3
# 0: Censoring; 1: Death during LoT; 2: Death out of LoT
modeldata$status_lot <- as.numeric((modeldata$lastLoT==1) & (modeldata$death_indicator=="Yes")) + 2*as.numeric(modeldata$lastLoT==0)
modeldata$status_lot[is.na(modeldata$status_lot)] <- 0
# LoT 4+
modeldata$status_lot[which(modeldata$linenumber==4)] <- ifelse(modeldata$death_indicator[which(modeldata$linenumber==4)]=="Yes", 1, 0)


# LoT 1
datalot1 <- modeldata[which(modeldata$linenumber == 1),]
# LoT 2
datalot2 <- modeldata[which(modeldata$linenumber == 2),]
# LoT 3
datalot3 <- modeldata[which(modeldata$linenumber == 3),]
# LoT 4+
datalot4 <- modeldata[which(modeldata$linenumber == 4),]


# Overall survival and time-to-event (competing risks)
datalot1$os <- datalot1$surv <- datalot1$dateofdeath_imp - datalot1$start_lot_date + 1
datalot1$os[is.na(datalot1$os)] <- datalot1$last_trt_date[is.na(datalot1$os)] - datalot1$start_lot_date[is.na(datalot1$os)] + 1
datalot1$surv[datalot1$status_lot == 0 | datalot1$status_lot == 2] <- datalot1$end_lot_date[datalot1$status_lot == 0 | datalot1$status_lot == 2] - datalot1$start_lot_date[datalot1$status_lot == 0 | datalot1$status_lot == 2] + 1
datalot1$status_os <- 0
datalot1$status_os[datalot1$status_lot == 1 | datalot1$status_lot == 2] <- 1

datalot2$os <- datalot2$surv <- datalot2$dateofdeath_imp - datalot2$start_lot_date + 1
datalot2$os[is.na(datalot2$os)] <- datalot2$last_trt_date[is.na(datalot2$os)] - datalot2$start_lot_date[is.na(datalot2$os)] + 1
datalot2$surv[datalot2$status_lot == 0 | datalot2$status_lot == 2] <- datalot2$end_lot_date[datalot2$status_lot == 0 | datalot2$status_lot == 2] - datalot2$start_lot_date[datalot2$status_lot == 0 | datalot2$status_lot == 2] + 1
datalot2$status_os <- 0
datalot2$status_os[datalot2$status_lot == 1 | datalot2$status_lot == 2] <- 1

datalot3$os <- datalot3$surv <- datalot3$dateofdeath_imp - datalot3$start_lot_date + 1
datalot3$os[is.na(datalot3$os)] <- datalot3$last_trt_date[is.na(datalot3$os)] - datalot3$start_lot_date[is.na(datalot3$os)] + 1
datalot3$surv[datalot3$status_lot == 0 | datalot3$status_lot == 2] <- datalot3$end_lot_date[datalot3$status_lot == 0 | datalot3$status_lot == 2] - datalot3$start_lot_date[datalot3$status_lot == 0 | datalot3$status_lot == 2] + 1
datalot3$status_os <- 0
datalot3$status_os[datalot3$status_lot == 1 | datalot3$status_lot == 2] <- 1

datalot4$os <- datalot4$dateofdeath_imp - datalot4$start_lot_date + 1
datalot4$os[is.na(datalot4$os)] <- datalot4$last_trt_date[is.na(datalot4$os)] - datalot4$start_lot_date[is.na(datalot4$os)] + 1
datalot4$surv <- datalot4$os
datalot4$status_os <- datalot4$status_lot


# Biomarker observation times since initiation of LoT
datalot1$timefromlot <- datalot1$testdate - datalot1$start_lot_date
datalot2$timefromlot <- datalot2$testdate - datalot2$start_lot_date
datalot3$timefromlot <- datalot3$testdate - datalot3$start_lot_date
datalot4$timefromlot <- datalot4$testdate - datalot4$start_lot_date


rm(list=setdiff(ls(), c("datalot1","datalot2","datalot3","datalot4")))
