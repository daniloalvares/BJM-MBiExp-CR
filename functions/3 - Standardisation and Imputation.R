vars <- c("age","albumin_bl","beta2_microglobulin_bl","creatinine_bl","hemoglobin_bl","ldh_bl","lymphocyte_bl",
          "neutrophil_count_bl","platelet_count_bl","immunoglobulin_a_bl","immunoglobulin_g_bl","immunoglobulin_m_bl")


# LOG SCALE
# Training set
trainlot1[,vars] <- log(trainlot1[,vars] + 0.1)
trainlot2[,vars] <- log(trainlot2[,vars] + 0.1)
trainlot3[,vars] <- log(trainlot3[,vars] + 0.1)
trainlot4[,vars] <- log(trainlot4[,vars] + 0.1)
# Test set
testlot1[,vars] <- log(testlot1[,vars] + 0.1)
testlot2[,vars] <- log(testlot2[,vars] + 0.1)
testlot3[,vars] <- log(testlot3[,vars] + 0.1)
testlot4[,vars] <- log(testlot4[,vars] + 0.1)


for(vv in vars){
    # Mean and standard deviation based on training set in LoT 1
    mean_bl <- mean(trainlot1[!duplicated(trainlot1$patientid),vv], na.rm = TRUE)
    sd_bl <- sd(trainlot1[!duplicated(trainlot1$patientid),vv], na.rm = TRUE)
    
    # STANDARDISATION
    # Training set
    trainlot1[,vv] <- (trainlot1[,vv] - mean_bl) / sd_bl
    trainlot2[,vv] <- (trainlot2[,vv] - mean_bl) / sd_bl
    trainlot3[,vv] <- (trainlot3[,vv] - mean_bl) / sd_bl
    trainlot4[,vv] <- (trainlot4[,vv] - mean_bl) / sd_bl
    # Test set
    testlot1[,vv] <- (testlot1[,vv] - mean_bl) / sd_bl
    testlot2[,vv] <- (testlot2[,vv] - mean_bl) / sd_bl
    testlot3[,vv] <- (testlot3[,vv] - mean_bl) / sd_bl
    testlot4[,vv] <- (testlot4[,vv] - mean_bl) / sd_bl
    
    # IMPUTATION
    # Training set
    trainlot1[is.na(trainlot1[,vv]),vv] <- 0
    trainlot2[is.na(trainlot2[,vv]),vv] <- 0
    trainlot3[is.na(trainlot3[,vv]),vv] <- 0
    trainlot4[is.na(trainlot4[,vv]),vv] <- 0
    # Test set
    testlot1[is.na(testlot1[,vv]),vv] <- 0
    testlot2[is.na(testlot2[,vv]),vv] <- 0
    testlot3[is.na(testlot3[,vv]),vv] <- 0
    testlot4[is.na(testlot4[,vv]),vv] <- 0
}


rm(list=setdiff(ls(), c("testlot1","testlot2","testlot3","testlot4","trainlot1","trainlot2","trainlot3","trainlot4")))
