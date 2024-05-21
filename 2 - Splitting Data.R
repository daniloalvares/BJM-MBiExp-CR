u_ID1 <- unique(datalot1$patientid)
u_ID2 <- unique(datalot2$patientid)
u_ID3 <- unique(datalot3$patientid)
u_ID4 <- unique(datalot4$patientid)

# Selecting patients for test set from LoT 4
set.seed(4); sel_test_lot4 <- sample(u_ID4, size = prop * length(u_ID4))

# Selecting patients for test set from LoT 3
u_ID3wt4 <- u_ID3[!(u_ID3 %in% u_ID4)]
size3 <- round(length(u_ID3) * prop - length(sel_test_lot4))
set.seed(3); sel_test_lot3 <- sort(c(sample(u_ID3wt4, size = size3), sel_test_lot4))

# Selecting patients for test set from LoT 2
u_ID2wt34 <- u_ID2[!(u_ID2 %in% u_ID3)]
size2 <- round(length(u_ID2) * prop - length(sel_test_lot3))
set.seed(2); sel_test_lot2 <- sort(c(sample(u_ID2wt34, size = size2), sel_test_lot3))

# Selecting patients for test set from LoT 1
u_ID1wt234 <- u_ID1[!(u_ID1 %in% u_ID2)]
size1 <- round(length(u_ID1) * prop - length(sel_test_lot2))
set.seed(1); sel_test_lot1 <- sort(c(sample(u_ID1wt234, size = size1), sel_test_lot2))


# Test sets
testlot1 <- datalot1[(datalot1$patientid %in% sel_test_lot1),]
testlot2 <- datalot2[(datalot2$patientid %in% sel_test_lot2),]
testlot3 <- datalot3[(datalot3$patientid %in% sel_test_lot3),]
testlot4 <- datalot4[(datalot4$patientid %in% sel_test_lot4),]

# Training sets
trainlot1 <- datalot1[!(datalot1$patientid %in% sel_test_lot1),]
trainlot2 <- datalot2[!(datalot2$patientid %in% sel_test_lot2),]
trainlot3 <- datalot3[!(datalot3$patientid %in% sel_test_lot3),]
trainlot4 <- datalot4[!(datalot4$patientid %in% sel_test_lot4),]


# Reference ID
testlot1 <- as.data.frame(testlot1 %>% group_by(patientid) %>% dplyr::mutate(ID1 = cur_group_id()))
testlot2 <- as.data.frame(testlot2 %>% group_by(patientid) %>% dplyr::mutate(ID1 = cur_group_id()))
testlot3 <- as.data.frame(testlot3 %>% group_by(patientid) %>% dplyr::mutate(ID1 = cur_group_id()))
testlot4 <- as.data.frame(testlot4 %>% group_by(patientid) %>% dplyr::mutate(ID1 = cur_group_id()))
trainlot1 <- as.data.frame(trainlot1 %>% group_by(patientid) %>% dplyr::mutate(ID1 = cur_group_id()))
trainlot2 <- as.data.frame(trainlot2 %>% group_by(patientid) %>% dplyr::mutate(ID1 = cur_group_id()))
trainlot3 <- as.data.frame(trainlot3 %>% group_by(patientid) %>% dplyr::mutate(ID1 = cur_group_id()))
trainlot4 <- as.data.frame(trainlot4 %>% group_by(patientid) %>% dplyr::mutate(ID1 = cur_group_id()))


# LoT 2 IDs that are also on LoT 1
u_test_LoT2 <- unique(testlot2$patientid)
u_train_LoT2 <- unique(trainlot2$patientid)
testlot1$ID2 <- NA
trainlot1$ID2 <- NA
    
aux_lot1 <- testlot1[(testlot1$patientid %in% u_test_LoT2),]
aux_lot1 <- (aux_lot1 %>% group_by(patientid) %>% dplyr::mutate(ID2 = cur_group_id()))$ID2
testlot1$ID2[(testlot1$patientid %in% u_test_LoT2)] <- aux_lot1

aux_lot1 <- trainlot1[(trainlot1$patientid %in% u_train_LoT2),]
aux_lot1 <- (aux_lot1 %>% group_by(patientid) %>% dplyr::mutate(ID2 = cur_group_id()))$ID2
trainlot1$ID2[(trainlot1$patientid %in% u_train_LoT2)] <- aux_lot1

# LoT 3 IDs that are also on LoT 2
u_test_LoT3 <- unique(testlot3$patientid)
u_train_LoT3 <- unique(trainlot3$patientid)
testlot2$ID2 <- NA
trainlot2$ID2 <- NA

aux_lot2 <- testlot2[(testlot2$patientid %in% u_test_LoT3),]
aux_lot2 <- (aux_lot2 %>% group_by(patientid) %>% dplyr::mutate(ID2 = cur_group_id()))$ID2
testlot2$ID2[(testlot2$patientid %in% u_test_LoT3)] <- aux_lot2

aux_lot2 <- trainlot2[(trainlot2$patientid %in% u_train_LoT3),]
aux_lot2 <- (aux_lot2 %>% group_by(patientid) %>% dplyr::mutate(ID2 = cur_group_id()))$ID2
trainlot2$ID2[(trainlot2$patientid %in% u_train_LoT3)] <- aux_lot2

# LoT 4 IDs that are also on LoT 3
u_test_LoT4 <- unique(testlot4$patientid)
u_train_LoT4 <- unique(trainlot4$patientid)
testlot3$ID2 <- NA
trainlot3$ID2 <- NA

aux_lot3 <- testlot3[(testlot3$patientid %in% u_test_LoT4),]
aux_lot3 <- (aux_lot3 %>% group_by(patientid) %>% dplyr::mutate(ID2 = cur_group_id()))$ID2
testlot3$ID2[(testlot3$patientid %in% u_test_LoT4)] <- aux_lot3

aux_lot3 <- trainlot3[(trainlot3$patientid %in% u_train_LoT4),]
aux_lot3 <- (aux_lot3 %>% group_by(patientid) %>% dplyr::mutate(ID2 = cur_group_id()))$ID2
trainlot3$ID2[(trainlot3$patientid %in% u_train_LoT4)] <- aux_lot3

testlot4$ID2 <- NA
trainlot4$ID2 <- NA


rm(list=setdiff(ls(), c("testlot1","testlot2","testlot3","testlot4","trainlot1","trainlot2","trainlot3","trainlot4")))
