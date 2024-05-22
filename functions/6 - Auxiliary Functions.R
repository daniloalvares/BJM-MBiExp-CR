# Fitting joint models using the joint specification
fit_jm <- function(data=data, LoT, iter=5000, warmup=1000, chains=3){
  
  vars <- c("sex","race_1","race_2","race_3",
            "ecog_1","ecog_2","ecog_3","iss_1","iss_2","iss_3",
            "age","albumin","beta2_microglobulin","creatinine",
            "hemoglobin","ldh","lymphocyte","neutrophil","platelet",
            "immunoglobulin_a","immunoglobulin_g","immunoglobulin_m")
  
  if(LoT==1){
    X <- data$Short[,vars]
  }else{
    X <- data$Short[,c(vars, "duration")]
  }
  
  # Model data
  N_M <- nrow(data$Long$M_Spike)
  N_F <- nrow(data$Long$FLC)
  n_M <- length(unique(data$Long$M_Spike$patientid))
  n_F <- length(unique(data$Long$FLC$patientid))
  y_M <- data$Long$M_Spike$y
  y_F <- data$Long$FLC$y
  times_M <- data$Long$M_Spike$time
  times_F <- data$Long$FLC$time
  n <- nrow(data$Short)
  nbeta <- ncol(X)
  ID_M1 <- (data$Long$M_Spike %>% group_by(patientid) %>% dplyr::mutate(ID = cur_group_id()))$ID
  ID_F1 <- (data$Long$FLC %>% group_by(patientid) %>% dplyr::mutate(ID = cur_group_id()))$ID
  ID_M2 <- data$Long$M_Spike$ID1
  ID_F2 <- data$Long$FLC$ID1
  Time <- data$Short$surv
  status_1 <- as.numeric((data$Short$status==0 & data$Short$lastLoT==0) | (data$Short$status==2 & data$Short$lastLoT==0))
  status_2 <- as.numeric(data$Short$status==1 | (data$Short$status==2 & data$Short$lastLoT==1))

  if(LoT != 4){
    # Setting initial values
    init_fun <- function(...){ 
      list(theta_M=c(0,0,0), theta_F=c(0,0,0), sigma2_M=1, sigma2_F=1,
           Omega_M=diag(rep(1,3)), Omega_F=diag(rep(1,3)),
           beta_1=rep(0,ncol(X)), beta_2=rep(0,ncol(X)),
           alphaB_1=c(0,0), alphaG_1=c(0,0), alphaD_1=c(0,0),
           alphaB_2=c(0,0), alphaG_2=c(0,0), alphaD_2=c(0,0),
           gamma_1=0, gamma_2=0, phi_1=1, phi_2=1,
           bi_M=matrix(0,nrow=n_M,ncol=3), bi_F=matrix(0,nrow=n_F,ncol=3))
    }
    
    options(mc.cores = parallel::detectCores())
    rstan_options(auto_write = TRUE)
    
    i.time <- Sys.time()
    fit <- suppressMessages(suppressWarnings( stan(model_code = jointmodel, 
                                                   data   = list(N_M=N_M, N_F=N_F, n_M=n_M, n_F=n_F, n=n, 
                                                                 nbeta=nbeta, y_M=log(y_M+1), y_F=log(y_F+1),
                                                                 ID_M1=ID_M1, ID_M2=ID_M2, ID_F1=ID_F1, ID_F2=ID_F2, 
                                                                 times_M=times_M, times_F=times_F, X=X,
                                                                 Time=Time, status_1=status_1, status_2=status_2),
                                                   init   = init_fun,
                                                   iter   = iter,
                                                   warmup = warmup,                 
                                                   chains = chains,
                                                   thin   = 1,
                                                   seed   = 123,
                                                   cores  = getOption("mc.cores",chains)) ))
    e.time <- Sys.time()
  }else{
    # Setting initial values
    init_fun <- function(...){ 
      list(theta_M=c(0,0,0), theta_F=c(0,0,0), sigma2_M=1, sigma2_F=1,
           Omega_M=diag(rep(1,3)), Omega_F=diag(rep(1,3)),
           beta_2=rep(0,ncol(X)),
           alphaB_2=c(0,0), alphaG_2=c(0,0), alphaD_2=c(0,0),
           gamma_2=0, phi_2=1,
           bi_M=matrix(0,nrow=n_M,ncol=3), bi_F=matrix(0,nrow=n_F,ncol=3))
    }
    
    options(mc.cores = parallel::detectCores())
    rstan_options(auto_write = TRUE)
    
    i.time <- Sys.time()
    fit <- suppressMessages(suppressWarnings( stan(model_code = jointmodelsurv, 
                                                   data   = list(N_M=N_M, N_F=N_F, n_M=n_M, n_F=n_F, n=n, 
                                                                 nbeta=nbeta, y_M=log(y_M+1), y_F=log(y_F+1),
                                                                 ID_M1=ID_M1, ID_M2=ID_M2, ID_F1=ID_F1, ID_F2=ID_F2, 
                                                                 times_M=times_M, times_F=times_F, X=X,
                                                                 Time=Time, status_2=data$Short$status),
                                                   init   = init_fun,
                                                   iter   = iter,
                                                   warmup = warmup,                 
                                                   chains = chains,
                                                   thin   = 1,
                                                   seed   = 123,
                                                   cores  = getOption("mc.cores",chains)) ))
    e.time <- Sys.time()
  }
  
  return( list(fit=fit, time=e.time-i.time) )
  
}


# Fitting joint models using the corrected two-stage specification
# Longitudinal submodels
fit_biexp <- function(data, biomarker="M-spike", iter=5000, warmup=1000, chains=3){

  # Biomarker information
  if(biomarker == "M-spike"){
    biom_data <- data$Long$M_Spike
    start <- data$Short$start_M
    stop <- data$Short$stop_M
  }{
    biom_data <- data$Long$FLC
    start <- data$Short$start_F
    stop <- data$Short$stop_F    
  }
  patientid <- data$Short$patientid
  
  # Model data
  N <- nrow(biom_data)
  n <- length(unique(biom_data$patientid))
  y <- biom_data$y
  tt <- biom_data$time
  start <- start[(patientid %in% unique(biom_data$patientid))]
  stop <- stop[(patientid %in% unique(biom_data$patientid))]
  ID <- (biom_data %>% group_by(patientid) %>% dplyr::mutate(ID = cur_group_id()))$ID
  
  # Setting initial values
  init_fun <- function(...){ 
    list(theta=c(0,0,0), sigma2=1, Omega=diag(rep(1,3)), bi=matrix(0,nrow=n,ncol=3))
  }
  
  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = TRUE)
  
  i.time <- Sys.time()
  fit <- suppressMessages(suppressWarnings( stan(model_code = BiExp, 
                                                 data   = list(N=N, n=n, y=log(y+1), ID=ID, times=tt, start=start, stop=stop),
                                                 init   = init_fun,
                                                 iter   = iter,
                                                 warmup = warmup,                 
                                                 chains = chains,
                                                 thin   = 1,
                                                 seed   = 123,
                                                 cores  = getOption("mc.cores",chains)) ))
  e.time <- Sys.time()
  
  return( list(fit=fit, time=e.time-i.time) )

}


# Survival submodels
fit_surv <- function(data=data, fit_LONG_M, fit_LONG_F, LoT, iter=1000, warmup=500, chains=3){
  
  vars <- c("sex","race_1","race_2","race_3",
            "ecog_1","ecog_2","ecog_3","iss_1","iss_2","iss_3",
            "age","albumin","beta2_microglobulin","creatinine",
            "hemoglobin","ldh","lymphocyte","neutrophil","platelet",
            "immunoglobulin_a","immunoglobulin_g","immunoglobulin_m")
  
  if(LoT==1){
    X <- data$Short[,vars]
  }else{
    X <- data$Short[,c(vars, "duration")]
  }
  
  # Model data
  N_M <- nrow(data$Long$M_Spike)
  N_F <- nrow(data$Long$FLC)
  n_M <- length(unique(data$Long$M_Spike$patientid))
  n_F <- length(unique(data$Long$FLC$patientid))
  y_M <- data$Long$M_Spike$y
  y_F <- data$Long$FLC$y
  times_M <- data$Long$M_Spike$time
  times_F <- data$Long$FLC$time
  n <- nrow(data$Short)
  nbeta <- ncol(X)
  ID_M1 <- (data$Long$M_Spike %>% group_by(patientid) %>% dplyr::mutate(ID = cur_group_id()))$ID
  ID_F1 <- (data$Long$FLC %>% group_by(patientid) %>% dplyr::mutate(ID = cur_group_id()))$ID
  ID_M2 <- data$Long$M_Spike$ID1
  ID_F2 <- data$Long$FLC$ID1
  Time <- data$Short$surv
  status_1 <- as.numeric((data$Short$status==0 & data$Short$lastLoT==0) | (data$Short$status==2 & data$Short$lastLoT==0))
  status_2 <- as.numeric(data$Short$status==1 | (data$Short$status==2 & data$Short$lastLoT==1))
  
  theta_M <- apply(extract(fit_LONG_M, "theta")$theta, 2, MAP_fc)
  theta_F <- apply(extract(fit_LONG_F, "theta")$theta, 2, MAP_fc)
  sigma2_M <- MAP_fc(extract(fit_LONG_M, "sigma2")$sigma2)
  sigma2_F <- MAP_fc(extract(fit_LONG_F, "sigma2")$sigma2)
  Omega_M <- apply(extract(fit_LONG_M, "Omega")$Omega, c(2,3), MAP_fc)
  Omega_F <- apply(extract(fit_LONG_F, "Omega")$Omega, c(2,3), MAP_fc)
  
  if(LoT != 4){
    # Setting initial values
    init_fun <- function(...){ 
      list(beta_1=rep(0,ncol(X)), beta_2=rep(0,ncol(X)),
           alphaB_1=c(0,0), alphaG_1=c(0,0), alphaD_1=c(0,0),
           alphaB_2=c(0,0), alphaG_2=c(0,0), alphaD_2=c(0,0),
           gamma_1=0, gamma_2=0, phi_1=1, phi_2=1,
           bi_M=matrix(0,nrow=n_M,ncol=3), bi_F=matrix(0,nrow=n_F,ncol=3))
    }
    
    options(mc.cores = parallel::detectCores())
    rstan_options(auto_write = TRUE)
    
    i.time <- Sys.time()
    fit <- suppressMessages(suppressWarnings( stan(model_code = WeibPHCompRisk, 
                                                   data   = list(N_M=N_M, N_F=N_F, n_M=n_M, n_F=n_F, n=n, 
                                                                 nbeta=nbeta, y_M=log(y_M+1), y_F=log(y_F+1),
                                                                 ID_M1=ID_M1, ID_M2=ID_M2, ID_F1=ID_F1, ID_F2=ID_F2, 
                                                                 times_M=times_M, times_F=times_F, X=X,
                                                                 Time=Time, status_1=status_1, status_2=status_2,
                                                                 theta_M=theta_M, theta_F=theta_F, sigma2_M=sigma2_M, 
                                                                 sigma2_F=sigma2_F, Omega_M=Omega_M, Omega_F=Omega_F),
                                                   init   = init_fun,
                                                   iter   = iter,
                                                   warmup = warmup,                 
                                                   chains = chains,
                                                   thin   = 1,
                                                   seed   = 123,
                                                   cores  = getOption("mc.cores",chains)) ))
    e.time <- Sys.time()
  }else{
    # Setting initial values
    init_fun <- function(...){ 
      list(beta_2=rep(0,ncol(X)),
           alphaB_2=c(0,0), alphaG_2=c(0,0), alphaD_2=c(0,0),
           gamma_2=0, phi_2=1,
           bi_M=matrix(0,nrow=n_M,ncol=3), bi_F=matrix(0,nrow=n_F,ncol=3))
    }
    
    options(mc.cores = parallel::detectCores())
    rstan_options(auto_write = TRUE)
    
    i.time <- Sys.time()
    fit <- suppressMessages(suppressWarnings( stan(model_code = WeibPHSurv, 
                                                   data   = list(N_M=N_M, N_F=N_F, n_M=n_M, n_F=n_F, n=n, 
                                                                 nbeta=nbeta, y_M=log(y_M+1), y_F=log(y_F+1),
                                                                 ID_M1=ID_M1, ID_M2=ID_M2, ID_F1=ID_F1, ID_F2=ID_F2, 
                                                                 times_M=times_M, times_F=times_F, X=X,
                                                                 Time=Time, status_2=data$Short$status,
                                                                 theta_M=theta_M, theta_F=theta_F, sigma_2_M=sigma2_M, 
                                                                 sigma2_F=sigma2_F, Omega_M=Omega_M, Omega_F=Omega_F),
                                                   init   = init_fun,
                                                   iter   = iter,
                                                   warmup = warmup,                 
                                                   chains = chains,
                                                   thin   = 1,
                                                   seed   = 123,
                                                   cores  = getOption("mc.cores",chains)) ))
    e.time <- Sys.time()
  }
  
  return( list(fit=fit, time=e.time-i.time) )
  
}


# Function for MAP
MAP_fc <- function(x){
    
  lim.inf <- min(x)-1
  lim.sup <- max(x)+1
  s <- density(x,from=lim.inf,to=lim.sup,bw=0.2)
  n <- length(s$y)
  v1 <- s$y[1:(n-2)]
  v2 <- s$y[2:(n-1)]
  v3 <- s$y[3:n]
  ix <- 1+which((v1<v2)&(v2>v3))
  out <- s$x[which(s$y==max(s$y))]
  
  return( out )
}


# Checking efficiency and convergence for joint specification and two-stage approaches
check_jm_ess_rhat <- function(fit_1, fit_2, fit_3, fit_4){
  
  pars <- c("theta_M","theta_F","sigma2_M","sigma2_F","Omega_M","Omega_F",
            "beta_1","beta_2","alphaB_1","alphaG_1","alphaD_1","alphaB_2",
            "alphaG_2","alphaD_2","gamma_1","gamma_2","phi_1","phi_2")
  
  # LoT 1
  ess_rhat_1 <- summary(fit_1, pars=pars)$summary
  # LoT 2
  ess_rhat_2 <- summary(fit_2, pars=pars)$summary
  # LoT 3
  ess_rhat_3 <- summary(fit_3, pars=pars)$summary
  # LoT 4
  ess_rhat_4 <- summary(fit_4, pars=c("theta_M","theta_F","sigma2_M","sigma2_F",
                                      "Omega_M","Omega_F","beta_2","alphaB_2",
                                      "alphaG_2","alphaD_2","gamma_2","phi_2"))$summary
  
  ess_all <- c(paste0(round(summary(ess_rhat_1[,9])[1]),"-",round(summary(ess_rhat_1[,9])[6])),
               paste0(round(summary(ess_rhat_2[,9])[1]),"-",round(summary(ess_rhat_2[,9])[6])),
               paste0(round(summary(ess_rhat_3[,9])[1]),"-",round(summary(ess_rhat_3[,9])[6])),
               paste0(round(summary(ess_rhat_4[,9])[1]),"-",round(summary(ess_rhat_4[,9])[6])))
  
  rhat_all <- c(paste0(round(summary(ess_rhat_1[,10])[1],3),"-",round(summary(ess_rhat_1[,10])[6],3)),
                paste0(round(summary(ess_rhat_2[,10])[1],3),"-",round(summary(ess_rhat_2[,10])[6],3)),
                paste0(round(summary(ess_rhat_3[,10])[1],3),"-",round(summary(ess_rhat_3[,10])[6],3)),
                paste0(round(summary(ess_rhat_4[,10])[1],3),"-",round(summary(ess_rhat_4[,10])[6],3)))
  
  out <- matrix(cbind(ess_all, rhat_all), ncol=2)
  colnames(out) <- c("ESS range", "R-hat range")
  rownames(out) <- c("LoT 1", "LoT 2", "LoT 3", "LoT 4")
  
  print( out )
}


check_long_ess_rhat <- function(fit_1_M, fit_1_F, fit_2_M, fit_2_F, fit_3_M, fit_3_F, fit_4_M, fit_4_F){

  pars <- c("theta","sigma2","Omega")
  
  # LoT 1
  ess_rhat_1_M <- summary(fit_1_M, pars=pars)$summary
  ess_rhat_1_F <- summary(fit_1_F, pars=pars)$summary
  # LoT 2
  ess_rhat_2_M <- summary(fit_2_M, pars=pars)$summary
  ess_rhat_2_F <- summary(fit_2_F, pars=pars)$summary
  # LoT 3
  ess_rhat_3_M <- summary(fit_3_M, pars=pars)$summary
  ess_rhat_3_F <- summary(fit_3_F, pars=pars)$summary
  # LoT 4
  ess_rhat_4_M <- summary(fit_4_M, pars=pars)$summary
  ess_rhat_4_F <- summary(fit_4_F, pars=pars)$summary
  
  ess_all <- c(paste0(round(summary(ess_rhat_1_M[,9])[1]),"-",round(summary(ess_rhat_1_M[,9])[6])),
               paste0(round(summary(ess_rhat_1_F[,9])[1]),"-",round(summary(ess_rhat_1_F[,9])[6])),
               paste0(round(summary(ess_rhat_2_M[,9])[1]),"-",round(summary(ess_rhat_2_M[,9])[6])),
               paste0(round(summary(ess_rhat_2_F[,9])[1]),"-",round(summary(ess_rhat_2_F[,9])[6])),
               paste0(round(summary(ess_rhat_3_M[,9])[1]),"-",round(summary(ess_rhat_3_M[,9])[6])),
               paste0(round(summary(ess_rhat_3_F[,9])[1]),"-",round(summary(ess_rhat_3_F[,9])[6])),
               paste0(round(summary(ess_rhat_4_M[,9])[1]),"-",round(summary(ess_rhat_4_M[,9])[6])),
               paste0(round(summary(ess_rhat_4_F[,9])[1]),"-",round(summary(ess_rhat_4_F[,9])[6])))
  
  rhat_all <- c(paste0(round(summary(ess_rhat_1_M[,10])[1],3),"-",round(summary(ess_rhat_1_M[,10])[6],3)),
                paste0(round(summary(ess_rhat_1_F[,10])[1],3),"-",round(summary(ess_rhat_1_F[,10])[6],3)),
                paste0(round(summary(ess_rhat_2_M[,10])[1],3),"-",round(summary(ess_rhat_2_M[,10])[6],3)),
                paste0(round(summary(ess_rhat_2_F[,10])[1],3),"-",round(summary(ess_rhat_2_F[,10])[6],3)),
                paste0(round(summary(ess_rhat_3_M[,10])[1],3),"-",round(summary(ess_rhat_3_M[,10])[6],3)),
                paste0(round(summary(ess_rhat_3_F[,10])[1],3),"-",round(summary(ess_rhat_3_F[,10])[6],3)),
                paste0(round(summary(ess_rhat_4_M[,10])[1],3),"-",round(summary(ess_rhat_4_M[,10])[6],3)),
                paste0(round(summary(ess_rhat_4_F[,10])[1],3),"-",round(summary(ess_rhat_4_F[,10])[6],3)))
  
  out <- matrix(cbind(ess_all, rhat_all),ncol=2)
  colnames(out) <- c("ESS range", "R-hat range")
  rownames(out) <- c("LoT 1: M-spike", "LoT 1: FLC", "LoT 2: M-spike", "LoT 2: FLC",
                     "LoT 3: M-spike", "LoT 3: FLC", "LoT 4: M-spike", "LoT 4: FLC")
  
  print( out )
  
}


check_surv_ess_rhat <- function(fit_1, fit_2, fit_3, fit_4){
    
  pars <- c("beta_1","beta_2","alphaB_1","alphaG_1","alphaD_1","alphaB_2",
            "alphaG_2","alphaD_2","gamma_1","gamma_2","phi_1","phi_2")
  
  # LoT 1
  ess_rhat_1 <- summary(fit_1, pars=pars)$summary
  # LoT 2
  ess_rhat_2 <- summary(fit_2, pars=pars)$summary
  # LoT 3
  ess_rhat_3 <- summary(fit_3, pars=pars)$summary
  # LoT 4
  ess_rhat_4 <- summary(fit_4, pars=c("beta_2","alphaB_2","alphaG_2",
                                      "alphaD_2","gamma_2","phi_2"))$summary

  ess_all <- c(paste0(round(summary(ess_rhat_1[,9])[1]),"-",round(summary(ess_rhat_1[,9])[6])),
               paste0(round(summary(ess_rhat_2[,9])[1]),"-",round(summary(ess_rhat_2[,9])[6])),
               paste0(round(summary(ess_rhat_3[,9])[1]),"-",round(summary(ess_rhat_3[,9])[6])),
               paste0(round(summary(ess_rhat_4[,9])[1]),"-",round(summary(ess_rhat_4[,9])[6])))
  
  rhat_all <- c(paste0(round(summary(ess_rhat_1[,10])[1],3),"-",round(summary(ess_rhat_1[,10])[6],3)),
                paste0(round(summary(ess_rhat_2[,10])[1],3),"-",round(summary(ess_rhat_2[,10])[6],3)),
                paste0(round(summary(ess_rhat_3[,10])[1],3),"-",round(summary(ess_rhat_3[,10])[6],3)),
                paste0(round(summary(ess_rhat_4[,10])[1],3),"-",round(summary(ess_rhat_4[,10])[6],3)))
    
  out <- matrix(cbind(ess_all, rhat_all), ncol=2)
  colnames(out) <- c("ESS range", "R-hat range")
  rownames(out) <- c("LoT 1", "LoT 2", "LoT 3", "LoT 4")
    
  print( out )
}


# Posterior summary for joint specification and two-stage approaches
posterior_summary_jm <- function(fit, type="long", LoT){
  
  out <- NULL
  if(type=="long"){
    out <- summary(fit, pars=c("theta_M","sigma2_M","Omega_M",
                               "theta_F","sigma2_F","Omega_F"))$summary[,c(1,4,8)]
  }else{
    if(LoT!=4){
      out <- summary(fit, pars=c("beta_1","alphaB_1","alphaG_1","alphaD_1","gamma_1","phi_1",
                                 "beta_2","alphaB_2","alphaG_2","alphaD_2","gamma_2","phi_2"))$summary[,c(1,4,8)]
    }else{
      out <- summary(fit, pars=c("beta_2","alphaB_2","alphaG_2",
                                 "alphaD_2","gamma_2","phi_2"))$summary[,c(1,4,8)]
    }
  }
  
  return( out )
  
}


posterior_summary_ts <- function(fit, type="long", LoT){
  
  out <- NULL
  if(type=="long"){
    out <- summary(fit, pars=c("theta","sigma2","Omega"))$summary[,c(1,4,8)]
  }else{
    if(LoT!=4){
      out <- summary(fit, pars=c("beta_1","beta_2","alphaB_1","alphaG_1",
                                 "alphaD_1","alphaB_2","alphaG_2","alphaD_2",
                                 "gamma_1","gamma_2","phi_1","phi_2"))$summary[,c(1,4,8)]
    }else{
      out <- summary(fit, pars=c("beta_2","alphaB_2","alphaG_2",
                                 "alphaD_2","gamma_2","phi_2"))$summary[,c(1,4,8)]
    }
  }
  
  return( out )
  
}

