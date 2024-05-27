# Fitting joint models using the joint specification
fit_jm <- function(data=data, LoT, iter=5000, warmup=1000, chains=3){
  
  vars <- c("sex","race_1","race_2","race_3",
            "ecog_1","ecog_2","ecog_3","iss_1","iss_2","iss_3",
            "age","albumin","beta2_microglobulin","creatinine",
            "hemoglobin","ldh","lymphocyte","neutrophil","platelet",
            "immunoglobulin_a","immunoglobulin_g","immunoglobulin_m")
  
  if(LoT == 1){
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
  
  options(mc.cores=parallel::detectCores())
  rstan_options(auto_write=TRUE)
  
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
    
    options(mc.cores=parallel::detectCores())
    rstan_options(auto_write=TRUE)
    
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
    
    options(mc.cores=parallel::detectCores())
    rstan_options(auto_write=TRUE)
    
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
  if(type == "long"){
    out <- summary(fit, pars=c("theta","sigma2","Omega"))$summary[,c(1,4,8)]
  }else{
    if(LoT != 4){
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


# Fitted average longitudinal trajectory
avg_trajectory_fc <- function(fit1, fit2, fit3, fit4, fit5=NULL, fit6=NULL, fit7=NULL, fit8=NULL, approach="JE"){
  
  # M-spike
  if(approach == "JE"){
    B10_M <- exp(extract(fit1, "theta_M")$theta_M[,1])
    G10_M <- exp(extract(fit1, "theta_M")$theta_M[,2])
    D10_M <- exp(extract(fit1, "theta_M")$theta_M[,3])
    
    B20_M <- exp(extract(fit2, "theta_M")$theta_M[,1])
    G20_M <- exp(extract(fit2, "theta_M")$theta_M[,2])
    D20_M <- exp(extract(fit2, "theta_M")$theta_M[,3])
    
    B30_M <- exp(extract(fit3, "theta_M")$theta_M[,1])
    G30_M <- exp(extract(fit3, "theta_M")$theta_M[,2])
    D30_M <- exp(extract(fit3, "theta_M")$theta_M[,3])
    
    B40_M <- exp(extract(fit4, "theta_M")$theta_M[,1])
    G40_M <- exp(extract(fit4, "theta_M")$theta_M[,2])
    D40_M <- exp(extract(fit4, "theta_M")$theta_M[,3])
  }else{
    B10_M <- exp(extract(fit1, "theta")$theta[,1])
    G10_M <- exp(extract(fit1, "theta")$theta[,2])
    D10_M <- exp(extract(fit1, "theta")$theta[,3])
    
    B20_M <- exp(extract(fit2, "theta")$theta[,1])
    G20_M <- exp(extract(fit2, "theta")$theta[,2])
    D20_M <- exp(extract(fit2, "theta")$theta[,3])
    
    B30_M <- exp(extract(fit3, "theta")$theta[,1])
    G30_M <- exp(extract(fit3, "theta")$theta[,2])
    D30_M <- exp(extract(fit3, "theta")$theta[,3])
    
    B40_M <- exp(extract(fit4, "theta")$theta[,1])
    G40_M <- exp(extract(fit4, "theta")$theta[,2])
    D40_M <- exp(extract(fit4, "theta")$theta[,3])
  }
  
  # FLC
  if(approach == "JE"){
    B10_F <- exp(extract(fit1, "theta_F")$theta_F[,1])
    G10_F <- exp(extract(fit1, "theta_F")$theta_F[,2])
    D10_F <- exp(extract(fit1, "theta_F")$theta_F[,3])
    
    B20_F <- exp(extract(fit2, "theta_F")$theta_F[,1])
    G20_F <- exp(extract(fit2, "theta_F")$theta_F[,2])
    D20_F <- exp(extract(fit2, "theta_F")$theta_F[,3])
    
    B30_F <- exp(extract(fit3, "theta_F")$theta_F[,1])
    G30_F <- exp(extract(fit3, "theta_F")$theta_F[,2])
    D30_F <- exp(extract(fit3, "theta_F")$theta_F[,3])
    
    B40_F <- exp(extract(fit4, "theta_F")$theta_F[,1])
    G40_F <- exp(extract(fit4, "theta_F")$theta_F[,2])
    D40_F <- exp(extract(fit4, "theta_F")$theta_F[,3])
  }else{
    B10_F <- exp(extract(fit5, "theta")$theta[,1])
    G10_F <- exp(extract(fit5, "theta")$theta[,2])
    D10_F <- exp(extract(fit5, "theta")$theta[,3])
    
    B20_F <- exp(extract(fit6, "theta")$theta[,1])
    G20_F <- exp(extract(fit6, "theta")$theta[,2])
    D20_F <- exp(extract(fit6, "theta")$theta[,3])
    
    B30_F <- exp(extract(fit7, "theta")$theta[,1])
    G30_F <- exp(extract(fit7, "theta")$theta[,2])
    D30_F <- exp(extract(fit7, "theta")$theta[,3])
    
    B40_F <- exp(extract(fit8, "theta")$theta[,1])
    G40_F <- exp(extract(fit8, "theta")$theta[,2])
    D40_F <- exp(extract(fit8, "theta")$theta[,3])
  }
  
  tt <- seq(0, 7, len=100)
  long11 <- long12 <- long13 <- long14 <- rep(NA, length(tt))
  long21 <- long22 <- long23 <- long24 <- rep(NA, length(tt))
  long31 <- long32 <- long33 <- long34 <- rep(NA, length(tt))
  long11L <- long12L <- long13L <- long14L <- rep(NA, length(tt))
  long21L <- long22L <- long23L <- long24L <- rep(NA, length(tt))
  long31L <- long32L <- long33L <- long34L <- rep(NA, length(tt))
  long11U <- long12U <- long13U <- long14U <- rep(NA, length(tt))
  long21U <- long22U <- long23U <- long24U <- rep(NA, length(tt))
  long31U <- long32U <- long33U <- long34U <- rep(NA, length(tt))
  
  for(l in 1:length(tt)){
    # M-spike
    longaux11 <- B10_M * (exp(G10_M * tt[l]) + exp(-D10_M * tt[l]) - 1)
    longaux12 <- B20_M * (exp(G20_M * tt[l]) + exp(-D20_M * tt[l]) - 1)
    longaux13 <- B30_M * (exp(G30_M * tt[l]) + exp(-D30_M * tt[l]) - 1)
    longaux14 <- B40_M * (exp(G40_M * tt[l]) + exp(-D40_M * tt[l]) - 1)
    
    long11[l] <- mean(longaux11)
    long12[l] <- mean(longaux12)
    long13[l] <- mean(longaux13)
    long14[l] <- mean(longaux14)
    
    long11L[l] <- quantile(longaux11, probs=0.025)
    long12L[l] <- quantile(longaux12, probs=0.025)
    long13L[l] <- quantile(longaux13, probs=0.025)
    long14L[l] <- quantile(longaux14, probs=0.025)
    
    long11U[l] <- quantile(longaux11, probs=0.975)
    long12U[l] <- quantile(longaux12, probs=0.975)
    long13U[l] <- quantile(longaux13, probs=0.975)
    long14U[l] <- quantile(longaux14, probs=0.975)
    
    # FLC
    longaux21 <- B10_F * (exp(G10_F * tt[l]) + exp(-D10_F * tt[l]) - 1)
    longaux22 <- B20_F * (exp(G20_F * tt[l]) + exp(-D20_F * tt[l]) - 1)
    longaux23 <- B30_F * (exp(G30_F * tt[l]) + exp(-D30_F * tt[l]) - 1)
    longaux24 <- B40_F * (exp(G40_F * tt[l]) + exp(-D40_F * tt[l]) - 1)
    
    long21[l] <- mean(longaux21)
    long22[l] <- mean(longaux22)
    long23[l] <- mean(longaux23)
    long24[l] <- mean(longaux24)
    
    long21L[l] <- quantile(longaux21, probs=0.025)
    long22L[l] <- quantile(longaux22, probs=0.025)
    long23L[l] <- quantile(longaux23, probs=0.025)
    long24L[l] <- quantile(longaux24, probs=0.025)
    
    long21U[l] <- quantile(longaux21, probs=0.975)
    long22U[l] <- quantile(longaux22, probs=0.975)
    long23U[l] <- quantile(longaux23, probs=0.975)
    long24U[l] <- quantile(longaux24, probs=0.975)
  }
  
  out <- data.frame(time=rep(tt,8),
                    long=c(long11,long12,long13,long14,long21,long22,long23,long24),
                    longL=c(long11L,long12L,long13L,long14L,long21L,long22L,long23L,long24L),
                    longU=c(long11U,long12U,long13U,long14U,long21U,long22U,long23U,long24U),
                    biom=factor(c(rep("M-spike",4*length(tt)),rep("FLC",4*length(tt)))),
                    LoT=factor(rep(rep(c("LoT 1","LoT 2","LoT 3","LoT 4"),each=length(tt)),2)))
  out$biom <- factor(out$biom, levels=c("M-spike","FLC"))
  
  return( out )
}


# Updating random effects (test set data and dynamic predictions)
update_bi <- function(data, fit1, fit2=NULL, fit3=NULL, approach="JE", LoT, iter=100){

  # Data
  vars <- c("sex","race_1","race_2","race_3",
            "ecog_1","ecog_2","ecog_3","iss_1","iss_2","iss_3",
            "age","albumin","beta2_microglobulin","creatinine",
            "hemoglobin","ldh","lymphocyte","neutrophil","platelet",
            "immunoglobulin_a","immunoglobulin_g","immunoglobulin_m")
  
  if(LoT == 1){
    X <- data$Short[,vars]
  }else{
    X <- data$Short[,c(vars, "duration")]
  }
  
  y_M <- log(data$Long$M_Spike$y+1)
  y_F <- log(data$Long$FLC$y+1)
  times_M <- data$Long$M_Spike$time
  times_F <- data$Long$FLC$time
  Time <- data$Short$surv
  status_1 <- as.numeric((data$Short$status==0 & data$Short$lastLoT==0) | (data$Short$status==2 & data$Short$lastLoT==0))
  if(LoT != 4){
    status_2 <- as.numeric(data$Short$status==1 | (data$Short$status==2 & data$Short$lastLoT==1))
  }else{
    status_2 <-  data$Short$status
  }
  
  # Parameters
  if(approach == "JE"){
    theta_M <- extract(fit1, "theta_M")$theta_M 
    theta_F <- extract(fit1, "theta_F")$theta_F
    sigma2_M <- extract(fit1, "sigma2_M")$sigma2_M
    sigma2_F <- extract(fit1, "sigma2_F")$sigma2_F
    Omega_M <- extract(fit1, "Omega_M")$Omega_M 
    Omega_F <- extract(fit1, "Omega_F")$Omega_F
    beta_2 <- extract(fit1, "beta_2")$beta_2
    alphaB_2 <- extract(fit1, "alphaB_2")$alphaB_2
    alphaG_2 <- extract(fit1, "alphaG_2")$alphaG_2
    alphaD_2 <- extract(fit1, "alphaD_2")$alphaD_2
    gamma_2 <- extract(fit1, "gamma_2")$gamma_2
    phi_2 <- extract(fit1, "phi_2")$phi_2
    if(LoT != 4){
      beta_1 <- extract(fit1, "beta_1")$beta_1
      alphaB_1 <- extract(fit1, "alphaB_1")$alphaB_1
      alphaG_1 <- extract(fit1, "alphaG_1")$alphaG_1
      alphaD_1 <- extract(fit1, "alphaD_1")$alphaD_1
      gamma_1 <- extract(fit1, "gamma_1")$gamma_1
      phi_1 <- extract(fit1, "phi_1")$phi_1
    }
  }else{
    theta_M <- extract(fit1, "theta")$theta 
    theta_F <- extract(fit2, "theta")$theta
    sigma2_M <- extract(fit1, "sigma2")$sigma2
    sigma2_F <- extract(fit2, "sigma2")$sigma2
    Omega_M <- extract(fit1, "Omega")$Omega 
    Omega_F <- extract(fit2, "Omega")$Omega
    beta_2 <- extract(fit3, "beta_2")$beta_2
    alphaB_2 <- extract(fit3, "alphaB_2")$alphaB_2
    alphaG_2 <- extract(fit3, "alphaG_2")$alphaG_2
    alphaD_2 <- extract(fit3, "alphaD_2")$alphaD_2
    gamma_2 <- extract(fit3, "gamma_2")$gamma_2
    phi_2 <- extract(fit3, "phi_2")$phi_2
    if(LoT != 4){
      beta_1 <- extract(fit3, "beta_1")$beta_1
      alphaB_1 <- extract(fit3, "alphaB_1")$alphaB_1
      alphaG_1 <- extract(fit3, "alphaG_1")$alphaG_1
      alphaD_1 <- extract(fit3, "alphaD_1")$alphaD_1
      gamma_1 <- extract(fit3, "gamma_1")$gamma_1
      phi_1 <- extract(fit3, "phi_1")$phi_1
    }
  }
  
  if(approach == "JE"){
    prop_theta_M <- mvrnorm(iter, mu=apply(theta_M, 2, MAP_fc), Sigma=diag(apply(theta_M, 2, var))) 
    prop_theta_F <- mvrnorm(iter, mu=apply(theta_F, 2, MAP_fc), Sigma=diag(apply(theta_F, 2, var))) 
    prop_sigma_M <- sqrt(exp(rnorm(iter, mean=MAP_fc(log(sigma2_M)), sd=sd(log(sigma2_M)))))
    prop_sigma_F <- sqrt(exp(rnorm(iter, mean=MAP_fc(log(sigma2_F)), sd=sd(log(sigma2_F)))))
    prop_Omega_M <- exp(mvrnorm(iter, mu=suppressWarnings(diag(log(apply(Omega_M, c(2,3), MAP_fc)))), Sigma=diag(suppressWarnings(diag(apply(log(Omega_M), c(2,3), var))))))
    prop_Omega_F <- exp(mvrnorm(iter, mu=suppressWarnings(diag(log(apply(Omega_F, c(2,3), MAP_fc)))), Sigma=diag(suppressWarnings(diag(apply(log(Omega_F), c(2,3), var))))))
    prop_beta_2 <- mvrnorm(iter, mu=apply(beta_2, 2, MAP_fc), Sigma=diag(apply(beta_2, 2, var)))
    prop_alphaB_2 <- mvrnorm(iter, mu=apply(alphaB_2, 2, MAP_fc), Sigma=diag(apply(alphaB_2, 2, var)))
    prop_alphaG_2 <- mvrnorm(iter, mu=apply(alphaG_2, 2, MAP_fc), Sigma=diag(apply(alphaG_2, 2, var)))
    prop_alphaD_2 <- mvrnorm(iter, mu=apply(alphaD_2, 2, MAP_fc), Sigma=diag(apply(alphaD_2, 2, var)))
    prop_gamma_2 <- rnorm(iter, mean=MAP_fc(gamma_2), sd=sd(gamma_2))
    prop_phi_2 <- exp(rnorm(iter, mean=MAP_fc(log(phi_2)), sd=sd(log(phi_2))))
    if(LoT != 4){
      prop_beta_1 <- mvrnorm(iter, mu=apply(beta_1, 2, MAP_fc), Sigma=diag(apply(beta_1, 2, var)))
      prop_alphaB_1 <- mvrnorm(iter, mu=apply(alphaB_1, 2, MAP_fc), Sigma=diag(apply(alphaB_1, 2, var)))
      prop_alphaG_1 <- mvrnorm(iter, mu=apply(alphaG_1, 2, MAP_fc), Sigma=diag(apply(alphaG_1, 2, var)))
      prop_alphaD_1 <- mvrnorm(iter, mu=apply(alphaD_1, 2, MAP_fc), Sigma=diag(apply(alphaD_1, 2, var)))
      prop_gamma_1 <- rnorm(iter, mean=MAP_fc(gamma_1), sd=sd(gamma_1))
      prop_phi_1 <- exp(rnorm(iter, mean=MAP_fc(log(phi_1)), sd=sd(log(phi_1))))
    }
  }else{
    prop_theta_M <- matrix(apply(theta_M, 2, MAP_fc), nrow=iter, ncol=length(apply(theta_M, 2, MAP_fc)), byrow=TRUE)
    prop_theta_F <- matrix(apply(theta_F, 2, MAP_fc), nrow=iter, ncol=length(apply(theta_F, 2, MAP_fc)), byrow=TRUE)
    prop_sigma_M <- sqrt(rep(MAP_fc(sigma2_M), iter))
    prop_sigma_F <- sqrt(rep(MAP_fc(sigma2_F), iter))
    prop_Omega_M <- matrix(diag(apply(Omega_M, c(2,3), MAP_fc)), nrow=iter, ncol=length(diag(apply(Omega_M, c(2,3), MAP_fc))), byrow=TRUE)
    prop_Omega_F <- matrix(diag(apply(Omega_F, c(2,3), MAP_fc)), nrow=iter, ncol=length(diag(apply(Omega_F, c(2,3), MAP_fc))), byrow=TRUE)
    prop_beta_2 <- matrix(apply(beta_2, 2, MAP_fc), nrow=iter, ncol=length(apply(beta_2, 2, MAP_fc)), byrow=TRUE)
    prop_alphaB_2 <- matrix(apply(alphaB_2, 2, MAP_fc), nrow=iter, ncol=length(apply(alphaB_2, 2, MAP_fc)), byrow=TRUE)
    prop_alphaG_2 <- matrix(apply(alphaG_2, 2, MAP_fc), nrow=iter, ncol=length(apply(alphaG_2, 2, MAP_fc)), byrow=TRUE)
    prop_alphaD_2 <- matrix(apply(alphaD_2, 2, MAP_fc), nrow=iter, ncol=length(apply(alphaD_2, 2, MAP_fc)), byrow=TRUE)
    prop_gamma_2 <- rep(MAP_fc(gamma_2), iter)
    prop_phi_2 <- rep(MAP_fc(phi_2), iter)
    if(LoT != 4){
      prop_beta_1 <- matrix(apply(beta_1, 2, MAP_fc), nrow=iter, ncol=length(apply(beta_1, 2, MAP_fc)), byrow=TRUE)
      prop_alphaB_1 <- matrix(apply(alphaB_1, 2, MAP_fc), nrow=iter, ncol=length(apply(alphaB_1, 2, MAP_fc)), byrow=TRUE)
      prop_alphaG_1 <- matrix(apply(alphaG_1, 2, MAP_fc), nrow=iter, ncol=length(apply(alphaG_1, 2, MAP_fc)), byrow=TRUE)
      prop_alphaD_1 <- matrix(apply(alphaD_1, 2, MAP_fc), nrow=iter, ncol=length(apply(alphaD_1, 2, MAP_fc)), byrow=TRUE)
      prop_gamma_1 <- rep(MAP_fc(gamma_1), iter)
      prop_phi_1 <- rep(MAP_fc(phi_1), iter)
    }
  }
  
  MAP_Omega <- diag(c(diag(apply(Omega_M, c(2,3), MAP_fc)),diag(apply(Omega_F, c(2,3), MAP_fc))))
  bi <- matrix(0, nrow=iter, ncol=6)
  bi_prev <- rep(0, 6)
  
  for(j in 1:iter){
    
    p.log <- function(b){
      
      # M-spike
      prop_logBi_M <- as.numeric(prop_theta_M[j,1] + b[1])
      prop_Gi_M <- as.numeric(exp(prop_theta_M[j,2] + b[2]))
      prop_Di_M <- as.numeric(exp(prop_theta_M[j,3] + b[3]))
      
      mu_prop_M <- prop_logBi_M + log( exp(prop_Gi_M * times_M) + exp(-prop_Di_M * times_M) - 1 )
      log_like_prop_M <- sum(dnorm(y_M, mean=mu_prop_M, sd=prop_sigma_M[j], log=TRUE))
      log_re_prop_M <- dmvnorm(b[1:3], mu=c(0,0,0), Sigma=diag(prop_Omega_M[j,]), log=TRUE)
      
      # FLC
      prop_logBi_F <- as.numeric(prop_theta_F[j,1] + b[4])
      prop_Gi_F <- as.numeric(exp(prop_theta_F[j,2] + b[5]))
      prop_Di_F <- as.numeric(exp(prop_theta_F[j,3] + b[6]))
      
      mu_prop_F <- prop_logBi_F + log( exp(prop_Gi_F * times_F) + exp(-prop_Di_F * times_F) - 1 )
      log_like_prop_F <- sum(dnorm(y_F, mean=mu_prop_F, sd=prop_sigma_F[j], log=TRUE))
      log_re_prop_F <- dmvnorm(b[4:6], mu=c(0,0,0), Sigma=diag(prop_Omega_F[j,]), log=TRUE)
      
      # Survival
      Xbeta_2 <- prop_gamma_2[j] + sum(X*prop_beta_2[j,])
      lBialpha_2 <- sum(c(prop_logBi_M,prop_logBi_F)*prop_alphaB_2[j,])
      lGialpha_2 <- sum(log(c(prop_Gi_M,prop_Gi_F))*prop_alphaG_2[j,])
      lDialpha_2 <- sum(log(c(prop_Di_M,prop_Di_F))*prop_alphaD_2[j,])
      
      loghaz_2 <- log(prop_phi_2[j]) + (prop_phi_2[j]-1)*log(Time) + Xbeta_2 + lBialpha_2 + lGialpha_2 + lDialpha_2
      cumHaz_2 <- Time^(prop_phi_2[j]) * exp( prop_gamma_2[j] + Xbeta_2 + lBialpha_2 + lGialpha_2 + lDialpha_2 )
      
      if(LoT != 4){
        Xbeta_1 <- prop_gamma_1[j] + sum(X*prop_beta_1[j,])
        lBialpha_1 <- sum(c(prop_logBi_M,prop_logBi_F)*prop_alphaB_1[j,])
        lGialpha_1 <- sum(log(c(prop_Gi_M,prop_Gi_F))*prop_alphaG_1[j,])
        lDialpha_1 <- sum(log(c(prop_Di_M,prop_Di_F))*prop_alphaD_1[j,])
        
        loghaz_1 <- log(prop_phi_1[j]) + (prop_phi_1[j]-1)*log(Time) + Xbeta_1 + lBialpha_1 + lGialpha_1 + lDialpha_1
        cumHaz_1 <- Time^(prop_phi_1[j]) * exp( prop_gamma_1[j] + Xbeta_1 + lBialpha_1 + lGialpha_1 + lDialpha_1 )
        
        log_like_surv <- status_1*loghaz_1 + status_2*loghaz_2 - cumHaz_1 - cumHaz_2
      }else{
        log_like_surv <- status_2*loghaz_2 - cumHaz_2
      }

      # Joint likelihood
      return( log_like_prop_M + log_re_prop_M + log_like_prop_F + log_re_prop_F + log_like_surv )
    }

    # Adaptive MCMC
    invisible(capture.output( bb <- MCMC(p.log, n=100, init=bi_prev, scale=MAP_Omega,
                                         adapt=TRUE, acc.rate=0.234, showProgressBar=FALSE) ))
    bi[j,] <- apply(bb$samples[-c(1:50),], 2, median)
    bi_prev <- bi[j,] 
  }
  
  if(LoT == 4){
    prop_beta_1 <- prop_alphaB_1 <- prop_alphaG_1 <- prop_alphaD_1 <- prop_gamma_1 <- prop_phi_1 <- NA
  }
  
  out <- data.frame(theta_M=prop_theta_M, sigma_M=prop_sigma_M, Omega_M=prop_Omega_M, bi_M=bi[,1:3], 
                    theta_F=prop_theta_F, sigma_M=prop_sigma_F, Omega_F=prop_Omega_F, bi_F=bi[,4:6],
                    beta_1=prop_beta_1, alphaB_1=prop_alphaB_1, alphaG_1=prop_alphaG_1,
                    alphaD_1=prop_alphaD_1, gamma_1=prop_gamma_1, phi_1=prop_phi_1,
                    beta_2=prop_beta_2, alphaB_2=prop_alphaB_2, alphaG_2=prop_alphaG_2,
                    alphaD_2=prop_alphaD_2, gamma_2=prop_gamma_2, phi_2=prop_phi_2)
  
  return( out )
  
}

update_RE <- function(data, fit_JM, fit_LONG_M, fit_LONG_F, fit_CR_Surv, LoT, progress=FALSE){
  
  uniqueID <- unique(data$Short$patientid)
  len_uniqueID <- length(uniqueID)
  update_JM <- update_TS <- list()
  
  for(i in 1:len_uniqueID){
    if(progress){ print(paste0(i,"/",len_uniqueID)) }
    dd <- data
    dd$Short <- dd1$Short[dd$Short$patientid == uniqueID[i],]
    dd$Long$M_Spike <- dd$Long$M_Spike[dd$Long$M_Spike$patientid == uniqueID[i],]
    dd$Long$FLC <- dd$Long$FLC[dd$Long$FLC$patientid == uniqueID[i],]
    update_JM[[i]] <- update_bi(data=dd, fit1=fit_JM$fit, fit2=NULL, fit3=NULL, approach="JE", LoT=LoT, iter=100)
    update_TS[[i]] <- update_bi(data=dd, fit1=fit_LONG_M$fit, fit2=fit_LONG_F$fit, fit3=fit_CR_Surv$fit, approach="TS", LoT=LoT, iter=100)
  }
  
  return( list(update_JM=update_JM, update_TS=update_TS) )
}
