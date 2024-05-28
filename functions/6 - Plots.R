# FITTED AVERAGE LONGITUDINAL TRAJECTORY
plot_avg_trajectory <- function(fit_JM1, fit_JM2, fit_JM3, fit_JM4,
                                fit_LONG1_M, fit_LONG2_M, fit_LONG3_M, fit_LONG4_M,
                                fit_LONG1_F, fit_LONG2_F, fit_LONG3_F, fit_LONG4_F){
  
  dta_JM <- avg_trajectory_fc(fit1=fit_JM1$fit, fit2=fit_JM2$fit, fit3=fit_JM3$fit, fit4=fit_JM4$fit, approach="JE")
  dta_TS <- avg_trajectory_fc(fit1=fit_LONG1_M$fit, fit2=fit_LONG2_M$fit, fit3=fit_LONG3_M$fit, fit4=fit_LONG4_M$fit, 
                           fit5=fit_LONG1_F$fit, fit6=fit_LONG2_F$fit, fit7=fit_LONG3_F$fit, fit8=fit_LONG4_F$fit, approach="TS")
  
  cols <- c("JE"="black", "TS"="gray")
  ppp <- ggplot(data=dta_TS, aes(x=time, y=long)) + geom_line(linewidth=0.2, linetype="dashed") + 
    geom_ribbon(aes(ymin=longL, ymax=longU, fill="TS"), alpha=0.3) + 
    xlab("Time (in years)") + ylab("Biomarkers") + theme_bw() +
    theme(legend.position="top", panel.grid.major=element_blank(), 
          panel.grid.minor=element_blank(), panel.background=element_blank()) +
    facet_grid(vars(biom),vars(LoT),scales="free") + labs(fill="Approach") +
    scale_x_continuous(breaks=seq(0,7,1)) + scale_fill_manual(values=cols)
  
  ppp + geom_line(data=dta_JM, aes(x=time, y=long), linewidth=0.2) + 
    geom_ribbon(data=dta_JM, aes(ymin=longL, ymax=longU, fill="JE"), alpha=0.2)
  
  # ggsave("fig.png", units="in", width=11, height=8, dpi=300)
}


# INDIVIDUAL WEIGHTED RESIDUALS (IWRES)
plot_iwres <- function(data1, data2, data3, data4,
                       fit_JM1, RE_JM1, fit_JM2, RE_JM2, fit_JM3, RE_JM3, fit_JM4, RE_JM4,
                       fit_LONG1_M, fit_LONG1_F, RE_TS1, fit_LONG2_M, fit_LONG2_F, RE_TS2,
                       fit_LONG3_M, fit_LONG3_F, RE_TS3, fit_LONG4_M, fit_LONG4_F, RE_TS4){

  # LoT 1
  pm1 <- posterior_mean_fc(fit_JM1$fit, RE_JM1, fit_LONG1_M$fit, fit_LONG1_F$fit, RE_TS1)
  # M-spike
  iwres_1_M_JM <- iwres_fc(data=data1$Long$M_Spike, start=data1$Short$start_M, stop=data1$Short$stop_M, 
                           lB=pm1$JM$theta_M[1] + pm1$JM$b1_M, lG=pm1$JM$theta_M[2] + pm1$JM$b2_M, lD=pm1$JM$theta_M[3] + pm1$JM$b3_M, 
                           sigma=mean(sqrt(extract(fit_JM1$fit, "sigma2_M")$sigma2_M)))
  iwres_1_M_TS <- iwres_fc(data=data1$Long$M_Spike, start=data1$Short$start_M, stop=data1$Short$stop_M, 
                           lB=pm1$TS$theta_M[1] + pm1$TS$b1_M, lG=pm1$TS$theta_M[2] + pm1$TS$b2_M, lD=pm1$TS$theta_M[3] + pm1$TS$b3_M, 
                           sigma=mean(sqrt(extract(fit_LONG1_M$fit, "sigma2")$sigma2)))
  # FLC
  iwres_1_F_JM <- iwres_fc(data=data1$Long$FLC, start=data1$Short$start_F, stop=data1$Short$stop_F, 
                           lB=pm1$JM$theta_F[1] + pm1$JM$b1_F, lG=pm1$JM$theta_F[2] + pm1$JM$b2_F, lD=pm1$JM$theta_F[3] + pm1$JM$b3_F, 
                           sigma=mean(sqrt(extract(fit_JM1$fit, "sigma2_F")$sigma2_F)))
  iwres_1_F_TS <- iwres_fc(data=data1$Long$FLC, start=data1$Short$start_F, stop=data1$Short$stop_F, 
                           lB=pm1$TS$theta_F[1] + pm1$TS$b1_F, lG=pm1$TS$theta_F[2] + pm1$TS$b2_F, lD=pm1$TS$theta_F[3] + pm1$TS$b3_F, 
                           sigma=mean(sqrt(extract(fit_LONG1_F$fit, "sigma2")$sigma2)))
  
  # LoT 2
  pm2 <- posterior_mean_fc(fit_JM2$fit, RE_JM2, fit_LONG2_M$fit, fit_LONG2_F$fit, RE_TS2)
  # M-spike
  iwres_2_M_JM <- iwres_fc(data=data2$Long$M_Spike, start=data2$Short$start_M, stop=data2$Short$stop_M, 
                           lB=pm2$JM$theta_M[1] + pm2$JM$b1_M, lG=pm2$JM$theta_M[2] + pm2$JM$b2_M, lD=pm2$JM$theta_M[3] + pm2$JM$b3_M, 
                           sigma=mean(sqrt(extract(fit_JM2$fit, "sigma2_M")$sigma2_M)))
  iwres_2_M_TS <- iwres_fc(data=data2$Long$M_Spike, start=data2$Short$start_M, stop=data2$Short$stop_M, 
                           lB=pm2$TS$theta_M[1] + pm2$TS$b1_M, lG=pm2$TS$theta_M[2] + pm2$TS$b2_M, lD=pm2$TS$theta_M[3] + pm2$TS$b3_M, 
                           sigma=mean(sqrt(extract(fit_LONG2_M$fit, "sigma2")$sigma2)))
  # FLC
  iwres_2_F_JM <- iwres_fc(data=data2$Long$FLC, start=data2$Short$start_F, stop=data2$Short$stop_F, 
                           lB=pm2$JM$theta_F[1] + pm2$JM$b1_F, lG=pm2$JM$theta_F[2] + pm2$JM$b2_F, lD=pm2$JM$theta_F[3] + pm2$JM$b3_F, 
                           sigma=mean(sqrt(extract(fit_JM2$fit, "sigma2_F")$sigma2_F)))
  iwres_2_F_TS <- iwres_fc(data=data2$Long$FLC, start=data2$Short$start_F, stop=data2$Short$stop_F, 
                           lB=pm2$TS$theta_F[1] + pm2$TS$b1_F, lG=pm2$TS$theta_F[2] + pm2$TS$b2_F, lD=pm2$TS$theta_F[3] + pm2$TS$b3_F, 
                           sigma=mean(sqrt(extract(fit_LONG2_F$fit, "sigma2")$sigma2)))

  # LoT 3
  pm3 <- posterior_mean_fc(fit_JM3$fit, RE_JM3, fit_LONG3_M$fit, fit_LONG3_F$fit, RE_TS3)
  # M-spike
  iwres_3_M_JM <- iwres_fc(data=data3$Long$M_Spike, start=data3$Short$start_M, stop=data3$Short$stop_M, 
                           lB=pm3$JM$theta_M[1] + pm3$JM$b1_M, lG=pm3$JM$theta_M[2] + pm3$JM$b2_M, lD=pm3$JM$theta_M[3] + pm3$JM$b3_M, 
                           sigma=mean(sqrt(extract(fit_JM3$fit, "sigma2_M")$sigma2_M)))
  iwres_3_M_TS <- iwres_fc(data=data3$Long$M_Spike, start=data3$Short$start_M, stop=data3$Short$stop_M, 
                           lB=pm3$TS$theta_M[1] + pm3$TS$b1_M, lG=pm3$TS$theta_M[2] + pm3$TS$b2_M, lD=pm3$TS$theta_M[3] + pm3$TS$b3_M, 
                           sigma=mean(sqrt(extract(fit_LONG3_M$fit, "sigma2")$sigma2)))
  # FLC
  iwres_3_F_JM <- iwres_fc(data=data3$Long$FLC, start=data3$Short$start_F, stop=data3$Short$stop_F, 
                           lB=pm3$JM$theta_F[1] + pm3$JM$b1_F, lG=pm3$JM$theta_F[2] + pm3$JM$b2_F, lD=pm3$JM$theta_F[3] + pm3$JM$b3_F, 
                           sigma=mean(sqrt(extract(fit_JM3$fit, "sigma2_F")$sigma2_F)))
  iwres_3_F_TS <- iwres_fc(data=data3$Long$FLC, start=data3$Short$start_F, stop=data3$Short$stop_F, 
                           lB=pm3$TS$theta_F[1] + pm3$TS$b1_F, lG=pm3$TS$theta_F[2] + pm3$TS$b2_F, lD=pm3$TS$theta_F[3] + pm3$TS$b3_F, 
                           sigma=mean(sqrt(extract(fit_LONG3_F$fit, "sigma2")$sigma2)))
  
  # LoT 4
  pm4 <- posterior_mean_fc(fit_JM4$fit, RE_JM4, fit_LONG4_M$fit, fit_LONG4_F$fit, RE_TS4)
  # M-spike
  iwres_4_M_JM <- iwres_fc(data=data4$Long$M_Spike, start=data4$Short$start_M, stop=data4$Short$stop_M, 
                           lB=pm4$JM$theta_M[1] + pm4$JM$b1_M, lG=pm4$JM$theta_M[2] + pm4$JM$b2_M, lD=pm4$JM$theta_M[3] + pm4$JM$b3_M, 
                           sigma=mean(sqrt(extract(fit_JM4$fit, "sigma2_M")$sigma2_M)))
  iwres_4_M_TS <- iwres_fc(data=data4$Long$M_Spike, start=data4$Short$start_M, stop=data4$Short$stop_M, 
                           lB=pm4$TS$theta_M[1] + pm4$TS$b1_M, lG=pm4$TS$theta_M[2] + pm4$TS$b2_M, lD=pm4$TS$theta_M[3] + pm4$TS$b3_M, 
                           sigma=mean(sqrt(extract(fit_LONG4_M$fit, "sigma2")$sigma2)))
  # FLC
  iwres_4_F_JM <- iwres_fc(data=data4$Long$FLC, start=data4$Short$start_F, stop=data4$Short$stop_F, 
                           lB=pm4$JM$theta_F[1] + pm4$JM$b1_F, lG=pm4$JM$theta_F[2] + pm4$JM$b2_F, lD=pm4$JM$theta_F[3] + pm4$JM$b3_F, 
                           sigma=mean(sqrt(extract(fit_JM4$fit, "sigma2_F")$sigma2_F)))
  iwres_4_F_TS <- iwres_fc(data=data4$Long$FLC, start=data4$Short$start_F, stop=data4$Short$stop_F, 
                           lB=pm4$TS$theta_F[1] + pm4$TS$b1_F, lG=pm4$TS$theta_F[2] + pm4$TS$b2_F, lD=pm4$TS$theta_F[3] + pm4$TS$b3_F, 
                           sigma=mean(sqrt(extract(fit_LONG4_F$fit, "sigma2")$sigma2)))
  
  dta <- data.frame(time=c(iwres_1_M_JM$time,iwres_1_M_TS$time,iwres_1_F_JM$time,iwres_1_F_JM$time,
                           iwres_2_M_JM$time,iwres_2_M_TS$time,iwres_2_F_JM$time,iwres_2_F_JM$time,
                           iwres_3_M_JM$time,iwres_3_M_TS$time,iwres_3_F_JM$time,iwres_3_F_JM$time,
                           iwres_4_M_JM$time,iwres_4_M_TS$time,iwres_4_F_JM$time,iwres_4_F_JM$time),
                    col=c(rep(data1$Long$M_Spike$patientid,2),rep(data1$Long$FLC$patientid,2),
                          rep(data2$Long$M_Spike$patientid,2),rep(data2$Long$FLC$patientid,2),
                          rep(data3$Long$M_Spike$patientid,2),rep(data3$Long$FLC$patientid,2),
                          rep(data4$Long$M_Spike$patientid,2),rep(data4$Long$FLC$patientid,2)),
                    type=c(rep("JE",nrow(iwres_1_M_JM)),rep("TS",nrow(iwres_1_M_TS)),
                           rep("JE",nrow(iwres_1_F_JM)),rep("TS",nrow(iwres_1_F_TS)),
                           rep("JE",nrow(iwres_2_M_JM)),rep("TS",nrow(iwres_2_M_TS)),
                           rep("JE",nrow(iwres_2_F_JM)),rep("TS",nrow(iwres_2_F_TS)),
                           rep("JE",nrow(iwres_3_M_JM)),rep("TS",nrow(iwres_3_M_TS)),
                           rep("JE",nrow(iwres_3_F_JM)),rep("TS",nrow(iwres_3_F_TS)),
                           rep("JE",nrow(iwres_4_M_JM)),rep("TS",nrow(iwres_4_M_TS)),
                           rep("JE",nrow(iwres_4_F_JM)),rep("TS",nrow(iwres_4_F_TS))),
                    lot=c(rep("LoT 1",nrow(iwres_1_M_JM)+nrow(iwres_1_M_TS)+nrow(iwres_1_F_JM)+nrow(iwres_1_F_TS)),
                          rep("LoT 2",nrow(iwres_2_M_JM)+nrow(iwres_2_M_TS)+nrow(iwres_2_F_JM)+nrow(iwres_2_F_TS)),
                          rep("LoT 3",nrow(iwres_3_M_JM)+nrow(iwres_3_M_TS)+nrow(iwres_3_F_JM)+nrow(iwres_3_F_TS)),
                          rep("LoT 4",nrow(iwres_4_M_JM)+nrow(iwres_4_M_TS)+nrow(iwres_4_F_JM)+nrow(iwres_4_F_TS))),
                    biom=c(rep("M-spike",nrow(iwres_1_M_JM)+nrow(iwres_1_M_TS)),rep("FLC",nrow(iwres_1_F_JM)+nrow(iwres_1_F_TS)),
                           rep("M-spike",nrow(iwres_2_M_JM)+nrow(iwres_2_M_TS)),rep("FLC",nrow(iwres_2_F_JM)+nrow(iwres_2_F_TS)),
                           rep("M-spike",nrow(iwres_3_M_JM)+nrow(iwres_3_M_TS)),rep("FLC",nrow(iwres_3_F_JM)+nrow(iwres_3_F_TS)),
                           rep("M-spike",nrow(iwres_4_M_JM)+nrow(iwres_4_M_TS)),rep("FLC",nrow(iwres_4_F_JM)+nrow(iwres_4_F_TS))),
                    res=c(iwres_1_M_JM$iwres,iwres_1_M_TS$iwres,iwres_1_F_JM$iwres,iwres_1_F_TS$iwres,
                          iwres_2_M_JM$iwres,iwres_2_M_TS$iwres,iwres_2_F_JM$iwres,iwres_2_F_TS$iwres,
                          iwres_3_M_JM$iwres,iwres_3_M_TS$iwres,iwres_3_F_JM$iwres,iwres_3_F_TS$iwres,
                          iwres_4_M_JM$iwres,iwres_4_M_TS$iwres,iwres_4_F_JM$iwres,iwres_4_F_TS$iwres))
  dta$biom <- factor(dta$biom, levels=c("M-spike","FLC"))

  ggplot(data=dta, aes(x=time, y=res)) + geom_point(aes(x=time, y=res), size=0.5, color="black") +
    geom_hline(yintercept=0, linetype="dashed", color="red") + xlab("Time (in years)") + ylab("IWRES") + theme_bw() +
    theme(legend.position="none", panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank()) +
    facet_nested(biom~lot+type, scales="free_y")

  # ggsave("fig.png", units="in", width=11, height=8, dpi=300)
}

