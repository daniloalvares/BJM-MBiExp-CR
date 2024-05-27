# FITTED AVERAGE LONGITUDINAL TRAJECTORY
plot_avg_trajectory <- function(fit_JM1, fit_JM2, fit_JM3, fit_JM4,
                                fit_LONG1_M, fit_LONG2_M, fit_LONG3_M, fit_LONG4_M,
                                fit_LONG1_F, fit_LONG2_F, fit_LONG3_F, fit_LONG4_F){
  
  dta_JM <- avg_trajectory(fit1=fit_JM1$fit, fit2=fit_JM2$fit, fit3=fit_JM3$fit, fit4=fit_JM4$fit, approach="JE")
  dta_TS <- avg_trajectory(fit1=fit_LONG1_M$fit, fit2=fit_LONG2_M$fit, fit3=fit_LONG3_M$fit, fit4=fit_LONG4_M$fit, 
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



