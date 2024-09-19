# A Bayesian joint model of multiple nonlinear longitudinal and competing risks outcomes for dynamic prediction in multiple myeloma: joint estimation and corrected two-stage approaches

The Main.R file includes all the commands (including comments) to run the joint estimation and corrected two-stage approaches for the following joint model:
$$
\begin{aligned}
  y_{kli}(t) &= B_{kli} \Big[ \exp\left\{ G_{kli} t \right\} + \exp\left\{ -D_{kli} t \right\} - 1 \Big] + \epsilon_{kli}(t), \\
  h_{lvi}(t) &= h_{lv0}(t)\exp\left\{{\bm{X}}_{lvi}^{\top}{\bm{\beta}}_{lv} + \sum_{k=1}^{2}\big(\alpha_{klv1}B_{kli}^{\ast} + \alpha_{klv2}G_{kli}^{\ast} + \alpha_{klv3}D_{kli}^{\ast}\big) \right\},
\end{aligned}
$$
