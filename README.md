# A Bayesian joint model of multiple nonlinear longitudinal and competing risks outcomes for dynamic prediction in multiple myeloma: joint estimation and corrected two-stage approaches

The **Main.R** file includes all the commands (including comments) to run the joint estimation and corrected two-stage approaches for the following joint model:

#### Longitudinal submodels

We specify the longitudinal processes that model M-spike ($k=1$) and free light chains ($k=2$) biomarkers through a bi-exponential model. Mathematically, such a model is given by

$$
\begin{align}
y_{kli}(t) &=& B_{kli} \left[ \exp\\{ G_{kli} t \\} + \exp\\{ -D_{kli} t \\} - 1 \right] + \epsilon_{kli}(t),
\end{align}
$$

where $y_{kli}(t)$ represents the observed value of biomarker $k=1,2$ in line of therapy (LoT) $l=1,2,3,4$ for patient $i=1,...,n_{l}$ at time $t$ ($t=0$ indicates the therapy start time). $B_{kli}$, $G_{kli}$, and $D_{kli}$ are parameters that take only positive values and represent baseline (the biomarker value at $t=0$), growth rate, and decay rate, respectively, which are characteristics associated with the biomarker's longitudinal trajectory. The residual errors, $\epsilon_{kl1}(t),...,\epsilon_{kln_{l}}(t)$, are assumed additive, independent and identically distributed as $\epsilon_{kli}(t) \sim \mbox{Normal}(0,\sigma_{kl}^{2})$. We redefine the three parameters of the longitudinal submodel as $B_{kli} = \exp\\{\theta_{1kl} + b_{kli1}\\}$, $G_{kli} = \exp\\{\theta_{2kl} + b_{kli2}\\}$, and $D_{kli} = \exp\\{\theta_{3kl} + b_{kli3}\\}$, where $\theta_{kl}=(\theta_{1kl},\theta_{2kl},\theta_{3kl})^{\top}$ are population parameters while $b_{kli}=(b_{kli1},b_{kli2},b_{kli3})^{\top}$ are random effects. In addition, we assume that $b_{kli} \sim \mbox{Normal}(0,\Omega_{kl})$, where $\Omega_{kl}$ is an unstructured variance-covariance matrix.

#### Competing risks submodels

We model time-to-death ($v=1$) and time-to-next-LoT ($v=2$) through a competing risks model, via a proportional cause-specific hazard specification. We denote $T_{lvi}$ as the time from the start of LoT $l$ to the occurrence of event $v$ for patient $i$; $C_{li}$ indicates the censoring time for patient $i$ in LoT $l$; $\delta_{li}=0,1,2$ is an event indicator, where $\delta_{li}=0$ represents censoring for both events in LoT $l$, $\delta_{li}=1$ indicates that patient $i$ died in LoT $l$, and $\delta_{li}=2$ that patient $i$ transitioned to LoT $l+1$; and $T_{li} = \min\\{T_{l1i},T_{l2i},C_{li}\\}$ represents the observed event time for patient $i$ in LoT $l$. For a LoT $l$, we specify the hazard function of patient $i$ for event $v$ at time $t$ given by

$$
\begin{align}
h_{lvi}(t) &=& h_{lv0}(t) \exp\left[ X_{lvi}^{\top}\beta_{lv} + \sum_{k=1}^{2}\left(\alpha_{klv1}B_{kli}^{\ast} + \alpha_{klv2}G_{kli}^{\ast} + \alpha_{klv3}D_{kli}^{\ast}\right) \right],
\end{align}
$$

where $h_{lv0}(t)$ represents a baseline hazard function and is defined throughout this work as a Weibull hazard given by $h_{lv0}(t) = \phi_{lv} t^{\phi_{lv}-1}\exp\\{\beta_{lv0}\\}$, where $\phi_{lv}$ and $\beta_{lv0}$ are shape and log-scale parameters; $X_{lvi}$ is a covariate vector with coefficients $\beta_{lv}$; $B_{kli}^{\ast} = \log(B_{kli})$, $G_{kli}^{\ast} = \log(G_{kli})$, and $D_{kli}^{\ast} = \log(D_{kli})$ are the baseline, growth rate, and decay rate (in log scale) of biomarker $k$ in LoT $l$ for patient $i$, shared from the longitudinal submodel, where $\alpha_{klv1}$, $\alpha_{klv2}$, and $\alpha_{klv3}$ have the role of measuring the strength of association between each characteristic of the biomarker trajectory and the risk for event $v$.
