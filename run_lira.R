library("lira")
library("readr")

obs <- read_csv("noras_simple.csv")

###accesses individual dataframe columns, returning them as vectors
kt_obs <- log10(obs$kt_obs/6.0)
l_obs <- log10(obs$l_obs/1.0e45)
kt_err <- obs$kt_obs_err
l_err <- obs$l_obs_err
covarxy <- obs$covarxy
y_thresh <- log10(obs$y_thresh/1.0e45)
z <- obs$z

###run lira
mcmc <- lira(x=kt_obs, y=l_obs, delta.x=kt_err, delta.y=l_err, covariance.xy=covarxy, 
y.threshold=y_thresh, z=z, time.factor="Ez", distance="luminosity", Omega.M0=0.3, 
Omega.L0=0.7, n.chains = 8, n.adapt = 1000, n.iter =10000, gamma.YIZ=0.42, 
gamma.mu.Z.Fz=0, gamma.sigma.Z.D="dt")