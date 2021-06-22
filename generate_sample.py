import pandas as pd
import numpy as np
from astropy.cosmology import FlatLambdaCDM

def Ez(z, Om, Ol):
	return np.sqrt(Om*(1+z**3) + Ol)

def calc_y(x, A, B, gamma, Ez, xnorm, ynorm):
	return ynorm * A * Ez**gamma * (x/xnorm)**B

np.random.seed(9)

H0 = 70
h = H0/100
cosmo = FlatLambdaCDM(H0=70.0, Om0=0.3, Tcmb0=2.725)

data = pd.read_csv('high_mass.csv')

data['Ez'] = Ez(data.z, 0.3, 0.7)
data['m'] = np.power(10, data.logm)

###calculate temperature from masses
#tm relation from picacs: a = 1.07+/-0.04, b=0.59+/-0.04
data['kt'] = calc_y(data.m, 1.07, 0.59, 2.0/3.0, data['Ez'], 5.0e14, 5.0)

###calculate luminosities from temperatures
#luminosities come from Pratt+09 and are bolometric (0.01-100 keV) and core excised
#lt relation from picacs: a=0.82+/-0.05, b=2.87+/-0.13, gamma=0.42+/-0.09
#adjusted normalisation from 5.0e44 to 8.0e44 in an effort to produce higher L (and
#therefore flux) clusters
data['l'] = calc_y(data.kt, 0.82, 2.87, 0.42, data['Ez'], 5.0, 8.0e44)

#covariance matrix from picacs - Table 6
#covariance measured in log10 space, there are errors
#correlation coefficient for TL, 0.37+/-0.31
#	T L
#T
#L
cov = np.array([[1.8e-3, 2.1e-3], [2.1e-3, 1.5e-2]])

###add intrinsic scatter
int_scatter = np.array([np.random.multivariate_normal(ii, cov) for ii in np.array([np.log10(data.kt), np.log10(data.l)]).T]).T
data['kt_int'] = np.power(10, int_scatter[0])
data['l_int'] = np.power(10, int_scatter[1])

###add observational scatter
cov_obs = np.array([[2.0e-04, 7.0e-05], [7.0e-05, 6.7e-05]])
obs_scatter = np.array([np.random.multivariate_normal(ii, cov_obs) for ii in np.array([np.log10(data.kt_int), np.log10(data.l_int)]).T]).T
data['kt_obs'] = np.power(10, obs_scatter[0])
data['l_obs'] = np.power(10, obs_scatter[1])
data['kt_obs_err'] = [np.sqrt(cov_obs[0][0]) for ii in range(len(data))]
data['l_obs_err'] = [np.sqrt(cov_obs[1][1]) for ii in range(len(data))]
data['covarxy'] = [cov_obs[0][1] for ii in range(len(data))]
data['y_thresh'] = 2.0e44

data_threshed = data[data.l_obs > 2.0e44]

noras = data_threshed.sample(301, random_state=9)

noras.to_csv('noras_simple.csv', index=False)