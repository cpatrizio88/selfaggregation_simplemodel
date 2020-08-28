from absorb import ABSORB_L93
import site
site.addsitedir('/Users/cpatrizio/Dropbox/research/code/thermlib/')
import numpy as np
import matplotlib.pyplot as plt


To = 50 #initial temperature of atmosphere (K)
N =  40 #number layers
deltap = 1e5/N #thickness of each atmospheric layer (Pa)
g = 9.81 #gravity (m s^-2)
cp = 1004 #heat capacity of air (J kg^-1 K^-1)
Lv = 2.260e6
m = deltap/g #mass of layer per unit area (kg m^-2)
So = 1363 #solar flux density (W m^-2)
sig = 5.67e-8 #stefan-boltzmann constant (W m^-2 K^-4)
alpha = 0.3 #average Earth albedo
S = (So/4)*(1-alpha) #downwelling solar radation averaged over Earth (W m^-2)
R = 287. #gas constant for dry air (J kg^-1 K^-1)
Rv = 461.5 #gas constant for water  vapor (J/(kg K))
ep = R/Rv
gamma_m = 6.5 #moist adiabatic lapse rate
gamma_d = 9.8 #dry adiabatic lapse rate
rho = 1.25
Ce = 1e-3 #slab ocean exchange coeff.
V = 10 #surface wind
cw = 4218 #heat capacity of water
rhow = 1000 #density of water

#pressure levels
plevs = deltap*np.arange(0, N+1)
plevs = plevs[::-1]
p_s = plevs[0]

h_a=7 #scale height (km)
gamma_PH = 7 #lapse rate (K km^-1)
zeta_T = 14.5 #km

T_pierrehumb = np.zeros(plevs.shape)

def findT(T_s, pl):
    
    zeta = -h_a*np.log(pl/p_s)
    
    if (zeta < zeta_T):
      T = T_s - gamma_PH*zeta
    else:
      T = T_s - gamma_PH*zeta_T
    
    return T

zlevs = -h_a*np.log(plevs/p_s)*1e3
delz = np.diff(zlevs)

p_BL = 950e2
q_BLbar = 0.017
q_FA = 0.001

qv = np.zeros(N+1)
qv[plevs > p_BL] = q_BLbar 
qv[plevs <= p_BL] = q_FA


beta_a1 = np.zeros(delz.size)
beta_a2 = np.zeros(delz.size)
beta_a3 = np.zeros(delz.size)
beta_a4 = np.zeros(delz.size)
beta_a5 = np.zeros(delz.size)
Tmids = np.zeros(delz.size)
zmids = np.zeros(delz.size)
pmids=np.zeros(delz.size)
qvmids=np.zeros(delz.size)


f1=183 #viewing frequency (GHz), corresponds to 1.3 mm (not sure what band this is)
f2=700#corresponds to  0.43 mm (infrared band)
f3=1000
f4=20000 #corresponds to 15 microns (infrared band)
f5=40000


CWC = 0
T_s= 302

for i, dz in enumerate(delz):
    zmid= (zlevs[i+1] + zlevs[i])/2.
    zmids[i] = zmid
    pmids[i] = p_s*np.exp(-zmid/(h_a*1e3))
    #pmid = interpP(hmid)
    Tmids[i] = findT(T_s, pmids[i])
    #qvmids[i] = (qv[i+1] + qv[i])/2.
    #print hmid,pmid,Tmid
    beta_a1[i] = ABSORB_L93(f1, Tmids[i], pmids[i]/100., qv[i], CWC)
    beta_a2[i] = ABSORB_L93(f2, Tmids[i], pmids[i]/100., qv[i], CWC)
    beta_a3[i] = ABSORB_L93(f3, Tmids[i], pmids[i]/100., qv[i], CWC)
    beta_a4[i] = ABSORB_L93(f4, Tmids[i], pmids[i]/100., qv[i], CWC)
    beta_a5[i] = ABSORB_L93(f5, Tmids[i], pmids[i]/100., qv[i], CWC)
    
    taulay1  = beta_a1*(delz/1000.)
    taulay2  = beta_a2*(delz/1000.)
    taulay3  = beta_a3*(delz/1000.)
    taulay4  = beta_a4*(delz/1000.)
    taulay5  = beta_a5*(delz/1000.)

theta=0
mu = np.cos((theta*np.pi)/180.)

#longwave emissivity 
epslw1 = np.zeros(N+1)
epslw2 = np.zeros(N+1)
epslw3 = np.zeros(N+1)
epslw4 = np.zeros(N+1)
epslw5 = np.zeros(N+1)

epslw1[0] = 1
epslw2[0] = 1
epslw3[0] = 1
epslw4[0] = 1
epslw5[0] = 1

epslw1[1:] = 1 - np.exp(-taulay1/mu)
epslw2[1:] = 1 - np.exp(-taulay2/mu)
epslw3[1:] = 1 - np.exp(-taulay3/mu)
epslw4[1:] = 1 - np.exp(-taulay4/mu)
epslw5[1:] = 1 - np.exp(-taulay5/mu)

epslw1[-1] = epslw1[-2]
epslw2[-1] = epslw2[-2]
epslw3[-1] = epslw3[-2]
epslw4[-1] = epslw4[-2]
epslw5[-1] = epslw5[-2]

plt.figure(1)
plt.plot(epslw1, plevs/100., 'x-', label='f={0}'.format(f1))
plt.plot(epslw2, plevs/100., 'x-', label='f={0}'.format(f2))
plt.plot(epslw3, plevs/100., 'x-', label='f={0}'.format(f3))
plt.plot(epslw4, plevs/100., 'x-', label='f={0}'.format(f4))
plt.plot(epslw5, plevs/100., 'x-', label='f={0}'.format(f5))
plt.xlabel('eps_lw')
plt.ylabel('p (hPa)')
plt.legend()
plt.title('N = {0}, T_s = {1} K, q_FA = {2}, q_BLbar = {3}, p_BL = {4} hPa'.format(N, T_s, q_FA, q_BLbar, p_BL/100))
plt.gca().invert_yaxis()
plt.show()

    
