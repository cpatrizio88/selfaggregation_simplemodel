import numpy as np
from absorb import ABSORB_L93
import os
import glob
import matplotlib.pyplot as plt

#plt.close('all')

cwd = os.getcwd()
#fpath = glob.glob(cwd + '.\\Dropbox\\A622\\Std-atmos.txt')[0]

fpath = glob.glob(cwd + '/Dropbox/A622/Std-atmos.txt')[0]

atmprof = np.loadtxt(fpath, skiprows=3)

z = atmprof[:,0] #height km
z=z*1e3
p = atmprof[:,1] #pressure mb
#p = p*1e2
T = atmprof[:,2] #temperature  K
rho = atmprof[:,3] #dry air density (g/m^3)
rho = rho*1e-3
h2o = atmprof[:,4] #h2o mixing ratio (g/m^3)
h2o = h2o*1e-3
h2onew = (0.10)*h2o

w = h2onew/rho #water vapor mixing ratio (kg/kg)

df = 5
f = np.arange(1, 300+df, df)
CWC = 0

numlevslow=1e4
numlevshigh=1e5
lowh = 2e3
zzlow = np.linspace(0, lowh, numlevslow)
zzhigh = np.linspace(lowh, z[-1], numlevshigh)
zz = np.concatenate([zzlow, zzhigh])

delz = np.diff(zz)

H = 500 #height (m), to calculate upwelling brightness temperature 
        #set at boundary layer height for now

H_index = np.where(zz > H)[0][0]

q_BLbar = 0.01

qv = np.zeros(zz.shape)

#idealized water vapor profile (FROM MY SIMPLE MODEL)
qv[0:H_index] = q_BLbar
qv[H_index:-1] = 0.0005

#calculate absorption coefficient, ABSORB_L93 takes 
#frequency (GHz), temerature (K), pressure (hPa), water vapor mixing ratio (kg/kg) 
#and cloud liquid water content (g/m^3)

#index at at h = 40 km
hi = np.where(np.abs(z - 10e3) < 0.1)[0][0]
lo = 1

beta_alo = np.zeros(f.size)
for i, freq in enumerate(f):
    beta_alo[i] = ABSORB_L93(freq, T[lo], p[lo], w[lo], CWC)
    
beta_ahi = np.zeros(f.size)
for i, freq in enumerate(f):
    beta_ahi[i] = ABSORB_L93(freq, T[hi], p[hi], w[hi], CWC)
    
params = {'legend.fontsize': 12}
plt.rcParams.update(params)
plt.figure(1, figsize=(12, 9)) 
plt.plot(f, beta_alo, 'r', label='near-surface layer (h=1km, p=899 mb)')
plt.plot(f, beta_ahi, 'g', label='high layer (z={0} km, p = {1} mb)'.format(z[hi]/1000., p[hi]))
plt.xlabel('frequency (GHz)', fontsize=11)
plt.ylabel('volumetric absorption coefficient (km^-1)', fontsize=11)
ax = plt.gcf().gca()
freqslabel = np.arange(0,320,20)
freqs_list = ["{0}".format(t) for t in freqslabel]
ax = plt.gcf().gca()
ax.set_xticks(freqslabel)
ax.set_xticklabels(freqs_list)
ax.spines["top"].set_visible(False)  
ax.spines["right"].set_visible(False)
ax.get_xaxis().tick_bottom()  
ax.get_yaxis().tick_left()  
plt.xlim(0,300)
plt.title('volume absorption coefficient in a cloud-free U.S. standard atmosphere')
plt.legend()
plt.show()

Tsurf = input('Enter surface temperature:')
#cosmic radiation background effective temperature
Tsurf = float(Tsurf)
Tcmb=2.7

#f=183 #viewing frequency (GHz), corresponds to 1.3 mm (not sure what band this is)
f=700 #corresponds to  0.43 mm (infrared band)
#f=20000 #corresponds to 15 microns (infrared band)
theta=input('Enter viewing angle:') #viewing zenith angle 
theta = float(theta)
epss = input('Enter surface emissivity:')
print()
epss = float(epss)
print('Calculating TOA upwelling and surface downwelling brightness temperatures...')
print()

h=7 #scale height (km)
p_s = 1000e2

pp=p_s*np.exp(-zz/(h*1e3))

gamma_PH = 6.2 #lapse rate (K km^-1)
zeta_T = 10.5 #km
T = np.zeros(zz.shape)
T_s = Tsurf

def findT(T_s, pl):
    
    zeta = -h*np.log(pl/p_s)
    
    if (zeta < zeta_T):
      T = T_s - gamma_PH*zeta
    else:
      T = T_s - gamma_PH*zeta_T
    
    return T
    
for i, plev in enumerate(delz):
    T[i] = findT(T_s, plev)

interpT = lambda x: np.interp(x, zz, T)
interpP = lambda x: np.interp(x, zz, pp)
interpw = lambda x: np.interp(x, zz, qv)

beta_a = np.zeros(delz.size)
Tmids = np.zeros(delz.size)
zmids = np.zeros(delz.size)
pmids=np.zeros(delz.size)

for i, dz in enumerate(delz):
    zmid= (zz[i+1] + zz[i])/2.
    zmids[i] = zmid
    pmids[i] = p_s*np.exp(-zmid/(h*1e3))
    #pmid = interpP(hmid)
    Tmids[i] = findT(T_s, pmids[i])
    wmid = (qv[i+1] + qv[i])/2.
    #print hmid,pmid,Tmid
    beta_a[i] = ABSORB_L93(f, Tmids[i], pmids[i]/100., wmid, CWC)
    #print 'Beta_a', beta_a[i]
    
taulay  = beta_a*(delz/1000.)

tau_zerotoz = np.cumsum(taulay)

mu = np.cos((theta*np.pi)/180.)

taustar = tau_zerotoz[-1]

tstar = np.exp(-taustar/mu)

epslay = 1 - np.exp(-taulay/mu)

Wdown = epslay*np.exp(-tau_zerotoz/mu)

Tbdown = Tcmb*tstar + np.sum(Tmids*Wdown)

print 'sum of tstar and sum(Wdown) should be 1. It equals: ', np.sum(Wdown) + tstar
print

#now code to find Tbup.. only difference should be tau_ztoinf
tau_ztoinf = np.cumsum(taulay[::-1])
taustarH = tau_ztoinf[H_index]
tstarH = np.exp(-taustarH/mu)
tau_ztoinf = tau_ztoinf[::-1]

Wup = epslay*np.exp(-tau_ztoinf/mu)

print 'sum of tstar and sum(Wup) should be 1. It equals: ', np.sum(Wup) + tstar
Tbup = ((1-epss)*Tbdown+epss*Tsurf)*tstar + np.sum(Tmids*Wup)

TbupH = ((1-epss)*Tbdown+epss*Tsurf)*tstarH + np.sum(Tmids[0:H_index]*Wup[0:H_index])

print
print 'Tbup', Tbup
print 'Tbdown', Tbdown
print 'Tbup at H = {0} m. is: {1}'.format(H, TbupH)


plt.figure(2)
ax = plt.gcf().gca()
ax.set_xscale('log')
plt.plot(Wdown, zmids/1000., label='downwelling')
plt.plot(Wup, zmids/1000., label='upwelling')
plt.xlim(0, 0.001)
plt.title('weighting functions for Schwarzschild Equation with Pierre-humbert temp. profile and idealized water vapor profile, f = {0} GHZ, zenith view angle = {1} degrees'.format(f, theta))
plt.xlabel('weight')
plt.ylabel('height (km)')
plt.legend()
plt.show()

plt.figure(3)
ax = plt.gcf().gca()
plt.plot(h2o, p, label = 'U.S. standard')
plt.plot(qv, pp/100., 'g', label='idealized profile')
ax.invert_yaxis()
plt.title('water vapor profiles')
plt.ylabel('pressure (hPa)')
plt.xlabel('water vapor mixing ratio (g/m^3)')
plt.legend()
plt.show()

plt.figure(4)
ax = plt.gcf().gca()
plt.plot(Tmids, pmids/100., label='Pierrehumbert profile')
plt.ylabel('pressure (hPa)')
plt.xlabel('temperature (K)')
ax.invert_yaxis()
plt.legend()
plt.show()













    
    


    

    


    





