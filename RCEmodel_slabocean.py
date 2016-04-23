import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker
import matplotlib
from absorb import ABSORB_L93
import site
site.addsitedir('\\Users\\Casey\\Dropbox\\research\\code\\thermlib')
#site.addsitedir('/Users/cpatrizio/Dropbox/research/code/thermlib/')
from tempfile import TemporaryFile
from esat import esat
import sys
import os 
import glob
from wsat import wsat
#import seaborn as sns

#sns.set_style("darkgrid")
#sns.despine()
#sns.set_palette(sns.color_palette())

matplotlib.style.use('ggplot')

To = 50 #initial temperature of atmosphere (K)
N = 20 #number layers
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
Ch = 1e-3 #slab sensible heat ocean exchange coeff.
Ce = 1e-3 #latent heat ocean exchange coeff.
V = 10 #surface wind
cw = 4218 #heat capacity of water
rhow = 1000 #density of water

cwd = os.getcwd()
#fpath = glob.glob(cwd + '.\\Dropbox\\A622\\Std-atmos.txt')[0]

fpath = glob.glob(cwd + '/Dropbox/A622/Std-atmos.txt')[0]

atmprof = np.loadtxt(fpath, skiprows=3)

zstd = atmprof[:,0] #height km
zstd=zstd*1e3
pstd = atmprof[:,1] #pressure mb
#p = p*1e2
Tstd = atmprof[:,2] #temperature  K


#Tinit = np.loadtxt("T0profile")
#Tinit = To*np.ones(N+1)
#T=Tinit
Tinit = np.zeros(N+1)


#timestep for evolution of temperature profile (sec)
deltat = 3600

month = 30*24*deltat #(a month in seconds)
tend = 60*month


#temperature profile, keeps the last convectively adjusted profile
Told = np.ones(N+1)
#net radiation profile
Fnet = np.zeros(N+1)
#pressure levels
plevs = deltap*np.arange(0, N+1)
plevs = plevs[::-1]
p_s = plevs[0]



h_a=7 #scale height (km)
gamma_PH = 7 #lapse rate (K km^-1)
zeta_T = 10.5 #km
T_pierrehumb = np.zeros(plevs.shape)

def findT(T_s, pl):
    
    zeta = -h_a*np.log(pl/p_s)
    
    if (zeta < zeta_T):
      T = T_s - gamma_PH*zeta
    else:
      T = T_s - gamma_PH*zeta_T
    
    return T

#initialize the profile    
T_sinit = 290
for i, plev in enumerate(plevs):
    Tinit[i] = findT(T_sinit, plev)

Tinit = To*np.ones(N+1)

T=Tinit

zlevs = -h_a*np.log(plevs/p_s)*1e3
delz = np.diff(zlevs)


p_BL = 950e2
delz_BL = -h_a*np.log(p_BL/p_s)*1e3
q_BLbar = 0.017
q_FA = 0.001

qv = np.zeros(N+1)
BL_index = np.where(delz_BL < zlevs)[0][0]
qv[plevs >= p_BL] = q_BLbar 
qv[plevs < p_BL] = q_FA

beta_a = np.zeros(delz.size)
Tmids = np.zeros(delz.size)
zmids = np.zeros(delz.size)
pmids=np.zeros(delz.size)
qvmids = np.zeros(delz.size)

#f=183 #viewing frequency (GHz), corresponds to 1.3 mm (not sure what band this is)
f=700 #corresponds to  0.43 mm (infrared band)
#f=20000 #corresponds to 15 microns (infrared band)

CWC = 0

for i, dz in enumerate(delz):
    zmid= (zlevs[i+1] + zlevs[i])/2.
    zmids[i] = zmid
    pmids[i] = p_s*np.exp(-zmid/(h_a*1e3))
    #pmid = interpP(hmid)
    Tmids[i] = findT(T_sinit, pmids[i])
    #print hmid,pmid,Tmid
    beta_a[i] = ABSORB_L93(f, Tmids[i], pmids[i]/100., qv[i], CWC)

taulay  = beta_a*(delz/1000.)

#theta=0
#mu = np.cos((theta*np.pi)/180.)
mubar=0.6

#longwave emissivity 
epslw = np.zeros(N+1)
epslw[0] = 1
epslw[1:] = 1 - np.exp(-taulay/mubar)
epslw[-1] = epslw[-2]

#longwave emissivity 
#epslw = np.zeros(N+1)
#epslw[plevs > 800e2] = 0.3
#epslw[np.bitwise_and(plevs > 400e2, plevs <= 800e2)] = 0.2
#epslw[plevs <= 400e2] = 0.1
#set surface layer longwave emissitvity to 1
#epslw[0] = 1
#epslw = epslw+0.02

#depth of slab ocean
h = 5

Tevo = Tinit

#epssw = np.zeros(N+1)
#shortwave emissivity (with ozone)
epssw = np.ones(N+1)*0.05
for i in range(N):
 epssw[N-i-1] = 0.5*epssw[N-i]
epssw[plevs > 200e2] = 0

#timesteps in sec
time = np.arange(0, tend, deltat)

tt = time/(1.*month)

pgrid, tgrid = np.meshgrid(tt, plevs)

#gamma_crit = 6.5 #critical lapse rate (K/km)

eps = 0.0001 #convergence threshold
Trdiff = (np.abs(Told - T))/Told #relative difference between old and new temperature profiles


for t in time:
#while np.any(Trdiff > eps):
    Trdiff = (np.abs(Told - T))/Told
    Told = T
    #initialize radiation profiles
    Fnet = np.zeros(N+1)
    LWupb = np.zeros(N+1)
    LWupt = np.zeros(N+1)
    LWdownb= np.zeros(N+1)
    LWdownt = np.zeros(N+1)
    SWdownb = np.zeros(N+1)
    SWdownt = np.zeros(N+1)
    
    for k, p in enumerate(plevs):
       
        #find net radiation at each layer
        #N - k is the downwelling index (starts from TOA)
        #k is the upwelling index (starts from surface)
        if k == 0:  #upwelling at surface, downwelling at TOA
            LWupb[k] = 0
            LWupt[k] = epslw[k]*sig*T[k]**4
            LWdownt[N-k] = 0
            LWdownb[N-k] = epslw[N-k]*sig*T[N-k]**4
            SWdownt[N-k] = S
            SWdownb[N-k] = (1-epssw[N-k])*SWdownt[N-k] 
        elif k == N: #downwelling at surface, upwelling at TOA
            LWupb[k] = LWupt[k-1]
            LWupt[k] = epslw[k]*sig*T[k]**4 + (1-epslw[k])*LWupb[k]
            LWdownt[N-k] = LWdownb[N-k+1]
            LWdownb[N-k] = 0
            SWdownt[N-k] = S
            SWdownb[N-k] = 0
        else:
            LWupb[k] = LWupt[k-1]
            LWupt[k] = epslw[k]*sig*T[k]**4 + (1-epslw[k])*LWupb[k]
            LWdownt[N-k] = LWdownb[N-k+1]  
            LWdownb[N-k] = epslw[N-k]*sig*T[N-k]**4 + (1-epslw[N-k])*LWdownt[N-k]
            SWdownt[N-k] = SWdownb[N-k+1]
            SWdownb[N-k] = (1-epssw[N-k])*SWdownt[N-k]
             
    Tslab = T[0]
    Tatm = T[1]     
    #T_BL = np.mean(T[1:BL_index])   
    SH = rho*cp*Ce*V*(Tslab - Tatm) 
    LE = Lv*rho*Ce*V*(wsat(T[0], p_s) - qv[0])
    Q_ocean = 25 #horizontal heat transport by ocean (W m^-2)
    Fnet = LWupb - LWupt + LWdownt - LWdownb + SWdownt - SWdownb 
    Fnet[0] = Fnet[0] - SH - Q_ocean - LE
    Fnet[1] = Fnet[1] + SH + LE
    #distribute the sensible and latent heat evenly (and simultaneously) throughout the boundary layer. 
    #for i in range(1, BL_index):
    #    Fnet[i] = Fnet[i] + (SH)*(delz[i]/delz_BL)
    #    T[i] = T[i] + (deltat*Fnet[i])/(m*cp)
    #T[BL_index:] = T[BL_index:] + (deltat*Fnet[BL_index:])/(m*cp)
    T[1:] = T[1:] + (deltat*Fnet[1:])/(m*cp)
    T[0] = T[0] + (deltat*Fnet[0])/(h*cw*rhow)
    
    

    #T = T + (deltat*Fnet)/(m*cp)
    
    #print 'radiative temp. profile', T
    
    #calculate layer thicknesses
    gamma_crit = np.zeros(plevs.size-1)
    
    ##convective adjustment
    for k in np.arange(1, plevs.size-2):
       
        wsatup = wsat(plevs[k+1], T[k+1])
        wsatlow = wsat(plevs[k], T[k])
        Tdiff = T[k+1] - T[k]
        dwsatdT = (wsatup - wsatlow)/Tdiff
        Tbar = (T[k+1] + T[k])/2.
        pratio = (plevs[k]/(1.*plevs[k+1]))
        delzz =  ((R*Tbar)/g)*np.log(pratio)
        delzz=delzz/1000.
        gamma = -Tdiff/delzz
        
        #gamma_crit[k] = gamma_d/(1+(Lv/cp)*dwsatdT)
        gamma_crit[k] = 6.5
       
    
        #convective adjustment
        #the layer is convectively unstable if gamma > gamma_crit
        if gamma > gamma_crit[k]:
           
            T[k+1] = (T[k]+T[k+1] - gamma_crit[k]*delzz)/2.
            T[k] = T[k+1] + gamma_crit[k]*delzz
            #T[k] = T[k+1] + gamma_crit*delz
    #Tinit = np.loadtxt('T0profile')        
    Tevo = np.column_stack((Tevo, T-Tinit))
    
    print 'after convective adjustment, T', T

np.savetxt("T0profile_slabocean", T)

T_s = T[0]
for i, plev in enumerate(plevs):
    T_pierrehumb[i] = findT(T_s, plev)

Tevo = Tevo[:,1:]

#calculate potential temperature
psurf = 1000*100. #ref. pressure
theta = T[:-1]*(psurf/plevs[:-1])**(R/cp)

#calculate lapse rate
Tbar = (T[:-1] + T[1:])/2.
pratio = (plevs[:-1]/(1.*plevs[1:]))
delz =  ((R*Tbar)/g)*np.log(pratio)
delz=delz/1000.

gamma = -(T[1:]-T[:-1])/delz

#Tinit = np.loadtxt('T0profile')

plt.figure(1)
#plt.subplot(1,2,1)
params = {'legend.fontsize': 12}
plt.rcParams.update(params)
#plt.plot(T, plevs/100., color = 'g', label='+0.02 longwave emissivity')
plt.plot(T, plevs/100., color= 'k', label='RCE profile, fixed SST = {0}, N = {1} layers, q_FA = {2}, q_BLbar = {3}'.format(T[0], N, q_FA, q_BLbar))
plt.plot(T_pierrehumb, plevs/100., color='r', label='Pierre humbert profile, zeta_T = {0} km, gamma_PH = {1} K km^-1'.format(zeta_T, gamma_PH))
#plt.plot(Tstd, pstd/100., color='g', label='Standard atmos.')
ax = plt.gcf().gca()
ax.set_yscale('log')
yticks = np.arange(0,1100,100)
ax.yaxis.set_ticks(yticks)
ax.invert_yaxis()
ax.yaxis.set_ticklabels( ['%1.0f' % i for i in yticks] )
plt.xlabel('temperature (K)', fontsize=12)
plt.ylabel('pressure (hPa)', fontsize=12)
plt.title('radiative-convective equilibrium temperature profile w/ slab-ocean') 
plt.legend()
plt.show()

plt.figure(6)
plt.plot(T_pierrehumb, plevs/100., color='r', label='Pierre humbert profile')
plt.plot(Tstd, pstd, color='g', label='Standard atmos.')
ax = plt.gcf().gca()
ax.yaxis.set_ticks(yticks)
ax.invert_yaxis()
plt.legend()
plt.show()

#plt.figure(2)
#plt.contourf(pgrid, tgrid/100., Tevo, 25, cmap=plt.cm.RdBu_r)
#cb = plt.colorbar()
#cb.set_label('temperarture anomaly (K)')
#ax = plt.gcf().gca()
##ax.set_yscale('log')
##yticks = np.arange(0,1100,100)
##ax.yaxis.set_ticks(yticks)
#ax.invert_yaxis()
##ax.yaxis.set_ticklabels( ['%1.0f' % i for i in yticks] )
#plt.xlabel('time (months)')
#plt.ylabel('pressure (hPa)')
#plt.title('evolution of temperature anomaly in radiative-equilibirium slab ocean climate model after +0.02 LW emissivity switch-on')
#plt.show()

#plt.subplot(1,2,2)
#plt.plot(gamma, plevs[:-1]/100., color='g', label = 'radiative equilibrium')
#plt.axvline(gamma_m, color = 'b', label='moist adiabatic', alpha = 0.2)
#plt.axvline(gamma_d, color = 'r', label='dry adiabatic', alpha=0.2)
#ax = plt.gcf().gca()
#ax.set_yscale('log')
#yticks = np.arange(0,1100,100)T
#ax.yaxis.set_ticks(yticks)
#ax.invert_yaxis()
#ax.yaxis.set_ticklabels( ['%1.0f' % i for i in yticks] )
#plt.xlabel('lapse rate (K/km)', fontsize=12)
#plt.ylabel('pressure (hPa)', fontsize=12)
##plt.title('equilibrium temperature profiles w/ ozone' )
#plt.legend()
#plt.show()

#Tsurf = Tevo[0,:] 
#np.savetxt("Tsurfh1", Tsurf)

#plt.figure(3)
#plt.plot(time/(1.*month), Tsurf, label='h = {0} m'.format(h))
#plt.xlabel('time (months)')
#plt.ylabel('ocean surface temperature anomaly (K)')
#plt.legend()
#plt.draw()
#plt.show()

#Tsurfh1 = np.loadtxt('Tsurfh1')
#Tsurfh2 = np.loadtxt('Tsurfh2')
#Tsurfh5 = np.loadtxt('Tsurfh5')
#Tsurfh20 = np.loadtxt('Tsurfh20')
#Tsurfh50 = np.loadtxt('Tsurfh50')
#Tsurfh100 = np.loadtxt('Tsurfh100')
#Tsurfh150 = np.loadtxt('Tsurfh150')

#sns.set_palette(sns.cubehelix_palette(7, start=.5, rot=-.75), n_colors=7)
#sns.set_palette(sns.color_palette("GnBu_d", n_colors=7), n_colors=7)

#tend = 200*month
#time = np.arange(0, tend, deltat)
#
#plt.figure(3)
#plt.plot(time/(1.*month), Tsurfh1, label='h = 1 m')
#plt.plot(time/(1.*month), Tsurfh2, label='h = 2 m')
#plt.plot(time/(1.*month), Tsurfh5, label='h = 5 m')
#plt.plot(time/(1.*month), Tsurfh20, label='h = 20 m')
#plt.plot(time/(1.*month), Tsurfh50, label='h = 50 m')
#plt.plot(time/(1.*month), Tsurfh100, label='h = 100 m')
#plt.plot(time/(1.*month), Tsurfh150, label='h = 150 m')
#plt.axhline(Tsurfh1[-1]*np.exp(-1), label='e-folding temp.', alpha=0.4)
#plt.xlabel('time (months)')
#plt.ylabel('ocean surface temperature anomaly (K)')
#plt.title('evolution of ocean surface temperature anomaly after +0.02 longwave emissivity switch-on')
#plt.legend(loc='best')
#plt.draw()
#plt.show()
#
#efoldh1i = np.where(Tsurfh1 < Tsurfh1[-1]*np.exp(-1))[0][-1]
#efoldh2i = np.where(Tsurfh2 < Tsurfh2[-1]*np.exp(-1))[0][-1]
#efoldh5i = np.where(Tsurfh5 < Tsurfh5[-1]*np.exp(-1))[0][-1]
#efoldh20i = np.where(Tsurfh20 < Tsurfh20[-1]*np.exp(-1))[0][-1]
#efoldh50i = np.where(Tsurfh50 < Tsurfh50[-1]*np.exp(-1))[0][-1]
#efoldh100i = np.where(Tsurfh100 < Tsurfh50[-1]*np.exp(-1))[0][-1]
#efoldh150i = np.where(Tsurfh150 < Tsurfh50[-1]*np.exp(-1))[0][-1]
#
#tefoldh1 = time[efoldh1i]/(1.*month)
#tefoldh2 = time[efoldh2i]/(1.*month)
#tefoldh5 = time[efoldh5i]/(1.*month)
#tefoldh20 = time[efoldh20i]/(1.*month)
#tefoldh50 = time[efoldh50i]/(1.*month)
#tefoldh100 = time[efoldh100i]/(1.*month)
#tefoldh150 = time[efoldh150i]/(1.*month)
#
#print tefoldh1
#print tefoldh2
#print tefoldh5
#print tefoldh20
#print tefoldh50
#print tefoldh100
#print tefoldh150
#
#efoldtimes = np.array([tefoldh1, tefoldh2, tefoldh5, tefoldh20, tefoldh50, tefoldh100, tefoldh150])
#depths = np.array([1, 2, 5, 20, 50, 100, 150])
#
#sns.set_palette(sns.color_palette())
#
#plt.figure(4)
#plt.scatter(depths, efoldtimes)
#plt.xlabel('ocean depth (m)')
#plt.ylabel('e-folding time (months)')
#plt.title('e-folding time (with respect to the equilbrium temperature after +0.02 longwave emissivity switch-on) for varying ocean depth')
#plt.show()
#
#    
#
#    
##function to calculate saturation specific humidity derivative, dq*/dT
##do in two ways: calculate q* at two temperatures, discrete derivative
##                q* = (eps*es)/(p - (1-eps)*es)
##                where eps = Rd/Rv
##                calculate directly using dq*/dT ~ dw*/dT = (eps*Lv*es)/(p*Rv*T**2)
##                the latter way assumes constant p, and q* ~ (eps*w*)/p






        
        







