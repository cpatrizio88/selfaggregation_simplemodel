import site
import sys
site.addsitedir('\\Users\\Casey\\Dropbox\\research\\code\\thermlib')
#site.addsitedir('/Users/cpatrizio/Dropbox/research/code/thermlib/')
from wsat import wsat
from esat import esat
from Tdfind import Tdfind
from findLCL0 import findLCL0
from constants import constants as c
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import numpy as np
import matplotlib
import thermo

matplotlib.style.use('ggplot')

sig = 5.67e-8
lv = 2.257e6
g = 9.81
#surface latent heat flux exchange coefficient
c_E = 1e-3
rho_w = 1000 #density of water


d = 100000e3 #length of domain (m)
x = np.linspace(0, d, 1e7)


#T_BL = 301 #boundary layer temperature
eps_BL = 0.5 #emissivity of boundary layer

S0 = 650.83 #solar constant
theta = 50.5 #zenith angle 
S = np.cos((theta*np.pi/180))*S0 

p_s = 1000e2 #surface pressure (Pa)
p_t = 200e2 #tropopause (Pa)
p_BL = 950e2 #boundary layer top (Pa)

delp_BL = p_s - p_BL #thickness of boundary layer (Pa)

#fs = np.linspace(0.80, 0.92, 50)

#for f in fs:
#    p_lcl, T_lcl = findLCL0(f*q_sat, p_s, T_BL)
#    if p_lcl > p_BL:
#      f_max = f
#      break

fs = np.linspace(290, 302, 100)
omega_ms = np.zeros(fs.shape)
omega_BLs = np.zeros(fs.shape)
l_ms = np.zeros(fs.shape)
l_ds = np.zeros(fs.shape)   
Ls = np.zeros(fs.shape) 
p_lcls = np.zeros(fs.shape)
W_ds = np.zeros(fs.shape)
W_BLs = np.zeros(fs.shape)
W_ms = np.zeros(fs.shape)
W_os = np.zeros(fs.shape)
q_BLms = np.zeros(fs.shape)
T_BLs = np.zeros(fs.shape)
T_as = np.zeros(fs.shape)
Es = np.zeros(fs.shape)
o_unbalance = np.zeros(fs.shape)
RH = np.zeros(fs.shape)

delss = np.zeros(fs.shape)
delq_ms = np.zeros(fs.shape)

Ebalance_err = np.zeros(fs.shape)

f=0.8

#the following is the temperature profile used in Pierrehumbert 1995 
#It is a very accurate approximation to a RCE profile in the tropics.
h=7 #scale height (km)
gamma_PH = 6.2 #lapse rate (K km^-1)
zeta_T = 16 #km

def findT(T_s, p):
    
    zeta = -h*np.log(p/p_s)
    
    if (zeta < zeta_T):
      T = T_s - gamma_PH*zeta
    else:
      T = T_s - gamma_PH*zeta_T
    
    return T
    

for i, T_s in enumerate(fs):
    
    q_sat = wsat(T_s, p_s) #mixing ratio above sea surface (100% saturated)
    thetae0 = thermo.theta_e(p_s, T_s, q_sat, 0) #theta_e in moist region 
                                                        #use surface temperature to get moist adiabat
    T_t = findT(T_s, p_t) #temperature of tropopause (outflow region)
    T_BLtop = findT(T_s, p_BL) #temperature of boundary layer top 
    T_BL = (T_s + T_BLtop)/2. #temperature of boundary layer, consistent with well-mixed assumption (linear mixing)
    q_BLsat = wsat(T_BL, (p_s + p_BL)/2.)
    q_BLtopsat = wsat(T_BLtop, p_BL)
    
    p_lcl, T_lcl = findLCL0(f*q_sat, (p_s+p_BL)/2., T_BL)
    if p_lcl > p_BL:
        p_lcl = p_BL
    
    M_trop = (p_BL - p_t)/g #mass of troposphere in kg m^-2
    M_BL = (p_s - p_BL)/g #mass of boundary layer in kg m^-2
            
    #delz_BL = ((c.Rd*T_BL)/g)*np.log(p_s/p_BL) #boundary layer thickness (m)
    rho_BLtop = p_BL/(c.Rd*T_BLtop)
    rho_s = p_s/(c.Rd*T_s)
    rho_BL = (rho_BLtop + rho_s)/2. #density of boundary layer 
    delz_BL = M_BL/rho_BL
        
    q_FA = wsat(T_t, p_t) #free troposphere water vapor mixing ratio
    T_c = T_t #cloud top temperature
        
    thickness = lambda p: (c.Rd/(p*g))*findT(T_s, p)
    delz_trop = integrate.quad(thickness, p_t, p_BL)[0]
        
    dels_trop = c.cpd*(T_BLtop - T_t) + g*(-delz_trop) #difference in dry static energy between boundary layer top and tropopause (J)
        
    delq = q_sat - q_FA #difference in mixing ratio between sea surface and tropopause
        
    s_BLtop = c.cpd*T_BLtop + g*delz_BL #dry static energy at top of boundary layer
    s_surf = c.cpd*T_s #dry static energy at surface
        
    s_BL = (s_BLtop + s_surf)/2. #dry static energy of BL (well-mixed)
    dels_BL = s_BL - s_BLtop #difference in dry static energy between BL and right above BL
        
    K = dels_BL/dels_trop #constant for calculating T_a
        
    eps_BL = 0.6 #boundary layer emissivity (CALCULATE THIS)
    
    #effective emission temperature of troposphere (K)
    #this was derived using the constraint that omega_BL driven by cooling in troposphere must be
    #equal to omega_BL driven by cooling in boundary layer
    T_a = ((eps_BL*T_BL**4*(K+2) + T_s**4*(K*(1-eps_BL) - eps_BL))/(eps_BL + 2*K))**0.25
    #T_a = 270
        
    #radiative warming rate in dry region interior (K s^-1)
    W_d = (eps_BL*sig*T_BL**4 + (1 - eps_BL)*sig*T_s**4 - 2*sig*T_a**4)/(c.cpd*M_trop)
        
    #radiative warming rate of moist region interior (K s^-1)
    W_m = (S - sig*T_c**4)/(c.cpd*M_trop)
        
    #radiative warming rate of dry region boundary layer (K s^-1)
    W_BL = (eps_BL*sig*T_a**4 + eps_BL*sig*T_s**4 - 2*eps_BL*sig*T_BL**4)/(c.cpd*M_BL)
        
    cpw = 4185.5 #heat capacity of water
    rho_w = 1000 
    delz_o = 5
    M_o = rho_w*delz_o
        
    #radiative warming rate of ocean surface (K s^-1)
    W_o = ((S + eps_BL*sig*T_BL**4) + (1-eps_BL)*sig*T_a**4 - sig*T_s**4)/(cpw*M_o)
        
    #vertical velocity in dry region
    omega_BL = (g*W_d*c.cpd*M_trop)/(dels_trop)
    w_BL = omega_BL/(-rho_BL*g)
    
    #horizontal velocity in dry region boundary layer
    u_BL = (omega_BL/delp_BL)*x
        
    #re-moistening length as defined in Wing & Cronin (2015)
    L_m = delz_BL/c_E
        
    #moisture in boundary layer
    q_BL = q_sat - (L_m*delq)*(1 - np.exp(-x/L_m))/x
    
    q_BLm = f*q_sat #moisture in moist region

    #difference in moisture between lifting condensation level and tropopause 
    delq_m = q_BLm - q_FA
    
    #vertical velocity in moist region
    omega_m = (g*W_m*c.cpd*M_trop)/(dels_trop + lv*delq_m)
    w_m = omega_m/(-rho_BL*g)
    
    #find where q_BL == f*q_sat
    hit = q_BL >= q_BLm
    
    hit_indices = np.where(hit)[0]
    moistedge_index = hit_indices[0]
    
    qBL_bar = np.mean(q_BL[1:moistedge_index])
    
    l_d = x[moistedge_index]
    
    E = ((omega_BL*delq)/g)*(l_d + L_m*(np.exp(-l_d/L_m) - 1))
    
    #l_m = -(l_d*W_d*c.cpd*M_trop + lv*E)/(W_m*c.cpd*M_trop)
    l_m = -(omega_BL*l_d)/omega_m
    
    P = -(omega_m/g)*l_m*delq_m
    
    #column relative humidity calculation
    PW = (1./(g*rho_w))*(p_BL*q_BL[1:moistedge_index] + (p_BL - p_t)*q_FA)
    
    #saturation precipitable water
    SPW_fn = lambda p: (1./(g*rho_w)*wsat(findT(T_s, p), p))
    SPW = integrate.quad(SPW_fn, p_t, p_BL)[0] + (1./(g*rho_w))*wsat(T_BL, (p_BL + p_s)/2.)*p_BL
    
    CRH = PW/SPW
    
    print 'T_s = {0}'.format(T_s)
    print 'check moisture balance:'
    print 'E = {0} (kg s^-1 m^-1)'.format(E)
    print 'P = {0} (kg s^-1 m^-1)'.format(P)
    print 'error = {0}', np.abs(E-P)/P
    print ''
    print 'check mass balance:'
    massdown = omega_BL*l_d
    massup = -omega_m*l_m
    print 'omega_BL * l_d = {0} (Pa m s^-1) '.format(massdown)
    print 'omega_m * l_m = {0} (Pa m s^-1)'.format(massup)
    print 'error = {0}'.format(np.abs((massdown - massup)/massdown))
    print ''
    evap_cool = lv*E
    o_warm =cpw*W_o*l_d*M_o
    o_unbalance[i] = o_warm - evap_cool
    print 'check energy balance:'
    print 'evaporative cooling = {0} (J m^-1 s^-1)'.format(evap_cool)
    print 'ocean radiative warming = {0} (J m^-1 s^-1)'.format(o_warm) 
    print 'surface energy unbalance = {0}'.format(o_warm - evap_cool)
    print ''
    moist_warming = c.cpd*W_m*l_m*M_trop 
    dry_cooling = -c.cpd*W_d*l_d*M_trop - lv*E
    print 'troposphere radiative cooling plus evaporative cooling in dry region = {0} (J m^-1 s^-1)'.format(dry_cooling)
    print 'moist radiative warming = {0} (J m^-1 s^-1)'.format(moist_warming)
    print 'error = {0}'.format(np.abs((moist_warming - dry_cooling)/dry_cooling))
    Ebalance_err[i] = np.abs((dry_cooling - moist_warming)/dry_cooling)
    print 'qBL_Bar = {0}, q_FA = {1}'.format(qBL_bar, q_FA)
    print 'T_a = {0} (K)', T_a
    #print 'radiative cooling in dry region + BL = {0} (J m^-1 s^-1)'.format(dry_cooling)
    #print 'radiative warming in moist region {0} (J m^-1 s^-1)'.format(moist_warming)
    #print 'error = {0}'.format(np.abs((dry_cooling - moist_warming)/dry_cooling))
    print '------------------'
    
    L = l_m + l_d
    
    hit_m = np.bitwise_and(x > l_d, x < L)
    
    omega_ms[i] = omega_m
    omega_BLs[i] = omega_BL
    T_BLs[i] = T_BL
    T_as[i] = T_a
    l_ms[i] = l_m
    l_ds[i] = l_d
    Ls[i] = L
    p_lcls[i] = p_lcl
    W_ds[i] = W_d
    W_ms[i] = W_m
    W_os[i] = W_o
    W_BLs[i] = W_BL
    q_BLms[i] = q_BLm
    delss[i] = dels_trop
    delq_ms[i] = delq_m
    Es[i] = E
    q_satBL = wsat(T_BL, (p_s+p_BL)/2.)
    RH[i] = q_BLm/q_satBL

#hit = q_BLms < q_BLtopsat 
#RH = q_BLms/q_BLsat
#p_BLeqLCLindex = np.where(hit)[0][-1]
#p_BLeqLCL = RH[p_BLeqLCLindex]

omega_mminindex = np.where(omega_ms == np.min(omega_ms))[0][0]

Tsomega_mmin = fs[omega_mminindex]

#lcl_below_BL = np.where(p_lcls > p_BL)[0]
#RH_cutoff = RH[lcl_below_BL[0]]
#omega_m_pos = np.where(omega_ms > 0)[0]
#if omega_m_pos:
#    omega_m_switch = RH[omega_m_pos[0]]
#else: 
#    omega_m_switch = RH_cutoff

plt.figure(1)
fig, axarr = plt.subplots(3,2, figsize=(12,9))
fig.set_size_inches(12,9)
axarr[0,0].plot(fs, omega_ms/(-rho_BL*g))
axarr[0,0].set_xlabel('T_s')
axarr[0,0].set_ylabel('w_m (m/s)')
#axarr[0,0].axvline(x = p_BLeqLCL, color='k', label='p_LCL = p_BL')
axarr[0,0].axvline(x = Tsomega_mmin, color='k', alpha=0.4, label='max w_m')
#axarr[0,0].axvline(x = RH_cutoff, color='r')

axarr[0,1].plot(fs, omega_BLs/(-rho_BL*g))
axarr[0,1].set_xlabel('T_s')
axarr[0,1].set_ylabel('w_BL (m/s)')

axarr[1,0].plot(fs, l_ms/1000.)
axarr[1,0].set_xlabel('T_s')
axarr[1,0].set_ylabel('l_m (km)')
#axarr[1,0].axvline(x = p_BLeqLCL, color='k', label='p_LCL = p_BL')
axarr[1,0].axvline(x = Tsomega_mmin, color='k', alpha=0.4, label='max w_m')
#axarr[1,0].axvline(x = omega_m_switch, color='g')
#axarr[1,0].axvline(x = RH_cutoff, color='r')

axarr[1,1].plot(fs, l_ds/1000.)
axarr[1,1].set_xlabel('T_s')
axarr[1,1].set_ylabel('l_d (km)')
#axarr[1,1].axvline(x = p_BLeqLCL, color='k', label='p_LCL = p_BL')
axarr[1,1].axvline(x = Tsomega_mmin, color='k', alpha=0.4, label='max w_m')
#axarr[1,1].axvline(x = omega_m_switch, color='g')
#axarr[1,1].axvline(x = RH_cutoff, color='r')

axarr[2,0].plot(fs, Ls/1000.)
axarr[2,0].set_xlabel('T_s')
#axarr[2,0].axvline(x = p_BLeqLCL, color='k', label='p_LCL = p_BL')
axarr[2,0].set_ylabel('L (km)')

axarr[2,1].plot(fs, l_ds/l_ms)
#axarr[2,1].axvline(x = p_BLeqLCL, color='k', label='p_LCL = p_BL')
axarr[2,1].axvline(x = Tsomega_mmin, color='k', alpha=0.4, label='max w_m')
axarr[2,1].set_xlabel('T_s')
axarr[2,1].set_ylabel('l_d/l_m')
plt.suptitle('f = {} , BL top = {} hPa, delz_BL = {:4.1f} m., eps_BL = {}, latent heat exchange coeff. = {}'.format(f, p_BL/100., delz_BL, eps_BL, c_E))
plt.legend()
plt.draw()
plt.show()
fig.savefig('fig1_lengths_vTs.pdf')
plt.close()



plt.figure(2)
fig, axarr = plt.subplots(3, 1, figsize=(12,9))
#axarr[0,].plot(RH, p_lcls/100.)
#axarr[0,].set_xlabel('RH')
#axarr[0,].set_ylabel('p_lcl')
#axarr[0,0].axvline(x = omega_m_switch, color='g')
#axarr[0,0].axvline(x = RH_cutoff, color='r')

o_warmings = np.multiply(l_ds, cpw*W_os*M_o)
drytrop_coolings = np.multiply(l_ds, c.cpd*W_ds*M_trop)
dryBL_coolings = np.multiply(l_ds, c.cpd*W_BLs*M_BL)
E_coolings = lv*Es
moist_warmings = np.multiply(l_ms, c.cpd*W_ms*M_trop)
#    
#plt.figure(6)
axarr[0,].plot(fs, E_coolings, 'b', label='E cooling')
axarr[0,].plot(fs, drytrop_coolings, 'r', label='dry trop. rad. cooling')
axarr[0,].plot(fs, moist_warmings, 'g', label='moist rad. warming')
#axarr[0,].axvline(x = p_BLeqLCL, color='k', label='p_LCL = p_BL')
#axarr[1,0].axvline(x = omega_m_switch, color='g')
#axarr[1,0].axvline(x = RH_cutoff, color='r')
axarr[0,].set_xlabel('T_s')
axarr[0,].set_ylabel('heating/cooling (J s^-1 m^-1)')
#plt.show()

net_warmingerror = moist_warmings - drytrop_coolings - E_coolings
net_oceanwarming = np.divide(o_unbalance/(cpw*M_o), l_ds)

#plt.figure(7)
axarr[1,].plot(fs, net_oceanwarming*(3600*24))
axarr[1,].set_xlabel('T_s')
axarr[1,].set_ylabel('rad ocean_warming - evap_cooling (K day^-1)')
#axarr[1,].axvline(x = p_BLeqLCL, color='k', label='p_LCL = p_BL')
#plt.show()

#plt.figure(8)
axarr[2,].plot(fs, W_ms*(3600*24), label='W_m', color ='g')
axarr[2,].plot(fs, W_ds*(3600*24), label='W_d', color = 'r')
axarr[2,].plot(fs, W_BLs*(3600*24), label='W_BL', color = 'y')
axarr[2,].plot(fs, W_os*(3600*24), label='W_o', color='b')
#axarr[2,].axvline(x = p_BLeqLCL, color='k', label='p_LCL = p_BL')
axarr[2,].set_ylabel('radiative warming (K day^-1)')
axarr[2,].set_xlabel('RH')
plt.suptitle('f = {} , BL top = {} hPa, delz_BL = {:4.1f} m., eps_BL = {}, latent heat exchange coeff. = {}'.format(f, p_BL/100., delz_BL, eps_BL, c_E))
plt.legend()
plt.show()
fig.savefig('fig2_radiation_vTs.pdf')
plt.close()
##plt.show()
#
#plt.figure(9)

fig = plt.figure(3, figsize=(12,9))
plt.plot(fs, net_warmingerror)
#plt.axvline(x = RHomega_mmin, color= 'k', alpha=0.4, label='max w_m')
plt.ylabel('moist rad. warming - dry trop. rad cooling - E_cooling (J s^-1 m^-1)')
plt.xlabel('T_s')
plt.legend()
fig.savefig('fig3_netraderror_vTs.pdf')
plt.show()
plt.close()


fig = plt.figure(4, figsize=(12,9))
plt.plot(fs, Ebalance_err)
#plt.axvline(x = RH_cutoff, color='r')
#plt.axvline(x = p_BLeqLCL, color='k', label='p_LCL = p_BL')
plt.ylabel('energy balance fractional error: (moist warming - dry_trop_cooling - evap_cooling)/total cooling')
plt.xlabel('T_s')
plt.legend()
fig.savefig('fig4_fracraderror_vTs.pdf')
plt.show()
plt.close()

hit = np.where(lv*delq_ms + delss > 0)[0][0]


fig = plt.figure(6, figsize=(12,9))
plt.title('dels_trop + lv*delq_m = MSE difference between boundary layer top and tropopause')
plt.plot(fs, -delss, color='r', label='-dels_trop')
plt.plot(fs, lv*delq_ms, color='b', label='lv*delq_m')
plt.axvline(x = fs[hit], color='k')
plt.ylabel('energy (J)')
plt.text(fs[hit] + 5, 55000, 'dh/dz > 0')
plt.text(fs[hit] - 5, 55000, 'dh/dz < 0')
plt.xlabel('T_s (K)')
plt.legend()
plt.suptitle('f = {} , BL top = {} hPa, delz_BL = {:4.1f} m., eps_BL = {}, latent heat exchange coeff. = {}'.format(f, p_BL/100., delz_BL, eps_BL, c_E))
fig.savefig('omega_m_negative.pdf')
plt.show()
plt.close()

fig = plt.figure(7, figsize=(12,9))
plt.plot(fs, RH)
plt.xlabel('T_s (K)')
plt.ylabel('RH')
plt.axvline(x = Tsomega_mmin, color='k', alpha=0.4, label='max w_m')
plt.suptitle('f = {} , BL top = {} hPa, delz_BL = {:4.1f} m., eps_BL = {}, latent heat exchange coeff. = {}'.format(f, p_BL/100., delz_BL, eps_BL, c_E))
plt.savefig('RHvsT_s.pdf')
plt.show()
plt.close()

fig = plt.figure(8, figsize=(12,9))
plt.plot(fs, p_lcls/100., 'b--', label='p_lcl')
plt.axhline(p_BL/100., label='p_BL')
plt.xlabel('T_s (K)')
plt.legend()
plt.ylabel('p_lcl (hPa)')
plt.show()
plt.savefig('LCLvsT_s.pdf')
plt.close()






























