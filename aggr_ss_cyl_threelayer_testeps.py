import site
import sys
#site.addsitedir('\\Users\\Casey\\Dropbox\\research\\code\\thermlib')
#site.addsitedir('\\Desktop\\Casey Patrizio\\repos\\research\\thermolib\\')
site.addsitedir('/Users/cpatrizio/repos/thermolib/')
from wsat import wsat
import findTmoist
from findTmoist import findTmoist
from findLCL0 import findLCL0
from constants import constants as c
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import numpy as np
import matplotlib
import thermo
from constants import constants
from scipy.optimize import fsolve, brentq



g = 9.81
sig = 5.67e-8
T_s = 302 
L_v = 2.257*1e6

S0 = 413.98 #solar constant
theta = 50.5 #zenith angle 
S = np.cos((theta*np.pi/180))*S0 

c_E = 1e-4
rho_w = 1000 #density of water


d = 100000e3 #

p_s = 1000e2 #surface pressure (Pa)
p_t = 200e2 #tropopause (Pa)
p_mid = 500e2
p_BL = 950e2 #boundary layer top (Pa)

domsizes = [768, 1536, 3072]
colors = ['k', 'r', 'g']

q_sat = wsat(T_s, p_s) #mixing ratio above sea surface (100% saturated)    
thetae0 = thermo.theta_e(T_s, p_s, q_sat, 0) #theta_e in moist region 
                                                #use surface temperature to get moist adiabat
                                                
T_a = np.mean(findTmoist(thetae0, np.linspace(p_t, p_mid, 100)))
T_mid = np.mean(findTmoist(thetae0, np.linspace(p_mid, p_BL, 100)))
T_BLtop = findTmoist(thetae0, p_BL) #temperature of boundary layer top 
T_t = findTmoist(thetae0, p_t) #temperature of tropopause (outflow region)
T_midtop = findTmoist(thetae0, p_mid) #temperature at top of mid layer
T_BL = (T_s + T_BLtop)/2. #temperature of boundary layer, consistent with well-mixed assumption (linear mixing) 

q_BLsat = wsat(T_BL, (p_s + p_BL)/2.)
q_BLtopsat = wsat(T_BLtop, p_BL)
    
q_FA = wsat(T_t, p_t) #free troposphere water vapor mixing ratio

thickness = lambda p: (c.Rd/(p*g))*findTmoist(thetae0, p)
delz_mid = integrate.quad(thickness, p_mid, p_BL)[0]
delz_FA = integrate.quad(thickness, p_t, p_mid)[0]

rho_BLtop = p_BL/(c.Rd*T_BLtop)
rho_midtop = p_BL
rho_s = p_s/(c.Rd*T_s)
rho_BL = (rho_BLtop + rho_s)/2. #density of boundary layer 
M_BL = (p_s - p_BL)/g #mass of bo+ boundary layer in kg m^-2
delz_BL = M_BL/rho_BL


dels_mid = c.cpd*(T_BLtop - T_midtop) + g*(-delz_mid) #difference in dry static energy between boundary layer top and tropopause (J)
dels_FA = c.cpd*(T_midtop - T_t) + g*(-delz_FA)
        
s_BLtop = c.cpd*T_BLtop + g*delz_BL #dry static energy at top of boundary layer
s_surf = c.cpd*T_s #dry static energy at surface
s_BL = (s_BLtop + s_surf)/2. #dry static energy of BL (well-mixed)
dels_BL = s_BL - s_BLtop #difference in dry static energy between BL and right above BL

#eps_mid = 1
#eps_a = 0.5

eps_mids = np.linspace(0, 1, 50)
eps_as = np.linspace(0, 1, 50)

omega_dFA = np.zeros((eps_mids.size, eps_as.size))
omega_dBL = np.zeros((eps_mids.size, eps_as.size))
eps_BL = np.zeros((eps_mids.size, eps_as.size))

for i, eps_mid in enumerate(eps_mids):
    
    for k, eps_a in enumerate(eps_as):
    
        BLcontinuity = lambda eps_BL: (eps_mid*sig*T_a**4 + (1-eps_BL)*eps_mid*sig*T_s**4 + eps_BL*eps_mid*sig*T_BL**4)/dels_mid \
                                    - (eps_BL*sig*T_s**4 + eps_BL*eps_mid*sig*T_mid**4 + eps_BL*(1-eps_BL)*sig*T_a**4 - 2*eps_BL*sig*T_BL**4)/dels_BL
                                    
        eps_BL[i,k] = fsolve(BLcontinuity, 0.5) 
        
        omega_dFA[i,k] = g*((1-eps_BL[i,k])*(1-eps_mid)*sig*T_s**4 + eps_BL[i,k]*(1-eps_mid)*sig*T_BL**4 + eps_mid*sig*T_mid**4 - 2*eps_a*sig*T_a**4)/(dels_FA)
        omega_dBL[i,k] = g*(eps_mid*eps_a*sig*T_a**4 + (1-eps_BL[i,k])*eps_mid*sig*T_s**4 + eps_BL[i,k]*eps_mid*sig*T_BL**4 - 2*eps_mid*sig*T_mid**4)/dels_mid

#print 'omega_dBL', omega_dBL
#print 'omega_dFA', omega_dFA

epsmidplot, epsaplot = np.meshgrid(eps_mids, eps_as)

maxomegaFA = np.max(omega_dFA)

maxomegaBL = np.max(omega_dBL)


omegaFA_vals = np.linspace(-maxomegaBL, maxomegaBL, 50)

omegaBL_vals = np.linspace(-maxomegaBL, maxomegaBL, 50)

eps_vals = np.linspace(0, 1, 50)

f, axarr = plt.subplots(3,1)
im = axarr[0,].contourf(epsmidplot, epsaplot, np.transpose(omega_dFA), omegaFA_vals, cmap='RdBu_r')
axarr[0,].set_xlabel(r'$\epsilon_{mid}$')
axarr[0,].set_ylabel(r'$\epsilon_{a}$')
#axarr[0,].set_title(r'$\omega_{d,FA}$')
f.colorbar(im, ax=axarr[0,], label=r'$\omega_{d,FA}$ (Pa/s)')

im = axarr[1,].contourf(epsmidplot, epsaplot, np.transpose(omega_dBL), omegaBL_vals, cmap='RdBu_r')
axarr[1,].set_xlabel(r'$\epsilon_{mid}$')
axarr[1,].set_ylabel(r'$\epsilon_{a}$')
#axarr[1,].set_title(r'$\omega_{d,BL}$')
f.colorbar(im, ax=axarr[1,], label=r'$\omega_{d,BL}$ (Pa/s)')

im = axarr[2,].contourf(epsmidplot, epsaplot, np.transpose(eps_BL), eps_vals, cmap='winter_r')
axarr[2,].set_xlabel(r'$\epsilon_{mid}$')
axarr[2,].set_ylabel(r'$\epsilon_{a}$')
#axarr[2,].set_title(r'$\epsilon_{BL}$')
f.colorbar(im, ax=axarr[2,], label=r'$\epsilon_{BL}$')

plt.suptitle(r'$p_{{t}}$ = {:3.0f} hPa, $p_{{mid}}$ = {:3.0f} hPa, $p_{{BL}}$ = {:3.0f} hPa'.format(p_t/1e2, p_mid/1e2, p_BL/1e2))
plt.show()





#for i, domsize in enumerate(domsizes):
#
#    l_d = (1/np.sqrt(2))*domsize*1000
#    
#    r = np.linspace(0, l_d, 1e6)
#
#    #q_FA = 0.01
#    
#    q_FAd = (q_sat + q_FA)/2.
#    q_FAd = q_FA
#    
#    
#    
#    #q_FA = 0.001
#    
#    #p_lcl, T_lcl = findLCL0(f*q_sat, p_s, T_BL)
#    
#    M_trop = (p_BL - p_t)/g #mass of troposphere in kg m^-2
#    M_BL = (p_s - p_BL)/g #mass of boundary layer in kg m^-2
#            
#    #delz_BL = ((c.Rd*T_BL)/g)*np.log(p_s/p_BL) #boundary layer thickness (m)
#    rho_BLtop = p_BL/(c.Rd*T_BLtop)
#    rho_s = p_s/(c.Rd*T_s)
#    rho_BL = (rho_BLtop + rho_s)/2. #density of boundary layer 
#    delz_BL = M_BL/rho_BL
#    
#    zhat = delz_BL/c_E
#    
#    
#    q_BL = q_sat + 2*zhat*(q_FAd - q_sat)/(l_d**2 - r**2)*(r + zhat - (zhat + l_d)*np.exp((r - l_d)/zhat))
#    
#
#    T_c = T_t #cloud top temperature
#        
#    thickness = lambda p: (c.Rd/(p*g))*findTmoist(thetae0, p)
#    delz_trop = integrate.quad(thickness, p_t, p_BL)[0]
#        
#    dels_trop = c.cpd*(T_BLtop - T_t) + g*(-delz_trop) #difference in dry static energy between boundary layer top and tropopause (J)
#        
#    delq = q_sat - q_FA #difference in mixing ratio between sea surface and tropopause
#        
#    s_BLtop = c.cpd*T_BLtop + g*delz_BL #dry static energy at top of boundary layer
#    s_surf = c.cpd*T_s #dry static energy at surface
#        
#    s_BL = (s_BLtop + s_surf)/2. #dry static energy of BL (well-mixed)
#    dels_BL = s_BL - s_BLtop #difference in dry static energy between BL and right above BL
#        
#    K = dels_BL/dels_trop #constant for calculating T_a
#        
#    
#    #effective emission temperature of troposphere (K)
#    #this was derived using the constraint that omega_BL driven by cooling in troposphere must be
#    #equal to omega_BL driven by cooling in boundary layer
#    T_a = ((eps_BL*T_BL**4*(K+2) + T_s**4*(K*(1-eps_BL) - eps_BL))/(eps_BL + 2*K))**0.25
#    #T_a = 270
#        
#    #radiative warming rate in dry region interior (K s^-1)
#    W_d = (eps_BL*sig*T_BL**4 + (1 - eps_BL)*sig*T_s**4 - 2*sig*T_a**4)/(c.cpd*M_trop)
#        
#    #radiative warming rate of moist region interior (K s^-1)
#    W_m = (S - sig*T_c**4)/(c.cpd*M_trop)
#        
#    #radiative warming rate of dry region boundary layer (K s^-1)
#    W_BL = (eps_BL*sig*T_a**4 + eps_BL*sig*T_s**4 - 2*eps_BL*sig*T_BL**4)/(c.cpd*M_BL)
#        
#        
#    #vertical velocity in dry region
#    omega_BL = (g*W_d*c.cpd*M_trop)/(dels_trop)
#    
#    w_BL = omega_BL/(-rho_BL*g) 
#    
#    delp_BL = p_s - p_BL
#    
#    u_BL = -(omega_BL/(2*delp_BL))*((l_d**2 - r**2))/r
#    
#    p_d = (rho_BL/8)*(omega_BL/delp_BL)**2*(r**4 - l_d**4)/r**2
#    
#    p_m = -(rho_BL/8)*(omega_BL/delp_BL)**2*(r**2)
#    
#    #f = 0.7
#    
#    
#    #solve nonlinear equation for f and l_m 
#    
#    def equations(p):
#        l_m, f = p
#        return (l_m - ((omega_BL*l_d**2)/(omega_BL - g*(S - sig*T_c**4)/(dels_trop + L_v*(f*q_sat -  q_FA))))**(0.5), \
#               (f-1)*q_sat - (2*zhat*(q_FAd - q_sat))/(l_d**2 - l_m**2)*(l_m + zhat - (zhat + l_d)*np.exp((l_m - l_d)/zhat)))
#    
#    #l_mfn = lambda x: l_m - ((omega_BL*l_d**2)/(omega_BL - (S - sig*T_c**4)/(dels_trop + L_v*(x*q_sat -  q_FAd))))**(0.5)
#    
#    #func = lambda f: (f-1)*q_sat - (2*zhat*(q_FAd - q_sat))/(l_d**2 - (l_mfn(f))**2)*(l_mfn(f) + zhat - (zhat + l_d)*np.exp((l_mfn(f) - l_d)/zhat))
#    
#    #fguess = 0
#    #fleft = 0
#    #fright = 2
#    
#    #fleft = 10
#    #fright = 500
#    #guess = 230
#    
#    #f = brentq(func, fleft, fright)
#    
#    l_m, f = fsolve(equations, [100e3, 0.5])
#    
#    #l_m = l_mfn(f)
#    
#    l_mi = np.where(r > l_m)[0]
#    
#    omega_m = g*(S - sig*T_c**4)/(dels_trop + L_v*(f*q_sat - q_FA))
#    w_m = -omega_m/(rho_BL*g)
#    
#    P = -(1./g)*omega_m*np.pi*l_m**2*(f*q_sat - q_FA)
#    
#    E = (1./g)*(2*np.pi*omega_BL)*(q_sat - q_FAd)*((l_d**2 - l_m**2)/2. + zhat*(l_d - l_m) - zhat*(zhat + l_d)*(1 - np.exp((l_m - l_d)/zhat)))
#    
#    u_m = omega_m/(2*delp_BL)*r
#    
#    q_BLm = f*q_sat
#    
#    cpw = 4185.5 #heat capacity of water
#    rho_w = 1000 
#    delz_o = 5
#    M_o = rho_w*delz_o
#        
#    #radiative warming rate of ocean surface (K s^-1)
#    W_o = ((S + eps_BL*sig*T_BL**4) + (1-eps_BL)*sig*T_a**4 - sig*T_s**4)/(cpw*M_o)
#    
#    
#    print '-----------------------------------'
#    print 'balances:'
#    print 'E - P:', E - P
#    print 'mass balance:', omega_m*l_m**2 + omega_BL*(l_d**2 - l_m**2)
#    
#    #check that the l_m is consistent with f*q_sat in moist region BL
#    
#    LHS = f*q_sat
#    
#    RHS= q_sat + 2*zhat*(q_FAd - q_sat)/(l_d**2 - l_m**2)*(l_m + zhat - (zhat + l_d)*np.exp((l_m - l_d)/zhat))
#    
#    omega_mP = g*(S - sig*T_c**4 + (L_v*P)/(np.pi*l_m**2))/(dels_trop)
#    
#    print 'RH in moist region BL = {0}'.format(q_BLm/q_BLsat)
#    print 'water vapor given by f:', RHS
#    print 'water vapor given by l_m:', LHS
#    print 'T_s = {0}'.format(T_s)
#    print 'check moisture balance:'
#    print 'E = {0} (kg s^-1)'.format(E)
#    print 'P/(pi*l_m**2*rho_w) = {0} (mm/day)'.format((P*3600*24*1e3)/(rho_w*np.pi*l_m**2))
#    print 'P = {0} (kg s^-1)'.format(P)
#    #print 'qBL_Bar = {0}, q_FA = {1}'.format(qBL_bar, q_FA)
#    print 'qsat = {0}'.format(q_sat)
#    print 'error = {0}', np.abs(E-P)/P
#    print ''
#    print 'check mass balance:'
#    massdown = omega_BL*np.pi*(l_d**2 - l_m**2)
#    massup = -omega_m*np.pi*l_m**2
#    print 'omega_BL * pi * l_d^2 = {0} (Pa m^2 s^-1) '.format(massdown)
#    print 'omega_m * pi * l_m^2 = {0} (Pa m^2 s^-1)'.format(massup)
#    print 'omega_m given by Precip = {0} (Pa s^-1)'.format(omega_mP)
#    print 'omega_m = {0} (Pa s^-1)'.format(omega_m)
#    print 'omega_BL = {0} (Pa s^-1)'.format(omega_BL)
#    print 'w_BL = {0} m s^-1'.format(w_BL)
#    print 'w_m = {0} m s^-1'.format(w_m)
#    print 'l_m = {0} km'.format(l_m/1e3)
#    print 'l_d = {0} km'.format(l_d/1e3)
#    print 'convective fractional area = {0}'.format((l_m/l_d)**2)
#    print 'M_BL = {0} kg/m^2, M_trop = {1} kg/m^2'.format(M_BL, M_trop)
#    print 'error = {0}'.format(np.abs((massdown - massup)/massdown))
#    print ''
#    evap_cool = L_v*E
#    o_warm =cpw*W_o*np.pi*(l_d**2 - l_m**2)*M_o
#    print 'check energy balance:'
#    print 'radiative warming in convective region free troposphere (K/day)', W_m*3600*24
#    print 'radiative cooling in dry region free troposphere (K/day)', W_d*3600*24
#    print 'evaporative cooling = {0} (J  s^-1)'.format(evap_cool)
#    print 'ocean radiative warming = {0} (J  s^-1)'.format(o_warm) 
#    print 'surface energy unbalance = {0}'.format(o_warm - evap_cool)
#    print ''
#    moist_warming = g*(S - sig*T_c**4)*l_m**2 + (L_v*g*P)/(np.pi)
#    dry_cooling = -g*(eps_BL*sig*T_BL**4 + (1 - eps_BL)*sig*T_s**4 - 2*sig*T_a**4)*(l_d**2 - l_m**2)
#    print 'troposphere radiative warming + condensation warming in moist region  = {0} (J  s^-1)'.format(moist_warming)
#    print 'dry region radiative cooling in free troposphere = {0} (J m^-1 s^-1)'.format(dry_cooling)
#    print 'error = {0}'.format(np.abs((moist_warming - dry_cooling)/dry_cooling))
#    print 'T_a = {0} (K)', T_a
#    
#    
#    plt.figure(1)
#    plt.plot(r/(1e3*domsize), u_BL, color = colors[i])
#    plt.ylim(0, 20)
#    plt.xlabel(r'$\hat{r}$')
#    plt.ylabel('radial wind velocity in dry region BL (m/s)')
#    plt.show()
#    
#    plt.figure(2)
#    plt.plot(r/(1e3*domsize), u_m, color = colors[i])
#    plt.ylim(0, 1)
#    plt.xlabel(r'$\hat{r}$')
#    plt.ylabel('radial wind velocity in moist region BL (m/s)')
#    plt.show()
#    
#    plt.figure(3)
#    plt.plot(r/(1e3*domsize), q_BL*1e3, color = colors[i])
#    plt.xlabel(r'$\hat{r}$')
#    plt.ylabel('BL water vapor mixing ratio (g/kg)')
#    plt.show()
#    
#    plt.figure(4)
#    plt.plot(r/(1e3*domsize), p_d, color = colors[i])
#    plt.xlabel(r'$\hat{r}$')
#    #plt.xlim(600, 768)
#    plt.ylim(-5,1)
#    plt.ylabel('surface pressure perturbation in dry region (Pa)')
#    plt.show()
#    
#    plt.figure(5)
#    plt.plot(r/(1e3*domsize), p_m, color = colors[i])
#    plt.xlabel(r'$\hat{r}$')
#    plt.ylim(-2,0)
#    plt.ylabel('surface pressure perturbation in moist region (Pa)')
#    plt.show()

