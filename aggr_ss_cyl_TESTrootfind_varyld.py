import site
import sys
#site.addsitedir('\\Users\\Casey\\Dropbox\\research\\code\\thermlib')
site.addsitedir('/Users/cpatrizio/repos/thermolib/')
from wsat import wsat
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

fout = '/Users/cpatrizio/Google Drive/figures/simplemodel/'


g = 9.81
sig = 5.67e-8
T_s = 302 
L_v = 2.257*1e6

eps_a = 0.5
eps_strat = 0.2
eps_BL = 0.9

#S0 = 413.98 #solar constant
#theta = 50.5 #zenith angle 
#S = np.cos((theta*np.pi/180))*S0 

c_E = .001
rho_w = 1000 #density of water


d = 100000e3 #

p_s = 1000e2 #surface pressure (Pa)
p_BL = 950e2 #boundary layer top (Pa)
p_out = 200e2 #tropopause (Pa)
p_t = 150e2 #tropopause (Pa)

domsizes = np.linspace(500, 12000, 40)
#domsizes = []
colors = ['k', 'r', 'g']

w_ms = np.zeros(len(domsizes))
l_ms = np.zeros(len(domsizes))
l_ds = np.zeros(len(domsizes))
massbalance = np.zeros(len(domsizes))
waterbalance = np.zeros(len(domsizes))
Ebalance = np.zeros(len(domsizes))
RHs = np.zeros(len(domsizes))
Ps = np.zeros(len(domsizes))
delhs = np.zeros(len(domsizes))
p_LCLs = np.zeros(len(domsizes))

q_sat = wsat(T_s, p_s) #mixing ratio above sea surface (100% saturated)    
thetae0 = thermo.theta_e(T_s, p_s, q_sat, 0) #theta_e in moist region 
                                            #use surface temperature to get moist adiabat
T_a = np.mean(findTmoist(thetae0, np.linspace(p_s, p_t)))
T_out = findTmoist(thetae0, p_out)
T_BL = T_s
M_BL = (p_s - p_BL)/g #mass of boundary layer in kg m^-2
rho_BL = p_s/(c.Rd*T_BL)
s_surf = c.cpd*T_s #dry static energy at surface

                                            
                                            
#def T_BLfn(p_BL):
#    T_BLtop = findTmoist(thetae0, p_BL) #temperature of boundary layer top 
#    T_BL = np.add(T_s,T_BLtop)/2. #temperature of boundary layer, consistent with well-mixed assumption (linear mixing)
#    return T_BL
    
#def T_afn(p_BL):
#    T_a = np.mean(findTmoist(thetae0, np.linspace(p_BL, p_out)))
#    return T_a
    
def dels_tropfn(p_BL):
    T_BLtop = findTmoist(thetae0, p_BL) #temperature of boundary layer top 
    thickness = lambda p: (c.Rd/(p*g))*findTmoist(thetae0, p)
    delz_trop = integrate.quad(thickness, p_out, p_BL)[0]
    dels_trop = c.cpd*np.subtract(T_BLtop, T_out) + g*(-delz_trop) #difference in dry static energy between boundary layer top and tropopause (J)
    return dels_trop
    
def dels_BLfn(p_BL): 
    T_BLtop = findTmoist(thetae0, p_BL) #temperature of boundary layer top 
    delz_BL = np.subtract(p_s, p_BL)/(g*rho_BL)
    s_BLtop = c.cpd*np.array(T_BLtop) + g*np.array(delz_BL) #dry static energy at top of boundary layer
    s_BL = np.add(s_BLtop,s_surf)/2. #dry static energy of BL (well-mixed)
    dels_BL = np.subtract(s_BL, s_BLtop) #difference in dry static energy between BL and right above BL
    return dels_BL
    
pBL_fn = lambda p_BL: eps_a*(eps_BL*sig*T_BL**4 + (1-eps_BL)*sig*T_s**4 - 2*sig*T_a**4)/dels_tropfn(p_BL) - eps_BL*(eps_a*sig*T_a**4 + sig*T_s**4 - 2*sig*T_BL**4)/dels_BLfn(p_BL)

p_BL = fsolve(pBL_fn, 950*1e2)

T_BLtop = findTmoist(thetae0, p_BL) #temperature of boundary layer top 


q_BLsat = wsat(T_BL, (p_s + p_BL)/2.)
q_BLtopsat = wsat(T_BLtop, p_BL)
                                            
T_a = np.mean(findTmoist(thetae0, np.linspace(p_BL, p_out)))

thickness = lambda p: (c.Rd/(p*g))*findTmoist(thetae0, p)
delz_trop = integrate.quad(thickness, p_out, p_BL)[0]
    
dels_trop = c.cpd*(T_BLtop - T_out) + g*(-delz_trop) #difference in dry static energy between boundary layer top and tropopause (J)

M_trop = (p_BL - p_out)/g #mass of troposphere in kg m^-2
M_BL = (p_s - p_BL)/g #mass of boundary layer in kg m^-2
        
#delz_BL = ((c.Rd*T_BL)/g)*np.log(p_s/p_BL) #boundary layer thickness (m)
delz_BL = M_BL/rho_BL

s_BLtop = c.cpd*T_BLtop + g*delz_BL #dry static energy at top of boundary layer
s_surf = c.cpd*T_s #dry static energy at surface
    
s_BL = (s_BLtop + s_surf)/2. #dry static energy of BL (well-mixed)
dels_BL = s_BL - s_BLtop #difference in dry static energy between BL and right above BL

q_BLsat = wsat(T_BL, (p_s + p_BL)/2.)
q_BLtopsat = wsat(T_BLtop, p_BL)
    
                                         
#for j, domsize in enumerate(domsizes):
#
#    l_d = (1/np.sqrt(2))*domsize*1000
#
#    r = np.linspace(0, l_d, 1e6)
#    
#    
#    #T_t = findTmoist(thetae0, p_out) #temperature of tropopause (outflow region)
#    
#    #T_strat = findTmoist(thetae0, p_out)
#    
#    q_FA = wsat(T_out, p_out) #free troposphere water vapor mixing ratio
#    #q_FA = 0.01
#    
#    q_FAd = q_FA
#    
#    #q_FA = 0.001
#    
#    
#    zhat = delz_BL/c_E
#    
#    
#    q_BL = q_sat + 2*zhat*(q_FAd - q_sat)/(l_d**2 - r**2)*(r + zhat - (zhat + l_d)*np.exp((r - l_d)/zhat))
#    
#    
#    q_FA = wsat(T_out, p_out) #free troposphere water vapor mixing ratio
#    #q_FA = 0.001
#    T_c = T_out #cloud top temperature
#    T_strat = T_out
#        
#        
#    delq = q_sat - q_FA #difference in mixing ratio between sea surface and tropopause
#
#        
#        
#    #radiative warming rate in dry region interior (K s^-1)
#    #W_d = (eps_BL*sig*T_BL**4 + (1 - eps_BL)*sig*T_s**4 - 2*sig*T_a**4)/(c.cpd*M_trop)
#    
#    Q_dtrop = eps_a*(eps_strat*sig*T_strat**4 + eps_BL*sig*T_BL**4 + (1 - eps_BL)*sig*T_s**4 - 2*sig*T_a**4)
#        
#    #radiative warming rate of moist region interior (K s^-1)
#    #W_m = (S - sig*T_c**4)/(c.cpd*M_trop)
#    #Q_mtrop = (S - sig*T_c**4)
#    
#    mheatrate=-0.2 #prescribed heating rate in moist region interior K/day
#    Q_mtrop = c.cpd*M_trop*mheatrate/(86400)
#    
#    #radiative warming rate of dry region boundary layer (K s^-1)
#    #W_BL = (eps_BL*sig*T_a**4 + eps_BL*sig*T_s**4 - 2*eps_BL*sig*T_BL**4)/(c.cpd*M_BL)
#    Q_dBL = (eps_BL*sig**eps_a*T_a**4 + eps_BL*sig*T_s**4 - 2*eps_BL*sig*T_BL**4 + (1 - eps_a)*eps_strat*sig*T_strat**4)
#        
#        
#    #vertical velocity in dry region
#    omega_BL = (g*Q_dtrop)/(dels_trop)
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
#    #solve nonlinear equation for f and l_m 
#    
#    def equations(p):
#        l_m, f = p
#        p_LCL, T_LCL = findLCL0(f*q_sat, p_s, T_BL)
#        T_LCL = findTmoist(thetae0, p_LCL)
#        delz_LCL = integrate.quad(thickness, p_out, p_LCL)[0]
#        dels_LCL = c.cpd*(T_LCL-T_out) + g*(-delz_LCL)
#return (l_m - ((omega_BL*l_d**2)/(omega_BL - (g*Q_mtrop)/(dels_LCL + L_v*(f*q_sat -  q_FA))))**(0.5), \
#               (f-1)*q_sat - (2*zhat*(q_FAd - q_sat))/(l_d**2 - l_m**2)*(l_m + zhat - (zhat + l_d)*np.exp((l_m - l_d)/zhat)))
#    
#    #l_mfn = lambda x: l_m - ((omega_BL*l_d**2)/(omega_BL - g*(S - sig*T_c**4)/(dels_trop + g*L_v*(x*q_sat -  q_FAd))))**(0.5)
#    
#    #RH_func = lambda f: (f-1)*q_sat - (2*zhat*(q_FAd - q_sat))/(l_d**2 - (l_mfn(f))**2)*(l_mfn(f) + zhat - (zhat + l_d)*np.exp((l_mfn(f) - l_d)/zhat))
#    
#    #l_ms = np.linspace(0, 500e3)
#    
#
#    
#    #for j, l_m in enumerate(l_ms):
#    #    
#    #    RH_func = lambda RH: q_sat - (2*zhat*(q_FAd - q_sat))/(l_d**2 - l_m**2)*(l_m + zhat - (zhat + l_d)*np.exp((l_m - l_d)/zhat))
#    #    
#    #    f = RH_func(l_m)
#    
#        #fguess = 0
#        #fleft = 0
#        #fright = 2
#        
#        #fleft = 10
#        #fright = 500
#        #guess = 230
#        
#        #f = brentq(func, fleft, fright)
#        
#    l_m, f = fsolve(equations, [100e3, 0.5], xtol=1e-13, maxfev=100000)
#    
#    p_LCL, T_LCL = findLCL0(f*q_sat, p_s, T_BL)
#    T_LCL = findTmoist(thetae0, p_LCL)
#    delz_LCL = integrate.quad(thickness, p_t, p_LCL)[0]
#    dels_LCL = c.cpd*(T_LCL-T_t) + g*(-delz_LCL)
#    
#    M_LCL = (p_LCL - p_t)/g  #mass of convective region (from LCL to tropopause) in kg m^-2
#    
#    W_m = Q_mtrop/(c.cpd*M_LCL)
#    W_d = Q_dtrop/(c.cpd*M_trop)
#    
#    #l_m = l_mfn(f)
#        
#    l_mi = np.where(r > l_m)[0][0]
#        
#    omega_m = g*(Q_mtrop)/(dels_LCL + L_v*(f*q_sat - q_FA))
#    w_m = -omega_m/(rho_BL*g)
#    
#    P = -(1./g)*omega_m*np.pi*l_m**2*(f*q_sat - q_FA)
#    
#    E = (1./g)*(2*np.pi*omega_BL)*(q_sat - q_FA)*((l_d**2 - l_m**2)/2. + zhat*(l_d - l_m) - zhat*(zhat + l_d)*(1 - np.exp((l_m - l_d)/zhat)))
#    
#    u_m = omega_m/(2*delp_BL)*r
#    
#    q_BLm = f*q_sat
#    
#    q_BL[:l_mi] = q_BLm
#    u_BL[:l_mi] = u_m[:l_mi]
#    
#    cpw = 4185.5 #heat capacity of water
#    rho_w = 1000 
#    delz_o = 5
#    M_o = rho_w*delz_o
#        
#    #radiative warming rate of ocean surface (K s^-1)
#    W_o = ((S + eps_BL*sig*T_BL**4) + (1-eps_BL)*sig*T_a**4 - sig*T_s**4)/(cpw*M_o)
#    
#    moist_warming = g*Q_mtrop*l_m**2 + (g*L_v*P)/np.pi
#    dry_cooling = -g*Q_dtrop*(l_d**2 - l_m**2)*(dels_LCL/dels_trop)
#    
#    waterbalance[j] = (E - P)*(1000*86400)/(rho_w*np.pi*(l_d**2))
#    massbalance[j] =  omega_m*l_m**2 + omega_BL*(l_d**2 - l_m**2)
#    Ebalance[j] = moist_warming - dry_cooling
#    l_ms[j] = l_m
#    l_ds[j] = l_d
#    w_ms[j] = w_m
#    Ps[j] = (P*3600*24*1000)/(rho_w*np.pi*l_m**2)
#    RHs[j] = q_BLm/q_BLsat
#    delhs[j] = dels_LCL + L_v*(f*q_sat -  q_FA)
#    p_LCLs[j] = p_LCL
#    
#    
#    
#    
#    #print '-----------------------------------'
#    #print 'balances:'
#    #print 'E - P:', E - P
#    #print 'mass balance:', omega_m*l_m**2 + omega_BL*(l_d**2 - l_m**2)
#    #print 'p_LCL (hPa)', p_LCL/1e2
#    #
#    #
#    ##check that the l_m is consistent with f*q_sat in moist region BL
#    #
#    #LHS = f*q_sat
#    #
#    #RHS= q_sat + 2*zhat*(q_FAd - q_sat)/(l_d**2 - l_m**2)*(l_m + zhat - (zhat + l_d)*np.exp((l_m - l_d)/zhat))
#    #
#    #omega_mP = g*(S - sig*T_c**4 + (L_v*P/(np.pi*l_m**2)))/(dels_trop)
#    #
#    #print 'RH in moist region BL = {0}'.format(q_BLm/q_BLsat)
#    #print 'water vapor given by f:', RHS
#    #print 'water vapor given by l_m:', LHS
#    #print 'T_s = {0}'.format(T_s)
#    #print 'check moisture balance:'
#    #print 'E = {0} (kg s^-1)'.format(E)
#    #print 'P/(pi*l_m**2*rho_w) = {0} (mm/day)'.format((P*3600*24*1000)/(rho_w*np.pi*l_m**2))
#    #print 'P = {0} (kg s^-1)'.format(P)
#    ##print 'qBL_Bar = {0}, q_FA = {1}'.format(qBL_bar, q_FA)
#    #print 'qsat = {0}'.format(q_sat)
#    #print 'error = {0}', np.abs(E-P)/P
#    #print ''
#    #print 'check mass balance:'
#    #massdown = omega_BL*np.pi*(l_d**2 - l_m**2)
#    #massup = -omega_m*np.pi*l_m**2
#    #print 'omega_BL * pi * l_d^2 = {0} (Pa m^2 s^-1) '.format(massdown)
#    #print 'omega_m * pi * l_m^2 = {0} (Pa m^2 s^-1)'.format(massup)
#    #print 'omega_m given by Precip = {0} (Pa s^-1)'.format(omega_mP)
#    #print 'omega_m = {0} (Pa s^-1)'.format(omega_m)
#    #print 'w_BL = {0} m s^-1'.format(w_BL)
#    #print 'w_m = {0} m s^-1'.format(w_m)
#    #print 'l_m = {0} km'.format(l_m/1e3)
#    #print 'l_d = {0} km'.format(l_d/1e3)
#    #print 'convective fractional area = {0}'.format((l_m/l_d)**2)
#    #print 'M_BL = {0} kg/m^2, M_trop = {1} kg/m^2, M_LCL = {2} km/m^2'.format(M_BL, M_trop, M_LCL)
#    #print 'error = {0}'.format(np.abs((massdown - massup)/massdown))
#    #print ''
#    #evap_cool = L_v*E
#    #o_warm =cpw*W_o*np.pi*(l_d**2 - l_m**2)*M_o
#    #print 'check energy balance:'
#    #print 'radiative warming in convective region free troposphere (K/day)', W_m*3600*24
#    #print 'radiative cooling in dry region free troposphere (K/day)', W_d*3600*24
#    #print 'evaporative cooling = {0} (J  s^-1)'.format(evap_cool)
#    #print 'ocean radiative warming = {0} (J  s^-1)'.format(o_warm) 
#    #print 'surface energy unbalance = {0}'.format(o_warm - evap_cool)
#    #print ''
#    #moist_warming = g*Q_mtrop*l_m**2 + (g*L_v*P)/np.pi
#    #dry_cooling = -g*Q_dtrop*(l_d**2 - l_m**2)*(dels_LCL/dels_trop)
#    #print 'troposphere radiative warming + condensation warming in moist region  = {0} (J  s^-1)'.format(moist_warming)
#    #print 'dry region radiative cooling in free troposphere = {0} (J m^-1 s^-1)'.format(dry_cooling)
#    #print 'error = {0}'.format(np.abs((moist_warming - dry_cooling)/dry_cooling))
#    #print 'T_a = {0} (K)', T_a
#    #
#    #
#    #
#    #
#    #plt.figure(1)
#    #plt.plot(r/(1e3), u_BL, color = colors[i])
#    #plt.ylim(-20, 0)
#    ##plt.xlabel(r'$\hat{r}$')
#    #plt.xlabel('r (km)')
#    #plt.ylabel('radial wind velocity in BL (m/s)')
#    #plt.show()
#    #
#    ##plt.figure(2)
#    ##plt.plot(r/(1e3), u_m, color = colors[i])
#    ##plt.ylim(-10, 0)
#    ###plt.xlabel(r'$\hat{r}$')
#    ##plt.xlabel('r (km)')
#    ##plt.ylabel('radial wind velocity in moist region BL (m/s)')
#    ##plt.show()
#    #
#    #plt.figure(3)
#    #plt.plot(r/(1e3), q_BL*1e3, color = colors[i])
#    ##plt.xlabel(r'$\hat{r}$')
#    #plt.xlabel('r (km)')
#    #plt.ylabel('BL water vapor mixing ratio (g/kg)')
#    #plt.show()
#    #
#    ##plt.figure(4)
#    ##plt.plot(r/(1e3), p_d, color = colors[i])
#    ###plt.xlabel(r'$\hat{r}$')
#    ##plt.xlabel('r (km)')
#    ###plt.xlim(600, 768)
#    ##plt.ylim(-5,1)
#    ##plt.ylabel('surface pressure perturbation in dry region (Pa)')
#    ##plt.show()
#    ##
#    ##plt.figure(5)
#    ##plt.plot(r/(1e3), p_m, color = colors[i])
#    ###plt.xlabel(r'$\hat{r}$')
#    ##plt.xlabel('r (km)')
#    ##plt.ylim(-2,0)
#    ##plt.ylabel('surface pressure perturbation in moist region (Pa)')
#    ##plt.show()
#        
plt.figure(1)
plt.plot(l_ds/1e3, massbalance, 'o-')
plt.xlabel(r'l$_d$ (km)')
plt.ylabel(r'mass down - mass up (Pa m$^{2}$ s$^{-1}$)')
plt.savefig(fout + 'massbalance_varyld.pdf')
plt.close()

plt.figure(2)
plt.plot(l_ds/1e3, waterbalance, 'o-')
plt.xlabel(r'l$_d$ (km)')
plt.ylabel(r'E - P (mm/day)')
plt.savefig(fout + 'waterbalance_varyld.pdf')
plt.close()

plt.figure(3)
plt.plot(l_ds/1e3, Ebalance, 'o-')
plt.xlabel(r'l$_d$ (km)')
plt.ylabel(r'condensation warming - radiative cooling (J s$^{-1}$)')
plt.savefig(fout + 'Ebalance_varyld.pdf')
plt.close()

plt.figure(4)
plt.plot(l_ds/1e3, l_ms/1e3, 'o-')
plt.xlabel(r'l$_d$ (km)')
plt.ylabel(r'l$_m$ (km)')
plt.savefig(fout + 'lm_varyld.pdf')
plt.close()

plt.figure(5)
plt.plot(l_ds/1e3, (l_ms**2)/(l_ds**2), 'o-')
plt.xlabel(r'l$_d$ (km)')
plt.ylabel(r'$\sigma$')
plt.savefig(fout + 'sigma_varyld.pdf')
plt.close()

plt.figure(6)
plt.plot(l_ds/1e3, w_ms, 'o-')
plt.xlabel(r'l$_d$ (km)')
plt.ylabel(r'w$_m$ (m/s)')
plt.savefig(fout + 'wm_varyld.pdf')
plt.close()

plt.figure(7)
plt.plot(l_ds/1e3, RHs, 'o-')
plt.ylim(0,1)
plt.xlabel(r'l$_d$ (km)')
plt.ylabel(r'RH$_m$')
plt.savefig(fout + 'RHm_varyld.pdf')
plt.close()

plt.figure(8)
plt.plot(l_ds/1e3, Ps, 'o-')
plt.xlabel(r'l$_d$ (km)')
plt.ylabel(r'P (mm/day)')
plt.savefig(fout + 'P_varyld.pdf')
plt.close()

plt.figure(9)
plt.plot(l_ds/1e3, delhs, 'o-')
plt.xlabel(r'l$_d$ (km)')
plt.ylabel(r'$\Delta h_{LCL}$ (J)')
plt.savefig(fout + 'deltah_varyld.pdf')
plt.close()

plt.figure(10)
plt.plot(l_ds/1e3, p_LCLs/1e2, 'o-')
plt.xlabel(r'l$_d$ (km)')
plt.ylabel(r'$p_{LCL}$ (hPa)')
plt.savefig(fout + 'pLCL_varyld.pdf')
plt.close()


        

