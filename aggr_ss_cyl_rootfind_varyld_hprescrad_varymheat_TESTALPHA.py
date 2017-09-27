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

S0 = 413.98 #solar constant
theta = 50.5 #zenith angle 
S = np.cos((theta*np.pi/180))*S0 

c_E = 0.001
rho_w = 1000 #density of water

alpha = 0.05

d = 100000e3 #

p_s = 1000e2 #surface pressure (Pa)
p_t = 200e2 #tropopause (Pa)
p_BL = 900e2 #boundary layer top (Pa)

#eps_a = 1
#eps_BL = 0.6

#domsizes = [500, 1000, 1500, 3000]
domsizes = np.linspace(7000,35e3, 20)
#colors = ['k', 'r', 'g']



mheatrates = np.array([-0.1, -0.2, -0.3])
end = len(mheatrates)
        

for v, mheatrate in enumerate(mheatrates):
    
    print mheatrate
    
    w_ms = np.zeros(len(domsizes))
    l_ms = np.zeros(len(domsizes))
    l_ds = np.zeros(len(domsizes))
    massbalance = np.zeros(len(domsizes))
    waterbalance = np.zeros(len(domsizes))
    #BLcontinuity = np.zeros(len(domsizes))
    w_BLpluss = np.zeros(len(domsizes))
    w_BLminuss = np.zeros(len(domsizes))
    Ebalance = np.zeros(len(domsizes))
    RHs = np.zeros(len(domsizes))
    Ps = np.zeros(len(domsizes))
    delhs = np.zeros(len(domsizes))
    p_LCLs = np.zeros(len(domsizes))


    for j, domsize in enumerate(domsizes):
        
    
        l_d = domsize*1000
        print l_d/1e3
        
        r = np.linspace(0, l_d, 1e6)
        
        
        q_sat = wsat(T_s, p_s) #mixing ratio above sea surface (100% saturated)    
        thetae0 = thermo.theta_e(T_s, p_s, q_sat, 0) #theta_e in moist region 
                                                    #use surface temperature to get moist adiaba
        
        T_BLtop = findTmoist(thetae0, p_BL) #temperature of boundary layer top 
        T_t = findTmoist(thetae0, p_t) #temperature of tropopause (outflow region)
        #T_BL = (T_s + T_BLtop)/2. #temperature of boundary layer, consistent with well-mixed assumption (linear mixing)
        T_BL = T_s
        q_BLsat = wsat(T_BL, (p_s + p_BL)/2.)
        q_BLtopsat = wsat(T_BLtop, p_BL)
        
        q_FA = wsat(T_t, p_t) #free troposphere water vapor mixing ratio
        #q_FA = 0.01
        
        q_FAd = q_FA
        
        
        #q_FA = 0.001
        
        M_trop = (p_BL - p_t)/g #mass of troposphere in kg m^-2
        M_BL = (p_s - p_BL)/g #mass of boundary layer in kg m^-2
                
        #delz_BL = ((c.Rd*T_BL)/g)*np.log(p_s/p_BL) #boundary layer thickness (m)
        rho_BLtop = p_BL/(c.Rd*T_BLtop)
        rho_s = p_s/(c.Rd*T_s)
        rho_BL = (rho_BLtop + rho_s)/2. #density of boundary layer 
        delz_BL = M_BL/rho_BL
        
        zhat = delz_BL/c_E
        
        
        q_BL = q_sat + 2*zhat*(q_FAd - q_sat)/(l_d**2 - r**2)*(r + zhat - (zhat + l_d)*np.exp((r - l_d)/zhat))
        
        
        q_FA = wsat(T_t, p_t) #free troposphere water vapor mixing ratio
        #q_FA = 0.001
        T_c = T_t #cloud top temperature
            
        thickness = lambda p: (c.Rd/(p*g))*findTmoist(thetae0, p)
        delz_trop = integrate.quad(thickness, p_t, p_BL)[0]
            
        dels_trop = c.cpd*(T_BLtop - T_t) + g*(-delz_trop) #difference in dry static energy between boundary layer top and tropopause (J)
            
        delq = q_sat - q_FA #difference in mixing ratio between sea surface and tropopause
            
        s_BLtop = c.cpd*T_BLtop + g*delz_BL #dry static energy at top of boundary layer
        s_surf = c.cpd*T_s #dry static energy at surface
            
        s_BL = (s_BLtop + s_surf)/2. #dry static energy of BL (well-mixed)
        dels_BL = s_BL - s_BLtop #difference in dry static energy between BL and right above BL
            
        #K = dels_BL/dels_trop #constant for calculating T_a
        
        #T_a = np.mean(findTmoist(thetae0, np.linspace(p_t, p_BL, 100)))
        
        #Ta_fn = lambda T_a: eps_a*(eps_BL*sig*T_BL**4 + (1-eps_BL)*sig*T_s**4 - 2*sig*T_a**4)/dels_trop - eps_BL*(eps_a*sig*T_a**4 + sig*T_s**4 - 2*sig*T_BL**4)/(dels_BL)
        
        #T_a = fsolve(Ta_fn, 270)
    
        #T_a = 280
            
        #eps_BL = 0.6 #boundary layer emissivity (CALCULATE THIS)
        
        #effective emission temperature of troposphere (K)
        #this was derived using the constraint that omega_BL driven by cooling in troposphere must be
        #equal to omega_BL driven by cooling in boundary layer
        #FIX THIS
        #T_a = ((eps_BL*T_BL**4*(K+2) + T_s**4*(K*(1-eps_BL) - eps_BL))/(eps_BL + 2*K))**0.25
        #T_a = 270
        
        #eps_BL = eps_a*(2*T_a**4 - T_s**4)/(T_BL**4*(eps_a - 2*K) - T_s**4*(eps_a + K) - eps_a*K*T_a**4)
            
        #radiative warming rate in dry region interior (K s^-1)
        #W_d = (eps_BL*sig*T_BL**4 + (1 - eps_BL)*sig*T_s**4 - 2*sig*T_a**4)/(c.cpd*M_trop)
        
        #Q_dtrop = eps_a*(eps_BL*sig*T_BL**4 + (1 - eps_BL)*sig*T_s**4 - 2*sig*T_a**4)
        
        dheatrate=-1 #prescribed heating rate in moist region interior K/day
        Q_dtrop = (c.cpd*M_trop*dheatrate)/(86400)
            
        #radiative warming rate of moist region interior (K s^-1)
        #W_m = (S - sig*T_c**4)/(c.cpd*M_trop)
        #Q_mtrop = (S - sig*T_c**4)
        #mheatrate=-0.2 #prescribed heating rate in moist region interior K/day
        Q_mtrop = (c.cpd*M_trop*mheatrate)/(86400)
        
        #radiative warming rate of dry region boundary layer (K s^-1)
        #W_BL = (eps_BL*sig*T_a**4 + eps_BL*sig*T_s**4 - 2*eps_BL*sig*T_BL**4)/(c.cpd*M_BL)
        #Q_dBL = (eps_BL*sig**eps_a*T_a**4 + eps_BL*sig*T_s**4 - 2*eps_BL*sig*T_BL**4)
        
        #vertical velocity in dry region
        omega_BL = (g*Q_dtrop)/(dels_trop)
        
        w_BL = omega_BL/(-rho_BL*g) 
        
        delp_BL = p_s - p_BL
        
        u_BL = -(omega_BL/(2*delp_BL))*((l_d**2 - r**2))/r
        
        p_d = (rho_BL/8)*(omega_BL/delp_BL)**2*(r**4 - l_d**4)/r**2
        
        p_m = -(rho_BL/8)*(omega_BL/delp_BL)**2*(r**2)
        
        #f = 0.7
        
        #solve nonlinear equation for f and l_m 
        
        def equations(p):
            l_m, f = p
            p_LCL, T_LCL = findLCL0(f*q_sat, p_s, T_BL)
            T_LCL = findTmoist(thetae0, p_LCL)
            delz_LCL = integrate.quad(thickness, p_t, p_LCL)[0]
            dels_LCL = c.cpd*(T_LCL-T_t) + g*(-delz_LCL)
            return (l_m - ((omega_BL*l_d**2)/(omega_BL - (g*Q_mtrop)/(dels_LCL + (1+ alpha/(1-alpha))*L_v*(f*q_sat -  q_FA))))**(0.5), \
                (f-1)*q_sat - (2*zhat*(q_FAd - q_sat))/(l_d**2 - l_m**2)*(l_m + zhat - (zhat + l_d)*np.exp((l_m - l_d)/zhat)))
        
    
        l_m, f = fsolve(equations, [100e3, 0.9], maxfev=1000000, xtol=1e-16)
        
        p_LCL, T_LCL = findLCL0(f*q_sat, p_s, T_BL)
        T_LCL = findTmoist(thetae0, p_LCL)
        delz_LCL = integrate.quad(thickness, p_t, p_LCL)[0]
        dels_LCL = c.cpd*(T_LCL-T_t) + g*(-delz_LCL)
        
        #l_m = l_mfn(f)
            
        l_mi = np.where(r > l_m)[0][0]
        
        M_LCL = (p_LCL - p_t)/g  #mass of convective region (from LCL to tropopause) in kg m^-2
        
        #W_m = Q_mtrop/(c.cpd*M_LCL)
        #W_d = Q_dtrop/(c.cpd*M_trop)
            
        omega_m = (g*Q_mtrop)/(dels_LCL + (1+ alpha/(1-alpha))*L_v*(f*q_sat - q_FA))
        w_m = -omega_m/(rho_BLtop*g)
        
        w_BL = -omega_BL/(rho_BLtop*g)
        
        #omega_BLplus = (g*Q_dtrop)/dels_trop
        #omega_BLminus = (g*Q_dBL)/dels_BL
        
        #w_BLplus = -omega_BLplus/(rho_BLtop*g)
        #w_BLminus = -omega_BLminus/(rho_BLtop*g)
        
        
        Po = -(1./g)*omega_m*np.pi*l_m**2*(f*q_sat - q_FA)
        P = Po/(1-alpha)
        
        E = (1./g)*(2*np.pi*omega_BL)*(q_sat - q_FA)*((l_d**2 - l_m**2)/2. + zhat*(l_d - l_m) - zhat*(zhat + l_d)*(1 - np.exp((l_m - l_d)/zhat))) + alpha*P
        
        
        u_m = omega_m/(2*delp_BL)*r
        
        q_BLm = f*q_sat
        
        q_BL[:l_mi] = q_BLm
        u_BL[:l_mi] = u_m[:l_mi]
        
        cpw = 4185.5 #heat capacity of water
        rho_w = 1000 
        delz_o = 5
        M_o = rho_w*delz_o
        
        moist_warming = g*Q_mtrop*l_m**2 + (g*L_v*P)/np.pi
        dry_cooling = -g*Q_dtrop*(l_d**2 - l_m**2)*(dels_LCL/dels_trop)
            
        #radiative warming rate of ocean surface (K s^-1)
        #W_o = ((S + eps_BL*sig*T_BL**4) + (1-eps_BL)*sig*T_a**4 - sig*T_s**4)/(cpw*M_o)
        
        waterbalance[j] = (E - P)*(1000*86400)/(rho_w*np.pi*(l_d**2))
        massbalance[j] =  omega_m*l_m**2 + omega_BL*(l_d**2 - l_m**2)
        #w_BLpluss[j] = w_BLplus 
        #w_BLminuss[j] = w_BLminus
        Ebalance[j] = moist_warming - dry_cooling
        l_ms[j] = l_m
        l_ds[j] = l_d
        w_ms[j] = w_m
        Ps[j] = (P*3600*24*1000)/(rho_w*np.pi*l_m**2)
        RHs[j] = q_BLm/q_BLsat
        delhs[j] = dels_LCL + L_v*(f*q_sat -  q_FA)
        p_LCLs[j] = p_LCL
    
        
    plt.figure(1)
    plt.plot(l_ds/1e3, massbalance, 'o-', label=r'$Q_{{m,trop}}$ = {:2.3f} K/day'.format(mheatrate))
    plt.xlabel(r'$l_d$ (km)')
    plt.ylabel(r'mass down - mass up (Pa m$^{2}$ s$^{-1}$)')
    plt.title(r'$Q_{{d,trop}}$ = {:2.1f} K/day, $\alpha$ = {:2.3f}'.format(dheatrate, alpha))
    if v == end-1:
        plt.legend()
        plt.savefig(fout + 'massbalance_varyld_varyQm_alpha{:2.3f}.pdf'.format(alpha))
        plt.close()
    else:
        plt.savefig(fout + 'massbalance_varyld_varyQm_alpha{:2.3f}.pdf'.format(alpha))

    
    plt.figure(2)
    plt.plot(l_ds/1e3, waterbalance, 'o-', label=r'$Q_{{m,trop}}$ = {:2.3f} K/day'.format(mheatrate))
    plt.xlabel(r'$l_d$ (km)')
    plt.ylabel(r'E - P (mm/day)')
    plt.title(r'$Q_{{d,trop}}$ = {:2.1f} K/day, $\alpha$ = {:2.3f}'.format(dheatrate, alpha))
    if v == end-1:
        plt.legend()
        plt.savefig(fout + 'waterbalance_varyld_varyQm_alpha{:2.3f}.pdf'.format(alpha))
        plt.close()
    else:
        plt.savefig(fout + 'waterbalance_varyld_varyQm_alpha{:2.3f}.pdf'.format(alpha))

    
    plt.figure(3)
    plt.plot(l_ds/1e3, Ebalance, 'o-', label=r'$Q_{{m,trop}}$ = {:2.3f} K/day'.format(mheatrate))
    plt.xlabel(r'$l_d$ (km)')
    plt.ylabel(r'condensation warming - radiative cooling (J s$^{-1}$)')
    plt.title(r'$Q_{{d,trop}}$ = {:2.1f} K/day, $\alpha$ = {:2.3f}'.format(dheatrate, alpha))
    if v == end-1:
        plt.legend(loc='best')
        plt.savefig(fout + 'Ebalance_varyld_varyQm_alpha{:2.3f}.pdf'.format(alpha))
        plt.close()
    else:
        plt.savefig(fout + 'Ebalance_varyld_varyQm_alpha{:2.3f}.pdf'.format(alpha))

    
    plt.figure(4)
    plt.plot(l_ds/1e3, l_ms/1e3, 'o-', label=r'$Q_{{m,trop}}$ = {:2.3f} K/day'.format(mheatrate))
    plt.xlabel(r'$l_d$ (km)')
    plt.ylabel(r'$l_m$ (km)')
    plt.title(r'$Q_{{d,trop}}$ = {:2.1f} K/day, $\alpha$ = {:2.3f}'.format(dheatrate, alpha))
    if v == end-1:
        plt.legend(loc='best')
        plt.savefig(fout + 'lm_varyld_varyQm_alpha{:2.3f}.pdf'.format(alpha))
        plt.close()
    else:
        plt.savefig(fout + 'lm_varyld_varyQm_alpha{:2.3f}.pdf'.format(alpha))

    
    plt.figure(5)
    plt.plot(l_ds/1e3, (l_ms**2/l_ds**2), 'o-', label=r'$Q_{{m,trop}}$ = {:2.3f} K/day'.format(mheatrate))
    plt.xlabel(r'$l_d$ (km)')
    plt.ylabel(r'$\sigma$')
    plt.title(r'$Q_{{d,trop}}$ = {:2.1f} K/day, $\alpha$ = {:2.3f}'.format(dheatrate, alpha))
    if v == end-1:
        plt.legend(loc='best')
        plt.savefig(fout + 'sigma_varyld_varyQm_alpha{:2.3f}.pdf'.format(alpha))
        plt.close()
    else:
        plt.savefig(fout + 'sigma_varyld_varyQm_alpha{:2.3f}.pdf'.format(alpha))

    
    plt.figure(6)
    plt.plot(l_ds/1e3, w_ms, 'o-', label=r'$Q_{{m,trop}}$ = {:2.3f} K/day'.format(mheatrate))
    plt.xlabel(r'$l_d$ (km)')
    plt.ylabel(r'w$_m$ (m/s)')
    plt.title(r'$Q_{{d,trop}}$ = {:2.1f} K/day, $\alpha$ = {:2.3f}'.format(dheatrate, alpha))
    if v == end-1:
        plt.legend(loc='best')
        plt.savefig(fout + 'wm_varyld_varyQm_alpha{:2.3f}.pdf'.format(alpha))
        plt.close()
    else:
        plt.savefig(fout + 'wm_varyld_varyQm_alpha{:2.3f}.pdf'.format(alpha))

    plt.figure(7)
    plt.plot(l_ds/1e3, RHs, 'o-', label=r'$Q_{{m,trop}}$ = {:2.3f} K/day'.format(mheatrate))
    plt.ylim(0,1)
    plt.xlabel(r'$l_d$ (km)')
    plt.ylabel(r'RH$_m$')
    plt.title(r'$Q_{{d,trop}}$ = {:2.1f} K/day, $\alpha$ = {:2.3f}'.format(dheatrate, alpha))
    if v == end-1:
        plt.legend(loc='best')
        plt.savefig(fout + 'RHm_varyld_varyQm_alpha{:2.3f}.pdf'.format(alpha))
        plt.close()
    else:
        plt.savefig(fout + 'RHm_varyld_varyQm_alpha{:2.3f}.pdf'.format(alpha))

    plt.figure(8)
    plt.plot(l_ds/1e3, Ps, 'o-', label=r'$Q_{{m,trop}}$ = {:2.3f} K/day'.format(mheatrate))
    plt.xlabel(r'$l_d$ (km)')
    plt.ylabel(r'P (mm/day)')
    plt.title(r'$Q_{{d,trop}}$ = {:2.1f} K/day, $\alpha$ = {:2.3f}'.format(dheatrate, alpha))
    if v == end-1:
        plt.legend(loc='best')
        plt.savefig(fout + 'P_varyld_varyQm_alpha{:2.3f}.pdf'.format(alpha))
        plt.close()
    else:
        plt.savefig(fout + 'P_varyld_varyQm_alpha{:2.3f}.pdf'.format(alpha))

    
    plt.figure(9)
    plt.plot(l_ds/1e3, -delhs, 'o-', label=r'$Q_{{m,trop}}$ = {:2.3f} K/day'.format(mheatrate))
    plt.xlabel(r'$l_d$ (km)')
    plt.ylabel(r'$h_{trop} - h_{LCL}$ (J)')
    plt.title(r'$Q_{{d,trop}}$ = {:2.1f} K/day, $\alpha$ = {:2.3f}'.format(dheatrate, alpha))
    if v == end-1:
        plt.legend(loc='best')
        plt.savefig(fout + 'deltah_varyld_varyQm_alpha{:2.3f}.pdf'.format(alpha))
        plt.close()
    else:
        plt.savefig(fout + 'deltah_varyld_varyQm_alpha{:2.3f}.pdf'.format(alpha))

    
    plt.figure(10)
    plt.plot(l_ds/1e3, p_LCLs/1e2, 'o-', label=r'$Q_{{m,trop}}$ = {:2.3f} K/day'.format(mheatrate))
    plt.xlabel(r'$l_d$ (km)')
    plt.ylabel(r'$p_{LCL}$ (hPa)')
    plt.title(r'$Q_{{d,trop}}$ = {:2.1f} K/day, $\alpha$ = {:2.3f}'.format(dheatrate, alpha))
    if v == end-1:
        plt.legend(loc='best')
        plt.savefig(fout + 'pLCL_varyld_varyQm_alpha{:2.3f}.pdf'.format(alpha))
        plt.close()
    else:
        plt.savefig(fout + 'pLCL_varyld_varyQm_alpha{:2.3f}.pdf'.format(alpha))

    


#plt.figure(10)
#plt.plot(l_ds/1e3, w_BLpluss, 'o-', label=r'$w_{BL+}$')
#plt.plot(l_ds/1e3, w_BLminuss, 'o-', label=r'$w_{BL-}$')
#plt.xlabel(r'l$_d$ (km)')
#plt.ylabel(r'$w$ (m/s)')
#plt.legend()
#plt.savefig(fout + 'BLcont_varyld.pdf')
#plt.close()
        
        
        
        




        

