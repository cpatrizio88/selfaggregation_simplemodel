import site
import sys
#site.addsitedir('\\Users\\Casey\\Dropbox\\research\\code\\thermlib')
site.addsitedir('/Users/cpatrizio/repos/thermolib/')
from wsatnew import wsat
from findTmoist_new import findTmoist
from findLCL0 import findLCL0
from constants import constants as c
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import numpy as np
import matplotlib
import thermo
from constants import constants
from matplotlib.ticker import FormatStrFormatter
from scipy.optimize import fsolve, brentq
fout = '/Users/cpatrizio/figures_arc/simple_model/'

#plt.style.use('seaborn')
matplotlib.rcParams.update({'font.size': 28})
matplotlib.rcParams.update({'figure.figsize': (16, 10)})
matplotlib.rcParams.update({'lines.linewidth': 4})
matplotlib.rcParams.update({'legend.fontsize': 20})
matplotlib.rcParams.update({'mathtext.fontset': 'cm'})

g = 9.81
sig = 5.67e-8
T_s = 302
L_v = 2.257*1e6

S0 = 413.98 #solar constant
theta = 50.5 #zenith angle 
S = np.cos((theta*np.pi/180))*S0 

c_E = 0.001
rho_w = 1000 #density of water


d = 100000e3 #

p_s = 1000e2 #surface pressure (Pa)
p_t = 150e2 #tropopause (Pa)
p_BL = 950e2 #boundary layer top (Pa)

#eps_a = 1
#eps_BL = 0.6

#domsizes = [500, 1000, 1500, 3000]r
domsizes = np.linspace(500, 5000, 50)
#colors = ['k', 'r', 'g']

#domsizes = np.linspace(200, 7000, 30)


#alphas = np.array([0, 0.05, 0.08, 0.10])
#alphas = np.array([0.12,0.08,0.04,0.02])
alphas = np.array([0.12,0.06,0.02])
#alphas = np.array([0])

#alphas = np.array([0.15])

#alphas = np.array([0.05])

#mheatrates = np.array([-0.1, -0.2, -0.3])
end = len(alphas)

#cplot = [(1.0,1.0,1.0), (0.7,0.7,0.7), (0.4,0.4,0.4)]
flag = True  
hs=[]
labels=[]      

for v, alpha in enumerate(alphas):
    
    cplot =  v/(1.0*len(alphas))
    if alpha < 0.1: 
        cplot = str(cplot + 0.1)
    else:
        cplot = str(cplot)
    
    print alpha
    
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
    delhseff = np.zeros(len(domsizes))
    p_LCLs = np.zeros(len(domsizes))
    q_BLms = np.zeros(len(domsizes))
    omega_ms = np.zeros(len(domsizes))
    T_BLcs = np.zeros(len(domsizes))
    p_BLcs = np.zeros(len(domsizes))
    theta_BLcs = np.zeros(len(domsizes))
        
        
    q_sat = wsat(T_s, p_s) #mixing ratio above sea surface (100% saturated)    
    thetae0 = thermo.theta_e(T_s, p_s, q_sat, 0) #theta_e in moist region 
                                                #use surface temperature to get moist adiaba
    
    T_BLtop = findTmoist(thetae0, p_BL) #temperature of boundary layer top 
    T_t = findTmoist(thetae0, p_t) #temperature of tropopause (outflow region)
    #T_BL = (T_s + T_BLtop)/2. #temperature of boundary layer, consistent with well-mixed assumption (linear mixing)
    T_BL = T_s
    #T_BLtop = T_BL
    q_BLsat = wsat(T_BL, (p_s + p_BL)/2.)
    q_BLtopsat = wsat(T_BLtop, p_BL)
    
    q_FA = wsat(T_t, p_t) #free troposphere water vapor mixing ratio
    #q_FA = 0.01
    
    #q_FAd = q_FA
    
    
    #q_FA = 0.001
    
    M_trop = (p_BL - p_t)/g #mass of troposphere in kg m^-2
    M_BL = (p_s - p_BL)/g #mass of boundary layer in kg m^-2
            
    #delz_BL = ((c.Rd*T_BL)/g)*np.log(p_s/p_BL) #boundary layer thickness (m)
    rho_BL = p_BL/(c.Rd*T_BL)
    rho_s = p_s/(c.Rd*T_s)
    #rho_BL = (rho_BLtop + rho_s)/2. #density of boundary layer 
    delz_BL = M_BL/rho_BL
    
    zhat = delz_BL/c_E

    
    #q_BL = q_sat + 2*zhat*(q_FA - q_sat)/(l_d**2 - r**2)*(r + zhat - (zhat + l_d)*np.exp((r - l_d)/zhat))
    
    
    q_FA = wsat(T_t, p_t) #free troposphere water vapor mixing ratio
    #q_FA = 0.001
    T_c = T_t #cloud top temperature
        
    thickness = lambda p: (c.Rd/(p*g))*findTmoist(thetae0, p)
    delz_trop = integrate.quad(thickness, p_t, p_BL)[0]
    
    s_trop = c.cpd*T_t + g*(delz_BL + delz_trop)
    
    s_BLtop = c.cpd*T_BLtop + g*delz_BL #dry static energy at top of boundary layer
    s_surf = c.cpd*T_s #dry static energy at surface
        
    s_BL = (s_BLtop + s_surf)/2.
    
    s_BL = c.cpd*T_BL + (g*delz_BL)/2.
        
    dels_trop = s_BL - s_trop #difference in dry static energy between boundary layer top and tropopause (J)
        
    delq = q_sat - q_FA #difference in mixing ratio between sea surface and tropopause
        
#dry static energy of BL (well-mixed)
    #dels_BL = s_BL - s_BLtop #difference in dry static energy between BL and right above BL
        
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
    mheatrate=-0.6#prescribed heating rate in moist region interior K/day
    Q_mtrop = (c.cpd*M_trop*mheatrate)/(86400)
    
    #radiative warming rate of dry region boundary layer (K s^-1)
    #W_BL = (eps_BL*sig*T_a**4 + eps_BL*sig*T_s**4 - 2*eps_BL*sig*T_BL**4)/(c.cpd*M_BL)
    #Q_dBL = (eps_BL*sig**eps_a*T_a**4 + eps_BL*sig*T_s**4 - 2*eps_BL*sig*T_BL**4)
    
    #vertical velocity in dry region
    omega_BL = (g*Q_dtrop)/(dels_trop)
    
    w_BL = omega_BL/(-rho_BL*g) 
    
    delp_BL = p_s - p_BL
    
    
    #p_d = (rho_BL/8)*(omega_BL/delp_BL)**2*(r**4 - l_d**4)/r**2
    
    #p_m = -(rho_BL/8)*(omega_BL/delp_BL)**2*(r**2)
    
    c_D = 0.001
    coef = (omega_BL/(2*delp_BL))**2    
    
    PI_s = c.cpd*1.0
    PI_BL = c.cpd*(p_BL/p_s)**(c.Rd/c.cpd)   
                
    dels_tropd = c.cpd*(T_s-T_t) + g*(-delz_trop)
    omega_BL = (g*Q_dtrop)/dels_tropd
    #rs = np.linspace(r, l_d, 1e5)
    
    def theta_v(r, l_d):
        rs = np.linspace(r, l_d, 1e5)
        coef = (omega_BL/(2*delp_BL))**2    
        dpdr = -rho_s*coef*((l_d**4 - rs**4)/rs**3 + (c_D/delz_BL)*((l_d**2 - rs**2)**2)/rs**2)
        dTdr = -(2*T_s)/(g*delz_BL)*(1/rho_s)*dpdr
        #dlnTdr = -(2)/(g*delz_BL)*(1/rho_s)*dpdr
        Ts = T_s + np.cumsum(dTdr[::-1][:-1]*np.diff(rs))
        #lnTs = np.log(T_s) + np.cumsum(dlnTdr[::-1][:-1]*np.diff(rs))
        #Ts = np.exp(lnTs)
        return Ts[-1]
        
    def p(r, l_d):
        rs = np.linspace(r, l_d, 1e5)
        coef = (omega_BL/(2*delp_BL))**2    
        #rs = np.linspace(r, l_d, 1e5)
        dpdr = -rho_s*coef*((l_d**4 - rs**4)/rs**3 + (c_D/delz_BL)*((l_d**2 - rs**2)**2)/rs**2)
        #dTdr = -(2/PI_s - PI_BL)*(1/rho_s)*dpdr
        #ps = p_s + np.cumsum(dpdr[::-1][:-1]*np.diff(rs))
        ps = p_s + np.cumsum(dpdr[::-1][:-1]*np.diff(rs))
        #Ts = (ps/(rho_s*c.Rd))
        return ps[-1]
        
    #def dpdr_fn(r):
    #    rs = np.linspace(r, l_d, 1e5)
    #    coef = (omega_BL/(2*delp_BL))**2    
    #    #rs = np.linspace(r, l_d, 1e5)
    #    dpdr = -rho_s*coef*((l_d**4 - rs**4)/rs**3 + (c_D/delz_BL)*((l_d**2 - rs**2)**2)/rs**2)
    #    return dpdr
        
    #def dTdr_fn(r):
    #    #rs = np.linspace(r, l_d, 1e5)
    #    #dTdr = -(2/(PI_s - PI_BL))*(1/rho_s)*dpdr_fn(r)
    #    dlnTdr = ((-g*delz_BL)/2.)*(1/rho_s)*dpdr_fn(r)
    #    return dlnTdr
        
    #def T_test(r):
    #    rs = np.linspace(r, l_d, 1e5)
    #    dTdr = dTdr_fn(r)
    #    return p_s + np.cumsum(dTdr[::-1][:-1]*np.diff(rs))[-1]
        
    #def p_test(r):
    #    rs = np.linspace(r, l_d, 1e5)
    #    dpdr = dpdr_fn(r)
    #    return T_s + np.cumsum(dpdr[::-1][:-1]*np.diff(rs))[-1]
        
    def thetav_to_T(thetav,q):
        T = thetav/(1+ 0.61*q)
        return T
        


    for j, domsize in enumerate(domsizes):
        
        l_d = domsize*1000
        print l_d/1e3
        
        r = np.linspace(0, l_d, 1e6)
     
        
        
        #f = 0.90
        
        #solve nonlinear equation for f and l_m 
        
        def equations(p):
            l_m, f, T_BLc = p
            #p_LCL, T_LCL = findLCL0(f*q_sat, p_s, T_BL)
            #p_LCL = p_BL
            #T_BLc = T_s
            #q_sat = wsat(T_BLc, p_BL)
            dels_tropc = c.cpd*(T_BLc-T_t) + g*(-delz_trop)
            #coef = (omega_BL/(2*delp_BL))**2  
            l_m_test = ((omega_BL*l_d**2)/(omega_BL - (g*Q_mtrop)/(dels_tropc + (1+ alpha/(1-alpha))*L_v*(f*q_sat -  q_FA))))**(0.5)
            thetav_test = theta_v(l_m, l_d)
            q_BLc = f*q_sat
            T_BLc_test = thetav_to_T(thetav_test,q_BLc)
            #q_BLm_test = q_sat + (2*zhat*(q_FA - q_sat))/(l_do*2 - l_m**2)*(l_m + zhat - (zhat + l_d)*np.exp((l_m - l_d)/zhat))
            #T_LCL = T_BLtop
            #T_LCL = findTmoist(thetae0, p_LCL)
            #delz_LCL = integrate.quad(thickness, p_t, p_LCL)[0]
            #dels_LCL = c.cpd*(T_LCL-T_t) + g*(-delz_LCL)
            return (l_m - l_m_test, \
                (f-1)*q_sat - (2*zhat*(q_FA - q_sat))/(l_d**2 - l_m**2)*(l_m + zhat - (zhat + l_d)*np.exp((l_m - l_d)/zhat)), \
                T_BLc - T_BLc_test)
            
    
        #l_m, f, T_BLc = fsolve(equations, [200e3, 0.9, 300], maxfev=100000000, xtol=1e-10)
        l_m, f, T_BLc = fsolve(equations, [1500e3, 0.8, 306], maxfev=100000000, xtol=1e-10)
        
        
        print 'l_m', l_m/1000
        print 'T_BLc', T_BLc
        print 'q_BLc', f*q_sat
        

        
        p_0 = p_s
        c_D = 0.001
        
        p_LCL = p_BL
        
    
        #T_BLc = T_s
        #PI_s = (p_s/p_0)**(c.Rd/c.cpd)
        #PI_BL = (p_BL/p_0)**(c.Rd/c.cpd)
        #p_LCL, T_LCL = findLCL0(f*q_sat, p_s, T_BL)
        #T_LCL = findTmoist(thetae0, p_LCL)
        #delz_trop = integrate.quad(thickness, p_t, p_BL)[0]
        #dels_trop = c.cpd*(T_BLc-T_t) + g*(-delz_trop)
        dels_tropc = c.cpd*(T_BLc-T_t) + g*(-delz_trop)
        #omega_BL = (g*Q_dtrop)/dels_trop
        #coef = (omega_BL/(2*delp_BL))**2  
        #exner = 2/(PI_s - PI_BL)
      
        q_BLm = f*q_sat
        
        #theta_BLc = T_BLc*(1+0.61*q_BLm)

        #f = q_BLm/wsat(T_BLc, p_BL)
        
        #T_BLc = T_s + coef*(-0.5*(2*l_d**2 - (l_d**4 + l_m**4)/l_m**2) + (c_D/delz_BL)*((-8/3)*l_d**3 + (l_d**4 + 2*l_m**2*l_d**2 - (1/3)*l_m**4)/l_m))
        #delz_LCL = integrate.quad(thickness, p_t, p_LCL)[0]
        #dels_LCL = c.cpd*(T_LCL-T_t) + g*(-delz_LCL)
        
        #print 'T_BLc', T_BLc
        q_BL = q_sat + 2*zhat*(q_FA - q_sat)/(l_d**2 - r**2)*(r + zhat - (zhat + l_d)*np.exp((r - l_d)/zhat))
        u_BL = -(omega_BL/(2*delp_BL))*((l_d**2 - r**2))/r
        
        #l_m = l_mfn(f)
    
        if r.max() > l_m:    
            l_mi = np.where(r > l_m)[0][0]
        else:
            l_mi= 0
        
        #M_LCL = (p_LCL - p_t)/g  #mass of convective region (from LCL to tropopause) in kg m^-2
        
        #W_m = Q_mtrop/(c.cpd*M_LCL)
        #W_d = Q_dtrop/(c.cpd*M_trop)
        
        delqv_trop = f*q_sat - q_FA
            
        omega_m = (g*Q_mtrop)/(dels_tropc + (1+ alpha/(1-alpha))*L_v*(delqv_trop))
        w_m = -omega_m/(rho_BL*g)
        
        w_BL = -omega_BL/(rho_BL*g)
        
        #omega_BLplus = (g*Q_dtrop)/dels_trop
        #omega_BLminus = (g*Q_dBL)/dels_BL
        
        #w_BLplus = -omega_BLplus/(rho_BLtop*g)
        #w_BLminus = -omega_BLminus/(rho_BLtop*g)
        
        
        Po = -(1./g)*omega_m*(f*q_sat - q_FA)
        P = Po/(1-alpha)
        
        Eo = (1./g)*(2*np.pi*omega_BL)*(q_sat - q_FA)*((l_d**2 - l_m**2)/2. + zhat*(l_d - l_m) - zhat*(zhat + l_d)*(1 - np.exp((l_m - l_d)/zhat))) 
        E = alpha*P*np.pi*l_m**2 + Eo
        
        
        u_m = omega_m/(2*delp_BL)*r
        
        #q_BLm = f*q_sat
        

        
        q_BL[:l_mi] = q_BLm
        u_BL[:l_mi] = u_m[:l_mi]
        
        cpw = 4185.5 #heat capacity of water
        rho_w = 1000 
        delz_o = 5
        M_o = rho_w*delz_o
        
        moist_warming = (Q_mtrop + L_v*P)*l_m**2
        dry_cooling = -Q_dtrop*(l_d**2 - l_m**2)*(dels_tropc/dels_tropd)
        
        
            
        #radiative warming rate of ocean surface (K s^-1)
        #W_o = ((S + eps_BL*sig*T_BL**4) + (1-eps_BL)*sig*T_a**4 - sig*T_s**4)/(cpw*M_o)
        print 'E - P', E - P*np.pi*l_m**2
        
        waterbalance[j] = E - P*np.pi*l_m**2
        massbalance[j] =  omega_m*l_m**2 + omega_BL*(l_d**2 - l_m**2)
        #w_BLpluss[j] = w_BLplus 
        #w_BLminuss[j] = w_BLminus
        Ebalance[j] = moist_warming - dry_cooling
        l_ms[j] = l_m
        l_ds[j] = l_d
        q_BLms[j] = q_BLm
        w_ms[j] = w_m
        T_BLcs[j] = T_BLc
        theta_BLcs[j] = T_BLc*(1+0.61*q_BLm)
        p_BLcs[j] = p(l_m, l_d)
        omega_ms[j] = omega_m
        Ps[j] = (P*3600*24*1000)/rho_w
        RHs[j] = q_BLm/wsat(T_BLc, p_BLcs[j])
        delhs[j] = dels_tropc + L_v*delqv_trop
        delhseff[j] = delhs[j] + ((L_v*alpha)/(1-alpha))*delqv_trop
        p_LCLs[j] = p_LCL
        
    
    toplot = np.abs(waterbalance) < 1
    
    l_ds = l_ds[toplot]
    waterbalance = waterbalance[toplot]
    massbalance = massbalance[toplot]
    Ebalance = Ebalance[toplot]
    l_ms = l_ms[toplot]
    w_ms = w_ms[toplot]
    Ps = Ps[toplot]
    RHs = RHs[toplot]
    delhs = delhs[toplot]
    delhseff = delhseff[toplot]
    q_BLms = q_BLms[toplot]
    omega_ms = omega_ms[toplot]
    T_BLcs = T_BLcs[toplot]
    p_BLcs = p_BLcs[toplot]
    theta_BLcs = theta_BLcs[toplot]
    
    delhs = -delhs
    delhseff = -delhseff
    
    negGMS = delhs <= 0
    #negGMS = delhs >= 0
    
    
    fmt = '-' if alpha > 0 else '--'
    #cplot = str((v+1.0)/(1.0*len(alphas))) if alpha > 0 else 'k' 
    
    plt.figure(1)
    plt.plot(l_ds/1e3, massbalance, fmt, color=cplot, label=r'$\alpha$ = {:2.2f}'.format(alpha))
    plt.xlabel(r'$l_d$ (km)')
    plt.ylabel(r'mass down - mass up (Pa m$^{2}$ s$^{-1}$)')
    #plt.title(r'$Q_{{d,net,trop}}$ = {:2.1f} K/day, $Q_{{c,net,trop}}$ = {:2.3f}'.format(dheatrate, mheatrate))
    if v == end-1:
        plt.legend()
        plt.xlim(0, domsizes[-1])
        plt.savefig(fout + 'massbalance_varyalpha_Qm{:2.3}.pdf'.format(mheatrate))
        plt.close()
    else:
        plt.savefig(fout + 'massbalance_varyalpha_Qm{:2.3}.pdf'.format(mheatrate))

    
    plt.figure(2)
    plt.plot(l_ds/1e3, waterbalance, fmt, color=cplot, label=r'$\alpha$ = {:2.2f}'.format(alpha))
    plt.xlabel(r'$l_d$ (km)')
    plt.ylabel(r'E - P (mm/day)')
    #plt.title(r'$Q_{{d,net,trop}}$ = {:2.1f} K/day, $Q_{{c,net,trop}}$ = {:2.3f} K/day'.format(dheatrate, mheatrate))
    if v == end-1:
        plt.legend()
        plt.xlim(0, domsizes[-1])
        plt.savefig(fout + 'waterbalance_varyalpha_Qm{:2.3}.pdf'.format(mheatrate))
        plt.close()
    else:
        plt.savefig(fout + 'waterbalance_varyalpha_Qm{:2.3}.pdf'.format(mheatrate))

    
    plt.figure(3)
    plt.plot(l_ds/1e3, Ebalance, fmt,  color=cplot, label=r'$\alpha$ = {:2.2f}'.format(alpha))
    plt.xlabel(r'$l_d$ (km)')
    plt.ylabel(r'condensation warming - radiative cooling (J s$^{-1}$)')
    #plt.title(r'$Q_{{d,net,trop}}$ = {:2.1f} K/day, $Q_{{c,net,trop}}$ = {:2.3f} K/day'.format(dheatrate, mheatrate))
    if v == end-1:
        plt.legend(loc='best')
        plt.xlim(0, domsizes[-1])
        plt.savefig(fout + 'Ebalance_varyalpha_Qm{:2.3}.pdf'.format(mheatrate))
        plt.close()
    else:
        plt.savefig(fout + 'Ebalance_varyalpha_Qm{:2.3}.pdf'.format(mheatrate))
    
    domsize_sim = np.array([768, 1536, 3072])    
    sigma_sim = np.array([0.051, 0.062, 0.082])
    wc_sim = np.array([0.069, 0.054, 0.038])
    omegac_sim = -rho_BL*g*wc_sim
    q_BLc_sim = np.array([17.0, 17.6, 18.1])
    RH_c_sim = np.array([89.3, 90.2, 90.3])
    TBL_c_sim = np.array([0.26, 0.46, 0.53])
    l_c_sim = np.sqrt((sigma_sim*domsize_sim**2)/np.pi)
    delh_sim = np.array([6581, 6338, 5452])
    P_sim = np.array([81.2, 58.9, 41.2])
    alpha_sim = np.array([0.05, 0.07, 0.09])
    colors_sim = ['k', 'r', 'g']
    delheff_sim = delh_sim + (-L_v*q_BLc_sim*1e-3)*(alpha_sim)/(1-alpha_sim)
    T_BLc_sim = np.array([0.26, 0.46, 0.53])

    
    plt.figure(4)
    plt.plot(l_ds/1e3, l_ms/1e3, fmt, color=cplot,  label=r'$\alpha$ = {:2.2f}'.format(alpha))
    #plt.plot(l_ds[negGMS]/1e3, l_ms[negGMS]/1e3, fmt, color='b', label=r'$\Delta h_{trop} < 0$'if v == end-2 else '')
    plt.plot(l_ds[negGMS]/1e3, l_ms[negGMS]/1e3, fmt, color='b')
    #plt.scatter(domsize_sim*(1./np.sqrt(2)), l_c_sim, 60, marker='x', c=colors_sim)
    if alpha > 0:
        plt.plot(l_ds[negGMS][0]/1e3, l_ms[negGMS][0]/1e3, '.', color='b', markersize=10, mew=3)
    plt.xlabel(r'$l_d$ (km)', fontsize=34)
    plt.ylabel(r'$l_c$ (km)', fontsize=34)
    #plt.title(r'$Q_{{d,net,trop}}$ = {:2.1f} K/day, $Q_{{c,net,trop}}$ = {:2.3f} K/day'.format(dheatrate, mheatrate))
    if v == end-1:
        plt.plot(domsize_sim[0]*(1./np.sqrt(2)), l_c_sim[0], 'x', color=colors_sim[0], markersize=10, mew=3, label='{:d} km'.format(domsize_sim[0]))
        plt.plot(domsize_sim[1]*(1./np.sqrt(2)), l_c_sim[1], 'x', color=colors_sim[1], markersize=10, mew=3, label='{:d} km'.format(domsize_sim[1]))
        plt.plot(domsize_sim[2]*(1./np.sqrt(2)), l_c_sim[2], 'x', color=colors_sim[2], markersize=10, mew=3, label='{:d} km'.format(domsize_sim[2]))
        plt.legend(loc='best')
        plt.xlim(0, domsizes[-1])
        plt.savefig(fout + 'lm_varyalpha_Qm{:2.3}new.pdf'.format(mheatrate))
        plt.close()
    else:
        plt.savefig(fout + 'lm_varyalpha_Qm{:2.3}new.pdf'.format(mheatrate))
        
    
    sigma = l_ms**2/l_ds**2

    
    plt.figure(5)
    plt.plot(l_ds/1e3, sigma, fmt,  color=cplot, label=r'$\alpha$ = {:2.2f}'.format(alpha))
    #plt.plot(l_ds[negGMS]/1e3, sigma[negGMS], fmt,  color='b', label=r'$\Delta h_{trop} < 0$'if v == end - 1 else '')
    plt.plot(l_ds[negGMS]/1e3, sigma[negGMS], fmt,  color='b')
    #plt.scatter(domsize_sim*(1./np.sqrt(2)), sigma_sim, 60, marker='x', c=colors_sim)
    if alpha > 0:
        plt.plot(l_ds[negGMS][0]/1e3, sigma[negGMS][0], '.', color ='b', markersize = 20, mew=3)
    plt.xlabel(r'$l_d$ (km)', fontsize=34)
    plt.ylabel(r'$\sigma$', fontsize=34)
    #plt.title(r'$Q_{{d,net,trop}}$ = {:2.1f} K/day, $Q_{{c,net,trop}}$ = {:2.3f} K/day'.format(dheatrate, mheatrate))
    if v == end-1:
        plt.plot(domsize_sim[0]*(1./np.sqrt(2)), sigma_sim[0], 'x', color=colors_sim[0], markersize=10, mew=3, label='{:d} km'.format(domsize_sim[0]))
        plt.plot(domsize_sim[1]*(1./np.sqrt(2)), sigma_sim[1], 'x', color=colors_sim[1], markersize=10, mew=3, label='{:d} km'.format(domsize_sim[1]))
        plt.plot(domsize_sim[2]*(1./np.sqrt(2)), sigma_sim[2], 'x', color=colors_sim[2], markersize=10, mew=3, label='{:d} km'.format(domsize_sim[2]))
        plt.legend(loc='best')
        plt.xlim(0, domsizes[-1])
        plt.savefig(fout + 'sigma_varyalpha_Qm{:2.3}new.pdf'.format(mheatrate))
        plt.close()
    else:
        plt.savefig(fout + 'sigma_varyalpha_Qm{:2.3}new.pdf'.format(mheatrate))

    
    plt.figure(6)
    plt.plot(l_ds/1e3, w_ms, fmt,  color=cplot, label=r'$\alpha$ = {:2.2f}'.format(alpha))
    #plt.plot(l_ds[negGMS]/1e3, w_ms[negGMS], fmt,  color='b', label=r'$\Delta h_{trop} < 0$'if v == end - 1 else '')
    plt.plot(l_ds[negGMS]/1e3, w_ms[negGMS], fmt,  color='b')
    #plt.scatter(domsize_sim*(1./np.sqrt(2)), wc_sim, 60, marker='x', c=colors_sim)
    if alpha > 0:
        plt.plot(l_ds[negGMS][0]/1e3, w_ms[negGMS][0], '.', color ='b', markersize = 20, mew=3)
    plt.xlabel(r'$l_d$ (km)', fontsize=34)
    plt.ylabel(r'$w_c$ (m/s)', fontsize=34)
    #plt.title(r'$Q_{{d,net,trop}}$ = {:2.1f} K/day, $Q_{{c,net,trop}}$ = {:2.3f} K/day'.format(dheatrate, mheatrate))
    if v == end-1:
        plt.plot(domsize_sim[0]*(1./np.sqrt(2)), wc_sim[0], 'x', color=colors_sim[0], markersize=10, mew=3, label='{:d} km'.format(domsize_sim[0]))
        plt.plot(domsize_sim[1]*(1./np.sqrt(2)), wc_sim[1], 'x', color=colors_sim[1], markersize=10, mew=3, label='{:d} km'.format(domsize_sim[1]))
        plt.plot(domsize_sim[2]*(1./np.sqrt(2)), wc_sim[2], 'x', color=colors_sim[2], markersize=10, mew=3, label='{:d} km'.format(domsize_sim[2]))
        plt.legend(loc='best')
        plt.xlim(0, domsizes[-1])
        plt.savefig(fout + 'wm_varyalpha_Qm{:2.3}new.pdf'.format(mheatrate))
        plt.close()
    else:
        plt.savefig(fout + 'wm_varyalpha_Qm{:2.3}new.pdf'.format(mheatrate))
        
    plt.figure(7)
    plt.plot(l_ds/1e3, T_BLcs - T_s, fmt,  color=cplot, label=r'$\alpha$ = {:2.2f}'.format(alpha))
    #plt.plot(l_ds[negGMS]/1e3, w_ms[negGMS], fmt,  color='b', label=r'$\Delta h_{trop} < 0$'if v == end - 1 else '')
    plt.plot(l_ds[negGMS]/1e3, T_BLcs[negGMS] - T_s, fmt,  color='b')
    #plt.scatter(domsize_sim*(1./np.sqrt(2)), wc_sim, 60, marker='x', c=colors_sim)
    if alpha > 0:
        plt.plot(l_ds[negGMS][0]/1e3, T_BLcs[negGMS][0] - T_s, '.', color ='b', markersize = 20, mew=3)
    plt.xlabel(r'$l_d$ (km)', fontsize=34)
    plt.ylabel(r'$\Delta T_{BL,c}$ (K)', fontsize=34)
    #plt.title(r'$Q_{{d,net,trop}}$ = {:2.1f} K/day, $Q_{{c,net,trop}}$ = {:2.3f} K/day'.format(dheatrate, mheatrate))
    if v == end-1:
        plt.plot(domsize_sim[0]*(1./np.sqrt(2)), TBL_c_sim[0], 'x', color=colors_sim[0], markersize=10, mew=3, label='{:d} km'.format(domsize_sim[0]))
        plt.plot(domsize_sim[1]*(1./np.sqrt(2)), TBL_c_sim[1], 'x', color=colors_sim[1], markersize=10, mew=3, label='{:d} km'.format(domsize_sim[1]))
        plt.plot(domsize_sim[2]*(1./np.sqrt(2)), TBL_c_sim[2], 'x', color=colors_sim[2], markersize=10, mew=3, label='{:d} km'.format(domsize_sim[2]))
        plt.legend(loc='best')
        plt.xlim(0, domsizes[-1])
        plt.savefig(fout + 'TBLc_varyalpha{:2.3}new.pdf'.format(mheatrate))
        plt.close()
    else:
        plt.savefig(fout + 'TBLc_varyalpha{:2.3}new.pdf'.format(mheatrate))
        
    plt.figure(8)
    plt.plot(l_ds/1e3, p_BLcs - p_s, fmt,  color=cplot, label=r'$\alpha$ = {:2.2f}'.format(alpha))
    #plt.plot(l_ds[negGMS]/1e3, w_ms[negGMS], fmt,  color='b', label=r'$\Delta h_{trop} < 0$'if v == end - 1 else '')
    plt.plot(l_ds[negGMS]/1e3, p_BLcs[negGMS] - p_s, fmt,  color='b')
    #plt.scatter(domsize_sim*(1./np.sqrt(2)), wc_sim, 60, marker='x', c=colors_sim)
    if alpha > 0:
        plt.plot(l_ds[negGMS][0]/1e3, p_BLcs[negGMS][0] - p_s, '.', color ='b', markersize = 20, mew=3)
    plt.xlabel(r'$l_d$ (km)', fontsize=34)
    plt.ylabel(r'$\Delta p_{BL,c}$ (Pa)', fontsize=34)
    #plt.title(r'$Q_{{d,net,trop}}$ = {:2.1f} K/day, $Q_{{c,net,trop}}$ = {:2.3f} K/day'.format(dheatrate, mheatrate))
    if v == end-1:
        #plt.plot(domsize_sim[0]*(1./np.sqrt(2)), wc_sim[0], 'x', color=colors_sim[0], markersize=10, mew=3, label='{:d} km'.format(domsize_sim[0]))
        #plt.plot(domsize_sim[1]*(1./np.sqrt(2)), wc_sim[1], 'x', color=colors_sim[1], markersize=10, mew=3, label='{:d} km'.format(domsize_sim[1]))
        #plt.plot(domsize_sim[2]*(1./np.sqrt(2)), wc_sim[2], 'x', color=colors_sim[2], markersize=10, mew=3, label='{:d} km'.format(domsize_sim[2]))
        #plt.legend(loc='best')
        plt.xlim(0, domsizes[-1])
        plt.savefig(fout + 'pBLc_varyalpha{:2.3}new.pdf'.format(mheatrate))
        plt.close()
    else:
        plt.savefig(fout + 'pBLc_varyalpha{:2.3}new.pdf'.format(mheatrate))
        
        


    fig=plt.figure(9)
    ax1= fig.gca()
    #if flag:
    #    ax1 = fig.gca()
    #    ax2 = ax1.twinx()
    #    flag = False
    h1, = ax1.plot(l_ds/1e3, q_BLms*1e3, fmt,  color=cplot) 
    #def update_ax2(ax1): 
    #    y1, y2 = ax1.get_ylim() 
    #    ax2.set_ylim(0.1*(y1/q_sat), 0.1*(y2/q_sat)) 
    #ax1.callbacks.connect("ylim_changed", update_ax2) 
    h2, =ax1.plot(l_ds[negGMS]/1e3, q_BLms[negGMS]*1e3, fmt,  color='b')
    #plt.scatter(domsize_sim*(1./np.sqrt(2)), q_BLc_sim, 60, marker='x', c=colors_sim)
    if alpha > 0:
        ax1.plot(l_ds[negGMS][0]/1e3, q_BLms[negGMS][0]*1e3, '.', color='b', markersize = 20, mew=3)
    ax1.set_xlabel(r'$l_d$ (km)', fontsize=34)
    #ax1.set_ylim(20, 24)
    #ax2.set_ylim(80, 100)
    ax1.set_ylabel(r'$q_{BL,c}$ (g/kg)', fontsize=34)
    #ax2.set_ylabel(r'$RH_c$ (%)')
    #plt.legend(loc='best')
    #plt.title(r'$Q_{{d,net,trop}}$ = {:2.1f} K/day, $Q_{{c,net,trop}}$ = {:2.3f} K/day'.format(dheatrate, mheatrate))
    label = r'$\alpha$ = {:2.2f}'.format(alpha)
    labels = labels + [label]
    hs = hs + [h1]
    if v == end-1:
        #label=r'$\Delta h_{trop} < 0$'
        #hs = hs + [h2]
        #labels = labels + [label]
        label='{:d} km'.format(domsize_sim[0])
        labels = labels + [label]
        label='{:d} km'.format(domsize_sim[1])
        labels = labels + [label]
        label='{:d} km'.format(domsize_sim[2])
        labels = labels + [label]
        h1, = plt.plot(domsize_sim[0]*(1./np.sqrt(2)), q_BLc_sim[0], 'x', color=colors_sim[0], markersize=10, mew=3)
        h2, = plt.plot(domsize_sim[1]*(1./np.sqrt(2)), q_BLc_sim[1], 'x', color=colors_sim[1], markersize=10, mew=3)
        h3, = plt.plot(domsize_sim[2]*(1./np.sqrt(2)), q_BLc_sim[2], 'x', color=colors_sim[2], markersize=10, mew=3)
        hs = hs + [h1, h2, h3]
        ax1.legend(hs, labels, loc=0)
        plt.xlim(0, domsizes[-1])
        plt.savefig(fout + 'qBLm_varyalpha_Qm{:2.3}new.pdf'.format(mheatrate))
        plt.close()
    else:
        plt.savefig(fout + 'qBLm_varyalpha_Qm{:2.3}new.pdf'.format(mheatrate))

    plt.figure(10)
    plt.plot(l_ds/1e3, Ps, fmt,  color=cplot, label=r'$\alpha$ = {:2.2f}'.format(alpha))
    #plt.plot(l_ds[negGMS]/1e3, Ps[negGMS], fmt,  color='b', label=r'$\Delta h_{trop} < 0$'if v == end - 1 else '')
    plt.plot(l_ds[negGMS]/1e3, Ps[negGMS], fmt,  color='b')
    if alpha > 0:
        plt.plot(l_ds[negGMS][0]/1e3, Ps[negGMS][0], '.', color ='b', markersize = 20, mew=3)
    plt.xlabel(r'$l_d$ (km)', fontsize=34)
    plt.ylabel(r'$P$ (mm/day)', fontsize=34)
    #plt.title(r'$Q_{{d,net,trop}}$ = {:2.1f} K/day, $Q_{{c,net,trop}}$ = {:2.3f} K/day'.format(dheatrate, mheatrate))
    if v == end-1:
        #label=r'$\Delta h_{trop} < 0$'
        #hs = hs + [h2]
        #labels = labels + [label]
        label='{:d} km'.format(domsize_sim[0])
        labels = labels + [label]
        label='{:d} km'.format(domsize_sim[1])
        labels = labels + [label]
        label='{:d} km'.format(domsize_sim[2])
        labels = labels + [label]
        h1, = plt.plot(domsize_sim[0]*(1./np.sqrt(2)), P_sim[0], 'x', color=colors_sim[0], markersize=10, mew=3)
        h2, = plt.plot(domsize_sim[1]*(1./np.sqrt(2)), P_sim[1], 'x', color=colors_sim[1], markersize=10, mew=3)
        h3, = plt.plot(domsize_sim[2]*(1./np.sqrt(2)), P_sim[2], 'x', color=colors_sim[2], markersize=10, mew=3)
        hs = hs + [h1, h2, h3]
        ax1.legend(hs, labels, loc=0)
        plt.xlim(0, domsizes[-1])
        plt.savefig(fout + 'P_varyalpha_Qm{:2.3}new.pdf'.format(mheatrate))
        plt.close()

    else:
        plt.savefig(fout + 'P_varyalpha_Qm{:2.3}new.pdf'.format(mheatrate))

    
    plt.figure(11)
    plt.plot(l_ds/1e3, delhs, fmt,  color=cplot, label=r'$\alpha$ = {:2.2f}'.format(alpha))
    #plt.plot(l_ds[negGMS]/1e3, delhs[negGMS], fmt,  color='b', label=r'$\Delta h_{trop} < 0$'if v == end - 1 else '')
    plt.plot(l_ds[negGMS]/1e3, delhs[negGMS], fmt,  color='b')
    #plt.scatter(domsize_sim*(1./np.sqrt(2)), delh_sim, 60, marker='x', c=colors_sim)
    if alpha > 0:
        plt.plot(l_ds[negGMS][0]/1e3, delhs[negGMS][0], '.', color ='b', markersize = 20, mew=3)
    plt.xlabel(r'$l_d$ (km)', fontsize=34)
    plt.ylabel(r'$\Delta h_{trop}$ (J)', fontsize=34)
    #plt.title(r'$Q_{{d,net,trop}}$ = {:2.1f} K/day, $Q_{{c,net,trop}}$ = {:2.3f} K/day'.format(dheatrate, mheatrate))
    if v == end-1:
        plt.plot(domsize_sim[0]*(1./np.sqrt(2)), delh_sim[0], 'x', color=colors_sim[0], markersize=10, mew=3, label='{:d} km'.format(domsize_sim[0]))
        plt.plot(domsize_sim[1]*(1./np.sqrt(2)), delh_sim[1], 'x', color=colors_sim[1], markersize=10, mew=3, label='{:d} km'.format(domsize_sim[1]))
        plt.plot(domsize_sim[2]*(1./np.sqrt(2)), delh_sim[2], 'x', color=colors_sim[2], markersize=10, mew=3, label='{:d} km'.format(domsize_sim[2]))
        plt.legend(loc='best')
        plt.axhline(0, color='k', linewidth=1)
        plt.xlim(0, domsizes[-1])
        plt.savefig(fout + 'deltah_varyalpha_Qm{:2.3}new.pdf'.format(mheatrate))
        plt.close()
    else:
        plt.savefig(fout + 'deltah_varyalpha_Qm{:2.3}new.pdf'.format(mheatrate))
        
    plt.figure(12)
    plt.plot(l_ds/1e3, delhseff, fmt,  color=cplot, label=r'$\alpha$ = {:2.2f}'.format(alpha))
    #plt.plot(l_ds[negGMS]/1e3, delhseff[negGMS], fmt,  color='b', label=r'$\Delta h_{trop} < 0$'if v == end - 1 else '')
    plt.plot(l_ds[negGMS]/1e3, delhseff[negGMS], fmt,  color='b')
    if alpha > 0:
        plt.plot(l_ds[negGMS][0]/1e3, delhseff[negGMS][0], '.', color ='b', markersize = 20, mew=3)
    plt.xlabel(r'$l_d$ (km)', fontsize=34)
    plt.ylabel(r'$\Delta h_{trop,eff}$ (J)', fontsize=34)
    #plt.title(r'$Q_{{d,net,trop}}$ = {:2.1f} K/day, $Q_{{c,net,trop}}$ = {:2.3f} K/day'.format(dheatrate, mheatrate))
    if v == end-1:
        plt.plot(domsize_sim[0]*(1./np.sqrt(2)), delheff_sim[0], 'x', color=colors_sim[0], markersize=10, mew=3, label='{:d} km'.format(domsize_sim[0]))
        plt.plot(domsize_sim[1]*(1./np.sqrt(2)), delheff_sim[1], 'x', color=colors_sim[1], markersize=10, mew=3, label='{:d} km'.format(domsize_sim[1]))
        plt.plot(domsize_sim[2]*(1./np.sqrt(2)), delheff_sim[2], 'x', color=colors_sim[2], markersize=10, mew=3, label='{:d} km'.format(domsize_sim[2]))
        plt.legend(loc='best')
        plt.axhline(0, color='k', linewidth=1)
        plt.xlim(0, domsizes[-1])
        plt.savefig(fout + 'deltaheff_varyalpha_Qm{:2.3}new.pdf'.format(mheatrate))
        plt.close()

    else:
        plt.savefig(fout + 'deltaheff_varyalpha_Qm{:2.3}new.pdf'.format(mheatrate))
        
    plt.figure(13)
    plt.plot(l_ds/1e3, RHs*100, fmt,  color=cplot, label=r'$\alpha$ = {:2.2f}'.format(alpha))
    #plt.plot(l_ds[negGMS]/1e3, RHs[negGMS]*100, fmt,  color='b', label=r'$\Delta h_{trop} < 0$'if v == end - 1 else '')
    plt.plot(l_ds[negGMS]/1e3, RHs[negGMS]*100, fmt,  color='b')
    #plt.scatter(domsize_sim*(1./np.sqrt(2)), RH_c_sim, 60, marker='x', c=colors_sim)
    if alpha > 0:
        plt.plot(l_ds[negGMS][0]/1e3, RHs[negGMS][0]*100, '.', color ='b', markersize = 20, mew=3)
    plt.xlabel(r'$l_d$ (km)', fontsize=34)
    plt.ylabel(r'$RH_c$ (%)', fontsize=34)
    #plt.title(r'$Q_{{d,net,trop}}$ = {:2.1f} K/day, $Q_{{c,net,trop}}$ = {:2.3f} K/day'.format(dheatrate, mheatrate))
    if v == end-1:
        plt.plot(domsize_sim[0]*(1./np.sqrt(2)), RH_c_sim[0], 'x', color=colors_sim[0], markersize=10, mew=3, label='{:d} km'.format(domsize_sim[0]))
        plt.plot(domsize_sim[1]*(1./np.sqrt(2)), RH_c_sim[1], 'x', color=colors_sim[1], markersize=10, mew=3, label='{:d} km'.format(domsize_sim[1]))
        plt.plot(domsize_sim[2]*(1./np.sqrt(2)), RH_c_sim[2], 'x', color=colors_sim[2], markersize=10, mew=3, label='{:d} km'.format(domsize_sim[2]))
        plt.legend(loc='best')
        plt.xlim(0, domsizes[-1])
        plt.savefig(fout + 'RH_varyalpha_Qm{:2.3}new.pdf'.format(mheatrate))
        plt.close()
        
    else:
        plt.savefig(fout + 'RH_varyalpha_Qm{:2.3}new.pdf'.format(mheatrate))
        
    GMS = w_ms*delhs
        
    plt.figure(14)
    plt.plot(l_ds/1e3, GMS, fmt,  color=cplot, label=r'$\alpha$ = {:2.2f}'.format(alpha))
    #plt.plot(l_ds[negGMS]/1e3, delhs[negGMS], fmt,  color='b', label=r'$\Delta h_{trop} < 0$'if v == end - 1 else '')
    plt.plot(l_ds[negGMS]/1e3, GMS[negGMS], fmt,  color='b')
    #plt.scatter(domsize_sim*(1./np.sqrt(2)), delh_sim, 60, marker='x', c=colors_sim)
    if alpha > 0:
        plt.plot(l_ds[negGMS][0]/1e3, GMS[negGMS][0], '.', color ='b', markersize = 20, mew=3)
    plt.xlabel(r'$l_d$ (km)', fontsize=34)
    plt.ylabel(r'$M$ (W/m$^2$)', fontsize=34)
    #plt.title(r'$Q_{{d,net,trop}}$ = {:2.1f} K/day, $Q_{{c,net,trop}}$ = {:2.3f} K/day'.format(dheatrate, mheatrate))
    if v == end-1:
        # plt.plot(domsize_sim[0]*(1./np.sqrt(2)), delh_sim[0], 'x', color=colors_sim[0], markersize=10, mew=3, label='{:d} km'.format(domsize_sim[0]))
        # plt.plot(domsize_sim[1]*(1./np.sqrt(2)), delh_sim[1], 'x', color=colors_sim[1], markersize=10, mew=3, label='{:d} km'.format(domsize_sim[1]))
        # plt.plot(domsize_sim[2]*(1./np.sqrt(2)), delh_sim[2], 'x', color=colors_sim[2], markersize=10, mew=3, label='{:d} km'.format(domsize_sim[2]))
        plt.legend(loc='best')
        plt.axhline(0, color='k', linewidth=1)
        plt.xlim(0, domsizes[-1])
        plt.savefig(fout + 'GMS_varyalpha_Qm{:2.3}new.pdf'.format(mheatrate))
        plt.close()
    else:
        plt.savefig(fout + 'GMS_varyalpha_Qm{:2.3}new.pdf'.format(mheatrate))
        
    plt.figure(15)
    plt.plot(l_ds/1e3, theta_BLcs - T_s, fmt,  color=cplot, label=r'$\alpha$ = {:2.2f}'.format(alpha))
    #plt.plot(l_ds[negGMS]/1e3, delhs[negGMS], fmt,  color='b', label=r'$\Delta h_{trop} < 0$'if v == end - 1 else '')
    plt.plot(l_ds[negGMS]/1e3, theta_BLcs[negGMS] - T_s, fmt,  color='b')
    #plt.scatter(domsize_sim*(1./np.sqrt(2)), delh_sim, 60, marker='x', c=colors_sim)
    if alpha > 0:
        plt.plot(l_ds[negGMS][0]/1e3, theta_BLcs[negGMS][0] - T_s, '.', color ='b', markersize = 20, mew=3)
    plt.xlabel(r'$l_d$ (km)', fontsize=34)
    plt.ylabel(r'$\Delta \theta_{BL,c}$ (K)', fontsize=34)
    #plt.title(r'$Q_{{d,net,trop}}$ = {:2.1f} K/day, $Q_{{c,net,trop}}$ = {:2.3f} K/day'.format(dheatrate, mheatrate))
    if v == end-1:
        # plt.plot(domsize_sim[0]*(1./np.sqrt(2)), delh_sim[0], 'x', color=colors_sim[0], markersize=10, mew=3, label='{:d} km'.format(domsize_sim[0]))
        # plt.plot(domsize_sim[1]*(1./np.sqrt(2)), delh_sim[1], 'x', color=colors_sim[1], markersize=10, mew=3, label='{:d} km'.format(domsize_sim[1]))
        # plt.plot(domsize_sim[2]*(1./np.sqrt(2)), delh_sim[2], 'x', color=colors_sim[2], markersize=10, mew=3, label='{:d} km'.format(domsize_sim[2]))
        plt.legend(loc='best')
        #plt.axhline(0, color='b', alpha=0.9)
        plt.xlim(0, domsizes[-1])
        plt.savefig(fout + 'thetaBLc_varyalpha_Qm{:2.3}new.pdf'.format(mheatrate))
        plt.close()
    else:
        plt.savefig(fout + 'thetaBLc_varyalpha_Qm{:2.3}new.pdf'.format(mheatrate))
        

r_test = np.linspace(600e3,2000e3,100)   

ps = np.zeros(r_test.shape)
Ts = np.zeros(r_test.shape)

for j,r in enumerate(r_test):
    ps[j] = p(r, 2000e3)
    Ts[j] = theta_v(r, 2000e3)
    
    
fig=plt.figure(figsize=(16,14))
#fig.tight_layout()
plt.plot(r_test/1e3, (ps-p_s))
plt.xlabel('r (km)')
plt.ylabel('$\Delta p$ (Pa)')
plt.savefig(fout + 'deltap_vs_r_new.png')

fig=plt.figure(figsize=(16,14))
#fig.tight_layout()
plt.plot(r_test/1e3, Ts-T_s)
plt.xlabel('r (km)')
plt.ylabel('$\Delta T$ (K)')
plt.savefig(fout + 'deltaT_vs_r_new.png')        


    
    # plt.figure(11)
    # plt.plot(l_ds/1e3, p_LCLs/1e2, '.-',  color=cplot,  label=r'$\alpha$ = {:2.3f}'.format(alpha))
    # plt.xlabel(r'$l_d$ (km)')
    # plt.ylabel(r'$p_{LCL}$ (hPa)')
    # plt.title(r'$Q_{{d,trop}}$ = {:2.1f} K/day, $Q_{{c,net,trop}}$ = {:2.3f} K/day'.format(dheatrate, mheatrate))
    # if v == end-1:
    #     plt.legend(loc='best')
    #     plt.savefig(fout + 'pLCL_varyalpha_Qm{:2.3}.pdf'.format(mheatrate))
    #     plt.close()
    # else:
    #     plt.savefig(fout + 'pLCL_varyalpha_Qm{:2.3}.pdf'.format(mheatrate))

    


#plt.figure(10)
#plt.plot(l_ds/1e3, w_BLpluss, 'o-', label=r'$w_{BL+}$')
#plt.plot(l_ds/1e3, w_BLminuss, 'o-', label=r'$w_{BL-}$')
#plt.xlabel(r'l$_d$ (km)')
#plt.ylabel(r'$w$ (m/s)')
#plt.legend()
#plt.savefig(fout + 'BLcont_varyld.pdf')
#plt.close()
        
        
        
        




        

