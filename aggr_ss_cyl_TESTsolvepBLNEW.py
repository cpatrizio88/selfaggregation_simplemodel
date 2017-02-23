import site
import sys
#site.addsitedir('\\Users\\Casey\\Dropbox\\research\\code\\thermlib')
site.addsitedir('/Users/cpatrizio/repos/thermolib/')
from wsat import wsat
#from findTmoist import findTmoist
from findLCL0 import findLCL0
from constants import constants as c
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import numpy as np
import matplotlib
import thermo
from constants import constants
import scipy.optimize
from scipy.optimize import fsolve, brentq, newton
import matplotlib.cm as cm


fout = '/Users/cpatrizio/Google Drive/figures/simplemodel/'

p_s = 1000e2 #surface pressure (Pa)
#p_BL = 950e2 #boundary layer top (Pa)
#p_out = 225e2 #tropopause (Pa)
#p_t = 150e2 #tropopause (Pa)

h=7. #scale height (km)
gamma_PH = 6.2 #lapse rate (K km^-1)
zeta_T = 15.0 #km

def findT(T_s, p):
    
    zeta = -h*np.log(p/p_s)
    n = zeta.size
    T = np.zeros(n)
    if n > 1:
        T[zeta < zeta_T] = T_s - gamma_PH*zeta[zeta < zeta_T]
        T[zeta >= zeta_T] = T_s - gamma_PH*zeta_T
    else:
       if (zeta < zeta_T):
          T = T_s - gamma_PH*zeta
       else:
          T = T_s - gamma_PH*zeta_T
    return T
    
    
    
p_t = p_s*np.exp(-zeta_T/h)


g = 9.81
sig = 5.67e-8
T_s = 302 
#T_BL = 302
L_v = 2.257*1e6

#eps_a = 1
eps_strat = 0.1
#eps_BL = 0.6

#S0 = 413.98 #solar constant
#theta = 50.5 #zenith angle 
#S = np.cos((theta*np.pi/180))*S0 

c_E = .001
rho_w = 1000 #density of water


d = 100000e3 #

delp_out = 50*1e2

p_s = 1000e2 #surface pressure (Pa)
#p_BL = 950e2 #boundary layer top (Pa)
p_out = p_t + delp_out #tropopause (Pa)
#p_t = 150e2 #tropopause (Pa)

domsizes = np.linspace(500, 12000, 40)
#domsizes = []
colors = ['k', 'r', 'g']

#w_ms = np.zeros(len(domsizes))
#l_ms = np.zeros(len(domsizes))
#l_ds = np.zeros(len(domsizes))
#massbalance = np.zeros(len(domsizes))
#waterbalance = np.zeros(len(domsizes))
#Ebalance = np.zeros(len(domsizes))
#RHs = np.zeros(len(domsizes))
#Ps = np.zeros(len(domsizes))
#delhs = np.zeros(len(domsizes))
#p_LCLs = np.zeros(len(domsizes))

eps_BLs = np.linspace(0, 1, 20)
eps_as = np.linspace(0, 1, 20)


eps_BLss, eps_ass = np.meshgrid(eps_BLs, eps_as)

p_BLs = np.zeros(eps_BLss.shape)
#T_as = np.zeros(eps_BLss.shape)

q_sat = wsat(T_s, p_s) #mixing ratio above sea surface (100% saturated)    
thetae0 = thermo.theta_e(T_s, p_s, q_sat, 0) #theta_e in moist region 
                                            #use surface temperature to get moist adiabat

plevs = np.linspace(50e2, p_s, 1000)    
#Tadb = findTmoist(thetae0, plevs)
omega_BLs = np.zeros(eps_BLss.shape)
w_BLs = np.zeros(eps_BLss.shape)
w_ms = np.zeros(eps_BLss.shape)
l_ms = np.zeros(eps_BLss.shape)
T_BLs = np.zeros(eps_BLss.shape)
netprecip = np.zeros(eps_BLss.shape)
RH_ms = np.zeros(eps_BLss.shape)
delhs = np.zeros(eps_BLss.shape)
massbalance = np.zeros(eps_BLss.shape)
p_LCLs = np.zeros(eps_BLss.shape)

l_d = 3000e3
                                                                                                                                                 

for i, eps_BL in enumerate(eps_BLs):
    
    print 'eps_BL', eps_BL
    
    for j, eps_a in enumerate(eps_as):
        
        print 'eps_a', eps_a

        #T_a = np.mean(Tadb)
        #T_out = findTmoist(thetae0, p_out)
        T_out = findT(T_s, p_out)
        T_strat = T_out
        #T_strat = np.mean(findTmoist(thetae0, np.linspace(95,p_out,50)))
        #rho_BL = p_s/(c.Rd*T_BL)
        s_surf = c.cpd*T_s #dry static energy at surface
        
                                                    
                                                    
        def T_BLfn(p_BL):
            #T_BLtop = findTmoist(thetae0, p_BL) #temperature of boundary layer top 
            T_BLtop = findT(T_s, p_BL)
            T_BL = np.add(T_s, T_BLtop)/2. #temperature of boundary layer, consistent with well-mixed assumption (linear mixing)
            return T_BL
            
        def T_afn(p_BL):
            T_a = np.mean(findT(T_s, np.linspace(p_out, p_BL, 50)))
            return T_a
            
        def dels_tropfn(p_BL):
            T_BLtop = findT(T_s, p_BL) #temperature of boundary layer top 
            thickness = lambda p: (c.Rd/(p*g))*findT(T_s, p)
            delz_trop = integrate.quad(thickness, p_out, p_BL)[0]
            #delz_trop = np.zeros(len(p_BL))
            #for i,p_BLe in enumerate(p_BL):
            #    delz_trop[i] = integrate.quad(thickness, p_out, p_BLe)[0]
            dels_trop = c.cpd*np.subtract(T_BLtop, T_out) + g*(-delz_trop) #difference in dry static energy between boundary layer top and tropopause (J)
            return dels_trop
            
        def dels_BLfn(p_BL): 
            T_BLtop = findT(T_s, p_BL) #temperature of boundary layer top 
            rho_BLtop = p_BL/np.multiply(c.Rd, T_BLtop)
            rho_s = p_s/(c.Rd*T_s)
            rho_BL = (rho_BLtop + rho_s)/2.
            delz_BL = np.subtract(p_s, p_BL)/(g*rho_BL)
            s_BLtop = c.cpd*np.array(T_BLtop) + g*np.array(delz_BL) #dry static energy at top of boundary layer
            s_BL = np.add(s_BLtop,s_surf)/2. #dry static energy of BL (well-mixed)
            dels_BL = np.subtract(s_BL, s_BLtop) #difference in dry static energy between BL and right above BL
            return dels_BL
            
        #Ta_fn = lambda T_a: eps_a*(eps_strat*sig*T_strat**4 + eps_BL*sig*T_BL**4 + (1-eps_BL)*sig*T_s**4 - 2*sig*T_a**4)/dels_tropfn(p_BL) - eps_BL*((1-eps_a)*eps_strat*sig*T_strat**4 + eps_a*sig*T_a**4 + sig*T_s**4 - 2*sig*T_BL**4)/dels_BLfn(p_BL)    
        pBL_fn = lambda p_BL: eps_a*(eps_strat*sig*T_strat**4 + eps_BL*sig*T_BLfn(p_BL)**4 + (1-eps_BL)*sig*T_s**4 - 2*sig*T_afn(p_BL)**4)/dels_tropfn(p_BL) - eps_BL*((1-eps_a)*eps_strat*sig*T_strat**4 + eps_a*sig*T_afn(p_BL)**4 + sig*T_s**4 - 2*sig*T_BLfn(p_BL)**4)/dels_BLfn(p_BL)
        #pBL_fn = lambda p_BL: eps_a*(eps_BL*sig*T_BL**4 + (1-eps_BL)*sig*T_s**4 - 2*sig*T_a**4)/dels_tropfn(p_BL) - eps_BL*( eps_a*sig*T_a**4 + sig*T_s**4 - 2*sig*T_BL**4)/dels_BLfn(p_BL)
        
        
        #p_BL = newton(pBL_fn, 950e2, maxiter=100000)
        p_BL = fsolve(pBL_fn, 950e2)[0]
        
        p_LCL = p_BL
        
        p_BLs[i,j] = p_BL
        p_LCLs[i,j] = p_LCL
        T_BL = T_BLfn(p_BL) 
        T_a = T_afn(p_BL)
        T_BLs[i,j] = T_BL
                
        T_BLtop = findT(T_s, p_BL) #temperature of boundary layer top 
        rho_BLtop = p_BL/(c.Rd*T_BLtop)
        
        rho_BLtop = p_BL/np.multiply(c.Rd, T_BLtop)
        rho_s = p_s/(c.Rd*T_s)
        rho_BL = (rho_BLtop + rho_s)/2.
 
        dels_BL = dels_BLfn(p_BL)
        dels_trop = dels_tropfn(p_BL)
            
        #dels_trop = c.cpd*(T_BLtop - T_out) + g*(-delz_trop) #difference in dry static energy between boundary layer top and tropopause (J)
        
        M_trop = (p_BL - p_out)/g #mass of troposphere in kg m^-2
        M_BL = (p_s - p_BL)/g #mass of boundary layer in kg m^-2

        #radiative heating in dry region interior 
        Q_dtrop = eps_a*(eps_strat*sig*T_strat**4 + eps_BL*sig*T_BL**4 + (1 - eps_BL)*sig*T_s**4 - 2*sig*T_a**4)

        mheatrate=-0.1 #prescribed heating rate in moist region interior K/day
        Q_mtrop = (c.cpd*M_trop*mheatrate)/(86400)
        
        #radiative warming rate of dry region boundary layer 
        Q_dBL = (eps_BL*sig**eps_a*T_a**4 + eps_BL*sig*T_s**4 - 2*eps_BL*sig*T_BL**4 + (1 - eps_a)*eps_strat*sig*T_strat**4)
        
        omega_BLs[i,j] = (g*Q_dtrop)/dels_trop
        omega_BL = omega_BLs[i,j]
        w_BLs[i,j] = omega_BLs[i,j]/(-rho_BLtop*g)
        
        q_FA = wsat(T_out, p_out)
        
                
        delp_BL = p_s - p_BL
        delz_BL = delp_BL/(rho_BL*g)
        zhat = delz_BL/c_E
        
        #u_BL = -(omega_BL/(2*delp_BL))*((l_d**2 - r**2))/r
        
        def equations(p):
            l_m, RH_m = p
            #p_LCL, T_LCL = findLCL0(f*q_sat, p_s, T_BL)
            #T_LCL = findTmoist(thetae0, p_LCL)
            #delz_LCL = integrate.quad(thickness, p_out, p_LCL)[0]
            #dels_LCL = c.cpd*(T_LCL-T_out) + g*(-delz_LCL)
            dels_LCL = dels_trop
            omega_m = (g*Q_mtrop)/(dels_LCL + L_v*(RH_m*q_sat - q_FA))
            
            f1 = omega_m*l_m**2 + omega_BL*(l_d**2 - l_m**2)
            f2 = (RH_m-1)*q_sat - (2*zhat*(q_FA - q_sat))/(l_d**2 - l_m**2)*(l_m + zhat - (zhat + l_d)*np.exp((l_m - l_d)/zhat))
            
            return (f1, f2)
            
        l_m, RH_m = fsolve(equations, [100e3, 0.5], xtol=1e-13, maxfev=100000)
            
        #p_LCL, T_LCL = findLCL0(f*q_sat, p_s, T_BL)
        #T_LCL = findTmoist(thetae0, p_LCL)
        #delz_LCL = integrate.quad(thickness, p_t, p_LCL)[0]
        #dels_LCL = c.cpd*(T_LCL-T_t) + g*(-delz_LCL)
        M_LCL = (p_LCL - p_t)/g  #mass of convective region (from LCL to tropopause) in kg m^-2
        
        #W_m = Q_mtrop/(c.cpd*M_LCL)
        #W_d = Q_dtrop/(c.cpd*M_trop)
        
        #l_m = l_mfn(f)
            
        #l_mi = np.where(r > l_m)[0][0]
        dels_LCL = dels_trop
        delh = dels_LCL + L_v*(RH_m*q_sat - q_FA)

            
        omega_m = (g*Q_mtrop)/delh
        w_m = -omega_m/(rho_BLtop*g)
        
        P = -(1./g)*omega_m*np.pi*l_m**2*(RH_m*q_sat - q_FA)
        
        E = (1./g)*(2*np.pi*omega_BL)*(q_sat - q_FA)*((l_d**2 - l_m**2)/2. + zhat*(l_d - l_m) - zhat*(zhat + l_d)*(1 - np.exp((l_m - l_d)/zhat)))
        
        qBL_m = RH_m*q_sat
        qBL_sat = wsat(T_BL, (p_s + p_BL)/2.)
        
        netprecip[i,j] = P - E
        l_ms[i,j] = l_m
        RH_ms[i,j] = qBL_m/qBL_sat
        delhs[i,j] = delh
        w_ms[i,j] = w_m
        massbalance[i,j] = omega_BL*(l_d**2 - l_m**2) + omega_m*l_m**2
        
    #    print '-----------------------------------'
    #    print 'balances:'
    #    print 'E - P:', E - P
    #    print 'mass balance:', omega_m*l_m**2 + omega_BL*(l_d**2 - l_m**2)
    #    print 'p_LCL (hPa)', p_LCLs[i,j]/1e2
    #
    #
    ##check that the l_m is consistent with f*q_sat in moist region BL
    #
    #    LHS = RH_m*q_sat
    #    
    #    RHS= q_sat + 2*zhat*(q_FA - q_sat)/(l_d**2 - l_m**2)*(l_m + zhat - (zhat + l_d)*np.exp((l_m - l_d)/zhat))
    #    
    #    omega_mP = g*(Q_mtrop + (L_v*P)/(np.pi*l_m**2))/(dels_trop)
    #    
    #    print 'RH in moist region BL = {0}'.format(qBL_m/qBL_sat)
    #    print 'water vapor given by f:', RHS
    #    print 'water vapor given by l_m:', LHS
    #    print 'T_s = {0}'.format(T_s)
    #    print 'check moisture balance:'
    #    print 'E = {0} (kg s^-1)'.format(E)
    #    print 'P/(pi*l_m**2*rho_w) = {0} (mm/day)'.format((P*3600*24*1000)/(rho_w*np.pi*l_m**2))
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
    #    print 'w_BL = {0} m s^-1'.format(w_BLs[i,j])
    #    print 'w_m = {0} m s^-1'.format(w_m)
    #    print 'l_m = {0} km'.format(l_m/1e3)
    #    print 'l_d = {0} km'.format(l_d/1e3)
    #    print 'convective fractional area = {0}'.format((l_m/l_d)**2)
    #    print 'M_BL = {0} kg/m^2, M_trop = {1} kg/m^2, M_LCL = {2} km/m^2'.format(M_BL, M_trop, M_LCL)
    #    print 'error = {0}'.format(np.abs((massdown - massup)/massdown))
    #    print ''
    #    #evap_cool = L_v*E
    #    #o_warm =cpw*W_o*np.pi*(l_d**2 - l_m**2)*M_o
    #    #print 'check energy balance:'
    #    #print 'radiative warming in convective region free troposphere (K/day)', W_m*3600*24
    #    #print 'radiative cooling in dry region free troposphere (K/day)', W_d*3600*24
    #    #print 'evaporative cooling = {0} (J  s^-1)'.format(evap_cool)
    #    #print 'ocean radiative warming = {0} (J  s^-1)'.format(o_warm) 
    #    #print 'surface energy unbalance = {0}'.format(o_warm - evap_cool)
    #    print ''
    #    moist_warming = g*Q_mtrop*l_m**2 + (g*L_v*P)/np.pi
    #    dry_cooling = -g*Q_dtrop*(l_d**2 - l_m**2)*(dels_LCL/dels_trop)
    #    print 'troposphere radiative warming + condensation warming in moist region  = {0} (J  s^-1)'.format(moist_warming)
    #    print 'dry region radiative cooling in free troposphere = {0} (J m^-1 s^-1)'.format(dry_cooling)
    #    print 'error = {0}'.format(np.abs((moist_warming - dry_cooling)/dry_cooling))
    #    print 'T_a = {0} (K)', T_a
        
            
 
cpBL = np.linspace(0, 1000, 100)
clm = np.linspace(0, l_d/1e3, 100)
cdelh = np.linspace(-1e7, 1e7, 100)
cmass = np.linspace(-4e10, 4e10, 100)
cwBL = np.linspace(-0.01, 0, 100)
cwm = np.linspace(0,0.1,100)
cnetP = np.linspace(-1e8, 1e8, 100)


        
plt.figure(1)
cs=plt.contourf(eps_BLss, eps_ass, p_BLs/1e2, cpBL, cmap=cm.RdYlBu_r, extend='both')
cs.cmap.set_under('grey')
cs.cmap.set_over('white')
plt.xlabel(r'$\epsilon_{BL}$')
plt.ylabel(r'$\epsilon_{a}$')
plt.colorbar(label=r'$p_{BL}$ (hPa)')
plt.show()

plt.figure(2)
cs=plt.contourf(eps_BLss, eps_ass, w_BLs, cwBL, cmap=cm.RdBu_r, extend='both')
cs.cmap.set_under('grey')
cs.cmap.set_over('white')
plt.xlabel(r'$\epsilon_{BL}$')
plt.ylabel(r'$\epsilon_{a}$')
plt.colorbar(label=r'$w_{BL}$ (m/s)')
plt.show()

plt.figure(3)
cs=plt.contourf(eps_BLss, eps_ass, w_ms, cwm, cmap=cm.RdBu_r, extend='both')
cs.cmap.set_under('grey')
cs.cmap.set_over('white')
plt.xlabel(r'$\epsilon_{BL}$')
plt.ylabel(r'$\epsilon_{a}$')
plt.colorbar(label=r'$w_{m}$ (m/s)')
plt.show()

plt.figure(4)
cs=plt.contourf(eps_BLss, eps_ass, l_ms/1e3, clm, cmap=cm.RdYlBu_r, extend='both')
cs.cmap.set_under('grey')
cs.cmap.set_over('white')
plt.xlabel(r'$\epsilon_{BL}$')
plt.ylabel(r'$\epsilon_{a}$')
plt.colorbar(label=r'$l_{m}$ (km)')
plt.show()

plt.figure(5)
cs=plt.contourf(eps_BLss, eps_ass, massbalance, cmass, cmap=cm.RdBu_r, extend='both')
cs.cmap.set_under('grey')
cs.cmap.set_over('white')
plt.xlabel(r'$\epsilon_{BL}$')
plt.ylabel(r'$\epsilon_{a}$')
plt.colorbar(label=r'mass balance (kg/s)')
plt.show()

plt.figure(6)
cs=plt.contourf(eps_BLss, eps_ass, netprecip, cnetP, cmap=cm.RdBu_r, extend='both')
cs.cmap.set_under('grey')
cs.cmap.set_over('white')
plt.xlabel(r'$\epsilon_{BL}$')
plt.ylabel(r'$\epsilon_{a}$')
plt.colorbar(label=r'net precipitation (kg/s)')
plt.show()

plt.figure(7)
cs=plt.contourf(eps_BLss, eps_ass, delhs, cdelh, cmap=cm.RdBu_r, extend='both')
cs.cmap.set_under('grey')
cs.cmap.set_over('white')
plt.xlabel(r'$\epsilon_{BL}$')
plt.ylabel(r'$\epsilon_{a}$')
plt.colorbar(label=r'$\Delta h_{trop}$')
plt.show()

plt.figure(8)
cs=plt.contourf(eps_BLss, eps_ass, RH_ms, np.linspace(0,100,100), cmap=cm.RdYlBu_r, extend='both')
cs.cmap.set_under('grey')
cs.cmap.set_over('white')
plt.xlabel(r'$\epsilon_{BL}$')
plt.ylabel(r'$\epsilon_{a}$')
plt.colorbar(label=r'relative humidity')
plt.show()






#    
#    
#    
#    
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
#plt.figure(1)
#plt.plot(l_ds/1e3, massbalance, 'o-')
#plt.xlabel(r'l$_d$ (km)')
#plt.ylabel(r'mass down - mass up (Pa m$^{2}$ s$^{-1}$)')
#plt.savefig(fout + 'massbalance_varyld.pdf')
#plt.close()
#
#plt.figure(2)
#plt.plot(l_ds/1e3, waterbalance, 'o-')
#plt.xlabel(r'l$_d$ (km)')
#plt.ylabel(r'E - P (mm/day)')
#plt.savefig(fout + 'waterbalance_varyld.pdf')
#plt.close()
#
#plt.figure(3)
#plt.plot(l_ds/1e3, Ebalance, 'o-')
#plt.xlabel(r'l$_d$ (km)')
#plt.ylabel(r'condensation warming - radiative cooling (J s$^{-1}$)')
#plt.savefig(fout + 'Ebalance_varyld.pdf')
#plt.close()
#
#plt.figure(4)
#plt.plot(l_ds/1e3, l_ms/1e3, 'o-')
#plt.xlabel(r'l$_d$ (km)')
#plt.ylabel(r'l$_m$ (km)')
#plt.savefig(fout + 'lm_varyld.pdf')
#plt.close()
#
#plt.figure(5)
#plt.plot(l_ds/1e3, (l_ms**2)/(l_ds**2), 'o-')
#plt.xlabel(r'l$_d$ (km)')
#plt.ylabel(r'$\sigma$')
#plt.savefig(fout + 'sigma_varyld.pdf')
#plt.close()
#
#plt.figure(6)
#plt.plot(l_ds/1e3, w_ms, 'o-')
#plt.xlabel(r'l$_d$ (km)')
#plt.ylabel(r'w$_m$ (m/s)')
#plt.savefig(fout + 'wm_varyld.pdf')
#plt.close()
#
#plt.figure(7)
#plt.plot(l_ds/1e3, RHs, 'o-')
#plt.ylim(0,1)
#plt.xlabel(r'l$_d$ (km)')
#plt.ylabel(r'RH$_m$')
#plt.savefig(fout + 'RHm_varyld.pdf')
#plt.close()
#
#plt.figure(8)
#plt.plot(l_ds/1e3, Ps, 'o-')
#plt.xlabel(r'l$_d$ (km)')
#plt.ylabel(r'P (mm/day)')
#plt.savefig(fout + 'P_varyld.pdf')
#plt.close()
#
#plt.figure(9)
#plt.plot(l_ds/1e3, delhs, 'o-')
#plt.xlabel(r'l$_d$ (km)')
#plt.ylabel(r'$\Delta h_{LCL}$ (J)')
#plt.savefig(fout + 'deltah_varyld.pdf')
#plt.close()
#
#plt.figure(10)
#plt.plot(l_ds/1e3, p_LCLs/1e2, 'o-')
#plt.xlabel(r'l$_d$ (km)')
#plt.ylabel(r'$p_{LCL}$ (hPa)')
#plt.savefig(fout + 'pLCL_varyld.pdf')
#plt.close()


        

