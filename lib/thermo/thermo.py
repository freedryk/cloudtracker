#!/usr/bin/env python

import numpy
import rootfinder

# Definition of thermodynamic constant:
Cpd = 1005.7  # Heat capacity at constant pressure for dry air [J kg^-1 K^-1]
Cpv = 1870.0  # Heat capacity at constant pressure of water vapor [J kg^-1 K^-1]
Cl = 4190.0   # Heat capacity of liquid water [J kg^-1 K^-1]
Rv = 461.5    # Gas constant of water vapor [J kg^-1 K^-1]
Rd = 287.04   # Gas constant of dry air [J kg^-1 K^-1]
Lv0 = 2.501e6 # Latent heat of vaporization at 0 deg C [J kg^-1]
g = 9.80616   # Accelleration of gravity [m s^-2]
rho_l = 1000.0 # Density of water [kg m^-3]
p_0 = 100000.  # Reference pressure [Pa]
epsilon = Rd/Rv

def invert_Se(Se, z, r):
    # Given Moist Static Energy Se [J], mixing ratio r and height z [m],
    # return Temperature T [K]
    T0 = (Se - Lv0*r - g*z)/Cpd
    def Tfind_Se(T, Se, z, r): return Se - Se_r(T, z, r)
    T = rootfinder.fzero(Tfind_Se, T0, Se, z, r)
    return T

def Tfind_thetal(p, T, thetal, rt):
    r = r_star(p, T)
    rl = rt - r
    err = thetal - theta_l(P, T, r, rl)
    return err

def invert_theta_l(thetal, p, rt):
    # Given Liquid Water Potential Temperature thetal [K],
    # pressure p [Pa] and total mixing ratio rt,
    # return Temperature T [K]
    T = theta_to_T(thetal, p)
    rstar = r_star(p, T)
    if r_star < rt:
        T0 = array([230., thetal])
        T = rootfinder.fzero(Tfind_thetal, p, T0, thetal, rt)
    return T

#   temperature [K] from pressure P [Pa] and potential temperature theta [K]
def theta_to_T(theta, p): return theta*(p/p_0)**(0.2854)

#  denisty rho (kg m^-3) from pressure p [Pa], temperature T [K]
#  and mixing ratio r [kg kg^-1]
def rho(p, T, r):
    p_v = e(r, p)
    p_d = p - p_v
    return (p_d/Rd + p_v/Rv)/T

#   specific humidity q [kg/kg] from mixing ratio r [kg/kg]
def r_to_q(r): return r / (1. + r)

#   mixing ratio r [kg/kg] from specific humidity q [kg/kg]
def q_to_r(q): return q / (1. - q)

def p_T_RH_to_r(p, T, RH):
    #   mixing ratio [Kg/Kg]
    #   from pressure P [Pa], temperature T [K]
    #   and relative humidity RH (range 0.0-1.0)
    #   Note: RH=p_v/p_{vsat} and mix rat = eps*p_v/(p-p_v) 
    result = RH*epsilon*e_star(T)/(p - RH*e_star(T))
    result[e_star(T) >= 0.616*p] = 0.0
        
    return result

# vapor pressure $e$ Pa from mixing ratio, $r$ (Kg/Kg) and pressure, $P$ (Pa).
def e(r, p): return p*r/(epsilon + r)

#   Saturation vapor pressure over water (Pa) (Emanuel) from temperature T (K)
def e_star(T):
    return 100.*numpy.exp(53.67957 - 6743.769/T - 4.8451*numpy.log(T))

#   saturation mixing ratio rs (Kg/Kg) from temperature t (K) and pressure P (Pa)
def r_star(p, T):
    return epsilon*e_star(T)/(p - e_star(T))

def destar_dT(T):
    return e_star(T)*(6743.769/T/T - 4.8451/T)

def drstar_dT(p, T):
    return epsilon*p*destar_dT(T)/((p-e_star(T))**2)

def dqstar_dT(p, T):
    return drstar_dT(p, T)/((1+r_star(p, T))**2)

# dew-point temperature [K] from pressure P(Pa), temperature t(K), mixing ratio r(g/g)
#   by inverting Bolton's formula. Note: Centigrade +273.15 --> K
def T_d(r, p): return 243.5/((17.67/numpy.log(e(r, p)/611.2)) - 1.) + 273.15

#   Virtual temperature [K] from temperature T [K], mixing ratio r [Kg/Kg]
def Tv_r(T, r, rl): return T*(1. + r/epsilon)/(1. + r + rl)

#   Virtual temperature [K] from temperature T [K], specific humidity q [Kg/Kg]
def Tv_q(T, q, ql):
    r = q_to_r(q)
    rl = q_to_r(ql)
    return T*(1. + r/epsilon)/(1. + r + rl)

#   potential temperature [K] from pressure P [Pa] and temperature T [K]
def theta(p, T): return T*(p_0/p)**(Rd/Cpd)
    
def theta_v(p, T, r, rl):
    return Tv_r(theta(p, T), r, rl)
    
def theta_l(p, T, r, rl):
    # Liquid Water potential temperature [K]
    # from pressure p [Pa], temperature T [K], vapor mixing ratio r [kg/kg]
    # and water mixing ration rl [kg/kg]
    rt = r + rl
    chi = (Rd + rt*Rv)/(Cpd + rt*Cpv)
    gamma = (rt*Rv)/(Cpd + rt*Cpv)
    return T * (p_0/p)**chi * (1 - rl/(epsilon + rt))**chi * (1 - rl/rt)**(-gamma) * numpy.exp(-Lv(T)*rl/(Cpd + rt*Cpv)/T)

def theta_e(p, T, r, rl):
    #   Reversably defined equivalent potential temperature
    #   based on the defination in Emanuel 94 [K]
    #   from pressure p[Pa], temperature t[K], vapor mixing ratio r[kg/kg],
    #   and water mixing ratio rl[kg/kg]
    es = e_star(T)
    RH = e(r, p)/es
    pd = p - e(r, p)
    rt = r + rl
                
    LV = Lv0 - (Cl - Cpv)*(T - 273.15)
    the1 = T*(p_0/pd)**(Rd/(Cpd + rt*Cl))
    the2 = RH**(-r*Rv/(Cpd + Cl*rt))
    the3 = numpy.exp(LV*r/((Cpd + Cl*rt)*T))
    the = the1*the2*the3
    return the
    
# Computes and returns the relative humidity from temperature, $T$ (K) and the dewpoint temperature, $T_d$ (K).
def RH_ttd(T, T_d): return esatb(T_d) / esatb(T)

def q_star(T, p):
    #   saturation specific humidity qs (Kg/Kg) from temperature t (K) and pressure P (Pa)	
    rs = r_star(p, T)
    return r_to_q(rs)

def h_RH(RH, p, gz, T):
    #   liquid water static energy per unit of moist air (including total water) h(J/kg)
    #   note the difference with liquid-water static energy per unit of dry air
    #   from relative humidity RH, pressure p (Pa) and temperature t (K)
    #    z=z_stp(p)
    #    gz=g*z
    LV = Lv0 - CpvMCL*(T - 273.15)    
    rs = r_star(p, T)
    estar = p*rs/(0.622 + rs)    
    e = RH*estar
    qv = 0.622*e/(p - e*(1. - 0.622))
    CPN = Cpd*(1. - qv) + Cpv*qv
    return CPN*t + gz

def qt_RH(RH, p, T):
    #   specfic total water content qt (g/g) (per unit of moist air including total water)
    #   from relative humidity RH, pressure p (Pa) and temperature t (K)
    rs = r_star(p, T)
    estar = p*rs/(0.622 + rs)    
    e = RH*estar
    return 0.622*e/(p - e*(1. - 0.622))

def T_lcl(p, T, r):
    #   Lifting condensation level temperature [K]
    #   from pressure p [Pa], temperature t [K], r=mixing ratio [Kg/Kg]
    #   Ref.: Bolton, Eq. 21  (good to 0.1 K)
    epress = e(r, p)
    result = 55. + 2840./(3.5*numpy.log(T) - numpy.log(epress) - 4.805)
    result[epress < 1.E-30] = 1.0
    return result

def theta_ep_RH(p, t, RH):
    #   Pesudo-equivalent potential temperature [K]
    #   (set RH=1.0 if saturation e.p.t. wanted)
    #   from pressure p [Pa], temperature t [K], relative humidity RH 
    #   Ref.: Bolton, Eq. 43 (accuracy: 0.3K)
    #   it can also return potential temperature th [k] if unsaturated but you force RH=0
    #   but note the pesudo-equiv-pot-temp. is not equal to the poten.-temperature     
    r = p_T_RH_to_r(p, t, RH)
    ARG = (3376./T_lcl(p, t, r) - 2.54)*r*(1. + 0.81*r)
    thep = theta(p, t)
    thep[r > 0.0] = theta_v(p[r > 0.0], t[r > 0.0], r[r > 0.0], 0.)*numpy.exp(ARG[r > 0.0])
    
    return thep

def Se_r(T, z, r):
    # Moist Static Energy (J)
    return T*Cpd + g*z + Lv(T)*r

def Lwse_r(T, z, rl):
    # Liquid Water Static Energy (J)
    return T*Cpd + g*z - Lv(T)*rl
        
def h_tqt(t, qt, p, gz):
    qs = q_star(t, p)
    if qt > qs:
        LV = Lv0 - (Cl - Cpv)*(t - 273.15)
        CPN = Cpd*(1. - qs) + Cpv*qs
        h = CPN*t - LV*(qt - qs) + gz
    elif qt <= qs:
        CPN = Cpd*(1. - qt) + Cpv*qt
        h = CPN*t + gz
    return h

def h_rl(rl, p, t):
    #   liquid water static energy per unit of moist air (including total water)  h(J/kg) 
    #   from liquid water mixing ratio rl (g/g), pressure p(Pa), temperature t(K)
    z = z_stp(p)
    gz = g*z
    LV = Lv0 - (Cl - Cpv)*(t - 273.15)
    qs = q_star(t,p)
    CPN = Cpd*(1. - qs) + Cpv*qs
    ql = r_to_q(rl)
    return CPN*t - LV*ql + gz

def qt_rl(rl, p, t):
    #   specific total water content qt (g/g) (per unit of moist air including total water) 
    #   from liquid water mixing ratio rl (g/g), pressure p(Pa), temperature t(K)
    qs = q_star(t, p)
    ql = r_to_q(rl)
    return qs + ql

# Latent heat of vaporization
def Lv(T): return 2.5008E+6*(273.15/T)**(0.167 + 3.67E-4*T)

# Latent heat of fusion
def Lf(T): return 3.337E+5 + T*(2031. - 10.47*T)

# Latent heat of sublimation
def Lf(T): return 2.6332E+6 + T*(1727. - 3.625*T)

def tv_rrl(t, r, rl):
    #   Virtual temperature including water content loading [K]
    #   from temperature t [K], vapor mixing ratio r [g/g]
    #   and liquid water mixing ratio rl [g/g]
    
    #return t * ( 1.+ r / epsilon ) / ( 1.+ r +rl )
    return Tv_r(t, r) * (1.- rl/(1. + rl + r))

def T_rho(T, qt, p):
    #   density temperature T_rho(K) for both saturated and unsaturated air
    #   from temperature t(K), specific total water content qt (g/g)
    #   and pressure p (mb) 
    rs = r_star(p, T)
    qs = r_to_q(rs)
    if (qt <= qs):
        result = T*(1. + qt/0.622 - qt)
        #print "unsat",result
    else:
        rt = q_to_r(qt)
        result = t*(1. + rs/0.622)/(1. + rt)
        #print "sat",qt-qs
    return result

def tmu(h, qt, p, gz):
    #   temperature calculated assuming it is unsaturated air
    #   from liquid water static energy per unit of moist air h(J/kg)
    #   specific total water content qt (g/g), pressure p (mb) and 
    #   geopotential height gz (m2/s2)

    t = 300.
    count = 0
    while 1 > 0:
        tgess1 = t
        LV = Lv0 - (Cl - Cpv)*(tgess1 - 273.15)
        CPN = Cpd*(1. - qt) + Cpv*qt
        h1 = CPN*tgess1 + gz
        if abs(h1 - h) < 0.01:
            result = tgess1
            return result
        tgess2 = t - 1.
        LV = Lv0 - (Cl - Cpv)*(tgess2 - 273.15)
        CPN = Cpd*(1. - qt) + Cpv*qt
        h2 = CPN*tgess2 + gz 
        t = t + (h1 - h)/(h2 - h1)
        count = count + 1
        if count == 100:
            result = -1
            return result
            
def tms(h, qt, p, gz):
    #   temperature calculated assuming it is saturated air
    #   from liquid water static energy per unit of moist air h(J/kg)
    #   specific total water content qt (g/g), pressure p (mb) and 
    #   geopotential height gz (m2/s2)
	
    #   z=z_stp(p); gz=g*z
    t = 300.
    count = 0
    while 1 > 0:
        tgess1 = t
        LV = Lv0 - (Cl - Cpv)*(tgess1 - 273.15)
        qs = q_star(tgess1, p)
        CPN = Cpd*(1. - qs) + Cpv*qs
        h1 = CPN*tgess1 - LV*(qt - qs) + gz
        if abs(h1 - h) < 0.01:
            result = tgess1
            return result
        tgess2 = t - 1
        LV = Lv0 - (Cl - Cpv)*(tgess2 - 273.15)
        qs = q_star(tgess2, p)
        CPN = Cpd*(1. - qs) + Cpv*qs
        h2 = CPN*tgess2 - LV*(qt - qs) + gz
        t = t + (h1 - h)/(h2 - h1)
        count = count + 1
        if count == 100:
            return -1

def t_uos(h, qt, p, gz):
    #   temperature for both saturated and unsaturated air
    #   from liquid water static energy per unit of moist air h(J/kg)
    #   specific total water content qt (g/g), pressure p (mb) and 
    #   geopotential height gz (m2/s2)

    result = {}
    t1 = tmu(h, qt, p, gz)
    t2 = tms(h, qt, p, gz)        
    qs1 = q_star(t1, p)
    qs2 = q_star(t2, p)
    if (qt > qs1 and qt > qs2):
        result["T"] = t2
        result["QL"] = qt - qs2
        result["X"] = 1
    else:
        result["T"] = t1
        result["QL"] = 0.0
        result["X"] = 0
    return result

def all_uos(h, qt, p, gz):
    #   temperature and other thermodynamic variables 
    #   from liquid water static energy per unit of moist air h(J/kg)
    #   specific total water content qt (g/g), pressure p (mb) and 
    #   geopotential height gz (m2/s2)
    result = {}
    test = t_uos(h, qt, p, gz)
    t = test["T"]
    result["T"] = t
    qs = q_star(t, p)
    result["ES"] = esatb(t)
    if test["X"] == 1:
        result["Q"] = qs
        result["Ql"] = qt - qs
        result["r"] = q_to_r(qs)
        result["rl"] = q_to_r(qt - qs)
        result["E"] = result["ES"]
        result["RH"] = 100.
        result["TRO"] = T_rho(t, qt, p)
    else:
        result["Q"] = qt
        result["Ql"] = 0
        result["r"] = q_to_r(qt)
        result["rl"] = 0
        result["E"] = e(result["r"], p)
        result["RH"] = result["E"]/result["ES"]
        result["TRO"] = T_rho(t, qt, p)
    return result

def T_theta_e(the, p, t):
    #   Computes temperature t (K)
    #   from pseduo-equiv. pot. temp. (theta-e) and pressure p (Pa).
    #   input t is the first guess of temperature i.e 290 K

    #   ratio rgas / cp where gas constant rgas=287.04
    #   specific heat of air cp=1005.7, t is the initial guess
    tmp = (p_0/p)**Rd/Cpd
    count = 0
    tguess1 = t
    th     = tguess1*tmp
    qs     = q_star(tguess1, p)
    the1   = th*numpy.exp( 2670.*qs / tguess1 )
    
    while abs(the1-the) > 0.01:
        tgeuss2 = t - 1.
        th2    = tguess2*tmp
        qs2    = q_star(tguess2,p)
        the2   = th2*numpy.exp( 2670.*qs2 / tguess2 )
        t      = t + ( the1 - the ) / ( the2 - the1 )
        count  = count + 1
        if count > 100:
	        raise "error in function T_theta_e()" 
        tguess1 = t
        th     = tguess1*tmp
        qs     = q_star(tguess1,p)
        the1   = th*numpy.exp( 2670.*qs / tguess1 )
        
    return t
    
 
if __name__=="__main__":
    pass
#      print "example: python tdd.py -10 30 400 1000"
#      print "argv:", sys.argv[0], sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
#      tmin  = float(sys.argv[1]) # minimum temperature in the TDD 
#      tmax  = float(sys.argv[2]) # maximum temperature in the TDD
#      pmin  = float(sys.argv[3]) # minimum pressure in the TDD
#      pmax  = float(sys.argv[4]) # maximum pressure in the TDD
#      _main(tmin, tmax, pmin, pmax)   # main driver to plot the TDD and soundings
