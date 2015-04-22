#!/usr/bin/env python

from pylab import *

import numpy
import sys
import thermo
import rootfinder
    
def Tfind(Tguess, p, theta_l, rt):
    # due to the check in invert_theta_l, we can assume that r = r_star.
    # However, we don't know the exact value of r_star because we don't know T.
    r = thermo.r_star(p, Tguess)
    rl = rt - r
    return theta_l - thermo.theta_l(p, Tguess, r, rl)
    
def invert_theta_l(theta_l, rt, p):
    T = thermo.theta_to_T(theta_l, p)
    r_star = thermo.r_star(p, T)
    if r_star < rt:
        T0 = numpy.array([233.15, theta_l])
        T = rootfinder.fzero(Tfind, T0, p, theta_l, rt)
    return T

def make_rt_vs_theta_l(tlmin, tlmax, rtmin, rtmax, p):
    #make a blank skewT diagram
    clf()
    
    #get a dense range of p, t0 to contour
    theta_l_vals = linspace(tlmin, tlmax, 100)
    r_vals = linspace(rtmin, rtmax, 100)
    theta_l_vals, r_vals = meshgrid(theta_l_vals, r_vals)
    
    Tvals = zeros(theta_l_vals.shape, float)

    for i in range(100):
        for j in range(100):
            Tvals[j,i] = invert_theta_l(theta_l_vals[j,i], r_vals[j,i], p)
    
    #use the real (data) value to get the potential temperature
    r_star = thermo.r_star(p, Tvals)
    rl = r_vals - r_star
    rl[rl < 0.] = 0.
    r = r_vals - rl

    theta_v = thermo.theta_v(p, Tvals, r)
    
    contour(theta_l_vals, r_vals*1000., -theta_v, 20, colors='k', linestyles=':')
    
    axis([tlmin, tlmax, rtmin*1000., rtmax*1000.])
    title('(theta_l, rt) Conserved Variable Diagram')
    ylabel('rt (g kg-1)')
    xlabel('theta_l (K)')
    
def plot_rt_vs_theta_l(p, T, r, rl):
    rt = r + rl
    theta_l = thermo.theta_l(p, T, r, rl)
    plot(theta_l, rt*1000)
    
def get_sounding(filename):
    T = []
    p = []
    RH = []
    soundingfile = open(filename)
    for line in soundingfile:
        line = line.split()
        p.append( float(line[0])*100. )
        T.append( float(line[1]) + 273.15 )
        RH.append( float(line[2]) )
        
    result = {}
    result['T'] = numpy.array(T)
    result['p'] = numpy.array(p)
    result['RH'] = numpy.array(RH)
    return result

def main(tmin, tmax, rmin, rmax):
    result = get_sounding('sounding.txt')
    T = result['T']
    p = result['p']
    RH = result['RH']
    rt = thermo.p_T_RH_to_r(p, T, RH)
        
    make_rt_vs_theta_l(270, 310, 0./1000., 6./1000., 84000.)
    
    r_star = thermo.r_star(p, T)
    rl = rt - r_star
    rl[rl < 0] = 0.
    r = rt - rl
        
    plot_rt_vs_theta_l(p, T, r, rl)
    axis([270,310,0,6])
    
    show()
 
if __name__=="__main__":
      #print "example: python tdd.py -10 30 400 1000"
      #print "argv:", sys.argv[0], sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
      #tmin  = float(sys.argv[1]) # minimum temperature in the TDD 
      #tmax  = float(sys.argv[2]) # maximum temperature in the TDD
      #pmin  = float(sys.argv[3]) # minimum pressure in the TDD
      #pmax  = float(sys.argv[4]) # maximum pressure in the TDD
      main(283, 294, 5, 18)   # main driver to plot the TDD and soundings
