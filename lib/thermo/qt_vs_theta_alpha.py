#!/usr/bin/env python

import numpy
from pylab import subplot
import sys
import thermo
    
def plot_rt_vs_theta_alpha(p, T, r, rl, alpha):
    theta_l = thermo.theta_l(p, T, r, rl)
    theta_e = thermo.theta_e(p, T, r, rl)
    theta_alpha = (1. - alpha)*theta_l + alpha*theta_e
    rt = r + rl
    
    ax = subplot(1,1,1)
    ax.plot(theta_alpha, rt*1000)
    
    yl = ax.get_ylim()
    if yl[0] < yl[1]: ax.set_ylim([yl[1], yl[0]])
    
def _get_sounding(filename):
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

def _main(tmin, tmax, rmin, rmax):
    result = _get_sounding('sounding.txt')
    T = result['T']
    p = result['p']
    RH = result['RH']
    rt = thermo.p_T_RH_to_r(p, T, RH)
    
    r_star = thermo.r_star(p, T)
    rl = rt - r_star
    rl[rl < 0] = 0.
    r = rt - rl
        
    plot_rt_vs_theta_alpha(p, T, r, rl, 1.)    
    show()
 
if __name__=="__main__":
      #print "example: python tdd.py -10 30 400 1000"
      #print "argv:", sys.argv[0], sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
      #tmin  = float(sys.argv[1]) # minimum temperature in the TDD 
      #tmax  = float(sys.argv[2]) # maximum temperature in the TDD
      #pmin  = float(sys.argv[3]) # minimum pressure in the TDD
      #pmax  = float(sys.argv[4]) # maximum pressure in the TDD
      _main(283, 294, 5, 18)   # main driver to plot the TDD and soundings
