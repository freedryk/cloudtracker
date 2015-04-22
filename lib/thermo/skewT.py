#!/usr/bin/env python

from pylab import *
import matplotlib as mpl
import sys
import thermo
import scipy.optimize.zeros as scizeros

def make_skewT(tmin, tmax, pmax, pmin, skew=30.):
    #make a blank skewT diagram
    clf()
    #get a dense range of p, t0 to contour
    yplot = linspace(1050, 100, 100)
    xplot = linspace(-50, 50, 100)
    xplot, yplot = meshgrid(xplot, yplot)
    
    #lay down a reference grid that labels xplot,yplot points 
    #in the new (skewT-lnP) coordinate system .
    # Each value of the temp matrix holds the actual (data) temperature
    # label (in deg C)  of the xplot, yplot coordinate pairs
    #note that we don't have to transform the y coordinate
    #it's still the pressure

    #use the real (data) value to get the potential temperature
    T = xplot + skew*log(0.001*yplot)
    Tk = T + 273.15 #convert from C to K for use in thermo functios
    p = yplot*100. #convert from hPa to Pa

    th = thermo.theta(p, Tk) #theta labels

    #add the mixing ratio
    rstar = thermo.r_star(p, Tk)  #wsat labels

    #saturated adiabat, so Tdew=Tk
    thetaeVals = thermo.theta_e(p, Tk, rstar, 0.)    
    
    tempLabels = arange(-140., 50., 10.)
    con1 = contour(xplot, yplot, T, levels = tempLabels, colors = 'k', linewidths=.5)
    ax = gca()
    ax.set_yscale('log')
    lines = arange(100., 1100., 100.)
    yticks(lines, ['100','200','300','400','500','600','700','800','900','1000'])
    for line in lines:
        axhline(line, ls=':', color='k', linewidth=.5)
    thetaLabels = arange(200., 380., 10.)
    con2 = contour(xplot, yplot, th, levels = thetaLabels, colors='b', linewidths=.5)
    rsLabels = [.1,.2,.4,.6, 1, 2, 3, 4, 5, 6, 8, 10, 15, 20, 25, 30, 40]
    con3 = contour(xplot, yplot, rstar*1.e3, levels=rsLabels, colors='g', linewidths=.5)
    thetaeLabels = linspace(200,400,21)
    con4 = contour(xplot, yplot, thetaeVals, levels = thetaeLabels, colors='r', linewidths=.5)
    axis([tmin, tmax, pmax, pmin])
    clabel(con1, inline = False, fmt = '%1.0f')
    clabel(con2, inline = False, fmt = '%1.0f')
    clabel(con3, inline = False, fmt = '%1.1f')
    clabel(con4, inline = False, fmt = '%1.0f')
    title('skew T - lnp chart')
    ylabel('pressure (hPa)')
    xlabel('temperature (black, degrees C)')
            
def skewIt(T, p, skew=30.):
    pz = p*0.01
    tz = T - 273.15
    skewedTemp = tz - skew*log(0.001*pz)
    return (skewedTemp, pz)
    
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
    result['T'] = array(T)
    result['p'] = array(p)
    result['RH'] = array(RH)
    return result

def main(tmin, tmax, pmin, pmax):
    result = get_sounding('sounding.txt')
    T = result['T']
    p = result['p']
    RH = result['RH']
    r = thermo.p_T_RH_to_r(p, T, RH)
    T_dew = thermo.T_d(r, p)
    make_skewT(tmin, tmax, pmax, pmin)
    tee, pee = skewIt(T, p, 30.)
    plot(tee, pee, 'k-')
    tee, pee = skewIt(T_dew, p, 30.)
    plot(tee, pee, 'k-')
    axis([tmin, tmax, pmax, pmin])
    show()
 
if __name__=="__main__":
      #print "example: python tdd.py -10 30 400 1000"
      #print "argv:", sys.argv[0], sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
      #tmin  = float(sys.argv[1]) # minimum temperature in the TDD 
      #tmax  = float(sys.argv[2]) # maximum temperature in the TDD
      #pmin  = float(sys.argv[3]) # minimum pressure in the TDD
      #pmax  = float(sys.argv[4]) # maximum pressure in the TDD
      main(-50, 50, 100, 1000)   # main driver to plot the TDD and soundings
