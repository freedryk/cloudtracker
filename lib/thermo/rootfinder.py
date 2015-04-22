#!/usr/bin/env python

import numpy

def find_interval(f, x, *args):
    x1 = x
    x2 = x
    if x == 0.:
        dx = 1./50.
    else:
        dx = x/50.
        
    maxiter = 40
    twosqrt = numpy.sqrt(2)
    a = x
    fa = f(a, *args)
    b = x
    fb = f(b, *args)
    
    for i in range(maxiter):
        dx = dx*twosqrt
        a = x - dx
        fa = f(a, *args)
        b = x + dx
        fb = f(b, *args)
        if (fa*fb < 0.): return (a, b)
        
    raise "Couldn't find a suitable range."

# This function evaluates a new point, sets the y range,
# and tests for convergence
def get_y(x, f, eps, ymax, ymin, *args):
    y = f(x, *args)
    ymax = max(ymax, y)
    ymin = min(ymin, y)
    converged = (abs(y) < eps*(ymax-ymin))
    return (y, ymax, ymin, converged)

def fzero(f, x, *args, **parms):
    # f is the function we wish to find the zeros of
    # x is an initial guess of the zero location.
    #   Can be a float or a sequence of two floats specifying a range
    # *args contains any other parameters needed for f
    # **parms can be eps (allowable error) or maxiter (max number of iterations.)

    if isinstance(x, (int, float)):
        if (f(x, *args) == 0): return x
        x1, x3 = find_interval(f, x, *args)
    elif isinstance(x, (tuple, list, numpy.ndarray)):
        if len(x) == 2:
            x1 = x[0]
            x3 = x[1]
        else:
            raise "problem"
    eps = 1.e-14
    maxiter = 40
    if 'eps' in parms: eps = parms['eps']
    if 'maxiter' in parms: maxiter = int(parms['maxiter'])
    xeps = eps*abs(x3-x1) #convergence criterion in x
    ymin = 1e-300         # set up for criterion in y
    ymax = -ymin
    converged = False
    
    # evaluate function at limits
    y1, ymax, ymin, converged = get_y(x1, f, eps, ymax, ymin, *args)
    if (y1 == 0): return x1    # only exact zero accepted here
    
    y3, ymax, ymin, converged = get_y(x3, f, eps, ymax, ymin, *args)
    if (y3 == 0): return x3    # only exact zero accepted here
    # initial points must straddle root
    if (y1 * y3 > 0): raise  "problem"

    # begin main bisection loop
    for i in range(maxiter):
        x2 = (x1 + x3)/2. # bisect x range
        # If we're bisecting, we stop when the interval gets small
        if (abs(x3 - x1) < eps): return x2
             
        y2, ymax, ymin, converged = get_y(x2, f, eps, ymax, ymin, *args)
        # perhaps we nailed the root?
        if converged: return x2

        # relabel points to keep root between x1 and x2
        if (y1*y2 > 0.):
            x3, x1 = x1, x3
            y3, y1 = y1, y3

        # here's where we try parabolic interpolation.
        # there are two criteria for accepting the point.
        # if any one of them fails, just bisect again.

        # note that y21 and y31 cannot be zero here
        # but y32 might be.
        y21 = y2 - y1
        y32 = y3 - y2
        if y32 != 0.:
            y31 = y3 - y1
            b = (x2 - x1)/y21
            c = (y21 - y32)/(y32*y31)

            # Test for the RTMI condition
            if(y3*y31 >= 2.0*y2*y21):
                xm = x1 - b*y1*(1. - c*y2)
                ym, ymax, ymin, converged = get_y(xm, f, eps, ymax, ymin, *args)

                # perhaps we nailed the root?
                if converged: return xm
                             
                # Relabel the points as needed
                # to keep root between x1 and x2
                if (ym*y1 < 0.0):
                    x2 = xm
                    y2 = ym
                else:
                    x1 = xm
                    y1 = ym
                
        # we didn't do parabolic interpolation.
        # just relabel points and continue.
        x3 = x2
        y3 = y2
    # end of bisection loop. If we get here, we did not converge.
    print "Iterate: Failed to converge in imax tries."
    return x2
    
def testfunc(x):
    return numpy.sin(x)
     
 
if __name__=="__main__":
    f = testfunc
    x = 1.
    print fzero(f, x)
    print fzero(f, x, eps=1e-300, maxiter = 80.)
