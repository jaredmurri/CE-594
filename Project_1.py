import numpy as np
import scipy
import math
import scipy.integrate._quadrature as integrate

##### Givens  #####
density= 1153
kin_viscos= 0.008
L1=1
L2=4.5
L3=0.5
D3=4
D1=1
outer_width=0.07620
inner_width=0.07493
weight=10000


def diameter (x):
    if x <= L3:
        return D3
    elif x>=L3+L2:
        return D1
    else:
        return D3-((D3-D1)/L2)*(x)

def velocity(x):
    return 92/(diameter(x)**2)

def Reynolds_number (x):
    return (velocity(x)*diameter(x))/kin_viscos

def drag_coeff(x):
    return 0.027/(Reynolds_number(x))**(1/7)

def Force_per_unit_length(x):
    perimeter=outer_width*4
    return perimeter*.5*drag_coeff(x)*density(x)*(velocity(x))**2
