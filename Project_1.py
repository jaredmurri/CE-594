#####Dependencies#####
import numpy as np
import scipy
import math

##### Givens  #####
density= 1153 #kg/m^3
kin_viscos= 0.008 #m/s^2
L1=1 #m
L2=4.5 #m
L3=0.5 #m
Ltotal=L1+L2+L3
D3=4 #m
D1=1 #m
v1=92 #m/s
outer_width=0.07620 #m
inner_width=0.07493 #m
self_weight=10000 #N
E= 69000000000 #Pa
A= outer_width**2-inner_width**2 #m^2

#####functions####
def diameter (x):
    if x <= L3:
        return D3
    elif x>=L3+L2:
        return D1
    else:
        return D3-((D3-D1)/L2)*(x) #m/s

def velocity(x):
    return (v1*(D1**2))/(diameter(x)**2)

def Rey_num (x):
    return (velocity(x)*diameter(x))/kin_viscos

def drag_coeff(x):
    return 0.027/(Rey_num(x))**(1/7)

def Force_distrib(x):
    perimeter=outer_width*4
    return perimeter*.5*drag_coeff(x)*density*(velocity(x)**2) #N/m

#####FEA Stuff#####
ksi_der=1/2

def gaussquad_points(num_points):
    if num_points == 1:
        points= np.array([0])
    elif num_points == 2:
        points= np.array([1/np.sqrt(3),-1/np.sqrt(3)])
    elif num_points == 3:
        points= np.array([0,np.sqrt(3/5),-np.sqrt(3/5)])
    elif num_points == 4:
        points= np.array([np.sqrt((2/7)-(2/7)*np.sqrt(6/5)),-np.sqrt((2/7)-(2/7)*np.sqrt(6/5)),np.sqrt((2/7)+(2/7)*np.sqrt(6/5)),-np.sqrt((2/7)+(2/7)*np.sqrt(6/5))])
    else:
        points= np.array([0,-(1/3)*np.sqrt(5-2*np.sqrt(10/7)),-(1/3)*np.sqrt(5-2*np.sqrt(10/7)),(1/3)*np.sqrt(5+2*np.sqrt(10/7)),-(1/3)*np.sqrt(5+2*np.sqrt(10/7))])
    return points

def gauss_weights (num_points):
    if num_points == 1:
        weight= np.array([2])
    elif num_points == 2:
        weight= np.array([1,1])
    elif num_points == 3:
        weight= np.array([0,8/9,8/9])
    elif num_points == 4:
        weight= np.array([(18+np.sqrt(30))/36,(18+np.sqrt(30))/36,(18-np.sqrt(30))/36,(18-np.sqrt(30))/36])
    else:
        weight= np.array([128/225,(322+13*np.sqrt(70))/900, (322+13*np.sqrt(70))/900, (322-13*np.sqrt(70))/900,(322-13*np.sqrt(70))/900])
    return weight

def global_stiff (num_elements):
    nodes=num_elements #determines the number of nodes. This one has the same number of nodes as elements becasue the last one is a fixed condition
    stiff_matrix=np.zeros((nodes,nodes)) #creates an empty glabal matrix
    Length_element=Ltotal/num_elements #length of an element
    num_points = 5
    local_stiff=np.zeros((2,2))
    for i in range(1,3):
        for j in range(1,3):
            local_stiff[i-1,j-1]
            
