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
weight=10000 #N
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

def global_stiff (num_elements):
    nodes=num_elements+1 #determines the number of nodes
    stiff_matrix=np.zeros((nodes,nodes)) #creates an empty glabal matrix
    Length_element=Ltotal/num_elements #length of an element

    #fill the matrix with the global stiffness values
    for i in range (1,nodes):
        stiff_matrix[i,i]=2*E*A/Length_element
        stiff_matrix[i-1,i]=-E*A/Length_element
        stiff_matrix[i,i-1]=-E*A/Length_element
    stiff_matrix[0,0]=E*A/Length_element
    stiff_matrix[nodes-1,nodes-1]=0 #this is a boundary condition

    return stiff_matrix

def force_matrix_create (num_elements):
    nodes=num_elements+1 #determines the number of nodes
    Length_element=Ltotal/num_elements #length of an element
    force_matrix=np.zeros((nodes,1)) #creates an empty force matrix
 
    for i in range(0,nodes): #loads the force values into the matrix
        lower=i*Length_element
        upper=(i+1)*Length_element
        force= scipy.integrate.quadrature(Force_distrib(nodes,lower,upper))
        force_matrix[i,1]= force
    return force_matrix

   
    
