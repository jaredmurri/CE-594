#####Dependencies#####
import numpy as np
import scipy
import math
import matplotlib.pyplot as plt

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
E= 69000000000 #Pa=(N/m)
A= outer_width**2-inner_width**2 #m^2

#####functions####
def diameter (x):
    if x <= L3:
        return D3
    elif x>=L3+L2:
        return D1
    else:
        return D3-((D3-D1)/L2)*(x) #m

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
        points= np.array([0,(-1/3)*np.sqrt(5-2*np.sqrt(10/7)),(1/3)*np.sqrt(5-2*np.sqrt(10/7)),(-1/3)*np.sqrt(5+2*np.sqrt(10/7)),(1/3)*np.sqrt(5+2*np.sqrt(10/7))])
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
    nodes=num_elements+1 #determines the number of nodes. 
    stiffness_matrix=np.zeros((nodes,nodes)) #creates an empty glabal matrix
    Length_element=Ltotal/num_elements #length of an element
    num_points=5
    local_stiff=np.zeros((2,2))
    gauss_value =0 #initializes variable
    map=2/Length_element #initializes variable
    for m in range(nodes-1): #loops through each element, filling the global matrix
        for i in range(0,2): #loops through the rows of the local stiffness matrix
            for j in range(0,2): #loops through the columsn of the local stiffness matrix
                for k in range(num_points): #loops through the gaussian points 
                    gauss_value=gauss_value+(gauss_weights(num_points)[k]*(ksi_der*((-1)**(j+1)))*(ksi_der*((-1)**(i+1))*(1/map)))
                local_stiff[i,j]=gauss_value #loads the local stiffness matrix with the results from each of the gaussian loops
                gauss_value=0 #reinitializes variable to zero for the next go round
        stiffness_matrix[m,m]=stiffness_matrix[m,m]+local_stiff[0,0]
        stiffness_matrix[m,m+1]=stiffness_matrix[m,m+1]+local_stiff[0,1]
        stiffness_matrix[m+1,m]=stiffness_matrix[m,m+1] #this make the A,-B and -A,B the same value
        stiffness_matrix[m+1,m+1]=stiffness_matrix[m+1,m+1]+local_stiff[1,1]
    stiffness_matrix[nodes-1,nodes-1]=0 #this is a boundary condition

    return stiffness_matrix[0:nodes-1,0:nodes-1]

def global_force(num_elements):
    nodes=num_elements+1 #determines the number of nodes.
    Length_element=Ltotal/num_elements #length of an element
    force_matrix=np.zeros(nodes) #creates an empty matrix
    force_matrix[0]=-self_weight
    for i in range(1,nodes-1):
        x_pos=i*Length_element
        x_pos1=(i+1)*Length_element
        force_matrix[i]= ((-1)*ksi_der*Length_element*(1/2)*(Force_distrib(x_pos)+Force_distrib(x_pos1)))+force_matrix[1-1]
    force_matrix[nodes-1]=-np.sum(force_matrix)
    return force_matrix[0:nodes-1]

def displacement_matrix(num_elements):
    return np.matmul(np.linalg.inv(global_stiff(num_elements)),global_force(num_elements))

def stress_matrix(num_elements):
    Length_element=Ltotal/num_elements #length of an element
    stress_matrix=np.zeros(num_elements) #creates an empty matrix
    displacement=displacement_matrix(num_elements) #brings in the displacements
    for i in range(0,num_elements-1):
        stress_matrix[i]= ((displacement[i+1]-displacement[i])/Length_element)/A #this is the derivative of the displacement per cross sectional area
    return stress_matrix

def strain_matrix(num_elements):
    strain_matrix=np.zeros(num_elements) #creates an empty matrix
    stress=stress_matrix(num_elements) #brings in the strain matrix
    for i in range(0,num_elements-1):
        strain_matrix[i]= stress[i]/E
    return strain_matrix

def max_stress(num_elements):
    return np.ndarray.max(stress_matrix(num_elements))

#####Plotting#####
def plot_stress(num_elements):
    x=np.linspace(0,Ltotal,num_elements)
    y= stress_matrix(num_elements)
    plt.plot(x,y, color="black")
    plt.show()

def plot_strain(num_elements):
    x=np.linspace(0,Ltotal,num_elements)
    y= strain_matrix(num_elements)
    plt.plot(x,y, color="black")
    plt.show()

def plot_displacement(num_elements):
    x=np.linspace(0,Ltotal,num_elements)
    y= displacement_matrix(num_elements)
    plt.plot(x,y, color="black")
    plt.show()