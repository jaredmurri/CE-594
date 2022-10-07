#####Dependencies#####
import numpy as np
import scipy
import math 
 
def transformation_matrix(theta):
    u= 1/np.sqrt(3)
    transformmatrix= np.array([[np.cos(theta)+(u**2)*(1-np.cos(theta)),(u**2)*(1-np.cos(theta))-u*np.sin(theta),(u**2)*(1-np.cos(theta))+u*np.sin(theta)],[(u**2)*(1-np.cos(theta))+u*np.sin(theta),np.cos(theta)+(u**2)*(1-np.cos(theta)),(u**2)*(1-np.cos(theta))-u*np.sin(theta)],[(u**2)*(1-np.cos(theta))-u*np.sin(theta),(u**2)*(1-np.cos(theta))+u*np.sin(theta),np.cos(theta)+(u**2)*(1-np.cos(theta))]])
    origonalmatrix=np.array([[1,-3,np.sqrt(2)],[-3,1,-np.sqrt(2)],[np.sqrt(2),-np.sqrt(2),4]])
    return transformmatrix, origonalmatrix, np.matmul(transformmatrix,origonalmatrix)

    


 #variables for the member matrix#
#c2=(np.cos(theta))**2
#cs=np.cos(theta)*np.sin(theta)
#s2=(np.sin(theta))**2

    #returns the member matrix
#member_matrix= (E*A/L_element)*np.array([[c2,cs,-c2,-cs],[cs,s2,-cs,-s2],[-c2,-cs,c2,cs],[-cs,-s2,cs,s2]])
#member_matrix_short= np.delete(member_matrix,1,1) #deletes column 2
#member_matrix_short= np.delete(member_matrix_short,2,1) #deletes what was column 4 but is now column 3
#member_matrix_short=np.delete(member_matrix_short,1,0) #deletes row 2
#member_matrix_short=np.delete(member_matrix_short,2,0) #deletes what was row 4 but is now row 3

    # adds the member matrix into the global matrix
#for k in range(0,2):
#    for m in range (0,2):
#       stiff_matrix [i,j]=stiff_matrix[i,j]+member_matrix_short[k,m]



