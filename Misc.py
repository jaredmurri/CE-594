#####Dependencies#####
import numpy as np
import scipy
import math 
 
def transformation_matrix(theta):
    u= 1/np.sqrt(3)
    theta=math.radians(theta)
    transformmatrix= np.array([[np.cos(theta)+(u**2)*(1-np.cos(theta)),(u**2)*(1-np.cos(theta))-u*np.sin(theta),(u**2)*(1-np.cos(theta))+u*np.sin(theta)],[(u**2)*(1-np.cos(theta))+u*np.sin(theta),np.cos(theta)+(u**2)*(1-np.cos(theta)),(u**2)*(1-np.cos(theta))-u*np.sin(theta)],[(u**2)*(1-np.cos(theta))-u*np.sin(theta),(u**2)*(1-np.cos(theta))+u*np.sin(theta),np.cos(theta)+(u**2)*(1-np.cos(theta))]])
    origonalmatrix=np.array([[1,-3,np.sqrt(2)],[-3,1,-np.sqrt(2)],[np.sqrt(2),-np.sqrt(2),4]])
    matrix= np.matmul(transformmatrix,origonalmatrix)
    return np.matmul(matrix,np.transpose(transformmatrix))