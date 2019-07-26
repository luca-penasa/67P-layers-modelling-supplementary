import numpy as np 

#SL parameters for the reference solution
centre = np.array([1.06, -0.346, 0.01])
scales = np.array([ 1.,0.76, 0.704])
angles = np.array([ 28.1, -11.2, -7.3])

#BL parameters from the reference solution
centre = np.array([-0.473, 0.33, -0.17])
scales = np.array([ 1.,0.805 ,0.544])
angles = np.array([ 44.8 , 15.0,  66.3])


def euler2mat(angles):    
    '''
    :param angles: array or list of length 3. in order the alpha,
    beta, gamma angles in z-y'-x'' intrinsic Tait–Bryan angles
    returns the correspondent rotation matrix as numpy array.
    from https://en.wikipedia.org/wiki/Euler_angles#Rotation_matrix
    '''
    c0, c1, c2 = [np.cos(angle) for angle in angles]
    s0, s1, s2 = [np.sin(angle) for angle in angles]
    
    R = np.array([[ c0*c1, c0*s1*s2 - c2*s0, s0*s2 + c0*c2*s1 ],
                  [ c1*s0, c0*c2 + s0*s1*s2, c2*s0*s1 - c0*s2 ],
                  [  -s1,        c1*s2,            c1*c2      ]])
    return R

def estimateEllipsoidalModel(p, centre, R, axes):
    '''
    This method returns the predicted lambda and grandient value 
    for a given Ellipsoidal Model (EM)
    
    :param p: 3D coords of the point for which to evaluate the EM
    :param centre: 3D coords of the centre of the ellipsoidal model
    :param R: the roataion matrix representing the EM orientation 
    :param axes: the 3 lengths of the axes (first element must be 1.0)
    '''
    a,b,c = axes
    pt = R.T.dot(p - centre)    
    tmp = np.array([b*c, c, b])
    t = np.sqrt(np.dot(tmp**2, pt**2))
    lam = 1/(b*c) * t
    grad = R.dot((tmp * pt) / (t*axes))
    
    return lam, grad

# a test point
ptest= np.array([10,22,1.0])

R =  euler2mat(angles/ 180*np.pi) # rot matrix from the euler angles
lam, grad = estimateEllipsoidalModel(ptest, centre, R, scales)

print ("lambda -RES value- {} km".format(lam))
print ("gradient {}".format(grad))   
