import numpy as np
import matplotlib.pyplot as plt
from harmonic_poly import basis_upto_order

# this is a helper function to get the total number of basis functions (upto order l)
def get_no_of_poly(l):
    sum = 0
    sum_arry = [sum]
    for i in range(l+1):
        sum += 2*(i+1) + 1
        sum_arry.append(sum)
    return sum_arry[-1]

# order of harmonic poly 
order = 2

# define harmonic polynomial coefficients 
harmonic_coeff = np.zeros(get_no_of_poly(order))
hc_non_zero = np.array([0.001,1,0.002,0.0004,0.036,0.003,0.042,0.023,0.0023,0.0045,0.0056,0.0034,0.0073,0.0056,0.0017])
harmonic_coeff[:len(hc_non_zero)] = hc_non_zero
harmonic_coeff

# get the harmonic polynomials upto order l 
Phi_x,Phi_y,Phi_z = basis_upto_order(order) 

# this evaluate basis functions(upto order l) on numpy meshgrid 
def get_Phi_upto(l,xv,yv,zv):
    Phi_xx = Phi_x(xv,yv,zv)
    Phi_yy = Phi_y(xv,yv,zv)
    Phi_zz = Phi_z(xv,yv,zv)
    
    m = get_no_of_poly(l)
    n = xv.flatten().shape[0]
    Phi = np.zeros((3,n,m))
    for i in range(m):
        Phi[0,:,i] = np.array(Phi_xx[i]).flatten()
        Phi[1,:,i] = np.array(Phi_yy[i]).flatten()
        Phi[2,:,i] = np.array(Phi_zz[i]).flatten()
    return Phi[0],Phi[1],Phi[2]


# get the magentic field (upto order l) on a meshgrid 
def get_mag_field(l,xv,yv,zv,harmonic_coeff=harmonic_coeff):
    Phi_xx,Phi_yy,Phi_zz = get_Phi_upto(l,xv,yv,zv)
    Bxx = Phi_xx @ harmonic_coeff
    Bzz = Phi_zz @ harmonic_coeff
    Byy = Phi_yy @ harmonic_coeff
    return np.array([Bxx,Byy,Bzz]).T 
