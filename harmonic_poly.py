import numpy as np
import sympy as sp 
from sympy import assoc_legendre
import matplotlib.pyplot as plt

l,m = sp.symbols('l,m',real=True)
x,y,z = sp.symbols('x,y,z',real=True,positive=True)
r,theta,phi = sp.symbols('r,theta,phi')

# function to convert spherical coords to cartesian 
def SphericalToCartesian(f):
  return f.subs([(r,sp.sqrt(x**2 + y**2 + z**2)),
                 (theta,sp.acos(z/sp.sqrt(x**2 + y**2 + z**2))),
                 (phi,sp.asin(y/sp.sqrt(x**2 + y**2)))])

# function to generate basis function in cartisian coords 
def basisFunction(l,m):
  l+=1
  c_plus = sp.factorial(l-1) * (-2)**(sp.Abs(m)) /sp.factorial(l + m) * sp.cos(m * phi)
  c_negt = sp.factorial(l-1) * (-2)**(sp.Abs(m)) /sp.factorial(l + sp.Abs(m)) * sp.sin(sp.Abs(m) * phi)

  s_plus  = SphericalToCartesian(c_plus* r**l * assoc_legendre(l,m,sp.cos(theta)))
  s_minus = SphericalToCartesian(c_negt* r**l * assoc_legendre(l,sp.Abs(m),sp.cos(theta)))
  
  if m >= 0:
    return [sp.expand_trig(sp.diff(s_plus, axis)).simplify() for axis in [x,y,z]]
  else:
    return [sp.expand_trig(sp.diff(s_minus, axis)).simplify().expand() for axis in [x,y,z]]

####################################################################################
  # These functions evalue array of harmonic polynomials on array of coordinates 
####################################################################################
  
# get the array of basis fucntion upto order l as a function of x,y,z
def basis_upto_order(l):
  Phi= [[basisFunction(i,j)[k] for i in range(l+1) for j in range(-i-1,i+2)] for k in range(3)]
  Phi_lambda = [sp.lambdify((x,y,z),Phi[i]) for i in range(3)]
  Phi_X,Phi_Y,Phi_Z = Phi_lambda[0],Phi_lambda[1],Phi_lambda[2]
  return Phi_X,Phi_Y,Phi_Z

# all basis function at coordinate array, the order is Phi_z,Phi_X,Phi_Y
# coord_array = [[x_1,y_1,z_1],[x_2,y_2,z_2],....[x_n,y_n,z_n]]
def all_basis_at_coordarray(coord_list,Phi_X,Phi_Y,Phi_Z):
  C_x = np.array(list(map(Phi_X,*coord_list.T)))
  C_y = np.array(list(map(Phi_Y,*coord_list.T)))
  C_z = np.array(list(map(Phi_Z,*coord_list.T)))
  C   = np.array([C_z,C_x,C_y])
  return C.T

# basis functions of oder l evaluated at coordinate array 
def basis_at_coordarray(coord_list,order,Phi_X,Phi_Y,Phi_Z):
  C_x = np.array(list(map(Phi_X,*coord_list.T)))
  C_y = np.array(list(map(Phi_Y,*coord_list.T)))
  C_z = np.array(list(map(Phi_Z,*coord_list.T)))
  C   = np.array([C_z,C_x,C_y])
  start_index, end_index =order*(order+2),(order +1)*(order + 3)
  return C.T[start_index:end_index]

# get value of all basis function at a given point(probe location)
def all_basis_at_point(point,coord_list,Phi_X,Phi_Y,Phi_Z):
  return all_basis_at_coordarray(coord_list,Phi_X,Phi_Y,Phi_Z)[:,point,:]

# get the value of basis function of given order at a probe location 
# This may not be the most effcient way (look some code below)
def basis_at_point(point,coord_list,order,Phi_X,Phi_Y,Phi_Z):
  return basis_at_coordarray(coord_list,order,Phi_X,Phi_Y,Phi_Z)[:,point,:]

####################################################################################
  # These functions evalue a harmonic polynomials on one or many coordinate points
  # these are light weight functions (dont need to evalue poly of order l beforehand)
####################################################################################

# function to evaluate a given basis function (l,m) at given cartisian coordinate point(u,v,z)
def get_basis_fv(u,v,w,l,m):
    bf = basisFunction(l,m)
    bfv = sp.lambdify([x,y,z],bf)
    return bfv(u,v,w)

# function to evaluate a given harmonics on selected coordinate mesh 
def get_harmonics_on_selected(xv,yv,zv,bool_mask,l,m):
    xx,yy,zz = get_selected_coords(xv,yv,zv,bool_mask)
    fvx,fvy,fvz  = get_basis_fv(xx,yy,zz,l,m)
    return np.array([fvx,fvy,fvz]).T

# function to filter coordinates from a grid using boolean masking  
def get_selected_coords(xv,yv,zv,bool_mask):
    x_coords = xv[bool_mask]
    y_coords = yv[bool_mask]
    z_coords = zv[bool_mask] 
    return np.array([x_coords,y_coords,z_coords]).T

# helper function to get the index of array where a fucntion is max 
def get_idx_maxf(func, pos_bool):
    max_val = np.max(func[pos_bool])
    idx = np.where(func[pos_bool]==max_val)[0]
    return idx

# this finds the x,y and z coordinates where the function is max on given surface 
def get_coords_maxf(func,pos_bool,xv,yv,zv):
    idx = get_idx_maxf(func,pos_bool)
    x_coords = xv[pos_bool][idx]
    y_coords = yv[pos_bool][idx]
    z_coords = zv[pos_bool][idx]
    return np.array([x_coords,y_coords,z_coords]).T