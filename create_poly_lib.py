import sympy as sp
import numpy as np

# Define symbolic variables
x, y, z = sp.symbols('x y z')
mat = [x,y,z]

# first order polynomials in (x,y,z)
f_1 = [[mat[i]] for i in range(3)]

# second order polynomials in (x,y,z)
f_2 = [[mat[i]*mat[j]] for i in range(3) for j in range(i,3)]

# Third order polynomials in (x,y,z)
f_3 = [[mat[i]*mat[j]*mat[k]] for i in range(3) for j in range(i,3) for k in range(j,3)]

# forth order polynomials in (x,y,z)
f_4 = [[mat[i]*mat[j]*mat[k]*mat[l]] for i in range(3) for j in range(i,3) for k in range(j,3) for l in range(k,3)]

# all the terms upto 4th order 
poly = [[1]] + f_1 + f_2 + f_3 + f_4

# number of poly terms upto given order 
terms_0 = 1
terms_1 = terms_0 + len(f_1)
terms_2 = terms_1 + len(f_2)
terms_3 = terms_2 + len(f_3)
terms_4 = terms_3 + len(f_4)

# this select polynomial terms upto given order(<=4)
def get_polynomial_upto_order(order):
    if order == 0:
        return poly[:terms_0]
    if order == 1:
        return poly[:terms_1]
    if order == 2:
        return poly[:terms_2]
    if order == 3:
        return poly[:terms_3]
    if order == 4:
        return poly[:terms_4]
    else:
        return "order can't be more than 4"
    
# this generate polynomial library 
def get_poly_lib(order,coords):
    poly_f = get_polynomial_upto_order(order)
    poly_f = sp.lambdify([x,y,z],poly_f)
    poly_fv = np.array(list(map(poly_f,*coords.T)))
    poly_fv = np.squeeze(poly_fv,axis=-1)
    return poly_fv

'''
Bx_f = a_0*1 + a_1*x + a_2*y + a_3*z + a_4*x**2 + a_5*x*y + a_6*x*z + a_7*y**2 + a_8*y*z + a_9*z**2 + a_10*x**3 + a_11*x**2*y + a_12*x**2*z +\
       a_13*x*y**2 + a_14*x*y*z + a_15*x*z**2 + a_16*y**3 + a_17*y**2*z + a_18*y*z**2 + a_19*z**3 + a_20*x**4 + a_21*x**3*y + a_22*x**3*z + a_23*x**2*y**2 +\
       a_24*x**2*y*z + a_25*x**2*z**2 + a_26*x*y**3 + a_27*x*y**2*z + a_28*x*y*z**2 + a_29*x*z**3 + a_30*y**4 + a_31*y**3*z + a_32*y**2*z**2 + a_33*y*z**3 + a_34*z**4
       
By_f = 'b_0*1 + b_1*x + b_2*y + b_3*z + b_4*x**2 + b_5*x*y + b_6*x*z + b_7*y**2 + b_8*y*z + b_9*z**2 + b_10*x**3 + b_11*x**2*y + b_12*x**2*z +\
    b_13*x*y**2 + b_14*x*y*z + b_15*x*z**2 + b_16*y**3 + b_17*y**2*z + b_18*y*z**2 + b_19*z**3 + b_20*x**4 + b_21*x**3*y + b_22*x**3*z + b_23*x**2*y**2 +\
    b_24*x**2*y*z + b_25*x**2*z**2 + b_26*x*y**3 + b_27*x*y**2*z + b_28*x*y*z**2 + b_29*x*z**3 + b_30*y**4 + b_31*y**3*z + b_32*y**2*z**2 + b_33*y*z**3 + b_34*z**4'
    
Bz_f = 'c_0*1 + c_1*x + c_2*y + c_3*z + c_4*x**2 + c_5*x*y + c_6*x*z + c_7*y**2 + c_8*y*z + c_9*z**2 + c_10*x**3 + c_11*x**2*y + c_12*x**2*z +\
    c_13*x*y**2 + c_14*x*y*z + c_15*x*z**2 + c_16*y**3 + c_17*y**2*z + c_18*y*z**2 + c_19*z**3 + c_20*x**4 + c_21*x**3*y + c_22*x**3*z + c_23*x**2*y**2 +\
    c_24*x**2*y*z + c_25*x**2*z**2 + c_26*x*y**3 + c_27*x*y**2*z + c_28*x*y*z**2 + c_29*x*z**3 + c_30*y**4 + c_31*y**3*z + c_32*y**2*z**2 + c_33*y*z**3 + c_34*z**4'

'''