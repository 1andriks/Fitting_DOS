import numpy as np 

# a=np.array([1,2,8,1,1,0,1])

# b=np.array([0,0,0,1,2,3,0])
# mask=b !=0
# c=a[mask] -b[mask]
# print(c)

COEFF=np.load('coefficients.npz')
print(COEFF[24][:,1])