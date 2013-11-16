import numpy as np
import rtools_test

theta_s=.2
theta_b=.2
N=10

sc_r,Cs_r,sc_w,Cs_w=rtools_test.scoord(theta_s,theta_b,N)
#print sc_r,Cs_r,sc_w,Cs_w

h=np.zeros(1,'f')+1000
zeta=0*h
hc=10

zr,zw=rtools_test.s_levels(h,zeta,hc,theta_s,theta_b,N)

print zr
print zw
