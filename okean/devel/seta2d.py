import pylab as pl
import numpy as np
from okean import calc


x0,y0=0,0
u,v=1,0.5
#u,v=0.05,0.05

#v=-0.5
#u=-1


#w=s/5
#w=min(w,s/3.)


def seta(anchor='start'):
  # anchor: start or center

  t1=35*np.pi/180
  t2=80*np.pi/180 # t0..90

  w=0.1
  h=0.2

  s=(u**2+v**2)**0.5
  h=min(h,s/3)
  w=min(w,s/5)

  P0=np.array([0,0])
  P1=np.array([0,w/2])
  w1=np.tan(t1)*h-w/2
  h1=w1/np.tan(t2)
  P2=np.array([s-h+h1,w/2])
  P3=np.array([s-h,w/2+w1])
  P4=np.array([s,0])

  P5=np.array([P3[0],-P3[1]])
  P6=np.array([P2[0],-P2[1]])
  P7=np.array([P1[0],-P1[1]])

  x=P1[0],P2[0],P3[0],P4[0],P5[0],P6[0],P7[0],P1[0]
  y=P1[1],P2[1],P3[1],P4[1],P5[1],P6[1],P7[1],P1[1]
  x=np.asarray(x)
  y=np.asarray(y)
  ang=np.arctan2(v,u)
  x,y=calc.rot2d(x,y,-ang)
  pl.plot(x,y)

u,v=2,2.5
seta()

u,v=1,0.5
seta()

u,v=.5,.5
seta()

u,v=.1,.1
seta()
u,v=0.01,0.01
seta()
pl.axis('equal')

