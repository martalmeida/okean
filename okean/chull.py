import numpy as np
import pylab as pl
from okean import calc


def _distances(x,y):
  foundJ=[]
  dst=[]
  for i in range(len(x)):
    if i in foundJ: continue
    d=np.sqrt((x[i]-x)**2+(y[i]-y)**2)
    d[i]=1e30
    j=np.where(d==d.min())[0][0]
    foundJ+=[j]
    dst+=[d[j]]

  print('avg min dist: %f'%np.mean(dst))
  return np.asarray(dst)

def _distances_kdt(x,y):
  from scipy import spatial
  points=np.vstack((x,y)).T
  tree = spatial.cKDTree(points)

  foundJ=[]
  dst=[]
  for i in range(len(x)):
    if i in foundJ: continue
    dist,ind=tree.query((x[i],y[i]),2)
    foundJ+=[ind[1]]
    dst+=[dist[1]]

  print('avg min dist: %f'%np.mean(dst))
  return np.asarray(dst)


def distances(x,y,kdtree=False):
  if kdtree: return _distances(x,y)
  else: return _distances_kdt(x,y) # which is faster!?


def clockwise_angle(x,y,r=None,r0=None):
  # angle between two vectors defined by three points,
  #(x0,y0),(x1,y1),(x2,y2)

  x1=x[0]-x[1]
  x2=x[2]-x[1]

  y1=y[0]-y[1]
  y2=y[2]-y[1]

  dot = x1*x2 + y1*y2
  det = x1*y2 - y1*x2
  angle = np.arctan2(det, dot)

  if angle<0: angle=-angle
  else: angle=2*np.pi-angle

  if r and r==r0: # may not make sense if r~=r0
    L=np.sqrt(x2**2+y2**2)
    a=np.arccos(L/r)
    angle=angle-a

    L0=np.sqrt(x1**2+y1**2)
    if L0<r0:
      a0=np.arccos(L0/r0)
      angle=angle-a0

  return angle


def find_next(X,Y,x0,y0,xi,yi,x,y,rmax,rmax_prev,method='circ',**kargs):
  tree=kargs.get('kdt',None)
  if method=='kdt':
    dist,ind=tree.query((xi,yi),rmax)
    ind=ind[1:]
    cond=np.zeros(x.shape,'bool')
    x=x[ind]
    y=y[ind]
    cond[ind]=1
  else:
    r=np.sqrt((x-xi)**2+(y-yi)**2)
    cond=(r<=rmax)&(r>0) # exclude xi,yi
    x,y=x[cond],y[cond]

  if len(x)==0:
    return None,None,None

  cond2=None
  if len(X)>2: # check intersections:
    cond2=np.ones(x.size,'bool')
    for i in range(x.size):
      xm,ym=calc.meetpoint(np.asarray(X[:-1]),np.asarray(Y[:-1]),
                           np.asarray([xi,x[i]]),np.asarray([yi,y[i]]))

      # forget about intersections at nodes of X,Y:
      if len(xm):
        # xm, ym may includes xi, yi. Will include at leats in the last step!
        cnd=np.asarray([(xi,yi)==(xm[j],ym[j]) for j in range(len(xm))])
        xm=xm[~cnd]
        ym=ym[~cnd]
        if len(xm):
          if method=='kdt': # allow 1st element only (so that polygon can be closed!)
            isXY=[(xm[j],ym[j])==(X[0],Y[0]) for j in range(len(xm))]
          else:
            isXY=[xm[j]+ym[j]*1j in np.asarray(X[:-1])+np.asarray(Y[:-1])*1j for j in range(len(xm))]

          isXY=np.asarray(isXY)
          xm=xm[~isXY]
          ym=ym[~isXY]

      cond2[i]=len(xm)==0

    x,y=x[cond2],y[cond2]

  if len(x)==0:
    return None,None,None

  if method=='circ':
    ang=np.asarray([clockwise_angle((x0,xi,x[i]),(y0,yi,y[i]),rmax,rmax_prev) for i in range(x.size)])
  else:
    ang=np.asarray([clockwise_angle((x0,xi,x[i]),(y0,yi,y[i]),0) for i in range(x.size)])

  # must be carefull with angles very close to zeros (colinear points)
  # better to use an epsilon as min angle:
  eps=np.finfo(x.dtype).eps
  eps=eps*1e4
  ang[ang<eps]=2*np.pi
  # also  be carefull about very close angles!! again the colinear points ...
  # let us assume ang is a float 64
  if ang.dtype==np.dtype('float64'):
    ang=ang.astype('float32')
  else:
    print('warning chull: ang is not float 64 -- possible colinear points problem !')

  i=np.where(ang==ang.min())[0]
  if len(ang)==1:
    i=i[0]
  else: # use nearest point: may solve some problems with colinear points...
    d=(x[i]-xi)**2+(y[i]-yi)**2
    j=np.where(d==d.min())[0][0]
    i=i[j]

  if cond2 is None:
    ii=np.where(cond)[0][i]
  else:
    j=np.where(cond2)[0][i]
    ii=np.where(cond)[0][j]

  return x[i],y[i],ii


def _chull(x,y,rmax,maxr='auto',method='circ'):
  status=0
  X=[]
  Y=[]
  segs=[]

  if method=='kdt' and rmax>len(x):
    print(rmax,'Max k is %d'%len(x))
    return X,Y,1

  if method=='kdt':
    from scipy import spatial
    points=np.vstack((x,y)).T
    tree = spatial.cKDTree(points)
  else: tree=None

  # 1st element:
  i=np.where(y==y.max())[0][0]
  X+=[x[i]]
  Y+=[y[i]]


  stop=False
  for j in range(x.size*2):
    if j==0:
      x0,y0=x.min()-1,Y[-1]
      #rprev=0
      rprev=rmax
    else:
      x0,y0=X[-2],Y[-2]

    xi,yi=X[-1],Y[-1]
    xx,yy,i=find_next(X,Y,x0,y0,xi,yi,x,y,rmax=rmax,rmax_prev=rprev,method=method,kdt=tree)
    rprev=rmax
    rmax2=rmax

    if method!='kdt':
      inverted_path=(xx==x0)&(yy==y0)
      if inverted_path: print('inverted...')
      while (xx is None) or (maxr and inverted_path):
        print('- looking for next inside while',rmax2)
        rmax2=rmax2*1.1

        if rmax2<=maxr:
          xx,yy,i=find_next(X,Y,x0,y0,xi,yi,x,y,rmax=rmax2,rmax_prev=rprev,method=method)
          rprev=rmax2
          inverted_path=(xx==x0)&(yy==y0)
        else:
          break

    if xx is None:
      print('cannot locate next point!')
      status=1
      break


    X+=[xx]
    Y+=[yy]

    if len(X)>2:
      if method=='kdt':
        is_closed=(xx,yy)==(X[0],Y[0])
        if is_closed:
          X=X[:-1]
          Y=Y[:-1]

      else:
        is_closed=(xi,yi,xx,yy)==(X[0],Y[0],X[1],Y[1])
        if is_closed:
          X=X[:-2]
          Y=Y[:-2]

      if is_closed:
        print('DONE: polygon is closed!')
        print('AREA=',calc.poly_area(np.asarray(X),np.asarray(Y)))
        break

    # current segment:
    seg=X[-2],X[-1],Y[-2],Y[-1]
    if seg in segs:
      print('ERROR: repeating segment !')
      status=1
      break

    segs+=[seg]


  # close polygon
  if not status:
    X=np.hstack((X,X[0]))
    Y=np.hstack((Y,Y[0]))

  return X,Y,status

def chull(x,y,rmax,maxr='auto',method='circ'):

  if method=='kdt':
    if maxr=='auto': maxr=rmax*10
  else:
    if maxr=='auto':
      maxr=np.sqrt((x.max()-x.min())**2+(y.max()-y.min())**2)

  status=1
  while status:
    xx,yy,status=_chull(x,y,rmax,maxr,method)
    if method=='kdt': rmax=rmax+1
    else: rmax=rmax*1.05

    if status and rmax>maxr:
      print('maxr reached!')
      break

  return xx,yy, status

