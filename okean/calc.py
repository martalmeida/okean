import numpy as np
from matplotlib import delaunay


def isarray(v,nomask=False):
  '''
  True for numpy arrays or numpy masked arrays (if nomask is False)
  '''

  if nomask: return isinstance(v,np.ndarray)
  else: return isinstance(v,np.ndarray) or np.ma.isMA(v)

def ismarray(v):
  '''
  True for numpy numpy masked arrays
  '''
  return np.ma.isMA(v)


def isiterable(*args):
  '''
  True for sequences
  '''

  try:
    for a in args: iter(a)
    return True
  except: return False


def inpolygon(x,y,xp,yp):
  '''
  Points inside polygon test.

  Based on matplotlib nxutils for old matplotlib versions
  http://matplotlib.sourceforge.net/faq/howto_faq.html#test-whether-a-point-is-inside-a-polygon

  Can also use Path.contains_point, for recent versions, or
  the pnpoly extension - point in polyon test - by W R Franklin,
  http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html
  '''

  try:
    from matplotlib.nxutils import pnpoly, points_inside_poly
    _use_mpl = 1
  except ImportError:
    from matplotlib.path import Path
    _use_mpl = 2
  except:
    import pnpoly
    _use_mpl = False

  shape=False
  try:
    if x.ndim>1:
      shape=x.shape
      x=x.flat
      y=y.flat
  except: pass

  if _use_mpl==1:
    verts=np.array(zip(xp,yp))
    if isiterable(x,y):
      points=np.array(zip(x,y))
      res=points_inside_poly(points,verts)
    else:
      res=pnpoly(x,y,verts)==1

  elif _use_mpl==2:
    verts=np.array(zip(xp,yp))
    p=Path(verts)
    if not isiterable(x,y): x,y,itr=[x],[y],0
    else: itr=1
    res=[p.contains_point((x[i],y[i]), radius=0.) for i in range(len(x))]
    res=np.asarray(res,'bool')
    if not itr: res=res[0]

  else:
    if not isiterable(x,y): x,y,itr=[x],[y],0
    else: itr=1
    res=np.zeros(len(x))
    res=[pnpoly.pnpoly(x[i],y[i],xp,yp,len(xp)) for i in range(len(x))]
    res=np.asarray(res)==1
    if not itr: res=res[0]

  if not shape is False: res.shape=shape
  return res


def inline(P,Q,T):
  '''
  2d point-in-line test
  After A.W. Paeth, Graphics Gems

  Given two points P and Q and a test point T
  return:
    0 if T is not on the (infinite) line PQ
    1 if T is on the open ray P
    2 if T is within the line segment PQ
    3 if T is on the open ray Q
  '''

  if abs((Q[1]-P[1])*(T[0]-P[0])-(T[1]-P[1])*(Q[0]-P[0]))>=max(abs(Q[0]-P[0]),
                                                               abs(Q[1]-P[1])): return 0
  elif (Q[0]<P[0] and P[0]<T[0]) or (Q[1]<P[1] and P[1]<T[1]): return 1
  elif (T[0]<P[0] and P[0]<Q[0]) or (T[1]<P[1] and P[1]<Q[1]): return 1
  elif (P[0]<Q[0] and Q[0]<T[0]) or (P[1]<Q[1] and Q[1]<T[1]): return 3
  elif (T[0]<Q[0] and Q[0]<P[0]) or (T[1]<Q[1] and Q[1]<P[1]): return 3
  else: return 2


def poly_area(x,y):
  '''
  Signed area of polygon
  If positive, the polygon is counter-clockwise (direct)
  The polygon need not to be closed
  '''
  r=range(1,x.size)+[0]
  return 0.5*np.sum(x*y[r]-y*x[r])


def poly_centroid(x,y):
  '''
  Geometric centre (Centroid) of non-self-intersecting polygon
  The polygon doesn't need to be closed
  '''
  r=range(1,x.size)+[0]
  a=poly_area(x,y)
  xc=1./(6.*a)*np.sum((x+x[r])*(x*y[r]-y*x[r]))
  yc=1./(6.*a)*np.sum((y+y[r])*(x*y[r]-y*x[r]))
  return xc,yc


def griddata(x,y,v,xi,yi,**kargs):
  '''
  Interpolates scattered or gridded (2d, regular or irregular) data
  (x,y,v) to some set of scattered or gridded points (xi,yi),

  Extrapolation:
  If karg extrap is True, extrapolation is done outside convex hull
  (of original data without the masked points).

  If extrap is False (default):

    - all points outside convex hull are returned as mask.

    - in case of gridded original data points (2d x,y,v) all points
      outside the domain boundary are returned as mask.

    - original masked regions, if present, are kept. Some masked
      regions may be inside the non masked convex hull and
      interpolation is done for such points. To avoid this the
      original mask is interpolated on the target points and where
      this new mask > some value, the target data is masked.
      The default value is 0.5 and can be changed with the keyword
      argument keepMaskVal (0 is the most strick case).


  Based on delaunay triangulation provided by matplotlib.
  Whenever triangulation is not possible numpy nans are returned

  Example:
    import numpy as np
    import pylab
    x  = np.linspace(-2,2,20)
    y  = x
    x,y=np.meshgrid(x,y)
    v  = x*np.exp(-x**2-y**2)
    xi = np.linspace(-1,3,100)
    yi = xi
    xi,yi=np.meshgrid(xi,yi)
    vi=griddata(x,y,v,xi,yi)
    pylab.figure()
    pylab.pcolormesh(x,y,v,shading='faceted')
    pylab.pcolormesh(xi,yi,vi)

  mma  2010
  '''

  mask2d = False
  # When v.ndim>2 and the same mask is used for all the dim 0 levels
  # setting mask2d is mutch fater since the triangulation is done only
  # once

  extrap=False
  # extrapolate outside data points convex hull and in masked region

  # keepMask    = not extrap
  keepMaskVal = 0.5
  # Notice  that even without extrap, the convex hull of the
  # non masked data may include masked points. In that case
  # interpolation is also done for such points!!
  # The solution is to interpolate the mask (as int, 0 or 1 where masked)
  # and consider mask where this value is greater than 0 in the most
  # strict case; 0.5 seems the 'natural' case; to explain, in 1d, it
  # means target points between on non masked cell and one masked cell
  # are not to be considered mask until half distance between then.
  # In case of extrapolation it makes no sense to keepMask

  # forceBoundary = not extrap
  #
  # When original data points represent a 2d grid, it is possible to
  # know the boundary of the points, which usually differs from the
  # convex hull (except when the sides of the region are straight
  # lines). We can then calculate the target points outside that boundary
  # and set then as mask.
  # forceBoundary only makes sense when extrap is False

  # finalextrap = extrap
  #
  # Even with extrap True, final data may have NaNs, depending on
  # te number and distribution of the original data.
  # In such case a second extrapolation is done to fill the NaN points
  # This extrapolation is done with the target points, so if for
  # instance they represent a line, the triangulation is not possible
  # and thus extrapolation is not done
  # finalextrap only makes sense when extrap is True

  norm_xy=True
  norm_xy=False
  # delaunay may give VERY VERY BAD bad results when x and y data range
  # differ by orders of magnitude. In such case the solution is to normalize
  # the coordinates

  if 'mask2d'        in kargs.keys(): mask2d        = kargs['mask2d']
  if 'extrap'        in kargs.keys(): extrap        = kargs['extrap']

  if extrap:
    keepMask      = False
    forceBoundary = False
    finalextrap   = True
  else:
    keepMask      = True
    forceBoundary = True
    finalextrap   = False

  if 'keepMask'      in kargs.keys(): keepMask      = kargs['keepMask']
  if 'keepMaskVal'   in kargs.keys(): keepMaskVal   = kargs['keepMaskVal']
  if 'forceBoundary' in kargs.keys(): forceBoundary = kargs['forceBoundary']
  if 'finalextrap'   in kargs.keys(): finalextrap   = kargs['finalextrap']

  if 'norm_xy' in kargs.keys(): norm_xy=kargs['norm_xy']

  if norm_xy:
    rxy=(x.max()-x.min())/(y.max()-y.min())
    y=y*rxy
    yi=yi*rxy


  # interp/extrap:
  try:    res=_griddataz(x,y,v,xi,yi,mask2d,extrap)
  except: res=xi*np.nan

  # final extrap:
  if finalextrap:
    try: res=mask_extrap(xi,yi,np.ma.masked_where(np.isnan(res),res))
    except: pass

  maskCond=np.isnan(res)
  # keep original mask:
  if keepMask:
    imask=False
    # check if there is a mask!
    if mask2d is False:
      if np.ma.isMA(v) and np.ma.count_masked(v)>0: imask=v.mask.astype('int8')
    else: imask=mask2d.astype('int8')

    if not imask is False:
      km=_griddata(x,y,imask,xi,yi,extrap=False)[0]
      km=np.where(np.isnan(km),np.nan,km)
      maskCond=maskCond | (km>keepMaskVal)

  # force original boundary:
  if forceBoundary:
    if (x.ndim,y.ndim)==(2,2):
      # get original points boundary:
      xb=var_border(x)
      yb=var_border(y)

      # get target points outside:
      cond=~inpolygon(xi.flat,yi.flat,xb,yb)
      cond.shape=xi.shape
      maskCond=maskCond | cond

  return np.ma.masked_where(maskCond,res)


def _griddataz(x,y,v,xi,yi,mask2d,extrap):
  '''
  Use griddata instead
  '''
  if v.ndim==x.ndim+1:
    tmp,tri=_griddata(x,y,v[0,...],xi,yi,extrap,tri=False,mask=mask2d)
    res=np.zeros([v.shape[0]]+list(tmp.shape),dtype=v.dtype)
    res[0,...]=tmp
    for i in range(1,v.shape[0]):
      res[i,...]=_griddata(x,y,v[i,...],xi,yi,extrap,tri,mask=mask2d)[0]
  else:
    res=_griddata(x,y,v,xi,yi,extrap,mask=mask2d)[0]

  return res


def _griddata(x,y,v,xi,yi,extrap=True,tri=False,mask=False):
  '''
  Use griddata instead
  '''

  if mask is False:
    if np.ma.isMA(v) and np.ma.count_masked(v)>0: mask=v.mask
    else: mask=np.zeros(v.shape,'bool')

  if 0:
    if not tri:
      tri=delaunay.Triangulation(x[~mask],y[~mask])

    if extrap: u=tri.nn_extrapolator(v[~mask])(xi,yi)
    else:      u=tri.nn_interpolator(v[~mask])(xi,yi)
  else:
    # deal with repeated pairs (problem for nn_extrapolator)
    xy=x[~mask]+1j*y[~mask]
    xy,ii=np.unique(xy,1)
    if not tri:
      tri=delaunay.Triangulation(x[~mask][ii],y[~mask][ii])

    if extrap: u=tri.nn_extrapolator(v[~mask][ii])(xi,yi)
    else:      u=tri.nn_interpolator(v[~mask][ii])(xi,yi)

  return u, tri


def mask_extrap(x,y,v):
  '''
  Extrapolate numpy array at masked points.
  Based on delaunay triangulation provided by matplotlib.
  '''

  if np.ma.isMA(v) and v.size!=v.count(): mask=v.mask
  else: return v

  if 0:
    tri=delaunay.Triangulation(x[~mask],y[~mask])
    v[mask]=tri.nn_extrapolator(v[~mask])(x[mask],y[mask])
  else:
    # deal with repeated pairs (problem for nn_extrapolator)
    xy=x[~mask]+1j*y[~mask]
    xy,ii=np.unique(xy,1)

    tri=delaunay.Triangulation(x[~mask][ii],y[~mask][ii])
    v[mask]=tri.nn_extrapolator(v[~mask][ii])(x[mask],y[mask])

  if np.any(np.isnan(v)):
    mask=np.isnan(v)
    tri=delaunay.Triangulation(x[~mask],y[~mask])
    v[mask]=tri.nn_extrapolator(v[~mask])(x[mask],y[mask])

  return v


def cyclic_index(time,t,cycle):
  '''
  Index arround value in periodic sequence

  Inputs:
    time   Sequence
    t      Value
    cycle  Sequence period

  Outputs:
    Index of of time above and below t
    t - time below, time above - t

  Example:
    import numpy as np
    x=np.arange(15,365+30,30)
    xi=5
    cycle=365
    i,d = cyclic_index(x,xi,cycle) # i=[11 0], d=[25 10]

  MMA 16-07-2008
  Dep. Earth Physics, UFBA, Salvador, Bahia, Brasil
  '''

  # time must be <= cycle
  time=time[time<=cycle]

  t=np.mod(t,cycle)

  Inds=np.arange(len(time)+2)
  Inds[0]=len(time)
  Inds[-1]=1
  Inds=Inds-1

  Time=np.arange(len(time)+2)
  Time[0]=-(cycle-time[-1])
  Time[-1]=cycle+time[0]
  Time[1:-1]=time

  i2=np.where(Time>t)[0][0]
  i1=i2-1
  I2=Inds[i2]
  I1=Inds[i1]

  d1=t-Time[i1]
  d2=Time[i2]-t

  Ind=np.array([I1, I2])
  d  =np.array([d1, d2])

  return Ind,d


def distance(lon,lat):
  '''
  Distance along path on sphere

  Example:
    import numpy as np
    lon=np.zeros(5)
    lat=np.arange(5)
    distance(lon,lat)/1000 # 0,111.13, 222.27, ...

  mma
  '''

  d=spheric_dist(lat[1:],lat[:-1],lon[1:],lon[:-1])
  return np.append(0.,d.cumsum())


def spheric_dist(lat1,lat2,lon1,lon2):
  '''
  Compute distances for a simple spheric earth
  between two points

  mma
  Python version after P. Penven (ROMSTOOLS)
  '''

  R=6367442.76
  # Determine proper longitudinal shift.
  l=abs(lon2-lon1)
  l[l>=180]=360-l[l>=180]
  # Convert Decimal degrees to radians.
  deg2rad=np.pi/180.
  lat1=lat1*deg2rad
  lat2=lat2*deg2rad
  l=l*deg2rad
  # Compute the distances
  dist=R*np.arcsin(np.sqrt(((np.sin(l)*np.cos(lat2))**2)+(((np.sin(lat2)*\
       np.cos(lat1))-(np.sin(lat1)*np.cos(lat2)*np.cos(l)))**2)))

  return dist


def var_border(v,di=1,dj=1):
  '''
  Border of 2d numpy array
  di,dj is the interval between points along columns and lines
  Corner points are kept even with di and dj not 1
  '''

  j,i=v.shape
  if (di,dj)==(1,1):
    xb=np.arange(2*i+2*j,dtype=v.dtype)
    yb=np.arange(2*i+2*j,dtype=v.dtype)

    xb[0:j]       = v[:,0]
    xb[j:j+i]     = v[-1,:]
    xb[j+i:j+i+j] = np.flipud(v[:,-1])
    xb[j+i+j:]    = np.flipud(v[0,:])
  else:
    # ensure corner points are kept!!
    tmp1 = v[::dj,0]
    tmp2 = v[-1,::di]
    tmp3 = np.flipud(v[:,-1])[::dj]
    tmp4 = np.flipud(v[0,:])[::di]
    xb=np.concatenate((tmp1,tmp2,tmp3,tmp4))

  return xb


def mode(x):
  '''
  Most frequently occuring element of sequence or numpy array
  Also returns the number of repetitions of the element
  '''

  x=np.array(x).flatten()
  x.sort()
  dist=np.diff(np.append(x,x[-1]-1))
  idx=np.where(dist!=0)[0]
  num=np.arange(len(idx))
  num[0]=idx[0]+1
  num[1:]=np.diff(idx)
  n=np.max(num)
  return x[idx[num==n]],n


def rot2d(x,y,ang, inverse=False):
  '''
  2d rotation
  Rotate x,y by angle ang (rad)

  R =  cos f  -sin f
       sin f   cos f

  X =  x cos f  +  y sin f
  Y = -x sin f  +  y cos f

  inverse:
  R'

  X = x cos f  -  y sin f
  Y = x sin f  +  y cos f

  ang can be a scalar or an array with the shape of x and y
  '''

  return rot3d(x,y,0*x,ang,0.,inverse)[:2]


def rot3d(x,y,z,fi,teta,inverse=False):
  '''
  3D rotation with azimuth fi and inclination teta
  NB: teta is not the elevation from the reference plane
  like returned by matlab cart2sph, for instance

  Ry =  cos f  -sin f  0
        sin f   cos f  0
        0       0      1

  Rz =  cos t   0  sin t
        0       1  0
       -sin t   0  cos t

  [x y z]*Rz*Ry
  inverse:
  [x y z]*(Rz*Ry)'

  the angles (rad) can have have the shape of x,y,z or be a scalar
  '''

  if inverse:
    xx =  x*np.cos(fi)*np.cos(teta) -y*np.sin(fi) +z*np.cos(fi)*np.sin(teta)
    yy =  x*np.sin(fi)*np.cos(teta) +y*np.cos(fi) +z*np.sin(fi)*np.sin(teta)
    zz = -x*np.sin(teta)+z*np.cos(teta)
  else:
    xx =  x*np.cos(fi)*np.cos(teta)+y*np.sin(fi)-z*np.sin(teta)
    yy = -x*np.sin(fi)+y*np.cos(fi)
    zz =  x*np.cos(fi)*np.sin(teta)+y*np.sin(fi)*np.sin(teta)+z*np.cos(teta)

  return xx,yy,zz


def cart2sph(x,y,z):
  '''
  Cartesian coordinates to spherical polar coordinates
  Returns r,fi (azimuth),teta (inclination, not elevation, ie,
  not latitude)
  '''

  if not isarray(x):
    x=np.array([x])
    y=np.array([y])
    z=np.array([z])

  r=0.*x
  teta=0.*x
  fi=0.*x

  r=np.sqrt(x**2+y**2+z**2)
  teta[r!=0]=np.arccos(z[r!=0]/r[r!=0])

  rr = np.sqrt(x**2+y**2)
  c1=(y>=0) & (rr!=0)
  c2=(y<0) & (rr!=0)
  fi[c1] = np.arccos(x[c1]/rr[c1])         # y>=0 & rr!=0
  fi[c2] = 2*np.pi-np.arccos(x[c2]/rr[c2]) # y<0  & rr!=0

  return r,fi,teta


def sph2cart(r,fi,teta):
  '''
  Spherical polar coordinates to Cartesian coordinates
  Input angles are fi - azimuth,teta - inclination
  '''

  x=r*np.sin(teta)*np.cos(fi)
  y=r*np.sin(teta)*np.sin(fi)
  z=r*np.cos(teta)
  return x,y,z


def meetpoint(x1,y1,x2,y2):
  '''
  Returns the intesection points of the two curves x1,y1 and y1,y2
  Can deal with coincident and with vertical segments.

  mma IEO-A Coruna 2008
  '''

  from alg import meetpoint as meetp
  mask=999.
  xy=[]
  xyi,N=meetp(x1,y1,x2,y2,len(x1),len(x2))
  xi=xyi[0]
  yi=xyi[1]
  for i in range(N):
    if  (xi[i],yi[i]) not in xy and xi[i]!=mask:
       xy+=  [(xi[i],yi[i])]

  ty=x1.dtype
  return np.array([i for i,j in xy],ty),np.array([j for i,j in xy],ty)


def matrix_paths(x,y):
  '''
  Returns 1d arrays of the 2d arrays x,y with the odd columns and lines
  flipped, ie, returns xi,yi with the odd columns flipped and xj,yj with
  the odd lines flipped.
  See mcross_points for more information

  mma PONG-TAMU 2010
  '''

  def m2path(a,ij):
    b=a.copy()
    if   ij=='i': b[:,1::2]=np.flipud(a[:,1::2]) 
    elif ij=='j': b[1::2,:]=np.fliplr(a[1::2,:]) 
    return b

  xi=m2path(x,'i').T
  yi=m2path(y,'i').T

  xj=m2path(x,'j')
  yj=m2path(y,'j')

  return xi.flatten(),yi.flatten(),xj.flatten(),yj.flatten()

def mcross_points(xm,ym,xl,yl,addl=False,plt=False):
  '''
  Returns the intersection between a stright line (xl,yl two points)
  and a gridded 2d region defined by the 2d arrays xm, ym.
  If addl the line points are also added to the output x,y points.
  The output is ordered with the distance from first line point.

  The intersection points are the max number of meet points between
  the line and the two possible "matrix paths" (see matrix_paths).
  For this reason better results are obtained with regular or quasi
  regular gridded regions. If the region as a L shape, for instance,
  and the line is horizontal many points are obtained in the | but none
  in the _. The solution would be calculate the number of cells crossed
  by the line and then calculate one inetrsection... but would be slow!

  mma PONG-TAMU 2010
  '''

  xi,yi,xj,yj=matrix_paths(xm,ym)

  #xxi,yyi=meetpoint(xi,yi,xl,yl)
  #xxj,yyj=meetpoint(xj,yj,xl,yl)
  # faster:
  xxi,yyi=meetpoint(xl,yl,xi,yi)
  xxj,yyj=meetpoint(xl,yl,xj,yj)

  if len(xxi)>len(xxj): xx,yy=xxi,yyi
  else: xx,yy=xxj,yyj

  # add line points:
  if addl or 1:
    for i in range(len(xl)): # no need for the loop, only 2 points!
      if xl[i] not in xx or yl[i] not in yy:
        xx=np.concatenate((xx,[xl[i]]))
        yy=np.concatenate((yy,[yl[i]]))

  # sort with distance form first point:
  d=(xx-xl[0])**2+(yy-yl[0])**2
  D=d.copy()
  D.sort()
  inds=np.searchsorted(D,d)

  xx[inds]=xx.copy()
  yy[inds]=yy.copy()


  if plt:
    import pylab as plt
    pl.figure()
    ax1=pl.subplot(311)
    ax1.plot(xi,yi,'b')
    ax1.plot(xl,yl,'k')
    ax1.plot(xxi,yyi,'kx')

    ax2=pl.subplot(312)
    ax2.plot(xj,yj,'r')
    ax2.plot(xl,yl,'k')
    ax2.plot(xxj,yyj,'kx')

    ax3=pl.subplot(313)
    ax3.plot(xi,yi,'k')
    ax3.plot(xj,yj,'r')
    ax3.plot(xl,yl,'k')
    ax3.plot(xx,yy,'kx')

    #for i in range(len(xx)):
    #  pl.plot(xx[i],yy[i],'ro')
    #  raw_input()

  return xx,yy

def ij_limits(x,y,xlim,ylim,margin=0):
  '''
  ij limits of xy so that xlim and ylim box is inside xy[ij]
  Returns i0,i1,j0,j1, so for 1D x[i0:i1], y[j0:j1], and for 2D
  x[j0:j1,i0:i1], y[j0:j1,i0:i1] include xlim and ylim by margin+1
  lines/columns
  '''
  if x.ndim==1:
    L=y.size
    M=x.size
    i1,=np.where(x<xlim[0])
    if len(i1): i1=i1[-1]
    else: i1=0

    i2,=np.where(x>xlim[1])
    if len(i2): i2=i2[0]+1
    else: i2=M

    j1,=np.where(y<ylim[0])
    if len(j1): j1=j1[-1]
    else: j1=0

    j2,=np.where(y>ylim[1])
    if len(j2): j2=j2[0]+1
    else: j2=L

  elif x.ndim==2:
    L,M=x.shape

    xp=np.array([xlim[0],xlim[1],xlim[1],xlim[0]])
    yp=np.array([ylim[0],ylim[0],ylim[1],ylim[1]])

    inp=inpolygon(x,y,xp,yp)
    try:    j,i=np.where(inp)
    except: j,i=np.ma.where(inp)
    if i.size:
      i1,i2=i.min()-1,i.max()+1
      j1,j2=j.min()-1,j.max()+1
    else: return [False]*4


  i1=max(0,i1-margin)
  j1=max(0,j1-margin)
  i2=min(M,i2+margin)
  j2=min(L,j2+margin)

  return i1,i2,j1,j2


def angle_point_line(x,y,xp,yp,n=1):
  '''angle between point line.
     line segments are x[i-n] to x[i+1] for angle[i], i=n:len(x)-1
     n should be 1 or 0 !
  '''
  def calc_ang(i):
    u=x[i+1]-x[i-n],y[i+1]-y[i-n]
    v=xp-x[i],yp-y[i]
    dot=u[0]*v[0]+u[1]*v[1]
    return np.arccos(dot/(np.abs(u[0]+u[1]*1j)*np.abs(v[0]+v[1]*1j)))

  if n==1:
    return np.squeeze([calc_ang(i) for i in range(1,len(x)-1)])
  elif n==0:
    return np.squeeze([calc_ang(i) for i in range(len(x)-1)])

