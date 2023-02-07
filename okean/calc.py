import numpy as np

def isarray(v,nomask=False):
  'True for numpy arrays or numpy masked arrays (if nomask is False)'
  if nomask: return isinstance(v,np.ndarray)
  else: return isinstance(v,np.ndarray) or np.ma.isMA(v)

def ismarray(v):
  'True for numpy numpy masked arrays'
  return np.ma.isMA(v)


def isiterable(*args):
  'True for sequences'
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
      # if masked will generate a MaskedIterator... not ok as we want a flatiter
      # so:
      if np.ma.isMA(x): x=x.data.flat
      else: x=x.flat

      if np.ma.isMA(y): y=y.data.flat
      else: y=y.flat
  except: pass

  if _use_mpl==1:
    verts=np.array(list(zip(xp,yp)))
    if isiterable(x,y):
      points=np.array(list(zip(x,y)))
      res=points_inside_poly(points,verts)
    else:
      res=pnpoly(x,y,verts)==1

  elif _use_mpl==2:
    verts=np.array(list(zip(xp,yp)))
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
  r=list(range(1,x.size))+[0]
  return 0.5*np.sum(x*y[r]-y*x[r])


def poly_centroid(x,y):
  '''
  Geometric centre (Centroid) of non-self-intersecting polygon
  The polygon doesn't need to be closed
  '''
  r=list(range(1,x.size))+[0]
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

  mask2d=kargs.pop('mask2d',False)
  # When v.ndim>2 and the same mask is used for all the dim 0 levels
  # setting mask2d is mutch fater since the triangulation is done only
  # once

  extrap=kargs.pop('extrap',False)
  # extrapolate outside data points convex hull and in masked region

  # keepMask    = not extrap
  #
  keepMaskVal=kargs.pop('keepMaskVal',0.5)
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

  norm_xy=kargs.pop('norm_xy',False)
  # matplotlib.delaunay may give VERY VERY BAD bad results when x and y data
  # range differ by orders of magnitude. In such case the solution is to
  # normalize the coordinates


  if extrap:
    keepMask      = False
    forceBoundary = False
    finalextrap   = True
  else:
    keepMask      = True
    forceBoundary = True
    finalextrap   = False

  keepMask      = kargs.pop('keepMask',      keepMask)
  forceBoundary = kargs.pop('forceBoundary', forceBoundary)
  finalextrap   = kargs.pop('finalextrap',   finalextrap)


  if norm_xy:
    dx=1.*(x.max()-x.min())
    dy=1.*(y.max()-y.min())
    x=100*(x-x.min())/dx
    y=100*(y-y.min())/dy


  # interp/extrap:
  res=_griddataz(x,y,v,xi,yi,mask2d,extrap,**kargs)
#  try:    res=_griddataz(x,y,v,xi,yi,mask2d,extrap,**kargs)
#  except: res=xi*np.nan

  # final extrap:
  if finalextrap:
    try: mask_extrap(xi,yi,np.ma.masked_where(np.isnan(res),res),inplace=True,**kargs)
    except: pass

  maskCond=np.isnan(res)
  # keep original mask:
  if keepMask:
    imask=False
    # check if there is a mask!
    if mask2d is False:
      if np.ma.isMA(v) and np.ma.count_masked(v)>0:
        if v.ndim==3: imask=v[-1].mask.astype('int8')
        else: imask=v.mask.astype('int8')
    else: imask=mask2d.astype('int8')

    if not imask is False:
      #print 'griddata for mask:'
      km=_griddata(x,y,imask,xi,yi,extrap=False,**kargs)[0]
      km=np.where(np.isnan(km),1,km)
      # note that km ndim may be lower than maskCond if v is 3d!
      maskCond=maskCond | (km>keepMaskVal)

  # force original boundary:
  if forceBoundary:
    if (x.ndim,y.ndim)==(2,2):
      # get original points boundary:
      xb=var_border(x)
      yb=var_border(y)

      # get target points outside:
      cond=~inpolygon(xi,yi,xb,yb)
      maskCond=maskCond | cond

  return np.ma.masked_where(maskCond,res)


def _griddataz(x,y,v,xi,yi,mask2d,extrap,**kargs):
  '''
  Use griddata instead
  '''

  mpl_tri=kargs.get('mpl_tri',True)

  if not mpl_tri:
    # nn_extrapolator/nn_interpolator may have problems dealing with
    # regular grids, nans may be obtained. One simple solution is to
    # rotate the domains!
    x,y=rot2d(np.squeeze(x),np.squeeze(y),np.pi/3.)
    xi,yi=rot2d(xi,yi,np.pi/3.)

  if v.ndim==x.ndim+1:
    tmp,tri=_griddata(x,y,v[0,...],xi,yi,extrap,tri=False,mask=mask2d,**kargs)
    res=np.zeros([v.shape[0]]+list(tmp.shape),dtype=v.dtype)
    res[0,...]=tmp
    for i in range(1,v.shape[0]):
      res[i,...]=_griddata(x,y,v[i,...],xi,yi,extrap,tri,mask=mask2d,**kargs)[0]
  else:
    res=_griddata(x,y,v,xi,yi,extrap,mask=mask2d,**kargs)[0]

  return res


def _griddata(x,y,v,xi,yi,extrap=True,tri=False,mask=False,**kargs):
  '''
  Use griddata instead
  '''

  mpl_tri=kargs.get('mpl_tri',True)
  tri_type=kargs.get('tri_type','cubic') # cubic or linear
  tri_kind=kargs.get('tri_kind','geom') # min_E or geom (for type cubic only)
  min_area=kargs.get('min_area','auto1') # right value of first bin if only 1st and
                                         # last bins have values; otherwise
                                         # area.mean()/1e10 (auto2)


  # warning, if x.shape=n,1  x[~mask] will also have 2 dims!! Thus better just use ravel...
  if x.shape!=x.size or y.shape!=y.size or v.shape!=v.size:
    x=x.ravel()
    y=y.ravel()
    v=v.ravel()
    if not mask is False: mask=mask.ravel()

  if mask is False:
    if np.ma.isMA(v) and np.ma.count_masked(v)>0: mask=v.mask
    else: mask=np.zeros(v.shape,'bool')


  if not mpl_tri:
    from matplotlib import delaunay # deprecated in version 1.4
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

      if extrap: u=tri.nn_extrapolator(v[~mask][ii].data)(xi,yi)
      else:      u=tri.nn_interpolator(v[~mask][ii])(xi,yi)

  else:
    import  matplotlib.tri as mtri
    if not tri:
      if extrap:
        # check if new points are needed in order to ensure extrap,
        # ie, check if xi,yi are not inside x,y convex hull:
        from shapely.geometry import MultiPoint
        chull=MultiPoint(np.vstack((x[~mask],y[~mask])).T).convex_hull
        # chull coords are given by xc,yc=chull.exterior.coords.xy
        points=MultiPoint(np.vstack((xi.ravel(),yi.ravel())).T)
        if not chull.contains(points):
          # add 4 corners to x,y:
          dx=xi.max()-xi.min()
          dy=yi.max()-yi.min()

          xv=np.asarray([xi.min()-dx,xi.max()+dy,xi.max()+dx,xi.min()-dx])
          yv=np.asarray([yi.min()-dy,yi.min()-dy,yi.max()+dy,yi.max()+dy])
          vv=np.zeros(4,v.dtype)
          mv=np.zeros(4,'bool')

          for i in range(4):
            d=(x[~mask]-xv[i])**2+(y[~mask]-yv[i])**2
            j=np.where(d==d.min())[0][0]
            vv[i]=v[~mask][j]

          x=np.ma.hstack((x,xv))
          y=np.ma.hstack((y,yv))
          v=np.ma.hstack((v,vv))
          mask=np.hstack((mask,mv))


      tri=mtri.Triangulation(x[~mask],y[~mask])

      # remove very small triangles:
      area=np.abs(np.asarray([poly_area(tri.x[i],tri.y[i]) for i in tri.triangles]))

      if min_area=='auto1':
        # check if all zeros except 1st and last, but 1st bin must
        # be smaller than last as it should contain outliers. A couple
        # of large triangles may appear in concave shapes and in such case
        # there are no outliers
        hist_y,hist_x=np.histogram(area)
        if np.all(hist_y[1:-1]==0) and hist_y[-1]>hist_y[0]: min_area=hist_x[1]
        else: min_area='auto2'

      if min_area=='auto2':
        tri.set_mask(area<area.mean()/1e10)
      else:
        tri.set_mask(area<min_area)

    else: # tri is provided as input arg
      tri,v,mask=tri

    if tri_type=='cubic':
      u = mtri.CubicTriInterpolator(tri, v[~mask],kind=tri_kind)(xi,yi)
    elif tri_type=='linear':
      u = mtri.LinearTriInterpolator(tri, v[~mask])(xi,yi)


  return u, (tri,v,mask) # new v and mask may be needed if extrap;
                         # returning v and mask will avoid recalculate
                         # convex hull and add new data corners


def mask_extrap(x,y,v,inplace=False,norm_xy=False,mpl_tri=True):
  '''
  Extrapolate numpy array at masked points.
  Based on delaunay triangulation provided by matplotlib.
  '''

  #if np.ma.isMA(v) and v.size!=v.count(): mask=v.mask
  if np.ma.is_masked(v): mask=v.mask
  else: return v

  sh=v.shape
  x=x.ravel()
  y=y.ravel()
  v=v.ravel()
  mask=mask.ravel()

  if inplace: u=v
  else: u=v.copy()

  if norm_xy:
    rxy=(x.max()-x.min())/(y.max()-y.min())
    y=y*rxy

  if not mpl_tri:
    from matplotlib import delaunay # deprecated in version 1.4

    # nn_extrapolator may have problems dealing with regular grids,
    # nans may be obtained. One simple solution is to rotate the domain!
    x,y=rot2d(x,y,np.pi/3.)

    if 0:
      tri=delaunay.Triangulation(x[~mask],y[~mask])
      u[mask]=tri.nn_extrapolator(u[~mask])(x[mask],y[mask])
    else:
      # deal with repeated pairs (problem for nn_extrapolator)
      xy=x[~mask]+1j*y[~mask]
      xy,ii=np.unique(xy,1)

      tri=delaunay.Triangulation(x[~mask][ii],y[~mask][ii])
      u[mask]=tri.nn_extrapolator(u[~mask][ii])(x[mask],y[mask])

    if np.any(np.isnan(u)):
      mask=np.isnan(u)
      tri=delaunay.Triangulation(x[~mask],y[~mask])
      u[mask]=tri.nn_extrapolator(u[~mask])(x[mask],y[mask])

  else:
    import  matplotlib.tri as mtri

    # add corners:
    xv=np.asarray([x.min()-1,x.max()+1,x.max()+1,x.min()-1])
    yv=np.asarray([y.min()-1,y.min()-1,y.max()+1,y.max()+1])
    vv=np.zeros(4,v.dtype)
    mv=np.zeros(4,'bool')

    for i in range(4):
      d=(x[~mask]-xv[i])**2+(y[~mask]-yv[i])**2
      j=np.where(d==d.min())[0][0]
      vv[i]=u[~mask][j]

    #x=np.ma.hstack((x.flat,xv)) # use ravel at top instead!
    x=np.ma.hstack((x,xv))
    y=np.ma.hstack((y,yv))
    u=np.ma.hstack((u,vv))
    mask=np.hstack((mask,mv))

    tri=mtri.Triangulation(x[~mask],y[~mask])
    #print u.shape,x.shape,y.shape,mask.shape
    u[mask] = mtri.CubicTriInterpolator(tri, u[~mask])(x[mask],y[mask])
    if np.any(np.isnan(u)):
      mask=np.isnan(u)
      tri=mtri.Triangulation(x[~mask],y[~mask])
      u[mask]=mtri.CubicTriInterpolator(tri,u[~mask])(x[mask],y[mask])

    u=u[:-4]

  u.shape=sh

  if not inplace:
    if not np.ma.is_masked(u): u=u.data
    return u


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

def lat_radius(lat):
  '''earth radius at lat'''
  lat=lat*np.pi/180
  r1=6378137
  r2=6356752.3

  r=(((r1**2*np.cos(lat))**2+(r2**2*np.sin(lat))**2)/((r1*np.cos(lat))**2+(r2*np.sin(lat))**2))**.5
  return r

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

  from .alg import meetpoint as meetp
  xy=[]
  xyi,N=meetp(x1,y1,x2,y2,len(x1),len(x2))
  return xyi[0,:N],xyi[1,:N]


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


def ij_limits_OLD_REMOVE(x,y,xlim,ylim,margin=0):
  '''
  ij limits of xy so that xlim and ylim box is inside xy[ij]
  Returns i0,i1,j0,j1, so for 1D x[i0:i1], y[j0:j1], and for 2D
  x[j0:j1,i0:i1], y[j0:j1,i0:i1] include xlim and ylim by margin+1
  lines/columns

  Ex:
  x,y=np.meshgrid(np.arange(10),np.arange(8))
  xl=2.5,5.5
  yl=2.5,4.5
  i0,i1,j0,j1=calc.ij_limits(x,y,xl,yl)

  pl.plot(x,y,'bo');
  pl.plot(xl,yl,'r*')
  pl.plot(x[j0:j1,i0:i1],y[j0:j1,i0:i1],'rs')
  '''

  if x.ndim==1:
    L=y.size
    M=x.size
    i1,=np.where(x<xlim[0])
    if len(i1): i1=i1[-1]
    else: i1=0

    i2,=np.where(x>xlim[1])
    if len(i2): i2=i2[0]
    else: i2=M

    j1,=np.where(y<ylim[0])
    if len(j1): j1=j1[-1]
    else: j1=0

    j2,=np.where(y>ylim[1])
    if len(j2): j2=j2[0]
    else: j2=L

  elif x.ndim==2:
    L,M=x.shape

    xp=np.array([xlim[0],xlim[1],xlim[1],xlim[0]])
    yp=np.array([ylim[0],ylim[0],ylim[1],ylim[1]])

    dxp=xp.max()-xp.min()
    dyp=yp.max()-yp.min()

    if dxp==0 or dyp==0: # vertical or horizontal line
      dx0=np.diff(np.abs(x)).max()
      dx1=np.diff(np.abs(x).T).max()
      dy0=np.diff(np.abs(y)).max()
      dy1=np.diff(np.abs(y).T).max()
      dxy=np.max([dx0,dx1,dy0,dy1])

      if dxp==0:
        xp[:2]=xp[0]-dxy
        xp[-2:]=xp[-1]+dxy

      if dyp==0:
        yp[:2]=yp[0]-dxy
        yp[-2:]=yp[-1]+dxy


    inp=inpolygon(x,y,xp,yp)
    try:    j,i=np.where(inp)
    except: j,i=np.ma.where(inp)
    if i.size:
      i1,i2=i.min()-1,i.max()+1
      j1,j2=j.min()-1,j.max()+1
    else: return [False]*4


  i1=max(0,i1-margin)
  j1=max(0,j1-margin)
  i2=min(M-1,i2+margin)
  j2=min(L-1,j2+margin)

  return i1,i2+1,j1,j2+1


def ij_limits(x,y,xlim,ylim,margin=1,check=0):

  if x.ndim==1:
    L=y.size
    M=x.size
    i1,=np.where(x<xlim[0])
    if len(i1): i1=i1[-1]
    else: i1=0

    i2,=np.where(x>xlim[1])
    if len(i2): i2=i2[0]
    else: i2=M

    j1,=np.where(y<ylim[0])
    if len(j1): j1=j1[-1]
    else: j1=0

    j2,=np.where(y>ylim[1])
    if len(j2): j2=j2[0]
    else: j2=L

    I0,I1,J0,J1=i1,i2,j1,j2
    J0=J0-margin
    J1=J1+1+margin
    I0=I0-margin
    I1=I1+1+margin
  else:
    # lower left and upper right:
    d=((x-xlim[0])**2+(y-ylim[0])**2)**.5
    j0,i0=np.where(d==d.min())
    j0,i0=j0[0],i0[0]

    d=((x-xlim[1])**2+(y-ylim[1])**2)**.5
    j1,i1=np.where(d==d.min())
    j1,i1=j1[0],i1[0]

    # upper left and lower right:
    d=((x-xlim[0])**2+(y-ylim[1])**2)**.5
    j2,i2=np.where(d==d.min())
    j2,i2=j2[0],i2[0]

    d=((x-xlim[1])**2+(y-ylim[0])**2)**.5
    j3,i3=np.where(d==d.min())
    j3,i3=j3[0],i3[0]

    add=1 # to ensure all xlim,ylim points are inside region
    J0=min(j0,j1,j2,j3)-add-margin
    J1=max(j0,j1,j2,j3)+1+add+margin
    I0=min(i0,i1,i2,i3)-add-margin
    I1=max(i0,i1,i2,i3)+1+add+margin

    L,M=x.shape

  J0=max(0,J0)
  I0=max(0,I0)
  J1=min(L,J1)
  I1=min(M,I1)

  if check:
    import pylab as pl
    pl.plot(xlim,ylim,'r')
    if x.ndim==1:
      [pl.plot([i,i],[y.min(),y.max()],'k') for i in x[I0:I1]]
      [pl.plot([x.min(),x.max()],[i,i],'k') for i in y[J0:J1]]

    else:
      pl.plot(x[j0,i0],y[j0,i0],'g*')
      pl.plot(x[j1,i1],y[j1,i1],'g*')
      pl.plot(x[J0:J1,I0:I1],y[J0:J1,I0:I1],'k+')

  return I0,I1,J0,J1


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


def _bilin_aux(xi,yi,vi,maski,**kargs):
  '''Used by bilin and bilin_aux.
  (processes extrapolation options)
  '''

  extrap=kargs.get('extrap',False)
  extrap_method=kargs.get('extrap_method','easy') # easy or tri

  condIn=maski==-1
  condOut=maski==0

  if extrap_method=='easy':
    extrp=easy_extrap
    args=dict(inplace=True)
  elif extrap_method=='tri':
    extrp=mask_extrap
    args=dict(inplace=True,norm_xy=True)

  if extrap=='in':
    vi=np.ma.masked_where(condIn,vi)
    extrp(xi,yi,vi,**args)
    return np.ma.masked_where(condOut,vi)
  elif extrap=='out':
    vi=np.ma.masked_where(condOut,vi)
    extrp(xi,yi,vi,**args)
    return np.ma.masked_where(condIn,vi)
  elif extrap in (0,1):
    vi=np.ma.masked_where(condIn|condOut,vi)
    if extrap==1: extrp(xi,yi,vi,**args)

  return vi


def bilin(x,y,v,xi,yi,**kargs):
  '''Bilinear interpolation.
  Domains x,y and xi,yi do not need to be fully regular!

  Options
    extrap: True, False, 'in','out' (extrap everywhere, do not extrap,
      extrap in masked points inside original domain, extrap in
      points outside original domain

    extrap_method: 'easy','tri' (extrap with a fast/simple/ugly method
      or with triangulation)
  '''

  from .alg import bilin as alg_bilin
  if np.ma.isMA(v) and v.size!=v.count(): mask=v.mask
  else: mask=np.zeros_like(v,'bool')
  vi,maski=alg_bilin(x,y,v,xi,yi,mask)
  return _bilin_aux(xi,yi,vi,maski,**kargs)


def bilin_coefs(x,y,xi,yi,mask=False):
  '''Bilinear interpolation coefficients.
  Returns coefs to be used as bilin_fast(v,coefs,...).
  Use Bilin instead.
  '''
  from .alg import bilin_coefs as alg_bc
  if mask is False: mask=np.zeros_like(x,'bool')
  coefs,inds=alg_bc(x,y,xi,yi,mask)
  return xi,yi,coefs,inds


def  bilin_fast(v,coefs,**kargs):
  '''Fast bilinear interpolation.
  Uses coefficients returned by bilin_coefs, making
  interpolations of diferent variables, but same domains,
  very fast.
  Use Bilin instead.
  '''

  from .alg import bilin_fast as alg_bf
  xi,yi,inds,coefs=coefs
  vi,maski=alg_bf(v,coefs,inds)
  return _bilin_aux(xi,yi,vi,maski,**kargs)


def easy_extrap(x,y,v,inplace=False):
  '''Simple, fast and ugly extrapolations method.
  Masked point value is calculated based on distance
  to nearest cells in each direction only!
  Use only when mask_extrap is too slow or extrapolation is not ok.
  '''

  from .alg import extrap2 as alg_extrap2
  if np.ma.isMA(v) and v.size!=v.count():
    directions=[1,1,1,1] # use all (E-W,W-E,N-S,S-N)
    if inplace: v[:]=alg_extrap2(x,y,v,v.mask,dir=directions)
    else:return alg_extrap2(x,y,v,v.mask,dir=directions)
  else: return v


class Bilin():
  '''Fast bilinear interpolation.
  Usage:
    a=Bilin(x,y,xi,yi,mask=False)
    vi=a(v)
    wi=a(w)
    zi=z(z,**kargs)

    # is faster but the same as:
    zi=bilin(x,y,z,xi,yi,**kargs)
    ...
    # and:
    C=bilin_coefs(x,y,xi,yi,mask=False)
    zi=bilin_fast(z,C,**kargs)
    ...

  kargs:
    extrap: 0,1,'in','out'
    extrap_method='easy','tri'

  See bilin from additional help
  '''

  def __init__(self,x,y,xi,yi,mask=False):
    self.coefs=bilin_coefs(x,y,xi,yi,mask)
  def __call__(self,v,**kargs):
    return bilin_fast(v,self.coefs,**kargs)


class InterpCommon:
  def solve(self,M,B):
    if 0:
      import lusolver
      lusolver.solve(M,B)
      return B
    else:
      return np.linalg.solve(M,B)

  def interp(self,x):
    try: iter(x)
    except: x=np.asarray([x])

    res=np.zeros(x.size)
    for j in range(x.size):
      i=np.where(self.x>x[j])[0]
      if i.size:
        i=i[0]-1
        if i<0: i+=1
      else:
        i=-1

      if self.C.shape[1]==3:
        a,b,c=self.C[i]
        res[j]=a*x[j]**2+b*x[j]+c
      elif self.C.shape[1]==4:
        a,b,c,d=self.C[i]
        res[j]=a*x[j]**3+b*x[j]**2+c*x[j]+d

    return res

  def __call__(self,x): return self.interp(x)


class Interp1_quad(InterpCommon):
  def __init__(self,x,y,m0='auto'):
    # calc coefs:
    nseg=x.size-1
    C=np.zeros((nseg,3),'f')
    for i in range(nseg):
      if i==0:
        if m0=='auto': C[i]=self.__coefs0(x[i:i+3],y[i:i+3])
        else: C[i]=self.__coefs(x[i:i+2],y[i:i+2],m0)
      else: C[i]=self.__coefs(x[i:i+2],y[i:i+2],2*C[i-1][0]*x[i]+C[i-1][1])

    self.x=x
    self.y=y
    self.m0=m0
    self.C=C


  def __coefs0(self,X,Y):
    if 0:
      # using 2 points and derivate in 2nd (X[1])
      M=np.asarray([
       [X[0]**2,   X[0], 1.],
       [X[1]**2,   X[1], 1.],
       [ 2*X[1],     1., 0.]])
      ym=(Y[2]-Y[0])/(X[2]-X[0])
      B=np.asarray([Y[0],Y[1],ym])
    else:
      # using 3 points:
      M=np.asarray([
       [X[0]**2,   X[0], 1.],
       [X[1]**2,   X[1], 1.],
       [X[2]**2,   X[2], 1.]])
      B=Y[:3]

    return self.solve(M,B)


  def __coefs(self,X,Y,ym):
    M=np.asarray([
       [X[0]**2,   X[0], 1.],
       [X[1]**2,   X[1], 1.],
       [ 2*X[0],     1., 0.]])
    B=np.asarray([Y[0],Y[1],ym])
    return self.solve(M,B)


class Interp1_cub(InterpCommon):
  def __init__(self,x,y,m0='auto',seg3=True):
    # calc coefs:
    nseg=x.size-1
    C=np.zeros((nseg,4),'f')
    for i in range(nseg):
      if seg3 and i%2==1:
        C[i]=C[i-1]
        continue

      if i==0:
        if m0=='auto': C[i]=self.__coefs0(x[i:i+4],y[i:i+4])
        else: C[i]=self.__coefs(x[i:i+3],y[i:i+3],m0)
      else:
        if i==nseg-1:
          if seg3:
            ym=3*C[i-1][0]*x[i]**2+2*C[i-1][1]*x[i]+C[i-1][2]
            C[i]=self.__coefs(x[i-1:i+2],y[i-1:i+2],ym,1)
          else:
            C[i]=C[i-1]

        else:
          ym=3*C[i-1][0]*x[i]**2+2*C[i-1][1]*x[i]+C[i-1][2]
          C[i]=self.__coefs(x[i:i+3],y[i:i+3],ym)

    self.x=x
    self.y=y
    self.m0=m0
    self.C=C


  def __coefs0(self,X,Y):
    M=np.asarray([
      [X[0]**3,  X[0]**2,  X[0], 1.],
      [X[1]**3,  X[1]**2,  X[1], 1.],
      [X[2]**3,  X[2]**2,  X[2], 1.],
      [X[3]**3,  X[3]**2,  X[3], 1.]])
    B=Y[:4]
    return self.solve(M,B)

  def __coefs(self,X,Y,ym,iym=0):#,ym2):
    if 0:
      # using 1st and 2nd derivative: 
      M=np.asarray([
        [  X[0]**3,  X[0]**2, X[0],  1.],
        [  X[1]**3,  X[1]**2, X[1],  1.],
        [3*X[iym]**2,   2*X[iym],   1.,  0.],
        [   6*X[0],       2.,   0.,  0.]])
      B=np.asarray([Y[0],Y[1],ym,ym2])
    else:
      # using 3 points and derivative in point iym
      M=np.asarray([
        [  X[0]**3,  X[0]**2, X[0], 1.],
        [  X[1]**3,  X[1]**2, X[1], 1.],
        [  X[2]**3,  X[2]**2, X[2], 1.],
        [3*X[0]**2,   2*X[0],   1., 0.]])
      B=np.asarray([Y[0],Y[1],Y[2],ym])

    return self.solve(M,B)


class Interp1_spl(InterpCommon):
  def __init__(self,x,y,m0=None,m1=None,im0=1,im1=1):
    '''
    x is assumed to be strictly increasing (see ext/pppack.f90)
    '''
    tau=x
    n=x.size
    c=np.zeros((4,n),'f')
    c[0]=y

    ibcbeg=0
    ibcend=0
    if not m0 is None:
      ibcbeg=im0
      c[1,0]=m0
    if not m1 is None:
      ibcbeg=im1
      c[1,-1]=m1

    import pppack


    pppack.cubspl ( tau, c, ibcbeg, ibcend )

    self.x=x
    self.y=y
    self.tau=tau
    self.c=c

  def interp(self,x,jderiv=0):
    try: iter(x)
    except: x=np.asarray([x])
    m=x.size

    #L=n-1
    #K=4 # signifying a piecewise cubic polynomial

    res=np.zeros(m)
    import pppack
    for j in range(m):
      res[j]=pppack.ppvalu(self.tau, self.c[:,:-1],x[j],jderiv)

    return res

def leggauss_ab(N,a=-1,b=1):
  '''Gauss-Legendre quadrature from a to b
  Based on numpy.polynomial.legendre.leggauss
  '''
  x,w=np.polynomial.legendre.leggauss(N)
  # Linear map from -1,1 to a,b
  x=(a*(1-x)+b*(1+x))/2.
  w=w*(b-a)/2.
  return x,w


def legendre(n,x,norm=True):
  '''Associated Legendre functions'''

  if n==0:
    if norm:
      return 1/np.sqrt(2)*(1+0*x)
    else:
      return 1+0*x

  rootn = np.sqrt(range(2*n+1))
  s = np.sqrt(1-x**2)
  P = np.zeros((n+3,x.size),x.dtype)

  e=np.geterr()['divide']
  np.seterr(divide='ignore')
  twocot = -2*x/s
  np.seterr(divide=e)

  sn=(-s)**n
  tol=np.finfo(x.dtype).tiny**0.5
  ind = np.where((s>0)& (np.abs(sn)<=tol))[0]
  if ind.size:
    v = 9.2-np.log(tol)/(n*s[ind])
    w = 1/np.log(v)
    m1 = 1+n*s[ind]*v*w*(1.0058+ w*(3.819 - w*12.173))
    m1 = np.minimum(n, np.floor(m1))

    # Column-by-column recursion
    for k in range(len(m1)):
        mm1 = int(m1[k])-1
        col = ind[k]
        P[mm1:n+1,col] = 0

        # Start recursion with proper sign
        tstart = np.finfo(x.dtype).eps
        P[mm1,col] = np.sign(np.remainder(mm1+1,2)-0.5)*tstart
        if x[col] < 0:
            P[mm1,col] = np.sign(np.remainder(n+1,2)-0.5)*tstart

        # Recur from m1 to m = 0, accumulating normalizing factor.
        sumsq = tol
        for m in range(mm1-1,-1,-1):
          P[m+1-1,col] = ((m+1)*twocot[col]*P[m+2-1,col]-
                rootn[n+m+3-1]*rootn[n-m-1]*P[m+3-1,col])/(rootn[n+m+2-1]*rootn[n-m+1-1])

          sumsq = P[m+1-1,col]**2 + sumsq

        scale = 1/np.sqrt(2*sumsq - P[0,col]**2)
        P[:mm1+1,col] = scale*P[:mm1+1,col]

  # Find the values of x,s for which there is no underflow, and for
  # which twocot is not infinite (x~=1).

  nind = np.where((x!=1)&(np.abs(sn)>=tol))[0]

  if nind.size:
    # Produce normalization constant for the m = n function
    c=(1-1/np.arange(2.,2*n+1,2)).prod()

    # Use sn = (-s).^n (written above) to write the m = n function
    P[n,nind] = np.sqrt(c)*sn[nind]
    P[n-1,nind] = P[n,nind]*twocot[nind]*n/rootn[-1]

    # Recur downwards to m = 0
    for m in  range(n-2,-1,-1):
      P[m,nind] = (P[m+1,nind]*twocot[nind]*(m+1) \
            -P[m+2,nind]*rootn[n+m+2]*rootn[n-m-1])/ \
            (rootn[n+m+2-1]*rootn[n-m+1-1])

  y = P[:n+1]

  # Polar argument   (x = +-1)
  s0 = np.where(s==0)[0]
  if s0.size:
    y[0,s0] = x[s0]**n

  if not norm:
    for m in range(n-1):
      y[m+1]=rootn[n-m:n+m+2].prod()*y[m+1]

    y[n] = rootn[1:].prod()*y[n]
  else:
    y = y*(n+0.5)**0.5
    const1 = -1
    for r in range(n+1):
        const1 = -const1
        y[r] = const1*y[r]


  return y

