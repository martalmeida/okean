# -*- coding: utf-8 -*-

import os
import datetime
import numpy as np
from .. import netcdf, calc, dateu, cookbook as cb


def depthof(v,z,val):
  '''
  Depths where variable, increasing/dec with depth, has some value
  val can be 2d

  Output is masked where surface is higher/lower than value or where all
  water column is lower/higher than value

  How to get thickness?
  Ex, if v inc with depth:
  A=output.data
  A[A==9999]=-h
  A=A+zeta
  '''

  try: val.shape==v.shape[1:]
  except: val=np.ones(v.shape[1:],v.dtype)*val

  from . import rtools
  mask=(~v[0].mask).astype('i')
  zz=rtools.depthof(v,z,mask,val)
  # 999 --> surface higher than val
  # 9999 --> all whater column lower than val

  return np.ma.masked_where((zz==999)|(zz==9999),zz)


def slicez(v,maskv,h,zeta,sparams,level,surface_masked=True,spline=True):
  tts,ttb,hc,Nr,vt,vs=sparams

  from . import rtools

  try: level.shape==v.shape[1:]
  except: level=np.ones(v.shape[1:],v.dtype)*level

  N,Ny,Nx=v.shape
  res=rtools.roms_slicez(v,h,zeta, tts,ttb,hc,Nr,vt,vs,
                         level,surface_masked,spline,
                         N,Ny,Nx)

  mask=np.where(res==-99.,0,1)*maskv==0
  res=np.ma.masked_where(mask,res)
  return res


def s_levels(h,zeta,sparams,rw=False,nomask=True):
  '''
  returns z_r and/or z_w
  '''
  tts,ttb,hc,n,vt,vs=sparams

  if nomask:
    if np.ma.isMA(h): h=np.where(h.mask,999,h)
    if np.ma.isMA(zeta): zeta=np.where(zeta.mask,0,zeta)
    
  # to deal with WET_DRY (See Nonlinear/set_depth.F)
  # h cannot be zero, at least with vt==1
  h=np.where(h==0.,1e-14,h)

  from . import rtools
  zr,zw=rtools.s_levels(h,zeta,tts,ttb,hc,n,vt,vs)

  zr=np.squeeze(zr)
  zw=np.squeeze(zw)

  if not nomask: # add zeta|h mask to z levels
    if np.ma.is_masked(zeta): mzeta=zeta.mask
    else: mzeta=np.zeros_like(zeta,'bool')

    if np.ma.is_masked(h): mh=h.mask
    else: mh=np.zeros_like(h,'bool')

    mask=mh|mzeta
    if np.any(mask):
      maskr=np.tile(mask,(n,)+(1,)*mask.ndim)
      maskw=np.tile(mask,(n+1,)+(1,)*mask.ndim)
    
      zr=np.ma.masked_where(maskr,zr)
      zw=np.ma.masked_where(maskw,zw)

  if rw==False: return zr,zw
  elif rw=='w': return zw
  else: return zr



  # add zeta mask to z levels:
#  if np.ma.isMA(zeta):
#  if np.ma.is_masked(zeta):
#    if 1 in zr.shape:
#      maskr=np.tile(np.squeeze(zeta).mask,(zr.shape[0],1))
#      maskw=np.tile(np.squeeze(zeta).mask,(zw.shape[0],1))
#      zr=np.squeeze(zr)
#      zw=np.squeeze(zw)
#    else:
#      maskr=np.tile(zeta.mask,(zr.shape[0],1,1))
#      maskw=np.tile(zeta.mask,(zw.shape[0],1,1))
#
#
#    if rw==False:
#      zr=np.ma.masked_where(maskr,zr)
#      zw=np.ma.masked_where(maskw,zw)
#      return zr,zw
#    elif rw=='w':
#      return np.ma.masked_where(maskw,zw)
#    else: # r
#      return np.ma.masked_where(maskr,zr)
#
#  else: # zeta has no mask, like in clm or ini files
#    if rw==False: return zr,zw
#    elif rw=='w': return zw
#    else: return zr


#def ll_rho2uvp(lonr,latr,Type=None):
#  '''
#  Converts lon lat at rho to u, v or p points
#  '''
#  type_=Type
#  Type=Type[0]
#  if Type not in ('u','v','p'):
#    print(':: ll_rho2uvp, unknown Type ',type_)
#    return
#
#  if Type=='u':
#    lon=(lonr[:,:-1]+lonr[:,1:])/2.
#    lat=(latr[:,:-1]+latr[:,1:])/2.
#  elif Type=='v':
#    lon=(lonr[:-1,:]+lonr[1:,:])/2.
#    lat=(latr[:-1,:]+latr[1:,:])/2.
#  elif Type=='p':
#    lon=(lonr[:,:-1]+lonr[:,1:])/2.
#    lat=(latr[:,:-1]+latr[:,1:])/2.
#
#    lon=(lon[:-1,:]+lon[1:,:])/2.
#    lat=(lat[:-1,:]+lat[1:,:])/2.
#
#  return lon,lat


def rho2uvp(*pargs):
  '''
  Convert fields at rho to u,v or p points
  Last parg is the type, u,v,p
  ex: lonu,latu=rho2uvp(lonr,latr,'u')
  will be the same as ll_rho2uvp
  '''
  Type=pargs[-1][0]
  pargs=list(pargs[:-1])
  for i in range(len(pargs)):
    v=pargs[i]
    if Type=='u':
      v=(v[:,:-1]+v[:,1:])/2.
    elif Type=='v':
      v=(v[:-1,:]+v[1:,:])/2.
    elif Type=='p':
      v=(v[:,:-1]+v[:,1:])/2.
      v=(v[:-1,:]+v[1:,:])/2.

    pargs[i]=v

  if i==0: pargs=pargs[0]
  return pargs


def rho2uvp3d(*pargs):
  '''
  Convert 3d fielfs at rho to u,v or p points
  (last 2 dims are used)
  ex: u = rho2uvp3d(ur,'u')
  '''
  Type=pargs[-1][0]
  pargs=list(pargs[:-1])
  for i in range(len(pargs)):
    v=pargs[i]
    if Type=='u':
      v=(v[:,:,:-1]+v[:,:,1:])/2.
    elif Type=='v':
      v=(v[:,:-1,:]+v[:,1:,:])/2.
    elif Type=='p':
      v=(v[:,:,:-1]+v[:,:,1:])/2.
      v=(v[:,:-1,:]+v[:,1:,:])/2.

    pargs[i]=v

  if i==0: pargs=pargs[0]
  return pargs


def psi2rho(vp): return psi2uvr(vp,'r')


def psi2uvr(vp,pt):
  '''
  Transfert a field at psi points to u, v or rho points
  '''

  M,L=vp.shape
  if pt[0]=='r': M_,L_=M+1,L+1
  elif  pt=='u': M_,L_=M+1,L
  elif  pt=='v': M_,L_=M,L+1

  if np.ma.isMA(vp):
    vr=np.ma.zeros((M_,L_),vp.dtype)
  else:
    vr=np.zeros((M_,L_),vp.dtype)

  if pt[0]=='r':
    vr[1:-1,1:-1]=0.25*(vp[:-1,:-1]+vp[:-1,1:]+vp[1:,:-1]+vp[1:,1:])
    vr[0,:]  = vr[1,:]
    vr[-1,:] = vr[-1,:]
    vr[:,0]  = vr[:,1]
    vr[:,-1] = vr[:,-2]
  elif pt=='u':
    vr[1:-1,:]=0.5*(vp[1:,:]+vp[:-1,:])
    vr[0,:] = vr[1,:]
    vr[-1,:] = vr[-2,:]
  elif pt=='v':
    vr[:,1:-1]=0.5*(vp[:,1:]+vp[:,:-1])
    vr[:,0] = vr[:,1]
    vr[:,-1] = vr[:,-2]

  return vr


def uvp_mask(rfield):
  '''Compute the mask at u,v and psi points
     Returns u,v,pmask
     Based on Pierrick Penven uvp_mask
     (http://www.brest.ird.fr/Roms_tools)

     mma 2007, La Coruña
  '''
  Mp,Lp=rfield.shape
  M=Mp-1
  L=Lp-1

  ufield=rfield[:,:-1]*rfield[:,1:]
  vfield=rfield[:-1,:]*rfield[1:,:]
  pfield=ufield[:-1,:]*ufield[1:,:]

  return ufield, vfield, pfield


def s_params(nc,show=0):
  """
  Get s-coordinates parameters from ROMS file
  Gets from ROMS output files the data needed to calculate vertical
  s-coordinate levels.
  This function looks for the required values among variables, file
  dimensions and file attributes.

  mma 28-7-2007;
  vs and vt added aug 2013 (Guayaquil)
  """

  if cb.isstr(nc) or isinstance(nc,list) or isinstance(nc,tuple): nc=netcdf.ncopen(nc)

  theta_s=theta_b=hc=N=False
  theta_s_src=theta_b_src=hc_src=N_src=''

  vt=1
  vs=1
  if 'Vtransform'  in nc.varnames: vt=nc.vars['Vtransform'][:]
  if 'Vstretching' in nc.varnames: vs=nc.vars['Vstretching'][:]


  v='theta_s'
  if v in nc.atts:
    theta_s=nc.atts[v].value
    theta_s_source='file attribute'
  elif v in nc.varnames:
    theta_s=nc.vars[v][:]
    theta_s_source='variable'

  v='theta_b'
  if v in nc.atts:
    theta_b=nc.atts[v].value
    theta_b_source='file attribute'
  elif v in nc.varnames:
    theta_b=nc.vars[v][:]
    theta_b_source='variable'

  v='hc'
  if v in  nc.atts:
    hc=nc.atts[v].value
    hc_source='file attribute'
  elif v in nc.varnames:
    hc=nc.vars[v][:]
    hc_source='variable'
  else:
    # if vtransform=1 --> hc=min(Tcline,hmin)
    # else: hc=Tcline

    if vt==1:
      hmin=[]
      Tcline=[]

      v='Tcline'
      if v in  nc.atts:
        Tcline=nc.atts[v].value
      elif v in nc.varnames:
        Tcline=nc.vars[v][:]

      if 'h' in nc.varnames:
        hmin=nc.vars['h'][:].min()


      if hmin and Tcline:
        #hc=np.min([hmin,Tcline])
        hc=np.min([np.max([hmin,0]),Tcline]) # to deal with WET_DRY (see Utility/set_scoord.F)
        hc_source = 'min of hmin and Tcline';
      else:
        hc='not found'

    else: # vt==2
      hc=Tcline
      hc_source = 'Tcline'


  if 's_rho' in nc.dims:
    N=nc.dims['s_rho']
    N_source = 'file dimension s_rho';
  elif 'sc_r' in nc.atts:
    N=nc.atts['sc_r']['value']
    N_source = 'file attribute sc_r';
  elif 'sc_r' in nc.varnames:
    N=len(nc.vars['sc_r'][:])
    N_source = 'length of variable sc_r';

  if show: # show sources:
    print('%-10s : %4.2f  %-10s' % ('theta_s',theta_s,theta_s_source))
    print('%-10s : %4.2f  %-10s' % ('theta_b',theta_b,theta_b_source))
    print('%-10s : %4.2f  %-10s' % ('hc',hc,hc_source))
    print('%-10s : %4d  %-10s' % ('N',N,N_source))
    print('%-10s : %4d  %-10s' % ('vt',vt,'variable'))
    print('%-10s : %4d  %-10s' % ('vs',vs,'variable'))

  return theta_s, theta_b, hc, N,vt,vs


def roms_read_out(f):
  '''Parse ROMS output text file
     Returns Time, Ek,Ep,Etotal,Volume

     Ex: time,kin,pot,tot,vol=roms_read_out('roms.out')

     Martinho MA, 2012
  '''

  tag='   STEP   Day'
  nskip=2

  badFormatStr='******'

  fid=open(f)
  lines=fid.readlines()
  n=-1
  for l in lines:
    if l.find('time_ref')>=0: print(str(int(float(l.split()[0]))))
    if l.find('time_ref')>=0: time_ref=dateu.parse_date(str(int(float(l.split()[0]))))
    if l.find(tag)==0: break
    n+=1

  n=n+nskip+1

  time=[]
  kin=[]
  pot=[]
  tot=[]
  vol=[]

  c=-1
  for l in lines[n:]:
    c+=1
    if c==0: L=len(l)
    if l.find('_')!=-1 or l.find('=')!=-1: continue # DEF_HIS, etc
    if l.find(badFormatStr): l.replace(badFormatStr,'0')
    if len(l)==L:
      tmp=l.split()

      # time:
      days=int(tmp[1])
      hh,mm,ss=[int(j) for j in tmp[2].split(':')]
      time+=[time_ref+datetime.timedelta(days=days,hours=hh,minutes=mm,seconds=ss)]

      kin+=[float(tmp[3])]
      pot+=[float(tmp[4])]
      tot+=[float(tmp[5])]
      vol+=[float(tmp[6])]

  return np.asarray(time),np.asarray(kin),np.asarray(pot),np.asarray(tot),np.asarray(vol)



def roms_read_out_agrif(f):
  '''
  ROMS_READ_OUT   Parse ROMS-AGRIF output text file
  Old. Use roms_read_out instead!!

     Syntax:
       TIME,KIN,POT,TOT,VOL = ROMS_READ_OUT(FILE)

     Inputs:
        FILE   ROMS output file (roms.out)

     Outputs:
        TIME
        KIN   Kinetic energy
        POT   Potential energy
        TOT   Total energy
        VOL   Net volume

     Example:
        time,kin,pot,tot,vol = roms_read_out_agrif('roms.out')

     Martinho MA (mma@odyle.net) and Janini P (janini@usp.br)
     Dep. Earth Physics, UFBA, Salvador, Bahia, Brasil
     02-10-2008
  '''

  tag=' MAIN: started time-steping.'
  nskip=2
  badFormatStr='******'

  fid=open(f)
  lines=fid.readlines()
  n=-1
  for l in lines:
    if l.find(tag)==0: break
    n+=1

  n=n+nskip+1

  data=range(len(lines))

  i=0
  for l in lines[n:]:
    if i==0: L=len(l)
    if l.find('_')!=-1: continue
    if l.find(badFormatStr): l.replace(badFormatStr,'0')
    if len(l)==L:
      d=[float(j) for j in l.split()]
      data[i]=d
      i+=1

  data = array(data[:i])
  time = data[:,1]
  kin  = data[:,2]
  pot  = data[:,3]
  tot  = data[:,4]
  vol  = data[:,5]

  return time,kin,pot,tot,vol


def is_roms_out_ok(ro):
  if not os.path.isfile(ro): return False
  s=open(ro).read()

  # ROMS AGRIF
  res=s[-12:-2]=='MAIN: DONE'

  if not res:
    # ROMS
    res=s.find('ROMS/TOMS: DONE')>0 and s.find('Blowing-up:')==-1

  return res


def barotropic(var,zeta,h,sparams):
  '''uvbar only works for u and v 3d!!
  this func works for all vars but correct
  zeta and h must be provided, ie, for var=u
  must be provided zetau and hu
  '''

  zr,zw=s_levels(h,zeta,sparams)
  dz=np.diff(np.squeeze(zw),axis=0)
  return np.sum(var*dz,0)/np.sum(dz,0)


def uvbar(u,v,zeta,h,sparams):
  '''u and v are 3d !'''

  try:
    zeta.shape==h.shape
  except:
    zeta=tile(zeta,h.shape).astype(h.dtype)
  pass

  L,M=h.shape

  zr,zw=s_levels(h,zeta,sparams)

  zwu=(zw[:,:,:-1]+zw[:,:,1:])/2
  zwv=(zw[:,:-1,:]+zw[:,1:,:])/2

  dzu=np.diff(zwu,axis=0)
  dzv=np.diff(zwv,axis=0);

  ubar=np.sum(u*dzu,0)/np.sum(dzu,0)
  vbar=np.sum(v*dzv,0)/np.sum(dzv,0)
  return ubar, vbar


def grid_vicinity(grid,x,y,margin=5,rect=False,retinds=False):
  '''
  Returns True for x,y points inside roms grid boundary plus margin.
  Margin is the number of cells to add around the grid.

  if rect is True returns True for all points in the smallest 2d xy grid
  (usually a rectangle) around the grid.
  In this case,margin is the rectangle margin, ie, in units of xy, not
  in  units of grid

  if retinds, in the case of rect, the rectangle j1,j2 and i1,i2 are
  also returned (cond,inds=grid_vicinity(....); i1,i2,j1,j2=inds)

  mma, TAMU 2011
  '''
  try:
    xg=grid.lon
    yg=grid.lat
    isGrid=True
  except:
    xg=netcdf.use(grid,'lon_rho')
    yg=netcdf.use(grid,'lat_rho')
    isGrid=False

  xlim=xg.min(),xg.max()
  ylim=yg.min(),yg.max()

  if x.ndim==1 and y.ndim==1: x,y=np.meshgrid(x,y)

  if rect:
    out=np.zeros(x.shape,'bool')
    i1,i2,j1,j2=calc.ij_limits(x,y,xlim,ylim,margin)
    out[j1:j2,i1:i2]=1
  else:
    if not isGrid:
      from roms import Grid
      g=Grid(grid)
    else: g=grid

    xb,yb=g.border(margin=-margin)
    out=calc.inpolygon(x,y,xb,yb)

  if rect and retinds: return out,(i1,i2,j1,j2)
  else: return out
