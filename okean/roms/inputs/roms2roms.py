import os
import numpy as np
import pylab as pl

from okean import calc, dateu as dt, netcdf
from okean.roms import roms
import prognostic
import gennc

sett=prognostic.DataAccess(xdim='xi_rho',ydim='eta_rho',zdim='s_rho',tdim='ocean_time',
                           x_name='lon_rho',y_name='lat_rho',z_name=False,
                           time_name='ocean_time',ssh_name='zeta')


def roms2clmbry(f0,clm,grd,sparams,**kargs):
  '''
  kargs:
  bry, bry file to create (False by default, no bry is created)
  times, all by default, but can be a range of times (inds or datetimes)
  tunits, 'days since 1970-01-01'

  ctitle, title for clm file (default used)
  btitle, title for bry file (defaul used)
  obc, bry open boundaries (defaul used)
  quiet, False

  other kargs of roms2roms:
    grd0
    sparams0
    tunits0
  '''

  bry    = False
  times  = 'all'
  tunits = 'days since 1970-01-01'
  ctitle = False
  btitle = False
  obc    = False
  quiet  = False

  if 'bry'    in kargs.keys(): bry    = kargs['bry']
  if 'times'  in kargs.keys(): times  = kargs['times']
  if 'tunits' in kargs.keys(): tunits = kargs['tunits']
  if 'ctitle' in kargs.keys(): ctitle = kargs['ctitle']
  if 'btitle' in kargs.keys(): btitle = kargs['btitle']
  if 'obc'    in kargs.keys(): obc    = kargs['obc']
  if 'quiet'  in kargs.keys(): quiet  = kargs['quiet']

  create=True

  Cargs={}
  if ctitle: Cargs['title']  = ctitle
  Cargs['tunits'] = tunits
  Cargs['create'] = create

  Bargs={}
  if btitle: Bargs['title'] = btitle
  if obc:    Bargs['obc']   = obc
  Bargs['tunits'] = tunits
  Bargs['create'] = create


  if times=='all':
    times=range(netcdf.fdim(f0,'ocean_time'))

  for tind in times:
    data,datab=roms2roms(f0,grd,sparams,tind,**kargs)

    prognostic.make_clm(clm,data,grd,sparams,quiet=quiet,**Cargs)
    if bry:
      prognostic.make_bry(bry,datab,grd,sparams,quiet=quiet,**Bargs)

    if create:
      Cargs['create']=False
      Bargs['create']=False



def roms2roms(f0,grd,sparams,tind=0,**kargs):
  '''
  tind: ind or datetime

  **kargs:
  grd0
  sparams0
  tunits0
  quiet

  '''

  grd0     = False
  sparams0 = False
  quiet    = False
  tunits0  = False
  Z='auto'
  ij='i'

  for k in kargs.keys():
    if   k=='grd0':     grd0      = kargs[k]
    elif k=='sparams0': sparams0 = kargs['sparams0']
    elif k=='quiet':    quiet    = kargs['quiet']
    elif k=='tunits0':  tunits0  = kargs['tunits0']
    elif k=='Z':  Z  = kargs['Z']
    elif k=='ij':  ij  = kargs['ij']


  r0=roms.His(f0,grd=grd0)
  g0=r0.grid
  g1=roms.Grid(grd)
  if sparams0 is False: sparams0=r0.s_params

  # load data:
  F0={}
  for i in ('temp','salt','u','v','ssh','time','misc'): F0[i]=f0
  F0['xy']=g0.name

  if tind=='all':
    times=range(r0.TIME)
  else: times=[tind]

  xlim=g1.lon.min(),g1.lon.max()
  ylim=g1.lat.min(),g1.lat.max()
  inds={}

  outdata={}
  outdatab={}
  for tind in times:
    inds[sett.tdim]=tind

    data,msg=prognostic.load_data(F0,quiet=quiet,settings=sett,inds=inds,t_units=tunits0)

    # z3d:
    data['z3d']=r0.s_levels(tind)

    # u,v at rho:
    print '  u,v at rho ...'
    u=np.ma.zeros((data['NZ'],data['NY'],data['NX']),data['u'].dtype)
    v=np.ma.zeros((data['NZ'],data['NY'],data['NX']),data['v'].dtype)

    u[:,:,1:-1]=(data['u'][:,:,:-1]+data['u'][:,:,1:])/2.
    u[:,:,0]=data['u'][:,:,0]
    u[:,:,-1]=data['u'][:,:,-1]

    v[:,1:-1,:]=(data['v'][:,:-1,:]+data['v'][:,1:,:])/2.
    v[:,0,:]=data['v'][:,0,:]
    v[:,-1,:]=data['v'][:,-1,:]

    data['u']=u
    data['v']=v

    # simplify data:
    print '  simplify data...'
    i0,i1,j0,j1=calc.ij_limits(g0.lon,g0.lat,xlim,ylim,margin=3)
    for v in 'z3d','temp','salt','u','v': data[v]=data[v][:,j0:j1,i0:i1]
    for v in 'ssh','lon','lat': data[v]=data[v][j0:j1,i0:i1]
    data['NZ'],data['NY'],data['NX']=data['temp'].shape

    # interp depths:
    if Z is 'auto':
      h=-g0.h
      Z=np.concatenate((np.arange(data['ssh'].max(),-2,-.05),np.arange(-2,-5,.2),np.arange(-5,-20,-1),np.arange(-20,-100,-2),
      np.arange(-100,-500,-5),np.arange(-500,-1000,-20),np.arange(-1000,h.min()-100,-100)))
      Z=Z[::3]
      if Z[-1]>h.min(): Z[-1]=h.min()

    data['depth']=Z

    for v in 'temp','salt','u','v','ssh':
      print '  %-6s %6.3f %6.3f'%(v,data[v].min(), data[v].max())

    # to z levels:
    Data=prognostic.data2z(data,quiet=quiet,ij=ij)
    for v in 'temp','salt','u','v','ssh':
      print '  %-6s %6.3f %6.3f'%(v,Data[v].min(), Data[v].max())


    # clm:
    data,HA=prognostic.data2roms(Data,grd,sparams,quiet=quiet,horizAux=True,ij=ij)
    for v in 'temp','salt','u','v','zeta':
      print '  %-6s %6.3f %6.3f'%(v,data[v].min(), data[v].max())

    # bry:
    datab=prognostic.data2romsbry(Data,grd,sparams,quiet=quiet,horizAux=HA)

    outdata[tind]  = data
    outdatab[tind] = datab

  if len(times)==1:  return outdata[tind],outdatab[tind]
  else: return outdata,outdatab
