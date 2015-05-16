'''
Tools to create ROMS river forcing files

mma  feb 2011, Texas A&M
'''

import numpy as np
import pylab as pl
import datetime
import os

from okean import netcdf, cookbook as cb, dateu as dts


class RiverFrc:
  '''
  Model river forcing file generation

  Parameters
  ----------
  fname : netcdf file name
  grid: model netcdf grid file
  nrivers: number or rivers to add
  nz: number of vertical levels

  **kargs:
  ncversion: netcdf version of new file (3)
  perm: file creation permission ('t')
  tunits: time variable units ('days since 1970-01-01 00:00:00')
  type: global attribute ('ROMS river forcing')
  title: global attribute ('ROMS river forcing')
  attr: dict of extra global attributes to add ({})
        ex: attr={'river_names':'Amazonas and Rio de la Plata'}
  '''

  ncversion = 3
  perm   = 't'
  tunits ='days since 1970-01-01 00:00:00'
  type   = 'ROMS river forcing'
  title  = 'ROMS river forcing'
  attr   = {} # dict with  of global attributes to add
              # ex: attr={'river_names':'Amazonas and Rio de la Plata'}

  def __init__(self,fname,grid,nrivers,nz,**kargs):
    self.grid    = grid
    self.fname   = fname
    self.nrivers = nrivers
    self.nz      = nz

    for k in kargs.keys():
      if   k=='ncversion': self.ncversion = kargs[k]
      elif k=='perm':      self.perm      = kargs[k]
      elif k=='tunits':    self.tunits    = kargs[k]
      elif k=='type':      self.type      = kargs[k]
      elif k=='title':     self.title     = kargs[k]
      elif k=='attr':      self.attr      = kargs[k]


  def create(self):
    '''
    Creates model netcdf river forcing file
    '''
    nc=netcdf.Pync(self.fname,self.perm,version=self.ncversion)

    nx=netcdf.fdim(self.grid,'xi_rho')
    ny=netcdf.fdim(self.grid,'eta_rho')

    # Dimensions:
    nc.add_dim('s_rho',self.nz)
    nc.add_dim('river',self.nrivers)
    nc.add_dim('river_time',0)

    # Variables:
    v=nc.add_var('river',np.dtype('d'),('river',))
    v.add_att('long_name','river runoff identification number')

    v=nc.add_var('river_Xposition',np.dtype('d'),('river',))
    v.add_att('long_name','river XI-position at RHO-points')
    v.add_att('valid_min',1)
    v.add_att('valid_max',nx-1)

    v=nc.add_var('river_Eposition',np.dtype('d'),('river',))
    v.add_att('long_name','river ETA-position at RHO-points')
    v.add_att('valid_min',1)
    v.add_att('valid_max',ny-1)

    v=nc.add_var('river_direction',np.dtype('d'),('river',))
    v.add_att('long_name','river runoff direction')

    v=nc.add_var('river_Vshape',np.dtype('d'),('s_rho','river'))
    v.add_att('long_name','river runoff mass transport vertical profile')

    v=nc.add_var('river_time',np.dtype('d'),('river_time',))
    v.add_att('long_name','river runoff time')
    v.add_att('units',self.tunits)
    v.add_att('add_offset',0)

    v=nc.add_var('river_transport',np.dtype('d'),('river_time','river'))
    v.add_att('long_name','river runoff vertically integrated mass transport')
    v.add_att('units','metre3 second-1')
    v.add_att('time','river_time')

    v=nc.add_var('river_temp',np.dtype('d'),('river_time','s_rho','river'))
    v.add_att('long_name','river runoff potential temperature')
    v.add_att('units','Celsius')
    v.add_att('time','river_time')

    v=nc.add_var('river_salt',np.dtype('d'),('river_time','s_rho','river'))
    v.add_att('long_name','river runoff salinity')
    v.add_att('units','Celsius')
    v.add_att('time','river_time')

    # Global Attributes:
    nc.add_att('type',self.type)
    nc.add_att('title',self.title)
    nc.add_att('grd_file',os.path.realpath(self.grid))
    nc.add_att('date',dts.currday().isoformat(' '))
    nc.add_att('author',cb.username()[1]+', '+cb.username()[0]+'@'+cb.machinename())

    # extra attrs:
    for i in self.attr.keys(): nc.add_att(i,self.attr[i])

    nc.close()


  def fill(self,data,quiet=1):
    '''
    Fills model netcdf river forcing file

    The input data shoud have as keys the river names
    (with river runoff, temp, salt, vshape, xi,eta and direction)
    and time (datetime)

    Example:
      data['amazonas']=q,temp,salt,vshape,Xi,Eta,Dir
      data['time']=np.arange(20)
    '''

    # convert time to tunits:
    time=data['time']
    for i in range(len(time)): time[i]=netcdf.date2num(time[i],self.tunits)

    nc=netcdf.Pync(self.fname,'a')

    if not quiet: print ' -- filling time...'
    for i in range(len(time)): nc.vars['river_time'][i]=time[i]

    cont=-1
    for k in data.keys():
      if k=='time': continue
      cont+=1

      q,temp,salt,vshape,Xi,Eta,Dir=data[k]
      if not quiet: print ' -- filling river',cont
      nc.vars['river'][cont]        = cont
      nc.vars['river_Xposition'][cont]   = Xi
      nc.vars['river_Eposition'][cont]   = Eta
      nc.vars['river_direction'][cont]   = Dir
      nc.vars['river_Vshape'][:,cont]    = vshape
      nc.vars['river_transport'][:,cont] = q
      nc.vars['river_temp'][:,:,cont]    = temp
      nc.vars['river_salt'][:,:,cont]    = salt

    nc.close()


def gen_vshape(nlevels,type='uniform',rbot=0,rtop=1):
  '''
  Returs river vshape in linear, exponential or uniform scale, for
  nlevels vertical levels, starting by rbot to rtop. Final shape
  is normalised to have sum 1.
  '''
  if   type.startswith('unif'): vshape=1./nlevels*np.ones(nlevels)
  elif type.startswith('lin'):  vshape=np.linspace(rbot,rtop,nlevels)/np.linspace(rbot,rtop,nlevels).sum()
  elif type.startswith('exp'):  vshape=np.logspace(rbot,rtop,nlevels)/np.logspace(rbot,rtop,nlevels).sum()

  return vshape

