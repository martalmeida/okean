from os.path import realpath
import numpy as np
from time import ctime

from okean import netcdf
from okean import cookbook as cb
from okean import dateu as dts

class GenCommon:
  '''
  Base class for model input netcdf files generation

  Parameters
  ----------
  filename : netcdf file name
  grid : model netcdf grid file
  sparams : s-coordinates parameters, theta_s,theta_b, hc and NLevels,
                                      Vtransform, Vstretching

  **kargs:
  type : global attribute ('ROMS file')
  title : global attribute ('ROMS file')
  tunits : time units ('days since 1970-01-01')
  cycle : time variables cycle_length (False, attribute not created)
  ncversion : netcdf file format/version (3)
  attr: dict of extra global attributes to add ({})

  '''

  def __init__(self,filename,grid,sparams=False,**kargs):
    type   = 'ROMS file'
    title  = 'ROMS file'
    tunits = 'days since 1970-01-01'
    cycle=False
    ncversion= 3
    attr={}

    if 'type'   in kargs.keys(): type   = kargs['type']
    if 'title'  in kargs.keys(): title  = kargs['title']
    if 'tunits' in kargs.keys(): tunits = kargs['tunits']
    if 'cycle'    in kargs.keys(): cycle=kargs['cycle']
    if 'ncversion' in kargs.keys(): ncversion = kargs['ncversion']
    if 'attr'   in kargs.keys(): attr = kargs['attr']

    self.type   = type
    self.title  = title
    self.tunits = tunits
    self.cycle = cycle
    self.ncversion = ncversion
    self.attr   = attr

    self.filename  = filename
    self.grid      = grid
    self.sparams   = sparams


  def add_horiz_dims(self,nc):
    '''
    Adds horizontal dimensions to netcdf file

    '''
    grd_dims=netcdf.fdim(self.grid)
    gdims='xi_rho','xi_u','xi_v','eta_rho','eta_u','eta_v'
    for name in gdims: nc.add_dim(name,grd_dims[name])


  def add_vert_dims(self,nc,addw=False):
    '''
    Adds vertical dimension(s) to netcdf file

    '''
    nc.add_dim('s_rho',self.sparams[3])
    if addw:
      nc.add_dim('s_w',self.sparams[3]+1)


  def set_time(self,time=0,tunits=False):
    '''
    Set (force) time and date info to user defined.
    '''

    nc=netcdf.Pync(self.filename,'w')
    for v in nc.varnames:
      if v.endswith('time'):
        nc.vars[v][:]=time
        if tunits: nc.vars[v].atts['units'].update(tunits)

    nc.close()


class GenIni(GenCommon):
  '''
  Model intial file generation

  Parameters
  ----------
  filename : netcdf file name
  grid : model netcdf grid file
  sparams : s-coordinates parameters, theta_s,theta_b, hc and NLevels,
                                      Vtransform, Vstretching

  **kargs:
  type : global attribute ('ROMS Initial file')
  title : global attribute ('ROMS Initial file')
  date : date as datetime (or something convertible into datetime through
         date_tools.parse_date) or number (0)
  tunits : time units ('days since 1970-01-01')

  see GenCommon for additional kargs

  Example
  -------
  with time reference:
  ob=GenIni('ini.nc','grid.nc',(5,0,5,20),title='My initial file',
            date='2011-01-01',tunits='seconds since 2000-01-01')

  without time reference:
  ob=GenIni('ini.nc','grid.nc',(5,0,5,20),title='My initial file',
            date=30,tunits='days')

  ob.create() # created the netcdf file
  ob.fill(data) # fills the netcdf file with data dict (may be provided
                # by prognostic.data2roms

  '''

  def __init__(self,filename,grid,sparams,**kargs):
    if not 'type'  in kargs.keys(): kargs['type']  = 'ROMS Initial file'
    if not 'title' in kargs.keys(): kargs['title'] = 'ROMS Initial file'

    GenCommon.__init__(self,filename,grid,sparams,**kargs)

    date = 0
    if 'date'   in kargs.keys(): date   = kargs['date']
    self.date   = date

    # about time:
    try:
      self.date=dts.parse_date(self.date)
      self.time=netcdf.date2num(self.date,self.tunits)
    except:
      self.time=self.date # date as number!


  def create(self):
    '''
    Creates model netcdf initial file
    '''

    nc=netcdf.Pync(self.filename,'t',version=self.ncversion)

    # Dimensions:
    # horizontal and vertical dims:
    self.add_horiz_dims(nc)
    self.add_vert_dims(nc)
    # time dim:
    nc.add_dim('ocean_time',1)

    # Variables:
    v=nc.add_var('ocean_time',np.dtype('d'),('ocean_time',))
    v.add_att('long_name','initial time')
    v.add_att('units',self.tunits)
    v.add_att('field','time, scalar, series')

    v=nc.add_var('temp',np.dtype('d'),('ocean_time','s_rho','eta_rho','xi_rho'))
    v.add_att('long_name','potential temperature')
    v.add_att('units','Celsius')
    v.add_att('time','ocean_time')
    v.add_att('coordinates','lon_rho lat_rho s_rho ocean_time')
    v.add_att('field','temperature, scalar, series')

    v=nc.add_var('salt',np.dtype('d'),('ocean_time','s_rho','eta_rho','xi_rho'))
    v.add_att('long_name','salinity')
    v.add_att('units','PSU')
    v.add_att('time','ocean_time')
    v.add_att('coordinates','lon_rho lat_rho s_rho ocean_time')
    v.add_att('field','salinity, scalar, series')

    v=nc.add_var('u',np.dtype('d'),('ocean_time','s_rho','eta_u','xi_u'))
    v.add_att('long_name','u-momentum component')
    v.add_att('units','metre second-1')
    v.add_att('time','ocean_time')
    v.add_att('coordinates','lon_u lat_u s_rho ocean_time')
    v.add_att('field','u-velocity, scalar, series')

    v=nc.add_var('v',np.dtype('d'),('ocean_time','s_rho','eta_v','xi_v'))
    v.add_att('long_name','v-momentum component')
    v.add_att('units','metre second-1')
    v.add_att('time','ocean_time')
    v.add_att('coordinates','lon_v lat_v s_rho ocean_time')
    v.add_att('field','v-velocity, scalar, series')

    v=nc.add_var('ubar',np.dtype('d'),('ocean_time','eta_u','xi_u'))
    v.add_att('long_name','vertically integrated  u-momentum component')
    v.add_att('units','metre second-1')
    v.add_att('time','ocean_time')
    v.add_att('coordinates','lon_u lat_u ocean_time')
    v.add_att('field','ubar-velocity, scalar, series')

    v=nc.add_var('vbar',np.dtype('d'),('ocean_time','eta_v','xi_v'))
    v.add_att('long_name','vertically integrated  v-momentum component')
    v.add_att('units','metre second-1')
    v.add_att('time','ocean_time')
    v.add_att('coordinates','lon_v lat_v ocean_time')
    v.add_att('field','vbar-velocity, scalar, series')

    v=nc.add_var('zeta',np.dtype('d'),('ocean_time','eta_rho','xi_rho'))
    v.add_att('long_name','free-surface')
    v.add_att('units','metre')
    v.add_att('time','ocean_time')
    v.add_att('coordinates','lon_rho lat_rho ocean_time')
    v.add_att('field','free-surface, scalar, series')


    # sparams vars:
    v=nc.add_var('theta_s',np.dtype('d'),())
    v.add_att('long_name','S-coordinate surface control parameter')

    v=nc.add_var('theta_b',np.dtype('d'),())
    v.add_att('long_name','S-coordinate bottom control parameter')

    v=nc.add_var('Tcline',np.dtype('d'),())
    v.add_att('long_name','S-coordinate surface/bottom layer width')
    v.add_att('units','metre')

    v=nc.add_var('hc',np.dtype('d'),())
    v.add_att('long_name','S-coordinate parameter, critical depth')
    v.add_att('units','metre')

    v=nc.add_var('Vtransform',np.dtype('i'),())
    v.add_att('long_name','vertical terrain-following transformation equation')

    v=nc.add_var('Vstretching',np.dtype('i'),())
    v.add_att('long_name','vertical terrain-following stretching function')


    # Global Attributes:
    nc.add_att('type',self.type)
    nc.add_att('title',self.title)
    nc.add_att('grd_file',realpath(self.grid))
    nc.add_att('history','ROMS ini file, '+ctime())
    nc.add_att('author',cb.username()[1]+', '+cb.machinename())

    # extra attrs:
    for i in self.attr.keys(): nc.add_att(i,self.attr[i])

    # fill time:
    nc.vars['ocean_time'][:]=self.time

    # fill sparams:
    nc.vars['theta_s'][:]=self.sparams[0]
    nc.vars['theta_b'][:]=self.sparams[1]
    nc.vars['Tcline'][:]=self.sparams[2]
    nc.vars['hc'][:]=self.sparams[2]
    nc.vars['Vtransform'][:]=self.sparams[4]
    nc.vars['Vstretching'][:]=self.sparams[5]

    nc.close()

  def fill(self,data,quiet=1):
    '''
    Fills model netcdf initial file (data can be provided by
    prognostic.data2roms

    '''

    nc=netcdf.Pync(self.filename,'w')
    names='temp','salt','u','v','ubar','vbar','zeta'
    if not quiet: print 'filling ini file %s' % self.filename
    for i in names:
      if not quiet: print '  %s' % i
      nc.vars[i][:]=data[i]

    nc.close()


class GenClm(GenCommon):
  '''
  Model climatology file generation

  Parameters
  ----------
  filename : netcdf file name
  grid : model netcdf grid file
  sparams : s-coordinates parameters, theta_s,theta_b, hc and NLevels,
                                      Vtransform, Vstretching

  **kargs:
  type : global attribute ('ROMS Climatology file')
  title : global attribute ('ROMS Climatology file')
  cycle : time variables cycle_length (False, attribute not created)
  tunits : time units ('days since 1970-01-01')
  one_time : if true (default), only one time variable is created
             (no tempo_time, zeta_time, ...)

  see GenCommon for additional kargs

  Example
  -------
  with time reference:
  ob=GenClm('clm.nc','grid.nc',(5,0,5,20),title='My clim file',
            tunits='seconds since 2000-01-01')

  without time reference (use cycle if climatological data):
  ... ,tunits='days',cycle=365)

  ob.create() # created the netcdf file
  ob.fill(data) # fills the netcdf file with data dict (may be provided
                # by prognostic.data2roms

  '''

  def __init__(self,filename,grid,sparams,**kargs):
    if not 'type'  in kargs.keys(): kargs['type']  = 'ROMS Climatology file'
    if not 'title' in kargs.keys(): kargs['title'] = 'ROMS Climatology file'

    GenCommon.__init__(self,filename,grid,sparams,**kargs)

    one_time = True # dont create several time variables, like temp_time, etc
    if 'one_time' in kargs.keys(): one_time=kargs['one_time']
    self.one_time=one_time

  def create(self):
    '''
    Creates model netcdf climatology file

    '''
    nc=netcdf.Pync(self.filename,'t',version=self.ncversion)

    # Dimensions:
    # horizontal and vertical dims:
    self.add_horiz_dims(nc)
    self.add_vert_dims(nc)
    # time dim:
    nc.add_dim('time',0)

    # Variables:
    v=nc.add_var('clim_time',np.dtype('d'),('time',))
    v.add_att('long_name','time for climatology')
    v.add_att('units',self.tunits)
    v.add_att('field','clim_time, scalar, series')
    if self.cycle: v.add_att('cycle_length',self.cycle)

    if not self.one_time:
      v=nc.add_var('temp_time',np.dtype('d'),('time',))
      v.add_att('long_name','time for temperature climatology')
      v.add_att('units',self.tunits)
      v.add_att('field','temp_time, scalar, series')
      if self.cycle: v.add_att('cycle_length',self.cycle)

    if not self.one_time:
      v=nc.add_var('salt_time',np.dtype('d'),('time',))
      v.add_att('long_name','time for salinity climatology')
      v.add_att('units',self.tunits)
      v.add_att('field','salt_time, scalar, series')
      if self.cycle: v.add_att('cycle_length',self.cycle)

    if not self.one_time:
      v=nc.add_var('u3d_time',np.dtype('d'),('time',))
      v.add_att('long_name','time for 3D U-momentum climatology')
      v.add_att('units',self.tunits)
      v.add_att('field','u3d_time, scalar, series')
      if self.cycle: v.add_att('cycle_length',self.cycle)

    if not self.one_time:
      v=nc.add_var('v3d_time',np.dtype('d'),('time',))
      v.add_att('long_name','time for 3D V-momentum climatology')
      v.add_att('units',self.tunits)
      v.add_att('field','v3d_time, scalar, series')
      if self.cycle: v.add_att('cycle_length',self.cycle)

    if not self.one_time:
      v=nc.add_var('u2d_time',np.dtype('d'),('time',))
      v.add_att('long_name','time for 2D U-momentum climatology')
      v.add_att('units',self.tunits)
      v.add_att('field','u2d_time, scalar, series')
      if self.cycle: v.add_att('cycle_length',self.cycle)

    if not self.one_time:
      v=nc.add_var('v2d_time',np.dtype('d'),('time',))
      v.add_att('long_name','time for 2D V-momentum climatology')
      v.add_att('units',self.tunits)
      v.add_att('field','v2d_time, scalar, series')
      if self.cycle: v.add_att('cycle_length',self.cycle)

    if not self.one_time:
      v=nc.add_var('zeta_time',np.dtype('d'),('time',))
      v.add_att('long_name','time for free surface climatology')
      v.add_att('units',self.tunits)
      v.add_att('field','zeta_time, scalar, series')
      if self.cycle: v.add_att('cycle_length',self.cycle)

    v=nc.add_var('temp',np.dtype('d'),('time','s_rho','eta_rho','xi_rho'))
    v.add_att('long_name','potential temperature climatology')
    v.add_att('units','Celsius')
    if self.one_time: v.add_att('time','clim_time')
    else:             v.add_att('time','temp_time')
    v.add_att('coordinates','lon_rho lat_rho s_rho temp_time')
    v.add_att('field','temperature, scalar, series')

    v=nc.add_var('salt',np.dtype('d'),('time','s_rho','eta_rho','xi_rho'))
    v.add_att('long_name','salinity climatology')
    v.add_att('units','PSU')
    if self.one_time: v.add_att('time','clim_time')
    else:             v.add_att('time','salt_time')
    v.add_att('coordinates','lon_rho lat_rho s_rho salt_time')
    v.add_att('field','salinity, scalar, series')

    v=nc.add_var('u',np.dtype('d'),('time','s_rho','eta_u','xi_u'))
    v.add_att('long_name','u-momentum component climatology')
    v.add_att('units','metre second-1')
    if self.one_time: v.add_att('time','clim_time')
    else:             v.add_att('time','u3d_time')
    v.add_att('coordinates','lon_u lat_u s_rho u3d_time')
    v.add_att('field','u-velocity, scalar, series')

    v=nc.add_var('v',np.dtype('d'),('time','s_rho','eta_v','xi_v'))
    v.add_att('long_name','v-momentum component climatology')
    v.add_att('units','metre second-1')
    if self.one_time: v.add_att('time','clim_time')
    else:             v.add_att('time','v3d_time')
    v.add_att('coordinates','lon_v lat_v s_rho v3d_time')
    v.add_att('field','v-velocity, scalar, series')

    v=nc.add_var('ubar',np.dtype('d'),('time','eta_u','xi_u'))
    v.add_att('long_name','vertically integrated  u-momentum component climatology')
    v.add_att('units','metre second-1')
    if self.one_time: v.add_att('time','clim_time')
    else:             v.add_att('time','u2d_time')
    v.add_att('coordinates','lon_u lat_u u2d_time')
    v.add_att('field','ubar-velocity, scalar, series')

    v=nc.add_var('vbar',np.dtype('d'),('time','eta_v','xi_v'))
    v.add_att('long_name','vertically integrated  v-momentum component climatology')
    v.add_att('units','metre second-1')
    if self.one_time: v.add_att('time','clim_time')
    else:             v.add_att('time','v2d_time')
    v.add_att('coordinates','lon_v lat_v v2d_time')
    v.add_att('field','vbar-velocity, scalar, series')

    v=nc.add_var('zeta',np.dtype('d'),('time','eta_rho','xi_rho'))
    v.add_att('long_name','free-surface climatology')
    v.add_att('units','metre')
    if self.one_time: v.add_att('time','clim_time')
    else:             v.add_att('time','zeta_time')
    v.add_att('coordinates','lon_rho lat_rho zeta_time')
    v.add_att('field','free-surface, scalar, series')


    # sparams vars:
    v=nc.add_var('theta_s',np.dtype('d'),())
    v.add_att('long_name','S-coordinate surface control parameter')

    v=nc.add_var('theta_b',np.dtype('d'),())
    v.add_att('long_name','S-coordinate bottom control parameter')

    v=nc.add_var('Tcline',np.dtype('d'),())
    v.add_att('long_name','S-coordinate surface/bottom layer width')
    v.add_att('units','metre')

    v=nc.add_var('hc',np.dtype('d'),())
    v.add_att('long_name','S-coordinate parameter, critical depth')
    v.add_att('units','metre')

    v=nc.add_var('Vtransform',np.dtype('i'),())
    v.add_att('long_name','vertical terrain-following transformation equation')

    v=nc.add_var('Vstretching',np.dtype('i'),())
    v.add_att('long_name','vertical terrain-following stretching function')


    # Global Attributes:
    nc.add_att('type',self.type)
    nc.add_att('title',self.title)
    nc.add_att('grd_file',realpath(self.grid))
    nc.add_att('history','ROMS clm file, '+ctime())
    nc.add_att('author',cb.username()[1]+', '+cb.machinename())

    # extra attrs:
    for i in self.attr.keys(): nc.add_att(i,self.attr[i])

    # fill sparams:
    nc.vars['theta_s'][:]=self.sparams[0]
    nc.vars['theta_b'][:]=self.sparams[1]
    nc.vars['Tcline'][:]=self.sparams[2]
    nc.vars['hc'][:]=self.sparams[2]
    nc.vars['Vtransform'][:]=self.sparams[4]
    nc.vars['Vstretching'][:]=self.sparams[5]

    nc.close()


  def fill(self,data,tind='next',quiet=1):
    '''
    Fills model netcdf climatology file (data can be provided by
    prognostic.data2roms

    '''
    nc=netcdf.Pync(self.filename,'w')
    if tind=='next': tind=nc.dims['time']


    if not quiet: print 'filling clm file %s' % self.filename

    # about time:
    try:
      date=dts.parse_date(data['date'])
      time=netcdf.date2num(date,self.tunits)
    except:
      time=data['date'] # date as number!

    for i in nc.varnames:
      if i.endswith('time'):
        if not quiet: print '  -- %s tind=%d %f' % (i,tind,time)
        nc.vars[i][tind]=time


    names='temp','salt','u','v','ubar','vbar','zeta'
    for i in names:
      if not quiet: print '  %s' % i
      nc.vars[i][tind,...]=data[i]

    nc.close()


class GenBry(GenCommon):
  '''
  Model boundary conditions file generation

  Parameters
  ----------
  filename : netcdf file name
  grid : model netcdf grid file
  sparams : s-coordinates parameters, theta_s,theta_b, hc and NLevels,
                                      Vtransform, Vstretching

  **kargs:
  type : global attribute ('ROMS Boundary forcing file')
  title : global attribute ('ROMS Boundary forcing file')
  cycle : time variables cycle_length (False, attribute not created)
  tunits : time units ('days since 1970-01-01')
  obc : open boundaries, the 4 boundaries are used by default ('nsew',
        for north, south, east and west)
  addxz : include distance x depth variables in the netcdf file (True)

  see GenCommon for additional kargs

  Example
  -------
  with time reference:
  ob=GenBry('bry.nc','grid.nc',(5,0,5,20),title='My bry file',
            tunits='seconds since 2000-01-01',obc='sew')

  without time reference (use cycle if climatological data):
  ... ,tunits='days',cycle=365)

  ob.create() # created the netcdf file
  ob.fill(data) # fills the netcdf file with data dict (may be provided
                # by prognostic.data2romsbry

  '''

  def __init__(self,filename,grid,sparams,**kargs):
    if not 'type'  in kargs.keys(): kargs['type']  = 'ROMS Boundary forcing file'
    if not 'title' in kargs.keys(): kargs['title'] = 'ROMS Boundary forcing file'

    GenCommon.__init__(self,filename,grid,sparams,**kargs)

    cycle = False
    obc   = 'nsew'
    addxz = True
    if 'cycle' in kargs.keys(): cycle = kargs['cycle']
    if 'obc'   in kargs.keys(): obc   = kargs['obc']
    if 'addxz' in kargs.keys(): addxz = kargs['addxz']

    self.cycle = cycle
    self.obc   = obc
    self.addxz = addxz


  def create(self):
    '''
    Creates model netcdf boundary conditions file

    '''

    nc=netcdf.Pync(self.filename,'t',version=self.ncversion)

    # Dimensions:
    # horizontal and vertical dims:
    self.add_horiz_dims(nc)
    self.add_vert_dims(nc)
    # time dim:
    nc.add_dim('time',0)

    # Variables:
    v=nc.add_var('bry_time',np.dtype('d'),('time',))
    v.add_att('long_name','open boundary conditions time')
    v.add_att('units',self.tunits)
    v.add_att('field','bry_time, scalar, series')
    if self.cycle: v.add_att('cycle_length',self.cycle)

    v=nc.add_var('temp_time',np.dtype('d'),('time',))
    v.add_att('long_name','potential temperature time')
    v.add_att('units',self.tunits)
    v.add_att('field','temp_time, scalar, series')
    if self.cycle: v.add_att('cycle_length',self.cycle)

    v=nc.add_var('salt_time',np.dtype('d'),('time',))
    v.add_att('long_name','salinity time')
    v.add_att('units',self.tunits)
    v.add_att('field','salt_time, scalar, series')
    if self.cycle: v.add_att('cycle_length',self.cycle)

    v=nc.add_var('u3d_time',np.dtype('d'),('time',))
    v.add_att('long_name','3D momentum time')
    v.add_att('units',self.tunits)
    v.add_att('field','u3d_time, scalar, series')
    if self.cycle: v.add_att('cycle_length',self.cycle)

    #v=nc.add_var('v3d_time',np.dtype('d'),('time',))
    #v.add_att('long_name','3D V-momentum time')
    #v.add_att('units',self.tunits)
    #v.add_att('field','v3d_time, scalar, series')
    #if self.cycle: v.add_att('cycle_length',self.cycle)

    v=nc.add_var('u2d_time',np.dtype('d'),('time',))
    v.add_att('long_name','2D momentum time')
    v.add_att('units',self.tunits)
    v.add_att('field','u2d_time, scalar, series')
    if self.cycle: v.add_att('cycle_length',self.cycle)

    #v=nc.add_var('v2d_time',np.dtype('d'),('time',))
    #v.add_att('long_name','2D V-momentum time')
    #v.add_att('units',self.tunits)
    #v.add_att('field','v2d_time, scalar, series')
    #if self.cycle: v.add_att('cycle_length',self.cycle)

    v=nc.add_var('zeta_time',np.dtype('d'),('time',))
    v.add_att('long_name','free surface time')
    v.add_att('units',self.tunits)
    v.add_att('field','zeta_time, scalar, series')
    if self.cycle: v.add_att('cycle_length',self.cycle)

    names={'n':'north','s':'south','e':'east','w':'west'}
    Names={'n':'northern','s':'southern','e':'eastern','w':'western'}

    Tdims={'n':('time','s_rho','xi_rho'),'s':('time','s_rho','xi_rho'),
           'e':('time','s_rho','eta_rho'),'w':('time','s_rho','eta_rho')}

    Udims={'n':('time','s_rho','xi_u'),'s':('time','s_rho','xi_u'),
           'e':('time','s_rho','eta_u'),'w':('time','s_rho','eta_u')}

    Vdims={'n':('time','s_rho','xi_v'),'s':('time','s_rho','xi_v'),
           'e':('time','s_rho','eta_v'),'w':('time','s_rho','eta_v')}


    for o in 'nsew':
      if o in self.obc:
        v=nc.add_var('temp_'+names[o],np.dtype('d'),Tdims[o])
        v.add_att('long_name','potential temperature '+Names[o]+' boundary condition')
        v.add_att('units','Celsius')
        v.add_att('time','temp_time')
        v.add_att('field','temp_'+names[o]+', scalar, series')

        v=nc.add_var('salt_'+names[o],np.dtype('d'),Tdims[o])
        v.add_att('long_name','salinity '+Names[o]+' boundary condition')
        v.add_att('units','Celsius')
        v.add_att('time','salt_time')
        v.add_att('field','salt_'+names[o]+', scalar, series')

        v=nc.add_var('u_'+names[o],np.dtype('d'),Udims[o])
        v.add_att('long_name','3D u-momentum '+Names[o]+' boundary condition')
        v.add_att('units','metre second-1')
        v.add_att('time','u3d_time')
        v.add_att('field','u_'+names[o]+', scalar, series')

        v=nc.add_var('v_'+names[o],np.dtype('d'),Vdims[o])
        v.add_att('long_name','3D v-momentum '+Names[o]+' boundary condition')
        v.add_att('units','metre second-1')
        v.add_att('time','u3d_time')
        v.add_att('field','v_'+names[o]+', scalar, series')

        v=nc.add_var('ubar_'+names[o],np.dtype('d'),Udims[o][::2])
        v.add_att('long_name','2D u-momentum '+Names[o]+' boundary condition')
        v.add_att('units','metre second-1')
        v.add_att('time','u2d_time')
        v.add_att('field','ubar_'+names[o]+', scalar, series')

        v=nc.add_var('vbar_'+names[o],np.dtype('d'),Vdims[o][::2])
        v.add_att('long_name','2D v-momentum '+Names[o]+' boundary condition')
        v.add_att('units','metre second-1')
        v.add_att('time','u2d_time')
        v.add_att('field','vbar_'+names[o]+', scalar, series')

        v=nc.add_var('zeta_'+names[o],np.dtype('d'),Tdims[o][::2])
        v.add_att('long_name','free surface '+Names[o]+' boundary condition')
        v.add_att('units','metre')
        v.add_att('time','zeta_time')
        v.add_att('field','zeta_'+names[o]+', scalar, series')

        if self.addxz:
          v=nc.add_var('dist_'+names[o],np.dtype('f'),Tdims[o][1:])
          v.add_att('long_name','distance of '+Names[o]+' boundary at rho points')
          v.add_att('units','metre')
          v.add_att('field','dist_'+names[o]+', scalar, series')

          v=nc.add_var('distu_'+names[o],np.dtype('f'),Udims[o][1:])
          v.add_att('long_name','distance of '+Names[o]+' boundary at u points')
          v.add_att('units','metre')
          v.add_att('field','distu_'+names[o]+', scalar, series')

          v=nc.add_var('distv_'+names[o],np.dtype('f'),Vdims[o][1:])
          v.add_att('long_name','distance of '+Names[o]+' boundary at v points')
          v.add_att('units','metre')
          v.add_att('field','distv_'+names[o]+', scalar, series')

          v=nc.add_var('depth_'+names[o],np.dtype('f'),Tdims[o])
          v.add_att('long_name','depths of '+Names[o]+' boundary at rho points')
          v.add_att('units','metre')
          v.add_att('field','depth_'+names[o]+', scalar, series')

          v=nc.add_var('depthu_'+names[o],np.dtype('f'),Udims[o])
          v.add_att('long_name','depths of '+Names[o]+' boundary at u points')
          v.add_att('units','metre')
          v.add_att('field','depthu_'+names[o]+', scalar, series')

          v=nc.add_var('depthv_'+names[o],np.dtype('f'),Vdims[o])
          v.add_att('long_name','depths of '+Names[o]+' boundary at v points')
          v.add_att('units','metre')
          v.add_att('field','depthv_'+names[o]+', scalar, series')


    # sparams vars:
    v=nc.add_var('theta_s',np.dtype('d'),())
    v.add_att('long_name','S-coordinate surface control parameter')

    v=nc.add_var('theta_b',np.dtype('d'),())
    v.add_att('long_name','S-coordinate bottom control parameter')

    v=nc.add_var('Tcline',np.dtype('d'),())
    v.add_att('long_name','S-coordinate surface/bottom layer width')
    v.add_att('units','metre')

    v=nc.add_var('hc',np.dtype('d'),())
    v.add_att('long_name','S-coordinate parameter, critical depth')
    v.add_att('units','metre')

    v=nc.add_var('Vtransform',np.dtype('i'),())
    v.add_att('long_name','vertical terrain-following transformation equation')

    v=nc.add_var('Vstretching',np.dtype('i'),())
    v.add_att('long_name','vertical terrain-following stretching function')


    # Global Attributes:
    nc.add_att('type',self.type)
    nc.add_att('title',self.title)
    nc.add_att('grd_file',realpath(self.grid))
    nc.add_att('history','ROMS bry file, '+ctime())
    nc.add_att('obc',self.obc)
    nc.add_att('author',cb.username()[1]+', '+cb.machinename())

    # extra attrs:
    for i in self.attr.keys(): nc.add_att(i,self.attr[i])

    # fill sparams:
    nc.vars['theta_s'][:]=self.sparams[0]
    nc.vars['theta_b'][:]=self.sparams[1]
    nc.vars['Tcline'][:]=self.sparams[2]
    nc.vars['hc'][:]=self.sparams[2]
    nc.vars['Vtransform'][:]=self.sparams[4]
    nc.vars['Vstretching'][:]=self.sparams[5]

    nc.close()


  def fill(self,data,tind='next',quiet=1):
    '''
    Fills model netcdf boundary conditions (data can be provided by
    prognostic.data2romsbry

    '''
    nc=netcdf.Pync(self.filename,'w')
    if tind=='next': tind=nc.dims['time']


    if not quiet: print 'filling bry file %s' % self.filename

    # about time:
    try:
      date=dts.parse_date(data['date'])
      time=netcdf.date2num(date,self.tunits)
    except:
      time=data['date'] # date as number!

    for i in nc.varnames:
      if i.endswith('time'):
        if not quiet: print '  -- %s tind=%d %f' % (i,tind,time)
        nc.vars[i][tind]=time

    names='temp','salt','u','v','ubar','vbar','zeta'
    if self.addxz:
      names=list(names)+['dist','distu','distv','depth','depthu','depthv']

    for i in names:
      for j in 'north','south','east','west':
        vname=i+'_'+j
        if vname in nc.varnames:
          if vname.startswith('dist'):
            if tind==0:
              if not quiet: print '  %s' % vname
              nc.vars[vname][:]=data[vname] # not time dependent
          else:
            if not quiet: print '  %s' % vname
            nc.vars[vname][tind,...]=data[vname]

    nc.close()


class GenBlk(GenCommon):
  '''
  Model bulk forcing file generation

  Parameters
  ----------
  filename : netcdf file name
  grid : model netcdf grid file

  **kargs:
  type : global attribute ('ROMS Bulk forcing file')
  title : global attribute ('ROMS Bulk forcing file')
  cycle : time variables cycle_length (False, attribute not created)
  tunits : time units ('seconds since 2000-01-01')

  see GenCommon for additional kargs

  Example
  -------
  with time reference:
  ob=GenBlk('blk.nc','grid.nc',title='My bulk file',
            tunits='seconds since 2000-01-01')

  without time reference (use cycle if climatological data):
  ... ,tunits='days',cycle=365)

  ob.create() # created the netcdf file
  ob.fill(data) # fills the netcdf file with data dict (may be provided
                # by surface.data2romsblk

  '''

  def __init__(self,filename,grid,**kargs):
    if not 'type'  in kargs.keys(): kargs['type']  = 'ROMS Bulk forcing file'
    if not 'title' in kargs.keys(): kargs['title'] = 'ROMS Bulk forcing file'

    GenCommon.__init__(self,filename,grid,**kargs)


  def create(self,model='roms',original=False):
    '''
    Creates netcdf bulk forcing file
    for model 'roms' or 'roms-agrif'
    '''
    nc=netcdf.Pync(self.filename,'t',version=self.ncversion)

    # Dimensions:
    # horizontal dims:
    self.add_horiz_dims(nc)
    # time dim:
    nc.add_dim('time',0)

    # original data dims:
    addOriginal= not original is False
    if addOriginal:
      nc.add_dim('x_original',original[1])
      nc.add_dim('y_original',original[0])

    # Variables:
    if model=='roms': vname='time'
    elif model=='roms-agrif': vname='bulk_time'
    v=nc.add_var(vname,np.dtype('d'),('time',))
    v.add_att('long_name','bulk formulation atmospheric forcing time')
    v.add_att('units',self.tunits)
    if self.cycle: v.add_att('cycle_length',self.cycle)

    if addOriginal:
      nc.add_var('x_original',np.dtype('d'),('y_original', 'x_original'))
      nc.add_var('y_original',np.dtype('d'),('y_original', 'x_original'))

    if model=='roms': vname='Tair'
    elif model=='roms-agrif': vname='tair'
    v=nc.add_var(vname,np.dtype('d'),('time','eta_rho', 'xi_rho'))
    v.add_att('long_name','surface air temperature')
    v.add_att('units','Celsius')
    v.add_att('time','time')

    if addOriginal:
      v=nc.add_var('tair_original',np.dtype('d'),('time','y_original', 'x_original'))
      v.add_att('units','Celsius')

    if model=='roms': vname,vunits='Pair','millibar'
    elif model=='roms-agrif': vname,vunits='pres','Pascal'
    v=nc.add_var(vname,np.dtype('d'),('time','eta_rho', 'xi_rho'))
    v.add_att('long_name','surface pressure')
    v.add_att('units',vunits)
    v.add_att('time','time')

    if addOriginal:
      v=nc.add_var('pres_original',np.dtype('d'),('time','y_original', 'x_original'))
      v.add_att('units','Pascal')

    if model=='roms': vname,vunits='Qair','percentage'
    elif model=='roms-agrif': vname,vunits='rhum','fraction'
    v=nc.add_var(vname,np.dtype('d'),('time','eta_rho', 'xi_rho'))
    v.add_att('long_name','relative humidity')
    v.add_att('units',vunits)
    v.add_att('time','time')

    if addOriginal:
      v=nc.add_var('rhum_original',np.dtype('d'),('time','y_original', 'x_original'))
      v.add_att('units','fraction')

    if model=='roms': vname,vunits='rain','kg m-2 s-1'
    elif model=='roms-agrif': vname,vunits='prate','cm day-1'
    v=nc.add_var(vname,np.dtype('d'),('time','eta_rho', 'xi_rho'))
    v.add_att('long_name','precipitation rate')
    v.add_att('units',vunits)
    v.add_att('time','time')

    if addOriginal:
      v=nc.add_var('prate_original',np.dtype('d'),('time','y_original', 'x_original'))
      v.add_att('units','cm day-1')

    if model=='roms': vname,vlname='lwrad','net longwave radiation flux'
    elif model=='roms-agrif': vname,vlname='radlw','outgoing longwave radiation'
    v=nc.add_var(vname,np.dtype('d'),('time','eta_rho', 'xi_rho'))
    v.add_att('long_name',vlname)
    v.add_att('units','Watts metre-2')
    v.add_att('time','time')
    if model=='roms':
      v.add_att('positive_value','downward flux, heating')
      v.add_att('negative_value','upward flux, cooling')
    elif model=='roms-agrif': # opposite to ROMS !!!
      v.add_att('positive','upward flux, cooling water')

    if addOriginal:
      v=nc.add_var('radlw_original',np.dtype('d'),('time','y_original', 'x_original'))
      v.add_att('units','Watts metre-2')
      if model=='roms':
        v.add_att('negative_value','upward flux, cooling')
      elif model=='roms-agrif':
        v.add_att('positive','upward flux, cooling water')

    if model=='roms': vname='lwrad_down'
    elif model=='roms-agrif': vname='dlwrf'
    v=nc.add_var(vname,np.dtype('d'),('time','eta_rho', 'xi_rho'))
    v.add_att('long_name','downward longwave radiation')
    v.add_att('units','Watts metre-2')
    v.add_att('time','time')
    if model=='roms':
      v.add_att('positive_value','downward flux, heating')
      v.add_att('negative_value','upward flux, cooling')
    elif model=='roms-agrif': # opposite to ROMS !!!
      v.add_att('positive','upward flux, cooling water')

    if addOriginal:
      v=nc.add_var('dlwrf_original',np.dtype('d'),('time','y_original', 'x_original'))
      v.add_att('units','Watts metre-2')
      if model=='roms':
        v.add_att('negative_value','upward flux, cooling')
      elif model=='roms-agrif':
        v.add_att('positive','upward flux, cooling water')

    if model=='roms': vname='swrad'
    elif model=='roms-agrif': vname='radsw'
    v=nc.add_var(vname,np.dtype('d'),('time','eta_rho', 'xi_rho'))
    v.add_att('long_name','shortwave radiation')
    v.add_att('units','Watts metre-2')
    v.add_att('time','time')
    v.add_att('positive_value','downward flux, heating')
    v.add_att('negative_value','upward flux, cooling')

    if addOriginal:
      v=nc.add_var('radsw_original',np.dtype('d'),('time','y_original', 'x_original'))
      v.add_att('units','Watts metre-2')

    if model=='roms': vname='Uwind'
    elif model=='roms-agrif': vname='uwnd'
    v=nc.add_var(vname,np.dtype('d'),('time','eta_rho', 'xi_rho'))
    v.add_att('long_name','u-wind')
    v.add_att('units','metre second-1')
    v.add_att('time','time')

    if model=='roms': vname='Vwind'
    elif model=='roms-agrif': vname='vwnd'
    v=nc.add_var(vname,np.dtype('d'),('time','eta_rho', 'xi_rho'))
    v.add_att('long_name','v-wind')
    v.add_att('units','metre second-1')
    v.add_att('time','time')

    if addOriginal:
      v=nc.add_var('uwnd_original',np.dtype('d'),('time','y_original', 'x_original'))
      v.add_att('units','metre second-1')
      v=nc.add_var('vwnd_original',np.dtype('d'),('time','y_original', 'x_original'))
      v.add_att('units','metre second-1')

    # cloud cover (for roms):
    if model=='roms' or 1: # add it anyway...
      v=nc.add_var('cloud',np.dtype('d'),('time','eta_rho', 'xi_rho'))
      v.add_att('long_name','cloud fraction')
      v.add_att('units','nondimensional')
      v.add_att('time','time')

    if addOriginal:
      v=nc.add_var('cloud_original',np.dtype('d'),('time','y_original', 'x_original'))
      v.add_att('units','cloud fraction')

    #  next only required by roms-agrif (wspd,sustr,svstr)
    v=nc.add_var('wspd',np.dtype('d'),('time','eta_rho', 'xi_rho'))
    v.add_att('long_name','wind speed 10m')
    v.add_att('units','metre second-1')
    v.add_att('time','time')

    v=nc.add_var('sustr',np.dtype('d'),('time','eta_u', 'xi_u'))
    v.add_att('long_name','surface u-momentum stress')
    v.add_att('units','Newton metre-2')
    v.add_att('time','time')

    v=nc.add_var('svstr',np.dtype('d'),('time','eta_v', 'xi_v'))
    v.add_att('long_name','surface v-momentum stress')
    v.add_att('units','Newton metre-2')
    v.add_att('time','time')

    if addOriginal:
      v=nc.add_var('wspd_original',np.dtype('d'),('time','y_original', 'x_original'))
      v.add_att('units','metre second-1')
      v=nc.add_var('sustr_original',np.dtype('d'),('time','y_original', 'x_original'))
      v.add_att('units','Newton metre-2')
      v=nc.add_var('svstr_original',np.dtype('d'),('time','y_original', 'x_original'))
      v.add_att('units','Newton metre-2')


    # Global Attributes:
    nc.add_att('type',self.type)
    nc.add_att('title',self.title)
    nc.add_att('grd_file',realpath(self.grid))
    nc.add_att('history','ROMS blk file, '+ctime())
    nc.add_att('author',cb.username()[1]+', '+cb.machinename())

    # extra attrs:
    for i in self.attr.keys(): nc.add_att(i,self.attr[i])

    nc.close()


  def fill(self,data,tind='next',quiet=1):
    '''
    Fills model netcdf bulk forcing file

    '''
    nc=netcdf.Pync(self.filename,'w')
    if tind=='next': tind=nc.dims['time']


    if not quiet: print 'filling blk file %s' % self.filename

    # about time:
    try:
      date=dts.parse_date(data['date'])
      time=netcdf.date2num(date,self.tunits)
    except:
      time=data['date'] # date as number!

    for i in nc.varnames:
      if i.endswith('time'):
        if not quiet: print '  -- %s tind=%d %f' % (i,tind,time)
        nc.vars[i][tind]=time


    if 'Tair' in nc.varnames: # roms
      names=('Tair','tair'),('Pair','pres'),('Qair','rhum'),('rain','prate'),\
            ('swrad','radsw'),('lwrad','radlw'),('Uwind','uwnd'),('Vwind','vwnd'),\
            'sustr','svstr','wspd','cloud',('lwrad_down','dlwrf')
      if not 'tair' in data.keys(): # assuming data has roms (not agrif) varnames:
        names='Tair','Pair','Qair','rain','swrad','lwrad','Uwind','Vwind','sustr','svstr','wspd','cloud','lwrad_down'

    elif 'tair' in nc.varnames: # roms-agrif
      names='tair','pres','rhum','prate','radlw','radsw','dlwrf','uwnd',\
            'vwnd','wspd','sustr','svstr',\
            'cloud' # not used, but add it anyway

    for i in names:
      if isinstance(i,basestring): filev,datav=i,i
      else:  filev,datav=i

      if datav not in data.keys():
        if not quiet: print '  Warning: data key %s not present' % datav
      else:
        if not quiet: print '  %s (%s) min=%8.3f max=%8.3f' % (filev.ljust(7),datav.ljust(7),
                                                               data[datav].min(),data[datav].max())
        nc.vars[filev][tind,...]=data[datav]

      # fill original data:
      orig=datav+'_original'
      if orig in data.keys() and not orig in nc.varnames:
        if not quiet: print '  Warning: original data will not be written %s' % orig
      elif not orig in data.keys() and orig in nc.varnames:
        if not quiet: print '  Warning: original data not present %s' % orig
      elif orig in data.keys() and orig in nc.varnames:
        if not quiet: print '  %s  min=%8.3f max=%8.3f' % (orig.ljust(7+9),
                                                          data[orig].min(),data[orig].max())
        nc.vars[orig][tind,...]=data[orig]

    # fill original x,y:
    if tind==0 and 'x_original' in data.keys() and 'x_original' in nc.varnames:
      if not quiet: print '  filling x,y original'
      nc.vars['x_original'][:]=data['x_original']
      nc.vars['y_original'][:]=data['y_original']

    nc.close()


  def update_wind(self,data,keepOld=True,quiet=0):
    # store old at first iteration


    nc=netcdf.Pync(self.filename,'w')
    time=netcdf.nctime(self.filename,'time')
    tind,=np.where(time==data['date'])

    if tind.size: tind=tind[0]
    else:
      if not quiet: print 'cannot find time=%s' % data['date'].isoformat(' ')
      return


    if 'Uwind' in nc.varnames: # roms
      names=('Uwind','uwnd'),('Vwind','vwnd'),'sustr','svstr','wspd'
    else: # roms-agrif
      names='uwnd','vwnd','sustr','svstr','wspd'

    # store old wind data:
    if False and tind==0: # not using for now.... too slow!
                          # TODO: remove keepOld form argument !?
      if 'x_wind' in data.keys():
        if not quiet: print 'adding x_wind, y_wind'
        nc.add_dim('x_wind_original',data['x_wind'].shape[1])
        nc.add_dim('y_wind_original',data['x_wind'].shape[0])

        nc.add_var('x_wind_original',data['x_wind'].dtype,
                   ('y_wind_original','x_wind_original'))
        nc.add_var('y_wind_original',data['y_wind'].dtype,
                   ('y_wind_original','x_wind_original'))

        nc.vars['x_wind_original'][...]=data['x_wind']
        nc.vars['y_wind_original'][...]=data['y_wind']


      for i in names:
        if isinstance(i,basestring): filev,datav=i,i
        else:  filev,datav=i

        newName = filev+'_old'
        newType = nc.vars[filev].dtype()
        newDims = nc.vars[filev].dims.keys()

        if  not quiet: print '  creating var %s' % newName
        nc.add_var(newName,newType,newDims)

        if  not quiet: print '    filling var %s' % newName
        nc.vars[newName][...]=nc.vars[filev][...]

        # also store old original data, if present:
        orig=datav+'_original'
        if orig in nc.varnames:
          newName = orig+'_old'
          newType = nc.vars[orig].dtype()
          newDims = nc.vars[orig].dims.keys()

          if  not quiet: print '  creating var %s' % newName
          nc.add_var(newName,newType,newDims)

          if  not quiet: print '    filling var %s' % newName
          nc.vars[newName][...]=nc.vars[orig][...]


    # store new wind data:
    for i in names:
      if isinstance(i,basestring): filev,datav=i,i 
      else:  filev,datav=i

      if not quiet: print '  %s (%s) min=%8.3f max=%8.3f' % (filev.ljust(7),datav.ljust(7),
                                                               data[datav].min(),data[datav].max())
      nc.vars[filev][tind,...]=data[datav]

      # also store new wind original data, if present:
      orig=datav+'_original'
      if not orig in data.keys(): continue
      if not quiet: print '  %s  min=%8.3f max=%8.3f' % (orig.ljust(7+9),
                                                        data[orig].min(),data[orig].max())
      nc.vars[orig][tind,...]=data[datav]


    # adding global atts:
    if 'attr' in data.keys():
      for i in data['attr'].keys(): nc.add_att(i,data['attr'][i])


    nc.close()
