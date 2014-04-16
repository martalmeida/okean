from okean import cookbook as cbt
import os
import numpy as np
import nctypes

# the default intrefaces by be given by the environment variable
# PY_NETCDF_INTERFACES, if not, will be used:
default_interfaces    = 'netcdf4','pycdf','pupynere','scientific'

maxint=cbt.maxint()

def __load_pycdf():
  pnc=False
  try: import pycdf as pnc
  except: pass
  return pnc

def __load_scientific():
  pnc=False
  try: from Scientific.IO import NetCDF as pnc
  except: pass
  return pnc

def __load_netcdf4():
  pnc=False
  try: import netCDF4 as pnc
  except: pass
  return pnc

def __load_pupynere():
  pnc=False
  try: import pupynere as pnc
  except: pass
  return pnc

def load_interface():
  DefaultInterfaces=os.environ.get('PY_NETCDF_INTERFACES')
  if DefaultInterfaces: default=DefaultInterfaces
  else: default=default_interfaces

  for i in default_interfaces:
    pnc=eval('__load_'+i+'()')
    if pnc: return pnc,i

  return False,False


pnc,Interface=load_interface()

def slice2str(elem):
  if not isinstance(elem,tuple): elem=(elem,)
  e=''
  for i in elem:
    if isinstance(i,int): e+=str(i)+','
    elif isinstance(i,slice):
      i0=i.start
      if i0 is None or i0==0: i0=''
      i1=i.stop
      if i1 is None or i1==maxint: i1=''
      i2=i.step
      if i2 is None: i2=''
      e+=str(i0) +':'+ str(i1) +':' +str(i2)+','



class Common:
  def __init__(self):
    # setting atts and attnames:
    self.atts=self.__atts()
    self.attnames=self.atts.keys()

  #def 

  def add_att(self,attname,attvalue,atttype=False):
    '''
    Creates file or variable attribute
    atttype can be a netcdf type name, numpy type name, dtype or
    numeric typecode
    '''
    if self._interface=='pycdf':
      if atttype: Type=nctypes.type_2pycdf(atttype)
      else:       Type=nctypes.type_var2pycdf(attvalue)
      a=self._nc.attr(attname)
      a.put(Type,attvalue)

    elif self._interface in ('scientific','netcdf4'):
      setattr(self._nc,attname,attvalue)

    newatt=Pyncatt(self._nc,attname,attvalue,
                   atttype=atttype,ncversion=self.ncversion,
                   interface=self._interface)

    # update self.atts and self.attnames:
    self.atts[attname]=newatt

    # update ncdump_info:
    if not isinstance(self,Pyncvar):
      if not self._ncdump_info is False:
        self._ncdump_info['attributes'][attname]={}

    return newatt


  def __atts(self):
    '''Returns netcdf ordered attributes, with types and values'''
    att=cbt.odict()
    if self._interface=='pycdf':
      for i in range(self._nc.inq_natts()):
        a=self._nc.attr(i)
        # ps: a.get() produces a segmentation fault for large
        # strings, so:
        L= a.inq_len()
        if L>=250: val='PYCDF ERROR: CANNOT READ !!!!'
        else: attvalue=a.get()

        # find nctype:
        nctype=nctypes.type_pycdf2nc(a.inq_type())

        att[a.inq_name()]=Pyncatt(self._nc,a.inq_name(),attvalue,
                                  nctype=nctype,ncversion=self.ncversion,
                                  interface=self._interface)

    elif self._interface=='scientific':
      a=self._nc.__dict__
      # get ordered atts:
      if not self._ncdump_info is False:
        oatts=self._ncdump_info['vattributes'][self.varname]
      else: oatts=False

      if oatts: keys=oatts
      else: keys= a.keys()

      for k in keys:
        attvalue=a[k]
        att[k]=Pyncatt(self._nc,k,attvalue,ncversion=self.ncversion,
                       interface=self._interface)

    elif self._interface=='netcdf4':
      aname=self._nc.ncattrs() # already ordered!
      for a in aname:
        attvalue=getattr(self._nc,a)
        att[a]=Pyncatt(self._nc,a,attvalue,ncversion=self.ncversion,
                       interface=self._interface)
    return att


  def att_rename(self,oldname,newname):
    '''Renames file or variable attribute'''
    if self._interface=='pycdf':
      self._nc.attr(oldname).rename(newname)
    elif self._interface in ('scientific','netcdf4'):
      setattr(self._nc,newname,getattr(self._nc,oldname))
      delattr(self._nc,oldname)

    # update self.atts:
    self.atts.rename_key(oldname,newname)

    # update ncdump info:
    if not self._ncdump_info is False:
      self._ncdump_info['attributes'].rename_key(oldname,newname)


class Pyncatt:
  def __init__(self,_nc,name,value,ncversion,
               nctype=False,dtype=False,atttype=False,interface=None):

    self._nc       = _nc
    self.name      = name
    self.value     = value
    self.ncversion = ncversion
    self.nctype    = nctype
    self.dtype     = dtype
    self.atttype   = atttype
    self._interface = interface

    self.set_types()

  def set_types(self):

    # find nctype and dtype:
    if self.atttype:
      # nctype:
      self.nctype = nctypes.type_2nc(self.atttype,ncver=self.ncversion)

      # dtype:
      isstr=isinstance(self.value,basestring)
      isbool=isinstance(self.value,bool)
      if isstr: strlen=len(str)
      else: strlen=None
      self.dtype  = nctypes.type_2numpy(self.atttype,isstr=isstr,strlen=strlen,
                               isbool=isbool)
    else:
      if not self.nctype: self.nctype = nctypes.type_var2nc(self.value,self.ncversion)
      if not self.dtype: self.dtype = nctypes.type_var2numpy(self.value)


  def __repr__(self):
    s=''
    s+='nctype: '+str(self.nctype)
    s+='\ndtype:  '+str(self.dtype)
    try:
      s+='\nvalue:  '+str(self.value)
    except: s+='\nvalue:  '+unicode(self.value)

    s='%s   %s  %s' % (self.name,str(self.dtype),str(self.value))
    return s

  def update(self,value,atttype=False):
      self.value=value
      self.atttype=atttype
      self.set_types()

      if self._interface=='pycdf':
        if self.atttype: Type=nctypes.type_2pycdf(self.atttype)
        else:            Type=nctypes.type_var2pycdf(self.value)
        a=self._nc.attr(self.name)
        a.put(Type,self.value)

      elif self._interface in ('scientific','netcdf4'):
        setattr(self._nc,self.name,value)



class Pyncvar(Common):
  def __init__(self,parent,varname):
    self._interface = parent._interface
    self.varname=varname
    self.parent=parent
    self.ncversion=self.parent.ncversion

    if  self._interface == 'pycdf':
      self._nc= parent._nc.var(self.varname)
    elif self._interface in ('scientific','netcdf4'):
      self._nc= parent._nc.variables[self.varname]

    Common.__init__(self)
    self.dims=self.__dims()
    self.dimnames=self.dims.keys()


  def __setitem__(self,elem,data):
    if  self._interface in ('pycdf','netcdf4'):
      self._nc.__setitem__(elem,data)
    elif self._interface == 'scientific':
      cmd='self._nc['+slice2str(elem)+']=data'
      exec(cmd)


  def __getitem__(self,elem):
    if  self._interface in ('pycdf','netcdf4'):
      # In some pycdf is required rebuild elem slice.
      # In such case, uncomment next two lines:
      # if isinstance(elem,slice) and elem.stop==maxint:
      #   elem=slice(elem.start,None,elem.step)

      data=self._nc.__getitem__(elem)

    elif self._interface == 'scientific':
      cmd='data=self._nc['+slice2str(elem)+']'
      try:    exec(cmd)
      except: data=self._nc.getValue()
      data=np.array(data)

    return data


  def dtype(self):
    '''Returns variable numpy dtype'''
    if self._interface=='netcdf4':
      return self._nc.dtype
    else:
      return nctypes.type_nc2numpy(self.type())
      # isstr and isbool are not important here !!


  def nctype(self):
    '''Returns variable NetCDF type'''
    if self._interface == 'pycdf':
      return nctypes.pycdf2nc(self._nc.inq_type())
    elif self._interface == 'scientific':
      return nctypes.type_numeric2nc(self._nc.typecode(),ncver=3)
      # should use self.ncversion but interface 2 only supports netcdf3!
    elif self._interface == 'netcdf4':
      return nctypes.type_numpy2nc(self._nc.dtype,ncver=self.ncversion)


  def range(self):
    '''Return variable range=(min, max)'''
    if self._interface in ('pycdf','netcdf4'):
      return self._nc[:].min(),self._nc[:].max()
    elif self._interface == 'scientific':
      return  np.min(self._nc[:]),np.max(self._nc[:])


  def shape(self):
    '''Returns variable shape tuple'''
    if self._interface in ('pycdf','netcdf4'):
      return self._nc.shape
    elif self._interface == 'scientific':
      return self._nc[:].shape


  def __dims(self):
    '''Returns variable dimensions dict'''
    d=cbt.odict()
    if   self._interface=='pycdf': names = self._nc.dimensions()
    elif self._interface=='scientific': names = self._nc.dimensions
    elif self._interface=='netcdf4': names = self._nc.dimensions

    if self._interface in ('pycdf','scientific'):
      for n in names: d[n]=self.parent.dims[n]

    # In netcdf 4, dimensions are scoped such that they can be seen
    # in all descendant groups. That is, dimensions can be shared
    # between variables in different groups, if they are defined in
    # a parent group.

    elif self._interface == 'netcdf4':
      for n in names:
        p=self.parent
        while True:
          if p.dims.has_key(n):
            d[n]=p.dims[n]
            break
          if p.parent: p=p.parent
          else: break

    return d


  def ndim(self):
    '''Returns variable number of dimensions'''
    return len(self.dims)


class Pyncgroup(Common):
  def __init__(self,nc,parent=False,name='root',info=False):

    self._nc=nc
    self.parent    = parent
    self.groupname = name

    # setting root:
    if name=='root':
      self.root=False
    else:
      if parent.root: self.root=parent.root
      else: self.root=parent

    # setting root info:
    if self.root:
      self._interface   = self.root._interface
      self.filename     = self.root.filename
      self.ncversion    = self.root.ncversion
      if not self.parent._ncdump_info is False:
        self._ncdump_info = self.parent._ncdump_info['groups'][self.groupname]
      else: self._ncdump_info=False

    else:
      self._interface   = info['interface']
      self.filename     = info['filename']
      self.ncversion    = info['version']
      self._ncdump_info = info['ncdump_info']


    # init Common:
    Common.__init__(self)

    # setting structure:
    if not self.root: self.structure=[self.groupname]
    else: self.structure=self.parent.structure+[self.groupname]

    # setting dims and dimnames:
    self.dims=self.__dims()
    self.dimnames=self.dims.keys()

    # setting groups and goupnames:
    self.groups = self.__groups()
    self.groupnames=self.groups.keys()

    # setting vars and varnames:
    self.vars=self.__vars()
    self.varnames=self.vars.keys()


  def add_dim(self,dimname,dimvalue):
    '''Creates file dimensions'''
    if self._interface=='pycdf':
      self._nc.def_dim(dimname,dimvalue)
    elif self._interface=='scientific':
      self._nc.createDimension(dimname,dimvalue)
    elif self._interface=='netcdf4':
      if dimvalue==0: dimvalue=None
      self._nc.createDimension(dimname,dimvalue)

    # update self.dims and self.dimnames:
    self.dims[dimname]=dimvalue

    # update ncdump_info:
    if not self._ncdump_info is False:
      self._ncdump_info['dimensions'][dimname]=str(dimvalue)


  def __dims(self):
    '''Returns file dimensions dict'''
    dims=cbt.odict()
    if self._interface=='pycdf':
      for i in range(self._nc.inq_ndims()):
        dims[self._nc.dim(i).inq()[0]]=self._nc.dim(i).inq()[1]
    elif self._interface=='scientific':
      # get ordered dims:
      if not self._ncdump_info is False:
        odims=self._ncdump_info['dimensions']
      else: odims=False

      if not odims: odims=self._nc.dimensions.keys()
      for k in odims: dims[k]=self._nc.dimensions[k]
    elif self._interface=='netcdf4':
      # get ordered dims:
      if not self._ncdump_info is False:
        odims=self._ncdump_info['dimensions']
      else: odims=False

      if not odims: odims=self._nc.dimensions.keys()
      for k in odims: dims[k]=len(self._nc.dimensions[k])

    return dims


  def dim_rename(self,oldname,newname):
    '''Renames file dimension'''
    if self._interface=='pycdf':
      self._nc.dim(oldname).rename(newname)
    elif self._interface=='scientific':
      print ':: dim_rename not implemented in ',self._interface
    elif self._interface=='netcdf4':
     self._nc.renameDimension(oldname,newname)

    # update self.dims:
    self.dims.rename_key(oldname,newname)

    # update ncdump info:
    if not self._ncdump_info is False:
      self._ncdump_info['dimensions'].rename_key(oldname,newname)


  def add_var(self,varname,vartype,dimnames,**kargs):
    '''
    Creates netcdf variable

    Ex: add_var('temp','float32',('lon','lat','z'))


    About compression:

    Compression kargs options (available for netcdf4 interface)
      zlib, default True, turn on compression
      lsd (least_significant_digit), default is None: lossy
        compression with lsd precision
      complevel, 1..9, copression level, default 4


    About vartype:

       i) netcdf type names:
         for netcdf3:
           byte,char,short,int,float,double

         for netcdf4:
           + int64, ushort, uint,uint64, string

       ii) Numeric type code: fd1silbwuc

       iii) Numpy dtype or type name

    kargs are required when vartype is a netcdf type name
    and the interface is scientific, see type_2numeric
    and when interface is netcdf4, see type_2dtype

    No need for kargs when:
     1. interface = pycdf
     2. interface is scientific and vartype is numpy type name or
        dtype, or vartipe is numeric typecode
     3. interface is netcdf4 and vartype is numpy type name or dtype
    '''

    if self._interface=='pycdf':
      vartype=nctypes.type_2pycdf(vartype)
      self._nc.def_var(varname,vartype,dimnames)
    elif self._interface=='scientific':
      vartype=type_2numeric(vartype,**kargs)
      self._nc.createVariable(varname,vartype,dimnames)
    elif self._interface=='netcdf4':
      vartype=nctypes.type_2numpy(vartype,**kargs)

      # NOTA: ocorre erro ao criar ex int64 em ncver 3. Pode fazer sentido
      # converter para int? ie, qd o tipo nao e suportado pela ncver!?

      zlib = kargs.get('zlib',True)
      lsd  = kargs.get('lsd',None)
      complevel  = kargs.get('complevel',4)
      fill_value = kargs.get('fill_value',None)

      self._nc.createVariable(varname,vartype,dimnames,zlib=zlib,
                              complevel=complevel,least_significant_digit=lsd,
                              fill_value=fill_value)

      newvar=Pyncvar(self,varname)

      # update self.vars and self.varnames:
      self.vars[varname]=newvar

      # update ncdump_info:
      if not self._ncdump_info is False:
        self._ncdump_info['variables'][varname]=cbt.odict()

    return newvar


  def __vars(self):
    '''Returns file variables'''
    var=cbt.odict()
    if self._interface=='pycdf':
      names=[self._nc.var(i).inq_name() for i in range(self._nc.inq_nvars())]
    elif self._interface in ('scientific','netcdf4'):
      # get ordered names:
      if not self._ncdump_info is False:
        onames=self._ncdump_info['variables'].keys()
      else: onames=False

      if onames: names=onames
      else: names=self._nc.variables.keys()

    for vname in names: var[vname]=Pyncvar(self,vname)
    return var


  def var_rename(self,oldname,newname):
    '''Renames variable'''
    if self._interface=='pycdf':
      self._nc.var(oldname).rename(newname)
    elif self._interface=='scientific':
      print ':: var_rename not implemented in ',self._interface
      return
    elif self._interface=='netcdf4':
      self._nc.renameVariable(oldname,newname)

    # update self.vars:
    self.vars.rename_key(oldname,newname)

    # update ncdump info:
    if not self._ncdump_info is False:
      self._ncdump_info['variables'].rename_key(oldname,newname)


  def add_group(self,groupname):
    '''
    Creates netcdf groups

    See:
    http://netcdf4-python.googlecode.com/svn/trunk/docs/netCDF4-module.html
    '''
    if self._interface=='netcdf4':
      self._nc.createGroup(groupname)
      newgroup=Pyncgroup(groupname,self)

      # update self.groups:
      self.groups[groupname]=newgroup

      # update ncdump info:
      if not self._ncdump_info is False:
        self._ncdump_info['groups'][groupname]=cbt.odict()

      return newgroup

    else:
      print ':: add_group only implemented in interface ','netcdf4'
      return False


  def __groups(self):
    out=cbt.odict()
    if self._interface=='netcdf4':
      gs=self._nc.groups
      try: gnames=gs.keys()
      except: gnames=()
      for k in gnames: out[k]=Pyncgroup(nc=gs[k],name=k,parent=self)
    else:
      print ':: groups only implemented in interface ','netcdf4'

    return out



class Pync(Pyncgroup):
  def __init__(self,filename,perm='r',interface=Interface,version=3,ncdump=False):

    if interface!='netcdf4': version=3

    if perm=='t': # trunc
      if os.path.isfile(filename): os.remove(filename)
      perm='w'

    if interface=='pycdf':
      if perm=='w':
        perm = pnc.NC.WRITE|pnc.NC.CREATE
      elif perm=='a':
        perm = pnc.NC.WRITE | pnc.NC.CREATE
      else:
        perm = pnc.NC.NOWRITE

      nc=pnc.CDF(filename,perm)
      nc.automode()

    elif interface=='scientific':
      nc=pnc.NetCDFFile(filename,perm)

    elif interface=='netcdf4':
      if perm=='w' and os.path.isfile(filename): perm='a'
      if not isinstance(version,basestring):
        format='NETCDF'+str(version)
        if format=='NETCDF3': format='NETCDF3_CLASSIC'
      else: format=version
      if not isinstance(filename,basestring) or filename.find('*')!=-1 or filename.find('?')!=-1:
        nc=pnc.MFDataset(filename)
      else:
        if os.path.isfile(filename):
          nc=pnc.Dataset(filename,perm)
        else:
          nc=pnc.Dataset(filename,perm,format=format)

    self.perm=perm

    # root info:
    info={}
    info['interface'] = interface
    info['filename']  = filename
    info['version']   = version
    if ncdump:
      from ncdump import ncdump_info
      info['ncdump_info']   = ncdump_info(filename)
    else: info['ncdump_info']=False

    # init group:
    Pyncgroup.__init__(self,nc,info=info)


  def sync(self):
    if self._interface=='netcdf4': self._nc.sync()
    else: print 'sync only in netcdf4'

  def close(self):
    self._nc.close()
