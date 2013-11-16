'''
 -------------------------------------------------------------------
 datatype  dtype.name name             nc4   description
 -------------------------------------------------------------------
 f4,f   float32      NC_FLOAT               32-bit floating point
 f8     float64      NC_DOUBLE              64-bit floating point

 i1     int8         NC_BYTE                 8-bit signed integer
 i2     int16        NC_SHORT               16-bit signed integer
 i4,i   int32        NC_INT or NC_LONG      32-bit signed integer
 i8     int64        NC_INT64           1   64-bit signed integer

 u1     uint8        NC_CHAR            2    8-bit unsigned integer
 u2     uint16       NC_USHORT          1   16-bit unsigned integer
 u4     uint32       NC_UINT            1   32-bit unsigned integer
 u8     uint64       NC_UINT64          1   64-bit unsigned integer

 S#     string8^#    NC_STRING          1   variable length character string

 b(i1)  bool         NC_BYTE
 --------------------------------------------------------------------

 1) Available only for netCDF-4 format files.
    All the unsigned ints (except NC_CHAR), the 64-bit ints, and
    string type are for netCDF-4 files only

 2) Char used in netcdf3 to represent strings!

 ------------------
 Numpy   Numeric
 ------------------
 f4      f
 f8      d

 i1      1
 i2      s
 i4      i
 i8      l

 u1      b
 u2      w
 u4      u
 u8      None
 S1      c

 b       1
 -------------------
'''

import numpy as np

nptypes={'float32':'NC_FLOAT','float64':'NC_DOUBLE','int8':'NC_BYTE',
         'int16':'NC_SHORT','int32':'NC_INT','int64':'NC_INT64',
         'uint8':'NC_CHAR','uint16':'NC_USHORT','uint32':'NC_UINT',
         'uint64':'NC_UINT64','stringN':'NC_STRING','bool':'NC_BYTE'}


np2numeric={'float32':'f','float64':'d','int8':'1','int16':'s',
            'int32':'i','int64':'l', 'uint8':'b','uint16':'w',
            'uint32':'u','uint64':False,'stringN':'c','bool':'1'}

pycdftypes=['byte','char','short','int','float','double']

numpynames=['float32','float64','int8','int16','int32','int64','uint8',
            'uint16','uint32','uint64','stringN','bool']

ncnames=['float','double','byte','short','int','long','int64','char',
         'ushort','uint','uint64','string']

def type_numpy2nc(type,ncver=4):
  '''
  Convert numpy typecode to netcdf type
  '''
  if isinstance(type,basestring):
    npname=np.dtype(type).name
  else: # is datype
    npname=type.name


  if ncver==3 and npname.lower().find('uint')==0: return 

  if npname.find('string')==0 or npname.find('unicode')==0:
    npname='stringN'
    if ncver==3: npname='uint8' # to return CHAR

  # down type case version 3:
  if ncver==3:
    if npname.lower().find('int')==0:
      sz=int(npname[3:])
      if sz>32: npname='int32'


  return nptypes[npname][3:]


def type_nc2numpy(type,strlen=1,isstr=False,isbool=False):
  '''
  Convert netcdf type to numpy dtype
  '''
  type=type.upper()
  if type=='STRING': return np.dtype('S'+str(strlen))

  if type=='CHAR' and isstr: return np.dtype('S1') # nc version 3

  if type=='BYTE':
    if isbool: return np.dtype('bool')
    else: return np.dtype('int8')

  for k in nptypes.keys():
    if nptypes[k][3:]==type: return np.dtype(k)


def type_numpy2numeric(type):
  '''
  Convert numpy dtype to numeric typecode
  '''
  if isinstance(type,basestring):
    npname=np.dtype(type).name
  else: # is datype
    npname=type.name

  if npname.find('string')==0:
    npname='stringN'

  return np2numeric[npname]


def type_numeric2numpy(type,strlen=1):
  '''
  Convert numeric typecode to numpy dtype
  '''
  if type=='c':
    return np.dtype('S'+str(strlen))

  for k in np2numeric.keys():
    if np2numeric[k]==type: return np.dtype(k)


def type_numeric2nc(type,ncver=4):
  '''
  Convert numeric typecode to netcdf type
  '''
  # strlen is not important here:
  tmp=type_numeric2numpy(type,strlen=1)
  return type_numpy2nc(tmp,ncver=ncver)


def type_nc2numeric(type,isstr=False):
  '''
  Convert netcdf type to numeric typecode

  isstr (default is False) required since netcdf 3 CHAR can be
  integers8 or strings
  '''
  # strlen and isbool not important
  tmp=type_nc2numpy(type,isstr=isstr)
  return type_numpy2numeric(tmp)


def type_nc2pycdf(type):
  '''
  Convert netcdf type to pycdf type number

  byte   pycdf.NC.BYTE   = 1
  char   pycdf.NC.CHAR   = 2
  short  pycdf.NC.SHORT  = 3
  int    pycdf.NC.INT    = 4
  float  pycdf.NC.FLOAT  = 5
  double pycdf.NC.DOUBLE = 6
  '''
  if pycdftypes.count(type.lower())==1:
    return pycdftypes.index(type.lower())+1

def type_pycdf2nc(num):
  '''Convert pycdf type numebr to netcdf type '''
  return pycdftypes[num-1]


def type_var2numpy(v):
  '''
  Numpy dtype from python data values
  Ex: type_var2numpy([1,2,3], type_var2numpy('a')
  '''
  return np.array(v).dtype


def type_var2numeric(v):
  '''
  Numeric typecode from python data values
  Ex: type_var2numeric([1,2,3], type_var2numeric('a')
  '''
  return type_numpy2numeric(type_var2numpy(v))


def type_var2nc(v,ncver=4):
  '''Netcdf type from python data values'''
  type=type_var2numpy(v)
  return type_numpy2nc(type,ncver=ncver)


def type_var2pycdf(v):
  '''Pycdf type code from python data values'''
  ncver=3
  type=type_var2nc(v,ncver=ncver)
  return type_nc2pycdf(type)


def type_2dtype(type,**kargs):
  '''Convert Numeric typecode or netcdf type to numpy dtype
     also supports numpy type names

  kargs:
    strlen: 1,when converting from Numeric character typecode.
    isstr:  False, netcdf type CHAR may be used as numpy S1 in
            netcdf 3.
    isbool: False, when converting nc to numpy, the type BYTE may be
            seen as numpy boolean.
  '''
  strlen=1
  isstr=False
  isbool=False
  if 'strlen' in kargs.keys(): strlen = kargs['strlen']
  if 'isstr'  in kargs.keys(): isstr  = kargs['isstr']
  if 'isbool' in kargs.keys(): isbool = kargs['isbool']

  if isinstance(type,basestring):
    if len(type)==1: # is a numeric typecode
      return type_numeric2numpy(type,strlen=strlen)
    elif type in numpynames:  # numpy type name
      return np.dtype(type)
    else: # is netcdf type name:
      return type_nc2numpy(type,strlen=strlen,isstr=isstr,isbool=isbool)

  elif isinstance(type,np.dtype): return type
  else: return False


def type_2numpy(type,**kargs):
  '''
  Same as type_2dtype
  '''
  return type_2dtype(type,**kargs)


def type_2nc(type,**kargs):
  '''
  Convert Numeric typecode, numpy dtype and type name
  to netcdf type

  if inputs is already a netcdf type its is returned as is

  kargs:
    ncver: 4, when converting from numpy name or dtype
  '''

  ncver=4
  if 'ncver'  in kargs.keys(): ncver  = kargs['ncver']

  if isinstance(type,basestring):
    if len(type)==1: # is a numeric typecode
      return type_numeric2nc(type,ncver=ncver)
    elif type in numpynames:
      return type_numpy2nc(type,ncver=ncver)
    else: # is netcdf type name
      return type

  elif isinstance(type,np.dtype):
    return type_numpy2nc(type,ncver=ncver)
  else: return False


def type_2pycdf(type):
  '''
  Convert Numeric typecode,  numpy dtype and type name or netcdf type
  to pycdf code
  '''

  strlen=1
  ncver=3
  type=type_2nc(type,strlen=strlen,ncver=ncver)
  return type_nc2pycdf(type)


def type_2numeric(type,**kargs):
  '''
  Convert numpy dtype and type name and netcdf typename to numeric
  typecode

  kargs:
    isstr, False, used when converting from netcdf type name
  '''

  isstr=False
  if 'isstr'  in kargs.keys(): isstr  = kargs['isstr']

  if isinstance(type,basestring):
    if len(type)==1: # numeric typecode
      return type
    elif type in numpynames:
      return type_numpy2numeric(type)
    else: # netcdf type name
      return type_nc2numeric(type,isstr=isstr)

  elif isinstance(type,nc.dtype):
    return type_numpy2numeric(type)
  else: return False

