'''
python 2x, 3x
'''

from collections import OrderedDict
from .. import cookbook as cb

def parse_dims(lines):
  i=-1
  i0=-1
  i1=-1
  while i<len(lines)-1:
    i+=1
    if lines[i].lstrip().find('dimensions:')==0: i0=i+1
    elif (lines[i].lstrip().find('variables:')==0 or\
           lines[i].lstrip().find('group:')==0 or\
           lines[i].lstrip().find('} // group ')==0) and i0>-1:
      i1=i
      break

  res=OrderedDict()
  for j in range(i0,i1):
      tmp=lines[j].split('=')
      dimname=tmp[0].strip()
      dimvalue=tmp[1][:-1].strip()
      res[dimname]=dimvalue

  return res


def parse_vars(lines):
  i=-1
  i0=-1
  i1=-1
  while i<len(lines)-1:
    i+=1
    if lines[i].lstrip().find('variables:')==0:
      i0=i+1
    elif (lines[i].lstrip().find('// global attributes:')==0 or\
         lines[i].lstrip().find('group:')==0 or\
         lines[i].lstrip().find('} // group ')==0 or\
         lines[i].strip()=='}') and i0>-1:
      i1=i
      break

  res=OrderedDict()
  for j in range(i0,i1):
    if lines[j].strip():
        if lines[j].find('=')==-1 and lines[j].strip()[0]!='"':
          varname=lines[j][:-2].strip().split()[1].split('(')[0]
          res[varname]=OrderedDict()
        else:
          if lines[j].strip()[0]!='"': # multiline !!
            name=lines[j].split(':')[0].strip()
            vatt=lines[j].split(':')[1].split('=')[0].strip()
            res[name][vatt]={}

  return res


def parse_atts(lines):
  i=-1
  i0=-1
  i1=-1
  while i<len(lines)-1:
    i+=1
    if lines[i].lstrip().find('// global attributes:')==0: i0=i+1
    elif (lines[i].lstrip().find('// global attributes:')==0 or\
         lines[i].lstrip().find('group:')==0 or\
         lines[i].lstrip().find('} // group ')==0) and i0>-1:
      i1=i
      break

  res=OrderedDict()
  for j in range(i0,i1):
    if lines[j].strip() and lines[j].strip()[0]==':':
      s=lines[j]
      res[s[s.index(':')+1:s.index('=')].strip()]={}

  return res


def parse_once(lines,root=False):
  name=lines[0].split()[1]

  res=OrderedDict()
  res['groups'] = OrderedDict()
  res['name']   = name
  res['dimensions'] = parse_dims(lines)
  res['variables']  = parse_vars(lines)
  res['attributes']  = parse_atts(lines)

  j0=False
  k=0
  for j in range(1,len(lines)):
    k+=1
    l=lines[j]
    if l.lstrip().find('group:')==0 and j0 is False:
      j0=j
      name=l.split()[1]

    if l.strip() == '} // group '+name  and not j0 is False:
      res['groups'][name]=parse_once(lines[j0:j+1])
      j0=False

  return res



def ncdump_info(f):
  '''
  List of dimensions, variables and attributes from ncdump
  Works with netcdf files version <= 4
  '''
  ncdump=cb.search('ncdump')
  if ncdump:
    #try :   out=cb.run(ncdump+' -h '+f+' 2>/dev/null')
    #except: out=cb.run(ncdump+' -h '+f)
    out=cb.run(ncdump,'-h',f)
  else: return

  if not out:
    res=OrderedDict()
    res['groups'] = OrderedDict()
    res['name']       = f
    res['dimensions'] = OrderedDict()
    res['variables']  = OrderedDict()
    res['attributes'] = OrderedDict()
  else:
    res=parse_once(out,True)

  return res


def ncdump_info3(f):
  '''
  List of dimensions, variables and attributes from ncdump
  Works with netcdf files version < 4
  '''

  ncdump=cb.search('ncdump')
  if ncdump:
    #try :   out=cb.run(ncdump+' -h '+f+' 2>/dev/null')
    #except: out=cb.run(ncdump+' -h '+f)
    out=cb.run(ncdump,'-h',f)
  else:
    return

  dims=[]
  vars=[]
  varatts={}
  atts=[]

  i1=i2=i3=-1;
  hasdims=hasvars=hasatts=False

  if len(out)>1:

    i1=i2=i3=-1;

    for n in range(len(out)):
      if out[n].find('dimensions')==0: i1=n+1
      if out[n].find('variables')==0:  i2=n+1
      if out[n].find('// global attributes')==0:  i3=n+1

    if i1!=-1: hasdims=True
    if hasdims:
      if i2==-1: i2=i3-1
      if i3==-1: i2=len(out)-1

    if i2!=-1: hasvars=True
    if hasvars:
      if i3==-1: i3=len(out)-1

    if i3!=-1: hasatts=True

    if hasdims:
      for i in range(i1,i2-1):
        dims+=[out[i].split('=')[0].strip()]

    if hasvars:
      for i in range(i2,i3-2):
        if out[i].find('=')==-1:
          vars+=[out[i][:-2].strip().split()[1].split('(')[0]]
          varatts[vars[-1]]=[]
        else:
          name=out[i].split(':')[0].strip()
          vatt=out[i].split(':')[1].split('=')[0].strip()
          varatts[name]+=[vatt]

    if hasatts:
      for i in range(i3,len(out)-1):
        s=out[i]
        if s.strip()[0]==':':
          atts+=[s[s.index(':')+1:s.index('=')].strip()]

  return {'dimensions':dims,'variables':vars,'attributes':atts,'vattributes':varatts}

