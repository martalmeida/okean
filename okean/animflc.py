"""creation and manipulation of FLI/FLC animations
MMA 11-2007, IEO La Corunha
"""

import os
import Image
from tempfile import mkstemp
from cookbook import unique

def run(s):
  out=os.popen(s)
  out=out.read()
  return out.split('\n')[:-1]

def ppm2fli(ppmfiles='*.ppm',animfile='anim.flc',fps=4,size='',options='',percent=1,convert=True,fig0=False,clean=False):
  """if convert, convert is used in resize (percent not 1)

  """

  if isinstance(ppmfiles,basestring):
    ppmfiles=run('ls '+ppmfiles)
    if not ppmfiles: return

  # convert fig0:
  if fig0:
    im=Image.open(fig0)
    if im.format not in ('PPM','PGM','PBM','FBM'):
      if convert:
        cmd='convert %s %s' % (fig0,fig0+'.ppm')
        fig0=fig0+'.ppm'
        run(cmd)
      else:
        im.save(fig0+'.ppm')
        fig0=fig0+'.ppm'

  files=[]
  istmp={}
  toRemove=[]
  isize=False
  for file0 in ppmfiles:
    file=file0
    im=Image.open(file0)
    if im.format not in ('PPM','PGM','PBM','FBM'):
      file=mkstemp('.ppm')[1]
      toRemove+=[file]
      istmp[file]=True
      if not convert:
        im.save(file)
      else:
        cmd='convert %s %s' % (file0,file)
        run(cmd)

    files+=[file]

  if not percent==1:
    Files=[]
    for f in files:
      if not convert:
        im=Image.open(f)
        sz=im.size[0]*percent,im.size[1]*percent
        isize=sz
        if not istmp[f]:
          f=mkstemp('.ppm')[1]
          toRemove+=[f]

        im=Image.open(f)
        im=im.resize(sz,resample=1)
        im.save(f)

      else:
        f0=f
        if not istmp[f]:
          f=mkstemp('.ppm')[1]
          toRemove+=[f]

        cmd='convert -geometry %%%d %s %s' % (percent*100,f0,f)
        run(cmd)

      Files+=[f]

    files=Files

  listFile = mkstemp('.ppmlist')[1]
  f=open(listFile,'w')
  for file in files:
    f.write(file+'\n')

  f.close()
  toRemove+=[listFile]

  s='-s '+str(1000./fps)

  if not size:
    im=Image.open(ppmfiles[0])
    size='-g %dx%d' % im.size
  else:
    size='-g '+size

  if fig0:
    cmd='ppm2fli -m '+fig0+' '+s+' '+size+' -N '+listFile+' '+options+' '+animfile;
  else:
    cmd='ppm2fli  '+s+' '+size+' -N '+listFile+' '+options+' '+animfile;

  print(cmd)
  a=run(cmd)
  if a:
    print('\ncreated '+animfile+' from:')
    i=0
    for f in files:
      i+=1
      print('%3d %s' % (i,f))

  # remove tmp files:
  if clean:
    toRemove=unique(toRemove)
    for f in toRemove:
      print('removing ',f)
      os.remove(f)
