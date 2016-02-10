# encoding: utf-8
''' 
OKEAN ocean modelling and analysis tools...

Based on the python/numpy/matplotlib scientific python suite. The toolkit
contains general modeling tools. Specific tools are also included for ROMS.
'''


__authors__ = 'Martinho Marta-Almeida <m.martalmeida@gmail.com> \
<Couto de Esteves, 3740-037, Portugal>',
__version__='2015-09-30 07:11:31.514064'


def doc(tag=None):
  '''
  Shows local doc files and opens ipynb documentation on github

  tag can be:
    'index' will open web doc index file
    ''      will open web doc main folder
    None    will show local files and print some help info

  ex: okean.doc('roms_glider')
  '''

  url='https://github.com/martalmeida/okean'
  urldoc=url+'/tree/master/okean/documentation/'
  urldoc0=url+'/blob/master/'

  import os
  import glob
  pdoc=os.path.join(__path__[0],'documentation')
  fdoc=os.path.join(pdoc,'*.ipynb')
  files=glob.glob(fdoc)
  F={}
  for f in files:
    name=os.path.splitext(os.path.basename(f))[0]
    name='_'.join(name.split('_')[1:])
    F[name]=f
    
  if tag==None:
    if not F:
      print ':: no local doc files found'
      print ':: folder %s'%pdoc
      print ':: visit %s'%urldoc
    else:
      print ':: available local doc files:'
      print ':: folder %s'%pdoc
      print ''
      print '     %-10s  %s'%('tag','filename')
      print '     ----------  ----------'
      for k in F.keys():
        print '   - %-10s %s'%(k,os.path.basename(F[k]))
      print ''
      print ':: doc files are also available in %s'%urldoc
      print ':: do okean(tag) to open web doc file with desired tag'
      print "    - ex: okean.doc('roms_glider')"
      print "    - ex: okean.doc('index') will open doc index file"
      print "    - ex: okean.doc('')      will open doc folder"

  else:
    import webbrowser
    pdoc=os.path.join(__path__[0],'documentation')
    if tag=='':
      dest=urldoc
    elif tag=='index':
      dest=urldoc0+'okean_documentation.ipynb'
    else:
      dest=urldoc+'okean_%s.ipynb'%tag

    webbrowser.open(dest)
