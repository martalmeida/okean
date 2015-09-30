# encoding: utf-8
''' 
OKEAN ocean modelling and analysis tools...

Based on the python/numpy/matplotlib scientific python suite. The toolkit
contains general modeling tools. Specific tools are also included for ROMS.
'''


__authors__ = 'Martinho Marta-Almeida <m.martalmeida@gmail.com> \
<Couto de Esteves, 3740-037, Portugal>',
__version__='2015-09-30 07:11:31.514064'


def doc(action='web',tag=''):
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
    name=name.split('_')[1:].join('_')
    F['name']=f

  if action=='local':
    if not F:
      print ':: no doc files found'
      print ':: (folder=%s)'%pdoc
      print ':: visit %s'%urldoc
    else:
      print ':: available doc files:'
      for k in F.keys():
        print '   - %-10s %s'%(k,F[k])

  elif action=='web':
    import webbrowser
    pdoc=os.path.join(__path__[0],'documentation')
    if tag=='':
      dest=urldoc
    elif tag=='index':
      dest=urldoc0+'okean_documentation.ipynb'
    else:
      dest=urldoc+'okean_%s.ipynb'%tag

    webbrowser.open(dest)



    
