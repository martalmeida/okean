import sys
import os
import site
import shutil
import glob
import subprocess


def go():
  here=os.path.dirname(__file__)

  dest=os.path.join(site.getsitepackages()[0],'okean')
  if os.path.isdir(dest):
    print('## folder exists: %s'%dest)
    shutil.rmtree(dest)

  os.makedirs(dest)

  # cp py files from main folder:
  src=os.path.join(here,'okean')
  files=glob.glob(os.path.join(src,'*.py'))
  for fsrc in files:
    print('  -- copying %-20s to %s'%(fsrc,dest))
    shutil.copy2(fsrc,dest)

  # copy contents of bin folder to bin:
  src=os.path.join(here,'okean','bin')
  dest_=os.path.dirname(sys.executable)
  files=glob.glob(os.path.join(src,'*'))
  for f in files:
    print('## copying %s to %s'%(f,dest_))
    shutil.copy2(f,dest_)

    # replace 1st line
    ff=os.path.join(dest_,os.path.basename(f))

    lines=open(ff).readlines()
    if lines[0].strip().startswith('#!'):
      print('   - fixing %s'%ff)
      lines[0]='#!'+sys.executable+'\n'
      open(ff,'w').writelines(lines)

  # copy many folders:
  dirs='data','datasets','misc','nc','roms','util','ext'
  for d in dirs:
    src=os.path.join(here,'okean',d)
    dest_=os.path.join(dest,d)
    print('## copying %s to %s'%(src,dest))
    shutil.copytree(src,dest_)

    # remove folders with upercase characters:
    for root, dirs, files in os.walk(dest_):
      for d in dirs:
         if d!=d.lower() or d=='__pycache__':# or d=='ext':
           D=os.path.join(root,d)
           print('   - removing %s'%D)
           shutil.rmtree(D)

  # extensions:
  # why use subprocess instead of f2py.compile? See:
  # https://numpy.org/doc/2.1/f2py/usage.html

  extDir=[] # folders to remove after compilation

  pwd=os.getcwd()
  soDir=dest
  extDir+=[os.path.join(soDir,'ext')]
  os.chdir(soDir)


  s=sys.executable,'-m','numpy.f2py','-c','ext/alg.f','-m','alg'
  subprocess.run(s)
  s=sys.executable,'-m','numpy.f2py','-c','ext/pnpoly.f','-m','pnpoly'
  subprocess.run(s)
  s=sys.executable,'-m','numpy.f2py','-c','ext/lu.f90','-m','lusolver'
  subprocess.run(s)

  soDir=os.path.join(dest,'roms')
  extDir+=[os.path.join(soDir,'ext')]
  os.chdir(soDir)

  s=sys.executable,'-m','numpy.f2py','-c','ext/rtools.f90','../ext/pppack.f90','-m','rtools'
  subprocess.run(s)

  # remove ext folders no longer needed:
  for e in extDir:
    print('## removing %s'%e)
    shutil.rmtree(e)

  os.chdir(pwd)

  # check if extensions created:
  print('\n## checking extensions')
  status=0
  for s in 'alg','lusolver','pnpoly',os.path.join('roms','rtools'):
    files=glob.glob(os.path.join(dest,s+'*.so'))
    if len(files)==1:
      f=files[0]#.replace(dest,'')
      print('   - file %s created'%f)
    else:
      print('#ERROR: %s*.so not created'%s)
      status=1

  if not status: print('\n   enjoy okean')



if __name__=='__main__':
  go()
