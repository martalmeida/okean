"""Okean ocean modelling and analysis tools ...

Requires:
    NumPy
    matplotlib with the Basemap toolkit
    netcdf interface for python (like netCDF4)

"""

classifiers = """\
Development Status :: 3 - Alpha
Environment :: Console
Intended Audience :: Science/Research
Intended Audience :: Developers
License :: OSI Approved :: European Union Public Licence 1.2 (EUPL 1.2)
Operating System :: OS Independent
Programming Language :: Python
Topic :: Scientific/Engineering
Topic :: Software Development :: Libraries :: Python Modules
"""

if 1:
  from numpy.distutils.core import Extension
  from numpy.distutils.command.install import install
  from numpy.distutils.core import setup
else:
  from setuptools import setup, Extension
  from setuptools.command.install import install

class my_install(install):
    def run(self):
        install.run(self)

        # post installation:


        if 0:
          # link pppack
          import os
          from distutils.sysconfig import get_python_lib
          p=get_python_lib(plat_specific=1) # installation folder
          src=os.path.join(p,'okean','roms','rtools.so')
          dest=os.path.join(p,'okean','pppack.so')
          #os.symlink(src,dest)

        print('''
        enjoy okean
        ''')


alg = Extension(name = 'alg',
                sources = ['okean/ext/alg.f'])

pnpoly = Extension(name = 'pnpoly',
                sources=['okean/ext/pnpoly.f'])

#rtools = Extension(name = 'roms.rtools',
#                sources=['okean/roms/ext/rtools.f90'])
rtools = Extension(name = 'roms.rtools',
                sources=['okean/roms/ext/rtools.f90','okean/ext/pppack.f90'])

pppack = Extension(name = 'pppack',
                sources=['okean/ext/pppack.f90'])

lu = Extension(name = 'lusolver',
                sources=['okean/ext/lu.f90'])

#rtools22 = Extension(name = 'roms.rtools_vs2vt2',
#                sources=['okean/roms/ext/rtools_vs2vt2.f'])
#
#rtools42 = Extension(name = 'roms.rtools_vs4vt2',
#                sources=['okean/roms/ext/rtools_vs4vt2.f'])
#
#rtools12 = Extension(name = 'roms.rtools_vs1vt2',
#                sources=['okean/roms/ext/rtools_vs1vt2.f'])

import glob
ncview_cm=glob.glob('okean/data/ncview_cmaps/*')
rgui_icons=glob.glob('okean/roms/gui/icons/*')
okean_doc=glob.glob('okean/documentation/*')

doclines = __doc__.split("\n")

def get_version():
  v='unknonw'
  lines=open('okean/__init__.py').readlines()
  for l in  lines:
    if l.startswith('__version__'): v=eval(l.split('=')[-1])

  return v
  

if __name__ == '__main__':

    setup(name = "okean",
          version = get_version(),
          description = doclines[0],
          long_description = "\n".join(doclines[2:]),
          author = "Martinho Marta-Almeida",
          author_email = "m.martalmeida@gmail.com",
          url = "https://github.com/martalmeida/okean",
          packages = ['okean',
                      'okean.roms',
                      'okean.roms.inputs',
                      'okean.roms.gui',
                      'okean.nc',
                      'okean.util',
                      'okean.datasets'],
#                      'okean.data'],
          license = 'EUPL',
          platforms = ["any"],
          ext_package='okean',
          ext_modules = [alg,pnpoly,rtools,pppack,lu],
          data_files = [('okean/roms/gui', ['okean/roms/gui/romsgui.derived',
                                            'okean/roms/gui/rgui.png']),
                         ('okean/roms/gui/icons', rgui_icons),
                         ('okean/data', ['okean/data/cities_world.txt',
                                         'okean/data/cities_more.txt']),
                         ('okean/misc', ['okean/misc/hull_code.tar.gz']),
                         ('okean/data/ncview_cmaps/', ncview_cm),
                         ('okean/documentation/',okean_doc),
#                         ('',['okean_documentation.ipynb']),
                         ('',['EUPL v1_2 PT.txt'])],
          classifiers = list(filter(None, classifiers.split("\n"))),
          scripts=['okean/bin/rgui','okean/bin/show_nctime','okean/bin/show','okean/bin/qstate',
                   'okean/bin/romsview','okean/bin/disp'],
          cmdclass={'install': my_install},
          )

