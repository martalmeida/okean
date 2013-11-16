"""Okean ocean modelling and analysis tools ...

Requires:
    NumPy
    matplotlib with the Basemap toolkit
    netcdf interface for python (like netCDF4)

"""

classifiers = """\
Development Status :: alpha
Environment :: Console
Intended Audience :: Science/Research
Intended Audience :: Developers
License :: European Union Public Licence - EUPL v.1.1
Operating System :: OS Independent
Programming Language :: Python
Topic :: Scientific/Engineering
Topic :: Software Development :: Libraries :: Python Modules
"""

from numpy.distutils.core import Extension
from numpy.distutils.command.install import install


class my_install(install):
    def run(self):
        install.run(self)

        # post installation:

        #installation folder:
        #from distutils.sysconfig import get_python_lib
        #get_python_lib(plat_specific=1)

        print '''
        enjoy okean
        '''
        '''
        ---------------------------- INFO -----------------------------
        The usage of some gshhs tools needs the env variable GSHHS_DATA
        gshhs data can be downloaded from
        http://www.ngdc.noaa.gov/mgg/shorelines/data/gshhs/

        Check folder misc for some utilities...


        enjoy okean
        ---------------------------------------------------------------
        '''


alg = Extension(name = 'alg',
                sources = ['okean/ext/alg.f'])

pnpoly = Extension(name = 'pnpoly',
                sources=['okean/ext/pnpoly.f'])

rtools = Extension(name = 'roms.rtools',
                sources=['okean/roms/ext/rtools.f'])

rtools22 = Extension(name = 'roms.rtools_vs2vt2',
                sources=['okean/roms/ext/rtools_vs2vt2.f'])

rtools42 = Extension(name = 'roms.rtools_vs4vt2',
                sources=['okean/roms/ext/rtools_vs4vt2.f'])


doclines = __doc__.split("\n")

def get_version():
  l=open('okean/version.py').readlines()
  s0=l[0].split('=')[1].strip()
  s1=l[1].split('=')[1].strip()
  return s0+'.'+s1

if __name__ == '__main__':
    from numpy.distutils.core import setup

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
                      'okean.datasets',
                      'okean.data'],
          license = 'EUPL',
          platforms = ["any"],
          ext_package='okean',
          ext_modules = [alg,pnpoly,rtools,rtools22,rtools42],
          data_files = [('okean/roms/gui', [#'okean/roms/gui/cities_world.txt',
                                            #'okean/roms/gui/cities_more.txt',
                                            'okean/roms/gui/romsgui.derived',
                                            'okean/roms/gui/rgui.png']),
                         ('okean/data', ['okean/data/cities_world.txt',
                                   'okean/data/cities_more.txt']),
                         ('okean/misc', ['okean/misc/hull_code.tar.gz']),
                         ('',['EUPL v.1.1 - licencia.pdf'])],
          classifiers = filter(None, classifiers.split("\n")),
          scripts=['okean/bin/rgui','okean/bin/show_nctime','okean/bin/show','okean/bin/qstate'],
          cmdclass={'install': my_install},
          )

