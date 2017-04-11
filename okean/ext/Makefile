# make FC=...
# for gfortran: make FC=gnu95

FC=
ifdef FC
  FComp=--fcompiler=$(FC)
else
  FComp=
endif

all: pnpoly alg

pnpoly:
	f2py -c pnpoly.f -m pnpoly $(FComp)

alg:
	f2py -c alg.f -m alg $(FComp)

check:
	f2py -c --help-fcompiler

lusolver:
	f2py -c lu.f90 -m lusolver $(FComp)

