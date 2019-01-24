import time
import numpy as np

as_consts={}
as_consts['rho_air']   = 1.22      # kg/m^3
as_consts['kappa']     = 0.4       # von Karman's constant
as_consts['emiss_lw']  = 0.985     # long-wave emissivity of ocean from Dickey et al
as_consts['sigmaSB']   = 5.6697e-8 # Stefan-Boltzmann constant [W/m^2/K^4]
as_consts['P_default'] = 1020      # default air pressure for Kinneret [mbars]


def qsat(Ta,Pa=as_consts['P_default']):
  '''
  Specific humidity at saturation

  q=QSAT(Ta) computes the specific humidity (kg/kg) at satuation at
  air temperature Ta (deg C). Dependence on air pressure, Pa, is small,
  but is included as an optional input.

  INPUT:   Ta - air temperature  [C]
           Pa - (optional) pressure [mb]

  OUTPUT:  q  - saturation specific humidity  [kg/kg]

  Version 1.0 used Tetens' formula for saturation vapor pressure_
  from Buck (1981), J. App. Meteor., 1527-1532.  This version_
  follows the saturation specific humidity computation in the COARE
  Fortran code v2.5b.  This results in an increase of ~5% in_
  latent heat flux compared to the calculation with version 1.0.


  from matlab air-sea tollbox at sea-mat
  '''

  # original code
  # a=(1.004.*6.112*0.6220)./Pa; 
  # q=a.*exp((17.502.*Ta)./(240.97+Ta))

  # as in Fortran code v2.5b for COARE
  ew = 6.1121*(1.0007+3.46e-6*Pa)*np.exp((17.502*Ta)/(240.97+Ta)) # in mb
  q  = 0.62197*(ew/(Pa-0.378*ew))                                 # mb -> kg/kg
  return q


def lwhf(Ts,dlw,dsw=False):
  '''
  LWHF: computes net longwave heat flux following Dickey et al (1994).
  qlw=LWHF(Ts,dlw) computes the net longwave heat flux into the ocean.
  Following Dickey et al (1994), J. Atmos. Oceanic Tech., 11, 1057-1078,
  the incident longwave flux can be corrected for sensor heating due to
  insolation if you are using Epply or Kipp & Zonen CG1 pyrgeometers.
  In this case, use qlw=LWHF(Ts,dlw,dsw). Epply is the default_
  pyrgeometer; change code for the Kipp & Zonen instrument.

    INPUT:  Ts  - sea surface temperature [K]
            dlw - (measured, positive) downward longwave flux [W/m^2]
            dsw - (measured) insolation [W/m^2] (needed for Eppley_
                  or Kipp & Zonen pyrgeometers)

    OUTPUT: qlw - net longwave heat flux [W/m^2]

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  3/8/97: version 1.0
  8/19/98: version 1.1 (revised for non-Epply pyrgeometers by RP)
  4/9/99: version 1.2 (included Kipp & Zonen CG1 pyrgeometers by AA)
  8/5/99: version 2.0
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  mma, IEO
  '''

  # correct dlw for sensor heating by insolation
  if dsw:
    # this line is for Epply pyrgeometers
    dlwc=dlw-0.036*dsw

    # this line is for Kipp & Zonen CG1 pyrgeometers
    # (the offset is specified as 25 W/m^2 at 1000 W/m^2)
    # dlwc=dlw-0.025.*dsw;
  else:
    dlwc=dlw

  # compute upward gray-body longwave flux
  lwup=-as_consts['emiss_lw']*as_consts['sigmaSB']*(Ts**4)

  # compute net flux
  qlw=lwup + as_consts['emiss_lw']*dlwc

  return qlw


def wind_stress(u,v,height=10.,rho_air=as_consts['rho_air']):
  '''Large and Pond (1981) u,v wind stress

     Inputs:
       u, v : East-West, North-West wind components (m/s)

     Outputs:
       taux, tauy: East-West, North-West wind stress components (Pa)
  '''

  speed=np.sqrt(u**2 + v**2)
  #Cd=(0.142 + 0.0764.*speed + 2.7./(speed+0.000001)).*0.001 # for height 10m
  Cd,sp10=cdnlp(speed,height)
  taux = rho_air*Cd*u*speed
  tauy = rho_air*Cd*v*speed

  return taux,tauy


def stresslp(speed,height,rho_air=as_consts['rho_air']):
  '''Large and Pond (1981) wind stress
  '''
  Cd,u10=cdnlp(speed,height)
  return rho_air*(Cd*u10**2)


def cdnlp(speed,height):
  '''Computes neutral drag coefficient following Large&Pond (1981),
     J. Phys. Oceanog., 11, 324-336

     Inputs:
       speed: wind speed (m/s)
       height: measurement height (m)


     Outputs:
       Cd, neutral drag coefficient at 10m
       speed10, wind speed at 10m (m/s)
  '''

  a=np.log(height/10.)/as_consts['kappa']  # log-layer correction factor
  tol=.001            # tolerance for iteration [m/s]

  u10o=np.zeros(speed.shape)
  cd=1.15e-3*np.ones(speed.shape)
  u10=speed/(1+a*np.sqrt(cd))

  ii=np.abs(u10-u10o)>tol;
  while np.any(ii):
    u10o=u10;
    cd=(4.9e-4+6.5e-5*u10o)   # compute cd(u10)
    # cd(u10o<10.15385)=1.15e-3;
    cd=np.where(u10o<10.15385,1.15e-3,cd)
    u10=speed/(1+a*np.sqrt(cd));   # next iteration
    ii=np.abs(u10-u10o)>tol      # keep going until iteration converges

  return cd,u10


def s2hms(secs):
  '''
  convert seconds to interger hour, minute, and seconds

  from air_sea toolbox (sea-mat)
  '''
  sec = np.round(secs)
  hr  = np.floor(sec/3600.)
  min = np.floor(np.remainder(sec,3600)/60)
  sec = np.round(np.remainder(sec,60))
  return hr,min,sec


def hms2h(h,m,s):
  '''
  convert hours, minutes, and seconds to hours

  from air_sea toolbox (sea-mat)
  '''

  return h+(m+s/60.)/60.;


def julianmd(y,m,d,h=0.,min=0.,sec=0.,startAtNoon=False):
  '''
  convert Gregorian calendar time to decimal Julian day
  converts Gregorian calendar dates to corresponding
  Julian day numbers.  Although the formal definition holds that Julian
  days start and end at noon, here Julian days start and end at midnight.
  In this convention, Julian day 2440000 began at 0000 UT, May 23, 1968.

  from air_sea toolbox (sea-mat)
  '''

  h=hms2h(h,min,sec)

  if not isinstance(y,np.ndarray):
    y=np.array([y])
    m=np.array([m])
    d=np.array([d])
    h=np.array([h])

  mo=m+9.
  yr=y-1.

  i=m>2
  mo[i]=m[i]-3
  yr[i]=y[i]
  c = np.floor(yr/100.)
  yr = yr - c*100.
  j = np.floor((146097.*c)/4.) + np.floor((1461.*yr)/4.) + \
      np.floor((153.*mo +2.)/5.) +d +1721119.

  if startAtNoon:
    j=j+(h-12)/24;
  else:
    j=j+h/24.

  return j


def greg2(yd,yr,startAtNoon=False):
  '''
  convert decimal yearday to standard Gregorian time
  In this convention, Julian day 2440000
  begins at 0000 UT May 23 1968.

  Returns year mo da hr mi sec

  from air_sea toolbox (sea-mat)
  '''
  js = julianmd(yr,1.,1.,startAtNoon=startAtNoon)
  julian = js + yd

  secs=np.remainder(julian,1)*24.*3600.
  j = np.floor(julian) - 1721119.
  In = 4*j -1;
  y = np.floor(In/146097.)
  j = In - 146097.*y
  In = np.floor(j/4.)
  In = 4.*In +3.
  j = np.floor(In/1461.)
  d = np.floor(((In - 1461.*j) +4.)/4.)
  In = 5.*d -3.
  m = np.floor(In/153.)
  d = np.floor(((In - 153.*m) +5.)/5.)
  y = y*100. +j
  mo=m-9.
  yr=y+1.
  i=m<10
  mo[i]=m[i]+3.
  yr[i]=y[i]
  hour,min,sec=s2hms(secs)

  return yr,mo,d,hour,min,sec


def yearday(y,m,d,h=0.,min=0.,sec=0.,y0=False):
  '''
  calender year, month and day, hour, min, sec into yearday
  y0 is the start/reference year
  yearday of yyyy-01-01 00h00 0secs is 0 (without y0)
  Ex:
   yearday(2005,1,1) --> 0
   yearday(2005,1,1,y0=2004) -> 366

  adapted (and enhanced) from sea_mat toolbox (sea-mat)
  '''

  if y0 is False: y0=y
  return julianmd(y,m,d,h,min,sec)-julianmd(y0,1,1)


#=====================================================================
# other tools, not in the air_sea original package:
#=====================================================================

C2K=273.15

def cloud_fraction(lw,Tsea,Tair,Rh,Wtype='net'):
  '''
  Cloud cover (fraction, 0..1)

  lw    - net longwave -- positive down -- or downward lw if Wtype='down'
  Tsea  - sst (C)
  Tair  - air temperature (C)
  Rh    - relative humudity (0..1)
  Wtype - lw type: 'down' or 'net' --> very low accuracy if is down

  If Wtype is down, net LW heat flux is computed  following Dickey et al (1994)

  Cloud cover is computed based on Berliand (1952) formula. The equation for
  saturation vapor pressure is from Gill (Atmosphere-Ocean Dynamics,
  pp 606). Here the coefficient in the cloud term is assumed constant,
  but it is a function of latitude varying from 1.0 at poles to 0.5 at
  the equator).

  mma, Texas A & M, aug 2011
  '''

  emmiss  = 0.97    # Infrared emmissivity
  StefBo  = 5.67e-8 # Stefan-Boltzmann constant (W/m2/K4)

  TairK=Tair+C2K
  TseaK=Tsea+C2K

  if Wtype is 'down':  LW=lwhf(TseaK,lw)
  elif Wtype is 'net': LW=lw

  cff  = (0.7859+0.03477*Tair)/(1.0+0.00412*Tair)
  sat  = 10.0**cff   # saturation vapor pressure (hPa or mbar)
  vap  = sat*Rh      # water vapor pressure (hPa or mbar)
  cff2 = TairK**3
  cff1 = TairK**4

  x=(-LW/(emmiss*StefBo) -cff2*4*(TseaK-TairK))/(cff1*(0.39-0.05*np.sqrt(vap)))
  try:
    x=np.min([x,1.])
    x=np.max([x,1.-0.6823])
  except: # numpy arrays:
    x[x>1.]=1.
    x[x<(1.-0.6823)]=1.-0.6823

  return np.sqrt((1-x)/0.6823)


def nlw(cloud,Tsea,Tair,Rh):
  '''
  Net longwave

  cloud (0..1)
  Tsea,Tair (C)
  Rh (0..1)

  Use Berliand (1952) formula to calculate net longwave radiation.
  The equation for saturation vapor pressure is from Gill (Atmosphere-
  Ocean Dynamics, pp 606). Here the coefficient in the cloud term
  is assumed constant, but it is a function of latitude varying from
  1.0 at poles to 0.5 at the equator).

  mma, Texas A & M, aug 2011
  '''

  emmiss  = 0.97    # Infrared emmissivity
  StefBo  = 5.67e-8 # Stefan-Boltzmann constant (W/m2/K4)

  TairK=Tair+C2K
  TseaK=Tsea+C2K

  cff  = (0.7859+0.03477*Tair)/(1.0+0.00412*Tair)
  sat  = 10.0**cff   # saturation vapor pressure (hPa or mbar)
  vap  = sat*Rh      # water vapor pressure (hPa or mbar)

  cff2 = TairK**3
  cff1 = TairK**4

  return -emmiss*StefBo*(cff1*(0.39-0.05*np.sqrt(vap))*(1.0-0.6823*cloud**2)+cff2*4.0*(TseaK-TairK))


def pvap(T):
  '''Saturation vapor pressure
  Using Augustt-Roche-Magnus (or Magnus-Tetens or Magnus) equation
    T in Celsius
    output in Pa

  Reference:
    Alduchov, O.A. and R.E. Eskridge, 1996.
    Improved Magnus Form Approximation of Saturation Vapor Pressure.
    J. Appl. Meteor., 35, 601–609
    https://doi.org/10.1175/1520-0450(1996)035<0601:IMFAOS>2.0.CO;2~

  Formulation:
    100*6.1094*np.exp(17.625*T/(T+243.04))

  Other formulations (see Table of reference), ex: Tetens (1930):
    s=100*6.11* 10**(7.5*T/(237.3+T))

  '''
  return 1000*0.61094*np.exp(17.625*T/(T+243.04))


def relative_humidity(T,Td):
  '''
  Relative humidity [0..1] from air temperature and dew point temperature
  T, Td: air and dew point temperature (C)

       saturation vapor pressure at Td (actual vapor pressure)
  RH = -------------------------------
       saturation vapor pressure at T

  '''
  return pvap(Td)/pvap(T)

