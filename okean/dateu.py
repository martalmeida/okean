# -*- coding: utf-8 -*-
'''dates utilities'''

try: basestring
except NameError: basestring=str # for python 3

import datetime
import operator

def fix_date(date):
  from calc import isiterable
  from air_sea import julianmd, greg2
  '''
  Support for "wrong" dates like (2011,13,-1)
  Input must be an iterable, like list or tuple, not string

  Example:
    fix_date((2011,13,-1))
    returns datetime.datetime(2011, 12, 30, 0, 0)
  '''

  if isiterable(date) and not isinstance(date,basestring):
    try:
      y,m,d=date
      hh,mm,ss=0.,0.,0.
    except: y,m,d,hh,mm,ss=date

    y,m,d,hh,mm,ss=greg2(julianmd(y,m,d,hh,mm,ss)-julianmd(y,1,1),y)
    return datetime.datetime(y,m,d,hh,mm,ss)

  else:
    return date


def parse_date(date,fmt='defaults',retfmt=False):
  '''
  Returns datetime object from string or iterable.
  If date cannot be parsed, returns False.
  If retfmt, the format used by date (as string) is returned.

  Example:
    '20110101', '201101011230', '2011-01-01 12:30', ...
    (2011,1,1), (2011,1,1,12,30), ...

  For other formats, ex: '01/01/2011 12h30min 30secs' specify the format
  fmt for date as string

  Example:
    '%m/%d/%Y %Hh%Mmin'

  >>> parse_date('01/01/2011 12h30min 50sec',fmt='%m/%d/%Y %Hh%Mmin %Ssec')
  datetime.datetime(2011, 1, 1, 12, 30,50)

  '''

  d,rfmt=False,''
  if isinstance(date,datetime.date):
    if isinstance(date,datetime.datetime): d=date
    else: d=datetime.datetime.fromordinal(date.toordinal())
  elif not isinstance(date,basestring): # iterable
    d=datetime.datetime(*date)
  else: # string

    if fmt=='defaults':
      if date.isdigit():
        fmt='%Y%m%d','%Y%m%d%H','%Y%m%d%H%M','%Y%m%d%H%M%S'
      else:
        fmt='%Y-%m-%d %H:%M:%S','%Y-%m-%d %H:%M','%Y-%m-%d %H', '%Y-%m-%d'

    if isinstance(fmt,basestring): fmt=[fmt]

    for i in fmt:
      try:
        d=datetime.datetime.strptime(date,i)
        rfmt=i
      except: pass

      if not d is False: break

  if retfmt: return d,rfmt
  else: return d


def next_date(date,n=1,samefmt=True):
  d,fmt=parse_date(date,retfmt=True)
  if samefmt and fmt:
    return (d+datetime.timedelta(days=n)).strftime(fmt)
  else: return d+datetime.timedelta(days=n)


def next_month(y,m,n=1):
  y1,m1=y+(m+n)//12,(m+n)%12
  if m1==0: y1,m1=y1-1,12
  return y1,m1


def drange(date0,date1,inclast=False,samefmt=True):
  date,fmt=parse_date(date0,retfmt=True)
  date1=parse_date(date1)
  res=[]
  if inclast: op=operator.le
  else: op=operator.lt
  while op(date,date1):
    if samefmt and fmt: res+=[date.strftime(fmt)]
    else: res+=[date]
    date=next_date(date,1)

  return res


def mrange(y0,m0,y1,m1,day=1,fmt=False):
  res=[]
  date0=datetime.datetime(y0,m0,day)
  date1=datetime.datetime(y1,m1,day)
  date=date0
  while date<=date1:
    res+=[date]
    y,m=next_month(date.year,date.month)
    date=datetime.datetime(y,m,day)

  if fmt: return [i.strftime(fmt) for i in res]
  else: return res


def mndays(y,m):
  y1,m1=next_month(y,m)
  return next_date((y1,m1,1),-1).day


def currday(local=False):
  if local: return datetime.datetime.now()
  else:     return datetime.datetime.utcnow()


def date_diff(date0,date1,asDays=False):
  dt=parse_date(date0)-parse_date(date1)
  if asDays: dt=dt.days+(dt.seconds+dt.microseconds/1000.)/86400.
  return dt


def yearday(date,y0=False):
  if y0:
    return date_diff(date,datetime.datetime(y0,1,1),True)
  else:
    date=parse_date(date)
    date0=datetime.datetime(date.year,1,1)
    return date_diff(date,date0,True)


def yearday2(yd,y0):
  dt=datetime.timedelta(days=yd)
  date0=datetime.datetime(y0,1,1)
  return date0+dt

def month_names(m,lang='en',abb=True):
  '''
  month names, in portuguese (lang 'pt'), spanish (lang 'es') and
  english (lang 'en')
  if abb, abbreviation name is returned
  m is the month number (Jan is 1)
  '''

  try:    pt_mar=unicode('Março','utf8') # python 2
  except: pt_mar='Março'

  mpt='Janeiro','Fevereiro',pt_mar,'Abril','Maio','Junho',\
      'Julho','Agosto','Setembro','Outubro','Novembro','Dezembro'

  men='January','February','March','April','May','June','July','August',\
      'September','October','November','December'

  mes='Enero','Febrero','Marzo','Abril','Mayo','Junio','Julio','Agosto',\
      'Septiembre','Noviembre','Diciembre'


  if   lang=='pt': mes=mpt
  elif lang=='es': mes=mes
  elif lang=='en': mes=men
  else: return 'TODO for %s' % lang

  if abb: return mes[m-1][:3]
  else: return mes[m-1]

