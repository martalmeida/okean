import Tkinter as tk
import tkFileDialog
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import numpy as np
import os
import pylab
from mpl_toolkits.basemap import Basemap
from mpl_toolkits. basemap import cm as basemap_cm

from okean import calc, netcdf, ticks, cookbook as cbk, pl_tools, locations
from okean.roms import roms
from okean.air_sea import greg2
from okean.dateu import month_names


pylab.rcParams['font.size']=8

import time as pytime


rGUIS=[]

def readrc(rcname='romsguirc'):
  import ConfigParser
  rc=ConfigParser.ConfigParser()
  rc.read(rcname)

  def parse(val):
#    v=val.split()
    v=val
    try: v=eval(v)
    except: pass
    return v
    

  res=cbk.odict()
  for sec in rc.sections():
    tmp=dict(rc.items(sec))
    for i in tmp.keys():
      res[sec.lower()+'.'+i]=parse(tmp[i])

  return res


#def get_cities(files='auto'):
#  if files=='auto':
#    from glob import glob
#    p=os.path.dirname(os.path.abspath(__file__))
#    Cities=glob(os.path.join(p,'cities*.txt'))
#  else: Cities=files
#
#  from string import join as sjoin
#
#  cities=[]
#  for c in Cities:
#    f=open(c).readlines()
#    for i in f:
#      if i.strip().startswith('#'): continue
#      tmp0=i.strip().split(',')
#      name=tmp0[0].strip()
#
#      tmp=tmp0[1].strip().split()
#      country=tmp[1]
#
#      for j in range(len(tmp)):
#        if tmp[j].isdigit(): break
#
#      country=sjoin(tmp[:j])
#
#      lat=float(tmp[j])+float(tmp[j+1])/60.
#      slat=tmp[j+2]
#      if slat=='S': lat=-lat
#      lon=float(tmp[j+3])+float(tmp[j+4])/60.
#      slon=tmp[j+5]
#      if slon=='W': lon=-lon
#
#      cities+=[{'name': name,'country':country,'lon':lon,'lat':lat}]
#
#  cities_lon=[]
#  cities_lat=[]
#  cities_name=[]
#  for c in cities:
#    cities_lon+=[c['lon']]
#    cities_lat+=[c['lat']]
#    cities_name+=[c['name']]
#
#  return cities_lon,cities_lat,cities_name #cities



try:
  locs=locations.Locations()
  cities_lon,cities_lat,cities_name,cities_c=locs.inside([-180,180],[-90,90])
except: cities_name=False

def get_ucmaps(name):
  from okean.pl_tools import ucmaps
  return getattr(ucmaps(),name)

def get_colormaps():
  matlabCmaps='jet hsv hot cool spring summer autumn winter gray bone'.split()+\
              'copper pink lines prism flag'.split()

  colorbrewer={
  'sequencial':'BuGn BuPu GnBu OrRd PuBu PuBuGn PuRd RdPu YlGn YlGnBu'.split()+\
               'YlOrBr YlOrRd Blues Greens Greys Oranges Purples Reds'.split(),
   'diverging': 'BrBG PiYG PRGn PuOr RdBu RdGy RdYlBu RdYlGn Spectral'.split(),
   'qualitative': 'Accent Dark2 Paired Pastel1 Pastel2 Set1 Set2 Set3'.split()}


  res=cbk.odict()
  res['matlab like']=[]
  res['gist']=[]
  res['colorbrewer']={'sequencial':[],'diverging':[],'qualitative':[]}
  res['other']=[]
  res['usersdef']=[]
  res['basemap']=[]

  from okean.pl_tools import ucmaps
  u=ucmaps()
  for k in u.cmaps:#.keys():
    res['usersdef']+=['_'+k]#u.cmaps[k].name]

  bcm=basemap_cm._cmapnames
  bcm.sort()
  res['basemap']=bcm

  try:    names=pylab.cm._cmapnames
  except: names=pylab.cm.cmapnames
  for n in names:
    for k in colorbrewer.keys():
      if n in colorbrewer[k]:
         res['colorbrewer'][k]+=[n]
         break
    else:
      if n.startswith('gist'): res['gist']+=[n]
      elif n in matlabCmaps:    res['matlab like']+=[n]
      else: res['other']+=[n]

  # sort:
  for k in res.keys():
    if isinstance(res[k],dict):
      for q in res[k].keys(): res[k][q].sort()
    else: res[k].sort()


  return res


def get_derived(derived_cfg='auto'):
  if derived_cfg=='auto':
    p=os.path.dirname(os.path.abspath(__file__))
    derived_cfg=os.path.join(p,'romsgui.derived')
    
  import ConfigParser
  c=ConfigParser.RawConfigParser()
  c.read(derived_cfg)

  res=cbk.odict()
  req=cbk.odict() # required vars
  vnames=c.sections()
  for v in vnames:
    tmp=dict(c.items(v))
    if tmp.has_key('req'):
      req[v]=tmp['req'].split(',')
      tmp.pop('req')

    res[v]=tmp

  return res, req


class gui_init:
  '''
  Deals with gui folders and configuration file. Also loads the
  configuration settings from the user's config file or, if not present,
  the default settings are used.
  '''

  _namerc     = 'romsguirc'
  _namedir    = '.romsgui'
  _nametmpdir = 'tmp'

  def __init__(self):
    self._instfolder  = os.path.join(os.path.dirname(__file__))
    self._default_tmp = os.path.join(self._instfolder,self._nametmpdir)
    self._default_rc  = os.path.join(self._instfolder,self._namerc)

    # set gui folders:
    self.folder    = self._instfolder
    self.tmpfolder = self._default_tmp
    self.__gui_folder()

    # set gui rc file:
    self.rcfile    = self._default_rc
    self.__gui_rc()


  def __gui_folder(self):
    '''
    Sets gui folders at user's HOME, creates if necessary. If not
    possible the default gui folders are used
    '''

    if 'HOME' in os.environ.keys():
      home=os.environ['HOME']
    elif 'TMPDIR' in os.environ.keys():
      home=os.environ['TMPDIR']
    else:
      home='.'
      print 'CANNOT find home folder!!'
      print 'rgui folder will be created in current location'

    rguiDir=os.path.join(home,self._namedir)
    rguiTmpDir=os.path.join(home,self._namedir,self._nametmpdir)
    if not os.path.isdir(rguiDir):
      print 'creating gui folder %s' % rguiDir
      try:
        os.mkdir(rguiDir)
      except:
        print 'unable to create gui folder'
        rguiDir=False

    if not os.path.isdir(rguiTmpDir):
      print 'creating gui tmp folder %s' % rguiTmpDir
      try:
        os.mkdir(rguiTmpDir)
      except:
        print 'unable to create gui tmp folder'
        rguiTmpDir=False

    if rguiDir!=False: self.folder=rguiDir
    if rguiTmpDir!=False: self.tmpfolder=rguiTmpDir


  def  __gui_rc(self):
    '''
    Sets gui configuration file, copies to gui folder if necessary.
    If not possible uses the default configuration file.
    '''

    if self.folder != self._instfolder:
      rcfile=os.path.join(self.folder,self._namerc)
      if not os.path.isfile(rcfile):
        try:
          print 'creating rc file %s' % rcfile
          import shutil
          shutil.copy2(self._default_rc,rcfile)
          self.rcfile=rcfile
        except:
          print 'unable to create rc file'
      else:
        self.rcfile=rcfile




class rgui:
  def __init__(self):
    TODO='''1. nao preciso ler os widgets editaveis, mais rapido usar uma tk var,
          2. corrigir o swapp para guardar nao o roms obj mas o nome do ficheiro
          pois quando se recalcula o tempo muda se o objecto, ie, como se
          toda a swapp fosse apagada,
          5. tirar para do titulo da swapp, pelo menos a parte do y orig
          6. refazer slice
          8. no menu tools colocar opcao de numero de slices na cache
          9. criar file de configuracao
          12. --> colocar o sliceuv na cache !!!!!!!!!!!!!!!!!
          '''

#    print TODO

    self.url='http://pong.tamu.edu/~mma'

    self.bg='#dddddd'
    self.fg='#333333'
    self.gui_size=640,660

    self.widgets={}

    self.active_tool = None

    global rGUIS

    if not rGUIS:
      self.root0=tk.Tk()
      rGUIS=[self.root0]
      self.root0.withdraw()
    else:
      self.root0=rGUIS[0]

    root=tk.Toplevel(self.root0)
    rGUIS+=[root]
    print rGUIS

    root.protocol("WM_DELETE_WINDOW", self.closegui)

    root.config(width=self.gui_size[0],height=self.gui_size[1],bg=self.bg)

    root.wm_title('romsgui')
    try:
      p=os.path.dirname(os.path.abspath(__file__))
      icon=os.path.join(p,'rgui.png')
      import ImageTk
      #    try:
      self.icon=ImageTk.PhotoImage(ImageTk.Image.open(icon))
      #    self.root0.tk.call('wm', 'iconphoto', self.root0._w, self.icon)
      root.tk.call('wm', 'iconphoto', root._w, self.icon)
      #    except: pass
    except: pass


    self.root=root


    self.swapp={}
    self.swapp['last_varname']=''
    self.swapp['last_proj']={}
    self.swapp['last_slices']=[]
    self.swapp['last_slices_data']=[]
    self.swapp['n_slices']=100
    self.swapp['caxis_mode']='auto'

    self.swapp['last_savename']=''

    self.options={}
    self.options['graf']='pcolor'


    self.options['draw_rivers']     = tk.IntVar()
    self.options['draw_countries']  = tk.IntVar()
    self.options['draw_scale']      = tk.StringVar()
    self.options['scale_length']    = tk.StringVar()
    self.options['draw_states']     = tk.IntVar()
    self.options['fill_continents'] = tk.IntVar()
    self.options['coastline_res']   = tk.StringVar()
    self.options['draw_grdborder']  = tk.IntVar()

    self.options['draw_rivers'].set(0)
    self.options['draw_countries'].set(0)
    self.options['draw_scale'].set('none')
    #self.options['scale_length'].set(500)
    self.options['scale_length'].set('auto')
    self.options['draw_states'].set(0)
    self.options['fill_continents'].set(1)
    self.options['coastline_res'].set('low')
    self.options['draw_grdborder'].set(1)

    self.options['projection']   = tk.StringVar()
    self.options['projection'].set('merc')

    self.options['colormap'] = tk.StringVar()
    self.options['colormap'].set('jet')
    self.options['colormap_r'] = tk.IntVar()
    self.options['colormap_r'].set(0)

    self.options['font_name']   = tk.StringVar()
    self.options['font_name'].set('helvetica')
    self.options['font_size']   = tk.StringVar()
    self.options['font_size'].set('7')

    self.options['zslice_SN']   = tk.IntVar()
    self.options['zslice_SN'].set(0)

    self.options['zslices_depths']=[]


    self.options['gui_size']=tk.StringVar()
    self.options['gui_size'].set(str(self.gui_size[0])+'x'+str(self.gui_size[1]))

    self.options['vslice_npts']   = tk.StringVar()
    self.options['vslice_npts'].set('auto')

    self.options['vslice_nlpts']   = tk.StringVar()
    self.options['vslice_nlpts'].set('auto')


    # cache:
    self.cache={}


    # menus:
    menubar=tk.Menu(self.root)

    # menu file:
    filemenu=tk.Menu(menubar,tearoff=0)
    filemenu.add_command(label='save as', command=self.savefig)
    filemenu.add_command(label='save', command=self.__savefig)
    filemenu.add_command(label='close',command=self.closegui)
    menubar.add_cascade(label='File',menu=filemenu)

    # menu options:
    optmenu=tk.Menu(menubar,tearoff=0)

    optmenu.add_checkbutton(label='show rivers',     command=self.show,variable=self.options['draw_rivers'])
    optmenu.add_checkbutton(label='show countries',  command=self.show,variable=self.options['draw_countries'])


    scale=tk.Menu(optmenu,tearoff=0)
    for i in 'none','NW','NE','SW','SE':
      scale.add_radiobutton(label=i, command=self.show,variable=self.options['draw_scale'])

    scaleL=tk.Menu(scale,tearoff=0)
    for i in ('auto',50,100,200,500,1000,2000):
      scaleL.add_radiobutton(label=str(i), command=self.show,variable=self.options['scale_length'])

    scale.add_separator()
    scale.add_cascade(label='length',menu=scaleL)
    optmenu.add_cascade(label='show scale',menu=scale)

    optmenu.add_checkbutton(label='show states',     command=self.show,variable=self.options['draw_states'])
    optmenu.add_checkbutton(label='show continents', command=self.show,variable=self.options['fill_continents'])
    optmenu.add_checkbutton(label='show grid border',command=self.show,variable=self.options['draw_grdborder'])
    optmenu.add_separator()

    cres=tk.Menu(optmenu,tearoff=0)
    cres.add_radiobutton(label='none',         command=self.show,variable=self.options['coastline_res'])
    cres.add_separator()
    cres.add_radiobutton(label='crude',        command=self.show,variable=self.options['coastline_res'])
    cres.add_radiobutton(label='low',          command=self.show,variable=self.options['coastline_res'])
    cres.add_radiobutton(label='intermediate', command=self.show,variable=self.options['coastline_res'])
    cres.add_radiobutton(label='high',         command=self.show,variable=self.options['coastline_res'])
    cres.add_radiobutton(label='full',         command=self.show,variable=self.options['coastline_res'])
    optmenu.add_cascade(label='coastline',menu=cres)


    def add_sub(main,label,options,cmd,var):
      tmp=tk.Menu(main,tearoff=0)
      for p in options: tmp.add_radiobutton(label=p, command=cmd,variable=var)
      main.add_cascade(label=label,menu=tmp)
      return tmp

    optmenu.add_separator()
    projs='merc','cyl','eqdc','mill','aeqd','poly','gnom','stere'
    prj=add_sub(optmenu,'projection',projs,self.show,self.options['projection'])


    optmenu.add_separator()
    colormaps=get_colormaps()
    cmaps=tk.Menu(optmenu,tearoff=0)
    cmaps.add_checkbutton(label='invert',command=self.show,variable=self.options['colormap_r'])
    cmaps.add_separator()
    for k in colormaps.keys():
      if isinstance(colormaps[k],dict):
        tmp=tk.Menu(cmaps,tearoff=0)
        for q in colormaps[k].keys(): add_sub(tmp,q,colormaps[k][q],self.show,self.options['colormap'])
        cmaps.add_cascade(label=k,menu=tmp)
      else:
        add_sub(cmaps,k,colormaps[k],self.show,self.options['colormap'])

    optmenu.add_cascade(label='colormap',menu=cmaps)


    optmenu.add_separator()
    optmenu.add_checkbutton(label='zslice surf NaNs',command=self.show,variable=self.options['zslice_SN'])


    # vertical slice number of points:
    optmenu.add_separator()

    vslc=tk.Menu(optmenu,tearoff=0)
    nslc=tk.Menu(optmenu,tearoff=0)
    for i in 'auto',50,100,200,500:
      nslc.add_radiobutton(label=str(i), variable=self.options['vslice_npts'])

    nlslc=tk.Menu(scale,tearoff=0)
    for i in 'auto',2,3:
      nlslc.add_radiobutton(label=str(i), variable=self.options['vslice_nlpts'])

    vslc.add_cascade(label='slice points',menu=nslc)
    vslc.add_cascade(label='line points',menu=nlslc)
    optmenu.add_cascade(label='vertical slice',menu=vslc)


    menubar.add_cascade(label='Options',menu=optmenu)


    # menu tools:
    toolsmenu=tk.Menu(menubar,tearoff=0)
    toolsmenu.add_command(label='clear cache',command=self.__clear_swapp)

    fonts=tk.Menu(toolsmenu,tearoff=0)
    fontn=tk.Menu(toolsmenu,tearoff=0)
    guisz=tk.Menu(toolsmenu,tearoff=0)
    for i in '56789':
      fonts.add_radiobutton(label=i,command=self.gui_chfont,variable=self.options['font_size'])

    for i in ('helvetica','times','courier','arial'):
      fontn.add_radiobutton(label=i,command=self.gui_chfont,variable=self.options['font_name'])

    for i in ('640x660','740x560'):
      guisz.add_radiobutton(label=i,command=self.gui_chsize,variable=self.options['gui_size'])

    toolsmenu.add_cascade(label='gui fontname', menu=fontn)
    toolsmenu.add_cascade(label='gui fontsize', menu=fonts)
    toolsmenu.add_cascade(label='gui size',     menu=guisz)


    menubar.add_cascade(label='Tools',menu=toolsmenu)

    # menu help:
    helpmenu=tk.Menu(menubar,tearoff=0)
    helpmenu.add_command(label='about',command=self.__about)
    helpmenu.add_command(label='homepage',command=self.__homepage)
    menubar.add_cascade(label='Help',menu=helpmenu)


    self.root.config(menu=menubar)
    self.widgets['menubar']=menubar
    self.widgets['menu_options']=optmenu
    self.widgets['menu_tools']=optmenu
    self.widgets['menu_help']=helpmenu
    self.widgets['menu_file']=filemenu

    # ------------------------------------

    Ht=.03
    Hm=.65
    Hb=1-(Ht+Hm)

    Ll=.15
    Lr=.15
    Lc=1-(Ll+Lr)

    # frame at top:
    ftop = tk.Frame(root, bg=self.bg,borderwidth=0)
    ftop.place(x=0,y=0,relwidth=1,relheight=Ht)
    self.widgets['ftop']=ftop

    # frame at left:
    fleft = tk.Frame(root, bg=self.bg,borderwidth=0)
    fleft.place(relx=0,rely=Ht,relwidth=Ll,relheight=Hm)
    #fleft.place(x=0,rely=Ht,width=100,relheight=Hm)
    self.widgets['fleft']=fleft

    # frame at center (Figure):
    f = Figure(frameon=True,facecolor=self.bg)
    canvas = FigureCanvasTkAgg(f,master=root)
    canvas.show()
    canvas.get_tk_widget().place(relx=Ll,rely=Ht,relwidth=Lc,relheight=Hm)
#####    canvas.get_tk_widget().pack()
    self.figure=f
    self.canvas=canvas

    # must let border 1, otherwise a strange white border appears at left and top !?
    self.canvas.get_tk_widget().config(highlightthickness=1,
                                       highlightcolor=self.bg,
                                       highlightbackground=self.bg)

    self.canvas.get_tk_widget().config(relief='flat',borderwidth=1,bg=self.bg)

    # frame at right:
    fright = tk.Frame(root, bg=self.bg,borderwidth=0)
    fright.place(relx=1-Lr,rely=Ht,relwidth=Lr,relheight=Hm)
    self.widgets['fright']=fright

    #frame at bottom:
    fbot = tk.Frame(root, bg=self.bg,borderwidth=0)
    fbot.place(x=0,rely=1-Hb,relwidth=1,relheight=Hb)
    self.widgets['fbot']=fbot

    #
    # contents of frame top:
    b = tk.Button(self.widgets['ftop'], text="roms files",command=self.select_file)
    b.place(relx=0,y=0,relw=.1,height=15)

    yorig=tk.Entry(self.widgets['ftop'])
    yorig.place(relx=1-.1,y=0,relw=.1,height=15)
    yorig.insert(0, "Y orig")
    yorig.bind('<Return>',self.show)
    self.widgets['yorig']=yorig

    ob=tk.Label(self.widgets['ftop'],text='',bg=self.bg)
    ob.place(relx=.1,y=0,relw=.8,h=15)
    self.widgets['finfo']=ob


#    clearSwapp=tk.Button(self.widgets['ftop'],text='rm swapp',command=self.__clear_swapp)
#    clearSwapp.place(relx=1-.205,y=1,relw=.1,height=15)
#    self.widgets['clearSwapp']=clearSwapp

    # contents of frame left:

    # list box for file vnames:
    scroll=tk.Scrollbar(self.widgets['fleft'],orient=tk.VERTICAL)
    listbox = tk.Listbox(self.widgets['fleft'],yscrollcommand=scroll.set)
    scroll.config(command=listbox.yview)

    listbox.place(x=2,y=2,width=70,height=180)
    scroll.place(x=72,y=2,width=15,height=180)
    listbox.bind('<ButtonRelease-1>',self.__show_after_select_var)
    self.widgets['varsList']=listbox
    ##self.widgets['varsScroll']=scroll

    # list box for derived vnames:
    scroll=tk.Scrollbar(self.widgets['fleft'],orient=tk.VERTICAL)
    listbox = tk.Listbox(self.widgets['fleft'],yscrollcommand=scroll.set)
    scroll.config(command=listbox.yview)

    listbox.place(x=2,y=190,width=70,height=90)
    scroll.place(x=72,y=190,width=15,height=90)
    listbox.bind('<ButtonRelease-1>',self.__show_after_select_Dvar)
    self.widgets['DvarsList']=listbox


#    hvals = tk.Entry(self.widgets['fleft'])
#    hvals.insert(0, "auto")##########200 1000")
#    hvals.place(x=2,y=184,width=80,height=20)
#    hvals.bind('<Return>',self.show)
#    self.widgets['hvals']=hvals



##    zlev = tk.Entry(self.widgets['fleft'])
##    zlev.insert(0, "-10")
##    zlev.place(x=2,y=206,width=80,height=20)
##    zlev.bind('<Return>',self.show)
##    #zlev.bind('<FocusOut>',self.show)
##    self.widgets['zlev']=zlev

###    cax = tk.Entry(self.widgets['fleft'])
###    cax.insert(0, "auto")
###    cax.place(x=2,y=228,width=80,height=20)
###    cax.bind('<Return>',self.__show_after_cax)
###    self.widgets['cax']=cax

#    # checkbutton:
#    var=tk.IntVar()
#    ck=tk.Checkbutton(self.widgets['fleft'],variable=var)
#    ck.place(x=70,y=228,width=12,height=20)
#    self.widgets['cax_ck']=var

#    # slicell
#    slc=tk.Button(self.widgets['fleft'], text="slice",command=self.draw_seg)
#    slc.place(x=2,y=250,width=50,height=20)

###    # contourf:
###    ob=tk.Button(self.widgets['fleft'], text="contourf",command=self.__draw_contourf)
###    ob.place(x=2,y=270,width=40,height=20)
###    # pcolor:
###    ob=tk.Button(self.widgets['fleft'], text="pcolor",command=self.__draw_pcolor)
###    ob.place(x=42,y=270,width=40,height=20)
###    # contour:
###    ob=tk.Button(self.widgets['fleft'], text="contour",command=self.__draw_contour)
###    ob.place(x=2,y=290,width=40,height=20)


    # ----------------------------------------------------------------
    # time
    # ----------------------------------------------------------------
    time = tk.Entry(self.widgets['fbot'])
    time.insert(0, "0")
    relw=.05
    x0=Ll+Lc/2.-relw/2
    y=2
    h=20
    time.place(relx=x0,y=y,relw=relw,height=h)
    time.bind('<Return>',self.show)
    self.widgets['time']=time

    # increase, decrease:
    Relw=.03
    x=x0-Relw
    timeL=tk.Button(self.widgets['fbot'], text="<",command=self.time_dec)
    timeL.place(relx=x,y=y,relw=Relw,h=h)
    x=x0+relw
    timeP=tk.Button(self.widgets['fbot'], text=">",command=self.time_inc)
    timeP.place(relx=x,y=y,relw=Relw,h=h)

    # update:
    ttime=tk.Button(self.widgets['fbot'],text='0',command=self.time_show)
    x=x0+relw+Relw
    ttime.place(relx=x,y=y,relw=relw,height=h)
    self.widgets['ttime']=ttime

    # labels:
    lab=tk.Label(self.widgets['fbot'],text='time',bg=self.bg,fg=self.fg)
    x=x0-Relw
    w=relw+Relw*2
    y=y+h
    lab.place(relx=x,y=y,relwidth=w,height=h)

    lab=tk.Label(self.widgets['fbot'],text='updt',bg=self.bg,fg=self.fg)
    x=x0+relw+Relw
    w=relw+Relw*2
    lab.place(relx=x,y=y,relw=relw,height=h)


    # ----------------------------------------------------------------
    # xlim, ylim
    # ----------------------------------------------------------------
    # lon1:
    lon1  = tk.Entry(self.widgets['fbot'])
    lon1.bind('<Return>',self.show)
    lon1L = tk.Button(self.widgets['fbot'],text="<",command=self.lon1_dec)
    lon1P = tk.Button(self.widgets['fbot'],text=">",command=self.lon1_inc)
    lon1.insert(0, "-180")
    relw=.07
    Relw=.03
    x0=Ll+Relw
    y=2
    h=20
    lon1.place(relx=x0,y=y,relw=relw,height=h)
    x=x0-Relw
    lon1L.place(relx=x,y=y,relw=Relw,height=h)
    x=x0+relw
    lon1P.place(relx=x,y=y,relw=Relw,height=h)
    self.widgets['lon1']=lon1

    # labels:
    lab=tk.Label(self.widgets['fbot'],text='xlim[0]',bg=self.bg,fg=self.fg)
    y1=y+h
    w=relw+Relw*2
    lab.place(relx=x0,y=y1,relw=relw,height=h)

    # lon2:
    lon2  = tk.Entry(self.widgets['fbot'])
    lon2.bind('<Return>',self.show)
    lon2L = tk.Button(self.widgets['fbot'],text="<",command=self.lon2_dec)
    lon2P = tk.Button(self.widgets['fbot'],text=">",command=self.lon2_inc)
    lon2.insert(0, "+180")
    x0=Ll+Lc-2*Relw-relw
    lon2.place(relx=x0,y=y,relw=relw,height=h)
    x=x0-Relw
    lon2L.place(relx=x,y=y,relw=Relw,height=h)
    x=x0+relw
    lon2P.place(relx=x,y=y,relw=Relw,height=h)
    self.widgets['lon2']=lon2

    # labels:
    lab=tk.Label(self.widgets['fbot'],text='xlim[1]',bg=self.bg,fg=self.fg)
    lab.place(relx=x0,y=y1,relw=relw,height=h)


    # lat1:
    relw=.5
    Relw=.25
    lat1  = tk.Entry(self.widgets['fright'])
    lat1.bind('<Return>',self.show)
    lat1L = tk.Button(self.widgets['fright'],text="v",command=self.lat1_dec)
    lat1P = tk.Button(self.widgets['fright'],text="^",command=self.lat1_inc)
    lat1.insert(0, "-90")
    x=0
    rh=.05
    y=1-rh
    lat1L.place(relx=x+relw/2-Relw/2,rely=1-rh,relw=Relw,relheight=rh)
    lat1.place(x=x,rely=1-2*rh,relw=relw,relheight=rh)
    lat1P.place(relx=x+relw/2-Relw/2,rely=1-3*rh,relw=Relw,relheight=rh)
    self.widgets['lat1']=lat1

    # lat2:
    lat2  = tk.Entry(self.widgets['fright'])
    lat2.bind('<Return>',self.show)
    lat2L = tk.Button(self.widgets['fright'],text="v",command=self.lat2_dec)
    lat2P = tk.Button(self.widgets['fright'],text="^",command=self.lat2_inc)
    lat2.insert(0, "+90")
    x=0
    y=0
    lat2L.place(relx=x+relw/2-Relw/2,y=y+2*h,relw=Relw,height=h)
    lat2.place(relx=x,y=h,relw=relw,height=h)
    lat2P.place(relx=x+relw/2-Relw/2,y=y,relw=Relw,height=h)
    self.widgets['lat2']=lat2

    # labels:
    #


    # z or s level:
    zlev = tk.Entry(self.widgets['fright'])
    zlev.insert(0, "-10")
    x=0
    y=5*h
    zlev.place(relx=x,y=y+h,relwidth=relw,height=h)
    zlev.bind('<Return>',self.show)
    #zlev.bind('<FocusOut>',self.show)
    self.widgets['zlev']=zlev

    # z or s increase or decrease:
    zlevL = tk.Button(self.widgets['fright'],text="v",command=self.zlev_dec)
    zlevP = tk.Button(self.widgets['fright'],text="^",command=self.zlev_inc)
    zlevL.place(relx=x+relw/2-Relw/2,y=y+2*h,relw=Relw,height=h)
    zlevP.place(relx=x+relw/2-Relw/2,y=y,relw=Relw,height=h)

    # zlev label:
    zlevLab=tk.Label(self.widgets['fright'],text='s or z',bg=self.bg,fg=self.fg)
    zlevLab.place(relx=x,y=y+3*h,relwidth=relw,height=h)


    # reset lims:
    ll_reset = tk.Button(self.widgets['fbot'],text="x",command=self.reset_lims)
    relx=1-Lr
    y=0
    Relw=.03
    ll_reset.place(relx=relx,y=y,relw=Relw,height=h)
    # zoom out:
    ll_out = tk.Button(self.widgets['fbot'],text="-",command=self.lims_out)
    relx=1-Lr+Relw
    ll_out.place(relx=relx,y=y,relw=Relw,height=h)
    # zoom in:
    ll_in = tk.Button(self.widgets['fbot'],text="+",command=self.lims_in)
    relx=1-Lr+Relw*2
    ll_in.place(relx=relx,y=y,relw=Relw,height=h)



###    # separate plot:
###    sep = tk.Button(self.widgets['fright'],text="separate",command=self.__show_newfig)
###    sep.place(x=0,rely=.5,relw=.7,height=h)

    # ----------------------------------------------------------------
    # vectors:
    # ----------------------------------------------------------------
    # dvectors [::d]
    dvec  = tk.Entry(self.widgets['fbot'])
    dvec.bind('<Return>',self.show)
    dvecL = tk.Button(self.widgets['fbot'],text="<",command=self.dvec_dec)
    dvecP = tk.Button(self.widgets['fbot'],text=">",command=self.dvec_inc)
    dvec.insert(0, "0")
    relw=.04
    Relw=.03
    x=Relw
    y=h*2.5
    yb0=y
    h=20
    dvec.place(relx=x,y=y,relw=relw,height=h)
    x=0
    dvecL.place(relx=x,y=y,relw=Relw,height=h)
    x=Relw+relw
    dvecP.place(relx=x,y=y,relw=Relw,height=h)
    self.widgets['dvec']=dvec

    # svectors (scale)
    svec  = tk.Entry(self.widgets['fbot'])
    svec.bind('<Return>',self.show)
    svecL = tk.Button(self.widgets['fbot'],text="<",command=self.svec_dec)
    svecP = tk.Button(self.widgets['fbot'],text=">",command=self.svec_inc)
    svec.insert(0, "1")
    x0=3.1*Relw+relw
    svec.place(relx=x0,y=y,relw=relw,height=h)
    x=x0-Relw
    svecL.place(relx=x,y=y,relw=Relw,height=h)
    x=x0+relw
    svecP.place(relx=x,y=y,relw=Relw,height=h)
    self.widgets['svec']=svec

    # labels:
    dvecLab=tk.Label(self.widgets['fbot'],text='vec ::n',bg=self.bg,fg=self.fg)
    y=y+h
    x=0
    w=2*Relw+relw
    dvecLab.place(relx=x,y=y,relw=w,height=h)

    svecLab=tk.Label(self.widgets['fbot'],text='vec scale',bg=self.bg,fg=self.fg)
    x=x0-Relw
    svecLab.place(relx=x,y=y,relw=w,height=h)


    # ----------------------------------------------------------------
    # isobaths:
    # ----------------------------------------------------------------
    x=.005
    y=y+h+h/4
    w=.175
    hvals = tk.Entry(self.widgets['fbot'])
    hvals.insert(0, "auto")##########200 1000")
    hvals.place(relx=x,y=y,relwidth=w,height=h)
    hvals.bind('<Return>',self.show)
    self.widgets['hvals']=hvals

    # isobaths label:
    lab=tk.Label(self.widgets['fbot'],text='isobaths',bg=self.bg,fg=self.fg)
    y=y+h
    lab.place(relx=x,y=y,relwidth=w,height=h)


    # ----------------------------------------------------------------
    # plot type:
    # ----------------------------------------------------------------
    # contourf:
    ob=tk.Button(self.widgets['fbot'], text="contourf",command=self.__draw_contourf)
    y=yb0
    x=6*relw
    x02=x
    w=2*relw
    ob.place(relx=x,y=y,relw=w,height=h)

    # pcolor:
    ob=tk.Button(self.widgets['fbot'], text="pcolor",command=self.__draw_pcolor)
    y=y+h
    ob.place(relx=x,y=y,relw=w,height=h)

    # contour:
    ob=tk.Button(self.widgets['fbot'], text="contour",command=self.__draw_contour)
    y=y+h
    ob.place(relx=x,y=y,relw=w,height=h)

    # label:
    lab=tk.Label(self.widgets['fbot'],text='plot type',bg=self.bg,fg=self.fg)
    y=y+h
    lab.place(relx=x,y=y,relwidth=w,height=h)

    # ----------------------------------------------------------------
    # vertical slice, etc:
    # ----------------------------------------------------------------
    # z slice
    ob=tk.Button(self.widgets['fbot'], text="vertical slc",command=self.interactive_vslice)
    x=x02+w*1.1
    y=yb0
    w=2.5*relw
    ob.place(relx=x,y=y,relw=w,height=h)
    ob.bind('<Button-3>', self.__vslice)
    self.widgets['vslice']=ob

    ob=tk.Button(self.widgets['fbot'], text="profile",command=self.interactive_profile)
    y=y+h
    ob.place(relx=x,y=y,relw=w,height=h)
    ob.bind('<Button-3>', self.__profile)
    self.widgets['profile']=ob

    ob=tk.Button(self.widgets['fbot'], text="time series",command=self.interactive_time_series)
    y=y+h
    ob.place(relx=x,y=y,relw=w,height=h)
    ob.bind('<Button-3>', self.__time_series)
    self.widgets['time_series']=ob

    ob=tk.Button(self.widgets['fbot'], text="hovmuller",command=self.interactive_hovmuller)
    y=y+h
    ob.place(relx=x,y=y,relw=w,height=h)
    ob.bind('<Button-3>', self.__hovmuller)
    self.widgets['hovmuller']=ob

    # ----------------------------------------------------------------
    # zoom, separate, printe, etc:
    # ----------------------------------------------------------------

    # zoom:
    ob = tk.Button(self.widgets['fbot'],text="zoom",command=self.interactive_zoom)
    y=yb0+3.2*h
    x=x02+w*2.2
    w=2.25*relw
    ob.place(relx=x,y=y,relw=w,height=h)
    self.widgets['zoom']=ob

    # separate plot:
    ob = tk.Button(self.widgets['fbot'],text="separate",command=self.__show_newfig)
    y=y+h
    ob.place(relx=x,y=y,relw=w,height=h)

    # print (save as):
    ob = tk.Button(self.widgets['fbot'],text="print",command=self.savefig)
    y=y+h
    ob.place(relx=x,y=y,relw=w,height=h)


    # ----------------------------------------------------------------
    # clim:
    # ----------------------------------------------------------------
    cax = tk.Entry(self.widgets['fbot'])
    cax.insert(0, "auto")
    #x=x02+w*2.2
    y=yb0
    ww=2.5*relw*2
    cax.place(relx=x,y=y,relwidth=ww,height=h)
    cax.bind('<Return>',self.__show_after_cax)
    self.widgets['cax']=cax

    # label:
    lab=tk.Label(self.widgets['fbot'],text='clim',bg=self.bg,fg=self.fg)
    y=y+h
    lab.place(relx=x,y=y,relwidth=ww,height=h)


    # ----------------------------------------------------------------
    # var info:
    # ----------------------------------------------------------------
    vinfo=tk.Listbox(self.widgets['fbot'],bg=self.bg,fg=self.fg,
                     borderwidth=0,
                     highlightthickness=1,
                     highlightbackground='#c4c4c4',
                     highlightcolor='#c4c4c4')
    vinfo.place(relx=.6,rely=.7,relw=.5,relh=.3)
    self.widgets['vinfo']=vinfo

    # messages:
    msg=tk.Label(self.widgets['fbot'],justify=tk.LEFT,fg='red',bg=self.bg)
    msg.place(relx=0,rely=.9,relw=.3,relh=.1)
    self.widgets['msg']=msg

########
    if 0:
      b = tk.Button(root, text="go!",command=self.show)
      b.place(x=10,y=400,width=20,height=20)
      self.widgets['go']=b

#    '''
#    f = Figure(figsize=(5,4), dpi=100)
#    a = f.add_subplot(111)
#    canvas = FigureCanvasTkAgg(f, master=root)
#    canvas.show()
##    canvas.get_tk_widget().place(x=figpos[0],y=figpos[1],width=figpos[2],height=figpos[3])
#    canvas.get_tk_widget().place(relx=.15,y=24,relwidth=0.7,relheight=.6)
#    '''

    self.root.bind('m',self.time_inc)
    self.root.bind('l',self.time_dec)

    for k in 'tsuvhzUV':
      self.root.bind(k,self.__set_varname)

    self.root.bind('S',self.__set_zlevel)
    self.root.bind('B',self.__set_zlevel)

#    self.root.bind('s',self.__set_varname)
#    self.root.bind('u',self.__set_varname)
#    self.root.bind('v',self.__set_varname)
#    self.root.bind('h',self.__set_varname)
#    self.root.bind('z',self.__set_varname)
#    self.root.bind('<KeyPress-U>',self.__set_varname)
#    self.root.bind('<KeyPress-V>',self.__set_varname)

    self.files={}
#    self.axes=a
#    self.canvas=canvas


#    self.select_file('/home/mma/roms_his.nc')
    #self.select_file('/home/mma/longRun_down/spring_2002/gen_ini/roms_ini_spr02.nc')

    # se font for all widgets:
    self.gui_chfont()

#    self.select_file('/home/mma/works_remo/iousp/op_figs_tmp/roms_his_20100324_a_n0.nc')

    #root.mainloop()

#  def __create_varsList(self,**kargs):
#    try:
#      self.widgets['varsList'].destroy()
#      self.widgets['varsScroll'].destroy()
#    except: pass
#
#    scroll=tk.Scrollbar(self.widgets['fleft'],orient=tk.VERTICAL)
#    listbox = tk.Listbox(self.widgets['fleft'],yscrollcommand=scroll.set,**kargs)
#    scroll.config(command=listbox.yview)
#    self.widgets['vscroll']=scroll
#
#    listbox.place(x=2,y=2,width=70,height=180)
#    scroll.place(x=72,y=2,width=15,height=180)
#    listbox.bind('<ButtonRelease-1>',self.show)
#    self.widgets['varsList']=listbox
#    self.widgets['varsScroll']=scroll


  def gui_chsize(self):
    sz=self.options['gui_size'].get().split('x')
    sz=int(sz[0]),int(sz[1])
    self.root.config(width=sz[0],height=sz[1])

  def gui_chfont(self,font=False):
    if not font:
      font=self.options['font_name'].get()+' '+\
           self.options['font_size'].get()
    
    print font
    def set_font(ch):
      for k in ch.children.keys():
        set_font(ch.children[k])
        try:
          ch.children[k].config(font=font)
        except: pass
    set_font(self.root)

    # changing font does not update size of listbox with place manager
    # so replace it!!
    v=self.widgets['varsList']
    p=v.place_info()
    pp=p.copy()
    for i in p.keys():
      if isinstance(p[i],basestring) and p[i].isdigit() and int(p[i])>0:
        pp[i]=int(p[i])+1

    self.widgets['varsList'].place(**pp)
    self.root.update_idletasks()
    self.widgets['varsList'].place(**p)


  def zlev_dec(self): self.chzlev('dec')
  def zlev_inc(self): self.chzlev('inc')

  def chzlev(self,w,show=True):
    lev=self.widgets['zlev']
    z=float(lev.get())
    zd=self.options['zslices_depths']
    if len(zd)==0: return

    zd=zd.tolist()

    if self.files.has_key('His'):
      ns=range(self.files['His'].S_RHO)
    else:
      ns=[]

    if w=='dec':
      tmp=ns[:]
      tmp.reverse()
      r=tmp+zd
    else:
      tmp=zd[:]
      tmp.reverse()
      r=tmp+[ns[-1]]

    r=np.array(r)
    print r,lev

    if w=='dec':
      i=np.where(r<z)[0]
    else:
      i=np.where(r>z)[0]

    if len(i):
      lev.delete(0,tk.END)
      lev.insert(0, '%d' % r[i[0]])

      if show: self.show()


  def dvec_dec(self): self.chvec('d',-1)
  def dvec_inc(self): self.chvec('d',+1)
  def svec_dec(self): self.chvec('s',-.2)
  def svec_inc(self): self.chvec('s',+.2)

  def chvec(self,what,add=1,show=True):
    vec=self.widgets[what+'vec']
    n=float(vec.get())
    vec.delete(0,tk.END)
    if what=='d':
      vec.insert(0, '%d' % max(0,int(n)+add))
    else:
      vec.insert(0, '%.1f' % max(0,n+add))

    if show: self.show()

  def lims_out(self,xy0=False):
    self.chlon(1,-1,xy0,show=False)
    self.chlon(2, 1,xy0,show=False)
    self.chlat(1,-1,xy0,show=False)
    self.chlat(2, 1,xy0,show=True)
  def lims_in(self):
    self.chlon(1, 1,xy0,show=False)
    self.chlon(2,-1,xy0,show=False)
    self.chlat(1, 1,xy0,show=False)
    self.chlat(2,-1,xy0,show=True)

  def chlon(self,a,add,xy0=False,show=True):
    if a==1:
      lon=self.widgets['lon1']
    else:
      lon=self.widgets['lon2']

    lon1=self.widgets['lon1']
    Lon1=float(lon1.get())
    lon2=self.widgets['lon2']
    Lon2=float(lon2.get())
    dlon=(Lon2-Lon1)/5.


    lon.delete(0,tk.END)
    if a==1: Lon=Lon1+add*dlon
    else:    Lon=Lon2+add*dlon
    lon.insert(0, '%6.2f' % Lon)

    if show: self.show()



  def lon1_dec(self): self.chlon(1,-1)
  def lon1_inc(self): self.chlon(1,+1)
  def lon2_dec(self): self.chlon(2,-1)
  def lon2_inc(self): self.chlon(2,+1)


  def chlat(self,a,add,show=True):
    if a==1:
      lat=self.widgets['lat1']
    else:
      lat=self.widgets['lat2']

    lat1=self.widgets['lat1']
    Lat1=float(lat1.get())
    lat2=self.widgets['lat2']
    Lat2=float(lat2.get())
    dlat=(Lat2-Lat1)/5.

    lat.delete(0,tk.END)
    if a==1: Lat=Lat1+add*dlat
    else:    Lat=Lat2+add*dlat
    lat.insert(0, '%6.2f' % Lat)

    if show: self.show()

  def lat1_dec(self): self.chlat(1,-1)
  def lat1_inc(self): self.chlat(1,+1)
  def lat2_dec(self): self.chlat(2,-1)
  def lat2_inc(self): self.chlat(2,+1)


  def ch_axis(self,axis=False,xy0=False,zoom='on',**args):
    zoomFact=1.5

    if 'zoomFact' in args.keys(): zoomFact=args['zoomFact']

    if zoom=='off': zoomFact=1./zoomFact

    xlim0,ylim0=self.__get_grid_lims()
    dx=xlim0[1]-xlim0[0]
    dy=ylim0[1]-ylim0[0]

    if xy0 is False:
      xy0=0.5*(xlim0[1]+xlim0[0]),0.5*(ylim0[1]+ylim0[0])

    if axis is False: # use zoom
      xlim=[xy0[0]-dx/(2*zoomFact), xy0[0]+dx/(2*zoomFact)]
      ylim=[xy0[1]-dy/(2*zoomFact), xy0[1]+dy/(2*zoomFact)]

    else: # use lims as axis
     xlim=axis[:2]
     ylim=axis[2:]

    self.__set_grid_lims(show=True,lims=xlim+ylim)


  def chtime(self,add,show=True):
    time=self.widgets['time']
    t=int(time.get())
    time.delete(0,tk.END)
    try:    tmax=self.files['His'].TIME-1
    except: tmax=False
    tnew=t+add
    if not tmax is False: tnew=min(tnew,tmax)
    tnew=max(tnew,-1)
    time.insert(0, str(tnew))

    if show: self.show()

  def time_inc(self,evt=False): self.chtime(1)
  def time_dec(self,evt=False): self.chtime(-1)
  def time_show(self):
     if self.files.has_key('His'):
       self.reset_file('his')
       self.widgets['ttime'].config(text=self.files['His'].TIME-1)

  def reset_lims(self):  self.__set_grid_lims(show=True)

  def __set_grid_lims(self,show=False,lims=False):
    if lims is False:
      grd=self.files['Grid']
      lonlims=grd.lon.min(),grd.lon.max()
      latlims=grd.lat.min(),grd.lat.max()
    else:
      lonlims=lims[:2]
      latlims=lims[2:]

    lonlims=list(lonlims)
    latlims=list(latlims)
    lonlims[0]=max(lonlims[0],-180.)
    lonlims[1]=min(lonlims[1],180.)
    latlims[0]=max(latlims[0],-90.)
    latlims[1]=min(latlims[1],90.)


    nDecimals=3
    fmt='%.'+str(nDecimals)+'f'

    x1=fmt % lonlims[0]
    x2=fmt % lonlims[1]
    y1=fmt % latlims[0]
    y2=fmt % latlims[1]
##    x1='%7.3f' % lonlims[0]
##    x2='%7.3f' % lonlims[1]
##    y1='%7.3f' % latlims[0]
##    y2='%7.3f' % latlims[1]

    if float(x2)>float(x1):###int(float(x2)*1e3)>int(float(x1)*1e3):
      self.widgets['lon1'].delete(0,tk.END)
      self.widgets['lon2'].delete(0,tk.END)
      self.widgets['lon1'].insert(0,x1)
      self.widgets['lon2'].insert(0,x2)

    if float(y2)>float(y1): ###int(float(y2)*1e3)>int(float(y1)*1e3):
      self.widgets['lat1'].delete(0,tk.END)
      self.widgets['lat2'].delete(0,tk.END)
      self.widgets['lat1'].insert(0,y1)
      self.widgets['lat2'].insert(0,y2)

    if show: self.show()

  def __get_grid_lims(self):
    lonlim1=float(self.widgets['lon1'].get())
    lonlim2=float(self.widgets['lon2'].get())
    latlim1=float(self.widgets['lat1'].get())
    latlim2=float(self.widgets['lat2'].get())
    return (lonlim1,lonlim2), (latlim1,latlim2)

  def __calc_xyticks(self):
    lonlims,latlims=self.__get_grid_lims()
    return ticks.loose_label(lonlims[0],lonlims[1]),ticks.loose_label(latlims[0],latlims[1])

  def select_file(self,f=False,show=True):
    def finfo(f):
     atts=netcdf.fatt(f)
     if 'title' in atts.keys(): s=netcdf.fatt(f,'title')
     elif 'type' in atts.keys(): s=netcdf.fatt(f,'type')
     else: s=''

     if len(f)>25+1:
       return s+' <+'+f[-25:]+'>'
     else:  return s+' <'+f+'>'

    if not f:
      f=tkFileDialog.askopenfilename(parent=self.root,title='Choose a file',
         initialdir='/home/mma/longRun_down/spring_2001/',
         defaultextension = ".nc", filetypes=[("All Types", ".*"),
         ("NC", ".nc")])
    try:
      atype=netcdf.fatt(f,'type').lower()
    except:
      try: atype=netcdf.fatt(f,'title').split()[0].lower() # agrif child grids
      except: atype='grid'

    print atype
    if atype.find('grid')>=0 or atype.find('grd')>=0:
      self.files['grid']=f
      self.files['grid_finfo']=finfo(f)

      #######self.files['Grid']=roms.Grid(f)
      self.reset_file('grid')
      self.__set_grid_lims()
      self.__init_depths()#hauto_hvals()
      ################################################self.show_vars('grid')
      # finfo:
      self.widgets['finfo'].config(text=self.files['grid_finfo'])

      self.show(what='grid')

    elif atype.find('history')>=0 or atype.find('initial')>=0 or\
         atype.find('average')>=0 or atype.find('restart')>=0 or\
         atype.find('climatology')>=0:
      self.files['his']=f
      self.files['his_finfo']=finfo(f)

      if netcdf.fatt(f).has_key('grd_file') and os.path.isfile(netcdf.fatt(f,'grd_file')):
        grd=netcdf.fatt(f,'grd_file')
      elif self.files.has_key('grid') and os.path.isfile(self.files['grid']):
        grd=self.files['grid']
      else:
        grd=''

      if grd:
        self.files['grid']=grd
        ########  self.files['Grid']=roms.Grid(grd)
        self.reset_file('grid')
        #        print f, grd
        #        self.files['His']=roms.Grid(f,grd)
      else:
        self.files['grid']=f
        self.reset_file('grid')

      #######self.files['His']=roms.His(f)
      self.reset_file('his')
      self.__set_grid_lims()
      self.__init_depths()#auto_hvals()
      self.time_show()

      self.show_vars()

      # finfo:
      self.widgets['finfo'].config(text=self.files['his_finfo'])

      # if hvals is auto, set them:
      # TODO

      if show: self.show()

#####    hvals.insert(0, "auto")##########200 1000")

  def reset_file(self,type):
    Type=type[0].upper()+type[1:]
    ob=eval('roms.'+Type)
    if type!='grid' and self.files.has_key('grid'):
      self.files[Type]=ob(self.files[type],grd=self.files['grid'])
    else:
      self.files[Type]=ob(self.files[type])

  def show_vars(self,type='his',show_grid_vars=False,**opts):
    listbox=self.widgets['varsList']

    notShow='lon+','lat+','pm','pn','time_step'
    vOrder='zeta','temp','salt','u','v','ubar','vbar','h'
    # vars will appear in this order is they exist...
    # then all the others

    Vars=[]
    self.variables={}
    vars,nc=netcdf.var(self.files[type])
    for k in vars.keys():
      if vars[k].ndim()>=2:
        Vars+=[k]
        self.variables[k]=self.files[type]
    nc.close()

    if type!='grid' and show_grid_vars:
      vars,nc=netcdf.var(self.files['grid'])
      for k in vars.keys():
        if vars[k].ndim()>=2 and k not in Vars:
          Vars+=[k]
          self.variables[k]=self.files['grid']
      nc.close()

    # remove some:
    for v in notShow:
      if v[-1]=='+':
        for vh in Vars[:]:
          if vh.find(v[:-1])==0:
            Vars.remove(vh)
            self.variables.pop(vh)

      elif v in Vars:
        Vars.remove(v)
        self.variables.pop(v)

    # order vars:
    vars=[]
    for i in vOrder:
      if i in Vars: vars+=[i]

    for i in Vars:
      if i not in vars: vars+=[i]

    # delete listbox entries:
    listbox.delete(0,tk.END)

    for v in vars:
        listbox.insert(tk.END, v)

    listbox.select_set(0)

    # store vinfo:
    self.swapp['vinfo']={}
    self.swapp['ftitle']={}
    for v in vars:
      f=self.variables[v]
      if not self.swapp['ftitle'].has_key(f):
        print f,v
        try:
          ftitle=netcdf.fatt(f,'title')
        except: ftitle='NO TITLE'
        self.swapp['ftitle'][f]=ftitle

      tmp=netcdf.vatt(f,v)
      atts={}
      for k in tmp.keys(): atts[k]=tmp[k].value
      self.swapp['vinfo'][v]=atts#netcdf.vatt(f,v)

    # show derived vars:
    Dlistbox=self.widgets['DvarsList']
    # delete listbox entries:
    Dlistbox.delete(0,tk.END)

    Dvars=[]
    DVars,DReq=get_derived()
    for k in DVars.keys():
      # required vars:
      if all([i in vars for i in DReq[k]]):
        Dvars+=['*'+k]
        self.swapp['vinfo']['*'+k]=DVars[k]


#####    if ('u' in vars and 'v' in vars) or ('ubar' in vars and 'vbar' in vars):
#####      Dvars+=['*speed','*ke']

    for v in Dvars:
        Dlistbox.insert(tk.END, v)
#####        self.swapp['vinfo'][v]={'long_name': 'some name','units':'some units'}
#####        print '='*50,'     -------> sacar isto de romsgui.derived !!'



  def clear_fig(self):
    if 0: # slower
      self.figure.clf()
      t0=pytime.time()
      self.axes = self.figure.add_axes((.12,.15,.8,.78))
      self.cbar=self.figure.add_axes((.1,.06,.8,.03))
    else:
      if hasattr(self,'axes'): self.axes.clear()
      else: self.axes = self.figure.add_axes((.1,.19,.85,.75))
      if hasattr(self,'cbar'): self.cbar.clear()
      else: self.cbar=self.figure.add_axes((.1,.1,.85,.03))

  def __get_varname(self,derived=False):
    if derived:
      vars=self.widgets['DvarsList']
    else:
      vars=self.widgets['varsList']

    i=vars.curselection()
    if i:
      i=int(i[0])
      return vars.get(i)
    else: return False

  def __set_varname(self,evt=False,varname=False,show=True):
    if evt:
      print evt.keysym
      if   evt.keysym=='t': varname='temp'
      elif evt.keysym=='s': varname='salt'
      elif evt.keysym=='u': varname='u'
      elif evt.keysym=='v': varname='v'
      elif evt.keysym=='z': varname='zeta'
      elif evt.keysym=='h': varname='h'
      elif evt.keysym=='U': varname='ubar'
      elif evt.keysym=='V': varname='vbar'

    if varname:
      self.select_var(varname)

    if show: self.__show_after_select_var()


  def __get_itime(self): return int(self.widgets['time'].get())

  def __get_zlevel(self): return float(self.widgets['zlev'].get())

  def __set_zlevel(self,evt=False,lev=False,show=True):
    if evt:
      if   evt.keysym=='B': lev=0
      elif evt.keysym=='S':
        var=self.__get_varname()
        ob=self.__get_romsobj(var)
        if 's_w' in netcdf.vdim(ob.name,var).keys():
          lev=ob.S_W-1
        else: lev=ob.S_RHO-1

    if not lev is False:
      z=self.widgets['zlev']
      z.delete(0,tk.END)
      z.insert(0,lev)

    if show: self.show()


  def __get_caxis(self):
    if self.swapp['caxis_mode']=='manual':
#    if self.widgets['cax_ck'].get():
      try: return [float(i) for i in self.widgets['cax'].get().split()]
      except: return False
    else:
      return False

  def __init_depths(self):#auto_hvals(self):
    h=self.files['Grid'].use('h')

    # set isobaths values:
    if  self.widgets['hvals'].get()=='auto':
      if h.min()!=h.max():
        tmp,stmp=ticks.loose_label(h.min(), h.max(),labels=1)
        if len(stmp)>2: s=stmp[1]+' '+stmp[2]
        else: s=stmp[0]+' '+stmp[1]
      else: s=''

      self.widgets['hvals'].delete(0,tk.END)
      self.widgets['hvals'].insert(0,s)

    # set depths for z slices:
    if h.min()!=h.max():
      zd=-ticks.tight(h.min(), h.max(),30)
    else: zd=[0]
    # want to have -10 at start:
    if -10. not in zd:
      zd=np.concatenate(([-10],zd))
      zd.sort()
      zd=zd[::-1]

    self.options['zslices_depths']=zd


  def __get_hvals(self):
    vals=self.widgets['hvals'].get()
    if vals.find(':')>=0:
      v=[float(i) for i in vals.split(':')]
      cvals=np.arange(v[0],v[1],v[2])
    else:
      cvals=[float(i) for i in vals.split()]

    return cvals


  def __set_caxis(self,cax):
    #if self.widgets['cax_ck'].get()==0:
    if self.swapp['caxis_mode']=='auto':
      self.widgets['cax'].delete(0,tk.END)
      self.widgets['cax'].insert(0, '%.2f %.2f' % (cax[0],cax[1]))


  def __get_vec(self):
    # dvec:
    dvec=self.widgets['dvec']
    d=int(dvec.get())

    # scale:
    svec=self.widgets['svec']
    s=float(svec.get())

    return d,s

  def select_var(self,varname):
    if not varname: return

    if varname[0]=='*':
      v=self.widgets['DvarsList']
    else:
      v=self.widgets['varsList']

    names=v.get(0,tk.END)
    if varname in names:
      v.selection_clear(0,tk.END)
      ind=list(names).index(varname)
      v.selection_set(ind)


  def __show_vinfo(self,varname=False):
    if not varname:  varname = self.__get_varname()
    ########f=self.variables[varname] will not work cos of the derived vars.... maybe fix this later !!! TODO
    f=self.files['his']
    

    info=self.swapp['vinfo'][varname]
    ftitle=self.swapp['ftitle'][f]

#    info=netcdf.vatt(f,varname)
#    ftitle=netcdf.fatt(f,'title')['value']

    # delete listbox entries:
    vinfo=self.widgets['vinfo']
    vinfo.delete(0,tk.END)

    vinfo.insert(tk.END, varname+' <'+ftitle+'>')
    for k in info.keys():
      #######print '======',varname, k, info[k]
      ########print '===>', k, info[k]####.value
      vinfo.insert(tk.END, k+': '+str(info[k]))#####.value))

  def __show_msg(self,msg,msgType=None):
    if   msgType=='error':   color='red'
    elif msgType=='warning': color='blue'
    else: color='black'

    w=self.widgets['msg']
    w.config(text=msg,fg=color)

  def __get_romsobj(self,varname):
    f=self.variables[varname]
    for k in self.files.keys():
      if self.files[k]==f: return self.files[k[0].upper()+k[1:]]

    return False

  def __start_proj(self):
    # current proj settings:
    lonlims,latlims=self.__get_grid_lims()
    resolution=self.options['coastline_res'].get()[0]
    if resolution=='n': resolution='l' # when using none

    curr_proj={}
    curr_proj['proj']=self.options['projection'].get()
    curr_proj['resolution']=resolution
    curr_proj['lonlims']=lonlims
    curr_proj['latlims']=latlims

    if self.swapp['last_proj']!=curr_proj:
      self.swapp['last_proj']=curr_proj

      if True:#curr_proj['proj']=='merc':
        map = Basemap(projection=curr_proj['proj'], lat_ts=0.0,
                      resolution=curr_proj['resolution'],
                      urcrnrlon=lonlims[1], urcrnrlat=latlims[1],
                      llcrnrlon=lonlims[0], llcrnrlat=latlims[0],
                      lon_0=0.5*(lonlims[0]+lonlims[1]),
                      lat_0=0.5*(latlims[0]+latlims[1]))
      else: map=False
      self.map=map

  def __get_slice(self,romsobj,varname,time,zlev,surf_nans,currents=False):
    curr_slice={}
    curr_slice['varname']   = varname
    curr_slice['romsobj']   = romsobj
    curr_slice['time']      = time
    curr_slice['zlev']      = zlev

    if zlev<0:   curr_slice['zslice_SN'] = surf_nans
    if currents: curr_slice['currents']  = True

    for i in range(len(self.swapp['last_slices'])):
      if self.swapp['last_slices'][i]==curr_slice:
         print 'USING SWAPP...'
         return self.swapp['last_slices_data'][i]

    # spherical:
    if romsobj.grid.use('spherical') in (0,'F'): spherical=False
    else: spherical=True

    # slice:
    if varname[0]=='*':
      x,y,z,v=romsobj.slice_derived(varname[1:],zlev,time)
      msg=''
    else:
      if zlev>=0:
        zlev=int(zlev)
        x,y,z,v,msg=romsobj.slicek(varname,zlev,time,msg=True)
      else:
        x,y,z,v,msg=romsobj.slicez(varname,zlev,time,msg=True,surf_nans=surf_nans)

    if currents:
      if varname in ('zeta','ubar','vbar'): vel_zlev='bar'
      else: vel_zlev=zlev
      x_vel,y_vel,z_vel,u_vel,v_vel=romsobj.sliceuv(vel_zlev,time)

      if spherical:
        # rotate to proj angles:
        u_vel,v_vel=self.map.rotate_vector(u_vel,v_vel,x_vel,y_vel)

      currents={'x':x_vel,'y':y_vel,'u':u_vel,'v':v_vel}

    if msg: # slice error!!
      return msg

    data=x,y,z,v, currents,spherical

    print 'STORING NEW SLICE...'
    # store slice:
    if len(self.swapp['last_slices'])<self.swapp['n_slices']:
      self.swapp['last_slices']+=[curr_slice]
      self.swapp['last_slices_data']+=[data]
    else:
      self.swapp['last_slices']=self.swapp['last_slices'][1:]+[curr_slice]
      self.swapp['last_slices_data']=self.swapp['last_slices_data'][1:]+[data]

    return data


  def __gen_title(self,romsobj,varname,zlev,time,currents):
    z_title=''
    t_title=''
    v_title=''
    if varname[0]=='*':  return '' # derived var

    # about var and currents:
    v_title=varname
    if currents:
      if varname in ('zeta','ubar','vbar'): v_title+='+bar currents'
      else: v_title+='+currents'

    # about depth:
    if romsobj.hasz(varname):
      if zlev>=0: z_title=' k='+str(zlev)+' '
      else: z_title=' z='+str(zlev)+'m '

    # about date:
    if romsobj.hast(varname):
      if not romsobj.datetime is False:
        try:
          t_title=' '+romsobj.datetime[time].isoformat(' ')
        except:
           # some old dates have not isoformat !!!
           t_title=' '+str(romsobj.datetime[time])
      else:
        try: yorig=int(self.widgets['yorig'].get())
        except: yorig=1
        ano,mes,dia,hora=greg2(romsobj.tdays[time],yorig)[:4]
        mes=int(mes)
        if yorig==1:
          date='%02d-%s %02dh' % (dia,month_names(mes),hora)
        else:
          date='%02d-%s-%04d %02dh' % (dia,month_names(mes),ano,hora)

        t_title=' tdays=%.3f %s ' % (romsobj.tdays[time],date)

    return v_title+z_title+t_title


  def __clear_swapp(self):
    self.swapp['last_slices']=[]
    self.swapp['last_slices_data']=[]

  def __draw_contourf(self):
     self.options['graf']='contourf'
     self.show()
  def __draw_contour(self):
     self.options['graf']='contour'
     self.show()
  def __draw_pcolor(self):
     self.options['graf']='pcolor'
     self.show()

  def __show_newfig(self):
    self.show(newfig=True)


  def __show_after_select_Dvar(self,evt=False):
    self.swapp['caxis_mode']='auto'
    self.show(derived=True)

  def __show_after_select_var(self,evt=False):
    self.swapp['caxis_mode']='auto'
    self.show()

  def __show_after_cax(self,evt=False):
    print 'setting caxis manual'
    self.swapp['caxis_mode']='manual'
    self.show()

  def show(self,evt=False,newfig=False,derived=False,**kargs):

    try: what=kargs['what']
    except: what=False

    if not self.files.has_key('his') and  self.files.has_key('grid'): what='grid'

    ######graf='contourf'
    ######if kargs.has_key('graf'): graf=kargs['graf']

    varname=''
    currents=False
    t0=pytime.time()

    if what !='grid':
      # varname, itime, zlevel:
      varname = self.__get_varname(derived=derived)

      if not varname:
        varname=self.swapp['last_varname']
        self.select_var(varname)

      if not varname: return
      self.swapp['last_varname']=varname

      #q=self.__get_romsobj(varname) isto j n pode ser usado devido as variaveis derived !!!
      q=self.files['His']

      if not q:
        self.__show_msg('roms obj not found %s' % varname,msgType='error')
        return

      time    = self.__get_itime()
      zlev    = self.__get_zlevel()
      print '0===',pytime.time()-t0

      t0=pytime.time()
      # show v info:
      self.__show_vinfo(varname)
      print '1 vinfo===',pytime.time()-t0
      t0=pytime.time()

      # surface nans options for zslice:
      surf_nans=int(self.options['zslice_SN'].get())

      # slice
      dvec,svec=self.__get_vec()
      if dvec!=0: currents=True
      data=self.__get_slice(q,varname,time,zlev,surf_nans,currents)

      # slice title:
      stitle=self.__gen_title(q,varname,zlev,time,currents)

      if isinstance(data,basestring): # error msg
        self.__show_msg(data,msgType='error')
        return
      else:
        self.__show_msg('')
        x,y,z,v,currents,spherical=data

      print '3===',pytime.time()-t0
      t0=pytime.time()

    else:

     try:
       print 'loading grid vars..........'
       x,y,v,mask=self.files['Grid'].vars()
       stitle=self.files['grid']
     except:
       self.__show_msg('Cannot use grid',msgType='error')
       return


    print 'get/set caxis'
    # get/set caxis:
    cax=self.__get_caxis()
    ntick=12
    ticksAuto=True
    if cax:
      try:    vmin,vmax=cax
      except:
        vmin,vmax,ntick=cax
        ticksAuto=False
    else:
      vmin,vmax=v.min(),v.max()
      self.__set_caxis((vmin,vmax))

    print '4===',pytime.time()-t0

    t0=pytime.time()
    # start plotting:
    pylab.ioff()
    print 'plotting...'
    if newfig:
      fig=pylab.figure()
      ax=pylab.axes()
      cbarOrientation='vertical'
      #cbar=fig.add_axes((.1,.06,.8,.03))#colorbar()
      #cbar=pylab.colorbar(orientation=cbarOrientation).ax
      cbar=None
    else:
      self.clear_fig()
      ax=self.axes
      cbar=self.cbar
      cbarOrientation='horizontal'

    print '4.5===',pytime.time()-t0
    t0=pytime.time()

    print 'here 0'
    # use projection:
    spherical=not self.files['Grid'].use('spherical') in (0,'F')

    if spherical:
      self.__start_proj()

      print 'here 1'

      print '5===',pytime.time()-t0
      t0=pytime.time()

      if self.options['coastline_res'].get()!='none':
        self.map.drawcoastlines(linewidth=.5,color='#999999',ax=ax)

      if self.options['draw_rivers'].get():
        self.map.drawrivers(linewidth=.5,color='blue',ax=ax)

      if self.options['draw_countries'].get():
        self.map.drawcountries(linewidth=.5,color='#999999',ax=ax)

      if self.options['draw_states'].get():
        self.map.drawstates(linewidth=.5,color='green',ax=ax)

      if self.options['fill_continents'].get():
        try: # may fail for high zooms
          self.map.fillcontinents(color='#ebebeb',ax=ax)
        except: pass

      scale=self.options['draw_scale'].get()
      scaleL=self.options['scale_length'].get()
      if scale!='none':
        lonlims,latlims=self.__get_grid_lims()
        if scaleL=='auto':
          grid_dist=calc.distance(np.array(lonlims),np.array(latlims))
          scaleL=ticks.nicenum(grid_dist[-1]/4000.,1)

        dlon=lonlims[1]-lonlims[0]
        dlat=latlims[1]-latlims[0]
        rx=5.
        ry=15.
        if scale=='NW':
          lon=lonlims[0]+dlon/rx
          lat=latlims[1]-dlat/ry
        if scale=='NE':
          lon=lonlims[1]-dlon/rx
          lat=latlims[1]-dlat/ry
        if scale=='SW':
          lon=lonlims[0]+dlon/rx
          lat=latlims[0]+dlat/ry
        if scale=='SE':
          lon=lonlims[1]-dlon/rx
          lat=latlims[0]+dlat/ry

        lon0,lat0=lon,lat
        self.map.drawmapscale(lon=lon,lat=lat,lon0=lon0,lat0=lat0,
                              length=scaleL,ax=ax,barstyle='fancy',
                              fontsize=7)

      meridians,paralels=self.__calc_xyticks()
      self.map.drawparallels(paralels,  labels=[1,0,0,0],linewidth=.5,
                        color='#4d4d4d',ax=ax,dashes=[1,4])
      self.map.drawmeridians(meridians, labels=[0,0,0,1],linewidth=.5,
                        color='#4d4d4d',ax=ax,dashes=[1,4])


    print '6===',pytime.time()-t0
    t0=pytime.time()

    # pcolor/contourf:
    if vmin==vmax:
      vals=np.array([vmin-1,vmin,vmin+1])
    elif ticksAuto:
      vals=ticks.loose_label_n(vmin,vmax,ntick=ntick)
    else:
      vals=np.linspace(vmin,vmax,ntick+1)

    cmap=self.options['colormap'].get()
    cmapr=self.options['colormap_r'].get()
    if cmapr: cmap=cmap+'_r'

    if cmap.startswith('_'): # usersdef
      cmap=get_ucmaps(cmap[1:])
    else:
      try:
        cmap=getattr(pylab.cm,cmap)
      except:
        cmap=getattr(basemap_cm,cmap)


    if 0: #varname=='h':
      from matplotlib.colors import rgb2hex
      def calc_colors(var,cmap=pylab.cm.jet,clim=False):
        mm=pylab.cm.ScalarMappable(cmap=cmap)
        if clim is False: clim=var.min(),var.max()
        mm.set_clim(clim)
        return mm.to_rgba(var)


      from mpl_toolkits.mplot3d import Axes3D
#      self.clear_fig()
      self.figure.clf()

      ax=Axes3D(self.figure)
      #handle=ax.plot_surface(x,y,v,rstride=5, cstride=5, cmap=pylab.cm.jet,linewidth=0, antialiased=False)
      if 1:
        handle=ax.plot_wireframe(x,y,-v,rstride=2, cstride=2,linewidth=.1)
        self.xx=x
        self.yy=y
        self.hh=v
      self.ax3=ax

      xx,yy,zz,vv=self.files['His'].slicez('temp',-2000,0)
      zz=np.ma.masked_where(vv.mask,zz)
#      handle=ax.plot_surface(xx,yy,zz,rstride=5, cstride=5,facecolors=vv)#,vmin=vmin,vmax=vmax)
      colors=calc_colors(vv[::3,::3])
      self.colors=colors
      print colors.shape, xx[::3].shape
      handle=ax.plot_surface(xx[::3,::3],yy[::3,::3],zz[::3,::3],rstride=1, cstride=1,facecolors=colors,
                             linewidth=0,shade=0,vmin=vmin,vmax=vmax)
      #handle=ax.plot_surface(xx[::5,::5],yy[::5,::5],zz[::5,::5],rstride=1, cstride=1,vmin=vmin,vmax=vmax,facecolors=colors)
      self.canvas.show()
      return
    else:
      if spherical: x,y=self.map(x,y)
      if self.options['graf']=='contourf':
        handle=ax.contourf(x,y,v,vals,cmap=cmap)
      elif self.options['graf']=='contour':
        handle=ax.contour(x,y,v,vals,cmap=cmap)
      elif self.options['graf']=='pcolor':
        handle=ax.pcolormesh(x,y,v,vmin=vmin,vmax=vmax,cmap=cmap)

      ax.set_title(stitle)


    print '7===',pytime.time()-t0
    t0=pytime.time()

    # vectors:
    if currents:
      #x,y,z,u,v=q.sliceuv(zlev,time)
      x=currents['x'][::dvec,::dvec]
      y=currents['y'][::dvec,::dvec]
      u=currents['u'][::dvec,::dvec]
      v=currents['v'][::dvec,::dvec]

      if spherical: xm,ym=self.map(x,y)
      else: xm,ym=x,y

      if svec==0: scale=None
      else:
        dx=xm.max()-xm.min()
        DX=x.max()-x.min()
        scale=1./svec*DX/dx

      qv=ax.quiver(xm,ym,u,v,units='x',scale=scale)
      #qv=self.map.quiver(xm,ym,u,v,units='x',scale=scale)
                   #headwidth=4,headlength=8,linewidth=.2
#TESTAR self.map.quiver !!!!!!!!!!

      print '======',qv.scale
      val=2*ticks.nicenum(np.sqrt(u**2+v**2).mean(),1)
      #qk = ax.quiverkey(qv, .05, -.1, val, str(val)+' m/s',labelpos='E',coordinates='axes')
      qk = ax.quiverkey(qv, .1, .15, val, str(val)+' m/s',labelpos='W',coordinates='figure')


  #    mostrar o valor com n+1 casas decimais, em q n=int(np.max((-np.floor(np.log10(svec)),0)))
  #    --> ver nicenum !!!

    if varname!='h':
      # colorbar:
      print 'adding cbar...',handle,cbar,cbarOrientation
      if self.options['graf']=='pcolor':cbarDrawEdges=False
      else: cbarDrawEdges=True
      cb=pylab.colorbar(handle,cax=cbar,orientation=cbarOrientation,drawedges=cbarDrawEdges)


      # bathy contour and border:
      hvals=self.__get_hvals()
      x,y=self.files['Grid'].lon,self.files['Grid'].lat
      if spherical: x,y=self.map(x,y)
      if hvals:
        ax.contour(x,y,self.files['Grid'].h,hvals,colors='k',linewidths=.5)

      if self.options['draw_grdborder'].get():
        xb,yb=self.files['Grid'].border()
        if spherical: xb,yb=self.map(xb,yb)
        ax.plot(xb,yb,lw=.5,color='#cccccc',ls='-')


      # cities: ----------------------------------------------
      if spherical and not cities_name is False:
        clon,clat=self.map(cities_lon,cities_lat)
        ax.plot(clon,clat,'b.')

        dx=(self.map.xmax-self.map.xmin)/100.
        dy=(self.map.ymax-self.map.ymin)/100.

        xp=np.array([self.map.xmin, self.map.xmax,self.map.xmax,self.map.xmin])
        yp=np.array([self.map.ymin, self.map.ymin,self.map.ymax,self.map.ymax])
        isin=calc.inpolygon(np.array(clon),np.array(clat),xp,yp)

        clon   = np.array(clon)[isin]
        clat   = np.array(clat)[isin]
        cnames = np.array(cities_name)[isin]

        if len(cnames) > 15:
          if 0:
            # remove some percentage of points:
            tmp=np.random.random(clon.shape)
            p=0.7
            tmp=tmp>p
            clon=clon[tmp]
            clat=clat[tmp]
            cnames=cnames[tmp]
          else:
            cnames=[] # show none !!

        if len(cnames):
          ax.plot(clon,clat,'bo')
          for i in range(len(cnames)):
            s=cnames[i].decode('UTF-8')
            ax.text(clon[i]+dx,clat[i]+dy,s,fontsize=8)

      #  ----------------------------------------------


      print '8===',pytime.time()-t0
      t0=pytime.time()

      if spherical: ax.axis([self.map.xmin, self.map.xmax, self.map.ymin, self.map.ymax])

    print 'now will draw!'
    if newfig:
      fig.show()
    else:
      self.canvas.show()

    pylab.ion()

    print '9===',pytime.time()-t0

    # activate tools like zoom, slices, etc:
    if   self.active_tool == 'zoom':        self.interactive_zoom(chstate=False)
    elif self.active_tool == 'vslice':      self.interactive_vslice(chstate=False)
    elif self.active_tool == 'profile':     self.interactive_profile(chstate=False)
    elif self.active_tool == 'time_series': self.interactive_time_series(chstate=False)
    elif self.active_tool == 'hovmuller':   self.interactive_hovmuller(chstate=False)

    # put focus on the main window, ie, remove focus from entry widgets
    # to allow use binds, like time change +-
    self.root.focus()


  def savefig(self,fname=False):
    if not fname:
      import tkFileDialog
      Formats = (
      ('PNG','*.png'),
      ('Encapsulated PostScript','*.eps'),
      ('PostScript','*.ps'),
      ('PDF','*.pdf'))

      fname = tkFileDialog.asksaveasfilename(
      parent=self.root,filetypes=Formats ,
      title="Save figure as...",
      defaultextension='.png')

    if fname:
      self.figure.savefig(fname,papertype='a4')
      self.swapp['last_savename']=fname

  def __savefig(self):
    self.savefig(fname=self.swapp['last_savename'])

  def closegui(self):
    import tkMessageBox
    res=tkMessageBox.askokcancel(title='close romsgui?',
                             message='are you sure you want to close?')

    if res:
      #self.root.destroy()
      self.root.quit()#destroy()
      global rGUIS
      print rGUIS
      rGUIS.remove(self.root)
      if len(rGUIS)==1:
        self.root0.destroy()
        self.root0=False
        rGUIS=[]

      print rGUIS

  def __about(self):
    import tkMessageBox
    msg='''
    Python romsgui
    Visualization interface for ocean model ROMS

    author: mma
    m.martalmeida@gmail.com
    '''
    tkMessageBox.showinfo('python romsgui',msg)

  def __homepage(self):
    import webbrowser
    webbrowser.open(self.url)


#  def draw_seg_OLDDD(self):
#    if not self.files.has_key('His'): return
#
#    try:
#      self.sliceline0=self.sliceline
#    except: pass
#
#    line, = self.axes.plot([0], [0])
#    self.axes.axis([self.map.xmin, self.map.xmax, self.map.ymin, self.map.ymax])
#    #self.__set_cursor('tcross')
#    self.__set_cursor('pencil')
#
#    try: prev=self.sliceline
#    except: pass
#
#    if 1:
#      from devel import line_builder
#      self.sliceline=self.linebuilder= line_builder.InteractiveLine(ax=self.axes,
#                                       cmd=self.__show_slicell,
#                                       axis=(self.map.xmin, self.map.xmax, self.map.ymin, self.map.ymax),
#                                       nmax=2)
#
#
#    else:
#      self.linebuilder= LineBuilder(line,self.__show_slicell)
#

  def interactive_zoom(self,chstate=True):
    '''chstate to be used with clicking with mouse on button zoom
    chstate not use when calling after show, to continue zooming
    '''
    if not self.files.has_key('His'): return

    ob=self.widgets['zoom']
    relief=ob['relief']

    if chstate and relief=='sunken':
      ob.config(relief='raised')
      self.active_tool=None
      self.zoom.cmd=False
      self.zoom.stop(False) # stop the zoom events that may be active
      self.__set_cursor()

    else:
      ob.config(relief='sunken')

      self.__set_cursor('tcross')

      self.zoom=pl_tools.InteractiveRect(ax=self.axes,axis=(self.map.xmin, self.map.xmax, self.map.ymin, self.map.ymax),
                                             cmd=self.__zoom)

      self.active_tool='zoom'

  def __zoom(self):
    ob=self.zoom # PUT THIS IN CACHE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    x,y=self.map(ob.x,ob.y,inverse=1)

    x0,x1=np.min(x),np.max(x)
    y0,y1=np.min(y),np.max(y)

    # check if zoom is too small... (maybe it was almost a click!)
    xl,yl=self._rgui__get_grid_lims()
    if (x1-x0)<(xl[1]-xl[0])/20. or (y1-y0)<(yl[1]-yl[0])/20. :
      if ob.button==3: zoomD='off'
      elif ob.button==1: zoomD='on'
      self.ch_axis(xy0=[x[0],y[0]],zoom=zoomD)
    else:
      lims=x0,x1,y0,y1
      self.ch_axis(axis=lims)



###  def draw_npoints(self,npoints,cmd=None):
###    if not self.files.has_key('His'): return
###
###    obj=pl_tools.InteractiveLine(ax=self.axes,cmd=cmd,
###                        axis=(self.map.xmin, self.map.xmax, self.map.ymin, self.map.ymax),
###                        npoints=npoints)
###
###    return obj


  def interactive_vslice(self,chstate=True):      self.__interactive_tool('vslice',chstate)
  def interactive_profile(self,chstate=True):     self.__interactive_tool('profile',chstate)
  def interactive_time_series(self,chstate=True): self.__interactive_tool('time_series',chstate)
  def interactive_hovmuller(self,chstate=True):   self.__interactive_tool('hovmuller',chstate)

  def __interactive_tool(self,tool,chstate):
    if not self.files.has_key('His'): return

    # stop any active tool:
    if self.active_tool and self.active_tool!=tool:  eval('self.interactive_'+self.active_tool)()

    # init tool cache data:
    if not self.cache.has_key(tool):
      self.cache[tool]={'obj':False}

    # button relief:
    button=self.widgets[tool]
    relief=button['relief']

    if chstate and relief=='sunken':
      button.config(relief='raised')

      # reset tool:
      self.active_tool=None

      # stop the events that may be active
      if self.cache[tool]['obj']: self.cache[tool]['obj'].stop(False)

      # reset cursor:
      self.__set_cursor()

    else:
      button.config(relief='sunken')

      # number of segment points:
      if tool in ('vslice','hovmuller'):
        npts=self.options['vslice_nlpts'].get()
        if npts=='auto': npts=-1
        else: npts=int(npts)
      else:
        npts=1

      # change cursor:
      self.__set_cursor('pencil')

      # start interactive tool:
      cmd=eval('self._rgui__'+tool)
      try:
        axis=self.map.xmin, self.map.xmax, self.map.ymin, self.map.ymax
      except:
        axis=self.axes.axis()

      self.cache[tool]['obj']=pl_tools.InteractiveLine(ax=self.axes,cmd=cmd,
                        axis=axis,
                        nmax=npts,type='blin')

      # set active tool
      self.active_tool=tool


  def __profile(self,prev=False): self.__show_profile('profile')

  def __time_series(self,prev=False): self.__show_profile('time_series')

  def __hovmuller(self,prev=False): self.__show_profile('hovmuller')


  def __show_profile(self,type):
    ob=self.cache[type]['obj']

    if not ob: return
    try:
      x,y=self.map(ob.x,ob.y,inverse=1)
    except: x,y=ob.x,ob.y

    # varname and time:
    var  = self.__get_varname()
    if not var:
      var=self.swapp['last_varname']
      self.select_var(var)

    q=self.__get_romsobj(var)
    itime = self.__get_itime()

    if type in ('profile','time_series'):
      x=x[0]
      y=y[0]

    if type=='profile':
        x,y,z,v=q.slicell(var,np.array(x),np.array(y),itime,dist=0)
    elif type=='time_series':
        zlev    = self.__get_zlevel()
        print var,x,y,zlev
        t,z,v=q.time_series(var,x,y,depth=zlev)


    # draw:
    if v.ndim>0:
      pylab.ioff()
      fig=pylab.figure()
    else: fig=False

    if v.ndim==1:
      if type=='profile':
        print 'HERE 4 prof'
        pylab.plot(v,z)
      elif type=='time_series':
        print 'HERE 4 ts'
        import datetime as dt
        try: yorig=int(self.widgets['yorig'].get())
        except: yorig=1
        
        dates=[dt.datetime(yorig,1,1)+dt.timedelta(days=i) for i in q.tdays]
        dates=pylab.date2num(dates)
#
        pylab.plot_date(dates,v,'-')
        pylab.grid()
#        print dates
#        pylab.plot(v)
        #pylab.plot(t,v) ## TODO plotar tempo como datas...
        #pylab.xlabel(var)

    elif v.ndim==2:
      pylab.pcolormesh(v)
      pylab.xlabel(var)
    else: print 'V=',v

    print 'HERE 5 fig=',fig
    if fig:
      fig.show()

    self.__set_cursor()
    self.select_var(var)
    pylab.ion()



  def __vslice(self,prev=False):
    if not prev:
      # vslice line:
      try:
        ob=self.cache['vslice']['obj']
      except:
        print 'cannot file vslice obj'
        return

      # covert to lon lat... after zoom the values change!
      try:
        ob.lon,ob.lat=self.map(ob.x,ob.y,inverse=1)
      except:
        ob.lon,ob.lat=ob.x,ob.y

      # store previous slice:
      self.cache['vslice']['obj_prev']=ob

    else:
      try:
        ob=self.cache['vslice']['obj_prev']
      except:
        print 'cannot file previous vslice obj'
        return

      # plot slice line:
      self.axes.plot(ob.x,ob.y,'k--x')
      self.canvas.draw()


    x,y=ob.lon,ob.lat
    if len(x)<2:
      print 'at least 2 points needed for vslice'

      # keep vslice active:
      if self.active_tool=='vslice':
        self.interactive_vslice(chstate=False)

      return

    # varname and time:
    var  = self.__get_varname()
    if not var:
      var=self.swapp['last_varname']
      self.select_var(var)

    q=self.__get_romsobj(var)
    itime = self.__get_itime()

    # number of points to use:
    npts=self.options['vslice_npts'].get()
    if npts=='auto':
      xm,ym=q.grid.vars(var)[:2]
      if len(x)==2:
        X,Y=calc.mcross_points(xm,ym,x,y) ######################--> talvez isto possa ser feito dentro de roms !?
      elif len(x)>2:
        X=np.array([])
        Y=np.array([])
        for i in range(len(x)-1):
          X_,Y_=calc.mcross_points(xm,ym,x[i:i+2],y[i:i+2])
          X=np.append(X,X_[:-1])
          Y=np.append(Y,Y_[:-1])

        X=np.append(X,X_[-1])
        Y=np.append(Y,Y_[-1])

    else:
      npts=int(npts)
      n=npts/len(x)
      X,Y=[],[]
      for i in range(len(x)-1):
        if i<len(x)-2: endpt=False
        else: endpt=True
        X.extend(np.linspace(x[i],x[i+1],n,endpoint=endpt).tolist())
        Y.extend(np.linspace(y[i],y[i+1],n,endpoint=endpt).tolist())

      X=np.array(X)
      Y=np.array(Y)


    # slicell:
    try:
      self.map.xmin
      d,z,v=q.slicell(var,X,Y,itime,dist=True)
    except:
      x,y,z,v=q.slicell(var,X,Y,itime,dist=False)
      x=x[0,...]
      y=y[0,...]
      #print x.shape, y.shape
      x=np.hstack((0,np.diff(x))).cumsum()
      y=np.hstack((0,np.diff(y))).cumsum()
      d=np.sqrt(x**2+y**2)

    # draw:
    pylab.ioff()
    f=pylab.figure()
    if v.ndim==2:
      pylab.pcolor(d/1000.,z,v)
      pylab.colorbar()
      pylab.ylabel('Depth [m]')
    elif v.ndim==1:
      pylab.plot(np.squeeze(d)/1000.,np.squeeze(v))

    pylab.xlabel('Position along the section [km]')

    if 0: # do not remove line !!
      try:
        self.sliceline0.remove_line()
      except: pass

    f.show()
    pylab.ion()

    self.__set_cursor()
    self.select_var(var)


    # keep vslice active:
    if self.active_tool=='vslice':
      self.interactive_vslice(chstate=False)


  def __set_cursor(self,type=''):
    self.canvas.get_tk_widget().config(cursor=type)


#class LineBuilder:
#    def __init__(self, line,cmd):
#        self.line = line
#        self.cmd=cmd
#        self.xs = list(line.get_xdata())
#        self.ys = list(line.get_ydata())
#        self.cid = line.figure.canvas.mpl_connect('button_press_event', self)
#
#    def __call__(self, event):
#        if event.button==1:
#          if event.inaxes!=self.line.axes: return
#          print event.xdata,event.ydata
#          if self.xs==[0] and self.ys==[0]:
#            self.xs=[event.xdata]
#            self.ys=[event.ydata]
#          else:
#            self.xs.append(event.xdata)
#            self.ys.append(event.ydata)
#
#          self.line.set_data(self.xs, self.ys)
#          self.line.figure.canvas.draw()
#        elif  event.button==3:
#          self.line.figure.canvas.mpl_disconnect(self.cid)
#          return self.cmd()
#'''
#fig = plt.figure()
#ax = fig.add_subplot(111)
#ax.set_title('click to build line segments')
#line, = ax.plot([0], [0])  # empty line
#linebuilder = LineBuilder(line)
#'''
#
#
#
#if __name__=='__main__':
#  r=rgui()
#  r.root.mainloop()
#
#
#
