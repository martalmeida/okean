import Tkinter as tk


class bButton():
  def __init__(self,master,image,state=0,toggle=True,**kargs):
    self.w=tk.Label(master,**kargs)
    if not toggle: state=0
    self.toggle=toggle
    self.state=state
    self.command=None

    import os
    p=os.path.dirname(os.path.abspath(__file__))

    self.icon_on=tk.PhotoImage(file=os.path.join(p,'icons',image+'_on.gif'))
    self.icon_off=tk.PhotoImage(file=os.path.join(p,'icons',image+'_off.gif'))
    self.icon_over=tk.PhotoImage(file=os.path.join(p,'icons',image+'_over.gif'))
    if state: s='on'
    else: s='off'
    icon=eval('self.icon_%s'%s)
    self.w.config(image=icon,width="16",height="16")

    self.w.bind("<Enter>", self.mouse_over)
    self.w.bind("<Leave>", self.mouse_out)
    self.w.bind("<Button-1>", self.mouse_click)
    self.w.bind("<ButtonRelease-1>", self.mouse_over)

  def mouse_over(self,evt):
    self.w.config(image=self.icon_over,cursor='hand1')

  def mouse_out(self,evt):
    if self.toggle:
      if self.state:
        self.w.config(image=self.icon_on)
      else:
        self.w.config(image=self.icon_off)
    else:
      self.w.config(image=self.icon_off)


  def mouse_click(self,evt):
    self.w.config(image=self.icon_on)
    self.state=int(not self.state)
    #if self.command:
    self.command()

def callback():
    print("click!")

def rgui_tools(master):
  out={}

  h=tk.Label(master,width=40,height=120)#,bg='green')#,**kargs)
  out['base']=h

  b = bButton(h, state=0,image='zoom',toggle=1)#, command=callback)
  b.w.grid(row=0, column=0)
  out['zoom']=b

  b = bButton(h, state=0,image='nozoom',toggle=0)#, command=callback)
  b.w.grid(row=0, column=1)
  out['nozoom']=b

  b = bButton(h, state=0,image='path',toggle=0)#, command=callback)
  b.w.grid(row=1, column=0)
  out['path']=b

  b = bButton(h, state=0,image='pathm',toggle=0)#, command=callback)
  b.w.grid(row=1, column=1)
  out['pathm']=b

  b = bButton(h, state=0,image='prof',toggle=1)#, command=callback)
  b.w.grid(row=2, column=0)
  out['prof']=b

  b = bButton(h, state=0,image='hov',toggle=0)#, command=callback)
  b.w.grid(row=2, column=1)
  out['hov']=b

  b = bButton(h, state=0,image='ts',toggle=0)#, command=callback)
  b.w.grid(row=3, column=0)
  out['ts']=b

  return out

#master = tk.Tk()
#o=rgui_tools(master)
#o['base'].pack()
#tk.mainloop()
