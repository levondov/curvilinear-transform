import wx
import sys
import threading
import time
import numpy as np
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as Navbar
from matplotlib.backend_bases import NavigationToolbar2 as Nb
from coord_transform import coord_transform

class ctGUI:

    # wx stuff
    ct = 0
    size = (1300,550)
    frame = 0
    panel = 0

    # buttons, checkboxes, and other
    btn1 = 0
    btn2 = 0
    btn3 = 0
    checkbox1 = 0
    checkbox2 = 0
    # checkboxes for categories of points
    precheckgp = 0
    precheckst = 0
    prechecktp = 0
    prechecktsp = 0
    precheckother = 0
    postcheckgp = 0
    postcheckst = 0
    postchecktp = 0
    postchecktsp = 0
    postcheckother = 0

    #matplotlib stuff
    fig1 = 0
    fig2 = 0
    ax1 = 0
    ax2 = 0
    canvas1 = 0
    canvas2 = 0
    toolbar1 = 0
    toolbar2 = 0
    preColor = False
    postColor = False
    prePoints = [True,True,True,True,True]
    postPoints = [True,True,True,True,True]

    # data
    xdata = 0
    xdataO = 0 #xdata, but organized, i.e. the categories put in order -> [gp points, tp points, sp points... etc]
    ydata = 0
    ydataO = 0
    zdata = 0
    zdataO = 0
    locations = ''

    # coord transform class variables
    centerlinegp = 10 #default line used as centerline is grid point 10m locations
    gpCord = 0
    stCord = 0
    tpCord = 0
    tspCord = 0
    otherCord = 0
    # coord for all data point transform
    s = 0
    n = 0
    l = 0
    # coord for each category
    s1 = 0; n1 = 0; l1 = 0;
    s2 = 0; n2 = 0; l2 = 0;
    s3 = 0; n3 = 0; l3 = 0;
    s4 = 0; n4 = 0; l4 = 0;
    s5 = 0; n5 = 0; l5 = 0;


    def __init__(self):

        # start wx GUI
        self.app = wx.App()

        # init frame
        self.frame = wx.Frame(None, -1, 'coord transform v1.2', style=wx.DEFAULT_FRAME_STYLE ^ wx.RESIZE_BORDER);
        self.frame.SetSize(self.size);
        self.frame.Centre(); # frame in the center of the screen

        # init panel
        self.panel = wx.Panel(self.frame, -1);

        # initialize our coord transform program
        self.xdata = np.genfromtxt("corrected_grid_sensor_ptsX_20131007.txt", skiprows=2, dtype=None)
        self.ydata = np.genfromtxt("corrected_grid_sensor_ptsY_20131007.txt", skiprows=2, dtype=None)
        self.zdata = np.genfromtxt("corrected_grid_sensor_ptsZ_20131007.txt", skiprows=2, dtype=None)
        self.locations = np.genfromtxt("corrected_grid_sensor_locations_20131007.txt", skiprows=2, dtype=None)

        # create object
        self.ct = coord_transform(self.xdata, self.ydata, self.zdata, self.locations) #give the object all the coordinate/location data

        # organize everything into their categories
        self.ct.organize()

        # get coordinate for all different type of points
        self.gpCord = self.ct.getgp()
        self.stCord = self.ct.getst()
        self.tspCord = self.ct.gettsp()
        self.tpCord = self.ct.gettp()
        self.otherCord = self.ct.getother()

        #decide which grid point lines to use as center, 0,5,10,15, or 20 meters
        self.ct.setCenterLine(self.centerlinegp)

        #create the centerline points through linear interpolation
        self.ct.createCenterline()

        # organize the data we have now and do the transformation
        self.transform()

        # create GUI stuff
        self.initialize()

        self.frame.Show(True) #show the GUI on the screen

        # main loop for wx GUI
        self.app.MainLoop()

    def initialize(self):

        # create canvas and toolbar objects
        self.fig1 = plt.figure()
        plt.grid(True)
        plt.title('Pre Transformation'); plt.xlabel('meters'); plt.ylabel('meters');
        self.fig2 = plt.figure()
        plt.grid(True)
        plt.title('Post Transformation'); plt.xlabel('meters'); plt.ylabel('meters');
        self.ax1 = self.fig1.add_subplot(111)#, projection='3d')
        self.ax2 = self.fig2.add_subplot(111)#, projection='3d')
        self.canvas1 = FigureCanvas(self.panel, -1, self.fig1)  #use matplotlib backend to grab canvas
        self.canvas2 = FigureCanvas(self.panel, -1, self.fig2)
        self.toolbar1 = Navbar(self.canvas1) #use matplotlib backend to grab toolbar
        self.toolbar2 = Navbar(self.canvas2)


        # default settings (upon launch) is just transform all data points and use a single color
        self.prePlot(self.preColor,self.prePoints)
        self.postPlot(self.postColor,self.postPoints)

        # create sizers
        vbox = wx.BoxSizer(wx.VERTICAL)
        hbox = wx.BoxSizer(wx.HORIZONTAL)
        hbox1 = wx.BoxSizer(wx.HORIZONTAL)
        hbox2 = wx.BoxSizer(wx.HORIZONTAL)

        hbox.Add(self.toolbar1,border=5)

        self.btn1 = wx.Button(self.panel, label='launch 3D', size=(80,30))
        hbox.Add(self.btn1, flag=wx.CENTER, border=5)

        self.checkbox1 = wx.CheckBox(self.panel, label='Display Colors')
        hbox.Add(self.checkbox1, flag=wx.CENTER, border=5)

        hbox.Add((315,-1))
        hbox.Add(self.toolbar2, border=5)

        self.btn2 = wx.Button(self.panel, label='launch 3D', size=(80,30))
        hbox.Add(self.btn2, flag=wx.CENTER, border=5)

        self.checkbox2 = wx.CheckBox(self.panel, label='Display Colors')
        hbox.Add(self.checkbox2, flag=wx.CENTER, border=5)

        hbox1.Add(self.canvas1, border=5)
        hbox1.Add((20,-1))
        hbox1.Add(self.canvas2, border=5)

        line = wx.StaticLine(self.panel)

        self.precheckgp = wx.CheckBox(self.panel, label='grid points')
        self.precheckgp.SetValue(True)
        hbox2.Add(self.precheckgp, flag=wx.CENTER,border=5)
        hbox2.Add((5,-1))
        self.precheckst = wx.CheckBox(self.panel, label='soil tubes')
        self.precheckst.SetValue(True)
        hbox2.Add(self.precheckst, flag=wx.CENTER,border=5)
        hbox2.Add((5,-1))
        self.prechecktp = wx.CheckBox(self.panel, label='temperature probes')
        self.prechecktp.SetValue(True)
        hbox2.Add(self.prechecktp, flag=wx.CENTER,border=5)
        hbox2.Add((5,-1))
        self.prechecktsp = wx.CheckBox(self.panel, label='total station points')
        self.prechecktsp.SetValue(True)
        hbox2.Add(self.prechecktsp, flag=wx.CENTER,border=5)
        hbox2.Add((5,-1))
        self.precheckother = wx.CheckBox(self.panel, label='other')
        self.precheckother.SetValue(True)
        hbox2.Add(self.precheckother, flag=wx.CENTER,border=5)
        hbox2.Add((215,-1))

        self.postcheckgp = wx.CheckBox(self.panel, label='grid points')
        self.postcheckgp.SetValue(True)
        hbox2.Add(self.postcheckgp, flag=wx.CENTER,border=5)
        hbox2.Add((5,-1))
        self.postcheckst = wx.CheckBox(self.panel, label='soil tubes')
        self.postcheckst.SetValue(True)
        hbox2.Add(self.postcheckst, flag=wx.CENTER,border=5)
        hbox2.Add((5,-1))
        self.postchecktp = wx.CheckBox(self.panel, label='temperature probes')
        self.postchecktp.SetValue(True)
        hbox2.Add(self.postchecktp, flag=wx.CENTER,border=5)
        hbox2.Add((5,-1))
        self.postchecktsp = wx.CheckBox(self.panel, label='total station points')
        self.postchecktsp.SetValue(True)
        hbox2.Add(self.postchecktsp, flag=wx.CENTER,border=5)
        hbox2.Add((5,-1))
        self.postcheckother = wx.CheckBox(self.panel, label='other')
        self.postcheckother.SetValue(True)
        hbox2.Add(self.postcheckother, flag=wx.CENTER,border=5)

        vbox.Add(hbox,border=5)
        vbox.Add(line,flag=wx.EXPAND,border=5)
        vbox.Add(hbox2,border=5)
        vbox.Add(hbox1,border=5)


        # event handlers
        self.btn1.Bind(wx.EVT_BUTTON, self.plot3dpre)
        self.btn2.Bind(wx.EVT_BUTTON, self.plot3dpost)
        self.checkbox1.Bind(wx.EVT_CHECKBOX, self.preChangeColor)
        self.checkbox2.Bind(wx.EVT_CHECKBOX, self.postChangeColor)
        self.precheckgp.Bind(wx.EVT_CHECKBOX, self.preChangePoints)
        self.precheckst.Bind(wx.EVT_CHECKBOX, self.preChangePoints)
        self.prechecktp.Bind(wx.EVT_CHECKBOX, self.preChangePoints)
        self.prechecktsp.Bind(wx.EVT_CHECKBOX, self.preChangePoints)
        self.precheckother.Bind(wx.EVT_CHECKBOX, self.preChangePoints)
        self.postcheckgp.Bind(wx.EVT_CHECKBOX, self.postChangePoints)
        self.postcheckst.Bind(wx.EVT_CHECKBOX, self.postChangePoints)
        self.postchecktp.Bind(wx.EVT_CHECKBOX, self.postChangePoints)
        self.postchecktsp.Bind(wx.EVT_CHECKBOX, self.postChangePoints)
        self.postcheckother.Bind(wx.EVT_CHECKBOX, self.postChangePoints)


        self.panel.SetSizer(vbox,wx.CENTER)

    def transform(self):

        # put data points in order before transforming
        # the order will be gpPoints -> stPoints -> tpPoints -> tspPoints -> otherPoints
        self.xdataO = np.append(self.gpCord[:,0],self.stCord[:,0])
        self.xdataO = np.append(self.xdataO,self.tpCord[:,0])
        self.xdataO = np.append(self.xdataO,self.tspCord[:,0])
        self.xdataO = np.append(self.xdataO,self.otherCord[:,0])
        self.ydataO = np.append(self.gpCord[:,1],self.stCord[:,1])
        self.ydataO = np.append(self.ydataO,self.tpCord[:,1])
        self.ydataO = np.append(self.ydataO,self.tspCord[:,1])
        self.ydataO = np.append(self.ydataO,self.otherCord[:,1])
        self.zdataO = np.append(self.gpCord[:,2],self.stCord[:,2])
        self.zdataO = np.append(self.zdataO,self.tpCord[:,2])
        self.zdataO = np.append(self.zdataO,self.tspCord[:,2])
        self.zdataO = np.append(self.zdataO,self.otherCord[:,2])

        # transform all data points put in order
        self.s,self.n,self.l = self.ct.transform(self.xdataO.astype(np.float),self.ydataO.astype(np.float))

        # separate the coordinates for each category after transformation
        self.s1 = self.s[0:len(self.gpCord[:,0])]
        self.n1 = self.n[0:len(self.gpCord[:,0])]
        self.s2 = self.s[len(self.s1):len(self.s1)+len(self.stCord[:,0])]
        self.n2 = self.n[len(self.s1):len(self.s1)+len(self.stCord[:,0])]
        self.s3 = self.s[len(self.s1)+len(self.s2):len(self.s1)+len(self.s2)+len(self.tpCord[:,0])]
        self.n3 = self.n[len(self.s1)+len(self.s2):len(self.s1)+len(self.s2)+len(self.tpCord[:,0])]
        self.s4 = self.s[len(self.s1)+len(self.s2)+len(self.s3):len(self.s1)+len(self.s2)+len(self.s3)+len(self.tspCord[:,0])]
        self.n4 = self.n[len(self.s1)+len(self.s2)+len(self.s3):len(self.s1)+len(self.s2)+len(self.s3)+len(self.tspCord[:,0])]
        self.s5 = self.s[len(self.s1)+len(self.s2)+len(self.s3)+len(self.s4):len(self.s1)+len(self.s2)+len(self.s3)+len(self.s4)+len(self.otherCord[:,0])]
        self.n5 = self.n[len(self.s1)+len(self.s2)+len(self.s3)+len(self.s4):len(self.s1)+len(self.s2)+len(self.s3)+len(self.s4)+len(self.otherCord[:,0])]

        # Now we can plot only certain categories only,
        # e.g. if you only want to see soil tubes, simply plot s2,n2 points

    def preChangeColor(self,event):

        if self.checkbox1.IsChecked(): #color on
            self.ax1.cla()
            self.preColor = True
            self.prePlot(self.preColor,self.prePoints)
            self.ax1.legend()
            self.ax1.set_title('Pre Transformation'); self.ax1.set_xlabel('meters'); self.ax1.set_ylabel('meters');
            self.canvas1.draw()
        else: # color off
            self.ax1.cla()
            self.preColor = False
            self.prePlot(self.preColor,self.prePoints)
            self.ax1.legend_ = None # turn off legend
            self.ax1.set_title('Pre Transformation'); self.ax1.set_xlabel('meters'); self.ax1.set_ylabel('meters');
            self.canvas1.draw()

    def postChangeColor(self,event):

        if self.checkbox2.IsChecked(): #color on
            self.ax2.cla()
            self.postColor = True
            self.postPlot(self.postColor,self.postPoints)
            self.ax2.legend()
            self.ax2.set_title('Post Transformation'); self.ax2.set_xlabel('meters'); self.ax2.set_ylabel('meters');
            self.canvas2.draw()
        else: # color off
            self.ax2.cla()
            self.postColor = False
            self.postPlot(self.postColor,self.postPoints)
            self.ax2.legend_ = None # turn off legend
            self.ax2.set_title('Post Transformation'); self.ax2.set_xlabel('meters'); self.ax2.set_ylabel('meters');
            self.canvas2.draw()

    def preChangePoints(self,event):
        checks = [self.precheckgp,self.precheckst,self.prechecktp,self.prechecktsp,self.precheckother]
        self.prePoints = [False,False,False,False,False]
        for i in range(0,5):
            if checks[i].IsChecked():
                self.prePoints[i] = True

        self.preChangeColor(event)

    def postChangePoints(self,event):
        checks = [self.postcheckgp,self.postcheckst,self.postchecktp,self.postchecktsp,self.postcheckother]
        self.postPoints = [False,False,False,False,False]
        for i in range(0,5):
            if checks[i].IsChecked():
                self.postPoints[i] = True

        self.postChangeColor(event)

    def prePlot(self,color,points):
        self.ax1.grid(True)
        if color: #if color
            # check each category individually
            if points[0]:
                self.ax1.scatter(self.gpCord[:,0].astype(np.float),self.gpCord[:,1].astype(np.float),color='k',label='grid points')
            if points[1]:
                self.ax1.scatter(self.stCord[:,0].astype(np.float),self.stCord[:,1].astype(np.float),color='b',label='soil tubes')
            if points[2]:
                self.ax1.scatter(self.tpCord[:,0].astype(np.float),self.tpCord[:,1].astype(np.float),color='r',label='temperature probes')
            if points[3]:
                self.ax1.scatter(self.tspCord[:,0].astype(np.float),self.tspCord[:,1].astype(np.float),color='g',label='total station points')
            if points[4]:
                self.ax1.scatter(self.otherCord[:,0].astype(np.float),self.otherCord[:,1].astype(np.float),color='c',label='other')
        else: # no color
            # still check each category
            if points[0]:
                self.ax1.scatter(self.gpCord[:,0].astype(np.float),self.gpCord[:,1].astype(np.float))
            if points[1]:
                self.ax1.scatter(self.stCord[:,0].astype(np.float),self.stCord[:,1].astype(np.float))
            if points[2]:
                self.ax1.scatter(self.tpCord[:,0].astype(np.float),self.tpCord[:,1].astype(np.float))
            if points[3]:
                self.ax1.scatter(self.tspCord[:,0].astype(np.float),self.tspCord[:,1].astype(np.float))
            if points[4]:
                self.ax1.scatter(self.otherCord[:,0].astype(np.float),self.otherCord[:,1].astype(np.float))
    def postPlot(self,color,points):
        self.ax2.grid(True)
        if color: #if color
            # check each category individually
            if points[0]:
                self.ax2.scatter(self.s1*self.l,self.n1,color='k',label='grid points')
            if points[1]:
                self.ax2.scatter(self.s2*self.l,self.n2,color='b',label='soil tubes')
            if points[2]:
                self.ax2.scatter(self.s3*self.l,self.n3,color='r',label='temperature probes')
            if points[3]:
                self.ax2.scatter(self.s4*self.l,self.n4,color='g',label='total station points')
            if points[4]:
                self.ax2.scatter(self.s5*self.l,self.n5,color='c',label='other')
        else: # no color
            # still check each category individually
            if points[0]:
                self.ax2.scatter(self.s1*self.l,self.n1)
            if points[1]:
                self.ax2.scatter(self.s2*self.l,self.n2)
            if points[2]:
                self.ax2.scatter(self.s3*self.l,self.n3)
            if points[3]:
                self.ax2.scatter(self.s4*self.l,self.n4)
            if points[4]:
                self.ax2.scatter(self.s5*self.l,self.n5)

    # show the current plots in 3d in a separate window
    def plot3dpre(self,event):
        fig = plt.figure()
        ax = fig.add_subplot((111), projection='3d')
        ax.scatter(self.xdata,self.ydata,self.zdata)
        plt.show()

    def plot3dpost(self,event):
        fig = plt.figure()
        ax = fig.add_subplot((111), projection='3d')
        ax.scatter(self.s*self.l,self.n,self.zdataO.astype(np.float))
        plt.show()