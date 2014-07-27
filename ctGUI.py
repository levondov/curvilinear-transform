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
    ydata = 0
    zdata = 0
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
        self.frame = wx.Frame(None, -1, 'coord transform v1', style=wx.DEFAULT_FRAME_STYLE ^ wx.RESIZE_BORDER);
        self.frame.SetSize(self.size);
        self.frame.Centre(); # frame in the center of the screen

        # init panel
        self.panel = wx.Panel(self.frame, -1);

        # initialize our coord transform program
        self.xdata = np.genfromtxt("corrected_grid_sensor_ptsX_20131007.txt", skiprows=2, dtype=None)
        self.ydata = np.genfromtxt("corrected_grid_sensor_ptsY_20131007.txt", skiprows=2, dtype=None)
        self.zdata = np.genfromtxt("corrected_grid_sensor_ptsZ_20131007.txt", skiprows=2, dtype=None)
        self.locations = np.genfromtxt("corrected_grid_sensor_locations_20131007.txt", skiprows=2, dtype=None)
        self.ct = coord_transform(self.xdata, self.ydata, self.zdata, self.locations) #give the object all the coordinate/location data
        self.ct.organize()
        #decide which grid point lines to use as center, 0,5,10,15, or 20 meters
        self.ct.setCenterLine(self.centerlinegp)
        #this method creates the centerline points through linear interpolation
        self.ct.createCenterline()

        # get coordinate for all different type of points
        self.gpCord = self.ct.getgp()
        self.stCord = self.ct.getst()
        self.tspCord = self.ct.gettsp()
        self.tpCord = self.ct.gettp()
        self.otherCord = self.ct.getother()

        # do all the transformations possible before application launches
        # transform each point separately
        self.s1,self.n1,self.l1 = self.ct.transform(self.gpCord[:,0].astype(np.float),self.gpCord[:,1].astype(np.float))
        self.s2,self.n2,self.l2 = self.ct.transform(self.stCord[:,0].astype(np.float),self.stCord[:,1].astype(np.float))
        self.s3,self.n3,self.l3 = self.ct.transform(self.tpCord[:,0].astype(np.float),self.tpCord[:,1].astype(np.float))
        self.s4,self.n4,self.l4 = self.ct.transform(self.tspCord[:,0].astype(np.float),self.tspCord[:,1].astype(np.float))
        self.s5,self.n5,self.l5 = self.ct.transform(self.otherCord[:,0].astype(np.float),self.otherCord[:,1].astype(np.float))
        # all data points
        self.s,self.n,self.l = self.ct.transform(self.xdata,self.ydata)

        self.initialize()
        self.frame.Show(True) #show the GUI on the screen

        # main loop for wx GUI
        self.app.MainLoop()

    def initialize(self):

        # create canvas and toolbar objects
        self.fig1 = plt.figure()
        plt.grid(True)
        self.fig2 = plt.figure()
        plt.grid(True)
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

        hbox.Add(self.toolbar1,border=5)

        self.btn1 = wx.Button(self.panel, label='launch 3D', size=(80,30))
        hbox.Add(self.btn1, flag=wx.CENTER, border=5)

        self.checkbox1 = wx.CheckBox(self.panel, label='Display Colors')
        hbox.Add(self.checkbox1, flag=wx.CENTER, border=5)

        hbox.Add((300,-1))
        hbox.Add(self.toolbar2, border=5)

        self.btn2 = wx.Button(self.panel, label='launch 3D', size=(80,30))
        hbox.Add(self.btn2, flag=wx.CENTER, border=5)

        self.checkbox2 = wx.CheckBox(self.panel, label='Display Colors')
        hbox.Add(self.checkbox2, flag=wx.CENTER, border=5)

        hbox1.Add(self.canvas1, border=5)
        hbox1.Add((20,-1))
        hbox1.Add(self.canvas2, border=5)

        vbox.Add(hbox,border=5)
        vbox.Add(hbox1,border=5)

        # event handlers
        self.btn1.Bind(wx.EVT_BUTTON, self.plot3dpre)
        self.btn2.Bind(wx.EVT_BUTTON, self.plot3dpost)
        self.checkbox1.Bind(wx.EVT_CHECKBOX, self.preChangeColor)
        self.checkbox2.Bind(wx.EVT_CHECKBOX, self.postChangeColor)


        self.panel.SetSizer(vbox,wx.CENTER)

    def preChangeColor(self,event):

        if self.checkbox1.IsChecked(): #color on
            self.ax1.cla()
            self.preColor = True
            self.prePlot(self.preColor,self.prePoints)
            self.ax1.legend()
            self.canvas1.draw()
        else: # color off
            self.ax1.cla()
            self.preColor = False
            self.prePlot(self.preColor,self.prePoints)
            self.ax1.legend_ = None # turn off legend
            self.canvas1.draw()

    def postChangeColor(self,event):

        if self.checkbox2.IsChecked(): #color on
            self.ax2.cla()
            self.postColor = True
            self.postPlot(self.postColor,self.postPoints)
            self.ax2.legend()
            self.canvas2.draw()
        else: # color off
            self.ax2.cla()
            self.postColor = False
            self.postPlot(self.postColor,self.postPoints)
            self.ax2.legend_ = None # turn off legend
            self.canvas2.draw()

    def prePlot(self,color,points):
        self.ax1.grid(True)
        if all(points): #if plotting all points

            if color: #if color
                self.ax1.scatter(self.gpCord[:,0].astype(np.float),self.gpCord[:,1].astype(np.float),color='k',label='grid points')
                self.ax1.scatter(self.stCord[:,0].astype(np.float),self.stCord[:,1].astype(np.float),color='b', label='soil tubes')
                self.ax1.scatter(self.tpCord[:,0].astype(np.float),self.tpCord[:,1].astype(np.float),color='r', label='temperature probes')
                self.ax1.scatter(self.tspCord[:,0].astype(np.float),self.tspCord[:,1].astype(np.float),color='g', label='total station points')
                self.ax1.scatter(self.otherCord[:,0].astype(np.float),self.otherCord[:,1].astype(np.float),color='c', label='other')

            else: # no color
                self.ax1.scatter(self.xdata,self.ydata)


    def postPlot(self,color,points):
        self.ax2.grid(True)
        if all(points): #if plotting all points

            if color: #if color
                self.ax2.scatter(self.s1*self.l1,self.n1,color='k',label='grid points')
                self.ax2.scatter(self.s2*self.l2,self.n2,color='b', label='soil tubes')
                self.ax2.scatter(self.s3*self.l3,self.n3,color='r', label='temperature probes')
                self.ax2.scatter(self.s4*self.l4,self.n4,color='g', label='total station points')
                self.ax2.scatter(self.s5*self.l5,self.n5,color='c', label='other')
            else: # no color
                self.ax2.scatter(self.s*self.l,self.n)

    # show the current plots in 3d in a separate window
    def plot3dpre(self,event):
        fig = plt.figure()
        ax = fig.add_subplot((111), projection='3d')
        ax.scatter(self.xdata,self.ydata,self.zdata)
        plt.show()

    def plot3dpost(self,event):
        fig = plt.figure()
        ax = fig.add_subplot((111), projection='3d')
        ax.scatter(self.s*self.l,self.n,self.zdata)
        plt.show()