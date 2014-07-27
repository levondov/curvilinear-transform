import wx
import sys
import threading
import time
import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as Navbar
from matplotlib.backends.backend_wx import NavigationToolbar2Wx as Navbar2
from matplotlib.backend_bases import NavigationToolbar2 as Nb
from coord_transform import coord_transform

class ctGUI:

    ct = 0
    size = (1300,550)
    frame = 0
    panel = 0

    #matplotlib stuff
    fig1 = 0
    fig2 = 0
    ax1 = 0
    ax2 = 0
    canvas1 = 0
    canvas2 = 0
    toolbar1 = 0
    toolbar2 = 0

    # data
    xdata = 0
    ydata = 0
    zdata = 0
    locations = ''

    # coord transform class variables
    centerlinegp = 10 #default line used as centerline is grid point 10m locations

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
        self.toolbar1 = Navbar2(self.canvas1) #use matplotlib backend to grab toolbar
        self.toolbar2 = Navbar2(self.canvas2)

        cid = self.fig1.canvas.mpl_connect('button_press_event', self.on_press)

        # default settings (upon launch) is just transform all data points
        self.ct.setCenterLine(self.centerlinegp) #decide which grid point lines to use as center, 0,5,10,15, or 20 meters
        self.ct.createCenterline() #this method creates the centerline points through linear interpolation
        s,n,l = self.ct.transform(self.xdata,self.ydata)

        # default plotting when launching

        self.ax1.scatter(self.xdata,self.ydata)
        self.ax2.scatter(s*l,n)

        # create sizers
        vbox = wx.BoxSizer(wx.VERTICAL)
        hbox = wx.BoxSizer(wx.HORIZONTAL)
        hbox1 = wx.BoxSizer(wx.HORIZONTAL)

        hbox.Add(self.toolbar1,border=5)
        hbox.Add((370,-1))
        hbox.Add(self.toolbar2, border=5)
        hbox1.Add(self.canvas1, border=5)
        hbox1.Add((20,-1))
        hbox1.Add(self.canvas2, border=5)

        vbox.Add(hbox,border=5)
        vbox.Add(hbox1,border=5)

        self.panel.SetSizer(vbox,wx.CENTER)

    def on_press(self,event):
        print('geee')