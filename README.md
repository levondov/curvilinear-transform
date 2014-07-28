curvilinear-transform
=====================

2d curvilinear coordinate transformation

This program does a 2d curvilinear coordinate transformation on a 
set of data given a 'centerline' to work with. The program is built
in a way to work with a specific set of data, but can easily be 
modified to work with different data. There is also a built in GUI
using wxpython that displays plots pre and post transformation. The
program was originally written in Matlab, but was rewritten in python. 
The Matlab code is available apon request.

The following python packages are needed to be able to run:
numpy,matplotlib,wxpython

To start the program, simply run/compile the file ctMain.py

v1.2
- compile time improved by organizing data manually and not with ct.transform()
- added support for displaying certain categories in the plot
- 3d plots will now display according to the same settings as the 2d plots (i.e. with colors and categories)


v1.1
- added color option for data points
- added option to show figures in 3d
