
import numpy as np  # for array operations
import numpy.matlib as npm
import copy




class coord_transform:
    # coordinate points
    xdata = 0
    ydata = 0
    zdata = 0
    # coord types/locations
    locdata = ''

    # lists containing x,y,z, and loc(names) of each point
    # make first row a empty row with zeros that will be removed later
    stCord = np.zeros((1, 4))  # soil tubes
    tpCord = np.zeros((1, 4))  # temperature probes
    gpCord = np.zeros((1, 4))  # grid points
    tspCord = np.zeros((1, 4))  # total station points
    otherCord = np.zeros((1, 4))  # all other points, C, S, P, etc...

    # grid points
    gp00 = np.zeros((1, 4))
    gp05 = np.zeros((1, 4))
    gp10 = np.zeros((1, 4))
    gp15 = np.zeros((1, 4))
    gp20 = np.zeros((1, 4))

    #centerline coordinates, see line 135 and below for more info
    Cup = np.zeros((1, 2))
    Cdown = np.zeros((1, 2))
    centerlineX = np.array([])
    centerlineY = np.array([])

    def __init__(self, x, y, z, loc):
        # define variables
        self.xdata = x
        self.ydata = y
        self.zdata = z
        self.locdata = loc
        np.set_printoptions(threshold='nan')
        self.organize()

    # creates nicely formatted lists from initial data
    # after the lists for each object will look like -> [x,y,z,name]
    def organize(self):

        # organize all the data points into their corresponding lists
        for i in range(0, len(self.xdata)):
            if self.locdata[i][0] == 'H':
                self.stCord = np.append(self.stCord, [[self.xdata[i], self.ydata[i], self.zdata[i], self.locdata[i]]],
                                        axis=0)
            elif self.locdata[i][0] == 'T':
                self.tpCord = np.append(self.tpCord, [[self.xdata[i], self.ydata[i], self.zdata[i], self.locdata[i]]],
                                        axis=0)
            elif self.locdata[i][0] == 'G':
                self.gpCord = np.append(self.gpCord, [[self.xdata[i], self.ydata[i], self.zdata[i], self.locdata[i]]],
                                        axis=0)
            elif self.locdata[i][0] == 'O' or self.locdata[i][0] == 'N' or self.locdata[i][0] == 'o':
                self.tspCord = np.append(self.tspCord, [[self.xdata[i], self.ydata[i], self.zdata[i], self.locdata[i]]],
                                         axis=0)
            else:
                self.otherCord = np.append(self.otherCord,
                                           [[self.xdata[i], self.ydata[i], self.zdata[i], self.locdata[i]]], axis=0)
        # remove first row of empty zeros
        self.stCord = np.delete(self.stCord, 0, 0)
        self.tpCord = np.delete(self.tpCord, 0, 0)
        self.gpCord = np.delete(self.gpCord, 0, 0)
        self.tspCord = np.delete(self.tspCord, 0, 0)
        self.otherCord = np.delete(self.otherCord, 0, 0)

        # put the grid point list in order by increasing grid size, i.e. G000.00 --> G640.20
        order = np.argsort(self.gpCord[:, 3])
        temp = copy.deepcopy(self.gpCord) # needed otherwise temp&self.gpCord would point to the same list
        for i in range(0,len(self.gpCord)):
            self.gpCord[i][0] = temp[order[i]][0]
            self.gpCord[i][1] = temp[order[i]][1]
            self.gpCord[i][2] = temp[order[i]][2]
            self.gpCord[i][3] = temp[order[i]][3]

        # organize grid points into lists by location, 00, 05, 10, 15, 20
        for i in range(0, len(self.gpCord)):
            if '.00' in self.gpCord[i][3]:
                self.gp00 = np.append(self.gp00, [[self.gpCord[i][0], self.gpCord[i][1], self.gpCord[i][2], self.gpCord[i][3]]], axis=0)
            elif '.05' in self.gpCord[i][3]:
                self.gp05 = np.append(self.gp05, [[self.gpCord[i][0], self.gpCord[i][1], self.gpCord[i][2], self.gpCord[i][3]]], axis=0)
            elif '.10' in self.gpCord[i][3]:
                self.gp10 = np.append(self.gp10, [[self.gpCord[i][0], self.gpCord[i][1], self.gpCord[i][2], self.gpCord[i][3]]], axis=0)
            elif '.15' in self.gpCord[i][3]:
                self.gp15 = np.append(self.gp15, [[self.gpCord[i][0], self.gpCord[i][1], self.gpCord[i][2], self.gpCord[i][3]]], axis=0)
            elif '.20' in self.gpCord[i][3]:
                self.gp20 = np.append(self.gp20, [[self.gpCord[i][0], self.gpCord[i][1], self.gpCord[i][2], self.gpCord[i][3]]], axis=0)

        # remove first row of zeros
        self.gp00 = np.delete(self.gp00, 0, 0)
        self.gp05 = np.delete(self.gp05, 0, 0)
        self.gp10 = np.delete(self.gp10, 0, 0)
        self.gp15 = np.delete(self.gp15, 0, 0)
        self.gp20 = np.delete(self.gp20, 0, 0)

    # set the grid points to be used as the centerline for the interpolation, default value is 10
    def setCenterLine(self,gpNum=10):

        # xCL -> x centerline point, yCL -> y centerline point
        if gpNum == 0:
            xCL = self.gp00[:,0].astype(np.float)
            yCL = self.gp00[:,1].astype(np.float)
        elif gpNum == 5:
            xCL = self.gp05[:,0].astype(np.float)
            yCL = self.gp05[:,1].astype(np.float)
        elif gpNum == 10:
            xCL = self.gp10[:,0].astype(np.float)
            yCL = self.gp10[:,1].astype(np.float)
        elif gpNum == 15:
            xCL = self.gp15[:,0].astype(np.float)
            yCL = self.gp15[:,1].astype(np.float)
        elif gpNum == 20:
            xCL = self.gp20[:,0].astype(np.float)
            yCL = self.gp20[:,1].astype(np.float)
        else:
            print('Error! Something went wrong, make sure you input an integer value of 0,5,10,15, or 20')

        # remove duplicates
        xCL = self.unique(xCL)
        yCL = self.unique(yCL)

        # organize from increasing to decreasing order,
        order = np.argsort(xCL)
        tempx = copy.deepcopy(xCL)
        tempy = copy.deepcopy(yCL)
        for i in range(0,len(xCL)):
            xCL[i] = tempx[order[i]]
            yCL[i] = tempy[order[i]]

        # Separate transect into two groups, points below/above -200Y,
        # include middle point in both to connect lines.
        for i in range(0,len(xCL)):
            if yCL[i] >= -200:
                self.Cup = np.append(self.Cup, [[xCL[i],yCL[i]]], axis=0)
            if yCL[i] < -195:
                self.Cdown = np.append(self.Cdown, [[xCL[i],yCL[i]]], axis=0)

        # remove first row of zeros
        self.Cup = np.delete(self.Cup,0,0)
        self.Cdown = np.delete(self.Cdown,0,0)

        # reverse order of Cup for points go from North -> South
        temp = copy.deepcopy(self.Cup)
        for i in range(0,len(self.Cup)):
            self.Cup[i][0] = temp[len(self.Cup)-1-i][0]
            self.Cup[i][1] = temp[len(self.Cup)-1-i][1]

    # use linear (affine) interpolation to create the centerline used for transformation
    def createCenterline(self):

        # define lists to be used for linear interpolation
        coeffup = np.zeros((len(self.Cup)-1,2))
        coeffdown = np.zeros((len(self.Cdown)-1,2))

        # solve least squares equation between every 2 points, [[1,xk],[1,xk+1]] = [yk,yk+1]
        for k in range(0,len(self.Cup)-1):
            A = np.array([[self.Cup[k][0],1],[self.Cup[k+1][0],1]])
            b = np.array([self.Cup[k][1],self.Cup[k+1][1]])
            x = np.linalg.solve(A, b)
            coeffup[k][0] = x[0]
            coeffup[k][1] = x[1]

        for k in range(0,len(self.Cdown)-1):
            A = np.array([[self.Cdown[k][0],1],[self.Cdown[k+1][0],1]])
            b = np.array([self.Cdown[k][1],self.Cdown[k+1][1]])
            x = np.linalg.solve(A, b)
            coeffdown[k][0] = x[0]
            coeffdown[k][1] = x[1]

        # create centerline coordinates. 10 evenly spaces points between each xk,xk+1
        for k in range(0,len(self.Cup)-1):
            x = np.linspace(self.Cup[k][0],self.Cup[k+1][0],10)
            y = np.polyval([coeffup[k][0],coeffup[k][1]],x)
            self.centerlineX = np.append(self.centerlineX,x)
            self.centerlineY = np.append(self.centerlineY,y)

        for k in range(0,len(self.Cdown)-1):
            x = np.linspace(self.Cdown[k][0],self.Cdown[k+1][0],10)
            y = np.polyval([coeffdown[k][0],coeffdown[k][1]],x)
            self.centerlineX = np.append(self.centerlineX,x)
            self.centerlineY = np.append(self.centerlineY,y)

    def getCenterline(self):
        return self.centerlineX,self.centerlineY

    def transform(self,x,y):
        l,junk = self.arclength()

        # output P is the closest point on the centerline to each point.
        # output n is the actual distance to point from centerline
        # output s is the fractional arc lengths(ds) for the centerline
        P,n,s = self.distance2curve(x,y)

        # these arrays contain the closest point on the centerline to each
        # corresponding point being mapped
        PX = P[0,:]
        PY = P[1,:]

        # the fractional arclength is defined as ds, and ds = sqrt(dx/dt^2 + dy/dt^2)
        # given fractional arc lengths, and centerline, Ct gives you the parametric derivatives dx/dt and dy/dt
        # Again, using linear method
        Ct = self.interparc(s)

        # Calculate normal vectors on centerline curve using parametric equations
        Cn = np.array([Ct[1,:], -1*Ct[0,:]])

        # Calculate angles between points and normal points
        angleMP = np.arctan2(PY-y, PX-x)

        # Calculate angles of normal vectors
        angleTn = np.arctan2(Cn[1,:], Cn[0,:])

        Nsign = np.abs(angleTn-angleMP) < (np.pi/2)
        for i in range(0,len(n)):
            if ~Nsign[i]:
                n[i] = n[i]*-1

        return s,n,l

    # compute chordal linear arclengths
    def arclength(self):
        # compute the sqrt( (xk+1-xk)^2 + (yk+1-yk) )
        # NOTE: there is probably a simpler way to do this

        # take the difference between xk and xk+1, yk and yk+1
        diffx = np.diff(self.centerlineX)
        diffy = np.diff(self.centerlineY)

        # square each of these differences
        squaredx = np.square(diffx)
        squaredy = np.square(diffy)

        # sum the two differences
        sum = np.zeros(len(squaredx))
        for i in range(0,len(squaredx)):
            sum[i] = squaredx[i] + squaredy[i]

        seglen = np.sqrt(sum)
        arclen = np.sum(seglen)

        return arclen, seglen

    # compute the distances from a two dimensional point to a curve (our centerline)
    def distance2curve(self,x,y):

        # calculate the scaled chordal linear arclength from 0 -> 1
        l,seglen = self.arclength()
        t = np.array([0])
        t = np.append(t,np.cumsum(seglen)/l)
        breaks = copy.deepcopy(t)

        # We need to build some parametric splines.
        # compute the splines, storing the polynomials in one 3-d array
        # the breaks for the splines will be t0, unless spline got fancy
        ppsegsX = np.array([[0],[0]])
        ppsegsY = np.array([[0],[0]])
        ppsegsX = ppsegsX.astype(np.float)
        ppsegsY = ppsegsY.astype(np.float)

        dt = np.diff(t)
        for i in range(0,len(self.centerlineX)-1):
            # when ever we try to divide by 0, set that value in array = 0 (was having problems with NaN)
            if dt[i] != 0:
                ppsegsX = np.insert(ppsegsX,i,[(self.centerlineX[i+1] - self.centerlineX[i])/dt[i],self.centerlineX[i]],axis=1)
                ppsegsY = np.insert(ppsegsY,i,[(self.centerlineY[i+1] - self.centerlineY[i])/dt[i],self.centerlineY[i]],axis=1)
            else:
                ppsegsX = np.insert(ppsegsX,i,[0,self.centerlineX[i]],axis=1)
                ppsegsY = np.insert(ppsegsY,i,[0,self.centerlineY[i]],axis=1)
        temp,length = np.shape(ppsegsX)
        ppsegsX = np.delete(ppsegsX,length-1,1)
        ppsegsY = np.delete(ppsegsY,length-1,1)

        # for each point in xdata,ydata, find the closest point to those
        # in centerlinex, centerliney.
        rowindex,columnindex,values = self.ipdm(x,y)

        # initialize the return variables, using the closest point
        # found in the centerlineX,Y.
        xy = np.array([[0],[0]])
        xy = xy.astype(np.float)

        for i in range(0,len(columnindex)):
            xy = np.insert(xy,i,[self.centerlineX[columnindex[i]],self.centerlineY[columnindex[i]]],axis=1)
        temp,length = np.shape(xy)
        xy = np.delete(xy,length-1,1)

        distance = values
        t = t[columnindex]
        #
        # Now do some checking
        #
        m = len(x)
        n = len(self.centerlineX)
        #segments, when there are many segments, but not many points to map.
        if n >= 5*m:
            # many segments, so loop over the points in x,y
            for i in range(0,m):
                xyi = np.array([x[i],y[i]])
                tnum = np.zeros([1,len(breaks) - 1])
                tden = copy.deepcopy(tnum)

                tden = tden + np.square(ppsegsX[0,:]) + np.square(ppsegsY[0,:])
                tnum = tnum + ppsegsX[0,:]*(xyi[0] - ppsegsX[1,:]) + ppsegsY[0,:]*(xyi[1] - ppsegsY[1,:])
                tmin = tnum/tden

                # toss out any element of tmin that is less than or equal to
                # zero, or or is greater than dt for that segment.
                tmin[tmin<=0] = np.nan
                tmin[tmin >= np.diff(breaks)] = np.nan

                # for any segments with a valid minimum distance inside the
                # segment itself, compute that distance.
                dmin = np.zeros([1,len(breaks)-1])
                dmin = dmin + np.square(ppsegsX[0,:]*tmin + ppsegsX[1,:] - xyi[0])
                dmin = dmin + np.square(ppsegsY[0,:]*tmin + ppsegsX[1,:] - xyi[0])
                dmin = np.sqrt(dmin)

                # get minimum distance between these segments
                mindist = np.nanmin(dmin)
                if ~np.isnan(mindist):
                    minindex = np.nanargmin(dmin)
                if ~np.isnan(mindist) and (distance[i] > mindist):
                    # there is a best segment, better than the
                    # closest point from centerlineX,Y.
                    distance[i] = mindist
                    t[i] = tmin[0][minindex] + breaks[minindex]

                    xy[0,i] = ppsegsX[0,minindex]*tmin[0][minindex] + ppsegsX[1,minindex]
                    xy[1,i] = ppsegsY[0,minindex]*tmin[0][minindex] + ppsegsY[1,minindex]
        else:
            for i in range(0,n-1):
                # the i'th segment of the curve
                t1 = breaks[i]
                t2 = breaks[i+1]

                # Compute the location (in t) of the minimal distance
                # point to x,y, for all points.

                tden = np.square(ppsegsX[0,i]) + np.square(ppsegsY[0,i])
                tnum = ppsegsX[0,i]*(x[:] - ppsegsX[1,i]) + ppsegsY[0,i]*(y[:] - ppsegsY[1,i])
                tmin = tnum/tden
                # We only care about those points for this segment where there
                # is a minimal distance to the segment that is internal to the segment.
                k = np.array([])
                for l in range(0,len(tmin)):
                    if tmin[l] > 0 and tmin[l] < (t2-t1):
                        k = np.insert(k,0,l)
                nk =len(k)
                if nk > 0:
                    # for any points with a valid minimum distance inside the
                    # segment itself, compute that distance.
                    k = k.astype(np.int)
                    #xymin = ppsegsX[0,i]*tmin[k] + ppsegsX[1,i]
                    dmin = np.square((ppsegsX[0,i]*tmin[k] + ppsegsX[1,i] - x[k])) + np.square((ppsegsY[0,i]*tmin[k] + ppsegsY[1,i] - y[k]))
                    dmin = np.sqrt(dmin)

                    L = dmin < distance[k]
                    # this segment has a closer point
                    # closest point from curvexy.
                    if L.any():
                        distance[k[L]] = dmin[L]
                        t[k[L]] = tmin[k[L]] + breaks[i]
                        xy[0,k[L]] = ppsegsX[0,i]*tmin[k[L]] + ppsegsX[1,i]
                        xy[1,k[L]] = ppsegsY[0,i]*tmin[k[L]] + ppsegsY[1,i]

        # for the linear case, t is identical to the fractional arc length
        # along the curve.
        t_a = t;
        return xy,distance,t_a

    # find nearest neighbor, used by the method 'distance2curve'
    def ipdm(self,x,y):
        # get the lengths of data and centerline points
        n1 = len(x)
        n2 = len(self.centerlineX)

        # sqrt((xk-xk-1)^2 + (yk-yk-1)^2), basically distance formula from coord points to centerline
        dx = np.square(np.tile(np.array([x]).T,(1,n2)) - np.tile(self.centerlineX,(n1,1)))
        dy = np.square(np.tile(np.array([y]).T,(1,n2)) - np.tile(self.centerlineY,(n1,1)))
        dist = np.sqrt(dx + dy)

        # get the minimum values, as well as their indices, row wise
        val = np.amin(dist,axis=1)
        j = np.argmin(dist,axis=1)

        # return row index, column index, and values
        rowindex = np.linspace(0,len(x)-1,len(x))
        return rowindex,j,val

    # Interpolates new points at any fractional point along
    # the curve defined by a list of points in 2
    # dimensions. The curve may be defined by any sequence
    # of non-replicated points.
    # uses linear chordal approximation
    # returns parametric derivatives (dx/dt, dy/dt, dz/dt, ...) as an array.
    def interparc(self,t):

        # number of points to be interpolated
        nt = len(t)
        n = len(self.centerlineX)

        # Compute the chordal (linear) arclength of each segment.
        chordlen = np.sqrt(np.square(np.diff(self.centerlineX)) + np.square(np.diff(self.centerlineY)))

        # Normalize the arclengths to a unit total
        chordlen = chordlen/np.sum(chordlen)

        cumarc = np.array([0])
        cumarc = np.append(cumarc,np.cumsum(chordlen))

        tbins = np.digitize(t,cumarc)

        # fix problems at the ends
        for i in range(0,len(tbins)):
            if tbins[i] <= 0 or t[i] <= 0:
                tbins[i] = 1;
            elif tbins[i] >= n or t[i] >= 1:
                tbins[i] = n-1

        # interpolate
        s = (t - cumarc[tbins-1])/chordlen[tbins-1]

        centerlineXY = np.array([self.centerlineX,self.centerlineY])
        pt = centerlineXY[:,tbins-1] + (centerlineXY[:,tbins] - centerlineXY[:,tbins-1])*npm.repmat(s,2,1)

        # finally compute partial derivatives
        dudt = (centerlineXY[:,tbins] - centerlineXY[:,tbins-1])/npm.repmat(chordlen[tbins-1],2,1)

        return dudt

    # used as a fast way of removing duplicates in a list
    def unique(self, seq):
        # order preserving
        noDupes = []
        [noDupes.append(i) for i in seq if not noDupes.count(i)]
        return noDupes

    def getst(self):
        return self.stCord

    def gettp(self):
        return self.tpCord

    def gettsp(self):
        return self.tspCord

    def getother(self):
        return self.otherCord

    def getgp(self):
        return self.gpCord


