import numpy as np
import os
import csv

import pickle
import sys

import matplotlib
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from matplotlib import cm

from mpl_toolkits.mplot3d import Axes3D

import scipy.interpolate

from pyevtk.hl import pointsToVTK 

import pyvtk


DATA_FORMAT_VERSION = 2

var_params = []

#path = "output/restart_002"
# record_times = [0, 50, 100]  # in Myrs
# record_depths = [31e3, 95e3, 45e3]

PLOT = 21
PLOT_TIME = 0
PLOT_DEPTH_BOTTOM = 5
PLOT_DEPTH = 0
PLOT_TIME1 = 0
PLOT_TIME2 = 1
PLOT_DEPTH_BELOW_THRUST = 0
ADD_INITIAL_DEPTH = True
# 1 = non-dim near-bottom temperature, at time PLOT_TIME
# 2 = non-dim change in temperature at PLOT_DEPTH_BELOW_THRUST, from PLOT_TIME1 to PLOT_TIME2


#path = "output/steadystate_003"
#record_times = [1800, 1900]
#record_depths = [35e3, 75e3, 125e3]

path = "output/001-ini"
record_times = np.array([10,180,200])
record_depths = np.array([10e3,99e3,20e3])

depth_adjust = record_depths * 0.0

ndatapoints = 0

rec_values = []
rec_times = []
rec_depths = []
rec_params = []
rec_ranks = []

for d in os.walk(path):
    meta_read = False
    params_read = False

    start_time = None
    irec = 0

    if ("params" in d[2]):
        params = []
        rparamfile = open(d[0] + "/" + "params")
        for row in rparamfile:
            params.append(eval(row))
        rparamfile.close()
        if PLOT == 20 or PLOT == 21:
            if ADD_INITIAL_DEPTH:
                depth_adjust = record_depths * 0.0 
                depth_adjust = depth_adjust + params[0]
            else:
                depth_adjust = record_depths * 0.0
        else:
            depth_adjust = record_depths * 0.0
    else:
        # skipped
        continue

    if ("meta" in d[2]):
        rmetafile = open(d[0] + "/" + "meta")
        csvmetar = csv.reader(rmetafile, delimiter=",", quotechar='"')
        irow = 0
        alltimesread = False
        for row in csvmetar:
            irow = irow + 1
            if irow == 1:
                if int(row[0]) != DATA_FORMAT_VERSION:
                    raise Exception("Incompatible data format")
                else:
                    continue
            elif irow == 2:
                mpi_rank = int(row[0])
                meta_read = True
                continue
            elif irow == 3:
                # just the header
                pass
            
            if irow == 4:
                # we've got data!
                start_time = float(row[1])
                rec_values.append([])
                rec_times.append([])
                rec_depths.append([])
                rec_params.append([])

                for p in params:
                    rec_params[-1].append(p)
                rec_ranks.append(mpi_rank)
                ndatapoints = ndatapoints + 1

            if irow > 3:
                if (float(row[1])-start_time) >= record_times[irec]:
                    irec = irec + 1
                    rec_values[-1].append([])
                    rec_depths[-1].append([])
                    time = float(row[1]) - start_time
                    rec_times[-1].append(time)
                    rnodefile = open(d[0] + "/" + "nodes." + str(row[0]), "rb")
                    csvnoder = csv.reader(rnodefile, delimiter=",", quotechar='"')
                    inoderow = 0
                    idepthrec = 0
                    skip_header = True
                    DEPTHCOL = 2
                    TEMPCOL = 3
                    alldepthsread = False
                    lastdepth = -np.inf
                    for noderow in csvnoder:
                        if skip_header:
                            skip_header = False
                            continue
                        if inoderow == 0:
                            for i in range(len(record_depths)):
                                if float(noderow[DEPTHCOL]) > record_depths[i] + depth_adjust[i]:
                                    # first depth recorder is larger than a requested recording depth
                                    idepthrec = i
                                    rec_values[-1][-1].append(np.nan)
                                    rec_depths[-1][-1].append(np.nan)
                                    lastdepth = float(noderow[DEPTHCOL])
                                    inoderow = inoderow + 1
                                    continue
                        #if idepthrec >= len(record_depths):
                        #    # we have already recorded more than requested amount of depths
                        #    rec_depths[-1][-1].append(-999)
                        #    rec_values[-1][-1].append(-999)
                        #    idepthrec = idepthrec + 1
                        #    lastdepth = float(noderow[DEPTHCOL])
                        elif float(noderow[DEPTHCOL]) > record_depths[idepthrec] + depth_adjust[idepthrec]:
                            idepthrec = idepthrec + 1
                            rec_values[-1][-1].append(float(noderow[TEMPCOL]))
                            rec_depths[-1][-1].append(float(noderow[DEPTHCOL]))
                            lastdepth = float(noderow[DEPTHCOL])
                        if idepthrec >= len(record_depths):
                            alldepthsread = True
                            break
                        inoderow = inoderow + 1
                    rnodefile.close()
                    while len(rec_depths[-1][-1]) < len(record_depths):
                        # fill unfound depths with NaNs
                        rec_depths[-1][-1].append(np.nan)
                        rec_values[-1][-1].append(np.nan)

                    if irec >= len(record_times):
                        alltimesread = True
                        break
        # all times checked
        # fill unfound times with NaNs
        while len(rec_values[-1]) < len(record_times):
            rec_values[-1].append([])
            rec_depths[-1].append([])
            rec_times[-1].append(np.nan)
            while len(rec_depths[-1][-1]) < len(record_depths):
                rec_depths[-1][-1].append(np.nan)
                rec_values[-1][-1].append(np.nan)


        rmetafile.close()

######
## SECONDARY HANDLING
######

print ndatapoints

np_rec_values = np.array(rec_values)
print "---"
print np_rec_values.shape

np_rec_ranks = np.array(rec_ranks)
np_rec_times = np.array(rec_times)
np_rec_depths = np.array(rec_depths)

np.set_printoptions(threshold=np.nan)

#Tscale = np_rec_values[:,0,1]
#dT = (np_rec_values[:,1,1] - np_rec_values[:,0,1])/Tscale

if PLOT == 1:
    plottitle = "bottom temp change from time zero to t=" + str(record_times[PLOT_TIME])
    Tscale = np_rec_values[:,0,PLOT_DEPTH_BOTTOM]
    dT = np_rec_values[:,PLOT_TIME,PLOT_DEPTH_BOTTOM]/Tscale
    p1 = np.array(rec_params)[:,0]
    p2 = np.array(rec_params)[:,1]
elif PLOT == 2:
    plottitle = "temp change from t=" + str(record_times[PLOT_TIME1]) + " to t=" + str(record_times[PLOT_TIME2]) + " at z = " + str(record_depths[PLOT_DEPTH_BELOW_THRUST])
    Tscale = np_rec_values[:,0,PLOT_DEPTH_BOTTOM]
    dT = np_rec_values[:,PLOT_TIME2,PLOT_DEPTH_BELOW_THRUST] - np_rec_values[:,PLOT_TIME1,PLOT_DEPTH_BELOW_THRUST]
    dT = dT / Tscale
    p1 = np.array(rec_params)[:,0]
    p2 = np.array(rec_params)[:,1]
elif PLOT == 20:
    plottitle = "what are we plotting here? " + str(record_times[PLOT_TIME])
    Tscale = 1350.0
    print np_rec_values.shape
    print np_rec_times.shape
    pickedT = np_rec_values[:,PLOT_TIME,PLOT_DEPTH]
    ip0 = []
    ip1 = []
    ip2 = []
    iy0 = []
    for i in range(0,len(rec_params)):
        if rec_params[i][2] > -1:
            ip0.append(rec_params[i][0])   # thickness of the overthrust sheet
            ip1.append(rec_params[i][1][1]-rec_params[i][1][0]) # original extent of the lithosphere
            ip2.append(rec_params[i][2]) # erosion speed
            iy0.append(pickedT[i]) # value itself
    ip0 = np.array(ip0)
    ip1 = np.array(ip1)
    ip2 = np.array(ip2)
    iy0 = np.array(iy0)

    p0 = ip0
    p0lab = "overthrust sheet thickness"
    p1 = ip1
    p1lab = "orig. lithosphere thickness"
    p2 = ip2
    p2lab = "erosion speed"
    y0 = iy0
    y0lab = "temperature"

    print ip0
    print ip1
    print ip2
    print iy0

elif PLOT == 21:
    plottitle = "what are we plotting here? " + str(record_times[PLOT_TIME])
    Tscale = 1350.0
    print np_rec_values.shape
    print np_rec_times.shape
    pickedT = np_rec_values[:,PLOT_TIME,PLOT_DEPTH]
    ip0 = []
    ip1 = []
    ip2 = []
    iy0 = []
    for i in range(0,len(rec_params)):
        if rec_params[i][2] > -1: 
            ip0.append(rec_params[i][0])   # thickness of the overthrust sheet
            ip1.append(rec_params[i][1]) # bottom temp
            ip2.append(rec_params[i][2]) # erosion speed
            iy0.append(pickedT[i]) # value itself
    ip0 = np.array(ip0)
    ip1 = np.array(ip1)
    ip2 = np.array(ip2)
    iy0 = np.array(iy0)

    p0 = ip0 
    p0lab = "overthrust sheet thickness"
    p1 = ip1 
    p1lab = "lithosphere bottom temp"
    p2 = ip2 
    p2lab = "erosion speed"
    y0 = iy0 
    y0lab = "temperature"

    print ip0 
    print ip1 
    print ip2 
    print iy0 

x0i = np.linspace(min(p0), max(p0))
x1i = np.linspace(min(p1), max(p1))
x2i = np.linspace(min(p2), max(p2))
x0, x1, x2 = np.meshgrid(x0i, x1i, x2i)

p1 = np.ma.array(p1, mask=np.isnan(y0))
p2 = np.ma.array(p2, mask=np.isnan(y0))
p0 = np.ma.array(p0, mask=np.isnan(y0))
y0 = np.ma.array(y0, mask=np.isnan(y0))
print "----"
print len(p0)
print len(p1)
print len(p2)
print len(y0)
print np.ma.vstack((p0,p1,p2)).T.shape
orig_x = np.ma.vstack((p0,p1,p2)).T
print "---"
##print len(x0), len(x0[0]), len(x0[0][0])
queryx = np.vstack((x0.flatten(), x1.flatten(), x2.flatten())).T
##print queryx.shape
ipolate = scipy.interpolate.LinearNDInterpolator(points=orig_x, values=y0)
zi = ipolate(queryx)
##print len(queryx)
##print len(zi)
#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#ax.scatter(queryx[:,0], queryx[:,1], queryx[:,2])
#plt.show()
#zi = scipy.interpolate.griddata((p1,p2), y0, (xi[None,:],yi[:,None]), method='cubic')
NN = len(queryx[:,0])
picks = np.arange(0,NN,1)
pointsToVTK("./points", queryx[picks,0], queryx[picks,1], queryx[picks,2], data = {"temp" : zi[picks]})
print queryx.shape
print zi.shape
#sys.exit()

do3d = False

for x in np.unique(queryx[picks,2]):
    idx = np.where(queryx[picks,2] == x)
    fig = plt.figure()
    if do3d:
        ax = fig.add_subplot(111,projection='3d')
        ax.scatter((queryx[picks,0])[idx], (queryx[picks,1])[idx], (queryx[picks,2])[idx], c=(zi[picks])[idx], cmap=plt.hot(), s=5, edgecolors='none')
    else:
        ax = fig.add_subplot(111)
        sc = ax.scatter((queryx[picks,0])[idx], (queryx[picks,1])[idx], c=(zi[picks])[idx], cmap=plt.hot(), edgecolors='none')
        plt.colorbar(sc)
        #cs = ax.contour((queryx[picks,0])[idx], (queryx[picks,1])[idx], (zi[picks])[idx])
        #ax.clabel(cs)
        ax.set_xlabel(p0lab)
        ax.set_ylabel(p1lab)
    plt.show()

#contour(xi, yi, zi)

#plt.figure()
#CS = plt.contour(X, Y, Z)
#plt.show()
        
