import sys
import numpy as np
import pyearth_sub as pe
from matplotlib import pyplot as plt
from scipy.interpolate import griddata
import csv
import os

sys.path.insert(0, ".")
from config import CONFIG

DATA_FORMAT_VERSION = 1   # increase according to compatibility

STATUS = { 'CONFIG' : CONFIG }

if CONFIG['OUTPUT_FILE']:
    Output_File_Path = CONFIG['OUTDIR'] + "/" + CONFIG['MODELNAME'] + "/"

    metafile = open(Output_File_Path + "meta", "rb")
    firstrow = metafile.readline()
    dataformat = int(firstrow.strip())
    metafile.close()

    if dataformat != DATA_FORMAT_VERSION:
        raise Exception("Incompatible data format")

    csvmetafile = open(Output_File_Path + "meta", "rb")
    csvmetareader = csv.reader(csvmetafile, delimiter=",", quotechar='"', quoting=csv.QUOTE_NONNUMERIC)
    # columns for metafile are: tstep, time_ma
else:
    raise Exception("No output according to input file")

tsteps = np.array((), dtype=int)
times = np.array(())
nnodes = np.array((), dtype=int)

minx = np.Inf
maxx = -np.Inf

rown = 0
for row in csvmetareader:
    rown += 1
    if rown < 3:
        continue
    tsteps = np.append(tsteps, int(row[0]))
    times = np.append(times, float(row[1]))
    csvfile = open(Output_File_Path + "nodes." + str(int(row[0])), "rb")
    csvreader = csv.reader(csvfile, delimiter=",", quotechar='"', quoting=csv.QUOTE_NONNUMERIC)
    nrow = 0
    for row in csvreader:
        nrow += 1
        if nrow > 1:
            maxx = max(maxx, row[2])
            minx = min(minx, row[2])
    csvfile.close()
    nrow -= 1 #header
    nnodes = np.append(nnodes, nrow)

print str(len(tsteps)) + " time steps"
print tsteps
print minx, maxx

csvmetafile.close()

maxnodes = int(max(nnodes))

T = np.zeros((len(tsteps), maxnodes)) * np.NaN
x = np.zeros((len(tsteps), maxnodes)) * np.NaN
t = np.zeros((len(tsteps), maxnodes)) * np.NaN

# data format
# "ix","t","x","T","k","cp","rho","H"

iit = -1

for it in tsteps:
    iit += 1
    csvfile = open(Output_File_Path + "nodes." + str(it), "rb")
    csvreader = csv.reader(csvfile, delimiter=",", quotechar='"', quoting=csv.QUOTE_NONNUMERIC)
    irow = -2
    for row in csvreader:
        irow += 1
        if irow < 0:
            continue  # skip header
        T[iit,irow] = float(row[3])
        x[iit,irow] = float(row[2])
        t[iit,irow] = float(row[1])
    csvfile.close()

    #print "nnodes: " + str(nnodes[iit])
    #print "irow: " + str(irow)

newxs = np.linspace(minx, maxx, 150)
Tinterp = np.zeros((len(newxs), len(times)))

for iit in range(len(times)):
    newT = newxs * 0.0
    pe.interpolate(x[iit,0:nnodes[iit]], T[iit,0:nnodes[iit]], newxs, newT, extrapolation=0.0)
    Tinterp[:,iit] = newT[:]

print times.shape
print newxs.shape
print Tinterp.shape
#CS = plt.contour(times, -newxs, Tinterp, 35, linewidth=0.5, colors='k')
CS = plt.contourf(times, -newxs, Tinterp, 35, cmap=plt.cm.jet)
plt.colorbar()
plt.show()

raise Exception("STOP")

xi = np.linspace(min(times),max(times),100)
yi = np.linspace(minx,maxx,100)
plot_x = np.reshape(t, np.size(t))
plot_y = np.reshape(x, np.size(x))
plot_z = np.reshape(T, np.size(T))


## to remove nans:
#idx = np.bitwise_or(np.bitwise_or(np.isnan(plot_x), np.isnan(plot_y)), np.isnan(plot_z))
#idx = np.bitwise_not(idx)
#plot_x = plot_x[idx]
#plot_y = plot_y[idx]
#plot_z = plot_z[idx]

zi = griddata((plot_x, plot_y), plot_z, (xi[None,:], yi[:,None]), method='linear', fill_value=0.0)

CS = plt.contour(xi,yi,zi,35,linewidths=0.5,colors='k')
CS = plt.contourf(xi,yi,zi,35,cmap=plt.cm.jet)
plt.colorbar() # draw colorbar
#plt.scatter(plot_x,plot_y,marker='o',c='b',s=5)
#  23 plt.xlim(-2,2)
plt.ylim(maxx,minx)
#  25 plt.title('griddata test (%d points)' % npts)
plt.show()