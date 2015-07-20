import numpy as np
import os
import csv

import matplotlib
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

import scipy.interpolate

DATA_FORMAT_VERSION = 2

var_params = []

#path = "output/restart_002"
# record_times = [0, 50, 100]  # in Myrs
# record_depths = [31e3, 95e3, 45e3]

PLOT = 1
PLOT_TIME = 4
PLOT_DEPTH_BOTTOM = 5
PLOT_TIME1 = 0
PLOT_TIME2 = 1
PLOT_DEPTH_BELOW_THRUST = 0
# 1 = non-dim near-bottom temperature, at time PLOT_TIME
# 2 = non-dim change in temperature at PLOT_DEPTH_BELOW_THRUST, from PLOT_TIME1 to PLOT_TIME2


#path = "output/steadystate_003"
#record_times = [1800, 1900]
#record_depths = [35e3, 75e3, 125e3]

path = "output/restart_003"
record_times = [0, 100, 200, 300, 400, 500]
record_depths = [31e3, 40e3, 60e3, 80e3, 100e3, 115e3]

rec_values = []
rec_ranks = []
rec_times = []
rec_depths = []
rec_params = []

for d in os.walk(path):
    meta_read = False
    params_read = False

    start_time = None
    irec = 0
    record_times_istep = [-1] * len(record_times)

    if ("params" in d[2]):
        params = []
        rparamfile = open(d[0] + "/" + "params")
	for row in rparamfile:
            params.append(float(row))
	rparamfile.close()
    else:
        # skipped
	continue

    if ("meta" in d[2]):
        rec_params.append([])
        for p in params:
            rec_params[-1].append(p)
        rmetafile = open(d[0] + "/" + "meta")
	csvmetar = csv.reader(rmetafile, delimiter=",", quotechar='"')
	rec_values.append([])
	rec_times.append([])
	rec_depths.append([])
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
		rec_ranks.append(mpi_rank)
		meta_read = True
		continue
	    elif irow == 3:
	        pass # just the header
            
	    if irow == 4:
	        start_time = float(row[1])

	    if irow > 3:
	        if (float(row[1])-start_time) >= record_times[irec]:
		    irec = irec + 1
		    record_times_istep[irec-1] = int(row[0])
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
		    for noderow in csvnoder:
		        if skip_header:
			    skip_header = False
			    continue
		        if inoderow == 0:
			    for i in range(len(record_depths)):
			        if float(noderow[DEPTHCOL]) > record_depths[i]:
				    # first depth recorder is larger than a requested recording depth
				    idepthrec = i+1
				    rec_values[-1][-1].append(np.nan)
				    rec_depths[-1][-1].append(np.nan)
		        inoderow = inoderow + 1
			if float(noderow[DEPTHCOL]) > record_depths[idepthrec]:
			    idepthrec = idepthrec + 1
			    rec_values[-1][-1].append(float(noderow[TEMPCOL]))
			    rec_depths[-1][-1].append(float(noderow[DEPTHCOL]))
			    if idepthrec >= len(record_depths):
			        alldepthsread = True
				break
		    rnodefile.close()

		    if irec >= len(record_times):
		        alltimesread = True
			break

	rmetafile.close()




######
## SECONDARY HANDLING
######

np_rec_values = np.array(rec_values)



#Tscale = np_rec_values[:,0,1]
#dT = (np_rec_values[:,1,1] - np_rec_values[:,0,1])/Tscale

if PLOT == 1:
    plottitle = "bottom temp change from time zero to t=" + str(record_times[PLOT_TIME])
    Tscale = np_rec_values[:,0,PLOT_DEPTH_BOTTOM]
    dT = np.array(rec_values)[:,PLOT_TIME,PLOT_DEPTH_BOTTOM]/Tscale
    p1 = np.array(rec_params)[:,0]
    p2 = np.array(rec_params)[:,1]
elif PLOT == 2:
    plottitle = "temp change from t=" + str(record_times[PLOT_TIME1]) + " to t=" + str(record_times[PLOT_TIME2]) + " at z = " + str(record_depths[PLOT_DEPTH_BELOW_THRUST])
    Tscale = np_rec_values[:,0,PLOT_DEPTH_BOTTOM]
    dT = np.array(rec_values)[:,PLOT_TIME2,PLOT_DEPTH_BELOW_THRUST] - np.array(rec_values)[:,PLOT_TIME1,PLOT_DEPTH_BELOW_THRUST]
    dT = dT / Tscale
    p1 = np.array(rec_params)[:,0]
    p2 = np.array(rec_params)[:,1]



xi = np.linspace(min(p1), max(p1))
yi = np.linspace(min(p2), max(p2))
zi = scipy.interpolate.griddata((p1,p2), dT, (xi[None,:],yi[:,None]), method='cubic')

plt.figure()
CS = plt.contourf(xi, yi, zi)
plt.colorbar()
plt.title(plottitle)
plt.scatter(p1,p2,marker='o',c='b',s=5)
plt.show()

#contour(xi, yi, zi)

#plt.figure()
#CS = plt.contour(X, Y, Z)
#plt.show()
        
