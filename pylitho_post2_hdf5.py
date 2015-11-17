import os
import numpy as np
import csv
import pandas as pd
import sys
try:
   import cPickle as pickle
except:
   import pickle
import matplotlib.pyplot as plt
import progressbar as pb
import h5py     # Using HDF5 we only need max_mem = size of one model
import scipy.interpolate
   
import interpolate as ip

def sameval(val1, val2, eps=1e-6):
    try:
        tmp = val1[0]
    except Exception as e:
        arr1 = False
    else:
        arr1 = True
    
    try:
        tmp = val2[0]
    except Exception as e:
        arr2 = False
    else:
        arr2 = True
        
    if arr1 and arr2:
        if len(val1) != len(val2):
            raise Exception("size of val1 != size of val2")
        return [abs((val2[i]-val1[i])/min(val1[i],val2[i])) < eps for i in range(len(val1))]
    elif arr1:
        return [abs((val2-val1[i])/min(val1[i],val2)) < eps for i in range(len(val1))]
    else:
        diff = abs((val2-val1)/min(val1,val2))
        if diff > eps:
            return False
        else:
            return True

def readData(path, outfilename, maxmodels=None, everytstep=1, exclcondition="False"):
    DATA_FORMAT_VERSION = 2 

    h5f = h5py.File(outfilename, 'w')
    
    nmodels = 0
    model_datafields = []       # ndatafields  ## NB! We assume all models have same datafields

    for d in os.walk(path):
        if maxmodels is not None and nmodels >= maxmodels:
            break
        
        sys.stdout.write("\n"+d[0]+": ")
        
        if "params" in d[2]:
            thismodel_params = []
            rparamfile = open(d[0] + "/" + "params")
            for row in rparamfile:
                thismodel_params.append(eval(row))
            rparamfile.close()
        else:
            # not an output dir
            sys.stdout.write("no output")
            continue
        
        if "meta" in d[2]:
            rmetafile = open(d[0] + "/" + "meta")
            csvmetar = csv.reader(rmetafile, delimiter=",", quotechar='"')
            irow = 0
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
                elif irow == 4:
                    # => we actually have data!
                    if eval(exclcondition):
                        irow = -1
                        sys.stdout.write("excluded")
                        break
                    nmodels = nmodels + 1
                    thismodel_times = []
                    thismodel_tsteps = []
                
                if irow > 3:
                    if (irow-4) % everytstep > 0:
                        continue

                    # actual metadata entries
                    thismodel_times.append(float(row[1]))
                    thismodel_tsteps.append(int(row[0]))
            rmetafile.close()
            if irow == -1:
                # this model was excluded by exclcondition
                continue
            thismodel_nsteps = len(thismodel_tsteps)
        else: 
            # no metafile inside dir d, not an output dir
            sys.stdout.write("no output")
            continue
        
        sys.stdout.write("\n")
        pbar = pb.ProgressBar(widgets=[pb.AnimatedMarker(), ' ', pb.Percentage(), ' ', pb.ETA()], maxval=len(thismodel_tsteps)).start()
        
        # params found, meta found => read data

        inidata = pd.read_csv(d[0] + "/" + "nodes_ini")
        ncol = inidata.shape[1]
        nrow = inidata.shape[0]
        thismodel_datafields = inidata.columns.tolist()
        if nmodels == 1:
            model_datafields = thismodel_datafields
        else:
            if not model_datafields == thismodel_datafields:
                raise Exception("Datafields not consistent")
        del inidata
        
        thismodel_data = np.zeros((thismodel_nsteps,nrow,ncol),dtype=np.float64)
        thismodel_data[:,:,:] = np.nan
        
        for it in range(thismodel_nsteps):
            pbar.update(it)
            nodedata = pd.read_csv(d[0] + "/" + "nodes." + str(thismodel_tsteps[it]))
            thismodel_thisstep_data = np.array(nodedata,dtype=np.float64)
            thisstep_nrow = thismodel_thisstep_data.shape[0]
            if thisstep_nrow > nrow:
                thismodel_data.resize((thismodel_nsteps, thisstep_nrow, ncol))
                thismodel_data[:,nrow:,:] = np.nan
                nrow = thisstep_nrow
            thismodel_data[it,:thisstep_nrow,:] = thismodel_thisstep_data[:,:]
        del nodedata
        del thismodel_thisstep_data
        
        pbar.finish()

        modelgrp = h5f.create_group("models/"+str(nmodels-1))
        modelgrp.attrs.create("mpi_rank", data=mpi_rank)
        modelgrp.attrs.create("params", data=thismodel_params,dtype=np.float64)
        modelgrp.attrs.create("times", data=thismodel_times,dtype=np.float64)
        modelgrp.create_dataset("modeldata", data=thismodel_data)
                              
    h5f.attrs.create("nmodels", data=nmodels)
    h5f.attrs.create("ndatafields", data=len(model_datafields))
    h5f.attrs.create("datafields", data=model_datafields)
    
    h5f.close()


def plot_allgeotherms(datafile, field="T", attime=0, reduceValue=None, reduceDepth=None):
    
    h5f = h5py.File(datafile, "r")
    
    N = h5f.attrs['nmodels']
    
    if reduceValue is None:
        reduceValue = np.zeros(N)
    else:
        reduceValue = [h5f['models'][str(i)].attrs['params'][reduceValue] for i in range(N)]
    
    if reduceDepth is None:
        reduceDepth = np.zeros(N)
    else:
        reduceDepth = [h5f['models'][str(i)].attrs['params'][reduceDepth] for i in range(N)]

    initimes = np.array([h5f['models'][str(i)].attrs['times'][0] for i in range(N)])
    lookattimes = initimes + attime
    print lookattimes
    lookatsteps = [ip.findNearest(h5f['models'][str(i)].attrs['times'], lookattimes[i], notFoundVal=np.nan) for i in range(N)]
    lookatsteps = [lookatsteps[i][lookatsteps[i][2]] for i in range(N)]
    print lookatsteps
    
    if field not in h5f.attrs['datafields']:
        raise Exception("Invalid fieldname")
   
    ifield = np.where(h5f.attrs['datafields']==field)[0][0]
    idepthf = np.where(h5f.attrs['datafields']=='x')[0][0]

    
    lines = [plt.plot(h5f['models'][str(i)]['modeldata'][lookatsteps[i],:,ifield]+reduceValue[i], -h5f['models'][str(i)]['modeldata'][lookatsteps[i],:,idepthf]+reduceDepth[i], '-') for i in range(N)]
    ret = [lines[i][0].set_label(h5f['models'][str(i)].attrs['params']) for i in range(N)]
    plt.title(field)
    plt.legend()
    
    h5f.close()
    
    plt.show()
    

def plot_field_at_depth(datafile, depth, field="T", reduceValue=None, fixedDepth=False, maxn=1e31):
    
    h5f = h5py.File(datafile, "r")
    
    N = h5f.attrs['nmodels']
    N = min(N, maxn)
    
    if reduceValue is None:
        reduceValue = np.zeros(N)
    else:
        reduceValue = [h5f['models'][str(imodel)].attrs['params'][reduceValue] for imodel in range(N)]
    
    ifield = np.where(h5f.attrs['datafields']==field)[0][0]
    idepthf = np.where(h5f.attrs['datafields']=='x')[0][0]

    pbar = pb.ProgressBar(widgets=[pb.AnimatedMarker(), ' ', pb.Percentage(), ' ', pb.ETA()], maxval=N).start()
    plt.figure()
    for imodel in range(N):
        ts = h5f['models'][str(imodel)].attrs['times']
        nsteps = len(ts)
        if not fixedDepth:
            idepth = [ip.findNearest(h5f['models'][str(imodel)]['modeldata'][i,:,idepthf], depth, notFoundVal=np.nan) for i in range(nsteps)]
        else:
            depthCompensate = [h5f['models'][str(imodel)]['modeldata'][i,0,idepthf] for i in range(nsteps)]
            idepth = [ip.findNearest(h5f['models'][str(imodel)]['modeldata'][i,:,idepthf], depth+depthCompensate[i], notFoundVal=np.nan) for i in range(nsteps)]
        idepth = [idepth[i][idepth[i][2]] for i in range(nsteps)]
        Ts = [h5f['models'][str(imodel)]['modeldata'][i,idepth[i],ifield]+reduceValue[imodel] for i in range(nsteps)]
        line, = plt.plot(ts, Ts, '-')
        line.set_label(h5f['models'][str(imodel)].attrs['params'])
        pbar.update(imodel)
    plt.title(field)
    plt.legend()
    pbar.finish()

    h5f.close()
    
    plt.show()
    
    
def get_isovalue(datafile, isovalue, field="T", fixedDepth=False, maxn=1e31, interpolate=False, aligndepths=False, progbar=True):
    h5f = h5py.File(datafile, "r")
    
    N = h5f.attrs['nmodels']
    N = min(N, maxn)
    
    ifield = np.where(h5f.attrs['datafields']==field)[0][0]
    idepthf = np.where(h5f.attrs['datafields']=='x')[0][0]

    if progbar:
        pbar = pb.ProgressBar(widgets=[pb.AnimatedMarker(), ' ', pb.Percentage(), ' ', pb.ETA()], maxval=N).start()
    
    timeseries = {}
    
    for imodel in range(N):
        ts = h5f['models'][str(imodel)].attrs['times']
        nsteps = len(ts)
        locs = [np.where(h5f['models'][str(imodel)]['modeldata'][istep,:,ifield] > isovalue) for istep in range(nsteps)]
        locs = [((~np.isnan(h5f['models'][str(imodel)]['modeldata'][istep,:,ifield])).sum()-1) if len(locs[istep][0]) == 0 else locs[istep][0][0] for istep in range(nsteps)]

        if interpolate:
            thmd = h5f['models'][str(imodel)]['modeldata']
            diffT = [thmd[istep,locs[istep],ifield]-thmd[istep,locs[istep]-1,ifield] for istep in range(nsteps)]
            diffD = [thmd[istep,locs[istep],idepthf] - thmd[istep,locs[istep]-1,idepthf] for istep in range(nsteps)]
            depths = [thmd[istep,locs[istep]-1,idepthf] + (isovalue-thmd[istep,locs[istep]-1,ifield])*diffD[istep]/diffT[istep] for istep in range(nsteps)]
        else:
            depths = [h5f['models'][str(imodel)]['modeldata'][istep,locs[istep],idepthf] for istep in range(nsteps)]
            
        if not fixedDepth:
            pass
        else:
            depths = [depths[istep] - h5f['models'][str(imodel)]['modeldata'][istep,0,idepthf] for istep in range(nsteps)]
        
        if aligndepths:
            depths = [depths[i] - depths[0] for i in range(len(depths))]

        timeseries[h5f['models'].keys()[imodel]] = {}
        timeseries[h5f['models'].keys()[imodel]]['t'] = np.array(ts)
        timeseries[h5f['models'].keys()[imodel]]['z'] = -np.array(depths)
        if progbar:
            pbar.update(imodel)
            
    h5f.close()
    if progbar:
        pbar.finish()
    
    return timeseries


def plot_isovalue(datafile, isovalue=None, field="T", fixedDepth=False, maxn=1e31, interpolate=False, aligndepths=False, progbar=True, usedata=None, filtercol=None, filterval=None, filtertol=1e-6, aligntimes=False, ):
    if usedata is None:
        if isovalue is None:
            raise Exception ("isovalue is required")
        data = get_isovalue(datafile, isovalue, field=field, fixedDepth=fixedDepth, maxn=maxn, interpolate=interpolate, aligndepths=aligndepths, progbar=progbar)
    else:
        data = usedata

    h5f = h5py.File(datafile, "r")
    
    if filtercol is not None and filterval is not None:
        try: tmp = filtercol[0]
        except Exception as e: 
            filtercol = [filtercol]
            filterval = [filterval]
            filtertol = [filtertol]
        filteridx = [True for imodel in data.keys()]
        
        for ifilter in range(len(filtercol)):
            thisfilter = [sameval(h5f['models'][imodel].attrs['params'][filtercol[ifilter]], filterval[ifilter], eps=filtertol[ifilter]) for imodel in data.keys()]
            filteridx = [filteridx[i] & thisfilter[i] for i in range(len(filteridx))]
            chosenvalues = np.array([h5f['models'][imodel].attrs['params'][filtercol[ifilter]] for imodel in data.keys()])[np.where(thisfilter)]
            print "Filter", ifilter, ", chosen values between", min(chosenvalues), " <> ", max(chosenvalues)
    else:
        filteridx = [True for imodel in data.keys()]
    filteridx = np.where(filteridx)[0]
    if len(filteridx) == 0:
        h5f.close()
        raise Exception("All models filtered out")
    print "Using ", len(filteridx), " models out of ", len(data.keys())


    plt.figure()
    for key in np.array(data.keys())[filteridx]:
        if not aligntimes:
            line, = plt.plot(data[key]['t'], data[key]['z'], '-')
        else:
            line, = plt.plot(data[key]['t']-data[key]['t'][0], data[key]['z'], '-')
        line.set_label(key)
    plt.title(field)
    plt.show(block=False)
    
    h5f.close()
    return plt

def plot_isovalue_dchange(datafile, isovalue=None, field="T", fixedDepth=False, maxn=1e31, interpolate=False, aligntimes=False, progbar=True, paramOnX=[0], filtercol=None, filterval=None, usedata=None, filtertol=1e-6):
    h5f = h5py.File(datafile, "r")
    
    # Assumes that provided 'data' has same number of models than the hdf5 file
    
    if filtercol in paramOnX:
        print "WARNING: filtercol used in paramOnX"

    if usedata is None:
        if isovalue is None:
            h5f.close()
            raise Exception ("isovalue is required")
        data = get_isovalue(datafile, isovalue, field=field, fixedDepth=fixedDepth, maxn=maxn, interpolate=interpolate, aligndepths=False, progbar=progbar)
    else:
        data = usedata
    
    keyidx = dict([(key, data.keys().index(key)) for key in data.keys()])
    
    maxvallocs = [np.where(data[key]['z']==np.max(data[key]['z']))[0] for key in data.keys()]
    maxvallocs = [np.nan if len(maxvallocs[i]) == 0 else maxvallocs[i][0] for i in range(len(maxvallocs))]
    maxvals = np.array([np.nan if np.isnan(maxvallocs[keyidx[key]]) else data[key]['z'][maxvallocs[keyidx[key]]] for key in data.keys()])
    startvals = np.array([data[key]['z'][0] for key in data.keys()])
    
    if filtercol is not None and filterval is not None:
        try: tmp = filtercol[0]
        except Exception as e: 
            filtercol = [filtercol]
            filterval = [filterval]
            filtertol = [filtertol]
        filteridx = [True for imodel in data.keys()]
        
        for ifilter in range(len(filtercol)):
            thisfilter = [sameval(h5f['models'][imodel].attrs['params'][filtercol[ifilter]], filterval[ifilter], eps=filtertol[ifilter]) for imodel in data.keys()]
            filteridx = [filteridx[i] & thisfilter[i] for i in range(len(filteridx))]
            chosenvalues = np.array([h5f['models'][imodel].attrs['params'][filtercol[ifilter]] for imodel in data.keys()])[np.where(thisfilter)]
            print "Filter", ifilter, ", chosen values between", min(chosenvalues), " <> ", max(chosenvalues)
    else:
        filteridx = [True for imodel in data.keys()]
    filteridx = np.where(filteridx)[0]
    if len(filteridx) == 0:
        h5f.close()
        raise Exception("All models filtered out")
    print "Using ", len(filteridx), " models out of ", len(data.keys())
        
    if len(paramOnX) == 1:
        plt.figure()
        param0s = np.array([h5f['models'][key].attrs['params'][paramOnX] for key in data.keys()])
        param0s = param0s[filteridx]
        plt.plot(param0s, (maxvals-startvals)[filteridx],'.')
        plt.show(block=False)
    elif len(paramOnX) == 2:
        plt.figure()
        param0s = np.array([h5f['models'][key].attrs['params'][paramOnX[0]] for key in data.keys()])
        param1s = np.array([h5f['models'][key].attrs['params'][paramOnX[1]] for key in data.keys()])
        param0s = param0s[filteridx]
        param1s = param1s[filteridx]
        
        xi = np.linspace(min(param0s), max(param0s), 100)
        yi = np.linspace(min(param1s), max(param1s), 100)
        
        zi = scipy.interpolate.griddata((param0s, param1s), (maxvals-startvals)[filteridx], (xi[None,:], yi[:,None]), method='linear', fill_value=np.nan)
        CS = plt.contourf(xi, yi, zi, 35, cmap=plt.cm.jet)
        cbar = plt.colorbar(CS)
        plt.show(block=False)
    elif len(paramOnX) == 3:
        return # this doesnt work yet
        param0s = np.array([h5f['models'][key].attrs['params'][paramOnX[0]] for key in data.keys()])
        param1s = np.array([h5f['models'][key].attrs['params'][paramOnX[1]] for key in data.keys()])
        param2s = np.array([h5f['models'][key].attrs['params'][paramOnX[2]] for key in data.keys()])
        param0s = param0s[filteridx]
        param1s = param1s[filteridx]
        param2s = param2s[filteridx]
        
        xi = np.linspace(min(param0s), max(param0s), 100)
        yi = np.linspace(min(param1s), max(param1s), 100)
        
        for p2 in param2s:
            plt.figure()
            idx = np.where(param2s == p2)[0]
            
            print idx, param0s, param1s, times
            print param0s[idx], param1s[idx], times[idx]
            zi = scipy.interpolate.griddata((param0s[idx], param1s[idx]), ((maxvals-startvals)[filteridx])[idx], (xi[None,:], yi[:,None]), method='linear', fill_value=np.nan)
            CS = plt.contourf(xi, yi, zi, 35, cmap=plt.cm.jet)
            plt.title(str(p2))
            cbar = plt.colorbar(CS)
            plt.show(block=False)
        
    else:
        h5f.close()
        raise Exception("invalid param value for paramOnX")
    h5f.close()
    return plt
    


def plot_isovalue_maxt(datafile, isovalue=None, field="T", fixedDepth=False, maxn=1e31, interpolate=False, aligntimes=False, progbar=True, paramOnX=[0], filtercol=None, filterval=None, usedata=None, filtertol=1e-6):
    h5f = h5py.File(datafile, "r")
    
    # Assumes that provided 'data' has same number of models than the hdf5 file
    
    if filtercol in paramOnX:
        print "WARNING: filtercol used in paramOnX"

    if usedata is None:
        if isovalue is None:
            h5f.close()
            raise Exception ("isovalue is required")
        data = get_isovalue(datafile, isovalue, field=field, fixedDepth=fixedDepth, maxn=maxn, interpolate=interpolate, aligndepths=False, progbar=progbar)
    else:
        data = usedata
    
    keyidx = dict([(key, data.keys().index(key)) for key in data.keys()])
    
    maxvallocs = [np.where(data[key]['z']==np.max(data[key]['z']))[0] for key in data.keys()]
    maxvallocs = [np.nan if len(maxvallocs[i]) == 0 else maxvallocs[i][0] for i in range(len(maxvallocs))]
    #maxvallocs = [maxvallocs[i] - maxvallocs[0] for i in range(len(maxvallocs))]
    if aligntimes:
        times = np.array([data[key]['t'][maxvallocs[keyidx[key]]]-data[key]['t'][0] for key in data.keys()])
    else:
        times = np.array([data[key]['t'][maxvallocs[keyidx[key]]] for key in data.keys()])
    
    if filtercol is not None and filterval is not None:
        try: tmp = filtercol[0]
        except Exception as e: 
            filtercol = [filtercol]
            filterval = [filterval]
            filtertol = [filtertol]
        filteridx = [True for imodel in data.keys()]
        
        for ifilter in range(len(filtercol)):
            thisfilter = [sameval(h5f['models'][imodel].attrs['params'][filtercol[ifilter]], filterval[ifilter], eps=filtertol[ifilter]) for imodel in data.keys()]
            filteridx = [filteridx[i] & thisfilter[i] for i in range(len(filteridx))]
            chosenvalues = np.array([h5f['models'][imodel].attrs['params'][filtercol[ifilter]] for imodel in data.keys()])[np.where(thisfilter)]
            print "Filter", ifilter, ", chosen values between", min(chosenvalues), " <> ", max(chosenvalues)
    else:
        filteridx = [True for imodel in data.keys()]
    filteridx = np.where(filteridx)[0]
    if len(filteridx) == 0:
        h5f.close()
        raise Exception("All models filtered out")
    print "Using ", len(filteridx), " models out of ", len(data.keys())
        
    times = times[filteridx]
    
    if len(paramOnX) == 1:
        plt.figure()
        param0s = np.array([h5f['models'][key].attrs['params'][paramOnX] for key in data.keys()])
        param0s = param0s[filteridx]
        plt.plot(param0s, times,'.')
        plt.show(block=False)
    elif len(paramOnX) == 2:
        plt.figure()
        param0s = np.array([h5f['models'][key].attrs['params'][paramOnX[0]] for key in data.keys()])
        param1s = np.array([h5f['models'][key].attrs['params'][paramOnX[1]] for key in data.keys()])
        param0s = param0s[filteridx]
        param1s = param1s[filteridx]
        
        xi = np.linspace(min(param0s), max(param0s), 100)
        yi = np.linspace(min(param1s), max(param1s), 100)
        
        zi = scipy.interpolate.griddata((param0s, param1s), times, (xi[None,:], yi[:,None]), method='linear', fill_value=np.nan)
        CS = plt.contourf(xi, yi, zi, 35, cmap=plt.cm.jet)
        cbar = plt.colorbar(CS)
        plt.show(block=False)
    elif len(paramOnX) == 3:
        return # this doesnt work yet
        param0s = np.array([h5f['models'][key].attrs['params'][paramOnX[0]] for key in data.keys()])
        param1s = np.array([h5f['models'][key].attrs['params'][paramOnX[1]] for key in data.keys()])
        param2s = np.array([h5f['models'][key].attrs['params'][paramOnX[2]] for key in data.keys()])
        param0s = param0s[filteridx]
        param1s = param1s[filteridx]
        param2s = param2s[filteridx]
        
        xi = np.linspace(min(param0s), max(param0s), 100)
        yi = np.linspace(min(param1s), max(param1s), 100)
        
        for p2 in param2s:
            plt.figure()
            idx = np.where(param2s == p2)[0]
            
            print idx, param0s, param1s, times
            print param0s[idx], param1s[idx], times[idx]
            zi = scipy.interpolate.griddata((param0s[idx], param1s[idx]), times[idx], (xi[None,:], yi[:,None]), method='linear', fill_value=np.nan)
            CS = plt.contourf(xi, yi, zi, 35, cmap=plt.cm.jet)
            plt.title(str(p2))
            cbar = plt.colorbar(CS)
            plt.show(block=False)
        
    else:
        h5f.close()
        raise Exception("invalid param value for paramOnX")
    h5f.close()
    return plt
    
