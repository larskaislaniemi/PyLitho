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
   
import interpolate as ip

def readData(path=None, outfilename=None, maxmodels=None, everytstep=1, exclcondition="False"):
    DATA_FORMAT_VERSION = 2 

    if path is None:
        #path = "/data/home/lkaislan/software/pylitho/output/002-ini/"
        path = "/globalscratch/lkaislan/pylitho/output/002-ini/"

    nmodels = 0
    modeldirs = []              # STR x nmodel
    model_mpi_ranks = []        # INT x nmodel
    model_params = []           # p x DBL x nmodel (p small)
    model_data = []             # HUUUGE
    model_times = []            # DBL x nsteps x nmodel (nsteps ~10^3)

    if maxmodels is None:
        maxmodels = 2**31

    for d in os.walk(path):
        if nmodels >= maxmodels:
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
                    modeldirs.append(d[0])
                    model_mpi_ranks.append(mpi_rank)
                    model_params.append(thismodel_params)
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
            model_times.append(thismodel_times)
        else: 
            # no metafile inside dir d, not an output dir
            sys.stdout.write("no output")
            continue
        
        sys.stdout.write("\n")
        pbar = pb.ProgressBar(widgets=[pb.AnimatedMarker(), ' ', pb.Percentage(), ' ', pb.ETA()], maxval=len(thismodel_tsteps)).start()
        
        # params found, meta found => read data
        thismodel_data = {}
        
        # find datafield names
        rnodefile = open(d[0] + "/" + "nodes." + str(thismodel_tsteps[0]), "rb")
        csvnoder = csv.reader(rnodefile, delimiter=",", quotechar='"')
        idatafieldname = 0
        datafieldnums = []
        datafieldnames = []
        for noderow in csvnoder:
            for datafieldname in noderow:
                thismodel_data[datafieldname] = []
                datafieldnums.append(idatafieldname)
                datafieldnames.append(datafieldname)
                idatafieldname = idatafieldname + 1
            break
        ndatafields = idatafieldname
        rnodefile.close()
        
        # actually read the data
        # use pandas module for sake of convenience...
        for it in range(thismodel_nsteps):
            pbar.update(it)
            nodedata = pd.read_csv(d[0] + "/" + "nodes." + str(thismodel_tsteps[it]))
            for datafieldname in datafieldnames:
                thismodel_data[datafieldname].append([]) # for new timestep
                thismodel_data[datafieldname][-1] = nodedata[datafieldname].tolist()
                
        pbar.finish()

        model_data.append(thismodel_data)
        
    
    alldata = {
        "nmodels" : nmodels,
        "modeldirs" : modeldirs,
        "model_mpi_ranks" : model_mpi_ranks, 
        "model_params" : model_params,
        "model_data" : model_data,
        "model_times" : model_times
    }
    
    if outfilename is not None:
        outfile = open(outfilename, 'wb')
        pickle.dump(alldata, outfile)
        outfile.close()
        
    return alldata


def plot_allgeotherms(data, field="T", attime=0, reduceValue=None, reduceDepth=None):
    N = data["nmodels"]
    
    if reduceValue is None:
        reduceValue = np.zeros(N)
    else:
        reduceValue = np.array([data["model_params"][i][reduceValue] for i in range(N)])
    
    if reduceDepth is None:
        reduceDepth = np.zeros(N)
    else:
        reduceDepth = np.array([data["model_params"][i][reduceDepth] for i in range(N)])
        
    initimes = np.array([data["model_times"][i][0] for i in range(N)])
    lookattimes = initimes + attime
    lookatsteps = [ip.findNearest(data["model_times"][i], lookattimes[i], notFoundVal=np.nan) for i in range(N)]
    lookatsteps = [lookatsteps[i][lookatsteps[i][2]] for i in range(N)]
    print lookatsteps
    

    lines = [plt.plot(np.array(data["model_data"][i][field][lookatsteps[i]])+reduceValue[i], -np.array(data["model_data"][i]["x"][lookatsteps[i]])+reduceDepth[i], '-') for i in range(N)]
    ret = [lines[i][0].set_label(data["model_params"][i]) for i in range(N)]
    plt.title(field)
    plt.legend()
    
    plt.show()
