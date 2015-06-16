from mpi4py import MPI
import sys
import numpy as np
import pyearth_sub as pe
#from matplotlib import pyplot as plt
import csv
import os
import math

sys.path.insert(0, ".")
from config import CONFIG

if CONFIG['MPI']:
    mpi_comm = MPI.COMM_WORLD
    mpi_size = mpi_comm.Get_size()
    mpi_rank = mpi_comm.Get_rank()
    mpi_file_postfix = "X" + ("{:0>3d}".format(mpi_rank))
    mpi_dir = "/x" + ("{:0>3d}".format(mpi_rank))
else:
    mpi_comm = None
    mpi_rank = 0
    mpi_size = 1
    mpi_file_postfix = ""
    mpi_dir = ""

kTfunc = pe.kT
cTfunc = pe.cT

DATA_FORMAT_VERSION = 2   # increase when output data format changes

STATUS = { 'CONFIG' : CONFIG }

if CONFIG['OUTPUT_FILE']:
    if CONFIG['MPI']:
        if mpi_rank == 0:
	    try:
	        os.mkdir(CONFIG['OUTDIR'] + "/" + CONFIG['MODELNAME'])
            except Exception as e:
	        if CONFIG['OUTPUT_OVERWRITE']:
		    print "Error in mkdir, directory exists? Don't care, OUTPUT_OVERWRITE == True"
		else:
		    raise e
	mpi_comm.Barrier()

        try:
	    os.mkdir(CONFIG['OUTDIR'] + "/" + CONFIG['MODELNAME'] + mpi_dir)
	except Exception as e:
	    if CONFIG['OUTPUT_OVERWRITE']:
	        print "Error in mkdir, directory exists? Don't care, OUTPUT_OVERWRITE == True"
	    else:
	        raise e
    else:
        try:
	    os.mkdir(CONFIG['OUTDIR'] + "/" + CONFIG['MODELNAME'])
        except Exception as e:
            if CONFIG['OUTPUT_OVERWRITE']:
                print "Error in mkdir, directory exists? Don't care, OUTPUT_OVERWRITE == True"
            else:
                raise e

    Output_File_Path = CONFIG['OUTDIR'] + "/" + CONFIG['MODELNAME'] + mpi_dir + "/"

    csvmetafile = open(Output_File_Path + "meta", "wb")
    csvmetafile.write(str(DATA_FORMAT_VERSION)+"\n")
    csvmetafile.write(str(mpi_rank)+"\n")
    csvmetafile.close()

    csvmetafile = open(Output_File_Path + "meta", "ab")
    csvmetawriter = csv.writer(csvmetafile, delimiter=",", quotechar='"', quoting=csv.QUOTE_NONNUMERIC)

    csvmetawriter.writerow(['tstep', 'time_ma'])

if CONFIG['RESTART']:
    Restart_In_Path = CONFIG['RESTART_INDIR'] + "/" + CONFIG['RESTART_MODELNAME'] + mpi_dir + "/"


STATUS['SECINYR'] = 60*60*24*365.25
SECINYR = STATUS['SECINYR']
STATUS['NX'] = CONFIG['NX']
STATUS['Erosion_Speed'] = CONFIG['EROSION_SPEED_M_MA'] / (1e6*SECINYR)
STATUS['Moho_Depth'] = CONFIG['MOHO_DEPTH_KM'] * 1e3
STATUS['L'] = np.array((CONFIG['L_KM'][0]*1e3, CONFIG['L_KM'][1]*1e3))
STATUS['MaxTime'] = CONFIG['MAXTIME_MA'] * SECINYR * 1e6
STATUS['MaxRunTime'] = CONFIG['MAXRUNTIME_MA'] * SECINYR * 1e6

if CONFIG['MPI'] and CONFIG['MPI_VARIATION_TYPE'] == 1:
    if CONFIG['KT_RELATION_TYPE'] != 1:
        raise Exception("Incompatible MPI_VARIATION_TYPE and KT_RELATION_TYPE")

    if mpi_size <= 1:
        raise Exception("MPI with one processor?")

    minval = CONFIG['MPI_VARIATION_PARAMS'][0]
    maxval = CONFIG['MPI_VARIATION_PARAMS'][1]
    CONFIG['KT_RELATION_PARAMS'][0] = minval + mpi_rank*(maxval-minval)/(mpi_size-1)
    print "MPI, #" + str(mpi_rank) + ": param is " + str(CONFIG['KT_RELATION_PARAMS'][0])

    paramfile = open(Output_File_Path + "params", "wb")
    for i in range(1):
        paramfile.write(str(CONFIG['KT_RELATION_PARAMS'][i]) + "\n")
    paramfile.close()

elif CONFIG['MPI'] and CONFIG['MPI_VARIATION_TYPE'] == 2:
    if CONFIG['KT_RELATION_TYPE'] != 2:
        raise Exception("Incompatible MPI_VARIATION_TYPE and KT_RELATION_TYPE")

    if mpi_size <= 1:
        raise Exception("MPI with one processor?")

    minval1 = CONFIG['MPI_VARIATION_PARAMS'][0]
    maxval1 = CONFIG['MPI_VARIATION_PARAMS'][1]
    minval2 = CONFIG['MPI_VARIATION_PARAMS'][2]
    maxval2 = CONFIG['MPI_VARIATION_PARAMS'][3]

    div1 = math.floor(math.sqrt(mpi_size))
    div2 = math.floor(mpi_size/div1)

    if mpi_rank >= div1*div2:
        print "MPI, #" + str(mpi_rank) + ": Skip\n"
	mpi_comm.Barrier()
        MPI.Finalize()
	sys.exit(0)

    
    CONFIG['KT_RELATION_PARAMS'][0] = minval1 + math.floor(mpi_rank/div2)*(maxval1-minval1)/(div1-1.0)
    CONFIG['KT_RELATION_PARAMS'][1] = minval2 + (mpi_rank - div2*math.floor(mpi_rank/div2))*(maxval2-minval2)/(div2-1.0)

    paramfile = open(Output_File_Path + "params", "wb")
    for i in range(2):
        paramfile.write(str(CONFIG['KT_RELATION_PARAMS'][i]) + "\n")
    paramfile.close()

    print "MPI, #" + str(mpi_rank) + ": params are " + str(CONFIG['KT_RELATION_PARAMS'][0]) + "," + str(CONFIG['KT_RELATION_PARAMS'][1]) + "\n"

# mesh
xs = np.linspace(STATUS['L'][0], STATUS['L'][1], num=STATUS['NX'])
# to demonstrate irreg grid:  + 4e3*(np.random.rand(NX)-0.5)


# initial fields

if CONFIG['RESTART']:
    STATUS['curTime'] = 0.0 #dummy

    Tin = np.array(())
    xsin = np.array(())

    rcsvmetafile = open(Restart_In_Path + "meta", "rb")
    csvmetareader = csv.reader(rcsvmetafile, delimiter=",", quotechar='"')

    firstRow = True
    rowfound = False
    mpirow = False
    for row in csvmetareader:
        if firstRow:
            firstRow = False
	    mpirow = True
            if int(row[0]) != DATA_FORMAT_VERSION:
                raise Exception("Incompatible data format version in restart")
            else:
                continue
        elif mpirow:
	    mpirow = False
	    if CONFIG['MPI']:
	        if int(row[0]) != mpi_rank:
		    raise Exception("Rank mismatch between running process and restart meta file")
		else:
		    continue

        if not row[0].isdigit():
            continue
        #print row[0], CONFIG['RESTART_TSTEP']
        if int(row[0]) == CONFIG['RESTART_TSTEP']:
            rowfound = True
            STATUS['curTime'] = float(row[1]) * SECINYR * 1e6

    if not rowfound:
        raise Exception("Error reading metadata file: tstep " + str(CONFIG['RESTART_TSTEP']) + " not found.")

    rcsvmetafile.close()

    rcsvfile = open(Restart_In_Path + "nodes." + str(CONFIG['RESTART_TSTEP']), "rb")
    csvreader = csv.reader(rcsvfile, delimiter=",", quotechar='"')

    skipFirst = True
    for row in csvreader:
        if skipFirst:
            skipFirst = False
        else:
            Tin = np.append(Tin, float(row[3]))
            xsin = np.append(xsin, float(row[2]))

    if len(xsin) == len(xs) and len(np.where(xsin == xs)) == STATUS['NX']:
        Tini = np.copy(Tin)
        # all nodes at same locations
    else:
        # nodes and/or extent in file read from restart file do not match those
        # configured in config file

        # changes in extent during restart not allowed
        if xsin[0] != xs[0] or xsin[-1] != xs[-1]:
            #print xsin
            #print xs
            raise Exception("Changes in extent during restart not allowed.")
        else:
            # match grid points, result must follow gridding configured in config file
            valueArrays = [Tin]
            pe.remesh(STATUS, xsin, xs, valueArrays)
            Tini = np.copy(Tin)

    rcsvfile.close()
    STATUS['ModelStartTime'] = STATUS['curTime']

else:
    STATUS['ModelStartTime'] = 0.0
    Tini = np.array([0.0] * STATUS['NX'])
    STATUS['curTime'] = 0.0

pe.initTemp(STATUS, Tini, xs)


if CONFIG['RHO0_TYPE'] == 0:
    # constant rho0
    rho = 0.0 * xs + 2800
elif CONFIG['RHO0_TYPE'] == 1:
    # different rho0 for crust/mantle
    rho = 0.0 * xs + 2800
    rho[xs > STATUS['Moho_Depth']] = 3300
else:
    raise Exception("Invalid RHO0_TYPE")

if CONFIG['K0_TYPE'] == 0:
    # constant k0
    K0val = 4.5
    k0 = 0.0 * xs + K0val
elif CONFIG['K0_TYPE'] == 1:
    K0val = 3.175
    k0 = 0.0 * xs + K0val
else:
    raise Exception("Invalid K0_TYPE")

if CONFIG['C0_TYPE'] == 0:
    #constant c0
    c0 = 0.0 * xs + 800
elif CONFIG['C0_TYPE'] == 1:
    #different c0 for crust/mantle
    c0 = 0.0 * xs + 800
    c0[xs >= STATUS['Moho_Depth']] = 700
elif CONFIG['C0_TYPE'] == 2:
    c0 = 0.0 * xs + 1100
else:
    raise Exception("Invalid C0_TYPE")

if CONFIG['H0_TYPE'] == 0:
    # no heat prod
    H = 0.0 * xs
elif CONFIG['H0_TYPE'] == 1:
    # constant heat prod in crust and mantle
    H = 0.0 * xs
    H[xs < STATUS['Moho_Depth']] = 2.5e-6 # W m^-3
    H[xs >= STATUS['Moho_Depth']] = 0.0 #1.6e-8
elif CONFIG['H0_TYPE'] == 2:
    # exponential decrease
    H_exp_param = 30e3
    H0val = 2.5e-6
    H = H0val * np.exp(-xs / H_exp_param)
else:
    raise Exception("Invalid H0_TYPE")


if CONFIG['RESTART_POST_MOD'] == 0:
    pass # nothing
elif CONFIG['RESTART_POST_MOD'] == 1:
    STATUS['L'][0] -= 30e3  # raise the surface by 30km => overthrust
    valueArrays = [Tini, H, c0, k0, rho]
    #print Tini
    pe.remesh(STATUS, xs, STATUS['L'], valueArrays, extrapolation=2)
    #print Tini
elif CONFIG['RESTART_POST_MOD'] == 2:
    H[xs < STATUS['Moho_Depth']] = H[xs < STATUS['Moho_Depth']] * 1.5
else:
    raise Exception("Invalid RESTART_POST_MOD")

# bnd cond
q0 = CONFIG['BND_BOT_HFLOW']
Tsurf = CONFIG['BND_TOP_TEMP']

T = np.copy(Tini)
new_T = np.copy(Tini)

it = 0

if not CONFIG['MPI']:
    plt.plot(Tini, -xs, 'ro-')
    plt.show()

csvfile = open(Output_File_Path + "nodes_ini", "wb")
csvfilewriter = csv.writer(csvfile, delimiter=",", quotechar='"', quoting=csv.QUOTE_NONNUMERIC)
csvfilewriter.writerow(["ix","t","x","T","k","cp","rho","H"])
for ix in range(len(xs)):
    csvfilewriter.writerow([ix, 1e-6*STATUS['curTime']/SECINYR, xs[ix], Tini[ix], k0[ix], c0[ix], rho[ix], H[ix]])
csvfile.close()

print "MPI, #" + str(mpi_rank) + ", start of time loop"
while (STATUS['curTime'] < STATUS['MaxTime']) and (STATUS['curTime']-STATUS['ModelStartTime'] < STATUS['MaxRunTime']):

    diffsChange = CONFIG['DIFFSCHANGE_ACCURACY'] + 1.0

    iterations = 0

    cp = cTfunc(CONFIG['CT_RELATION_TYPE'], c0, CONFIG['CT_RELATION_PARAMS'], T)
    k = kTfunc(CONFIG['KT_RELATION_TYPE'], k0, CONFIG['KT_RELATION_PARAMS'], T)
    idx = k <= 0
    if sum(idx) > 0:
        raise Exception("k == 0 at " + str(np.where(idx)))

    dtadjust = 1.0
    ndtadjust = 0

    while diffsChange > CONFIG['DIFFSCHANGE_ACCURACY']:
        iterations += 1
	if math.floor(math.log10(iterations)) > ndtadjust:
	    ndtadjust = ndtadjust + 1
	    dtadjust = dtadjust / 2.0
	    if dtadjust < 0.1:
	        print "NB! #" + str(mpi_rank) + ": dtadjust = " + str(dtadjust) + ", it = " + str(it) + "\n"

        dt = dtadjust * CONFIG['TSTEP_MULTI'] * pe.maxdt(k, rho, cp, xs)
        pe.diffstep(STATUS, T, new_T, xs, k, dt, Tsurf, q0, rho, cp, H)

        kold = np.copy(k)
        cpold = np.copy(cp)
        diffsOld = k / (cp * rho)
        cp = cTfunc(CONFIG['CT_RELATION_TYPE'], c0, CONFIG['CT_RELATION_PARAMS'], new_T)
        k = kTfunc(CONFIG['KT_RELATION_TYPE'], k0, CONFIG['KT_RELATION_PARAMS'], new_T)
	idx = k <= 0
	if sum(idx) > 0:
	    #k[idx] = kold[idx]
	    #print "#" + str(mpi_rank) + ", k == 0 at " + str(np.where(idx)) + ", it = " + str(it)
	    raise Exception("#" + str(mpi_rank) + ", k == 0 at " + str(np.where(idx)) + ", it = " + str(it) + ", new_T = " + str(new_T[idx])) 

        diffsChange = max(abs(k/(cp*rho)-diffsOld))

    T = np.copy(new_T)
    
    if it % CONFIG['OUTPUT_FILE_EVERY_TSTEP'] == 0 and not CONFIG['MPI']:
        print "curtime:", 1e-6 * STATUS['curTime'] / SECINYR, "\tdt:", 1e-6 * dt / SECINYR, "\t", "." * iterations

    if CONFIG['OUTPUT_FILE']:
        if it % CONFIG['OUTPUT_FILE_EVERY_TSTEP'] == 0:
            csvmetawriter.writerow([it, STATUS['curTime']/(SECINYR*1e6)])

            csvfile = open(Output_File_Path + "nodes." + str(it) , "wb")
            csvfilewriter = csv.writer(csvfile, delimiter=",", quotechar='"', quoting=csv.QUOTE_NONNUMERIC)
            csvfilewriter.writerow(["ix","t","x","T","k","cp","rho","H"])
            for ix in range(len(xs)):
                csvfilewriter.writerow([ix, 1e-6*STATUS['curTime']/SECINYR, xs[ix], T[ix], kold[ix], cpold[ix], rho[ix], H[ix]])
            csvfile.close()

    pe.getErosionSpeed(STATUS)

    STATUS['L'][0] += STATUS['Erosion_Speed'] * dt
    valueArrays = [T, new_T, H, c0, k0, rho]
    pe.remesh(STATUS, xs, STATUS['L'], valueArrays)

    STATUS['curTime'] += dt

    it += 1


if CONFIG['OUTPUT_FILE']:
    csvmetafile.close()

if CONFIG['MPI']:
    mpi_comm.Barrier()
    MPI.Finalize()
    sys.exit(0)

### analytical solutions:
# for constant H
# T_analytical = (1.0/C) * ((1+C*Tsurf)*np.exp((C/K0SINGLE)*(q0*xs-H0SINGLE*xs*xs/2.0))-1.0)
# for H decr exponentially:
if not CONFIG['MPI']:
    if CONFIG['KT_RELATION_TYPE'] == 1 and CONFIG['CT_RELATION_TYPE'] == 0 and \
                CONFIG['K0_TYPE'] == 0 and CONFIG['H0_TYPE'] == 2 and \
                STATUS['Erosion_Speed'] == 0 and not CONFIG['RESTART']:
        T_analytical = (1.0/CONFIG['KT_RELATION_PARAMS'][0]) * \
                   ((1+CONFIG['KT_RELATION_PARAMS'][0]*Tsurf)  *
                      np.exp((CONFIG['KT_RELATION_PARAMS'][0]/K0val)*(H0val*H_exp_param*H_exp_param*(1-np.exp(-xs/H_exp_param))-H0val*H_exp_param*xs+q0*xs))
                    - 1.0)
    else:
        T_analytical = xs * np.nan


    fig, (ax0, ax1) = plt.subplots(ncols=2, figsize=(8, 4))
    ax0.plot(T, -xs, 'ro-')
    ax0.plot(T_analytical, -xs, 'g-')
    ax0.plot(Tini, -xs, 'b-')
    ax1.plot(kold, -xs, 'ro-')
    plt.tight_layout()
    plt.show()

    print T

