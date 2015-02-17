import sys
import numpy as np
import pyearth_sub as pe
from matplotlib import pyplot as plt
import csv
import os

sys.path.insert(0, ".")
from config import CONFIG

#print sys.path

# **** CONFIGURATION ****

kTfunc = pe.kT
cTfunc = pe.cT

# **** END OF CONFIGURATION ****

STATUS = { 'CONFIG' : CONFIG }


if CONFIG['OUTPUT_FILE']:
    try:
        os.mkdir(CONFIG['OUTDIR'] + "/" + CONFIG['MODELNAME'])
    except Exception as e:
        if CONFIG['OUTPUT_OVERWRITE']:
            print "Error in mkdir, directory exists? Don't care, OUTPUT_OVERWRITE == True"
        else:
            raise e

    Output_File_Path = CONFIG['OUTDIR'] + "/" + CONFIG['MODELNAME'] + "/"

    csvmetafile = open(Output_File_Path + "meta", "wb")
    csvmetafile.close()

    csvmetafile = open(Output_File_Path + "meta", "ab")
    csvmetawriter = csv.writer(csvmetafile, delimiter=",", quotechar='"', quoting=csv.QUOTE_NONNUMERIC)

    csvmetawriter.writerow(['tstep', 'time_ma'])

if CONFIG['RESTART']:
    Restart_In_Path = CONFIG['RESTART_INDIR'] + "/" + CONFIG['RESTART_MODELNAME'] + "/"


STATUS['NX'] = CONFIG['NX']
STATUS['Erosion_Speed'] = CONFIG['EROSION_SPEED_M_MA'] / (1e6*365.25*24*60*60)
STATUS['Moho_Depth'] = CONFIG['MOHO_DEPTH_KM'] * 1e3
SECINYR = 60*60*24*365.25
STATUS['L'] = np.array((CONFIG['L_KM'][0]*1e3, CONFIG['L_KM'][1]*1e3))
STATUS['MaxTime'] = CONFIG['MAXTIME_MA'] * SECINYR * 1e6
STATUS['MaxRunTime'] = CONFIG['MAXRUNTIME_MA'] * SECINYR * 1e6

# mesh
xs = np.linspace(STATUS['L'][0], STATUS['L'][1], num=STATUS['NX'])
# to demonstrate irreg grid:  + 4e3*(np.random.rand(NX)-0.5)

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
else:
    raise Exception("Invalid K0_TYPE")

if CONFIG['C0_TYPE'] == 0:
    #constant c0
    c0 = 0.0 * xs + 800
elif CONFIG['C0_TYPE'] == 1:
    #different c0 for crust/mantle
    c0 = 0.0 * xs + 800
    c0[xs >= STATUS['Moho_Depth']] = 700
else:
    raise Exception("Invalid C0_TYPE")

if CONFIG['H0_TYPE'] == 0:
    # no heat prod
    H = 0.0 * xs
elif CONFIG['H0_TYPE'] == 1:
    # constant heat prod in crust and mantle
    H = 0.0 * xs
    H[xs < STATUS['Moho_Depth']] = 2.5e-6 # W m^-3
    H[xs >= STATUS['Moho_Depth']] = 3300 * 5e-12
elif CONFIG['H0_TYPE'] == 2:
    # exponential decrease
    H_exp_param = 30e3
    H0val = 2.5e-6
    H = H0val * np.exp(-xs / H_exp_param)
else:
    raise Exception("Invalid H0_TYPE")


# initial fields

if CONFIG['RESTART']:
    curtime = 0.0 #dummy

    Tin = np.array(())
    xsin = np.array(())

    rcsvmetafile = open(Restart_In_Path + "meta", "rb")
    csvmetareader = csv.reader(rcsvmetafile, delimiter=",", quotechar='"')

    rowfound = False
    for row in csvmetareader:
        if not row[0].isdigit():
            continue
        #print row[0], CONFIG['RESTART_TSTEP']
        if int(row[0]) == CONFIG['RESTART_TSTEP']:
            rowfound = True
            curtime = float(row[1]) * SECINYR * 1e6

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
            Tin = np.append(Tin, float(row[2]))
            xsin = np.append(xsin, float(row[1]))

    if len(xsin) == len(xs) and len(np.where(xsin == xs)) == STATUS['NX']:
        Tini = np.copy(Tin)
        # all nodes at same locations
    else:
        # nodes and/or extent in file read from restart file do not match those
        # configured in config file

        # changes in extent during restart not allowed
        if xsin[0] != xs[0] or xsin[-1] != xs[-1]:
            raise Exception("Changes in extent during restart not allowed.")
        else:
            # match grid points, result must follow gridding configured in config file
            valueArrays = [Tin]
            pe.remesh(STATUS, xsin, xs, valueArrays)
            Tini = np.copy(Tin)

    rcsvfile.close()
else:
    Tini = np.array([0.0] * CONFIG['NX'])
    curtime = 0.0

pe.initTemp(STATUS, Tini, xs)


if CONFIG['RESTART_POST_MOD'] == 0:
    pass # nothing
elif CONFIG['RESTART_POST_MOD'] == 1:
    STATUS['L'][0] -= 30e3  # raise the surface by 30km => overthrust
    valueArrays = [Tini, H, c0, k0, rho]
    #print Tini
    pe.remesh(STATUS, xs, STATUS['L'], valueArrays, extrapolation=2)
    #print Tini


# bnd cond
q0 = CONFIG['BND_BOT_HFLOW']  # 17e-3  # W m^-2
Tsurf = CONFIG['BND_TOP_TEMP'] # 20  # deg C

T = np.copy(Tini)
new_T = np.copy(Tini)

it = 0


plt.plot(Tini, -xs, 'ro-')
plt.show()

starttime = curtime

while (curtime < STATUS['MaxTime']) and (curtime-starttime < STATUS['MaxRunTime']):

    diffsChange = CONFIG['DIFFSCHANGE_ACCURACY'] + 1.0

    iterations = 0

    cp = cTfunc(CONFIG['CT_RELATION_TYPE'], c0, CONFIG['CT_RELATION_PARAMS'], T)
    k = kTfunc(CONFIG['KT_RELATION_TYPE'], k0, CONFIG['KT_RELATION_PARAMS'], T)

    while diffsChange > CONFIG['DIFFSCHANGE_ACCURACY']:
        iterations += 1
        dt = CONFIG['TSTEP_MULTI'] * pe.maxdt(k, rho, cp, xs)
        pe.diffstep(T, new_T, xs, k, dt, Tsurf, q0, rho, cp, H)

        kold = np.copy(k)
        cpold = np.copy(cp)
        diffsOld = k / (cp * rho)
        cp = cTfunc(CONFIG['CT_RELATION_TYPE'], c0, CONFIG['CT_RELATION_PARAMS'], new_T)
        k = kTfunc(CONFIG['KT_RELATION_TYPE'], k0, CONFIG['KT_RELATION_PARAMS'], new_T)

        diffsChange = max(abs(k/(cp*rho)-diffsOld))

    T = np.copy(new_T)
    
    STATUS['L'][0] += STATUS['Erosion_Speed'] * dt
    valueArrays = [T, H, c0, k0, rho]
    pe.remesh(STATUS, xs, STATUS['L'], valueArrays)

    curtime += dt

    if it % 1000 == 0:
        print "curtime:", 1e-6 * curtime / SECINYR, "\tdt:", 1e-6 * dt / SECINYR, "\t", "." * iterations

    if CONFIG['OUTPUT_FILE']:
        if it % CONFIG['OUTPUT_FILE_EVERY_TSTEP'] == 0:
            csvmetawriter.writerow([it, curtime/(SECINYR*1e6)])

            csvfile = open(Output_File_Path + "nodes." + str(it) , "wb")
            csvfilewriter = csv.writer(csvfile, delimiter=",", quotechar='"', quoting=csv.QUOTE_NONNUMERIC)
            csvfilewriter.writerow(["ix","x","T"])
            for ix in range(len(xs)):
                csvfilewriter.writerow([ix, xs[ix], T[ix]])
            csvfile.close()

    it += 1


if CONFIG['OUTPUT_FILE']:
    csvmetafile.close()

### analytical solutions:
# for constant H
# T_analytical = (1.0/C) * ((1+C*Tsurf)*np.exp((C/K0SINGLE)*(q0*xs-H0SINGLE*xs*xs/2.0))-1.0)
# for H decr exponentially:
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

