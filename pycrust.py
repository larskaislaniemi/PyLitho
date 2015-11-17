import numpy as np
import PrplxWrap
import sys
import pprint as pp
import matplotlib.pyplot as plt
import interpolate as ip
import pylitho_exceptions as pyle
import copy
import Grid
import Crust
from PhysicalConstants import *

# reload classes that might have been modified during developing
# and interactive sessions like IPython
reload(ip)
reload(PrplxWrap)
reload(pyle)
reload(Grid)
reload(Crust)


### MAIN PROGRAM:

config = { 
    'DT_MULTIPLIER' : 0.5,
    'T_ITER_THRESHOLD' : 1e-1,   # converge criteria for dT=f(dH) iteration (Kelvins)
    'GRID_NX' : 50,     # num of grid points (= num of elements + 1)
    'GRID_H' : 10000.0,  # original height of the model, meters
    'GRID_STATIC_POINT' : 0,   # which grid point is kept static (x doesn't change) during volume changes
    'BND_UPPER_PRESSURE' : 1370.0,           # pressure at the upper boundary, bar
    'BND_UPPER_TEMPERATURE' : 273.0 + 200.0, # temperature at the upper boundary, Kelvin
    'BND_LOWER_TEMPERATURE' : 273.0 + 600.0, # temperature at the lower boundary, Kelvin
    'RECORD_TYPE' : 1,
#    'EVENT_TYPE' : [2, 3], 
#    'EVENT_TIMING' : [1e3 * SECINYR, 500e3 * SECINYR],
    'EVENT_TYPE' : [1],
    'EVENT_TIMING' : [1e3 * SECINYR],
    'INI_FIELD_TYPE' : 1,
}


cr = Crust.Crust("crustmod", config=config)
cr.updatePerplex()
cr.updateCpRho()
cr.updateGridToRho()
cr.updatePressure()
cr.initRecord()

eventDone = [False] * len(cr.config['EVENT_TYPE'])

while cr.time < 1e6 *SECINYR:
    cr.diffuseT()
    cr.updateCpRho()
    cr.updateGridToRho()
    cr.updatePressure()
    cr.doRecord()

    plt.plot(cr.T_e, -cr.xs_e)
    plt.hold(True)
    
    for ievent in range(len(eventDone)):
        if (not eventDone[ievent]) and cr.time > cr.config['EVENT_TIMING'][ievent]:
            eventType = cr.config['EVENT_TYPE'][ievent]
            print " *** Event ", ievent, ", type ", eventType
            eventDone[ievent] = True
            if eventType == 1:
                ne_add = 3  # num of elems to add (and of grid points to add)
                h_add = 500.  # total length to add
                loc_add = int(cr.nx/2) # grid point to add to 
                print "loc_add:", loc_add
                new_xs = np.zeros(cr.nx+ne_add)
                print new_xs
                new_xs[0:loc_add] = cr.xs[0:loc_add]
                new_xs[(loc_add+ne_add):(cr.nx+ne_add)] = cr.xs[loc_add:cr.nx] + h_add
                new_xs[loc_add:(loc_add+ne_add)] = cr.xs[loc_add] + np.arange(ne_add) * h_add / float(ne_add)
                new_nx = len(new_xs)
                new_xs_e = (new_xs[1:new_nx]+new_xs[0:(new_nx-1)]) * 0.5
                new_ne = len(new_xs_e)

                new_T_val = 1200. # K
                cr.T_e = np.resize(cr.T_e, new_ne)
                print loc_add, ne_add, cr.ne, new_ne
                cr.T_e[(loc_add+ne_add):] = cr.T_e[loc_add:cr.ne]  # copy old values forward
                cr.T_e[loc_add:(loc_add+ne_add)] = new_T_val    # fill the gap

                new_k_val = 3.5
                cr.k_e = np.resize(cr.k_e, new_ne)
                cr.k_e[(loc_add+ne_add):] = cr.k_e[loc_add:cr.ne]
                cr.k_e[loc_add:(loc_add+ne_add)] = new_k_val

                new_compo_val = np.array([0.0327, 0.0248, 0.1540, 0.6662, 0.0280, 0.0359, 0.0064, 0.0010, 0.0504, 0.0100])
                cr.C_e = np.resize(cr.C_e, (new_ne, cr.ncomponents))
                for i in range(cr.ne-1, loc_add, -1):
                    cr.C_e[i+ne_add][:] = cr.C_e[i][:] 
                for i in range(ne_add):
                    cr.C_e[loc_add+i][:] = new_compo_val[:]

                cr.xs = np.resize(cr.xs, new_nx)
                cr.xs[:] = new_xs[:]
                cr.xs_e = np.resize(cr.xs_e, new_ne)
                cr.xs_e[:] = new_xs_e[:]
                cr.nx = new_nx
                cr.ne = new_ne

                # only resize needed, no copying of values
                cr.Cp_e = np.resize(cr.Cp_e, new_ne)
                cr.rho_e = np.resize(cr.rho_e, new_ne)
                cr.pres_e = np.resize(cr.pres_e, new_ne)
                cr.pres = np.resize(cr.pres, new_nx)
                cr.mass_e = np.resize(cr.mass_e, new_ne)

                print "pressure in elements was:"
                print cr.pres_e
                print cr.xs_e
                cr.updatePerplex(initial=True)
                cr.updateGridToRho()
                cr.updatePressure()
                print "pressure in elements is:"
                print cr.pres_e
                print cr.xs_e
            elif eventType == 2:
                cr.bnd_lower_temperature = 273. + 900.
            elif eventType == 3:
                cr.bnd_lower_temperature = 273. + 600.
            else:
                raise Exception("Unknown event type")


#plt.show()

#allrecdpt = np.array(cr.record_depths)
#allrecval = np.array(cr.record_values)
#allreclab = cr.record_valuesources.keys()
#allrectime = np.array(cr.record_times)

uniquephases = []
for recrec in cr.recdata:
    for i in range(recrec.grid.nx):
        for phase in recrec.data["!ph"][i]:
            if phase in uniquephases:
                pass
            else:
                uniquephases.append(phase)
uniquephases.sort()

# TODO account for immiscibility (Fsp, etc )
# currently just mixed
stablephases = []
for rec in cr.recdata:
    stablephases.append(Grid.GriddedData(rec.grid.nx, rec.grid.xs))
    stablephases[-1].addMetaData("time", rec.metadata["time"])
    for aphase in uniquephases:
        stablephases[-1].addData(aphase, [0] * rec.grid.nx)
    for idepth in range(rec.grid.nx):
        for aphase in uniquephases:
            phaselocs = [i for i, j in enumerate(rec.data["!ph"][idepth]) if j == aphase]
            for i in phaselocs:
                # add together similar phases --> immiscibility disappears and are treated
                # as mixed phases
                stablephases[-1].data[aphase][idepth] += rec.data["!wp"][idepth][i]

#### plot single value fields
plt.close('all')

plt.figure()
depthres = 60
times = np.zeros(len(stablephases))
xs = np.linspace(0,35e3,depthres)
recfields = ['T', 'rho', 'Vs', 'Vp']
plotsqr = int(len(recfields)**0.5-1e-12)+1
TDW = np.zeros((len(recfields), depthres, len(cr.recdata)))
for i in range(len(cr.recdata)):
    sys.stdout.write(" " + str(i) + "\r")
    sys.stdout.flush()
    for ifield in range(len(recfields)):
        fieldvals = np.array(ip.interpolate(cr.recdata[i].grid.xs, cr.recdata[i].data[recfields[ifield]], xs))
        TDW[ifield,:,i] = fieldvals
    times[i] = cr.recdata[i].metadata["time"]

XS, YS = np.meshgrid(times/(SECINYR*1e3), -xs)
for i in range(len(recfields)):
    plt.subplot(plotsqr,plotsqr,i)
    CS = plt.contourf(XS, YS, TDW[i,:,:], 30)
    cbar = plt.colorbar(CS)
    plt.title(recfields[i])

#### plot phase assemblages and weight percentages 

plt.figure()
depthres = 60
TDA = np.zeros((depthres, len(stablephases)))
times = np.zeros(len(stablephases))
xs = np.linspace(0,35e3,depthres)
considerphases = ['Bio', 'Gt', 'Melt', 'q']
plotsqr = int(len(considerphases)**0.5-1e-12)+1
TDW = np.zeros((len(considerphases), depthres, len(stablephases)))
for i in range(len(stablephases)):
    sys.stdout.write(" " + str(i) + "\r")
    sys.stdout.flush()
    iasm = 1
    phaseex = xs * 0.0
    iphase = 0
    for aphase in considerphases:
        phasewts = np.array(ip.interpolate(stablephases[i].grid.xs, stablephases[i].data[aphase], xs, notFoundVal=0.0))
        TDW[iphase,:,i] = phasewts
        phaseex = phaseex + (phasewts > 0) * float(iasm)
        iasm = iasm * 2
        iphase = iphase + 1
    TDA[:,i] = phaseex
    times[i] = stablephases[i].metadata["time"]

XS, YS = np.meshgrid(times/(60.*60.*24.*365.*1e3), -xs)
#CS = plt.contourf(XS, YS, TDA)
#cbar = plt.colorbar(CS)
#plt.show()

for i in range(len(considerphases)):
    plt.subplot(plotsqr,plotsqr,i)
    CS = plt.contourf(XS, YS, TDW[i,:,:])
    cbar = plt.colorbar(CS)
    plt.title(considerphases[i])


plt.show()
#considerphases = ['Melt', 'Bio', 'Gt']
#assemblages = []
#for rec in stablephases:
#    assemblages.append(Grid.GriddedData(rec.grid.nx, rec.grid.xs))
#    for aphase in considerphases:
#        assemblages[-1].addData(aphase, [0] * rec.grid.nx)
#
#    for idepth in range(rec.grid.nx):
#        for aphase in considerphases:
#            if aphase not in rec.data.keys()
#                raise Exception("requested a phase that is nowhere stable")
#            if rec.data[aphase][idepth] > 0.0:
#                assemblages
#

print "You can do plt.show()"
x = "."
while len(x) == 0 or x[0] != "!":
    x=raw_input(":")
    try:
        exec(x)
    except Exception as e:
        print e

mineral = 'Melt'
for it in range(0,len(cr.recdata),10):
    x=cr.recdata[it].grid.xs
    y1=cr.recdata[it].data["T"]
    y2=stablephases[it].data[mineral]
    plt.plot(x,y1)
    plt.twinx()
    plt.plot(x,y2)
    plt.show()

#plt.close('all')
#f, (ax1, ax2) = plt.subplots(1,2,sharey=True)
#ax1.plot(cr.T_e, -cr.xs_e)
#ax2.plot(cr.Cp_e, -cr.xs_e)
#plt.show()


