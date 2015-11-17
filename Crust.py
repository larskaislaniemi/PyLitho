import numpy as np
import PrplxWrap
import sys
import pprint as pp
import interpolate as ip
import pylitho_exceptions as pyle
import copy
import Grid
from PhysicalConstants import *

reload(ip)
reload(PrplxWrap)
reload(pyle)
reload(Grid)

class Container(object):
    @property
    def T(self):
        return self._T
    @T.setter
    def T(self, value):
        self._T = value
        self.perplexResult = None
    @T.deleter
    def T(self):
        self._T = None
        self.perplexResult = None

    @property
    def P(self):
        return self._P
    @P.setter
    def P(self, value):
        self._P = value
        self.perplexResult = None
    @P.deleter
    def P(self):
        self._P = None
        self.perplexResult = None
    

    def __init__(self, copyfrom=None, compos=None, masses=None, perplex=None):
        if copyfrom is None:
            if compos is None:
                self.components = []
                self.masses = []
                self.ncomponents = 0
                self._T = None
                self._P = None
                self._H = None
            else:
                self.components = compos[:]
                self.ncomponents = len(self.components)
                if masses is None:
                    self.masses = [0] * self.ncomponents
                else:
                    if self.ncomponents != len(masses):
                        raise pyle.PyLerr_Undefined("len(components) != len(masses)")
                    self.masses = masses[:]
                self._T = None
                self._P = None
                self._H = None
        else:
            if not isinstance(copyfrom, Container):
                raise pyle.PyLerr_TypeError("Argument is not instance of Container")
            self.components = copy.components[:]
            self.masses = copyfrom.masses[:]
            self.ncomponents = copyfrom.ncomponents
            self._T = copyfrom._T
            self._P = copyfrom._P
            self._H = copyfrom._H
        self.perplex = perplex
        self.perplexResult = None

    #T = property(get_T, set_T, del_T, "Temperature")
    #P = property(get_P, set_P, del_P, "Pressure")
    #self.H = property(get_H, set_H, del_H, "Enthalpy")

    def updatePerplex(self, perplex=None):
        if perplex is None:
            if self.perplex is None:
                raise pyle.PyLerr_Undefined("No valid perplex module")
            else:
                perplex = self.perplex
        else:
            self.perplex = perplex

        if set(self.components) != set(perplex.callerComponents) or len(self.components) != len(perplex.callerComponents):
            raise pyle.PyLerr_Undefined("component set mismatch")

        if self._P is None or self._T is None:
            raise pyle.PyLerr_Undefined("P or T not set")

        idxs = [self.components.index(a) for a in perplex.callerComponents]

        masses = np.array(self.masses, dtype=float)
        masses = masses/np.sum(masses)
        masses = masses[idxs]

        self.perplexResult = None
        self.perplexResult = perplex.phaseq(self._P, self._T, masses.tolist(), debug=False)

    def getVolume(self, perplex=None):
        if self.perplexResult is None:
            self.updatePerplex(perplex=perplex)
        rho = self.perplexResult['SYSPROP'][self.perplex.syspropnum['rho']]
        return np.sum(self.masses)/rho

    def move_isentropic(self, dP, perplex=None):
        if self.perplexResult is None:
            self.updatePerplex(perplex=perplex)

        if self._P is None or self._T is None:
            raise pyle.PyLerr_Undefined("T/P not defined")

        maxerr = 0.01
        curerr = 1.1*maxerr

        T0 = self.T
        P0 = self.P
        V0 = self.perplexResult['SYSPROP'][self.perplex.syspropnum['V']]
        Cp0 = self.perplexResult['SYSPROP'][self.perplex.syspropnum['Cp']]
        S0 = self.perplexResult['SYSPROP'][self.perplex.syspropnum['S']]
        H0 = self.perplexResult['SYSPROP'][self.perplex.syspropnum['H']]

        dH = V0 * dP  # q=0, only expansion work

        dT = dH / Cp0
        self.T = T0 + dT
        self.P = P0 + dP
        self.updatePerplex(perplex=perplex)
        S1 = self.perplexResult['SYSPROP'][self.perplex.syspropnum['S']]
        Cp1 = self.perplexResult['SYSPROP'][self.perplex.syspropnum['Cp']]

        dSph = S1-S0  # entropy mismatch (entropy  change of phase changes)
        dHadj = (-dSph) * self.T  # to gain it back, modify enthalpy
        dTadj = dHadj / Cp1     # calculate corresponding temp change
        self.T = self.T + dTadj
        self.updatePerplex(perplex=perplex)
        S1 = self.perplexResult['SYSPROP'][self.perplex.syspropnum['S']]
        Cp1 = self.perplexResult['SYSPROP'][self.perplex.syspropnum['Cp']]

    def move_adiab(self, dPtot, dP=3e6*1e-5, perplex=None):
        # Assume adiabatic AND isentropic process
        #
        # dS = 0
        # dH = V dp
        # 
        # 1. Read current H (H_1), S (S_1) and V (V_1) at T_1, P_1
        # 2. Calc dH = V dP
        # 3. Calc H_2 = H_1 + dH
        # 4. Find T_2, P_2==P_1 where H==H_2
        # 5. Check that S_2 == S_1
        #
        # 1. Estimate dH: dH = V dP
        # 2. Estimate dT: dT = dH / Cp 
        # 3. 

        method = 2

        if self.perplexResult is None:
            self.updatePerplex(perplex=perplex)

        if self._P is None or self._T is None:
            raise pyle.PyLerr_Undefined("T/P not defined")

        maxdT = 5.
        if dPtot < 0:
            dPdir = -1.
            if dP > 0:
                dP = -dP
        else:
            dPdir = 1.

        dPsum = 0.0
        dTsum = 0.0

        T0 = self._T
        rho0 = self.perplexResult['SYSPROP'][self.perplex.syspropnum['rho']] # kg/m3
        Cp0 = 1000. * self.perplexResult['SYSPROP'][self.perplex.syspropnum['Cp']] / self.perplexResult['SYSPROP'][self.perplex.syspropnum['N']] # J K^-1 kg^-1
        alpha0 = self.perplexResult['SYSPROP'][self.perplex.syspropnum['alpha']]
        #V0 = sum(self.masses) / rho0  
        V0 = 1000. * self.perplexResult['SYSPROP'][self.perplex.syspropnum['V']] / self.perplexResult['SYSPROP'][self.perplex.syspropnum['N']] # J/bar/kg
        H0 = 1000. * self.perplexResult['SYSPROP'][self.perplex.syspropnum['H']] / self.perplexResult['SYSPROP'][self.perplex.syspropnum['N']] # J/kg


        origP = self.P
        origT = T0
        origV = V0
        origH = H0
        origCp = Cp0

        while abs(dPtot - dPsum) > dPdir*dPtot/1e3:
            dP0 = dPdir * min(dPdir * dP, dPdir * dPtot - dPdir * dPsum)
            if method == 1:
                dT0 = T0 * alpha0 * dP0*1e5 / (Cp0 * rho0) 
            elif method == 2:
                dH0 = dP0 * V0
                dT0 = dH0 / Cp0
            self.T = self.T + dT0
            self.P = self.P + dP0
            dPsum = dPsum + dP0
            print dP0, dPsum, dT0, self.T, self.P
            self.updatePerplex(perplex=perplex)
            T0 = self.T
            Cp0 = 1000. * self.perplexResult['SYSPROP'][self.perplex.syspropnum['Cp']] / self.perplexResult['SYSPROP'][self.perplex.syspropnum['N']]
            alpha0 = self.perplexResult['SYSPROP'][self.perplex.syspropnum['alpha']]
            rho0 = self.perplexResult['SYSPROP'][self.perplex.syspropnum['rho']]
            V0 = 1000. * self.perplexResult['SYSPROP'][self.perplex.syspropnum['V']] / self.perplexResult['SYSPROP'][self.perplex.syspropnum['N']]
            if method == 2:
                H1 = 1000. * self.perplexResult['SYSPROP'][self.perplex.syspropnum['H']] / self.perplexResult['SYSPROP'][self.perplex.syspropnum['N']]
                print "H mismatch: ", (H1-H0-dH0), (H1-H0), dH0
                H0 = H1

            print self.perplexResult['NAMEPHASES']
            if "melt(HP)" in self.perplexResult['NAMEPHASES']:
                print "Melt!"

        finalV = 1000. * self.perplexResult['SYSPROP'][self.perplex.syspropnum['V']] / self.perplexResult['SYSPROP'][self.perplex.syspropnum['N']]
        finalH = 1000. * self.perplexResult['SYSPROP'][self.perplex.syspropnum['H']] / self.perplexResult['SYSPROP'][self.perplex.syspropnum['N']]
        finalCp = 1000. * self.perplexResult['SYSPROP'][self.perplex.syspropnum['Cp']] / self.perplexResult['SYSPROP'][self.perplex.syspropnum['N']]
        print "dPtot:", dPtot
        print "P, T, V, H, Cp:"
        print "orig: ", origP, origT, origV, origH, origCp
        print "final:", self.P, self.T, finalV, finalH, finalCp

        

    def addComponent(self, component, mass, overwrite=False):
        if component in self.components:
            idx = self.components.index(component)
            if overwrite:
                masses[idx] = mass
            else:
                masses[idx] = masses[idx] + mass
        else:
            self.components.append(component)
            self.components.append(mass)
        self.perplexResult = None

class Crust:

    def __init__(self, datfile, config=None):
        if config is None:
            self.config = {
                'DT_MULTIPLIER' : 0.5,
                'T_ITER_THRESHOLD' : 1e-1,   # converge criteria for dT=f(dH) iteration (Kelvins)
                'GRID_NX' : 35,     # num of grid points (= num of elements + 1)
                'GRID_H' : 35000.0,  # original height of the model, meters
                'GRID_STATIC_POINT' : 0,   # which grid point is kept static (x doesn't change) during volume changes
                'BND_UPPER_PRESSURE' : 3000.0,           # pressure at the upper boundary, bar
                'BND_UPPER_TEMPERATURE' : 273.0 + 300.0, # temperature at the upper boundary, Kelvin
                'BND_LOWER_TEMPERATURE' : 273.0 + 600.0, # temperature at the lower boundary, Kelvin
                'RECORD_TYPE' : 1,
                'EVENT_TYPE' : [2, 3],
                'EVENT_TIMING' : [1e3 * SECINYR, 500e3 * SECINYR],
                'INI_FIELD_TYPE' : 1,
            }
        else:
            self.config = config

        self.time = 0.0    # model time
        self.last_dt = 0.0 # latest time step taken
        self.timestep = 0  # counter for time steps

        # (initial) grid specification
        self.nx = self.config['GRID_NX']  # num of grid points
        self.ne = self.nx-1               # num of elements
        self.xs = np.linspace(0.0,self.config['GRID_H'],self.nx)        # xs at grid points
        self.xs_e = 0.5 * (self.xs[0:self.nx-1] + self.xs[1:self.nx])   # xs at elements

        # boundary condition stuff
        self.bnd_upper_pressure = self.config['BND_UPPER_PRESSURE']
        self.bnd_upper_temperature = self.config['BND_UPPER_TEMPERATURE']
        self.bnd_lower_temperature = self.config['BND_LOWER_TEMPERATURE']

        # physical constants
        self.accel_g = 9.81

        # (initial) field specification
        # "_e" refers to elements (grid points in the middle of elements);
        # without '_e': the main grid points, defining (e.g.) upper and lower boundaries
        # of the model and placed also in between lithological units.
        # Num of main grid points (self.nx) = num of element grid points (self.ne) + 1
        # Most of the data is stored in the elements.
        self.T_e = self.getIniField("temp", self.xs_e) # K
        self.C_e = self.getIniField("compo", self.xs_e)           # wt% (0..1)
        self.k_e = self.getIniField("conductivity", self.xs_e)    #
        self.Cp_e = self.getIniField("heatcapacity", self.xs_e)   # J/kgK
        self.rho_e = self.getIniField("density", self.xs_e)       # kg/m3
        self.pres_e = self.getIniField("pressure", self.xs_e)     # bar
        self.pres = self.getIniField("pressure", self.xs)
        self.mass_e = self.getIniField("mass", self.xs_e)         # kg
        #self.T = self.getIniField("temp", self.xs)               # K # seems not to be needed....

        # initiate perplex
        self.components = self.getIniField("componames")
        self.perplex = PrplxWrap.perplex(datfile, callerCompnames=self.components)
        self.ncomponents = len(self.components)

        # latest perplex results
        self.perplexResult = None

        # flag to indicate whether perplex data calculated
        # for the grid points is up to date
        self.perplexOK = False

        # generate initial info by perplex
        self.updatePerplex(initial=True)  # initial call also adjusts Cp/Rho and pressure and mass

        # update grid according to new density (=volume) values from PerpleX
        self.updateGridToRho()  # in principle, after initial call to updatePerpleX(), this shouldn't do anything


    def output(self):
        pass

    def updatePressure(self):
        # NB! and TODO:
        # Pressure changes involve change in enthalpy (dH=VdP)
        # which is not currently calculated
        self.pres[:] = self.bnd_upper_pressure
        self.pres_e[:] = self.bnd_upper_pressure
        for i in range(1,self.ne):
            self.pres[i] += sum(self.mass_e[0:i]) * self.accel_g * 1e-5
            self.pres_e[i] = self.pres[i]   # simplification: the pressure of the element is the pressure
                                            # at the upper surface of the element

    def updatePerplex(self, ielem=-1, initial=False):
        if ielem >= 0 and initial:
            raise Exception("ielem >=0 and initial == True in updatePerplex()")
        if ielem < 0:
            self.perplexResult = []
        self.updatePressure()
        for i in range(self.ne):
            if ielem < 0:
                self.perplexResult.append(self.perplex.phaseq(self.pres_e[i], self.T_e[i], self.C_e[i]))
                testval = np.sum(np.sum(np.array(self.perplexResult[-1]['WTPHASES'])))
                if (self.perplexResult[-1])['RETVAL'] != 0 or testval != testval:
                    pp.pprint(self.perplexResult[-1], stream=sys.stderr)
                    raise Exception("Ooo..ps")
                self.perplexOK = True
                if initial:
                    self.updateCpRho(i)
                    self.mass_e[i] = self.rho_e[i] * (self.xs[i+1] - self.xs[i])
                    self.updatePressure()
            elif ielem == i:
                self.perplexResult[ielem] = self.perplex.phaseq(self.pres_e[i], self.T_e[i], self.C_e[i])
                testval = np.sum(np.sum(np.array(self.perplexResult[ielem]['WTPHASES'])))
                if (self.perplexResult[ielem])['RETVAL'] != 0 or testval != testval:
                    pp.pprint(self.perplexResult[ielem], stream=sys.stderr)
                    raise Exception("Ooo..ps")

    def updateCpRho(self, ielem=-1):
        #if not self.perplexOK:
        #    self.updatePerplex()
        if ielem >= 0:
            doElems = [ielem]
        else:
            doElems = range(self.ne)
        for i in doElems:
            self.rho_e[i] = self.perplexResult[i]['SYSPROP'][self.perplex.syspropnum['rho']]
            self.Cp_e[i] = G2KG * self.perplexResult[i]['SYSPROP'][self.perplex.syspropnum['Cp']] / self.perplexResult[i]['SYSPROP'][self.perplex.syspropnum['N']]

    def updateGridToRho(self):
        volumes = self.mass_e / self.rho_e  # in 1D this is directly dx
        new_x = self.xs[:] * 0.0
        for i in range(0,self.config['GRID_STATIC_POINT']):
            new_x[i] = self.xs[self.config['GRID_STATIC_POINT']] - sum(volumes[i:self.config['GRID_STATIC_POINT']])
        for i in range(self.config['GRID_STATIC_POINT'],self.nx):
            new_x[i] = self.xs[self.config['GRID_STATIC_POINT']] + sum(volumes[self.config['GRID_STATIC_POINT']:i])
        self.xs[:] = new_x[:]
        print " * updateGridToRho(): New extents are ", new_x[0], new_x[-1]

    def maxdt(self):
        # return maximum time step for diffusion,
        # assume that information for k, rho and Cp is up-to-date
        dx = self.xs[1:self.nx] - self.xs[0:(self.nx-1)]
        diff = self.k_e / (self.rho_e * self.Cp_e)
        maxdt = min(0.5 * dx * dx / diff)
        return maxdt

    def addEnthalpy(self, ielem, dH):
        # input:
        #  ielem: element to which the enthalpy is added
        #  dH:    amount of added enthalpy [J]

        #if not self.perplexOK:
        #    self.updatePerplex()  # to make sure pressure field is OK.
                                  # probably unnecessary.

        doNewIteration = True

        # dH is given in Joules (total), transform to J/kg
        dH = dH / (self.rho_e[ielem] * (self.xs[ielem+1] - self.xs[ielem]))

        Cp0 = G2KG * self.perplexResult[ielem]['SYSPROP'][self.perplex.syspropnum['Cp']] / self.perplexResult[ielem]['SYSPROP'][self.perplex.syspropnum['N']]
        T0 = self.T_e[ielem]
        H0 = G2KG * self.perplexResult[ielem]['SYSPROP'][self.perplex.syspropnum['H']] / self.perplexResult[ielem]['SYSPROP'][self.perplex.syspropnum['N']]
        Tini = T0
        Hini = H0
        Cpini = Cp0

        niter = 0
        stepDecr = 1e10   # decrease in solution between successive steps
        stepMultip = 1.0

        while doNewIteration:
            T1 = T0 - stepMultip * (H0-Hini-dH) / Cp0
            self.T_e[ielem] = T1
            self.updatePerplex(ielem)
            testval = self.perplexResult[ielem]['SYSPROP'][self.perplex.syspropnum['Cp']]
            if testval != testval:
                pp.pprint(self.perplexResult[ielem], stream=sys.stderr)
                sys.exit(0)
            H1 = G2KG * self.perplexResult[ielem]['SYSPROP'][self.perplex.syspropnum['H']] / self.perplexResult[ielem]['SYSPROP'][self.perplex.syspropnum['N']]
            Cp1 = G2KG * self.perplexResult[ielem]['SYSPROP'][self.perplex.syspropnum['Cp']] / self.perplexResult[ielem]['SYSPROP'][self.perplex.syspropnum['N']]

            if abs(T1-T0) < self.config['T_ITER_THRESHOLD']:
                doNewIteration = False
            else:
                if abs(T1-T0) >= stepDecr:
                    # halve the step
                    stepMultip = stepMultip / 2.0
                stepDecr = abs(T1-T0)
                niter = niter + 1
                if niter > 10:
                    sys.stdout.write(" WARN niter = " + str(niter) + "(T1-T0 =" + str(T1-T0) + ")\r")
                T0 = T1
                H0 = H1
                Cp0 = Cp1

        #if niter > 0:
        #    print "Needed", niter, "iterations"
        T_react = T1 - Tini - dH/Cpini
        print "dH = ", dH, "; dT = ", (T1-Tini), "; Cp = ", Cp1, ", Treact = ", T_react

        return T1-Tini


    def diffuseT(self, dt = 0.0):
        self.timestep = self.timestep + 1
        print " * Diffusion, time step", self.timestep, ", time = ", (self.time / SECINYR), " yrs"

        T_e = self.T_e

        if dt <= 0:
            # estimate largest time step, use a fraction of that
            dt = self.config['DT_MULTIPLIER'] * self.maxdt()
            self.last_dt = dt

        # calculate conductivity for main grid points by averaging from elements (wgtd avg)
        # boundaries excluded
        k_g =       0.5 * (self.xs[2:self.nx] - self.xs[1:(self.nx-1)]) * self.k_e[1:self.ne]
        k_g = k_g + 0.5 * (self.xs[1:(self.nx-1)] - self.xs[0:(self.nx-2)]) * self.k_e[0:(self.ne-1)]
        k_g = k_g / (0.5 * (self.xs[2:self.nx] - self.xs[0:(self.nx-2)]))
        n_kg = k_g.size

        ## heat diffusion for
        # ... internal grid points.
        # At elemental grid points we define T, rho, Cp, k, d2T/dz2;
        # at main grid points we define dT/dz.
        dT1 = k_g[1:n_kg]     * (T_e[2:self.ne]     - T_e[1:(self.ne-1)]) / (self.xs_e[2:self.ne]     - self.xs_e[1:(self.ne-1)])
        dT2 = k_g[0:(n_kg-1)] * (T_e[1:(self.ne-1)] - T_e[0:(self.ne-2)]) / (self.xs_e[1:(self.ne-1)] - self.xs_e[0:(self.ne-2)])
        d2T = (dT1 - dT2) / (0.5 * (self.xs_e[2:self.ne] - self.xs_e[0:(self.ne-2)]))
        #DTinternal = d2T * dt / (self.rho_e[1:(self.ne-1)] * self.Cp_e[1:(self.ne-1)])
        DHinternal = d2T * dt   # J/m3
        DHinternal = DHinternal * (self.xs[2:(self.nx-1)] - self.xs[1:(self.nx-2)])  # = Cp * m * dT = dH

        # ... uppermost grid point
        dT1 = k_g[0]     * (T_e[1]     - T_e[0]) / (self.xs_e[1]     - self.xs_e[0])
        dT2 = self.k_e[0] * (T_e[0] - self.bnd_upper_temperature) / (self.xs_e[0]-self.xs[0])
        d2T = (dT1 - dT2) / (0.5 * (self.xs_e[1] - (-self.xs_e[0] + self.xs[0])))
        #DTupper = d2T * dt / (self.rho_e[0] * self.Cp_e[0])
        DHupper = d2T * dt
        DHupper = DHupper * (self.xs[1] - self.xs[0])

        # ... lowermost grid point
        dT1 = self.k_e[self.ne-1] * (self.bnd_lower_temperature - T_e[self.ne-1]) / (self.xs[self.nx-1]-self.xs_e[self.ne-1])
        dT2 = k_g[n_kg-1] * (T_e[self.ne-1] - T_e[self.ne-2]) / (self.xs_e[self.ne-1] - self.xs_e[self.ne-2])
        d2T = (dT1 - dT2) / (0.5 * (self.xs[self.nx-1] + (self.xs[self.nx-1]-self.xs_e[self.ne-1]) - self.xs_e[self.ne-2]))
        #DTlower = d2T * dt / (self.rho_e[self.ne-1] * self.Cp_e[self.ne-1])
        DHlower = d2T * dt
        DHlower = DHlower * (self.xs[self.nx-1]-self.xs[self.nx-2])

        self.addEnthalpy(0,DHupper)
        for ielem in range(1,self.ne-1):
            if (DHinternal[ielem-1] != 0.0):
                self.addEnthalpy(ielem, DHinternal[ielem-1])
        self.addEnthalpy(self.ne-1,DHlower)
        self.perplexOK = True  # ... since we just called addEnthalpy on all elements

        self.time = self.time + self.last_dt
        print "\tdt used: ", self.last_dt/SECINYR, "yrs"


    def getIniField(self, field, xs=[0]):
        fieldType = self.config['INI_FIELD_TYPE']

        if field == "temp":
            if fieldType in [0,1]:
                #retField = xs * 0.0 + 1000.0 + 273.0
                retField = self.bnd_upper_temperature + (self.bnd_lower_temperature-self.bnd_upper_temperature) * (xs-xs[0]) / (xs[-1]-xs[0])
            else:
                raise Exception("Invalid fieldType")
        elif field == "componames":
            if fieldType in [0,1]:
                retField = ['NA2O',
                            'MGO',
                            'AL2O3',
                            'SIO2',
                            'K2O',
                            'CAO',
                            'TIO2',
                            'MNO',
                            'FEO',
                            'H2O']
            else:
                raise Exception("Invalid fieldType")
        elif field == "compo":
            if fieldType in [0,1]:
                retField = []
                for i in range(len(xs)):
                    retField.append([[]] * 10)
                    retField[i][0] = 0.0327
                    retField[i][1] = 0.0248
                    retField[i][2] = 0.1540
                    retField[i][3] = 0.6662
                    retField[i][4] = 0.0280
                    retField[i][5] = 0.0359
                    retField[i][6] = 0.0064
                    retField[i][7] = 0.0010
                    retField[i][8] = 0.0504
                    retField[i][9] = 0.0100
                    if fieldType == 1 and xs[i] > 0.9*xs[-1]:
                        retField[i][9] = 0.0000 
            else:
                raise Exception("Invalid fieldType")
        elif field == "conductivity":
            if fieldType in [0,1]:
                retField = 0.0 * xs + 3.5
            else:
                raise Exception("Invalid fieldType")
        elif field == "heatcapacity":
            retField = 0.0 * xs
            # (this is a dummy, need perplex to do properly)
        elif field == "mass":
            retField = 0.0 * xs
            # (this is a dummy, need perplex to do properly)
        elif field == "density":
            retField = 0.0 * xs
            # (this is a dummy, need perplex to do properly)
        elif field == "pressure":
            retField = 0.0 * xs
            # (this is a dummy, need perplex to do properly)
        else:
            raise Exception("Invalid field")

        return retField

    def initRecord(self):
        if self.config['RECORD_TYPE'] == 1:
            self.recdata = []
            self.recdata_depths = [] # these will be update on every call to doRecord()
                                     # (for RECORD_TYPE 1)
            self.recdata_sources = { 'T'   : "self.T_e",
                                     'Cp'  : "self.Cp_e",
                                     'rho' : "self.rho_e",
                                     'Vp'  : "[dic['SYSPROP'][self.perplex.syspropnum['Vp']] for dic in self.perplexResult[:]]",
                                     'Vs'  : "[dic['SYSPROP'][self.perplex.syspropnum['Vs']] for dic in self.perplexResult[:]]",
                                     '!ph' : "[dic['NAMEPHASES'] for dic in self.perplexResult[:]]",
                                     '!wp' : "[dic['WTPHASES'] for dic in self.perplexResult[:]]",
                                   }  # N.B. these must be values at the elements, *_e
        else: 
            raise Exception("Invalid RECORD_TYPE")

    def doRecord(self):
        if self.config['RECORD_TYPE'] == 1:
            # we'll always record the data at elements
            self.recdata_depths = list(self.xs_e)

            # add new record grid for this timestep
            self.recdata.append(Grid.GriddedData(self.ne, self.xs_e))
            self.recdata[-1].addMetaData("time", self.time)

            # add data to the grid
            for dataname in self.recdata_sources.keys():
                if dataname[0] != "!":
                    # interpolatable data
                    self.recdata[-1].addData(dataname, ip.interpolate(self.xs_e, eval(self.recdata_sources[dataname]), self.recdata_depths))
                else:
                    # non-interpolatable data
                    self.recdata[-1].addData(dataname, eval(self.recdata_sources[dataname]))


