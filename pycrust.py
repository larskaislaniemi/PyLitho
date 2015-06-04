import numpy as np
import PrplxWrap

G2KG = 1000.0
SECINYR = 60.0*60.0*24.0*365.25

class Crust:
    perplex = PrplxWrap.perplex("in26klb_new")

    def __init__(self):
        self.config = {
            'DT_MULTIPLIER' : 0.5,
            'DT_MULTIPLIER_MAX' : 0.9,
            'CP_ITER_THRESHOLD' : 1e-4,  # converge criteria for Cp iteration,
                                         # (Cp_new - Cp_old)/Cp_new has to be smaller than this
            'GRID_NX' : 20,     # num of grid points (= num of elements + 1)
            'GRID_H' : 5000.0,  # height of the model, meters
            'BND_UPPER_PRESSURE' : 20000.0,           # pressure at the upper boundary, bar
            'BND_UPPER_TEMPERATURE' : 273.0 + 1000.0, # temperature at the upper boundary, Kelvin
            'BND_LOWER_TEMPERATURE' : 273.0 + 1125.0, # temperature at the lower boundary, Kelvin
        }

        self.components = self.perplex.components
        self.ncomponents = len(self.components)

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

        # latest perplex results
        self.perplexResult = None

        # physical constants
        self.accel_g = 9.81

        # flag to indicate whether perplex data calculated
        # for the grid points is up to date
        self.perplexOK = False

        # (initial) field specification
        # "_e" refers to elements (grid points in the middle of elements);
        # without '_e': the main grid points, defining (e.g.) upper and lower boundaries
        # of the model and placed also in between lithological units.
        # Num of main grid points (self.nx) = num of element grid points (self.ne) + 1
        # Most of the data is stored in the elements.
        self.T_e = self.getIniField("temp", self.xs_e)            # K
        self.C_e = self.getIniField("compo", self.xs_e)           # wt% (0..1)
        self.k_e = self.getIniField("conductivity", self.xs_e)    #
        self.Cp_e = self.getIniField("heatcapacity", self.xs_e)   # J/kgK
        self.rho_e = self.getIniField("density", self.xs_e)       # kg/m3
        self.pres_e = self.getIniField("pressure", self.xs_e)     # bar
        self.pres = self.getIniField("pressure", self.xs)
        self.mass_e = self.getIniField("mass", self.xs_e)         # kg
        #self.T = self.getIniField("temp", self.xs)                # K # seems not to be needed....

        # generate initial info by perplex
        self.updatePerplex()
        self.updateCpRho()

        # generate the mass information
        for i in range(self.ne):
            self.mass_e = self.rho_e * (self.xs[i+1]-self.xs[i])

    def updatePressure(self):
        for i in range(self.ne):
            self.pres[i] = self.bnd_upper_pressure + sum(self.rho_e[0:i] * self.accel_g * (self.xs[i+1]-self.xs[i])) * 1e-5
            self.pres_e[i] = self.pres[i]   # simplification: the pressure of the element is the pressure
                                            # at the upper surface of the element

    def updatePerplex(self):
        #self.updatePressure()
        self.perplexResult = []
        for i in range(self.ne):
            self.pres[i] = self.bnd_upper_pressure + sum(self.rho_e[0:i] * self.accel_g * (self.xs[i+1]-self.xs[i])) * 1e-5
            self.pres_e[i] = self.pres[i]
            self.perplexResult.append(self.perplex.phaseq(self.pres_e[i], self.T_e[i], self.C_e[i]))
        self.perplexOK = True

    def updateCpRho(self):
        if not self.perplexOK:
            self.updatePerplex()
        for i in range(self.ne):
            self.rho_e[i] = self.perplexResult[i]['SYSPROP'][self.perplex.syspropnum['rho']]
            self.Cp_e[i] = G2KG * self.perplexResult[i]['SYSPROP'][self.perplex.syspropnum['Cp']] / self.perplexResult[i]['SYSPROP'][self.perplex.syspropnum['N']]

    def maxdt(self):
        # return maximum time step for diffusion,
        # assume that information for k, rho and Cp is up-to-date
        dx = self.xs[1:self.nx] - self.xs[0:(self.nx-1)]
        diff = self.k_e / (self.rho_e * self.Cp_e)
        maxdt = min(0.5 * dx * dx / diff)
        return maxdt

    def diffuseT(self, dt = 0.0):
        self.timestep = self.timestep + 1
        print " * Diffusion, time step", self.timestep, ", time = ", (self.time / SECINYR), " yrs"
        self.updateCpRho()
        self.updatePressure()

        doNewIteration = True
        niterations = 0

        self.T_e_old = np.copy(self.T_e)
        T_e = self.T_e

        if dt == 0:
            # estimate largest time step, use a fraction of that
            dt = self.config['DT_MULTIPLIER'] * self.maxdt()
            self.last_dt = dt

        firstIteration = True

        while doNewIteration:
            # calculate conductivity for main grid points by averaging from elements
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
            DTinternal = d2T * dt / (self.rho_e[1:(self.ne-1)] * self.Cp_e[1:(self.ne-1)])

            # ... uppermost grid point
            dT1 = k_g[0]     * (T_e[1]     - T_e[0]) / (self.xs_e[1]     - self.xs_e[0])
            dT2 = self.k_e[0] * (T_e[0] - self.bnd_upper_temperature) / (self.xs_e[0]-self.xs[0])
            d2T = (dT1 - dT2) / (0.5 * (self.xs_e[1] - (-self.xs_e[0] + self.xs[0])))
            DTupper = d2T * dt / (self.rho_e[0] * self.Cp_e[0])

            # ... lowermost grid point
            dT1 = self.k_e[self.ne-1] * (self.bnd_lower_temperature - T_e[self.ne-1]) / (self.xs[self.nx-1]-self.xs_e[self.ne-1])
            dT2 = k_g[n_kg-1] * (T_e[self.ne-1] - T_e[self.ne-2]) / (self.xs_e[self.ne-1] - self.xs_e[self.ne-2])
            d2T = (dT1 - dT2) / (0.5 * (self.xs[self.nx-1] + (self.xs[self.nx-1]-self.xs_e[self.ne-1]) - self.xs_e[self.ne-2]))
            DTlower = d2T * dt / (self.rho_e[self.ne-1] * self.Cp_e[self.ne-1])

            self.T_e[1:(self.ne-1)] = self.T_e[1:(self.ne-1)] + DTinternal
            self.T_e[0] = self.T_e[0] + DTupper
            self.T_e[self.ne-1] = self.T_e[self.ne-1] + DTlower

            self.perplexOK = False  # temperature has been changed, perplex calculations out of date

            self.Cp_e_old = np.copy(self.Cp_e)
            self.rho_e_old = np.copy(self.rho_e)
            self.updateCpRho()
            self.updatePressure()


            #print self.Cp_e-self.Cp_e_old
            #print self.T_e-self.T_e_old

            if (max(abs((self.Cp_e-self.Cp_e_old)/self.Cp_e)) > self.config['CP_ITER_THRESHOLD']) or firstIteration:
                doNewIteration = True
                firstIteration = False
                self.T_e[0:self.ne] = self.T_e_old[0:self.ne]
                niterations = niterations + 1

                # relaxation
                self.Cp_e = 0.5 * self.Cp_e_old + 0.5 * self.Cp_e
                self.rho_e = 0.5 * self.rho_e_old + 0.5 * self.rho_e

                newdt = self.maxdt()
                if dt > self.config['DT_MULTIPLIER_MAX'] * newdt:
                    # Adjusting Cp/rho has lead to situation where timestep to be used is too close
                    # to the maximum possible timestep. Adjust timestep and restart iteration.
                    firstIteration = True
                    print "\tAdjusting dt from", self.last_dt/SECINYR, "to", newdt
                    self.last_dt = newdt
                    dt = self.last_dt
            else:
                doNewIteration = False

        self.time = self.time + self.last_dt
        print "\tIterations needed: ", niterations
        print "\tdt used: ", self.last_dt/SECINYR, "yrs"


    def getIniField(self, field, xs, fieldType = 0):
        if field == "temp":
            if fieldType == 0:
                retField = xs * 0.0 + 1000.0 + 273.0
            else:
                raise Exception("Invalid fieldType")
        elif field == "compo":
            if fieldType == 0:
                retField = []
                for i in range(len(xs)):
                    retField.append([[]] * self.ncomponents)
                    retField[i][0] = 0.4448
                    retField[i][1] = 0.0359
                    retField[i][2] = 0.0810
                    retField[i][3] = 0.3922
                    retField[i][4] = 0.0344
                    retField[i][5] = 0.0030
            else:
                raise Exception("Invalid fieldType")
        elif field == "conductivity":
            if fieldType == 0:
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





### MAIN PROGRAM:


cr = Crust()

while cr.time < 1e5 * SECINYR:
    cr.diffuseT()

print cr.T_e
print cr.Cp_e
print cr.rho_e









