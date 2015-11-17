import numpy as np


def initTemp(STATUS, T, x):
    # return values for initial temperature field

    if not isinstance(x, np.ndarray):
        raise Exception("x must be numpy.ndarray")

    if not isinstance(T, np.ndarray):
        raise Exception("T must be numpy.ndarray")

    if len(T) != len(x):
        raise Exception("Incompatible x and T")

    if STATUS['CONFIG']['TINI_TYPE'] == -1:
        pass # T read at restart and no modification needed
    elif STATUS['CONFIG']['TINI_TYPE'] == 0:
        T[:] = 1.0
    elif STATUS['CONFIG']['TINI_TYPE'] == 1:
        T[:] = 1519.0 * x[:]/STATUS['L'][1] + 10.0
    elif STATUS['CONFIG']['TINI_TYPE'] == 10:
        # extra 30km overthrust
        idx = x < 30e3
        T[idx] = 750.0 * x[idx]/30e3

        idx = (x >= 30e3) & (x < 30e3+STATUS['Moho_Depth'])
        T[idx] = 750.0 * (x[idx]-30e3)/STATUS['Moho_Depth']

        idx = x >= 30e3+STATUS['Moho_Depth']
        T[idx] = 750.0 + (1250.0-750.0) * (x[idx]-(30e3+STATUS['Moho_Depth']))/(STATUS['L'][1]-(30e3+STATUS['Moho_Depth']))
    elif STATUS['CONFIG']['TINI_TYPE'] == 101:
        # T read at restart but we want to modify it
        pass
    else:
        raise Exception("Invalid TINI_TYPE " + str(STATUS['CONFIG']['TINI_TYPE']))


def prod(x):
    ret = 1
    for i in x:
        ret = ret * i
    return ret


def kT(reltype, k0, params, T):

    if reltype == 0:
        retk = np.copy(k0)

    elif reltype == 1:
        retk = k0/(1+params[0]*T)

    elif reltype == 2:
        MAXK = 1e2
        MINK = 0
        retk = k0/(1+params[0]*T) + params[1]*T**3
        retk[retk < MINK] = MINK
        retk[retk > MAXK] = MAXK

    else:
        raise Exception("Invalid reltype")

    return retk


def cT(reltype, c0, params, T):
    if reltype == 0:
        retc = np.copy(c0)
        pass
    elif reltype == 1:
        # From Waples and Waples (2004)
        Tt = 20 # measuring temp
        cpnt1 = params[0] * Tt**3 + params[1] * Tt**2 + params[2]*Tt + params[3]
        D = params[4] # Debye temp of a silicate

        Tlimited = np.copy(T)
        Tlimited[Tlimited > D] = D

        cpnt2 = params[0] * Tlimited**3 + params[1] * Tlimited**2 + params[2]*Tlimited + params[3]

        retc = c0 * cpnt2 / cpnt1
    else:
        raise Exception("Invalid reltype")

    return retc


def maxdt(k, rho, cp, xs):
    NX = np.size(k)
    dx = xs[1:NX] - xs[0:(NX-1)]
    diff = k / (rho*cp)
    Khalf = 0.5 * (diff[1:NX] + diff[0:(NX-1)])
    return min(0.5 * dx * dx / Khalf)


def diffstep(STATUS, T, new_T, xs, k, dt, Tsurf, q0, rho, cp, H):
    NX = np.size(T)

    if NX != np.size(new_T):
        raise Exception("dimension mismatch")

    if NX != np.size(k):
        raise Exception("dimension mismatch")

    if NX != np.size(xs):
        raise Exception("dimension mismatch")

    if NX != np.size(rho):
        raise Exception("dimension mismatch")

    if NX != np.size(cp):
        raise Exception("dimension mismatch")
        
    khalf = 0.5 * (k[1:NX] + k[0:(NX-1)])
    dx = xs[1:NX] - xs[0:(NX-1)]
    dxhalf = 0.5 * (dx[1:(NX-1)] + dx[0:(NX-2)])

    # upper bnd
    new_T[0] = Tsurf
    
    # inner points
    new_T[1:(NX-1)] = khalf[1:(NX-1)] * (T[2:NX]-T[1:(NX-1)]) / dx[1:(NX-1)]   -   khalf[0:(NX-2)] * (T[1:(NX-1)] - T[0:(NX-2)]) / dx[0:(NX-2)]
    new_T[1:(NX-1)] = 2.0 * new_T[1:(NX-1)] / dxhalf[0:(NX-2)]
    new_T[1:(NX-1)] = new_T[1:(NX-1)] * dt / (rho[1:(NX-1)] * cp[1:(NX-1)])   +   T[1:(NX-1)]
    new_T[1:(NX-1)] = new_T[1:(NX-1)]   +   H[1:(NX-1)] * dt / (rho[1:(NX-1)] * cp[1:(NX-1)])

    # W m-3 => J m-3 => J m-3 J-1 kg K => K kg m-3 => K

    #for ix in range(1,NX-1):
    #    new_T[ix] = khalf[ix] * (T[ix+1] - T[ix]) / dx[ix] - khalf[ix-1] * (T[ix] - T[ix-1]) / dx[ix-1]
    #    new_T[ix] = 2.0 * new_T[ix] / dxhalf[ix-1]
    #    new_T[ix] = new_T[ix] * dt / (rho[ix] * cp[ix]) + T[ix]
    #    new_T[ix] = new_T[ix] + H[ix] * dt / (rho[ix] * cp[ix])
    
    # lower bnd
    if STATUS['CONFIG']['BND_BOT_TYPE'] == 0:
        new_T[NX-1] = STATUS['CONFIG']['BND_BOT_TEMP']
    elif STATUS['CONFIG']['BND_BOT_TYPE'] == 1:
        new_T[NX-1] = q0 * dx[NX-2] / khalf[NX-2] + T[NX-2]
    elif STATUS['CONFIG']['BND_BOT_TYPE'] == 9:
        TL = 1250.0
        try:
            botloc = min(np.where(new_T >= TL)[0])
        except ValueError:
            botloc = NX-1

        botq = (rho[botloc]**2 * cp[botloc] * 1.0**4 * 60.0**4 * 9.81 * 3.5e-5 * (xs[botloc]-xs[0])**2) \
               / (1e19*new_T[botloc]**2)

        for ix in range(botloc, NX):
            new_T[ix] = botq * dx[ix-1] / khalf[ix-1] + T[ix-1]

        #print botq, botloc
    else:
        raise Exception("Invalid BND_BOT_TYPE")


def getErosionSpeed(STATUS):
    if STATUS['CONFIG']['EROSION_SPEED_TYPE'] == 0:
        STATUS['Erosion_Speed'] = STATUS['CONFIG']['EROSION_SPEED_M_MA'] / (1e6*STATUS['SECINYR'])

    elif STATUS['CONFIG']['EROSION_SPEED_TYPE'] == 1:
        if (STATUS['curTime'] - STATUS['ModelStartTime']) > STATUS['MaxTimeToErode']:
            STATUS['Erosion_Speed'] = 0.0
        else:
            STATUS['Erosion_Speed'] = STATUS['CONFIG']['EROSION_SPEED_M_MA'] / (1e6*STATUS['SECINYR'])
    elif STATUS['CONFIG']['EROSION_SPEED_TYPE'] == 10:
        # only erode the original overthrust sheet
	if STATUS['CONFIG']['RESTART_POST_MOD'] != 2:
	    raise Exception("EROSION_SPEED_TYPE == 10 requires RESTART_POST_MOD == 2")

        if (STATUS['curTime'] - STATUS['ModelStartTime']) > STATUS['CONFIG']['RESTART_POST_MOD_PARAMS'][0] / STATUS['CONFIG']['EROSION_SPEED_M_MA']:
	    STATUS['Erosion_Speed'] = 0.0
	else:
	    STATUS['Erosion_Speed'] = STATUS['CONFIG']['EROSION_SPEED_M_MA'] / (1e6*STATUS['SECINYR'])

    else:
        raise Exception("Invalid EROSION_SPEED_TYPE")


def remesh(STATUS, curxs, newExt, arrays, extrapolation=0):
    # modify curxs to newExt and interpolate values within the arrays
    # newExt is either the extent (min,max)
    # or an array of new xs

    NX = 0

    if len(newExt) == 2:
        if newExt[0] == curxs[0] and newExt[1] == curxs[-1]:
            # we're done here
            return
        NX = STATUS['NX']
        if newExt[1] == curxs[-1]:
            # bottom location hasn't changed

            if (newExt[0] > curxs[0]) and (newExt[0] < curxs[1]):
                # only the uppermost grid point "has" moved
                #print newExt[0], curxs[0], curxs[1], curxs[2]
                if (curxs[1] - newExt[0]) < 0.5 * (curxs[2] - curxs[1]):
                    # grid spacing getting too tight, remove uppermost grid point
                    newxs = np.zeros(NX-1)
                    newxs[:] = curxs[1:]
                    newxs[0] = newExt[0]
                    NX = NX-1
                    STATUS['NX'] = NX
                    #print NX
                else:
                    newxs = np.copy(curxs)
                    newxs[0] = newExt[0]
            else:
                # TODO:
                # make this better so that a jump over one grid point
                # is also properly handled without extra interpolations
                newxs = np.linspace(newExt[0], newExt[1], num=NX)
        else:
            newxs = np.linspace(newExt[0], newExt[1], num=NX)

    else:
        # newExt is an array of xs
        if len(np.where(newExt == curxs)) == len(curxs):
            # we're done here
            return
        NX = len(newExt)
        STATUS['NX'] = NX
        newxs = np.copy(newExt)

    addGridPoints = NX - len(curxs)
    if addGridPoints != 0:
        #print "remesh(): Adding " + str(addGridPoints) + " grid points"
        pass

    # apply the new grid to the value arrays
    for i in range(len(arrays)):
        newarr = np.zeros(NX)
        interpolate(curxs, arrays[i], newxs, newarr, extrapolation)
        if addGridPoints != 0:
            arrays[i].resize(NX, refcheck=False)
        arrays[i][0:NX] = newarr[0:NX]

    if addGridPoints != 0:
        curxs.resize(NX, refcheck=False)

    curxs[0:NX] = newxs[0:NX]

    return None

def posMin(arr):
    return np.where(arr == min(arr[arr > 0]))

def findPoint(xs, val):
    # returns the first point in xs that is larger than loc
    # raises ValueError is not found
    # returns the last occurrence
    dist = xs - val
    i = posMin(dist)

    if len(i[0]) != 1:
        raise Exception("Error in findPoint()")
    return i[0][0]


def interpolate(xs1, arr1, xs2, arr2, extrapolation):
    # interpolate from arr1 to arr2

    for i in range(len(xs2)):
        if xs2[i] < xs1[0]:
            if extrapolation <= 0:
                arr2[i] = -extrapolation
            elif extrapolation == 1:
                arr2[i] = np.NaN
            elif extrapolation == 2:
                relloc = xs2[i] - xs2[0]
                j = findPoint(xs1-xs1[0], relloc)
                arr2[i] = arr1[j-1] + (relloc+xs1[0] - xs1[j-1]) * (arr1[j] - arr1[j-1]) / (xs1[j] - xs1[j-1])
            else:
                raise Exception("Don't know how to extrapolate")
        elif xs2[i] > xs1[-1]:
            if extrapolation == 1:
                arr2[i] = np.NaN
            elif extrapolation == 2:
                raise Exception("Extrapolation: Copying at the lower bnd not supported yet")
            else:
                raise Exception("Don't know how to extrapolate")
        else:
            idx = xs2[i] == xs1
            whereIdx = np.where(idx)
            nSameNodes = len(whereIdx[0])
            if nSameNodes == 0:
                # interpolate, linear
                idx2 = np.where(xs1 > xs2[i])[0][0]
                idx1 = idx2 - 1
                arr2[i] = arr1[idx1] + (xs2[i] - xs1[idx1]) * (arr1[idx2] - arr1[idx1]) / (xs1[idx2] - xs1[idx1])
            elif nSameNodes == 1:
                # copy value directly
                arr2[i] = arr1[whereIdx]
            elif nSameNodes > 1:
                raise Exception("Invalid grid xs1: multiple same values")
