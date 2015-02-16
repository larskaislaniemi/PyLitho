import numpy as np


def initTemp(STATUS, T, x):
    # return values for initial temperature field

    if not isinstance(x, np.ndarray):
        raise Exception("x must be numpy.ndarray")

    if not isinstance(T, np.ndarray):
        raise Exception("T must be numpy.ndarray")

    if len(T) != len(x):
        raise Exception("Incompatible x and T")

    if STATUS['CONFIG']['TINI_TYPE'] == 0:
        T[:] = 0.0
    elif STATUS['CONFIG']['TINI_TYPE'] == 1:
        T[:] = 1519.0 * x[:]/STATUS['L'][1]
    elif STATUS['CONFIG']['TINI_TYPE'] == 10:
        # extra 30km overthrust
        idx = x < 30e3
        T[idx] = 750.0 * x[idx]/30e3

        idx = (x >= 30e3) & (x < 30e3+STATUS['Moho_Depth'])
        T[idx] = 750.0 * (x[idx]-30e3)/STATUS['Moho_Depth']

        idx = x >= 30e3+STATUS['Moho_Depth']
        T[idx] = 750.0 + (1250.0-750.0) * (x[idx]-(30e3+STATUS['Moho_Depth']))/(STATUS['L'][1]-(30e3+STATUS['Moho_Depth']))
    else:
        raise Exception("Invalid TINI_TYPE")


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


def diffstep(T, new_T, xs, k, dt, Tsurf, q0, rho, cp, H):
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
    new_T[NX-1] = q0 * dx[NX-2] / khalf[NX-2] + T[NX-2]
    #new_T[NX-1] = 1350


def remesh(curxs, newExt, arrays):
    # NB : currently ignores new extent lower bnd

    if newExt[0] == curxs[0] and newExt[1] == curxs[-1]:
        # we're done here
        return

    NX = len(curxs)
    #curL = (min(curxs), max(curxs))

    #if curL[1] != newExt[1]:
    #    raise Exception("Only surface level adjustment supported at the moment")

    #newxs = np.copy(curxs)
    #
    #jumpedNodes = sum(curxs < newExt[0]) - 1
    #
    #if jumpedNodes == 0:
    #    # surface level still within the first element
    #    lowerAdjIndex = 2
    #    newxs[0:(lowerAdjIndex+1)] = np.linspace(newExt[0], curxs[lowerAdjIndex], num=lowerAdjIndex+1)
    #else:
    #    lowerAdjIndex = jumpedNodes*2   # num of grid points to adjust is double of those actually jumped over
    #    newxs[0:(lowerAdjIndex+1)] = np.linspace(newExt[0], curxs[lowerAdjIndex], num=lowerAdjIndex+1)

    newxs = np.linspace(newExt[0], curxs[-1], num=NX)

    # apply the new grid to the value arrays
    for i in range(len(arrays)):
        newarr = np.copy(arrays[i])
        interpolate(curxs, arrays[i], newxs, newarr)
        arrays[i][0:NX] = newarr[0:NX]

    curxs[0:NX] = newxs[0:NX]

    return None


def interpolate(xs1, arr1, xs2, arr2):
    # interpolate from arr1 to arr2

    if min(xs1) > min(xs2) or max(xs1) < max(xs2):
        raise Exception("This is for interpolation, not extrapolation!")

    for i in range(len(xs2)):
        idx = xs2[i] == xs1
        nSameNodes = len(np.where(idx)[0])
        #nSameNodes = sum(idx)
        if nSameNodes == 0:
            # interpolate, linear
            idx2 = np.where(xs1 > xs2[i])[0][0]
            idx1 = idx2 - 1
            arr2[i] = arr1[idx1] + (xs2[i] - xs1[idx1]) * (arr1[idx2] - arr1[idx1]) / (xs1[idx2] - xs1[idx1])
        elif nSameNodes == 1:
            # copy value directly
            arr2[i] = arr1[np.where(idx)]
        elif nSameNodes > 1:
            raise Exception("Invalid grid xs1: multiple same values")
