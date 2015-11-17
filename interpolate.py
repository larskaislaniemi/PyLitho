import pylitho_exceptions as pyle 
import numpy as np

def interpolate(locarr, valarr, locarrnew, notFoundVal=float('NaN')):
    # find values from valarr at locations locarrnew

    valarrnew = []

    for loc in locarrnew:
        if loc < locarr[0] or loc > locarr[-1]:
            valarrnew.append(notFoundVal)
            continue

        try:
            (ilt, igt, which) = findNearest(locarr, loc)
        except pyle.PyLerr_OutsideRange as e:
            valarrnew.append(notFoundVal)
            continue
        except Exception as e:
            raise e

        if ilt == igt:
            newval = valarr[ilt]
        elif which >= 0:
            newval = valarr[ilt] + (valarr[igt]-valarr[ilt])*(loc-locarr[ilt])/(locarr[igt]-locarr[ilt])
        else:
            raise pyle.PyLerr_Undefined ("interpolate(): I don't understand findNearest()")

        valarrnew.append(newval)

    return valarrnew


def findNearest(locarr, loc, notFoundVal=None):
    # returns nearest index from locarr (minimize abs(loc-locarr[?]))
    # return values: idxleft, idxright, idxclosest (idxclosesti is either 0 or 1)
    #
    # About NaNs: Currently can handle a situation where there are NaNs at
    # either end of the array. Intermittent NaNs cause unknonwn behaviour.

    if loc < locarr[0] or loc > locarr[-1]:
        if notFoundVal is None:
            raise pyle.PyLerr_OutsideRange (loc)
        else:
            return notFoundVal
    gtearr = [loc >= val for val in locarr]
    ltearr = [loc <= val for val in locarr]
    andarr = [gtearr[i] and ltearr[i] for i in range(len(gtearr))]
    xorarr = [gtearr[i]  ^  ltearr[i] for i in range(len(gtearr))]
    if np.nansum(andarr) == 1:  # found exact location
        idx = andarr.index(True)
        return (idx, idx, 0)  #  0 could as well be 1
    elif np.nansum(andarr) == 0 and np.nansum(xorarr) == len(gtearr)-np.isnan(locarr).sum():
        ilt = ltearr.index(True)-1
        igt = [True if locarr[i]!=locarr[i] else gtearr[i] for i in range(len(gtearr))].index(False)
        # igt = gtearr.index(False) # <-- this is the idea, but this can't handle Falses in between, caused by nans
        if ilt != igt-1:
            raise pyle.PyLerr_Undefined ("ilt != igt-1")
        if abs(locarr[ilt]-loc) <= abs(locarr[igt]-loc):
            return (ilt, igt, 0)
        else:
            return (ilt, igt, 1)
    else:
        raise pyle.PyLerr_Undefined ("Strange locarr, unordered?")
