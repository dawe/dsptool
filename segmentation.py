# !/usr/bin/env python3.7
import numpy as np
import pyBigWig as bw
from scipy import ndimage
from UliEngineering.SignalProcessing.Utils import zero_crossings
# -----------function Definitions-----------------------------------------------------------------------
def Sobel_filters(data):
    Kx = np.array([-1, 0, 1], np.float32)
    Ix = ndimage.filters.convolve(data, Kx)
    return(Ix)
#-----------------------------------------------
def Maxima(data):
    x, y, z = [],[],[]
    CrossingCoordinates=zero_crossings(data)
    for i in CrossingCoordinates:
        if data[i] < data[i+1]:
            x.append(i)
            if data[i] == 0:
                z.append(i)
        else:
            y.append(i)
    return(x, y, z)
#----------------------------------------------
def Merge(max,flat):
    dic={}
    for i in max:
        if i in flat:
            dic[i-1] = i
        else:
            dic[i] = i+1
    return(dic)

#---------------------------------------------
def Callregions (chr,coordinates,sp,st,header,value,shift):
    # check if the output file cannot be created because of the output folder does not exist or permissions an error message with clean exit will occur
    try:
        output = bw.open("Peaks.bw", "w")
        output.addHeader(header)
    except IOError:
        sys.stderr.write("\n No permission to write. \n")
        exit()
    j,k=0,0
    for i in coordinates:
        if st >= i-j:
            pass
        else:
            k = i*st+shift
            intensity = value[i]
            output.addEntries(chr, k, values=[intensity], span=sp, step=st)
            j = i
    output.close
