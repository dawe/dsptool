sg_step=sg_span=50
# !/usr/bin/env python3.7
import gc
import numpy as np
import pyBigWig as bw
from scipy import ndimage
from UliEngineering.SignalProcessing.Utils import zero_crossings
import argparse, sys, time
gc.enable()

# -----------function Definitions-----------------------------------------------------------------------
def Sobel_filters(data):
    Kx = np.array([-1, 0, 1], np.float32)
    Ix = ndimage.filters.convolve(data, Kx)
    return(Ix)
#-----------------------------------------------
def Maxima(data):
    x, y, z = [],[],[]
    CrossingCoordinates = zero_crossings(data)
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
#----------------------------------------------
def segmentation(chr,str,end,step):
    print(' - Peak finding chromosome',chr,'from',str,'to',end)
    sg_value = sg_input.values(chr, str, end, numpy=True)[::step]
    ZeroCrossList = Sobel_filters(sg_value)     #using ndimage
    ListofMaxima, ListofMinima, ListofPlateau = Maxima(ZeroCrossList)
    MergedList = Merge(ListofMaxima,ListofPlateau)
    j,k=0,0
    if MergedList:
        for i in MergedList:
            if sg_step > i-j:
                pass
            else:
                k = i*step
                intensity =float(sg_value[i])
                j = i
                sg_output.addEntries(chr, k, values=[intensity], span=step, step=step)
#----------------------------------------------
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=("""
##############################################################################
################ dsptool is a tool for ChIP-seq data analysis ################
##############################################################################
This module could find the peaks from the denoised BigWig files. You could define the name of the output file, step size and the region of intreset. Program will scan entire input file if no region defined by user.
"""))
parser.add_argument("-i", "--input", help="Input file for the denoising module.")
parser.add_argument("-o", "--output", help="Define the name and location of the BigWig file as input. e.g. -i path/file.bw ",required=True)
parser.add_argument("-p", "--peak", help="Define the name and the location of the output BigWig file. If no output file defined, the result will be stored in Peak.bw file in the running folder.", default="Peaks.bw")
parser.add_argument("-S", "--step", type=int, help="The step size that denoising process applied by, the default step size is 50bps.", default=50)
parser.add_argument("-r", "--region", help="A single section of input file could be defined as chromosomename:startindex:endindex")
parser.add_argument("-seg", "--segmentation", action='store_true', help="Run segmentation")
args = parser.parse_args()
input = args.output
output = args.peak
sg_region = args.region
sg_step = int(args.step)

sg_input = bw.open(input)
sg_output = bw.open(output, "w")
sg_output.addHeader(list(sg_input.chroms().items()))
chrlist=list(sg_input.chroms())
if sg_region:
    try:
        # Characters before column-sign should be a chromosome name
        _splited = sg_region.strip().split(':')
        sg_name = _splited[0]
        # The phrase after column-sign consists of two coordinates that separated by a dash-sign
        _splited = _splited[1].strip().split('-')
        sg_start = int(_splited[0])
        sg_end = int(_splited[1])
    except:
        sys.stderr.write('The region of interest must follow the standard pattern ChromosomeName:StartBaseIndex-EndBaseIndex\n')
        exit()
    if (int(sg_input.chroms(sg_name)) < sg_end) :
        sys.stderr.write('Interval definition is incorrect (Larger than the length of the chromosome).\n')
        exit()
    segmentation(sg_name, sg_start, sg_end, sg_step)
else:
    sg_start=1
    for line in chrlist:
        sg_end=sg_input.chroms(line)
        segmentation(line, sg_start, sg_end, sg_step)

sg_input.close()
sg_output.close()
