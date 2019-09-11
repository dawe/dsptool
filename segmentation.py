# !/usr/bin/env python3.7
import numpy as np
import pyBigWig as bw
from scipy import ndimage
from UliEngineering.SignalProcessing.Utils import zero_crossings
import argparse, sys, pybedtools, gc
from pybedtools import BedTool
from skimage import data, filters
gc.enable()
global save
save=[]
# -----------function Definitions-----------------------------------------------------------------------
def Sobel_filters(data):
    Kx = np.array([-1, 0, 1], np.float32)
    Ix = ndimage.filters.convolve(data, Kx)
    return(Ix)
#-----------------------------------------------
def Canny(data,M):
    low = 0.005*M
    high = 0.1*M
    hyst = filters.apply_hysteresis_threshold(data, low, high)
    return(hyst)
#-----------------------------------------------
def Maxima(data,value):
    x, y, b = [],[],{}
    CrossingCoordinates = zero_crossings(data)
    xp, yp, j = 0, 0, 0
    for i in CrossingCoordinates:
        if data[i] < j:
            x.append(i)
            if abs((value[xp]-value[i])/(value[xp]+value[i]+1e-9)) < 0.1 and abs((value[yp]-value[i])/(value[yp]+value[i]+1e-9)) < 0.1 :
                b[xp] = i
            else:
                b[i] = i
                xp = i
        else:
            y.append(i)
            yp=i
        j=data[i]
    return(x, y, b)
#----------------------------------------------
def Boundaries(value, Max, chr, start, end, step, Trues):
    end=int((end-start)/step)
    dxU=np.gradient(value)
    dxxU = np.gradient(dxU)
    s = np.sign(dxxU)
    signchangeco = np.where(np.diff(np.sign(dxxU)))[0]
    matrix = [[0]*5]
    for i in Max:
        if Trues[i]:
            #print('i=',i,'Max[i]',Max[i])
            sco = [x for x in signchangeco if x < i]
            sco = (sco[-1] if sco else 0)*step+start
            eco = [x for x in signchangeco if x > Max[i]]
            eco = (eco[0] if eco else end)*step+start
            matrix = np.vstack([matrix, (chr,sco,eco,i,Max[i])])
    matrix = np.delete(matrix, (0), axis=0)
    return(matrix)
#----------------------------------------------
def Write2BED(matrix,st,start):

    for i in matrix:
        a=i[0]
        b=int(i[1])
        c=(i[2])
        s=("%s\t%s\t%s" % (a,b,c))
        save.append((s))
    pybedtools.BedTool(save).saveas(s_boundary)
    #pybedtools.BedTool(save()).saveas('counted.bed')
#----------------------------------------------
def segmentation(chr,start,end,step):
    print('  - Peak finding chromosome',chr,'from',start,'to',end)
    sg_value = sg_input.values(chr, start, end, numpy=True)[::step]
    Med=(np.median(sg_value))
    ZeroCrossList = Sobel_filters(sg_value)     #using ndimage
    Trues = Canny(sg_value,Med)
    ListofMaxima, ListofMinima, CoMax = Maxima(ZeroCrossList,sg_value)
    Edges = Boundaries(sg_value, CoMax, chr, start, end, step, Trues)
    Write2BED(Edges,step,start)
    j,k=0,0
    if ListofMaxima:
        for i in ListofMaxima:
            if Trues[i]:
                if sg_step > i-j:
                    pass
                else:
                    k = i*step+start
                    intensity =float(sg_value[i])
                    j = i
                    sg_output.addEntries(chr, k, values=[intensity], span=step, step=step)
            else:
                pass
#----------------------------------------------
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=("""
##############################################################################
################ dsptool is a tool for ChIP-seq data analysis ################
##############################################################################
This module could find the peaks from the denoised BigWig files. You could define the name of the output file, step size and the region of intreset. Program will scan entire input file if no region defined by user.
"""))
parser.add_argument("-i", "--input", help="Input file for the denoising module.", required=True)
parser.add_argument("-p", "--peak", help="Define the name and the location of the output BigWig file. If no output file defined, the result will be stored in Peak.bw file in the running folder.", default="Peaks.bw")
parser.add_argument("-b", "--boundary", help="Define the name and the location of the output BigWig file. If no output file defined, the result will be stored in Peak.bw file in the running folder.", default="boundary.bed")
parser.add_argument("-S", "--step", type=int, help="The step size that denoising process applied by, the default step size is 50bps.", default=50)
parser.add_argument("-r", "--region", help="A single section of input file could be defined as chromosomename:startindex:endindex")
parser.add_argument("-seg", "--segmentation", action='store_true', help="Run segmentation")
parser.add_argument("-l", "--interval", help="A list of regions could be defined as a BED file format")
parser.add_argument("-H", "--highresolution", action='store_true', help="Force the program to define step and span equal to 1 to have the highest resolution as possible (One base resolution)")

args = parser.parse_args()
s_input = args.input
s_output = args.peak
s_boundary = args.boundary
sg_region = args.region
sg_interval = args.interval
sg_span = sg_step = int(args.step)
if args.highresolution:
    sg_step, sg_span = 1, 1
sg_input = bw.open(s_input)
sg_output = bw.open(s_output, "w")
sg_output.addHeader(list(sg_input.chroms().items()))
chrlist=list(sg_input.chroms())
if sg_region:
    try:
        # Characters before column-sign should be a chromosome name
        _splited = sg_region.strip().split(':')
        sg_name = _splited[0]
        chrlist.index(sg_name)
        # The phrase after column-sign consists of two coordinates that separated by a dash-sign
        _splited = _splited[1].strip().split('-')
        sg_start = int(_splited[0])
        sg_end = int(_splited[1])
    except:
        sys.stderr.write('The region of interest must follow the standard pattern ChromosomeName:StartBaseIndex-EndBaseIndex\n')
        exit()
    if (int(sg_input.chroms(sg_name)) < sg_end) :
        sys.stderr.write('Interval definition is incorrect (Larger than the length of the chromosome).\n')
        question = input('Do you want to continue? (Y/N)').lower()
        if question == 'y':
            sg_end = int(sg_input.chroms(sg_name))
        else:
            print('Please correct the end coordinate and try again.')
            exit()


    segmentation(sg_name, sg_start, sg_end, sg_step)
elif sg_interval:
    import  tempfile, pybedtools, pathlib
    bed_file = pathlib.Path(sg_interval)
    # define a temperory file and put the list of the chromosome name inside it, because the Bedtools does not accept the list as a array
    tmp = tempfile.NamedTemporaryFile()
    with open(tmp.name, 'w') as _tmp:
        for line in chrlist:
            _tmp.write(str(line+"\n"))
    # Call temprory file
    with open(tmp.name) as _tmp:
        # Check if the Bed file present, continue
        if bed_file.is_file():
            try:
                # Interval data derived from Bed file sorted and merged if their steps overlapped
                sites = list(pybedtools.BedTool(sg_interval).sort(g = _tmp.name).merge(d=sg_step-1))
            except:
                # for bedtools sort Version: v2.26.0
                sites = list(pybedtools.BedTool(sg_interval).sort(faidx = _tmp.name).merge(d=sg_step-1))
            # For each interval the value reads and stored in a varible
            for line in sites:
                _temp=str(line)
                L = _temp.strip().split()
                print(L)
                segmentation(L[0], int(L[1]), int(L[2]), sg_step)
        else:
            sys.stderr.write("File \"" + sg_interval + "\" is not exist. Check the BED file and it\'s path. \n")
            exit()
else:
    sg_start=1
    for line in chrlist:
        sg_end=sg_input.chroms(line)
        segmentation(line, sg_start, sg_end, sg_step)
sg_input.close()
sg_output.close()
