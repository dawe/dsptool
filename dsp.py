#!/usr/bin/env python

from __future__ import with_statement
__author__ = ('Shahram Saghaei (saghaei.shahram@hsr.it), '
              'Davide Cittaro (cittaro.davide@hsr.it)')
__version__ = '0.1.00'
__date__ = '12 Jan 2020'

import time, sys, os
dsp_time = time.time()
sys.stdout.write('\r Please wait ...')
os.system('clear')


# def parseArgs(*args, **kwargs):
#   print(f"parsing args, \n args: {args}\n kwargs: {kwargs}")

# def one():
#   print("process one")

# if __name__ == "__main__":
#   parseArgs(*sys.argv)
#   one()

# print("this is sourced")


# if __name__ == "__main__":


try:
    import pyBigWig, sys, argparse, re, os, tempfile, pybedtools, pathlib, gc ,skimage
    import scipy.signal
    from pathlib import Path
    from itertools import groupby
    from operator import itemgetter
    import numpy as np
    import sklearn.mixture
    import matplotlib.pyplot as plt
    from package import excluder as EX
    from package import segment as SG
except ImportError:
    raise ImportError("""\n#\n# Excluder error: The following packages are required, please install them.:

package                    |            build
---------------------------|-----------------
pyBigWig                   |           0.3.17
python-3.7.6               |       h359304d_2
------------------------------------------------------------

Also you can use dsp.yml file to create a conda enviroment and activate it, use following commands:\n -  conda env create -f dsp.yml \n and then \n -  conda activate dsptool\n""")
    sys.exit()

d_warning=False
gc.enable()
# -----------available filters list-----------------------------------------------------------------
# Define a class for the optional parameter "--list-filters"
class _ListAction(argparse.Action):
    def __init__(self,option_strings,dest=argparse.SUPPRESS,default=argparse.SUPPRESS,help=None):
        super(_ListAction, self).__init__(option_strings=option_strings,dest=dest,default=default,nargs=0,help=help)
    def __call__(self, parser, namespace, values, option_string=None):
        sys.stdout.write("\n Avalable filters for Scipy Signal: \n -  hann, hamming, blackman, blackmanharris, triang, parzen, cosine, nuttall, bohman, barthann, bartlett \n -  The Blackman filter is selected by defult.\n ")
        parser.exit()

# -----------input definitions----------------------------------------------------------------------
# Required and optional parameters definition usnig "argparse" package
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=("""
##############################################################################
################ dsptool is a tool for ChIP-seq data analysis ################
##############################################################################
Version:   v1.0
About:     developed in the HSR
Code:      https://github.com/dawe/dsptool
Mail:
Usage:     dsptool -i input -o output -r region [options]
The dsptool options for denoise include:
[ Main Parameters ]
    input file         Input files must be in BigWig file format.
    output file        A BigWig file that could be explored by IGV.
[ Interval Definitions ]
    region              Apply filter on a specific region defined by chromosome name and location.
    interval            Cover a list of regions using BED files.
    entire              Apply filter on the entire BigWig input file.
[ Optional Parameters ]
    filter              Scipy signal filter selection changes the pattern of applied signal (Blackman by default selected).
    size                Window size of the filter in the basepair unit.
    step                The distance between the start of one region and the end of the previous region.
    span                The fixed distance from the start position which defines a selected region.
[ General help ]
    help                Print this help menu.
    version             What version of bedtools are you using?.
    filter list         Print all the available scipy signal filter name that could be applied on the input signal.
"""))
parser.add_argument("-i", "--input", help="Input file for the denoising module.",required=True)
parser.add_argument("-d", "--denoised", help="Define the name and location of the BigWig file as input. e.g. -i path/file.bw ",required=True)
parser.add_argument("-S", "--step", type=int, help="The step size that denoising process applied by, the default step size is 50bps.", default=50)
parser.add_argument("-r", "--region", help="A single section of input file could be defined as chromosomename:startindex:endindex")
parser.add_argument("-seg", "--segmentation", action='store_true', help="Run segmentation exclusively")
parser.add_argument("-l", "--interval", help="A list of regions could be defined as a BED file format")
parser.add_argument("-H", "--highresolution", action='store_true', help="Force the program to define step and span equal to 1 to have the highest resolution as possible (One base resolution)")
parser.add_argument("-s", "--windowsize", type=int, help="Define the window size of applying filter", required=False, default="131072")
parser.add_argument("-V", "--version", action='version', version="%(prog)s (Version 1.0)")
parser.add_argument("-e", "--exclude", help="Define the name and location of the BigWig file as BlackList_file. e.g. -i data/File.bw ")
parser.add_argument("-temp", "--tempfolder", help="The address of the temporary folder.", default="temp/")
parser.add_argument("-f", "--filter", help="Insert the filter name of Scipy Signal filter list. Blackman filter selected by default", required=False, default="blackman")
parser.add_argument("-list", "--filterlist", action=_ListAction, help="Print all the possible signal filters")
parser.add_argument("-fold", "--foldtimes", type=int, help="Number of folds of low resolution window size in compare to high resolution denoising window size.", required=False, default="3")
parser.add_argument("-b", "--boundaryfile", help="Define the name of the BED file as output. e.g. -b file.bed", required=False, default="Boundaries.bed")
parser.add_argument("-save", "--saveinterval", help="Export interval after exclusion of zero-filled-region and defined blacklist as a BED file. e.g. -save Interval.bed")
parser.add_argument("-t", "--threshold", action='store_true', help="Turns on the Auto-Thresholding using the parameters of a Gaussian mixture distribution.")

# Use "args" instead of the subcommand "parse_args"
args = parser.parse_args()
# Store the inserted values for switches on the variables
source_data = args.input
denoised_data = args.denoised
d_filter = args.filter
userdefined_windowsize = int(args.windowsize)
Region_of_interest = args.region
Span_size = Step_size = int(args.step)
user_defined_interval = args.interval
BlackList_file = args.exclude
TEMPFolder = args.tempfolder
foldtimes = args.foldtimes
boundaryfile = args.boundaryfile
SaveInterval_file = args.saveinterval
GoToSegmentation = args.segmentation

if TEMPFolder[-1] != '/':
    TEMPFolder=str(TEMPFolder)+'/'
try:
    os.makedirs(TEMPFolder)
except:
    # print('Failed to make Dir and passed\n')
    pass

# -----------input/step validation----------------------
# If the high resolution option selected by user, it is useful for the short peaks that needs more resulations. it means all the bases of defined region will be called and processed.
if args.highresolution:
    Step_size, Span_size = 1, 1

# find log base 2 of user defined window size and sum it by the foldtime which represents the windows size of the low resolution smoothed signal
primewindowsize = 2 << (int(np.log2(userdefined_windowsize))+foldtimes-1)
Exwindowsize = 2 << (int(np.log2(userdefined_windowsize))+foldtimes)
# Window size of prime signal
PWS = primewindowsize >> int(np.log2(Step_size))
# ExWS = Exwindowsize >> int(np.log2(Step_size))                <<<<<< For add extremely low resolution
# Window size of the signal resized regaurd the input file
WSize = userdefined_windowsize >> int(np.log2(Step_size))
# Location of the input file stores on "my_file" variable
my_file = pathlib.Path(source_data)
# Check whether the path/file presence or not
if my_file.is_file():
    pass
else:
    sys.stderr.write('\n# Filter error: File \"'+ str(my_file) +'\" is not exist. Check the file name and relative path.\n')
    sys.exit()
# current_folder = (os.getcwd())                            <<<<<<< Remove
# print(current_folder)                                     <<<<<<< Remove
# -----------Excluder / Blacklist-----------------------

# check if the output file cannot be created because of the output folder does not exist or permissions, an error message with clean exit will occur
if GoToSegmentation:
    try:
        Open_file = pyBigWig.open(denoised_data)
        ChrLIST = Open_file.chroms()
        Header = list(Open_file.chroms().items())
    except IOError:
        sys.stderr.write("\n# Filter error:  The PATH/FILE `"+ denoised_data +"` not valid, please check the path and try again. \n - \n")
        sys.exit()
else:
    try:
        Open_file = pyBigWig.open(source_data)
        ChrLIST = Open_file.chroms()
        Header = list(Open_file.chroms().items())
    # This module could be applied on denoised signals before segmentation. Denoised file will load if exist
    except:
        sys.stderr.write('\n# Excluder error: Input file not defined.\n')
        sys.exit()
    try:
        OUT_file = pyBigWig.open(denoised_data, "w")
        OUT_file.addHeader(Header)
        L = denoised_data[:-3:]
        low_resolution_output = pyBigWig.open(L+'L.bw', "w")
        low_resolution_output.addHeader(Header)
    except IOError:
        sys.stderr.write("\n# Filter error:  The PATH/FILE `"+ denoised_data +"` has not write permission, please check the path and try again. \n - \n")
        sys.exit()

if primewindowsize > int(Header[0][1]):
    sys.stderr.write("\n# Filter error:  Window Size for low_resolution denoising is greater than chromosme size. \n - \n")
    sys.exit()

sys.stdout.write('\r Excluding ...')

if BlackList_file:
    try:
        # Import blacklist file, and merge the intervals if they have overlap
        BL = pybedtools.BedTool(BlackList_file).sort().merge()
    except:
        sys.stderr.write('\n# Excluder error: Blacklist file not found.\n')
        sys.exit()
else:
    BL =[]

if Region_of_interest:
    target='Region'
    zone=Region_of_interest
elif user_defined_interval:
    target='Interval'
    zone=user_defined_interval
else:
    target='Entrie_genome'
    zone=''
INPint = EX.intervalconverter(BL, target, ChrLIST, Step_size, zone) # list of intervals
INint = pybedtools.BedTool(INPint)
#---------------Zero function-------------
Z=[]
if len(INPint) == 0:
    sys.stderr.write("\r# Excluder error: The interval file is empty.\n\n")
    exit()
for line in INPint:
    sys.stdout.write('\r Excluding '+line.chrom+'           ')
    try:
        openvalue = Open_file.values(line[0], int(line[1]), int(line[2]), numpy=True)
    except:
        sys.stderr.write("\n# Excluder error: pyBigWig package is not installed completely, please install it using conda.\n To assure it is installed properly, import pyBigWig in python, if you get '0' for 'pyBigWig.np' command, please install it using conda as following:\nconda install -c bioconda pybigwig\n")
        exit()
    openvalue[np.isnan(openvalue)] = 0
        # INPcoord.append(line)                                 # <<<<<<< Remove
    #INPvalue.append(openvalue)                                 # <<<<<<< Remove


    # Find the index list of non-zero values within the input signal
    indices_nonzeros = np.nonzero(openvalue != 0)[0]
    # Find stretches of zero by at least window_size length
    if len(indices_nonzeros) > 0:   # If any nonzero region identified, all the zeros around it will be stored in an array

        if indices_nonzeros[0] != 0:    # If the region of interest started by 0 values, it will be reported in the zero regions
            initial = [0,indices_nonzeros[0]-1]
        else:
            initial = False

        if indices_nonzeros[-1] != len(openvalue)-1:   # If the region of interest ends by 0 values, it will be reported in the zero regions
            # zeros = np.insert(zeros, len(zeros), [int(indices_nonzeros[-1]),int(len(Mat))], axis=0)           <<<<<<<<<< Remove
            final = [indices_nonzeros[-1]+1, len(openvalue)]
        else:
            final = None
        RZ = EX.ranges(indices_nonzeros, userdefined_windowsize, initial, final)   # array of relative cooridnates
        print(RZ)
        zeros = [(line[0], str(i[0]+int(line[1])), str(i[1]+int(line[1]))) for i in RZ]
    # If no signal found in the region of interest, entrie of that will be reported as zero region
    else:
        zeros = [tuple(line)]                           # nonzero value could not find, all the interval defined as zero filled
    Z.extend(zeros)                                                              # array of averaged value arrays
del openvalue
Zint = pybedtools.BedTool(list(Z))                     # Zero-fill coordinates in bed format

if BL != []:
    BLRemoved = INint.subtract(BL)
else:
    BLRemoved = INint

INTERVAL = BLRemoved.subtract(Zint)                    # List of coordinates to be processed
if SaveInterval_file:
    pybedtools.BedTool(INTERVAL).saveas(SaveInterval_file)
    print('Interval exported : ',SaveInterval_file)

sys.stdout.write('\r Exclusion completed.           \n')

# -----------End of excluder / blacklist----------------


# --------------------Denoise---------------------------
# -----------------Filter checker-----------------------
# Store filter signal curve with the specified size within an array
try:
    ScipySignal=eval('scipy.signal.%s(%s)' % (d_filter,WSize)) # Cross the scaling factor
    low_resolution_Signal=eval('scipy.signal.%s(%s)' % (d_filter,PWS))
    # Ex_resolution_Signal=eval('scipy.signal.%s(%s)' % (d_filter,ExWS))
    # Apply different scales on the signal
except AttributeError:
    sys.stderr.write("\n# Filter error: The filter name "+ d_filter +" for the following command is not valid!\n -  scipy.signal." + d_filter)
    sys.exit()
# # --------------scalingfactor-------------------------
# if ScaFa:
#     pass
# else:
#     ScaFa = 1 # One scaling Factor for entire the genome must be calculated here
# --------------output validation and premission check -------------------------

InHighValue, InLowValue, InExtreamValue = [], [], []
if GoToSegmentation:
    for line in INTERVAL:
        convolved = Open_file.values(line.chrom, line.start, line.end, numpy=True)[::Step_size]
        low_resolution_convolve = Open_file2.values(line.chrom, line.start, line.end, numpy=True)[::Step_size]
else:
    for line in INTERVAL:
        openvalue = Open_file.values(line.chrom, line.start, line.end, numpy=True)
        openvalue[np.isnan(openvalue)] = 0
        Average_value=EX.AverageofRegion(openvalue, Step_size)                               # array of mean over stepsize values
        if len(openvalue) < Step_size*20:
            d_warning=True
        # Combine signals derived from the input file for each step's value by filter with specific windows size
        convolved = scipy.signal.fftconvolve(Average_value, ScipySignal, mode="same")
        # convolved=[i if (i > 0) else 0.0 for i in convolved]
        # This process may takes time, so one message informs the process stage
        # sys.stdout.flush()
        sys.stdout.write('\r '+line.chrom+'                     ')
        # Write the combined values for each interval to the output file
        OUT_file.addEntries(line.chrom, line.start, values=convolved, span=Span_size, step=Step_size)
        InHighValue.append(convolved)
        # Low resolution
        low_resolution_convolve = scipy.signal.fftconvolve(Average_value, low_resolution_Signal, mode="same")
        InLowValue.append(low_resolution_convolve)
        # Ex_resolution_convolve = scipy.signal.fftconvolve(Average_value, Ex_resolution_Signal, mode="same")
        # InExtreamValue.append(Ex_resolution_convolve)
        low_resolution_output.addEntries(line.chrom, line.start, values=low_resolution_convolve, span=Span_size, step=Step_size)
    sys.stdout.write('\r Saving bigWig file...')
    OUT_file.close()
    low_resolution_output.close()
    sys.stdout.write('\r Denoising completed. \n')# in',"%.5f seconds" % (time.time()-s_time))

# -----------------WARNING--------------------
# If the ratio of region/interval size and step size is not ideal, a warning message prints for the user
if d_warning==True:
    sys.stderr.write('\n# Warning: The defined Span=%s is too large. It could cause an error in the results because of the number of samples. It is recommended to use the smaller step.\n - ' % Span_size)

# --------------End of denoise----------------


# ----------------Segmentaion-----------------

# For each interval the value reads and stored in a varible
Wazowski, Salivan, Henry, counter = [],[],[], 0
for line in INTERVAL:
    sys.stdout.write('\r '+line.chrom+'       ')
    Sub_region = SG.high_resolution(line.chrom, line.start, line.end, Step_size, InHighValue[counter])
    Region = SG.low_resolution(line.chrom, line.start, line.end, Step_size, InLowValue[counter], TEMPFolder)
    # Extream_region = SG.low_resolution(line.chrom, line.start, line.end, Step_size, InExtreamValue[counter], TEMPFolder)
    if Sub_region is not None:
        Wazowski.extend(Sub_region)
    if Region is not None:
        Salivan.extend(Region)
    # if Extream_region is not None:
    #     Henry.extend(Extream_region)
    counter=counter+1
pybedtools.BedTool(list(Wazowski)).saveas(str(TEMPFolder)+'Wazowski.bed')
pybedtools.BedTool(list(Salivan)).saveas(str(TEMPFolder)+'Salivan.bed')
# pybedtools.BedTool(list(Henry)).saveas(str(TEMPFolder)+'Henry.bed')


# # find address of source file
# saveto = os.path.dirname(denoised_data)
# if saveto and saveto[-1] != '/':
#     saveto=str(saveto)+'/'

sys.stdout.write("\r Execution times %.5f seconds\n" % (time.time()-dsp_time))
