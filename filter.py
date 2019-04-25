# !/usr/bin/env python3.7
valid=0                                                                                                # Variable definitions
read_window_size=50                                                                                    # Windows size that applied on the BW files
# ------------code starts from here ---------------------------------------------------------------------
import pyBigWig, scipy.signal, sys, argparse, re, os, tempfile      
from pathlib import Path
import numpy as np
from pybedtools import BedTool
import pybedtools as bt
# -----------library import-----------------------------------------------------------------------------
def AVERAGE(segnumber, step, value, zero):                                                            # Improve Picking Accuracy
    d_matrix = np.matrix(value)                                                                       # Define a matrix with the values of BigWig file
    d_matrix = np.pad(d_matrix,[(0,0),(0,zero)], mode='constant', constant_values=0.000000001)        # fill the free cells of the matrix with a small and unique values
    d_matrix = np.reshape(d_matrix,(segnumber, step))                                                 # reshape the matrix from the one dimensional to step size rows and multiple columns
    masked = np.ma.masked_equal(d_matrix, 0.000000001)                                                # Create a mask and remove the value of the added cells
    # print(masked)
    Result=masked.mean(axis=1)                                                                        # Mean of the rows
    return(Result)                                                                                    # Returns the array of the means
# -----------function Definitions-----------------------------------------------------------------------
class _ListAction(argparse.Action):                                                                   # Define a class for the optional parameter "--list-filters" 
    def __init__(self,option_strings,dest=argparse.SUPPRESS,default=argparse.SUPPRESS,help=None):
        super(_ListAction, self).__init__(option_strings=option_strings,dest=dest,default=default,nargs=0,help=help)
    def __call__(self, parser, namespace, values, option_string=None):
        print("Avalable filters for Scipy Signal: \n hann, hamming, blackman")
        parser.exit()
    # -----------list available filters-----------------------------------------------------------------
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=("""
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # #  dsptool is a tool for ChIP-seq data analysis # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

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
    entire              Apply filter on the entire BigWig input file (One base resolution).

[ Optional Parameters ]
    filter              Scipy signal filter selection change the pattern of applied signal (Blackman by default selected)
    size                Window size of the filter in the basepair unit.
    step                The distance between the start of one region and the end of the previous region.
    span                The fixed distance from the start position which defines a selected region.

[ General help ]
    help                Print this help menu.
    version             What version of bedtools are you using?.
    filter list         Print all the available scipy signal filter name that could be apply on input signal. 

"""))
parser.add_argument("-i", "--input", help="Define the name and location of the BigWig file as input. e.g. -i data/File.bw ",required=True)
parser.add_argument("-o", "--output", help="Define the name and the location of the output BigWig file", required=True)
parser.add_argument("-f", "--filter", help="Insert the filter name of Scipy Signal filter list. Blackman filter selected by default", required=False, default="blackman")
parser.add_argument("-s", "--size", type=int, help="Define the window size of applying filter", required=False, default="65536")
parser.add_argument("-r", "--region", help="A single section of input file could be defined as chromosomename:startindex:endindex")
parser.add_argument("-l", "--interval", help="A list of regions could be defined as a BED file format")
parser.add_argument("-S", "--step", type=int, help="The distance between the end point of each region from the relevant start point as basepair", default=50)
parser.add_argument("--span", type=int, help="The intended length from the beginning of each step", default=50)
parser.add_argument('-e', "--entire", action='store_true', help="Force program to read all indexes within the selected region")
parser.add_argument("-L", "--list-filters", action=_ListAction, help="Print all the possible signal filters")
parser.add_argument("-V", "--version", action='version', version="%(prog)s (Version 1.0)")
args = parser.parse_args()          # Use "args" instead of the subcommand "parse_args" 
d_inputfile = args.input            # Store the results of switches on the variables
d_outputfile = args.output
d_filter = args.filter
d_size = int(args.size)
d_region = str(args.region)
d_step = int(args.step)
d_span = int(args.span)
d_interval = args.interval
# -----------input definitions----------------------
if args.entire == True:             # If user call "--entire", resolution increases to maximum (one base) 
    d_step=1
    d_span=1
my_file = Path(d_inputfile)         # Location of the input file stores on "my_file" variable
if my_file.is_file():               # Check the path/file presence or not
    pass
else:
    print("File \"" + d_inputfile+ "\" is not exist. Check the file name and it\'s path.")  # --------->>>>>>>>>>>>>> STD
    exit()
# -----------input/step validation----------------------
d_size = d_size >> int(np.log2(read_window_size))       # Window size of the signal resized regaurd the input file
# --------------size of window correction-------------
if d_filter=='hanning':              # Prevent to continue with old filter name which arise the internal error from scipy pakage
    print("`hanning` is deprecated, use `scipy.signal.windows.hann` instead!")
    exit()
try:                                 # Apply the filter name and windows size, if there is any problem with them, a clean exit with a message that points to the filter name will show
    d_signal=eval('scipy.signal.%s(%s)' % (d_filter,d_size))
except AttributeError:
    print("The filter name '%s' for the following command is not valid!" % d_filter)
    print("scipy.signal.%s " % d_filter)
    exit()
# --------------filterchecker------------------------
d_open = pyBigWig.open(d_inputfile)     # define the call the input bigwig file by pyBigWig library function on a variable
try:                                    # check if the output file cannot be created because of the output folder does not exist or permissions an error message with clean exit will occur
    d_output = pyBigWig.open(d_outputfile, "w")
    d_output.addHeader(list(d_open.chroms().items()))
except IOError:
    print ("Error: can\'t find path or write data")
    exit()
else:
    print ("File %s created successfully" % d_outputfile)
# --------------output validation-------------------------
templist=list(d_open.chroms())          # Define the temp list by all the chromosme names that present with in the input file (BigWig)
tmp = tempfile.NamedTemporaryFile()     # define a temperory file and put the list of the chromosome name inside it, because the Bedtools does not accept the list as a variable
with open(tmp.name, 'w') as _tmp:
    for line in templist:
        _tmp.write(str(line+"\n")) 
# --------------names.txt------------------------
d_arglist = sys.argv                    # Save all the arguments that are imported by user within a variable
if '--interval' in d_arglist:           # If call the interval, the locations will read from the Bed file 
    valid=2
elif '-l' in d_arglist:
    valid=2
if '--region' in d_arglist:             # If the used address an specific region, the interval will be ignored
    valid=1
elif '-r' in d_arglist:
    valid=1
del argparse
# --------------optional switiches------------------
if valid==1:
    if re.search('\w{1,12}:\d{1,12}-\d{1,12}$', d_region) == None:              # Extract the chromosome name, start and end lication from the user input
        print("The region of interest must follow the standard structure ChromosomeName:StartBaseIndex-EndBaseIndex")
        exit()

    _x = re.match('\w{1,10}:', d_region)
    _xx = re.search(':\d{1,10}-', d_region)
    _xxx = re.search('-\d{1,14}$', d_region)
    d_regname=_x.group(0)
    d_regname=d_regname [:-1]
    d_regs=_xx.group(0)
    d_regs=d_regs [1:-1]
    d_rege=_xxx.group(0)
    d_rege=d_rege [1:]
    d_rege=int(d_rege)+1            # The end location could have value which should keep in mind
    if d_regname in templist:       # If the request chromosome name present in the input file process starts
        if (int(d_open.chroms(d_regname)) >= (d_rege)) :        # If the region of interest present in the input file goes ahead
            d_openvalue = d_open.values(d_regname, int(d_regs), d_rege, numpy=True)# [::d_step]     # ????????????
            d_length=len(d_openvalue)
            d_segnum=(d_length/d_step)
            if d_segnum < 20:               # If the size of step is too large for the region of interest, in a way that less that 20 steps called, the warning appears and program continues 
                print('Warning: The defined "span" is too large for this interval, it could cause an error in the results because of the number of samples. It is recommended to use the smaller step.')
            if d_length % d_step == 0:      # Find the number of required cells to fill the matrix by AVERAGE function 
                d_dif=0
                d_segnum=int(d_segnum)
            else:
                d_segnum=int(d_segnum)+1
                d_dif=(d_segnum*d_step)-d_length

            d_averages=AVERAGE(d_segnum,d_step,d_openvalue,d_dif)           # Calclate the average of values relative to the each step and store the result as an array
            # --------------length / Step------------------------
            d_convolve = scipy.signal.fftconvolve(d_averages, d_signal, mode="same")        # Combine signals derived from the step's value by filter with specific windows size
            print(d_regname, int(d_regs), int(d_rege), d_span, d_step)      # Write the combined values for the region of intreset to the output file 
            d_output.addEntries(d_regname, int(d_regs), ends=int(d_rege), values=d_convolve, span=d_span, step=d_step)      # Write the combined values for the region of intreset to the output file 
            d_output.close()
        else:
            print('Interval definition is incorrect. The length of the chromosome ' + str(d_regname) + ' is ' + str(d_open.chroms(d_regname)) + '. Please correct it and try again.') # =====>>> STD
            exit()
    else:
        print('The chromosome \"'+ d_regname +'\" is not present in the input file.')
        exit()
# --------------region refinement-------------------
elif valid==2:                      # If the user defined an Bed file these codes run
    bed_file = Path(d_interval)     # Store the file name and location in one variable
    with open(tmp.name) as _tmp:    # Call temprory file 
        if bed_file.is_file():      # Check if the Bed file present, continue
            sys.stdout.write('')    # Print a standard output
            sys.stdout.flush()
            if args.entire == True:         # If user use the "-e" switch, Step changes to 1 
                sites = list(bt.BedTool(d_interval).sort(g = _tmp.name))        # Interval data derived from Bed file sorts by name
            else:
                sites = list(bt.BedTool(d_interval).sort(g = _tmp.name).merge(d=d_step-1))          # Interval data derived from Bed file sorted and merged if their steps ooverlapped
            for line in sites:                  # For each interval the value reads and stored in a varible
                _temp=str(line)
                L = _temp.strip().split()
                d_openvalue = d_open.values(L[0], int(L[1]), int(L[2]), numpy=True)# [::d_step]        # ????
                d_length=len(d_openvalue)
                d_segnum=(d_length/d_step)      # The number of columns of matrix is equal to the step-size, so the number of rows are the length devided by step-size
                if d_length % d_step == 0:      # If divmod is zero it means the length of the string is fit to the matrix
                    d_dif=0
                    d_segnum=int(d_segnum)
                else:
                    d_segnum=int(d_segnum)+1    
                    d_dif=(d_segnum*d_step)-d_length
                d_averages=AVERAGE(d_segnum,d_step,d_openvalue,d_dif)
                d_convolve = scipy.signal.fftconvolve(d_averages, d_signal, mode="same")[::d_step]      # Combine signals derived from the input file for each step's value by filter with specific windows size
                sys.stdout.write('\r'+'Working on '+L[0])       # Update the standard output
                sys.stdout.flush()
                d_output.addEntries(str(L[0]), int(L[1]), values=d_convolve, span=d_span, step=d_step)  # Write the combined values for each interval to the output file 
                # --------------output processings---------------------
            sys.stdout.write('\r'+'Action completed                           ')        # Update the standard output
            print('\n')
            d_output.close()
            del (d_open, d_convolve, sys, re, os)
        else:
            print("File \"" + d_interval + "\" is not exist. Check the BED file and it\'s path.")          # ===>>> STD
            exit()
# -------------- interval refinement ------------
else:
    for line in templist:     # For each interval the value reads and stored in a varible for the entire input data
        d_openvalue = d_open.values(line, 1, d_open.chroms(line), numpy=True)[::d_step]
        d_convolve = scipy.signal.fftconvolve(d_openvalue, d_signal, mode="same")
        print("- analyzing chromosome "+line+" from begining to "+str(d_open.chroms(line)))
        d_output.addEntries(line, 1, ends= d_open.chroms(line), values=d_convolve, span=d_span, step=d_step)
# -------------- entire chromosome ------------
# -----------------------------------------------END OF DENOISING PART---------------------------------
