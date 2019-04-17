#!/usr/bin/env python
valid=0
read_window_size=50
#------------code starts from here -----------
import pyBigWig, scipy.signal, sys, argparse, re, os
from pathlib import Path
import numpy as np
from pybedtools import BedTool
import pybedtools as bt

#-----------library import---------------------
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=("""\
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
    input         Input files must be in BigWig file format.
    output        A BigWig file that could be explored by IGV.

[ Interval Definitions ]
    region        Apply filter on a specific region defined by chromosome name and location.
    interval      Cover a list of regions using BED files.
    entire        Apply filter on the entire BigWig input file.

[ Optional Parameters ]
    filter        Scipy signal filter selection change the pattern of applied signal (Blackman by default selected)
    size          Window size of the filter in the basepair unit.
    step          The distance between the start of one region and the end of the previous region.
    span          The fixed distance from the start position which defines a selected region.

[ General help ]
    help        Print this help menu.
    version     What version of bedtools are you using?.

"""))
parser.add_argument("-i", "--input", help="Define the name and location of the BigWig file as input. e.g. -i data/File.bw ",required=True)
parser.add_argument("-o", "--output", help="Define the name and the location of the output BigWig file.", required=True)
parser.add_argument("-f", "--filter", help="Insert the filter name of Scipy Signal filter list. Blackman filter selected by default", required=False, default="blackman")
parser.add_argument("-s", "--size", type=int, help="Define the window size of applying filter.", required=False, default="65536")
parser.add_argument("-r", "--region", help="A single section of input file could be defined as chromosomename:startindex:endindex")
parser.add_argument("-l", "--interval", help="A list of regions could be defined as a BED file format")
parser.add_argument("-S", "--step", type=int, help="The distance between the end point of each region from the relevant start point as basepair", default=50)
parser.add_argument("-e", "--entire", help="Force program to read all indexes within the selected region")
parser.add_argument("--list-filters", help="Print all the possible signal filters.")
parser.add_argument('-V', '--version', action='version', version="%(prog)s (Version 1.0)")

args = parser.parse_args()
d_inputfile = args.input
d_outputfile = args.output
d_filter = args.filter
d_size = int(args.size)
d_region = str(args.region)
d_step = int(args.step)
d_interval = args.interval
#-----------input definitions----------------------
my_file = Path(d_inputfile)
if my_file.is_file():
    pass
else:
    print("File \"" + d_inputfile+ "\" is not exist. Check the file name and it\'s path.")
    exit()
#-----------input validation----------------------
def WinfuncChecker():
    try:
        tt=eval('scipy.signal.%s' % d_filter)
    except AttributeError:
        print("The filter name '%s' for the following command is not valid!" % d_filter)
        print("scipy.signal.%s " % d_filter)
        exit()
WinfuncChecker()
#--------------filterchecker------------------------
d_size = d_size >> int(np.log2(read_window_size))
#--------------size of window corrector-------------
d_open = pyBigWig.open(d_inputfile)
d_output = pyBigWig.open(d_outputfile, "w")
d_output.addHeader(list(d_open.chroms().items()))
d_signal = eval('scipy.signal.%s(%s)' % (d_filter,d_size))
#--------------common part-------------------------
templist=list(d_open.chroms())
with open("names.txt", 'w') as file_handler:
    for item in templist:
        file_handler.write("{}\n".format(item))
#--------------names.txt------------------------
d_arglist = sys.argv
if '--interval' in d_arglist:
    valid=2
elif '-l' in d_arglist:
    valid=2
if '--region' in d_arglist:
    valid=1
elif '-r' in d_arglist:
    valid=1
del argparse

#--------------optional switiches------------------

if valid==1:
    if re.search('\w{1,12}:\d{1,12}:\d{1,12}$', d_region) == None:
        print("The region of interest must follow the standard structure ChromosomeName:StartBaseIndex:EndBaseIndex")
        exit()

    x = re.match('\w{1,10}:', d_region)
    xx = re.search(':\d{1,10}:', d_region)
    xxx = re.search(':\d{1,14}$', d_region)
    d_regname=x.group(0)
    d_regname=d_regname [:-1]
    d_regs=xx.group(0)
    d_regs=d_regs [1:-1]
    d_rege=xxx.group(0)
    d_rege=d_rege [1:]
    if d_regname in templist:
        if (int(d_open.chroms(d_regname)) >= int(d_rege)) :
            d_openvalue = eval('d_open.values("%s", %s, %s, numpy=True)[::d_step]' % (d_regname,d_regs,d_rege))
            d_convolve = scipy.signal.fftconvolve(d_openvalue, d_signal, mode="same")
            #d_output.addEntries(d_regname, int(d_regs), values=d_convolve, span=50, step=d_step)
            d_output.addEntries(d_regname, int(d_regs), ends=int(d_rege), values=d_convolve, span=50, step=d_step)
            d_output.close()
        else:
            print('Interval definition is incorrect. The length of the chromosome ' + str(d_regname) + ' is ' + str(d_open.chroms(d_regname)) + '. Please correct it and try again.')
            exit()
    else:
        print('The chromosome \"'+ d_regname +'\" is not present in the input file.')
        exit()

#--------------region refinement-------------------
elif valid==2:
    if d_filter!='hanning':
        bed_file = Path(d_interval)
        if bed_file.is_file():
            sites = bt.BedTool(d_interval).sort(g = 'names.txt').merge(d=d_step-1, c=5, o='sum').saveas('temp.bed')
            os.remove('names.txt')
            del (sites)
            with open('temp.bed')as d_intline:
                for line in d_intline:
                        L = line.strip().split()
                        d_openvalue = eval('d_open.values("%s", %s, %s, numpy=True)[::d_step]' % (L[0],L[1],L[2]))
                        d_convolve = scipy.signal.fftconvolve(d_openvalue, d_signal, mode="same")
                        d_output.addEntries(str(L[0]), int(L[1]), values=d_convolve, span=50, step=d_step)
                        #--------------output processings---------------------
                print('\n')
                os.remove('temp.bed')
                d_output.close()
                del (d_open, d_convolve, sys, re, os, time)
        else:
            print("File \"" + d_interval + "\" is not exist. Check the BED file and it\'s path.")
            exit()
    else:
        exit()
#-------------- interval refinement ------------
else:
    for line in templist:
        d_openvalue = d_open.values(line, 1, d_open.chroms(line), numpy=True)[::d_step]
        d_convolve = scipy.signal.fftconvolve(d_openvalue, d_signal, mode="same")
        print("- analyzing chromosome "+line+" from begining to "+str(d_open.chroms(line)))
        d_output.addEntries(line, 1, ends= d_open.chroms(line), values=d_convolve, span=50, step=d_step)
#-------------- entire chromosome ------------
