#!/usr/bin/env python

from __future__ import with_statement
__authors__ = ('Shahram Saghaei (saghaei.shahram@hsr.it), '
               'Davide Cittaro (cittaro.davide@hsr.it)')
__version__ = '1.4.18'
__date__ = '9 March 2021'

import time, sys, os
dsp_time = time.time()
sys.stdout.write('\n Please wait ...       ')

# Color codes which are useing in the message print
RED   = "\033[1;31m"
GREEN = "\033[0;32m"
RESET = "\033[0;0m"

# -----------modules check------------------------------------------------------
# Check the availability of all the packages required for running the program
try:
    import pyBigWig, sys, argparse, os, pathlib
    from pybedtools import BedTool
    import scipy.signal
    from pathlib import Path
    from itertools import groupby
    import numpy as np
    import pandas as pd
    import sklearn.mixture
    import matplotlib.pyplot as plt
    from package import excluder as EX
    from package import segment as SG
    from package import gmm as GMM
    from package import kmean as KM
    from package import otsu as OTSU


except ImportError:
    raise ImportError(sys.exit("""\n# The following packages are required, please install them:

\033[1;36m-------------------------------------------------------------------------------
|\033[0;0m package                    \033[1;36m|\033[0;0m        Version       \033[1;36m|\033[0;0m          Build          \033[1;36m|
|----------------------------|----------------------|-------------------------|
|\033[0;0m pyBigWig                   \033[1;36m|\033[0;0m        0.3.17        \033[1;36m|\033[0;0m       py36h8c929e3_2    \033[1;36m|
|\033[0;0m pybedtools                 \033[1;36m|\033[0;0m        0.8.1         \033[1;36m|\033[0;0m       py36h0c56d2d_2    \033[1;36m|
|\033[0;0m sklearn                    \033[1;36m|\033[0;0m        0.24.0        \033[1;36m|\033[0;0m       py36hed11f80_0    \033[1;36m|
|\033[0;0m pathlib                    \033[1;36m|\033[0;0m        1.0.1         \033[1;36m|\033[0;0m       py36h9f0ad1d_3    \033[1;36m|
|\033[0;0m matplotlib                 \033[1;36m|\033[0;0m        3.3.3         \033[1;36m|\033[0;0m       py36h79c6626_0    \033[1;36m|
|\033[0;0m argparse                   \033[1;36m|\033[0;0m        1.4.0         \033[1;36m|\033[0;0m       py36_0            \033[1;36m|
-------------------------------------------------------------------------------\033[0;0m\n
DSP-tool should have a package folder containing excluder.py, segment.py, gmm.py, OTSU.py, and kmean.py.
in order to preventing the package versions conflict, it is suggested to create a conda virtual environment. e.g.
conda create -n DSP-tool -c bioconda pybigwig pybedtools scikit-learn pathlib matplotlib argparse seaborn
"""))

d_warning=False

# -----------available filters list---------------------------------------------
# Define a class for the optional parameter "--list-filters"
class _ListAction(argparse.Action):
    def __init__(self,option_strings,dest=argparse.SUPPRESS,default=argparse.SUPPRESS,help=None):
        super(_ListAction, self).__init__(option_strings=option_strings,dest=dest,default=default,nargs=0,help=help)
    def __call__(self, parser, namespace, values, option_string=None):
        sys.stdout.write("\n Avalable filters Signal for denoising process: \n - hann, hamming, blackman, blackmanharris, triang, parzen, cosine, nuttall, bohman, barthann, bartlett, morlet, and flattop. \n - The Blackman filter is selected by defult.\n \n")
        parser.exit()

# -----------input definitions--------------------------------------------------
# Required and optional parameters definition usnig "argparse" package
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=("""
##############################################################################
################ DSP-tool is a tool for ChIP-seq data analysis ################
##############################################################################
Version:   v1.4.18
About:     developed in the HSR.it
Code:      https://github.com/dawe/DSP-tool
Mail:      saghaei.shahram@hsr.it

The minimum needed arguments:     DSP-tool -i input -o output -r region [options]
This program is made with python 3.6.

\033[1;36m[ Input Parameters ]\033[0;0m
    input file          Input files with bigWig file format as a raw signal (non-denoised).
    denoised file       A previously denoised bigWig file as input to bypass the denoising step.
\033[1;36m[ Output Parameters ]\033[0;0m
    output file         Either a '*.bw' file for denoised signal or a '.bed' file for enriched region.
    temp folder         The location of a folder with read and write permission.
\033[1;36m[ Interval Definitions ]\033[0;0m
    region              Select a specific address by chromosome name, start and end coordinates.
    interval            Inspect a list of regions stored in a BED files.
    entire genome       Do not select any interval to check the entire input file.
    exclude list        Define a Bed file as black/gray list file to exclude for inspection.
    save interval       Save the inspecting interval in a '.bed' file.
\033[1;36m[ Optional Parameters ]\033[0;0m
    filter              Name of filter-signal that will be used for denoising the input signal.
    window size         The window size of the convolution process in the basepair unit.
    step size           The length of the reading interval from the input file.
    fold Size           The number of times to increase the window size in the low-resolution denoising process, with respect to the high-resolution.
    clustering method   Select on of the clustering methods Gaussian Mixture Models, KMeans or OTSU algorithm.
    components number   The number of components user suggests for clustering. It can be automatically computed by software.
\033[1;36m[ General help ]\033[0;0m
    help                Print this help menu.
    version             Verion of the DSP-tool.
    filter list         Print the name of all the available filter-signal in DSP-tool.
"""))
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument("-i", "--input", help="A non-processed bigWig file to input to the denoising process. Possible to run the whole program or denoising part exclusively. e.g. '-i path/file.bw'")
group.add_argument("-d", "--denoised", help="A previously denoised bigWig file as input for running the segmentation process exclusively e.g. '-d path/file.bw'")

parser.add_argument("-o", "--output", help="Either a bigWig file for denoised signal or a Bed file for final identified boundaries could be the output of the program. It depends on which processing step(s) the user applies e.g. '-o path/file.bw' or '-o path/file.bed' ",required=True)
parser.add_argument("-temp", "--tempfolder", help="The location of a folder to store runtime temporary files with reading and writing permission.", default="")
parser.add_argument("-e", "--exclude", help="Define a Bed file as a BlackList file e.g. '-e path/blacklist.bed'")
parser.add_argument("-save", "--saveinterval", help="Save the inspection interval in a '.bed' file. e.g. '-e path/investigatedintervals.bed'")

parser.add_argument("-s", "--step", type=int, help="The step size that the denoising process applied by, the default step size is 50 bps.", default=50)
parser.add_argument("-H", "--highresolution", action='store_true', help="Force the program to define step and span equal to 1 to have the highest resolution as possible (One base resolution).")

parser.add_argument("-r", "--region", help="Select a specific address by chromosome name, start and end coordinates. '-r chr1:1-10000000'")
parser.add_argument("-l", "--interval", help="A bed file contains the list of regions with three columns: chromosome name, start position, and end position.")
parser.add_argument("-w", "--windowsize", type=int, help="The window size of the convolution process in the basepair unit. The default window size is 16 Kb", required=False, default="16384")
parser.add_argument("-f", "--filter", help="Name of filter-signal that will be used for denoising the input signal. 'blackmanharris' filter selected by default", required=False, default="blackmanharris")
parser.add_argument("-list", "--filterlist", action=_ListAction, help="Print the name of all the available filter-signal in DSP-tool.")
parser.add_argument("-fold", "--foldtimes", type=int, help="Number of folds of low resolution window size in compare to high resolution denoising window size.", required=False, default="2")
parser.add_argument("-c", "--componentsnumber", help="The number of components use in GMM Thresholding.", required=False, default="auto")
parser.add_argument("-m", "--method", help="The thresholding method which could be KMean, GMM or OTSU.", required=False, default="kmean")
parser.add_argument("-v", "--version", action='version', version="%(prog)s (\033[1;36mVersion v1.4.18\033[0;0m)")


# Use "args" instead of the subcommand "parse_args"
args = parser.parse_args()
# Store the inserted values for switches on the variables
source_data = args.input
denoised_data = args.denoised
d_filter = args.filter.lower()
userdefined_windowsize = int(args.windowsize)
Region_of_interest = args.region
Span_size = Step_size = int(args.step)
user_defined_interval = args.interval
BlackList_file = args.exclude
TEMPFolder = args.tempfolder
foldtimes = args.foldtimes
cn = args.componentsnumber
output = args.output
SaveInterval_file = args.saveinterval
METHOD = args.method.lower()

# -----------end of parameters--------------------------------------------------

# Save directory checking and make it if it is not present
saveloc = os.path.dirname(os.path.realpath(output))
output = os.path.realpath(output)
try:
    os.makedirs(saveloc)
except:
    pass

if METHOD not in ["kmean", "gmm", "otsu", "none"]:
    exit(RED+'\n* Error: The defined thresholding method "'+METHOD+'" is not valid. Avalable methods are KMean, GMM or OTSU. Please check and try again.\n'+RESET)

# -----------tempdir------------------------------------------------------------

# Temp directory checking and make it if it is not present
if TEMPFolder != "":
    if TEMPFolder[-1] != '/':
        TEMPFolder=str(TEMPFolder)+'/'
    if TEMPFolder[0] == '/':
        TEMPFolder=str(TEMPFolder[1:])

else:
    TEMPFolder = saveloc+"/temp/"
try:
    os.makedirs(TEMPFolder)
except:
    pass

# -----------running modules decision-------------------------------------------

# Result file could be either a smoothed signal or an interval file. If the
# smoothed file ".bw" requested by the user, the "excluder" module and
# "denoising" modules would be applied on the input file. While if interval file
# ".bed" is requested, two other module "segmentation" and "clustering" will be
# applied after denoising process.
if output.endswith(".bw"):
    # if a bigWig file is the output file
    DoutH = output
    save_path_file = os.path.abspath(DoutH)[:-3:]
    filename = save_path_file.split("/")[-1]

    # if a bigWig file is the input file
    if source_data:
        Din = os.path.realpath(source_data)
        # input
        try:
            Open_file = pyBigWig.open(Din)
            ChrLIST = Open_file.chroms()
            Header = list(Open_file.chroms().items())
        except:
            sys.stderr.write(RED+'\n* Error: Input file \"' + Din + '\" was not found or corrupted.\n'+RESET)
            sys.exit()

        # output
        try:
            DoutH = output
            DoutL = save_path_file + '_L.bw'
            OUT_file = pyBigWig.open(DoutH, "w")
            OUT_file.addHeader(Header)
            OUT_file_L = pyBigWig.open(DoutL, "w")
            OUT_file_L.addHeader(Header)
            # no bed needed
        except IOError:
            sys.stderr.write(RED+"\n* Error:  The PATH/FILE \""+ DoutH +"\" has not write permission, please check the path and try again. \n"+RESET)
            sys.exit()

    elif denoised_data:
        exit(RED+"\n* Error: The output file of the segmentation process is a bed file. Please define an appropriate file name for the output.\nYou defined "+denoised_data +RESET)

elif output.endswith(".bed"):
    # if a bed file is the output file
    save_path_file = os.path.abspath(output)[:-4:]
    filename = save_path_file.split("/")[-1]

    # if a non-processed bigWig file is the input file
    if source_data:
        Din = os.path.realpath(source_data)
        Pin = Din
        try:
            Open_file = pyBigWig.open(Din)
            ChrLIST = Open_file.chroms()
            Header = list(Open_file.chroms().items())
        except:
            sys.stderr.write(RED+'\n* Error: Input file \"' + Din + '\" was not found or corrupted.\n'+RESET)
            sys.exit()
        try:
            DoutH = TEMPFolder+filename+ '_denoised_H.bw'
            DoutL = TEMPFolder+filename+ '_denoised_L.bw'
            OUT_file = pyBigWig.open(DoutH, "w")
            OUT_file.addHeader(Header)
            OUT_file_L = pyBigWig.open(DoutL, "w")
            OUT_file_L.addHeader(Header)
        except IOError:
            sys.stderr.write(RED+"\n* Error:  The PATH/FILE \""+ filename +"\" has not write permission, please check the path and try again. \n"+RESET)
            sys.exit()
        SinH, SinL, Sout = DoutH, DoutL, output

    # if a denoised-bigWig file is the input file
    elif denoised_data:
        if denoised_data.endswith(".bw"):
            SinH = denoised_data
            Pin = SinH
            Sout = output
            # Denoising process will be skipped (Din, Dout), because two
            # denoised files are used as input files
            if denoised_data.endswith("_H.bw"):
                Sin = os.path.abspath(SinH)[:-5:]
            else:
                Sin = os.path.abspath(SinH)[:-3:]
            SinL = Sin + "_L.bw"
            # input
            try:
                Open_file = pyBigWig.open(SinH)
                ChrLIST = Open_file.chroms()
                Header = list(Open_file.chroms().items())
            except:
                sys.stderr.write(RED+'\n* Error: Input file \"' + SinH + '\" was not found or corrupted.\n'+RESET)
                sys.exit()
            if os.path.exists(SinL):
                pass
            else:
                sys.exit(RED+'\n* Error: Input file \"' + SinL + '\" was not found or corrupted.\n'+RESET)
        else:
            sys.exit(RED+'\n* Error: Input file \"' + denoised_data + '\" is not ba bigWig file.\n'+RESET)


else:
    sys.exit(RED+'\n* Error: The output file of the system is either a bed file or a bw file. Please define an appropriate file name for the output.'+RESET)


# -----------window size/step size----------------------------------------------
# If the high-resolution option selected by the user. It is useful for the short
# peaks that need more resolutions. it means all the bases of the defined region
# will be called and processed.
if args.highresolution:
    Step_size, Span_size = 1, 1

# Window size and fold size
WSize = userdefined_windowsize
PWS = 2 ** (int(np.log2(userdefined_windowsize))+foldtimes)

if source_data:
    if Step_size > int(Header[0][1]):
        sys.stderr.write(RED+"\n* Error:  Window Size for low_resolution denoising is greater than chromosme size. \n "+RESET)
        sys.exit()


# -----------excluder/blacklist-------------------------------------------------
sys.stdout.write('\r Excluding ...')

if BlackList_file:
    try:
        # Import blacklist file, and merge the blacklist intervals if they have overlap
        BL = BedTool(BlackList_file).sort().merge()
    except:
        sys.stderr.write(RED+'\n* Error: Blacklist file was not found or corrupted.\n'+RESET)
        sys.exit()
else:
    BL =[]

# Highlighting the region of inspection
if Region_of_interest:
    target='Region'
    zone=Region_of_interest
elif user_defined_interval:
    target='Interval'
    zone=user_defined_interval
else:
    target='Entrie_genome'
    zone=''
# The "excluder" module converts the region of inspection to a list of intervals.
INPint = EX.intervalconverter(BL, target, ChrLIST, Step_size, zone) # list of intervals
# Convert the interval list to the Bedtool format. Thus Bedtool could be applied to it.
INint = BedTool(INPint)

#---------------Zero function---------------------------------------------------
# Very long zero-valued regions after convolution creates noisy regions which
# consume require computational power to discard them. Therefore, in order to
# prevent wasting time and computational power, it is better to exclude these
# regions from the inspecting regions.

Zetta=[]
if len(INPint) == 0:
    sys.stderr.write(RED+'\n* Error: The interval file is empty.\n'+RESET)
    exit()
for line in INPint:
    sys.stdout.write('\r\x1b[K Excluding '+line.chrom)
    try:
        openvalue = Open_file.values(line[0], int(line[1]), int(line[2]), numpy=True)
        openvalue=np.nan_to_num(openvalue)
    except:
        sys.stderr.write(RED+"\n* Error: pyBigWig package is not installed properly, please install it using conda.\n To assure it is installed properly, import pyBigWig in python, if you get '0' for 'pyBigWig.np' command, please install it using conda as following:\nconda install -c bioconda pybigwig\n"+RESET)
        exit()
    openvalue[np.isnan(openvalue)] = 0


    # Find the index list of non-zero values within the input signal
    indices_nonzeros = np.nonzero(openvalue != 0)[0]
    # Find stretches of zero by at least window-size length
    if len(indices_nonzeros) > 0:
        # If the region of interest started by 0 values, it will be reported in the zero regions
        if indices_nonzeros[0] != 0:
            initial = [0,indices_nonzeros[0]-1]
        else:
            initial = False

        # If the region of interest ends by 0 values, it will be reported in the zero regions
        if indices_nonzeros[-1] != len(openvalue)-1:
            final = [indices_nonzeros[-1]+1, len(openvalue)]
        else:
            final = False

        # identification of all the long gaps between the scaffolds and store them in an array.
        RZ = EX.ranges(indices_nonzeros, userdefined_windowsize, initial, final)
        try:
            zeros = [(line[0], str(i[0]+int(line[1])), str(i[1]+int(line[1]))) for i in RZ]
        except:
            zeros = [[line[0], str(int(RZ[0])+int(line[1])), str(int(RZ[1])+int(line[1]))]]

    # If no signal found in the region of interest, entire of that will be reported as a zero region
    else:
        zeros = [tuple(line)]
    # An array of zero-valued arrays
    Zetta.extend(zeros)
del openvalue

# Convert zero-fill regions into the bed format
Zettaint = BedTool(list(Zetta))

# Dissection of blacklist intervals from the inspecting intervals
if BL != []:
    BLRemoved = INint.subtract(BL)
else:
    BLRemoved = INint

# Dissection of zero-fill regions from the inspecting intervals
INTERVAL = BLRemoved.subtract(Zettaint)

sys.stdout.write(GREEN+'\r\x1b[K Exclusion operation completed.\n'+RESET)

# If the export of the inspection intervals is requested by the user, a ".bed"
# will be created in the selected location.
if SaveInterval_file:
    BedTool(INTERVAL).saveas(SaveInterval_file)
    print(' Interval exported: ',SaveInterval_file)

# -----------End-of-the-excluder/blacklist--------------------------------------


# --------------------Denoise---------------------------------------------------
# In order to remove the white-noise of the input signal, two parallel denoising
# process with two size of a filter-signal applies on input signal.

# -----------------Filter checker-----------------------------------------------
# At the same time evaluates the validity of the selected filter name and store
# filter-signal curve with the specified window size within an array.

try:
    # Cross the scaling factor
    ScipySignal=eval('scipy.signal.%s(%s)' % (d_filter,WSize))
    low_resolution_Signal=eval('scipy.signal.%s(%s)' % (d_filter,PWS))
except AttributeError:
    sys.stderr.write(RED+"\n* Error: The filter name \""+ d_filter +'\" is not valid! PLease check filter name.\n'+RESET)
    sys.exit()

# If the input files are already smoothed ".bw" files program will skip and
# starts the segmentation procedure.
InHighValue, InLowValue, = [], []
if denoised_data:
    pass
else:
    # Interval file reads line by line. Each line carries a chromosome name, and coordinates
    # of two sides of a fragment. The value of the corresponding region of this coordinate in
    # the bigWig file called and stored in an intensity array.
    for line in INTERVAL:
        sys.stdout.write('\r\x1b[K Denoising '+line.chrom)
        openvalue = Open_file.values(line.chrom, line.start, line.end, numpy=True)
        openvalue=np.nan_to_num(openvalue)

        # Downsampling of the intensity array according to the step size. The average of one partition
        # of the array with the length of step size will be representative of that fragment.
        Average_value=EX.AverageofRegion(openvalue, Step_size)
        if len(openvalue) < Step_size*20:
            d_warning=True

        # Fast Fourier transform convolution of the intensity array with a filter-signal makes it smooth.
        # The level of this denoising effect depends on the window size and shape of the filter-signal.
        # high resolution
        convolved = scipy.signal.fftconvolve(Average_value, ScipySignal, mode="same")

        # If a signal with a negative value uses as the filter-signal, the convolved signal could become
        # negative in some regions with a low level of amplitude.
        if (convolved.shape[0] > 0) and (min(convolved) < 0):
                # The numerical value of these negative points is very small and is corrected by subtracting it
                # from the whole signal.
                convolved = convolved-min(convolved)

        values0 = np.array(convolved, dtype=np.float64)
        # Write the combined values of each interval into the output files.
        try:
            OUT_file.addEntries(line.chrom, line.start, ends=line.end, values=values0, span=Span_size, step=Step_size)
        except:
            sys.stderr.write(RED+'\n* Error: The defined step size = '+ str(Step_size) +' is too large for the inspecting fragment.\n'+RESET)
            sys.exit()

        # Store the combined values of each interval into an array for the next processing step.
        InHighValue.append(convolved)

        # Low resolution
        low_resolution_convolve = scipy.signal.fftconvolve(Average_value, low_resolution_Signal, mode="same")
        InLowValue.append(low_resolution_convolve)


        OUT_file_L.addEntries(line.chrom, line.start, ends=line.end, values=low_resolution_convolve, span=Span_size, step=Step_size)
    sys.stdout.write('\r\x1b[K Saving bigWig file...')

    # Close the output ".bw" files for two denoised signals.
    OUT_file.close()
    OUT_file_L.close()

    # Create a config file stores all the settings used to run the program. This file will
    # be created in the temp folder.
    configfile = TEMPFolder+filename+'.cfg'
    conf = open(str(configfile),"w+")
    conf.write('Input_filename= '+str(Din))
    conf.write('\nDenoised_filename= '+str(DoutH))
    conf.write('\nBlackList_file= '+str(BL))
    conf.write('\nFoldtimes= '+str(foldtimes))
    conf.write('\nStep_size= '+str(Step_size))
    conf.write('\nWindowsize= '+str(userdefined_windowsize))
    conf.write('\nTempfolder= '+str(TEMPFolder))
    conf.close()

    sys.stdout.write(GREEN+'\r\x1b[K Denoising operation completed. \n'+RESET)

    if output.endswith(".bw"):
        sys.stdout.write("\r Execution times %.5f seconds\n" % (time.time()-dsp_time))
        exit()

# -----------------WARNING------------------------------------------------------
# # If the ratio of region/interval size and step size is not ideal, a warning message prints for the user
# if d_warning==True:
#     sys.stderr.write('\n# Warning: The defined Span=%s is too large. It could cause an \n* Error in the results because of the number of samples. It is recommended to use the smaller step.\n ' % Span_size)

# --------------End of denoise--------------------------------------------------


# ----------------Segmentaion---------------------------------------------------
# A coupled process including the extremum recognition and boundary determination
# of each peak which leads to peak calling of the signal.

sys.stdout.write('\r Segmenting ... ')

# Load the denoised bigWig files from the temporary folder or the defined folder
# by the user.
try:
    Load_bw_file = pyBigWig.open(SinH)
    Load_bw_file2 = pyBigWig.open(SinL)
except IOError:
    sys.stderr.write(RED+ "\n* Error:  The PATH/FILE \"" + SinH + "\" or \"" + SinL + "\" was not found, please check the path and try again. \n"+RESET)
    sys.exit()

# if the input files are denoised signals, the desired intervals of them load to
# the two arrays.
if denoised_data:
    InHighValue, InLowValue = [], []
    for line in INTERVAL:
        sys.stdout.write('\r\x1b[K Segmenting '+str(line.chrom))
        convolved = Load_bw_file.values(line.chrom, line.start, line.end, numpy=True)[::Step_size]
        convolved=np.nan_to_num(convolved)
        InHighValue.append(convolved)
        low_resolution_convolve = Load_bw_file2.values(line.chrom, line.start, line.end, numpy=True)[::Step_size]
        low_resolution_convolve=np.nan_to_num(low_resolution_convolve)
        InLowValue.append(low_resolution_convolve)

# Boundary tracing for each peak applies to denoised signals and all the
# candidate regions stores in two arrays.
Wazowski, Salivan, counter = [],[], 0
for line in INTERVAL:
    # Candidate regions in the signal smoothed with smaller window-size are more in numbers and
    # have shorter lengths. These regions considered as the sub-regions of the identified regions.
    Sub_region = SG.boundary_finder(line.chrom, line.start, line.end, Step_size, InHighValue[counter])

    # Candidate regions in the signal smoothed with larger window-size are less in numbers and
    # have longer lengths.
    Region = SG.boundary_finder(line.chrom, line.start, line.end, Step_size, InLowValue[counter])

    if Sub_region is not None:
        Wazowski.extend(Sub_region)
    if Region is not None:
        Salivan.extend(Region)
    counter=counter+1
sys.stdout.write('\r\x1b[K Segmentaion Saving ... ')

# List of the candidate regions stores in two bedfiles which allows bedtools to
# manipulate them.
Waz=TEMPFolder+filename+'_highresolution.bed'
Sal=TEMPFolder+filename+'_lowresolution.bed'

# List of the candidate regions exports in the temporary folder.
BedTool(list(Wazowski)).sort().saveas(Waz)
BedTool(list(Salivan)).sort().saveas(Sal)
sys.stdout.write(GREEN+'\r\x1b[K Segmentaion operation completed. \n'+RESET)

# --------------End of segmentation---------------------------------------------

# --------------Clustering------------------------------------------------------
#  Although a large amount of the white noise removed from the input signal, but
# still candidate peaks called from the signal have the source of the foreground
# and background. In order to discard the background derived peaks, a clustering
# step required.

sys.stdout.write(' Gaussian Mixture Model Clustering...')

gmmfile = TEMPFolder+filename+'_modified.bed'
CLT=TEMPFolder+'plt'

# On the candidate list with more than 10 regions, a clustering method applies.
if len(Wazowski) >= 10:

    # clustering using the Gaussian Mixture Model with the certain number of components
    # devides all intervals into two (foreground/background) groups. This action applies
    # only on the sub-regions. If user uses "none" for the clustering method, it will be
    # skipped.
    if METHOD == "gmm":
        GMM.evaluate(Load_bw_file, BedTool(Waz), gmmfile, CLT, cn)
        sys.stdout.write(GREEN+'\r\x1b[K Gaussian Mixture Model applied. \n'+RESET)
    elif METHOD == "otsu":
        OTSU.evaluate(Load_bw_file, BedTool(Waz), gmmfile, CLT)
        sys.stdout.write(GREEN+'\r\x1b[K OTSU thresholding method applied. \n'+RESET)
    elif METHOD == "kmean":
        KM.evaluate(Load_bw_file, BedTool(Waz), gmmfile, CLT, cn)
        sys.stdout.write(GREEN+'\r\x1b[K KMean thresholding method applied. \n'+RESET)
    elif METHOD == "none":
        print(GREEN+'\r\x1b[K Clustering skipped by the user. \n'+RESET)
        gmmfile = Waz
    A1, A2 = BedTool(gmmfile), BedTool(Sal)
else:
    print(GREEN+'\r\x1b[K Clustering skipped due to the low number of the candidate sub-regions. \n'+RESET)
    A1, A2 = BedTool(Waz), BedTool(Sal)

if len(A2) == 0:
    sys.exit(RED+"\n* Error: In the low-resolution examination, no peak found. Please re-execute the denoising procedure and decrease either Window size and/or fold size. Current parameter sizes are too large for the defined interval.\n"+RESET+" foldtimes = "+str(foldtimes) + "\n windowsize = "+str(userdefined_windowsize))

# The intersection of a region and sub-regions found in its corresponding
# coordination suggests that this region is likely to be an enriched region.
try:
    # First step of the intersection
    intersect1=A2.intersect(A1, wa=True, wb=True).groupby(c=(4,5,6), o=("distinct","min","max"))
    intersection = True
except:
    intersection = False
    sys.stdout.write(RED+'\n* Error: Intersection failed.\n'+RESET)
    output="\tfailed"

if intersection:
    section=[]
    for i in intersect1:
        section.append([i.name,i.score,i.strand])
    A3 = BedTool(section)
    b3 = TEMPFolder+'bed3.bed'
    # Second step of the intersection. The result of the intersection saves as a
    # Bed3 file.
    bed3=A3.intersect(A1, wa=True, wb=True).groupby(c=(5,6), o=("collapse","collapse")).saveas(b3)

    # A bed12 file is capable to show all the regions and their relative
    # sub-regions in a single file. This file could be used in a differencial
    # enrichment analysis of the ChIP-seq file.
    LINE=[]
    for i in bed3:

        LthickStart=[int(x) for x in (i.name).split(",")]
        LthickEnd=[int(x) for x in (i.score).split(",")]
        sortindex=SG.SortIdx(LthickStart)
        SthickStart=[LthickStart[i] for i in sortindex]
        SthickEnd=[LthickEnd[i] for i in sortindex]
        thickStart=",".join([str(j) for j in SthickStart])
        thickEnd=",".join([str(j) for j in SthickEnd])
        IblockSizes=list(map(int.__sub__, SthickEnd, SthickStart))
        blockSizes=",".join([str(j) for j in IblockSizes])
        Mstart=min(SthickStart)
        IblockStarts = list(map(lambda x: x - Mstart, SthickStart))
        blockStarts=",".join([str(j) for j in IblockStarts])
        LINE.append([i.chrom, i.start, i.end, Mstart, 1000, ".", i.start, i.end, "128,0,128", len(SthickStart), blockSizes, blockStarts])

    BedTool(LINE).saveas(output)

sys.stdout.write("\r Execution times %.5f seconds\n" % (time.time()-dsp_time))

# Reporting of all the files created by the run and their location in a table at
# the end of the run.
if source_data:
    print("""
----------------------------------------------------------------------------------
|                               Created Files
----------------------------------------------------------------------------------
|       Name        |                       Saved Location
----------------------------------------------------------------------------------
| input             | %s
| Enriched regions  | %s
| Enriched in Bed3  | %s
| Hi-res boundaries | %s
| Hi-res signal     | %s
| Lo-res boundaries | %s
| Lo-res signal     | %s
| Config file       | %s
| Plot folder       | %s
----------------------------------------------------------------------------------

"""% (Pin, output, b3, Waz, SinH, Sal, SinL, configfile, CLT))
