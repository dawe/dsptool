#!/usr/bin/env python

from __future__ import with_statement
__author__ = ('Shahram Saghaei (saghaei.shahram@hsr.it), '
              'Davide Cittaro (cittaro.davide@hsr.it)')
__version__ = '0.1.00'
__date__ = '12 Jan 2020'

import time, sys
dsp_time = time.time()

##############filter
if __name__ == "__main__":
    try:
        import pyBigWig, sys, argparse, re, os, tempfile, numpy, pybedtools, pathlib, gc ,skimage, numpy
        import scipy.signal
        from pathlib import Path
        from itertools import groupby
        from operator import itemgetter
        import package.filter as fl
        import package.segmentation as sg
    except ImportError:
        raise ImportError("\n#\n# Excluder error: Some required packages missing, please install them.\nRequired pakages:\n -  pyBigWig=0.3.17\n -  sys\n -  argparse\n -  re\n -  os\n -  tempfile\n -  numpy=1.17.4\n -  pybedtools=0.8.1\n -  pathlib\n -  gc \n -  scikit-image=0.15.0\n -  scipy=1.3.1\n -  pathlib \n -  itertools\n -  operator\n -  uliengineering=0.3.4 \n Also you can use dsp.yml file to create a conda enviroment and activate it, use following commands:\n -  conda env create -f dsp.yml \n and then \n -  conda activate dsptool\n")
        sys.exit()

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
    parser.add_argument("-p", "--peak", help="Define the name and the location of the output BigWig file. If no output file defined, the result will be stored in Peak.bw file in the running folder.")
    parser.add_argument("-S", "--step", type=int, help="The step size that denoising process applied by, the default step size is 50bps.", default=50)
    parser.add_argument("-r", "--region", help="A single section of input file could be defined as chromosomename:startindex:endindex")
    parser.add_argument("-seg", "--segmentation", action='store_true', help="Run segmentation")
    parser.add_argument("-l", "--interval", help="A list of regions could be defined as a BED file format")
    parser.add_argument("-H", "--highresolution", action='store_true', help="Force the program to define step and span equal to 1 to have the highest resolution as possible (One base resolution)")
    parser.add_argument("-s", "--windowsize", type=int, help="Define the window size of applying filter", required=False, default="65536")
    parser.add_argument("-V", "--version", action='version', version="%(prog)s (Version 1.0)")
    parser.add_argument("-e", "--exclude", help="Define the name and location of the BigWig file as BlackList_file. e.g. -i data/File.bw ")
    # parser.add_argument("-scale", "--scalingfactor", type=int, help="Define the size of the scaling factor which computes by the software by default. e.g. -scale 247 ")
    parser.add_argument("-t", "--tempfolder", help="The address of the temporary folder.", default="temp/")
    parser.add_argument("-low", "--lowresolution", help="Export low resolution denoised file in defined location")
    parser.add_argument("-f", "--filter", help="Insert the filter name of Scipy Signal filter list. Blackman filter selected by default", required=False, default="blackman")
    parser.add_argument("-list", "--filterlist", action=_ListAction, help="Print all the possible signal filters")
    parser.add_argument("-fold", "--foldtimes", type=int, help="Number of folds of low resolution window size in compare to high resolution denoising window size.", required=False, default="4")

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
    SaveLowReso = args.lowresolution
    foldtimes = args.foldtimes
    if args.peak:
        s_output = args.peak
        prime_output = 'peak_low_resolution_'+args.peak
    else:
        s_output = None


    if TEMPFolder[-1] != '/':
        TEMPFolder=str(TEMPFolder)+'/'
    try:
        os.makedirs(TEMPFolder)
    except:
        # print('Failed to make Dir and passed\n')
        pass
    # -----------input/step validation----------------------
    d_span=Step_size
    # If the high resolution option selected by user, it is useful for the short peaks that needs more resulations. it means all the bases of defined region will be called and processed.
    if args.highresolution:
        Step_size, d_span = 1, 1

    # find log base 2 of user defined window size and sum it by the foldtime which represents the windows size of the low resolution smoothed signal
    primewindowsize = 2 << (int(numpy.log2(userdefined_windowsize))+foldtimes-1)
    PWS = primewindowsize >> int(numpy.log2(Step_size))
    # Window size of the signal resized regaurd the input file
    WSize = userdefined_windowsize >> int(numpy.log2(Step_size))

    # Location of the input file stores on "my_file" variable
    my_file = pathlib.Path(source_data)
    current_folder = (os.getcwd())
    # print(current_folder)
    # -----------Excluder / Blacklist-----------------------

    # input file validates and import if exist
    try:
        Open_file = pyBigWig.open(source_data)
        chrlist=list(Open_file.chroms())
    # This module could be applied on denoised signals before segmentation. Denoised file will load if exist
    except:
        sys.stderr.write('\n# Excluder error: Input file not defined.\n')
        sys.exit()


    if BlackList_file:
        try:
            # Import blacklist file, and merge the intervals if they have overlap
            a = pybedtools.BedTool(BlackList_file)
            BL = a.sort().merge()
        except:
            sys.stderr.write('\n# Excluder error: Blacklist file not found.\n')
            sys.exit()
    else:
        BL =[]

    sys.stdout.write('\r Please wait ...')
    if Region_of_interest:
        target='Region'
        Region=Region_of_interest
    elif user_defined_interval:
        target='Interval'
        Region=user_defined_interval
    else:
        target='Entrie_genome'
        Region=''
    fl.intervalconverter(BL, target, chrlist, Open_file, TEMPFolder, WSize, Step_size, Region)
    fl.Excluder(BL, (TEMPFolder+'intervals.bed'), TEMPFolder, Open_file, Step_size)

    sys.stdout.write('\r Exclusion completed.\n')# in',"%.5f seconds" % (time.time()-s_time))

    # -----------End of excluder / blacklist----------------

    # Check whether the path/file presence or not
    if my_file.is_file():
        pass
    else:
        sys.stderr.write('\n# Filter error: File \"'+ my_file +'\" is not exist. Check the file name and relative path.\n')
        sys.exit()


    # --------------filterchecker------------------------
    # Store filter signal curve with the specified size within an array
    try:
        d_signal=eval('scipy.signal.%s(%s)' % (d_filter,WSize)) # Cross the scaling factor
        low_resolution_signal=eval('scipy.signal.%s(%s)' % (d_filter,PWS))
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
    # Call the input bigwig file by pyBigWig library function and store in a array
    d_open = pyBigWig.open(source_data)
    # check if the output file cannot be created because of the output folder does not exist or permissions, an error message with clean exit will occur
    try:
        d_output = pyBigWig.open(denoised_data, "w")
    except IOError:
        sys.stderr.write("\n# Filter error:  The PATH/FILE `"+ denoised_data +"` has not write permission, please check the path and try again. \n - \n")
        sys.exit()
    if SaveLowReso:
        try:
            low_resolution_output = pyBigWig.open(SaveLowReso, "w")
        except IOError:
            sys.stderr.write("\n# Filter error:  The PATH/FILE `"+ SaveLowReso +"` has not write permission, please check the path and try again. \n - \n")
            sys.exit()
    # If file created successfully, set the header of output file similar as input
    else:
        name=denoised_data.strip().split('/')[-1]
        low_resolution_output = pyBigWig.open(str(TEMPFolder)+"lowresolution_"+name, "w")

    d_output.addHeader(list(d_open.chroms().items()))
    low_resolution_output.addHeader(list(d_open.chroms().items()))

    # --------------names.txt------------------------
    # Define the temp list by all the chromosome names that present within the input file (BigWig)
    chrlist=list(d_open.chroms())
    Header=list(d_open.chroms().items())
    if primewindowsize > int(Header[2][1]):
        sys.stderr.write("\n# Filter error:  Window Size for low_resolution denoising is greater than chromosme size. \n - \n")
        sys.exit()

    my_interval = pybedtools.BedTool(str(TEMPFolder)+'intervals.bed')
    for line in my_interval:
        d_openvalue = d_open.values(line.chrom, line.start, line.end, numpy=True)
        d_averages=fl.Average(d_openvalue, Step_size)
        if len(d_openvalue) < Step_size*20:
            d_warning=True
        try:
            # Combine signals derived from the input file for each step's value by filter with specific windows size
            d_convolve = scipy.signal.fftconvolve(d_averages, d_signal, mode="same")
            # print(d_convolve) #<<<< remove me
            # This process may takes time, so one message informs the process stage
            sys.stdout.flush()
            sys.stdout.write('\r '+line.chrom+'                     ')
            # Update the standard output
            sys.stdout.flush()
            # Write the combined values for each interval to the output file
            d_output.addEntries(line.chrom, line.start, values=d_convolve, span=d_span, step=Step_size)
        except:
            sys.stderr.write("\n# Filter error: Please check the Bed file \"" + user_defined_interval + "\" .\n")
            sys.exit()
        low_resolution_convolve = scipy.signal.fftconvolve(d_averages, low_resolution_signal, mode="same")
        low_resolution_output.addEntries(line.chrom, line.start, values=low_resolution_convolve, span=d_span, step=Step_size)
    sys.stdout.write('\r Saving bigWig file...')

    d_output.close()
    low_resolution_output.close()
    # --------------  WARNING  ------------
    # If the ratio of region/interval size and step size is not ideal, a warning message prints for the user
    # if d_warning==True:
    #     sys.stderr.write('\n# Warning: The defined Span=%s is too large. It could cause an error in the results because of the number of samples. It is recommended to use the smaller step.\n - ' % d_span)
    sys.stdout.write('\r Denoising completed. \n')# in',"%.5f seconds" % (time.time()-s_time))

    ######################End of filter

    saveme=[]
    # ------------ load data ---------------------

    if SaveLowReso:
        prime_input = SaveLowReso
    else:
        name=denoised_data.strip().split('/')[-1]
        prime_input = str(TEMPFolder)+"lowresolution_"+name
    try:
        sg_input = pyBigWig.open(denoised_data)
    except:
        print(denoised_data+' is not found.')
        exit()
    try:
        prime_sg_input = pyBigWig.open(prime_input)
    except:
        print(prime_input+' is not found.')
        exit()
    if s_output is not None:
        sg_output = pyBigWig.open(s_output, "w")
        sg_output.addHeader(list(sg_input.chroms().items()))
        chrlist=list(sg_input.chroms())
        prime_sg_output = pyBigWig.open(prime_output, "w")
        prime_sg_output.addHeader(list(sg_input.chroms().items()))

    # ---------------- import BED file --------------
    try:
        data = pybedtools.BedTool(str(TEMPFolder)+'intervals.bed')
    except:
        sys.stderr.write("File \"" + str(TEMPFolder)+'intervals.bed' + "\" is not exist. Check the BED file and it\'s path. \n")
        exit()
    # ---------------- segmentaion ------------------
    # For each interval the value reads and stored in a varible
    for line in data:
        sys.stdout.write('\r '+line.chrom+'       ')
        Sub_regions=sg.high_resolution(line.chrom, line.start, line.end, Step_size, sg_input, s_output, prime_sg_input)
        sg.low_resolution(line.chrom, line.start, line.end, Step_size, Sub_regions, prime_sg_input, TEMPFolder, s_output, saveme)
    saveto = os.path.dirname(denoised_data)
    pybedtools.BedTool(saveme).saveas(saveto+'/Boundaries.bed')

    if s_output:
        sg_input.close()
        sg_output.close()
        prime_sg_input.close()
        prime_sg_output.close()

sys.stdout.write("\r Execution times %.5f seconds\n" % (time.time()-dsp_time))
