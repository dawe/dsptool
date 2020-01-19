# !/usr/bin/env python3.7
from scipy import ndimage
from UliEngineering.SignalProcessing.Utils import zero_crossings
import argparse, sys, pybedtools, gc
from pybedtools import BedTool
from skimage import data, filters
import scipy.signal, time, pybedtools

#----------------------------------------------
def boundaries(intensity, LoCoMax, LoCoMin, chr, start, end, Step_size=50):
    """boundaries function uses the list_of_maxima_coordinates, list_of_minima_coordinates and inflection points (second derivative
    of intensity) and returns a matrix of boundaries related to each Maximum point.

    Keyword arguments:
    intensity -- A list of input value regarding the defined interval
    LoCoMax -- A coordinates list of all the Maximum that found in the called input corresponding to the relative start coordinates 
    LoCoMin -- A coordinates list of all the Minimum that found in the called input corresponding to the relative start coordinates 
    chr -- chromosome name of the input values
    start -- start coordinate of the input interval corresponding to the start of the chromosome
    end -- end coordinate of the input interval corresponding to the start of the chromosome
    Step_size -- The step size of the reading the data from the bigWig file (default 50)

    Function description:
    The nearest inflection points to the minimum point of the curve in two sides of a maximum point report in a matrix with 3 columns.
    The first column, chromosome name, second column start coordinate of the boundary that relates to the Maximum, third column end
    coordinate of the boundary that relates to the Maximum.
    If at least 2 indeces are present in the input data, for each maxima within it, one boundary will be report. If there is any minimum
    for the input sample, and the maximum coordinate is between two minimum,

    """
    import numpy

    eco=0
    if len(intensity) >= 2 :
        dxU=numpy.gradient(intensity)
        dxxU = numpy.gradient(dxU)
        LoCo_infl = numpy.where(numpy.diff(numpy.sign(dxxU)))[0]
        matrix = [[0]*3]
        sco = start
        for i in LoCoMax:
            if len(LoCoMin) > 0: # any minimum
                insertion_index=numpy.searchsorted(LoCoMin, i) # if 0 or length of list, it means out of the list
                if insertion_index == 0:    # if there is no minimum before the maximum point
                    sco = LoCo_infl[0]*Step_size+start
                    Rmin = LoCoMin[0]
                    Rinflection_index=numpy.searchsorted(LoCo_infl, Rmin)
                    eco=LoCo_infl[Rinflection_index-1]*Step_size+start
                elif len(LoCoMin) == insertion_index: # if there is no minimum after the maximum point
                    Lmin = LoCoMin[-1]
                    Linflection_index=numpy.searchsorted(LoCo_infl, Lmin)
                    sco=LoCo_infl[Linflection_index]*Step_size+start
                    eco = LoCo_infl[-1]*Step_size+start
                elif insertion_index < len(LoCoMin):
                    Lmin = LoCoMin[insertion_index-1]
                    Rmin = LoCoMin[insertion_index]
                    Linflection_index=numpy.searchsorted(LoCo_infl, Lmin)
                    sco=LoCo_infl[Linflection_index]*Step_size+start
                    Rinflection_index=numpy.searchsorted(LoCo_infl, Rmin)
                    eco=LoCo_infl[Rinflection_index-1]*Step_size+start

            elif len(LoCoMin) == 0 and len(LoCoMax) >= 1:
                if len(LoCo_infl) > 0:
                    sco = LoCo_infl[0]*Step_size+start
                    eco = LoCo_infl[-1]*Step_size+start

            if int(matrix[-1][1])!=int(sco) and int(sco) < int(eco):
                 matrix = numpy.vstack([matrix, (chr,sco,eco)])
        matrix = numpy.delete(matrix, (0), axis=0)
        return(matrix)

#----------------------------------------------
def write2BED(Sub_regions, Regions, TEMPFolder):
    """Take two matrices as input and create two bed12 files.

    Keyword arguments:
    blocks -- A matrix contain boundaries of the high_resolution segmentaion im four columns created by boundaries function
    edges -- A matrix contain boundaries of the low_resolution segmentaion im four columns created by boundaries function
    Step_size -- The step size of the reading the data from the bigWig file (default 50)
    start -- start coordinate of the input interval corresponding to the start of the chromosome

    Error notifications:
    if the start coorindate is higher or equal to the end cooridnate, program will stop the program and warn it.

    Function description:
    Each bed12 file is containg chromosome name, The starting and end positions of the feature in the chromosome, name of the
    feature equal to the starting position, and other informations
    """
    save1, save3 = [], []
    for i in Sub_regions:
        a, b, c=i[0], int(i[1]), int(i[2])
        # in the case of the problem, start must have lower value in compare to the end. If wrong data found, program will warn and close.
        if c-b <= 0:
            print('During the bed file writing, one error occurs. The start position ',b,'is higher than end position',c,'at chromosome',a)
            exit()
        string=("%s\t%s\t%s\tH") % (a,b,c)
        save1.append((string))
        save3.append((string))
    pybedtools.BedTool(save1).saveas(str(TEMPFolder)+'blocks.bed')
    H = '%sblocks.bed' % (TEMPFolder)
    save2=[]
    for i in Regions:
        a, b, c=i[0], int(i[1]), int(i[2])
        # in the case of the problem, start must have lower value in compare to the end. If wrong data found, program will warn and close.
        if c-b < 0:
            print('During the bed file writing, one error occurs. The start position ',b,'is higher than end position',c,'at chromosome',a)
            exit()
        string=("%s\t%s\t%s\tL") % (a,b,c)
        save2.append((string))      #could be deleted
        save3.append((string))
    # pybedtools.BedTool(save2).saveas('%slow_resolution_edges_%s.bed' % (str(TEMPFolder),denoised_data)) # could be deleted
    # L = '%slow_resolution_edges_%s.bed' % (str(TEMPFolder),denoised_data)
    pybedtools.BedTool(save3).sort().merge().saveas(str(TEMPFolder)+'sum.bed') # could be deleted
    return(H)
# ------------- peak finding at the smaller window size ------------------
def high_resolution(chr, start, end, Step_size, sg_input, s_output, prime_sg_input):
    """high_resolution function reads the data from the input file (with defined window_size) for a specific region, find all the maxima and minima and relative boundary
    for each maximum, and returns maxima by creating a bigWig file.

    Keyword arguments:
    chr -- chromosome name of the input values
    start -- start coordinate of the input interval corresponding to the start of the chromosome
    end -- end coordinate of the input interval corresponding to the start of the chromosome
    Step_size -- The step size of the reading the data from the bigWig file (default 50)

    Function description:
    For each interval, scipy.signal.find_peaks identifies all the extremum of the input and returns the Edgeblocks for each identified maximum boundary.
    """
    # print(' - Peak finding chromosome',chr,'from',start,'to',end)
    sg_value = sg_input.values(chr, start, end, numpy=True)[::Step_size]
    # --------------- Find the maxima of the signal ----------------------
    Maxs=(scipy.signal.find_peaks(sg_value)[0])
    # --------------- invert the signal in order to detect the minima ----
    inv_value=sg_value*(-1)
    Mins=(scipy.signal.find_peaks(inv_value)[0])
    # --------------- Identify the boundary of each peak and report it as a matrix --------
    Sub_regions = boundaries(sg_value, Maxs, Mins, chr, start, end, Step_size)
    # --------------- Export the peaks with relative intensity as a bigWig file -----------
    if s_output:
        j,k=0,0
        for i in Maxs:
            k = i*Step_size+start
            intensity =float(sg_value[i])
            j = i
            sg_output.addEntries(chr, k, values=[intensity], Step_size=Step_size, step=Step_size)
    return(Sub_regions)

# ------------- peak finding at the larger window size ------------------
def low_resolution(chr,start,end,Step_size,Sub_regions, prime_sg_input, TEMPFolder, s_output, saveme):
    """low_resolution function reads the data from the input file (with higher window_size) for a specific region, find all the maxima and minima and relative boundary
    for each maximum, run the write2BED function for the boundaries in both low_resolution and high_resolution function and returns maxima by creating a bigWig file.

    Keyword arguments:
    chr -- chromosome name of the input values
    start -- start coordinate of the input interval corresponding to the start of the chromosome
    end -- end coordinate of the input interval corresponding to the start of the chromosome
    Step_size -- The step size of the reading the data from the bigWig file (default 50)

    Function description:
    For each interval, scipy.signal.find_peaks identifies all the extremum of the input and returns the Edgeblocks for each identified maximum boundary. The boundaries 
    identified using the boundaries function send to the write2BED function to create a Bed12 file.
    """
    prime_sg_value = prime_sg_input.values(chr, start, end, numpy=True)[::Step_size]
    # --------------- Find the maxima of the signal ----------------------
    PMaxs=(scipy.signal.find_peaks(prime_sg_value)[0])
    # --------------- invert the signal in order to detect the minima ----
    inv_value=prime_sg_value*(-1)
    PMins=(scipy.signal.find_peaks(inv_value)[0])
    # --------------- Identify the boundary of each peak and report it as a matrix --------
    Regions = boundaries(prime_sg_value, PMaxs, PMins, chr, start, end, Step_size)
    # -------------- Write the boundary matrixes as bed file ------------------------------
    if Sub_regions is not None:
        H = write2BED(Sub_regions,Regions, TEMPFolder)
        intersect=pybedtools.BedTool(str(TEMPFolder)+'sum.bed').intersect(H, C=True)
        sub_regions = pybedtools.BedTool(H)
        regions = pybedtools.BedTool(str(TEMPFolder)+'sum.bed')
        counts = [int(x[3]) for x in intersect if int(x[3]) != 0]
        blocksize, starts, c= '','',0
        for n in range (len(counts)):
            blocksize, starts = '',''
            for x in range (c,counts[n]+c):
                # print('************\n',starts,(sub_regions[x].start)-(regions[n].start))
                starts=("%s,%s") % (starts,(sub_regions[x].start)-(regions[n].start))
                blocksize=("%s,%s") % (blocksize,(sub_regions[x].end-sub_regions[x].start))
            starts=starts.lstrip(',')
            blocksize=blocksize.lstrip(',')
            try:
                string=("%s\t%s\t%s\t%s\t1000\t.\t%s\t%s\t255,50,50\t%d\t%s\t%s") % (regions[n].chrom,regions[n].start,regions[n].end,regions[n].start,sub_regions[n].start,sub_regions[counts[n]].end,counts[n],blocksize,starts)
            except:
                string=("%s\t%s\t%s\t%s\t1000\t.\t%s\t%s\t255,50,50\t%d\t%s\t%s") % (regions[n].chrom,regions[n].start,regions[n].end,regions[n].start,regions[n].start,regions[n].end,counts[n],blocksize,starts)
            c=c+counts[n]
            saveme.append((string))
            # print(string) #<<<<<< Remove me
    # --------------- Export the peaks with relative intensity as a bigWig file -----------
    if s_output:
        j,k=0,0
        for i in PMaxs:
            k = i*Step_size+start
            intensity =float(prime_sg_value[i])
            j = i
            prime_sg_output.addEntries(chr, k, values=[intensity], span=Step_size, step=Step_size)
    return()


# -------------- Write the boundary matrixes as bed file ------------------------------
# pybedtools.BedTool(bl).saveas('edges(%s).bed' % denoised_data)
# pybedtools.BedTool(eg).saveas('prime_edges(%s).bed' % denoised_data)
# ---------------- close the files --------------
