#!/usr/bin/env python
# -----------functions Definition-----------------------------------------------

# ---------------Ranges---------------------------------------------------------
def ranges(nonzerolist, WSize, init, fin):
    """
Range function identifies the connectivity of the indices of the zero-value
points within the input file and returns them as a paired-coordinate list if
they cover more than the window size.

Arguments:
nonzerolist: The list of indices of the nonzero-value bases.
WSize: Window size in base pairs unit.
init: A variable that carries a boolean value indicates the interval starts with
      a zero-valued base or not.
fin:  A variable that carries a boolean value indicates the interval ends with a
      zero-valued base or not.

Function description:
The input list is sorted, so the range function identifies all the gaps between
the elements (the value of the gap is equal to zero). If the length of these
gaps is more than window size (default 16384 bp) they will be reported as a list.
The list "Zetta" has two elements of start and end coordinates.
    """
    import numpy as np
    index_gaps = np.where(np.diff(nonzerolist) >= WSize)
    Zetta = list([[nonzerolist[i],nonzerolist[i+1]] for i in index_gaps[0]])
    if init and len(Zetta) > 0:
        Zetta = np.insert(Zetta, 0, init, axis=0)
    elif init and len(Zetta) == 0:
        Zetta = [init]
    if fin:
        Zetta = np.insert(Zetta, len(Zetta), fin, axis=0)
    return (Zetta)

# --------------size of window correction-------------

def AverageofRegion(VALUE, Step_size):
    """
When the step size is higher than 1 bp, downsampling with one basepair for each
step-size length is required. As the input signal is noisy and the value of
successive bases is fluctuating, this function increases the accuracy of
downsampling by taking the average of each step and report it as the
representative of that step, instead of taking one sample in each step.

Arguments:
VALUE: A list of intensities that loaded from the input file in the inspecting
       regions.
Step_size: The distance between two consecutive samples in the regions under
           investigation.

Function description:
If the step size is equal to 1, it means there is no need to downsampling of the
input signal. The input array converts to a two-dimensional matrix. The size of
columns in this matrix is equal to the step size. The empty cells of the matrix
at the end of the array fill with zero numbers. An array of the average of each
row provides the sample values of the input signal.
    """

    import numpy as np
    if Step_size==1:
        return(VALUE)
    else:

        d_matrix = [VALUE]
        # Define the number of the cells should be added on the tail of the 1D-matrix to fit a matrix.
        zero = Step_size - (len(VALUE) % Step_size)

        # Find the number of rows of the matrix.
        rownum = len(VALUE) // Step_size
        if zero > 0:
            rownum += 1

        # Fill the free cells of the matrix with 0 value.
        d_matrix = np.pad(d_matrix,[(0,0),(0,zero)], mode='constant', constant_values=0)

        # reshape the matrix from the one dimensional array to a matrix.
        d_matrix = np.reshape(d_matrix,(rownum, Step_size))

        # An array of rows' means of the matrix, provides the sample values of the input signal.
        Mat=d_matrix.mean(axis=1)
        return(Mat)

# ---------------Interval coverter--------------

def intervalconverter(DOMAIN, HEADER, Step_size, REGION):
    """
This function uniforms the different types of intervals requested by the user.
The user could define one coordinate for a region by using the "-r" switch, or a
list of the regions by defining a bed file by using the "-l" switch, or even an
entire input file without any switch definition. Function divides the standard
coordinate into relative information and if they are valid one BED file stores
the information. In the same way interval file after import validates and stores
in the temp folder. For the entire genome, all the chromosome names present in
the header of the BigWig file with relative chromosome size.

Arguments:
DOMAIN: A string contains the type of the input interval ("Region", "Interval",
        or "Entrie_genome").
HEADER: The header of the input bigWig file contains the chromosome names and
        their lengths.
Step_size: The distance between two consecutive samples in the regions under
           investigation.
REGION: A string containing chromosome name, start and end coordinates.

Function description:
If the target domain for the inspecting is a region of one single chromosome,
first program checks that the text pattern is correct. Then splits the domain to
the chromome name, start and the end coordinates. According the header of the
input file, if the region of interest presents in the input file, the interval
saves in an array with the format of a bed file.
    """

    import sys, pybedtools
    save = []
    # If switch "r" or "Region" defined by the user, the inserted value must follow
    # the standard pattern that contains chromosome name, start and end coordinates
    if DOMAIN=='Region' :
        try:
            # Characters before column-sign should be a chromosome name
            _splited = REGION.strip().split(':')
            chrname = _splited[0]

            # The phrase after column-sign consists of two coordinates that separated by a dash-sign
            _splited = _splited[1].strip().split('-')
            start = int(_splited[0])
            end = int(_splited[1])
        except:
            sys.stderr.write('\n* error: The region of interest must follow the standard pattern ChromosomeName:StartPosition-EndPosition.\n')
            exit()

        if chrname in HEADER.keys():
            # If the region of interest present in the input file, the program splits the input region
            # to chromosome name, start and end coordinate.
            if (int(HEADER[chrname]) < (end)) :
                end = int(HEADER[chrname])
                print('\n* Excluder: Defined end cordination is out of the chromosome length.\n\033[1;36m            The end coordinate corrected to \033[0;0m',HEADER[chrname])
            if (int(HEADER[chrname]) < (start)) :
                sys.stderr.write('\n* error: The start position is wrong.\n')
                exit()
        else:
            sys.stderr.write('\n* error: The chromosome \"'+ chrname +'\" is not present in the input file.\n')
            exit()
        INPUT = pybedtools.BedTool([(chrname, start, end)])

    # If the domain is a bedfile, program loads, sorts, and merges the ".bed" file using the pybedtools,
    elif DOMAIN=='Interval':
        import pathlib, pybedtools
        bed_file = pathlib.Path(REGION)
        if bed_file.is_file():
            try:
                INPUT = pybedtools.BedTool(REGION).sort().merge(d=Step_size-1)
            except:
                sys.stderr.write('\n* error: The interval file must have 3 columns, consist of chromosome Name, start position and end position.\n')
                exit()
        else:
            sys.stderr.write('\n* error: The interval BED file could not found. Please check the interval file name/location and try again.\n')
            exit()

    # If the entire of the input file is requested for inspection, all the chromosome names with their
    # length from the header of the input signal extracted and stored in a ".bed" file format.
    else:

        INPUT = pybedtools.BedTool([(c, 1, HEADER[c]) for c in HEADER.keys()])
    return(INPUT)

# -----------End of functions Definition-----------------------------------------------------------------------
