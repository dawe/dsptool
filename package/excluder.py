#!/usr/bin/env python
# -----------functions Definition-----------------------------------------------------------------------

# ---------------Ranges--------------
def ranges(nonzerolist, WSize, init, fin):
    """Range function identifies the connectivity of the indices of the zero-value points within the input file and return them
    as a paired-coordinate list if they cover more than 2500 base pairs.

    Keyword arguments:
    nonzerolist: list of the indices of the nonzero-value points within the 1D-matrix

    Function description:
    The input list is sorted, so range function identifies all the gaps between the elements (equal to zero-value indices) and
    if the length of these gaps are more than window size (131072) they will be reported as a list with two elements
    start and end positions.
    """
    import numpy as np
    gaps = [[s, e] for s, e in zip(nonzerolist, nonzerolist[1:]) if s+1 < e]
    zero = list([[i[0]+1,i[1]+1] for i in gaps if i[1]-i[0]>WSize])
    if init and len(zero) > 0:
        zero = np.insert(zero, 0, init, axis=0)
    elif init and len(zero) == 0:
        zero = [init]
    if fin and len(zero) > 1:
        zero = np.insert(zero, -1, fin, axis=0)
    return (zero)

# --------------size of window correction-------------

# The function increase the accuracy of sampling when Step is higher than 1. Two parameters are used in this fuction, VALUE is a list of intensities that loaded from input file related to the specific regions and STEP is the distance between each sample through out this regions.
def AverageofRegion(VALUE, Step_size):
    import numpy as np
    if Step_size==1:
        return(VALUE)
    else:                                                                                       ####
        # Perform average calculation using a 2-dimensional matrix, empty cells filled by zero     #
        d_matrix = [VALUE]                                                                         #
        # Define the number of the cells must be added on the tail of the 1D-matrix                #
        zero = Step_size - (len(VALUE) % Step_size)                                                #
        # Find the number of required cells to fill the matrix by AVERAGE function                 #
        segnum = len(VALUE) // Step_size                                                           #
        if zero > 0:                                                                               #  Average of the input signal
            segnum += 1                                                                            #
        # fill the free cells of the matrix with 0                                                 #
        d_matrix = np.pad(d_matrix,[(0,0),(0,zero)], mode='constant', constant_values=0)           #
        # reshape the matrix from the one dimensional to step size rows and multiple columns       #
        d_matrix = np.reshape(d_matrix,(segnum, Step_size))                                        #
        # Average of each raw reported as the intensity of that step                               #
        Mat=d_matrix.mean(axis=1)                                                               ####
        return(Mat)

# ---------------Interval coverter--------------

def intervalconverter(BL, target, HEADER, Step_size, Region):
    """ The interval converter function converts the region of interest or defined interval to an input.bed
    file, store within the temp folder.

    Keyword arguments:
    Step_size -- An integer number defined by the user, the size of the data reading from the bigWig file

    Function description:
    User could define one coordinate for a region by using "-r" switch, or list of the regions by define
    a bed file by using "-l" switch, or even entire input file without any switch definition.
    Function divides the standard coordinate to relative information and if they are valid one BED file
    stores the information. In the same way interval file after import validates and stores in the temp folder.
    For the entire genome, all the chromosome names present in the header of the BigWig file with relative
    chromosome size stores in the temporary file "input.bed".
    """
    import sys, pybedtools
    save = []
    # If switch "r" or "Region" defined by the user, the inserted value must follow the standard pattern that contains chromosome name, start and end coordinates
    if target=='Region' :
        try:
            # Characters before column-sign should be a chromosome name
            _splited = Region.strip().split(':')
            chrname = _splited[0]
            # The phrase after column-sign consists of two coordinates that separated by a dash-sign
            _splited = _splited[1].strip().split('-')
            start = int(_splited[0])
            end = int(_splited[1])
        except:
            sys.stderr.write('\n# Excluder error: The region of interest must follow the standard pattern ChromosomeName:StartPosition-EndPosition.\n')
            exit()
        if chrname in HEADER.keys():
            # If the region of interest present in the input file, the program goes ahead
            if (int(HEADER[chrname]) < (end)) :
                end = int(HEADER[chrname])
                print('\n# Excluder: The end position is higher than the chromosome length.\n            The end coordinate corrected to ',HEADER[chrname])
            if (int(HEADER[chrname]) < (start)) :
                sys.stderr.write('\n# Excluder error: The start position is wrong.\n')
                exit()
        else:
            sys.stderr.write('\n# Excluder error: The chromosome \"'+ chrname +'\" is not present in the input file.\n')
            exit()
        INPUT = pybedtools.BedTool([(chrname, start, end)])

    elif target=='Interval':
        import pathlib, pybedtools
        bed_file = pathlib.Path(Region)
        if bed_file.is_file():
            try:
                INPUT = pybedtools.BedTool(Region).sort().merge(d=Step_size-1)
            except:
                sys.stderr.write('\n# Excluder error: The interval file must have 3 columns, consist of chromosome Name, start position and end position.\n')
                exit()
        else:
            sys.stderr.write('\n# Excluder error: The interval BED file could not found. Please check the interval file name/location and try again.\n')
            exit()

    else:

        INPUT = pybedtools.BedTool([(c, 1, HEADER[c]) for c in HEADER.keys()])
    return(INPUT)

# -----------End of functions Definition-----------------------------------------------------------------------
