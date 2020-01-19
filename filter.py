# !/usr/bin/env python3.7

# -----------function Definitions-----------------------------------------------------------------------

def ranges(nonzerolist):
    """Range function identifies the connectivity of the indices of the zero-value points within the input file and return them
    as a paired-coordinate list if they cover more than 2500 base pairs.

    Keyword arguments:
    nonzerolist -- list of the indices of the nonzero-value points within the 1D-matrix

    Function description:
    The input list is sorted, so range function identifies all the gaps between the elements (equal to zero-value indices) and
    if the length of these gaps are more than 50 times of step size (50 x 50) they will be reported as a list with two elements
    start and end positions.
    """
    gaps = [[s, e] for s, e in zip(nonzerolist, nonzerolist[1:]) if s+1 < e]
    return (list([i for i in gaps if i[1]-i[0]>50]))


# The function increase the accuracy of sampling when Step is higher than 1. Two parameters are used in this fuction, VALUE is a list of intensities that loaded from input file related to the specific regions and STEP is the distance between each sample through out this regions.
def Average(VALUE, STEP):
    import numpy
    if STEP==1:
        return(VALUE)
    else:
        # Perform average calculation using a 2-dimensional matrix, empty cells filled by zero
        d_matrix = [VALUE]
        # Define the number of the cells must be added on the tail of the 1D-matrix
        zero = STEP - (len(VALUE) % STEP)
        # Find the number of required cells to fill the matrix by AVERAGE function
        segnum = len(VALUE) // STEP
        if zero > 0:
            segnum += 1
        # fill the free cells of the matrix with 0
        d_matrix = numpy.pad(d_matrix,[(0,0),(0,zero)], mode='constant', constant_values=0)
        # reshape the matrix from the one dimensional to step size rows and multiple columns
        d_matrix = numpy.reshape(d_matrix,(segnum, STEP))
        # Average of each raw reported as the intensity of that step
        Mat=d_matrix.mean(axis=1)
        return(Mat)

# ---------------Average and Ranges--------------
def ranges(nonzerolist, z):
    """Range function identifies the connectivity of the indices of the zero-value points within the input file and return them
    as a paired-coordinate list if they cover more than 2500 base pairs.

    Keyword arguments:
    nonzerolist: list of the indices of the nonzero-value points within the 1D-matrix

    Function description:
    The input list is sorted, so range function identifies all the gaps between the elements (equal to zero-value indices) and
    if the length of these gaps are more than 50 times of step size (50 x 50) they will be reported as a list with two elements
    start and end positions.
    """
    gaps = [[s, e] for s, e in zip(nonzerolist, nonzerolist[1:]) if s+1 < e]
    return (list([[i[0]+1,i[1]+1] for i in gaps if i[1]-i[0]>z]))

# --------------size of window correction-------------

# The function increase the accuracy of sampling when Step is higher than 1. Two parameters are used in this fuction, VALUE is a list of intensities that loaded from input file related to the specific regions and STEP is the distance between each sample through out this regions.
def AverageofRegion(VALUE, window_size, Step_size):
    import numpy
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
        d_matrix = numpy.pad(d_matrix,[(0,0),(0,zero)], mode='constant', constant_values=0)        #
        # reshape the matrix from the one dimensional to step size rows and multiple columns       #
        d_matrix = numpy.reshape(d_matrix,(segnum, Step_size))                                     #
        # Average of each raw reported as the intensity of that step                               #
        Mat=d_matrix.mean(axis=1)                                                               ####

        # Find the index list of non-zero values within the input signal
        indices_nonzeros = numpy.nonzero(Mat != 0)[0]
        # Find stretches of zero by at least window_size length
        if len(indices_nonzeros) == 0:
            indices_nonzeros=[0,(len(Mat)-1)]
        elif indices_nonzeros[0] != 0 and indices_nonzeros[-1] != 0:
            indices_nonzeros = numpy.insert(indices_nonzeros, 0, 0)
            indices_nonzeros = numpy.append(indices_nonzeros, (len(Mat)-1))
        elif indices_nonzeros[0] != 0:
            indices_nonzeros = numpy.insert(indices_nonzeros, 0, 0)
        elif indices_nonzeros[-1] != (len(Mat)):
            indices_nonzeros = numpy.append(indices_nonzeros, (len(Mat)-1))
        zeros = ranges(indices_nonzeros, window_size)
        # print(zeros) #<<<< remove me
        return(zeros)

# ---------------Interval coverter--------------
def intervalconverter(BL, target, chrlist, Open_file, TEMPFolder, WSize, Step_size, Region):
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
    import pybedtools, numpy, sys
    save = []
    # If switch "r" or "Region" defined by the user, the inserted value must follow the standard pattern that contains chromosome name, start and end coordinates
    if target=='Region' :
        try:
            # Characters before column-sign should be a chromosome name
            _splited = Region.strip().split(':')
            chrname = _splited[0]
            # The phrase after column-sign consists of two coordinates that separated by a dash-sign
            _splited = _splited[1].strip().split('-')
            startprime = int(_splited[0])
            endprime = int(_splited[1])
        except:
            sys.stderr.write('\n# Excluder error: The region of interest must follow the standard pattern ChromosomeName:StartPosition-EndPosition.\n')
            exit()
        if chrname in chrlist:
            # If the region of interest present in the input file, the program goes ahead
            if (int(Open_file.chroms(chrname)) < (endprime)) :
                endprime = int(Open_file.chroms(chrname))
                print('\n# Excluder: The end position is higher than the chromosome length.\n            The end coordinate corrected to ',Open_file.chroms(chrname))
            if (int(Open_file.chroms(chrname)) < (startprime)) :
                sys.stderr.write('\n# Excluder error: The start position is wrong.\n')
                exit()
        else:
            sys.stderr.write('\n# Excluder error: The chromosome \"'+ chrname +'\" is not present in the input file.\n')
            exit()
        string=("%s\t%s\t%s") % (chrname,startprime,endprime)
        save.append(string)

        pybedtools.BedTool(save).saveas(TEMPFolder+'input.bed')

    elif target=='Interval':
        import pathlib
        bed_file = pathlib.Path(Region)
        if bed_file.is_file():
            bed_file = pybedtools.BedTool(Region)
            # define a temporary file and put the list of the chromosome name inside it, because the Bedtools does not accept the list as a array
            temporaryfile = open(TEMPFolder+'genome.txt',"w")
            for i in Open_file.chroms().items():
                temporaryfile.write(i[0]+"\t"+str(i[1])+"\n")
            temporaryfile.close()
            tmp = open(TEMPFolder+"genome.txt")
            # Check if the Bed file present, continue
            try:
                # for bedtools sort Version: v2.26.0
                sites = (bed_file.sort(faidx = tmp.name).merge(d=Step_size-1))
            except:
                sys.stderr.write('\n# Excluder error: The interval file must have 3 columns, consist of chromosome Name, start position and end position.\n')
                exit()
            pybedtools.BedTool(sites).saveas(TEMPFolder+'input.bed')

        else:
            sys.stderr.write('\n# Excluder error: The interval BED file could not found. Please check the interval file name/location and try again.\n')
            exit()

    else:
        for CHR in chrlist:
            string=("%s\t%s\t%s") % (CHR, 1, Open_file.chroms(CHR))
            save.append(string)

        pybedtools.BedTool(save).saveas(TEMPFolder+'input.bed')
    save=[]
    ###########################Zero function
    my_interval = pybedtools.BedTool(TEMPFolder+'input.bed')
    if len(my_interval) == 0:
        sys.stderr.write("\r# Excluder error: The interval file is empty.\n\n")
        exit()
    for line in my_interval:
        try:
            d_openvalue = Open_file.values(line.chrom, line.start, line.end, numpy=True)
        except:
            sys.stderr.write("\n# Excluder error: pyBigWig package is not installed completely, please install it using conda.\n To assure it is installed properly, import pyBigWig in python, if you get '0' for 'pyBigWig.numpy' command, please install it using conda as following:\nconda install -c bioconda pybigwig\n")
            exit()
        Zeros=AverageofRegion(d_openvalue, WSize, Step_size)
        for i in Zeros:
            string=("%s\t%s\t%s") % (line.chrom,i[0]*Step_size+line.start,i[1]*Step_size+line.start)
            save.append(string)
        pybedtools.BedTool(save).saveas(TEMPFolder+'zeroz.bed')

# ---------------Excluder definition------------
# function that imports the blacklist merged file and exclude from the interval
def Excluder(BL, output, TEMPFolder, Open_file, Step_size):

    """ In the case of a blacklist presents, the Excluder function loads blacklist intervals and input.bed
    file, then store the intersect of them as a one-bed file.

    Keyword arguments:
    BL -- A class of interval(s) which the user wants to exclude.
    Step_size -- An integer number defined by the user, the size of the data reading from the bigWig file.
    input -- A folder/file address of the input file in Bed format, produced in interval converter function.
    output -- An string contains 'temp/intervals.bed', used as the output of the module.

    Function description:
    Both bed file Blacklist and input file loads and blacklist intervals exclude using the BedTool subtract 
    function, the result stored in a temporary folder ("temp/" by default).
    """
    import pybedtools
    save =[]
    bedfile = pybedtools.BedTool(TEMPFolder+'input.bed')

    # define a temporary file and put the list of the chromosome name inside it, because the Bedtools does not accept the list as a array
    temporaryfile = open(TEMPFolder+'genome.txt',"w")
    for i in Open_file.chroms().items():
        temporaryfile.write(i[0]+"\t"+str(i[1])+"\n")
    temporaryfile.close()
    tmp =  open(TEMPFolder+"genome.txt")
    # Check if the Bed file present, continue
    try:
        # Interval data derived from Bed file sorted and merged if their steps overlapped
        sites = (bedfile.sort(g = tmp.name).merge(d=Step_size-1))
    except:
        # for bedtools sort Version: v2.26.0
        sites = (bedfile.sort(faidx = tmp.name).merge(d=Step_size-1))
    if BL != []:
        save = sites.subtract(BL)
    else:
        save = sites
    Z = pybedtools.BedTool(TEMPFolder+'zeroz.bed')
    blocks = save.subtract(Z)
    pybedtools.BedTool(blocks).saveas(output)
