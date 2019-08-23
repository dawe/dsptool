# !/usr/bin/env python3.7
Target, d_warning, esc, read_window_size ='Entire', False, False, 50
PACKAGES=['pyBigWig', 'scipy', 'sys', 'argparse', 're', 'os', 'tempfile', 'numpy', 'pybedtools', 'pathlib', 'scipy', 'UliEngineering','gc','skimage']#, 'numba']
# -----------library import-----------------------------------------------------------------------------
# Check whether all the required packages are installed or not, one by one
for i in PACKAGES:
    try:
        # Import the package if exist
        module_obj = __import__(i)
        # create a global object containging our module
        globals()[i] = module_obj
    except:
        print('The required package ',i,' is not present.\n')
        question = input("Would you like to install it? (Y/N)").lower()
        if question == 'y':
            try:
                import pip, subprocess, sys
                # Install the missing package using pip, then import it
                subprocess.call([sys.executable, "-m", "pip", "install", i])
                module_obj = __import__(i)
                # create a global object containging our module
                globals()[i] = module_obj
            except:
                if i == 'pyBigWig':
                    os.system('conda install pybigwig -c bioconda')
                elif i == 'pybedtools':
                    os.system('conda install --channel conda-forge --channel bioconda pybedtools')
                elif i == 'skimage':
                    os.system('pip3 install opencv-python')
        else:
            print('Please intsall \"',i,'\" package and retry.')
            exit()
import scipy.signal
gc.enable()
# -----------function Definitions-----------------------------------------------------------------------
# The function increase the accuracy of sampling when Step is higher than 1. Two parameters are used in this fuction, VALUE is a list of intensities that loaded from input file related to the specific regions and STEP is the distance between each sample through out this regions.
def AverageofRegion(VALUE, STEP=50):
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
        Result=d_matrix.mean(axis=1)
        return(Result)

# -----------available filters list-----------------------------------------------------------------
# Define a class for the optional parameter "--list-filters"
class _ListAction(argparse.Action):
    def __init__(self,option_strings,dest=argparse.SUPPRESS,default=argparse.SUPPRESS,help=None):
        super(_ListAction, self).__init__(option_strings=option_strings,dest=dest,default=default,nargs=0,help=help)
    def __call__(self, parser, namespace, values, option_string=None):
        sys.stdout.write(" Avalable filters for Scipy Signal: \n hann, hamming, blackman \n The Blackman filter is selected by defult.\n")
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
parser.add_argument("-i", "--input", help="Define the name and location of the BigWig file as input. e.g. -i data/File.bw ",required=True)
parser.add_argument("-o", "--output", help="Define the name and the location of the output BigWig file", required=True)
parser.add_argument("-f", "--filter", help="Insert the filter name of Scipy Signal filter list. Blackman filter selected by default", required=False, default="blackman")
parser.add_argument("-s", "--size", type=int, help="Define the window size of applying filter", required=False, default="65536")
parser.add_argument("-r", "--region", help="A single section of input file could be defined as chromosomename:startindex:endindex")
parser.add_argument("-l", "--interval", help="A list of regions could be defined as a BED file format")
parser.add_argument("-S", "--step", type=int, help="The distance between the end point of each region from the relevant start point as basepair, the default step size is 50bps.", default=50)
parser.add_argument("--span", type=int, help="The intended length from the beginning of each step, the default span is 50bps.", default=50)
parser.add_argument("-L", "--list-filters", action=_ListAction, help="Print all the possible signal filters")
parser.add_argument("-H", "--highresolution", action='store_true', help="Force the program to define step and span equal to 1 to have the highest resolution as possible (One base resolution)")
parser.add_argument("-seg", "--segmentation", action='store_true', help="Run segmentation")
parser.add_argument("-p", "--peak", help="Name of the segmentation file.")
parser.add_argument("-V", "--version", action='version', version="%(prog)s (Version 1.0)")
# Use "args" instead of the subcommand "parse_args"
args = parser.parse_args()
# Store the inserted values for switches on the variables
d_inputfile = args.input
d_outputfile = args.output
d_filter = args.filter
d_size = int(args.size)
d_region = args.region
d_step = int(args.step)
d_span = int(args.span)
d_interval = args.interval

# -----------input/step validation----------------------
# If the high resolution option selected by user, it is useful for the short peaks that needs more resulations. it means all the bases of defined region will be called and processed.
if args.highresolution:
    d_step, d_span = 1, 1
elif d_step<d_span:
    # The Span size must be less or equal by Step size of sampling
    d_span=d_step
# Location of the input file stores on "my_file" variable
my_file = pathlib.Path(d_inputfile)
# Check whether the path/file presence or not
if my_file.is_file():
    pass
else:
    sys.stderr.write('File \"'+ d_inputfile +'\" is not exist. Check the file name and relative path.\n')
    exit()

# --------------size of window correction-------------
# Window size of the signal resized regaurd the input file
d_size = d_size >> int(numpy.log2(read_window_size))

# --------------filterchecker------------------------
# Store filter signal curve with the specified size within an array
try:
    d_signal=eval('scipy.signal.%s(%s)' % (d_filter,d_size))
except AttributeError:
    sys.stderr.write("The filter name "+ d_filter +" for the following command is not valid!\n scipy.signal." + d_filter)
    exit()

# --------------output validation and premission check -------------------------
# Call the input bigwig file by pyBigWig library function and store in a array
d_open = pyBigWig.open(d_inputfile)
# check if the output file cannot be created because of the output folder does not exist or permissions, an error message with clean exit will occur
try:
    d_output = pyBigWig.open(d_outputfile, "w")
except IOError:
    sys.stderr.write("\n The PATH/FILE `"+ d_outputfile +"` has not write permission, please check the path and try again. \n")
    exit()
# If file created successfully, set the header of output file similar as input
else:
    d_output.addHeader(list(d_open.chroms().items()))

# --------------names.txt------------------------
# Define the temp list by all the chromosome names that present within the input file (BigWig)
chrlist=list(d_open.chroms())
Header=list(d_open.chroms().items())

# --------------region refinement-------------------
# If switch "r" or "Region" defined by the user, the inserted value must follow the standard pattern that contains chromosome name, start and end coordinates
if d_region:
    Target = 'Region'
    try:
        # Characters before column-sign should be a chromosome name
        _splited = d_region.strip().split(':')
        d_regname = _splited[0]
        # The phrase after column-sign consists of two coordinates that separated by a dash-sign
        _splited = _splited[1].strip().split('-')
        d_regs = int(_splited[0])
        d_rege = int(_splited[1])
    except:
        sys.stderr.write('The region of interest must follow the standard pattern ChromosomeName:StartBaseIndex-EndBaseIndex\n')
        exit()
    # If the request chromosome name present in the input file process starts
    if d_regname in chrlist:
        # If the region of interest present in the input file, the program goes ahead
        if (int(d_open.chroms(d_regname)) < (d_rege)) :
            sys.stderr.write('Interval definition is incorrect. The length of the chromosome ' + str(d_regname) + ' is ' + str(d_open.chroms(d_regname)) + '.\n')
            question = input('Do you want to continue with EndCoordinate=' + str(d_open.chroms(d_regname)) + ' ? (Y/N)').lower()
            if question == 'y':
                d_rege = int(d_open.chroms(d_regname))
            else:
                print('Please correct the end coordinate and try again.')
                exit()
        try:
            # Sample once a step (defult each 50 bases)
            d_openvalue = d_open.values(d_regname, int(d_regs), d_rege, numpy=True)
            # It is important to define a proper step size, if large step size for short segment, results false negative
            if len(d_openvalue) < d_step*20:
                d_warning=True
                if len(d_openvalue) < d_step:
                    exit()
            # Instead of intensity of one point, the average of values relative to the each step used
            d_averages=AverageofRegion(d_openvalue,d_step)
            # Combine signals derived from the input value by filter with specific windows size
            d_convolve = scipy.signal.fftconvolve(d_averages, d_signal, mode="same")
            # Write the combined values for the region of intreset to the output file and close it
            d_output.addEntries(d_regname, int(d_regs), ends=int(d_rege), values=d_convolve, span=d_span, step=d_step)
            d_output.close()
        except:
            esc=True
            sys.stderr.write('Interval definition or step size is incorrect.\n')
            exit()

    else:
        sys.stderr.write('The chromosome \"'+ d_regname +'\" is not present in the input file.\n')
        exit()

# -------------- interval refinement ------------
# If switch "l" or "Interval" defined by the user, a BED file must be introduced
elif d_interval:
    Target = 'Region'
    # Store the file name and location in one array
    bed_file = pathlib.Path(d_interval)
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
                sites = list(pybedtools.BedTool(d_interval).sort(g = _tmp.name).merge(d=d_step-1))
            except:
                # for bedtools sort Version: v2.26.0
                sites = list(pybedtools.BedTool(d_interval).sort(faidx = _tmp.name).merge(d=d_step-1))
            # For each interval the value reads and stored in a varible
            for line in sites:
                _temp=str(line)
                L = _temp.strip().split()
                d_openvalue = d_open.values(L[0], int(L[1]), int(L[2]), numpy=True)
                d_regname,d_regs=L[0],int(L[1])
                # Instead of intensity of one point, the average of values relative to the each step used
                d_averages=AverageofRegion(d_openvalue, d_step)
                # It is important to define a proper step size, if large step size for short segment, results false negative
                if len(d_openvalue) < d_step*20:
                    d_warning=True
                try:
                    # Combine signals derived from the input file for each step's value by filter with specific windows size
                    d_convolve = scipy.signal.fftconvolve(d_averages, d_signal, mode="same")
                    # This process may takes time, so one message informs the process stage
                    sys.stdout.write('\r'+'Working on '+L[0])
                    # Update the standard output
                    sys.stdout.flush()
                    # Write the combined values for each interval to the output file
                    d_output.addEntries(str(L[0]), int(L[1]), values=d_convolve, span=d_span, step=d_step)
                except:
                    sys.stderr.write("Please check the Bed file \"" + d_interval + "\" .")
                    exit()
            sys.stdout.write('\r'+'Denoising completed                           \n')
            d_output.close()
        else:
            sys.stderr.write("File \"" + d_interval + "\" is not exist. Check the BED file and it\'s path. \n")
            exit()

# -------------- entire chromosome ------------
# If any specific region/interval did not define by the user, the entire inputfile value from the beginning to the end will be processed
else:
    # For each chromosome of the entire input data, the value reads and stored in a varible
    for line in chrlist:
        # This process may takes time, so one message informs the process stage
        sys.stdout.write("- Denoising chromosome "+line+" from begining to "+str(d_open.chroms(line))+"\n")
        # Update the standard output
        sys.stdout.flush()
        d_openvalue = d_open.values(line, 1, d_open.chroms(line), numpy=True)
        # Instead of intensity of one point, the average of values relative to the each step used
        d_openvalue = AverageofRegion(d_openvalue, d_step)
        # Combine signals derived from the input file for each step's value by filter with specific windows size
        try:
            d_convolve = scipy.signal.fftconvolve(d_openvalue, d_signal, mode="same")
        except:
            # In the case of memory problem
            sys.stderr.write('This step size cause the Memory Error, please increase the step size respect to the available memory.\n')
            exit()
        # Write the combined values for each interval to the output file
        d_output.addEntries(line, 1, ends= d_open.chroms(line), values=d_convolve, span=d_span, step=d_step)
    d_output.close()
# --------------  WARNING  ------------
# If the ratio of region/interval size and step size is not ideal, a warning message prints for the user
if d_warning==True:
    sys.stderr.write('Warning: The defined Span=%s is too large. It could cause an error in the results because of the number of samples. It is recommended to use the smaller step.\n' % d_span)
