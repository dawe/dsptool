#!/usr/bin/env python
# -----------functions Definition-----------------------------------------------

# ----------------clustering-evaluation-----------------------------------------

def evaluate(BW, Bed, cutfile, CLT):
    """
This function evaluates the candidate regions with a clustering method based on
the OTSU method in order to discard the background noise from the results.

Keyword arguments:
BW: An input file in the format of pyBigWig contains the header and intensity
    values for a smoothed signal.
Bed: An array of the intervals with the "bed" file format, contains chromosome
     name, start and end coordinates for each interval.
cutfile: The file name and the location of the output file with the "bed" file
        format.
CLT: The folder location to save the plotted png images.

Function description:
The intensity array of the candidate regions evaluates with the OTSU method and one threshold defines by this method. All the regions that have an intensity mean higher than this threshold will keep as the enriched regions and the regions that have intensity levels equal or lower than this threshold will be discarded.
    """

    import numpy as np
    import pybedtools, itertools, os, sys
    import matplotlib.pyplot as plt
    from skimage import filters

    # Assess the presence of the folder for plotting the images.
    try:
        os.mkdir(CLT)
    except:
        pass

    # Extract the list of chromosome names from input file and the interval regions.
    CHRbw ,CHRbed = list(BW.chroms()), [x.chrom for x in Bed]

    # List of the chromosome names requested for the investigation which are present
    # in the input file saves in an array.
    gen = [x for x in CHRbw if x in CHRbed]
    linenumber, truesignal = 0, []

    # Two arrays with intensity value means and length of the fragments created.
    for Chromosome in gen:
        sys.stdout.write('\r\x1b[K Thresholding '+str(Chromosome))
        array_of_intensity = np.array([np.mean(BW.values(i[0], int(i[1]), int(i[2]), numpy=True)[::50]) for i in Bed if i[0]==Chromosome])
        array_of_length = np.array([int(i[2])-int(i[1]) for i in Bed if i[0]==Chromosome])

        # A list of non-zero indices for intensity list
        nZ = np.nonzero(array_of_intensity)[0]

        # Creates an array of non-intensity values.
        ARRAYint=array_of_intensity[nZ]

        if len(ARRAYint)>1:

            # Apply the OTSU method to the array of intensities.
            otsu_thr_i = filters.threshold_otsu(np.sort(ARRAYint))

            # Define a threshold value for the foreground intensity level.
            T=otsu_thr_i

            # Filter the candidate regions with the threshold.
            FGI = np.where((array_of_intensity > T))
            truesignal.extend([(i+linenumber) for i in FGI[0]])

            # The candidate regions with intensities equal or less than the threshold
            # stores in an array for plot.
            NFGI = np.where((array_of_intensity <= T))


            # R=[array_of_intensity[i] for i in FGI[0]]

            # Plot the distribution of the candidates according to their length and
            # intensity levels before and after thresholding.
            plt.figure(figsize=(12, 10))
            plt.scatter(np.log10(array_of_length[FGI[0]]), np.log10(array_of_intensity[FGI[0]]), 13, color='darkorange')
            plt.scatter(np.log10(array_of_length[NFGI[0]]), np.log10(array_of_intensity[NFGI[0]]), 13, color='navy')
            plt.xlabel('log10 Length of fragment',size=18)
            plt.title('Clusters of the identified fragments within '+Chromosome+' using OTSU',size=18)
            plt.ylabel('log10 Mean of fragment intensities',size=18)
            plt.ylim(bottom=0,top=5)
            plt.xlim(left=0,right=6.5)
            sv = ("%s/%s") % (CLT, Chromosome)
            plt.savefig(sv+'.png', bbox_inches='tight', dpi=300)
            plt.close('all')

            plt.figure(figsize=(12, 10))
            plt.scatter(np.log10(array_of_length[FGI[0]]), np.log10(array_of_intensity[FGI[0]]), 13, color='darkorange')
            plt.xlabel('log10 Length of fragment',size=18)
            plt.title('Clusters of the identified fragments within '+Chromosome+' using OTSU',size=18)
            plt.ylabel('log10 Mean of fragment intensities',size=18)
            plt.ylim(bottom=0,top=5)
            plt.xlim(left=0,right=6.5)
            sv = ("%s/%s") % (CLT, Chromosome)
            plt.savefig(sv+'after.png', bbox_inches='tight', dpi=300)
            plt.close('all')
        linenumber = len(array_of_intensity)+linenumber

    # Creates a new BedTool array with only intervals at lines indices.
    TS=Bed.at(truesignal)
    pybedtools.BedTool(TS).saveas(cutfile)

# ----------------end-of-clustering-evaluation----------------------------------

# ------------------run-the-module-stand-alone----------------------------------
# This module could be used for the clustering the regions of a bigWig file if the
# required arguments provided by the user.

if __name__ == "__main__":
    import pyBigWig, pybedtools, sys
    try:
        CLT=str(sys.argv[4])
        cutfile=str(sys.argv[3])
        BW = pyBigWig.open(sys.argv[1])
        Bed = pybedtools.BedTool(sys.argv[2])
    except:
        exit("python cut.py input.bw input.bed output.bed plot_folder components_number")
    evaluate(BW, Bed, cutfile, CLT)
