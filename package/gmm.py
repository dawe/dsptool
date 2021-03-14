#!/usr/bin/env python
# -----------functions Definition-----------------------------------------------

# ----------------clustering-evaluation-----------------------------------------

def evaluate(BW, Bed, GMMFILE, CLT, CN):
    """
This function evaluates the candidate regions with a clustering method based on
the GMM in order to discard the background noise from the results.

Keyword arguments:
BW: An input file in the format of pyBigWig contains the header and intensity
    values for a smoothed signal.
Bed: An array of the intervals with the "bed" file format, contains chromosome
     name, start and end coordinates for each interval.
GMMFILE: The file name and the location of the output file with the "bed" file
        format.
CLT: The folder location to save the plotted png images.
CN: The number of the components used for clustering. This variable could
    contain the string "auto".

Function description:
Intensity array of the candidate regions clustered by the number of the
components defined by the user or automatically computed by the BIC method. The
clustering performed by the GMM method. The cluster that has the lowest mean
value discards from the candidate list and the rest of the candidate will be
saved in a bed file format. The distributions of the two arrays before and after
clustering plots in the scatter plots for each chromosome.
    """
    from sklearn import mixture
    import numpy as np
    import pybedtools, itertools, os, sys
    import matplotlib.pyplot as plt

    # Assess the presence of the folder for plotting the images.
    try:
        os.mkdir(CLT)
    except:
        pass
    color_iter = itertools.cycle(['green','darkorange', 'navy', 'red', 'cornflowerblue'])
    Tlist, Real_subregions=[],[]

    # Extract the list of chromosome names from input file and the interval regions.
    CHRbw ,CHRbed = list(BW.chroms()), [x.chrom for x in Bed]

    # List of the chromosome names requested for the investigation which are present
    # in the input file saves in an array.
    gen = [x for x in CHRbw if x in CHRbed]

    # Two arrays with intensity value means and length of the fragments created.
    for Chromosome in gen:
        sys.stdout.write('\r\x1b[K Clustering '+str(Chromosome))
        array_of_intensity = np.array([np.mean(BW.values(i[0], int(i[1]), int(i[2]), numpy=True)[::50]) for i in Bed if i[0]==Chromosome])
        array_of_length = np.array([int(i[2])-int(i[1]) for i in Bed if i[0]==Chromosome])

        # A list of intervals found in the input file.
        chrbed = [i for i in Bed if i[0]==Chromosome]

        # A list of non-zero indices for intensity list
        nZ = np.nonzero(array_of_intensity)[0]

        # Log10 scaling of the two arrays.
        ARRAYint=np.log10(array_of_intensity[nZ])
        ARRAYlen=np.log10(array_of_length[nZ])

        # An array filled by "1" to use as the second dimension of interval array.
        ARY=np.ones(len(nZ))

        if len(ARRAYint)>5:
            # Converts the 1D-array to 2D-array using vertical stack function.
            ARRAY2D = np.vstack([ARRAYlen,ARRAYint]).transpose()
            ARRAY = np.vstack([ARY,ARRAYint]).transpose()

            # If the number of the components did not specify by the user, there is
            # a function that finds the component number fits the distribution using
            # the BIC.
            if CN=="auto":
                Component_number = Components(np.sort(ARRAYint))
            else:
                Component_number = int(CN)

            # The main function of the clustering using the GMM method.
            gmm = mixture.GaussianMixture(Component_number, covariance_type='full', n_init=10)
            gmm.fit(ARRAY)

            # The cluster with the lowest mean value identified.
            Downleftcomponent = np.argmin(gmm.means_, axis=0)

            #assign every data point to one of the clusters (Predict the labels for the data samples in X using the trained model.)
            PP = gmm.predict(ARRAY)

            # The indices of the elements of the cluster with the lowest mean value stores in an array.
            Downleft = np.where(PP == Downleftcomponent[1])[0]

            # Plot the distribution of the candidates according to their length and
            # intensity levels.
            plt.figure(figsize=(12, 10))
            for i, (mean, color) in enumerate(zip(gmm.means_, color_iter)):
                plt.scatter(ARRAY2D[PP == i, 0], ARRAY2D[PP == i, 1], 13, color=color)
            plt.xlabel('log10 Length of fragment',size=18)
            plt.title(str(Component_number)+' clusters of the identified fragments within '+Chromosome+' using GMM',size=18)
            plt.ylabel('log10 Mean of fragment intensities',size=18)
            plt.ylim(bottom=0,top=5)
            plt.xlim(left=0,right=6.5)
            sv = ("%s/%s") % (CLT, Chromosome)
            plt.savefig(sv+'.png', bbox_inches='tight', dpi=300)

            # An array of the candidates passed the clustering procedure.
            intlist = [str(b) for b,t in zip(chrbed,nZ) if t not in Downleft]
            Real_subregions.extend(intlist)
            plt.close('all')

    # Save the foreground regions list in a "bed" file format.
    pybedtools.BedTool(Real_subregions).saveas(GMMFILE)

    # Plot the distribution of the foreground regions according to their length
    # and intensity levels.
    for Chromosome in gen:
        sys.stdout.write('\r\x1b[K Ploting '+str(Chromosome)+" ")
        BED = pybedtools.BedTool(GMMFILE)
        array_of_intensity = np.array([np.mean(BW.values(i[0], int(i[1]), int(i[2]), numpy=True)[::50]) for i in BED if i[0]==Chromosome])
        array_of_length = np.array([int(i[2])-int(i[1]) for i in BED if i[0]==Chromosome])
        nZ = np.nonzero(array_of_intensity)[0]
        if len(nZ)>0:
            ARRAY2D = np.log10(np.vstack([array_of_length[nZ],array_of_intensity[nZ]]).transpose())
            plt.figure(figsize=(12, 10))
            plt.scatter(ARRAY2D[:, 0], ARRAY2D[:, 1], 13, color="green")
            plt.title('Remaining fragments of '+str(Chromosome)+" after deleting the lowest cluster",size=18)
            plt.xlabel('log10 Length of fragment',size=18)
            plt.ylabel('log10 Mean of fragment intensities',size=18)
            plt.ylim(bottom=0,top=5)
            plt.xlim(left=0,right=6.5)
            sv = ("%s/%s") % (CLT, Chromosome)
            plt.savefig(sv+'after.png', bbox_inches='tight', dpi=300)
            plt.close('all')
    return()

# ----------------end-of-clustering-evaluation----------------------------------

# --------------components-number-determination---------------------------------

def Components(ARRAY):
    """
This function computes the best number of components using the bayesian
information criterion (BIC) based on the likelihood function. It is a criterion
for model selection among a finite set of models. The model with the lowest BIC
is preferred.

Keyword arguments:
ARRAY: A 1D-array contains the intensity levels used as the input.

Function description:
The GMM clutering applies on the input data for components number between 2 and
6, the BIC for each clutering profile will be calculated. Resulting BIC(Bayesian
Information Criterion) for different number of components creates a curve with
the BIC values (y-axis) across all models with different compnents (x-axis). The
elbow of the curve is the lowest number of component with the respect to the
lowest BIC value.
    """
    from sklearn import mixture
    from kneed import KneeLocator
    import numpy as np
    cn_full = []

    # It is exported to have less than 6 clusters in the input data.
    n_components_range = range(1, 7)
    for n_components in n_components_range:

        # Fit the input array with the Gaussian mixture model with 10 times of initiation
        # with all the possible cluster shapes.
        gmm = mixture.GaussianMixture(n_components=n_components, covariance_type='full', n_init=10)
        gmm.fit(ARRAY[: , None])
        cn_full.append(gmm.bic(ARRAY[: , None]))

    # Finds the knee position of the curve drawn by BIC numbers. This position gives the
    # lowest number of components whereby with the increase of the cluster numbers, the
    # structure of the clusters does not change at an equal rate.
    kn = KneeLocator(range(1, 7), cn_full, curve='convex', direction='decreasing')
    Component_number=range(1, 7)[kn.knee]

    return(Component_number)

# ------------end-of-components-number-determination----------------------------

# ------------------run-the-module-stand-alone----------------------------------
# This module could be used for the clustering the regions of a bigWig file if the
# required arguments provided by the user.

if __name__ == "__main__":
    import pyBigWig, pybedtools, sys
    from sklearn import mixture
    import numpy as np
    import pybedtools, itertools, os
    import matplotlib.pyplot as plt

    try:
        CLT=str(sys.argv[4])
        GMMFILE=str(sys.argv[3])
        BW = pyBigWig.open(sys.argv[1])
        Bed = pybedtools.BedTool(sys.argv[2])
    except:
        exit("python gmm.py input.bw input.bed output.bed plot_folder components_number")
    try:
        cn=int(sys.argv[5])
    except:
        cn="auto"
    evaluate(BW, Bed, GMMFILE, CLT, cn)
