import sklearn.mixture
import numpy as np
import pyBigWig, pybedtools, sys
# from skimage import filters
BW = pyBigWig.open(sys.argv[1])
Bed = pybedtools.BedTool(sys.argv[2])
Tlist = []
for k in BW.chroms().items():
    array, data  = [], []
    # data = np.array([(BW.values(i[0], int(i[1]), int(i[2]), numpy=True)[::50]) for i in Bed if i[0]==k[0]])
    array = np.array([np.mean(BW.values(i[0], int(i[1]), int(i[2]), numpy=True)[::50]) for i in Bed if i[0]==k[0]])
    if len(array) > 1:
        gmm = sklearn.mixture.GaussianMixture(2, n_init=10).fit(array[: , None])
        thr = np.argmax(gmm.means_, axis=0)
        tflist = np.argmax(gmm.predict_proba(array[: , None]), axis=1) == thr
        Tlist.extend(tflist)

Real_subregions = [str(a) for a,b in zip(Bed,Tlist) if b]
pybedtools.BedTool(Real_subregions).saveas(sys.argv[3])
