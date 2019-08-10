import sys, time
s_time = time.time()
arguments = sys.argv[1:]
import filter as fl

if '-seg' in arguments:
    import segmentation as SG

    if fl.Target == 'Region':
        ZeroCrossList=SG.Sobel_filters(fl.d_convolve)
        ListofMaxima, ListofMinima, ListofHighPlain = SG.Maxima(ZeroCrossList)
        MergedList=SG.Merge(ListofMaxima,ListofHighPlain)
        PeaksIntensity=SG.Callregions(fl.d_regname, MergedList, fl.d_span, fl.d_step,fl.Header,fl.d_convolve,fl.d_regs)
print('segmentation completed in',"%.2f seconds" % (time.time()-s_time))
