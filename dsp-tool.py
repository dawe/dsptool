import sys, time, os, pathlib
s_time = time.time()
folder, filename = os.path.split(os.path.abspath(__file__))
arguments = sys.argv[1:]
import filter as fl
if fl.esc:
    exit()
arg = ' '.join(arguments)
if fl.args.segmentation:
    del fl
    del sys.modules["filter"]
    segmentation=str('python '+ folder+ '/SG.py '+arg)
    os.system(segmentation)
else:
    pass
print('Completed in',"%.2f seconds" % (time.time()-s_time))
