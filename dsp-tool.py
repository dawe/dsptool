import sys, time, os
s_time = time.time()
arguments = sys.argv[1:]
import filter as fl
if fl.esc:
    exit()
arg = ' '.join(arguments)
if fl.args.segmentation:
    del fl
    del sys.modules["filter"]
    os.system('python SG.py '+arg)
else:
    pass
print('Completed in',"%.2f seconds" % (time.time()-s_time))
