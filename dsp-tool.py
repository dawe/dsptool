import sys, time, os, pathlib, re
s_time = time.time()
my_data = {}
# In the config file, some of the variables are T/F
def str2bool(v):
  return v.lower() in ("True", "true", "t", "T")
try:
    file = open("config.ini","r")
    for line in file:
        if not line.startswith('#')and line.strip():
            c_splited = line.strip().split(' ')
            if c_splited[2] == 'True' or  c_splited[2] == 'False':
                my_data[c_splited[0]] = str2bool(c_splited[2])
            else:
                my_data[c_splited[0]] = c_splited[2:]
    file.close()
except:
    my_data['Config.ini']=False
if my_data['Config.ini']:
    filter=my_data['filter']
    span=my_data['span']
    step=my_data['step']
    HeatMap=my_data['HeatMap']
    Segmen=my_data['Segmentation']
    Individuals=my_data['Individuals']
    Denoised=my_data['Denoised']
    Boundaries=my_data['Boundaries']
    read_window_size=my_data['read_window_size']
    Entire=my_data['Entire']
    Regions=my_data['Regions']
    Intervals=my_data['Intervals']
    Peaks=my_data['Peaks']
    Heatmapnames=my_data['Heatmapnames']
    if Entire == True:
        for n in range (len(Individuals)):
            print(' - Denoising ', Individuals[n])
            run='python filter.py -i '+str(Individuals[n])+' -o '+str(Denoised[n])+' -f '+filter[0]+' -s '+read_window_size[0]+' -S '+step[0]
            os.system(run)
    elif Intervals != False:
        for n in range (len(Individuals)):
            for m in range (len(Intervals)):
                print(' - Denoising ', Individuals[n], Intervals[m])
                run='python filter.py -i '+str(Individuals[n])+' -o '+str(Denoised[n])+' -f '+filter[0]+' -s '+read_window_size[0]+' -l '+Intervals[m]+' -S '+step[0]
                os.system(run)
    else:
        for n in range (len(Individuals)):
            for m in range (len(Regions)):
                print(' - Denoising ', Individuals[n], Regions[m])
                run='python filter.py -i '+str(Individuals[n])+' -o '+str(Denoised[n])+' -f '+filter[0]+' -s '+read_window_size[0]+' -r '+Regions[m]+' -S '+step[0]
                os.system(run)
    if Segmen:
        if Entire == True:
            for n in range (len(Individuals)):
                print(' - Peak finding ', Individuals[n])
                run='python segmentation.py -i '+str(Denoised[n])+' -p '+str(Peaks[n])+' -b '+Boundaries[n]+' -S '+step[0]
                os.system(run)
        elif Intervals != False:
            for n in range (len(Individuals)):
                for m in range (len(Intervals)):
                    #print('Peak finding ', Individuals[n], Intervals[m])
                    run='python segmentation.py -i '+str(Denoised[n])+' -p '+str(Peaks[n])+' -b '+Boundaries[n]+' -S '+step[0]+' -l '+Intervals[m]
                    os.system(run)
        else:
            for n in range (len(Individuals)):
                for m in range (len(Regions)):
                    #print('Peak finding ', Individuals[n], Regions[m])
                    run='python segmentation.py -i '+str(Denoised[n])+' -p '+str(Peaks[n])+' -b '+Boundaries[n]+' -S '+step[0]+' -r '+Regions[m]
                    os.system(run)
        if HeatMap:
            BWfiles= ' '.join(Denoised)
            Bedfiles= ' '.join(Boundaries)
            files = '%s %s' % (Bedfiles, BWfiles)
            print(' - Plotting HeatMap file')
            run='python heatmap.py '+files+' -S '+step[0]
            os.system(run)

    print(' Completed in',"%.2f seconds" % (time.time()-s_time))
    exit()


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
print(' Completed in',"%.2f seconds" % (time.time()-s_time))
