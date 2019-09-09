
import sys, re, os
import numpy as np
import pyBigWig as bw
BED,ID,step,dictionary=[],[],50,{}
import seaborn as sn
import pandas as pd
from pylab import savefig
arguments = sys.argv[1:]

# ---------------------------- Function Definitions --------------------------------
def findbed(w):
    return re.search('.bed$', w)

# FInd the bigwig file names from the arguments
def findbw(w):
    return re.search('.bw$', w)

# Read BED, split them and list all the coordinates within the dictionary with the key of chromosome name
def insertbed(name, dic):
    if os.stat(name).st_size == 0:
        return()
    c=[]
    with open(name)as f:
        for line in f:
            L = line.strip().split()
            c.append(L[1])
            c.append(L[2])
        dic[L[0]]=(c)
    return(dic)

# Matrix with row number = chromosomes, column number coordinates
def coordinate(item):
    x=(dictionary[item])
    x.sort(reverse = False,key=len)
    x = list(dict.fromkeys(x))
    return(x)

# Plot the heatmap from the matrix
def plot(Matrix, name):
    svm = sn.heatmap(Matrix, annot=False, cmap='coolwarm')#, linecolor='white', linewidths=1)
    figure = svm.get_figure()
    figure.savefig(name+'.png', dpi=600)
    del(svm,figure)
# --------------- End of Definitions ---------------------------------

for x in arguments:
    if findbed(x):
        BED.append(x)
    elif findbw(x):
        ID.append(x)
    elif x == "-S":
        y=(arguments.index(x))+1
        step=int(arguments[y])

if len(ID) != len(BED):
    sys.stderr.write('Please check the input files. Each individual must have one bed file and one bigwig file in correct order.\n')
    exit()
# Run the insert function and fill the dictionary
for x in BED:
    dictionary=insertbed(x,dictionary)

# For each chromosome one specific matrix created containing all the individuals relative intensities
try:
    for item in dictionary:
        Coor=coordinate(item)
        df = pd.DataFrame({'ID': Coor})
        for x in ID:
            h_input = bw.open(x)
            intensity=[]
            for co in Coor:
                end=(Coor.index(co))+1
                try:
                    value =  h_input.values(item, int(co), int(Coor[end]), numpy=True)[::step]
                except:
                    pass
                intensity.append(np.average(value))
            df[x] = intensity
        df=df.set_index(['ID'])
        #print(df)
        plot(df, str(item))
except:
    sys.stderr.write('Please check the input files. Each individual must have one bed file and one bigwig file in correct order.\n')
