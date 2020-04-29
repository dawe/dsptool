import sys


import matplotlib.pyplot as plt
def PLOT(ARRAY,SAVENAME,COLPOS,COLNAME,AXTITLE='',TITLE='Title',XLab='xlable'):
    fig1 = plt.figure(dpi=300)
    fig1.suptitle(TITLE, fontsize=14, fontweight='bold')
    ax1 = fig1.add_subplot(111)
    ax1.boxplot(ARRAY)
    plt.xticks(COLPOS,COLNAME)
    ax1.set_title(AXTITLE)
    ax1.set_xlabel(XLab)
    ax1.set_ylabel('ylabel')
    plt.savefig(SAVENAME, bbox_inches='tight')

S,R,G,D,P = sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5]
TITLE,AXTITLE,COLNAME=' '.join(sys.argv[6].split(",")),' '.join(sys.argv[7].split(",")),sys.argv[8].split(",")
COLPOS=list(range(1,len(COLNAME)+1))
# S,R,G,D=input[0],input[1],input[2],input[3] # Direct run in snakemake
OSM=[float(i.strip()) for i in open(S,"r")]
ORM=[float(i.strip()) for i in open(R,"r")]
OGM=[float(i.strip()) for i in open(G,"r")]
ODM=[float(i.strip()) for i in open(D,"r")]
xlab, xl='# of nonZero regions OSM=', []
xlab+=str(len([i for i in OSM if i > 0.0]))
xlab+=', ORM='
xlab+=str(len([i for i in ORM if i > 0.0]))
xlab+=', OGM='
xlab+=str(len([i for i in OGM if i > 0.0]))
xlab+=', ODM='
xlab+=str(len([i for i in ODM if i > 0.0]))

array =[OSM,ORM,OGM,ODM]
PLOT(array,P,COLPOS,COLNAME,AXTITLE,TITLE,xlab)
