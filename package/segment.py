#!/usr/bin/env python
# -----------functions Definition-----------------------------------------------

# -------------boundary-finder--------------------------------------------------

def boundary_finder(CHR, START, END, STEP, INPUT):
    """
In the beginning, the program examines the sequence under study to find non-zero
sequences and the distances between them. For each region, all the extrema points
of the signal identify. Each maximum has a boundary limited by inflection points
between that maximum and neighboring minima. Coordinations of the peak boundaries
converts to the absolute coordinations of the chromosome. Each block represents
a local peak within the signal. All these candidate intervals stores in an array
with the bed file format "chromosome name, initialtion coordination, ending
coordination". Function returns the array as the identified boundaries.

Keyword arguments:
CHR: The name of the inspecting chromosome.
START: The start position of the investigating region which is corresponding to
       the start of the chromosome.
END: The end position of the investigating region which is corresponding to the
     start of the chromosome.
STEP: The STEP size of the reading the input data.

Function description:
The function searches within the convolved signal for finite numbers and names
them "scaffolds". If there is any gap between finite numbers they report as the
"GAPS" region. Indices of the scaffolds stores in "IdxScaffolds". For each
scaffold, all the minima and maxima calculate using sicpy.signal module. If the
region has only one maximum without any minimum, the boundaries around that
maximum limited by two inflection points will be reported as the candidate peak.
If a region has only one minimum without any maximum, two ranges will be
considered. The range started from the beginning of the region and limits to the
inflection point before the minima as the first range, and the range starts from
the inflection point placed after the minima and the end of the region as the
second range. The other possible situation is more common. The region has both
maxima and minima. The position of each maximum within the minima list makes
finding the peak possible. Each peak range is limited by two inflection points
between each maximum and its neighboring minima.
    """

    from itertools import groupby
    import numpy as np, scipy.signal, sys
    sys.stdout.write('\r\x1b[K Boundary tracing '+str(CHR))

    # Finite numbers within the input region stores in an array.
    scaffolds = [list(v) for k,v in groupby(INPUT,np.isfinite) if k]

    # If there is a gap between the regions stores in another array
    GAPS = np.where(np.diff(np.isfinite(INPUT)))
    GAPS = [STEP*G+START for G in GAPS]
    GAPS = np.append(GAPS, END)
    if np.isfinite(INPUT[0]):
        GAPS = np.insert(GAPS, 0, START-1, axis=0)
    for i in range(0, len(GAPS), 2): GAPS[i]+=1
    IdxScaffolds = [GAPS[i:i + 2] for i in range(0, len(GAPS), 2)]


    LIST=[]

    for S,I in zip(scaffolds,IdxScaffolds):
            A,B = 0,0

            # "Maxs" is a list of maxima within the input region.
            Maxs=scipy.signal.find_peaks(S)[0]

            # "Mins" is a list of minima within the input region
            Mins=scipy.signal.find_peaks([ -x for x in S])[0]

            # If both the minima and maxima are found within the region.
            if Maxs.size > 0 and  Mins.size > 0:

                # Index list of all the positions of the maximum in the list of the minima,
                # combined by the index of the maximum in the "Max's" list.
                Max_in_Mins=[(np.searchsorted(Mins, M), M) for M in Maxs]
                for i in Max_in_Mins:
                    min,min1,min2,max,P,P1,P2,R,R1,R2=[],[],[],[],[],[],[],[],[],[]

                    # If the first extremum is maximum, means signal starts with an increasing line
                    # leads to a local maximum.
                    if i[0] == 0:
                        min=Mins[0]
                        max=Maxs[0]
                        R1=S[:max]
                        R2=S[max:min]

                        # If there is an inflection point before the maximum, it is the left boundary of
                        # the peak. If there is no inflection point (in the case of a line with constant
                        # slope), the middle point of the increasing line will be considered as the left
                        # border of the peak.
                        try:
                            P1=np.where(np.diff(R1)==np.max(np.diff(R1)))[0]
                            if P1.size == 0:
                                P1=np.array([int(len(R1)/2)])
                        except:
                            P1=np.array([int(len(R1)/2)])

                        # The position of the left border is not an absolute chromosomeal coordinate. "P1"
                        # is an index number while the coordinates in the chromosome are a one-base number.
                        sco = I[0]+(P1[0]+1)*STEP

                        # If there is an inflection point between the maximum and the following minimum,
                        # it is the right boundary of the peak. If there is no inflection point, the middle
                        # point of the decreasing line will be considered as the right border of the peak.
                        try:
                            P2=np.where(np.diff(R2)==np.min(np.diff(R2)))[0]
                            if P2.size == 0:
                                P2=np.array([int(len(R2)/2)])
                        except:
                            P2=np.array([int(len(R2)/2)])

                        # The position of the right border is not an absolute chromosomal coordinate. "P2"
                        # is an index number while the coordinates in the chromosome are a one-base number.
                        eco = I[0]+(P2[-1]+max+1)*STEP
                        if eco > END:
                            eco = END

                        # The boundaries of the peak and the chromosome name with the "bed" file format adds
                        # to the end of an array.
                        LIST.append([CHR,int(sco),int(eco)])

                    # If the last extremum is maximum, means signal ends with an decreasing line leads
                    # from a local maximum.
                    elif i[0] == len(Mins):
                        min=Mins[-1]
                        max=Maxs[-1]
                        R1=S[min:max]
                        R2=S[max:]

                         # If the distance between the maximum and the previous minimum is longer than 2
                         # positions, this part will run.
                        if len(R1) > 1:
                            # If there is an inflection point between the maximum and the previous minimum,
                            # it is the left boundary of the peak.
                            P1=np.where(np.diff(R1)==np.max(np.diff(R1)))[0]
                            P1=np.array(P1+min)

                        # Or if there is not any position between the maximum and its previous minimum, the
                        # minimum position will be considered as the left border of the peak boundary.
                        elif max-min == 1:
                            P1=np.array([min])

                        # If the distance between the maximum and the previous minimum is longer than 2
                        # positions and there is no inflection point, the middle point of the increasing
                        # line will be considered as the left border of the peak.
                        if P1.size == 0 and len(R1) >= 2:
                            P1 = np.array([int(len(R1)/2)+min])

                        # The position of the left border is not an absolute chromosomeal coordinate. "P1"
                        # is an index number while the coordinates in the chromosome are a one-base number.
                        sco = I[0]+(P1[0]+1)*STEP

                        # If the distance between the maximum and the following minimum is one position or
                        # more, the inflection between the maximum and the following minimum is the right
                        # border of the peak boundary.
                        if len(R2) >= 1:
                            P2=np.where(np.diff(R2)==np.min(np.diff(R2)))[0]
                            P2 = np.array([max+P2[-1]])
                        # If there is no inflection point, the middle point of the decreasing line will be
                        # considered as the right border of the peak.
                            if P2.size == 0:
                                P2=np.array([int(len(R2)/2)])
                                P2 = np.array([max1+P2[0]])
                        else:
                            P2 = I[1]
                        # The position of the right border is not an absolute chromosomal coordinate. "P2"
                        # is an index number while the coordinates in the chromosome are a one-base number.
                        eco = I[0]+(P2[0]+1)*STEP
                        if eco > END:
                            eco = END

                        # The boundaries of the peak and the chromosome name with the "bed" file format
                        # adds to the end of an array.
                        LIST.append([CHR,int(sco),int(eco)])


                    # If the maximum is surrounded by two or more minima;
                    elif i[0] > 0:
                        min1=Mins[i[0]-1]
                        min2=Mins[i[0]]
                        max1=i[1]

                        # P1
                        # If the distance between the previous minimum and the maximum is longer than 1,
                        # the inflection point in this region is the left border of the peak boundary.
                        if max1-min1 > 1:
                            R1=S[min1:max1]
                            P1=np.where(np.diff(R1)==np.max(np.diff(R1)))[0]

                            # If there is no inflection point, the middle point of the increasing line
                            # will be considered as the left border of the peak.
                            if P1.size == 0:
                                P1 = [int(len(R1)/2)]
                            P1 = np.array([min1+P1[0]])

                            # The position of the left border is not an absolute chromosomeal coordinate.
                            # "P1" is an index number while the coordinates in the chromosome are a
                            # one-base number.
                            sco = I[0]+(P1[0]+1)*STEP

                        # P2
                        # If the distance between the maximum and the next minimum is longer than 1, the
                        # inflection point in this region is the right border of the peak boundary.
                        if min2-max1 > 1:
                            R2=S[max1:min2]
                            P2=np.where(np.diff(R2)==np.min(np.diff(R2)))[0]

                            # If there is no inflection point, the middle point of the decreasing line
                            # will be considered as the right border of the peak.
                            if P2.size == 0:
                                P2 = [int(len(R2)/2)]
                            P2 = np.array([max1+P2[0]])

                            # The position of the right border is not an absolute chromosomal coordinate.
                            # "P2" is an index number while the coordinates in the chromosome are a
                            # one-base number.
                            eco = I[0]+(P2[0]+1)*STEP

                            # The boundaries of the peak and the chromosome name with the "bed" file
                            # format adds to the end of an array.
                            if eco > sco:
                                LIST.append([CHR,int(sco),int(eco)])

                        # If there is an increasing line after the last local minimum, if there is an
                        # inflection point in this region, the area between the inflection point and the
                        # end of the region will be considered as a peak boundary.
                        if i[0] == len(Mins)-1 & i[0]==Max_in_Mins[-1][0]:
                            R3=S[min2:]
                            P3=np.where(np.diff(R3)==np.max(np.diff(R3)))[0]
                            if P3.size > 0:
                                sco = I[0]+(P3[0]+min2)*STEP
                                eco = I[1]

                                # The boundaries of the peak and the chromosome name with the "bed" file
                                # format adds to the end of an array.
                                LIST.append([CHR,int(sco),int(eco)])

                    # If the signal starts with a decreasing line that leads to a local Minimum,
                    # the inflection point in this line will determine the right border of the
                    # peak boundary.
                    if Max_in_Mins[0][0] == 1:
                        min=Mins[0]
                        R=S[:min]
                        P=np.where(np.diff(R)==np.min(np.diff(R)))[0]

                        if P.size > 0:
                            sco = I[0]
                            eco = I[0]+(P[-1]+1)*STEP

                            # The boundaries of the peak and the chromosome name with the "bed" file
                            # format adds to the end of an array.
                            LIST.append([CHR,int(sco),int(eco)])

            # If the region has only one maximum without any minimum, the boundaries around
            # that maximum limited by two inflection points will be reported as the candidate
            # peak.
            elif Maxs.size > 0:
                for j in range(Maxs.size):
                    max=Maxs[j]
                    R1=S[:max]
                    R2=S[max:]

                    # If there is an inflection point before the maximum, it is the left border of
                    # the peak boundary. If there is no inflection point, the middle point of the
                    # increasing line will be considered as the left border of the peak.
                    if len(R1) > 2:
                        P1=np.where(np.diff(R1)==np.max(np.diff(R1)))[0]
                        if P1.size == 0:
                            indx=[int(len(R1)/2)]
                            P1 = np.where(S==S[indx])
                    else:
                        P1 = [max]

                    # If there is an inflection point after the maximum, it is the right border of
                    # the peak boundary. If there is no inflection point, the middle point of the
                    # increasing line will be considered as the left border of the peak.
                    if len(R2) > 2:
                        P2=np.where(np.diff(R2)==np.min(np.diff(R2)))[0]
                        if P2.size == 0:
                            P2=[int(len(R2)/2)]

                        #  The positions of the left and right borders are not the absolute chromosomal
                        # coordinate, but they are index numbers while the coordinates in the chromosome
                        # are a one-base number.
                        sco=I[0]+(P1[0]+1)*STEP
                        eco=I[0]+(P2[0]+1+max)*STEP

                        # The boundaries of the peak and the chromosome name with the "bed" file format
                        # adds to the end of an array.
                        if eco > sco:
                            LIST.append([CHR,int(sco),int(eco)])

            # If the region has only one minimum without any maximum, two ranges will be considered.
            # The range started from the beginning of the region and limits to the inflection point
            # before the minima as the first range, and the range starts from the inflection point
            # placed after the minima and the end of the region as the second range.
            elif Mins.size > 0:
                min=Mins[0]
                R1=S[:min]
                R2=S[min:]
                P1=np.where(np.diff(R1)==np.min(np.diff(R1)))[0]
                P2=np.where(np.diff(R2)==np.max(np.diff(R2)))[0]

                if P1.size > 0:
                    sco=I[0]
                    eco=I[0]+(P1[-1]+1)*STEP

                    # The boundaries of the peak and the chromosome name with the "bed" file format adds
                    # to the end of an array.
                    LIST.append([CHR,int(sco),int(eco)])

                if len(P2) > 1:
                    sco=I[0]+(min+P2[1]+1)*STEP
                    eco=I[1]
                    if eco > END:
                        eco = END

                    # The boundaries of the peak and the chromosome name with the "bed" file format adds
                    # to the end of an array.
                    LIST.append([CHR,int(sco),int(eco)])

    return(LIST)
# -------------The-end-of-boundary-finder---------------------------------------


# ---------------------Sort-index-----------------------------------------------
def SortIdx(LIST):
    """
The indices of a list after concatenation list could become disordered. This
function sorts the "LIST" array and returns the sorted and duplicate removed
items.

Keyword arguments:
LIST: An unsorted array of the starting coordinates for filling the "thickStart"
field of the bed12 file format.

Function description:
The Dict.__getitem__ is actually equivalent to lambda x: Dict[x]. This fuction
repeated for the number of the elements in the list and takes an array and using
the dictionary removes the duplicated values. in the end all the elements will
be sorted.
    """
    return (sorted(range(len(LIST)), key=LIST.__getitem__))
# ---------------The-end-of-Sort-index------------------------------------------
