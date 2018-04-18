# -*- coding: utf-8 -*-

'''
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
@description: bed-related methods.
'''

# DEPENDENCIES =================================================================

import numpy as np
import os
import pandas as pd
import pybedtools as pbt

# FUNCTIONS ====================================================================

def is_overlapping(bed):
    '''Check if a bed contains overlapping features.

    Args:
        bed (pbt.BedTools): parsed bed file.

    Returns:
        bool
    '''
    isect = bed.intersect(bed, wao = True)
    return(bed.count() != isect.count())

def calc_stats(bed):
    '''Calculate bin-wise statistics.

    Args:
        bed (pbt.BedTool): parsed bed file.

    Returns:
        pd.DataFrame: data frame with bin statistics.
    '''
    
    data = {}
    with open(bed.fn, "r+") as IH:
        for line in IH:
            i = line.strip().split("\t")
            ilabel = "\t".join(i[:3])
            if ilabel in data.keys():
                data[ilabel].append(float(i[4]))
            else:
                data[ilabel] = [float(i[4])]

    s = []
    for (k, v) in data.items():
        tmp = k.split("\t")
        tmp[1] = int(tmp[1])
        tmp[2] = int(tmp[2])
        tmp.extend([sum(v), np.mean(v), np.std(v), len(v)])
        s.append(tmp)

    df = pd.DataFrame(s)
    df.columns = ["chrom", "start", "end", "sum", "mean", "std", "count"]
    df["cond_nreads"] = sum(df['sum'].values)

    return(df)

def get_chr_size(bed, d = None):
    '''Extract chromosome size from a bed file.

    Args:
        bed (str): path to bed file, needs to be fully stored in memory.

    Returns:
        dict: (chrom, size) item couples.
    '''
    
    if type("") == type(bed):
        assert os.path.isfile(bed), "missing file: %s" % bed

    if type(None) == type(d): d = {}

    with open(bed, "r+") as IH:
        for line in IH:
            if not line.startswith("track"):
                i = line.strip().split("\t")
                if not i[0] in d.keys():
                    d[i[0]] = int(i[2])
                else:
                    if int(i[2]) >= d[i[0]]:
                        d[i[0]] = int(i[2])

    return(d)

def mk_windows(chr_sizes, bsize, bstep):
    '''Generate sub-chromosome windows.

    Args:
        chr_sizes (dict): (chrom, size) item tuples.
        bsize (int): bin size.
        bstep (int): bin step.

    Returns:
        pbt.BedTools: bin bed.
    '''

    assert bsize >= bstep, "bin size must be greater than or equal to bin step."

    s = ""
    for (c, e) in chr_sizes.items():
        for start in range(0, e, bstep):
            s += "%s\t%d\t%d\t\n" % (c, start, start + bsize)
    return(pbt.BedTool(s, from_string = True))

def normalize(normbed, bed):
    '''Normalize one bed over another. Discards empty intersections.
    Consider only the first intersection of any region of bed with normbed.

    Args:
        normbed (pbt.BedTool): normbed bed.
        bed (pbt.BedTool): parsed bed.

    Returns:
        pbt.BedTool: normalized bed.
    '''
    isect = normbed.intersect(bed, wb = True)
    
    encountered = set()
    
    s = ""
    with open(isect.fn, "r+") as IH:
        for line in IH:
            i = line.strip().split("\t")

            s += line
            if i[3] in encountered: continue
            encountered.add(i[3])


            # if 0 < int(i[4]):
            #     data = i[5:-1]
            #     data.append(int(i[9]) / float(i[4]))
            #     data.extend([i[4], i[9]])
            #     s += "%s\t%s\t%s\t%s\t%.2f\t%s\t%s\n" % tuple(data)

    return(pbt.BedTool(s, from_string = True))

def to_bins(bins, bed):
    '''Assign regions to bins. Each bin will appear once per each intersecting
    region, with the region value field appended. 

    If the bins are non-overlapping, for each region in the bed only the largest
    intersection is considered.

    Each region in the output will have the value 'row_XXX' in the name column.

    Args:
        bins (pbt.BedTool): bins bed.
        bed (pbt.BedTool): parsed bed.

    Returns:
        pbt.BedTool: grouped bed.
    '''

    # Enforce bins to BED3
    bins = bins.cut(range(3))

    # Perform intersection
    isect = bins.intersect(bed, wao = True)

    d = {}              # Output container
    bi = 1              # Region counter

    # Iterate intersection
    def parsegen(isect):
        with open(isect.fn, "r+") as IH:
            for line in IH:
                yield line.strip().split("\t")

    if not is_overlapping(bins): # Retain only largest intersections
        def d_update(d, i, bi):
            ''''''
            data = i[:3]
            data.append("row_%d" % bi)
            if float(i[7]) < 0: data.append("0")
            else: data.append(i[7])
            d[i[6]] = (int(i[8]), bi, data)
            return(d)

        for i in parsegen(isect):
            # Retain only largest intersection
            if i[6] in d.keys():
                if int(i[8]) > d[i[6]][0]: d = d_update(d, i, d[i[6]][1])
            else:
                d = d_update(d, i, bi)
                bi += 1

    else: # Retain all intersactions
        for i in parsegen(isect):
            data = i[:3]
            data.append("row_%d" % bi)
            if float(i[7]) < 0: data.append("0")
            else: data.append(i[7])
            d[bi] = (int(i[8]), bi, data)
            bi += 1


    # Assemble
    d = "\n".join(["\t".join(x[2]) for x in d.values()])

    # Format as bed file
    return(pbt.BedTool(d, from_string = True))

def to_combined_bins(bins, bed, fcomb = None):
    '''Groups UMI-uniqued reads from bed file to bins. For each region in bed,
    only the largest intersection is considered. Each region in the output will
    have the value 'row_XXX' in the name column.

    Args:
        bins (pbt.BedTool): bins bed.
        bed (pbt.BedTool): parsed bed.
        fcomb (fun): lambda(x,y) for combining.

    Returns:
        pbt.BedTool: grouped bed.
    '''

    # Enforce bins to BED3
    bins = bins.cut(range(3))

    # Default combination style: sum
    if type(None) == type(fcomb): fcomb = lambda x, y: x + y

    # Perform intersection
    isect = bins.intersect(bed, wao = True)

    d2 = {}
    bi = 1  # Region counter

    if not is_overlapping(bins): # Retain only largest intersection

        # Extract largest intersections ----------------------------------------

        def d_update(d, i):
            ''''''
            data = i[:3]
            if float(i[7]) < 0: data.append("0")
            else: data.append(i[7])
            d[i[6]] = (int(i[8]), data)
            return(d)

        d = {}
        with open(isect.fn, "r+") as IH:
            for line in IH:
                i = line.strip().split("\t")

                # Retain only largest intersection
                if i[6] in d.keys():
                    if int(i[8]) > d[i[6]][0]: d = d_update(d, i)
                else: d = d_update(d, i)

        # Combine intersections ------------------------------------------------

        for (isize, i) in d.values():
            i[-1] = float(i[-1])
            ilabel = " ".join(i[:3])

            # Combine
            if not ilabel in d2.keys():
                if i[-1] < 0: i[-1] = 0
                d2[ilabel] = [i[0], int(i[1]), int(i[2]),
                    "row_%d" % bi, i[-1]]
                bi += 1
            else:
                d2[ilabel][4] = fcomb(d2[ilabel][4], i[-1])
    
    else: # Retain all intersections

        with open(isect.fn, "r+") as IH:
            for line in IH:
                i = line.strip().split("\t")

                i[7] = float(i[7])
                ilabel = " ".join(i[:3])

                # Combine
                if not ilabel in d2.keys():
                    if i[7] < 0: i[7] = 0
                    d2[ilabel] = [i[0], int(i[1]), int(i[2]),
                        "row_%d" % bi, i[7]]
                    bi += 1
                else:
                    d2[ilabel][4] = fcomb(d2[ilabel][4], i[7])

    # Format as bed file
    s = "\n".join(["%s\t%d\t%d\t%s\t%d" % tuple(v) for v in d2.values()])
    return(pbt.BedTool(s, from_string = True))

# END ==========================================================================

################################################################################
