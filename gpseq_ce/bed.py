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
                data[ilabel].append(float(i[3]))
            else:
                data[ilabel] = [float(i[3])]

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
    '''Normalize one bed over another.
    Discards empty intersections.

    Args:
        normbed (pbt.BedTool): normbed bed.
        bed (pbt.BedTool): parsed bed.

    Returns:
        pbt.BedTool: normalized bed.
    '''
    isect = normbed.intersect(bed, wb = True)
    
    s = ""
    with open(isect.fn, "r+") as IH:
        for line in IH:
            i = line.strip().split("\t")
            if 0 < int(i[4]):
                data = i[5:-1]
                data.append(int(i[9]) / float(i[4]))
                s += "%s\t%s\t%s\t%s\t%.2f\n" % tuple(data)

    return(pbt.BedTool(s, from_string = True))

def read_or_check(bed):
    '''Checks if a bed file was already parsed, and parses it otherwise.

    Args:
        bed (str): path to bed file, needs to be fully stored in memory.
        bed (pbt.BedTool): parsed bed file content.

    Returns:
        pbt.BedTools: parsed bed file content.
    '''
    if type("") == type(bed):
        assert os.path.isfile(bed), "missing file: %s" % bed
        bed = pbt.BedTool(bed)

    assert_msg = "either path or BedTool() expected."
    assert type(pbt.BedTool()) == type(bed), assert_msg

    return(bed)

def sort_by(bed, fkey = None):
    '''Sort a bed file. By default, sorts by chromosome ID.
    It is possible to provide a custom key function for sorting.

    Args:
        bed (str): path to bed file, needs to be fully stored in memory.
        bed (pbt.BedTool): parsed bed file content.
        fkey (fun): function to extract sorting keys.

    Returns:
        list: list of intervals sorted by key.
    '''
    bed = bed_read_or_check(bed)
    if type(None) == type(fkey): fkey = lambda x: chr2chrid(x)
    return(sorted(bed, key = fkey))

def to_bins(bins, bed):
    '''Assign regions to bins. Each bin will appear once per each intersecting
    region, with the region value field appended.

    Args:
        bins (pbt.BedTool): bins bed.
        bed (pbt.BedTool): parsed bed.

    Returns:
        pbt.BedTool: grouped bed.
    '''

    # Perform intersection
    isect = bins.intersect(bed, wa = True, wb = True, loj = True)

    # Sum read counts
    d = []
    bi = 1
    with open(isect.fn, "r+") as IH:
        for line in IH:
            i = line.strip().split("\t")
            data = i[:3]
            if float(i[7]) < 0: data.append("0")
            else: data.append(i[7])
            d.append("\t".join(data))

    # Format as bed file
    return(pbt.BedTool("\n".join(d), from_string = True))

def to_combined_bins(bins, bed, fcomb = None):
    '''Groups UMI-uniqued reads from bed file to bins.

    Args:
        bins (pbt.BedTool): bins bed.
        bed (pbt.BedTool): parsed bed.
        fcomb (fun): lambda(x,y) for combining.

    Returns:
        pbt.BedTool: grouped bed.
    '''

    if type(None) == type(fcomb):
        fcomb = lambda x, y: x + y

    # Perform intersection
    isect = bins.intersect(bed, wa = True, wb = True, loj = True)

    # Sum read counts
    d = {}
    bi = 1
    with open(isect.fn, "r+") as IH:
        for line in IH:
            i = line.strip().split("\t")
            i[7] = float(i[7])
            ilabel = " ".join(i[:3])
            if not ilabel in d.keys():
                if i[7] < 0: i[7] = 0
                d[ilabel] = [i[0], int(i[1]), int(i[2]),
                    "row_%d" % bi, i[7]]
                bi += 1
            else:
                d[ilabel][4] = fcomb(d[ilabel][4], i[7])

    # Format as bed file
    s = "\n".join(["%s\t%d\t%d\t%s\t%d" % tuple(v) for v in d.values()])
    return(pbt.BedTool(s, from_string = True))

# END ==========================================================================

################################################################################
