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

def is_overlapping(bed):
    '''Check if a bed contains overlapping features.

    Args:
        bed (pbt.BedTools): parsed bed file.

    Returns:
        bool
    '''
    isect = bed.intersect(bed, loj = True)
    return(bed.count() != isect.count())

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
    '''Normalize one bed over another. Discards empty intersections. The two bed
    files are expected to be sorted and contain exactly the same regions.

    Args:
        normbed (pbt.BedTool): normbed bed.
        bed (pbt.BedTool): parsed bed.

    Returns:
        pbt.BedTool: normalized bed.
    '''

    assert_msg = "the two bed files are expected to have the same regions."
    assert normbed.count() == bed.count(), assert_msg

    # Sort bed files
    normbed = normbed.sort()
    bed = bed.sort()

    # Normalize
    s = ""
    for i in range(bed.count()):
        # Skip if no reads in the normbed
        b = float(normbed[i]['score'])
        if b == 0: continue

        tmp = list(bed[i][:4])
        a = float(bed[i]['score'])
        tmp.append("%.2f" % (a / b))
        s += "\t".join(tmp) + "\n"

    return(pbt.BedTool(s, from_string = True))

def to_bins(bins, bed, noValues = False, skipEmpty = True):
    '''Assign regions to bins. Each bin will appear once per each intersecting
    region, with the region value field appended. 

    If the bins are non-overlapping, for each region in the bed only the largest
    intersection is considered.

    Each region in the output will have the value 'row_XXX' in the name column.

    Args:
        bins (pbt.BedTool): bins bed.
        bed (pbt.BedTool): parsed bed.
        noValues (bool): do not output value column.
        skipEmpty (bool): do not output empty bins.

    Returns:
        pbt.BedTool: grouped bed.
    '''

    if not noValues:
        assert_msg = "missing score column, run with 'noValues = True'."
        assert bed.field_count() >= 5, assert_msg
        bed = bed.cut(range(5)).sort() # Force to BED5

    # Enforce bins to BED3
    bins = bins.cut(range(3)).sort()

    # Perform intersection
    isect = bins.intersect(bed, wao = True)

    d = {}              # Output container
    bi = 1              # Region counter

    # Iterate intersection
    def parsegen(isect):
        with open(isect.fn, "r+") as IH:
            for line in IH:
                yield line.strip().split("\t")
    gen = (i for i in parsegen(isect))
    if skipEmpty: gen = (i for i in gen if float(i[-2]) >= 0)

    if not is_overlapping(bins): # Retain only largest intersections
        def lab(i):
            return("%s:%s-%s" % tuple(i[3:6]))

        def d_update(d, i, bi):
            data = i[:3]
            data.append("row_%d" % bi)
            data.append(i[-2] if float(i[-2]) >= 0 else "0")
            
            if lab(i) == ".:-1--1": d["row_%d" % bi] = (int(i[-1]), bi, data)
            else: d[lab(i)] = (int(i[-1]), bi, data)

            return(d)

        for i in gen:
            # Retain only largest intersection
            if lab(i) in d.keys():
                if int(i[-1]) > d[lab(i)][0]: d = d_update(d, i, d[lab(i)][1])
            else:
                d = d_update(d, i, bi)
                bi += 1

    else: # Retain all intersections
        for i in gen:
            data = i[:3]
            data.append("row_%d" % bi)
            data.append(i[-2])
            d[bi] = (int(i[-1]), bi, data)
            bi += 1


    # Assemble
    d = "\n".join(["\t".join(x[2]) for x in d.values()])

    # Format as bed file
    d = pbt.BedTool(d, from_string = True)
    if not noValues: d = d.cut(range(5))

    return(d)

def to_combined_bins(bins, bed, fcomb = "sum"):
    '''Groups UMI-uniqued reads from bed file to bins. For each region in bed,
    only the largest intersection is considered. Each region in the output will
    have the value 'row_XXX' in the name column.

    Args:
        bins (pbt.BedTool): bins bed.
        bed (pbt.BedTool): parsed bed.
        fcomb (str): one of the bedtools merge methods:
                        sum, min, max, absmin, absmax,
                        mean, median, mode, antimode
                        stdev, sstdev

    Returns:
        pbt.BedTool: grouped bed.
    '''

    methods = ['sum', 'min', 'max', 'absmin', 'absmax', 'mean', 'median',
        'mode', 'antimode', 'stdev', 'sstdev']
    assert_msg = "unrecognized combination method, got '%s'." % fcomb
    assert_msg += " Available: %s" % methods
    assert fcomb in methods, assert_msg

    # Enforce bins to BED3
    bins = bins.cut(range(3))
    bed = bed.cut(range(5))

    # Perform intersection
    isect = bins.intersect(bed, wao = True)

    if not is_overlapping(bins): # Retain only largest intersection

        # Retain only largest intersections ------------------------------------

        d = {}
        with open(isect.fn, "r+") as IH:
            for line in IH:
                i = line.strip().split("\t")
                if float(i[-2]) < 0: continue

                # Retain only largest intersection
                if i[6] in d.keys():
                    if int(i[-1]) > d[i[6]][0]:
                        d[i[6]] = (int(i[-1]), line)
                else: d[i[6]] = (int(i[-1]), line)

        isect = "".join([v[1] for v in d.values()])
        isect = pbt.BedTool(isect, from_string = True)
    
    # Combine intersections ----------------------------------------------------

    # Merge and reset naming system
    isect = isect.merge(d = -50, c = "7,8", o = "first,%s" % fcomb)
    tmp = ""
    for i in range(len(isect)):
        line = isect[i]
        line[3] = "row_%d" % (i + 1)
        tmp += "\t".join(line) + "\n"

    return(pbt.BedTool(tmp, from_string = True))

# END ==========================================================================

################################################################################
