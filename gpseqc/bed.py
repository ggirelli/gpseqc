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

from gpseqc import stats

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

def parse_genome_size_file(ipath):
    ''''''
    
    assert os.path.isfile(ipath), "file not found: %s" % ipath

    chrom_sizes = {}
    line_counter = 0

    with open(ipath, "r+") as IH:
        for line in IH:
            line_counter += 1

            if "\t" in line:
                tmp_array = line.strip().split("\t")

                if tmp_array[1].isdecimal():
                    chrom_sizes[tmp_array[0]] = int(tmp_array[1])
                else:
                    print("Skipped line #%d, not decimal size." % line_counter)
            else:
                print("Skipped line #%d, single field found." % line_counter)

    return(chrom_sizes)

def get_chr_size(bed, di = None):
    '''Extract chromosome size from a bed file.

    Args:
        bed (str): path to bed file, needs to be fully stored in memory.

    Returns:
        dict: (chrom, size) item couples.
    '''
    
    if type("") == type(bed):
        assert os.path.isfile(bed), "missing file: %s" % bed

    if type(None) == type(di): di = {}
    do = {}

    with open(bed, "r+") as IH:
        for line in IH:
            if not line.startswith("track"):
                i = line.strip().split("\t")
                if not i[0] in do.keys():
                    if not i[0] in di.keys():
                        do[i[0]] = int(i[2])
                    else:
                        do[i[0]] = di[i[0]]
                else:
                    if int(i[2]) >= do[i[0]]:
                        do[i[0]] = int(i[2])

    return(do)

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
    return(pbt.BedTool(s, from_string = True).sort())

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

    return(pbt.BedTool(s, from_string = True).sort())

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
    isect = bins.intersect(bed, wao = True, sorted = True)

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
            
            if lab(i) == ".:-1--1":
                d["row_%d" % bi] = (int(i[-1]), bi, data)
            else:
                d[lab(i)] = (int(i[-1]), bi, data)

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
            data.append(i[-2] if float(i[-2]) >= 0 else "0")
            d[bi] = (int(i[-1]), bi, data)
            bi += 1

    # Assemble
    d = "\n".join(["\t".join(x[2]) for x in d.values()])

    # Format as bed file
    d = pbt.BedTool(d, from_string = True).sort()
    if not noValues: d = d.cut(range(5))

    return(d)

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
    bed = bed.cut(range(5))

    # Default combination style: sum
    if type(None) == type(fcomb): fcomb = lambda x, y: x + y

    # Perform intersection
    isect = bins.intersect(bed, wao = True, sorted = True)

    d2 = {}
    bi = 1  # Region counter

    if not is_overlapping(bins): # Retain only largest intersection

        # Extract largest intersections ----------------------------------------

        def d_update(d, i):
            data = i[:3]
            data.append(i[-2])
            d[i[6]] = (int(i[-1]), data)
            return(d)

        d = {}
        with open(isect.fn, "r+") as IH:
            for line in IH:
                i = line.strip().split("\t")
                if float(i[-2]) < 0: continue

                # Retain only largest intersection
                if i[6] in d.keys():
                    if int(i[-1]) > d[i[6]][0]: d = d_update(d, i)
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

                i[-2] = float(i[-2])
                if i[-2] < 0: continue

                ilabel = " ".join(i[:3])

                # Combine
                if not ilabel in d2.keys():
                    if i[-2] < 0: i[-2] = 0
                    d2[ilabel] = [i[0], int(i[1]), int(i[2]),
                        "row_%d" % bi, i[-2]]
                    bi += 1
                else:
                    d2[ilabel][4] = fcomb(d2[ilabel][4], i[-2])

    # Format as bed file
    s = "\n".join(["%s\t%d\t%d\t%s\t%d" % tuple(v) for v in d2.values()])
    return(pbt.BedTool(s, from_string = True).sort())

def identify_outliers(bed, stype, prob = .99, lim = 1.5):
    '''Select only the outliers from a bed file.

    Args:
        bed (pbt.BedTool).
        stype (str): outlier score type. "Z" calculates normal scores, "t"
            calculates t-Student scores, "chi2" gives chi-squared scores. "IQR"
            considers only values lower than or equal to the first quartile or
            greater than or equal to the third quartile, and sets the  score to
            the ratio between the distance to the closest quartile and the IQR.
            The score for values between the first and third quartiles is set to
            0. "MAD" gives the difference between each value and the median,
            divided by the median absolute deviation.
        prob (float): if set, the corresponding p-values instead of scores are
            given. If set to 1, p-values are returned. otherwise, a logical
            vector is formed, indicating which values are exceeding the
            specified threshold. "IQR" mode does not support probabilities, use
            the lim argument instead.
        lim (float): this value can be set for "IQR" scores, to form a logical
            vector for scores that exceed the limit (after module).
    Returns:
        pbt.BedTool: outliers bed.
    '''

    scores = np.array([float(region.score) for region in bed])

    if "IQR" == stype:
        outliers = stats.score_outliers(scores, stype, lim = lim)
    else:
        outliers = stats.score_outliers(scores, stype, prob = prob)

    outliers = np.where(outliers)[0]

    record_id = 0
    outliers_data = []
    with open(bed.fn, "r+") as IH:
        for record in IH:
            if not record.startswith("track"):
                if record_id in outliers:
                    outliers_data.append(record)
                record_id += 1

    return pbt.BedTool("".join(outliers_data), from_string = True).sort()

# END ==========================================================================

################################################################################