# -*- coding: utf-8 -*-

'''
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
@description: cutsite methods.
'''

# DEPENDENCIES =================================================================

import os
import pybedtools as pbt

from gpseqc import bed

# FUNCTIONS ====================================================================

def build(bedfiles, mode, csbed = None, groups = None):
    '''Build bed file with list of cutsites.
    The cutsite domain can be specified as follows:
    1 - all genomic cutsites (universe)
    2 - all cutsites restricted in the experiment (union)
    3 - all cutsites restricted in a condition (separate)
    4 - all cutsites restricted in all conditions (intersection)
    Default is 3 (separate). Also, note that if option 1 is selected, an
    additional argument -l is required.

    If option 1 (universe) is used, either a cutsite list or a group list is
    used as the cutsite domain. If both are provided, only the group list is
    considered, where each group is used as a cutsite.

    Args:
        bedfiles (list): list of bed files.
        mode (int): custite domain mode (see description).
        csbed (pbt.BedTools): parsed bed with all cutsite in the reference.
        groups (pbt.BedTools): parsed bed with groups, if grouping.
    '''

    assert mode in range(1, 5), "allowed modes: %s" % str(list(range(1, 5)))

    # Identify common cutsites/groups lists for domain preparation -------------

    if 1 == mode:   # Universe
        if type(None) == type(groups):
            # Use cutsite list
            assert_msg = "parsed cutsite bed required, got: %s" % type(csbed)
            assert type(pbt.BedTool()) == type(csbed), assert_msg
        else:
            # Use groups as cutsites
            assert_msg = "parsed groups bed required, got: %s" % type(groups)
            assert type(pbt.BedTool()) == type(groups), assert_msg
            return(groups)
    elif 2 == mode: # Union
        # Identify cutsites from all conditions
        print("> Performing union over beds...")
        csbed = bedfiles[0].cat(*bedfiles[1:], force_truncate = True)
    elif 3 == mode: # Separated
        csbed = None
    elif 4 == mode: # Intersection
        print("> Performing intersection over beds...")
        csbed = bedfiles[0]
        for i in range(len(bedfiles[1:])):
            csbed += bedfiles[i]
            csbed = csbed.cut(range(3))

    # Group domain if needed and save counts -----------------------------------

    if type(None) != type(csbed) and type(None) != type(groups):
        # Group cutsites
        csbed = bed.to_bins(groups, csbed, noValues = True)
        csbed = csbed.cut(range(3)).merge(d = -1)

    return(csbed)

def apply(bedfiles, csbed):
    '''Apply cutsite domain to bed files.

    Args:
        bedfiles (list):.
        csbed (pbt.BedTool):.
    '''

    if type(None) == type(csbed):
        print("Removing empty sites...")
        for i in range(len(bedfiles)):
            # Remove zero-loci or empty groups
            s = ""
            with open(bedfiles[i].fn, "r+") as IH:
                for line in IH:
                    if line.startswith("track"): continue
                    interval = line.strip().split("\t")
                    if float(interval[4]) != 0: s += line
            bedfiles[i] = pbt.BedTool(s, from_string = True)
    else:
        print("Applying cutsite domain...")
        for bid in range(len(bedfiles)):
            bedfiles[bid] = bed.to_bins(csbed, bedfiles[bid], skipEmpty = False)
            bedfiles[bid] = bedfiles[bid].merge(
                d = -1, o = "first,sum", c = "4,5")

    return(bedfiles)

# END ==========================================================================

################################################################################
