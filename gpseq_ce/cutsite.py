# -*- coding: utf-8 -*-

'''
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
@description: cutsite methods.
'''

# DEPENDENCIES =================================================================

import os
import pybedtools as pbt

# FUNCTIONS ====================================================================

def mk_domain(bedfiles, mode, csbed = None, groups = None):
    '''Build bed file with list of cutsites.
    The cutsite domain can be specified as follows:
    1 - all genomic cutsites (universe)
    2 - all cutsites restricted in the experiment (union)
    3 - all cutsites restricted in a condition (separate)
    4 - all cutsites restricted in all conditions (intersection)
    Default is 3 (separate). Also, note that if option 1 is selected, an
    additional argument -l is required.

    Args:
        bedfiles (list): list of bed files.
        mode (int): custite domain mode (see description).
        csbed (pbt.BedTools): parsed bed with all cutsite in the reference.
        groups (pbt.BedTools): parsed bed with groups, if grouping.
    '''

    assert mode in range(1, 5), "allowed modes: %s" % str(list(range(1, 5)))

    # Identify common cutsites/groups lists for domain preparation -------------

    csbed = None
    if 1 == mode:   # Universe
        if type(None) == type(groups):
            # Use cutsite list
            assert_msg = "cutsite bed required for mode 1."
            assert type("") == type(csbed[0]), assert_msg
            assert os.path.isfile(csbed[0]), "file not found: %s" % csbed[0]
            csbed = pbt.BedTool(csbed[0])
        else:
            # Use groups as cutsites
            assert_msg = "parsed bed tool required, got: %s" % type(groups)
            assert type(pbt.BedTool()) == type(groups), assert_msg
            csbed = groups
    elif 2 == mode: # Union
        # Identify cutsites from all conditions
        print("> Performing union over beds...")
        csbed = bedfiles[0].cat(*bedfiles[1:], force_truncate = True)
    elif 3 == mode: # Separated
        pass
    elif 4 == mode: # Intersection
        print("> Performing intersection over beds...")
        csbed = bedfiles[0]
        for i in range(len(bedfiles[1:])):
            csbed += bedfiles[i]

    # Group domain if needed and save counts -----------------------------------

    if type(None) != type(csbed):
        if type(None) != type(groups):
            # Group cutsites
            csbed = bed_to_combined_bins(groups, csbed)

    return(csbed)

def apply_domain(bedfiles, csbed, groups):
    '''Apply cutsite domain to bed files.

    Args:
        bedfiles (list):.
        csbed (pbt.BedTool):.
        groups (pbt.BedTool):.
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
        print("Applying cutsite domain... [TODO]")

    return(bedfiles)

# END ==========================================================================

################################################################################
