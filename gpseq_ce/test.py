#!/usr/bin/env python3

# ==============================================================================
# 
# 180416 - Gabriele Girelli
# Project: GPSeq centrality estimate
# 
# Aim:
# 	Testing pybedtools package.
#	
# ==============================================================================

# DEPENDENCIES =================================================================

import pybedtools as pbt

# PARAMETERS ===================================================================

# FUNCTIONS ====================================================================

# RUN ==========================================================================

#print(pbt.check_for_bedtools())

b1 = pbt.BedTool("/media/bicroserver2/sequencing_data_2/BICRO55/bed files/TK93_1min_GG__cutsiteLoc-umiCount.bed")
b2 = pbt.BedTool("/media/bicroserver2/sequencing_data_2/BICRO55/bed files/TK94_5min_GG__cutsiteLoc-umiCount.bed")

print(b1.intersect(b2))

# END ==========================================================================

################################################################################
