
# Adds information to the bed file.
#
# Args:
#	cnr (int): number of reads in the condition.
#	bfi (int): bed file ID.
#
# Returns:
#	bed-like table with added columns.
#	cnr as #3 column, the rest shifted of 1 and bfi as last.

BEGIN {
	OFS = FS = "\t";
}

{
	$4 = cnr OFS $4;
	print $0 OFS bfi;
}
