
# Add name column to bed file and shift columns of 1
# 
# Args:
#  prefix (string): name prefix.

BEGIN {
	OFS = FS = "\t";
}

{
	$4 = prefix NR OFS $4
	print $0;
}
