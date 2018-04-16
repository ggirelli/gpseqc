
# Add header to the file
# 
# Args:
# 	header (string): the header line.

BEGIN {
	OFS = FS = "\t";
	print header;
}

( 1 == 1 ) {
	print $0;
}
