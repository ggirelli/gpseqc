
# Produce a bed of bins

BEGIN {
	OFS = FS = "\t";
}

{
	for ( i = 0; i < $2; i += step ) {
		print $1 OFS i OFS i+size;
	}
}
