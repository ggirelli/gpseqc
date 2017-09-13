
# Merge two bed files based on first three columns

BEGIN {
	OFS = FS = "\t";
	sep = "~";
}

( FNR == NR ) {
	k = $1 sep $2 sep $3;
	a[k] = $0;
	next;
}

{
	k = $1 sep $2 sep $3;
	if ( k in a ) {
		print $0 OFS a[k];
	}
}
