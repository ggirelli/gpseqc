
# Pre-process for metric.
#
# Args:
#	type (string): 'p' for probability, 'ff' for Fano factor, and 'cv' for coefficient of variation

BEGIN {
	OFS = FS = "\t";
}

{
	switch(type) {
		case "p": {
			if ( 0 == $1 || 0 == $3 ) {
				print "nan";
			} else {
				print $2 / ($1 * $3);
			}
			break;
		}
		case "ff": {
			if ( 0 == $1 ) {
				print "nan";
			} else {
				print $2 ** 2 / $1;
			}
			break;
		}
		case "cv": {
			if ( 0 == $1 ) {
				print "nan";
			} else {
				print $2 / $1;
			}
			break;
		}
	}
}
