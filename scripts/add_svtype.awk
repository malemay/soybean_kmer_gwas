#!/usr/bin/gawk -f

# Changing the output field separator to \t
BEGIN {OFS = "\t"}

# Header lines are printed as is
/#/ {print $0}

# for other lines, the SVTYPE value is written depending on the allele lengths
!/^#/ {
    alt = length($5)
    ref = length($4)
    if(alt >= ref) {
	    svtype="INS"
    } else {
            svtype="DEL"
    }

    print $0 ";SVTYPE=" svtype
} 

