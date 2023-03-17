#!/usr/bin/gawk -f

BEGIN {OFS = "\t"}

/^#/ {print}

/^##contig=<ID=Gm..,/ {match($0, /Gm../, chr) ; match($0, /[0-9]{6,}/, N); chrom_length[chr[0]] = N[0]}

! /^#/ {if ( ($2 - length($4) > 150 && $2 - length($5) > 150 ) && ($2 + length($4) < (chrom_length[$1] - 150) && $2 + length($5) < (chrom_length[$1] - 150) ) ) {print}}

