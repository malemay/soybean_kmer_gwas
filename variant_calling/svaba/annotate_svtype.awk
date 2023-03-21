#!/usr/bin/gawk -f

BEGIN {OFS="\t"} 

/^#/ {print} 

/^##reference/ {print $0; print "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">"} 

!/^#/ {
	if(length($4) > length($5)) {
		svtype="DEL"
	} else {
	svtype="INS"
        }

	$8 = $8 ";SVTYPE=" svtype
	print
}

