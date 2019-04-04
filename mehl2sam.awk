#!/usr/bin/awk -f

BEGIN {
	FS="\t"
}

!/^#/ {
	flag=0
	if($10=="-") {
		flag=16
	}
	cigar=""
	split($14, tmp, ";")
	for(i=0; i<length(tmp)-1; i++) {
		match(tmp[i], /^([A-Z])([0-9]+)$/, c)
		cigar=cigar""c[2]""c[1]
	}
	tags=""
	print $2"\t"flag"\t"$13"\t"$11"\t255\t"cigar"\t*\t0\t"$12-$11+1"\t*\t*\t"tags
}
