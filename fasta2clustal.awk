#!/usr/bin/awk -f
BEGIN{print "CLUSTAL W"; print ""} /^>/{if(head){print head"\t"seq; seq=""}; head=substr($1, 2, length($1))} !/^>/{seq=seq""$1} END{print head"\t"seq}
