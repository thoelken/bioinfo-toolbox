#!/usr/bin/perl

use strict;
use warnings "all";

# Variables
our (@chrom, @chromStart, @chromEnd, @name, @score, @strand, @thickStart, @thickEnd, @itemRgb, @blockCount, @blockSizes, @blockStarts, @attributes);	# columns in BED
our ($seqid, $source, $type, $start, $end, $score, $strand, $phase, $attributes);	# columns in GFF
our @filter;

# Parameters
our @files = ();
foreach my $option (@ARGV) {
	if		($option !~ /^-/)				{ push(@files,$option); }
	elsif	($option =~ m/^--?f=(.+)/)		{ @filter = split(/,/, $1); }
	else	{ print_usage(); die "Invalid command line option: \'$option\'!\n"; }
}

check_Files();

my $line_count = 0;		# counts lines with data which are used from gff file 

# convert Files
foreach my $file (@files) {
	$line_count = 0;
	open(GFF,"$file") || die ("Could not open GFF '$file': $!");
	while(my $line = <GFF>) {
		# use only rows with data
		if ($line !~ m/^#/) {
			# fill each column in that row in a Variable
			($seqid, $source, $type, $start, $end, $score, $strand, $phase, $attributes) = split(/\t/, $line);
			# check for filters and use only rows in which type match a filter value
			if (scalar(@filter) > 0) {
				for (my $i = 0, my $n = scalar(@filter); $i < $n; $i++) {
					if ($filter[$i] eq $type) {
						write_col_arrays($line_count);
						$line_count++
					}
				}
			}
			# or use every row
			else {
				write_col_arrays($line_count);
				$line_count++
			}		
		}
	}
	close(GFF);
	# write the bad file with the stored data
	my $output = "$file.bed";
	open(BED, ">", "$output") || die ("Could not open or create '$output': $!");
	for (my $i = 0, my $n = scalar(@chrom); $i < $n; $i++) {
		print BED "$chrom[$i]\t$chromStart[$i]\t$chromEnd[$i]\t$name[$i]\t$score[$i]\t$strand[$i]\t$thickStart[$i]\t$thickEnd[$i]\t$itemRgb[$i]\t$blockCount[$i]\t$blockSizes[$i]\t$blockStarts[$i]\t$attributes[$i]";	
	}
	close(BED);
}

# subroutines

# writes data from gff variables in the column arrays for the bed file
sub write_col_arrays {
	$chrom[$_[0]] = $seqid;
	$start--;
	$chromStart[$_[0]] = $start;
	$chromEnd[$_[0]] = $end;
	$name[$_[0]] = $type;
	$score[$_[0]] = $score;
	$strand[$_[0]] = $strand;
	if ($type eq "CDS")	{
		if ($phase =~ m/[0-2]/) {
			if ($strand eq "+") {
				#$start = $start + $phase;
				$thickStart[$_[0]] = $start + $phase;
				$thickEnd[$_[0]] = $end;
			}
			else {
				#$end = $end - $phase;
				$thickStart[$_[0]] = $start;
				$thickEnd[$_[0]] = $end - $phase;
			}
		}
		else {
			$thickStart[$_[0]] = $start;
			$thickEnd[$_[0]] = $end;
		}
	}
	else {
		$thickStart[$_[0]] = $start;
		$thickEnd[$_[0]] = $start;
	}
	$itemRgb[$_[0]] = ".";
	$blockCount[$_[0]] = ".";
	$blockSizes[$_[0]] = ".";
	$blockStarts[$_[0]] = ".";
	if ($attributes =~ m/.*\s.*/) {
		$attributes =~ tr/ /_/;
	}
	$attributes[$_[0]] = $attributes;
}

sub check_Files {
	if (scalar(@files) < 1) {
		print_usage();
		Error("bad_gffs.pl needs at least 1 GFF File!");
	}
}

sub print_usage {
	print STDERR << "JUS";
Only for GFF3 (especially NCBI GFF Files)!
Usage: gff2bed [OPTIONS] GFF1 [GFF...]
Options:  -f=<keyword>  filter: only convert rows which contain <keyword> in third column
                        diffrent keywords may be given seperated by comma.
example: gff2bed -f=gene,CDS,rRNA NC_*.gff
JUS
}

sub Error {
	print STDERR "Error: ".$_[0]."\n";
	exit 0;
}
