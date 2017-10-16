#!/usr/bin/perl

#check command-line arguments
if (@ARGV != 2) {
	printf STDERR "aln2sto.pl <arg1> <arg2>\n";
	printf STDERR "<arg1>: clustalw alignment file\n";
	printf STDERR "<arg2>: stockholm1.0 alignment file\n";
	exit;
}
#check if clustalw alignment file exists
if (not -f $ARGV[0]) {
	printf STDERR "aln2sto.pl error: $ARGV[0] not found\n";
	exit;
}
#open clustalw alignment file
open(ALN, $ARGV[0]);
@alntxt = <ALN>;
$alnlen = @alntxt;
close(ALN);
#check the title of clustalw alignment file
if ($alntxt[0] !~ /CLUSTAL/) {
	printf STDERR "aln2sto.pl error: title line not found in alignment\n";
	exit;
}
#read in the clustalw alignment
for ($nbloc = 0, $i = 1; $i < $alnlen; $i++) {
	chomp $alntxt[$i];
	#multiple-alignment block starts with a non-blank & non-star line
	if (length($alntxt[$i]) != 0 and substr($alntxt[$i], 0, 1) ne " ") {
		for ($nseq = 0, $j = $i; $j < $alnlen; $j++) {
			chomp $alntxt[$j];
			#multiple-alignment block stops with a blank or star line
			if (length($alntxt[$j]) == 0 or substr($alntxt[$j], 0, 1) eq " ") {
				$i = $j;
				last;
			}
			#split the string
			@arr = split / +/, $alntxt[$j];
			#sequence name
			$seqname[$nseq] = $arr[0];
			#alignment segment
			$alnbloc[$nbloc][$nseq] = $arr[1];
			#sequence counter
			$nseq++;
		}
		#alignment block counter
		$nbloc++;
	}
}
#build stockholm1.0 alignment
for ($i = 0; $i < $nseq; $i++) {
	#reconstruct alignment
	$alignment[$i] = $alnbloc[0][$i];
	for ($j = 1; $j < $nbloc; $j++) {
		$alignment[$i] = $alignment[$i].$alnbloc[$j][$i];
	}
	#alignment length
	$strlen = length($alignment[$i]);
	for ($j = 0; $j < $strlen; $j++) {
		if (substr($alignment[$i], $j, 1) eq "-") {
			substr($alignment[$i], $j, 1) = ".";
		}
	}
}
#maximal length of sequence name
for ($maxnmlen = -100, $i = 0; $i < $nseq; $i++) {
	if (length($seqname[$i]) > $maxnmlen) {
		$maxnmlen = length($seqname[$i]);
	}
}
#output stockholm1.0 alignment
open(STO, ">$ARGV[1]");
printf STO "#stockholm1.0\n\n";
for ($i = 0; $i < $nseq; $i++) {
	if ($maxnmlen <= 5) {
		printf STO ("%-6s", $seqname[$i]);
	}
	elsif ($maxnmlen > 5 and $maxnmlen <= 10) {
		printf STO ("%-11s", $seqname[$i]);
	}
	elsif ($maxnmlen > 10 and $maxnmlen <= 15) {
		printf STO ("%-16s", $seqname[$i]);
	}
	elsif ($maxnmlen > 15 and $maxnmlen <= 20) {
		printf STO ("%-21s", $seqname[$i]);
	}
	elsif ($maxnmlen > 20 and $maxnmlen <= 25) {
		printf STO ("%-26s", $seqname[$i]);
	}
	elsif ($maxnmlen > 25 and $maxnmlen <= 30) {
		printf STO ("%-31s", $seqname[$i]);
	}
	elsif ($maxnmlen > 30 and $maxnmlen <= 35) {
		printf STO ("%-36s", $seqname[$i]);
	}
	else {
		printf STDERR "aln2sto.pl error: sequence name too long\n";
		exit;
	}
	printf STO "$alignment[$i]\n";
}
printf STO "//";
close(STO);
