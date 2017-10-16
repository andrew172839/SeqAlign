#!/usr/bin/perl

#cdrh3 length range
$minlen = 5;
$maxlen = 35;
#check the command-line argumnets
if (@ARGV != 2) {
	printf "usuage: get_all_cdr3.pl <arg1> <arg2>\n";
	printf "<arg1>: input nucleotide sequence file\n";
	printf "<arg2>: output nucleotide sequence file\n";
	exit;
}
#check if the nucleotide sequence file exists
if (not -f $ARGV[0]) {
	printf "get_all_cdr3.pl error: $ARGV[0] doesn't exist\n";
	exit;
}
#read in the nucleotide sequences
open(NUC, $ARGV[0]);
@nuctxt = <NUC>;
$nuclen = @nuctxt;
close(NUC);
for ($nseq1 = 0, $i = 0; $i < $nuclen; $i++) {
	if ($nuctxt[$i] =~ /^>/) {
		chomp $nuctxt[$i];
		$nuccom[$nseq1] = $nuctxt[$i];
		for ($ntseq = "", $j = $i + 1; $j < $nuclen; $j++) {
			last if ($nuctxt[$j] =~ /^>/);
			chomp $nuctxt[$j];
			$ntseq = $ntseq.$nuctxt[$j];
		}
		$nucseq[$nseq1] = $ntseq;
		$nseq1++;
	}
}
open(OUT, ">$ARGV[1]");
#output cdr3 sequences
for ($i = 0; $i < $nseq1; $i++) {
	#parse the comment line
	$strlen = length($nuccom[$i]);
	$strtmp = substr($nuccom[$i], 1, $strlen - 1);
	@arr1 = split / +/, $strtmp;
	$narr1 = @arr1;
	for ($j = 0; $j < $narr1; $j++) {
		#parameters for germline divergence
		if ($arr1[$j] =~ /cdrh3_seq=/ or $arr1[$j] =~ /cdrl3_seq=/) {
			@arr2 = split /=/, $arr1[$j];
			$cdrseq = $arr2[1];
		}
		if ($arr1[$j] =~ /cdrh3_len=/ or $arr1[$j] =~ /cdrl3_len=/) {
			@arr2 = split /=/, $arr1[$j];
			$cdrlen = $arr2[1];
		}
	}
	#skip cdr3 outside the range
	next if ($cdrlen < $minlen or $cdrlen > $maxlen);
	#print out cdr3 sequence
	printf OUT ("$nuccom[$i]\n");
	$strlen = length($cdrseq);
	for ($j = 0; $j < $strlen; $j++) {
		printf OUT ("%1s", substr($cdrseq, $j, 1));
		if (int(($j + 1) / 60.0) * 60 == $j + 1) {
			printf OUT ("\n");
		}
	}
	if (int($j / 60.0) * 60 != $j) {
		printf OUT ("\n");
	}
}
close(OUT);
