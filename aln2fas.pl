#!/usr/bin/perl

#check command-line arguments
if (@ARGV != 2) {
	printf STDERR "aln2fas.pl <arg1> <arg2>\n";
	printf STDERR "<arg1>: clustalw alignment file\n";
	printf STDERR "<arg2>: fasta-format sequence file\n";
	exit;
}
#check if clustalw alignment file exists
if (not -f $ARGV[0]) {
	printf STDERR "aln2fas.pl error: $ARGV[0] not found\n";
	exit;
}
#open clustalw alignment file
open(ALN, $ARGV[0]);
@alntxt = <ALN>;
$alnlen = @alntxt;
close(ALN);
#check the title of clustalw alignment file
if ($alntxt[0] !~ /CLUSTAL/) {
	printf STDERR "aln2fas.pl error: title line not found in alignment\n";
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
#build alignment array
for ($i = 0; $i < $nseq; $i++) {
	#reconstruct alignment
	$alignment[$i] = $alnbloc[0][$i];
	#attach pieces from dif ferent blocks
	for ($j = 1; $j < $nbloc; $j++) {
		$alignment[$i] = $alignment[$i].$alnbloc[$j][$i];
	}
}
#output fasta-format sequence file
open(NEW, ">$ARGV[1]");
for ($i = 0; $i < $nseq; $i++) {
	printf NEW (">$seqname[$i]\n");
	$strlen = length($alignment[$i]);
	for ($ncount = 0, $k = 0; $k < $strlen; $k++) {
		if (substr($alignment[$i], $k, 1) ne "-") {
			printf NEW ("%1s", substr($alignment[$i], $k, 1));
			$ncount++;
			if (int($ncount / 60.0) * 60 == $ncount) {
				printf NEW ("\n");
			}
		}
	}
	if (int($ncount / 60.0) * 60 != $ncount) {
		printf NEW ("\n");
	}
}
close(NEW);
