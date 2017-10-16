#!/usr/bin/perl

#set random number generator
srand;
#emboss transeq function
$transeq = "/usr/local/emboss/bin/transeq";
#check the command-line argumnets
if (@ARGV != 3) {
	printf "bclust_first_seq.pl <arg1> <arg2> <arg3>\n";
	printf "pick the first sequence of each cluster\n";
	printf "listed in blastclust output\n";
	printf "<arg1>: input sequence file used for clustering\n";
	printf "<arg2>: blastclust output file\n";
	printf "<arg3>: output sequence file\n";
	exit;
}
#check if the nucleotide sequence file exists
if (not -f $ARGV[0]) {
	printf "error: $ARGV[0] doesn't exist\n";
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
#pull out sequence indexes
for ($i = 0; $i < $nseq1; $i++) {
	#parse the comment line
	$strlen = length($nuccom[$i]);
	$strtmp = substr($nuccom[$i], 1, $strlen - 1);
	@arr1 = split / +/, $strtmp;
	$narr1 = @arr1;
	#store the first item of comment line
	$iseqidx[$i] = $arr1[0];
}
#check if the clustering results exist
if (not -f $ARGV[1]) {
	printf "error: $ARGV[1] doesn't exist!\n";
	exit;
}
#read in clustering results
open(CLS, $ARGV[1]);
@clstxt = <CLS>;
$clslen = @clstxt;
$nclust = @clstxt;
close(CLS);
$nclust = $clslen;
#store the clustering info into 2D array
for ($i = 0; $i < $clslen; $i++) {
	chomp $clstxt[$i];
	#parse each cluster
	@arr = split / +/, $clstxt[$i];
	$narr = @arr;
	#number of sequences in current cluster
	$cseqnum[$i] = $narr;
	#store cluster sequence index
	for ($j = 0; $j < $cseqnum[$i]; $j++) {
		$cseqidx[$i][$j] = $arr[$j];
	}
}
#open ouput sequence file
open(NEW, ">$ARGV[2]");
#go through each cluster
for ($ncount = 0, $i = 0; $i < $nclust; $i++) {
	#pick the #1 sequence of each cluster
	for ($j = 0; $j < 1; $j++) {
		#go through all the sequences
		for ($k = 0; $k < $nseq1; $k++) {
			#if sequence indexes (names) match
			if ($cseqidx[$i][$j] eq $iseqidx[$k]) {
				#output sequences
				printf NEW "$nuccom[$k]\n";
				$strlen = length($nucseq[$k]);
				for ($l = 0; $l < $strlen; $l++) {
					printf NEW ("%1s", substr($nucseq[$k], $l, 1));
					if (int(($l + 1) / 60.0) * 60 == $l + 1) {
						printf NEW ("\n");
					}
				}
				if (int($strlen / 60.0) * 60 != $strlen) {
					printf NEW ("\n");
				}
				#update the sequence counter
				$ncount++;
				#output to the screen
				printf "$i $ncount\n";
			}
		}
	}
}
#close the output sequence file
close(NEW);

#log 10 function
sub log10 {
	my $n = shift;
	return log($n) / log(10);
}
