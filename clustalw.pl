#!/usr/bin/perl

#clustalw2 executable
$clustalw2_exe = "clustalw2";
#check the command-line arguments
if (@ARGV != 2) {
	printf STDERR "clustalw.pl <template.fa> <target.fa>\n";
	printf STDERR "will output pir for mat alignment file\n";
	exit;
}
#check if template.fa exists
if (not -f $ARGV[0]) {
	printf STDERR "error: template - $ARGV[0] - doesn't exist\n";
	exit;
}
#basename of the template
$strlen = length($ARGV[0]);
$wantlen = index $ARGV[0], ".fa";
$tempnam = substr($ARGV[0], 0, $wantlen);
#read in the template sequence
&readfa($ARGV[0]);
$tempseq = $sequence;
#check if target.fa exists
if (not -f $ARGV[1]) {
	printf STDERR "error: target - $ARGV[1] - doesn't exist\n";
	exit;
}
#basename of the target
$strlen = length($ARGV[1]);
$wantlen = index $ARGV[1], ".fa";
$targnam = substr($ARGV[1], 0, $wantlen);
#read in the target sequence
&readfa($ARGV[1]);
$targseq = $sequence;
#seed the random number generator
srand;
#get a random number for file name
$randnum = int(rand 1000000);
#generate a temporary sequence file
open(SEQ, ">$randnum.fa");
printf SEQ ">$tempnam\n";
printf SEQ "$tempseq\n";
printf SEQ ">$targnam\n";
printf SEQ "$targseq\n";
close(SEQ);
#run sequence alignment using clustalw2
system "$clustalw2_exe -INFILE = $randnum.fa -OUTFILE = $targnam\_on_$tempnam.pir -OUTPUT = PIR >& /dev/null";
#clean up the scene
unlink "$randnum.fa", "$randnum.dnd";
#read in the alignment file
&readaln("$targnam\_on_$tempnam.pir");
#rewrite the alignment file
open(PIR, ">$targnam\_on_$tempnam.pir");
printf PIR ">p1; $tempnam\n";
printf PIR "structureX: $tempnam: 1:::: t1: t2: t3: \n";
printf PIR "$alignment[0]*\n";
printf PIR ">P1; $targnam\n";
printf PIR "sequence: $targnam: 1:::: q1: q2: q3: \n";
printf PIR "$alignment[1]*\n";
close(PIR);

#read in protein sequence and store in a global array
sub readfa() {
	my $ifa, $jfa, $sequence;
	my $strlen, $aa_tmp;
	my $aastr = "GgAaVvLlIiSsTtCcPpFfYyWwHhDdNnEeQqMmKkRrXx";
	open(FAS, "$_[0]");
	@fastxt = <FAS>;
	$faslen = @fastxt;
	close(FAS);
	$sequence = "";
	for ($if a = 0; $ifa < $faslen; $if a++) {
		chomp $fastxt[$ifa];
		$strlen = length($fastxt[$ifa]);
		next if ($fastxt[$ifa] =~ /^>/);
		for ($jfa = 0; $jfa < $strlen; $jfa++) {
			$aa_tmp = substr($fastxt[$ifa], $jfa, 1);
			if ($aastr =~ /$aa_tmp/) {
				$sequence = $sequence.$aa_tmp;
			}
		}
	}
};

#read in sequence alignment and store in a global array
sub readaln() {
	my $ialn, $jaln, $alnseq1, $alnseq2;
	my $strlen, $aa_tmp;
	open(ALN, "$_[0]");
	@alntxt = <ALN>;
	$alnlen = @alntxt;
	close(ALN);
	for ($idx = -1, $ialn = 0; $ialn < $alnlen; $ialn++) {
		chomp $alntxt[$ialn];
		$strlen = length($alntxt[$ialn]);
		#count when meet the ">" sign and skip
		if ($alntxt[$ialn] =~ /^>/) {
			$idx++; next;
		};
		#skip the empty line and *
		next if ($strlen == 0 || $alntxt[$ialn] =~ /^\*/);
		$alignment[$idx] = $alignment[$idx].$alntxt[$ialn];
	}
};
