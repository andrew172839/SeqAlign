#!/usr/bin/perl

#define the sequence in 3-letter and 1-letter formats
@amino3 = ("GLY", "ALA", "VAL", "LEU", "ILE", "SER", "THR", "CYS", "PRO", "PHE", "TYR", "TRP", "HIS", "ASP", "ASN", "GLU", "GLN", "MET", "LYS", "ARG");
@amino1_upper = ("G", "A", "V", "L", "I", "S", "T", "C", "P", "F", "Y", "W", "H", "D", "N", "E", "Q", "M", "K", "R");
@amino1_lower = ("g", "a", "v", "l", "i", "s", "t", "c", "p", "f", "y", "w", "h", "d", "n", "e", "q", "m", "k", "r");
$namino = @amino3;
#set up random number generator
srand;
#check the command-line argumnets
if (@ARGV != 2) {
	printf "shared_change.pl <arg1> <arg2>\n";
	printf "<arg1>: input sequences file\n";
	printf "<arg2>: reference sequence\n";
	exit;
}
#check if the input sequence file exists
if (not -f $ARGV[0]) {
	printf "error: $ARGV[0] doesn't exist\n";
	exit;
}
#read in the input sequences
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
#check if the germline sequence file exists
if (not -f $ARGV[1]) {
	printf "error: $ARGV[1] doesn't exist\n";
	exit;
}
#read in the germline sequences
open(NUC, $ARGV[1]);
@nuctxt = <NUC>;
$nuclen = @nuctxt;
close(NUC);
for ($i = 0; $i < $nuclen; $i++) {
	if ($nuctxt[$i] =~ /^>/) {
		chomp $nuctxt[$i];
		$grmcom = $nuctxt[$i];
		for ($ntseq = "", $j = $i + 1; $j < $nuclen; $j++) {
			last if ($nuctxt[$j] =~ /^>/);
			chomp $nuctxt[$j];
			$ntseq = $ntseq.$nuctxt[$j];
		}
		$grmseq = $ntseq;
	}
}
#random number for temporary file name
$rndnam = int(rand 100000);
open(TMP, ">$rndnam.fasta");
#output input sequences sequences
for ($j = 0; $j < $nseq1; $j++) {
	printf TMP ("$nuccom[$j]\n");
	$strlen = length($nucseq[$j]);
	for ($k = 0; $k < $strlen; $k++) {
		printf TMP ("%1s", substr($nucseq[$j], $k, 1));
		if (int(($k + 1) / 60.0) * 60 == $k + 1) {
			printf TMP ("\n");
		}
	}
	if (int($strlen / 60.0) * 60 != $strlen) {
		printf TMP ("\n");
	}
}
printf TMP (">reference\n");
$strlen = length($grmseq);
for ($k = 0; $k < $strlen; $k++) {
	printf TMP ("%1s", substr($grmseq, $k, 1));
	if (int(($k + 1) / 60.0) * 60 == $k + 1) {
		printf TMP ("\n");
	}
}
if (int($strlen / 60.0) * 60 != $strlen) {
	printf TMP ("\n");
}
close(TMP);
#generate multiple sequence alignment using clustalw
system "clustalw $rndnam.fasta >& clustalw.out";
#read in the multiple sequence alignment file
($sharemut) = &readmsa("$rndnam.aln", "reference");
#clean up the scene
printf "shared mutations with respect to the reference: $sharemut\n";

#read in clustalw2 multiple sequence alignment
sub readmsa() {
	my $imsa, $jmsa, $kmsa;
	my @arr_aln, $narr_aln;
	my @alntxt, $alnlen;
	my $nam_tmp, $aln_tmp;
	#initialization
	$nmsa = 0;
	$sharemut = 0;
	#open multiple sequence alignment (msa) file
	open(ALN, "$_[0]");
	@alntxt = <ALN>;
	$alnlen = @alntxt;
	close(ALN);
	#get the reference sequence index
	$tgtnam = $_[1];
	#read in the clustalw2 msa file
	for ($nmsa = 0, $imsa = 0; $imsa < $alnlen; $imsa++) {
		chomp $alntxt[$imsa];
		$strlen = length($alntxt[$imsa]);
		#skip the title line
		next if ($alntxt[$imsa] =~ /CLUSTAL/);
		#skip the empty line or consensus line
		next if ($strlen == 0 or substr($alntxt[$imsa], 0, 1) eq " ");
		#cut the line into two parts: name and sequence
		@arr_aln = split / +/, $alntxt[$imsa];
		$narr_aln = @arr_aln;
		#check the for mat, for clustal, it has to be 2
		if ($narr_aln != 2) {
			printf STDERR "alignment error: $alntxt[$imsa]\n";
			exit;
		}
		#get the sequence name and alignment fragment
		$nam_tmp = $arr_aln[0];
		$aln_tmp = $arr_aln[1];
		#insert current read into alignment list
		for ($find = 0, $jmsa = 0; $jmsa < $nmsa; $jmsa++) {
			if ($seq_nam_arr[$jmsa] eq $nam_tmp) {
				$seq_aln_arr[$jmsa] = $seq_aln_arr[$jmsa].$aln_tmp;
				$find = 1;
			}
		}
		if ($find == 0) {
			$seq_aln_arr[$nmsa] = $aln_tmp;
			$seq_nam_arr[$nmsa] = $nam_tmp;
			$nmsa++;
		}
	}
	#alignment range
	$range = length($seq_aln_arr[0]);
	#locate the reference sequence
	for ($tgtidx = -1, $imsa = 0; $imsa < $nmsa; $imsa++) {
		if ($seq_nam_arr[$imsa] eq $tgtnam) {
			$tgtidx = $imsa;
		}
	}
	#calculate the shared mutations
	for ($imsa = 0; $imsa < $range; $imsa++) {
		#get nt / aa of a site in the reference seq
		$char0 = substr($seq_aln_arr[$tgtidx], $imsa, 1);
		#skip if this site is a gap
		next if ($char0 eq "-");
		#initialize counting array
		for ($jmsa = 0; $jmsa <= $namino; $jmsa++) {
			$ncount[$jmsa] = 0;
		}
		#count nt / aa of this site in all input seqs
		for ($jmsa = 0; $jmsa < $nmsa; $jmsa++) {
			#skip the reference
			next if ($jmsa == $tgtidx);
			#get nt/aa of this site in an input seqs
			$char1 = substr($seq_aln_arr[$jmsa], $imsa, 1);
			for ($kmsa = 0; $kmsa < $namino; $kmsa++) {
				if ($char1 eq $amino1_upper[$kmsa]) {
					$ncount[$kmsa]++;
				}
			}
			#if it is a gap in input, add to the last
			if ($char1 eq "-") {
				$ncount[$namino]++;
			}
		}
		#ranking ...
		&sort2($namino + 1, @ncount);
		#shared mutations should be:
		#(1) the largest number of repeats = the number of input seqs
		#(2) the consensus nt/aa in input seqs = nt / aa in reference seq
		if ($ncount[$bakptr[$namino]] == $nmsa - 1 and $amino1_upper[$bakptr[$namino]] ne $char0) {
			$sharemut++;
		}
	}
	return ($sharemut);
}

#sort a series of numbers
sub sort2() {
	my($n, @list) = @_;
	my($i, $j, $k, $index, $keys, $lists, @key);
	#initialize index into the original ordering
	for ($i = 0; $i < $n; $i++) {
		$key[$i] = $i;
	}
	#perfor m the heapsort of the input list
	$k  =  $n / 2;
	$index  =  $n - 1;
	#begin the sorting
	while ($n > 0) {
		if ($k > 0) {
			$k = $k - 1;
			$lists = $list[$k];
			$keys = $key[$k];
		}
		else {
			$lists = $list[$index];
			$keys = $key[$index];
			$list[$index] = $list[0];
			$key[$index] = $key[0];
			$index = $index - 1;
			if ($index <= 0) {
				$list[0] = $lists;
				$key[0] = $keys;
				@bakptr = @key;
				return;
			}
		}
		$i = $k;
		$j = 2 * ($k + 1) - 1;
		while ($j <= $index) {
			if ($j < $index && $list[$j] < $list[$j + 1]) {
				$j++;
			}
			if ($lists < $list[$j]) {
				$list[$i] = $list[$j];
				$key[$i] = $key[$j];
				$i = $j;
				$j = 2 * ($i + 1) - 1;
			}
			else {
				$j = $index + 1;
			}
		}
		$list[$i] = $lists;
		$key[$i] = $keys;
	}
	@bakptr = @key;
	return;
}
