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
	printf "consensus.pl <arg1> <arg2>\n";
	printf "<arg1>: input sequences file\n";
	printf "<arg2>: output consensus sequence file\n";
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
close(TMP);
#generate multiple sequence alignment using clustalw
system "clustalw2 $rndnam.fasta >& clustalw.out";
#read in the multiple sequence alignment file
($consenq) = &readmsa("$rndnam.aln");
#clean up the scene
unlink "clustalw.out", "$rndnam.fasta", "$rndnam.aln", "$rndnam.dnd";
#prepare output consensus file
open(NEW1, ">$ARGV[1]");
printf NEW1 (">c$rndnam\n");
$strlen = length($conseq);
for ($i = 0; $i < $strlen; $i++) {
	printf NEW1 ("%1s", substr($conseq, $i, 1));
	if (int(($i + 1) / 60.0) * 60 == $i + 1) {
		printf NEW1 ("\n");
	}
}
if (int($strlen / 60.0) * 60 != $strlen) {
	printf NEW1 ("\n");
}
close(NEW1);

#read in clustalw2 multiple sequence alignment
sub readmsa() {
	my $imsa, $jmsa, $kmsa;
	my @arr_aln, $narr_aln;
	my @alntxt, $alnlen;
	my $nam_tmp, $aln_tmp;
	#initialization
	$nmsa = 0;
	$conseq = "";
	#open multiple sequence alignment (msa) file
	open(ALN, "$_[0]");
	@alntxt = <ALN>;
	$alnlen = @alntxt;
	close(ALN);
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
	#calculate the consensus
	for ($imsa = 0; $imsa < $range; $imsa++) {
		#initialize counting array
		for ($jmsa = 0; $jmsa <= $namino; $jmsa++) {
			$ncount[$jmsa] = 0;
		}
		#counting ...
		for ($jmsa = 0; $jmsa < $nmsa; $jmsa++) {
			$char = substr($seq_aln_arr[$jmsa], $imsa, 1);
			for ($kmsa = 0; $kmsa < $namino; $kmsa++) {
				if ($char eq $amino1_upper[$kmsa]) {
					$ncount[$kmsa]++;
				}
			}
			if ($char eq "-") {
				$ncount[$namino]++;
			}
		}
		#ranking ...
		&sort2($namino + 1, @ncount);
		#get the consensus
		if ($bakptr[$namino] != $namino) {
			$conseq = $conseq.$amino1_upper[$bakptr[$namino]];
		}
	}
	return ($conseq);
}

#sort a series of number
sub sort2() {
	my ($n, @list) = @_;
	my ($i, $j, $k, $index, $keys, $lists, @key);
	#initialize index into the original ordering
	for ($i = 0; $i < $n; $i++) {
		$key[$i] = $i;
	}
	#perfor m the heapsort of the input list
	$k = $n / 2;
	$index = $n - 1;
	#begin sorting
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
			if ($index <=  0) {
				$list[0] = $lists;
				$key[0] = $keys;
				@bakptr = @key;
				return;
			}
		}
		$i = $k;
		$j = 2 * ($k + 1) - 1;
		while ($j <=  $index) {
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
				$j = $index  +  1;
			}
		}
		$list[$i] = $lists;
		$key[$i] = $keys;
	}
	@bakptr = @key;
	return;
}
