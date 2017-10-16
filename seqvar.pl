#!/usr/bin/perl

#define the sequence in 3-letter and 1-letter formats
@amino3 = ("GLY", "ALA", "VAL", "LEU", "ILE", "SER", "THR", "CYS", "PRO", "PHE", "TYR", "TRP", "HIS", "ASP", "ASN", "GLU", "GLN", "MET", "LYS", "ARG");
@amino1_upper = ("G", "A", "V", "L", "I", "S", "T", "C", "P", "F", "Y", "W", "H", "D", "N", "E", "Q", "M", "K", "R");
@amino1_lower = ("g", "a", "v", "l", "i", "s", "t", "c", "p", "f", "y", "w", "h", "d", "n", "e", "q", "m", "k", "r");
$namino = @amino3;
#emboss transeq function
$transeq = "/usr/local/emboss/bin/transeq";
#random number generator
srand;
#check the command-line argumnets
if (@ARGV != 3) {
	printf "seqvar.pl <arg1> <arg2> <arg3>\n";
	printf "<arg1>: input nucleotide sequence file\n";
	printf "<arg2>: reference nucleotide sequence file\n";
	printf "<arg4>: basename for output file\n";
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
#check if the reference sequence exists
if (not -f $ARGV[1]) {
	printf " error: $ARGV[1] doesn't exist\n";
	exit
}
#read in the reference sequence (nt)
open(REF, $ARGV[1]);
@reftxt = <REF>;
$reflen = @reftxt;
close(REF);
for ($i = 0; $i < $reflen; $i++) {
	if ($reftxt[$i] =~ /^>/) {
		chomp $reftxt[$i];
		$refcom = $reftxt[$i];
		for ($ntseq = "", $j = $i + 1; $j < $reflen; $j++) {
			last if ($reftxt[$j] =~ /^>/);
			chomp $reftxt[$j];
			$ntseq = $ntseq.$reftxt[$j];
		}
		$refnuc = $ntseq;
	}
}
#align the 454 sequences to reference
#sequence to determine the variation
for ($i = 0; $i < $nseq1; $i++) {
	#sequence index
	$strlen = length($nuccom[$i]);
	$strtmp = substr($nuccom[$i], 1, $strlen - 1);
	@arr1 = split / +/, $strtmp;
	@arr2 = split /=/, $arr1[0];
	$pidx = $arr2[1];
	#output nucleotide sequence
	open(NUC, ">s$pidx.fasta");
	printf NUC "$nuccom[$i]\n";
	printf NUC "$nucseq[$i]\n";
	close(NUC);
	#output reference nucleotide sequence
	open(NUC, ">tref.fasta");
	printf NUC "$refcom\n";
	printf NUC "$refnuc\n";
	close(NUC);
	system "clustalw.pl tref.fasta s$pidx.fasta";
	#read in the sequence alignment
	&readpir("s$pidx\_on_tref.pir");
	#clean up the scene
	unlink "s$pidx.fasta", "tref.fasta", "s$pidx.fa", "s$pidx\_on_tref.pir";
	#output index, name, nmutate, ngap1(insertion) and ngap2 (deletion)
	printf ("%10d%10d%5d%5d%5d\n", $i, $pidx, $nmutate, $ngap1, $ngap2);
	#bin statistics
	$mutate_bin[$nmutate]++;
	$ngap1_bin[$ngap1]++;
	$ngap2_bin[$ngap2]++;
	#rms variation
	$mutate_rms += $nmutate * $nmutate;
	$ngap1_rms += $ngap1 * $ngap1;
	$ngap2_rms += $ngap2 * $ngap2;
}
#output file 1 - mutation
open(NEW1, ">$ARGV[2]_mut.dat");
for ($i = 0; $i < @mutate_bin; $i++) {
	printf NEW1 ("%10d%10d\n", $i, $mutate_bin[$i]);
}
close(NEW1);
#output file 2 - insertion
open(NEW1, ">$ARGV[2]_ins.dat");
for ($i = 0; $i < @ngap1_bin; $i++) {
	printf NEW1 ("%10d%10d\n", $i, $ngap1_bin[$i]);
}
close(NEW1);
#output file 3 - deletion
open(NEW1, ">$ARGV[2]_del.dat");
for ($i = 0; $i < @ngap2_bin; $i++) {
	printf NEW1 ("%10d%10d\n", $i, $ngap2_bin[$i]);
}
close(NEW1);
#calculate rms variation
$mutate_rms = sqrt($mutate_rms / $nseq1);
$ngap1_rms = sqrt($ngap1_rms / $nseq1);
$ngap2_rms = sqrt($ngap2_rms / $nseq1);
printf "mutation variation: $mutate_rms\n";
printf "insertion variation: $ngap1_rms\n";
printf "deletion variation: $ngap2_rms\n";

#read in .pir for mat sequence alignment file
sub readpir() {
	my $ipir, $jpir, $kpir, $iseq;
	my $tmpini, $tmpend, $qryini, $qryend;
	#initialization
	$modellable = 0;
	$mixseq = "";
	#open alignment file
	open(PIR, "$_[0]");
	@pirtxt = <PIR>;
	$pirlen = @pirtxt;
	close(PIR);
	#processing the alignment
	for ($iseq = 0, $ipir = 0; $ipir < $pirlen; $ipir++) {
		chomp $pirtxt[$ipir];
		#skip the structurex and sequence
		next if ($pirtxt[$ipir] =~ /structureX/ or $pirtxt[$ipir] =~ /sequence/);
		#start a new sequence
		if ($pirtxt[$ipir] =~ /^>/) {
			for ($naln = 0, $jpir = $ipir + 1; $jpir < $pirlen; $jpir++) {
				chomp $pirtxt[$jpir];
				#skip the structureX and sequence
				next if ($pirtxt[$jpir] =~ /structureX/ or $pirtxt[$jpir] =~ /sequence/);
				last if ($pirtxt[$jpir] =~ /^>/);
				#line length
				$strlen = length($pirtxt[$jpir]);
				#read in each character of the line
				for ($kpir = 0; $kpir < $strlen; $kpir++) {
					#skip the *
					next if (substr($pirtxt[$jpir], $kpir, 1) eq '*');
					#store the alignment
					$align[$iseq][$naln] = substr($pirtxt[$jpir], $kpir, 1);
					$naln++;
				}
			}
			$iseq++;
		}
	}
	#first residue in template
	for ($ipir = 0; $ipir < $naln; $ipir++) {
		if ($align[0][$ipir] ne "-") {
			$tmpini = $ipir;
			last;
		}
	}
	#last residue in template
	for ($ipir = $naln - 1; $ipir >= 0; $ipir--) {
		if ($align[0][$ipir] ne "-") {
			$tmpend = $ipir;
			last;
		}
	}
	#first residue in query
	for ($ipir = 0; $ipir < $naln; $ipir++) {
		if ($align[1][$ipir] ne "-") {
			$qryini = $ipir;
			last;
		}
	}
	#last residue in query
	for ($ipir = $naln - 1; $ipir >= 0; $ipir--) {
		if ($align[1][$ipir] ne "-") {
			$qryend = $ipir;
			last;
		}
	}
	#generating a mixed sequence
	$ini_shorter = $qryini;
	#query n-terminus is longer
	if ($tmpini > $qryini) {
		$ini_shorter = $tmpini;
		for ($ipir = 0; $ipir < $tmpini; $ipir++) {
			$align[2][$ipir] = "-";
		}
	}
	#template c-terminus is shorter
	$end_shorter = $qryend;
	if ($tmpend < $qryend) {
		$end_shorter = $tmpend;
		for ($ipir = $naln-1; $ipir > $tmpend; $ipir--) {
			$align[2][$ipir] = "-";
		}
	}
	#fill the sequence using query amino acids
	for ($ipir = $ini_shorter; $ipir <= $end_shorter; $ipir++) {
		$align[2][$ipir] = $align[1][$ipir];
	}
	#generate mixed sequence
	for ($ipir = 0; $ipir < $naln; $ipir++) {
		if ($align[2][$ipir] ne "-") {
			$mixseq = $mixseq.$align[2][$ipir];
		}
	}
	#sequence comparison statistics
	$nmutate = 0;
	$ngap1 = 0;
	$ngap2 = 0;
	for ($ipir = $ini_shorter; $ipir <= $end_shorter; $ipir++) {
		if ($align[0][$ipir] eq "-") {
			$ngap1++;
		}
		if ($align[1][$ipir] eq "-") {
			$ngap2++;
		}
		if ($align[0][$ipir] ne "-" and $align[1][$ipir] ne "-") {
			#count mutated nucleotides
			if ($align[0][$ipir] ne $align[1][$ipir]) {
				$nmutate++;
			}
		}
	}
}

#sort a series of number
sub sort2() {
	my($n, @list) = @_;
	my($i, $j, $k, $index, $keys, $lists, @key);
	#initialize index into the original ordering
	for ($i = 0; $i < $n; $i++) {
		$key[$i] = $i;
	}
	#perfor m the heapsort of the input list
	$k = $n / 2;
	$index = $n - 1;
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
			if ($index < = 0) {
				$list[0] = $lists;
				$key[0] = $keys;
				@bakptr = @key;
				return;
			}
		}
		$i = $k;
		$j = 2 * ($k + 1) - 1;
		while ($j < =  $index) {
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
