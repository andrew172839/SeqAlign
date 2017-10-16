#!/usr/bin/perl

#check the command line
if (@ARGV != 1) {
	print STDERR "seqid.pl usuage: seqid.pl <sequence_alignment_file_in_pir_format>\n";
	exit;
}
#check if the file exists
if (not -f "$ARGV[0]") {
	printf STDERR "seqid.pl error: alignment file doesn't exist\n";
	exit;
}
#read in the alignment file
open(ALI, "$ARGV[0]");
@alitxt = <ALI>;
$alilen = @alitxt;
close(ALI);
#locate the starting line
for ($ncount = 0, $i = 0; $i < $alilen; $i++) {
	if ($alitxt[$i] =~ /^>/) {
		$alistr[$ncount] = $i + 2;
		$ncount++;
	}
}
#locate the ending line
for ($ncount = 0, $i = 0; $i < $alilen; $i++) {
	if ($alitxt[$i] =~ /\*/) {
		$aliend[$ncount] = $i;
		$ncount++;
	}
}
#read in sequence 1
for ($nr1 = 0, $i = $alistr[0]; $i <= $aliend[0]; $i++) {
	$line = $alitxt[$i];
	chomp $line;
	if ($i != $aliend[0]) {
		$len = length($line);
	}
	else {
		$len = length($line) - 1;
	}
	for ($j = 0; $j < $len; $j++) {
		$sequence[0][$nr1] = substr($line, $j, 1);
		$nr1++;
	}
}
#read in sequence 2
for ($nr2 = 0, $i = $alistr[1]; $i <= $aliend[1]; $i++) {
	$line = $alitxt[$i];
	chomp $line;
	if ($i != $aliend[1]) {
		$len = length($line);
	}
	else {
		$len = length($line) - 1;
	}
	for ($j = 0; $j < $len; $j++) {
		$sequence[1][$nr2] = substr($line, $j, 1);
		$nr2++;
	}
}
if ($nr1 != $nr2) {
	printf "number of residues in sequence 1: $nr1\n";
	printf "number of residues in sequence 2: $nr2\n";
}
#calculate sequence identity
$nsame = 0;
$nqres = 0;
for ($i = 0; $i < $nr1; $i++) {
	if ($sequence[0][$i] eq $sequence[1][$i] and $sequence[0][$i] ne '-') {
		$nsame++;
	}
	if ($sequence[1][$i] ne '-') {
		$nqres++;
	}
}
$seqid = $nsame * 100.0 / $nqres;
printf "sequence identity = $seqid\n";
