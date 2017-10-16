#!/usr/bin/perl

#define the reasonable pdb atom name
#define a, u, g and c in rna
%index = ("ADE_P", 1, "ADE_O1P", 2, "ADE_O2P", 3, "ADE_O5'", 4, "ADE_C5'", 5, "ADE_C4'", 6,
	"ADE_O4'", 7, "ADE_C3'", 8, "ADE_O3'", 9, "ADE_C2'", 10, "ADE_O2'", 11, "ADE_C1'", 12,
	"ADE_N9", 13, "ADE_C8", 14, "ADE_N7", 15, "ADE_C5", 16, "ADE_C6", 17, "ADE_N6", 18,
	"ADE_N1", 19, "ADE_C2", 20, "ADE_N3", 21, "ADE_C4", 22,
	"URA_P", 23, "URA_O1P", 24, "URA_O2P", 25, "URA_O5'", 26, "URA_C5'", 27, "URA_C4'", 28,
	"URA_O4'", 29, "URA_C3'", 30, "URA_O3'", 31, "URA_C2'", 32, "URA_O2'", 33, "URA_C1'", 34,
	"URA_N1", 35, "URA_C2", 36, "URA_O2", 37, "URA_N3", 38, "URA_C4", 39, "URA_O4", 40,
	"URA_C5", 41, "URA_C6", 42,
	"GUA_P", 43, "GUA_O1P", 44, "GUA_O2P", 45, "GUA_O5'", 46, "GUA_C5'", 47, "GUA_C4'", 48,
	"GUA_O4'", 49, "GUA_C3'", 50, "GUA_O3'", 51, "GUA_C2'", 52, "GUA_O2'", 53, "GUA_C1'", 54,
	"GUA_N9", 55, "GUA_C8", 56, "GUA_N7", 57, "GUA_C5", 58, "GUA_C6", 59, "GUA_O6", 60,
	"GUA_N1", 61, "GUA_C2", 62, "GUA_N2", 63, "GUA_N3", 64, "GUA_C4", 65,
	"CYT_P", 66, "CYT_O1P", 67, "CYT_O2P", 68, "CYT_O5'", 69, "CYT_C5'", 70, "CYT_C4'", 71,
	"CYT_O4'", 72, "CYT_C3'", 73, "CYT_O3'", 74, "CYT_C2'", 75, "CYT_O2'", 76, "CYT_C1'", 77,
	"CYT_N1", 78, "CYT_C2", 79, "CYT_O2", 80, "CYT_N3", 81, "CYT_C4", 82, "CYT_N4", 83,
	"CYT_C5", 84, "CYT_C6", 85
);
@list = keys(%index);
$nindex = @list;
#define 4 nucleotide names
@base1 = ("A", "U", "G", "C");
@base3 = ("ADE", "URA", "GUA", "CYT");
#define number of atoms in nucleotides
@natm_base = (22, 20, 23, 20);
$nbase = @base1;
#check the commond line
if (@ARGV != 4) {
	printf STDERR "shiftseq_rna.pl <arg1> <arg2> <arg3> <arg4>\n";
	printf STDERR "<arg1>: input pdb file\n";
	printf STDERR "<arg2>: shift of residue number\n";
	printf STDERR "<arg3>: desired chain id, _ means blank\n";
	printf STDERR "<arg4>: output pdb file\n";
	exit;
}
#check model pdb file
if (not -f $ARGV[0]) {
	printf "shiftseq_rna.pl error: pdb file $ARGV[0] doesn't exist\n";
	exit;
}
else {
	#read model pdb
	$modpdb = $ARGV[0];
	&readpdb($modpdb);
}
#change residue number
for ($i = 0; $i < $natm; $i++) {
	$rseq[$i] = $rseq[$i] + $ARGV[1];
	if (substr($ARGV[2], 0, 1) ne '_') {
		$chid[$i] = substr($ARGV[2], 0, 1);
	}
	else {
		$chid[$i] = ' ';
	}
}
&writepdb($ARGV[3]);

#read in protein pdb file and store info in global arrays
sub readpdb() {
	my $ipdb, $jpdb, $pdblen, @pdbtxt;
	my $strtmp, @dual, $pdbnam, $find;
	open(PROPDB, "$_[0]");
	@pdbtxt = <PROPDB>;
	$pdblen = @pdbtxt;
	close(PROPDB);
	for ($natm = 0, $ipdb = 0; $ipdb < $pdblen; $ipdb++) {
		if ($pdbtxt[$ipdb] =~ /^ATOM/) {
			chomp $pdbtxt[$ipdb];
			$strtmp = substr($pdbtxt[$ipdb], 11, 10);
			@dual = split(/ +/, $strtmp);
			#convert op1 and op2 to o1p and 02p
			if ($dual[1] eq "OP1") {
				$dual[1] = "O1P";
			}
			if ($dual[1] eq "OP2") {
				$dual[1] = "O2P";
			}
			#1-letter base name
			if (length($dual[2]) == 1) {
				for ($jpdb = 0; $jpdb < $nbase; $jpdb++) {
					if ($dual[2] eq $base1[$jpdb]) {
						$dual[2] = $base3[$jpdb];
					}
				}
			}
			$pdbnam = $dual[2]."_".$dual[1];
			$find = 0;
			for ($jpdb = 0; $jpdb < $nindex; $jpdb++) {
				if ($pdbnam eq $list[$jpdb]) {
					$find = 1;
					last;
				}
			}
			next if ($find == 0);
			$anam[$natm] = $dual[1];
			$rnam[$natm] = $dual[2];
			$chid[$natm] = substr($pdbtxt[$ipdb], 21, 1);
			$rseq[$natm] = substr($pdbtxt[$ipdb], 22, 4);
			$xpdb[$natm] = substr($pdbtxt[$ipdb], 30, 8);
			$ypdb[$natm] = substr($pdbtxt[$ipdb], 38, 8);
			$zpdb[$natm] = substr($pdbtxt[$ipdb], 46, 8);
			$natm++;
		}
	}
	#number of residues
	$nres = 0;
	$rseq_prev = -1000;
	for ($ipdb = 0; $ipdb < $natm; $ipdb++) {
		if ($rseq[$ipdb] != $rseq_prev) {
			$rseq_prev = $rseq[$ipdb];
			$ratm_star[$nres] = $ipdb;
			$nres++;
		}
	}
}

#write a protein pdb file in standard format
sub writepdb() {
	my $ipdb;
	open(PROPDB, ">$_[0]");
	for ($ipdb = 0; $ipdb < $natm; $ipdb++) {
		printf PROPDB ("ATOM%7d  %-4s%-4s%1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n", $ipdb + 1, $anam[$ipdb], $rnam[$ipdb], $chid[$ipdb], $rseq[$ipdb], $xpdb[$ipdb], $ypdb[$ipdb], $zpdb[$ipdb], 1.00, 0.00);
	}
	close(PROPDB);
}
