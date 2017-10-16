#!/usr/bin/perl

#check the command-line arguments
if (@ARGV != 2) {
	printf "match_align.pl <arg1> <arg2>\n";
	printf "<arg1>: alignment of two sequences, in .pir format\n";
	printf "<arg2>: number of a particular residue in template\n";
	printf "output: corresponding residue in query sequence\n";
	exit;
}
#check the alignment file
if (not -f $ARGV[0]) {
	printf "match_align.pl error: cannot find the alignment file - $ARGV[0]\n";
	exit;
}
#check the file name
$idx = index($ARGV[0], '.pir');
if ($idx < 0) {
	printf "match_align.pl error: check the for mat of alignment file - $ARGV[0]\n";
	printf ".pir for mat is required\n";
	exit;
}
#read in the sequence alignment
&readpir($ARGV[0]);
#look for the corresponding residue
for ($idx1 = 0, $idx2 = 0, $i = 0; $i < $naln; $i++) {
	if ($align[0][$i] ne '-') {
		$idx1++;
	}
	if ($align[1][$i] ne '-') {
		$idx2++;
	}
	if ($idx1 == $ARGV[1]) {
		printf ("residue in template: %4d   %1s\n", $idx1, $align[0][$i]);
		printf ("residue in query: %4d %1s%1s%1s%1s%1s\n", $idx2, $align[1][$i - 2], $align[1][$i - 1], $align[1][$i], $align[1][$i + 1], $align[1][$i + 2]);
		exit;
	}
}

#read in .pir format sequence alignment file
sub readpir() {
	my $ipir, $jpir, $kpir, $iseq;
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
				#skip the structurex and sequence
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
}
