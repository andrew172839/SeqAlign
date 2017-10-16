#!/usr/bin/perl
#pudge dir
$seqid_bin = "~/seqid_seg.pl";
#check the command line
if (@ARGV != 5) {
	printf STDERR "calc_seqid.pl usuage: calc_seqid.pl <arg1> <arg2> <arg3> <arg4> <arg5>\n";
	printf STDERR "<arg1>: template list screen-copied from pudge output\n";
	printf STDERR "<arg2>: directory where the alignment files are stored\n";
	printf STDERR "<arg3>: pudge id, e.g., p05191\n";
	printf STDERR "<arg4>: first residue of segment\n";
	printf STDERR "<arg5>: second residue of segment\n";
	exit;
}
#check the list file
if (not -f $ARGV[0]) {
	printf STDERR "calc_seqid.pl error: can not find list file\n";
	exit;
}
open(LST, $ARGV[0]);
@tmptxt = <LST>;
$tmplen = @tmptxt;
close(LST);
#check the dir
if (not -d $ARGV[1]) {
	printf STDERR "calc_seqid.pl error: can not find alignment dir\n";
	exit;
}
#read the list and copy the alignment to local dir
for ($i = 0; $i < $tmplen; $i++) {
	chomp $tmptxt[$i];
	#get the basename
	@arr = split /\t+/, $tmptxt[$i];
	$basename = $arr[0];
	#get the alignment method
	@arr = split /\./, $basename;
	$almethod = $arr[1];
	if ($almethod eq "bl5") {
		$almethod = "nal";
	}
	#make up the full alignment name
	$algnname = $ARGV[2]."_on_".$basename.".".$almethod.".pir";
	#copy the file if it exists
	if (-f "$ARGV[1]/$algnname") {
		$line = `$seqid_bin $ARGV[1]/$algnname $ARGV[3] $ARGV[4]`;
		chomp $line;
		@arr = split /=/, $line;
		$sid = $arr[1];
		printf ("%s\t%8.3f\n", $tmptxt[$i], $sid);
	}
}
