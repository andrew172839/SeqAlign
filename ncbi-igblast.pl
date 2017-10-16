#!/usr/bin/perl

#igblast directory
$igblast_dir = "~/programs/IgBLAST/ncbi-igblast-1.0.0";
#igblast options
$igblast_opt1 = "-germline_db_V $igblast_dir/database/human_gl_V -germline_db_J $igblast_dir/database/human_gl_J -germline_db_D $igblast_dir/database/human_gl_D -organism human -domain_system kabat -auxiliary_data $igblast_dir/optional_file/human_gl.aux -num_alignments 1 -show_translation -outfmt 4";
$igblast_opt2 = "-germline_db_V $igblast_dir/database/human_gl_V -germline_db_J $igblast_dir/database/human_gl_J -germline_db_D $igblast_dir/database/human_gl_D -organism human -domain_system kabat -auxiliary_data $igblast_dir/optional_file/human_gl.aux -num_alignments 1 -show_translation -outfmt 7";
#check the command-line arguments
if (@ARGV ne 2) {
	printf STDERR "ncbi-igblast.pl <arg1> <arg2>\n";
	printf STDERR "<arg1>: input nucleotide sequence file\n";
	printf STDERR "<arg2>: igBLAST for mat option: 4 or 7\n";
	exit;
}
#check if the protein pdb exists
if (not -f $ARGV[0]) {
	printf STDERR "ncbi-igblast.pl error: $ARGV[0] doesn't exist\n";
	exit;
}
#run ncbi igblast
if ($ARGV[1] == 4) {
	system "$igblast_dir/bin/igblastn $igblast_opt1 -query $ARGV[0]";
}
elsif ($ARGV[1] == 7) {
	system "$igblast_dir/bin/igblastn $igblast_opt2 -query $ARGV[0]";
}
else {
	printf STDERR "ncbi-igblast.pl error: format option error\n";
}
