#!/usr/bin/perl
#check the command-line argumnets
if (@ARGV != 6) {
	printf "clustseq.pl <arg1> <arg2> <arg3> <arg4> <arg5> <arg6>\n";
	printf "<arg1>: input sequence file\n";
	printf "<arg2>: sequence identity cutoff, e.g., 0.75\n";
	printf "<arg3>: number of desired clusters, e.g., 10\n";
	printf "<arg4>: number of clustering cycles, e.g., 10\n";
	printf "<arg5>: scoring function tag, e.g., thrd1\n";
	printf "<arg6>: if reverse the score, e.g., 0: no; 1: yes\n";
	exit;
}
#check if the sequence file exists
if (not -f $ARGV[0]) {
	printf "error: sequence file - $ARGV[0] doesn't exist\n";
	exit;
}
#read in the sequences
open(INP, $ARGV[0]);
@inptxt = <INP>;
$inplen = @inptxt;
close(INP);
for ($nseq = 0, $i = 0; $i < $inplen; $i++) {
	if ($inptxt[$i] =~ /^>/) {
		chomp $inptxt[$i];
		$inpcom[$nseq] = $inptxt[$i];
		for ($seqtmp = "", $j = $i + 1; $j < $inplen; $j++) {
			last if ($inptxt[$j] =~ /^>/);
			chomp $inptxt[$j];
			$seqtmp = $seqtmp.$inptxt[$j];
		}
		$inpseq[$nseq] = $seqtmp;
		$nseq++;
	}
}
$nseg = $nseq;
#check sequence identity cutoff
if ($ARGV[1] > 1.0) {
	printf "error: sequence identity cutoff must be < 1.0\n";
	exit;
}
$cutoff = 0.1 / $ARGV[1];
#check number of clusters
if ($ARGV[2] > 500) {
	printf "error: too many clusters - $ARGV[2]\n";
	exit;
}
$nsite = $ARGV[2];
#check number of clustering cycles
if ($ARGV[3] > 20) {
	printf "error: too many cycles - $ARGV[3]\n";
	exit;
}
$ncycle = $ARGV[3];
#must be no less than 3 cycles
if ($ncycle < 3) {
	$ncycle = 3;
}
#pull out segment scores
for ($i = 0; $i < $nseg; $i++) {
	$strlen = length($inpcom[$i]);
	$strtmp = substr($inpcom[$i], 1, $strlen - 1);
	@arr1 = split / +/, $strtmp;
	$narr1 = @arr1;
	#search for the scoring function tag
	for ($j = 0; $j < $narr1; $j++) {
		if ($arr1[$j] =~ /$ARGV[4]=/) {
			@arr2 = split /=/, $arr1[$j];
			$narr2 = @arr2;
			$score = $arr2[$narr2 - 1];
			last;
		}
	}
	#search for the cdr-h3 / cdr-l3
	for ($j = 0; $j < $narr1; $j++) {
		if ($arr1[$j] =~ /cdrh3_seq=/ or $arr1[$j] =~ /cdrl3_seq=/) {
			@arr2 = split /=/, $arr1[$j];
			$narr2 = @arr2;
			$cdrseq[$i] = $arr2[$narr2 - 1];
			last;
		}
	}
	#reverse the score based on user input
	if ($ARGV[5] == 0) {
		$seg_scr[$i] = $score;
	}
	else {
		$seg_scr[$i] = 1.0 / $score;
	}
}
#perfor m the energy-driven clustering algorithm
&cluster();
#print out the representatives
open(CLS, ">clustseq.out");
for ($i = 0; $i < $nsite; $i++) {
	next if ($topptr[$i] < 0);
	printf CLS "cluster $i: $topnum[$i]\n";
	printf CLS "$inpcom[$topptr[$i]]\n";
	printf CLS "$inpseq[$topptr[$i]]\n";
}
close(CLS);
#output sequences
open(NEW, ">clustseq.fasta");
for ($i = 0; $i < $nsite; $i++) {
	for ($j = 0; $j < $nseg; $j++) {
		#check the cluster label
		next if ($seg_in_grp[$j] != $i);
		#output the sequence with cluster label
		printf NEW ("$inpcom[$j] clst = %d\n", $seg_in_grp[$j]);
		$strlen = length($inpseq[$j]);
		for ($k = 0; $k < $strlen; $k++) {
			printf NEW ("%1s", substr($inpseq[$j], $k, 1));
			if (int(($k + 1) / 60.0) * 60 == $k + 1) {
				printf NEW ("\n");
			}
		}
		if (int($strlen / 60.0) * 60 != $strlen) {
			printf NEW ("\n");
		}
	}
}
close(NEW);
#output blastclust-format file
open(NEW, ">bclust.out");
for ($i = 0; $i < $nsite; $i++) {
	for ($j = 0; $j < $nseg; $j++) {
		#check the cluster label
		next if ($seg_in_grp[$j] != $i);
		#search for the sequence index
		$strlen = length($inpcom[$j]);
		$strtmp = substr($inpcom[$j], 1, $strlen - 1);
		@arr1 = split / +/, $strtmp;
		$narr1 = @arr1;
		#search for the scoring function tag
		for ($k = 0; $k < $narr1; $k++) {
			if ($arr1[$k] =~ /index=/) {
				@arr2 = split /=/, $arr1[$k];
				$narr2 = @arr2;
				$seqidx = $arr2[$narr2 - 1];
				last;
			}
		}
		printf NEW ("%d ", $seqidx);
	}
	printf NEW ("\n");
}
close(NEW);

#divide segments into clusters
sub cluster() {
	my ($i, $j, $k, $l, $itmp, $imax, $enemax);
	my ($flag1, $flag2, $flag3);
	my ($replaced, $nreplace, $ireplace);
	my (@cls_scr, @cls_rep, @cls_rep_pre);
	my ($nclust_cur, $nclust_pre);
	my (@seg_bel);
	#set the cutoff for clustering
	#initialization before clustering
	for ($i = 0; $i < $nsite; $i++) {
		#energy of each cluster representative
		$cls_scr[$i] = 10000.0;
		#segment index of each cluster representative
		$cls_rep[$i] = -1;
		$cls_rep_pre[$i] = -1;
	}
	$nclust_pre = $nclust_cur = -1;
	#to which cluster each segment belongs
	for ($i = 0; $i < $nseg; $i++) {
		$seg_bel[$i] = -1;
	}
	$imax = 0;
	$enemax = 10000.0;
	#iterate until convergence
	LABEL0: for ($i = 0; $i < $ncycle; $i++) {
		printf ("cluster cycle: %3d\n", $i + 1);
		#visit all segments in the pool
		LABEL1: for ($j = 0; $j < $nseg; $j++) {
			printf ("processing %5d\r", $j + 1);
			#assign a segment to a cluster
			$replaced = 0;
			$nreplace = 0;
			for ($k = 0; $k < $nsite; $k++) {
				$flag1 = 0;
				if ($cls_rep[$k] == -1) {
					$flag1 = 1;
				}
				if ($flag1 == 0) {
					$rmsd = &getrmsd3($j, $cls_rep[$k]);
				}
				if ($flag1 == 1) {
					$rmsd = 5000.0;
				}
				#if rmsd less than the cutoff
				if ($rmsd < $cutoff) {
					#if energy is less favorite, go back
					if ($seg_scr[$j] > $cls_scr[$k]) {
						$seg_bel[$j] = $k;
						next LABEL1;
					}
					#if energy is better, mark it as will-be-replaced
					$ireplace[$nreplace] = $k;
					$nreplace++;
				}
			}
			#always replace the first in the should-be-replaced list
			if ($nreplace > 0) {
				$cls_scr[$ireplace[0]] = $seg_scr[$j];
				$cls_rep[$ireplace[0]] = $j;
				$seg_bel[$j] = $ireplace[0];
				$itmp = $ireplace[0];
				#change the flag here, the replacement done
				$replaced = 1;
			}
			#merge other clusters into the first in the replaced list
			for ($k = 1; $k < $nreplace; $k++) {
				$cls_scr[$ireplace[$k]] = 10000.0;
				$cls_rep[$ireplace[$k]] = -1;
			}
			#locate the empty-seat (max-ene one) for the next time fill-in
			for ($k = 1; $k<$nreplace; $k++) {
				#if one in the will-be-replaced list happens to be the max-ene one, but it has been merged into ireplace[0], use the cluster (ireplace[0]) as the max-ene one and then test if it is true by going through the whole cluster list
				if ($ireplace[$k] == $imax) {
					$enemax = $seg_scr[$j];
					$imax = $ireplace[0];
					for ($l = 0; $l < $nsite; $l++) {
						if ($enemax < $cls_scr[$l]) {
							$enemax = $cls_scr[$l];
							$imax = $l;
						}
					}
				}
			}
			#a new cluster should be filled into the empty-seat in the cluster list (marked with max-ene) and then relocate the max-ene one by going through the whole cluster list
			if ($replaced == 0) {
				$cls_scr[$imax] = $seg_scr[$j];
				$cls_rep[$imax] = $j;
				$seg_bel[$j] = $imax;
				$enemax = $seg_scr[$j];
				$itmp = $imax;
				for ($k = 0; $k < $nsite; $k++) {
					if ($enemax < $cls_scr[$k]) {
						$enemax = $cls_scr[$k];
						$imax = $k;
					}
				}
			}
		}
		#convergence criteria (0: stop; 1: continue)
		$flag3 = $flag4 = 0;
		#check if cluster centers remain the same
		for ($j = 0; $j < $nsite; $j++) {
			if ($cls_rep[$j] != $cls_rep_pre[$j]) {
				$flag3 = 1;
				last;
			}
		}
		#calculate the number of actual clusters
		for ($nclust_cur = 0, $j = 0; $j < $nsite; $j++) {
			if ($cls_rep[$j] != -1) {
				$nclust_cur++;
			}
		}
		if ($nclust_pre != $nclust_cur) {
			$flag4 = 1;
		}
		#printf out number of actual clusters
		printf ("found: %3d\n", $nclust_cur);
		printf ("flag3: %3d  flag4: %3d\n", $flag3, $flag4);
		#determine whether to terminate the clustering
		last LABEL0 if ($i > 1 and ($flag3 == 0 or $flag4 == 0));
		#update criteria
		for ($j = 0; $j < $nsite; $j++) {
			$cls_rep_pre[$j] = $cls_rep[$j];
		}
		$nclust_pre = $nclust_cur;
	}
	#sort the cluster representatives
	&sort2($nsite, @cls_scr);
	#pull out the top clusters after sorting
	for ($i = 0; $i < $nsite; $i++) {
		$topptr[$i] = $cls_rep[$bakptr[$i]];
		$topnum[$i] = 0;
		for ($j = 0; $j < $nseg; $j++) {
			if ($seg_bel[$j] == $bakptr[$i]) {
				$topnum[$i]++;
			}
		}
	}
	#assign each segment to a cluster(group)
	for ($i = 0; $i < $nseg; $i++) {
		$seg_in_grp[$i] = -1;
		for ($j = 0; $j < $nsite; $j++) {
			if ($seg_bel[$i] == $bakptr[$j]) {
				$seg_in_grp[$i] = $j;
			}
		}
	}
}

#calculate segment pairwise rmsd
sub getrmsd2() {
	my ($i, $j, $nwant1, $nwant2, $seqid, $rmsd);
	$nwant1 = $_[0];
	$nwant2 = $_[1];
	#initialization
	$rmsd = 0.0;
	#go back if comparing the same sequence
	if ($nwant1 == $nwant2) {
		return $rmsd;
	}
	#write seq file 1
	open(SEQ1, ">s$nwant1.fa");
	printf SEQ1 "$inpcom[$nwant1]\n";
	printf SEQ1 "$inpseq[$nwant1]\n";
	close(SEQ1);
	#write seq file 2
	open(SEQ2, ">s$nwant2.fa");
	printf SEQ2 "$inpcom[$nwant2]\n";
	printf SEQ2 "$inpseq[$nwant2]\n";
	close(SEQ2);
	#align two sequences
	system "clustalw.pl s$nwant2.fa s$nwant1.fa";
	#calculate sequence identity
	$output = `seqid.pl s$nwant1\_on_s$nwant2.pir`;
	chomp $output;
	@arr = split /=/, $output;
	$narr = @arr;
	$seqid = $arr[1] * 0.01;
	#set RMSD to be inversely weighted seqid
	$rmsd = 0.1 / $seqid;
	#clean up the scene
	unlink "s$nwant1.fa", "s$nwant2.fa";
	unlink "s$nwant1\_on_s$nwant2.pir";
	#return the rmsd value
	return $rmsd;
}

#calculate segment pairwise rmsd
sub getrmsd3() {
	my ($i, $j, $nwant1, $nwant2, $seqid, $rmsd);
	$nwant1 = $_[0];
	$nwant2 = $_[1];
	#initialization
	$rmsd = 0.0;
	#go back if comparing the same sequence
	if ($nwant1 == $nwant2) {
		return $rmsd;
	}
	#write seq file 1
	open(SEQ1, ">s$nwant1.fa");
	printf SEQ1 "$inpcom[$nwant1]\n";
	printf SEQ1 "$cdrseq[$nwant1]\n";
	close(SEQ1);
	#write seq file 2
	open(SEQ2, ">s$nwant2.fa");
	printf SEQ2 "$inpcom[$nwant2]\n";
	printf SEQ2 "$cdrseq[$nwant2]\n";
	close(SEQ2);
	#align sequence 1 onto 2
	system "clustalw.pl s$nwant2.fa s$nwant1.fa";
	#calculate sequence identity
	$output = `seqid.pl s$nwant1\_on_s$nwant2.pir`;
	chomp $output;
	@arr = split /=/, $output;
	$narr = @arr;
	$seqid = $arr[1] * 0.01;
	#align sequence 2 onto 1
	system "clustalw.pl s$nwant1.fa s$nwant2.fa";
	#calculate sequence identity
	$output = `seqid.pl s$nwant2\_on_s$nwant1.pir`;
	chomp $output;
	@arr = split /=/, $output;
	$narr = @arr;
	$seqid = ($seqid + $arr[1] * 0.01) / 2.0;
	#set RMSD to be inversely weighted seqid
	$rmsd = 0.1 / $seqid;
	#clean up the scene
	unlink "s$nwant1.fa", "s$nwant2.fa";
	unlink "s$nwant1\_on_s$nwant2.pir";
	unlink "s$nwant2\_on_s$nwant1.pir";
	#return the RMSD value
	return $rmsd;
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
			if ($index <=  0) {
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
				$j = $index  +  1;
			}
		}
		$list[$i] = $lists;
		$key[$i] = $keys;
	}
	@bakptr = @key;
	return;
}
