#!/usr/bin/perl

#$infile = "./ChatGPT_genefunction.md";
#$outfile = "./ChatGPT_genefunction_table.tsv";

$infile = "./ChatGPT_genefunction2b.md";
$outfile = "./ChatGPT_genefunction_table2b.tsv";

open(IN,"$infile");
@lines = <IN>;
close(IN);
$nlines = @lines;

open(OUT,">$outfile");
print OUT "Source\tGeneSetID\tGeneSetName\tGeneList\tn_Genes\tLLM Name\tLLM Analysis\tScore\n";
for($i = 0;$i < $nlines;$i++) {
	$l = $lines[$i];
	chomp $l;
	if ($l =~ /^\*\*Process\:.+\(([\d\.]+?)\)\*\*/) {
		$score = $1;
		$drid = '';
		$drid = $1 if ($prev2 =~ /\*\*(.+?)\*\*/);
	} elsif ($l =~ /^\*\*(\d+)\. (.+?)\*\*/) {
		$num = $1;
		$func = $2;
		$cluster = "Cluster\_$drid\_$num";
		$i++;
		$l2 = $lines[$i];
		undef @ita;
		$i++;
		$l3 = $lines[$i];
		chomp $l3;
		$j = $i;
		$l3t = '';
		while ($j != 0) {
#			while ($l3 =~ /[^\*]\*([^\*]+?)\*[^\*]/g) {
			while ($l3 =~ /(?<!\*)\*([^\*]+?)\*(?!\*)/g) {
				push @ita, $1;
			}
			$j2 = $j + 2;
			$l3t .= "$l3 ";
			if ($lines[$j2] !~ /^\*\*/) {
				$l3 = $lines[$j2];
				chomp $l3;
				$j = $j2;
			} else {
				$i = $j;
				$j = 0;
				$l3 = $l3t
			}
		}
		undef %gene;
		foreach $x (@ita) {
			@ita2 = split(/\, /,$x);
			foreach $y (@ita2) {
				$gene{$1} = 1 if ($y =~ /^([A-Z][A-Z0-9]+)/);
#					push @genes, $1;
			}
		}
		@genes = (sort keys %gene);
		$gn = @genes;
		$glist = join(' ',@genes);
		print OUT "NeST\t$cluster\t$cluster\t$glist\t$gn\t$func\t$l3\t$score\n";
		print "NeST\t$cluster\t$cluster\t$glist\t$gn\t$func\t$l3\t$score\n";
#		$aaa = <STDIN>;
	}
	$prev2 = $prev1;
	$prev1 = $l;
}
