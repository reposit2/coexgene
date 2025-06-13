#!/usr/local/bin/perl

# Description:
#   This script identifies genes that are likely to be regulated by a specified
#   nucleic acid-binding protein (NABP) based on their DeepLIFT-derived rank scores.
#   It selects genes where the NABP ranks within a specified top threshold,
#   indicating potential regulatory relevance.

# Inputs:
#   $infile  - Path to the NABP rank file (output from dlscore_ranksort.pl), containing ranked NABPs for each gene
#   $nabp    - Symbol name of the NABP of interest (e.g., AKAP8, PKM)
#   $rank    - Rank threshold: a gene is selected if the NABP ranks ≤ this value for that gene

# Outputs:
#   $outfile  - Tab-delimited file listing selected genes and their corresponding scores
#   $outfile2 - Plain text file containing only Ensembl gene IDs of selected genes

# Usage:
#   Update the variables below for your target NABP, cell type, and model timestamp.
#   Uncomment the appropriate configuration block.

# Example configurations:
#   AKAP8 rank 1 - K562 / Neuronal cells - model 2024-10-08_03-21-22
#   PKM rank 2   - K562 / Neuronal cells - model 2024-10-08_03-21-22
#                  HepG2 cells - model 2024-10-08_01-30-05

# Note:
#   - Rank scores are assumed to be sorted in descending order (higher scores = stronger importance).
#   - Genes are selected only if the NABP appears within the top $rank entries and has a non-zero score.

# AKAP8 rank 1 ; K562 cells : DNA_3.txt.gz ; model 2024-10-08_03-21-22
$infile = './3celltypes_k562/2024-10-08_03-21-22/DeepLIFT/DNA_3_rank.txt.gz';
$outfile = './geneidlist/K562_AKAP8_rank1_2024-10-08_03-21-22.txt';
$nabp = 'AKAP8';
$rank = 1;

# AKAP8 rank 1 ; neuronal cells : DNA_2.txt.gz ; model 2024-10-08_03-21-22
#$infile = './3celltypes_k562/2024-10-08_03-21-22/DeepLIFT/DNA_2_rank.txt.gz';
#$outfile = './geneidlist/neuronal_AKAP8_rank1_2024-10-08_03-21-22.txt';
#$nabp = 'AKAP8';
#$rank = 1;

# PKM rank 2 ; K562 cells ; PKM rank 2 ; model 2024-10-08_03-21-22
#$infile = './3celltypes_k562/2024-10-08_03-21-22/DeepLIFT/DNA_3_rank.txt.gz';
#$outfile = './geneidlist/K562_PKM_rank2_2024-10-08_03-21-22.txt';
#$nabp = 'PKM';
#$rank = 2;

# PKM rank 2 ; neuronal cells ; PKM rank 2 ; model 2024-10-08_03-21-22
#$infile = './3celltypes_k562/2024-10-08_03-21-22/DeepLIFT/DNA_3_rank.txt.gz';
#$outfile = './geneidlist/neuronal_PKM_rank2_2024-10-08_03-21-22.txt';
#$nabp = 'PKM';
#$rank = 2;


# For HepG2 cell analysis
# HepG2 cells : DNA_3.txt.gz ; model 2024-10-08_01-30-05
#$infile = "./3celltypes_hepg2/2024-10-08_01-30-05/DeepLIFT/DNA_3_rank.txt.gz";
# neuronal cells : DNA_2.txt.gz
#$infile = "./3celltypes_hepg2/2024-10-08_01-30-05/DeepLIFT/DNA_2_rank.txt.gz";


open(IN,"zcat $infile |");
while($l = <IN>) {
	chomp $l;
	@list = split(/\t/,$l);
	@list2 = split(/\,/,$list[1]);
	for($i = 0;$i < @list2;$i++) {			# Iterate in descending order (from high to low scores)
#	for($i = @list2 - 1;$i >= 0;$i--) {		# Iterate in ascending order (from low to high scores)
		@list3 = split(/\|/,$list2[$i]);
		if ($list3[1] != 0) {
#			$cou++;
			last if ($cou == $rank);
			$cou++;
#			next if ($cou <= 3);	 	# Skip the top 3 entries (1st, 2nd, and 3rd)
			if ($list3[0] eq $nabp) {
				$sc{$list[0]} = $list3[1];
				print "$list[0]\t$list3[0]\t$list3[1]\n";
				last;
			}
		}
	}
	$cou = 0;
}
close(IN);

$r = (keys %sc);
$outfile =~ s/\.txt/\_$r.txt/;
$outfile2 = $outfile;
$outfile2 =~ s/\.txt/\_gid\.txt/;
print "$outfile\n$outfile2\n";
#open(OUT,"| gzip -9c > $outfile");
open(OUT2,">$outfile2");
open(OUT,">$outfile");
foreach $x (sort {$sc{$b} <=> $sc{$a}} keys %sc) {
	print OUT "$x\t$sc{$x}\n";
	print OUT2 "$x\n";
}
close(OUT);
close(OUT2);

