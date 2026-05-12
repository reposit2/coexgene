#!/usr/local/bin/perl

# Description:
#   This script normalizes DeepLIFT scores for each gene in a given cell type,
#   and ranks nucleic acid-binding proteins (NABPs) based on the normalized scores.

# Inputs:
#   $infile1  - Path to the input DeepLIFT score matrix (genes × NABPs)
#   $infile2  - Path to the file containing NABP names (used for filtering and labeling)

# Outputs:
#   $outfile1 - Path to the output file containing normalized DeepLIFT scores
#   $outfile2 - Path to the output file containing ranked NABPs for each gene

# Usage:
#   Modify the input/output file paths below according to the target cell type and model timestamp.
#   Uncomment the appropriate block for the desired cell type.

# Example cell types and model versions:
#   K562     - model 2024-10-08_03-21-22 (DNA_3.txt.gz)
#   HepG2    - model 2024-10-08_01-30-05 (DNA_3.txt.gz)
#   Neuronal - model 2024-10-08_03-21-22 2024-10-08_01-30-05 (DNA_2.txt.gz)

# K562 cells : DNA_3.txt.gz ; model 2024-10-08_03-21-22
$infile1 = "./3celltypes_k562/2024-10-08_03-21-22/DeepLIFT/DNA_3.txt.gz";
$infile2 = "./3celltypes_k562/2024-10-08_03-21-22/DeepLIFT/promoter_importance_mean.txt";
$outfile1 = "./3celltypes_k562/2024-10-08_03-21-22/DeepLIFT/DNA_3_norm.txt.gz";
$outfile2 = "./3celltypes_k562/2024-10-08_03-21-22/DeepLIFT/DNA_3_rank.txt.gz";

# neuronal cells : DNA_2.txt.gz ; model 2024-10-08_03-21-22
# (uncomment to use)
#$infile1 = "./3celltypes_k562/2024-10-08_03-21-22/DeepLIFT/DNA_2.txt.gz";
#$infile2 = "./3celltypes_k562/2024-10-08_03-21-22/DeepLIFT/promoter_importance_mean.txt";
#$outfile1 = "./3celltypes_k562/2024-10-08_03-21-22/DeepLIFT/DNA_2_norm.txt.gz";
#$outfile2 = "./3celltypes_k562/2024-10-08_03-21-22/DeepLIFT/DNA_2_rank.txt.gz";

# HepG2 cells : DNA_3.txt.gz ; model 2024-10-08_01-30-05
# (uncomment to use)
#$infile1 = "./3celltypes_hepg2/2024-10-08_01-30-05/DeepLIFT/DNA_3.txt.gz";
#$infile2 = "./3celltypes_hepg2/2024-10-08_01-30-05/DeepLIFT/promoter_importance_mean.txt";
#$outfile1 = "./3celltypes_hepg2/2024-10-08_01-30-05/DeepLIFT/DNA_3_norm.txt.gz";
#$outfile2 = "./3celltypes_hepg2/2024-10-08_01-30-05/DeepLIFT/DNA_3_rank.txt.gz";

# neuronal cells : DNA_2.txt.gz ; model 2024-10-08_01-30-05
# (uncomment to use)
#$infile1 = "./3celltypes_hepg2/2024-10-08_01-30-05/DeepLIFT/DNA_2.txt.gz";
#$infile2 = "./3celltypes_hepg2/2024-10-08_01-30-05/DeepLIFT/promoter_importance_mean.txt";
#$outfile1 = "./3celltypes_hepg2/2024-10-08_01-30-05/DeepLIFT/DNA_2_norm.txt.gz";
#$outfile2 = "./3celltypes_hepg2/2024-10-08_01-30-05/DeepLIFT/DNA_2_rank.txt.gz";

open(IN,"$infile2");
while($l = <IN>) {
        chomp $l;
        @list = split(/\,/,$l);
        next if ($list[1] !~ /promoter_annot/);
        $tid[$list[2]] = (split(/[\_\.]/,$list[1]))[2];
}
close(IN);

for($i = 1;$i <= 1310;$i++) {
        if ($header eq '') {
                $header = $tid[$i];
        } else {
                $header .= "\,$tid[$i]";
        }
}

open(IN,"zcat $infile1 |");
while($l = <IN>) {
    chomp $l;
    @list = split(/\t/,$l);
	@list2 = split(/\,/,$list[2]);
	shift @list2;
    $nlist2 = @list2;
    $sum = 0;
    for($i = 0;$i < $nlist2;$i++) {
        $sum += abs($list2[$i]);
    }
	if ($sum == 0) {
		$c++;
		next;
	}
    for($i = 0;$i < $nlist2;$i++) {
		$t{$list[1]} .= ',' if ($t{$list[1]} ne '');
        $t{$list[1]} .= ($list2[$i] / $sum);
    }
}
close(IN);

open(OUT,"| gzip -9c >$outfile1");
print OUT "$header\n";
foreach $x (sort keys %t) {
	print OUT "$x\t$t{$x}\n";
}
close(OUT); 

#print "$c\n";


open(OUT,"| gzip -9c > $outfile2");
open(IN,"zcat $outfile1 |");
$l = <IN>;
chomp $l;
@list = split(/\t/,$l);
@tfsym = split(/\,/,$list[1]);
while($l = <IN>) {
        chomp $l;
        @list = split(/\t/,$l);
        @list2 = split(/\,/,$list[1]);
        for($i = 0;$i < @list2;$i++) {
                $sc{$i} = $list2[$i];
        }
        $sc2 = '';
        foreach $x (sort {$sc{$b} <=> $sc{$a}} keys %sc) {
                $sc2 .= ',' if ($sc2 ne '');
                $sc2 .= "$tfsym[$x]\|$sc{$x}";
        }
        print OUT "$list[0]\t$sc2\n";
        undef %sc;
}
close(IN);
close(OUT);

