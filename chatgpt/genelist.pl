#!/usr/bin/perl

#$infile = '/Volumes/Mac20220720/work/学会勉強会/RECOMB2025/mac_work_coexpress/PANTHERin/gid2/full_model_20240215_5b_k_20240812_20241008_rankl_PAX5_5_324_gid.txt';
#$infile = '/Volumes/Mac20220720/work/学会勉強会/RECOMB2025/mac_work_coexpress/PANTHERin/gid2/full_model_20240215_5b_k_20240812_20241008_rankl_AKAP8_1_gid.txt';
#$infile = '/Users/user/work/coexpress/PANTHERin/gid2/full_model_20240215_5b_k_20240812_20241008_rankl_PKM_2_332_gid.txt';
#$infile = '/Users/user/work/coexpress/PANTHERin/gid2/full_model_20240215_4b_hepg2_20241007_rankl_ZKSCAN5_l5wo12_gid.txt';
#$infile = '/Users/user/work/coexpress/PANTHERin/gid2/full_model_20240215_4b_hepg2_20241007_rankl_ZBTB34_3.txt.gz';
#$infile = '/Users/user/work/coexpress/PANTHERin/gid2/full_model_20240215_4b_hepg2_20241007_rankl_ZNF557_l2_831_gid.txt';
#$infile = '/Users/user/work/coexpress/PANTHERin/gid2/full_model_20240215_4b_hepg2_20241007_rankl_ZNF446_l5_336_gid.txt';
#$infile = '/Users/user/work/coexpress/PANTHERin/gid2/full_model_20240215_4b_hepg2_20241007_rankl_ZC3H4_l3_gid.txt';
#$infile = '/Users/user/work/coexpress/PANTHERin/gid2/full_model_20240215_5b_k_20240812_20241008_rankl_POLR2G_1.txt.gz';
#$infile = '/Users/user/work/coexpress/PANTHERin/gid2/full_model_20240215_5b_k_20240812_20241008_rankl_RPLP0_8_73_gid.txt';
#$infile = '/Users/user/work/coexpress/PANTHERin/gid2/full_model_20240215_4b_hepg2_20241007_rankl_SPEN_l5_162_gid.txt';
# 20250509
#$infile = '/Users/user/work/coexpress/PANTHERin/gid2/full_model_20240215_5b_k_20240812_20241008_DNA2_rankl_PKM_2_326_gid.txt';
$infile = 'PANTHERin/gid2/full_model_20240215_5b_k_20240812_20241008_rankl_PKM_2_332_gid.txt';

$infile2 = './gencode.v19.annotation.gtf.gz';

#$outfile = '/Volumes/Mac20220720/work/学会勉強会/RECOMB2025/mac_work_coexpress/PANTHERin/gid2/full_model_20240215_5b_k_20240812_20241008_rankl_PAX5_5_324_symbol2.txt';
#$outfile = '/Volumes/Mac20220720/work/学会勉強会/RECOMB2025/mac_work_coexpress/PANTHERin/gid2/full_model_20240215_5b_k_20240812_20241008_rankl_AKAP8_1_symbol.txt';
#$outfile = '/Users/user/work/coexpress/PANTHERin/gid2/full_model_20240215_5b_k_20240812_20241008_rankl_PKM_2_332_symbol.txt';
#$outfile = '/Users/user/work/coexpress/PANTHERin/gid2/full_model_20240215_4b_hepg2_20241007_rankl_ZKSCAN5_l5wo12_symbol.txt';
#$outfile = '/Users/user/work/coexpress/PANTHERin/gid2/full_model_20240215_4b_hepg2_20241007_rankl_ZBTB34_3_symbol.txt';
#$outfile = '/Users/user/work/coexpress/PANTHERin/gid2/full_model_20240215_4b_hepg2_20241007_rankl_ZNF557_l2_831_symbol.txt';
#$outfile = '/Users/user/work/coexpress/PANTHERin/gid2/full_model_20240215_4b_hepg2_20241007_rankl_ZNF446_l5_336_symbol.txt';
#$outfile = '/Users/user/work/coexpress/PANTHERin/gid2/full_model_20240215_4b_hepg2_20241007_rankl_ZC3H4_l3_336_symbol.txt';
#$outfile = '/Users/user/work/coexpress/PANTHERin/gid2/full_model_20240215_4b_hepg2_20241007_rankl_POLR2G_1_233_symbol.txt';
#$outfile = '/Users/user/work/coexpress/PANTHERin/gid2/full_model_20240215_4b_hepg2_20241007_rankl_RPLP0_8_73_symbol.txt';
#$outfile = '/Users/user/work/coexpress/PANTHERin/gid2/full_model_20240215_4b_hepg2_20241007_rankl_SPEN_l5_162_symbol.txt';
# 20250509
#$outfile = '/Users/user/work/coexpress/PANTHERin/gid2/full_model_20240215_5b_k_20240812_20241008_DNA2_rankl_PKM_2_326_symbol.txt';
$outfile = 'PANTHERin/gid2/full_model_20240215_5b_k_20240812_20241008_rankl_PKM_2_332_symbol2.txt';

print "$infile\n$outfile\n";

open(IN,"gzip -dc $infile2 |");
while($l = <IN>) {
	chomp $l;
	@list = split(/\t/,$l);
	if ($list[8] =~ /gene_id \"(.+?)\..+gene_name \"(.+?)\"/) {
		$gid = $1;
		$name = $2;
		$gids{$gid} = $name;
	}
}
close(IN);

$n = 0;
open(OUT,">$outfile");
($infile =~ /gz$/) ? open(IN,"gzip -dc $infile |") : open(IN,"$infile");
while($l = <IN>) {
	chomp $l;
	$l = (split(/\t/,$l))[0];
	print "\," if ($n != 0);
	print OUT "\," if ($n != 0);
	print "$gids{$l}" if ($gids{$l} ne '');
	print OUT "$gids{$l}" if ($gids{$l} ne '');
	$n++;
}
close(IN);
close(OUT);
print "\n";
print "$n\n";
