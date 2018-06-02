#!/usr/bin/perl
use strict;
use warnings;

my @pj = qw(brain breast colorectal hnsc kidney lung ov prad thca ucec);
my $ascat_dir = "/Volumes/areca42TB/tcga/CNA";
open (OUT,">purity_list.tsv");
print OUT "patient_id\tpurity\n";
foreach my $pj (@pj){
		my $file = "$ascat_dir/$pj/cel/annotate_ascat.tsv.gz";
		(-e $file) or die "ERROR::not exist $file!!\n";
		open(IN,"gunzip -c $file|");
		my $patient_id="";
		<IN>;
		while(<IN>){
				chomp;
				my @line = split(/\t/,);
				if($line[1] eq $patient_id){next;}
				$patient_id = $line[1];
				print OUT "$patient_id\t$line[12]\n";
		}
		close IN;
}
close OUT;
exit;

