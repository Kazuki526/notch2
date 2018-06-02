#!/usr/bin/perl
use warnings;
use strict;

#check nkf was installed??
my $nkfpath=`which nkf`;chomp $nkfpath;
($nkfpath and -e $nkfpath) or die "ERROR:nkf was not installed. please do\nbrew install nkf\n";

#check top driver 105genes bed file exist?
my $bed="$ENV{HOME}/git/notch2/notch2.json";
(-e $bed) or die "ERROR:not exist bed file:$bed\n";

#check existing bamslicing.pl
my $download_pl="$ENV{HOME}/git/notch2/bamslicing.pl";
(-e $download_pl) or die "ERROR::download script $download_pl is not exist!!\n";
#check reference of bam are exist?
my $ref="/Volumes/areca42TB/GRCh38.d1.vd1.fa";
(-e $ref)or die "ERROR:not exist ref fasta:$ref\n";

mkdir "norm_bam";
mkdir "tumor_bam";
#`curl --request POST --header \"Content-Type: application/json\" --data \@all_patient_bam.json 'https://gdc-api.nci.nih.gov/files'|nkf -Lu >response.tsv`;
my $response = "response_20180603.tsv";

# rad response and make norm & tumor manifest
my %data = ();
my($norm_manifest,$tumor_manifest) = ("norm_bam/norm_manifest.tsv","tumor_bam/tumor_manifest.tsv");
open(RES,"$response")or die "ERROR::cannot open $response\n";
open(OUTN,">$norm_manifest");
open(OUTT,">$tumor_manifest");
print OUTN "id\tfilename\n";print OUTT "id\tfilename\n";
my @header = split(/\t/,<RES>);chomp @header;
my %header_num=();
for(my $i=0;@header>$i;$i++){
		if($header[$i] eq "file_name"){$header_num{'file'}=$i;
		}elsif($header[$i] eq "cases.0.submitter_id"){$header_num{'case_id'}=$i;
		}elsif($header[$i] eq "cases.0.samples.0.sample_type"){$header_num{'sample_type'}=$i;
		}elsif($header[$i] eq "id"){$header_num{'file_id'}=$i;
		}elsif($header[$i] eq "cases.0.project.project_id"){$header_num{'project'}=$i;
		}elsif($header[$i] eq "cases.0.primary_site"){$header_num{'primary_site'}=$i;
		}else{die "ERROR::$response have wrong header what is $header[$i]\n";}
}
while(<RES>){
		chomp;
		my @line=split(/\t/,);
		my $pid = $line[$header_num{case_id}];
		if($line[$header_num{sample_type}] eq "Primary Tumor"){
				if(defined $data{$pid}{tumor_bam}){
						$data{$pid}{tumor_bam}.=";tumor_bam/$line[$header_num{file}]";
				}else{
						$data{$pid}{tumor_bam}="tumor_bam/$line[$header_num{file}]";
				}
				print OUTT "$line[$header_num{file_id}]\t$line[$header_num{file}]\n";
		}else{
				if(defined $data{$pid}{norm_bam}){
						$data{$pid}{norm_bam}.=";norm_bam/$line[$header_num{file}]";
				}else{
						$data{$pid}{norm_bam}="norm_bam/$line[$header_num{file}]";
				}
				print OUTN "$line[$header_num{file_id}]\t$line[$header_num{file}]\n";
		}
		$data{$pid}{project} = $line[$header_num{project}];
		$data{$pid}{primary_site} = $line[$header_num{primary_site}];
}
close RES;
close OUTN;
close OUTT;

#bamslicing
system("perl $download_pl $norm_manifest $bed norm_bam");
system("perl $download_pl $tumor_manifest $bed tumor_bam");

# read download error file
my %dl_err =();
open(ERN,"norm_bam/result_download.txt");
while(<ERN>){
		chomp;
		if($_ =~ /^\d+:(.+) download more than 10 times so this file cannot download/){
				$dl_err{"norm_bam/$1"} = "error";
		}
}
close ERN;
open(ERT,"tumor_bam/result_download.txt");
while(<ERT>){
		chomp;
		if($_ =~ /^\d+:(.+) download more than 10 times so this file cannot download/){
				$dl_err{"tumor_bam/$1"} = "error";
		}
}
close ERT;

#read purity file
open(PL,"purity_list.tsv") or die "ERROR::not exist purity file:purity_list.tsv\n";
<PL>;
while(<PL>){
		chomp;
		my @line = split(/\t/,);
		$data{$line[0]}{purity} = $line[1];
}
close PL;

#merge bam files to 1bam file when 1patient have many bam
print "merging bam files\n";
foreach my $pid(keys %data){
		if((!defined $data{$pid}{tumor_bam})||(!defined $data{$pid}{norm_bam})){$data{$pid}{focal}="no";next;}else{$data{$pid}{focal}="ok";}
		if($data{$pid}{norm_bam} =~ /;/){
				my @bam=split(/;/,$data{$pid}{norm_bam});
				my $focal = 0;
				foreach my $bam (@bam){if(defined$dl_err{$bam}){$focal++;}}
				if($focal == scalar(@bam)){$data{$pid}{focal}="no";next;}
				`samtools merge -f norm_bam/$pid.bam @bam`;
				`samtools index norm_bam/$pid.bam`;
				$data{$pid}{'file_norm'}="norm_bam/$pid.bam";
		}else{
				if(defined$dl_err{$data{$pid}{norm_bam}}){$data{$pid}{focal}="no";}
				`samtools index $data{$pid}{norm_bam}`;
				$data{$pid}{'file_norm'}=$data{$pid}{norm_bam};
		}
		if($data{$pid}{tumor_bam} =~ /;/){
				my @bam=split(/;/,$data{$pid}{tumor_bam});
				my $focal = 0;
				foreach my $bam (@bam){if(defined$dl_err{$bam}){$focal++;}}
				if($focal == scalar(@bam)){$data{$pid}{focal}="no";next;}
				`samtools merge -f tumor_bam/$pid.bam @bam`;
				`samtools index tumor_bam/$pid.bam`;
				$data{$pid}{'file_tumor'}="tumor_bam/$pid.bam";
		}else{
				if(defined$dl_err{$data{$pid}{norm_bam}}){$data{$pid}{focal}="no";}
				`samtools index $data{$pid}{tumor_bam}`;
				$data{$pid}{'file_tumor'}=$data{$pid}{tumor_bam};
		}
}
# varscan
print "doing varscan\n";
open(FPL,">focal_patient_list.tsv");
print FPL "patient_id\tproject\tprimary_site\tnorm_bam\ttumor_bam\n";
mkdir "varsca_vcf";
foreach my $pid(keys %data){
		if((!defined$data{$pid}{focal})||($data{$pid}{focal} eq "no")){print "$pid cannot download both tumor & normal files\n";next;}
		print FPL "$pid\t$data{$pid}{project}\t$data{$pid}{primary_site}\t$data{$pid}{file_norm}\t$data{$pid}{file_tumor}\n";
		my $purity =1;
		if(defined$data{$pid}{purity}){$purity = $data{$pid}{purity};}
		my $mpile = "samtools mpileup -q 10 -f $ref $data{$pid}{file_norm} $data{$pid}{file_tumor}";
		`zsh -c \"varscan somatic <\($mpile\) varscan_vcf/$pid --tumor-purity $purity --p-value 0.1 --output-vcf 1 --mpileup 1\"`;
}
close FPL;

# make depth file
mkdir "tdepth";
mkdir "ndepth";
my @norm_bam=();
my @tumor_bam=();
my @patients=();
my $file_num=0;
foreach my $pid(keys %data){
		if((!defined$data{$pid}{focal})||($data{$pid}{focal} eq "no")){next;}
		push(@patients,$pid);
		push(@norm_bam,$data{$pid}{file_norm});
		push(@tumor_bam,$data{$pid}{file_tumor});
		if(scalar(@patients)==50){
				$file_num++;
				my($tdepth,$ndepth)=("tdepth/tdepth$file_num.tsv","ndepth/ndepth$file_num.tsv");
				open(OUT,">$tdepth");print OUT "chr\tposition\t",join("\t",@patients)."\n";close OUT;
				open(OUT,">$ndepth");print OUT "chr\tposition\t",join("\t",@patients)."\n";close OUT;
				print "make depth $file_num\n";
				`samtools depth -q 10 @tumor_bam >>$tdepth`;
				`samtools depth -q 10 @norm_bam >>$ndepth`;
				@patients=();
				@tumor_bam=();
				@norm_bam=();
		}
}
$file_num++;
my($tdepth,$ndepth)=("tdepth/tdepth$file_num.tsv","ndepth/ndepth$file_num.tsv");
open(OUT,">$tdepth");print OUT "chr\tposition\t",join("\t",@patients)."\n";close OUT;
open(OUT,">$ndepth");print OUT "chr\tposition\t",join("\t",@patients)."\n";close OUT;
`samtools depth -q 10 @tumor_bam >>$tdepth`;
`samtools depth -q 10 @norm_bam >>$ndepth`;

exit;
