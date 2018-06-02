#!/usr/bin/perl
use strict;
use warnings;
use Parallel::ForkManager;

############   perl bamslicing.pl  MANIFEST BED_JSON OUT_DIR

#input check?
my ($manifest,$json,$out_dir)=@ARGV;
if( $manifest !~/manifest/|$json !~/json$/){die "ERROR::input files are correct??\n";}
(-e $out_dir) or die "ERROR::out dir is not exist!\n";

#perl bamslicing.pl MANIFEST.tsv out_DIR maxparallel
my $token_path=`ls ~/git/innanlab/gdc|grep 'gdc-user-token'`;
if(!$token_path){die "!!ERROR!!:token file not exitst!!";}
chomp $token_path;
$token_path="~/git/innanlab/gdc/$token_path";
my $token=`cat $token_path`;

($manifest and -e "$manifest") or die "ERROR: there arenot manifest file!!\n";
($json     and -e "$json"    ) or die "ERROR: there arenot json file!!\n";

open(MAN,"$manifest");
my $linen=0;
my $max_processes=10;
my $pm = new Parallel::ForkManager($max_processes);
my $header=<MAN>;chomp $header;
my @header = split(/\t/,$header);
if(($header[0] ne "id")||($header[1] ne "filename")){die "ERROR::manifest file colume is correct??\n";}

open (ERR,">$out_dir/result_download.txt");
while(<MAN>){
		$linen++;
		#forks and returns the pid for child
		my $pid = $pm->start and next;

#		if($start_line_num > $line_num){next;
#		}elsif($end_line_num < $line_num){last;}
		chomp;
		my @line=split(/\t/,);
		my @ls=`ls $out_dir`;chomp @ls;
		if(!grep{$_ eq "$line[1]"}@ls){
				print "$linen:curl $line[1] now\n";
				print ERR "$linen:curl $line[1] now\n";
				system("curl --header \"X-Auth-Token: $token\" --request POST https://api.gdc.cancer.gov/slicing/view/$line[0] --header \"Content-Type: application/json\" -d\@$json --output $out_dir/$line[1] > /dev/null 2>&1");
		}
		my $focal=0;
		my $time=0;
		while($focal==0){
				$time++;
				my $sam_head = `samtools view $out_dir/$line[1] 2>&1|head -n 1`;
				if($time >10){
						print "$linen:$line[1] download more than 10 times so this file cannot download?\n";
						print ERR "$linen:$line[1] download more than 10 times so this file cannot download?\n";
						$focal++;
				}elsif(($sam_head !~ /EOF\smarker\sis\sabsent/) && ($sam_head !~ /Parse\serror/)){
						$focal++;
						print "$linen:$line[1] downloaded is ok\n";
						print ERR "$linen:$line[1] downloaded is ok\n";
				}else{
						print "$linen:$line[1] download error. download again\n";
						print ERR "$linen:$line[1] download error. download again\n";
				}
				if($focal==0){
						system("curl --header \"X-Auth-Token: $token\" --request POST https://api.gdc.cancer.gov/slicing/view/$line[0] --header \"Content-Type: application/json\" -d\@$json --output $out_dir/$line[1] > /dev/null 2>&1");
				}
		}

		$pm->finish; #terminates the child process
}
close MAN;
close ERR;
exit;
