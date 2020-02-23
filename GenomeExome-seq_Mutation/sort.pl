#!/usr/bin/perl

my $strSampleList = "sampleSheet.txt";
open(LIST, "<".$strSampleList) or die "Can not open $strSampleList.\r\n";
while(<LIST>){
	my $strLine = $_;
	$strLine =~ s/^\s+//;
	$strLine =~ s/\s+$//;	

my $input = "./bam/".$strLine."_addRG.bam";
my $output =  "./bam/".$strLine."_addRG_sorted.bam";
my $job = $strLine."_sort.erro";

my $strCmd="bsub -e $job -n 6 'java -jar /work/home/luoqing/software/picard.jar SortSam INPUT=$input OUTPUT=$output SORT_ORDER=coordinate CREATE_INDEX=true'";
	system($strCmd);
}
close

