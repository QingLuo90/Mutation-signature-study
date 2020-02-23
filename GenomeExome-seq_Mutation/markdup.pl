#!/usr/bin/perl

my $strSampleList = "sampleSheet.txt";
open(LIST, "<".$strSampleList) or die "Can not open $strSampleList.\r\n";
while(<LIST>){
	my $strLine = $_;
	$strLine =~ s/^\s+//;
	$strLine =~ s/\s+$//;	

my $input = "./bam/".$strLine."_addRG_sorted.bam";
my $output =  "./bam/".$strLine."_addRG_sorted_markdup.bam";
my $metrics = "./bam/".$strLine."_metrics.txt";
my $job = $strLine."_markdup.erro";

my $strCmd="bsub -e $job -n 6 'java -jar /work/home/luoqing/software/picard.jar MarkDuplicates INPUT=$input OUTPUT=$output METRICS_FILE=$metrics CREATE_INDEX=true'";
	system($strCmd);
}
close

