#!/usr/bin/perl

my $strSampleList = "sampleSheet.txt";
open(LIST, "<".$strSampleList) or die "Can not open $strSampleList.\r\n";
while(<LIST>){
	my $strLine = $_;
	$strLine =~ s/^\s+//;
	$strLine =~ s/\s+$//;	

my $input = "./bam/".$strLine.".bam";
my $output =  "./bam/".$strLine."_addRG.bam";
my $RGID= "H0164.".$strLine ;
my $RGLB= "LB_".$strLine ;
my $RGSM= $strLine ;
my $job = $strLine."_addhead.erro";

my $strCmd="bsub -e $job -n 6 'java -jar /work/home/luoqing/software/picard.jar AddOrReplaceReadGroups INPUT=$input OUTPUT=$output RGID=$RGID RGLB=RGLB RGPL=illumina RGPU=unit2 RGSM=$RGSM'";
	system($strCmd);
}
close

