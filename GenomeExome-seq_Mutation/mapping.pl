#!/usr/bin/perl

my $strSampleList = "sampleSheet.txt";
open(LIST, "<".$strSampleList) or die "Can not open $strSampleList.\r\n";
while(<LIST>){
	my $strLine = $_;
	$strLine =~ s/^\s+//;
	$strLine =~ s/\s+$//;	
	my $strR1 = "./".$strLine."_R1.fq.gz";
	my $strR2 = "./".$strLine."_R2.fq.gz";
	my $strSam = "./bam/".$strLine.".sam";
      my $job = 	$strLine.".erro";
	my $strCmd="bsub -e $job -n 4 'bwa mem -M -t 8 /work/home/luoqing/project/mm10/mm10.fa $strR1 $strR2 > $strSam' ";
	system($strCmd);
}
close(LIST);