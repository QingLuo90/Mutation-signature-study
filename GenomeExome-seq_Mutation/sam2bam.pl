#!/usr/bin/perl


my $strSampleList = "sampleSheet.txt";
open(LIST, "<".$strSampleList) or die "Can not open $strSampleList.\r\n";
while(<LIST>){
	my $strLine = $_;
	$strLine =~ s/^\s+//;
	$strLine =~ s/\s+$//;	
	my $strSam = "./bam/".$strLine.".sam";
	my $strBam = "./bam/".$strLine.".bam";
        my $job = $strLine."_samtobam.erro";
	my $strCmd="bsub -n 6 'samtools import /work/home/luoqing/project/mm10/mm10.fa $strSam $strBam'";
	system($strCmd);
}
close(LIST);
