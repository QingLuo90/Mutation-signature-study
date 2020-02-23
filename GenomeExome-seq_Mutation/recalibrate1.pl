#!/usr/bin/perl

my $strSampleList = "sampleSheet.txt";
open(LIST, "<".$strSampleList) or die "Can not open $strSampleList.\r\n";
while(<LIST>){
	my $strLine = $_;
	$strLine =~ s/^\s+//;
	$strLine =~ s/\s+$//;	

my $input = "./bam/".$strLine."_addRG_sorted_markdup.bam";
my $output =  "./bam/".$strLine."_recal.table";
my $job = $strLine."_recalibrate1.erro";

my $strCmd="bsub -e $job -n 8 'java -jar /work/home/luoqing/software/GenomeAnalysisTK.jar -T BaseRecalibrator -R /work/home/luoqing/project/mm10/mm10.fa -I $input -knownSites /work/home/luoqing/project/ZhaoningLu/snp.vcf -knownSites /work/home/luoqing/project/ZhaoningLu/indel.vcf -o $output'";
	system($strCmd);
}
close

