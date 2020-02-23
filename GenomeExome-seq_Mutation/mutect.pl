#!/usr/bin/perl

my $strSampleList = "sampleSheet2.txt";
open(LIST, "<".$strSampleList) or die "Can not open $strSampleList.\r\n";
while(<LIST>){
	my $strLine = $_;
	$strLine =~ s/^\s+//;
	$strLine =~ s/\s+$//;	
      my @mouse = split(/_/,$strLine );
      my $mouseIndex = @mouse[0];
	my $input1 = "./".$strLine."_recal.bam.bam";
	my $input2 = "./".$mouseIndex ."_Tail.bam";
      my $output = "./".$strLine.".vcf";
	my $strCmd="bsub  -n 10 'java -jar /work/home/luoqing/software/GenomeAnalysisTK.jar -T MuTect2 -R /work/home/luoqing/project/mm10/mm10.fa -I:tumor $input1 -I:normal $input2 --dbsnp /work/home/luoqing/project/ZhaoningLu/snp.vcf -nct 80 -o $output'";
	system($strCmd);
}
close(LIST);


