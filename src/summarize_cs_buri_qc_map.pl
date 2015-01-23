#!/usr/bin/perl -w
$HOME=$ARGV[0];
$sourcedir=$ARGV[1];

@fastqclist =  ("Basic Statistics", "Per base sequence quality", "Per tile sequence quality", "Per sequence quality scores", "Per base sequence content", "Per sequence GC content", "Per base N content", "Sequence Length Distribution", "Sequence Duplication Levels", "Overrepresented sequences", "Adapter Content", "Kmer Content");

print "Run\tSampleName\tPercentMapped\tWolbachiaDepth\tWolbachiaBreadth\tPredictedInfectionStatus\t";
print "NumRead1\tNumBadRead1\tLengthRead1\tBasicStatisticsRead1\tPerBaseSequenceQualityRead1\tPerTileSequenceQualityRead1\tPerSequenceQualityScoresRead1\tPerBaseSequenceContentRead1\tPerSequenceGCContentRead1\tPerBaseNContentRead1\tSequenceLengthDistributionRead1\tSequenceDuplicationLevelsRead1\tOverrepresentedSequencesRead1\tAdapterContentRead1\tKmerContentRead1\t";
print "NumRead2\tNumBadRead2\tLengthRead2\tBasicStatisticsRead2\tPerBaseSequenceQualityRead2\tPerTileSequenceQualityRead2\tPerSequenceQualityScoresRead2\tPerBaseSequenceContentRead2\tPerSequenceGCContentRead2\tPerBaseNContentRead2\tSequenceLengthDistributionRead2\tSequenceDuplicationLevelsRead2\tOverrepresentedSequencesRead2\tAdapterContentRead2\tKmerContentRead2\t\n";

# set up list of sample names
@samples = ("TZ_GGACTCCT-TAGATCGC", "TP_TAGGCATG-TAGATCGC", "JC_AGGCAGAA-TAGATCGC", "HS_TCCTGAGC-TAGATCGC", "BS_CGTACTAG-TAGATCGC");

# submit fastQC, bowtie, bedtools and samtools jobs
# submit fastQC, bowtie, bedtools and samtools jobs
foreach $sample (@samples) {
	for ($lane=1; $lane<5; $lane++) {
		
		$run="$sample"."_L00"."$lane";

		# get % of mapped reads
		open FLAGSTAT, "$sourcedir/$run".".flagstat";
		while ($line = <FLAGSTAT>) {
			if ($line =~ /mapped\s\((\d+\.\d+)\%/) {
				$mapped=$1;
			}
		}

		# estimate depth and breadth of Wolbachia coverage
		open COVERAGE, "$sourcedir/$run".".coverage";
		$weightedSum=0;
		$zeroDepth=0;
		while ($line = <COVERAGE>) {
			if ($line =~ /chr\t(\d+)\t(\d+)/) {
				$depth=$1;
				$numSites=$2;
				$weightedSum=$weightedSum+($depth*$numSites);
				if ($depth == 0){
					$zeroDepth = $numSites/1267782;
				}
			}
		}
		
		# predict Wolbachia infection status
		$depth=$weightedSum/1267782;
		$breadth=1-$zeroDepth;
		$predictedStatus="uninfected";
		if (($breadth > 0.90) && ($depth > 1)) {
			$predictedStatus = "infected";
		}

		# get fastQC stats for read1 file
		open QCRUN1, "$sourcedir/$run"."_R1_001_fastqc/fastqc_data.txt";
		while ($line = <QCRUN1>) {
			if ($line =~ /Total Sequences\s(\d+)/) {
				$numseq1 = $1;
			}
			if ($line =~ /Sequences flagged as poor quality\s(\d+)/) {
				$numbad1 = $1;
			}
			if ($line =~ /Sequence length\s(\d+)/) {
				$length1 = $1;
			}
		}

		my %qc1;
		open SUMMARYRUN1, "$sourcedir/$run"."_R1_001_fastqc/summary.txt";
		while ($line = <SUMMARYRUN1>) {
			if ($line =~ /(\S+)\t(.+?)\t/) {
				$qc1{$2} = $1;
			}
		}

		# get fastQC stats for read2 file	
		open QCRUN2, "$sourcedir/$run"."_R2_001_fastqc/fastqc_data.txt";
        while ($line = <QCRUN2>) {
			if ($line =~ /Total Sequences\s(\d+)/) {
				$numseq2 = $1;
            }
			if ($line =~ /Sequences flagged as poor quality\s(\d+)/) {
				$numbad2 = $1;
            }
            if ($line =~ /Sequence length\s(\d+)/) {
            	$length2 = $1;
            }
        }

		my %qc2;
		open SUMMARYRUN2, "$sourcedir/$run"."_R2_001_fastqc/summary.txt";
		while ($line = <SUMMARYRUN2>) {
			if ($line =~ /(\S+)\t(.+?)\t/) {
				$qc2{$2} = $1;
			}
		}

		#print summary stats	
		print "$run\t$sample\t$mapped\t$depth\t$breadth\t$predictedStatus\t";
		print "$numseq1\t$numbad1\t$length1\t";
		foreach $stat (@fastqclist) {
			print "$qc1{$stat}\t";
		}

		print "$numseq2\t$numbad2\t$length2\t";
		foreach $stat (@fastqclist) {
			print "$qc2{$stat}\t";
		}
		
		print "\n";
	}
}
