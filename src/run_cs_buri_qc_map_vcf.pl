#!/usr/bin/perl -w
$HOME=$ARGV[0];
$sourcedir=$ARGV[1];

# get and unpack D. melanogaster reference sequence
`wget -q -P $sourcedir http://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.fa.gz`;
`gunzip -c $sourcedir/dm6.fa.gz > $sourcedir/dm6.fa`;

# get and reformat wolbachia reference sequence (change header to "chr" to be able to load in UCSC microbes browser)
`wget -P $sourcedir ftp://ftp.ncbi.nlm.nih.gov/genbank/genomes/Bacteria/Wolbachia_endosymbiont_of_Drosophila_melanogaster_uid272/AE017196.fna`;
`echo ">chr" > $sourcedir/w_pipientis.fa`;
`grep -v \\> $sourcedir/AE017196.fna >> $sourcedir/w_pipientis.fa`;

# make hologenome & clean up unused files
`cat $sourcedir/dm6.fa $sourcedir/w_pipientis.fa > $sourcedir/dm6_hologenome_raw.fa`;
`fasta_formatter -i $sourcedir/dm6_hologenome_raw.fa -o $sourcedir/dm6_hologenome.fa -w 50`;
`rm $sourcedir/dm6.fa $sourcedir/w_pipientis.fa $sourcedir/AE017196.fna $sourcedir/dm6_hologenome_raw.fa $sourcedir/dm6.fa.gz`;

# index reference hologenome
`bowtie2-build $sourcedir/dm6_hologenome.fa $sourcedir/dm6_hologenome.fa`;
`samtools faidx $sourcedir/dm6_hologenome.fa`;

# get sizes of hologenome sequences
`faCount $sourcedir/dm6_hologenome.fa | grep -v \\# | cut -f1,2 > $sourcedir/dm6_hologenome.sizes`;

# set up list of sample names
@samples = ("TZ_GGACTCCT-TAGATCGC", "TP_TAGGCATG-TAGATCGC", "JC_AGGCAGAA-TAGATCGC", "HS_TCCTGAGC-TAGATCGC");

# submit fastQC, bowtie, bedtools and samtools jobs
foreach $sample (@samples) {
	for ($lane=1; $lane<5; $lane++) {
		$run="$sample"."_L00"."$lane";
		$input1="$sample"."_L00"."$lane"."_R1_001.fastq.gz";
		$input1="$sample"."_L00"."$lane"."_R2_001.fastq.gz";
		system("qsub -l cores=6 -l mem=128 -V -b y -cwd -o $HOME/out -e $HOME/out -N fastq_1_$run \"fastqc --extract $sourcedir/$input1\"");	
		system("qsub -l cores=6 -l mem=128 -V -b y -cwd -o $HOME/out -e $HOME/out -N fastq_1_$run \"fastqc --extract $sourcedir/$input2\"");		
		system("qsub -l cores=6 -l mem=128 -V -b y -cwd -o $HOME/out -e $HOME/out -N map_$run \"bowtie2 -p 6 -x $sourcedir/dm6_hologenome.fa -1 $sourcedir/$input1 -2 $sourcedir/$input2 | samtools view -bS - | samtools sort - $sourcedir/$run\"");	
		$bamfile="$run".".bam";
		$coveragefile="$run".".coverage";
		system("qsub -l cores=6 -l mem=128 -V -b y -cwd -o $HOME/out -e $HOME/out -N bedtools_$run -hold_jid map_$run \"bedtools genomecov -ibam $sourcedir/$bamfile -g $sourcedir/dm6_hologenome.sizes > $sourcedir/$coveragefile\"");	
		$flagstatfile="$run".".flagstat";
		system("qsub -l cores=6 -l mem=128 -V -b y -cwd -o $HOME/out -e $HOME/out -N flagstat_$run -hold_jid map_$run \"samtools flagstat $sourcedir/$bamfile > $sourcedir/$flagstatfile\"");	
		$idxstatfile="$run".".idxstat";
		system("qsub -l cores=6 -l mem=128 -V -b y -cwd -o $HOME/out -e $HOME/out -N idx_$run -hold_jid map_$run \"samtools index $sourcedir/$bamfile\"");	
		system("qsub -l cores=6 -l mem=128 -V -b y -cwd -o $HOME/out -e $HOME/out -N idxstat_$run -hold_jid idx_$run \"samtools idxstats $sourcedir/$bamfile > $sourcedir/$idxstatfile\"");	
		$holdlist="map_$run".",";
	}
	$bamfiles="$sample"."_*.bam";
	system("qsub -l cores=6 -l mem=128 -V -b y -cwd -o $HOME/out -e $HOME/out -N merge_$sample -hold_jid $holdlist \"samtools merge $sourcedir/$sample.bam $sourcedir/$bamfiles\"");	 
	system("qsub -l cores=6 -l mem=128 -V -b y -cwd -o $HOME/out -e $HOME/out -N pileup_$sample -hold_jid merge_$sample \"samtools mpileup -d 100000 -ugf $sourcedir/dm6_hologenome.fa $sourcedir/$sample.bam | bcftools view -bvcg - > $sourcedir/$sample.bcf; bcftools view $sourcedir/$sample.bcf | vcfutils.pl varFilter -D 1000 > $sourcedir/$sample.vcf\"");	 
}