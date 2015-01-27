#!/usr/bin/perl -w
$inputdir=$ARGV[0];
$outputdir=$ARGV[1];
$paired_fastq_file="$outputdir"."/cs_buri_paired_fastq_sheet.tsv";
$fofn="$outputdir"."/cs_buri.fofn";

$instrument_model="Illumina HiSeq 2500";
$library_source="GENOMIC";
$library_selection="RANDOM";
$library_strategy="WGS";
$design_description="genome resequencing";
$library_construction_protocol="Nextera DNA Sample Preparation Kit (Product No. FC-121-1030)";
$insert_size="200";

# set up list of sample names
@samples = ("TZ_GGACTCCT-TAGATCGC", "TP_TAGGCATG-TAGATCGC", "JC_AGGCAGAA-TAGATCGC", "HS_TCCTGAGC-TAGATCGC", "BS_CGTACTAG-TAGATCGC");

open (CSV, ">$paired_fastq_file");
open (FOFN, ">$fofn");

print CSV "sample_alias\tinstrument_model\tlibrary_name\tlibrary_source\tlibrary_selection\tlibrary_strategy\tdesign_description\tlibrary_construction_protocol\tinsert_size\tforward_file_name\tforward_file_md5\treverse_file_name\treverse_file_md5\n";
# make md5 sums and print to .csv and .fofn files
foreach $sample (@samples) {
	for ($lane=1; $lane<5; $lane++) {
	
		$sample_alias = substr $sample, 0, 2;

		$forward_file_name="$sample"."_L00"."$lane"."_R1_001.fastq.gz";
		system("md5sum $inputdir/$forward_file_name > $outputdir/$forward_file_name.md5");	
		open FORMD5, "$outputdir/$forward_file_name.md5";
		$line = <FORMD5>;
		@elements = split (/\s+/, $line);
		$forward_file_md5 = $elements[0];

		$reverse_file_name="$sample"."_L00"."$lane"."_R2_001.fastq.gz";
		system("md5sum $inputdir/$reverse_file_name > $outputdir/$reverse_file_name.md5");	
		open REVMD5, "$outputdir/$reverse_file_name.md5";
		$line = <REVMD5>;
		@elements = split (/\s+/, $line);
		$reverse_file_md5 = $elements[0];

	print CSV "CS_$sample_alias\t$instrument_model\t$sample_alias\t$library_source\t$library_selection\t$library_strategy\t$design_description\t$library_construction_protocol\t$insert_size\t$forward_file_name\t$forward_file_md5\t$reverse_file_name\t$reverse_file_md5\n";
	print FOFN "$inputdir/$forward_file_name\n$inputdir/$reverse_file_name\n";
	}
}

close (CSV);
close (FOFN);