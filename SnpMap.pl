#!/usr/bin/perl

$name = "LEb6-1";

if(@ARGV != 1) { die "./CloudMap.pk $name\n"; }
$name = $ARGV[0];

$workdir = "/home/ji/data/Mao/$name";
`mkdir -p $workdir`;
`ln -s /home/ji/data/Mao/fastq/$name\_R1.fastq $workdir`;
`ln -s /home/ji/data/Mao/fastq/$name\_R2.fastq $workdir`;

$fastq = "$workdir/$name";

##############################################################################################################
$ref_gnm = "/home/ji/data/ref/ce10/ce10.fasta";
$picrd_dir = "/home/ji/tools/picard-tools/picard-tools-1.100";
$gatk = "/home/ji/tools/GATK/GenomeAnalysisTK.jar";
##############################################################################################################


### Mapping ### S5
$time = `date`;
chomp $time;
print "$time\tMapping...\n";
`bwa aln -n 0.04 -o 1 -e -1 -d 16 -i 5 -l 1 -k 2 -M 3 -O 11 -E 4 $ref_gnm $fastq\_R1.fastq > $fastq\_R1.sai`;
`bwa aln -n 0.04 -o 1 -e -1 -d 16 -i 5 -l 1 -k 2 -M 3 -O 11 -E 4 $ref_gnm $fastq\_R2.fastq > $fastq\_R2.sai`;
`bwa sampe -n 3 -N 10 -a 500 -o 100000 $ref_gnm $fastq\_R1.sai $fastq\_R2.sai $fastq\_R1.fastq $fastq\_R2.fastq > $fastq.sam`;

###  filter unmapped sam ### S7
`samtools view -b -F 4 -S -b $fastq.sam > $fastq.filter.bam`;

### Add or Replace Groups ### S9
$time = `date`;
chomp $time;
print "$time\tAdd or Replace Groups...\n";
`java -Xmx2g -jar $picrd_dir/AddOrReplaceReadGroups.jar ID=1 SM=rgSM LB=rgLB PL=illumina PU=rgPU INPUT=$fastq.filter.bam OUTPUT=$fastq.AddOrReplaceReadGroups.filter.bam`;
`samtools sort $fastq.AddOrReplaceReadGroups.filter.bam $fastq.AddOrReplaceReadGroups.filter.sorted`;
`samtools index $fastq.AddOrReplaceReadGroups.filter.sorted.bam`;

### Realigner Target Creator/Indel Realigner ### S10+S11
$time = `date`;
chomp $time;
print "$time\tRealigner Target Creator/Indel Realigner...\n";
`java -Xmx2g -jar $gatk -T RealignerTargetCreator -I $fastq.AddOrReplaceReadGroups.filter.sorted.bam -R $ref_gnm -o $fastq.IndelRealigner.intervals`;
`java -Xmx2g -jar $gatk -T IndelRealigner -I $fastq.AddOrReplaceReadGroups.filter.sorted.bam -R $ref_gnm -targetIntervals $fastq.IndelRealigner.intervals -o $fastq.realignedBam.bam`;

### Mark Duplicate ### S12
$time = `date`;
chomp $time;
print "$time\tMark Duplicate...\n";
`java -Xmx2g -jar $picrd_dir/MarkDuplicates.jar REMOVE_DUPLICATES=true ASSUME_SORTED=true READ_NAME_REGEX='[a-zA-Z0-9]+.*:[0-9]:([0-9]+):([0-9]+):([0-9]+)\$' OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 INPUT=$fastq.realignedBam.bam OUTPUT=$fastq.realignedBam.dedupped.bam M=$fastq.dupmetrics.dat`;

### Unified Genotyper ### S13
$time = `date`;
chomp $time;
print "$time\tUnified Genotyper...\n";
`samtools sort $fastq.realignedBam.dedupped.bam $fastq.realignedBam.dedupped.bam.sorted`;
`samtools index $fastq.realignedBam.dedupped.bam.sorted.bam`;
`java -Xmx2g -jar $gatk -T UnifiedGenotyper -R $ref_gnm -I $fastq.realignedBam.dedupped.bam.sorted.bam -glm BOTH -stand_call_conf 30 -stand_emit_conf 30 -hets 0.001 -pcr_error 0.0001 -mbq 17 -deletions 0.05 -maxAltAlleles 5 -minIndelCnt 5 -indelHeterozygosity 0.000125 -indelGCP 10 -indelGOP 45 -o $fastq.UnifiedGenotyper.vcf`;

### Create a BedGraph of genome coverage ### S14
$time = `date`;
chomp $time;
print "$time\tBedGraph of genome coverage...\n";
`bedtools genomecov -ibam -scale -bga -i $fastq.realignedBam.dedupped.bam.sorted.bam -g /home/ji/data/ref/ce10/ce10_chrom_info.txt > $fastq.S14.BedGraph`;

### SnpSift Filter ### S17+S18
$time = `date`;
chomp $time;
print "$time\tSnpEff...\n";
`cat $fastq.UnifiedGenotyper.vcf | java -jar /home/ji/tools/snpEff/SnpSift.jar filter "isHom( GEN[0] ) & isVariant( GEN[0] )" > $fastq.S17.vcf`;
`grep -w 0\$ $fastq.S14.BedGraph > $fastq.S18.filtered.BedGraph`;

### SnpEff ### S20+S21
`java -Xmx2g -jar /home/ji/tools/snpEff/snpEff.jar eff -i vcf -o vcf -upDownStreamLen 10000 -no None -s $fastq.snpEff_summary1.html WS241 $fastq.S17.vcf > $fastq.S20.vcf`;
`java -Xmx2g -jar /home/ji/tools/snpEff/snpEff.jar eff -i bed -o bed -upDownStreamLen 10000 -no None -s $fastq.snpEff_summary2.html WS241 $fastq.S18.filtered.BedGraph > $fastq.S21.bed`;

$time = `date`;
chomp $time;
print "$time\tFinish!\n";

