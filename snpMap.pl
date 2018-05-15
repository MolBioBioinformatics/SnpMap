#!/usr/bin/perl
use File::Basename;

### Read input ###
for($i=0;$i<@ARGV;$i++){
  if($ARGV[$i] eq "-f"){$func = $ARGV[$i+1];}
  if($ARGV[$i] eq "-config"){$config = $ARGV[$i+1];}
  if($ARGV[$i] eq "-t"){$ref = $ARGV[$i+1]}
  if($ARGV[$i] eq "-in"){$fq1 = $ARGV[$i+1];$fq2 = $ARGV[$i+2]}
  if($ARGV[$i] eq "-d"){$snp_dir = $ARGV[$i+1]}
}

if($func eq "align"){
  if($config eq "" | $ref eq "" | $fq1 eq ""){
    die "./SnpMap.pl -f align -t bwa/ref.fa -in input.fq -config config.txt\n";
  }
}else{
  if($func eq "freq"){
    if($snp_dir eq ""){die "./SnpMap.pl -f freq -d sample/snp\n";}
  }else{
    die "./SnpMap.pl -f align\n./SnpMap.pl -f freq\n";
  }
}

if($fq2 =~ /-/ | $fq2 eq ""){ # fq2 another flag?
  $mode = "single";
}else{
  $mode = "pair";
}

### Read Config ###
open IN, "<$config";
while(<IN>){
  chomp $_;
  @tmp = split(/\t/,$_);
  if($tmp[0] eq "BWA"){$bwa = $tmp[1]}
  if($tmp[0] eq "Picard"){$picrd_dir = $tmp[1]}
  if($tmp[0] eq "GATK"){$gatk = $tmp[1]}
  if($tmp[0] eq "snpEff"){$snpEff = $tmp[1]}
  if($tmp[0] eq "snpEff_ref_gnm"){$snpEff_ref = $tmp[1]}
}

$curdir = `pwd`;
chomp $curdir;

if($func eq "align"){goto ALIGN;}
if($func eq "freq"){goto FREQ;}

ALIGN:

$name = basename($fq1);
$name =~ s/.fq//g;
$name =~ s/.fastq//g;

$workdir = "$curdir/$name";
`mkdir -p $workdir`;

$ref_id = basename($ref);
$ref_file = basename($ref);
@t = split(/\./,$ref_id);
pop @t;
$ref_id = join(".",@t);

####################################################################
#$ref_gnm = "/data/ji/ref/ce10/ce10.fasta";
#$gnm_info = "/data/ji/ref/ce10/ce10_chrom_info.txt";
#$picrd_dir = "/data/ji/tools/picard-tools/picard-tools-1.100";
####################################################################
#
### Mapping ### S5
$time = `date`;
chomp $time;
print "$time\tMapping...\n";
`cp $ref* $workdir/`;
if($mode eq "pair"){
  `$bwa aln -n 0.04 -o 1 -e -1 -d 16 -i 5 -l 1 -k 2 -M 3 -O 11 -E 4 $ref $fq1 > $workdir/fq1.sai`;
  `$bwa aln -n 0.04 -o 1 -e -1 -d 16 -i 5 -l 1 -k 2 -M 3 -O 11 -E 4 $ref $fq2 > $workdir/fq2.sai`;
  `$bwa sampe -n 3 -N 10 -a 500 -o 100000 $ref $workdir/fq1.sai $workdir/fq2.sai $fq1 $fq2 > $workdir/$name.sam`;
}
if($mode eq "single"){
  `$bwa aln -n 0.04 -o 1 -e -1 -d 16 -i 5 -l 1 -k 2 -M 3 -O 11 -E 4 $ref $curdir/$fq1 > $workdir/fq1.sai`;
  `$bwa samse -n 100000 -f $workdir/$name.sam $ref $workdir/fq1.sai $fq1`;
}

`samtools view -S -H $workdir/$name.sam > $workdir/$ref_id.dict`;
###  filter unmapped sam ### S7
`samtools view -b -F 4 -S -b $workdir/$name.sam > $workdir/$name.filter.bam`;

### Add or Replace Groups ### S9
$time = `date`;
chomp $time;
print "$time\tAdd or Replace Groups...\n";
`java -Xmx2g -jar $picrd_dir/AddOrReplaceReadGroups.jar ID=1 SM=rgSM LB=rgLB PL=illumina PU=rgPU INPUT=$workdir/$name.filter.bam OUTPUT=$workdir/$name.AddOrReplaceReadGroups.filter.bam`;
`samtools sort -o $workdir/$name.AddOrReplaceReadGroups.filter.sorted.bam $workdir/$name.AddOrReplaceReadGroups.filter.bam`;
`samtools index $workdir/$name.AddOrReplaceReadGroups.filter.sorted.bam`;

### Realigner Target Creator/Indel Realigner ### S10+S11
$time = `date`;
chomp $time;
print "$time\tRealigner Target Creator/Indel Realigner...\n";
`java -Xmx2g -jar $gatk -T RealignerTargetCreator -I $workdir/$name.AddOrReplaceReadGroups.filter.sorted.bam -R $workdir/$ref_file -o $workdir/$name.IndelRealigner.intervals`;
#print "java -Xmx2g -jar $gatk -T RealignerTargetCreator -I $workdir/$name.AddOrReplaceReadGroups.filter.sorted.bam -R $workdir/$ref_file -o $workdir/$name.IndelRealigner.intervals\n";
`java -Xmx2g -jar $gatk -T IndelRealigner -I $workdir/$name.AddOrReplaceReadGroups.filter.sorted.bam -R $workdir/$ref_file -targetIntervals $workdir/$name.IndelRealigner.intervals -o $workdir/$name.realignedBam.bam`;

### Mark Duplicate ### S12
$time = `date`;
chomp $time;
print "$time\tMark Duplicate...\n";
`java -Xmx2g -jar $picrd_dir/MarkDuplicates.jar REMOVE_DUPLICATES=true ASSUME_SORTED=true READ_NAME_REGEX='[a-zA-Z0-9]+.*:[0-9]:([0-9]+):([0-9]+):([0-9]+)\$' OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 INPUT=$workdir/$name.realignedBam.bam OUTPUT=$workdir/$name.realignedBam.dedupped.bam M=$workdir/$name.dupmetrics.dat`;

### Unified Genotyper ### S13
$time = `date`;
chomp $time;
print "$time\tUnified Genotyper...\n";
`samtools sort -o $workdir/$name.realignedBam.dedupped.bam.sorted.bam $workdir/$name.realignedBam.dedupped.bam`;
`samtools index $workdir/$name.realignedBam.dedupped.bam.sorted.bam`;
`java -Xmx2g -jar $gatk -T UnifiedGenotyper -R $workdir/$ref_file -I $workdir/$name.realignedBam.dedupped.bam.sorted.bam -glm BOTH -stand_call_conf 30 -stand_emit_conf 30 -hets 0.001 -pcr_error 0.0001 -mbq 17 -deletions 0.05 -maxAltAlleles 5 -minIndelCnt 5 -indelHeterozygosity 0.000125 -indelGCP 10 -indelGOP 45 -o $workdir/$name.UnifiedGenotyper.vcf`;

## Create a BedGraph of genome coverage ### S14
$time = `date`;
chomp $time;
print "$time\tBedGraph of genome coverage...\n";
`samtools faidx $workdir/$ref_file`;
`cat $workdir/$ref_file.fai | cut -f1-2 > $workdir/$ref_id.chrom_info.txt`;
`bedtools genomecov -ibam -scale -bga -i $workdir/$name.realignedBam.dedupped.bam.sorted.bam -g $workdir/$ref_id.chrom_info.txt > $workdir/$name.S14.BedGraph`;

### SnpSift Filter ### S17+S18
$time = `date`;
chomp $time;
print "$time\tSnpEff...\n";
`cat $workdir/$name.UnifiedGenotyper.vcf | java -jar $snpEff/SnpSift.jar filter "isHom( GEN[0] ) & isVariant( GEN[0] )" > $workdir/$name.S17.vcf`;
`grep -w 0\$ $workdir/$name.S14.BedGraph > $workdir/$name.S18.filtered.BedGraph`;

### SnpEff ### S20+S21
`java -Xmx2g -jar $snpEff/snpEff.jar eff -no-downstream -no-intergenic -no-intron -no-upstream -i vcf -o vcf -upDownStreamLen 10000 -no None -s $workdir/$name.snpEff_summary1.html $snpEff_ref $workdir/$name.S17.vcf > $workdir/$name.S20.vcf`;
`java -Xmx2g -jar $snpEff/snpEff.jar eff -i bed -o bed -upDownStreamLen 10000 -no None -s $workdir/$name.snpEff_summary2.html $snpEff_ref $workdir/$name.S18.filtered.BedGraph > $workdir/$name.S21.bed`;

### vcf file to anno table ###
open IN, "<$workdir/$name.S20.vcf";
open OUT, ">$workdir/$name.snp.txt";
while(<IN>){
  if($_ =~ /^#/){
  }
  else{
    chomp $_;
    @line = split(/\t/,$_);
    for($i=0;$i<7;$i++){
      $chr = $line[0];
      $pos = $line[1];
      $ref = $line[3];
      $alt = $line[4];
      $qua = $line[5];
    }
    @score = split(/\;/,$line[7]);
    $ac = $score[0];
    $af = $score[1];
    $an = $score[2];
    $dp = $score[3];
    @info_array = split(/\,/,$score[13]);
    foreach $info (@info_array){
      @tmp = split(/\|/,$info);
      $header = $tmp[0];
      $annotation = $tmp[1];
      $Annotation_Impact = $tmp[2];
      $Gene_name = $tmp[3];
      $Gene_ID = $tmp[4];
      $Feature_ID = $tmp[6];
      $Transcript_BioType = $tmp[7];
      $Rank = $tmp[8];
      $HGVS_C = $tmp[9];
      $HGVS_P = $tmp[10];
      $cDNA_pos = $tmp[11];
      $AA_pos = $tmp[13];
      $Distance = $tmp[14];
      if($annotation eq "missense_variant" || $annotation eq "missense_variant&splice_region_variant" || $annotation eq "splice_region_variant&intron_variant" || $annotation eq "splice_region_variant&non_coding_exon_variant" || $annotation eq "splice_region_variant&synonymous_variant" || $annotation eq "stop_gained"){
        print OUT "$chr\t$pos\t$ref\t$alt\t$qua\t";
        print OUT "$dp\t";
        print OUT "$annotation\t$Gene_ID\t$Gene_name";
        print OUT "\n";
      }
    }
  }
}
`cat $workdir/$name.snp.txt | sort | uniq > $name.snp.txt`;

$time = `date`;
chomp $time;
print "$time\tFinish!\n";

exit;
FREQ:
@list = `ls $snp_dir/*.txt`;
my %cnt;
foreach $file (@list){
  open IN, "<$file";
  while(<IN>){
    chomp $_;
    @tmp = split(/\t/,$_);
    $snp = "$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[3]\t$tmp[6]\t$tmp[7]\t$tmp[8]";
    $cnt{$snp}++;
  }
}
open OUT, ">gene.freq.txt.tmp";
foreach $snp (keys %cnt){
  print OUT "$snp\t$cnt{$snp}\n";
}
`cat gene.freq.txt.tmp | sort -k8,8nr > gene.freq.txt`;
`rm gene.freq.txt.tmp`;

