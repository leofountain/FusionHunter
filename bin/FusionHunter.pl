#!/usr/bin/perl

##FusionHunter Pipeline
use POSIX;
use FindBin qw($Bin);
$usage = qq{

	Copy FusionHunter.cfg from FusionHunter directory to your working space and edit its variables
		Then type: DirtoFusionHunter/bin/FusionHunter.pl FusionHunter.cfg

};

if(scalar(@ARGV) == 0){
	print $usage;
	exit(0);
}
while(scalar(@ARGV)>0)
{
	my $this_arg = shift @ARGV;
	if($this_arg eq 'FusionHunter.cfg')
	{
		$in=$this_arg;
	}    
	else{print "Unexpected parameter";}
}

open($GAcfg,"<",$in);
while(<$GAcfg>){
	chomp;
	next unless $_;
	if(/#/){
		next;
	}
	@GAcmd=split(/\s*=\s*/,$_);
	if($GAcmd[0] eq 'CORE'){ $CORE = $GAcmd[1];}
	elsif($GAcmd[0] eq 'L'){ $L = $GAcmd[1];}
	elsif($GAcmd[0] eq 'R'){ $R = $GAcmd[1];}
	elsif($GAcmd[0] eq 'Reference'){ $Reference = $GAcmd[1];}
	elsif($GAcmd[0] eq 'BowtieIdx'){ $BowtieIdx = $GAcmd[1];}
	elsif($GAcmd[0] eq 'Gene_annotation'){ $Gene_annotation = $GAcmd[1];}
	elsif($GAcmd[0] eq 'Repeats'){ $Repeats = $GAcmd[1];}
	elsif($GAcmd[0] eq 'EST'){ $EST = $GAcmd[1];}
	elsif($GAcmd[0] eq 'SelfAlign'){ $SelfAlign = $GAcmd[1];}
	elsif($GAcmd[0] eq 'M1'){ $M1 = $GAcmd[1];}
	elsif($GAcmd[0] eq 'K1'){ $K1 = $GAcmd[1];}
	elsif($GAcmd[0] eq 'M2'){ $M2 = $GAcmd[1];}
	elsif($GAcmd[0] eq 'K2'){ $K2 = $GAcmd[1];}
	elsif($GAcmd[0] eq 'REAPTOVLP'){ $REAPTOVLP = $GAcmd[1];}
	elsif($GAcmd[0] eq 'CHAINNUM'){ $CHAINNUM = $GAcmd[1];}
	elsif($GAcmd[0] eq 'RPTOVLP'){ $RPTOVLP = $GAcmd[1];}
	elsif($GAcmd[0] eq 'CHAINOVP'){ $CHAINOVP = $GAcmd[1];}
	elsif($GAcmd[0] eq 'CHAINDIS'){ $CHAINDIS = $GAcmd[1];}
	elsif($GAcmd[0] eq 'READOVLP'){ $READOVLP = $GAcmd[1];}
	elsif($GAcmd[0] eq 'PAIRNUM'){ $PAIRNUM = $GAcmd[1];}
	elsif($GAcmd[0] eq 'MINSPAN'){ $MINSPAN = $GAcmd[1];}
	elsif($GAcmd[0] eq 'MINOVLP'){ $MINOVLP = $GAcmd[1];}
	elsif($GAcmd[0] eq 'segment_size'){ $SZ = $GAcmd[1];}
	elsif($GAcmd[0] eq 'fusion_output'){ $JK = $GAcmd[1];}
	elsif($GAcmd[0] eq 'MASK'){ $MASK = $GAcmd[1];}
	elsif($GAcmd[0] eq 'TILE'){ $TILE = $GAcmd[1];}
	elsif($GAcmd[0] eq 'MISMATCH'){ $V = $GAcmd[1];}
	elsif($GAcmd[0] eq 'readthrough_output'){ $RTname = $GAcmd[1];}
	elsif($GAcmd[0] eq 'IF_HUMAN'){ $IF_HUMAN = $GAcmd[1];}



}

$D="FusionHunter_temp";
$RTdir = "Readthrough_temp";
$RG="fusion.list";
$I="Z";
$X="X";
$LN="left";
$RN="right";
$now_string = localtime;
print "$now_string:   ";
print "Checking configurations\n"; 
system("echo $now_string:   Checking configurations > run.log");
if($L eq $R){die "You should check your configuration file, L read should not be the same as R read!\n";}
if(! -e $Reference){die "$Reference not found\n";}
if(! -e $Gene_annotation){die "$Gene_annotation not found\n";}
if($IF_HUMAN){
	if(! -e $Repeats){die "$Repeats not found\n";}
	if(! -e $SelfAlign){die "$SelfAlign not found\n";}
}
else{
	if(! -e $Repeats){print "FusionHunter will run without repeat annotation!\n";open I,">Repeats_fake";$Repeats = "Repeats_fake";}
	if(! -e $SelfAlign){print "FusionHunter will run without self-chain annotation!\n"; open I,">SelfAlign_fake";$SelfAlign = "SelfAlign_fake";}
	if(! -e $EST){
		system("$Bin/EST2p.pl $Gene_annotation > MY_EST");
		$EST = "MY_EST";
	}
}
if(! -e "$BowtieIdx.1.ebwt"){die "Bowtie index not found\n";}
if($IF_HUMAN){
	@t = split(/\/|\./,$Repeats);

	$hg = $t[$#t-1];
	$Refseq="$Bin/refseq/refseqExons.$hg";
	if(! -e $Refseq ){die "$Refseq not found\n";}
}
else{
	system("$Bin/EST2point.pl $Gene_annotation > MY_MODEL"); 
	$Refseq = "MY_MODEL";
}



if($TILE < 2){$TILE=2;}
if($V > 4){$V=4;}
if($MINSPAN < 1){$MINSPAN = 1;}
if($MINOVLP < $TILE){$MINOVLP = $TILE;}
#test paired-end reads
open($inL,"<",$L);
open($inR,"<",$R);
my @Llines,@Rlines;
for my $i (0 .. 7){
	my $Lin = <$inL>;
	my $Rin = <$inR>;
	push @Llines,$Lin;
	push @Rlines,$Rin;
}
close($inL);
close($inR);
chomp(@Llines);
chomp(@Rlines);
@Lline1=split(/\s+/,$Llines[0]);
@Lline2=split(/\s+/,$Llines[4]);
@Rline1=split(/\s+/,$Rlines[0]);
@Rline2=split(/\s+/,$Rlines[4]);
$read_length = length($Llines[1]);
#print $read_length,"\n";
if(substr ($Lline1[0],-1) == 1 && substr ($Lline2[0],-1) == 1 && substr ($Rline1[0],-1) == 2 && substr ($Rline2[0],-1) == 2){$NB = "";}
elsif($Lline1[0] eq $Rline1[0] && $Lline2[0] eq $Rline2[0]){
	$NB = -ncbi;
}
else{print "Check again your reads\n";exit 0;}
#if($IF_HUMAN){
	system("sed -e '/uc010axl.1/d' -e '/uc010tyt.1/d' $Gene_annotation > annotation");
#}
#else{
#	`cat $Gene_annotation > annotation`;
#}

#map_full:

$now_string = localtime;
print "$now_string:   ";
print "Mapping full length reads against $Reference\n";
system("echo $now_string:   Mapping full length reads against $Reference >> run.log");
system("$Bin/bowtie --suppress 6 -k $K2 -m $M2 -v 2 -p $CORE --un L.un.fq $BowtieIdx $L ${X}_$LN.full.bwt 2>>run.log");
system("$Bin/bowtie --suppress 6 -k $K2 -m $M2 -v 2 -p $CORE --un R.un.fq $BowtieIdx $R ${X}_$RN.full.bwt 2>>run.log");
##split:
#
system("$Bin/presentJunc.pl $Gene_annotation > SJS");
system("$Bin/reMapFromJunction SJS $BowtieIdx.fa $read_length > REMAP.fa");
system("$Bin/bowtie-build --quiet REMAP.fa REMAP; $Bin/bowtie --quiet -k 20 -v 3 -p $CORE --un L.REMAP.un.fq REMAP L.un.fq L.REMAP.bwt; mv L.REMAP.un.fq L.un.fq ");
system("$Bin/bowtie-build --quiet REMAP.fa REMAP; $Bin/bowtie --quiet -k 20 -v 3 -p $CORE --un R.REMAP.un.fq REMAP R.un.fq R.REMAP.bwt; mv R.REMAP.un.fq R.un.fq ");


$now_string = localtime;
print "$now_string:   ";
print "Split original reads to $SZ bp long partial reads\n";
system("echo $now_string:   Split original reads to $SZ bp long partial reads >> run.log");
system("$Bin/splitReads L.un.fq -r=$SZ $NB -s=0 > ${X}_$LN.part.fq");
system("$Bin/splitReads R.un.fq -r=$SZ $NB -s=1 > ${X}_$RN.part.fq");

##map_split:
$now_string = localtime;
print "$now_string:   ";
print "Mapping $SZ long partial reads against $Reference\n";
system("echo $now_string:   Mapping $SZ long partial reads against $Reference >> run.log");
system("$Bin/bowtie --suppress 6 -k $K1 -m $M1 -v 2 -p $CORE $BowtieIdx ${X}_$LN.part.fq ${X}_$LN.part.bwt 2>>run.log");
system("$Bin/bowtie --suppress 6 -k $K1 -m $M1 -v 2 -p $CORE $BowtieIdx ${X}_$RN.part.fq ${X}_$RN.part.bwt 2>>run.log");

##trim:
$now_string = localtime;
print "$now_string:   ";
print "Trimming bowtie results\n";
system("echo $now_string:   Trimming bowtie results >> run.log");
system("$Bin/trimOriginalBwt $NB ${X}_$LN.full.bwt > ${X}_$LN.full.bwt.tmp");
system("mv -f ${X}_$LN.full.bwt.tmp ${X}_$LN.full.bwt");
system("$Bin/trimOriginalBwt $NB ${X}_$RN.full.bwt > ${X}_$RN.full.bwt.tmp");
system("mv -f ${X}_$RN.full.bwt.tmp ${X}_$RN.full.bwt");
system("$Bin/trimOriginalBwt ${X}_$LN.part.bwt > ${X}_$LN.part.bwt.tmp");
system("mv -f ${X}_$LN.part.bwt.tmp ${X}_$LN.part.bwt");
system("$Bin/trimOriginalBwt ${X}_$RN.part.bwt > ${X}_$RN.part.bwt.tmp");
system("mv -f ${X}_$RN.part.bwt.tmp ${X}_$RN.part.bwt");
system("rm -f ${X}_$LN.part.fq ${X}_$RN.part.fq");

system("$Bin/knownExonIntron annotation > ${I}_known.exon.intron");
system("$Bin/reduceBwt $NB -repeatOverlap=$REAPTOVLP -chainNum=$CHAINNUM ${X}_$LN.full.bwt $Repeats $SelfAlign > ${I}_$LN.full.slim.bwt");
system("$Bin/reduceBwt $NB -repeatOverlap=$REAPTOVLP -chainNum=$CHAINNUM ${X}_$RN.full.bwt $Repeats $SelfAlign > ${I}_$RN.full.slim.bwt");
system("$Bin/reduceBwt $NB -repeatOverlap=$REAPTOVLP -chainNum=$CHAINNUM ${X}_$LN.part.bwt $Repeats $SelfAlign > ${I}_$LN.part.slim.bwt");
system("$Bin/reduceBwt $NB -repeatOverlap=$REAPTOVLP -chainNum=$CHAINNUM ${X}_$RN.part.bwt $Repeats $SelfAlign > ${I}_$RN.part.slim.bwt");

$now_string = localtime;
print "$now_string:   ";
print "Getting filtered pair alignments\n";
system("echo $now_string:   Getting filtered pair alignments >> run.log");
system("$Bin/leftRightOvlp -repeatOverlap=$RPTOVLP ${I}_$LN.full.slim.bwt ${I}_$RN.full.slim.bwt $Repeats $NB > ${I}_pair.mapping.init");
system("$Bin/postLeftRightOvlp -chainOverlap=$CHAINOVP -chainDist=$CHAINDIS ${I}_pair.mapping.init $SelfAlign > ${I}_pair.mapping.second");
system("$Bin/removePcrReplicate ${I}_pair.mapping.second > ${I}_pair.mapping.final");

$now_string = localtime;
print "$now_string:   ";
print "Getting regions covered by alignments and gene annotation\n";
system("echo $now_string:   Getting regions covered by alignments and gene annotation >> run.log");
system("$Bin/putativeExons ${X}_$LN.full.bwt ${X}_$RN.full.bwt ${X}_$LN.part.bwt ${X}_$RN.part.bwt > ${I}_exons.bed");
system("$Bin/constructRegions ${I}_exons.bed > ${I}_regions.bed.init");
system("$Bin/considerKnownGenes ${I}_regions.bed.init annotation > ${I}_regions.bed.final");
system("mv ${I}_regions.bed.final ${I}_regions.bed.all.final");
system("$Bin/onlyConsiderKnownGenes ${I}_regions.bed.all.final annotation > ${I}_regions.bed.final");

##pair:
$now_string = localtime;
print "$now_string:   ";
print "Getting list for candidate fusion genes\n";
system("echo $now_string:   Getting list for candidate fusion genes >> run.log");
system("$Bin/regionPairs -readOverlap=$READOVLP ${I}_regions.bed.final ${I}_pair.mapping.final > ${I}_regions.pairs.init");
system("$Bin/regionPairsList -minPair=$PAIRNUM ${I}_regions.pairs.init > ${I}_regions.pairs.second");
system("$Bin/filter.pl ${I}_regions.pairs.second > $RG ");

if(! -z $RG){

#dir:
	$now_string = localtime;
	print "$now_string:   ";
	print "Getting correlative pseudo reference for each pair of candidate fusion genes\n";
	system("echo $now_string:   Getting correlative pseudo reference for each pair of candidate fusion genes >> run.log");
	system("rm -rf $D/R*");
	system("$Bin/allMakeDirFasta $RG $Reference $D");

#reads:
	$now_string = localtime;
	print "$now_string:   ";
	print "Getting correlative reads for each pair of candidate fusion genes\n";
	system("echo $now_string:   Getting correlative reads for each pair of candidate fusion genes >> run.log");
	system("$Bin/allWriteReadsToDir $NB $RG L.un.fq R.un.fq ${I}_$LN.part.slim.bwt ${I}_$RN.part.slim.bwt $D");
##GappedAlignment:
	$now_string = localtime;
	print "$now_string:   ";
	print "Running gapped alignment for all candidate fusion spanning reads\n";
	system("echo $now_string:   Running gapped alignment for all candidate fusion spanning reads >> run.log");
	system("$Bin/allRunGappedAlignment.pl -rg $RG -h $SZ -v $V -t $TILE -o $MINOVLP -mask $MASK -d $D -p $CORE -e ${I}_known.exon.intron -b $Bin");
#Present:

	$now_string = localtime;
	print "$now_string:   ";
	print "Present gene fusion events in $JK\n";
	system("echo $now_string:   Present gene fusion events in $JK >> run.log");
	system("$Bin/allFusionJunction -minSpan=$MINSPAN -minOverlap=$MINOVLP $RG ${I}_known.exon.intron $D $JK $Bin");
}
else{print "No gene fusion candidates for your data\n";}

$now_string = localtime;
print "$now_string:   ";
print "Searching readthrough\n";
system("echo $now_string:   Searching readthrough >> run.log");
if($NB eq "-ncbi"){$NCBI=1;}
else{$NCBI="";}
system("$Bin/getReadThroughCandidate $PAIRNUM annotation ${I}_pair.mapping.second $Refseq > RG");
$now_string = localtime;
print "$now_string:   ";
print "Getting correlative pseudo reference for each pair of candidate readthrough genes\n";
system("echo $now_string:   Getting correlative pseudo reference for each pair of candidate readthrough genes >> run.log");
system("$Bin/allMakeDirFasta RG $Reference $RTdir");
$now_string = localtime;
print "$now_string:   ";
print "Getting correlative reads for each pair of candidate readthrough genes\n";
system("echo $now_string:   Getting correlative reads for each pair of candidate readthrough genes >> run.log");
system("$Bin/allWriteReadsToDir $NB RG L.un.fq R.un.fq ${I}_$LN.part.slim.bwt ${I}_$RN.part.slim.bwt $RTdir");
$now_string = localtime;
print "$now_string:   ";
print "Running gapped alignment for all candidate readthrough spanning reads\n";
system("echo $now_string:   Running gapped alignment for all candidate readthrough spanning reads >> run.log");
system("$Bin/allRunGappedAlignment.pl -rg RG -h $SZ -v $V -t $TILE -o $MINOVLP -mask $MASK -d $RTdir -p $CORE -e ${I}_known.exon.intron -b $Bin");
$now_string = localtime;
print "$now_string:   ";
print "Present readthrough events in $RTname\n";
system("echo $now_string:   Present gene readthrough events in $RTname >> run.log");
system("$Bin/allFusionJunction -minSpan=$MINSPAN -minOverlap=$MINOVLP RG ${I}_known.exon.intron $RTdir Temp $Bin");
system("$Bin/getReadthroughfromFusion.pl Temp $D");
system("$Bin/readthroughAnno.pl Temp $EST > $RTname");

#clean:

$now_string = localtime;
print "$now_string:   ";

print "Clean temp files\n";
	system("rm ${I}_* ${X}_* annotation RG Namelist fusion.list REMAP* L* R*");
	if(-e "Temp"){
               system("rm Temp");
       }



