#!/usr/bin/perl
##GappedAlignment
use Parallel::ForkManager;

while(scalar(@ARGV)>0){
        my $this_arg = shift @ARGV;
        if($this_arg eq '-d'){
                $D=shift @ARGV;
        }
        elsif($this_arg eq '-p'){
                $P=shift @ARGV;
        }
        elsif($this_arg eq '-r'){
                $Ref=shift @ARGV;
        }
	elsif($this_arg eq '-i'){
		$Input_fq=shift @ARGV;
	}
	elsif($this_arg eq '-e'){
		$Exon = shift @ARGV;
	}
	elsif($this_arg eq '-mask'){
		$MASK = shift @ARGV;
	}
	elsif($this_arg eq '-src'){
                $SRC = shift @ARGV;
        }
	elsif($this_arg eq '-t'){
                $TILE = shift @ARGV;
        }
	elsif($this_arg eq '-b'){
                $Bin = shift @ARGV;
        }
	elsif($this_arg eq '-o'){
                $MINOVLP = shift @ARGV;
        }
	elsif($this_arg eq '-v'){
                $totalMis = shift @ARGV;
        }
	elsif($this_arg eq '-h'){
                $halfLength = shift @ARGV;
        }
	elsif($this_arg eq '-rg'){
                $RG = shift @ARGV;
        }


        else{print "Unexpected Parameter";}
}
open($in,"<",$RG);
@FILE=<$in>;
$N=scalar @FILE;
$pm=new Parallel::ForkManager($P);
for ($i=1;$i<=$N;$i++){
        if(-e "$D/R$i"){
		if(! -z "$D/R$i/reads.fq"){
                my $pid = $pm->start and next;
                #print "Running singleGappedAlignment on $D/R$i\n";
		system("cd $D/R$i; $Bin/bowtie-build ref.fa ref 1> run.log");
		system("cd $D/R$i; $Bin/bowtie --suppress 6 --quiet -m 1 -v 2 -p 2 --un reads.un.fq --max reads.max.fq ref reads.fq reads.bwt");
		system("cd $D/R$i; $Bin/new_splitReads reads.un.fq -r=$halfLength -ncbi > split.fq");
		system("cd $D/R$i; $Bin/bowtie --suppress 6 --quiet -k 8 -m 20 -v 2 -p 2 ref split.fq split.bwt");
		system("cd $D/R$i; $Bin/singleGappedAlignment $halfLength $totalMis $MINOVLP $TILE $MASK ../../$Exon reads.bwt reads.un.fq ref.fa split.bwt>GappedAlignment.sam");
		system("cd $D/R$i; mkdir output/ 2>run.log;cp GappedAlignment.sam output/final.sam");
		$pm->finish;
		}
		else{
		system("cd $D/R$i; mkdir output/ 2>run.log; touch output/final.sam");
		}
             
        }
        else{
                last;
        }
}
$pm->wait_all_children;

