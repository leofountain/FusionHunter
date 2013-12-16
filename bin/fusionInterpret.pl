#!/usr/bin/perl

$FS_ori = shift @ARGV;
$ref = shift @ARGV;
$dir = shift @ARGV;
open($in,"<",$FS_ori);
@head = <$in>;
@h=split(/\s+/,$head[0]);
@m=($h[3]=~m/\d+/g);
$k=index($h[3],"[");
$flag=0;
$flag2=0;
for $i (1 .. $#head-2){
	if(/>/,$head[$i]){
	@h1=split(/\s+/,$head[$i]);
	if($h1[$#h1] eq "known" && $h1[$#h1-2] eq "known"){
		$flag=1;
		last;
	}
	if($h1[$#h1] eq "known" || $h1[$#h1-2] eq "known"){
                $flag2=1;
        }
	}
}
close($in);

@e=($head[$#head-1]=~m/\d+/g);
if($m[0] != $m[3]){
#	if(($flag == 0 && $e[0]<3) || ($flag == 0 && $flag2 == 0 && $e[0]<5)){
#		system("mv $FS_ori $dir/fusion.filteredout");
#		system("touch $FS_ori");
#	}
	
}
else{
	if(substr($h[3],$k+1,2) ne "++"){exit 0;}
	if($m[2]>$m[4]){exit (0);}
	@h=split(/\s+/,$head[1]);
	@s1=split(/:|-/,$h[4]);
	@s2=split(/:|-/,$h[5]);
	if(abs($s1[2]-$s2[1])<60000){
		system("mv $FS_ori $dir/readthrough");
		system("touch $FS_ori");
	}
	else{
		
		$chr="chr$m[0]";
		open($in2,"<",$ref);
		$end=$s1[2];
		$begin=$s2[1];
		while(<$in2>){
			@t=split(/\s+/,$_);
			if($t[4] eq $h[6]){
				if($t[3]>$end){$end = $t[3];}
			}
			if($t[4] eq $h[8]){
				if($t[2]<$begin){$begin = $t[2];}
			}
				
		}
		open($in2,"<",$ref);
                $flag=0;
		%exonHash=();
                while(<$in2>){
                        @t=split(/\s+/,$_);
                        if($t[1] eq $chr && $t[2] > $end && $t[3] < $begin && $t[4] ne $h[6] && $t[4] ne $h[8]){
				if($exonHash{$t[4]} eq ""){
					$flag++;
					$exonHash{$t[4]} = 1;
				}
			}
                }
		open($in2,"<",$ref);
                %geneHash=();;
                
                while(<$in2>){
                        @t=split(/\s+/,$_);
			if($exonHash{$t[4]}==1){
				if($t[2] < $end || $t[3] > $begin){
					$geneHash{$t[4]}=1;
				}
			}
		}
		$c=$flag-(scalar (keys %geneHash));
		if($flag<2 || abs($begin-$end)<30000 || $c < 2){
			system("mv $FS_ori $dir/readthrough");
                	system("touch $FS_ori");
		}
	}
}




