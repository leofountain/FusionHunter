#!/usr/bin/perl

$FS_ori = shift @ARGV;
$ref = shift @ARGV;
#$dir = shift @ARGV;
$junctionHash=();
open($in2,"<",$ref);

while(<$in2>){
        chomp;
	$junctionHash{$_}=1;
            
}
open($in,"<",$FS_ori);
@head = <$in>;
@h=split(/\s+/,$head[0]);
@m=($h[3]=~m/\d+/g);
$k=index($h[3],"[");
$flag=0;
for $i (0 .. $#head-1){
	@pp=split(/\s+/,$head[$i]);
	if($pp[0] eq "#" && $pp[1] eq "Fusion:"){
		@h1=split(/\s+/,$head[$i+1]);
		@d1=($h1[4]=~m/\d+/g);
		@d2=($h1[5]=~m/\d+/g);
		$chr = "chr$d1[0]";
		$t = $d2[1]-1;
		$junction = "$chr.$d1[2].$t";
		if($junctionHash{$junction} eq ""){
			$flag=1;
			$pp[1] = "Readthrough:";
			print join(" ",@pp),"\n";
			next;
		}
		else{
			$flag=0;
		}
	}
	if($flag == 1){
		print $head[$i];
	}

}
