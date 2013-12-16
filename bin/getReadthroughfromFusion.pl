#!/usr/bin/perl
$rt=shift @ARGV;
$ft=shift @ARGV;

$RTlist="Namelist";
%Name=();
open($in,"<",$RTlist);
while(<$in>){
	@m=split(/\s+/,$_);
	$Name{$m[0]}=1;
	$Name{$m[1]}=1;
}
for ($i=1;;$i++){
	if(-e "$ft/R$i"){
		if(-e "$ft/R$i/output/readthrough"){
			open($in2,"<","$ft/R$i/output/readthrough");
			@RT=<$in2>;
			@m=split(/\s+/,$RT[1]);
			if($Name{$m[6]} == 1 || $Name{$m[8]} == 1 ){next;}
			else{system("cat $ft/R$i/output/readthrough >> $rt");}
		}
	}
	else{last;}
}
			
		


