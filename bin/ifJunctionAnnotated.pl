#!/usr/bin/perl
%hash = ();
$sam = shift @ARGV;
$present = shift @ARGV;
open($in,"<",$sam)|| die "cannot open";
while(<$in>){
	@m=split(/\s+/,$_);
	@s=split(/:/,$m[12]);
	$hash{$m[0]} =$s[2];
}

open($in1,"<",$present);
while(<$in1>){
	chomp($_);
	if(/>/){
		@m=split(/\s+/,$_);
		print $_;
		if($hash{$m[1]} == 0){print "\tunknown x unknown\n";}
		if($hash{$m[1]} == 1){print "\tknown x unknown\n";}
		if($hash{$m[1]} == 2){print "\tunknown x known\n";}
		if($hash{$m[1]} == 3){print "\tknown x known\n";}
	}
	else{print "$_\n";}		
}
	
