#!/usr/bin/perl
%hash = ();
my $this_arg = shift @ARGV;
open($in,"<",$this_arg)|| die "cannot open";
while(<$in>){
	@m=split(/\s+/,$_);
	if($hash{$m[0]} eq ""){$hash{$m[0]} =1;}
	else {$hash{$m[0]}++;}
	if($hash{$m[1]} eq ""){$hash{$m[1]} =1;}
        else {$hash{$m[1]}++;}	
}
close($this_arg);
open($in,"<",$this_arg);
while(<$in>){
        @m=split(/\s+/,$_);
	$s = substr($m[4],1,-1);
        if( ( $hash{$m[0]} > 2 || $hash{$m[1]} > 2 ) && $s < 5){next;}
	print $_;
}
	
