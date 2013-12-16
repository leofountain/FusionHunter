#!/usr/bin/perl


open($in,"<",shift@ARGV);
$S = <$in>;
while(<$in>){
	@m = split(/\s+/);
	$chr = $m[1];
	$N = $m[7];
	@L = split(/,/,$m[8]);
	@S = split(/,/,$m[9]);
	
	for $i (0 .. $N-2){
		$hash{$chr}{$S[$i]}{$L[$i+1]} = 1;
	}
#@	print $chr,"\t",$S[$i]+$L[$i],"\t",$S[$i+1],"\n";
}

foreach $key (keys %hash){
	foreach $b (keys %{$hash{$key}}){
		foreach $e (keys %{$hash{$key}{$b}}){
			print $key,".",$b,".",$e,"\n";
		}
	}
}
	
