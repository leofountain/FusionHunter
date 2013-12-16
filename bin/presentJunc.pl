#!/usr/bin/perl


open($in,"<",shift@ARGV);
$S = <$in>;
while(<$in>){
	@m = split(/\s+/);
	$chr = $m[1];
	$N = $m[7];
	@L = split(/,/,$m[8]);
	@S = split(/,/,$m[9]);
	
	
	for $i (0 .. $N-1){

		if($i != $N-1){
			$M = $S[$i] + 1;
			$N = $L[$i+1] + 1;
			$hash{"$chr\t$M\t$N"} = $m[$#m];
		}
#		$exon{"$chr\t$L[$i]\t$S[$i]"} = $m[$#m];
	}

#@	print $chr,"\t",$S[$i]+$L[$i],"\t",$S[$i+1],"\n";
}

foreach $key (keys %hash){
	print $key,"\t",$hash{$key},"\n";
}

#foreach $key (keys %exon){
#	print "@\t",$key,"\t",$exon{$key},"\n";
#}

