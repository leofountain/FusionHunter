#!/bin/sh

# install commands for FusionHunter


cd lib/ ; make ; cd ../

cd src/ ; make ; cd ../



# check and install ForkManager package for parallel processes with Perl
# downloaded from 
# http://search.cpan.org/~dlux/Parallel-ForkManager-0.7.5/ForkManager.pm
#
#
# if you are not a sudoer in your server, ask someone who has the sudo permission to help install ForkManager package
#


CMD="make install"

cd Parallel-ForkManager-0.7.5/ ; perl Makefile.PL ; make ; make test
echo $CMD
if $CMD ; then
	echo "\nParallel::ForkManager has been successfully installed"
else
	echo "\ninstallation for Parallel::ForkManager failed, maybe you need a root privilege"
		exit 1
fi




