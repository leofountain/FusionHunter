ARCH ?= $(shell uname -m)
CC = gcc
CXX = g++ 
CDEBUG =#-ggdb -g #-pg
OPTM = -O
OMP = -fopenmp
WARN = -W -Wall
KLIB = ../lib/$(ARCH)
KINC = ../inc
CFLAGS = $(WARN) $(OPTM) -I. -I$(KINC)
CLIB = $(KLIB)/jkweb.a -lm
BIN = ../bin

RM = rm -rf
DATE = `date +'%F'`

ALLSRC = splitReads putativeExons readChains constructRegions regionPairs \
				 leftRightOvlp chainPairs postLeftRightOvlp considerKnownGenes \
				 knownExonIntron onlyConsiderKnownGenes regionPairsList reduceBwt \
				 trimOriginalBwt removePcrReplicate  \
				 allMakeDirFasta allWriteReadsToDir allFusionJunction \
				 singleMakeDirFasta  singleFusionJunction \
				  new_splitReads 

OBJ = myutil.o

all: $(OBJ) $(ALLSRC) getReadThroughCandidate singleGappedAlignment reMapFromJunction

%.o: %.c %.h
	$(CC) $(CDEBUG) $(CFLAGS) -c $(addsuffix .c, $(basename $@)) 

%: %.c $(OBJ) $(KLIB)/jkweb.a
	$(CC) $(CDEBUG) $(CFLAGS) $(OMP) $+ $(CLIB) -o $@
	mv $@ $(BIN); strip $(BIN)/$@

singleGappedAlignment:
	$(CXX) $(WARN) -o singleGappedAlignment singleGappedAlignment.cpp
	mv singleGappedAlignment $(BIN); strip $(BIN)/singleGappedAlignment

getReadThroughCandidate:
	$(CXX) $(WARN) -o getReadThroughCandidate getReadThroughCandidate.cpp
	mv getReadThroughCandidate $(BIN); strip $(BIN)/getReadThroughCandidate

reMapFromJunction:
	$(CXX) -o reMapFromJunction reMapFromJunction.cpp
	mv reMapFromJunction $(BIN); strip $(BIN)/reMapFromJunction

.PHONY: clean
clean:
	$(RM) $(OBJ) $(ALLSRC) *.dSYM singleGappedAlignment getReadThroughCandidate
	cd $(BIN); $(RM) $(ALLSRC) singleGappedAlignment getReadThroughCandidate reMapFromJunction
.PHONY: bak
bak:
	mkdir -p BAK/bak.$(DATE)
	cp *.c *.h Makefile BAK/bak.$(DATE)

