MACHTYPE ?= $(shell uname -m)

CC = gcc

CFLAGS = -w 

CFLAGS += -O
CFLAGS += -D_FILE_OFFSET_BITS=64
CFLAGS += -D_LARGEFILE_SOURCE
CFLAGS += -D_GNU_SOURCE
CFLAGS += -DMACHTYPE_${MACHTYPE}

CFLAGS += -I../inc

OBJ := $(patsubst %.c,%.o, $(wildcard *.c))

foo.a: $(OBJ)
	if [ ! -d $(MACHTYPE) ]; then mkdir $(MACHTYPE); fi
	ar rcus $(MACHTYPE)/jkweb.a $(OBJ)

clean:
	if [ -d $(MACHTYPE) ]; then rm -rf $(MACHTYPE); fi
	rm -f $(OBJ)

.PHONY: clean

