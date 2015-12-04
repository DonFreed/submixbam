CC       = gcc
CPPFLAGS = 
CFLAGS   = -g -Wall -O2
LDFLAGS  =
LIBS     =

INCLUDES = -I . -I $(HTSDIR)

HTSDIR = /usr/local/lib
HTSLIB = $(HTSDIR)/libhts.a

all: submixbam

.SUFFIXES: .c .o

%.o : %.c
	$(CC) -c $(CFLAGS) $(INCLUDES) $< -o $@

submixbam: submixbam.o bedidx.o $(HTSLIB)
	$(CC) -pthread $(LDFLAGS) -o $@ submixbam.o bedidx.o $(HTSLIB) $(LIBS) -lz

submixbam.o: submixbam.c
bedidx.o: bedidx.c $(HTSDIR)/htslib/ksort.h $(HTSDIR)/htslib/kseq.h $(HTSDIR)/htslib/khash.h

clean:
	-rm -f *.o submixbam
