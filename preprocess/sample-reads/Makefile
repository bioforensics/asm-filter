CC=gcc
CFLAGS=-g -Wall -O2 -Wno-unused-function

all:sample-reads

sample-reads:sample-reads.c pcg_basic.h ksort.h kseq.h
		$(CC) $(CFLAGS) sample-reads.c pcg_basic.c -o $@ -lz -lm

clean:
		rm -fr gmon.out *.o ext/*.o a.out sample-reads *~ *.a *.dSYM session*
