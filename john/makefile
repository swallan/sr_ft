qprobs.o : qprob.f
	gfortran -c qprob.f
	ld -r qprob.o -o qprobs.o
	rm qprob.o

# C source is also available.  The following is the makefile for this:
# CC=gcc
# qprobs.o : qprob.c pow_dd.c qprob.o pow_dd.o
#	gcc -c -munaligned-doubles qprob.c pow_dd.c
#	ld -r qprob.o pow_dd.o -lm -o qprobs.o
#	rm qprob.o pow_dd.o


