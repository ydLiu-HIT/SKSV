CC=	gcc
CFLAGS=     -g -Wc++-compat -Wall -O3 -Wno-unused-variable -Wno-unused-but-set-variable -Wno-unused-function
CPPFLAGS= 
INCLUDES= 
OBJS=		ktime.o bseq.o main.o var_aligner.o form.o index_loader.o bit_operation.o sdp.o lv.o
PROG=		SKSV-skeleton	
LIBS=		-lm -lz -lpthread -lgomp



.PHONY:all extra clean depend
.SUFFIXES:.c .o


.c.o:
		$(CC) -c $(CFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

SKSV-skeleton:$(OBJS)
	$(CC) $(CFLAGS) $(OBJS) $(LIBS) -o $@

clean:
		@rm -f *.o $(PROG)

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) -- *.c)

# # DO NOT DELETE
bseq.o: bseq.h kseq1.h
main.o: db_index.h var_aligner.h form.h desc.h
var_aligner.o: form.h var_aligner.h desc.h index_loader.h bit_operation.h lv.h sdp.h bseq.h ktime.h
form.o: form.h desc.h
index_loader.o: index_loader.h
bit_operation.o: bit_operation.h
sdp.o: sdp.h bit_operation.h lv.h
lv.o: bit_operation.h
ktime.o: ktime.h
