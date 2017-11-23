### snapshot options #######
#EXTRAS += -DLEN_MANUEL
EXTRAS += -DPRINT_XYZ

#CC
CC     := $(OMPP) gcc $(DOMPP)
DC     := -DNTHREADS=4
#DC     += -DLOCK
CFLAGS := -Wall -O3 -fopenmp -g
GSLL   := -lgsl -lgslcblas
LIBS   := -lm $(GSLL)

.PHONY : cleanall clean todo 

MAKEFILE := Makefile

OBJS := grid.o variables.o leesloanreal.o iden.o

HEADERS := $(patsubst %.o,$.h,$(OBJS))

EXEC := mendieta.x

todo: $(EXEC)

%.o: %.c %.h $(MAKEFILE)
	$(CC) $(EXTRAS) $(CFLAGS) $(DC) -c $<

mendieta.x: mendieta.c $(OBJS)
	$(CC) $(CFLAGS) $(EXTRAS) $^  -o $@ $(LIBS)

clean:
	rm -rf $(OBJS)
	rm -rf mendieta.o
	rm -rf $(EXEC)
