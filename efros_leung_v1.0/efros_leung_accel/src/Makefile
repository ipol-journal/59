CC = gcc
CFLAGS	= -ansi -Wall -Wextra -Werror -O3 -march=native -fopenmp 
LDFLAGS = -L/usr/lib -lpng -lm -lblas -llapack -fopenmp 

BIN	= efrosLeung_accel
INCDIR = ../inc/
CSRC	 = efros_leung_cmd.c efros_leung.c synthesis.c auxiliary.c \
         pca.c svd.c linear_algebra.c \
				 io_png.c mt.c \

OBJ	= $(CSRC:.c=.o)

all: $(BIN)


efros_leung_cmd.o: $(INCDIR)/efros_leung.h
efros_leung.o: $(INCDIR)/synthesis.h $(INCDIR)/auxiliary.h
synthesis.o: $(INCDIR)/pca.h $(INCDIR)/mt.h
auxiliary.o: $(INCDIR)/io_png.h
pca.o:$(INCDIR)/auxiliary.h $(INCDIR)/linear_algebra.h $(INCDIR)/svd.h


%.o	: %.c
	$(CC) $(CFLAGS) -I$(INCDIR) -c -o $@ $<

$(BIN) : $(OBJ)
	$(CC) -o $@  $^ $(LDFLAGS)

.PHONY: all clean test

clean	:
	$(RM) $(OBJ)


