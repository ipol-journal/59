CC = gcc
CFLAGS	= -ansi -Wall -Wextra -O3 -fopenmp 
LDFLAGS = -L/usr/lib -lpng -lm -fopenmp

BIN	= efros_leung
INCDIR = ../inc/
CSRC	 = efros_leung_cmd.c efros_leung.c synthesis.c auxiliary.c \
				 io_png.c mt.c \

OBJ	= $(CSRC:.c=.o)

all: $(BIN)


efros_leung_cmd.o: $(INCDIR)/efros_leung.h
efros_leung.o: $(INCDIR)/synthesis.h $(INCDIR)/auxiliary.h
synthesis.o: $(INCDIR)/mt.h
auxiliary.o: $(INCDIR)/io_png.h


%.o	: %.c
	$(CC) $(CFLAGS) -I$(INCDIR) -c -o $@ $<

$(BIN) : $(OBJ)
	$(CC) -o $@  $^ $(LDFLAGS)

.PHONY: all clean test

clean	:
	$(RM) $(OBJ)


