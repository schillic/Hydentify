CC=mex
LDFLAGS=-lgmp -lgmpxx -lppl -cxx
CFLAGS=-g -I'/usr/local/include/' -L'/usr/local/lib'
OBJ=ppl_mexif.o ppl_wrap.o

#mex -c ppl_mexif.cpp  ppl_wrap.cpp - -lgmp -lgmpxx -lppl 

all: pplmex

.PHONY: all clean

pplmex: $(OBJ)
	$(CC) $(LDFLAGS) -o $@ $^

%.o: %.cpp
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	rm -f pplmex.* *.o
