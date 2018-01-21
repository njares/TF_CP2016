CC=cc

CFLAGS=-Wall -Wextra -std=c99
#CFLAGS=-Wall -Wextra -std=c99 -O2 -ffast-math

LDFLAGS=-lm

all: genera-script

heat: heat.o netpbm.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

.PHONY: clean

clean:
	rm -rf *.o heat
