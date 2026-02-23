# Examples Makefile for libgeo

CC := cc
CFLAGS := -g -Wall -Wextra -I../include
LDFLAGS := -L../lib -L/usr/lib
LDLIBS := -lgeo -lqmap -lqsys -lxxhash

EXAMPLES := basic iteration persistence 3d_world

all: $(EXAMPLES)

%: %.c
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS) $(LDLIBS)

clean:
	rm -f $(EXAMPLES) *.o

.PHONY: all clean
