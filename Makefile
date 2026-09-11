all := libislet
LDLIBS-libislet := -lqsys -lqmap
CFLAGS := -g

-include ../mk/include.mk

test:
	$(MAKE) -C tests test

bench:
	$(MAKE) -C tests bench
