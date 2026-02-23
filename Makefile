all := libgeo
LDLIBS-libgeo := -lqsys -lqmap
CFLAGS := -g

test:
	$(MAKE) -C tests test

-include ../mk/include.mk
