LIB-LDLIBS := -lqsys -lqdb
LIB := geo
HEADERS := morton.h point.h
CFLAGS := -g

npm-lib := @tty-pt/qsys @tty-pt/qdb

-include ../mk/include.mk
