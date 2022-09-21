#
# makefile for testing geom module
# the boost library is assumed to be in default location, like 
# /usr/lib etc.
# we note that geom is based on general
#
GENERAL              = ../../general
NAME                 = test
CC		   	         = g++
CFLAGS       			= -Wall  -g -fopenmp
INCLUDE              = -I../include -I${GENERAL}/include -I${TBBROOT}/include
MACRO                = -DDEBUG 
OBJC                 = excep.o globalinfor.o xcfunc.o filerw.o textread.o \
							  parameterparsing.o element.o molecule.o zmatrix.o geomutil.o main.o
LIB                  = -lboost_filesystem -lboost_system -lboost_thread -ltbb -ltbbmalloc 

ALL: $(OBJC)
	$(CC) $(CFLAGS) -o $(NAME) $(OBJC) $(LIB)

globalinfor.o: ${GENERAL}/src/globalinfor.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

textread.o: ${GENERAL}/src/textread.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

xcfunc.o: ${GENERAL}/src/xcfunc.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

excep.o: ${GENERAL}/src/excep.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

filerw.o: ${GENERAL}/src/filerw.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

parameterparsing.o: ${GENERAL}/src/parameterparsing.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

element.o: ../src/element.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

molecule.o: ../src/molecule.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

zmatrix.o: ../src/zmatrix.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

geomutil.o: ../src/geomutil.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

main.o: main.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

.PHONY: clean
clean:
	@rm *.o
	@rm test.out
	@rm test

