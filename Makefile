# Makefile
# Nicolas ROUGON | Telecom SudParis / ARTEMIS Department - 15/06/2015

####### Compiler, tools and options

CC = gcc

CFLAGS  = -Wall
INCPATH =
LINK    = gcc
LIBS    = -lm

SOURCES = pdefilter.c utilities.c
OBJECTS = pdefilter.o utilities.o
TARGET  = pdefilter.exe

####### Implicit rules

.SUFFIXES: .c

.c.o:
	$(CC) -c $(CFLAGS) $(INCPATH) -o $@ $<

####### Build rules

all: Makefile $(TARGET)

$(TARGET): $(OBJECTS)
	$(LINK) -o $(TARGET) $(OBJECTS) $(LIBS)

clean:
	@ rm -f *% *.o core pdefilter resultat.???*
