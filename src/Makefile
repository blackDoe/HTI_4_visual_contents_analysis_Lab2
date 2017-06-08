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
TARGET  = ../bin/pdefilter

SOURCES2 = permal.c utilities.c
OBJECTS2 = permal.o utilities.o
TARGET2  = ../bin/permal

SOURCES3 = nordstrom.c utilities.c
OBJECTS3 = nordstrom.o utilities.o
TARGET3  = ../bin/nordstrom

SOURCES4 = nordstrom_entropique.c utilities.c
OBJECTS4 = nordstrom_entropique.o utilities.o
TARGET4  = ../bin/nordstrom_entropique

####### Implicit rules

.SUFFIXES: .c

.c.o:
	$(CC) -c $(CFLAGS) $(INCPATH) -o $@ $<

####### Build rules

all: Makefile $(TARGET) $(TARGET2) $(TARGET3) $(TARGET4) # clean

$(TARGET): $(OBJECTS)
	$(LINK) -o $(TARGET) $(OBJECTS) $(LIBS)

$(TARGET2): $(OBJECTS2)
	$(LINK) -o $(TARGET2) $(OBJECTS2) $(LIBS)

$(TARGET3): $(OBJECTS3)
	$(LINK) -o $(TARGET3) $(OBJECTS3) $(LIBS)

$(TARGET4): $(OBJECTS4)
	$(LINK) -o $(TARGET4) $(OBJECTS4) $(LIBS)

clean:
	@ rm -f *% *.o core pdefilter resultat.???*
