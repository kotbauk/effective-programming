BUILDDIR=build
SRCDIR=src
CC=gcc
OPT_LEVEL=
PROFILE_OPTIONS=
DEBUG_MODE=
CFLAGS=-std=c11 -Wall -fno-omit-frame-pointer $(OPT_LEVEL) $(PROFILE_OPTIONS)
LDFLAGS=-lm
HEADERS=my_matrix.h
MATRIX_SOURCES=my_matrix.c
SOURCES=main.c $(MATRIX_SOURCES)
OBJECTS=$(addprefix $(BUILDDIR)/,$(SOURCES:.c=.o))
EXECUTABLE=lab1

all: $(EXECUTABLE)

$(BUILDDIR):
	mkdir -p $@

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) -o $@ $(LDFLAGS)

$(BUILDDIR)/%.o: $(SRCDIR)/%.c  
	$(CC) $(CFLAGS) -I./$(SRCDIR) $(DEBUG_MODE) $^ -c -o $@ $(LDFLAGS) 
	
clean:
	rm -rf $(OBJECTS) $(EXECUTABLE)

clear: clean

rebuild: clean all
