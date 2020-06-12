BUILDDIR=build
SRCDIR=src
CC=clang
OPT_LEVEL=-O3
PROFILE_OPTIONS=
DEBUG_MODE=
PARALLEL_MODE=
CFLAGS=-std=c11 -Wall $(OPT_LEVEL) $(PROFILE_OPTIONS) $(PARALLEL_MODE)
LDFLAGS=-lm
SOURCES=lab3.c
OBJECTS=$(addprefix $(BUILDDIR)/,$(SOURCES:.c=.o))
EXECUTABLE=lab3

all: $(EXECUTABLE)

$(BUILDDIR):
	mkdir -p $@

$(EXECUTABLE): $(OBJECTS) $(BUILDDIR)
	$(CC) $(CFLAGS) $(OBJECTS) -o $@ $(LDFLAGS)

$(BUILDDIR)/%.o: $(SRCDIR)/%.c
	$(CC) $(CFLAGS) -I./$(SRCDIR) $(DEBUG_MODE) $^ -c -o $@
	
clean:
	rm -rf $(OBJECTS) $(EXECUTABLE)

clear: clean

rebuild: clean all
