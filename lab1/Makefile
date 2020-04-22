BUILDDIR=build
SRCDIR=src
CC=clang
OPT_LEVEL=
PROFILE_OPTIONS=
DEBUG_MODE=
CFLAGS=-std=c99 -Wall $(OPT_LEVEL) $(PROFILE_OPTIONS)
LDFLAGS=-lm
SOURCES=lab1.c
OBJECTS=$(addprefix $(BUILDDIR)/,$(SOURCES:.c=.o))
EXECUTABLE=lab1

all: $(EXECUTABLE)

$(BUILDDIR):
	mkdir -p $@

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) -o $@ $(LDFLAGS)

$(BUILDDIR)/%.o: $(SRCDIR)/%.c
	$(CC) $(CFLAGS) -I./$(SRCDIR) $(DEBUG_MODE) $^ -c -o $@
	
clean:
	rm -rf $(OBJECTS) $(EXECUTABLE)

clear: clean

rebuild: clean all
