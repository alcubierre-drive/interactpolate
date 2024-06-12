NAME_LIB := libinteractpolate.so

CC := gcc
LD := g++
CXX := g++
STRIP := strip

CXXFLAGS = $(FLAGS) -std=c++17 -fPIC -fopenmp -MMD -Wall -Wextra -pedantic $(DEFINES) $(INCLUDES)
CFLAGS = $(FLAGS) -std=c11 -fPIC -fopenmp -MMD -Wall -Wextra -pedantic $(DEFINES) $(INCLUDES)
LDFLAGS = -fPIC -fopenmp -Wall -Wextra -pedantic $(LIBS)

FLAGS = 
INCLUDES = 
DEFINES = 
LIBS = -lm

SRC_C := $(wildcard *.c)
SRC_CPP := $(wildcard *.cpp)

OBJ_C := $(patsubst %.c,%.c.o,$(SRC_C))
OBJ_CPP := $(patsubst %.cpp,%.cpp.o,$(SRC_CPP))

DEP_C := $(patsubst %.c,%.c.d,$(SRC_C))
DEP_CPP := $(patsubst %.cpp,%.cpp.d,$(SRC_CPP))

HEADERS := $(wildcard *.h *.hpp)

-include Makefile.inc
-include Makefile.local

.PHONY: lib clean types

lib: types $(NAME_LIB)

-include $(DEP_C) $(DEP_CPP)

$(NAME_LIB): $(OBJ_C) $(OBJ_CPP)
	$(LD) -shared $^ $(LDFLAGS) -o $@
	$(STRIP) $(NAME_LIB)

%.c.o: %.c %.h Makefile $(wildcard Makefile.local Makefile.inc)
	$(CC) -c $< -o $@ $(CFLAGS)
%.c.o: %.c Makefile $(wildcard Makefile.local Makefile.inc)
	$(CC) -c $< -o $@ $(CFLAGS)

%.cpp.o: %.cpp %.h Makefile $(wildcard Makefile.local Makefile.inc)
	$(CXX) -c $< -o $@ $(CXXFLAGS)
%.cpp.o: %.cpp Makefile $(wildcard Makefile.local Makefile.inc)
	$(CXX) -c $< -o $@ $(CXXFLAGS)

types: types.vim
types.vim: $(SRC_C) $(SRC_CPP) $(HEADERS) Makefile $(wildcard Makefile.local Makefile.inc)
	@echo "adding ctypes to vim..."
	@ctags --c-kinds=t -o- $(SRC_C) $(SRC_CPP) $(HEADERS) |\
		awk '{print "syntax keyword Type " $$1}' | uniq > $@

clean:
	-@$(RM) -rf $(NAME_LIB) types.vim
	-@$(RM) -rf $(shell find . -type f -name '*.o') $(shell find . -type f -name '*.d')

