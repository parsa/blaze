#==================================================================================================
#
#  Makefile for the dmatdmatmult module of the Blaze test suite
#
#  Copyright (C) 2011 Klaus Iglberger - All Rights Reserved - All Rights Reserved
#
#  This file is part of the Blaze library. This library is free software; you can redistribute
#  it and/or modify it under the terms of the GNU General Public License as published by the
#  Free Software Foundation; either version 3, or (at your option) any later version.
#
#  This library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
#  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#  See the GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License along with a special
#  exception for linking and compiling against the Blaze library, the so-called "runtime
#  exception"; see the file COPYING. If not, see http://www.gnu.org/licenses/.
#
#==================================================================================================


# Including the compiler and library settings
ifneq ($(MAKECMDGOALS),clean)
-include ../../Makeconfig
endif


# Setting the source, object and dependency files
SRC = $(wildcard ./*.cpp)
DEP = $(SRC:.cpp=.d)
OBJ = $(SRC:.cpp=.o)
BIN = $(SRC:.cpp=)


# Default rule
default: $(BIN)


# Build rules
M2x2aM2x2a: M2x2aM2x2a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M2x2aM2x2b: M2x2aM2x2b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M2x2aMDa: M2x2aMDa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M2x2aMDb: M2x2aMDb.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M2x2bM2x2a: M2x2bM2x2a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M2x2bM2x2b: M2x2bM2x2b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M2x2bMDa: M2x2bMDa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M2x2bMDb: M2x2bMDb.o
	@$(CXX) -o $@ $< $(LIBRARIES)

M3x3aM3x3a: M3x3aM3x3a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M3x3aM3x3b: M3x3aM3x3b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M3x3aMDa: M3x3aMDa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M3x3aMDb: M3x3aMDb.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M3x3bM3x3a: M3x3bM3x3a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M3x3bM3x3b: M3x3bM3x3b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M3x3bMDa: M3x3bMDa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M3x3bMDb: M3x3bMDb.o
	@$(CXX) -o $@ $< $(LIBRARIES)

M4x4aM4x4a: M4x4aM4x4a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M4x4aM4x4b: M4x4aM4x4b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M4x4aMDa: M4x4aMDa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M4x4aMDb: M4x4aMDb.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M4x4bM4x4a: M4x4bM4x4a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M4x4bM4x4b: M4x4bM4x4b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M4x4bMDa: M4x4bMDa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M4x4bMDb: M4x4bMDb.o
	@$(CXX) -o $@ $< $(LIBRARIES)

M5x5aM5x5a: M5x5aM5x5a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M5x5aM5x5b: M5x5aM5x5b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M5x5aMDa: M5x5aMDa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M5x5aMDb: M5x5aMDb.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M5x5bM5x5a: M5x5bM5x5a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M5x5bM5x5b: M5x5bM5x5b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M5x5bMDa: M5x5bMDa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M5x5bMDb: M5x5bMDb.o
	@$(CXX) -o $@ $< $(LIBRARIES)

M5x7aM7x5a: M5x7aM7x5a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M5x7aM7x5b: M5x7aM7x5b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M5x7aM7x9a: M5x7aM7x9a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M5x7aM7x9b: M5x7aM7x9b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M5x7aMDa: M5x7aMDa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M5x7aMDb: M5x7aMDb.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M5x7bM7x5a: M5x7bM7x5a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M5x7bM7x5b: M5x7bM7x5b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M5x7bM7x9a: M5x7bM7x9a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M5x7bM7x9b: M5x7bM7x9b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M5x7bMDa: M5x7bMDa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M5x7bMDb: M5x7bMDb.o
	@$(CXX) -o $@ $< $(LIBRARIES)

M6x6aM6x6a: M6x6aM6x6a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M6x6aM6x6b: M6x6aM6x6b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M6x6aMDa: M6x6aMDa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M6x6aMDb: M6x6aMDb.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M6x6bM6x6a: M6x6bM6x6a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M6x6bM6x6b: M6x6bM6x6b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M6x6bMDa: M6x6bMDa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M6x6bMDb: M6x6bMDb.o
	@$(CXX) -o $@ $< $(LIBRARIES)

M9x7aM7x5a: M9x7aM7x5a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M9x7aM7x5b: M9x7aM7x5b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M9x7aM7x9a: M9x7aM7x9a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M9x7aM7x9b: M9x7aM7x9b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M9x7aMDa: M9x7aMDa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M9x7aMDb: M9x7aMDb.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M9x7bM7x5a: M9x7bM7x5a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M9x7bM7x5b: M9x7bM7x5b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M9x7bM7x9a: M9x7bM7x9a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M9x7bM7x9b: M9x7bM7x9b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M9x7bMDa: M9x7bMDa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M9x7bMDb: M9x7bMDb.o
	@$(CXX) -o $@ $< $(LIBRARIES)

MDaM2x2a: MDaM2x2a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MDaM2x2b: MDaM2x2b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MDaM3x3a: MDaM3x3a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MDaM3x3b: MDaM3x3b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MDaM4x4a: MDaM4x4a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MDaM4x4b: MDaM4x4b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MDaM5x5a: MDaM5x5a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MDaM5x5b: MDaM5x5b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MDaM6x6a: MDaM6x6a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MDaM6x6b: MDaM6x6b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MDaM7x5a: MDaM7x5a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MDaM7x5b: MDaM7x5b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MDaM7x9a: MDaM7x9a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MDaM7x9b: MDaM7x9b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MDaMDa: MDaMDa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MDaMDb: MDaMDb.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MDbM2x2a: MDbM2x2a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MDbM2x2b: MDbM2x2b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MDbM3x3a: MDbM3x3a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MDbM3x3b: MDbM3x3b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MDbM4x4a: MDbM4x4a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MDbM4x4b: MDbM4x4b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MDbM5x5a: MDbM5x5a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MDbM5x5b: MDbM5x5b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MDbM6x6a: MDbM6x6a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MDbM6x6b: MDbM6x6b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MDbM7x5a: MDbM7x5a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MDbM7x5b: MDbM7x5b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MDbM7x9a: MDbM7x9a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MDbM7x9b: MDbM7x9b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MDbMDa: MDbMDa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MDbMDb: MDbMDb.o
	@$(CXX) -o $@ $< $(LIBRARIES)

AliasingTest: AliasingTest.o
	@$(CXX) -o $@ $< $(LIBRARIES)


# Cleanup
clean:
	@$(RM) $(DEP) $(OBJ) $(BIN)


# Makefile includes
ifneq ($(MAKECMDGOALS),clean)
-include $(DEP)
endif


# Makefile generation
%.d: %.cpp
	@$(CXX) -MM -MP -MT "$*.o $*.d" -MF $@ $(CXXFLAGS) $<


# Setting the independent commands
.PHONY: default clean
