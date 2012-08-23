#==================================================================================================
#
#  Makefile for the smatdmatsub module of the Blaze test suite
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
MCaM2x2a: MCaM2x2a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MCaM2x2b: MCaM2x2b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MCaM3x3a: MCaM3x3a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MCaM3x3b: MCaM3x3b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MCaM4x4a: MCaM4x4a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MCaM4x4b: MCaM4x4b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MCaM5x5a: MCaM5x5a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MCaM5x5b: MCaM5x5b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MCaM6x6a: MCaM6x6a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MCaM6x6b: MCaM6x6b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MCaM7x13a: MCaM7x13a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MCaM7x13b: MCaM7x13b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MCaM16x8a: MCaM16x8a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MCaM16x8b: MCaM16x8b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MCaMDa: MCaMDa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MCaMDb: MCaMDb.o
	@$(CXX) -o $@ $< $(LIBRARIES)

MCbM2x2a: MCbM2x2a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MCbM2x2b: MCbM2x2b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MCbM3x3a: MCbM3x3a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MCbM3x3b: MCbM3x3b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MCbM4x4a: MCbM4x4a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MCbM4x4b: MCbM4x4b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MCbM5x5a: MCbM5x5a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MCbM5x5b: MCbM5x5b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MCbM6x6a: MCbM6x6a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MCbM6x6b: MCbM6x6b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MCbM7x13a: MCbM7x13a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MCbM7x13b: MCbM7x13b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MCbM16x8a: MCbM16x8a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MCbM16x8b: MCbM16x8b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MCbMDa: MCbMDa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MCbMDb: MCbMDb.o
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
