#==================================================================================================
#
#  Makefile for the tsvecdmatmult module of the Blaze test suite
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
VCaM2x2a: VCaM2x2a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VCaM2x2b: VCaM2x2b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VCaM3x3a: VCaM3x3a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VCaM3x3b: VCaM3x3b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VCaM4x4a: VCaM4x4a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VCaM4x4b: VCaM4x4b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VCaM5x5a: VCaM5x5a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VCaM5x5b: VCaM5x5b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VCaM6x6a: VCaM6x6a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VCaM6x6b: VCaM6x6b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VCaM7x13a: VCaM7x13a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VCaM7x13b: VCaM7x13b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VCaM16x8a: VCaM16x8a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VCaM16x8b: VCaM16x8b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VCaMDa: VCaMDa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VCaMDb: VCaMDb.o
	@$(CXX) -o $@ $< $(LIBRARIES)

VCbM2x2a: VCbM2x2a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VCbM2x2b: VCbM2x2b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VCbM3x3a: VCbM3x3a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VCbM3x3b: VCbM3x3b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VCbM4x4a: VCbM4x4a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VCbM4x4b: VCbM4x4b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VCbM5x5a: VCbM5x5a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VCbM5x5b: VCbM5x5b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VCbM6x6a: VCbM6x6a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VCbM6x6b: VCbM6x6b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VCbM7x13a: VCbM7x13a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VCbM7x13b: VCbM7x13b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VCbM16x8a: VCbM16x8a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VCbM16x8b: VCbM16x8b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VCbMDa: VCbMDa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VCbMDb: VCbMDb.o
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
