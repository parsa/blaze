#==================================================================================================
#
#  Makefile for the dmatdvecmult module of the Blaze test suite
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
M2x2aV2a: M2x2aV2a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M2x2aV2b: M2x2aV2b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M2x2aVDa: M2x2aVDa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M2x2aVDb: M2x2aVDb.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M2x2bV2a: M2x2bV2a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M2x2bV2b: M2x2bV2b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M2x2bVDa: M2x2bVDa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M2x2bVDb: M2x2bVDb.o
	@$(CXX) -o $@ $< $(LIBRARIES)

M3x3aV3a: M3x3aV3a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M3x3aV3b: M3x3aV3b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M3x3aVDa: M3x3aVDa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M3x3aVDb: M3x3aVDb.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M3x3bV3a: M3x3bV3a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M3x3bV3b: M3x3bV3b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M3x3bVDa: M3x3bVDa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M3x3bVDb: M3x3bVDb.o
	@$(CXX) -o $@ $< $(LIBRARIES)

M4x4aV4a: M4x4aV4a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M4x4aV4b: M4x4aV4b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M4x4aVDa: M4x4aVDa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M4x4aVDb: M4x4aVDb.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M4x4bV4a: M4x4bV4a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M4x4bV4b: M4x4bV4b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M4x4bVDa: M4x4bVDa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M4x4bVDb: M4x4bVDb.o
	@$(CXX) -o $@ $< $(LIBRARIES)

M5x5aV5a: M5x5aV5a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M5x5aV5b: M5x5aV5b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M5x5aVDa: M5x5aVDa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M5x5aVDb: M5x5aVDb.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M5x5bV5a: M5x5bV5a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M5x5bV5b: M5x5bV5b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M5x5bVDa: M5x5bVDa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M5x5bVDb: M5x5bVDb.o
	@$(CXX) -o $@ $< $(LIBRARIES)

M6x6aV6a: M6x6aV6a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M6x6aV6b: M6x6aV6b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M6x6aVDa: M6x6aVDa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M6x6aVDb: M6x6aVDb.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M6x6bV6a: M6x6bV6a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M6x6bV6b: M6x6bV6b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M6x6bVDa: M6x6bVDa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M6x6bVDb: M6x6bVDb.o
	@$(CXX) -o $@ $< $(LIBRARIES)

M7x13aV13a: M7x13aV13a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M7x13aV13b: M7x13aV13b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M7x13aVDa: M7x13aVDa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M7x13aVDb: M7x13aVDb.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M7x13bV13a: M7x13bV13a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M7x13bV13b: M7x13bV13b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M7x13bVDa: M7x13bVDa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M7x13bVDb: M7x13bVDb.o
	@$(CXX) -o $@ $< $(LIBRARIES)

M16x8aV8a: M16x8aV8a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M16x8aV8b: M16x8aV8b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M16x8aVDa: M16x8aVDa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M16x8aVDb: M16x8aVDb.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M16x8bV8a: M16x8bV8a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M16x8bV8b: M16x8bV8b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M16x8bVDa: M16x8bVDa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M16x8bVDb: M16x8bVDb.o
	@$(CXX) -o $@ $< $(LIBRARIES)

MDaV2a: MDaV2a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MDaV2b: MDaV2b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MDaV3a: MDaV3a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MDaV3b: MDaV3b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MDaV4a: MDaV4a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MDaV4b: MDaV4b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MDaV5a: MDaV5a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MDaV5b: MDaV5b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MDaV6a: MDaV6a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MDaV6b: MDaV6b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MDaVDa: MDaVDa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MDaVDb: MDaVDb.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MDbV2a: MDbV2a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MDbV2b: MDbV2b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MDbV3a: MDbV3a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MDbV3b: MDbV3b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MDbV4a: MDbV4a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MDbV4b: MDbV4b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MDbV5a: MDbV5a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MDbV5b: MDbV5b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MDbV6a: MDbV6a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MDbV6b: MDbV6b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MDbVDa: MDbVDa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MDbVDb: MDbVDb.o
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
