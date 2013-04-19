#==================================================================================================
#
#  Makefile for the tdvecdmatmult module of the Blaze test suite
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
V2aM2x2a: V2aM2x2a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
V2aM2x2b: V2aM2x2b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
V2aMDa: V2aMDa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
V2aMDb: V2aMDb.o
	@$(CXX) -o $@ $< $(LIBRARIES)
V2bM2x2a: V2bM2x2a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
V2bM2x2b: V2bM2x2b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
V2bMDa: V2bMDa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
V2bMDb: V2bMDb.o
	@$(CXX) -o $@ $< $(LIBRARIES)

V3aM3x3a: V3aM3x3a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
V3aM3x3b: V3aM3x3b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
V3aMDa: V3aMDa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
V3aMDb: V3aMDb.o
	@$(CXX) -o $@ $< $(LIBRARIES)
V3bM3x3a: V3bM3x3a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
V3bM3x3b: V3bM3x3b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
V3bMDa: V3bMDa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
V3bMDb: V3bMDb.o
	@$(CXX) -o $@ $< $(LIBRARIES)

V4aM4x4a: V4aM4x4a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
V4aM4x4b: V4aM4x4b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
V4aMDa: V4aMDa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
V4aMDb: V4aMDb.o
	@$(CXX) -o $@ $< $(LIBRARIES)
V4bM4x4a: V4bM4x4a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
V4bM4x4b: V4bM4x4b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
V4bMDa: V4bMDa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
V4bMDb: V4bMDb.o
	@$(CXX) -o $@ $< $(LIBRARIES)

V5aM5x5a: V5aM5x5a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
V5aM5x5b: V5aM5x5b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
V5aMDa: V5aMDa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
V5aMDb: V5aMDb.o
	@$(CXX) -o $@ $< $(LIBRARIES)
V5bM5x5a: V5bM5x5a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
V5bM5x5b: V5bM5x5b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
V5bMDa: V5bMDa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
V5bMDb: V5bMDb.o
	@$(CXX) -o $@ $< $(LIBRARIES)

V6aM6x6a: V6aM6x6a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
V6aM6x6b: V6aM6x6b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
V6aMDa: V6aMDa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
V6aMDb: V6aMDb.o
	@$(CXX) -o $@ $< $(LIBRARIES)
V6bM6x6a: V6bM6x6a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
V6bM6x6b: V6bM6x6b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
V6bMDa: V6bMDa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
V6bMDb: V6bMDb.o
	@$(CXX) -o $@ $< $(LIBRARIES)

V7aM7x13a: V7aM7x13a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
V7aM7x13b: V7aM7x13b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
V7bM7x13a: V7bM7x13a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
V7bM7x13b: V7bM7x13b.o
	@$(CXX) -o $@ $< $(LIBRARIES)

V16aM16x8a: V16aM16x8a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
V16aM16x8b: V16aM16x8b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
V16bM16x8a: V16bM16x8a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
V16bM16x8b: V16bM16x8b.o
	@$(CXX) -o $@ $< $(LIBRARIES)

VDaM2x2a: VDaM2x2a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VDaM2x2b: VDaM2x2b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VDaM3x3a: VDaM3x3a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VDaM3x3b: VDaM3x3b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VDaM4x4a: VDaM4x4a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VDaM4x4b: VDaM4x4b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VDaM5x5a: VDaM5x5a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VDaM5x5b: VDaM5x5b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VDaM6x6a: VDaM6x6a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VDaM6x6b: VDaM6x6b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VDaM7x13a: VDaM7x13a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VDaM7x13b: VDaM7x13b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VDaM16x8a: VDaM16x8a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VDaM16x8b: VDaM16x8b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VDaMDa: VDaMDa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VDaMDb: VDaMDb.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VDbM2x2a: VDbM2x2a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VDbM2x2b: VDbM2x2b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VDbM3x3a: VDbM3x3a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VDbM3x3b: VDbM3x3b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VDbM4x4a: VDbM4x4a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VDbM4x4b: VDbM4x4b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VDbM5x5a: VDbM5x5a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VDbM5x5b: VDbM5x5b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VDbM6x6a: VDbM6x6a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VDbM6x6b: VDbM6x6b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VDbM7x13a: VDbM7x13a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VDbM7x13b: VDbM7x13b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VDbM16x8a: VDbM16x8a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VDbM16x8b: VDbM16x8b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VDbMDa: VDbMDa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VDbMDb: VDbMDb.o
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
