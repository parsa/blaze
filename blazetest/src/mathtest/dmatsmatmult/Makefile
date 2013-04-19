#==================================================================================================
#
#  Makefile for the dmatsmatmult module of the Blaze test suite
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
M2x2aMCa: M2x2aMCa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M2x2aMCb: M2x2aMCb.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M2x2bMCa: M2x2bMCa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M2x2bMCb: M2x2bMCb.o
	@$(CXX) -o $@ $< $(LIBRARIES)

M3x3aMCa: M3x3aMCa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M3x3aMCb: M3x3aMCb.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M3x3bMCa: M3x3bMCa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M3x3bMCb: M3x3bMCb.o
	@$(CXX) -o $@ $< $(LIBRARIES)

M4x4aMCa: M4x4aMCa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M4x4aMCb: M4x4aMCb.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M4x4bMCa: M4x4bMCa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M4x4bMCb: M4x4bMCb.o
	@$(CXX) -o $@ $< $(LIBRARIES)

M5x5aMCa: M5x5aMCa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M5x5aMCb: M5x5aMCb.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M5x5bMCa: M5x5bMCa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M5x5bMCb: M5x5bMCb.o
	@$(CXX) -o $@ $< $(LIBRARIES)

M6x6aMCa: M6x6aMCa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M6x6aMCb: M6x6aMCb.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M6x6bMCa: M6x6bMCa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M6x6bMCb: M6x6bMCb.o
	@$(CXX) -o $@ $< $(LIBRARIES)

M7x13aMCa: M7x13aMCa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M7x13aMCb: M7x13aMCb.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M7x13bMCa: M7x13bMCa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M7x13bMCb: M7x13bMCb.o
	@$(CXX) -o $@ $< $(LIBRARIES)

M16x8aMCa: M16x8aMCa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M16x8aMCb: M16x8aMCb.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M16x8bMCa: M16x8bMCa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M16x8bMCb: M16x8bMCb.o
	@$(CXX) -o $@ $< $(LIBRARIES)

MDaMCa: MDaMCa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MDaMCb: MDaMCb.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MDbMCa: MDbMCa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MDbMCb: MDbMCb.o
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
