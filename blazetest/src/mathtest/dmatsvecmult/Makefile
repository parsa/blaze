#==================================================================================================
#
#  Makefile for the dmatsvecmult module of the Blaze test suite
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
M2x2aVCa: M2x2aVCa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M2x2aVCb: M2x2aVCb.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M2x2bVCa: M2x2bVCa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M2x2bVCb: M2x2bVCb.o
	@$(CXX) -o $@ $< $(LIBRARIES)

M3x3aVCa: M3x3aVCa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M3x3aVCb: M3x3aVCb.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M3x3bVCa: M3x3bVCa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M3x3bVCb: M3x3bVCb.o
	@$(CXX) -o $@ $< $(LIBRARIES)

M4x4aVCa: M4x4aVCa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M4x4aVCb: M4x4aVCb.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M4x4bVCa: M4x4bVCa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M4x4bVCb: M4x4bVCb.o
	@$(CXX) -o $@ $< $(LIBRARIES)

M5x5aVCa: M5x5aVCa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M5x5aVCb: M5x5aVCb.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M5x5bVCa: M5x5bVCa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M5x5bVCb: M5x5bVCb.o
	@$(CXX) -o $@ $< $(LIBRARIES)

M6x6aVCa: M6x6aVCa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M6x6aVCb: M6x6aVCb.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M6x6bVCa: M6x6bVCa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M6x6bVCb: M6x6bVCb.o
	@$(CXX) -o $@ $< $(LIBRARIES)

M7x13aVCa: M7x13aVCa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M7x13aVCb: M7x13aVCb.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M7x13bVCa: M7x13bVCa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M7x13bVCb: M7x13bVCb.o
	@$(CXX) -o $@ $< $(LIBRARIES)

M16x8aVCa: M16x8aVCa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M16x8aVCb: M16x8aVCb.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M16x8bVCa: M16x8bVCa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
M16x8bVCb: M16x8bVCb.o
	@$(CXX) -o $@ $< $(LIBRARIES)

MDaVCa: MDaVCa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MDaVCb: MDaVCb.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MDbVCa: MDbVCa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MDbVCb: MDbVCb.o
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
