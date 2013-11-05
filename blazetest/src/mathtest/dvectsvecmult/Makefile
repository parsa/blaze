#==================================================================================================
#
#  Makefile for the dvectsvecmult module of the Blaze test suite
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
V2aVCa: V2aVCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V2aVCb: V2aVCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V2bVCa: V2bVCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V2bVCb: V2bVCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

V3aVCa: V3aVCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V3aVCb: V3aVCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V3bVCa: V3bVCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V3bVCb: V3bVCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

V4aVCa: V4aVCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V4aVCb: V4aVCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V4bVCa: V4bVCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V4bVCb: V4bVCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

V5aVCa: V5aVCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V5aVCb: V5aVCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V5bVCa: V5bVCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V5bVCb: V5bVCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

V6aVCa: V6aVCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V6aVCb: V6aVCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V6bVCa: V6bVCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V6bVCb: V6bVCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

VDaVCa: VDaVCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDaVCb: VDaVCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDbVCa: VDbVCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDbVCb: VDbVCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

AliasingTest: AliasingTest.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)


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
