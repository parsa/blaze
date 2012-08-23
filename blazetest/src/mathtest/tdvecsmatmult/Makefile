#==================================================================================================
#
#  Makefile for the tdvecsmatmult module of the Blaze test suite
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
V2aMCa: V2aMCa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
V2aMCb: V2aMCb.o
	@$(CXX) -o $@ $< $(LIBRARIES)
V2bMCa: V2bMCa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
V2bMCb: V2bMCb.o
	@$(CXX) -o $@ $< $(LIBRARIES)

V3aMCa: V3aMCa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
V3aMCb: V3aMCb.o
	@$(CXX) -o $@ $< $(LIBRARIES)
V3bMCa: V3bMCa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
V3bMCb: V3bMCb.o
	@$(CXX) -o $@ $< $(LIBRARIES)

V4aMCa: V4aMCa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
V4aMCb: V4aMCb.o
	@$(CXX) -o $@ $< $(LIBRARIES)
V4bMCa: V4bMCa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
V4bMCb: V4bMCb.o
	@$(CXX) -o $@ $< $(LIBRARIES)

V5aMCa: V5aMCa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
V5aMCb: V5aMCb.o
	@$(CXX) -o $@ $< $(LIBRARIES)
V5bMCa: V5bMCa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
V5bMCb: V5bMCb.o
	@$(CXX) -o $@ $< $(LIBRARIES)

V6aMCa: V6aMCa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
V6aMCb: V6aMCb.o
	@$(CXX) -o $@ $< $(LIBRARIES)
V6bMCa: V6bMCa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
V6bMCb: V6bMCb.o
	@$(CXX) -o $@ $< $(LIBRARIES)

VDaMCa: VDaMCa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VDaMCb: VDaMCb.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VDbMCa: VDbMCa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VDbMCb: VDbMCb.o
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
