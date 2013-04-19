#==================================================================================================
#
#  Makefile for the smatdvecmult module of the Blaze test suite
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
MCaV2a: MCaV2a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MCaV2b: MCaV2b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MCaV3a: MCaV3a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MCaV3b: MCaV3b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MCaV4a: MCaV4a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MCaV4b: MCaV4b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MCaV5a: MCaV5a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MCaV5b: MCaV5b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MCaV6a: MCaV6a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MCaV6b: MCaV6b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MCaVDa: MCaVDa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MCaVDb: MCaVDb.o
	@$(CXX) -o $@ $< $(LIBRARIES)

MCbV2a: MCbV2a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MCbV2b: MCbV2b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MCbV3a: MCbV3a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MCbV3b: MCbV3b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MCbV4a: MCbV4a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MCbV4b: MCbV4b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MCbV5a: MCbV5a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MCbV5b: MCbV5b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MCbV6a: MCbV6a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MCbV6b: MCbV6b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MCbVDa: MCbVDa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
MCbVDb: MCbVDb.o
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
