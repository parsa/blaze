#==================================================================================================
#
#  Makefile for the svectdvecmult module of the Blaze test suite
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
VCaV2a: VCaV2a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VCaV2b: VCaV2b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VCbV2a: VCbV2a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VCbV2b: VCbV2b.o
	@$(CXX) -o $@ $< $(LIBRARIES)

VCaV3a: VCaV3a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VCaV3b: VCaV3b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VCbV3a: VCbV3a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VCbV3b: VCbV3b.o
	@$(CXX) -o $@ $< $(LIBRARIES)

VCaV4a: VCaV4a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VCaV4b: VCaV4b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VCbV4a: VCbV4a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VCbV4b: VCbV4b.o
	@$(CXX) -o $@ $< $(LIBRARIES)

VCaV5a: VCaV5a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VCaV5b: VCaV5b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VCbV5a: VCbV5a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VCbV5b: VCbV5b.o
	@$(CXX) -o $@ $< $(LIBRARIES)

VCaV6a: VCaV6a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VCaV6b: VCaV6b.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VCbV6a: VCbV6a.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VCbV6b: VCbV6b.o
	@$(CXX) -o $@ $< $(LIBRARIES)

VCaVDa: VCaVDa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VCaVDb: VCaVDb.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VCbVDa: VCbVDa.o
	@$(CXX) -o $@ $< $(LIBRARIES)
VCbVDb: VCbVDb.o
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
