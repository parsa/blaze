#==================================================================================================
#
#  Makefile for the dmatdmatsub module of the Blaze test suite
#
#  Copyright (C) 2013 Klaus Iglberger - All Rights Reserved - All Rights Reserved
#
#  This file is part of the Blaze library. You can redistribute it and/or modify it under
#  the terms of the New (Revised) BSD License. Redistribution and use in source and binary
#  forms, with or without modification, are permitted provided that the following conditions
#  are met:
#
#  1. Redistributions of source code must retain the above copyright notice, this list of
#     conditions and the following disclaimer.
#  2. Redistributions in binary form must reproduce the above copyright notice, this list
#     of conditions and the following disclaimer in the documentation and/or other materials
#     provided with the distribution.
#  3. Neither the names of the Blaze development group nor the names of its contributors
#     may be used to endorse or promote products derived from this software without specific
#     prior written permission.
#
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
#  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
#  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
#  SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
#  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
#  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
#  BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
#  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
#  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
#  DAMAGE.
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


# General rules
default: all
all: $(BIN)
essential: M3x3aM3x3a MDaMDa AliasingTest


# Build rules
M2x2aM2x2a: M2x2aM2x2a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M2x2aM2x2b: M2x2aM2x2b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M2x2aMDa: M2x2aMDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M2x2aMDb: M2x2aMDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M2x2bM2x2a: M2x2bM2x2a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M2x2bM2x2b: M2x2bM2x2b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M2x2bMDa: M2x2bMDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M2x2bMDb: M2x2bMDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

M3x3aM3x3a: M3x3aM3x3a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M3x3aM3x3b: M3x3aM3x3b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M3x3aMDa: M3x3aMDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M3x3aMDb: M3x3aMDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M3x3bM3x3a: M3x3bM3x3a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M3x3bM3x3b: M3x3bM3x3b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M3x3bMDa: M3x3bMDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M3x3bMDb: M3x3bMDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

M4x4aM4x4a: M4x4aM4x4a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M4x4aM4x4b: M4x4aM4x4b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M4x4aMDa: M4x4aMDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M4x4aMDb: M4x4aMDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M4x4bM4x4a: M4x4bM4x4a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M4x4bM4x4b: M4x4bM4x4b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M4x4bMDa: M4x4bMDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M4x4bMDb: M4x4bMDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

M5x5aM5x5a: M5x5aM5x5a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M5x5aM5x5b: M5x5aM5x5b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M5x5aMDa: M5x5aMDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M5x5aMDb: M5x5aMDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M5x5bM5x5a: M5x5bM5x5a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M5x5bM5x5b: M5x5bM5x5b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M5x5bMDa: M5x5bMDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M5x5bMDb: M5x5bMDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

M6x6aM6x6a: M6x6aM6x6a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M6x6aM6x6b: M6x6aM6x6b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M6x6aMDa: M6x6aMDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M6x6aMDb: M6x6aMDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M6x6bM6x6a: M6x6bM6x6a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M6x6bM6x6b: M6x6bM6x6b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M6x6bMDa: M6x6bMDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M6x6bMDb: M6x6bMDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

M7x13aM7x13a: M7x13aM7x13a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M7x13aM7x13b: M7x13aM7x13b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M7x13aMDa: M7x13aMDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M7x13aMDb: M7x13aMDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M7x13bM7x13a: M7x13bM7x13a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M7x13bM7x13b: M7x13bM7x13b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M7x13bMDa: M7x13bMDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M7x13bMDb: M7x13bMDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

M16x8aM16x8a: M16x8aM16x8a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M16x8aM16x8b: M16x8aM16x8b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M16x8aMDa: M16x8aMDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M16x8aMDb: M16x8aMDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M16x8bM16x8a: M16x8bM16x8a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M16x8bM16x8b: M16x8bM16x8b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M16x8bMDa: M16x8bMDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M16x8bMDb: M16x8bMDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

MDaM2x2a: MDaM2x2a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDaM2x2b: MDaM2x2b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDaM3x3a: MDaM3x3a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDaM3x3b: MDaM3x3b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDaM4x4a: MDaM4x4a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDaM4x4b: MDaM4x4b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDaM5x5a: MDaM5x5a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDaM5x5b: MDaM5x5b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDaM6x6a: MDaM6x6a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDaM6x6b: MDaM6x6b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDaM7x13a: MDaM7x13a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDaM7x13b: MDaM7x13b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDaM16x8a: MDaM16x8a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDaM16x8b: MDaM16x8b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDaMDa: MDaMDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDaMDb: MDaMDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDbM2x2a: MDbM2x2a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDbM2x2b: MDbM2x2b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDbM3x3a: MDbM3x3a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDbM3x3b: MDbM3x3b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDbM4x4a: MDbM4x4a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDbM4x4b: MDbM4x4b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDbM5x5a: MDbM5x5a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDbM5x5b: MDbM5x5b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDbM6x6a: MDbM6x6a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDbM6x6b: MDbM6x6b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDbM7x13a: MDbM7x13a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDbM7x13b: MDbM7x13b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDbM16x8a: MDbM16x8a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDbM16x8b: MDbM16x8b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDbMDa: MDbMDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDbMDb: MDbMDb.o
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
.PHONY: default all essential clean
