#==================================================================================================
#
#  Makefile for the dmatdvecmult module of the Blaze test suite
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
essential: M3x3aV3a MDaVHa MDaVDa AliasingTest


# Build rules
M2x2aV2a: M2x2aV2a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M2x2aV2b: M2x2aV2b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M2x2aVDa: M2x2aVDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M2x2aVDb: M2x2aVDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M2x2aVHa: M2x2aVHa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M2x2aVHb: M2x2aVHb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M2x2bV2a: M2x2bV2a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M2x2bV2b: M2x2bV2b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M2x2bVDa: M2x2bVDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M2x2bVDb: M2x2bVDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M2x2bVHa: M2x2bVHa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M2x2bVHb: M2x2bVHb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

M3x3aV3a: M3x3aV3a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M3x3aV3b: M3x3aV3b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M3x3aVDa: M3x3aVDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M3x3aVDb: M3x3aVDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M3x3aVHa: M3x3aVHa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M3x3aVHb: M3x3aVHb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M3x3bV3a: M3x3bV3a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M3x3bV3b: M3x3bV3b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M3x3bVDa: M3x3bVDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M3x3bVDb: M3x3bVDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M3x3bVHa: M3x3bVHa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M3x3bVHb: M3x3bVHb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

M4x4aV4a: M4x4aV4a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M4x4aV4b: M4x4aV4b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M4x4aVDa: M4x4aVDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M4x4aVDb: M4x4aVDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M4x4aVHa: M4x4aVHa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M4x4aVHb: M4x4aVHb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M4x4bV4a: M4x4bV4a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M4x4bV4b: M4x4bV4b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M4x4bVDa: M4x4bVDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M4x4bVDb: M4x4bVDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M4x4bVHa: M4x4bVHa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M4x4bVHb: M4x4bVHb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

M5x5aV5a: M5x5aV5a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M5x5aV5b: M5x5aV5b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M5x5aVDa: M5x5aVDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M5x5aVDb: M5x5aVDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M5x5aVHa: M5x5aVHa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M5x5aVHb: M5x5aVHb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M5x5bV5a: M5x5bV5a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M5x5bV5b: M5x5bV5b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M5x5bVDa: M5x5bVDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M5x5bVDb: M5x5bVDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M5x5bVHa: M5x5bVHa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M5x5bVHb: M5x5bVHb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

M6x6aV6a: M6x6aV6a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M6x6aV6b: M6x6aV6b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M6x6aVDa: M6x6aVDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M6x6aVDb: M6x6aVDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M6x6aVHa: M6x6aVHa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M6x6aVHb: M6x6aVHb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M6x6bV6a: M6x6bV6a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M6x6bV6b: M6x6bV6b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M6x6bVDa: M6x6bVDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M6x6bVDb: M6x6bVDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M6x6bVHa: M6x6bVHa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M6x6bVHb: M6x6bVHb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

M7x13aV13a: M7x13aV13a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M7x13aV13b: M7x13aV13b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M7x13aVDa: M7x13aVDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M7x13aVDb: M7x13aVDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M7x13aVHa: M7x13aVHa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M7x13aVHb: M7x13aVHb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M7x13bV13a: M7x13bV13a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M7x13bV13b: M7x13bV13b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M7x13bVDa: M7x13bVDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M7x13bVDb: M7x13bVDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M7x13bVHa: M7x13bVHa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M7x13bVHb: M7x13bVHb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

M16x8aV8a: M16x8aV8a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M16x8aV8b: M16x8aV8b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M16x8aVDa: M16x8aVDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M16x8aVDb: M16x8aVDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M16x8aVHa: M16x8aVHa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M16x8aVHb: M16x8aVHb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M16x8bV8a: M16x8bV8a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M16x8bV8b: M16x8bV8b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M16x8bVDa: M16x8bVDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M16x8bVDb: M16x8bVDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M16x8bVHa: M16x8bVHa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M16x8bVHb: M16x8bVHb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

MDaV2a: MDaV2a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDaV2b: MDaV2b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDaV3a: MDaV3a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDaV3b: MDaV3b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDaV4a: MDaV4a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDaV4b: MDaV4b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDaV5a: MDaV5a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDaV5b: MDaV5b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDaV6a: MDaV6a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDaV6b: MDaV6b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDaVDa: MDaVDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDaVDb: MDaVDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDaVHa: MDaVHa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDaVHb: MDaVHb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDbV2a: MDbV2a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDbV2b: MDbV2b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDbV3a: MDbV3a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDbV3b: MDbV3b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDbV4a: MDbV4a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDbV4b: MDbV4b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDbV5a: MDbV5a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDbV5b: MDbV5b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDbV6a: MDbV6a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDbV6b: MDbV6b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDbVDa: MDbVDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDbVDb: MDbVDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDbVHa: MDbVHa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDbVHb: MDbVHb.o
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
