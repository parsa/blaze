#==================================================================================================
#
#  Makefile for the dmatdmatmult module of the Blaze test suite
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
default: M3x3aM3x3a M3x3aM3x3b M3x3bM3x3a M3x3bM3x3b \
         MHaMHa MHaMHb MHbMHa MHbMHb \
         MDaMDa MDaMDb MDbMDa MDbMDb \
         SDaSDa SDaSDb SDbSDa SDbSDb \
         AliasingTest
all: $(BIN)
essential: M3x3aM3x3a MHaMHa MDaMDa SDaSDa AliasingTest
single: MDaMDa


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

M5x7aM7x5a: M5x7aM7x5a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M5x7aM7x5b: M5x7aM7x5b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M5x7aM7x9a: M5x7aM7x9a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M5x7aM7x9b: M5x7aM7x9b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M5x7aMDa: M5x7aMDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M5x7aMDb: M5x7aMDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M5x7bM7x5a: M5x7bM7x5a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M5x7bM7x5b: M5x7bM7x5b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M5x7bM7x9a: M5x7bM7x9a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M5x7bM7x9b: M5x7bM7x9b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M5x7bMDa: M5x7bMDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M5x7bMDb: M5x7bMDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

M6x6aM6x6a: M6x6aM6x6a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M6x6aM6x6b: M6x6aM6x6b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M6x6aMDa: M6x6aMDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M6x6aMDb: M6x6aMDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M6x6aMHa: M6x6aMHa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M6x6aMHb: M6x6aMHb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M6x6bM6x6a: M6x6bM6x6a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M6x6bM6x6b: M6x6bM6x6b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M6x6bMDa: M6x6bMDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M6x6bMDb: M6x6bMDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M6x6bMHa: M6x6bMHa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M6x6bMHb: M6x6bMHb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

M9x7aM7x5a: M9x7aM7x5a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M9x7aM7x5b: M9x7aM7x5b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M9x7aM7x9a: M9x7aM7x9a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M9x7aM7x9b: M9x7aM7x9b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M9x7aMDa: M9x7aMDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M9x7aMDb: M9x7aMDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M9x7bM7x5a: M9x7bM7x5a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M9x7bM7x5b: M9x7bM7x5b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M9x7bM7x9a: M9x7bM7x9a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M9x7bM7x9b: M9x7bM7x9b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M9x7bMDa: M9x7bMDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M9x7bMDb: M9x7bMDb.o
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
MDaM7x5a: MDaM7x5a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDaM7x5b: MDaM7x5b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDaM7x9a: MDaM7x9a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDaM7x9b: MDaM7x9b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDaMDa: MDaMDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDaMDb: MDaMDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDaMHa: MDaMHa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDaMHb: MDaMHb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDaSDa: MDaSDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDaSDb: MDaSDb.o
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
MDbM7x5a: MDbM7x5a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDbM7x5b: MDbM7x5b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDbM7x9a: MDbM7x9a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDbM7x9b: MDbM7x9b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDbMDa: MDbMDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDbMDb: MDbMDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDbMHa: MDbMHa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDbMHb: MDbMHb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDbSDa: MDbSDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDbSDb: MDbSDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

MHaM6x6a: MHaM6x6a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MHaM6x6b: MHaM6x6b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MHaMDa: MHaMDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MHaMDb: MHaMDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MHaMHa: MHaMHa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MHaMHb: MHaMHb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MHbM6x6a: MHbM6x6a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MHbM6x6b: MHbM6x6b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MHbMDa: MHbMDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MHbMDb: MHbMDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MHbMHa: MHbMHa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MHbMHb: MHbMHb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

S3x3aS3x3a: S3x3aS3x3a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
S3x3aS3x3b: S3x3aS3x3b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
S3x3bS3x3a: S3x3bS3x3a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
S3x3bS3x3b: S3x3bS3x3b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

SDaMDa: SDaMDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SDaMDb: SDaMDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SDaSDa: SDaSDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SDaSDb: SDaSDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SDbMDa: SDbMDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SDbMDb: SDbMDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SDbSDa: SDbSDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SDbSDb: SDbSDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

SHaSHa: SHaSHa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SHaSHb: SHaSHb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SHbSHa: SHbSHa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SHbSHb: SHbSHb.o
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
.PHONY: default all essential single clean
