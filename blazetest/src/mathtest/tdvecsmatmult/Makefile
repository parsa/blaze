#==================================================================================================
#
#  Makefile for the tdvecsmatmult module of the Blaze test suite
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
essential: V3aMCa VHaMCa VDaMCa AliasingTest


# Build rules
V2aMCa: V2aMCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V2aMCb: V2aMCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V2bMCa: V2bMCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V2bMCb: V2bMCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

V3aMCa: V3aMCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V3aMCb: V3aMCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V3bMCa: V3bMCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V3bMCb: V3bMCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

V4aMCa: V4aMCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V4aMCb: V4aMCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V4bMCa: V4bMCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V4bMCb: V4bMCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

V5aMCa: V5aMCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V5aMCb: V5aMCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V5bMCa: V5bMCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V5bMCb: V5bMCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

V6aMCa: V6aMCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V6aMCb: V6aMCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V6bMCa: V6bMCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V6bMCb: V6bMCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

VDaMCa: VDaMCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDaMCb: VDaMCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDbMCa: VDbMCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDbMCb: VDbMCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

VHaMCa: VHaMCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VHaMCb: VHaMCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VHbMCa: VHbMCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VHbMCb: VHbMCb.o
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
