#==================================================================================================
#
#  Makefile for the traits module of the Blaze test suite
#
#  Copyright (C) 2012-2020 Klaus Iglberger - All Rights Reserved - All Rights Reserved
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
ifneq ($(MAKECMDGOALS),reset)
ifneq ($(MAKECMDGOALS),clean)
-include ../../Makeconfig
endif
endif


# Setting the source, object and dependency files
SRC = $(wildcard ./*.cpp)
DEP = $(SRC:.cpp=.d)
OBJ = $(SRC:.cpp=.o)
BIN = $(SRC:.cpp=)


# Build rules
default: all

all: $(BIN) addtrait crosstrait decldiagtrait declhermtrait decllowtrait declstrlowtrait \
     declstrupptrait declsymtrait declunilowtrait decluniupptrait declupptrait divtrait \
     expandtrait krontrait maptrait multtrait schurtrait solvetrait subtrait

essential: all

single: all

IncludeTest: IncludeTest.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)


# Internal rules
addtrait:
	@echo
	@echo "Building the AddTrait tests..."
	@$(MAKE) --no-print-directory -C ./addtrait $(MAKECMDGOALS)

crosstrait:
	@echo
	@echo "Building the CrossTrait tests..."
	@$(MAKE) --no-print-directory -C ./crosstrait $(MAKECMDGOALS)

decldiagtrait:
	@echo
	@echo "Building the DeclDiagTrait tests..."
	@$(MAKE) --no-print-directory -C ./decldiagtrait $(MAKECMDGOALS)

declhermtrait:
	@echo
	@echo "Building the DeclHermTrait tests..."
	@$(MAKE) --no-print-directory -C ./declhermtrait $(MAKECMDGOALS)

decllowtrait:
	@echo
	@echo "Building the DeclLowTrait tests..."
	@$(MAKE) --no-print-directory -C ./decllowtrait $(MAKECMDGOALS)

declstrlowtrait:
	@echo
	@echo "Building the DeclStrLowTrait tests..."
	@$(MAKE) --no-print-directory -C ./declstrlowtrait $(MAKECMDGOALS)

declstrupptrait:
	@echo
	@echo "Building the DeclStrUppTrait tests..."
	@$(MAKE) --no-print-directory -C ./declstrupptrait $(MAKECMDGOALS)

declsymtrait:
	@echo
	@echo "Building the DeclSymTrait tests..."
	@$(MAKE) --no-print-directory -C ./declsymtrait $(MAKECMDGOALS)

declunilowtrait:
	@echo
	@echo "Building the DeclUniLowTrait tests..."
	@$(MAKE) --no-print-directory -C ./declunilowtrait $(MAKECMDGOALS)

decluniupptrait:
	@echo
	@echo "Building the DeclUniUppTrait tests..."
	@$(MAKE) --no-print-directory -C ./decluniupptrait $(MAKECMDGOALS)

declupptrait:
	@echo
	@echo "Building the DeclUppTrait tests..."
	@$(MAKE) --no-print-directory -C ./declupptrait $(MAKECMDGOALS)

divtrait:
	@echo
	@echo "Building the DivTrait tests..."
	@$(MAKE) --no-print-directory -C ./divtrait $(MAKECMDGOALS)

expandtrait:
	@echo
	@echo "Building the ExpandTrait tests..."
	@$(MAKE) --no-print-directory -C ./expandtrait $(MAKECMDGOALS)

krontrait:
	@echo
	@echo "Building the KronTrait tests..."
	@$(MAKE) --no-print-directory -C ./krontrait $(MAKECMDGOALS)

maptrait:
	@echo
	@echo "Building the MapTrait tests..."
	@$(MAKE) --no-print-directory -C ./maptrait $(MAKECMDGOALS)

multtrait:
	@echo
	@echo "Building the MultTrait tests..."
	@$(MAKE) --no-print-directory -C ./multtrait $(MAKECMDGOALS)

schurtrait:
	@echo
	@echo "Building the SchurTrait tests..."
	@$(MAKE) --no-print-directory -C ./schurtrait $(MAKECMDGOALS)

solvetrait:
	@echo
	@echo "Building the SolveTrait tests..."
	@$(MAKE) --no-print-directory -C ./solvetrait $(MAKECMDGOALS)

subtrait:
	@echo
	@echo "Building the SubTrait tests..."
	@$(MAKE) --no-print-directory -C ./subtrait $(MAKECMDGOALS)


# Cleanup
reset:
	@$(RM) $(OBJ) $(BIN)
	@$(MAKE) --no-print-directory -C ./addtrait reset
	@$(MAKE) --no-print-directory -C ./crosstrait reset
	@$(MAKE) --no-print-directory -C ./decldiagtrait reset
	@$(MAKE) --no-print-directory -C ./declhermtrait reset
	@$(MAKE) --no-print-directory -C ./decllowtrait reset
	@$(MAKE) --no-print-directory -C ./declstrlowtrait reset
	@$(MAKE) --no-print-directory -C ./declstrupptrait reset
	@$(MAKE) --no-print-directory -C ./declsymtrait reset
	@$(MAKE) --no-print-directory -C ./declunilowtrait reset
	@$(MAKE) --no-print-directory -C ./decluniupptrait reset
	@$(MAKE) --no-print-directory -C ./declupptrait reset
	@$(MAKE) --no-print-directory -C ./divtrait reset
	@$(MAKE) --no-print-directory -C ./expandtrait reset
	@$(MAKE) --no-print-directory -C ./krontrait reset
	@$(MAKE) --no-print-directory -C ./maptrait reset
	@$(MAKE) --no-print-directory -C ./multtrait reset
	@$(MAKE) --no-print-directory -C ./schurtrait reset
	@$(MAKE) --no-print-directory -C ./solvetrait reset
	@$(MAKE) --no-print-directory -C ./subtrait reset

clean:
	@$(RM) $(OBJ) $(BIN) $(DEP)
	@$(MAKE) --no-print-directory -C ./addtrait clean
	@$(MAKE) --no-print-directory -C ./crosstrait clean
	@$(MAKE) --no-print-directory -C ./decldiagtrait clean
	@$(MAKE) --no-print-directory -C ./declhermtrait clean
	@$(MAKE) --no-print-directory -C ./decllowtrait clean
	@$(MAKE) --no-print-directory -C ./declstrlowtrait clean
	@$(MAKE) --no-print-directory -C ./declstrupptrait clean
	@$(MAKE) --no-print-directory -C ./declsymtrait clean
	@$(MAKE) --no-print-directory -C ./declunilowtrait clean
	@$(MAKE) --no-print-directory -C ./decluniupptrait clean
	@$(MAKE) --no-print-directory -C ./declupptrait clean
	@$(MAKE) --no-print-directory -C ./divtrait clean
	@$(MAKE) --no-print-directory -C ./expandtrait clean
	@$(MAKE) --no-print-directory -C ./krontrait clean
	@$(MAKE) --no-print-directory -C ./maptrait clean
	@$(MAKE) --no-print-directory -C ./multtrait clean
	@$(MAKE) --no-print-directory -C ./schurtrait clean
	@$(MAKE) --no-print-directory -C ./solvetrait clean
	@$(MAKE) --no-print-directory -C ./subtrait clean


# Setting the independent commands
.PHONY: default all essential single reset clean \
        addtrait crosstrait decldiagtrait declhermtrait decllowtrait declstrlowtrait \
        declstrupptrait declsymtrait declunilowtrait decluniupptrait declupptrait divtrait \
        expandtrait krontrait maptrait multtrait schurtrait solvetrait subtrait
