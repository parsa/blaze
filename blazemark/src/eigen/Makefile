#==================================================================================================
#
#  Makefile for the Eigen module of the Blaze benchmark suite
#
#  Copyright (C) 2011 Klaus Iglberger - All Rights Reserved
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


# Module
MODULE = EIGEN


# Includes
CXXFLAGS += -I$(INSTALL_PATH)
ifneq ($(BLAZE_INCLUDE_PATH),)
CXXFLAGS += -isystem $(BLAZE_INCLUDE_PATH)
endif
ifneq ($(BOOST_INCLUDE_PATH),)
CXXFLAGS += -isystem $(BOOST_INCLUDE_PATH)
endif
ifneq ($(BLAS_INCLUDE_PATH),)
CXXFLAGS += -isystem $(BLAS_INCLUDE_PATH)
endif
ifneq ($(EIGEN_INCLUDE_PATH),)
CXXFLAGS += -isystem $(EIGEN_INCLUDE_PATH)
endif


# Setting the source, object and dependency files
SRC = $(wildcard ./*.cpp)
OBJ = $(SRC:.cpp=.o)
DEP = $(SRC:.cpp=.d)


# Rules
default: $(OBJ)
	@$(foreach dir,$(OBJ),cp -p $(dir) $(OBJECT_PATH)/$(MODULE)_$(notdir $(dir));)

clean:
	@$(RM) $(OBJ) $(DEP)


# Makefile includes
ifneq ($(MAKECMDGOALS),clean)
-include $(DEP)
endif


# Makefile generation
%.d: %.cpp
	@$(CXX) -MM -MP -MT "$*.o $*.d" -MF $@ $(CXXFLAGS) $<


# Setting the independent commands
.PHONY: default clean
