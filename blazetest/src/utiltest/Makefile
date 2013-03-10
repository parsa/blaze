#==================================================================================================
#
#  Makefile for the utiltest module of the Blaze test suite
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


# Rules
default: uniqueptr uniquearray

uniqueptr:
	@echo
	@echo "Building the unique pointer tests..."
	@$(MAKE) --no-print-directory -C ./uniqueptr/

uniquearray:
	@echo
	@echo "Building the unique array tests..."
	@$(MAKE) --no-print-directory -C ./uniquearray/

clean:
	@$(MAKE) --no-print-directory -C ./uniqueptr clean
	@$(MAKE) --no-print-directory -C ./uniquearray clean
	@$(RM) $(OBJ) $(DEP)


# Setting the independent commands
.PHONY: default clean uniqueptr uniquearray
