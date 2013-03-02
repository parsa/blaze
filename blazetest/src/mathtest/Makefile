#==================================================================================================
#
#  Makefile for the mathtest module of the Blaze test suite
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
default: staticvector dynamicvector \
         dvecdvecadd dvecsvecadd svecdvecadd svecsvecadd \
         dvecdvecsub dvecsvecsub svecdvecsub svecsvecsub \
         dvecdvecmult dvecsvecmult svecdvecmult svecsvecmult \
         dvecdveccross dvecsveccross svecdveccross svecsveccross \
         tdvecdvecmult tdvecsvecmult tsvecdvecmult tsvecsvecmult \
         dvectdvecmult dvectsvecmult svectdvecmult svectsvecmult \
         dmatdvecmult dmatsvecmult smatdvecmult smatsvecmult \
         tdvecdmatmult tsvecdmatmult tdvecsmatmult tsvecsmatmult \
         dmatdmatadd dmatsmatadd smatdmatadd smatsmatadd \
         dmatdmatsub dmatsmatsub smatdmatsub smatsmatsub \
         dmatdmatmult dmatsmatmult smatdmatmult smatsmatmult

staticvector:
	@echo
	@echo "Building the StaticVector tests..."
	@$(MAKE) --no-print-directory -C ./staticvector/

dynamicvector:
	@echo
	@echo "Building the DynamicVector tests..."
	@$(MAKE) --no-print-directory -C ./dynamicvector/

dvecdvecadd:
	@echo
	@echo "Building the dense vector/dense vector addition tests..."
	@$(MAKE) --no-print-directory -C ./dvecdvecadd/

dvecsvecadd:
	@echo
	@echo "Building the dense vector/sparse vector addition tests..."
	@$(MAKE) --no-print-directory -C ./dvecsvecadd/

svecdvecadd:
	@echo
	@echo "Building the sparse vector/dense vector addition tests..."
	@$(MAKE) --no-print-directory -C ./svecdvecadd/

svecsvecadd:
	@echo
	@echo "Building the sparse vector/sparse vector addition tests..."
	@$(MAKE) --no-print-directory -C ./svecsvecadd/

dvecdvecsub:
	@echo
	@echo "Building the dense vector/dense vector subtraction tests..."
	@$(MAKE) --no-print-directory -C ./dvecdvecsub/

dvecsvecsub:
	@echo
	@echo "Building the dense vector/sparse vector subtraction tests..."
	@$(MAKE) --no-print-directory -C ./dvecsvecsub/

svecdvecsub:
	@echo
	@echo "Building the sparse vector/dense vector subtraction tests..."
	@$(MAKE) --no-print-directory -C ./svecdvecsub/

svecsvecsub:
	@echo
	@echo "Building the sparse vector/sparse vector subtraction tests..."
	@$(MAKE) --no-print-directory -C ./svecsvecsub/

dvecdvecmult:
	@echo
	@echo "Building the dense vector/dense vector multiplication tests..."
	@$(MAKE) --no-print-directory -C ./dvecdvecmult/

dvecsvecmult:
	@echo
	@echo "Building the dense vector/sparse vector multiplication tests..."
	@$(MAKE) --no-print-directory -C ./dvecsvecmult/

svecdvecmult:
	@echo
	@echo "Building the sparse vector/dense vector multiplication tests..."
	@$(MAKE) --no-print-directory -C ./svecdvecmult/

svecsvecmult:
	@echo
	@echo "Building the sparse vector/sparse vector multiplication tests..."
	@$(MAKE) --no-print-directory -C ./svecsvecmult/

dvecdveccross:
	@echo
	@echo "Building the dense vector/dense vector cross product tests..."
	@$(MAKE) --no-print-directory -C ./dvecdveccross/

dvecsveccross:
	@echo
	@echo "Building the dense vector/sparse vector cross product tests..."
	@$(MAKE) --no-print-directory -C ./dvecsveccross/

svecdveccross:
	@echo
	@echo "Building the sparse vector/dense vector cross product tests..."
	@$(MAKE) --no-print-directory -C ./svecdveccross/

svecsveccross:
	@echo
	@echo "Building the sparse vector/sparse vector cross product tests..."
	@$(MAKE) --no-print-directory -C ./svecsveccross/

tdvecdvecmult:
	@echo
	@echo "Building the dense vector/dense vector inner product tests..."
	@$(MAKE) --no-print-directory -C ./tdvecdvecmult/

tdvecsvecmult:
	@echo
	@echo "Building the dense vector/sparse vector inner product tests..."
	@$(MAKE) --no-print-directory -C ./tdvecsvecmult/

tsvecdvecmult:
	@echo
	@echo "Building the sparse vector/dense vector inner product tests..."
	@$(MAKE) --no-print-directory -C ./tsvecdvecmult/

tsvecsvecmult:
	@echo
	@echo "Building the sparse vector/sparse vector inner product tests..."
	@$(MAKE) --no-print-directory -C ./tsvecsvecmult/

dvectdvecmult:
	@echo
	@echo "Building the dense vector/dense vector outer product tests..."
	@$(MAKE) --no-print-directory -C ./dvectdvecmult/

dvectsvecmult:
	@echo
	@echo "Building the dense vector/sparse vector outer product tests..."
	@$(MAKE) --no-print-directory -C ./dvectsvecmult/

svectdvecmult:
	@echo
	@echo "Building the sparse vector/dense vector outer product tests..."
	@$(MAKE) --no-print-directory -C ./svectdvecmult/

svectsvecmult:
	@echo
	@echo "Building the sparse vector/sparse vector outer product tests..."
	@$(MAKE) --no-print-directory -C ./svectsvecmult/

dmatdvecmult:
	@echo
	@echo "Building the dense matrix/dense vector multiplication tests..."
	@$(MAKE) --no-print-directory -C ./dmatdvecmult/

dmatsvecmult:
	@echo
	@echo "Building the dense matrix/sparse vector multiplication tests..."
	@$(MAKE) --no-print-directory -C ./dmatsvecmult/

smatdvecmult:
	@echo
	@echo "Building the sparse matrix/dense vector multiplication tests..."
	@$(MAKE) --no-print-directory -C ./smatdvecmult/

smatsvecmult:
	@echo
	@echo "Building the sparse matrix/sparse vector multiplication tests..."
	@$(MAKE) --no-print-directory -C ./smatsvecmult/

tdvecdmatmult:
	@echo
	@echo "Building the dense vector/dense matrix multiplication tests..."
	@$(MAKE) --no-print-directory -C ./tdvecdmatmult/

tsvecdmatmult:
	@echo
	@echo "Building the sparse vector/dense matrix multiplication tests..."
	@$(MAKE) --no-print-directory -C ./tsvecdmatmult/

tdvecsmatmult:
	@echo
	@echo "Building the dense vector/sparse matrix multiplication tests..."
	@$(MAKE) --no-print-directory -C ./tdvecsmatmult/

tsvecsmatmult:
	@echo
	@echo "Building the sparse vector/sparse matrix multiplication tests..."
	@$(MAKE) --no-print-directory -C ./tsvecsmatmult/

dmatdmatadd:
	@echo
	@echo "Building the dense matrix/dense matrix addition tests..."
	@$(MAKE) --no-print-directory -C ./dmatdmatadd/

dmatsmatadd:
	@echo
	@echo "Building the dense matrix/sparse matrix addition tests..."
	@$(MAKE) --no-print-directory -C ./dmatsmatadd/

smatdmatadd:
	@echo
	@echo "Building the sparse matrix/dense matrix addition tests..."
	@$(MAKE) --no-print-directory -C ./smatdmatadd/

smatsmatadd:
	@echo
	@echo "Building the sparse matrix/sparse matrix addition tests..."
	@$(MAKE) --no-print-directory -C ./smatsmatadd/

dmatdmatsub:
	@echo
	@echo "Building the dense matrix/dense matrix subtraction tests..."
	@$(MAKE) --no-print-directory -C ./dmatdmatsub/

dmatsmatsub:
	@echo
	@echo "Building the dense matrix/sparse matrix subtraction tests..."
	@$(MAKE) --no-print-directory -C ./dmatsmatsub/

smatdmatsub:
	@echo
	@echo "Building the sparse matrix/dense matrix subtraction tests..."
	@$(MAKE) --no-print-directory -C ./smatdmatsub/

smatsmatsub:
	@echo
	@echo "Building the sparse matrix/sparse matrix subtraction tests..."
	@$(MAKE) --no-print-directory -C ./smatsmatsub/

dmatdmatmult:
	@echo
	@echo "Building the dense matrix/dense matrix multiplication tests..."
	@$(MAKE) --no-print-directory -C ./dmatdmatmult/

dmatsmatmult:
	@echo
	@echo "Building the dense matrix/sparse matrix multiplication tests..."
	@$(MAKE) --no-print-directory -C ./dmatsmatmult/

smatdmatmult:
	@echo
	@echo "Building the sparse matrix/dense matrix multiplication tests..."
	@$(MAKE) --no-print-directory -C ./smatdmatmult/

smatsmatmult:
	@echo
	@echo "Building the sparse matrix/sparse matrix multiplication tests..."
	@$(MAKE) --no-print-directory -C ./smatsmatmult/

clean:
	@$(MAKE) --no-print-directory -C ./staticvector clean
	@$(MAKE) --no-print-directory -C ./dvecdvecadd clean
	@$(MAKE) --no-print-directory -C ./dvecsvecadd clean
	@$(MAKE) --no-print-directory -C ./svecdvecadd clean
	@$(MAKE) --no-print-directory -C ./svecsvecadd clean
	@$(MAKE) --no-print-directory -C ./dvecdvecsub clean
	@$(MAKE) --no-print-directory -C ./dvecsvecsub clean
	@$(MAKE) --no-print-directory -C ./svecdvecsub clean
	@$(MAKE) --no-print-directory -C ./svecsvecsub clean
	@$(MAKE) --no-print-directory -C ./dvecdvecmult clean
	@$(MAKE) --no-print-directory -C ./dvecsvecmult clean
	@$(MAKE) --no-print-directory -C ./svecdvecmult clean
	@$(MAKE) --no-print-directory -C ./svecsvecmult clean
	@$(MAKE) --no-print-directory -C ./dvecdveccross clean
	@$(MAKE) --no-print-directory -C ./dvecsveccross clean
	@$(MAKE) --no-print-directory -C ./svecdveccross clean
	@$(MAKE) --no-print-directory -C ./svecsveccross clean
	@$(MAKE) --no-print-directory -C ./tdvecdvecmult clean
	@$(MAKE) --no-print-directory -C ./tdvecsvecmult clean
	@$(MAKE) --no-print-directory -C ./tsvecdvecmult clean
	@$(MAKE) --no-print-directory -C ./tsvecsvecmult clean
	@$(MAKE) --no-print-directory -C ./dvectdvecmult clean
	@$(MAKE) --no-print-directory -C ./dvectsvecmult clean
	@$(MAKE) --no-print-directory -C ./svectdvecmult clean
	@$(MAKE) --no-print-directory -C ./svectsvecmult clean
	@$(MAKE) --no-print-directory -C ./dmatdvecmult clean
	@$(MAKE) --no-print-directory -C ./dmatsvecmult clean
	@$(MAKE) --no-print-directory -C ./smatdvecmult clean
	@$(MAKE) --no-print-directory -C ./smatsvecmult clean
	@$(MAKE) --no-print-directory -C ./tdvecdmatmult clean
	@$(MAKE) --no-print-directory -C ./tsvecdmatmult clean
	@$(MAKE) --no-print-directory -C ./tdvecsmatmult clean
	@$(MAKE) --no-print-directory -C ./tsvecsmatmult clean
	@$(MAKE) --no-print-directory -C ./dmatdmatadd clean
	@$(MAKE) --no-print-directory -C ./dmatsmatadd clean
	@$(MAKE) --no-print-directory -C ./smatdmatadd clean
	@$(MAKE) --no-print-directory -C ./smatsmatadd clean
	@$(MAKE) --no-print-directory -C ./dmatdmatsub clean
	@$(MAKE) --no-print-directory -C ./dmatsmatsub clean
	@$(MAKE) --no-print-directory -C ./smatdmatsub clean
	@$(MAKE) --no-print-directory -C ./smatsmatsub clean
	@$(MAKE) --no-print-directory -C ./dmatdmatmult clean
	@$(MAKE) --no-print-directory -C ./dmatsmatmult clean
	@$(MAKE) --no-print-directory -C ./smatdmatmult clean
	@$(MAKE) --no-print-directory -C ./smatsmatmult clean
	@$(RM) $(OBJ) $(DEP)


# Setting the independent commands
.PHONY: default clean \
        staticvector \
        dvecdvecadd dvecsvecadd svecdvecadd svecsvecadd \
        dvecdvecsub dvecsvecsub svecdvecsub svecsvecsub \
        dvecdvecmult dvecsvecmult svecdvecmult svecsvecmult \
        dvecdveccross dvecsveccross svecdveccross svecsveccross \
        tdvecdvecmult tdvecsvecmult tsvecdvecmult tsvecsvecmult \
        dvectdvecmult dvectsvecmult svectdvecmult svectsvecmult \
        dmatdvecmult dmatsvecmult smatdvecmult smatsvecmult \
        tdvecdmatmult tsvecdmatmult tdvecsmatmult tsvecsmatmult \
        dmatdmatadd dmatsmatadd smatdmatadd smatsmatadd \
        dmatdmatsub dmatsmatsub smatdmatsub smatsmatsub \
        dmatdmatmult dmatsmatmult smatdmatmult smatsmatmult
