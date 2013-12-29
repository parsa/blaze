#==================================================================================================
#
#  Makefile for the mathtest module of the Blaze test suite
#
#  Copyright (C) 2013 Klaus Iglberger - All Rights Reserved
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


# General rules
default: all

all: functions intrinsics \
     staticvector hybridvector dynamicvector compressedvector \
     staticmatrix dynamicmatrix compressedmatrix \
     densesubvector sparsesubvector \
     densesubmatrix sparsesubmatrix \
     denserow densecolumn sparserow sparsecolumn \
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
     dmatdmatmult dmatsmatmult smatdmatmult smatsmatmult \
     vectorserializer matrixserializer

essential: all

functions:
	@echo
	@echo "Building the mathematical functions operation tests..."
	@$(MAKE) --no-print-directory -C ./functions $(MAKECMDGOALS)

intrinsics:
	@echo
	@echo "Building the intrinsics operation tests..."
	@$(MAKE) --no-print-directory -C ./intrinsics $(MAKECMDGOALS)

staticvector:
	@echo
	@echo "Building the StaticVector tests..."
	@$(MAKE) --no-print-directory -C ./staticvector $(MAKECMDGOALS)

hybridvector:
	@echo
	@echo "Building the HybridVector tests..."
	@$(MAKE) --no-print-directory -C ./hybridvector $(MAKECMDGOALS)

dynamicvector:
	@echo
	@echo "Building the DynamicVector tests..."
	@$(MAKE) --no-print-directory -C ./dynamicvector $(MAKECMDGOALS)

compressedvector:
	@echo
	@echo "Building the CompressedVector tests..."
	@$(MAKE) --no-print-directory -C ./compressedvector $(MAKECMDGOALS)

staticmatrix:
	@echo
	@echo "Building the StaticMatrix tests..."
	@$(MAKE) --no-print-directory -C ./staticmatrix $(MAKECMDGOALS)

dynamicmatrix:
	@echo
	@echo "Building the DynamicMatrix tests..."
	@$(MAKE) --no-print-directory -C ./dynamicmatrix $(MAKECMDGOALS)

compressedmatrix:
	@echo
	@echo "Building the CompressedMatrix tests..."
	@$(MAKE) --no-print-directory -C ./compressedmatrix $(MAKECMDGOALS)

densesubvector:
	@echo
	@echo "Building the DenseSubvector tests..."
	@$(MAKE) --no-print-directory -C ./densesubvector $(MAKECMDGOALS)

sparsesubvector:
	@echo
	@echo "Building the SparseSubvector tests..."
	@$(MAKE) --no-print-directory -C ./sparsesubvector $(MAKECMDGOALS)

densesubmatrix:
	@echo
	@echo "Building the DenseSubmatrix tests..."
	@$(MAKE) --no-print-directory -C ./densesubmatrix $(MAKECMDGOALS)

sparsesubmatrix:
	@echo
	@echo "Building the SparseSubmatrix tests..."
	@$(MAKE) --no-print-directory -C ./sparsesubmatrix $(MAKECMDGOALS)

denserow:
	@echo
	@echo "Building the DenseRow tests..."
	@$(MAKE) --no-print-directory -C ./denserow $(MAKECMDGOALS)

densecolumn:
	@echo
	@echo "Building the DenseColumn tests..."
	@$(MAKE) --no-print-directory -C ./densecolumn $(MAKECMDGOALS)

sparserow:
	@echo
	@echo "Building the SparseRow tests..."
	@$(MAKE) --no-print-directory -C ./sparserow $(MAKECMDGOALS)

sparsecolumn:
	@echo
	@echo "Building the SparseColumn tests..."
	@$(MAKE) --no-print-directory -C ./sparsecolumn $(MAKECMDGOALS)

dvecdvecadd:
	@echo
	@echo "Building the dense vector/dense vector addition tests..."
	@$(MAKE) --no-print-directory -C ./dvecdvecadd $(MAKECMDGOALS)

dvecsvecadd:
	@echo
	@echo "Building the dense vector/sparse vector addition tests..."
	@$(MAKE) --no-print-directory -C ./dvecsvecadd $(MAKECMDGOALS)

svecdvecadd:
	@echo
	@echo "Building the sparse vector/dense vector addition tests..."
	@$(MAKE) --no-print-directory -C ./svecdvecadd $(MAKECMDGOALS)

svecsvecadd:
	@echo
	@echo "Building the sparse vector/sparse vector addition tests..."
	@$(MAKE) --no-print-directory -C ./svecsvecadd $(MAKECMDGOALS)

dvecdvecsub:
	@echo
	@echo "Building the dense vector/dense vector subtraction tests..."
	@$(MAKE) --no-print-directory -C ./dvecdvecsub $(MAKECMDGOALS)

dvecsvecsub:
	@echo
	@echo "Building the dense vector/sparse vector subtraction tests..."
	@$(MAKE) --no-print-directory -C ./dvecsvecsub $(MAKECMDGOALS)

svecdvecsub:
	@echo
	@echo "Building the sparse vector/dense vector subtraction tests..."
	@$(MAKE) --no-print-directory -C ./svecdvecsub $(MAKECMDGOALS)

svecsvecsub:
	@echo
	@echo "Building the sparse vector/sparse vector subtraction tests..."
	@$(MAKE) --no-print-directory -C ./svecsvecsub $(MAKECMDGOALS)

dvecdvecmult:
	@echo
	@echo "Building the dense vector/dense vector multiplication tests..."
	@$(MAKE) --no-print-directory -C ./dvecdvecmult $(MAKECMDGOALS)

dvecsvecmult:
	@echo
	@echo "Building the dense vector/sparse vector multiplication tests..."
	@$(MAKE) --no-print-directory -C ./dvecsvecmult $(MAKECMDGOALS)

svecdvecmult:
	@echo
	@echo "Building the sparse vector/dense vector multiplication tests..."
	@$(MAKE) --no-print-directory -C ./svecdvecmult $(MAKECMDGOALS)

svecsvecmult:
	@echo
	@echo "Building the sparse vector/sparse vector multiplication tests..."
	@$(MAKE) --no-print-directory -C ./svecsvecmult $(MAKECMDGOALS)

dvecdveccross:
	@echo
	@echo "Building the dense vector/dense vector cross product tests..."
	@$(MAKE) --no-print-directory -C ./dvecdveccross $(MAKECMDGOALS)

dvecsveccross:
	@echo
	@echo "Building the dense vector/sparse vector cross product tests..."
	@$(MAKE) --no-print-directory -C ./dvecsveccross $(MAKECMDGOALS)

svecdveccross:
	@echo
	@echo "Building the sparse vector/dense vector cross product tests..."
	@$(MAKE) --no-print-directory -C ./svecdveccross $(MAKECMDGOALS)

svecsveccross:
	@echo
	@echo "Building the sparse vector/sparse vector cross product tests..."
	@$(MAKE) --no-print-directory -C ./svecsveccross $(MAKECMDGOALS)

tdvecdvecmult:
	@echo
	@echo "Building the dense vector/dense vector inner product tests..."
	@$(MAKE) --no-print-directory -C ./tdvecdvecmult $(MAKECMDGOALS)

tdvecsvecmult:
	@echo
	@echo "Building the dense vector/sparse vector inner product tests..."
	@$(MAKE) --no-print-directory -C ./tdvecsvecmult $(MAKECMDGOALS)

tsvecdvecmult:
	@echo
	@echo "Building the sparse vector/dense vector inner product tests..."
	@$(MAKE) --no-print-directory -C ./tsvecdvecmult $(MAKECMDGOALS)

tsvecsvecmult:
	@echo
	@echo "Building the sparse vector/sparse vector inner product tests..."
	@$(MAKE) --no-print-directory -C ./tsvecsvecmult $(MAKECMDGOALS)

dvectdvecmult:
	@echo
	@echo "Building the dense vector/dense vector outer product tests..."
	@$(MAKE) --no-print-directory -C ./dvectdvecmult $(MAKECMDGOALS)

dvectsvecmult:
	@echo
	@echo "Building the dense vector/sparse vector outer product tests..."
	@$(MAKE) --no-print-directory -C ./dvectsvecmult $(MAKECMDGOALS)

svectdvecmult:
	@echo
	@echo "Building the sparse vector/dense vector outer product tests..."
	@$(MAKE) --no-print-directory -C ./svectdvecmult $(MAKECMDGOALS)

svectsvecmult:
	@echo
	@echo "Building the sparse vector/sparse vector outer product tests..."
	@$(MAKE) --no-print-directory -C ./svectsvecmult $(MAKECMDGOALS)

dmatdvecmult:
	@echo
	@echo "Building the dense matrix/dense vector multiplication tests..."
	@$(MAKE) --no-print-directory -C ./dmatdvecmult $(MAKECMDGOALS)

dmatsvecmult:
	@echo
	@echo "Building the dense matrix/sparse vector multiplication tests..."
	@$(MAKE) --no-print-directory -C ./dmatsvecmult $(MAKECMDGOALS)

smatdvecmult:
	@echo
	@echo "Building the sparse matrix/dense vector multiplication tests..."
	@$(MAKE) --no-print-directory -C ./smatdvecmult $(MAKECMDGOALS)

smatsvecmult:
	@echo
	@echo "Building the sparse matrix/sparse vector multiplication tests..."
	@$(MAKE) --no-print-directory -C ./smatsvecmult $(MAKECMDGOALS)

tdvecdmatmult:
	@echo
	@echo "Building the dense vector/dense matrix multiplication tests..."
	@$(MAKE) --no-print-directory -C ./tdvecdmatmult $(MAKECMDGOALS)

tsvecdmatmult:
	@echo
	@echo "Building the sparse vector/dense matrix multiplication tests..."
	@$(MAKE) --no-print-directory -C ./tsvecdmatmult $(MAKECMDGOALS)

tdvecsmatmult:
	@echo
	@echo "Building the dense vector/sparse matrix multiplication tests..."
	@$(MAKE) --no-print-directory -C ./tdvecsmatmult $(MAKECMDGOALS)

tsvecsmatmult:
	@echo
	@echo "Building the sparse vector/sparse matrix multiplication tests..."
	@$(MAKE) --no-print-directory -C ./tsvecsmatmult $(MAKECMDGOALS)

dmatdmatadd:
	@echo
	@echo "Building the dense matrix/dense matrix addition tests..."
	@$(MAKE) --no-print-directory -C ./dmatdmatadd $(MAKECMDGOALS)

dmatsmatadd:
	@echo
	@echo "Building the dense matrix/sparse matrix addition tests..."
	@$(MAKE) --no-print-directory -C ./dmatsmatadd $(MAKECMDGOALS)

smatdmatadd:
	@echo
	@echo "Building the sparse matrix/dense matrix addition tests..."
	@$(MAKE) --no-print-directory -C ./smatdmatadd $(MAKECMDGOALS)

smatsmatadd:
	@echo
	@echo "Building the sparse matrix/sparse matrix addition tests..."
	@$(MAKE) --no-print-directory -C ./smatsmatadd $(MAKECMDGOALS)

dmatdmatsub:
	@echo
	@echo "Building the dense matrix/dense matrix subtraction tests..."
	@$(MAKE) --no-print-directory -C ./dmatdmatsub $(MAKECMDGOALS)

dmatsmatsub:
	@echo
	@echo "Building the dense matrix/sparse matrix subtraction tests..."
	@$(MAKE) --no-print-directory -C ./dmatsmatsub $(MAKECMDGOALS)

smatdmatsub:
	@echo
	@echo "Building the sparse matrix/dense matrix subtraction tests..."
	@$(MAKE) --no-print-directory -C ./smatdmatsub $(MAKECMDGOALS)

smatsmatsub:
	@echo
	@echo "Building the sparse matrix/sparse matrix subtraction tests..."
	@$(MAKE) --no-print-directory -C ./smatsmatsub $(MAKECMDGOALS)

dmatdmatmult:
	@echo
	@echo "Building the dense matrix/dense matrix multiplication tests..."
	@$(MAKE) --no-print-directory -C ./dmatdmatmult $(MAKECMDGOALS)

dmatsmatmult:
	@echo
	@echo "Building the dense matrix/sparse matrix multiplication tests..."
	@$(MAKE) --no-print-directory -C ./dmatsmatmult $(MAKECMDGOALS)

smatdmatmult:
	@echo
	@echo "Building the sparse matrix/dense matrix multiplication tests..."
	@$(MAKE) --no-print-directory -C ./smatdmatmult $(MAKECMDGOALS)

smatsmatmult:
	@echo
	@echo "Building the sparse matrix/sparse matrix multiplication tests..."
	@$(MAKE) --no-print-directory -C ./smatsmatmult $(MAKECMDGOALS)

vectorserializer:
	@echo
	@echo "Building the VectorSerializer class tests..."
	@$(MAKE) --no-print-directory -C ./vectorserializer $(MAKECMDGOALS)

matrixserializer:
	@echo
	@echo "Building the MatrixSerializer class tests..."
	@$(MAKE) --no-print-directory -C ./matrixserializer $(MAKECMDGOALS)


# Cleanup
clean:
	@$(MAKE) --no-print-directory -C ./functions clean
	@$(MAKE) --no-print-directory -C ./intrinsics clean
	@$(MAKE) --no-print-directory -C ./staticvector clean
	@$(MAKE) --no-print-directory -C ./hybridvector clean
	@$(MAKE) --no-print-directory -C ./dynamicvector clean
	@$(MAKE) --no-print-directory -C ./compressedvector clean
	@$(MAKE) --no-print-directory -C ./staticmatrix clean
	@$(MAKE) --no-print-directory -C ./dynamicmatrix clean
	@$(MAKE) --no-print-directory -C ./compressedmatrix clean
	@$(MAKE) --no-print-directory -C ./densesubvector clean
	@$(MAKE) --no-print-directory -C ./sparsesubvector clean
	@$(MAKE) --no-print-directory -C ./densesubmatrix clean
	@$(MAKE) --no-print-directory -C ./sparsesubmatrix clean
	@$(MAKE) --no-print-directory -C ./denserow clean
	@$(MAKE) --no-print-directory -C ./densecolumn clean
	@$(MAKE) --no-print-directory -C ./sparserow clean
	@$(MAKE) --no-print-directory -C ./sparsecolumn clean
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
	@$(MAKE) --no-print-directory -C ./vectorserializer clean
	@$(MAKE) --no-print-directory -C ./matrixserializer clean
	@$(RM) $(OBJ) $(DEP)


# Setting the independent commands
.PHONY: default all essential clean \
        functions intrinsics \
        staticvector hybridvector dynamicvector compressedvector \
        staticmatrix dynamicmatrix compressedmatrix \
        densesubvector sparsesubvector \
        densesubmatrix sparsesubmatrix \
        denserow densecolumn sparserow sparsecolumn \
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
        dmatdmatmult dmatsmatmult smatdmatmult smatsmatmult \
        vectorserializer matrixserializer
