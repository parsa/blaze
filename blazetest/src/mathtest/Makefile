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


# Build rules
default: all

all: functions simd blas lapack typetraits \
     densevector sparsevector densematrix sparsematrix \
     staticvector hybridvector dynamicvector customvector compressedvector \
     staticmatrix hybridmatrix dynamicmatrix custommatrix compressedmatrix \
     symmetricmatrix hermitianmatrix \
     lowermatrix unilowermatrix strictlylowermatrix \
     uppermatrix uniuppermatrix strictlyuppermatrix \
     diagonalmatrix \
     subvector submatrix row column \
     dvecdvecadd dvecsvecadd svecdvecadd svecsvecadd \
     dvecdvecsub dvecsvecsub svecdvecsub svecsvecsub \
     dvecdvecmult dvecsvecmult svecdvecmult svecsvecmult \
     dvecdveccross dvecsveccross svecdveccross svecsveccross \
     dvecdvecinner dvecsvecinner svecdvecinner svecsvecinner \
     dvectdvecmult dvectsvecmult svectdvecmult svectsvecmult \
     dvecdvecdiv svecdvecdiv \
     dmatdvecmult dmatsvecmult smatdvecmult smatsvecmult \
     tdvecdmatmult tsvecdmatmult tdvecsmatmult tsvecsmatmult \
     dmatdmatadd dmatsmatadd smatdmatadd smatsmatadd \
     dmatdmatsub dmatsmatsub smatdmatsub smatsmatsub \
     dmatdmatmult dmatsmatmult smatdmatmult smatsmatmult \
     determinant decomposition inversion \
     vectorserializer matrixserializer

essential: all

single: all

noop: functions simd blas lapack typetraits \
      densevector sparsevector densematrix sparsematrix \
      staticvector hybridvector dynamicvector customvector compressedvector \
      staticmatrix hybridmatrix dynamicmatrix custommatrix compressedmatrix \
      symmetricmatrix hermitianmatrix \
      lowermatrix unilowermatrix strictlylowermatrix \
      uppermatrix uniuppermatrix strictlyuppermatrix \
      diagonalmatrix \
      subvector submatrix row column \
      determinant decomposition inversion \
      vectorserializer matrixserializer


# Internal rules
functions:
	@echo
	@echo "Building the mathematical functions operation tests..."
	@$(MAKE) --no-print-directory -C ./functions $(MAKECMDGOALS)

simd:
	@echo
	@echo "Building the SIMD operation tests..."
	@$(MAKE) --no-print-directory -C ./simd $(MAKECMDGOALS)

blas:
	@echo
	@echo "Building the BLAS operation tests..."
	@$(MAKE) --no-print-directory -C ./blas $(MAKECMDGOALS)

lapack:
	@echo
	@echo "Building the LAPACK operation tests..."
	@$(MAKE) --no-print-directory -C ./lapack $(MAKECMDGOALS)

typetraits:
	@echo
	@echo "Building the typetraits operation tests..."
	@$(MAKE) --no-print-directory -C ./typetraits $(MAKECMDGOALS)

densevector:
	@echo
	@echo "Building the DenseVector operation tests..."
	@$(MAKE) --no-print-directory -C ./densevector $(MAKECMDGOALS)

sparsevector:
	@echo
	@echo "Building the SparseVector operation tests..."
	@$(MAKE) --no-print-directory -C ./sparsevector $(MAKECMDGOALS)

densematrix:
	@echo
	@echo "Building the DenseMatrix operation tests..."
	@$(MAKE) --no-print-directory -C ./densematrix $(MAKECMDGOALS)

sparsematrix:
	@echo
	@echo "Building the SparseMatrix operation tests..."
	@$(MAKE) --no-print-directory -C ./sparsematrix $(MAKECMDGOALS)

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

customvector:
	@echo
	@echo "Building the CustomVector tests..."
	@$(MAKE) --no-print-directory -C ./customvector $(MAKECMDGOALS)

compressedvector:
	@echo
	@echo "Building the CompressedVector tests..."
	@$(MAKE) --no-print-directory -C ./compressedvector $(MAKECMDGOALS)

staticmatrix:
	@echo
	@echo "Building the StaticMatrix tests..."
	@$(MAKE) --no-print-directory -C ./staticmatrix $(MAKECMDGOALS)

hybridmatrix:
	@echo
	@echo "Building the HybridMatrix tests..."
	@$(MAKE) --no-print-directory -C ./hybridmatrix $(MAKECMDGOALS)

dynamicmatrix:
	@echo
	@echo "Building the DynamicMatrix tests..."
	@$(MAKE) --no-print-directory -C ./dynamicmatrix $(MAKECMDGOALS)

custommatrix:
	@echo
	@echo "Building the CustomMatrix tests..."
	@$(MAKE) --no-print-directory -C ./custommatrix $(MAKECMDGOALS)

compressedmatrix:
	@echo
	@echo "Building the CompressedMatrix tests..."
	@$(MAKE) --no-print-directory -C ./compressedmatrix $(MAKECMDGOALS)

symmetricmatrix:
	@echo
	@echo "Building the SymmetricMatrix tests..."
	@$(MAKE) --no-print-directory -C ./symmetricmatrix $(MAKECMDGOALS)

hermitianmatrix:
	@echo
	@echo "Building the HermitianMatrix tests..."
	@$(MAKE) --no-print-directory -C ./hermitianmatrix $(MAKECMDGOALS)

lowermatrix:
	@echo
	@echo "Building the LowerMatrix tests..."
	@$(MAKE) --no-print-directory -C ./lowermatrix $(MAKECMDGOALS)

unilowermatrix:
	@echo
	@echo "Building the UniLowerMatrix tests..."
	@$(MAKE) --no-print-directory -C ./unilowermatrix $(MAKECMDGOALS)

strictlylowermatrix:
	@echo
	@echo "Building the StrictlyLowerMatrix tests..."
	@$(MAKE) --no-print-directory -C ./strictlylowermatrix $(MAKECMDGOALS)

uppermatrix:
	@echo
	@echo "Building the UpperMatrix tests..."
	@$(MAKE) --no-print-directory -C ./uppermatrix $(MAKECMDGOALS)

uniuppermatrix:
	@echo
	@echo "Building the UniUpperMatrix tests..."
	@$(MAKE) --no-print-directory -C ./uniuppermatrix $(MAKECMDGOALS)

strictlyuppermatrix:
	@echo
	@echo "Building the StrictlyUpperMatrix tests..."
	@$(MAKE) --no-print-directory -C ./strictlyuppermatrix $(MAKECMDGOALS)

diagonalmatrix:
	@echo
	@echo "Building the DiagonalMatrix tests..."
	@$(MAKE) --no-print-directory -C ./diagonalmatrix $(MAKECMDGOALS)

subvector:
	@echo
	@echo "Building the Subvector tests..."
	@$(MAKE) --no-print-directory -C ./subvector $(MAKECMDGOALS)

submatrix:
	@echo
	@echo "Building the Submatrix tests..."
	@$(MAKE) --no-print-directory -C ./submatrix $(MAKECMDGOALS)

row:
	@echo
	@echo "Building the Row tests..."
	@$(MAKE) --no-print-directory -C ./row $(MAKECMDGOALS)

column:
	@echo
	@echo "Building the Column tests..."
	@$(MAKE) --no-print-directory -C ./column $(MAKECMDGOALS)

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

dvecdvecinner:
	@echo
	@echo "Building the dense vector/dense vector inner product tests..."
	@$(MAKE) --no-print-directory -C ./dvecdvecinner $(MAKECMDGOALS)

dvecsvecinner:
	@echo
	@echo "Building the dense vector/sparse vector inner product tests..."
	@$(MAKE) --no-print-directory -C ./dvecsvecinner $(MAKECMDGOALS)

svecdvecinner:
	@echo
	@echo "Building the sparse vector/dense vector inner product tests..."
	@$(MAKE) --no-print-directory -C ./svecdvecinner $(MAKECMDGOALS)

svecsvecinner:
	@echo
	@echo "Building the sparse vector/sparse vector inner product tests..."
	@$(MAKE) --no-print-directory -C ./svecsvecinner $(MAKECMDGOALS)

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

dvecdvecdiv:
	@echo
	@echo "Building the dense vector/dense vector division tests..."
	@$(MAKE) --no-print-directory -C ./dvecdvecdiv $(MAKECMDGOALS)

svecdvecdiv:
	@echo
	@echo "Building the sparse vector/dense vector division tests..."
	@$(MAKE) --no-print-directory -C ./svecdvecdiv $(MAKECMDGOALS)

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

determinant:
	@echo
	@echo "Building the matrix determinant tests..."
	@$(MAKE) --no-print-directory -C ./determinant $(MAKECMDGOALS)

decomposition:
	@echo
	@echo "Building the matrix decomposition tests..."
	@$(MAKE) --no-print-directory -C ./decomposition $(MAKECMDGOALS)

inversion:
	@echo
	@echo "Building the matrix inversion tests..."
	@$(MAKE) --no-print-directory -C ./inversion $(MAKECMDGOALS)

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
	@$(MAKE) --no-print-directory -C ./simd clean
	@$(MAKE) --no-print-directory -C ./blas clean
	@$(MAKE) --no-print-directory -C ./lapack clean
	@$(MAKE) --no-print-directory -C ./typetraits clean
	@$(MAKE) --no-print-directory -C ./densevector clean
	@$(MAKE) --no-print-directory -C ./sparsevector clean
	@$(MAKE) --no-print-directory -C ./densematrix clean
	@$(MAKE) --no-print-directory -C ./sparsematrix clean
	@$(MAKE) --no-print-directory -C ./staticvector clean
	@$(MAKE) --no-print-directory -C ./hybridvector clean
	@$(MAKE) --no-print-directory -C ./dynamicvector clean
	@$(MAKE) --no-print-directory -C ./customvector clean
	@$(MAKE) --no-print-directory -C ./compressedvector clean
	@$(MAKE) --no-print-directory -C ./staticmatrix clean
	@$(MAKE) --no-print-directory -C ./hybridmatrix clean
	@$(MAKE) --no-print-directory -C ./dynamicmatrix clean
	@$(MAKE) --no-print-directory -C ./custommatrix clean
	@$(MAKE) --no-print-directory -C ./compressedmatrix clean
	@$(MAKE) --no-print-directory -C ./symmetricmatrix clean
	@$(MAKE) --no-print-directory -C ./hermitianmatrix clean
	@$(MAKE) --no-print-directory -C ./lowermatrix clean
	@$(MAKE) --no-print-directory -C ./unilowermatrix clean
	@$(MAKE) --no-print-directory -C ./strictlylowermatrix clean
	@$(MAKE) --no-print-directory -C ./uppermatrix clean
	@$(MAKE) --no-print-directory -C ./uniuppermatrix clean
	@$(MAKE) --no-print-directory -C ./strictlyuppermatrix clean
	@$(MAKE) --no-print-directory -C ./diagonalmatrix clean
	@$(MAKE) --no-print-directory -C ./subvector clean
	@$(MAKE) --no-print-directory -C ./submatrix clean
	@$(MAKE) --no-print-directory -C ./row clean
	@$(MAKE) --no-print-directory -C ./column clean
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
	@$(MAKE) --no-print-directory -C ./dvecdvecinner clean
	@$(MAKE) --no-print-directory -C ./dvecsvecinner clean
	@$(MAKE) --no-print-directory -C ./svecdvecinner clean
	@$(MAKE) --no-print-directory -C ./svecsvecinner clean
	@$(MAKE) --no-print-directory -C ./dvectdvecmult clean
	@$(MAKE) --no-print-directory -C ./dvectsvecmult clean
	@$(MAKE) --no-print-directory -C ./svectdvecmult clean
	@$(MAKE) --no-print-directory -C ./svectsvecmult clean
	@$(MAKE) --no-print-directory -C ./dvecdvecdiv clean
	@$(MAKE) --no-print-directory -C ./svecdvecdiv clean
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
	@$(MAKE) --no-print-directory -C ./determinant clean
	@$(MAKE) --no-print-directory -C ./decomposition clean
	@$(MAKE) --no-print-directory -C ./inversion clean
	@$(MAKE) --no-print-directory -C ./vectorserializer clean
	@$(MAKE) --no-print-directory -C ./matrixserializer clean
	@$(RM) $(OBJ) $(DEP)


# Setting the independent commands
.PHONY: default all essential single noop clean \
        functions simd blas lapack typetraits \
        densevector sparsevector densematrix sparsematrix \
        staticvector hybridvector dynamicvector customvector compressedvector \
        staticmatrix hybridmatrix dynamicmatrix custommatrix compressedmatrix \
        symmetricmatrix hermitianmatrix \
        lowermatrix unilowermatrix strictlylowermatrix \
        uppermatrix uniuppermatrix strictlyuppermatrix \
        diagonalmatrix \
        subvector submatrix row column \
        dvecdvecadd dvecsvecadd svecdvecadd svecsvecadd \
        dvecdvecsub dvecsvecsub svecdvecsub svecsvecsub \
        dvecdvecmult dvecsvecmult svecdvecmult svecsvecmult \
        dvecdveccross dvecsveccross svecdveccross svecsveccross \
        dvecdvecinner dvecsvecinner svecdvecinner svecsvecinner \
        dvectdvecmult dvectsvecmult svectdvecmult svectsvecmult \
        dvecdvecdiv svecdvecdiv \
        dmatdvecmult dmatsvecmult smatdvecmult smatsvecmult \
        tdvecdmatmult tsvecdmatmult tdvecsmatmult tsvecsmatmult \
        dmatdmatadd dmatsmatadd smatdmatadd smatsmatadd \
        dmatdmatsub dmatsmatsub smatdmatsub smatsmatsub \
        dmatdmatmult dmatsmatmult smatdmatmult smatsmatmult \
        determinant decomposition inversion \
        vectorserializer matrixserializer
