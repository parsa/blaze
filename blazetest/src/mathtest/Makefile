#==================================================================================================
#
#  Makefile for the mathtest module of the Blaze test suite
#
#  Copyright (C) 2012-2020 Klaus Iglberger - All Rights Reserved
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

all: shims simd blas lapack typetraits traits constraints functors \
     densevector sparsevector densematrix sparsematrix \
     staticvector hybridvector dynamicvector customvector initializervector compressedvector \
     staticmatrix hybridmatrix dynamicmatrix custommatrix initializermatrix compressedmatrix \
     symmetricmatrix hermitianmatrix \
     lowermatrix unilowermatrix strictlylowermatrix \
     uppermatrix uniuppermatrix strictlyuppermatrix \
     diagonalmatrix identitymatrix \
     uniformvector uniformmatrix \
     zerovector zeromatrix \
     subvector elements submatrix row rows column columns band \
     dvecdvecadd dvecsvecadd svecdvecadd svecsvecadd \
     dvecdvecsub dvecsvecsub svecdvecsub svecsvecsub \
     dvecdvecmult dvecsvecmult svecdvecmult svecsvecmult \
     dvecdveccross dvecsveccross svecdveccross svecsveccross \
     dvecdvecinner dvecsvecinner svecdvecinner svecsvecinner \
     dvecdvecouter dvecsvecouter svecdvecouter svecsvecouter \
     dvecdveckron dvecsveckron svecdveckron svecsveckron \
     dvecdvecdiv svecdvecdiv \
     dvecdvecmin dvecdvecmax \
     dvecreduce svecreduce \
     dvecexpand svecexpand \
     dmatdvecmult dmatsvecmult smatdvecmult smatsvecmult \
     tdvecdmatmult tsvecdmatmult tdvecsmatmult tsvecsmatmult \
     dmatdmatadd dmatsmatadd smatdmatadd smatsmatadd \
     dmatdmatsub dmatsmatsub smatdmatsub smatsmatsub \
     dmatdmatschur dmatsmatschur smatdmatschur smatsmatschur \
     dmatdmatmult dmatsmatmult smatdmatmult smatsmatmult \
     dmatdmatkron dmatsmatkron smatdmatkron smatsmatkron \
     dmatdmatmin dmatdmatmax \
     dmatreduce smatreduce \
     determinant lu llh qr rq ql lq inversion lse exponential eigen svd \
     vectorserializer matrixserializer

essential: all

single: all

noop: shims simd blas lapack typetraits traits constraints functors \
      densevector sparsevector densematrix sparsematrix \
      staticvector hybridvector dynamicvector customvector initializervector compressedvector \
      staticmatrix hybridmatrix dynamicmatrix custommatrix initializermatrix compressedmatrix \
      symmetricmatrix hermitianmatrix \
      lowermatrix unilowermatrix strictlylowermatrix \
      uppermatrix uniuppermatrix strictlyuppermatrix \
      diagonalmatrix identitymatrix \
      uniformvector uniformmatrix \
      zerovector zeromatrix \
      subvector elements submatrix row rows column columns band \
      determinant lu llh qr rq ql lq inversion lse exponential eigen svd \
      vectorserializer matrixserializer


# Internal rules
shims:
	@echo
	@echo "Building the mathematical shims operation tests..."
	@$(MAKE) --no-print-directory -C ./shims $(MAKECMDGOALS)

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

traits:
	@$(MAKE) --no-print-directory -C ./traits $(MAKECMDGOALS)

constraints:
	@echo
	@echo "Building the constraints operation tests..."
	@$(MAKE) --no-print-directory -C ./constraints $(MAKECMDGOALS)

functors:
	@echo
	@echo "Building the functors operation tests..."
	@$(MAKE) --no-print-directory -C ./functors $(MAKECMDGOALS)

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

initializervector:
	@echo
	@echo "Building the InitializerVector tests..."
	@$(MAKE) --no-print-directory -C ./initializervector $(MAKECMDGOALS)

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

initializermatrix:
	@echo
	@echo "Building the InitializerMatrix tests..."
	@$(MAKE) --no-print-directory -C ./initializermatrix $(MAKECMDGOALS)

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

identitymatrix:
	@echo
	@echo "Building the IdentityMatrix tests..."
	@$(MAKE) --no-print-directory -C ./identitymatrix $(MAKECMDGOALS)

uniformvector:
	@echo
	@echo "Building the UniformVector tests..."
	@$(MAKE) --no-print-directory -C ./uniformvector $(MAKECMDGOALS)

uniformmatrix:
	@echo
	@echo "Building the UniformMatrix tests..."
	@$(MAKE) --no-print-directory -C ./uniformmatrix $(MAKECMDGOALS)

zerovector:
	@echo
	@echo "Building the ZeroVector tests..."
	@$(MAKE) --no-print-directory -C ./zerovector $(MAKECMDGOALS)

zeromatrix:
	@echo
	@echo "Building the ZeroMatrix tests..."
	@$(MAKE) --no-print-directory -C ./zeromatrix $(MAKECMDGOALS)

subvector:
	@echo
	@echo "Building the Subvector tests..."
	@$(MAKE) --no-print-directory -C ./subvector $(MAKECMDGOALS)

elements:
	@echo
	@echo "Building the Elements tests..."
	@$(MAKE) --no-print-directory -C ./elements $(MAKECMDGOALS)

submatrix:
	@echo
	@echo "Building the Submatrix tests..."
	@$(MAKE) --no-print-directory -C ./submatrix $(MAKECMDGOALS)

row:
	@echo
	@echo "Building the Row tests..."
	@$(MAKE) --no-print-directory -C ./row $(MAKECMDGOALS)

rows:
	@echo
	@echo "Building the Rows tests..."
	@$(MAKE) --no-print-directory -C ./rows $(MAKECMDGOALS)

column:
	@echo
	@echo "Building the Column tests..."
	@$(MAKE) --no-print-directory -C ./column $(MAKECMDGOALS)

columns:
	@echo
	@echo "Building the Columns tests..."
	@$(MAKE) --no-print-directory -C ./columns $(MAKECMDGOALS)

band:
	@echo
	@echo "Building the Band tests..."
	@$(MAKE) --no-print-directory -C ./band $(MAKECMDGOALS)

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

dvecdvecouter:
	@echo
	@echo "Building the dense vector/dense vector outer product tests..."
	@$(MAKE) --no-print-directory -C ./dvecdvecouter $(MAKECMDGOALS)

dvecsvecouter:
	@echo
	@echo "Building the dense vector/sparse vector outer product tests..."
	@$(MAKE) --no-print-directory -C ./dvecsvecouter $(MAKECMDGOALS)

svecdvecouter:
	@echo
	@echo "Building the sparse vector/dense vector outer product tests..."
	@$(MAKE) --no-print-directory -C ./svecdvecouter $(MAKECMDGOALS)

svecsvecouter:
	@echo
	@echo "Building the sparse vector/sparse vector outer product tests..."
	@$(MAKE) --no-print-directory -C ./svecsvecouter $(MAKECMDGOALS)

dvecdveckron:
	@echo
	@echo "Building the dense vector/dense vector Kronecker product tests..."
	@$(MAKE) --no-print-directory -C ./dvecdveckron $(MAKECMDGOALS)

dvecsveckron:
	@echo
	@echo "Building the dense vector/sparse vector Kronecker product tests..."
	@$(MAKE) --no-print-directory -C ./dvecsveckron $(MAKECMDGOALS)

svecdveckron:
	@echo
	@echo "Building the sparse vector/dense vector Kronecker product tests..."
	@$(MAKE) --no-print-directory -C ./svecdveckron $(MAKECMDGOALS)

svecsveckron:
	@echo
	@echo "Building the sparse vector/sparse vector Kronecker product tests..."
	@$(MAKE) --no-print-directory -C ./svecsveckron $(MAKECMDGOALS)

dvecdvecdiv:
	@echo
	@echo "Building the dense vector/dense vector division tests..."
	@$(MAKE) --no-print-directory -C ./dvecdvecdiv $(MAKECMDGOALS)

svecdvecdiv:
	@echo
	@echo "Building the sparse vector/dense vector division tests..."
	@$(MAKE) --no-print-directory -C ./svecdvecdiv $(MAKECMDGOALS)

dvecdvecmin:
	@echo
	@echo "Building the dense vector/dense vector minimum tests..."
	@$(MAKE) --no-print-directory -C ./dvecdvecmin $(MAKECMDGOALS)

dvecdvecmax:
	@echo
	@echo "Building the dense vector/dense vector maximum tests..."
	@$(MAKE) --no-print-directory -C ./dvecdvecmax $(MAKECMDGOALS)

dvecreduce:
	@echo
	@echo "Building the dense vector reduction tests..."
	@$(MAKE) --no-print-directory -C ./dvecreduce $(MAKECMDGOALS)

svecreduce:
	@echo
	@echo "Building the sparse vector reduction tests..."
	@$(MAKE) --no-print-directory -C ./svecreduce $(MAKECMDGOALS)

dvecexpand:
	@echo
	@echo "Building the dense vector expansion tests..."
	@$(MAKE) --no-print-directory -C ./dvecexpand $(MAKECMDGOALS)

svecexpand:
	@echo
	@echo "Building the sparse vector expansion tests..."
	@$(MAKE) --no-print-directory -C ./svecexpand $(MAKECMDGOALS)

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

dmatdmatschur:
	@echo
	@echo "Building the dense matrix/dense matrix Schur product tests..."
	@$(MAKE) --no-print-directory -C ./dmatdmatschur $(MAKECMDGOALS)

dmatsmatschur:
	@echo
	@echo "Building the dense matrix/sparse matrix Schur product tests..."
	@$(MAKE) --no-print-directory -C ./dmatsmatschur $(MAKECMDGOALS)

smatdmatschur:
	@echo
	@echo "Building the sparse matrix/dense matrix Schur product tests..."
	@$(MAKE) --no-print-directory -C ./smatdmatschur $(MAKECMDGOALS)

smatsmatschur:
	@echo
	@echo "Building the sparse matrix/sparse matrix Schur product tests..."
	@$(MAKE) --no-print-directory -C ./smatsmatschur $(MAKECMDGOALS)

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

dmatdmatkron:
	@echo
	@echo "Building the dense matrix/dense matrix Kronecker product tests..."
	@$(MAKE) --no-print-directory -C ./dmatdmatkron $(MAKECMDGOALS)

dmatsmatkron:
	@echo
	@echo "Building the dense matrix/sparse matrix Kronecker product tests..."
	@$(MAKE) --no-print-directory -C ./dmatsmatkron $(MAKECMDGOALS)

smatdmatkron:
	@echo
	@echo "Building the sparse matrix/dense matrix Kronecker product tests..."
	@$(MAKE) --no-print-directory -C ./smatdmatkron $(MAKECMDGOALS)

smatsmatkron:
	@echo
	@echo "Building the sparse matrix/sparse matrix Kronecker product tests..."
	@$(MAKE) --no-print-directory -C ./smatsmatkron $(MAKECMDGOALS)

dmatdmatmin:
	@echo
	@echo "Building the dense matrix/dense matrix minimum tests..."
	@$(MAKE) --no-print-directory -C ./dmatdmatmin $(MAKECMDGOALS)

dmatdmatmax:
	@echo
	@echo "Building the dense matrix/dense matrix maximum tests..."
	@$(MAKE) --no-print-directory -C ./dmatdmatmax $(MAKECMDGOALS)

dmatreduce:
	@$(MAKE) --no-print-directory -C ./dmatreduce $(MAKECMDGOALS)

smatreduce:
	@$(MAKE) --no-print-directory -C ./smatreduce $(MAKECMDGOALS)

determinant:
	@echo
	@echo "Building the matrix determinant tests..."
	@$(MAKE) --no-print-directory -C ./determinant $(MAKECMDGOALS)

lu:
	@echo
	@echo "Building the matrix LU decomposition tests..."
	@$(MAKE) --no-print-directory -C ./lu $(MAKECMDGOALS)

llh:
	@echo
	@echo "Building the matrix LLH decomposition tests..."
	@$(MAKE) --no-print-directory -C ./llh $(MAKECMDGOALS)

qr:
	@echo
	@echo "Building the matrix QR decomposition tests..."
	@$(MAKE) --no-print-directory -C ./qr $(MAKECMDGOALS)

rq:
	@echo
	@echo "Building the matrix RQ decomposition tests..."
	@$(MAKE) --no-print-directory -C ./rq $(MAKECMDGOALS)

ql:
	@echo
	@echo "Building the matrix QL decomposition tests..."
	@$(MAKE) --no-print-directory -C ./ql $(MAKECMDGOALS)

lq:
	@echo
	@echo "Building the matrix LQ decomposition tests..."
	@$(MAKE) --no-print-directory -C ./lq $(MAKECMDGOALS)

inversion:
	@echo
	@echo "Building the matrix inversion tests..."
	@$(MAKE) --no-print-directory -C ./inversion $(MAKECMDGOALS)

lse:
	@echo
	@echo "Building the matrix LSE tests..."
	@$(MAKE) --no-print-directory -C ./lse $(MAKECMDGOALS)

exponential:
	@echo
	@echo "Building the matrix exponential tests..."
	@$(MAKE) --no-print-directory -C ./exponential $(MAKECMDGOALS)

eigen:
	@echo
	@echo "Building the matrix eigenvalue/eigenvector tests..."
	@$(MAKE) --no-print-directory -C ./eigen $(MAKECMDGOALS)

svd:
	@echo
	@echo "Building the matrix singular value/vector tests..."
	@$(MAKE) --no-print-directory -C ./svd $(MAKECMDGOALS)

vectorserializer:
	@echo
	@echo "Building the VectorSerializer class tests..."
	@$(MAKE) --no-print-directory -C ./vectorserializer $(MAKECMDGOALS)

matrixserializer:
	@echo
	@echo "Building the MatrixSerializer class tests..."
	@$(MAKE) --no-print-directory -C ./matrixserializer $(MAKECMDGOALS)


# Cleanup
reset:
	@$(MAKE) --no-print-directory -C ./shims reset
	@$(MAKE) --no-print-directory -C ./simd reset
	@$(MAKE) --no-print-directory -C ./blas reset
	@$(MAKE) --no-print-directory -C ./lapack reset
	@$(MAKE) --no-print-directory -C ./typetraits reset
	@$(MAKE) --no-print-directory -C ./traits reset
	@$(MAKE) --no-print-directory -C ./constraints reset
	@$(MAKE) --no-print-directory -C ./functors reset
	@$(MAKE) --no-print-directory -C ./densevector reset
	@$(MAKE) --no-print-directory -C ./sparsevector reset
	@$(MAKE) --no-print-directory -C ./densematrix reset
	@$(MAKE) --no-print-directory -C ./sparsematrix reset
	@$(MAKE) --no-print-directory -C ./staticvector reset
	@$(MAKE) --no-print-directory -C ./hybridvector reset
	@$(MAKE) --no-print-directory -C ./dynamicvector reset
	@$(MAKE) --no-print-directory -C ./customvector reset
	@$(MAKE) --no-print-directory -C ./initializervector reset
	@$(MAKE) --no-print-directory -C ./compressedvector reset
	@$(MAKE) --no-print-directory -C ./staticmatrix reset
	@$(MAKE) --no-print-directory -C ./hybridmatrix reset
	@$(MAKE) --no-print-directory -C ./dynamicmatrix reset
	@$(MAKE) --no-print-directory -C ./custommatrix reset
	@$(MAKE) --no-print-directory -C ./initializermatrix reset
	@$(MAKE) --no-print-directory -C ./compressedmatrix reset
	@$(MAKE) --no-print-directory -C ./symmetricmatrix reset
	@$(MAKE) --no-print-directory -C ./hermitianmatrix reset
	@$(MAKE) --no-print-directory -C ./lowermatrix reset
	@$(MAKE) --no-print-directory -C ./unilowermatrix reset
	@$(MAKE) --no-print-directory -C ./strictlylowermatrix reset
	@$(MAKE) --no-print-directory -C ./uppermatrix reset
	@$(MAKE) --no-print-directory -C ./uniuppermatrix reset
	@$(MAKE) --no-print-directory -C ./strictlyuppermatrix reset
	@$(MAKE) --no-print-directory -C ./diagonalmatrix reset
	@$(MAKE) --no-print-directory -C ./identitymatrix reset
	@$(MAKE) --no-print-directory -C ./uniformvector reset
	@$(MAKE) --no-print-directory -C ./uniformmatrix reset
	@$(MAKE) --no-print-directory -C ./zerovector reset
	@$(MAKE) --no-print-directory -C ./zeromatrix reset
	@$(MAKE) --no-print-directory -C ./subvector reset
	@$(MAKE) --no-print-directory -C ./elements reset
	@$(MAKE) --no-print-directory -C ./submatrix reset
	@$(MAKE) --no-print-directory -C ./row reset
	@$(MAKE) --no-print-directory -C ./rows reset
	@$(MAKE) --no-print-directory -C ./column reset
	@$(MAKE) --no-print-directory -C ./columns reset
	@$(MAKE) --no-print-directory -C ./band reset
	@$(MAKE) --no-print-directory -C ./dvecdvecadd reset
	@$(MAKE) --no-print-directory -C ./dvecsvecadd reset
	@$(MAKE) --no-print-directory -C ./svecdvecadd reset
	@$(MAKE) --no-print-directory -C ./svecsvecadd reset
	@$(MAKE) --no-print-directory -C ./dvecdvecsub reset
	@$(MAKE) --no-print-directory -C ./dvecsvecsub reset
	@$(MAKE) --no-print-directory -C ./svecdvecsub reset
	@$(MAKE) --no-print-directory -C ./svecsvecsub reset
	@$(MAKE) --no-print-directory -C ./dvecdvecmult reset
	@$(MAKE) --no-print-directory -C ./dvecsvecmult reset
	@$(MAKE) --no-print-directory -C ./svecdvecmult reset
	@$(MAKE) --no-print-directory -C ./svecsvecmult reset
	@$(MAKE) --no-print-directory -C ./dvecdveccross reset
	@$(MAKE) --no-print-directory -C ./dvecsveccross reset
	@$(MAKE) --no-print-directory -C ./svecdveccross reset
	@$(MAKE) --no-print-directory -C ./svecsveccross reset
	@$(MAKE) --no-print-directory -C ./dvecdvecinner reset
	@$(MAKE) --no-print-directory -C ./dvecsvecinner reset
	@$(MAKE) --no-print-directory -C ./svecdvecinner reset
	@$(MAKE) --no-print-directory -C ./svecsvecinner reset
	@$(MAKE) --no-print-directory -C ./dvecdvecouter reset
	@$(MAKE) --no-print-directory -C ./dvecsvecouter reset
	@$(MAKE) --no-print-directory -C ./svecdvecouter reset
	@$(MAKE) --no-print-directory -C ./svecsvecouter reset
	@$(MAKE) --no-print-directory -C ./dvecdveckron reset
	@$(MAKE) --no-print-directory -C ./dvecsveckron reset
	@$(MAKE) --no-print-directory -C ./svecdveckron reset
	@$(MAKE) --no-print-directory -C ./svecsveckron reset
	@$(MAKE) --no-print-directory -C ./dvecdvecdiv reset
	@$(MAKE) --no-print-directory -C ./svecdvecdiv reset
	@$(MAKE) --no-print-directory -C ./dvecdvecmin reset
	@$(MAKE) --no-print-directory -C ./dvecdvecmax reset
	@$(MAKE) --no-print-directory -C ./dvecreduce reset
	@$(MAKE) --no-print-directory -C ./svecreduce reset
	@$(MAKE) --no-print-directory -C ./dvecexpand reset
	@$(MAKE) --no-print-directory -C ./svecexpand reset
	@$(MAKE) --no-print-directory -C ./dmatdvecmult reset
	@$(MAKE) --no-print-directory -C ./dmatsvecmult reset
	@$(MAKE) --no-print-directory -C ./smatdvecmult reset
	@$(MAKE) --no-print-directory -C ./smatsvecmult reset
	@$(MAKE) --no-print-directory -C ./tdvecdmatmult reset
	@$(MAKE) --no-print-directory -C ./tsvecdmatmult reset
	@$(MAKE) --no-print-directory -C ./tdvecsmatmult reset
	@$(MAKE) --no-print-directory -C ./tsvecsmatmult reset
	@$(MAKE) --no-print-directory -C ./dmatdmatadd reset
	@$(MAKE) --no-print-directory -C ./dmatsmatadd reset
	@$(MAKE) --no-print-directory -C ./smatdmatadd reset
	@$(MAKE) --no-print-directory -C ./smatsmatadd reset
	@$(MAKE) --no-print-directory -C ./dmatdmatsub reset
	@$(MAKE) --no-print-directory -C ./dmatsmatsub reset
	@$(MAKE) --no-print-directory -C ./smatdmatsub reset
	@$(MAKE) --no-print-directory -C ./smatsmatsub reset
	@$(MAKE) --no-print-directory -C ./dmatdmatschur reset
	@$(MAKE) --no-print-directory -C ./dmatsmatschur reset
	@$(MAKE) --no-print-directory -C ./smatdmatschur reset
	@$(MAKE) --no-print-directory -C ./smatsmatschur reset
	@$(MAKE) --no-print-directory -C ./dmatdmatmult reset
	@$(MAKE) --no-print-directory -C ./dmatsmatmult reset
	@$(MAKE) --no-print-directory -C ./smatdmatmult reset
	@$(MAKE) --no-print-directory -C ./smatsmatmult reset
	@$(MAKE) --no-print-directory -C ./dmatdmatkron reset
	@$(MAKE) --no-print-directory -C ./dmatsmatkron reset
	@$(MAKE) --no-print-directory -C ./smatdmatkron reset
	@$(MAKE) --no-print-directory -C ./smatsmatkron reset
	@$(MAKE) --no-print-directory -C ./dmatdmatmin reset
	@$(MAKE) --no-print-directory -C ./dmatdmatmax reset
	@$(MAKE) --no-print-directory -C ./dmatreduce reset
	@$(MAKE) --no-print-directory -C ./smatreduce reset
	@$(MAKE) --no-print-directory -C ./determinant reset
	@$(MAKE) --no-print-directory -C ./lu reset
	@$(MAKE) --no-print-directory -C ./llh reset
	@$(MAKE) --no-print-directory -C ./qr reset
	@$(MAKE) --no-print-directory -C ./rq reset
	@$(MAKE) --no-print-directory -C ./ql reset
	@$(MAKE) --no-print-directory -C ./lq reset
	@$(MAKE) --no-print-directory -C ./inversion reset
	@$(MAKE) --no-print-directory -C ./lse reset
	@$(MAKE) --no-print-directory -C ./exponential reset
	@$(MAKE) --no-print-directory -C ./eigen reset
	@$(MAKE) --no-print-directory -C ./svd reset
	@$(MAKE) --no-print-directory -C ./vectorserializer reset
	@$(MAKE) --no-print-directory -C ./matrixserializer reset

clean:
	@$(MAKE) --no-print-directory -C ./shims clean
	@$(MAKE) --no-print-directory -C ./simd clean
	@$(MAKE) --no-print-directory -C ./blas clean
	@$(MAKE) --no-print-directory -C ./lapack clean
	@$(MAKE) --no-print-directory -C ./typetraits clean
	@$(MAKE) --no-print-directory -C ./traits clean
	@$(MAKE) --no-print-directory -C ./constraints clean
	@$(MAKE) --no-print-directory -C ./functors clean
	@$(MAKE) --no-print-directory -C ./densevector clean
	@$(MAKE) --no-print-directory -C ./sparsevector clean
	@$(MAKE) --no-print-directory -C ./densematrix clean
	@$(MAKE) --no-print-directory -C ./sparsematrix clean
	@$(MAKE) --no-print-directory -C ./staticvector clean
	@$(MAKE) --no-print-directory -C ./hybridvector clean
	@$(MAKE) --no-print-directory -C ./dynamicvector clean
	@$(MAKE) --no-print-directory -C ./customvector clean
	@$(MAKE) --no-print-directory -C ./initializervector clean
	@$(MAKE) --no-print-directory -C ./compressedvector clean
	@$(MAKE) --no-print-directory -C ./staticmatrix clean
	@$(MAKE) --no-print-directory -C ./hybridmatrix clean
	@$(MAKE) --no-print-directory -C ./dynamicmatrix clean
	@$(MAKE) --no-print-directory -C ./custommatrix clean
	@$(MAKE) --no-print-directory -C ./initializermatrix clean
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
	@$(MAKE) --no-print-directory -C ./identitymatrix clean
	@$(MAKE) --no-print-directory -C ./uniformvector clean
	@$(MAKE) --no-print-directory -C ./uniformmatrix clean
	@$(MAKE) --no-print-directory -C ./zerovector clean
	@$(MAKE) --no-print-directory -C ./zeromatrix clean
	@$(MAKE) --no-print-directory -C ./subvector clean
	@$(MAKE) --no-print-directory -C ./elements clean
	@$(MAKE) --no-print-directory -C ./submatrix clean
	@$(MAKE) --no-print-directory -C ./row clean
	@$(MAKE) --no-print-directory -C ./rows clean
	@$(MAKE) --no-print-directory -C ./column clean
	@$(MAKE) --no-print-directory -C ./columns clean
	@$(MAKE) --no-print-directory -C ./band clean
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
	@$(MAKE) --no-print-directory -C ./dvecdvecouter clean
	@$(MAKE) --no-print-directory -C ./dvecsvecouter clean
	@$(MAKE) --no-print-directory -C ./svecdvecouter clean
	@$(MAKE) --no-print-directory -C ./svecsvecouter clean
	@$(MAKE) --no-print-directory -C ./dvecdveckron clean
	@$(MAKE) --no-print-directory -C ./dvecsveckron clean
	@$(MAKE) --no-print-directory -C ./svecdveckron clean
	@$(MAKE) --no-print-directory -C ./svecsveckron clean
	@$(MAKE) --no-print-directory -C ./dvecdvecdiv clean
	@$(MAKE) --no-print-directory -C ./svecdvecdiv clean
	@$(MAKE) --no-print-directory -C ./dvecdvecmin clean
	@$(MAKE) --no-print-directory -C ./dvecdvecmax clean
	@$(MAKE) --no-print-directory -C ./dvecreduce clean
	@$(MAKE) --no-print-directory -C ./svecreduce clean
	@$(MAKE) --no-print-directory -C ./dvecexpand clean
	@$(MAKE) --no-print-directory -C ./svecexpand clean
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
	@$(MAKE) --no-print-directory -C ./dmatdmatschur clean
	@$(MAKE) --no-print-directory -C ./dmatsmatschur clean
	@$(MAKE) --no-print-directory -C ./smatdmatschur clean
	@$(MAKE) --no-print-directory -C ./smatsmatschur clean
	@$(MAKE) --no-print-directory -C ./dmatdmatmult clean
	@$(MAKE) --no-print-directory -C ./dmatsmatmult clean
	@$(MAKE) --no-print-directory -C ./smatdmatmult clean
	@$(MAKE) --no-print-directory -C ./smatsmatmult clean
	@$(MAKE) --no-print-directory -C ./dmatdmatkron clean
	@$(MAKE) --no-print-directory -C ./dmatsmatkron clean
	@$(MAKE) --no-print-directory -C ./smatdmatkron clean
	@$(MAKE) --no-print-directory -C ./smatsmatkron clean
	@$(MAKE) --no-print-directory -C ./dmatdmatmin clean
	@$(MAKE) --no-print-directory -C ./dmatdmatmax clean
	@$(MAKE) --no-print-directory -C ./dmatreduce clean
	@$(MAKE) --no-print-directory -C ./smatreduce clean
	@$(MAKE) --no-print-directory -C ./determinant clean
	@$(MAKE) --no-print-directory -C ./lu clean
	@$(MAKE) --no-print-directory -C ./llh clean
	@$(MAKE) --no-print-directory -C ./qr clean
	@$(MAKE) --no-print-directory -C ./rq clean
	@$(MAKE) --no-print-directory -C ./ql clean
	@$(MAKE) --no-print-directory -C ./lq clean
	@$(MAKE) --no-print-directory -C ./inversion clean
	@$(MAKE) --no-print-directory -C ./lse clean
	@$(MAKE) --no-print-directory -C ./exponential clean
	@$(MAKE) --no-print-directory -C ./eigen clean
	@$(MAKE) --no-print-directory -C ./svd clean
	@$(MAKE) --no-print-directory -C ./vectorserializer clean
	@$(MAKE) --no-print-directory -C ./matrixserializer clean


# Setting the independent commands
.PHONY: default all essential single noop reset clean \
        shims simd blas lapack typetraits traits constraints functors \
        densevector sparsevector densematrix sparsematrix \
        staticvector hybridvector dynamicvector customvector initializervector compressedvector \
        staticmatrix hybridmatrix dynamicmatrix custommatrix initializermatrix compressedmatrix \
        symmetricmatrix hermitianmatrix \
        lowermatrix unilowermatrix strictlylowermatrix \
        uppermatrix uniuppermatrix strictlyuppermatrix \
        diagonalmatrix identitymatrix \
        uniformvector uniformmatrix \
        zerovector zeromatrix \
        subvector elements submatrix row rows column columns band \
        dvecdvecadd dvecsvecadd svecdvecadd svecsvecadd \
        dvecdvecsub dvecsvecsub svecdvecsub svecsvecsub \
        dvecdvecmult dvecsvecmult svecdvecmult svecsvecmult \
        dvecdveccross dvecsveccross svecdveccross svecsveccross \
        dvecdvecinner dvecsvecinner svecdvecinner svecsvecinner \
        dvecdvecouter dvecsvecouter svecdvecouter svecsvecouter \
        dvecdveckron dvecsveckron svecdveckron svecsveckron \
        dvecdvecdiv svecdvecdiv \
        dvecdvecmin dvecdvecmax \
        dvecreduce svecreduce \
        dvecexpand svecexpand \
        dmatdvecmult dmatsvecmult smatdvecmult smatsvecmult \
        tdvecdmatmult tsvecdmatmult tdvecsmatmult tsvecsmatmult \
        dmatdmatadd dmatsmatadd smatdmatadd smatsmatadd \
        dmatdmatsub dmatsmatsub smatdmatsub smatsmatsub \
        dmatdmatschur dmatsmatschur smatdmatschur smatsmatschur \
        dmatdmatmult dmatsmatmult smatdmatmult smatsmatmult \
        dmatdmatkron dmatsmatkron smatdmatkron smatsmatkron \
        dmatdmatmin dmatdmatmax \
        dmatreduce smatreduce \
        determinant lu llh qr rq ql lq inversion lse exponential eigen svd \
        vectorserializer matrixserializer
