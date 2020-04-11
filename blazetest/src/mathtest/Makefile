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
     identitymatrix \
     uniformvector uniformmatrix \
     zerovector zeromatrix \
     views \
     adaptors \
     operations \
     vectorserializer matrixserializer

essential: all

single: all

noop: shims simd blas lapack typetraits traits constraints functors \
      densevector sparsevector densematrix sparsematrix \
      staticvector hybridvector dynamicvector customvector initializervector compressedvector \
      staticmatrix hybridmatrix dynamicmatrix custommatrix initializermatrix compressedmatrix \
      identitymatrix \
      uniformvector uniformmatrix \
      zerovector zeromatrix \
      views \
      adaptors \
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

views:
	@$(MAKE) --no-print-directory -C ./views $(MAKECMDGOALS)

adaptors:
	@$(MAKE) --no-print-directory -C ./adaptors $(MAKECMDGOALS)

operations:
	@$(MAKE) --no-print-directory -C ./operations $(MAKECMDGOALS)

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
	@$(MAKE) --no-print-directory -C ./identitymatrix reset
	@$(MAKE) --no-print-directory -C ./uniformvector reset
	@$(MAKE) --no-print-directory -C ./uniformmatrix reset
	@$(MAKE) --no-print-directory -C ./zerovector reset
	@$(MAKE) --no-print-directory -C ./zeromatrix reset
	@$(MAKE) --no-print-directory -C ./views reset
	@$(MAKE) --no-print-directory -C ./adaptors reset
	@$(MAKE) --no-print-directory -C ./operations reset
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
	@$(MAKE) --no-print-directory -C ./identitymatrix clean
	@$(MAKE) --no-print-directory -C ./uniformvector clean
	@$(MAKE) --no-print-directory -C ./uniformmatrix clean
	@$(MAKE) --no-print-directory -C ./zerovector clean
	@$(MAKE) --no-print-directory -C ./zeromatrix clean
	@$(MAKE) --no-print-directory -C ./views clean
	@$(MAKE) --no-print-directory -C ./adaptors clean
	@$(MAKE) --no-print-directory -C ./operations clean
	@$(MAKE) --no-print-directory -C ./vectorserializer clean
	@$(MAKE) --no-print-directory -C ./matrixserializer clean


# Setting the independent commands
.PHONY: default all essential single noop reset clean \
        shims simd blas lapack typetraits traits constraints functors \
        densevector sparsevector densematrix sparsematrix \
        staticvector hybridvector dynamicvector customvector initializervector compressedvector \
        staticmatrix hybridmatrix dynamicmatrix custommatrix initializermatrix compressedmatrix \
        identitymatrix \
        uniformvector uniformmatrix \
        zerovector zeromatrix \
        views \
        adaptors \
        operations \
        vectorserializer matrixserializer
