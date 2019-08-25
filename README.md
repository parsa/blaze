![blaze300x150.jpg](https://bitbucket.org/blaze-lib/blaze/wiki/images/blaze300x150.jpg)

**Blaze** is an open-source, high-performance C++ math library for dense and sparse arithmetic. With its state-of-the-art *Smart Expression Template* implementation **Blaze** combines the elegance and ease of use of a domain-specific language with HPC-grade performance, making it one of the most intuitive and fastest C++ math libraries available.

The **Blaze** library offers ...

  * ... **high performance** through the integration of BLAS libraries and manually tuned HPC math kernels
  * ... **vectorization** by SSE, SSE2, SSE3, SSSE3, SSE4, AVX, AVX2, AVX-512, FMA, and SVML
  * ... **parallel execution** by OpenMP, HPX, C++11 threads and Boost threads
  * ... the **intuitive** and **easy to use** API of a domain specific language
  * ... **unified arithmetic** with dense and sparse vectors and matrices
  * ... **thoroughly tested** matrix and vector arithmetic
  * ... completely **portable**, **high quality** C++ source code

Get an impression of the clear but powerful syntax of **Blaze** in the [Getting Started](https://bitbucket.org/blaze-lib/blaze/wiki/Getting_Started) tutorial and of the impressive performance in the [Benchmarks](https://bitbucket.org/blaze-lib/blaze/wiki/Benchmarks) section.

----

## Download ##

![white20x120.jpg](https://bitbucket.org/blaze-lib/blaze/wiki/images/white20x120.jpg)
[![blaze-3.6.jpg](https://bitbucket.org/blaze-lib/blaze/wiki/images/blaze-3.6.jpg)](https://bitbucket.org/blaze-lib/blaze/downloads/blaze-3.6.tar.gz)
![white40x120.jpg](https://bitbucket.org/blaze-lib/blaze/wiki/images/white40x120.jpg)
[![blaze-docu-3.6.jpg](https://bitbucket.org/blaze-lib/blaze/wiki/images/blaze-docu-3.6.jpg)](https://bitbucket.org/blaze-lib/blaze/downloads/blaze-docu-3.6.tar.gz)

Older releases of **Blaze** can be found in the [downloads](https://bitbucket.org/blaze-lib/blaze/downloads) section or in our [release archive](https://bitbucket.org/blaze-lib/blaze/wiki/Release Archive).

----

## Blaze Projects ##

[Blaze CUDA](https://github.com/STEllAR-GROUP/blaze_cuda): Add CUDA capabilities to the **Blaze** library (Jules Pénuchot)

[blaze_tensor](https://github.com/STEllAR-GROUP/blaze_tensor): An implementation of 3D tensors for the **Blaze** library (Stellar Group)

[BlazeIterative](https://github.com/tjolsen/BlazeIterative): A collection of iterative solvers (CG, BiCGSTAB, ...) for the **Blaze** library (Tyler Olsen)

[RcppBlaze](https://github.com/ChingChuan-Chen/RcppBlaze): A **Blaze** port for the R language (ChingChuan Chen)

----

## News ##

**25.8.2019**: On time for [CppCon 2019](https://cppcon.org) and [SC19](https://sc19.supercomputing.org) we release the next **Blaze** milestone. **Blaze** 3.6 comes with a multitude of new features including the Kronecker product for [vectors](https://bitbucket.org/blaze-lib/blaze/wiki/Vector-Vector Multiplication#!kronecker-product) and [matrices](https://bitbucket.org/blaze-lib/blaze/wiki/Matrix-Matrix Multiplication#!kronecker-product), the [`mean()`](https://bitbucket.org/blaze-lib/blaze/wiki/Vector Operations#!mean), [`var()`](https://bitbucket.org/blaze-lib/blaze/wiki/Vector Operations#!var) and [`stddev()`](https://bitbucket.org/blaze-lib/blaze/wiki/Vector Operations#!stddev) statistic functions for all kinds of vectors and matrices, [scalar additions](https://bitbucket.org/blaze-lib/blaze/wiki/Addition#!scalar_addition), [scalar subtractions](https://bitbucket.org/blaze-lib/blaze/wiki/Addition#!scalar_subtraction), and [scalar expansion](https://bitbucket.org/blaze-lib/blaze/wiki/Vector Operations#!scalar-expansion). Furthermore we have integrated various bitwise operations for dense vectors and dense matrices ([bitwise shift](https://bitbucket.org/blaze-lib/blaze/wiki/Bitwise Shift), [AND](https://bitbucket.org/blaze-lib/blaze/wiki/Bitwise AND), [OR](https://bitbucket.org/blaze-lib/blaze/wiki/Bitwise OR), and [XOR](https://bitbucket.org/blaze-lib/blaze/wiki/Bitwise XOR)) and logical operations ([logical NOT](https://bitbucket.org/blaze-lib/blaze/wiki/Logical NOT), [AND](https://bitbucket.org/blaze-lib/blaze/wiki/Logical AND), and [OR](https://bitbucket.org/blaze-lib/blaze/wiki/Logical OR)). We hope you enjoy this amazing release of **Blaze**!

**26.2.2019**: Today we present the next evolution of the **Blaze** library, **Blaze** 3.5. This new release introduces several new, requested features:

* New vector and matrix types, specifically [```UniformVector```](https://bitbucket.org/blaze-lib/blaze/wiki/Vector%20Types#!uniformvector), [```UniformMatrix```](https://bitbucket.org/blaze-lib/blaze/wiki/Matrix%20Types#!uniformmatrix), [```ZeroVector```](https://bitbucket.org/blaze-lib/blaze/wiki/Vector%20Types#!zerovector), and [```ZeroMatrix```](https://bitbucket.org/blaze-lib/blaze/wiki/Matrix%20Types#!zeromatrix):

```
#!c++
blaze::UniformVector<int> u( 5UL );          // Creating a 5D uniform vector
blaze::UniformMatrix<double> U( 4UL, 6UL );  // Creating a 4x6 uniform matrix

blaze::ZeroVector<float> z( 4UL );           // Creating a 4D zero vector
blaze::ZeroMatrix<double> Z( 3UL, 7UL );     // Creating a 3x7 zero matrix
```

* More flexible element selections, row selections, and column selections:

```
#!c++
blaze::DynamicVector<double,blaze::rowVector> x{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };
blaze::DynamicMatrix<double,blaze::rowMajor> A( 9UL, 9UL );

// Selecting all even elements of the vector, i.e. selecting (1,3,5,7,9)
auto e = elements( x, []( size_t i ){ return i*2UL; }, 5UL );

// Selecting all odd rows of the matrix, i.e. selecting the rows 1, 3, 5, and 7
auto rs = rows( A, []( size_t i ){ return i*2UL+1UL; }, 4UL );

// Reversing the columns of the matrix, i.e. selecting the columns 8, 7, 6, 5, 4, 3, 2, 1, and 0
auto cs = columns( A, [max=A.columns()-1UL]( size_t i ){ return max-i; }, 9UL );
```

* Vector expansion via the [```expand()```](https://bitbucket.org/blaze-lib/blaze/wiki/Vector%20Operations#!vector-expansion) function:

```
#!c++
blaze::DynamicVector<int,columnVector> a{ 1, 2, 3 };
blaze::CompressedVector<int,rowVector> b{ 1, 0, 3, 0, 5 };

// Expand the dense column vector ( 1 2 3 ) into a dense 3x5 column-major matrix
//
//   ( 1 1 1 1 1 )
//   ( 2 2 2 2 2 )
//   ( 3 3 3 3 3 )
//
expand( a, 5 );  // Runtime parameter
expand<5>( a );  // Compile time parameter

// Expand the sparse row vector ( 1 0 3 0 5 ) into a sparse 3x5 row-major matrix
//
//   ( 1 0 3 0 5 )
//   ( 1 0 3 0 5 )
//   ( 1 0 3 0 5 )
//
expand( b, 3 );  // Runtime parameter
expand<3>( b );  // Compile time parameter
```

With the release of **Blaze** 3.5 we also officially deprecate the **Blazemark**, which means that we will eventually remove it in an upcoming release. We hope that you enjoy this new release!

----

## Wiki: Table of Contents ##

* [Configuration and Installation](https://bitbucket.org/blaze-lib/blaze/wiki/Configuration and Installation)
* [Getting Started](https://bitbucket.org/blaze-lib/blaze/wiki/Getting Started)
* [Vectors](https://bitbucket.org/blaze-lib/blaze/wiki/Vectors)
    * [Vector Types](https://bitbucket.org/blaze-lib/blaze/wiki/Vector Types)
    * [Vector Operations](https://bitbucket.org/blaze-lib/blaze/wiki/Vector Operations)
* [Matrices](https://bitbucket.org/blaze-lib/blaze/wiki/Matrices)
    * [Matrix Types](https://bitbucket.org/blaze-lib/blaze/wiki/Matrix Types)
    * [Matrix Operations](https://bitbucket.org/blaze-lib/blaze/wiki/Matrix Operations)
* [Adaptors](https://bitbucket.org/blaze-lib/blaze/wiki/Adaptors)
    * [Symmetric Matrices](https://bitbucket.org/blaze-lib/blaze/wiki/Symmetric Matrices)
    * [Hermitian Matrices](https://bitbucket.org/blaze-lib/blaze/wiki/Hermitian Matrices)
    * [Triangular Matrices](https://bitbucket.org/blaze-lib/blaze/wiki/Triangular Matrices)
* [Views](https://bitbucket.org/blaze-lib/blaze/wiki/Views)
    * [Subvectors](https://bitbucket.org/blaze-lib/blaze/wiki/Subvectors)
    * [Element Selections](https://bitbucket.org/blaze-lib/blaze/wiki/Element Selections)
    * [Submatrices](https://bitbucket.org/blaze-lib/blaze/wiki/Submatrices)
    * [Rows](https://bitbucket.org/blaze-lib/blaze/wiki/Rows)
    * [Row Selections](https://bitbucket.org/blaze-lib/blaze/wiki/Row Selections)
    * [Columns](https://bitbucket.org/blaze-lib/blaze/wiki/Columns)
    * [Column Selections](https://bitbucket.org/blaze-lib/blaze/wiki/Column Selections)
    * [Bands](https://bitbucket.org/blaze-lib/blaze/wiki/Bands)
* [Arithmetic Operations](https://bitbucket.org/blaze-lib/blaze/wiki/Arithmetic Operations)
    * [Addition](https://bitbucket.org/blaze-lib/blaze/wiki/Addition)
    * [Subtraction](https://bitbucket.org/blaze-lib/blaze/wiki/Subtraction)
    * [Scalar Multiplication](https://bitbucket.org/blaze-lib/blaze/wiki/Scalar Multiplication)
    * [Vector/Vector Multiplication](https://bitbucket.org/blaze-lib/blaze/wiki/Vector-Vector Multiplication)
        * [Componentwise Multiplication](https://bitbucket.org/blaze-lib/blaze/wiki/Vector-Vector Multiplication#!componentwise-multiplication)
        * [Inner Product / Scalar Product / Dot Product](https://bitbucket.org/blaze-lib/blaze/wiki/Vector-Vector Multiplication#!inner-product-scalar-product-dot-product)
        * [Outer Product](https://bitbucket.org/blaze-lib/blaze/wiki/Vector-Vector Multiplication#!outer-product)
        * [Cross Product](https://bitbucket.org/blaze-lib/blaze/wiki/Vector-Vector Multiplication#!cross-product)
        * [Kronecker Product](https://bitbucket.org/blaze-lib/blaze/wiki/Vector-Vector Multiplication#!kronecker-product)
    * [Vector/Vector Division](https://bitbucket.org/blaze-lib/blaze/wiki/Vector-Vector Division)
    * [Matrix/Vector Multiplication](https://bitbucket.org/blaze-lib/blaze/wiki/Matrix-Vector Multiplication)
    * [Matrix/Matrix Multiplication](https://bitbucket.org/blaze-lib/blaze/wiki/Matrix-Matrix Multiplication)
        * [Schur Product](https://bitbucket.org/blaze-lib/blaze/wiki/Matrix-Matrix Multiplication#!componentwise-multiplication-schur-product)
        * [Matrix Product](https://bitbucket.org/blaze-lib/blaze/wiki/Matrix-Matrix Multiplication#!matrix-product)
        * [Kronecker Product](https://bitbucket.org/blaze-lib/blaze/wiki/Matrix-Matrix Multiplication#!kronecker-product)
* [Bitwise Operations](https://bitbucket.org/blaze-lib/blaze/wiki/Arithmetic Operations)
    * [Bitwise Shift](https://bitbucket.org/blaze-lib/blaze/wiki/Bitwise Shift)
    * [Bitwise AND](https://bitbucket.org/blaze-lib/blaze/wiki/Bitwise AND)
    * [Bitwise OR](https://bitbucket.org/blaze-lib/blaze/wiki/Bitwise OR)
    * [Bitwise XOR](https://bitbucket.org/blaze-lib/blaze/wiki/Bitwise XOR)
* [Logical Operations](https://bitbucket.org/blaze-lib/blaze/wiki/Arithmetic Operations)
    * [Logical NOT](https://bitbucket.org/blaze-lib/blaze/wiki/Logical NOT)
    * [Logical AND](https://bitbucket.org/blaze-lib/blaze/wiki/Logical AND)
    * [Logical OR](https://bitbucket.org/blaze-lib/blaze/wiki/Logical OR)
* [Shared-Memory Parallelization](https://bitbucket.org/blaze-lib/blaze/wiki/Shared Memory Parallelization)
    * [HPX Parallelization](https://bitbucket.org/blaze-lib/blaze/wiki/HPX Parallelization)
    * [C++11 Thread Parallelization](https://bitbucket.org/blaze-lib/blaze/wiki/Cpp Thread Parallelization)
    * [Boost Thread Parallelization](https://bitbucket.org/blaze-lib/blaze/wiki/Boost Thread Parallelization)
    * [OpenMP Parallelization](https://bitbucket.org/blaze-lib/blaze/wiki/OpenMP Parallelization)
    * [Serial Execution](https://bitbucket.org/blaze-lib/blaze/wiki/Serial Execution)
* [Serialization](https://bitbucket.org/blaze-lib/blaze/wiki/Serialization)
    * [Vector Serialization](https://bitbucket.org/blaze-lib/blaze/wiki/Vector Serialization)
    * [Matrix Serialization](https://bitbucket.org/blaze-lib/blaze/wiki/Matrix Serialization)
* [Customization](https://bitbucket.org/blaze-lib/blaze/wiki/Customization)
    * [Configuration Files](https://bitbucket.org/blaze-lib/blaze/wiki/Configuration Files)
    * [Vector and Matrix Customization](https://bitbucket.org/blaze-lib/blaze/wiki/Vector and Matrix Customization)
        * [Custom Data Members](https://bitbucket.org/blaze-lib/blaze/wiki/Vector and Matrix Customization#!custom-data-members)
        * [Custom Operations](https://bitbucket.org/blaze-lib/blaze/wiki/Vector and Matrix Customization#!custom-operations)
        * [Custom Data Types](https://bitbucket.org/blaze-lib/blaze/wiki/Vector and Matrix Customization#!custom-data-types)
    * [Error Reporting Customization](https://bitbucket.org/blaze-lib/blaze/wiki/Error Reporting Customization)
* [BLAS Functions](https://bitbucket.org/blaze-lib/blaze/wiki/BLAS Functions)
* [LAPACK Functions](https://bitbucket.org/blaze-lib/blaze/wiki/LAPACK Functions)
* [Block Vectors and Matrices](https://bitbucket.org/blaze-lib/blaze/wiki/Block Vectors and Matrices)
* [Intra-Statement Optimization](https://bitbucket.org/blaze-lib/blaze/wiki/Intra-Statement Optimization)
* [Frequently Asked Questions (FAQ)](https://bitbucket.org/blaze-lib/blaze/wiki/FAQ)
* [Issue Creation Guidelines](https://bitbucket.org/blaze-lib/blaze/wiki/Issue Creation Guidelines)
* [Blaze References](https://bitbucket.org/blaze-lib/blaze/wiki/Blaze References)
* [Blazemark: The Blaze Benchmark Suite](https://bitbucket.org/blaze-lib/blaze/wiki/Blazemark)
* [Benchmarks/Performance Results](https://bitbucket.org/blaze-lib/blaze/wiki/Benchmarks)

----

## License ##

The **Blaze** library is licensed under the New (Revised) BSD license. Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

  * Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
  * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
  * Neither the names of the **Blaze** development group nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

----

## Compiler Compatibility ##

**Blaze** supports the C++14 standard and is compatible with a wide range of C++ compilers. In fact, **Blaze** is constantly tested with the GNU compiler collection (version 6.0 through 9.1), the Clang compiler (version 5.0 through 8.0), and Visual C++ 2017 (Win64 only). Other compilers are not explicitly tested, but might work with a high probability.

If you are looking for a C++98 compatible math library you might consider using an older release of **Blaze**. Until the release 2.6 **Blaze** was written in C++-98 and constantly tested with the GNU compiler collection (version 4.5 through 5.0), the Intel C++ compiler (12.1, 13.1, 14.0, 15.0), the Clang compiler (version 3.4 through 3.7), and Visual C++ 2010, 2012, 2013, and 2015 (Win64 only).

----

## Publications ##

* K. Iglberger, G. Hager, J. Treibig, and U. Rüde: **Expression Templates Revisited: A Performance Analysis of Current Methodologies** ([Download](http://epubs.siam.org/sisc/resource/1/sjoce3/v34/i2/pC42_s1)). SIAM Journal on Scientific Computing, 34(2): C42--C69, 2012
* K. Iglberger, G. Hager, J. Treibig, and U. Rüde: **High Performance Smart Expression Template Math Libraries** ([Download](http://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=06266939)). Proceedings of the 2nd International Workshop on New Algorithms and Programming Models for the Manycore Era (APMM 2012) at HPCS 2012

----

## Contributions ##

[Klaus Iglberger](https://www.linkedin.com/in/klaus-iglberger-2133694/) -- Project initiator and main developer

[Georg Hager](http://www.rrze.uni-erlangen.de/wir-ueber-uns/organigramm/mitarbeiter/index.shtml/georg-hager.shtml) -- Performance analysis and optimization

[Christian Godenschwager](http://www10.informatik.uni-erlangen.de/~godenschwager/) -- Visual Studio 2010/2012/2013/2015 bug fixes and testing

Tobias Scharpff -- Sparse matrix multiplication algorithms

byzhang -- Bug fixes

Emerson Ferreira -- Bug fixes

Fabien Péan -- CMake support

Denis Demidov -- Export CMake package configuration

Jannik Schürg -- AVX-512 support and cache size detection for macOS in CMake

Marcin Copik -- CMake fixes

Hartmut Kaiser -- HPX backend

[Patrick Diehl](http://www.diehlpk.de/) -- Integration of HPX to the Blazemark and maintainer of the **Blaze** Fedora package

Mario Emmenlauer -- Blazemark extensions

Jeff Pollock -- CMake extensions

Darcy Beurle -- Integration of **Blaze** into the Compiler Explorer

Robert Schumacher -- CMake fixes

Jan Rudolph -- CMake fixes

Mikhail Katliar -- LAPACK extensions
