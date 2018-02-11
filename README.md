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
[![blaze-3.3.jpg](https://bitbucket.org/blaze-lib/blaze/wiki/images/blaze-3.3.jpg)](https://bitbucket.org/blaze-lib/blaze/downloads/blaze-3.3.tar.gz)
![white40x120.jpg](https://bitbucket.org/blaze-lib/blaze/wiki/images/white40x120.jpg)
[![blaze-docu-3.3.jpg](https://bitbucket.org/blaze-lib/blaze/wiki/images/blaze-docu-3.3.jpg)](https://bitbucket.org/blaze-lib/blaze/downloads/blaze-docu-3.3.tar.gz)

Older releases of **Blaze** can be found in the [downloads](https://bitbucket.org/blaze-lib/blaze/downloads) section or in our [release archive](https://bitbucket.org/blaze-lib/blaze/wiki/Release Archive).

----

## Blaze Projects ##

[BlazeIterative](https://github.com/tjolsen/BlazeIterative): A collection of iterative solvers (CG, BiCGSTAB, ...) for the **Blaze** library (Tyler Olsen)
[RcppBlaze](https://github.com/ChingChuan-Chen/RcppBlaze): A **Blaze** port for the R language (ChingChuan Chen)

----

## News ##

**11.02.2018**: Today we proudly present the next release of the **Blaze** library! In **Blaze** 3.3 we have focused on views. First of all, we have upgraded all existing views. Subvectors, submatrices, rows and columns have been introduced in the **Blaze** 1.x releases when we were still limited by C++03. Now we have completely reworked them with C++14 and have made them more flexible and more powerful than ever before. Most remarkable, it is now possible to configure views at compile time:

```
#!c++
blaze::DynamicVector<double,blaze::rowVector> x;
// ... Resizing and initialization

// Create a subvector from index 4 with a size of 12 (i.e. in the range [4..15]) (compile time arguments)
auto sv1 = subvector<4UL,12UL>( x );

// Create a subvector from index 8 with a size of 16 (i.e. in the range [8..23]) (runtime arguments)
auto sv2 = subvector( x, 8UL, 16UL );
```

```
#!c++
blaze::DynamicMatrix<double,blaze::rowMajor> A;
// ... Resizing and initialization

// Creating a dense submatrix of size 4x8, starting in row 3 and column 0 (compile time arguments)
auto sm1 = submatrix<3UL,0UL,4UL,8UL>( A );

// Creating a dense submatrix of size 8x16, starting in row 0 and column 4 (runtime arguments)
auto sm2 = submatrix( A, 0UL, 4UL, 8UL, 16UL );
```

Furthermore, we have introduced [element selections](https://bitbucket.org/blaze-lib/blaze/wiki/Element%20Selections), [row selections](https://bitbucket.org/blaze-lib/blaze/wiki/Row%20Selections), [column selections](https://bitbucket.org/blaze-lib/blaze/wiki/Column%20Selections) and [bands](https://bitbucket.org/blaze-lib/blaze/wiki/Bands) to enable almost every possible view on vectors and matrices:

```
#!c++
blaze::DynamicVector<double,blaze::rowVector> x;
blaze::DynamicMatrix<double,blaze::rowMajor> A;
// ... Resizing and initialization

// Selecting the elements 4, 6, 8, and 10 (compile time arguments)
auto e1 = elements<4UL,6UL,8UL,10UL>( x );

// Selecting the elements 3, 2, and 1 (runtime arguments via an initializer list)
auto e2 = elements( x, { 3UL, 2UL, 1UL } );

// Selecting the rows 4, 6, 8, and 10 (compile time arguments)
auto rs1 = rows<4UL,6UL,8UL,10UL>( A );

// Selecting the rows 3, 2, and 1 (runtime arguments via an initializer list)
auto rs2 = rows( A, { 3UL, 2UL, 1UL } );

// Selecting the columns 4, 6, 8, and 10 (compile time arguments)
auto cs1 = columns<4UL,6UL,8UL,10UL>( A );

// Selecting the columns 3, 2, and 1 (runtime arguments via an initializer list)
auto cs2 = columns( A, { 3UL, 2UL, 1UL } );

// Creating a reference to the 1st lower band of matrix A (compile time index)
auto b1 = band<-1L>( A );

// Creating a reference to the 2nd upper band of matrix A (runtime index)
auto b2 = band( A, 2L );
```

Also noteworthy are our new [HPX shared-memory parallelization backend](https://bitbucket.org/blaze-lib/blaze/wiki/HPX%20Parallelization) and the introduction of [vector norms](https://bitbucket.org/blaze-lib/blaze/wiki/Vector%20Operations#!norms) and [matrix norms](https://bitbucket.org/blaze-lib/blaze/wiki/Matrix%20Operations#!norms). We hope you enjoy this new release!

**18.8.2017**: Today, after nearly six month of hard work, we officially release **Blaze** 3.2! This version is dedicated to several of the most anticipated features: **Blaze** finally provides [CMake support](https://bitbucket.org/blaze-lib/blaze/wiki/Configuration%20and%20Installation) and an [advanced configuration system](https://bitbucket.org/blaze-lib/blaze/wiki/Configuration Files), which allows you to configure each single detail of **Blaze** from the command line. Additionally, **Blaze** finally provides complete support of AVX-512 and introduces the [`IdentityMatrix`](https://bitbucket.org/blaze-lib/blaze/wiki/Matrix Types#!identitymatrix) class. Furthermore, **Blaze** finally features [binary custom operations](https://bitbucket.org/blaze-lib/blaze/wiki/Vector and Matrix Customization#!custom-operations) and the [componentwise matrix multiplication (Schur Product)](https://bitbucket.org/blaze-lib/blaze/wiki/Matrix-Matrix Multiplication#!componentwise-multiplication-schur-product). Of course we have also spent time on a lot of smaller features and tweaked countless little details. We hope you enjoy this new release and the ton of new features.

We don't want to miss the opportunity to thank our many contributors: Thanks a lot for your efforts to make **Blaze** a better library!

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
    * [Vector/Vector Division](https://bitbucket.org/blaze-lib/blaze/wiki/Vector-Vector Division)
    * [Matrix/Vector Multiplication](https://bitbucket.org/blaze-lib/blaze/wiki/Matrix-Vector Multiplication)
    * [Matrix/Matrix Multiplication](https://bitbucket.org/blaze-lib/blaze/wiki/Matrix-Matrix Multiplication)
        * [Schur Product](https://bitbucket.org/blaze-lib/blaze/wiki/Matrix-Matrix Multiplication#!componentwise-multiplication-schur-product)
        * [Matrix Product](https://bitbucket.org/blaze-lib/blaze/wiki/Matrix-Matrix Multiplication#!matrix-product)
* [Shared-Memory Parallelization](https://bitbucket.org/blaze-lib/blaze/wiki/Shared Memory Parallelization)
    * [OpenMP Parallelization](https://bitbucket.org/blaze-lib/blaze/wiki/OpenMP Parallelization)
    * [C++11 Thread Parallelization](https://bitbucket.org/blaze-lib/blaze/wiki/Cpp Thread Parallelization)
    * [Boost Thread Parallelization](https://bitbucket.org/blaze-lib/blaze/wiki/Boost Thread Parallelization)
    * [HPX Parallelization](https://bitbucket.org/blaze-lib/blaze/wiki/HPX Parallelization)
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

**Blaze** supports the C++14 standard and is compatible with a wide range of C++ compilers. In fact, **Blaze** is constantly tested with the GNU compiler collection (version 5.0 through 7.2), the Intel C++ compiler (16.0 through 17.4), the Clang compiler (version 3.9 through 5.0), and Visual C++ 2015 and 2017 (Win64 only). Other compilers are not explicitly tested, but might work with a high probability.

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

Patrick Diehl -- Integration of HPX to the Blazemark
