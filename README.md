![blaze300x150.jpg](https://bitbucket.org/blaze-lib/blaze/wiki/images/blaze300x150.jpg)

**Blaze** is an open-source, high-performance C++ math library for dense and sparse arithmetic. With its state-of-the-art *Smart Expression Template* implementation **Blaze** combines the elegance and ease of use of a domain-specific language with HPC-grade performance, making it one of the most intuitive and fastest C++ math libraries available.

The **Blaze** library offers ...

  * ... **high performance** through the integration of BLAS libraries and manually tuned HPC math kernels
  * ... **vectorization** by SSE, SSE2, SSE3, SSSE3, SSE4, AVX, AVX2, AVX-512, FMA, and SVML
  * ... **parallel execution** by OpenMP, C++11 threads and Boost threads
  * ... the **intuitive** and **easy to use** API of a domain specific language
  * ... **unified arithmetic** with dense and sparse vectors and matrices
  * ... **thoroughly tested** matrix and vector arithmetic
  * ... completely **portable**, **high quality** C++ source code

Get an impression of the clear but powerful syntax of **Blaze** in the [Getting Started](https://bitbucket.org/blaze-lib/blaze/wiki/Getting_Started) tutorial and of the impressive performance in the [Benchmarks](https://bitbucket.org/blaze-lib/blaze/wiki/Benchmarks) section.

----

## Download ##

![white20x120.jpg](https://bitbucket.org/blaze-lib/blaze/wiki/images/white20x120.jpg)
[![blaze-3.1.jpg](https://bitbucket.org/blaze-lib/blaze/wiki/images/blaze-3.1.jpg)](https://bitbucket.org/blaze-lib/blaze/downloads/blaze-3.1.tar.gz)
![white40x120.jpg](https://bitbucket.org/blaze-lib/blaze/wiki/images/white40x120.jpg)
[![blaze-docu-3.1.jpg](https://bitbucket.org/blaze-lib/blaze/wiki/images/blaze-docu-3.1.jpg)](https://bitbucket.org/blaze-lib/blaze/downloads/blaze-docu-3.1.tar.gz)

Older releases of **Blaze** can be found in the [downloads](https://bitbucket.org/blaze-lib/blaze/downloads) section or in our [release archive](https://bitbucket.org/blaze-lib/blaze/wiki/Release Archive).

----

## Blaze Projects ##

[BlazeIterative](https://github.com/tjolsen/BlazeIterative): A collection of iterative solvers (CG, BiCGSTAB, ...) for the **Blaze** library (Tyler Olsen)

----

## News ##

**3.7.2017**: We are proud to announce the first of hopefully many **Blaze** projects: [BlazeIterative](https://github.com/tjolsen/BlazeIterative). Check out this collection of iterative solves, which neatly integrate with the **Blaze** library.

**20.2.2017**: We are happy to announce that there is now a port of the **Blaze** library for the R language available: [RcppBlaze](https://github.com/ChingChuan-Chen/RcppBlaze).

**18.2.2017**: Rejoice, **Blaze** 3.1 is finally online! And the waiting time was worthwhile: This new version comes with a bunch of new features. Amongst others **Blaze** 3.1 introduces the functionality to compute eigenvalues and eigenvectors and to perform a singular value decomposition. Most remarkable, however, are the new [`declsym()`](https://bitbucket.org/blaze-lib/blaze/wiki/Matrix%20Operations#!declsym), [`declherm()`](https://bitbucket.org/blaze-lib/blaze/wiki/Matrix%20Operations#!declherm), [`decllow()`](https://bitbucket.org/blaze-lib/blaze/wiki/Matrix%20Operations#!decllow), [`declupp()`](https://bitbucket.org/blaze-lib/blaze/wiki/Matrix%20Operations#!declupp), and [`decldiag()`](https://bitbucket.org/blaze-lib/blaze/wiki/Matrix%20Operations#!decldiag) functions. Via these functions it is now possible to explicitly declare a matrix or matrix expression as either symmetric, Hermitian, lower or upper triangular, or diagonal. And this can have significant effects on the performance. Consider for instance the following example:

```
#!c++
blaze::DynamicMatrix<double> A, B;
// ... Resizing and initialization

B = A * trans(A);             // The result of the matrix multiplication is symmetric, but a
                              // full matrix multiplication is performed
B = declsym( A * trans(A) );  // Declare the result of the matrix multiplication as symmetric
                              // i.e. perform an optimized matrix multiplication
```

In this example, the [`declsym()`](https://bitbucket.org/blaze-lib/blaze/wiki/Matrix%20Operations#!declsym) function unlocks the potential to perform an optimized matrix multiplication, exploiting the compile time information that the result of the multiplication will be a symmetric matrix. Especially large matrix multiplications can benefit tremendously.

**24.8.2016**: Exactly four years after the release of **Blaze** 1.0 we are very proud to release **Blaze** 3.0! For four years **Blaze** has been supporting C++98, and with that essentially every possible compiler and system. Also, within these four years we have managed to never introduce a breaking change in the released interface of **Blaze**. However, we believe it is time to upgrade to the new capabilities of the C++ standard. Now, for the first time we introduce a breaking change to the **Blaze** library by committing to C++14.

Our new release comes with a huge variety of new features. **Blaze** 3.0 provides support for fused multiply-add (FMA), the Intel SVML, and improved and extended support for AVX-512. Also, it introduces [custom operations](https://bitbucket.org/blaze-lib/blaze/wiki/Custom Operations) via the [`forEach()`](https://bitbucket.org/blaze-lib/blaze/wiki/Matrix%20Operations#!foreach) function. With this comes a huge variety of new componentwise operations: [`floor()`](https://bitbucket.org/blaze-lib/blaze/wiki/Matrix%20Operations#!floor-ceil), [`ceil()`](https://bitbucket.org/blaze-lib/blaze/wiki/Matrix%20Operations#!floor-ceil), [`sqrt()`](https://bitbucket.org/blaze-lib/blaze/wiki/Matrix%20Operations#!sqrt-invsqrt), [`invsqrt()`](https://bitbucket.org/blaze-lib/blaze/wiki/Matrix%20Operations#!sqrt-invsqrt), [`cbrt()`](https://bitbucket.org/blaze-lib/blaze/wiki/Matrix%20Operations#!cbrt-invcbrt), [`invcbrt()`](https://bitbucket.org/blaze-lib/blaze/wiki/Matrix%20Operations#!cbrt-invcbrt), [`clip()`](https://bitbucket.org/blaze-lib/blaze/wiki/Matrix%20Operations#!clip), [`pow()`](https://bitbucket.org/blaze-lib/blaze/wiki/Matrix%20Operations#!pow), [`exp()`](https://bitbucket.org/blaze-lib/blaze/wiki/Matrix%20Operations#!exp), [`log()`](https://bitbucket.org/blaze-lib/blaze/wiki/Matrix%20Operations#!log-log10), [`log10()`](https://bitbucket.org/blaze-lib/blaze/wiki/Matrix%20Operations#!log-log10), [`sin()`](https://bitbucket.org/blaze-lib/blaze/wiki/Matrix%20Operations#!sin-cos-tan-asin-acos-atan), [`asin()`](https://bitbucket.org/blaze-lib/blaze/wiki/Matrix%20Operations#!sin-cos-tan-asin-acos-atan), [`cos()`](https://bitbucket.org/blaze-lib/blaze/wiki/Matrix%20Operations#!sin-cos-tan-asin-acos-atan), [`acos()`](https://bitbucket.org/blaze-lib/blaze/wiki/Matrix%20Operations#!sin-cos-tan-asin-acos-atan), [`tan()`](https://bitbucket.org/blaze-lib/blaze/wiki/Matrix%20Operations#!sin-cos-tan-asin-acos-atan), [`atan()`](https://bitbucket.org/blaze-lib/blaze/wiki/Matrix%20Operations#!sin-cos-tan-asin-acos-atan), [`sinh()`](https://bitbucket.org/blaze-lib/blaze/wiki/Matrix%20Operations#!sinh-cosh-tanh-asinh-acosh-atanh), [`asinh()`](https://bitbucket.org/blaze-lib/blaze/wiki/Matrix%20Operations#!sinh-cosh-tanh-asinh-acosh-atanh), [`cosh()`](https://bitbucket.org/blaze-lib/blaze/wiki/Matrix%20Operations#!sinh-cosh-tanh-asinh-acosh-atanh), [`acosh()`](https://bitbucket.org/blaze-lib/blaze/wiki/Matrix%20Operations#!sinh-cosh-tanh-asinh-acosh-atanh), [`tanh()`](https://bitbucket.org/blaze-lib/blaze/wiki/Matrix%20Operations#!sinh-cosh-tanh-asinh-acosh-atanh), [`atanh()`](https://bitbucket.org/blaze-lib/blaze/wiki/Matrix%20Operations#!sinh-cosh-tanh-asinh-acosh-atanh), [`erf()`](https://bitbucket.org/blaze-lib/blaze/wiki/Matrix%20Operations#!sinh-cosh-tanh-asinh-acosh-atanh), [`erfc()`](https://bitbucket.org/blaze-lib/blaze/wiki/Matrix%20Operations#!sinh-cosh-tanh-asinh-acosh-atanh). Last but not least it introduces [vector-vector divisions](https://bitbucket.org/blaze-lib/blaze/wiki/Vector-Vector Division).

Along with the upgrade to C++14, **Blaze** 3.0 implements several additional, smaller changes to the interface: First, the direct initialization constructors for `StaticVector` and `StaticMatrix` have been removed. Instead, it is now possible to use initializer lists for the initialization of or the assignment to all dense vectors and matrices. Second, the `Dense` and `Sparse` prefixes of all views have been removed, which for instance means that `DenseSubvector` and `SparseSubvector` have been merged into the `Subvector` class template. And third, the `byDefault` inversion flag has been removed.

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
    * [Submatrices](https://bitbucket.org/blaze-lib/blaze/wiki/Submatrices)
    * [Rows](https://bitbucket.org/blaze-lib/blaze/wiki/Rows)
    * [Columns](https://bitbucket.org/blaze-lib/blaze/wiki/Columns)
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
* [Custom Operations](https://bitbucket.org/blaze-lib/blaze/wiki/Custom Operations)
* [Shared-Memory Parallelization](https://bitbucket.org/blaze-lib/blaze/wiki/Shared Memory Parallelization)
    * [OpenMP Parallelization](https://bitbucket.org/blaze-lib/blaze/wiki/OpenMP Parallelization)
    * [C++11 Thread Parallelization](https://bitbucket.org/blaze-lib/blaze/wiki/Cpp Thread Parallelization)
    * [Boost Thread Parallelization](https://bitbucket.org/blaze-lib/blaze/wiki/Boost Thread Parallelization)
    * [Serial Execution](https://bitbucket.org/blaze-lib/blaze/wiki/Serial Execution)
* [Serialization](https://bitbucket.org/blaze-lib/blaze/wiki/Serialization)
    * [Vector Serialization](https://bitbucket.org/blaze-lib/blaze/wiki/Vector Serialization)
    * [Matrix Serialization](https://bitbucket.org/blaze-lib/blaze/wiki/Matrix Serialization)
* [BLAS Functions](https://bitbucket.org/blaze-lib/blaze/wiki/BLAS Functions)
* [LAPACK Functions](https://bitbucket.org/blaze-lib/blaze/wiki/LAPACK Functions)
* [Configuration Files](https://bitbucket.org/blaze-lib/blaze/wiki/Configuration Files)
* [Block Vectors and Matrices](https://bitbucket.org/blaze-lib/blaze/wiki/Block Vectors and Matrices)
* [Custom Data Types](https://bitbucket.org/blaze-lib/blaze/wiki/Custom Data Types)
* [Error Reporting Customization](https://bitbucket.org/blaze-lib/blaze/wiki/Error Reporting Customization)
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

**Blaze** supports the C++14 standard and is compatible with a wide range of C++ compilers. In fact, **Blaze** is constantly tested with the GNU compiler collection (version 4.9 through 6.3), the Intel C++ compiler (16.0), the Clang compiler (version 3.7 through 4.0), and Visual C++ 2015 (Win64 only). Other compilers are not explicitly tested, but might work with a high probability.

If you are looking for a C++98 compatible math library you might consider using an older release of **Blaze**. Until the release 2.6 **Blaze** was written in C++-98 and constantly tested with the GNU compiler collection (version 4.5 through 5.0), the Intel C++ compiler (12.1, 13.1, 14.0, 15.0), the Clang compiler (version 3.4 through 3.7), and Visual C++ 2010, 2012, 2013, and 2015 (Win64 only).

----

## Publications ##

  * K. Iglberger, G. Hager, J. Treibig, and U. Rüde: **Expression Templates Revisited: A Performance Analysis of Current Methodologies** ([Download](http://epubs.siam.org/sisc/resource/1/sjoce3/v34/i2/pC42_s1)). SIAM Journal on Scientific Computing, 34(2): C42--C69, 2012
  * K. Iglberger, G. Hager, J. Treibig, and U. Rüde: **High Performance Smart Expression Template Math Libraries** ([Download](http://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=06266939)). Proceedings of the 2nd International Workshop on New Algorithms and Programming Models for the Manycore Era (APMM 2012) at HPCS 2012

----

## Contributions ##

[Klaus Iglberger](http://www.xing.com/profile/Klaus_Iglberger) -- Project initiator and main developer

[Georg Hager](http://www.rrze.uni-erlangen.de/wir-ueber-uns/organigramm/mitarbeiter/index.shtml/georg-hager.shtml) -- Performance analysis and optimization

[Christian Godenschwager](http://www10.informatik.uni-erlangen.de/~godenschwager/) -- Visual Studio 2010/2012/2013/2015 bug fixes and testing

Tobias Scharpff -- Sparse matrix multiplication algorithms

byzhang -- Bug fixes

Emerson Ferreira -- Bug fixes

Fabien Péan -- CMake support

Denis Demidov -- Export CMake package configuration

Jannik Schürg -- Cache size detection for macOS in CMake
