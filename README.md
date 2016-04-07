![blaze300x150.jpg](https://bitbucket.org/blaze-lib/blaze/wiki/images/blaze300x150.jpg)

**Blaze** is an open-source, high-performance C++ math library for dense and sparse arithmetic. With its state-of-the-art *Smart Expression Template* implementation **Blaze** combines the elegance and ease of use of a domain-specific language with HPC-grade performance, making it one of the most intuitive and fastest C++ math libraries available.

The **Blaze** library offers ...

  * ... **high performance** through the integration of BLAS libraries and manually tuned HPC math kernels
  * ... **parallel execution** by OpenMP, C++11 threads and Boost threads
  * ... the **intuitive** and **easy to use** API of a domain specific language
  * ... **unified arithmetic** with dense and sparse vectors and matrices
  * ... **thoroughly tested** matrix and vector arithmetic
  * ... completely **portable**, **high quality** C++ source code

Get an impression of the clear but powerful syntax of **Blaze** in the [Getting Started](https://bitbucket.org/blaze-lib/blaze/wiki/Getting_Started) tutorial and of the impressive performance in the [Benchmarks](https://bitbucket.org/blaze-lib/blaze/wiki/Benchmarks) section.

----

## Download ##

![white20x120.jpg](https://bitbucket.org/blaze-lib/blaze/wiki/images/white20x120.jpg)
[![blaze-2.6.jpg](https://bitbucket.org/blaze-lib/blaze/wiki/images/blaze-2.6.jpg)](https://bitbucket.org/blaze-lib/blaze/downloads/blaze-2.6.tar.gz)
![white40x120.jpg](https://bitbucket.org/blaze-lib/blaze/wiki/images/white40x120.jpg)
[![blaze-docu-2.6.jpg](https://bitbucket.org/blaze-lib/blaze/wiki/images/blaze-docu-2.6.jpg)](https://bitbucket.org/blaze-lib/blaze/downloads/blaze-docu-2.6.tar.gz)

Older releases of **Blaze** can be found in the [downloads](https://bitbucket.org/blaze-lib/blaze/downloads) section or in our [release archive](https://bitbucket.org/blaze-lib/blaze/wiki/Release Archive).

----

## News ##

**16.2.2016**: Release 2.6 represents a major step in the history of the **Blaze** library. Within the last four and a half month we have implemented and brought to perfection several major features that have been on the wishlist of many **Blaze** users. Among these features are the LAPACK-based [LU](https://bitbucket.org/blaze-lib/blaze/wiki/Matrix Operations#!lu-decomposition), [Cholesky](https://bitbucket.org/blaze-lib/blaze/wiki/Matrix Operations#!cholesky-decomposition), [QR](https://bitbucket.org/blaze-lib/blaze/wiki/Matrix Operations#!qr-decomposition), [RQ](https://bitbucket.org/blaze-lib/blaze/wiki/Matrix Operations#!rq-decomposition), [QL](https://bitbucket.org/blaze-lib/blaze/wiki/Matrix Operations#!ql-decomposition), and [LQ](https://bitbucket.org/blaze-lib/blaze/wiki/Matrix Operations#!lq-decomposition) decomposition of dense matrices as well as the [inversion of dense matrices](https://bitbucket.org/blaze-lib/blaze/wiki/Matrix Operations#!matrix_inversion). Along with these come a large number of [LAPACK wrapper functions](https://bitbucket.org/blaze-lib/blaze/wiki/LAPACK Functions) that considerably simplify the use of LAPACK. Furthermore, the introduction of [`CustomVector`](https://bitbucket.org/blaze-lib/blaze/wiki/Vector Types#!customvector) and [`CustomMatrix`](https://bitbucket.org/blaze-lib/blaze/wiki/Matrix Types#!custommatrix) now enables custom memory management for both dense vector and matrices and allows to represent manually acquired memory as native **Blaze** data structure.

A major thank you goes to [CD-adapco](http://www.cd-adapco.com), a leader in the field of computational fluid dynamics (CFD). Since their decision to use **Blaze** within their flagship product [Star-CCM+](http://www.cd-adapco.com/products/star-ccm®) to drive linear algebra computations we have received considerable support and thus were able to completely focus on this release for the past few month. Thanks a lot for the ongoing help and support!

----

**1.10.2015**: After some month of hard work we finally release **Blaze** 2.5! This version comes with a fair amount of new features, many of them focused on computations with complex numbers: [Hermitian matrices](https://bitbucket.org/blaze-lib/blaze/issues/15/introduce-hermitian-matrices-into-blaze), the [`conj()`](https://bitbucket.org/blaze-lib/blaze/issues/7/provide-support-for-conjugate-complex), [`ctrans()`](https://bitbucket.org/blaze-lib/blaze/issues/8/provide-support-for-conjugate-transpose), [`real()`](https://bitbucket.org/blaze-lib/blaze/issues/17/provide-a-real-function-for-vectors-and), and [`imag()`](https://bitbucket.org/blaze-lib/blaze/issues/18/provide-an-imag-function-for-vectors-and) operations, and the [vectorization of integral complex numbers](https://bitbucket.org/blaze-lib/blaze/issues/21/support-the-vectorization-of-integral).

In addition, we are particularly proud about the performance improvements of several [dense matrix/sparse matrix multiplication](https://bitbucket.org/blaze-lib/blaze/issues/9/improve-the-performance-of-dense-matrix) and [sparse matrix/dense matrix multiplication](https://bitbucket.org/blaze-lib/blaze/issues/11/improve-the-performance-of-sparse-matrix) kernels. Furthermore, to extend the support for special environments, we have enabled the [customization of the error reporting mechanism](https://bitbucket.org/blaze-lib/blaze/issues/19/provide-means-to-adapt-the-error-reporting). Enjoy!

----

**4.7.2015**: We are proud to announce the release of **Blaze** 2.4! With this release, we say goodbye to the [GoogleCode](http://code.google.com/p/blaze-lib) platform and move to our new home:

![white250x100.jpg](https://bitbucket.org/blaze-lib/blaze/wiki/images/white250x100.jpg)
[![bitbucket300x100.jpg](https://bitbucket.org/blaze-lib/blaze/wiki/images/bitbucket300x100.jpg)](http://bitbucket.org/blaze-lib/blaze)

Since GoogleCode is shutting down, **Blaze** 2.4 will be the last release on [GoogleCode](http://code.google.com/p/blaze-lib). From now on we will focus entirely on our [Bitbucket](http://bitbucket.org/blaze-lib/blaze) repositories. Please update your bookmarks accordingly!

In **Blaze** 2.2 we introduced the first adaptor, `SymmetricMatrix`. In **Blaze** 2.3 we continued in this direction and introduced adaptors for lower and upper triangular matrices. Now, **Blaze** 2.4 concludes this line of development with unitriangular, strictly triangular and diagonal matrices. This now enables a total of 36 different matrix types, which provide the unique possibility to adapt the type of matrices exactly to the problem at hand and to achieve maximum performance due to specifically optimized kernels.

Additionally, in **Blaze** 2.4 we have removed a couple of hidden, unofficial features that caused trouble due to negligence and lack of proper testing: All solvers, `RotationMatrix` and `Quaternion` are not part of the release anymore. However, we plan to reintroduce both `RotationMatrix` and `Quaternion` in an upcoming release.

----

**24.3.2015**: You might have heard already: Google has officially announced to shut down the Google Code platform. This is of course bad news for us since this means that we have to move the **Blaze** website and repository to a new platform. However, since there have been several requests from the **Blaze** community to move **Blaze** to a different platform we have been thinking and working in this direction for some time anyway. Given the current situation we have decided to officially move our repository with the next release of **Blaze**: Version 2.4 will be the last release on Google Code and the first release on our new home. We are currently wrapping up the 2.4 release and at the same time preparing our move to the new platform. Stay tuned, we will inform you as soon as possible about the location of our new home.

----

**11.3.2015**: You were asking for them, here they are! The new big feature of **Blaze** 2.3 are lower and upper triangular matrices. And they come with full support within the **Blaze** library: Specifically adapted and optimized compute kernels and full support for parallelization via OpenMP, C++11 and Boost threads. Also noteworthy: The implementation is thoroughly checked by 748 new test cases. See the **Blaze** [tutorial](https://bitbucket.org/blaze-lib/blaze/wiki/Triangular_Matrices) for how the new `LowerMatrix` and `UpperMatrix` adaptors work!

----

**3.12.2014**: After a total of five and a half months, a little late for [SC'14](http://sc14.supercomputing.org), but right on time for [Meeting C++](http://meetingcpp.com), we finally release **Blaze** 2.2! But the waiting time was worthwhile! This release comes with several bug fixes and hundreds of improvements, many based on your hints, suggestions and ideas. Thank you very much for your support and help to make the **Blaze** library even better!

The big new feature of **Blaze** 2.2 is symmetric matrices. And this is not just any implementation of symmetric matrices, but one of the most complete and powerful implementations available. See the **Blaze** [tutorial](https://bitbucket.org/blaze-lib/blaze/wiki/Symmetric_Matrices) to get an idea of how symmetric matrices work and how they can help you prevent some inadvertent pessimizations of your code.

----

**20.6.2014**: After three month, just in time for [ISC'14](http://www.isc-events.com/isc14/), **Blaze** 2.1 has finally been released! The focus of this release is the improvement and extension of the shared memory parallelization capabilities. In addition to the OpenMP parallelization, **Blaze** 2.1 now also enables shared memory parallelization based on C++11 and Boost threads. With that, shared memory parallel execution is now available for virtually every platform and compiler. Additionally, the underlying architecture for the parallel execution and automatic, context-depending optimization has been significantly improved. Therefore **Blaze** 2.1 should prove to be the fastest and most flexible **Blaze** release yet.

----

**23.3.2014**: **Blaze** 2.0 goes parallel!

One of the main motivations of the **Blaze** 1.x releases was to provide maximum performance on a single CPU core for all possible operations. However, the [free lunch is over](http://www.gotw.ca/publications/concurrency-ddj.htm) and today's CPUs are not single core anymore. In order to fully utilize the performance potential of a multicore CPU, computations have to be parallelized across all available cores of a CPU. Therefore, starting with **Blaze** 2.0, the **Blaze** library provides automated shared memory parallelization. Get an impression of the possible performance boost in the [Benchmarks](https://bitbucket.org/blaze-lib/blaze/wiki/Benchmarks) section.

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
    * [Matrix/Vector Multiplication](https://bitbucket.org/blaze-lib/blaze/wiki/Matrix-Vector Multiplication)
    * [Matrix/Matrix Multiplication](https://bitbucket.org/blaze-lib/blaze/wiki/Matrix-Matrix Multiplication)
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

**Blaze** is written in C++-98 and therefore compatible with a wide range of C++ compilers. In fact, **Blaze** is constantly tested with the GNU compiler collection (version 4.5 through 5.0), the Intel C++ compiler (12.1, 13.1, 14.0, 15.0), the Clang compiler (version 3.4 through 3.7), and Visual C++ 2010, 2012, 2013, and 2015 (Win64 only). Other compilers are not explicitly tested, but might work with a high probability.

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
