//=================================================================================================
/*!
//  \file blaze/Blaze.h
//  \brief Primary include file of the Blaze library
//
//  Copyright (C) 2013 Klaus Iglberger - All Rights Reserved
//
//  This file is part of the Blaze library. You can redistribute it and/or modify it under
//  the terms of the New (Revised) BSD License. Redistribution and use in source and binary
//  forms, with or without modification, are permitted provided that the following conditions
//  are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this list of
//     conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright notice, this list
//     of conditions and the following disclaimer in the documentation and/or other materials
//     provided with the distribution.
//  3. Neither the names of the Blaze development group nor the names of its contributors
//     may be used to endorse or promote products derived from this software without specific
//     prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
//  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
//  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
//  SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
//  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
//  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
//  BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
//  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
//  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
//  DAMAGE.
*/
//=================================================================================================

#ifndef _BLAZE_BLAZE_H_
#define _BLAZE_BLAZE_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/Math.h>
#include <blaze/Util.h>




//=================================================================================================
//
//  DOXYGEN DOCUMENTATION
//
//=================================================================================================

//*************************************************************************************************
//! Namespace of the \b Blaze C++ math library.
namespace blaze {}
//*************************************************************************************************


//**Mainpage***************************************************************************************
/*!\mainpage
//
// \image html blaze300x150.jpg
//
// This is the API for the \b Blaze high performance C++ math library. It gives a complete
// overview of the individual features and sublibraries of \b Blaze. To get a first impression
// on \b Blaze, the short \ref getting_started tutorial is a good place to start. Afterwards,
// the following long tutorial covers the most important aspects of the \b Blaze math library.
// The tabs at the top of the page allow a direct access to the individual modules, namespaces,
// classes, and files of the \b Blaze library.\n\n
//
// \section table_of_content Table of Contents
//
// <ul>
//    <li> \ref configuration_and_installation </li>
//    <li> \ref getting_started </li>
//    <li> Tutorial\n
//       <ul>
//          <li> Vectors
//             <ul>
//                <li> \ref vector_types </li>
//                <li> \ref vector_operations </li>
//             </ul>
//          </li>
//          <li> Matrices
//             <ul>
//                <li> \ref matrix_types </li>
//                <li> \ref matrix_operations </li>
//             </ul>
//          </li>
//          <li> Views
//             <ul>
//                <li> \ref views_subvectors </li>
//                <li> \ref views_submatrices </li>
//                <li> \ref views_rows </li>
//                <li> \ref views_columns </li>
//             </ul>
//          </li>
//          <li> Arithmetic Operations
//             <ul>
//                <li> \ref addition </li>
//                <li> \ref subtraction </li>
//                <li> \ref scalar_multiplication </li>
//                <li> \ref vector_vector_multiplication
//                   <ul>
//                      <li> \ref componentwise_multiplication </li>
//                      <li> \ref inner_product </li>
//                      <li> \ref outer_product </li>
//                      <li> \ref cross_product </li>
//                   </ul>
//                </li>
//                <li> \ref matrix_vector_multiplication </li>
//                <li> \ref matrix_matrix_multiplication </li>
//             </ul>
//          </li>
//          <li> Shared-Memory Parallelization
//             <ul>
//                <li> \ref openmp_parallelization </li>
//                <li> \ref serial_execution </li>
//             </ul>
//          </li>
//          <li> Serialization
//             <ul>
//                <li> \ref vector_serialization </li>
//                <li> \ref matrix_serialization </li>
//             </ul>
//          </li>
//       </ul>
//    </li>
//    <li> \ref intra_statement_optimization </li>
//    <li> \ref configuration_files </li>
// </ul>
*/
//*************************************************************************************************


//**Configuration and Installation*****************************************************************
/*!\page configuration_and_installation Configuration and Installation
//
// <center> Next: \ref getting_started </center> \n
//
// Setting up the \b Blaze library on a particular system is a fairly easy two step process. Since
// \b Blaze is a template library and therefore mainly consists of header files no compilation is
// required. In the following, this two step process is explained in detail, preceded only by a
// short summary of the requirements.
//
//
// \n \section requirements Requirements
// <hr>
//
// In order for \b Blaze to work properly, the Boost library must be installed on the system. It
// is recommended to use the newest Boost library available, but \b Blaze requires at minimum the
// Boost version 1.39. If you don't have Boost installed on your system, you can download it for
// free from 'http://www.boost.org'.
//
// Additionally, for maximum performance \b Blaze expects you to have a BLAS library installed
// (<a href="http://software.intel.com/en-us/articles/intel-mkl/">Intel MKL</a>,
// <a href="http://developer.amd.com/libraries/acml/">ACML</a>,
// <a href="http://math-atlas.sourceforge.net">Atlas</a>,
// <a href="http://www.tacc.utexas.edu/tacc-projects/gotoblas2">Goto</a>, ...). If you don't
// have a BLAS library installed on your system, \b Blaze will still work and will not be reduced
// in functionality, but performance may be severely limited. Thus it is strongly recommended to
// install a BLAS library.
//
//
// \n \section step_1_configuration Step 1: Configuration
// <hr>
//
// \subsection step_1_configuration_unix Linux/MacOSX User
//
// The first step is to adapt the \c Configfile in the \b Blaze home directory to the local
// configuration. Any text editor can be used for this task:

   \code
   vi ./Configfile
   \endcode

// In the \c Configfile, the kind of installation (debug or release), the library types (static
// and/or dynamic), the compiler including compiler flags, and several include paths have to be
// specified. Afterwards, the \c configure script can be run, which uses the \c Configfile to
// update and create several files:

   \code
   ./configure
   \endcode

// This step can also be omitted, but results in a default configuration that does not guarantee
// the highest performance for all operations. For instance, without running the \c configure
// script, \b Blaze assumes that no BLAS library is installed on the system and cannot use BLAS
// functionality for instance for the matrix/matrix multiplication.
//
// In order to further customize the \b Blaze library the header files in the <em>./blaze/config/</em>
// subdirectory can be adapted. See section \ref configuration_files for more details.
//
// \n \subsection step_1_configuration_windows Windows User
//
// Unfortunately, for Windows users there is no \c configure script available (yet). Therefore
// Windows user have to manually configure the \b Blaze library. Most configuration headers are
// located in the <em>./blaze/config/</em> subdirectory. The one exception is the \c BLAS.h
// header in the <em>./blaze/system/</em> subdirectory that contains the configuration of the
// BLAS functionality. Note that in case the \c BLAZE_BLAS_MODE symbol is set to 1, the correct
// BLAS header file has to be specified!
//
//
// \n \section step_2_installation Step 2: Installation
// <hr>
//
// \subsection step_2_configuration_unix Linux/MacOSX User
//
// The second step is the installation of the header files. Since \b Blaze mainly consists of
// header files, the <em>./blaze</em> subdirectory can be simply copied to a standard include
// directory (note that this requires root privileges):

   \code
   cp -r ./blaze /usr/local/include
   \endcode

// Alternatively, on Unix-based machines (which includes Linux and Mac OS X) the
// \c CPLUS_INCLUDE_PATH environment variable can be set. The specified directory will be
// searched after any directories specified on the command line with the option \c -I and
// before the standard default directories (such as \c /usr/local/include and \c /usr/include).
// Assuming a user misterX, the environment variable can be set as follows:

   \code
   CPLUS_INCLUDE_PATH=/usr/home/misterX/blaze
   export CPLUS_INCLUDE_PATH
   \endcode

// Last but not least, the <em>./blaze</em> subdirectory can be explicitly specified on the
// command line. The following example demonstrates this by means of the GNU C++ compiler:

   \code
   g++ -I/usr/home/misterX/blaze -o BlazeTest BlazeTest.cpp
   \endcode

// \n \subsection step_2_configuration_windows Windows User
//
// Windows doesn't have a standard include directory. Therefore the \b Blaze header files can be
// copied to any other directory or simply left in the default \b Blaze directory. However, the
// chosen include directory has to be explicitly specified as include path. In Visual Studio,
// this is done via the project property pages, configuration properties, C/C++, General settings.
// Here the additional include directories can be specified. Note that there are small differences
// between VS2008 and VS2010:
// <a href="http://blogs.msdn.com/b/vsproject/archive/2009/07/07/vc-directories.aspx">VC++ Directories</a>.
//
//
// \n \section step_3_compilation Step 3 (Optional): Compilation
// <hr>
//
// \subsection step_3_configuration_unix Linux/MacOSX User
//
// Next to the math library, \b Blaze also contains a small number of additional (sub-)libraries.
// If these libraries, such as blaze::ThreadPool or the blaze::logging functionality, are required
// it is necessary to create the \b Blaze library files. For that purpose, the \c configure script
// has created a \c Makefile that can be used for the compilation process:

   \code
   make
   \endcode

// Afterwards, the \c libblaze.so and/or \c libblaze.a libraries are contained in the \a lib
// subdirectory and can be copied to a standard library directory (note that this requires
// root privilages). However, this step can be omitted if only the \b Blaze math library is
// required.

   \code
   cp ./lib/* /usr/local/lib
   \endcode

// Alternatively, on Unix-based systems the \c LD_LIBRARY_PATH environment variable can be
// extended to also consider the \b Blaze \a lib directory:

   \code
   LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/usr/home/misterX/blaze/lib
   export LD_LIBRARY_PATH
   \endcode

// \n \subsection step_3_configuration_windows Windows User
//
// For Windows users, a comfortable compilation of the extended \b Blaze features is not (yet)
// supported.
//
// \n <center> Next: \ref getting_started </center>
*/
//*************************************************************************************************


//**Getting Started********************************************************************************
/*!\page getting_started Getting Started
//
// <center> Previous: \ref configuration_and_installation &nbsp; &nbsp; Next: \ref vector_types </center> \n
//
// This short tutorial serves the purpose to give a quick overview of the way mathematical
// expressions have to be formulated in \b Blaze. Starting with \ref vector_types, the following
// long tutorial covers the most important aspects of the \b Blaze math library.
//
//
// \n \section getting_started_vector_example A First Example
//
// \b Blaze is written such that using mathematical expressions is as close to mathematical
// textbooks as possible and therefore as intuitive as possible. In nearly all cases the seemingly
// easiest solution is the right solution and most users experience no problems when trying to
// use \b Blaze in the most natural way. The following example gives a first impression of the
// formulation of a vector addition in \b Blaze:

   \code
   #include <iostream>
   #include <blaze/Math.h>

   using blaze::StaticVector;
   using blaze::DynamicVector;

   // Instantiation of a static 3D column vector. The vector is directly initialized as
   //   ( 4 -2  5 )
   StaticVector<int,3UL> a( 4, -2, 5 );

   // Instantiation of a dynamic 3D column vector. Via the subscript operator the values are set to
   //   ( 2  5 -3 )
   DynamicVector<int> b( 3UL );
   b[0] = 2;
   b[1] = 5;
   b[2] = -3;

   // Adding the vectors a and b
   DynamicVector<int> c = a + b;

   // Printing the result of the vector addition
   std::cout << "c =\n" << c << "\n";
   \endcode

// Note that the entire \b Blaze math library can be included via the \c blaze/Math.h header
// file. Alternatively, the entire \b Blaze library, including both the math and the entire
// utility module, can be included via the \c blaze/Blaze.h header file. Also note that all
// classes and functions of \b Blaze are contained in the blaze namespace.\n\n
//
// Assuming that this program resides in a source file called \c FirstExample.cpp, it can be
// compiled for instance via the GNU C++ compiler:

   \code
   g++ -ansi -O3 -DNDEBUG -mavx -o FirstExample FirstExample.cpp
   \endcode

// Note the definition of the \c NDEBUG preprocessor symbol. In order to achieve maximum
// performance, it is necessary to compile the program in release mode, which deactivates
// all debugging functionality inside \b Blaze. It is also strongly recommended to specify
// the available architecture specific instruction set (as for instance the AVX instruction
// set, which if available can be activated via the \c -mavx flag). This allows \b Blaze
// to optimize computations via vectorization.\n\n
//
// When running the resulting executable \c FirstExample, the output of the last line of
// this small program is

   \code
   c =
   6
   3
   2
   \endcode

// \n \section getting_started_matrix_example An Example Involving Matrices
//
// Similarly easy and intuitive are expressions involving matrices:

   \code
   #include <blaze/Math.h>

   using namespace blaze;

   // Instantiating a dynamic 3D column vector
   DynamicVector<int> x( 3UL );
   x[0] =  4;
   x[1] = -1;
   x[2] =  3;

   // Instantiating a dynamic 2x3 row-major matrix, preinitialized with 0. Via the function call
   // operator three values of the matrix are explicitly set to get the matrix
   //   ( 1  0  4 )
   //   ( 0 -2  0 )
   DynamicMatrix<int> A( 2UL, 3UL, 0 );
   A(0,0) =  1;
   A(0,2) =  4;
   A(1,1) = -2;

   // Performing a matrix/vector multiplication
   DynamicVector<int> y = A * x;

   // Printing the resulting vector
   std::cout << "y =\n" << y << "\n";

   // Instantiating a static column-major matrix. The matrix is directly initialized as
   //   (  3 -1 )
   //   (  0  2 )
   //   ( -1  0 )
   StaticMatrix<int,3UL,2UL,columnMajor> B( 3, 0, -1, -1, 2, 0 );

   // Performing a matrix/matrix multiplication
   DynamicMatrix<int> C = A * B;

   // Printing the resulting matrix
   std::cout << "C =\n" << C << "\n";
   \endcode

// The output of this program is

   \code
   y =
   16
   2

   C =
   ( -1 -1 )
   (  0  4 )
   \endcode

// \n \section getting_started_complex_example A Complex Example
//
// The following example is much more sophisticated. It shows the implementation of the Conjugate
// Gradient (CG) algorithm (http://en.wikipedia.org/wiki/Conjugate_gradient) by means of the
// \b Blaze library:
//
// \image html cg.jpg
//
// In this example it is not important to understand the CG algorithm itself, but to see the
// advantage of the API of the \b Blaze library. In the \b Blaze implementation we will use a
// sparse matrix/dense vector multiplication for a 2D Poisson equation using \f$ N \times N \f$
// unknowns. It becomes apparent that the core of the algorithm is very close to the mathematical
// formulation and therefore has huge advantages in terms of readability and maintainability,
// while the performance of the code is close to the expected theoretical peak performance:

   \code
   const size_t NN( N*N );

   blaze::CompressedMatrix<double,rowMajor> A( NN, NN );
   blaze::DynamicVector<double,columnVector> x( NN, 1.0 ), b( NN, 0.0 ), r( NN ), p( NN ), Ap( NN );
   double alpha, beta, delta;

   // ... Initializing the sparse matrix A

   // Performing the CG algorithm
   r = b - A * x;
   p = r;
   delta = (r,r);

   for( size_t iteration=0UL; iteration<iterations; ++iteration )
   {
      Ap = A * p;
      alpha = delta / (p,Ap);
      x += alpha * p;
      r -= alpha * Ap;
      beta = (r,r);
      if( std::sqrt( beta ) < 1E-8 ) break;
      p = r + ( beta / delta ) * p;
      delta = beta;
   }
   \endcode

// \n Hopefully this short tutorial gives a good first impression of how mathematical expressions
// are formulated with \b Blaze. The following long tutorial, starting with \ref vector_types,
// will cover all aspects of the \b Blaze math library, i.e. it will introduce all vector and
// matrix types, all possible operations on vectors and matrices, and of course all possible
// mathematical expressions.
//
// \n <center> Previous: \ref configuration_and_installation &nbsp; &nbsp; Next: \ref vector_types </center>
*/
//*************************************************************************************************


//**Vector Types***********************************************************************************
/*!\page vector_types Vector Types
//
// <center> Previous: \ref getting_started &nbsp; &nbsp; Next: \ref vector_operations </center> \n
//
//
// \tableofcontents
//
//
// The \b Blaze library currently offers three dense vector types (\ref vector_types_static_vector,
// \ref vector_types_dynamic_vector and \ref vector_types_hybrid_vector) and one sparse vector type
// (\ref vector_types_compressed_vector). All vectors can be specified as either column vectors

                          \f$\left(\begin{array}{*{1}{c}}
                          1 \\
                          2 \\
                          3 \\
                          \end{array}\right)\f$

// or row vectors

                          \f$\left(\begin{array}{*{3}{c}}
                          1 & 2 & 3 \\
                          \end{array}\right).\f$

// Per default, all vectors in \b Blaze are column vectors.
//
//
// \n \section vector_types_static_vector StaticVector
// <hr>
//
// The blaze::StaticVector class template is the representation of a fixed-size vector with
// statically allocated elements of arbitrary type. It can be included via the header file

   \code
   #include <blaze/math/StaticVector.h>
   \endcode

// The type of the elements, the number of elements, and the transpose flag of the vector can
// be specified via the three template parameters:

   \code
   template< typename Type, size_t N, bool TF >
   class StaticVector;
   \endcode

//  - \c Type: specifies the type of the vector elements. StaticVector can be used with any
//             non-cv-qualified, non-reference, non-pointer element type.
//  - \c N   : specifies the total number of vector elements. It is expected that StaticVector is
//             only used for tiny and small vectors.
//  - \c TF  : specifies whether the vector is a row vector (\c blaze::rowVector) or a column
//             vector (\c blaze::columnVector). The default value is \c blaze::columnVector.
//
//
// \n \section vector_types_dynamic_vector DynamicVector
// <hr>
//
// The blaze::DynamicVector class template is the representation of an arbitrary sized vector
// with dynamically allocated elements of arbitrary type. It can be included via the header file

   \code
   #include <blaze/math/DynamicVector.h>
   \endcode

// The type of the elements and the transpose flag of the vector can be specified via the two
// template parameters:

   \code
   template< typename Type, bool TF >
   class DynamicVector;
   \endcode

//  - \c Type: specifies the type of the vector elements. DynamicVector can be used with any
//             non-cv-qualified, non-reference, non-pointer element type.
//  - \c TF  : specifies whether the vector is a row vector (\c blaze::rowVector) or a column
//             vector (\c blaze::columnVector). The default value is \c blaze::columnVector.
//
//
// \n \section vector_types_hybrid_vector HybridVector
// <hr>
//
// The blaze::HybridVector class template combines the advantages of the blaze::StaticVector and
// the blaze::DynamicVector class templates. It represents a fixed-size vector with statically
// allocated elements, but still can be dynamically resized (within the bounds of the available
// memory). It can be included via the header file

   \code
   #include <blaze/math/HybridVector.h>
   \endcode

// The type of the elements, the number of elements, and the transpose flag of the vector can
// be specified via the three template parameters:

   \code
   template< typename Type, size_t N, bool TF >
   class HybridVector;
   \endcode

//  - \c Type: specifies the type of the vector elements. HybridVector can be used with any
//             non-cv-qualified, non-reference, non-pointer element type.
//  - \c N   : specifies the maximum number of vector elements. It is expected that HybridVector
//             is only used for tiny and small vectors.
//  - \c TF  : specifies whether the vector is a row vector (\c blaze::rowVector) or a column
//             vector (\c blaze::columnVector). The default value is \c blaze::columnVector.
//
//
// \n \section vector_types_compressed_vector CompressedVector
// <hr>
//
// The blaze::CompressedVector class is the representation of an arbitrarily sized sparse
// vector, which stores only non-zero elements of arbitrary type. It can be included via the
// header file

   \code
   #include <blaze/math/CompressedVector.h>
   \endcode

// The type of the elements and the transpose flag of the vector can be specified via the two
// template parameters:

   \code
   template< typename Type, bool TF >
   class CompressedVector;
   \endcode

//  - \c Type: specifies the type of the vector elements. CompressedVector can be used with any
//             non-cv-qualified, non-reference, non-pointer element type.
//  - \c TF  : specifies whether the vector is a row vector (\c blaze::rowVector) or a column
//             vector (\c blaze::columnVector). The default value is \c blaze::columnVector.
//
//
// \n <center> Previous: \ref getting_started &nbsp; &nbsp; Next: \ref vector_operations </center>
*/
//*************************************************************************************************


//**Vector Operations******************************************************************************
/*!\page vector_operations Vector Operations
//
// <center> Previous: \ref vector_types &nbsp; &nbsp; Next: \ref matrix_types </center>
//
//
// \tableofcontents
//
//
// \n \section vector_operations_constructors Constructors
// <hr>
//
// Instantiating and setting up a vector is very easy and intuitive. However, there are a few
// rules to take care of:
//  - In case the last template parameter (the transpose flag) is omitted, the vector is per
//    default a column vector.
//  - The elements of a StaticVector are default initialized (i.e. built-in data types are
//    initialized to 0, class types are initialized via the default constructor).
//  - Newly allocated elements of a DynamicVector or CompressedVector remain uninitialized
//    if they are of built-in type and are default constructed if they are of class type.
//
// \n \subsection vector_operations_default_construction Default Construction

   \code
   using blaze::StaticVector;
   using blaze::DynamicVector;
   using blaze::CompressedVector;

   // All vectors can be default constructed. Whereas the size
   // of StaticVectors is fixed via the second template parameter,
   // the initial size of a default constructed DynamicVector or
   // CompressedVector is 0.
   StaticVector<int,2UL> v1;                // Instantiation of a 2D integer column vector.
                                            // All elements are initialized to 0.
   StaticVector<long,3UL,columnVector> v2;  // Instantiation of a 3D long integer column vector.
                                            // Again, all elements are initialized to 0L.
   DynamicVector<float> v3;                 // Instantiation of a dynamic single precision column
                                            // vector of size 0.
   DynamicVector<double,rowVector> v4;      // Instantiation of a dynamic double precision row
                                            // vector of size 0.
   CompressedVector<int> v5;                // Instantiation of a compressed integer column
                                            // vector of size 0.
   CompressedVector<double,rowVector> v6;   // Instantiation of a compressed double precision row
                                            // vector of size 0.
   \endcode

// \n \subsection vector_operations_size_construction Construction with Specific Size
//
// The DynamicVector and CompressedVector classes offer a constructor that allows to immediately
// give the vector the required size. Whereas DynamicVector uses this information to allocate
// memory for all vector elements, CompressedVector merely acquires the size but remains empty.

   \code
   DynamicVector<int,columnVector> v7( 9UL );   // Instantiation of an integer dynamic column vector
                                                // of size 9. The elements are NOT initialized!
   DynamicVector< complex<float> > v8( 2UL );   // Instantiation of a column vector with two single
                                                // precision complex values. The elements are
                                                // default constructed.
   CompressedVector<int,rowVector> v9( 10UL );  // Instantiation of a compressed row vector with
                                                // size 10. Initially, the vector provides no
                                                // capacity for non-zero elements.
   \endcode

// \n \subsection vector_operations_initialization_constructors Initialization Constructors
//
// All dense vector classes offer a constructor that allows for a direct, homogeneous initialization
// of all vector elements. In contrast, for sparse vectors the predicted number of non-zero elements
// can be specified

   \code
   StaticVector<int,3UL,rowVector> v10( 2 );            // Instantiation of a 3D integer row vector.
                                                        // All elements are initialized to 2.
   DynamicVector<float> v11( 3UL, 7.0F );               // Instantiation of a dynamic single precision
                                                        // column vector of size 3. All elements are
                                                        // set to 7.0F.
   CompressedVector<float,rowVector> v12( 15UL, 3UL );  // Instantiation of a single precision column
                                                        // vector of size 15, which provides enough
                                                        // space for at least 3 non-zero elements.
   \endcode

// The StaticVector class offers a special initialization constructor. For StaticVectors of up
// to 6 elements (i.e. 6D vectors) the vector elements can be individually specified in the
// constructor:

   \code
   using blaze::StaticVector;

   StaticVector<int,1UL>                v13( 4 );
   StaticVector<long,2UL>               v14( 1L, -2L );
   StaticVector<float,3UL,columnVector> v15( -0.1F, 4.2F, -7.1F );
   StaticVector<double,4UL,rowVector>   v16( 1.3, -0.4, 8.3, -1.2 );
   StaticVector<size_t,5UL>             v17( 3UL, 4UL, 1UL, 9UL, 4UL );
   StaticVector<long,6UL>               v18( 1L, 3L, -2L, 9L, 4L, -3L );
   \endcode

// \n \subsection vector_operations_array_construction Array Construction
//
// Alternatively, all dense vector classes offer a constructor for an initialization with a dynamic
// or static array. If the vector is initialized from a dynamic array, the constructor expects the
// actual size of the array as first argument, the array as second argument. In case of a static
// array, the fixed size of the array is used:

   \code
   const double array1* = new double[2];
   // ... Initialization of the dynamic array

   float array2[4] = { 1.0F, 2.0F, 3.0F, 4.0F };

   blaze::StaticVector<double,2UL> v1( 2UL, array1 );
   blaze::DynamicVector<float>     v2( array2 );

   delete[] array1;
   \endcode

// \n \subsection vector_operations_copy_construction Copy Construction
//
// All dense and sparse vectors can be created as the copy of any other dense or sparse vector
// with the same transpose flag (i.e. blaze::rowVector or blaze::columnVector).

   \code
   StaticVector<int,9UL,columnVector> v19( v7 );  // Instantiation of the dense column vector v19
                                                  // as copy of the dense column vector v7.
   DynamicVector<int,rowVector> v20( v9 );        // Instantiation of the dense row vector v20 as
                                                  // copy of the sparse row vector v9.
   CompressedVector<int,columnVector> v21( v1 );  // Instantiation of the sparse column vector v21
                                                  // as copy of the dense column vector v1.
   CompressedVector<float,rowVector> v22( v12 );  // Instantiation of the sparse row vector v22 as
                                                  // copy of the row vector v12.
   \endcode

// Note that it is not possible to create a StaticVector as a copy of a vector with a different
// size:

   \code
   StaticVector<int,5UL,columnVector> v23( v7 );  // Runtime error: Size does not match!
   StaticVector<int,4UL,rowVector> v24( v10 );    // Compile time error: Size does not match!
   \endcode

// \n \section vector_operations_assignment Assignment
// <hr>
//
// There are several types of assignment to dense and sparse vectors:
// \ref vector_operations_homogeneous_assignment, \ref vector_operations_array_assignment,
// \ref vector_operations_copy_assignment, and \ref vector_operations_compound_assignment.
//
// \n \subsection vector_operations_homogeneous_assignment Homogeneous Assignment
//
// Sometimes it may be necessary to assign the same value to all elements of a dense vector.
// For this purpose, the assignment operator can be used:

   \code
   blaze::StaticVector<int,3UL> v1;
   blaze::DynamicVector<double> v2;

   // Setting all integer elements of the StaticVector to 2
   v1 = 2;

   // Setting all double precision elements of the DynamicVector to 5.0
   v2 = 5.0;
   \endcode

// \n \subsection vector_operations_array_assignment Array Assignment
//
// Dense vectors can also be assigned a static array:

   \code
   blaze::StaticVector<float,2UL> v1;
   blaze::DynamicVector<double,rowVector> v2;

   float  array1[2] = { 1.0F, 2.0F };
   double array2[5] = { 2.1, 4.0, -1.7, 8.6, -7.2 };

   v1 = array1;
   v2 = array2;
   \endcode

// \n \subsection vector_operations_copy_assignment Copy Assignment
//
// For all vector types it is generally possible to assign another vector with the same transpose
// flag (i.e. blaze::columnVector or blaze::rowVector). Note that in case of StaticVectors, the
// assigned vector is required to have the same size as the StaticVector since the size of a
// StaticVector cannot be adapted!

   \code
   blaze::StaticVector<int,3UL,columnVector> v1;
   blaze::DynamicVector<int,columnVector>    v2( 3UL );
   blaze::DynamicVector<float,columnVector>  v3( 5UL );
   blaze::CompressedVector<int,columnVector> v4( 3UL );
   blaze::CompressedVector<float,rowVector>  v5( 3UL );

   // ... Initialization of the vectors

   v1 = v2;  // OK: Assignment of a 3D dense column vector to another 3D dense column vector
   v1 = v4;  // OK: Assignment of a 3D sparse column vector to a 3D dense column vector
   v1 = v3;  // Runtime error: Cannot assign a 5D vector to a 3D static vector
   v1 = v5;  // Compilation error: Cannot assign a row vector to a column vector
   \endcode

// \n \subsection vector_operations_compound_assignment Compound Assignment
//
// Next to plain assignment, it is also possible to use addition assignment, subtraction
// assignment, and multiplication assignment. Note however, that in contrast to plain assignment
// the size and the transpose flag of the vectors has be to equal in order to able to perform a
// compound assignment.

   \code
   blaze::StaticVector<int,5UL,columnVector>   v1;
   blaze::DynamicVector<int,columnVector>      v2( 5UL );
   blaze::CompressedVector<float,columnVector> v3( 7UL );
   blaze::DynamicVector<float,rowVector>       v4( 7UL );
   blaze::CompressedVector<float,rowVector>    v5( 7UL );

   // ... Initialization of the vectors

   v1 += v2;  // OK: Addition assignment between two column vectors of the same size
   v1 += v3;  // Runtime error: No compound assignment between vectors of different size
   v1 -= v4;  // Compilation error: No compound assignment between vectors of different transpose flag
   v4 *= v5;  // OK: Multiplication assignment between two row vectors of the same size
   \endcode

// \n \section vector_operations_common_vector_operations Common Vector Operations
// <hr>
//
// \subsection vector_operations_size Size of a Vector
//
// Via the \c size() function, the current size of a vector can be queried:

   \code
   // Instantiating a dynamic vector with size 10
   blaze::DynamicVector<int> v1( 10UL );
   v1.size();  // Returns 10

   // Instantiating a compressed vector with size 12 and capacity for 3 non-zero elements
   blaze::CompressedVector<double> v2( 12UL, 3UL );
   v2.size();  // Returns 12
   \endcode

// \n \subsection vector_operations_capacity Capacity of a Vector
//
// Via the \c capacity() function the internal capacity of a DynamicVector or CompressedVector
// can be queried. Note that the capacity of a vector doesn't have to be equal to the size
// of a vector. In case of a dense vector the capacity will always be greater or equal than
// the size of the vector, in case of a sparse vector the capacity may even be less than
// the size.

   \code
   blaze::DynamicVector<int> v1( 10UL );
   v1.capacity();  // returns at least 10
   \endcode

// \n \subsection vector_operations_nonzeros Number of Non-Zero Elements
//
// For both dense and sparse vectors the number of non-zero elements can be determined via the
// \c nonZeros() function. Sparse vectors directly return their number of non-zero elements,
// dense vectors traverse their elements and count the number of non-zero elements.

   \code
   blaze::DynamicVector<int> v1( 10UL );
   blaze::CompressedVector<double> v2( 20UL );

   // ... Initializing the vectors

   v1.nonZeros();  // Returns the number of non-zero elements in the dense vector
   v2.nonZeros();  // Returns the number of non-zero elements in the sparse vector
   \endcode

// \n \section vector_operations_resize_reserve Resize/Reserve
// <hr>
//
// The size of a StaticVector is fixed by the second template parameter. In contrast, the size
// of DynamicVectors as well as CompressedVectors can be changed via the \c resize() function:

   \code
   using blaze::DynamicVector;
   using blaze::CompressedVector;

   DynamicVector<int,columnVector> v1;
   CompressedVector<int,rowVector> v2( 4 );
   v2[1] = -2;
   v2[3] = 11;

   // Adapting the size of the dynamic and compressed vectors. The (optional) second parameter
   // specifies whether the existing elements should be preserved. Per default, the existing
   // elements are not preserved.
   v1.resize( 5UL );         // Resizing vector v1 to 5 elements. Elements of built-in type remain
                             // uninitialized, elements of class type are default constructed.
   v1.resize( 3UL, false );  // Resizing vector v1 to 3 elements. The old elements are lost, the
                             // new elements are NOT initialized!
   v2.resize( 8UL, true );   // Resizing vector v2 to 8 elements. The old elements are preserved.
   v2.resize( 5UL, false );  // Resizing vector v2 to 5 elements. The old elements are lost.
   \endcode

// Note that resizing a vector invalidates all existing views (see e.g. \ref views_subvectors)
// on the vector:

   \code
   typedef blaze::DynamicVector<int,rowVector>  VectorType;
   typedef blaze::DenseSubvector<VectorType>    SubvectorType;

   VectorType v1( 10UL );                         // Creating a dynamic vector of size 10
   SubvectorType sv = subvector( v1, 2UL, 5UL );  // Creating a view on the range [2..6]
   v1.resize( 6UL );                              // Resizing the vector invalidates the view
   \endcode

// When the internal capacity of a vector is no longer sufficient, the allocation of a larger
// junk of memory is triggered. In order to avoid frequent reallocations, the \c reserve()
// function can be used up front to set the internal capacity:

   \code
   blaze::DynamicVector<int> v1;
   v1.reserve( 100 );
   v1.size();      // Returns 0
   v1.capacity();  // Returns at least 100
   \endcode
//
// Note that the size of the vector remains unchanged, but only the internal capacity is set
// according to the specified value!
//
//
// \n \section vector_operations_element_access Element Access
// <hr>
//
// The easiest and most intuitive way to access a dense or sparse vector is via the subscript
// operator. The indices to access a vector are zero-based:

   \code
   blaze::DynamicVector<int> v1( 5UL );
   v1[0] = 1;
   v1[1] = 3;
   // ...

   blaze::CompressedVector<float> v2( 5UL );
   v2[2] = 7.3F;
   v2[4] = -1.4F;
   \endcode

// Whereas using the subscript operator on a dense vector only accesses the already existing
// element, accessing an element of a sparse vector via the subscript operator potentially
// inserts the element into the vector and may therefore be more expensive. Consider the
// following example:

   \code
   blaze::CompressedVector<int> v1( 10UL );

   for( size_t i=0UL; i<v1.size(); ++i ) {
      ... = v1[i];
   }
   \endcode

// Although the compressed vector is only used for read access within the for loop, using the
// subscript operator temporarily inserts 10 non-zero elements into the vector. Therefore, all
// vectors (sparse as well as dense) offer an alternate way via the \c begin() and \c end()
// functions to traverse only the currently contained elements by iterators. In case of
// non-const vectors, \c begin() and \c end() return an Iterator, which allows a manipulation
// of the non-zero value, in case of a constant vector a ConstIterator is returned:

   \code
   using blaze::CompressedVector;

   CompressedVector<int> v1( 10UL );

   // ... Initialization of the vector

   // Traversing the vector by Iterator
   for( CompressedVector<int>::Iterator it=v1.begin(); it!=v1.end(); ++it ) {
      it->value() = ...;  // OK: Write access to the value of the non-zero element.
      ... = it->value();  // OK: Read access to the value of the non-zero element.
      it->index() = ...;  // Compilation error: The index of a non-zero element cannot be changed.
      ... = it->index();  // OK: Read access to the index of the non-zero element.
   }

   // Traversing the vector by ConstIterator
   for( CompressedVector<int>::ConstIterator it=v1.begin(); it!=v1.end(); ++it ) {
      it->value() = ...;  // Compilation error: Assignment to the value via a ConstIterator is invalid.
      ... = it->value();  // OK: Read access to the value of the non-zero element.
      it->index() = ...;  // Compilation error: The index of a non-zero element cannot be changed.
      ... = it->index();  // OK: Read access to the index of the non-zero element.
   }
   \endcode

// \n \section vector_operations_element_insertion Element Insertion
// <hr>
//
// In contrast to dense vectors, that store all elements independent of their value and that
// offer direct access to all elements, spares vectors only store the non-zero elements contained
// in the vector. Therefore it is necessary to explicitly add elements to the vector. The first
// option to add elements to a sparse vector is the subscript operator:

   \code
   using blaze::CompressedVector;

   CompressedVector<int> v1( 3UL );
   v1[1] = 2;
   \endcode

// In case the element at the given index is not yet contained in the vector, it is automatically
// inserted. Otherwise the old value is replaced by the new value 2. The operator returns a
// reference to the sparse vector element.\n
// However, insertion of elements can be better controlled via the \c insert() function. In contrast
// to the subscript operator it emits an exception in case the element is already contained in
// the matrix. In order to check for this case, the \c find() function can be used:

   \code
   // In case the element at index 4 is not yet contained in the matrix it is inserted
   // with a value of 6.
   if( v1.find( 4 ) == v1.end() )
      v1.insert( 4, 6 );
   \endcode

// Although the \c insert() function is very flexible, due to performance reasons it is not suited
// for the setup of large sparse vectors. A very efficient, yet also very low-level way to fill
// a sparse vector is the \c append() function. It requires the sparse vector to provide enough
// capacity to insert a new element. Additionally, the index of the new element must be larger
// than the index of the previous element. Violating these conditions results in undefined
// behavior!

   \code
   v1.reserve( 5 );     // Reserving space for 5 non-zero elements
   v1.append( 5, -2 );  // Appending the element -2 at index 5
   v1.append( 6,  4 );  // Appending the element 4 at index 6
   // ...
   \endcode

// \n \section vector_operations_reset_clear Reset/Clear
// <hr>
//
// In order to reset all elements of a vector, the \c reset() function can be used:

   \code
   // Setup of a single precision column vector, whose elements are initialized with 2.0F.
   blaze::DynamicVector<float> v1( 3UL, 2.0F );

   // Resetting all elements to 0.0F. Only the elements are reset, the size of the vector is unchanged.
   reset( v1 );  // Resetting all elements
   v1.size();    // Returns 3: size and capacity remain unchanged
   \endcode

// In order to return a vector to its default state (i.e. the state of a default constructed
// vector), the \c clear() function can be used:

   \code
   // Setup of a single precision column vector, whose elements are initialized with -1.0F.
   blaze::DynamicVector<float> v1( 5, -1.0F );

   // Resetting the entire vector.
   clear( v1 );  // Resetting the entire vector
   v1.size();    // Returns 0: size is reset, but capacity remains unchanged
   \endcode

// Note that resetting or clearing both dense and sparse vectors does not change the capacity
// of the vectors.
//
//
// \n \section vector_operations_vector_transpose Vector Transpose
// <hr>
//
// As already mentioned, vectors can either be column vectors (blaze::columnVector) or row vectors
// (blaze::rowVector). A column vector cannot be assigned to a row vector and vice versa. However,
// vectors can be transposed via the \c trans() function:

   \code
   blaze::DynamicVector<int,columnVector> v1( 4UL );
   blaze::CompressedVector<int,rowVector> v2( 4UL );

   v1 = v2;            // Compilation error: Cannot assign a row vector to a column vector
   v1 = trans( v2 );   // OK: Transposing the row vector to a column vector and assigning it
                       //     to the column vector v1
   v2 = trans( v1 );   // OK: Transposing the column vector v1 and assigning it to the row vector v2
   v1 += trans( v2 );  // OK: Addition assignment of two column vectors
   \endcode

// \n \section vector_operations_length Vector Length
// <hr>
//
// In order to calculate the length of a vector, both the \c length() and \c sqrLength() function
// can be used:

   \code
   blaze::StaticVector<float,3UL,rowVector> v( -1.2F, 2.7F, -2.3F );

   const float len    = length   ( v );  // Computes the current length of the vector
   const float sqrlen = sqrLength( v );  // Computes the square length of the vector
   \endcode

// Note that both functions can only be used for vectors with built-in or complex element type!
//
// \n \section vector_operators_abs Absolute Values
// <hr>
//
// The \c abs() function can be used to compute the absolute values of each element of a vector.
// For instance, the following computation

   \code
   blaze::StaticVector<int,3UL,rowVector> a( -1, 2, -3 );
   blaze::StaticVector<int,3UL,rowVector> b( abs( a ) );
   \endcode

// results in the vector

                          \f$ b = \left(\begin{array}{*{1}{c}}
                          1 \\
                          2 \\
                          3 \\
                          \end{array}\right)\f$

// \n \section vector_operations_normalize Normalize
// <hr>
//
// The \c normalize() function can be used to scale any non-zero vector to a length of 1. In
// case the vector does not contain a single non-zero element (i.e. is a zero vector), the
// \c normalize() function returns a zero vector.

   \code
   blaze::DynamicVector<float,columnVector>     v1( 10UL );
   blaze::CompressedVector<double,columnVector> v2( 12UL );

   v1 = normalize( v1 );  // Normalizing the dense vector v1
   length( v1 );          // Returns 1 (or 0 in case of a zero vector)
   v1 = normalize( v2 );  // Assigning v1 the normalized vector v2
   length( v1 );          // Returns 1 (or 0 in case of a zero vector)
   \endcode

// Note that the \c normalize() function only works for floating point vectors. The attempt to
// use it for an integral vector results in a compile time error.
//
// \n \section vector_operations_swap Swap
// <hr>
//
// Via the \c swap() function it is possible to completely swap the contents of two vectors of
// the same type:

   \code
   blaze::DynamicVector<int,columnVector> v1( 10UL );
   blaze::DynamicVector<int,columnVector> v2( 20UL );

   swap( v1, v2 );  // Swapping the contents of v1 and v2
   \endcode

// \n <center> Previous: \ref vector_types &nbsp; &nbsp; Next: \ref matrix_types </center>
*/
//*************************************************************************************************


//**Matrix Types***********************************************************************************
/*!\page matrix_types Matrix Types
//
// <center> Previous: \ref vector_operations &nbsp; &nbsp; Next: \ref matrix_operations </center> \n
//
//
// \tableofcontents
//
//
// The \b Blaze library currently offers two dense matrix types (\ref matrix_types_static_matrix
// and \ref matrix_types_dynamic_matrix) and one sparse matrix type (\ref matrix_types_compressed_matrix).
// All matrices can either be stored as row-major matrices or column-major matrices. Per default,
// all matrices in \b Blaze are row-major matrices.
//
//
// \n \section matrix_types_static_matrix StaticMatrix
// <hr>
//
// The blaze::StaticMatrix class template is the representation of a fixed-size matrix with
// statically allocated elements of arbitrary type. It can be included via the header file

   \code
   #include <blaze/math/StaticMatrix.h>
   \endcode

// The type of the elements, the number of rows and columns, and the storage order of the matrix
// can be specified via the four template parameters:

   \code
   template< typename Type, size_t M, size_t N, bool SO >
   class StaticMatrix;
   \endcode

//  - \c Type: specifies the type of the matrix elements. StaticMatrix can be used with any
//             non-cv-qualified, non-reference element type.
//  - \c M   : specifies the total number of rows of the matrix.
//  - \c N   : specifies the total number of columns of the matrix. Note that it is expected
//             that StaticMatrix is only used for tiny and small matrices.
//  - \c SO  : specifies the storage order (blaze::rowMajor, blaze::columnMajor) of the matrix.
//             The default value is blaze::rowMajor.
//
//
// \n \section matrix_types_dynamic_matrix DynamicMatrix
// <hr>
//
// The blaze::DynamicMatrix class template is the representation of an arbitrary sized matrix
// with \f$ M \cdot N \f$ dynamically allocated elements of arbitrary type. It can be included
// via the header file

   \code
   #include <blaze/math/DynamicMatrix.h>
   \endcode

// The type of the elements and the storage order of the matrix can be specified via the two
// template parameters:

   \code
   template< typename Type, bool SO >
   class DynamicMatrix;
   \endcode

//  - \c Type: specifies the type of the matrix elements. DynamicMatrix can be used with any
//             non-cv-qualified, non-reference element type.
//  - \c SO  : specifies the storage order (blaze::rowMajor, blaze::columnMajor) of the matrix.
//             The default value is blaze::rowMajor.
//
//
// \n \section matrix_types_compressed_matrix CompressedMatrix
// <hr>
//
// The blaze::CompressedMatrix class template is the representation of an arbitrary sized sparse
// matrix with \f$ M \cdot N \f$ dynamically allocated elements of arbitrary type. It can be
// included via the header file

   \code
   #include <blaze/math/CompressedMatrix.h>
   \endcode

// The type of the elements and the storage order of the matrix can be specified via the two
// template parameters:

   \code
   template< typename Type, bool SO >
   class CompressedMatrix;
   \endcode

//  - \c Type: specifies the type of the matrix elements. CompressedMatrix can be used with
//             any non-cv-qualified, non-reference, non-pointer element type.
//  - \c SO  : specifies the storage order (blaze::rowMajor, blaze::columnMajor) of the matrix.
//             The default value is blaze::rowMajor.
//
//
// \n <center> Previous: \ref vector_operations &nbsp; &nbsp; Next: \ref matrix_operations </center>
*/
//*************************************************************************************************


//**Matrix Operations******************************************************************************
/*!\page matrix_operations Matrix Operations
//
// <center> Previous: \ref matrix_types &nbsp; &nbsp; Next: \ref views_subvectors </center>
//
//
// \tableofcontents
//
//
// \n \section matrix_operations_constructors Constructors
// <hr>
//
// Matrices are just as easy and intuitive to create as vectors. Still, there are a few rules
// to be aware of:
//  - In case the last template parameter (the storage order) is omitted, the matrix is per
//    default stored in row-major order.
//  - The elements of a StaticMatrix are default initialized (i.e. built-in data types are
//    initialized to 0, class types are initialized via the default constructor).
//  - Newly allocated elements of a DynamicMatrix or CompressedMatrix remain uninitialized
//    if they are of built-in type and are default constructed if they are of class type.
//
// \n \subsection matrix_operations_default_construction Default Construction

   \code
   using blaze::StaticMatrix;
   using blaze::DynamicMatrix;
   using blaze::CompressedMatrix;

   // All matrices can be default constructed. Whereas the size of
   // a StaticMatrix is fixed via the second and third template
   // parameter, the initial size of a constructed DynamicMatrix
   // or CompressedMatrix is 0.
   StaticMatrix<int,2UL,2UL> M1;           // Instantiation of a 2x2 single precision row-major
                                           // matrix. All elements are initialized to 0.
   DynamicMatrix<float> M2;                // Instantiation of a single precision dynamic
                                           // row-major matrix with 0 rows and 0 columns.
   DynamicMatrix<double,columnMajor> M3;   // Instantiation of a double precision dynamic
                                           // column-major matrix with 0 rows and 0 columns.
   CompressedMatrix<int> M4;               // Instantiation of a compressed integer
                                           // row-major matrix of size 0x0.
   CompressedMatrix<double,rowVector> M5;  // Instantiation of a compressed double precision
                                           // column-major matrix of size 0x0.
   \endcode

// \n \subsection matrix_operations_size_construction Construction with Specific Size
//
// The DynamicMatrix and CompressedMatrix classes offer a constructor that allows to immediately
// give the matrices a specific number of rows and columns:

   \code
   DynamicMatrix<int> M6( 5UL, 4UL );                 // Instantiation of a 5x4 dynamic row-major
                                                      // matrix. The elements are not initialized.
   DynamicMatrix<double,columnMajor> M7( 3UL, 7UL );  // Instantiation of a 3x7 dynamic column-major
                                                      // matrix. The elements are not initialized.
   CompressedMatrix<float,rowMajor> M8( 8UL, 6UL );   // Instantiation of a 8x6 compressed row-major
                                                      // matrix. The elements are not initialized.
   \endcode

// Note that dense matrices (in this case DynamicMatrix) immediately allocate enough capacity
// for all matrix elements. Sparse matrices on the other hand (in this example CompressedMatrix)
// merely acquire the size, but don't necessarily allocate memory.
//
//
// \n \subsection matrix_operations_initialization_constructors Initialization Constructors
//
// All dense matrix classes offer a constructor for a direct, homogeneous initialization of all
// matrix elements. In contrast, for sparse matrices the predicted number of non-zero elements
// can be specified.

   \code
   StaticMatrix<int,4UL,3UL,columnMajor> M9( 7 );  // Instantiation of a 4x3 integer row-major
                                                    // matrix. All elements are initialized to 7.
   DynamicMatrix<float> M10( 2UL, 5UL, 2.0F );      // Instantiation of a 2x5 single precision row-major
                                                    // matrix. All elements are initialized to 2.0F.
   CompressedMatrix<int> M11( 3UL, 4UL, 4 );        // Instantiation of a 3x4 integer column-major
                                                    // matrix. All elements are initialized to 4.
   \endcode

// The StaticMatrix class offers a special initialization constructor. For StaticMatrix of up
// to 10 elements the vector elements can be individually specified in the constructor:

   \code
   using blaze::StaticMatrix;

   StaticMatrix<int,3UL,1UL>               M12( 2, 5, -1 );
   StaticMatrix<float,2UL,3UL,columnMajor> M13( -0.1F, 4.2F, -7.1F,
                                                -0.8F, 1.3F,  4.2F );
   StaticMatrix<double,3UL,3UL,rowVector>  M14( 1.3, -0.4,  8.3,
                                                0.2, -1.5, -2.6,
                                                1.3,  9.3, -7.1 );
   \endcode

// \n \subsection matrix_operations_array_construction Array Construction
//
// Alternatively, all dense matrix classes offer a constructor for an initialization with a
// dynamic or static array. If the matrix is initialized from a dynamic array, the constructor
// expects the dimensions of values provided by the array as first and second argument, the
// array as third argument. In case of a static array, the fixed size of the array is used:

   \code
   const double array1* = new double[6];
   // ... Initialization of the dynamic array

   float array2[3][2] = { { 3.1F, 6.4F }, { -0.9F, -1.2F }, { 4.8F, 0.6F } };

   blaze::StaticMatrix<double,2UL,3UL> v1( 2UL, 3UL, array1 );
   blaze::DynamicMatrix<float>         v2( array2 );

   delete[] array1;
   \endcode

// \n \subsection matrix_operations_copy_construction Copy Construction
//
// All dense and sparse matrices can be created as a copy of another dense or sparse matrix.

   \code
   StaticMatrix<int,5UL,4UL,rowMajor> M15( M6 );  // Instantiation of the dense row-major matrix M15
                                                  // as copy of the dense row-major matrix M7.
   DynamicMatrix<int,columnMajor> M16( M8 );      // Instantiation of the dense column-major matrix M16
                                                  // as copy of the sparse row-major matrix M9.
   CompressedMatrix<double,rowMajor> M17( M7 );   // Instantiation of the compressed row-major matrix
                                                  // M17 as copy of the dense column-major matrix M8.
   CompressedMatrix<float,rowMajor> M18( M8 );    // Instantiation of the compressed row-major matrix
                                                  // M18 as copy of the compressed row-major matrix M9.
   \endcode

// Note that it is not possible to create a StaticMatrix as a copy of a matrix with a different
// number of rows and/or columns:

   \code
   StaticMatrix<int,4UL,5UL,rowMajor> M19( M6 );     // Runtime error: Number of rows and columns
                                                     // does not match!
   StaticMatrix<int,4UL,4UL,columnMajor> M20( M9 );  // Compile time error: Number of columns does
                                                     // not match!
   \endcode

// \n \section matrix_operations_assignment Assignment
// <hr>
//
// There are several types of assignment to dense and sparse matrices:
// \ref matrix_operations_homogeneous_assignment, \ref matrix_operations_array_assignment,
// \ref matrix_operations_copy_assignment, and \ref matrix_operations_compound_assignment.
//
// \n \subsection matrix_operations_homogeneous_assignment Homogeneous Assignment
//
// It is possible to assign the same value to all elements of a dense matrix. All dense matrix
// classes provide an according assignment operator:

   \code
   blaze::StaticMatrix<int,3UL,2UL> M1;
   blaze::DynamicMatrix<double> M2;

   // Setting all integer elements of the StaticMatrix to 4
   M1 = 4;

   // Setting all double precision elements of the DynamicMatrix to 3.5
   M2 = 3.5
   \endcode

// \n \subsection matrix_operations_array_assignment Array Assignment
//
// Dense matrices can also be assigned a static array:

   \code
   blaze::StaticMatrix<int,2UL,2UL,rowMajor> M1;
   blaze::StaticMatrix<int,2UL,2UL,columnMajor> M2;
   blaze::DynamicMatrix<double> M3;

   int array1[2][2] = { { 1, 2 }, { 3, 4 } };
   double array2[3][2] = { { 3.1, 6.4 }, { -0.9, -1.2 }, { 4.8, 0.6 } };

   M1 = array1;
   M2 = array1;
   M3 = array2;
   \endcode

// Note that due to the different storage order, the matrix M1 is initialized differently than
// matrix M2:

                          \f$ M1 = \left(\begin{array}{*{2}{c}}
                          1 & 2 \\
                          3 & 4 \\
                          \end{array}\right),\quad
                          M2 = \left(\begin{array}{*{2}{c}}
                          1 & 3 \\
                          2 & 4 \\
                          \end{array}\right)\f$

// Also note that the dimensions of the static array have to match the size of a StaticMatrix,
// whereas a DynamicMatrix is resized according to the array dimensions:

                          \f$ M1 = \left(\begin{array}{*{2}{c}}
                           3.1 &  6.4 \\
                          -0.9 & -1.2 \\
                           4.8 &  0.6 \\
                          \end{array}\right)\f$

// \n \subsection matrix_operations_copy_assignment Copy Assignment
//
// All kinds of matrices can be assigned to each other. The only restriction is that since a
// StaticMatrix cannot change its size, the assigned matrix must match both in the number of
// rows and in the number of columns.

   \code
   blaze::StaticMatrix<int,3UL,2UL,rowMajor>  M1;
   blaze::DynamicMatrix<int,rowMajor>         M2( 3UL, 2UL );
   blaze::DynamicMatrix<float,rowMajor>       M3( 5UL, 2UL );
   blaze::CompressedMatrix<int,rowMajor>      M4( 3UL, 2UL );
   blaze::CompressedMatrix<float,columnMajor> M5( 3UL, 2UL );

   // ... Initialization of the matrices

   M1 = M2;  // OK: Assignment of a 3x2 dense row-major matrix to another 3x2 dense row-major matrix
   M1 = M4;  // OK: Assignment of a 3x2 sparse row-major matrix to a 3x2 dense row-major matrix
   M1 = M3;  // Runtime error: Cannot assign a 5x2 matrix to a 3x2 static matrix
   M1 = M5;  // OK: Assignment of a 3x2 sparse column-major matrix to a 3x2 dense row-major matrix
   \endcode

// \n \subsection matrix_operations_compound_assignment Compound Assignment
//
// Compound assignment is also available for matrices: addition assignment, subtraction assignment,
// and multiplication assignment. In contrast to plain assignment, however, the number of rows
// and columns of the two operands have to match according to the arithmetic operation.

   \code
   blaze::StaticMatrix<int,2UL,3UL,rowMajor>    M1;
   blaze::DynamicMatrix<int,rowMajor>           M2( 2UL, 3UL );
   blaze::CompressedMatrix<float,columnMajor>   M3( 2UL, 3UL );
   blaze::CompressedMatrix<float,rowMajor>      M4( 2UL, 4UL );
   blaze::StaticMatrix<float,2UL,4UL,rowMajor>  M5;
   blaze::CompressedMatrix<float,rowMajor>      M6( 3UL, 2UL );

   // ... Initialization of the matrices

   M1 += M2;  // OK: Addition assignment between two row-major matrices of the same dimensions
   M1 -= M3;  // OK: Subtraction assignment between between a row-major and a column-major matrix
   M1 += M4;  // Runtime error: No compound assignment between matrices of different size
   M1 -= M5;  // Compilation error: No compound assignment between matrices of different size
   M2 *= M6;  // OK: Multiplication assignment between two row-major matrices
   \endcode

// Note that the multiplication assignment potentially changes the number of columns of the
// target matrix:

                          \f$\left(\begin{array}{*{3}{c}}
                          2 & 0 & 1 \\
                          0 & 3 & 2 \\
                          \end{array}\right) \times
                          \left(\begin{array}{*{2}{c}}
                          4 & 0 \\
                          1 & 0 \\
                          0 & 3 \\
                          \end{array}\right) =
                          \left(\begin{array}{*{2}{c}}
                          8 & 3 \\
                          3 & 6 \\
                          \end{array}\right)\f$

// Since a StaticMatrix cannot change its size, only a quadratic StaticMatrix can be used in a
// multiplication assignment with other quadratic matrices of the same dimensions.
//
//
// \n \section matrix_operations_common_matrix_operations Common Matrix Operations
// <hr>
//
// \subsection matrix_operations_rows Number of Rows of a Matrix
//
// The current number of rows of a matrix can be acquired via the \c rows() function:

   \code
   // Instantiating a dynamic matrix with 10 rows and 8 columns
   blaze::DynamicMatrix<int> M1( 10UL, 8UL );
   M1.rows();  // Returns 10

   // Instantiating a compressed matrix with 5 rows and 12 columns
   blaze::CompressedMatrix<double> M2( 5UL, 12UL );
   M2.rows();  // Returns 5
   \endcode

// \n \subsection matrix_operations_columns Number of Columns of a Matrix
//
// The current number of columns of a matrix can be acquired via the \c columns() function:

   \code
   // Instantiating a dynamic matrix with 6 rows and 8 columns
   blaze::DynamicMatrix<int> M1( 6UL, 8UL );
   M1.columns();  // Returns 8

   // Instantiating a compressed matrix with 4 rows and 7 columns
   blaze::CompressedMatrix<double> M2( 4UL, 7UL );
   M2.columns();  // Returns 7
   \endcode

// \n \subsection matrix_operations_capacity Capacity of a Matrix
//
// The \c capacity() function returns the internal capacity of a DynamicMatrix or CompressedMatrix.
// Note that the capacity of a matrix doesn't have to be equal to the size of a matrix. In case of
// a dense matrix the capacity will always be greater or equal than the total number of elements
// of the matrix. In case of a sparse matrix, the capacity will usually be much less than the
// total number of elements.

   \code
   blaze::DynamicMatrix<float> M1( 5UL, 7UL );
   M1.capacity();  // Returns at least 35
   \endcode

// \n \subsection matrix_operations_nonzeros Number of Non-Zero Elements
//
// For both dense and sparse matrices the current number of non-zero elements can be queried via
// the \c nonZeros() function. In case of matrices there are two flavors of the nonZeros() function:
// One returns the total number of non-zero elements in the matrix, the second returns the number
// of non-zero elements in a specific row (in case of a row-major matrix) or column (in case of
// a column-major matrix). Sparse matrices directly return their number of non-zero elements,
// dense matrices traverse their elements and count the number of non-zero elements.

   \code
   blaze::DynamicMatrix<int,rowMajor> M1( 3UL, 5UL );
   blaze::CompressedMatrix<double,columnMajor> M2( 4UL, 7UL );

   // ... Initializing the matrices

   M1.nonZeros();     // Returns the total number of non-zero elements in the dense matrix
   M1.nonZeros( 2 );  // Returns the number of non-zero elements in row 2

   M2.nonZeros();     // Returns the total number of non-zero elements in the sparse matrix
   M2.nonZeros( 3 );  // Returns the number of non-zero elements in column 3
   \endcode

// \n \section matrix_operations_resize_reserve Resize/Reserve
// <hr>
//
// The dimensions of a StaticMatrix are fixed at compile time by the second and third template
// parameter. In contrast, the number or rows and/or columns of DynamicMatrix and CompressedMatrix
// can be changed at runtime:

   \code
   using blaze::DynamicMatrix;
   using blaze::CompressedMatrix;

   DynamicMatrix<int,rowMajor> M1;
   CompressedMatrix<int,columnMajor> M2( 3UL, 2UL );

   // Adapting the number of rows and columns via the resize() function. The (optional)
   // third parameter specifies whether the existing elements should be preserved.
   M1.resize( 2UL, 2UL );         // Resizing matrix M1 to 2x2 elements. Elements of built-in type
                                  // remain uninitialized, elements of class type are default
                                  // constructed.
   M1.resize( 3UL, 1UL, false );  // Resizing M1 to 3x1 elements. The old elements are lost, the
                                  // new elements are NOT initialized!
   M2.resize( 5UL, 7UL, true );   // Resizing M2 to 5x7 elements. The old elements are preserved.
   M2.resize( 3UL, 2UL, false );  // Resizing M2 to 3x2 elements. The old elements are lost.
   \endcode

// Note that resizing a matrix invalidates all existing views (see e.g. \ref views_submatrices)
// on the matrix:

   \code
   typedef blaze::DynamicMatrix<int,rowMajor>  MatrixType;
   typedef blaze::DenseRow<MatrixType>         RowType;

   MatrixType M1( 10UL, 20UL );    // Creating a 10x20 matrix
   RowType row8 = row( M1, 8UL );  // Creating a view on the 8th row of the matrix
   M1.resize( 6UL, 20UL );         // Resizing the matrix invalidates the view
   \endcode

// When the internal capacity of a matrix is no longer sufficient, the allocation of a larger
// junk of memory is triggered. In order to avoid frequent reallocations, the \c reserve()
// function can be used up front to set the internal capacity:

   \code
   blaze::DynamicMatrix<int> M1;
   M1.reserve( 100 );
   M1.rows();      // Returns 0
   M1.capacity();  // Returns at least 100
   \endcode

// Additionally it is possible to reserve memory in a specific row (for a row-major matrix) or
// column (for a column-major matrix):

   \code
   blaze::CompressedMatrix<int> M1( 4UL, 6UL );
   M1.reserve( 1, 4 );  // Reserving enough space for four non-zero elements in row 4
   \endcode

// \n \section matrix_operations_element_access Element Access
// <hr>
//
// The easiest way to access a specific dense or sparse matrix element is via the function call
// operator. The indices to access a matrix are zero-based:

   \code
   blaze::DynamicMatrix<int> M1( 4UL, 6UL );
   M1(0,0) = 1;
   M1(0,1) = 3;
   // ...

   blaze::CompressedMatrix<double> M2( 5UL, 3UL );
   M2(0,2) =  4.1;
   M2(1,1) = -6.3;
   \endcode

// Since dense matrices allocate enough memory for all contained elements, using the function
// call operator on a dense matrix directly returns a reference to the accessed value. In case
// of a sparse matrix, if the accessed value is currently not contained in the matrix, the
// value is inserted into the matrix prior to returning a reference to the value, which can
// be much more expensive than the direct access to a dense matrix. Consider the following
// example:

   \code
   blaze::CompressedMatrix<int> M1( 4UL, 4UL );

   for( size_t i=0UL; i<M1.rows(); ++i ) {
      for( size_t j=0UL; j<M1.columns(); ++j ) {
         ... = M1(i,j);
      }
   }
   \endcode

// Although the compressed matrix is only used for read access within the for loop, using the
// function call operator temporarily inserts 16 non-zero elements into the matrix. Therefore,
// all matrices (sparse as well as dense) offer an alternate way via the \c begin() and \c end()
// functions to traverse all contained elements by iterator. Note that it is not possible to
// traverse all elements of the matrix, but that it is only possible to traverse elements in a
// row/column-wise fashion. In case of a non-const matrix, \c begin() and \c end() return an
// Iterator, which allows a manipulation of the non-zero value, in case of a constant matrix a
// ConstIterator is returned:

   \code
   using blaze::CompressedMatrix;

   CompressedMatrix<int,rowMajor> M1( 4UL, 6UL );

   // Traversing the matrix by Iterator
   for( size_t i=0UL; i<A.rows(); ++i ) {
      for( CompressedMatrix<int,rowMajor>::Iterator it=A.begin(i); it!=A.end(i); ++it ) {
         it->value() = ...;  // OK: Write access to the value of the non-zero element.
         ... = it->value();  // OK: Read access to the value of the non-zero element.
         it->index() = ...;  // Compilation error: The index of a non-zero element cannot be changed.
         ... = it->index();  // OK: Read access to the index of the non-zero element.
      }
   }

   // Traversing the matrix by ConstIterator
   for( size_t i=0UL; i<A.rows(); ++i ) {
      for( CompressedMatrix<int,rowMajor>::ConstIterator it=A.begin(i); it!=A.end(i); ++it ) {
         it->value() = ...;  // Compilation error: Assignment to the value via a ConstIterator is invalid.
         ... = it->value();  // OK: Read access to the value of the non-zero element.
         it->index() = ...;  // Compilation error: The index of a non-zero element cannot be changed.
         ... = it->index();  // OK: Read access to the index of the non-zero element.
      }
   }
   \endcode

// \n \section matrix_operations_element_insertion Element Insertion
// <hr>
//
// Whereas a dense matrix always provides enough capacity to store all matrix elements, a sparse
// matrix only stores the non-zero elements. Therefore it is necessary to explicitly add elements
// to the matrix. The first possibility to add elements to a sparse matrix is the function call
// operator:

   \code
   using blaze::CompressedMatrix;

   CompressedMatrix<int> M1( 3UL, 4UL );
   M1(1,2) = 9;
   \endcode

// In case the element at the given position is not yet contained in the sparse matrix, it is
// automatically inserted. Otherwise the old value is replaced by the new value 2. The operator
// returns a reference to the sparse vector element.\n
// However, insertion of elements can be better controlled via the \c insert() function. In
// contrast to the function call operator it emits an exception in case the element is already
// contained in the matrix. In order to check for this case, the \c find() function can be used:

   \code
   // In case the element at position (2,3) is not yet contained in the matrix it is inserted
   // with a value of 4.
   if( M1.find( 2, 3 ) == M1.end( 2 ) )
      M1.insert( 2, 3, 4 );
   \endcode

// Although the \c insert() function is very flexible, due to performance reasons it is not
// suited for the setup of large sparse matrices. A very efficient, yet also very low-level
// way to fill a sparse matrix is the \c append() function. It requires the sparse matrix to
// provide enough capacity to insert a new element in the specified row. Additionally, the
// index of the new element must be larger than the index of the previous element in the same
// row. Violating these conditions results in undefined behavior!

   \code
   M1.reserve( 0, 3 );     // Reserving space for three non-zero elements in row 0
   M1.append( 0, 1,  2 );  // Appending the element 2 in row 0 at column index 1
   M1.append( 0, 2, -4 );  // Appending the element -4 in row 0 at column index 2
   // ...
   \endcode

// The most efficient way to fill a sparse matrix with elements, however, is a combination of
// \c reserve(), \c append(), and the \c finalize() function:

   \code
   blaze::CompressedMatrix<int> M1( 3UL, 5UL );
   M1.reserve( 3 );       // Reserving enough space for 3 non-zero elements
   M1.append( 0, 1, 1 );  // Appending the value 1 in row 0 with column index 1
   M1.finalize( 0 );      // Finalizing row 0
   M1.append( 1, 1, 2 );  // Appending the value 2 in row 1 with column index 1
   M1.finalize( 1 );      // Finalizing row 1
   M1.append( 2, 0, 3 );  // Appending the value 3 in row 2 with column index 0
   M1.finalize( 2 );      // Finalizing row 2
   \endcode

// \n \section matrix_operations_reset_clear Reset/Clear
// <hr>
//
// In order to reset all elements of a dense or sparse matrix, the \c reset() function can be
// used. The number of rows and columns of the matrix are preserved:

   \code
   // Setting up a single precision row-major matrix, whose elements are initialized with 2.0F.
   blaze::DynamicMatrix<float> M1( 4UL, 5UL, 2.0F );

   // Resetting all elements to 0.0F.
   reset( M1 );  // Resetting all elements
   M1.rows();    // Returns 4: size and capacity remain unchanged
   \endcode

// In order to return a matrix to its default state (i.e. the state of a default constructed
// matrix), the \c clear() function can be used:

   \code
   // Setting up a single precision row-major matrix, whose elements are initialized with 2.0F.
   blaze::DynamicMatrix<float> M1( 4UL, 5UL, 2.0F );

   // Resetting all elements to 0.0F.
   clear( M1 );  // Resetting the entire matrix
   M1.rows();    // Returns 0: size is reset, but capacity remains unchanged
   \endcode

// \n \section matrix_operations_matrix_transpose Matrix Transpose
// <hr>
//
// Matrices can be transposed via the \c trans() function. Row-major matrices are transposed into
// a column-major matrix and vice versa:

   \code
   blaze::DynamicMatrix<int,rowMajor> M1( 5UL, 2UL );
   blaze::CompressedMatrix<int,columnMajor> M2( 3UL, 7UL );

   M1 = M2;            // Assigning a column-major matrix to a row-major matrix
   M1 = trans( M2 );   // Assigning the transpose of M2 (i.e. a row-major matrix) to M1
   M1 += trans( M2 );  // Addition assignment of two row-major matrices
   \endcode

// \n \section matrix_operators_abs Absolute Values
// <hr>
//
// The \c abs() function can be used to compute the absolute values of each element of a matrix.
// For instance, the following computation

   \code
   blaze::StaticMatrix<int,2UL,3UL,rowMajor> A( -1, 2, -3, 4, -5, 6 );
   blaze::StaticMatrix<int,2UL,3UL,rowMajor> B( abs( A ) );
   \endcode

// results in the matrix

                          \f$ B = \left(\begin{array}{*{3}{c}}
                          1 & 2 & 3 \\
                          4 & 5 & 6 \\
                          \end{array}\right)\f$

// \n \section matrix_operations_swap Swap
//
// Via the \c \c swap() function it is possible to completely swap the contents of two matrices
// of the same type:

   \code
   blaze::DynamicMatrix<int,rowMajor> M1( 10UL, 15UL );
   blaze::DynamicMatrix<int,rowMajor> M2( 20UL, 10UL );

   swap( M1, M2 );  // Swapping the contents of M1 and M2
   \endcode

// \n <center> Previous: \ref matrix_types &nbsp; &nbsp; Next: \ref views_subvectors </center>
*/
//*************************************************************************************************


//**Subvectors*************************************************************************************
/*!\page views_subvectors Subvectors
//
// <center> Previous: \ref matrix_operations &nbsp; &nbsp; Next: \ref views_submatrices </center> \n
//
//
// \tableofcontents
//
//
// Subvectors provide views on a specific part of a dense or sparse vector. As such, subvectors
// act as a reference to a specific range within a vector. This reference is valid and can be
// used in every way any other dense or sparse vector can be used as long as the vector containing
// the subvector is not resized or entirely destroyed. The subvector also acts as an alias to the
// vector elements in the specified range: Changes made to the elements (e.g. modifying values,
// inserting or erasing elements) are immediately visible in the vector and changes made via the
// vector are immediately visible in the subvector. \b Blaze provides two subvector types:
// \ref views_dense_subvector and \ref views_sparse_subvector.
//
//
// \n \section views_dense_subvector DenseSubvector
// <hr>
//
// The blaze::DenseSubvector template represents a view on a specific subvector of a dense vector
// primitive. It can be included via the header file

   \code
   #include <blaze/math/DenseSubvector.h>
   \endcode

// The type of the dense vector is specified two template parameters:

   \code
   template< typename VT, bool AF >
   class DenseSubvector;
   \endcode

//  - \c VT: specifies the type of the dense vector primitive. DenseSubvector can be used with
//           every dense vector primitive or view, but does not work with any vector expression
//           type.
//  - \c AF: the alignment flag specifies whether the subvector is aligned (blaze::aligned) or
//           unaligned (blaze::unaligned). The default value is blaze::unaligned.
//
//
// \n \section views_sparse_subvector SparseSubvector
// <hr>
//
// The blaze::SparseSubvector template represents a view on a specific subvector of a sparse
// vector primitive. It can be included via the header file

   \code
   #include <blaze/math/SparseSubvector.h>
   \endcode

// The type of the sparse vector is specified via two template parameters:

   \code
   template< typename VT, bool AF >
   class SparseSubvector;
   \endcode

//  - \c VT: specifies the type of the sparse vector primitive. As in case of DenseSubvector, a
//           SparseSubvector can be used with every sparse vector primitive or view, but does not
//           work with any vector expression type.
//  - \c AF: the alignment flag specifies whether the subvector is aligned (blaze::aligned) or
//           unaligned (blaze::unaligned). The default value is blaze::unaligned.
//
//
// \n \section views_subvectors_setup Setup of Subvectors
// <hr>
//
// A view on a dense or sparse subvector can be created very conveniently via the \c subvector()
// function. This view can be treated as any other vector, i.e. it can be assigned to, it can
// be copied from, and it can be used in arithmetic operations. A subvector created from a row
// vector can be used as any other row vector, a subvector created from a column vector can be
// used as any other column vector. The view can also be used on both sides of an assignment:
// The subvector can either be used as an alias to grant write access to a specific subvector
// of a dense vector primitive on the left-hand side of an assignment or to grant read-access
// to a specific subvector of a vector primitive or expression on the right-hand side of an
// assignment. The following example demonstrates this in detail:

   \code
   typedef blaze::DynamicVector<double,blaze::rowVector>  DenseVectorType;
   typedef blaze::CompressedVector<int,blaze::rowVector>  SparseVectorType;

   DenseVectorType  d1, d2;
   SparseVectorType s1, s2;
   // ... Resizing and initialization

   // Creating a view on the first ten elements of the dense vector d1
   blaze::DenseSubvector<DenseVectorType> dsv = subvector( d1, 0UL, 10UL );

   // Creating a view on the second ten elements of the sparse vector s1
   blaze::SparseSubvector<SparseVectorType> ssv = subvector( s1, 10UL, 10UL );

   // Creating a view on the addition of d2 and s2
   dsv = subvector( d2 + s2, 5UL, 10UL );

   // Creating a view on the multiplication of d2 and s2
   ssv = subvector( d2 * s2, 2UL, 10UL );
   \endcode

// The \c subvector() function can be used on any dense or sparse vector, including expressions,
// as demonstrated in the example. Note however that a \ref views_dense_subvector or
// \ref views_sparse_subvector can only be instantiated with a dense or sparse vector primitive,
// respectively, i.e. with types that can be written, and not with an expression type.
//
//
// \n \section views_subvectors_common_operations Common Operations
// <hr>
//
// A subvector view can be used like any other dense or sparse vector. For instance, the current
// number of elements can be obtained via the \c size() function, the current capacity via the
// \c capacity() function, and the number of non-zero elements via the \c nonZeros() function.
// However, since subvectors are references to a specific range of a vector, several operations
// are not possible on views, such as resizing and swapping. The following example shows this by
// means of a dense subvector view:

   \code
   typedef blaze::DynamicVector<int,blaze::rowVector>  VectorType;
   typedef blaze::DenseSubvector<VectorType>           SubvectorType;

   VectorType v( 42UL );
   // ... Resizing and initialization

   // Creating a view on the range [5..15] of vector v
   SubvectorType sv = subvector( v, 5UL, 10UL );

   sv.size();          // Returns the number of elements in the subvector
   sv.capacity();      // Returns the capacity of the subvector
   sv.nonZeros();      // Returns the number of non-zero elements contained in the subvector

   sv.resize( 84UL );  // Compilation error: Cannot resize a subvector of a vector

   SubvectorType sv2 = subvector( v, 15UL, 10UL );
   swap( sv, sv2 );   // Compilation error: Swap operation not allowed
   \endcode

// \n \section views_subvectors_element_access Element Access
// <hr>
//
// The elements of a subvector can be directly accessed via the subscript operator:

   \code
   typedef blaze::DynamicVector<double,blaze::rowVector>  VectorType;
   VectorType v;
   // ... Resizing and initialization

   // Creating an 8-dimensional subvector, starting from index 4
   blaze::DenseSubvector<VectorType> sv = subvector( v, 4UL, 8UL );

   // Setting the 1st element of the subvector, which corresponds to
   // the element at index 5 in vector v
   sv[1] = 2.0;
   \endcode

   \code
   typedef blaze::CompressedVector<double,blaze::rowVector>  VectorType;
   VectorType v;
   // ... Resizing and initialization

   // Creating an 8-dimensional subvector, starting from index 4
   blaze::SparseSubvector<VectorType> sv = subvector( v, 4UL, 8UL );

   // Setting the 1st element of the subvector, which corresponds to
   // the element at index 5 in vector v
   sv[1] = 2.0;
   \endcode

// The numbering of the subvector elements is

                             \f[\left(\begin{array}{*{5}{c}}
                             0 & 1 & 2 & \cdots & N-1 \\
                             \end{array}\right),\f]

// where N is the specified size of the subvector. Alternatively, the elements of a subvector can
// be traversed via iterators. Just as with vectors, in case of non-const subvectors, \c begin()
// and \c end() return an Iterator, which allows a manipulation of the non-zero values, in case
// of constant subvectors a ConstIterator is returned:

   \code
   typedef blaze::DynamicVector<int,blaze::rowVector>  VectorType;
   typedef blaze::DenseSubvector<VectorType>           SubvectorType;

   VectorType v( 256UL );
   // ... Resizing and initialization

   // Creating a reference to a specific subvector of the dense vector v
   SubvectorType sv = subvector( v, 16UL, 64UL );

   for( SubvectorType::Iterator it=sv.begin(); it!=sv.end(); ++it ) {
      *it = ...;  // OK: Write access to the dense subvector value.
      ... = *it;  // OK: Read access to the dense subvector value.
   }

   for( SubvectorType::ConstIterator it=sv.begin(); it!=sv.end(); ++it ) {
      *it = ...;  // Compilation error: Assignment to the value via a ConstIterator is invalid.
      ... = *it;  // OK: Read access to the dense subvector value.
   }
   \endcode

   \code
   typedef blaze::CompressedVector<int,blaze::rowVector>  VectorType;
   typedef blaze::SparseSubvector<VectorType>             SubvectorType;

   VectorType v( 256UL );
   // ... Resizing and initialization

   // Creating a reference to a specific subvector of the sparse vector v
   SubvectorType sv = subvector( v, 16UL, 64UL );

   for( SubvectorType::Iterator it=sv.begin(); it!=sv.end(); ++it ) {
      it->value() = ...;  // OK: Write access to the value of the non-zero element.
      ... = it->value();  // OK: Read access to the value of the non-zero element.
      it->index() = ...;  // Compilation error: The index of a non-zero element cannot be changed.
      ... = it->index();  // OK: Read access to the index of the sparse element.
   }

   for( SubvectorType::ConstIterator it=sv.begin(); it!=sv.end(); ++it ) {
      it->value() = ...;  // Compilation error: Assignment to the value via a ConstIterator is invalid.
      ... = it->value();  // OK: Read access to the value of the non-zero element.
      it->index() = ...;  // Compilation error: The index of a non-zero element cannot be changed.
      ... = it->index();  // OK: Read access to the index of the sparse element.
   }
   \endcode

// \n \section views_subvectors_element_insertion Element Insertion
// <hr>
//
// Inserting/accessing elements in a sparse subvector can be done by several alternative functions.
// The following example demonstrates all options:

   \code
   typedef blaze::CompressedVector<double,blaze::rowVector>  VectorType;
   VectorType v( 256UL );  // Non-initialized vector of size 256

   typedef blaze::SparseSubvector<VectorType>  SubvectorType;
   SubvectorType sv( subvector( v, 10UL, 60UL ) );  // View on the range [10..69] of v

   // The subscript operator provides access to all possible elements of the sparse subvector,
   // including the zero elements. In case the subscript operator is used to access an element
   // that is currently not stored in the sparse subvector, the element is inserted into the
   // subvector.
   sv[42] = 2.0;

   // An alternative for inserting elements into the subvector is the insert() function. However,
   // it inserts the element only in case the element is not already contained in the subvector.
   sv.insert( 50UL, 3.7 );

   // Just as in case of vectors, elements can also be inserted via the append() function. In
   // case of subvectors, append() also requires that the appended element's index is strictly
   // larger than the currently largest non-zero index of the subvector and that the subvector's
   // capacity is large enough to hold the new element. Note however that due to the nature of
   // a subvector, which may be an alias to the middle of a sparse vector, the append() function
   // does not work as efficiently for a subvector as it does for a vector.
   sv.reserve( 10UL );
   sv.append( 51UL, -2.1 );
   \endcode

// \n \section views_subvectors_arithmetic_operations Arithmetic Operations
// <hr>
//
// Both dense and sparse subvectors can be used in all arithmetic operations that any other dense
// or sparse vector can be used in. The following example gives an impression of the use of dense
// subvectors within arithmetic operations. All operations (addition, subtraction, multiplication,
// scaling, ...) can be performed on all possible combinations of dense and sparse subvectors with
// fitting element types:

   \code
   typedef blaze::DynamicVector<double,blaze::rowVector>     DenseVectorType;
   typedef blaze::CompressedVector<double,blaze::rowVector>  SparseVectorType;
   DenseVectorType d1, d2, d3;
   SparseVectorType s1, s2;

   // ... Resizing and initialization

   typedef blaze::DynamicMatrix<double,blaze::rowMajor>  DenseMatrixType;
   DenseMatrixType A;

   typedef blaze::DenseSubvector<DenseVectorType>  SubvectorType;
   SubvectorType dsv( subvector( d1, 0UL, 10UL ) );  // View on the range [0..9] of vector d1

   dsv = d2;                          // Dense vector initialization of the range [0..9]
   subvector( d1, 10UL, 10UL ) = s1;  // Sparse vector initialization of the range [10..19]

   d3 = dsv + d2;                           // Dense vector/dense vector addition
   s2 = s1 + subvector( d1, 10UL, 10UL );   // Sparse vector/dense vector addition
   d2 = dsv * subvector( d1, 20UL, 10UL );  // Component-wise vector multiplication

   subvector( d1, 3UL, 4UL ) *= 2.0;      // In-place scaling of the range [3..6]
   d2 = subvector( d1, 7UL, 3UL ) * 2.0;  // Scaling of the range [7..9]
   d2 = 2.0 * subvector( d1, 7UL, 3UL );  // Scaling of the range [7..9]

   subvector( d1, 0UL , 10UL ) += d2;   // Addition assignment
   subvector( d1, 10UL, 10UL ) -= s2;   // Subtraction assignment
   subvector( d1, 20UL, 10UL ) *= dsv;  // Multiplication assignment

   double scalar = subvector( d1, 5UL, 10UL ) * trans( s1 );  // Scalar/dot/inner product between two vectors

   A = trans( s1 ) * subvector( d1, 4UL, 16UL );  // Outer product between two vectors
   \endcode

// \n \section views_aligned_subvectors Aligned Subvectors
// <hr>
//
// Usually subvectors can be defined anywhere within a vector. They may start at any position and
// may have an arbitrary size (only restricted by the size of the underlying vector). However, in
// contrast to vectors themselves, which are always properly aligned in memory and therefore can
// provide maximum performance, this means that subvectors in general have to be considered to be
// unaligned. This can be made explicit by the blaze::unaligned flag:

   \code
   using blaze::unaligned;

   typedef blaze::DynamicVector<double,blaze::rowVector>  DenseVectorType;

   DenseVectorType x;
   // ... Resizing and initialization

   // Identical creations of an unaligned subvector in the range [8..23]
   blaze::DenseSubvector<DenseVectorType>           sv1 = subvector           ( x, 8UL, 16UL );
   blaze::DenseSubvector<DenseVectorType>           sv2 = subvector<unaligned>( x, 8UL, 16UL );
   blaze::DenseSubvector<DenseVectorType,unaligned> sv3 = subvector           ( x, 8UL, 16UL );
   blaze::DenseSubvector<DenseVectorType,unaligned> sv4 = subvector<unaligned>( x, 8UL, 16UL );
   \endcode

// All of these calls to the \c subvector() function are identical. Whether the alignment flag is
// explicitly specified or not, it always returns an unaligned subvector. Whereas this may provide
// full flexibility in the creation of subvectors, this might result in performance disadvantages
// in comparison to vector primitives (even in case the specified subvector could be aligned).
// Whereas vector primitives are guaranteed to be properly aligned and therefore provide maximum
// performance in all operations, a general view on a vector might not be properly aligned. This
// may cause a performance penalty on some platforms and/or for some operations.
//
// However, it is also possible to create aligned subvectors. Aligned subvectors are identical to
// unaligned subvectors in all aspects, except that they may pose additional alignment restrictions
// and therefore have less flexibility during creation, but don't suffer from performance penalties
// and provide the same performance as the underlying vector. Aligned subvectors are created by
// explicitly specifying the blaze::aligned flag:

   \code
   using blaze::aligned;

   // Creating an aligned dense subvector in the range [8..23]
   blaze::DenseSubvector<DenseVectorType,aligned> sv = subvector<aligned>( x, 8UL, 16UL );
   \endcode

// The alignment restrictions refer to system dependent address restrictions for the used element
// type and the available vectorization mode (SSE, AVX, ...). The following source code gives some
// examples for a double precision dense vector, assuming that AVX is available, which packs 4
// \c double values into an intrinsic vector:

   \code
   using blaze::columnVector;

   typedef blaze::DynamicVector<double,columnVector>  VectorType;
   typedef blaze::DenseSubvector<VectorType,aligned>  SubvectorType;

   VectorType d( 17UL );
   // ... Resizing and initialization

   // OK: Starts at the beginning and the size is a multiple of 4
   SubvectorType dsv1 = subvector<aligned>( d, 0UL, 12UL );

   // OK: Start index and the size are both a multiple of 4
   SubvectorType dsv2 = subvector<aligned>( d, 4UL, 8UL );

   // OK: The start index is a multiple of 4 and the subvector includes the last element
   SubvectorType dsv3 = subvector<aligned>( d, 8UL, 9UL );

   // Error: Start index is not a multiple of 4
   SubvectorType dsv4 = subvector<aligned>( d, 5UL, 8UL );

   // Error: Size is not a multiple of 4 and the subvector does not include the last element
   SubvectorType dsv5 = subvector<aligned>( d, 8UL, 5UL );
   \endcode

// Note that the discussed alignment restrictions are only valid for aligned dense subvectors.
// In contrast, aligned sparse subvectors at this time don't pose any additional restrictions.
// Therefore aligned and unaligned sparse subvectors are truly fully identical. Still, in case
// the blaze::aligned flag is specified during setup, an aligned subvector is created:

   \code
   using blaze::aligned;

   typedef blaze::CompressedVector<double,blaze::rowVector>  SparseVectorType;

   SparseVectorType x;
   // ... Resizing and initialization

   // Creating an aligned subvector in the range [8..23]
   blaze::SparseSubvector<SparseVectorType,aligned> sv = subvector<aligned>( x, 8UL, 16UL );
   \endcode

// \n \section views_subvectors_on_subvectors Subvectors on Subvectors
// <hr>
//
// It is also possible to create a subvector view on another subvector. In this context it is
// important to remember that the type returned by the \c subvector() function is the same type
// as the type of the given subvector, not a nested subvector type, since the view on a subvector
// is just another view on the underlying vector:

   \code
   typedef blaze::DynamicVector<double,blaze::rowVector>  VectorType;
   typedef blaze::DenseSubvector<VectorType>              SubvectorType;

   VectorType d1;

   // ... Resizing and initialization

   // Creating a subvector view on the dense vector d1
   SubvectorType sv1 = subvector( d1, 5UL, 10UL );

   // Creating a subvector view on the dense subvector sv1
   SubvectorType sv2 = subvector( sv1, 1UL, 5UL );
   \endcode

// \n <center> Previous: \ref matrix_operations &nbsp; &nbsp; Next: \ref views_submatrices </center>
*/
//*************************************************************************************************


//**Submatrices************************************************************************************
/*!\page views_submatrices Submatrices
//
// <center> Previous: \ref views_subvectors &nbsp; &nbsp; Next: \ref views_rows </center> \n
//
//
// \tableofcontents
//
//
// Submatrices provide views on a specific part of a dense or sparse matrix just as subvectors
// provide views on specific parts of vectors. As such, submatrices act as a reference to a
// specific block within a matrix. This reference is valid and can be used in evary way any
// other dense or sparse matrix can be used as long as the matrix containing the submatrix is
// not resized or entirely destroyed. The submatrix also acts as an alias to the matrix elements
// in the specified block: Changes made to the elements (e.g. modifying values, inserting or
// erasing elements) are immediately visible in the matrix and changes made via the matrix are
// immediately visible in the submatrix. \b Blaze provides two submatrix types:
// \ref views_dense_submatrix and \ref views_sparse_submatrix.
//
//
// \n \section views_dense_submatrix DenseSubmatrix
// <hr>
//
// The blaze::DenseSubmatrix template represents a view on a specific submatrix of a dense matrix
// primitive. It can be included via the header file

   \code
   #include <blaze/math/DenseSubmatrix.h>
   \endcode

// The type of the dense matrix is specified via two template parameters:

   \code
   template< typename MT, bool AF >
   class DenseSubmatrix;
   \endcode

//  - \c MT: specifies the type of the dense matrix primitive. DenseSubmatrix can be used with
//           every dense matrix primitive, but does not work with any matrix expression type.
//  - \c AF: the alignment flag specifies whether the submatrix is aligned (blaze::aligned) or
//        unaligned (blaze::unaligned). The default value is blaze::unaligned.
//
//
// \n \section views_sparse_submatrix SparseSubmatrix
// <hr>
//
// The blaze::SparseSubmatrix template represents a view on a specific submatrix of a sparse
// matrix primitive. It can be included via the header file

   \code
   #include <blaze/math/SparseSubmatrix.h>
   \endcode

// The type of the sparse matrix is specified via two template parameters:

   \code
   template< typename MT, bool AF >
   class SparseSubmatrix;
   \endcode

//  - \c MT: specifies the type of the sparse matrix primitive. SparseSubmatrix can be used with
//           every sparse matrix primitive, but does not work with any matrix expression type.
//  - \c AF: the alignment flag specifies whether the submatrix is aligned (blaze::aligned) or
//           unaligned (blaze::unaligned). The default value is blaze::unaligned.
//
//
// \n \section views_submatrices_setup Setup of Submatrices
// <hr>
//
// A view on a submatrix can be created very conveniently via the \c submatrix() function.
// This view can be treated as any other matrix, i.e. it can be assigned to, it can be copied
// from, and it can be used in arithmetic operations. A submatrix created from a row-major
// matrix will itself be a row-major matrix, a submatrix created from a column-major matrix
// will be a column-major matrix. The view can also be used on both sides of an assignment:
// The submatrix can either be used as an alias to grant write access to a specific submatrix
// of a dense matrix primitive on the left-hand side of an assignment or to grant read-access
// to a specific submatrix of a matrix primitive or expression on the right-hand side of an
// assignment. The following example demonstrates this in detail:

   \code
   typedef blaze::DynamicMatrix<double,blaze::rowMajor>     DenseMatrixType;
   typedef blaze::CompressedVector<int,blaze::columnMajor>  SparseMatrixType;

   DenseMatrixType  D1, D2;
   SparseMatrixType S1, S2;
   // ... Resizing and initialization

   // Creating a view on the first 8x16 block of the dense matrix D1
   blaze::DenseSubmatrix<DenseMatrixType> dsm = submatrix( D1, 0UL, 0UL, 8UL, 16UL );

   // Creating a view on the second 8x16 block of the sparse matrix S1
   blaze::SparseSubmatrix<SparseMatrixType> ssm = submatrix( S1, 0UL, 16UL, 8UL, 16UL );

   // Creating a view on the addition of D2 and S2
   dsm = submatrix( D2 + S2, 5UL, 10UL, 8UL, 16UL );

   // Creating a view on the multiplication of D2 and S2
   ssm = submatrix( D2 * S2, 7UL, 13UL, 8UL, 16UL );
   \endcode
//
//
// \n \section views_submatrices_common_operations Common Operations
// <hr>
//
// The current size of the matrix, i.e. the number of rows or columns can be obtained via the
// \c rows() and \c columns() functions, the current total capacity via the \c capacity() function,
// and the number of non-zero elements via the \c nonZeros() function. However, since submatrices
// are views on a specific submatrix of a matrix, several operations are not possible on views,
// such as resizing and swapping:

   \code
   typedef blaze::DynamicMatrix<int,blaze::rowMajor>  MatrixType;
   typedef blaze::DenseSubmatrix<MatrixType>          SubmatrixType;

   MatrixType A;
   // ... Resizing and initialization

   // Creating a view on the a 8x12 submatrix of matrix A
   SubmatrixType sm = submatrix( A, 0UL, 0UL, 8UL, 12UL );

   sm.rows();      // Returns the number of rows of the submatrix
   sm.columns();   // Returns the number of columns of the submatrix
   sm.capacity();  // Returns the capacity of the submatrix
   sm.nonZeros();  // Returns the number of non-zero elements contained in the submatrix

   sm.resize( 10UL, 8UL );  // Compilation error: Cannot resize a submatrix of a matrix

   SubmatrixType sm2 = submatrix( A, 8UL, 0UL, 12UL, 8UL );
   swap( sm, sm2 );  // Compilation error: Swap operation not allowed
   \endcode

// \n \section views_submatrices_element_access Element Access
// <hr>
//
// The elements of a submatrix can be directly accessed with the function call operator:

   \code
   typedef blaze::DynamicMatrix<double,blaze::rowMajor>  MatrixType;
   MatrixType A;
   // ... Resizing and initialization

   // Creating a 8x8 submatrix, starting from position (4,4)
   blaze::DenseSubmatrix<MatrixType> sm = submatrix( A, 4UL, 4UL, 8UL, 8UL );

   // Setting the element (0,0) of the submatrix, which corresponds to
   // the element at position (4,4) in matrix A
   sm(0,0) = 2.0;
   \endcode

   \code
   typedef blaze::CompressedMatrix<double,blaze::rowMajor>  MatrixType;
   MatrixType A;
   // ... Resizing and initialization

   // Creating a 8x8 submatrix, starting from position (4,4)
   blaze::SparseSubmatrix<MatrixType> sm = submatrix( A, 4UL, 4UL, 8UL, 8UL );

   // Setting the element (0,0) of the submatrix, which corresponds to
   // the element at position (4,4) in matrix A
   sm(0,0) = 2.0;
   \endcode

// Alternatively, the elements of a submatrix can be traversed via (const) iterators. Just as
// with matrices, in case of non-const submatrices, \c begin() and \c end() return an Iterator,
// which allows a manipulation of the non-zero values, in case of constant submatrices a
// ConstIterator is returned:

   \code
   typedef blaze::DynamicMatrix<int,blaze::rowMajor>  MatrixType;
   typedef blaze::DenseSubmatrix<MatrixType>          SubmatrixType;

   MatrixType A( 256UL, 512UL );
   // ... Resizing and initialization

   // Creating a reference to a specific submatrix of the dense matrix A
   SubmatrixType sm = submatrix( A, 16UL, 16UL, 64UL, 128UL );

   // Traversing the elements of the 0th row via iterators to non-const elements
   for( SubmatrixType::Iterator it=sm.begin(0); it!=sm.end(0); ++it ) {
      *it = ...;  // OK: Write access to the dense submatrix value.
      ... = *it;  // OK: Read access to the dense submatrix value.
   }

   // Traversing the elements of the 1st row via iterators to const elements
   for( SubmatrixType::ConstIterator it=sm.begin(1); it!=sm.end(1); ++it ) {
      *it = ...;  // Compilation error: Assignment to the value via a ConstIterator is invalid.
      ... = *it;  // OK: Read access to the dense submatrix value.
   }
   \endcode

   \code
   typedef blaze::CompressedMatrix<int,blaze::rowMajor>  MatrixType;
   typedef blaze::SparseSubmatrix<MatrixType>            SubmatrixType;

   MatrixType A( 256UL, 512UL );
   // ... Resizing and initialization

   // Creating a reference to a specific submatrix of the sparse matrix A
   SubmatrixType sm = submatrix( A, 16UL, 16UL, 64UL, 128UL );

   // Traversing the elements of the 0th row via iterators to non-const elements
   for( SubmatrixType::Iterator it=sm.begin(0); it!=sm.end(0); ++it ) {
      it->value() = ...;  // OK: Write access to the value of the non-zero element.
      ... = it->value();  // OK: Read access to the value of the non-zero element.
      it->index() = ...;  // Compilation error: The index of a non-zero element cannot be changed.
      ... = it->index();  // OK: Read access to the index of the sparse element.
   }

   // Traversing the elements of the 1st row via iterators to const elements
   for( SubmatrixType::ConstIterator it=sm.begin(1); it!=sm.end(1); ++it ) {
      it->value() = ...;  // Compilation error: Assignment to the value via a ConstIterator is invalid.
      ... = it->value();  // OK: Read access to the value of the non-zero element.
      it->index() = ...;  // Compilation error: The index of a non-zero element cannot be changed.
      ... = it->index();  // OK: Read access to the index of the sparse element.
   }
   \endcode

// \n \section views_submatrices_element_insertion Element Insertion
// <hr>
//
// Inserting/accessing elements in a sparse submatrix can be done by several alternative functions.
// The following example demonstrates all options:

   \code
   typedef blaze::CompressedMatrix<double,blaze::rowMajor>  MatrixType;
   MatrixType A( 256UL, 512UL );  // Non-initialized matrix of size 256x512

   typedef blaze::SparseSubmatrix<MatrixType>  SubmatrixType;
   SubmatrixType sm = submatrix( A, 10UL, 10UL, 16UL, 16UL );  // View on a 16x16 submatrix of A

   // The function call operator provides access to all possible elements of the sparse submatrix,
   // including the zero elements. In case the subscript operator is used to access an element
   // that is currently not stored in the sparse submatrix, the element is inserted into the
   // submatrix.
   sm(2,4) = 2.0;

   // An alternative for inserting elements into the submatrix is the insert() function. However,
   // it inserts the element only in case the element is not already contained in the submatrix.
   sm.insert( 2UL, 6UL, 3.7 );

   // Just as in case of sparse matrices, elements can also be inserted via the append() function.
   // In case of submatrices, append() also requires that the appended element's index is strictly
   // larger than the currently largest non-zero index in the according row or column of the
   // submatrix and that the according row's or column's capacity is large enough to hold the new
   // element. Note however that due to the nature of a submatrix, which may be an alias to the
   // middle of a sparse matrix, the append() function does not work as efficiently for a
   // submatrix as it does for a matrix.
   sm.reserve( 2UL, 10UL );
   sm.append( 2UL, 10UL, -2.1 );
   \endcode

// \n \section views_submatrices_arithmetic_operations Arithmetic Operations
// <hr>
//
// Both dense and sparse submatrices can be used in all arithmetic operations that any other dense
// or sparse matrix can be used in. The following example gives an impression of the use of dense
// submatrices within arithmetic operations. All operations (addition, subtraction, multiplication,
// scaling, ...) can be performed on all possible combinations of dense and sparse matrices with
// fitting element types:

   \code
   typedef blaze::DynamicMatrix<double,blaze::rowMajor>     DenseMatrixType;
   typedef blaze::CompressedMatrix<double,blaze::rowMajor>  SparseMatrixType;
   DenseMatrixType D1, D2, D3;
   SparseMatrixType S1, S2;

   typedef blaze::CompressedVector<double,blaze::columnVector>  SparseVectorType;
   SparseVectorType a, b;

   // ... Resizing and initialization

   typedef DenseSubmatrix<DenseMatrixType>  SubmatrixType;
   SubmatrixType sm = submatrix( D1, 0UL, 0UL, 8UL, 8UL );  // View on the 8x8 submatrix of matrix D1
                                                            // starting from row 0 and column 0

   submatrix( D1, 0UL, 8UL, 8UL, 8UL ) = D2;  // Dense matrix initialization of the 8x8 submatrix
                                              // starting in row 0 and column 8
   sm = S1;                                   // Sparse matrix initialization of the second 8x8 submatrix

   D3 = sm + D2;                                    // Dense matrix/dense matrix addition
   S2 = S1  - submatrix( D1, 8UL, 0UL, 8UL, 8UL );  // Sparse matrix/dense matrix subtraction
   D2 = sm * submatrix( D1, 8UL, 8UL, 8UL, 8UL );   // Dense matrix/dense matrix multiplication

   submatrix( D1, 8UL, 0UL, 8UL, 8UL ) *= 2.0;      // In-place scaling of a submatrix of D1
   D2 = submatrix( D1, 8UL, 8UL, 8UL, 8UL ) * 2.0;  // Scaling of the a submatrix of D1
   D2 = 2.0 * sm;                                   // Scaling of the a submatrix of D1

   submatrix( D1, 0UL, 8UL, 8UL, 8UL ) += D2;  // Addition assignment
   submatrix( D1, 8UL, 0UL, 8UL, 8UL ) -= S1;  // Subtraction assignment
   submatrix( D1, 8UL, 8UL, 8UL, 8UL ) *= sm;  // Multiplication assignment

   a = submatrix( D1, 4UL, 4UL, 8UL, 8UL ) * b;  // Dense matrix/sparse vector multiplication
   \endcode

// \n \section views_aligned_submatrices Aligned Submatrices
// <hr>
//
// Usually submatrices can be defined anywhere within a matrix. They may start at any position and
// may have an arbitrary extension (only restricted by the extension of the underlying matrix).
// However, in contrast to matrices themselves, which are always properly aligned in memory and
// therefore can provide maximum performance, this means that submatrices in general have to be
// considered to be unaligned. This can be made explicit by the blaze::unaligned flag:

   \code
   using blaze::unaligned;

   typedef blaze::DynamicMatrix<double,blaze::rowMajor>  DenseMatrixType;

   DenseMatrixType A;
   // ... Resizing and initialization

   // Identical creations of an unaligned submatrix of size 8x8, starting in row 0 and column 0
   blaze::DenseSubmatrix<DenseMatrixType>           sm1 = submatrix           ( A, 0UL, 0UL, 8UL, 8UL );
   blaze::DenseSubmatrix<DenseMatrixType>           sm2 = submatrix<unaligned>( A, 0UL, 0UL, 8UL, 8UL );
   blaze::DenseSubmatrix<DenseMatrixType,unaligned> sm3 = submatrix           ( A, 0UL, 0UL, 8UL, 8UL );
   blaze::DenseSubmatrix<DenseMatrixType,unaligned> sm4 = submatrix<unaligned>( A, 0UL, 0UL, 8UL, 8UL );
   \endcode

// All of these calls to the \c submatrix() function are identical. Whether the alignment flag is
// explicitly specified or not, it always returns an unaligned submatrix. Whereas this may provide
// full flexibility in the creation of submatrices, this might result in performance disadvantages
// in comparison to matrix primitives (even in case the specified submatrix could be aligned).
// Whereas matrix primitives are guaranteed to be properly aligned and therefore provide maximum
// performance in all operations, a general view on a matrix might not be properly aligned. This
// may cause a performance penalty on some platforms and/or for some operations.
//
// However, it is also possible to create aligned submatrices. Aligned submatrices are identical to
// unaligned submatrices in all aspects, except that they may pose additional alignment restrictions
// and therefore have less flexibility during creation, but don't suffer from performance penalties
// and provide the same performance as the underlying matrix. Aligned submatrices are created by
// explicitly specifying the blaze::aligned flag:

   \code
   using blaze::aligned;

   // Creating an aligned submatrix of size 8x8, starting in row 0 and column 0
   blaze::DenseSubmatrix<DenseMatrixType,aligned> sv = submatrix<aligned>( A, 0UL, 0UL, 8UL, 8UL );
   \endcode

// The alignment restrictions refer to system dependent address restrictions for the used element
// type and the available vectorization mode (SSE, AVX, ...). The following source code gives some
// examples for a double precision dense matrix, assuming that AVX is available, which packs 4
// \c double values into an intrinsic vector:

   \code
   using blaze::rowMajor;

   typedef blaze::DynamicMatrix<double,rowMajor>      MatrixType;
   typedef blaze::DenseSubmatrix<MatrixType,aligned>  SubmatrixType;

   MatrixType D( 13UL, 17UL );
   // ... Resizing and initialization

   // OK: Starts at position (0,0) and the number of rows and columns are a multiple of 4
   SubmatrixType dsm1 = submatrix<aligned>( D, 0UL, 0UL, 8UL, 12UL );

   // OK: First row and column and the number of rows and columns are all a multiple of 4
   SubmatrixType dsm2 = submatrix<aligned>( D, 4UL, 12UL, 8UL, 16UL );

   // OK: First row and column are a multiple of 4 and the submatrix includes the last row and column
   SubmatrixType dsm3 = submatrix<aligned>( D, 4UL, 0UL, 9UL, 17UL );

   // Error: First row is not a multiple of 4
   SubmatrixType dsm4 = submatrix<aligned>( D, 2UL, 4UL, 12UL, 12UL );

   // Error: First column is not a multiple of 4
   SubmatrixType dsm5 = submatrix<aligned>( D, 0UL, 2UL, 8UL, 8UL );

   // Error: The number of rows is not a multiple of 4 and the submatrix does not include the last row
   SubmatrixType dsm6 = submatrix<aligned>( D, 0UL, 0UL, 7UL, 8UL );

   // Error: The number of columns is not a multiple of 4 and the submatrix does not include the last column
   SubmatrixType dsm6 = submatrix<aligned>( D, 0UL, 0UL, 8UL, 11UL );
   \endcode

// Note that the discussed alignment restrictions are only valid for aligned dense submatrices.
// In contrast, aligned sparse submatrices at this time don't pose any additional restrictions.
// Therefore aligned and unaligned sparse submatrices are truly fully identical. Still, in case
// the blaze::aligned flag is specified during setup, an aligned submatrix is created:

   \code
   using blaze::aligned;

   typedef blaze::CompressedMatrix<double,blaze::rowMajor>  SparseMatrixType;

   SparseMatrixType A;
   // ... Resizing and initialization

   // Creating an aligned submatrix of size 8x8, starting in row 0 and column 0
   blaze::SparseSubmatrix<SparseMatrixType,aligned> sv = submatrix<aligned>( A, 0UL, 0UL, 8UL, 8UL );
   \endcode

// \n \section views_submatrices_on_submatrices Submatrices on Submatrices
// <hr>
//
// It is also possible to create a submatrix view on another submatrix. In this context it is
// important to remember that the type returned by the \c submatrix() function is the same type
// as the type of the given submatrix, since the view on a submatrix is just another view on the
// underlying matrix:

   \code
   typedef blaze::DynamicMatrix<double,blaze::rowMajor>  MatrixType;
   typedef blaze::DenseSubmatrix<MatrixType>             SubmatrixType;

   MatrixType D1;

   // ... Resizing and initialization

   // Creating a submatrix view on the dense matrix D1
   SubmatrixType sm1 = submatrix( D1, 4UL, 4UL, 8UL, 16UL );

   // Creating a submatrix view on the dense submatrix sm1
   SubmatrixType sm2 = submatrix( sm1, 1UL, 1UL, 4UL, 8UL );
   \endcode

// \n <center> Previous: \ref views_subvectors &nbsp; &nbsp; Next: \ref views_rows </center>
*/
//*************************************************************************************************


//**Rows*******************************************************************************************
/*!\page views_rows Rows
//
// <center> Previous: \ref views_submatrices &nbsp; &nbsp; Next: \ref views_columns </center> \n
//
//
// \tableofcontents
//
//
// Rows provide views on a specific row of a dense or sparse matrix. As such, rows act as a
// reference to a specific row. This reference is valid and can be used in every way any other
// row vector can be used as long as the matrix containing the row is not resized or entirely
// destroyed. The row also acts as an alias to the row elements: Changes made to the elements
// (e.g. modifying values, inserting or erasing elements) are immediately visible in the matrix
// and changes made via the matrix are immediately visible in the row. \b Blaze provides two
// row types: \ref views_dense_row and \ref views_sparse_row.
//
//
// \n \section views_dense_row DenseRow
// <hr>
//
// The blaze::DenseRow class template represents a reference to a specific row of a dense matrix
// primitive. It can be included via the header file

   \code
   #include <blaze/math/DenseRow.h>
   \endcode

// The type of the dense matrix is specified via template parameter:

   \code
   template< typename MT >
   class DenseRow;
   \endcode

// \c MT specifies the type of the dense matrix primitive. DenseRow can be used with every dense
// matrix primitive, but does not work with any matrix expression type.
//
//
// \n \section views_sparse_row SparseRow
// <hr>
//
// The blaze::SparseRow class template represents a reference to a specific row of a sparse matrix
// primitive. It can be included via the header file

   \code
   #include <blaze/math/SparseRow.h>
   \endcode

// The type of the sparse matrix is specified via template parameter:

   \code
   template< typename MT >
   class SparseRow;
   \endcode

// \c MT specifies the type of the sparse matrix primitive. SparseRow can be used with every
// sparse matrix primitive, but does not work with any matrix expression type.
//
//
// \n \section views_rows_setup Setup of Rows
// <hr>
//
// A reference to a dense or sparse row can be created very conveniently via the \c row() function.
// This reference can be treated as any other row vector, i.e. it can be assigned to, it can be
// copied from, and it can be used in arithmetic operations. The reference can also be used on
// both sides of an assignment: The row can either be used as an alias to grant write access to a
// specific row of a matrix primitive on the left-hand side of an assignment or to grant read-access
// to a specific row of a matrix primitive or expression on the right-hand side of an assignment.
// The following two examples demonstrate this for dense and sparse matrices:

   \code
   typedef blaze::DynamicVector<double,rowVector>     DenseVectorType;
   typedef blaze::CompressedVector<double,rowVector>  SparseVectorType;
   typedef blaze::DynamicMatrix<double,rowMajor>      DenseMatrixType;
   typedef blaze::CompressedMatrix<double,rowMajor>   SparseMatrixType;

   DenseVectorType  x;
   SparseVectorType y;
   DenseMatrixType  A, B;
   SparseMatrixType C, D;
   // ... Resizing and initialization

   // Setting the 2nd row of matrix A to x
   blaze::DenseRow<DenseMatrixType> row2 = row( A, 2UL );
   row2 = x;

   // Setting the 3rd row of matrix B to y
   row( B, 3UL ) = y;

   // Setting x to the 4th row of the result of the matrix multiplication
   x = row( A * B, 4UL );

   // Setting y to the 2nd row of the result of the sparse matrix multiplication
   y = row( C * D, 2UL );
   \endcode

// The \c row() function can be used on any dense or sparse matrix, including expressions, as
// illustrated by the source code example. However, both \ref views_dense_row and
// \ref views_sparse_row cannot be instantiated for expression types, but only for dense and
// sparse matrix primitives, respectively, i.e. for matrix types that offer write access.
//
//
// \n \section views_rows_common_operations Common Operations
// <hr>
//
// A row view can be used like any other row vector. For instance, the current number of elements
// can be obtained via the \c size() function, the current capacity via the \c capacity() function,
// and the number of non-zero elements via the \c nonZeros() function. However, since rows are
// references to specific rows of a matrix, several operations are not possible on views, such
// as resizing and swapping. The following example shows this by means of a dense row view:

   \code
   typedef blaze::DynamicMatrix<int,rowMajor>  MatrixType;
   typedef blaze::DenseRow<MatrixType>         RowType;

   MatrixType A( 42UL, 42UL );
   // ... Resizing and initialization

   // Creating a reference to the 2nd row of matrix A
   RowType row2 = row( A, 2UL );

   row2.size();          // Returns the number of elements in the row
   row2.capacity();      // Returns the capacity of the row
   row2.nonZeros();      // Returns the number of non-zero elements contained in the row

   row2.resize( 84UL );  // Compilation error: Cannot resize a single row of a matrix

   RowType row3 = row( A, 3UL );
   swap( row2, row3 );   // Compilation error: Swap operation not allowed
   \endcode

// \n \section views_rows_element_access Element Access
// <hr>
//
// The elements of the row can be directly accessed with the subscript operator. The numbering
// of the row elements is

                             \f[\left(\begin{array}{*{5}{c}}
                             0 & 1 & 2 & \cdots & N-1 \\
                             \end{array}\right),\f]

// where N is the number of columns of the referenced matrix. Alternatively, the elements of
// a row can be traversed via iterators. Just as with vectors, in case of non-const rows,
// \c begin() and \c end() return an Iterator, which allows a manipulation of the non-zero
// value, in case of a constant row a ConstIterator is returned:

   \code
   typedef blaze::DynamicMatrix<int,rowMajor>  MatrixType;
   typedef blaze::DenseRow<MatrixType>         RowType;

   MatrixType A( 128UL, 256UL );
   // ... Resizing and initialization

   // Creating a reference to the 31st row of matrix A
   RowType row31 = row( A, 31UL );

   for( RowType::Iterator it=row31.begin(); it!=row31.end(); ++it ) {
      *it = ...;  // OK; Write access to the dense row value
      ... = *it;  // OK: Read access to the dense row value.
   }

   for( RowType::ConstIterator it=row31.begin(); it!=row31.end(); ++it ) {
      *it = ...;  // Compilation error: Assignment to the value via a ConstIterator is invalid.
      ... = *it;  // OK: Read access to the dense row value.
   }
   \endcode

   \code
   typedef blaze::CompressedMatrix<int,rowMajor>  MatrixType;
   typedef blaze::SparseRow<MatrixType>           RowType;

   MatrixType A( 128UL, 256UL );
   // ... Resizing and initialization

   // Creating a reference to the 31st row of matrix A
   RowType row31 = row( A, 31UL );

   for( RowType::Iterator it=row31.begin(); it!=row31.end(); ++it ) {
      it->value() = ...;  // OK: Write access to the value of the non-zero element.
      ... = it->value();  // OK: Read access to the value of the non-zero element.
      it->index() = ...;  // Compilation error: The index of a non-zero element cannot be changed.
      ... = it->index();  // OK: Read access to the index of the sparse element.
   }

   for( RowType::Iterator it=row31.begin(); it!=row31.end(); ++it ) {
      it->value() = ...;  // Compilation error: Assignment to the value via a ConstIterator is invalid.
      ... = it->value();  // OK: Read access to the value of the non-zero element.
      it->index() = ...;  // Compilation error: The index of a non-zero element cannot be changed.
      ... = it->index();  // OK: Read access to the index of the sparse element.
   }
   \endcode

// \n \section views_rows_element_insertion Element Insertion
// <hr>
//
// Inserting/accessing elements in a sparse row can be done by several alternative functions.
// The following example demonstrates all options:

   \code
   typedef blaze::CompressedMatrix<double,blaze::rowMajor>  MatrixType;
   MatrixType A( 10UL, 100UL );  // Non-initialized 10x100 matrix

   typedef blaze::SparseRow<MatrixType>  RowType;
   RowType row0( row( A, 0UL ) );  // Reference to the 0th row of A

   // The subscript operator provides access to all possible elements of the sparse row,
   // including the zero elements. In case the subscript operator is used to access an element
   // that is currently not stored in the sparse row, the element is inserted into the row.
   row0[42] = 2.0;

   // An alternative for inserting elements into the row is the insert() function. However,
   // it inserts the element only in case the element is not already contained in the row.
   row0.insert( 50UL, 3.7 );

   // A very efficient way to add new elements to a sparse row is the append() function.
   // Note that append() requires that the appended element's index is strictly larger than
   // the currently largest non-zero index of the row and that the row's capacity is large
   // enough to hold the new element.
   row0.reserve( 10UL );
   row0.append( 51UL, -2.1 );
   \endcode

// \n \section views_rows_arithmetic_operations Arithmetic Operations
// <hr>
//
// Both dense and sparse rows can be used in all arithmetic operations that any other dense or
// sparse row vector can be used in. The following example gives an impression of the use of
// dense rows within arithmetic operations. All operations (addition, subtraction, multiplication,
// scaling, ...) can be performed on all possible combinations of dense and sparse rows with
// fitting element types:

   \code
   blaze::DynamicVector<double,blaze::rowVector> a( 2UL, 2.0 ), b;
   blaze::CompressedVector<double,blaze::rowVector> c( 2UL );
   c[1] = 3.0;

   typedef blaze::DynamicMatrix<double,blaze::rowMajor>  DenseMatrix;
   DenseMatrix A( 4UL, 2UL );  // Non-initialized 4x2 matrix

   typedef blaze::DenseRow<DenseMatrix>  RowType;
   RowType row0( row( A, 0UL ) );  // Reference to the 0th row of A

   row0[0] = 0.0;        // Manual initialization of the 0th row of A
   row0[1] = 0.0;
   row( A, 1UL ) = 1.0;  // Homogeneous initialization of the 1st row of A
   row( A, 2UL ) = a;    // Dense vector initialization of the 2nd row of A
   row( A, 3UL ) = c;    // Sparse vector initialization of the 3rd row of A

   b = row0 + a;              // Dense vector/dense vector addition
   b = c + row( A, 1UL );     // Sparse vector/dense vector addition
   b = row0 * row( A, 2UL );  // Component-wise vector multiplication

   row( A, 1UL ) *= 2.0;     // In-place scaling of the 1st row
   b = row( A, 1UL ) * 2.0;  // Scaling of the 1st row
   b = 2.0 * row( A, 1UL );  // Scaling of the 1st row

   row( A, 2UL ) += a;              // Addition assignment
   row( A, 2UL ) -= c;              // Subtraction assignment
   row( A, 2UL ) *= row( A, 0UL );  // Multiplication assignment

   double scalar = row( A, 1UL ) * trans( c );  // Scalar/dot/inner product between two vectors

   A = trans( c ) * row( A, 1UL );  // Outer product between two vectors
   \endcode

// \n \section views_rows_non_fitting_storage_order Views on Matrices with Non-Fitting Storage Order
// <hr>
//
// Especially noteworthy is that row views can be created for both row-major and column-major
// matrices. Whereas the interface of a row-major matrix only allows to traverse a row directly
// and the interface of a column-major matrix only allows to traverse a column, via views it is
// possible to traverse a row of a column-major matrix or a column of a row-major matrix. For
// instance:

   \code
   typedef blaze::CompressedMatrix<int,columnMajor>  MatrixType;
   typedef blaze::SparseRow<MatrixType>              RowType;

   MatrixType A( 64UL, 32UL );
   // ... Resizing and initialization

   // Creating a reference to the 31st row of a column-major matrix A
   RowType row1 = row( A, 1UL );

   for( RowType::Iterator it=row1.begin(); it!=row1.end(); ++it ) {
      // ...
   }
   \endcode

// However, please note that creating a row view on a matrix stored in a column-major fashion
// can result in a considerable performance decrease in comparison to a view on a matrix with
// a fitting storage orientation. This is due to the non-contiguous storage of the matrix
// elements. Therefore care has to be taken in the choice of the most suitable storage order:

   \code
   // Setup of two column-major matrices
   CompressedMatrix<double,columnMajor> A( 128UL, 128UL );
   CompressedMatrix<double,columnMajor> B( 128UL, 128UL );
   // ... Resizing and initialization

   // The computation of the 15th row of the multiplication between A and B ...
   CompressedVector<double,rowVector> x = row( A * B, 15UL );

   // ... is essentially the same as the following computation, which multiplies
   // the 15th row of the column-major matrix A with B.
   CompressedVector<double,rowVector> x = row( A, 15UL ) * B;
   \endcode

// Although \b Blaze performs the resulting vector/matrix multiplication as efficiently as possible
// using a row-major storage order for matrix A would result in a more efficient evaluation.
//
// \n <center> Previous: \ref views_submatrices &nbsp; &nbsp; Next: \ref views_columns </center>
*/
//*************************************************************************************************


//**Columns****************************************************************************************
/*!\page views_columns Columns
//
// <center> Previous: \ref views_rows &nbsp; &nbsp; Next: \ref addition </center> \n
//
//
// \tableofcontents
//
//
// Just as rows provide a view on a specific row of a matrix, columns provide views on a specific
// column of a dense or sparse matrix. As such, columns act as a reference to a specific column.
// This reference is valid an can be used in every way any other column vector can be used as long
// as the matrix containing the column is not resized or entirely destroyed. Changes made to the
// elements (e.g. modifying values, inserting or erasing elements) are immediately visible in the
// matrix and changes made via the matrix are immediately visible in the column. \b Blaze provides
// two column types: \ref views_dense_column and \ref views_sparse_column.
//
//
// \n \section views_dense_column DenseColumn
// <hr>
//
// The blaze::DenseColumn class template represents a reference to a specific column of a dense
// matrix primitive. It can be included via the header file

   \code
   #include <blaze/math/DenseColumn.h>
   \endcode

// The type of the dense matrix is specified via template parameter:

   \code
   template< typename MT >
   class DenseColumn;
   \endcode

// \c MT specifies the type of the dense matrix primitive. DenseColumn can be used with every
// dense matrix primitive, but does not work with any matrix expression type.
//
//
// \n \section views_sparse_column SparseColumn
// <hr>
//
// The blaze::SparseColumn class template represents a reference to a specific column of a sparse
// matrix primitive. It can be included via the header file

   \code
   #include <blaze/math/SparseColumn.h>
   \endcode

// The type of the sparse matrix is specified via template parameter:

   \code
   template< typename MT >
   class SparseColumn;
   \endcode

// \c MT specifies the type of the sparse matrix primitive. SparseColumn can be used with every
// sparse matrix primitive, but does not work with any matrix expression type.
//
//
// \n \section views_colums_setup Setup of Columns
// <hr>
//
// Similar to the setup of a row, a reference to a dense or sparse column can be created very
// conveniently via the \c column() function. This reference can be treated as any other column
// vector, i.e. it can be assigned to, copied from, and be used in arithmetic operations. The
// column can either be used as an alias to grant write access to a specific column of a matrix
// primitive on the left-hand side of an assignment or to grant read-access to a specific column
// of a matrix primitive or expression on the right-hand side of an assignment. The following
// two examples demonstrate this for dense and sparse matrices:

   \code
   typedef blaze::DynamicVector<double,columnVector>     DenseVectorType;
   typedef blaze::CompressedVector<double,columnVector>  SparseVectorType;
   typedef blaze::DynamicMatrix<double,columnMajor>      DenseMatrixType;
   typedef blaze::CompressedMatrix<double,columnMajor>   SparseMatrixType;

   DenseVectorType  x;
   SparseVectorType y;
   DenseMatrixType  A, B;
   SparseMatrixType C, D;
   // ... Resizing and initialization

   // Setting the 1st column of matrix A to x
   blaze::DenseColumn<DenseMatrixType> col1 = column( A, 1UL );
   col1 = x;

   // Setting the 4th column of matrix B to y
   column( B, 4UL ) = y;

   // Setting x to the 2nd column of the result of the matrix multiplication
   x = column( A * B, 2UL );

   // Setting y to the 2nd column of the result of the sparse matrix multiplication
   y = column( C * D, 2UL );
   \endcode

// The \c column() function can be used on any dense or sparse matrix, including expressions,
// as illustrated by the source code example. However, both \ref views_dense_column and
// \ref views_sparse_column cannot be instantiated for expression types, but only for dense
// and sparse matrix primitives, respectively, i.e. for matrix types that offer write access.
//
//
// \n \section views_columns_common_operations Common Operations
// <hr>
//
// A column view can be used like any other column vector. For instance, the current number of
// elements can be obtained via the \c size() function, the current capacity via the \c capacity()
// function, and the number of non-zero elements via the \c nonZeros() function. However, since
// columns are references to specific columns of a matrix, several operations are not possible on
// views, such as resizing and swapping. The following example shows this by means of a dense
// column view:

   \code
   typedef blaze::DynamicMatrix<int,columnMajor>  MatrixType;
   typedef blaze::DenseColumn<MatrixType>         ColumnType;

   MatrixType A( 42UL, 42UL );
   // ... Resizing and initialization

   // Creating a reference to the 2nd column of matrix A
   ColumnType col2 = column( A, 2UL );

   col2.size();          // Returns the number of elements in the column
   col2.capacity();      // Returns the capacity of the column
   col2.nonZeros();      // Returns the number of non-zero elements contained in the column

   col2.resize( 84UL );  // Compilation error: Cannot resize a single column of a matrix

   ColumnType col3 = column( A, 3UL );
   swap( col2, col3 );   // Compilation error: Swap operation not allowed
   \endcode

// \n \section views_columns_element_access Element Access
// <hr>
//
// The elements of the column can be directly accessed with the subscript operator. The numbering
// of the column elements is

                             \f[\left(\begin{array}{*{5}{c}}
                             0 & 1 & 2 & \cdots & N-1 \\
                             \end{array}\right),\f]

// where N is the number of rows of the referenced matrix. Alternatively, the elements of
// a column can be traversed via iterators. Just as with vectors, in case of non-const columns,
// \c begin() and \c end() return an Iterator, which allows a manipulation of the non-zero
// value, in case of a constant column a ConstIterator is returned:

   \code
   typedef blaze::DynamicMatrix<int,columnMajor>  MatrixType;
   typedef blaze::DenseColumn<MatrixType>         ColumnType;

   MatrixType A( 128UL, 256UL );
   // ... Resizing and initialization

   // Creating a reference to the 31st column of matrix A
   ColumnType col31 = column( A, 31UL );

   for( ColumnType::Iterator it=col31.begin(); it!=col31.end(); ++it ) {
      *it = ...;  // OK; Write access to the dense column value
      ... = *it;  // OK: Read access to the dense column value.
   }

   for( ColumnType::ConstIterator it=col31.begin(); it!=col31.end(); ++it ) {
      *it = ...;  // Compilation error: Assignment to the value via a ConstIterator is invalid.
      ... = *it;  // OK: Read access to the dense column value.
   }
   \endcode

   \code
   typedef blaze::CompressedMatrix<int,columnMajor>  MatrixType;
   typedef blaze::SparseColumn<MatrixType>           ColumnType;

   MatrixType A( 128UL, 256UL );
   // ... Resizing and initialization

   // Creating a reference to the 31st column of matrix A
   ColumnType col31 = column( A, 31UL );

   for( ColumnType::Iterator it=col31.begin(); it!=col31.end(); ++it ) {
      it->value() = ...;  // OK: Write access to the value of the non-zero element.
      ... = it->value();  // OK: Read access to the value of the non-zero element.
      it->index() = ...;  // Compilation error: The index of a non-zero element cannot be changed.
      ... = it->index();  // OK: Read access to the index of the sparse element.
   }

   for( ColumnType::Iterator it=col31.begin(); it!=col31.end(); ++it ) {
      it->value() = ...;  // Compilation error: Assignment to the value via a ConstIterator is invalid.
      ... = it->value();  // OK: Read access to the value of the non-zero element.
      it->index() = ...;  // Compilation error: The index of a non-zero element cannot be changed.
      ... = it->index();  // OK: Read access to the index of the sparse element.
   }
   \endcode

// \n \section views_columns_element_insertion Element Insertion
// <hr>
//
// Inserting/accessing elements in a sparse column can be done by several alternative functions.
// The following example demonstrates all options:

   \code
   typedef blaze::CompressedMatrix<double,blaze::columnMajor>  MatrixType;
   MatrixType A( 100UL, 10UL );  // Non-initialized 10x100 matrix

   typedef blaze::SparseColumn<MatrixType>  ColumnType;
   ColumnType col0( column( A, 0UL ) );  // Reference to the 0th column of A

   // The subscript operator provides access to all possible elements of the sparse column,
   // including the zero elements. In case the subscript operator is used to access an element
   // that is currently not stored in the sparse column, the element is inserted into the column.
   col0[42] = 2.0;

   // An alternative for inserting elements into the column is the insert() function. However,
   // it inserts the element only in case the element is not already contained in the column.
   col0.insert( 50UL, 3.7 );

   // A very efficient way to add new elements to a sparse column is the append() function.
   // Note that append() requires that the appended element's index is strictly larger than
   // the currently largest non-zero index of the column and that the column's capacity is
   // large enough to hold the new element.
   col0.reserve( 10UL );
   col0.append( 51UL, -2.1 );
   \endcode

// \n \section views_columns_arithmetic_operations Arithmetic Operations
// <hr>
//
// Both dense and sparse columns can be used in all arithmetic operations that any other dense or
// sparse column vector can be used in. The following example gives an impression of the use of
// dense columns within arithmetic operations. All operations (addition, subtraction, multiplication,
// scaling, ...) can be performed on all possible combinations of dense and sparse columns with
// fitting element types:

   \code
   blaze::DynamicVector<double,blaze::columnVector> a( 2UL, 2.0 ), b;
   blaze::CompressedVector<double,blaze::columnVector> c( 2UL );
   c[1] = 3.0;

   typedef blaze::DynamicMatrix<double,blaze::columnMajor>  MatrixType;
   MatrixType A( 2UL, 4UL );  // Non-initialized 2x4 matrix

   typedef blaze::DenseColumn<DenseMatrix>  RowType;
   RowType col0( column( A, 0UL ) );  // Reference to the 0th column of A

   col0[0] = 0.0;           // Manual initialization of the 0th column of A
   col0[1] = 0.0;
   column( A, 1UL ) = 1.0;  // Homogeneous initialization of the 1st column of A
   column( A, 2UL ) = a;    // Dense vector initialization of the 2nd column of A
   column( A, 3UL ) = c;    // Sparse vector initialization of the 3rd column of A

   b = col0 + a;                 // Dense vector/dense vector addition
   b = c + column( A, 1UL );     // Sparse vector/dense vector addition
   b = col0 * column( A, 2UL );  // Component-wise vector multiplication

   column( A, 1UL ) *= 2.0;     // In-place scaling of the 1st column
   b = column( A, 1UL ) * 2.0;  // Scaling of the 1st column
   b = 2.0 * column( A, 1UL );  // Scaling of the 1st column

   column( A, 2UL ) += a;                 // Addition assignment
   column( A, 2UL ) -= c;                 // Subtraction assignment
   column( A, 2UL ) *= column( A, 0UL );  // Multiplication assignment

   double scalar = trans( c ) * column( A, 1UL );  // Scalar/dot/inner product between two vectors

   A = column( A, 1UL ) * trans( c );  // Outer product between two vectors
   \endcode

// \n \section views_columns_non_fitting_storage_order Views on Matrices with Non-Fitting Storage Order
// <hr>
//
// Especially noteworthy is that column views can be created for both row-major and column-major
// matrices. Whereas the interface of a row-major matrix only allows to traverse a row directly
// and the interface of a column-major matrix only allows to traverse a column, via views it is
// possible to traverse a row of a column-major matrix or a column of a row-major matrix. For
// instance:

   \code
   typedef blaze::CompressedMatrix<int,rowMajor>  MatrixType;
   typedef blaze::SparseColumn<MatrixType>        ColumnType;

   MatrixType A( 64UL, 32UL );
   // ... Resizing and initialization

   // Creating a reference to the 31st column of a row-major matrix A
   ColumnType col1 = column( A, 1UL );

   for( ColumnType::Iterator it=col1.begin(); it!=col1.end(); ++it ) {
      // ...
   }
   \endcode

// However, please note that creating a column view on a matrix stored in a row-major fashion
// can result in a considerable performance decrease in comparison to a view on a matrix with
// a fitting storage orientation. This is due to the non-contiguous storage of the matrix
// elements. Therefore care has to be taken in the choice of the most suitable storage order:

   \code
   // Setup of two row-major matrices
   CompressedMatrix<double,rowMajor> A( 128UL, 128UL );
   CompressedMatrix<double,rowMajor> B( 128UL, 128UL );
   // ... Resizing and initialization

   // The computation of the 15th column of the multiplication between A and B ...
   CompressedVector<double,columnVector> x = column( A * B, 15UL );

   // ... is essentially the same as the following computation, which multiplies
   // the 15th column of the row-major matrix B with A.
   CompressedVector<double,columnVector> x = A * column( B, 15UL );
   \endcode

// Although \b Blaze performs the resulting matrix/vector multiplication as efficiently as possible
// using a column-major storage order for matrix B would result in a more efficient evaluation.
//
// \n <center> Previous: \ref views_rows &nbsp; &nbsp; Next: \ref addition </center>
*/
//*************************************************************************************************


//**Addition***************************************************************************************
/*!\page addition Addition
//
// <center> Previous: \ref views_columns &nbsp; &nbsp; Next: \ref subtraction </center> \n
//
// The addition of vectors and matrices is as intuitive as the addition of scalar values. For both
// the vector addition as well as the matrix addition the addition operator can be used. It even
// enables the addition of dense and sparse vectors as well as the addition of dense and sparse
// matrices:

   \code
   blaze::DynamicVector<int>      v1( 5UL ), v3;
   blaze::CompressedVector<float> v2( 5UL );

   // ... Initializing the vectors

   v3 = v1 + v2;  // Addition of a two column vectors of different data type
   \endcode

   \code
   blaze::DynamicMatrix<float,rowMajor>        M1( 7UL, 3UL );
   blaze::CompressedMatrix<size_t,columnMajor> M2( 7UL, 3UL ), M3;

   // ... Initializing the matrices

   M3 = M1 + M2;  // Addition of a row-major and a column-major matrix of different data type
   \endcode

// Note that it is necessary that both operands have exactly the same dimensions. Violating this
// precondition results in an exception. Also note that in case of vectors it is only possible to
// add vectors with the same transpose flag:

   \code
   blaze::DynamicVector<int,columnVector>   v1( 5UL );
   blaze::CompressedVector<float,rowVector> v2( 5UL );

   v1 + v2;           // Compilation error: Cannot add a column vector and a row vector
   v1 + trans( v2 );  // OK: Addition of two column vectors
   \endcode

// In case of matrices, however, it is possible to add row-major and column-major matrices. Note
// however that in favor of performance the addition of two matrices with the same storage order
// is favorable. The same argument holds for the element type: In case two vectors or matrices
// with the same element type are added, the performance can be much higher due to vectorization
// of the operation.

   \code
   blaze::DynamicVector<double>v1( 100UL ), v2( 100UL ), v3;

   // ... Initialization of the vectors

   v3 = v1 + v2;  // Vectorized addition of two double precision vectors
   \endcode

   \code
   blaze::DynamicMatrix<float> M1( 50UL, 70UL ), M2( 50UL, 70UL ), M3;

   // ... Initialization of the matrices

   M3 = M1 + M2;  // Vectorized addition of two row-major, single precision dense matrices
   \endcode

// \n <center> Previous: \ref views_columns &nbsp; &nbsp; Next: \ref subtraction </center>
*/
//*************************************************************************************************


//**Subtraction************************************************************************************
/*!\page subtraction Subtraction
//
// <center> Previous: \ref addition &nbsp; &nbsp; Next: \ref scalar_multiplication </center> \n
//
// The subtraction of vectors and matrices works exactly as intuitive as the addition, but with
// the subtraction operator. For both the vector subtraction as well as the matrix subtraction
// the subtraction operator can be used. It also enables the subtraction of dense and sparse
// vectors as well as the subtraction of dense and sparse matrices:

   \code
   blaze::DynamicVector<int>      v1( 5UL ), v3;
   blaze::CompressedVector<float> v2( 5UL );

   // ... Initializing the vectors

   v3 = v1 - v2;  // Subtraction of a two column vectors of different data type


   blaze::DynamicMatrix<float,rowMajor>        M1( 7UL, 3UL );
   blaze::CompressedMatrix<size_t,columnMajor> M2( 7UL, 3UL ), M3;

   // ... Initializing the matrices

   M3 = M1 - M2;  // Subtraction of a row-major and a column-major matrix of different data type
   \endcode

// Note that it is necessary that both operands have exactly the same dimensions. Violating this
// precondition results in an exception. Also note that in case of vectors it is only possible to
// subtract vectors with the same transpose flag:

   \code
   blaze::DynamicVector<int,columnVector>   v1( 5UL );
   blaze::CompressedVector<float,rowVector> v2( 5UL );

   v1 - v2;           // Compilation error: Cannot subtract a row vector from a column vector
   v1 - trans( v2 );  // OK: Subtraction of two column vectors
   \endcode

// In case of matrices, however, it is possible to subtract row-major and column-major matrices.
// Note however that in favor of performance the subtraction of two matrices with the same storage
// order is favorable. The same argument holds for the element type: In case two vectors or matrices
// with the same element type are added, the performance can be much higher due to vectorization
// of the operation.

   \code
   blaze::DynamicVector<double>v1( 100UL ), v2( 100UL ), v3;

   // ... Initialization of the vectors

   v3 = v1 - v2;  // Vectorized subtraction of two double precision vectors


   blaze::DynamicMatrix<float> M1( 50UL, 70UL ), M2( 50UL, 70UL ), M3;

   // ... Initialization of the matrices

   M3 = M1 - M2;  // Vectorized subtraction of two row-major, single precision dense matrices
   \endcode

// \n <center> Previous: \ref addition &nbsp; &nbsp; Next: \ref scalar_multiplication </center>
*/
//*************************************************************************************************


//**Scalar Multiplication**************************************************************************
/*!\page scalar_multiplication Scalar Multiplication
//
// <center> Previous: \ref subtraction &nbsp; &nbsp; Next: \ref componentwise_multiplication </center> \n
//
// The scalar multiplication is the multiplication of a scalar value with a vector or a matrix.
// In \b Blaze it is possible to use all built-in/fundamental data types except bool as scalar
// values. Additionally, it is possible to use std::complex values with the same built-in data
// types as element type.

   \code
   blaze::StaticVector<int,3UL> v1( 1, 2, 3 );

   blaze::DynamicVector<double>   v2 = v1 * 1.2;
   blaze::CompressedVector<float> v3 = -0.3F * v1;
   \endcode

   \code
   blaze::StaticMatrix<int,3UL,2UL> M1( 1, 2, 3, 4, 5, 6 );

   blaze::DynamicMatrix<double>   M2 = M1 * 1.2;
   blaze::CompressedMatrix<float> M3 = -0.3F * M1;
   \endcode

// Vectors and matrices cannot be used for as scalar value for scalar multiplications (see the
// following example). However, each vector and matrix provides the \c scale() function, which
// can be used to scale a vector or matrix element-wise with arbitrary scalar data types:

   \code
   blaze::CompressedMatrix< blaze::StaticMatrix<int,3UL,3UL> > M1;
   blaze::StaticMatrix<int,3UL,3UL> scalar;

   M1 * scalar;  // No scalar multiplication, but matrix/matrix multiplication

   M1.scale( scalar );  // Scalar multiplication
   \endcode

// \n <center> Previous: \ref subtraction &nbsp; &nbsp; Next: \ref componentwise_multiplication </center>
*/
//*************************************************************************************************


//**Vector/Vector Multiplication*******************************************************************
/*!\page vector_vector_multiplication Vector/Vector Multiplication
//
// <center> Previous: \ref scalar_multiplication &nbsp; &nbsp; Next: \ref matrix_vector_multiplication </center> \n
//
// \n \section componentwise_multiplication Componentwise Multiplication
// <hr>
//
// Multiplying two vectors with the same transpose flag (i.e. either blaze::columnVector or
// blaze::rowVector) via the multiplication operator results in a componentwise multiplication
// of the two vectors:

   \code
   using blaze::DynamicVector;
   using blaze::CompressedVector;

   CompressedVector<int,columnVector> v1( 17UL );
   DynamicVector<int,columnVector>    v2( 17UL );

   StaticVector<double,10UL,rowVector> v3;
   DynamicVector<double,rowVector>     v4( 10UL );

   // ... Initialization of the vectors

   CompressedVector<int,columnVector> v5( v1 * v2 );  // Componentwise multiplication of a sparse and
                                                      // a dense column vector. The result is a sparse
                                                      // column vector.
   DynamicVector<double,rowVector>    v6( v3 * v4 );  // Componentwise multiplication of two dense row
                                                      // vectors. The result is a dense row vector.
   \endcode

// \n \section inner_product Inner Product / Scalar Product / Dot Product
// <hr>
//
// The multiplication between a row vector and a column vector results in an inner product between
// the two vectors:

   \code
   blaze::StaticVector<int,3UL,rowVector> v1( 2, 5, -1 );

   blaze::DynamicVector<int,columnVector> v2( 3UL );
   v2[0] = -1;
   v2[1] = 3;
   v2[2] = -2;

   int result = v1 * v2;  // Results in the value 15
   \endcode

// The \c trans() function can be used to transpose a vector as necessary:

   \code
   blaze::StaticVector<int,3UL,rowVector> v1(  2, 5, -1 );
   blaze::StaticVector<int,3UL,rowVector> v2( -1, 3, -2 );

   int result = v1 * trans( v2 );  // Also results in the value 15
   \endcode

// Alternatively, the comma operator can used for any combination of vectors (row or column vectors)
// to perform an inner product:

   \code
   blaze::StaticVector<int,3UL,rowVector> v1(  2, 5, -1 );
   blaze::StaticVector<int,3UL,rowVector> v2( -1, 3, -2 );

   int result = (v1,v2);  // Inner product between two row vectors
   \endcode

// Please note the brackets embracing the inner product expression. Due to the low precedence
// of the comma operator (lower even than the assignment operator) these brackets are strictly
// required for a correct evaluation of the inner product.
//
//
// \n \section outer_product Outer Product
// <hr>
//
// The multiplication between a column vector and a row vector results in the outer product of
// the two vectors:

   \code
   blaze::StaticVector<int,3UL,columnVector> v1( 2, 5, -1 );

   blaze::DynamicVector<int,rowVector> v2( 3UL );
   v2[0] = -1;
   v2[1] = 3;
   v2[2] = -2;

   StaticMatrix<int,3UL,3UL> M1 = v1 * v2;
   \endcode

// The \c trans() function can be used to transpose a vector as necessary:

   \code
   blaze::StaticVector<int,3UL,rowVector> v1(  2, 5, -1 );
   blaze::StaticVector<int,3UL,rowVector> v2( -1, 3, -2 );

   int result = trans( v1 ) * v2;
   \endcode

// \n \section cross_product Cross Product
// <hr>
//
// Two column vectors can be multiplied via the cross product. The cross product between two
// vectors \f$ a \f$ and \f$ b \f$ is defined as

   \f[
   \left(\begin{array}{*{1}{c}}
   c_0 \\
   c_1 \\
   c_2 \\
   \end{array}\right)
   =
   \left(\begin{array}{*{1}{c}}
   a_1 b_2 - a_2 b_1 \\
   a_2 b_0 - a_0 b_2 \\
   a_0 b_1 - a_1 b_0 \\
   \end{array}\right).
   \f]

// Due to the absence of a \f$ \times \f$ operator in the C++ language, the cross product is
// realized via the modulo operator (i.e. \c operator%):

   \code
   blaze::StaticVector<int,3UL,columnVector> v1( 2, 5, -1 );

   blaze::DynamicVector<int,rowVector> v2( 3UL );
   v2[0] = -1;
   v2[1] = 3;
   v2[2] = -2;

   blaze::StaticVector<int,3UL,columnVector> v3( v1 % v2 );
   \endcode

// Please note that the cross product is restricted to three dimensional (dense and sparse)
// vectors.
//
// \n <center> Previous: \ref scalar_multiplication &nbsp; &nbsp; Next: \ref matrix_vector_multiplication </center>
*/
//*************************************************************************************************


//**Matrix/Vector Multiplication*******************************************************************
/*!\page matrix_vector_multiplication Matrix/Vector Multiplication
//
// <center> Previous: \ref outer_product &nbsp; &nbsp; Next: \ref matrix_matrix_multiplication </center> \n
//
// In \b Blaze matrix/vector multiplications can be as intuitively formulated as in mathematical
// textbooks. Just as in textbooks there are two different multiplications between a matrix and
// a vector: a matrix/column vector multiplication and a row vector/matrix multiplication:

   \code
   using blaze::StaticVector;
   using blaze::DynamicVector;
   using blaze::DynamicMatrix;

   DynamicMatrix<int>                  M1( 39UL, 12UL );
   StaticVector<int,12UL,columnVector> v1;

   // ... Initialization of the matrix and the vector

   DynamicVector<int,columnVector> v2 = M1 * v1;           // Matrix/column vector multiplication
   DynamicVector<int,rowVector>    v3 = trans( v1 ) * M1;  // Row vector/matrix multiplication
   \endcode

// Note that the storage order of the matrix poses no restrictions on the operation. Also note,
// that the highest performance for a multiplication between a dense matrix and a dense vector can
// be achieved if both the matrix and the vector have the same scalar element type.
//
// \n <center> Previous: \ref outer_product &nbsp; &nbsp; Next: \ref matrix_matrix_multiplication </center>
*/
//*************************************************************************************************


//**Matrix/Matrix Multiplication*******************************************************************
/*!\page matrix_matrix_multiplication Matrix/Matrix Multiplication
//
// <center> Previous: \ref matrix_vector_multiplication &nbsp; &nbsp; Next: \ref openmp_parallelization </center> \n
//
// The matrix/matrix multiplication can be formulated exactly as in mathematical textbooks:

   \code
   using blaze::DynamicMatrix;
   using blaze::CompressedMatrix;

   DynamicMatrix<double>   M1( 45UL, 85UL );
   CompressedMatrix<float> M2( 85UL, 37UL );

   // ... Initialization of the matrices

   DynamicMatrix<double> M3 = M1 * M2;
   \endcode

// The storage order of the two matrices poses no restrictions on the operation, all variations
// are possible. Note however that the highest performance for a multiplication between two dense
// matrices can be expected for two matrices with the same scalar element type.
//
// \n <center> Previous: \ref matrix_vector_multiplication &nbsp; &nbsp; Next: \ref openmp_parallelization </center>
*/
//*************************************************************************************************


//**OpenMP Parallelization*************************************************************************
/*!\page openmp_parallelization OpenMP Parallelization
//
// <center> Previous: \ref matrix_matrix_multiplication &nbsp; &nbsp; Next: \ref serial_execution </center> \n
//
// One of the main motivations of the \b Blaze 1.x releases was to achieve maximum performance
// on a single CPU core for all possible operations. However, today's CPUs are not single core
// anymore, but provide several (homogeneous or heterogeneous) compute cores. In order to fully
// exploit the performance potential of a multicore CPU, computations have to be parallelized
// across all available cores of a CPU. Therefore, starting with \b Blaze 2.0, the \b Blaze
// library provides shared memory parallelization with OpenMP.
//
//
// \n \section openmp_setup OpenMP Setup
// <hr>
//
// To enable OpenMP-based parallelization, all that needs to be done is to explicitly specify
// the use of OpenMP on the command line:

   \code
   -fopenmp   // GNU C++ compiler
   -openmp    // Intel C++ compiler
   /openmp    // Visual Studio
   \endcode

// This simple action will cause the \b Blaze library to automatically try to run all operations
// in parallel with the specified number of threads.
//
// As common for OpenMP, the number of threads can be specified either via an environment variable

   \code
   export OMP_NUM_THREADS=4
   \endcode

// or via an explicit call to the \c omp_set_num_threads() function:

   \code
   omp_set_num_threads( 4 );
   \endcode

// Either way, the best performance can be expected if the specified number of threads matches
// the available number of cores.
//
//
// \n \section openmp_configuration OpenMP Configuration
// <hr>
//
// Note that \b Blaze is not unconditionally running an operation in parallel. In case \b Blaze
// deems the parallel execution as counterproductive for the overall performance, the operation
// is executed serially. One of the main reasons for not executing an operation in parallel is
// the size of the operands. For instance, a vector addition is only executed in parallel if the
// size of both vector operands exceeds a certain threshold. Otherwise, the performance could
// seriously decrease due to the overhead caused by the thread setup. However, in order to be
// able to adjust the \b Blaze library to a specific system, it is possible to configure these
// thresholds manually. All OpenMP thresholds are contained within the configuration file
// <em>./blaze/config/Thresholds.h</em>.
//
//
// \n \section openmp_first_touch First Touch Policy
// <hr>
//
// So far the \b Blaze library does not (yet) automatically initialize dynamic memory according
// to the first touch principle. Consider for instance the following vector triad example:

   \code
   using blaze::columnVector;

   const size_t N( 1000000UL );

   blaze::DynamicVector<double,columnVector> a( N ), b( N ), c( N ), d( N );

   // Initialization of the vectors b, c, and d
   for( size_t i=0UL; i<N; ++i ) {
      b[i] = rand<double>();
      c[i] = rand<double>();
      d[i] = rand<double>();
   }

   // Performing a vector triad
   a = b + c * d;
   \endcode

// If this code, which is prototypical for many OpenMP applications that have not been optimized
// for ccNUMA architectures, is run across several locality domains (LD), it will not scale
// beyond the maximum performance achievable on a single LD if the working set does not fit into
// the cache. This is because the initialization loop is executed by a single thread, writing to
// \c b, \c c, and \c d for the first time. Hence, all memory pages belonging to those arrays will
// be mapped into a single LD.
//
// As mentioned above, this problem can be solved by performing vector initialization in parallel:

   \code
   // ...

   // Initialization of the vectors b, c, and d
   #pragma omp parallel for
   for( size_t i=0UL; i<N; ++i ) {
      b[i] = rand<double>();
      c[i] = rand<double>();
      d[i] = rand<double>();
   }

   // ...
   \endcode

// This simple modification makes a huge difference on ccNUMA in memory-bound situations (as for
// instance in all BLAS level 1 operations and partially BLAS level 2 operations). Therefore, in
// order to achieve the maximum possible performance, it is imperative to initialize the memory
// according to the later use of the data structures.
//
// \n <center> Previous: \ref matrix_matrix_multiplication &nbsp; &nbsp; Next: \ref serial_execution </center>
*/
//*************************************************************************************************


//**Serial Execution*******************************************************************************
/*!\page serial_execution Serial Execution
//
// <center> Previous: \ref openmp_parallelization &nbsp; &nbsp; Next: \ref vector_serialization </center> \n
//
// Sometimes it may be necessary to enforce the serial execution of specific operations. For this
// purpose, the \b Blaze library offers two possible options.
//
// The first option is the temporary and local enforcement of a serial execution via the
// \c BLAZE_SERIAL_SECTION:

   \code
   using blaze::rowMajor;
   using blaze::columnVector;

   blaze::DynamicMatrix<double,rowMajor> A;
   blaze::DynamicVector<double,columnVector> b, c, d, x, y, z;

   // ... Resizing and initialization

   // Parallel execution
   // If possible and beneficial for performance the following operation is executed in parallel.
   x = A * b;

   // Serial execution
   // All operations executed within the serial section are guaranteed to be executed in
   // serial (even if a parallel execution would be possible and/or beneficial).
   BLAZE_SERIAL_SECTION
   {
      y = A * c;
      z = A * d;
   }

   // Parallel execution continued
   // ...
   \endcode

// Within the scope of the \c BLAZE_SERIAL_SECTION, all operations are guaranteed to run in serial.
// Outside the scope of the serial section, all operations are run in parallel (if beneficial for
// the performance).
//
// The second option is the general deactivation of the parallel execution (even in case OpenMP
// is enabled on the command line). This can be achieved via the <em>./blaze/config/SMP.h</em>
// configuration file: In case the \c BLAZE_ENABLE_PARALLEL_EXECUTION switch is set to 0, the
// shared-memory parallelization is deactivated altogether.
//
// \n <center> Previous: \ref openmp_parallelization &nbsp; &nbsp; Next: \ref vector_serialization </center>
*/
//*************************************************************************************************


//**Vector Serialization***************************************************************************
/*!\page vector_serialization Vector Serialization
//
// <center> Previous: \ref matrix_matrix_multiplication &nbsp; &nbsp; Next: \ref matrix_serialization </center> \n
//
// Sometimes it is necessary to store vector and/or matrices on disk, for instance for storing
// results or for sharing specific setups with other people. The \b Blaze math serialization
// module provides the according functionality to create platform independent, portable, binary
// representations of vectors and matrices that can be used to store the \b Blaze data structures
// without loss of precision and to reliably transfer them from one machine to another.
//
// The following example demonstrates the (de-)serialization of dense and sparse vectors:

   \code
   using blaze::columnVector;
   using blaze::rowVector;

   // Serialization of both vectors
   {
      blaze::StaticVector<double,5UL,rowVector> d;
      blaze::CompressedVector<int,columnVector> s;

      // ... Resizing and initialization

      // Creating an archive that writes into a the file "vectors.blaze"
      blaze::Archive<std::ofstream> archive( "vectors.blaze" );

      // Serialization of both vectors into the same archive. Note that d lies before s!
      archive << d << s;
   }

   // Reconstitution of both vectors
   {
      blaze::DynamicVector<double,rowVector> d1;
      blaze::DynamicVector<int,rowVector> d2;

      // Creating an archive that reads from the file "vectors.blaze"
      blaze::Archive<std::ifstream> archive( "vectors.blaze" );

      // Reconstituting the former d vector into d1. Note that it is possible to reconstitute
      // the vector into a differrent kind of vector (StaticVector -> DynamicVector), but that
      // the type of elements has to be the same.
      archive >> d1;

      // Reconstituting the former s vector into d2. Note that is is even possible to reconstitute
      // a sparse vector as a dense vector (also the reverse is possible) and that a column vector
      // can be reconstituted as row vector (and vice versa). Note however that also in this case
      // the type of elements is the same!
      archive >> d2
   }
   \endcode

// The (de-)serialization of vectors is not restricted to vectors of built-in data type, but can
// also be used for vectors with vector or matrix element type:

   \code
   // Serialization
   {
      blaze::CompressedVector< blaze::DynamicVector< blaze::complex<double> > > vec;

      // ... Resizing and initialization

      // Creating an archive that writes into a the file "vector.blaze"
      blaze::Archive<std::ofstream> archive( "vector.blaze" );

      // Serialization of the vector into the archive
      archive << vec;
   }

   // Deserialization
   {
      blaze::CompressedVector< blaze::DynamicVector< blaze::complex<double> > > vec;

      // Creating an archive that reads from the file "vector.blaze"
      blaze::Archive<std::ifstream> archive( "vector.blaze" );

      // Reconstitution of the vector from the archive
      archive >> vec;
   }
   \endcode

// As the examples demonstrates, the vector serialization offers an enormous flexibility. However,
// several actions result in errors:
//
//  - vectors cannot be reconstituted as matrices (and vice versa)
//  - the element type of the serialized and reconstituted vector must match, which means
//    that on the source and destination platform the general type (signed/unsigned integral
//    or floating point) and the size of the type must be exactly the same
//  - when reconstituting a StaticVector, its size must match the size of the serialized vector
//
// In case an error is encountered during (de-)serialization, a \a std::runtime_exception is
// thrown.
//
// \n <center> Previous: \ref matrix_matrix_multiplication &nbsp; &nbsp; Next: \ref matrix_serialization </center>
*/
//*************************************************************************************************


//**Matrix Serialization***************************************************************************
/*!\page matrix_serialization Matrix Serialization
//
// <center> Previous: \ref vector_serialization &nbsp; &nbsp; Next: \ref intra_statement_optimization </center> \n
//
// The serialization of matrices works in the same manner as the serialization of vectors. The
// following example demonstrates the (de-)serialization of dense and sparse matrices:

   \code
   using blaze::rowMajor;
   using blaze::columnMajor;

   // Serialization of both matrices
   {
      blaze::StaticMatrix<double,3UL,5UL,rowMajor> D;
      blaze::CompressedMatrix<int,columnMajor> S;

      // ... Resizing and initialization

      // Creating an archive that writes into a the file "matrices.blaze"
      blaze::Archive<std::ofstream> archive( "matrices.blaze" );

      // Serialization of both matrices into the same archive. Note that D lies before S!
      archive << D << S;
   }

   // Reconstitution of both matrices
   {
      blaze::DynamicMatrix<double,rowMajor> D1;
      blaze::DynamicMatrix<int,rowMajor> D2;

      // Creating an archive that reads from the file "matrices.blaze"
      blaze::Archive<std::ifstream> archive( "matrices.blaze" );

      // Reconstituting the former D matrix into D1. Note that it is possible to reconstitute
      // the matrix into a differrent kind of matrix (StaticMatrix -> DynamicMatrix), but that
      // the type of elements has to be the same.
      archive >> D1;

      // Reconstituting the former S matrix into D2. Note that is is even possible to reconstitute
      // a sparse matrix as a dense matrix (also the reverse is possible) and that a column-major
      // matrix can be reconstituted as row-major matrix (and vice versa). Note however that also
      // in this case the type of elements is the same!
      archive >> D2
   }
   \endcode

// Note that also in case of matrices it is possible to (de-)serialize matrices with vector or
// matrix elements:

   \code
   // Serialization
   {
      blaze::CompressedMatrix< blaze::DynamicMatrix< blaze::complex<double> > > mat;

      // ... Resizing and initialization

      // Creating an archive that writes into a the file "matrix.blaze"
      blaze::Archive<std::ofstream> archive( "matrix.blaze" );

      // Serialization of the matrix into the archive
      archive << mat;
   }

   // Deserialization
   {
      blaze::CompressedMatrix< blaze::DynamicMatrix< blaze::complex<double> > > mat;

      // Creating an archive that reads from the file "matrix.blaze"
      blaze::Archive<std::ifstream> archive( "matrix.blaze" );

      // Reconstitution of the matrix from the archive
      archive >> mat;
   }
   \endcode

// Note that just as the vector serialization, the matrix serialization is restricted by a
// few important rules:
//
//  - matrices cannot be reconstituted as vectors (and vice versa)
//  - the element type of the serialized and reconstituted matrix must match, which means
//    that on the source and destination platform the general type (signed/unsigned integral
//    or floating point) and the size of the type must be exactly the same
//  - when reconstituting a StaticMatrix, the number of rows and columns must match those of
//    the serialized matrix
//
// In case an error is encountered during (de-)serialization, a \a std::runtime_exception is
// thrown.
//
// \n <center> Previous: \ref vector_serialization &nbsp; &nbsp; Next: \ref intra_statement_optimization </center> \n
*/
//*************************************************************************************************


//**Intra-Statement Optimization*******************************************************************
/*!\page intra_statement_optimization Intra-Statement Optimization
//
// <center> Previous: \ref matrix_serialization &nbsp; &nbsp; Next: \ref configuration_files </center> \n
//
//
// One of the prime features of the \b Blaze library is the automatic intra-statement optimization.
// In order to optimize the overall performance of every single statement \b Blaze attempts to
// rearrange the operands based on their types. For instance, the following addition of dense and
// sparse vectors

   \code
   blaze::DynamicVector<double> d1, d2, d3;
   blaze::CompressedVector<double> s1;

   // ... Resizing and initialization

   d3 = d1 + s1 + d2;
   \endcode

// is automatically rearranged and evaluated as

   \code
   // ...
   d3 = d1 + d2 + s1;  // <- Note that s1 and d2 have been rearranged
   \endcode

// This order of operands is highly favorable for the overall performance since the addition of
// the two dense vectors \c d1 and \c d2 can be handled much more efficiently in a vectorized
// fashion.
//
// This intra-statement optimization can have a tremendous effect on the performance of a statement.
// Consider for instance the following computation:

   \code
   blaze::DynamicMatrix<double> A, B;
   blaze::DynamicVector<double> x, y;

   // ... Resizing and initialization

   y = A * B * x;
   \endcode

// Since multiplications are evaluated from left to right, this statement would result in a
// matrix/matrix multiplication, followed by a matrix/vector multiplication. However, if the
// right subexpression is evaluated first, the performance can be dramatically improved since the
// matrix/matrix multiplication can be avoided in favor of a second matrix/vector multiplication.
// The \b Blaze library exploits this by automatically restructuring the expression such that the
// right multiplication is evaluated first:

   \code
   // ...
   y = A * ( B * x );
   \endcode

// Note however that although this intra-statement optimization may result in a measurable or
// even significant performance improvement, this behavior may be undesirable for several reasons,
// for instance because of numerical stability. Therefore, in case the order of evaluation matters,
// the best solution is to be explicit and to separate a statement into several statements:

   \code
   blaze::DynamicVector<double> d1, d2, d3;
   blaze::CompressedVector<double> s1;

   // ... Resizing and initialization

   d3  = d1 + s1;  // Compute the dense vector/sparse vector addition first ...
   d3 += d2;       // ... and afterwards add the second dense vector
   \endcode

   \code
   // ...
   blaze::DynamicMatrix<double> A, B, C;
   blaze::DynamicVector<double> x, y;

   // ... Resizing and initialization

   C = A * B;  // Compute the left-hand side matrix-matrix multiplication first ...
   y = C * x;  // ... before the right-hand side matrix-vector multiplication
   \endcode

// \n <center> Previous: \ref matrix_serialization &nbsp; &nbsp; Next: \ref configuration_files </center> \n
*/
//*************************************************************************************************


//**Configuration Files****************************************************************************
/*!\page configuration_files Configuration Files
//
// <center> Previous: \ref intra_statement_optimization </center> \n
//
//
// \tableofcontents
//
//
// Sometimes it might necessary to adapt \b Blaze to specific requirements. For this purpose
// \b Blaze provides several configuration files in the <em>./blaze/config/</em> subdirectory,
// which provide ample opportunity to customize internal settings, behavior, and thresholds.
// This chapter explains the most important of these configuration files.
//
//
// \n \section transpose_flag Default Vector Storage
// <hr>
//
// The \b Blaze default is that all vectors are created as column vectors (if not specified
// explicitly):

   \code
   blaze::StaticVector<double,3UL> x;  // Creates a 3-dimensional static column vector
   \endcode

// The header file <em>./blaze/config/TransposeFlag.h</em> allows the configuration of the default
// vector storage (i.e. the default transpose flag of the vectors). Via the \a defaultTransposeFlag
// value the default transpose flag for all vector of the \b Blaze library can be specified:

   \code
   const bool defaultTransposeFlag = columnVector;
   \endcode

// Valid settings for the \a defaultTransposeFlag are blaze::rowVector and blaze::columnVector.
//
//
// \n \section storage_order Default Matrix Storage
// <hr>
//
// Matrices are by default created as row-major matrices:

   \code
   blaze::StaticMatrix<double,3UL,3UL>  A;  // Creates a 3x3 row-major matrix
   \endcode

// The header file <em>./blaze/config/StorageOrder.h</em> allows the configuration of the default
// matrix storage order. Via the \a defaultStorageOrder value the default storage order for all
// matrices of the \b Blaze library can be specified.

   \code
   const bool defaultStorageOrder = rowMajor;
   \endcode

// Valid settings for the \a defaultStorageOrder are blaze::rowMajor and blaze::columnMajor.
//
//
// \n \section vectorization Vectorization
//
// In order to achieve maximum performance and to exploit the compute power of a target platform
// the \b Blaze library attempts to vectorize all linear algebra operations by SSE, AVX, and/or
// MIC intrinsics, depending on which instruction set is available. However, it is possible to
// disable the vectorization entirely by the compile time switch in the configuration file
// <em>./blaze/config/Vectorization.h</em>:

   \code
   #define BLAZE_USE_VECTORIZATION 1
   \endcode

// In case the switch is set to 1, vectorization is enabled and the \b Blaze library is allowed
// to use intrinsics to speed up computations. In case the switch is set to 0, vectorization is
// disabled entirely and the \b Blaze library chooses default, non-vectorized functionality for
// the operations. Note that deactivating the vectorization may pose a severe performance
// limitation for a large number of operations!
//
//
// \n \section thresholds Thresholds
//
// \b Blaze provides several thresholds that can be adapted to the characteristics of the target
// platform. For instance, the \a DMATDVECMULT_THRESHOLD specifies the threshold between the
// application of the custom Blaze kernels for small dense matrix/dense vector multiplications
// and the BLAS kernels for large multiplications. All thresholds, including the thresholds for
// the OpenMP-based parallelization, are contained within the configuration file
// <em>./blaze/config/Thresholds.h</em>.
//
//
// \n \section streaming Streaming (Non-Temporal Stores)
//
// For vectors and matrices that don't fit into the cache anymore non-temporal stores can provide
// a significant performance advantage of about 20%. However, this advantage is only in effect in
// case the memory bandwidth of the target architecture is maxed out. If the target architecture's
// memory bandwidth cannot be exhausted the use of non-temporal stores can decrease performance
// instead of increasing it.
//
// The configuration file <em>./blaze/config/Streaming.h</em> provides a compile time switch that
// can be used to (de-)activate streaming:

   \code
   const bool useStreaming = true;
   \endcode

// If \a useStreaming is set to \a true streaming is enabled, if it is set to \a false streaming
// is disabled. It is recommended to consult the target architecture's white papers to decide
// whether streaming is beneficial or hurtful for performance.
//
//
// \n <center> Previous: \ref intra_statement_optimization </center>
*/
//*************************************************************************************************

#endif
