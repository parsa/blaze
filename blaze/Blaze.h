//=================================================================================================
/*!
//  \file blaze/Blaze.h
//  \brief Primary include file of the Blaze library
//
//  Copyright (C) 2011 Klaus Iglberger - All Rights Reserved
//
//  This file is part of the Blaze library. This library is free software; you can redistribute
//  it and/or modify it under the terms of the GNU General Public License as published by the
//  Free Software Foundation; either version 3, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
//  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//  See the GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License along with a special
//  exception for linking and compiling against the Blaze library, the so-called "runtime
//  exception"; see the file COPYING. If not, see http://www.gnu.org/licenses/.
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
//! Namespace of the Blaze C++ math library.
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
//                <li> \ref view_types </li>
//                <li> \ref view_operations </li>
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
//          <li> Serialization
//             <ul>
//                <li> \ref vector_serialization </li>
//                <li> \ref matrix_serialization </li>
//             </ul>
//          </li>
//       </ul>
//    </li>
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
// have a BLAS library installed on your system, *Blaze* will still work and will not be reduced
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
// subdirectory can be adapted. For instance, in the header file <em>./blaze/config/StorageOrder.h</em>
// the default matrix storage order (i.e. row-major or column-major) can be specified.
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
// long tutorial covers all aspects of the \b Blaze math library.
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
// classes and functions of Blaze are contained in the blaze namespace.\n\n
//
// The output of the last line of this small program is

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
// The \b Blaze library currently offers two dense vector types (\ref vector_types_static_vector
// and \ref vector_types_dynamic_vector) and one sparse vector type (\ref vector_types_compressed_vector).
// All vectors can be specified as either column vectors

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
// statically allocated elements of arbitrary type. The type of the elements, the number of
// elements, and the transpose flag of the vector can be specified via the three template
// parameters:

   \code
   template< typename Type, size_t N, bool TF >
   class StaticVector;
   \endcode

//  - Type: specifies the type of the vector elements. StaticVector can be used with any
//          non-cv-qualified, non-reference, non-pointer element type.
//  - N   : specifies the total number of vector elements. It is expected that StaticVector is
//          only used for tiny and small vectors.
//  - TF  : specifies whether the vector is a row vector (\c blaze::rowVector) or a column
//          vector (\c blaze::columnVector). The default value is \c blaze::columnVector.
//
//
// \n \section vector_types_dynamic_vector DynamicVector
// <hr>
//
// The blaze::DynamicVector class template is the representation of an arbitrary sized vector
// with dynamically allocated elements of arbitrary type. The type of the elements and the
// transpose flag of the vector can be specified via the two template parameters:

   \code
   template< typename Type, bool TF >
   class DynamicVector;
   \endcode

//  - Type: specifies the type of the vector elements. DynamicVector can be used with any
//          non-cv-qualified, non-reference, non-pointer element type.
//  - TF  : specifies whether the vector is a row vector (\c blaze::rowVector) or a column
//          vector (\c blaze::columnVector). The default value is \c blaze::columnVector.
//
//
// \n \section vector_types_compressed_vector CompressedVector
// <hr>
//
// The blaze::CompressedVector class is the representation of an arbitrarily sized sparse
// vector, which stores only non-zero elements of arbitrary type. The type of the elements
// and the transpose flag of the vector can be specified via the two template parameters:

   \code
   template< typename Type, bool TF >
   class CompressedVector;
   \endcode

//  - Type: specifies the type of the vector elements. CompressedVector can be used with any
//          non-cv-qualified, non-reference, non-pointer element type.
//  - TF  : specifies whether the vector is a row vector (\c blaze::rowVector) or a column
//          vector (\c blaze::columnVector). The default value is \c blaze::columnVector.
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
// As already mentioned, vectors can be either column vectors (blaze::columnVector) or row vectors
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
// statically allocated elements of arbitrary type. The type of the elements, the number of
// rows and columns, and the storage order of the matrix can be specified via the four template
// parameters:

   \code
   template< typename Type, size_t M, size_t N, bool SO >
   class StaticMatrix;
   \endcode

//  - Type: specifies the type of the matrix elements. StaticMatrix can be used with any
//          non-cv-qualified, non-reference element type.
//  - M   : specifies the total number of rows of the matrix.
//  - N   : specifies the total number of columns of the matrix. Note that it is expected
//          that StaticMatrix is only used for tiny and small matrices.
//  - SO  : specifies the storage order (blaze::rowMajor, blaze::columnMajor) of the matrix.
//          The default value is blaze::rowMajor.
//
//
// \n \section matrix_types_dynamic_matrix DynamicMatrix
// <hr>
//
// The blaze::DynamicMatrix class template is the representation of an arbitrary sized matrix
// with \f$ M \cdot N \f$ dynamically allocated elements of arbitrary type. The type of the
// elements and the storage order of the matrix can be specified via the two template parameters:

   \code
   template< typename Type, bool SO >
   class DynamicMatrix;
   \endcode

//  - Type: specifies the type of the matrix elements. DynamicMatrix can be used with any
//          non-cv-qualified, non-reference element type.
//  - SO  : specifies the storage order (blaze::rowMajor, blaze::columnMajor) of the matrix.
//          The default value is blaze::rowMajor.
//
//
// \n \section matrix_types_compressed_matrix CompressedMatrix
// <hr>
//
// The blaze::CompressedMatrix class template is the representation of an arbitrary sized sparse
// matrix with \f$ M \cdot N \f$ dynamically allocated elements of arbitrary type. The type of the
// elements and the storage order of the matrix can be specified via the two template parameters:

   \code
   template< typename Type, bool SO >
   class CompressedMatrix;
   \endcode

//  - Type: specifies the type of the matrix elements. CompressedMatrix can be used with
//          any non-cv-qualified, non-reference, non-pointer element type.
//  - SO  : specifies the storage order (blaze::rowMajor, blaze::columnMajor) of the matrix.
//          The default value is blaze::rowMajor.
//
//
// \n <center> Previous: \ref vector_operations &nbsp; &nbsp; Next: \ref matrix_operations </center>
*/
//*************************************************************************************************


//**Matrix Operations******************************************************************************
/*!\page matrix_operations Matrix Operations
//
// <center> Previous: \ref matrix_types &nbsp; &nbsp; Next: \ref view_types </center>
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
   blaze::DynamicMatrix<double,3UL,1UL> M3;

   int    array1[4] = { 1, 2, 3, 4 };
   double array2[3] = { 3.1, 6.4, -0.9 };

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

// Note that resizing a matrix invalidates all existing views (see \ref view_types) on the matrix:

   \code
   typedef blaze::DynamicMatrix<int,rowMajor>  MatrixType;

   MatrixType M1( 10UL, 20UL );                 // Creating a 10x20 matrix
   DenseRow<MatrixType> row8 = row( M1, 8UL );  // Creating a view on the 8th row of the matrix
   M1.resize( 6UL, 20UL );                      // Resizing the matrix invalidates the view
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

// \n \section matrix_operations_swap Swap
//
// Via the \c \c swap() function it is possible to completely swap the contents of two matrices
// of the same type:

   \code
   blaze::DynamicMatrix<int,rowMajor> M1( 10UL, 15UL );
   blaze::DynamicMatrix<int,rowMajor> M2( 20UL, 10UL );

   swap( M1, M2 );  // Swapping the contents of M1 and M2
   \endcode

// \n <center> Previous: \ref matrix_types &nbsp; &nbsp; Next: \ref view_types </center>
*/
//*************************************************************************************************


//**View Types*************************************************************************************
/*!\page view_types View Types
//
// <center> Previous: \ref matrix_operations &nbsp; &nbsp; Next: \ref view_operations </center> \n
//
//
// \tableofcontents
//
//
// Views are a very powerful feature to select a specific row or column of a matrix. The Blaze
// library currently offers two different views on dense matrices (\ref view_types_dense_row and
// \ref view_types_dense_column) and two views on sparse matrices (\ref view_types_sparse_row and
// \ref view_types_sparse_column).
//
//
// \n \section view_types_dense_row DenseRow
// <hr>
//
// The blaze::DenseRow class template represents a reference to a specific row of a dense matrix
// primitive. The type of the dense matrix is specified via template parameter:

   \code
   template< typename MT >
   class DenseRow;
   \endcode

//  - MT: specifies the type of the dense matrix primitive. DenseRow can be used with any dense
//        matrix primitive, but does not work with any matrix expression type.
//
//
// \n \section view_types_dense_column DenseColumn
// <hr>
//
// The blaze::DenseColumn class template represents a reference to a specific column of a dense
// matrix primitive. The type of the dense matrix is specified via template parameter:

   \code
   template< typename MT >
   class DenseColumn;
   \endcode

//  - MT: specifies the type of the dense matrix primitive. DenseColumn can be used with any
//        dense matrix primitive, but does not work with any matrix expression type.
//
//
// \n \section view_types_sparse_row SparseRow
// <hr>
//
// The blaze::SparseRow class template represents a reference to a specific row of a sparse matrix
// primitive. The type of the sparse matrix is specified via template parameter:

   \code
   template< typename MT >
   class SparseRow;
   \endcode

//  - MT: specifies the type of the sparse matrix primitive. SparseRow can be used with any sparse
//        matrix primitive, but does not work with any matrix expression type.
//
//
// \n \section view_types_sparse_column SparseColumn
// <hr>
//
// The blaze::SparseColumn template represents a reference to a specific column of a sparse
// matrix primitive. The type of the sparse matrix is specified via template parameter:

   \code
   template< typename MT >
   class SparseColumn;
   \endcode

//  - MT: specifies the type of the sparse matrix primitive. SparseColumn can be used with any
//        sparse matrix primitive, but does not work with any matrix expression type.
//
//
// \n <center> Previous: \ref matrix_operations &nbsp; &nbsp; Next: \ref view_operations </center>
*/
//*************************************************************************************************


//**View Operations********************************************************************************
/*!\page view_operations View Operations
//
// <center> Previous: \ref view_types &nbsp; &nbsp; Next: \ref addition </center>
//
//
// \tableofcontents
//
//
// \n \section view_operations_row_setup Setup of Rows
// <hr>
//
// A reference to a dense or sparse row can very conveniently be created via the \c row() function.
// This reference can be treated as any other row vector, i.e. it can be assigned to, it can be
// copied from, and it can be used in arithmetic operations. The reference can also be used on
// both sides of an assignment: The row can be either used as an alias to grant write access to a
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

// \n \section view_operations_column_setup Setup of Columns
// <hr>
//
// Similar to the setup of a row, a reference to a dense or sparse column can very conveniently
// be created via the \c column() function. This reference can be treated as any other column
// vector, i.e. it can be assigned to, copied from, and be used in arithmetic operations. The
// column can be either used as an alias to grant write access to a specific column of a matrix
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

// \n \section view_operations_common_operations Common Operations
// <hr>
//
// A row view can be used like any other row vector and a column view can be used like any other
// column vector. For instance, the current number of elements can be obtained via the \c size()
// function, the current capacity via the \c capacity() function, and the number of non-zero
// elements via the \c nonZeros() function. However, since rows and columns are references to
// specific rows and columns of a matrix, several operations are not possible on views, such as
// resizing and swapping. The following example shows this by means of a row view:

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

// \n \section view_operations_element_access Element Access
// <hr>
//
// The elements of the row and column can be directly accessed with the subscript operator. The
// numbering of the row/column elements is

                             \f[\left(\begin{array}{*{5}{c}}
                             0 & 1 & 2 & \cdots & N-1 \\
                             \end{array}\right),\f]

// where N is the number of columns/rows of the referenced matrix. Alternatively, the elements of
// a row or column can be traversed via iterators. Just as with vectors, in case of non-const rows
// or columns, \c begin() and \c end() return an Iterator, which allows a manipulation of the
// non-zero value, in case of a constant rows or columns a ConstIterator is returned:

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

// \n \section view_operations_non_fitting_storage_order Views on Matrices with Non-Fitting Storage Order
// <hr>
//
// Especially noteworthy is that row and column views can be created for both row-major and
// column-major matrices. Whereas the interface of a row-major matrix only allows to traverse
// a row directly and the interface of a column-major matrix only allows to traverse a column,
// via views it is possible to traverse a row of a column-major matrix or a column of a row-major
// matrix. For instance:

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
// or a column-view on a row-major matrix can result in a considerable performance decrease in
// comparison to a view on a matrix with a fitting storage orientation. This is due to the
// non-contiguous storage of the matrix elements. Therefore care has to be taken in the choice
// of the most suitable storage order:

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

// Although Blaze performs the resulting vector/matrix multiplication as efficiently as possible
// using a row-major storage order for matrix A would result in a more efficient evaluation.
//
//
// <center> Previous: \ref view_types &nbsp; &nbsp; Next: \ref addition </center>
*/
//*************************************************************************************************


//**Addition***************************************************************************************
/*!\page addition Addition
//
// <center> Previous: \ref view_operations &nbsp; &nbsp; Next: \ref subtraction </center> \n
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

// \n <center> Previous: \ref view_operations &nbsp; &nbsp; Next: \ref subtraction </center>
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
// <center> Previous: \ref matrix_vector_multiplication &nbsp; &nbsp; Next: \ref vector_serialization </center> \n
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
// <center> Previous: \ref matrix_vector_multiplication &nbsp; &nbsp; Next: \ref vector_serialization </center>
*/
//*************************************************************************************************


//**Vector serialization***************************************************************************
/*!\page vector_serialization Vector Serialization
//
// <center> Previous: \ref matrix_matrix_multiplication &nbsp; &nbsp; Next: \ref matrix_serialization </center> \n
//
// Sometimes it is necessary to store vector and/or matrices on disk, for instance for storing
// results or for sharing specific setups with other people. The \b Blaze math serialization
// module provides the according functionality to create platform independent, portable, binary
// representations of vectors and matrices that can be used to store the \a Blaze data structures
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
      blaze::Archive<std::ofstream> archive( "vectors.blaze" );

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
      blaze::Archive<std::ofstream> archive( "vector.blaze" );

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
// <center> Previous: \ref matrix_matrix_multiplication &nbsp; &nbsp; Next: \ref matrix_serialization </center>
*/
//*************************************************************************************************


//**Matrix serialization***************************************************************************
/*!\page matrix_serialization Matrix Serialization
//
// <center> Previous: \ref vector_serialization </center> \n
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
      blaze::Archive<std::ofstream> archive( "matrices.blaze" );

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
      blaze::Archive<std::ofstream> archive( "matrix.blaze" );

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
// <center> Previous: \ref vector_serialization </center>
*/
//*************************************************************************************************

#endif
