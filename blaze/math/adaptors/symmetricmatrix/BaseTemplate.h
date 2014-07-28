//=================================================================================================
/*!
//  \file blaze/math/adaptors/symmetricmatrix/BaseTemplate.h
//  \brief Header file for the implementation of the base template of the SymmetricMatrix
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

#ifndef _BLAZE_MATH_ADAPTORS_SYMMETRICMATRIX_BASETEMPLATE_H_
#define _BLAZE_MATH_ADAPTORS_SYMMETRICMATRIX_BASETEMPLATE_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/typetraits/IsColumnMajorMatrix.h>
#include <blaze/math/typetraits/IsDenseMatrix.h>
#include <blaze/util/typetraits/IsNumeric.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\defgroup symmetric_matrix SymmetricMatrix
// \ingroup adaptors
*/
/*!\brief Matrix adapter for symmetric \f$ N \times N \f$ matrices.
// \ingroup symmetric_matrix
//
// \section symmetricmatrix_general General
//
// The SymmetricMatrix class template is an adapter for existing dense and sparse matrix types.
// It inherits the properties and the interface of the given matrix type \a MT and extends it
// by enforcing the additional invariant of symmetry (i.e. the matrix is always equal to its
// transpose \f$ A = A^T \f$). The type of the adapted matrix can be specified via the first
// template parameter:

   \code
   template< typename MT, bool DF, bool NF >
   class SymmetricMatrix;
   \endcode

//  - MT: specifies the type of the matrix to be adapted. SymmetricMatrix can be used with any
//        non-cv-qualified, non-reference, non-pointer, non-expression dense or sparse matrix
//        type. Note that the given matrix type must be either resizable (as for instance
//        HybridMatrix or DynamicMatrix) or must be square at compile time (as for instance
//        StaticMatrix).
//  - DF: specifies whether the given matrix type is a dense or sparse matrix type. This template
//        parameter doesn't have to be defined explicitly, it is automatically derived from the
//        first template parameter. Defining the parameter explicitly may result in a compilation
//        error!
//  - NF: determines how the elements of the matrix are handled internally. This template parameter
//        must \b NOT be defined explicitly, it is automatically derived from the first template
//        parameter. Defining the parameter explicitly may result in a compilation error!
//
// The following example give an impression of several possible symmetric matrices:

   \code
   // Definition of a 3x3 row-major dense symmetric matrix with static memory
   blaze::SymmetricMatrix< blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> > A;

   // Definition of a resizable column-major dense symmetric matrix based on HybridMatrix
   blaze::SymmetricMatrix< blaze::HybridMatrix<float,4UL,4UL,blaze::columnMajor> B;

   // Definition of a resizable row-major dense symmetric matrix based on DynamicMatrix
   blaze::SymmetricMatrix< blaze::DynamicMatrix<double,blaze::rowMajor> > C;

   // Definition of a compressed row-major single precision symmetric matrix
   blaze::SymmetricMatrix< blaze::CompressedMatrix<float,blaze::rowMajor> > D;
   \endcode

// The storage order of a symmetric matrix is depending on the storage order of the adapted matrix
// type \a MT. In case the adapted matrix is stored in a row-wise fashion (i.e. is specified as
// blaze::rowMajor), the symmetric matrix will also be a row-major matrix. Otherwise, if the
// adapted matrix is column-major (i.e. is specified as blaze::columnMajor), the symmetric matrix
// will also be a column-major matrix.
//
//
// \n \section symmetricmatrix_special_properties Special Properties of Symmetric Matrices
//
// A symmetric matrix is used exactly like a matrix of the underlying, adapted matrix type \a MT.
// It also provides (nearly) the same interface as the underlying matrix type. However, there are
// some important exceptions resulting from the symmetry constraint:
//
//  -# <b>\ref symmetricmatrix_square</b>
//  -# <b>\ref symmetricmatrix_symmetry</b>
//  -# <b>\ref symmetricmatrix_initialization</b>
//
// \n \subsection symmetricmatrix_square Symmetric Matrices Must Always be Square!
//
// In case a resizable matrix is used (as for instance blaze::HybridMatrix, blaze::DynamicMatrix,
// or blaze::CompressedMatrix), this means that the according constructors, the \c resize() and
// the \c extend() functions only expect a single parameter, which specifies both the number of
// rows and columns, instead of two (one for the number of rows and one for the number of columns):

   \code
   using blaze::DynamicMatrix;
   using blaze::SymmetricMatrix;
   using blaze::rowMajor;

   // Default constructed, default initialized, row-major 3x3 symmetric dynamic matrix
   SymmetricMatrix< DynamicMatrix<double,rowMajor> > A( 3 );

   // Resizing the matrix to 5x5
   A.resize( 5 );

   // Extending the number of rows and columns by 2, resulting in a 7x7 matrix
   A.extend( 2 );
   \endcode

// In case a matrix with a fixed size is used (as for instance blaze::StaticMatrix), the number
// of rows and number of columns must be specified equally:

   \code
   using blaze::StaticMatrix;
   using blaze::SymmetricMatrix;
   using blaze::columnMajor;

   // Correct setup of a fixed size column-major 3x3 symmetric static matrix
   SymmetricMatrix< StaticMatrix<int,3UL,3UL,columnMajor> > A;

   // Compilation error: the provided matrix type is not a square matrix type
   SymmetricMatrix< StaticMatrix<int,3UL,4UL,columnMajor> > B;
   \endcode

// \n \subsection symmetricmatrix_symmetry The Symmetric Property is Always Enforced!
//
// This means that modifying the element \f$ a_{ij} \f$ of a symmetric matrix also modifies its
// counterpart element \f$ a_{ji} \f$. Also, it is only possible to assign matrices that are
// symmetric themselves:

   \code
   using blaze::CompressedMatrix;
   using blaze::DynamicMatrix;
   using blaze::StaticMatrix;
   using blaze::SymmetricMatrix;
   using blaze::rowMajor;

   // Default constructed, row-major 3x3 symmetric compressed matrix
   SymmetricMatrix< CompressedMatrix<double,rowMajor> > A( 3 );

   // Initializing three elements via the function call operator
   A(0,0) = 1.0;  // Initialization of the diagonal element (0,0)
   A(0,2) = 2.0;  // Initialization of the elements (0,2) and (2,0)

   // Inserting three more elements via the insert() function
   A.insert( 1, 1, 3.0 );  // Inserting the diagonal element (1,1)
   A.insert( 1, 2, 4.0 );  // Inserting the elements (1,2) and (2,1)

   // Access via a non-const iterator
   *A.begin(1UL) = 10.0;  // Modifies both elements (1,0) and (0,1)

   // Erasing elements via the erase() function
   A.erase( 0, 0 );  // Erasing the diagonal element (0,0)
   A.erase( 0, 2 );  // Erasing the elements (0,2) and (2,0)

   // Construction from a symmetric dense matrix
   StaticMatrix<double,3UL,3UL> B(  3.0,  8.0, -2.0,
                                    8.0,  0.0, -1.0,
                                   -2.0, -1.0,  4.0 );

   SymmetricMatrix< DynamicMatrix<double,rowMajor> > C( B );  // OK

   // Assignment of a non-symmetric dense matrix
   StaticMatrix<double,3UL,3UL> D(  3.0,  8.0, -2.0,
                                    8.0,  0.0, -1.0,
                                   -2.0, -1.0,  4.0 );

   C = D;  // Throws an exception; symmetric invariant would be violated!
   \endcode

// The symmetry property is also enforced for views (rows, columns, submatrices, ...) on the
// symmetric matrix. The following example demonstrates that modifying the elements of an entire
// row of the symmetric matrix also affects the counterpart elements in the according column of
// the matrix:

   \code
   using blaze::DynamicMatrix;
   using blaze::SymmetricMatrix;

   // Setup of the symmetric matrix
   //
   //       ( 0 1 0 2 )
   //   A = ( 1 3 4 0 )
   //       ( 0 4 0 5 )
   //       ( 2 0 5 0 )
   //
   SymmetricMatrix< DynamicMatrix<int> > A( 4 );
   A(0,1) = 1;
   A(0,3) = 2;
   A(1,1) = 3;
   A(1,2) = 4;
   A(2,3) = 5;

   // Setting all elements in the 1st row to 0 results in the matrix
   //
   //       ( 0 0 0 2 )
   //   A = ( 0 0 0 0 )
   //       ( 0 0 0 5 )
   //       ( 2 0 5 0 )
   //
   row( A, 1 ) = 0;
   \endcode

// The same restriction also applies to the \c append() function for sparse matrices: Appending
// the element \f$ a_{ij} \f$ additionally inserts the element \f$ a_{ji} \f$ into the matrix.
// Despite the additional insertion, the \c append() function still provides the most efficient
// way to set up a symmetric sparse matrix. In order to achieve the maximum efficiency, the
// capacity of the individual rows/columns of the matrix should to be specifically prepared with
// \c reserve() calls:

   \code
   using blaze::CompressedMatrix;
   using blaze::SymmetricMatrix;
   using blaze::rowMajor;

   // Setup of the symmetric matrix
   //
   //       ( 0 1 3 )
   //   A = ( 1 2 0 )
   //       ( 3 0 0 )

   SymmetricMatrix< CompressedMatrix<double,rowMajor> > A( 3 );

   A.reserve( 5 );         // Reserving enough space for 5 non-zero elements
   A.reserve( 0, 2 );      // Reserving two non-zero elements in the first row
   A.reserve( 1, 2 );      // Reserving two non-zero elements in the second row
   A.reserve( 2, 1 );      // Reserving a single non-zero element in the third row
   A.append( 0, 1, 1.0 );  // Appending the value 1 at position (0,1) and (1,0)
   A.append( 1, 1, 2.0 );  // Appending the value 2 at position (1,1)
   A.append( 2, 0, 3.0 );  // Appending the value 3 at position (2,0) and (0,2)
   \endcode

// \n \subsection symmetricmatrix_initialization The Elements of a Dense Symmetric Matrix are Always Default Initialized!
//
// Although this results in a small loss of efficiency (especially in case all default values are
// overridden afterwards), this property is important since otherwise the symmetric property of
// dense symmetric matrices could not be guaranteed:

   \code
   using blaze::DynamicMatrix;
   using blaze::SymmetricMatrix;

   // Uninitialized, 5x5 row-major dynamic matrix
   DynamicMatrix<int,rowMajor> A( 5, 5 );

   // Default initialized, 5x5 row-major symmetric dynamic matrix
   SymmetricMatrix< DynamicMatrix<int,rowMajor> > B( 5 );
   \endcode

// \n \section symmetricmatrix_arithmetic_operations Arithmetic Operations
//
// A SymmetricMatrix matrix can participate in numerical operations in any way any other dense
// or sparse matrix can participate. It can also be combined with any other dense or sparse vector
// or matrix. The following code example gives an impression of the use of SymmetricMatrix within
// arithmetic operations:

   \code
   using blaze::SymmetricMatrix;
   using blaze::DynamicMatrix;
   using blaze::StaticMatrix;
   using blaze::CompressedMatrix;
   using blaze::rowMajor;
   using blaze::columnMajor;

   CompressedMatrix<float> E( 3, 3 );  // Empty row-major sparse single precision 3x3 matrix

   SymmetricMatrix< HybridMatrix<float,3UL,3UL,rowMajor> > F;
   SymmetricMatrix< StaticMatrix<float,3UL,3UL,columnMajor> > G;

   F = A + B;     // Matrix addition and assignment to a row-major symmetric matrix
   G = A - C;     // Matrix subtraction and assignment to a column-major symmetric matrix
   G = A * E;     // Matrix multiplication between a dense and a sparse matrix

   A *= 2.0;      // In-place scaling of matrix A
   F  = 2.0 * B;  // Scaling of matrix B
   G  = E * 2.0;  // Scaling of matrix E

   F += A - B;    // Addition assignment
   G -= A + C;    // Subtraction assignment
   G *= A * E;    // Multiplication assignment
   \endcode

// Note that it is possible to use have block-structured symmetric matrices:

   \code
   // Definition of a block-structured row-major dense symmetric matrix based on CompressedMatrix
   blaze::SymmetricMatrix< blaze::CompressedMatrix< blaze::StaticMatrix<double,3UL> > > A;
   \endcode

// Also in this case, the SymmetricMatrix class template enforces the invariant of symmetry and
// guarantees that a modifications of element \f$ a_{ij} \f$ of the adapted matrix is also
// applied to element \f$ a_{ji} \f$.
*/
template< typename MT                                             // Type of the adapted matrix
        , bool DF = IsDenseMatrix<MT>::value                      // Density flag
        , bool NF = IsNumeric<typename MT::ElementType>::value >  // Numeric flag
class SymmetricMatrix;
//*************************************************************************************************

} // namespace blaze

#endif
