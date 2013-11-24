//=================================================================================================
/*!
//  \file blazetest/config/MathTest.h
//  \brief General configuration file for the math tests of the blaze test suite
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


//*************************************************************************************************
/*!\brief First data type for the math tests.
//
// The math tests are always conducted with two data types to test mixed type arithmetic
// expressions. This type definition defines the first of these two data types.
*/
typedef int  TypeA;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Second data type for the math tests.
//
// The math tests are always conducted with two data types to test mixed type arithmetic
// expressions. This type definition defines the second of these two data types.
*/
typedef double  TypeB;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compilation switch for the basic tests.
//
// This compilation switch triggers the basic tests for all test scenarios. In case the
// basic tests are activated, each operation is tested separately and without any other
// combined operation. The following example demonstrates this by means of the vector
// addition:

   \code
   blaze::DynamicVector<double> a, b, c;
   c = a + b;  // Basic vector addition
   \endcode

// The following settings are possible:
//
//   - 0: The basic tests are not included in the compilation process and not executed
//   - 1: The basic tests are included in the compilation process, but not executed
//   - 2: The basic tests are included in the compilation process and executed
*/
#define BLAZETEST_MATHTEST_TEST_BASIC_OPERATION 2
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compilation switch for the negation tests.
//
// This compilation switch triggers the negation tests for all test scenarios. In case the
// negation tests are activated, each operation is tested in combination with a negation.
// The following example demonstrates this by means of the vector addition:

   \code
   blaze::DynamicVector<double> a, b, c;
   c = -( a + b );  // Negated vector addition
   \endcode

// The following settings are possible:
//
//   - 0: The negation tests are not included in the compilation process and not executed
//   - 1: The negation tests are included in the compilation process, but not executed
//   - 2: The negation tests are included in the compilation process and executed
*/
#define BLAZETEST_MATHTEST_TEST_NEGATED_OPERATION 2
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compilation switch for the scaling tests.
//
// This compilation switch triggers the scaling tests for all test scenarios. In case the
// scaling tests are activated, each operation is tested in combination with a scalar
// multiplication or division. The following example demonstrates this by means of the
// vector addition:

   \code
   blaze::DynamicVector<double> a, b, c;
   c = 1.1 * ( a + b );  // Left-multiplied vector addition
   c = ( a + b ) * 1.1;  // Right-multiplied vector addition
   c = ( a + b ) / 1.1;  // Scalar-divided vector addition
   \endcode

// The following settings are possible:
//
//   - 0: The scaling tests are not included in the compilation process and not executed
//   - 1: The scaling tests are included in the compilation process, but not executed
//   - 2: The scaling tests are included in the compilation process and executed
*/
#define BLAZETEST_MATHTEST_TEST_SCALED_OPERATION 2
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compilation switch for the transpose tests.
//
// This compilation switch triggers the transpose tests for all test scenarios. In case the
// transpose tests are activated, each operation is tested in combination with a transpose
// operation. The following example demonstrates this by means of the vector addition:

   \code
   blaze::DynamicVector<double,false> a, b, ;
   blaze::DynamicVector<double,true> c;
   c = trans( a + b );  // Transpose vector addition
   \endcode

// The following settings are possible:
//
//   - 0: The transpose tests are not included in the compilation process and not executed
//   - 1: The transpose tests are included in the compilation process, but not executed
//   - 2: The transpose tests are included in the compilation process and executed
*/
#define BLAZETEST_MATHTEST_TEST_TRANSPOSE_OPERATION 2
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compilation switch for the abs tests.
//
// This compilation switch triggers the abs tests for all test scenarios. In case the abs
// tests are activated, each operation is tested in combination with an abs operation. The
// following example demonstrates this by means of the vector addition:

   \code
   blaze::DynamicVector<double> a, b, c;
   c = abs( a + b );  // Absolute value vector addition
   \endcode

// The following settings are possible:
//
//   - 0: The abs tests are not included in the compilation process and not executed
//   - 1: The abs tests are included in the compilation process, but not executed
//   - 2: The abs tests are included in the compilation process and executed
*/
#define BLAZETEST_MATHTEST_TEST_ABS_OPERATION 2
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compilation switch for the subvector tests.
//
// This compilation switch triggers the subvector tests for all test scenarios. In case the
// subvector tests are activated, all operations resulting in vectors are tested in combination
// with a subvector operation. The following example gives an impression by means of the vector
// addition:

   \code
   blaze::DynamicVector<double> a, b, c;
   for( ... )
      subvector( c, ... ) = subvector( a + b, ... );  // Subvector-wise vector addition
   \endcode

// The following settings are possible:
//
//   - 0: The subvector tests are not included in the compilation process and not executed
//   - 1: The subvector tests are included in the compilation process, but not executed
//   - 2: The subvector tests are included in the compilation process and executed
*/
#define BLAZETEST_MATHTEST_TEST_SUBVECTOR_OPERATION 2
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compilation switch for the submatrix tests.
//
// This compilation switch triggers the submatrix tests for all test scenarios. In case the
// submatrix tests are activated, all operations resulting in matrices are tested in combination
// with a submatrix operation. The following example gives an impression by means of the matrix
// addition:

   \code
   blaze::DynamicMatrix<double> A, B, C;
   for( ... )
      for( ... )
         submatrix( C, ... ) = submatrix( A + B, ... );  // Submatrix-wise matrix addition
   \endcode

// The following settings are possible:
//
//   - 0: The submatrix tests are not included in the compilation process and not executed
//   - 1: The submatrix tests are included in the compilation process, but not executed
//   - 2: The submatrix tests are included in the compilation process and executed
*/
#define BLAZETEST_MATHTEST_TEST_SUBMATRIX_OPERATION 2
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compilation switch for the row tests.
//
// This compilation switch triggers the row tests for all test scenarios. In case the row
// tests are activated, all operations resulting in matrices are tested in combination with
// a row operation. The following example demonstrates this by means of the matrix addition:

   \code
   blaze::DynamicMatrix<double> A, B, C;
   for( size_t i=0UL; i<A.rows(); ++i )
      row( C, i ) = row( A + B, i );  // Row-wise matrix addition
   \endcode

// The following settings are possible:
//
//   - 0: The row tests are not included in the compilation process and not executed
//   - 1: The row tests are included in the compilation process, but not executed
//   - 2: The row tests are included in the compilation process and executed
*/
#define BLAZETEST_MATHTEST_TEST_ROW_OPERATION 2
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compilation switch for the column tests.
//
// This compilation switch triggers the column tests for all test scenarios. In case the column
// tests are activated, all operations resulting in matrices are tested in combination with
// a column operation. The following example demonstrates this by means of the matrix addition:

   \code
   blaze::DynamicMatrix<double> A, B, C;
   for( size_t i=0UL; i<A.column(); ++i )
      column( C, i ) = column( A + B, i );  // Column-wise matrix addition
   \endcode

// The following settings are possible:
//
//   - 0: The column tests are not included in the compilation process and not executed
//   - 1: The column tests are included in the compilation process, but not executed
//   - 2: The column tests are included in the compilation process and executed
*/
#define BLAZETEST_MATHTEST_TEST_COLUMN_OPERATION 2
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Number of repetitions for a single test case.
//
// The \a repetitions value specifies the number of repetitions for each single test case. In
// each repetition the test case is run with random input values. Therefore a higher number of
// repetitions increases the likelihood of detecting implementationn errors. On the downside
// the execution of all tests takes longer.
*/
const size_t repetitions = 3;
//*************************************************************************************************
