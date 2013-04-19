//=================================================================================================
/*!
//  \file blazetest/mathtest/tsvecsmatmult/AliasingTest.h
//  \brief Header file for the sparse vector/sparse matrix multiplication aliasing test
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

#ifndef _BLAZETEST_MATHTEST_TSVECSMATMULT_ALIASINGTEST_H_
#define _BLAZETEST_MATHTEST_TSVECSMATMULT_ALIASINGTEST_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <sstream>
#include <stdexcept>
#include <string>
#include <blaze/math/CompressedMatrix.h>
#include <blaze/math/CompressedVector.h>
#include <blaze/math/DynamicVector.h>
#include <blaze/math/StaticVector.h>


namespace blazetest {

namespace mathtest {

namespace tsvecsmatmult {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Auxiliary class template for the sparse vector/sparse matrix multiplication aliasing test.
//
// This class represents a test suite for all sparse vector/sparse matrix multiplication aliasing
// tests. It performs a series of runtime tests to assure that all mathematical operations work
// correctly even in the presence of aliasing.
*/
class AliasingTest
{
 private:
   //**Type definitions****************************************************************************
   typedef blaze::DynamicVector<int,blaze::rowVector>       TDVec;  //!< Dense row vector type.
   typedef blaze::CompressedVector<int,blaze::rowVector>    TSVec;  //!< Sparse row vector type.
   typedef blaze::CompressedMatrix<int,blaze::rowMajor>     SMat;   //!< Row-major sparse matrix type.
   typedef blaze::CompressedMatrix<int,blaze::columnMajor>  TSMat;  //!< Column-major sparse matrix type.
   typedef blaze::StaticVector<int,3UL,blaze::rowVector>    TRVec;  //!< Result row vector type.
   //**********************************************************************************************

 public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit AliasingTest();
   // No explicitly declared copy constructor.
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   // No explicitly declared destructor.
   //**********************************************************************************************

 private:
   //**Test functions******************************************************************************
   /*!\name Test functions */
   //@{
   void testTSVecSMatMult ();
   void testTSVecTSMatMult();

   template< typename T1, typename T2 >
   void checkResult( const T1& computedResult, const T2& expectedResult );
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   void initialize();
   //@}
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   SMat sA4x3_;    //!< The first row-major sparse matrix.
                   /*!< The \f$ 3 \times 4 \f$ matrix is initialized as
                        \f[\left(\begin{array}{*{3}{c}}
                        -1 & 0 & -2 \\
                         0 & 2 & -3 \\
                         0 & 1 &  2 \\
                         1 & 0 & -2 \\
                        \end{array}\right)\f]. */
   SMat sB3x3_;    //!< The second row-major sparse matrix.
                   /*!< The \f$ 3 \times 3 \f$ matrix is initialized as
                        \f[\left(\begin{array}{*{3}{c}}
                        0 & -1 &  0 \\
                        1 & -2 &  2 \\
                        0 &  0 & -3 \\
                        \end{array}\right)\f]. */
   TSMat tsA4x3_;  //!< The first column-major sparse matrix.
                   /*!< The \f$ 4 \times 3 \f$ matrix is initialized as
                        \f[\left(\begin{array}{*{3}{c}}
                        -1 & 0 & -2 \\
                         0 & 2 & -3 \\
                         0 & 1 &  2 \\
                         1 & 0 & -2 \\
                        \end{array}\right)\f]. */
   TSMat tsB3x3_;  //!< The second column-major sparse matrix.
                   /*!< The \f$ 3 \times 3 \f$ matrix is initialized as
                        \f[\left(\begin{array}{*{3}{c}}
                        0 & -1 &  0 \\
                        1 & -2 &  2 \\
                        0 &  0 & -3 \\
                        \end{array}\right)\f]. */
   TSVec tsa4_;    //!< The first sparse row vector.
                   /*!< The 4-dimensional vector is initialized as
                        \f[\left(\begin{array}{*{1}{c}}
                        -1 \\
                         0 \\
                        -3 \\
                         2 \\
                        \end{array}\right)\f]. */
   TSVec tsb4_;    //!< The second sparse row vector.
                   /*!< The 4-dimensional vector is initialized as
                        \f[\left(\begin{array}{*{1}{c}}
                         0 \\
                         1 \\
                         2 \\
                        -1 \\
                        \end{array}\right)\f]. */
   TSVec tsc3_;    //!< The third sparse row vector.
                   /*!< The 3-dimensional vector is initialized as
                        \f[\left(\begin{array}{*{1}{c}}
                        1 \\
                        2 \\
                        3 \\
                        \end{array}\right)\f]. */
   TSVec tsd3_;    //!< The fourth sparse row vector.
                   /*!< The 3-dimensional vector is initialized as
                        \f[\left(\begin{array}{*{1}{c}}
                        0 \\
                        2 \\
                        1 \\
                        \end{array}\right)\f]. */
   TDVec tda4_;    //!< The first dense row vector.
                   /*!< The 4-dimensional vector is initialized as
                        \f[\left(\begin{array}{*{1}{c}}
                        -1 \\
                         0 \\
                        -3 \\
                         2 \\
                        \end{array}\right)\f]. */
   TDVec tdb3_;    //!< The second dense row vector.
                   /*!< The 3-dimensional vector is initialized as
                        \f[\left(\begin{array}{*{1}{c}}
                        0 \\
                        2 \\
                        1 \\
                        \end{array}\right)\f]. */
   TRVec result_;  //!< The dense vector for the reference result.

   std::string test_;  //!< Label of the currently performed test.
   //@}
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Checking and comparing the computed result.
//
// \param computedResult The computed result.
// \param expectedResult The expected result.
// \return void
// \exception std::runtime_error Incorrect result detected.
//
// This function is called after each test case to check and compare the computed result.
// In case the computed and the expected result differ in any way, a \a std::runtime_error
// exception is thrown.
*/
template< typename T1    // Vector type of the computed result
        , typename T2 >  // Vector type of the expected result
void AliasingTest::checkResult( const T1& computedResult, const T2& expectedResult )
{
   if( computedResult != expectedResult ) {
      std::ostringstream oss;
      oss.precision( 20 );
      oss << " Test : " << test_ << "\n"
          << " Error: Incorrect result detected\n"
          << " Details:\n"
          << "   Computed result:\n" << computedResult << "\n"
          << "   Expected result:\n" << expectedResult << "\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Testing the sparse vector/sparse matrix multiplication in the presence of aliasing.
//
// \return void
*/
void runTest()
{
   AliasingTest();
}
//*************************************************************************************************




//=================================================================================================
//
//  MACRO DEFINITIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Macro for the execution of the sparse vector/sparse matrix multiplication aliasing test.
*/
#define RUN_TSVECSMATMULT_ALIASING_TEST \
   blazetest::mathtest::tsvecsmatmult::runTest()
/*! \endcond */
//*************************************************************************************************

} // namespace tsvecsmatmult

} // namespace mathtest

} // namespace blazetest

#endif
