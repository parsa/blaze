//=================================================================================================
/*!
//  \file blazetest/mathtest/dmatdmatadd/OperationTest.h
//  \brief Header file for the dense matrix/dense matrix addition operation test
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

#ifndef _BLAZETEST_MATHTEST_DMATDMATADD_OPERATIONTEST_H_
#define _BLAZETEST_MATHTEST_DMATDMATADD_OPERATIONTEST_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <sstream>
#include <stdexcept>
#include <string>
#include <typeinfo>
#include <blaze/math/CompressedMatrix.h>
#include <blaze/math/constraints/Computation.h>
#include <blaze/math/constraints/DenseMatrix.h>
#include <blaze/math/constraints/MatMatAddExpr.h>
#include <blaze/math/constraints/SparseMatrix.h>
#include <blaze/math/constraints/StorageOrder.h>
#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/shims/Equal.h>
#include <blaze/math/StaticMatrix.h>
#include <blaze/math/traits/AddExprTrait.h>
#include <blaze/math/traits/AddTrait.h>
#include <blaze/math/typetraits/BaseElementType.h>
#include <blaze/math/typetraits/IsRowMajorMatrix.h>
#include <blaze/math/Views.h>
#include <blaze/util/constraints/Numeric.h>
#include <blaze/util/constraints/SameType.h>
#include <blazetest/system/MathTest.h>
#include <blazetest/mathtest/Creator.h>
#include <blazetest/mathtest/IsEqual.h>
#include <blazetest/mathtest/RandomMaximum.h>
#include <blazetest/mathtest/RandomMinimum.h>


namespace blazetest {

namespace mathtest {

namespace dmatdmatadd {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Auxiliary class template for the dense matrix/dense matrix addition operation test.
//
// This class template represents one particular matrix addition test between two matrices of
// a particular type. The two template arguments \a MT1 and \a MT2 represent the types of the
// left-hand side and right-hand side matrix, respectively.
*/
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
class OperationTest
{
 private:
   //**Type definitions****************************************************************************
   typedef typename MT1::OppositeType               OMT1;  //!< Matrix type 1 with opposite storage order
   typedef typename MT2::OppositeType               OMT2;  //!< Matrix type 2 with opposite storage order
   typedef typename MT1::TransposeType              TMT1;  //!< Transpose matrix type 1
   typedef typename MT2::TransposeType              TMT2;  //!< Transpose matrix type 2
   typedef typename blaze::AddTrait<MT1,MT2>::Type  RE;    //!< Default result type
   typedef typename RE::OppositeType                ORE;   //!< Default result type with opposite storage order
   typedef typename RE::TransposeType               TRE;   //!< Transpose default result type
   typedef typename ORE::TransposeType              TORE;  //!< Transpose default result type with opposite storage order
   //**********************************************************************************************

   //**Type definitions****************************************************************************
   typedef typename MT1::ElementType           ET1;     //!< Element type 1
   typedef typename MT2::ElementType           ET2;     //!< Element type 2
   typedef typename RE::ElementType            RET;     //!< Resulting element type
   typedef blaze::DynamicMatrix<ET1,false>     RT1;     //!< Reference type 1
   typedef blaze::CompressedMatrix<ET2,false>  RT2;     //!< Reference type 2
   typedef blaze::DynamicMatrix<RET,false>     DRRE;    //!< Dense reference result type
   typedef blaze::CompressedMatrix<RET,false>  SRRE;    //!< Sparse reference result type
   typedef typename DRRE::OppositeType         ODRRE;   //!< Dense reference result type with opposite storage order
   typedef typename SRRE::OppositeType         OSRRE;   //!< Sparse reference result type with opposite storage order
   typedef typename DRRE::TransposeType        TDRRE;   //!< Transpose dense reference result type
   typedef typename SRRE::TransposeType        TSRRE;   //!< Transpose sparse reference result type
   typedef typename ODRRE::TransposeType       TODRRE;  //!< Transpose dense reference result type with opposite storage order
   typedef typename OSRRE::TransposeType       TOSRRE;  //!< Transpose sparse reference result type with opposite storage order
   typedef RE                                  DRE;     //!< Dense result type
   typedef SRRE                                SRE;     //!< Sparse result type
   typedef ORE                                 ODRE;    //!< Dense result type with opposite storage order
   typedef OSRRE                               OSRE;    //!< Sparse result type with opposite storage order
   typedef TRE                                 TDRE;    //!< Transpose dense result type
   typedef TSRRE                               TSRE;    //!< Transpose sparse result type
   typedef TORE                                TODRE;   //!< Transpose dense result type with opposite storage order
   typedef TOSRRE                              TOSRE;   //!< Transpose sparse result type with opposite storage order

   //! Type of the matrix/matrix addition expression
   typedef typename blaze::AddExprTrait<MT1,MT2>::Type  MatMatAddExprType;

   //! Type of the matrix/transpose matrix addition expression
   typedef typename blaze::AddExprTrait<MT1,OMT2>::Type  MatTMatAddExprType;

   //! Type of the transpose matrix/matrix addition expression
   typedef typename blaze::AddExprTrait<OMT1,MT2>::Type  TMatMatAddExprType;

   //! Type of the transpose matrix/transpose matrix addition expression
   typedef typename blaze::AddExprTrait<OMT1,OMT2>::Type  TMatTMatAddExprType;
   //**********************************************************************************************

 public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit OperationTest( const Creator<MT1>& creator1, const Creator<MT2>& creator2 );
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
                          void testInitialStatus     ();
                          void testAssignment        ();
                          void testElementAccess     ();
                          void testBasicOperation    ();
                          void testNegatedOperation  ();
   template< typename T > void testScaledOperation   ( T scalar );
                          void testTransposeOperation();
                          void testAbsOperation      ();
                          void testSubmatrixOperation();
                          void testRowOperation      ();
                          void testColumnOperation   ();
   //@}
   //**********************************************************************************************

   //**Error detection functions*******************************************************************
   /*!\name Error detection functions */
   //@{
   template< typename LT, typename RT > void checkResults();
   template< typename LT, typename RT > void checkTransposeResults();
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   void initResults();
   void initTransposeResults();
   template< typename LT, typename RT > void convertException( const std::exception& ex );
   //@}
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   MT1   lhs_;     //!< The left-hand side dense matrix.
   MT2   rhs_;     //!< The right-hand side dense matrix.
   OMT1  olhs_;    //!< The left-hand side dense matrix with opposite storage order.
   OMT2  orhs_;    //!< The right-hand side dense matrix with opposite storage order.
   DRE   dres_;    //!< The dense result matrix.
   SRE   sres_;    //!< The sparse result matrix.
   ODRE  odres_;   //!< The dense result matrix with opposite storage order.
   OSRE  osres_;   //!< The sparse result matrix with opposite storage order.
   TDRE  tdres_;   //!< The transpose dense result matrix.
   TSRE  tsres_;   //!< The transpose sparse result matrix.
   TODRE todres_;  //!< The transpose dense result matrix with opposite storage order.
   TOSRE tosres_;  //!< The transpose sparse result matrix with opposite storage order.
   RT1   reflhs_;  //!< The reference left-hand side matrix.
   RT2   refrhs_;  //!< The reference right-hand side matrix.
   DRRE  refres_;  //!< The reference result.

   std::string test_;   //!< Label of the currently performed test.
   std::string error_;  //!< Description of the current error type.
   //@}
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE ( MT1    );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE ( MT2    );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE ( OMT1   );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE ( OMT2   );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE ( TMT1   );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE ( TMT2   );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE ( RT1    );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( RT2    );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE ( DRRE   );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( SRRE   );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE ( ODRRE  );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( OSRRE  );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE ( TDRRE  );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( TSRRE  );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE ( TODRRE );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( TOSRRE );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE ( DRE    );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( SRE    );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE ( ODRE   );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( OSRE   );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE ( TDRE   );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( TSRE   );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE ( TODRE  );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( TOSRE  );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE   ( MT1    );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE   ( MT2    );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( OMT1   );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( OMT2   );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( TMT1   );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( TMT1   );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE   ( RT1    );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE   ( RT2    );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE   ( DRRE   );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE   ( SRRE   );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( ODRRE  );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( OSRRE  );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( TDRRE  );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( TSRRE  );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE   ( TODRRE );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE   ( TOSRRE );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE   ( DRE    );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE   ( SRE    );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( ODRE   );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( OSRE   );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( TDRE   );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( TSRE   );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE   ( TODRE  );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE   ( TOSRE  );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ET1, typename OMT1::ElementType   );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ET2, typename OMT2::ElementType   );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ET1, typename TMT1::ElementType   );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ET2, typename TMT2::ElementType   );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MT1, typename OMT1::OppositeType  );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MT2, typename OMT2::OppositeType  );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MT1, typename TMT1::TransposeType );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MT2, typename TMT2::TransposeType );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( RE , typename ORE::OppositeType   );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( RE , typename TRE::TransposeType  );
   BLAZE_CONSTRAINT_MUST_BE_MATMATADDEXPR_TYPE( MatMatAddExprType   );
   BLAZE_CONSTRAINT_MUST_BE_MATMATADDEXPR_TYPE( MatTMatAddExprType  );
   BLAZE_CONSTRAINT_MUST_BE_MATMATADDEXPR_TYPE( TMatMatAddExprType  );
   BLAZE_CONSTRAINT_MUST_BE_MATMATADDEXPR_TYPE( TMatTMatAddExprType );
   BLAZE_CONSTRAINT_MUST_BE_COMPUTATION_TYPE( MatMatAddExprType   );
   BLAZE_CONSTRAINT_MUST_BE_COMPUTATION_TYPE( MatTMatAddExprType  );
   BLAZE_CONSTRAINT_MUST_BE_COMPUTATION_TYPE( TMatMatAddExprType  );
   BLAZE_CONSTRAINT_MUST_BE_COMPUTATION_TYPE( TMatTMatAddExprType );
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the dense matrix/dense matrix addition operation test.
//
// \param creator1 The creator for the left-hand side dense matrix of the matrix addition.
// \param creator2 The creator for the right-hand side dense matrix of the matrix addition.
// \exception std::runtime_error Operation error detected.
*/
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
OperationTest<MT1,MT2>::OperationTest( const Creator<MT1>& creator1, const Creator<MT2>& creator2 )
   : lhs_( creator1() )  // The left-hand side dense matrix
   , rhs_( creator2() )  // The right-hand side dense matrix
   , olhs_( lhs_ )       // The left-hand side dense matrix with opposite storage order
   , orhs_( rhs_ )       // The right-hand side dense matrix with opposite storage order
   , dres_()             // The dense result matrix
   , sres_()             // The sparse result matrix
   , odres_()            // The dense result matrix with opposite storage order
   , osres_()            // The sparse result matrix with opposite storage order
   , tdres_()            // The transpose dense result matrix
   , tsres_()            // The transpose sparse result matrix
   , todres_()           // The transpose dense result matrix with opposite storage order
   , tosres_()           // The transpose sparse result matrix with opposite storage order
   , reflhs_( lhs_ )     // The reference left-hand side matrix
   , refrhs_( rhs_ )     // The reference right-hand side matrix
   , refres_()           // The reference result
   , test_()             // Label of the currently performed test
   , error_()            // Description of the current error type
{
   testInitialStatus();
   testAssignment();
   testElementAccess();
   testBasicOperation();
   testNegatedOperation();
   testScaledOperation( 2 );
   testScaledOperation( 2UL );
   testScaledOperation( 2.0F );
   testScaledOperation( 2.0 );
   testTransposeOperation();
   testAbsOperation();
   testSubmatrixOperation();
   testRowOperation();
   testColumnOperation();
}
//*************************************************************************************************




//=================================================================================================
//
//  TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Tests on the initial status of the matrices.
//
// \return void
// \exception std::runtime_error Initialization error detected.
//
// This function runs tests on the initial status of the matrices. In case any initialization
// error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
void OperationTest<MT1,MT2>::testInitialStatus()
{
   //=====================================================================================
   // Performing initial tests with the row-major types
   //=====================================================================================

   // Checking the number of rows of the left-hand side operand
   if( lhs_.rows() != reflhs_.rows() ) {
      std::ostringstream oss;
      oss << " Test: Initial size comparison of left-hand side row-major dense operand\n"
          << " Error: Invalid number of rows\n"
          << " Details:\n"
          << "   Row-major dense matrix type:\n"
          << "     " << typeid( MT1 ).name() << "\n"
          << "   Detected number of rows = " << lhs_.rows() << "\n"
          << "   Expected number of rows = " << reflhs_.rows() << "\n";
      throw std::runtime_error( oss.str() );
   }

   // Checking the number of columns of the left-hand side operand
   if( lhs_.columns() != reflhs_.columns() ) {
      std::ostringstream oss;
      oss << " Test: Initial size comparison of left-hand side row-major dense operand\n"
          << " Error: Invalid number of columns\n"
          << " Details:\n"
          << "   Row-major dense matrix type:\n"
          << "     " << typeid( MT1 ).name() << "\n"
          << "   Detected number of columns = " << lhs_.columns() << "\n"
          << "   Expected number of columns = " << reflhs_.columns() << "\n";
      throw std::runtime_error( oss.str() );
   }

   // Checking the number of rows of the right-hand side operand
   if( rhs_.rows() != refrhs_.rows() ) {
      std::ostringstream oss;
      oss << " Test: Initial size comparison of right-hand side row-major dense operand\n"
          << " Error: Invalid number of rows\n"
          << " Details:\n"
          << "   Row-major dense matrix type:\n"
          << "     " << typeid( MT2 ).name() << "\n"
          << "   Detected number of rows = " << rhs_.rows() << "\n"
          << "   Expected number of rows = " << refrhs_.rows() << "\n";
      throw std::runtime_error( oss.str() );
   }

   // Checking the number of columns of the right-hand side operand
   if( rhs_.columns() != refrhs_.columns() ) {
      std::ostringstream oss;
      oss << " Test: Initial size comparison of right-hand side row-major dense operand\n"
          << " Error: Invalid number of columns\n"
          << " Details:\n"
          << "   Row-major dense matrix type:\n"
          << "     " << typeid( MT2 ).name() << "\n"
          << "   Detected number of columns = " << rhs_.columns() << "\n"
          << "   Expected number of columns = " << refrhs_.columns() << "\n";
      throw std::runtime_error( oss.str() );
   }

   // Checking the initialization of the left-hand side operand
   if( !isEqual( lhs_, reflhs_ ) ) {
      std::ostringstream oss;
      oss << " Test: Initial test of initialization of left-hand side row-major dense operand\n"
          << " Error: Invalid matrix initialization\n"
          << " Details:\n"
          << "   Row-major dense matrix type:\n"
          << "     " << typeid( MT1 ).name() << "\n"
          << "   Current initialization:\n" << lhs_ << "\n"
          << "   Expected initialization:\n" << reflhs_ << "\n";
      throw std::runtime_error( oss.str() );
   }

   // Checking the initialization of the right-hand side operand
   if( !isEqual( rhs_, refrhs_ ) ) {
      std::ostringstream oss;
      oss << " Test: Initial test of initialization of right-hand side row-major dense operand\n"
          << " Error: Invalid matrix initialization\n"
          << " Details:\n"
          << "   Row-major dense matrix type:\n"
          << "     " << typeid( MT2 ).name() << "\n"
          << "   Current initialization:\n" << rhs_ << "\n"
          << "   Expected initialization:\n" << refrhs_ << "\n";
      throw std::runtime_error( oss.str() );
   }


   //=====================================================================================
   // Performing initial tests with the column-major types
   //=====================================================================================

   // Checking the number of rows of the left-hand side operand
   if( olhs_.rows() != reflhs_.rows() ) {
      std::ostringstream oss;
      oss << " Test: Initial size comparison of left-hand side column-major dense operand\n"
          << " Error: Invalid number of rows\n"
          << " Details:\n"
          << "   Column-major dense matrix type:\n"
          << "     " << typeid( OMT1 ).name() << "\n"
          << "   Detected number of rows = " << olhs_.rows() << "\n"
          << "   Expected number of rows = " << reflhs_.rows() << "\n";
      throw std::runtime_error( oss.str() );
   }

   // Checking the number of columns of the left-hand side operand
   if( olhs_.columns() != reflhs_.columns() ) {
      std::ostringstream oss;
      oss << " Test: Initial size comparison of left-hand side column-major dense operand\n"
          << " Error: Invalid number of columns\n"
          << " Details:\n"
          << "   Column-major dense matrix type:\n"
          << "     " << typeid( OMT1 ).name() << "\n"
          << "   Detected number of columns = " << olhs_.columns() << "\n"
          << "   Expected number of columns = " << reflhs_.columns() << "\n";
      throw std::runtime_error( oss.str() );
   }

   // Checking the number of rows of the right-hand side operand
   if( orhs_.rows() != refrhs_.rows() ) {
      std::ostringstream oss;
      oss << " Test: Initial size comparison of right-hand side column-major dense operand\n"
          << " Error: Invalid number of rows\n"
          << " Details:\n"
          << "   Column-major dense matrix type:\n"
          << "     " << typeid( OMT2 ).name() << "\n"
          << "   Detected number of rows = " << orhs_.rows() << "\n"
          << "   Expected number of rows = " << refrhs_.rows() << "\n";
      throw std::runtime_error( oss.str() );
   }

   // Checking the number of columns of the right-hand side operand
   if( orhs_.columns() != refrhs_.columns() ) {
      std::ostringstream oss;
      oss << " Test: Initial size comparison of right-hand side column-major dense operand\n"
          << " Error: Invalid number of columns\n"
          << " Details:\n"
          << "   Column-major dense matrix type:\n"
          << "     " << typeid( OMT2 ).name() << "\n"
          << "   Detected number of columns = " << orhs_.columns() << "\n"
          << "   Expected number of columns = " << refrhs_.columns() << "\n";
      throw std::runtime_error( oss.str() );
   }

   // Checking the initialization of the left-hand side operand
   if( !isEqual( olhs_, reflhs_ ) ) {
      std::ostringstream oss;
      oss << " Test: Initial test of initialization of left-hand side column-major dense operand\n"
          << " Error: Invalid matrix initialization\n"
          << " Details:\n"
          << "   Column-major dense matrix type:\n"
          << "     " << typeid( OMT1 ).name() << "\n"
          << "   Current initialization:\n" << olhs_ << "\n"
          << "   Expected initialization:\n" << reflhs_ << "\n";
      throw std::runtime_error( oss.str() );
   }

   // Checking the initialization of the right-hand side operand
   if( !isEqual( orhs_, refrhs_ ) ) {
      std::ostringstream oss;
      oss << " Test: Initial test of initialization of right-hand side column-major dense operand\n"
          << " Error: Invalid matrix initialization\n"
          << " Details:\n"
          << "   Column-major dense matrix type:\n"
          << "     " << typeid( OMT2 ).name() << "\n"
          << "   Current initialization:\n" << orhs_ << "\n"
          << "   Expected initialization:\n" << refrhs_ << "\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the matrix assignment.
//
// \return void
// \exception std::runtime_error Assignment error detected.
//
// This function tests the matrix assignment. In case any error is detected, a
// \a std::runtime_error exception is thrown.
*/
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
void OperationTest<MT1,MT2>::testAssignment()
{
   //=====================================================================================
   // // Performing an assignment with the row-major types
   //=====================================================================================

   try {
      lhs_ = reflhs_;
      rhs_ = refrhs_;
   }
   catch( std::exception& ex ) {
      std::ostringstream oss;
      oss << " Test: Assignment with the row-major types\n"
          << " Error: Failed assignment\n"
          << " Details:\n"
          << "   Left-hand side row-major dense matrix type:\n"
          << "     " << typeid( MT1 ).name() << "\n"
          << "   Right-hand side row-major dense matrix type:\n"
          << "     " << typeid( MT2 ).name() << "\n"
          << "   Error message: " << ex.what() << "\n";
      throw std::runtime_error( oss.str() );
   }

   if( !isEqual( lhs_, reflhs_ ) ) {
      std::ostringstream oss;
      oss << " Test: Checking the assignment result of left-hand side row-major dense operand\n"
          << " Error: Invalid matrix initialization\n"
          << " Details:\n"
          << "   Row-major dense matrix type:\n"
          << "     " << typeid( MT1 ).name() << "\n"
          << "   Current initialization:\n" << lhs_ << "\n"
          << "   Expected initialization:\n" << reflhs_ << "\n";
      throw std::runtime_error( oss.str() );
   }

   if( !isEqual( rhs_, refrhs_ ) ) {
      std::ostringstream oss;
      oss << " Test: Checking the assignment result of right-hand side row-major dense operand\n"
          << " Error: Invalid matrix initialization\n"
          << " Details:\n"
          << "   Row-major dense matrix type:\n"
          << "     " << typeid( MT2 ).name() << "\n"
          << "   Current initialization:\n" << rhs_ << "\n"
          << "   Expected initialization:\n" << refrhs_ << "\n";
      throw std::runtime_error( oss.str() );
   }


   //=====================================================================================
   // Performing an assignment with the column-major types
   //=====================================================================================

   try {
      olhs_ = reflhs_;
      orhs_ = refrhs_;
   }
   catch( std::exception& ex ) {
      std::ostringstream oss;
      oss << " Test: Assignment with the column-major types\n"
          << " Error: Failed assignment\n"
          << " Details:\n"
          << "   Left-hand side column-major dense matrix type:\n"
          << "     " << typeid( OMT1 ).name() << "\n"
          << "   Right-hand side column-major dense matrix type:\n"
          << "     "  << typeid( OMT2 ).name() << "\n"
          << "   Error message: " << ex.what() << "\n";
      throw std::runtime_error( oss.str() );
   }

   if( !isEqual( olhs_, reflhs_ ) ) {
      std::ostringstream oss;
      oss << " Test: Checking the assignment result of left-hand side column-major dense operand\n"
          << " Error: Invalid matrix initialization\n"
          << " Details:\n"
          << "   Column-major dense matrix type:\n"
          << "     " << typeid( OMT1 ).name() << "\n"
          << "   Current initialization:\n" << olhs_ << "\n"
          << "   Expected initialization:\n" << reflhs_ << "\n";
      throw std::runtime_error( oss.str() );
   }

   if( !isEqual( orhs_, refrhs_ ) ) {
      std::ostringstream oss;
      oss << " Test: Checking the assignment result of right-hand side column-major dense operand\n"
          << " Error: Invalid matrix initialization\n"
          << " Details:\n"
          << "   Column-major dense matrix type:\n"
          << "     " << typeid( OMT2 ).name() << "\n"
          << "   Current initialization:\n" << orhs_ << "\n"
          << "   Expected initialization:\n" << refrhs_ << "\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the matrix element access.
//
// \return void
// \exception std::runtime_error Element access error detected.
//
// This function tests the element access via the subscript operator. In case any
// error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
void OperationTest<MT1,MT2>::testElementAccess()
{
   using blaze::equal;


   //=====================================================================================
   // Testing the element access with two row-major matrices
   //=====================================================================================

   if( lhs_.rows() > 0UL && lhs_.columns() > 0UL )
   {
      if( !equal( ( lhs_ + rhs_ )(0UL,0UL), ( reflhs_ + refrhs_ )(0UL,0UL) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of addition expression\n"
             << " Error: Unequal resulting elements at element (0,0) detected\n"
             << " Details:\n"
             << "   Left-hand side row-major dense matrix type:\n"
             << "     " << typeid( MT1 ).name() << "\n"
             << "   Right-hand side row-major dense matrix type:\n"
             << "     " << typeid( MT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( ( lhs_ + eval( rhs_ ) )(0UL,0UL), ( reflhs_ + eval( refrhs_ ) )(0UL,0UL) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of right evaluated addition expression\n"
             << " Error: Unequal resulting elements at element (0,0) detected\n"
             << " Details:\n"
             << "   Left-hand side row-major dense matrix type:\n"
             << "     " << typeid( MT1 ).name() << "\n"
             << "   Right-hand side row-major dense matrix type:\n"
             << "     " << typeid( MT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( ( eval( lhs_ ) + rhs_ )(0UL,0UL), ( eval( reflhs_ ) + refrhs_ )(0UL,0UL) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of left evaluated addition expression\n"
             << " Error: Unequal resulting elements at element (0,0) detected\n"
             << " Details:\n"
             << "   Left-hand side row-major dense matrix type:\n"
             << "     " << typeid( MT1 ).name() << "\n"
             << "   Right-hand side row-major dense matrix type:\n"
             << "     " << typeid( MT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( ( eval( lhs_ ) + eval( rhs_ ) )(0UL,0UL), ( eval( reflhs_ ) + eval( refrhs_ ) )(0UL,0UL) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of fully evaluated addition expression\n"
             << " Error: Unequal resulting elements at element (0,0) detected\n"
             << " Details:\n"
             << "   Left-hand side row-major dense matrix type:\n"
             << "     " << typeid( MT1 ).name() << "\n"
             << "   Right-hand side row-major dense matrix type:\n"
             << "     " << typeid( MT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Testing the element access with a row-major matrix and a column-major matrix
   //=====================================================================================

   if( lhs_.rows() > 0UL && lhs_.columns() > 0UL )
   {
      if( !equal( ( lhs_ + orhs_ )(0UL,0UL), ( reflhs_ + refrhs_ )(0UL,0UL) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of addition expression\n"
             << " Error: Unequal resulting elements at element (0,0) detected\n"
             << " Details:\n"
             << "   Left-hand side row-major dense matrix type:\n"
             << "     " << typeid( MT1 ).name() << "\n"
             << "   Right-hand side column-major dense matrix type:\n"
             << "     " << typeid( OMT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( ( lhs_ + eval( orhs_ ) )(0UL,0UL), ( reflhs_ + eval( refrhs_ ) )(0UL,0UL) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of right evaluated addition expression\n"
             << " Error: Unequal resulting elements at element (0,0) detected\n"
             << " Details:\n"
             << "   Left-hand side row-major dense matrix type:\n"
             << "     " << typeid( MT1 ).name() << "\n"
             << "   Right-hand side column-major dense matrix type:\n"
             << "     " << typeid( OMT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( ( eval( lhs_ ) + orhs_ )(0UL,0UL), ( eval( reflhs_ ) + refrhs_ )(0UL,0UL) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of left evaluated addition expression\n"
             << " Error: Unequal resulting elements at element (0,0) detected\n"
             << " Details:\n"
             << "   Left-hand side row-major dense matrix type:\n"
             << "     " << typeid( MT1 ).name() << "\n"
             << "   Right-hand side column-major dense matrix type:\n"
             << "     " << typeid( OMT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( ( eval( lhs_ ) + eval( orhs_ ) )(0UL,0UL), ( eval( reflhs_ ) + eval( refrhs_ ) )(0UL,0UL) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of fully evaluated addition expression\n"
             << " Error: Unequal resulting elements at element (0,0) detected\n"
             << " Details:\n"
             << "   Left-hand side row-major dense matrix type:\n"
             << "     " << typeid( MT1 ).name() << "\n"
             << "   Right-hand side column-major dense matrix type:\n"
             << "     " << typeid( OMT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Testing the element access with a column-major matrix and a row-major matrix
   //=====================================================================================

   if( lhs_.rows() > 0UL && lhs_.columns() > 0UL )
   {
      if( !equal( ( olhs_ + rhs_ )(0UL,0UL), ( reflhs_ + refrhs_ )(0UL,0UL) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of addition expression\n"
             << " Error: Unequal resulting elements at element (0,0) detected\n"
             << " Details:\n"
             << "   Left-hand side column-major dense matrix type:\n"
             << "     " << typeid( OMT1 ).name() << "\n"
             << "   Right-hand side row-major dense matrix type:\n"
             << "     " << typeid( MT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( ( olhs_ + eval( rhs_ ) )(0UL,0UL), ( reflhs_ + eval( refrhs_ ) )(0UL,0UL) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of right evaluated addition expression\n"
             << " Error: Unequal resulting elements at element (0,0) detected\n"
             << " Details:\n"
             << "   Left-hand side column-major dense matrix type:\n"
             << "     " << typeid( OMT1 ).name() << "\n"
             << "   Right-hand side row-major dense matrix type:\n"
             << "     " << typeid( MT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( ( eval( olhs_ ) + rhs_ )(0UL,0UL), ( eval( reflhs_ ) + refrhs_ )(0UL,0UL) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of left evaluated addition expression\n"
             << " Error: Unequal resulting elements at element (0,0) detected\n"
             << " Details:\n"
             << "   Left-hand side column-major dense matrix type:\n"
             << "     " << typeid( OMT1 ).name() << "\n"
             << "   Right-hand side row-major dense matrix type:\n"
             << "     " << typeid( MT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( ( eval( olhs_ ) + eval( rhs_ ) )(0UL,0UL), ( eval( reflhs_ ) + eval( refrhs_ ) )(0UL,0UL) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of fully evaluated addition expression\n"
             << " Error: Unequal resulting elements at element (0,0) detected\n"
             << " Details:\n"
             << "   Left-hand side column-major dense matrix type:\n"
             << "     " << typeid( OMT1 ).name() << "\n"
             << "   Right-hand side row-major dense matrix type:\n"
             << "     " << typeid( MT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Testing the element access with two column-major matrices
   //=====================================================================================

   if( olhs_.rows() > 0UL && olhs_.columns() > 0UL )
   {
      if( !equal( ( olhs_ + orhs_ )(0UL,0UL), ( reflhs_ + refrhs_ )(0UL,0UL) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of transpose addition expression\n"
             << " Error: Unequal resulting elements at element (0,0) detected\n"
             << " Details:\n"
             << "   Transpose left-hand side dense matrix type:\n"
             << "     " << typeid( OMT1 ).name() << "\n"
             << "   Transpose right-hand side dense matrix type:\n"
             << "     " << typeid( OMT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( ( olhs_ + eval( orhs_ ) )(0UL,0UL), ( reflhs_ + eval( refrhs_ ) )(0UL,0UL) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of right evaluated transpose addition expression\n"
             << " Error: Unequal resulting elements at element (0,0) detected\n"
             << " Details:\n"
             << "   Transpose left-hand side dense matrix type:\n"
             << "     " << typeid( OMT1 ).name() << "\n"
             << "   Transpose right-hand side dense matrix type:\n"
             << "     " << typeid( OMT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( ( eval( olhs_ ) + orhs_ )(0UL,0UL), ( eval( reflhs_ ) + refrhs_ )(0UL,0UL) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of left evaluated transpose addition expression\n"
             << " Error: Unequal resulting elements at element (0,0) detected\n"
             << " Details:\n"
             << "   Transpose left-hand side dense matrix type:\n"
             << "     " << typeid( OMT1 ).name() << "\n"
             << "   Transpose right-hand side dense matrix type:\n"
             << "     " << typeid( OMT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( ( eval( olhs_ ) + eval( orhs_ ) )(0UL,0UL), ( eval( reflhs_ ) + eval( refrhs_ ) )(0UL,0UL) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of fully evaluated transpose addition expression\n"
             << " Error: Unequal resulting elements at element (0,0) detected\n"
             << " Details:\n"
             << "   Transpose left-hand side dense matrix type:\n"
             << "     " << typeid( OMT1 ).name() << "\n"
             << "   Transpose right-hand side dense matrix type:\n"
             << "     " << typeid( OMT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the plain dense matrix/dense matrix addition.
//
// \return void
// \exception std::runtime_error Addition error detected.
//
// This function tests the plain matrix addition with plain assignment, addition assignment,
// and subtraction assignment. In case any error resulting from the addition or the subsequent
// assignment is detected, a \a std::runtime_error exception is thrown.
*/
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
void OperationTest<MT1,MT2>::testBasicOperation()
{
#if BLAZETEST_MATHTEST_TEST_BASIC_OPERATION
   if( BLAZETEST_MATHTEST_TEST_BASIC_OPERATION > 1 )
   {
      //=====================================================================================
      // Addition
      //=====================================================================================

      // Addition with the given matrices
      {
         test_  = "Addition with the given matrices";
         error_ = "Failed addition operation";

         try {
            initResults();
            dres_   = lhs_ + rhs_;
            odres_  = lhs_ + rhs_;
            sres_   = lhs_ + rhs_;
            osres_  = lhs_ + rhs_;
            refres_ = reflhs_ + refrhs_;
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   = lhs_ + orhs_;
            odres_  = lhs_ + orhs_;
            sres_   = lhs_ + orhs_;
            osres_  = lhs_ + orhs_;
            refres_ = reflhs_ + refrhs_;
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   = olhs_ + rhs_;
            odres_  = olhs_ + rhs_;
            sres_   = olhs_ + rhs_;
            osres_  = olhs_ + rhs_;
            refres_ = reflhs_ + refrhs_;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   = olhs_ + orhs_;
            odres_  = olhs_ + orhs_;
            sres_   = olhs_ + orhs_;
            osres_  = olhs_ + orhs_;
            refres_ = reflhs_ + refrhs_;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Addition with evaluated matrices
      {
         test_  = "Addition with evaluated matrices";
         error_ = "Failed addition operation";

         try {
            initResults();
            dres_   = eval( lhs_ ) + eval( rhs_ );
            odres_  = eval( lhs_ ) + eval( rhs_ );
            sres_   = eval( lhs_ ) + eval( rhs_ );
            osres_  = eval( lhs_ ) + eval( rhs_ );
            refres_ = eval( reflhs_ ) + eval( refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   = eval( lhs_ ) + eval( orhs_ );
            odres_  = eval( lhs_ ) + eval( orhs_ );
            sres_   = eval( lhs_ ) + eval( orhs_ );
            osres_  = eval( lhs_ ) + eval( orhs_ );
            refres_ = eval( reflhs_ ) + eval( refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   = eval( olhs_ ) + eval( rhs_ );
            odres_  = eval( olhs_ ) + eval( rhs_ );
            sres_   = eval( olhs_ ) + eval( rhs_ );
            osres_  = eval( olhs_ ) + eval( rhs_ );
            refres_ = eval( reflhs_ ) + eval( refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   = eval( olhs_ ) + eval( orhs_ );
            odres_  = eval( olhs_ ) + eval( orhs_ );
            sres_   = eval( olhs_ ) + eval( orhs_ );
            osres_  = eval( olhs_ ) + eval( orhs_ );
            refres_ = eval( reflhs_ ) + eval( refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Addition with addition assignment
      //=====================================================================================

      // Addition with addition assignment with the given matrices
      {
         test_  = "Addition with addition assignment with the given matrices";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += lhs_ + rhs_;
            odres_  += lhs_ + rhs_;
            sres_   += lhs_ + rhs_;
            osres_  += lhs_ + rhs_;
            refres_ += reflhs_ + refrhs_;
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   += lhs_ + orhs_;
            odres_  += lhs_ + orhs_;
            sres_   += lhs_ + orhs_;
            osres_  += lhs_ + orhs_;
            refres_ += reflhs_ + refrhs_;
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   += olhs_ + rhs_;
            odres_  += olhs_ + rhs_;
            sres_   += olhs_ + rhs_;
            osres_  += olhs_ + rhs_;
            refres_ += reflhs_ + refrhs_;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   += olhs_ + orhs_;
            odres_  += olhs_ + orhs_;
            sres_   += olhs_ + orhs_;
            osres_  += olhs_ + orhs_;
            refres_ += reflhs_ + refrhs_;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Addition with addition assignment with evaluated matrices
      {
         test_  = "Addition with addition assignment with evaluated matrices";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += eval( lhs_ ) + eval( rhs_ );
            odres_  += eval( lhs_ ) + eval( rhs_ );
            sres_   += eval( lhs_ ) + eval( rhs_ );
            osres_  += eval( lhs_ ) + eval( rhs_ );
            refres_ += eval( reflhs_ ) + eval( refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   += eval( lhs_ ) + eval( orhs_ );
            odres_  += eval( lhs_ ) + eval( orhs_ );
            sres_   += eval( lhs_ ) + eval( orhs_ );
            osres_  += eval( lhs_ ) + eval( orhs_ );
            refres_ += eval( reflhs_ ) + eval( refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   += eval( olhs_ ) + eval( rhs_ );
            odres_  += eval( olhs_ ) + eval( rhs_ );
            sres_   += eval( olhs_ ) + eval( rhs_ );
            osres_  += eval( olhs_ ) + eval( rhs_ );
            refres_ += eval( reflhs_ ) + eval( refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   += eval( olhs_ ) + eval( orhs_ );
            odres_  += eval( olhs_ ) + eval( orhs_ );
            sres_   += eval( olhs_ ) + eval( orhs_ );
            osres_  += eval( olhs_ ) + eval( orhs_ );
            refres_ += eval( reflhs_ ) + eval( refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Addition with subtraction assignment with the given matrices
      //=====================================================================================

      // Addition with subtraction assignment with the given matrices
      {
         test_  = "Addition with subtraction assignment with the given matrices";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= lhs_ + rhs_;
            odres_  -= lhs_ + rhs_;
            sres_   -= lhs_ + rhs_;
            osres_  -= lhs_ + rhs_;
            refres_ -= reflhs_ + refrhs_;
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   -= lhs_ + orhs_;
            odres_  -= lhs_ + orhs_;
            sres_   -= lhs_ + orhs_;
            osres_  -= lhs_ + orhs_;
            refres_ -= reflhs_ + refrhs_;
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   -= olhs_ + rhs_;
            odres_  -= olhs_ + rhs_;
            sres_   -= olhs_ + rhs_;
            osres_  -= olhs_ + rhs_;
            refres_ -= reflhs_ + refrhs_;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   -= olhs_ + orhs_;
            odres_  -= olhs_ + orhs_;
            sres_   -= olhs_ + orhs_;
            osres_  -= olhs_ + orhs_;
            refres_ -= reflhs_ + refrhs_;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Addition with subtraction assignment with evaluated matrices
      {
         test_  = "Addition with subtraction assignment with evaluated matrices";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= eval( lhs_ ) + eval( rhs_ );
            odres_  -= eval( lhs_ ) + eval( rhs_ );
            sres_   -= eval( lhs_ ) + eval( rhs_ );
            osres_  -= eval( lhs_ ) + eval( rhs_ );
            refres_ -= eval( reflhs_ ) + eval( refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   -= eval( lhs_ ) + eval( orhs_ );
            odres_  -= eval( lhs_ ) + eval( orhs_ );
            sres_   -= eval( lhs_ ) + eval( orhs_ );
            osres_  -= eval( lhs_ ) + eval( orhs_ );
            refres_ -= eval( reflhs_ ) + eval( refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   -= eval( olhs_ ) + eval( rhs_ );
            odres_  -= eval( olhs_ ) + eval( rhs_ );
            sres_   -= eval( olhs_ ) + eval( rhs_ );
            osres_  -= eval( olhs_ ) + eval( rhs_ );
            refres_ -= eval( reflhs_ ) + eval( refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   -= eval( olhs_ ) + eval( orhs_ );
            odres_  -= eval( olhs_ ) + eval( orhs_ );
            sres_   -= eval( olhs_ ) + eval( orhs_ );
            osres_  -= eval( olhs_ ) + eval( orhs_ );
            refres_ -= eval( reflhs_ ) + eval( refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the negated dense matrix/dense matrix addition.
//
// \return void
// \exception std::runtime_error Addition error detected.
//
// This function tests the negated matrix addition with plain assignment, addition assignment,
// and subtraction assignment. In case any error resulting from the addition or the subsequent
// assignment is detected, a \a std::runtime_error exception is thrown.
*/
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
void OperationTest<MT1,MT2>::testNegatedOperation()
{
#if BLAZETEST_MATHTEST_TEST_NEGATED_OPERATION
   if( BLAZETEST_MATHTEST_TEST_NEGATED_OPERATION > 1 )
   {
      //=====================================================================================
      // Negated addition
      //=====================================================================================

      // Negated addition with the given matrices
      {
         test_  = "Negated addition with the given matrices";
         error_ = "Failed addition operation";

         try {
            initResults();
            dres_   = -( lhs_ + rhs_ );
            odres_  = -( lhs_ + rhs_ );
            sres_   = -( lhs_ + rhs_ );
            osres_  = -( lhs_ + rhs_ );
            refres_ = -( reflhs_ + refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   = -( lhs_ + orhs_ );
            odres_  = -( lhs_ + orhs_ );
            sres_   = -( lhs_ + orhs_ );
            osres_  = -( lhs_ + orhs_ );
            refres_ = -( reflhs_ + refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   = -( olhs_ + rhs_ );
            odres_  = -( olhs_ + rhs_ );
            sres_   = -( olhs_ + rhs_ );
            osres_  = -( olhs_ + rhs_ );
            refres_ = -( reflhs_ + refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   = -( olhs_ + orhs_ );
            odres_  = -( olhs_ + orhs_ );
            sres_   = -( olhs_ + orhs_ );
            osres_  = -( olhs_ + orhs_ );
            refres_ = -( reflhs_ + refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Negated addition with evaluated matrices
      {
         test_  = "Negated addition with evaluated matrices";
         error_ = "Failed addition operation";

         try {
            initResults();
            dres_   = -( eval( lhs_ ) + eval( rhs_ ) );
            odres_  = -( eval( lhs_ ) + eval( rhs_ ) );
            sres_   = -( eval( lhs_ ) + eval( rhs_ ) );
            osres_  = -( eval( lhs_ ) + eval( rhs_ ) );
            refres_ = -( eval( reflhs_ ) + eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   = -( eval( lhs_ ) + eval( orhs_ ) );
            odres_  = -( eval( lhs_ ) + eval( orhs_ ) );
            sres_   = -( eval( lhs_ ) + eval( orhs_ ) );
            osres_  = -( eval( lhs_ ) + eval( orhs_ ) );
            refres_ = -( eval( reflhs_ ) + eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   = -( eval( olhs_ ) + eval( rhs_ ) );
            odres_  = -( eval( olhs_ ) + eval( rhs_ ) );
            sres_   = -( eval( olhs_ ) + eval( rhs_ ) );
            osres_  = -( eval( olhs_ ) + eval( rhs_ ) );
            refres_ = -( eval( reflhs_ ) + eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   = -( eval( olhs_ ) + eval( orhs_ ) );
            odres_  = -( eval( olhs_ ) + eval( orhs_ ) );
            sres_   = -( eval( olhs_ ) + eval( orhs_ ) );
            osres_  = -( eval( olhs_ ) + eval( orhs_ ) );
            refres_ = -( eval( reflhs_ ) + eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Negated addition with addition assignment
      //=====================================================================================

      // Negated addition with addition assignment with the given matrices
      {
         test_  = "Negated addition with addition assignment with the given matrices";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += -( lhs_ + rhs_ );
            odres_  += -( lhs_ + rhs_ );
            sres_   += -( lhs_ + rhs_ );
            osres_  += -( lhs_ + rhs_ );
            refres_ += -( reflhs_ + refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   += -( lhs_ + orhs_ );
            odres_  += -( lhs_ + orhs_ );
            sres_   += -( lhs_ + orhs_ );
            osres_  += -( lhs_ + orhs_ );
            refres_ += -( reflhs_ + refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   += -( olhs_ + rhs_ );
            odres_  += -( olhs_ + rhs_ );
            sres_   += -( olhs_ + rhs_ );
            osres_  += -( olhs_ + rhs_ );
            refres_ += -( reflhs_ + refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   += -( olhs_ + orhs_ );
            odres_  += -( olhs_ + orhs_ );
            sres_   += -( olhs_ + orhs_ );
            osres_  += -( olhs_ + orhs_ );
            refres_ += -( reflhs_ + refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Negated addition with addition assignment with the given matrices
      {
         test_  = "Negated addition with addition assignment with evaluated matrices";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += -( eval( lhs_ ) + eval( rhs_ ) );
            odres_  += -( eval( lhs_ ) + eval( rhs_ ) );
            sres_   += -( eval( lhs_ ) + eval( rhs_ ) );
            osres_  += -( eval( lhs_ ) + eval( rhs_ ) );
            refres_ += -( eval( reflhs_ ) + eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   += -( eval( lhs_ ) + eval( orhs_ ) );
            odres_  += -( eval( lhs_ ) + eval( orhs_ ) );
            sres_   += -( eval( lhs_ ) + eval( orhs_ ) );
            osres_  += -( eval( lhs_ ) + eval( orhs_ ) );
            refres_ += -( eval( reflhs_ ) + eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   += -( eval( olhs_ ) + eval( rhs_ ) );
            odres_  += -( eval( olhs_ ) + eval( rhs_ ) );
            sres_   += -( eval( olhs_ ) + eval( rhs_ ) );
            osres_  += -( eval( olhs_ ) + eval( rhs_ ) );
            refres_ += -( eval( reflhs_ ) + eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   += -( eval( olhs_ ) + eval( orhs_ ) );
            odres_  += -( eval( olhs_ ) + eval( orhs_ ) );
            sres_   += -( eval( olhs_ ) + eval( orhs_ ) );
            osres_  += -( eval( olhs_ ) + eval( orhs_ ) );
            refres_ += -( eval( reflhs_ ) + eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Negated addition with subtraction assignment
      //=====================================================================================

      // Negated addition with subtraction assignment with the given matrices
      {
         test_  = "Negated addition with subtraction assignment with the given matrices";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= -( lhs_ + rhs_ );
            odres_  -= -( lhs_ + rhs_ );
            sres_   -= -( lhs_ + rhs_ );
            osres_  -= -( lhs_ + rhs_ );
            refres_ -= -( reflhs_ + refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   -= -( lhs_ + orhs_ );
            odres_  -= -( lhs_ + orhs_ );
            sres_   -= -( lhs_ + orhs_ );
            osres_  -= -( lhs_ + orhs_ );
            refres_ -= -( reflhs_ + refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   -= -( olhs_ + rhs_ );
            odres_  -= -( olhs_ + rhs_ );
            sres_   -= -( olhs_ + rhs_ );
            osres_  -= -( olhs_ + rhs_ );
            refres_ -= -( reflhs_ + refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   -= -( olhs_ + orhs_ );
            odres_  -= -( olhs_ + orhs_ );
            sres_   -= -( olhs_ + orhs_ );
            osres_  -= -( olhs_ + orhs_ );
            refres_ -= -( reflhs_ + refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Negated addition with subtraction assignment with evaluated matrices
      {
         test_  = "Negated addition with subtraction assignment with evaluated matrices";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= -( eval( lhs_ ) + eval( rhs_ ) );
            odres_  -= -( eval( lhs_ ) + eval( rhs_ ) );
            sres_   -= -( eval( lhs_ ) + eval( rhs_ ) );
            osres_  -= -( eval( lhs_ ) + eval( rhs_ ) );
            refres_ -= -( eval( reflhs_ ) + eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   -= -( eval( lhs_ ) + eval( orhs_ ) );
            odres_  -= -( eval( lhs_ ) + eval( orhs_ ) );
            sres_   -= -( eval( lhs_ ) + eval( orhs_ ) );
            osres_  -= -( eval( lhs_ ) + eval( orhs_ ) );
            refres_ -= -( eval( reflhs_ ) + eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   -= -( eval( olhs_ ) + eval( rhs_ ) );
            odres_  -= -( eval( olhs_ ) + eval( rhs_ ) );
            sres_   -= -( eval( olhs_ ) + eval( rhs_ ) );
            osres_  -= -( eval( olhs_ ) + eval( rhs_ ) );
            refres_ -= -( eval( reflhs_ ) + eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   -= -( eval( olhs_ ) + eval( orhs_ ) );
            odres_  -= -( eval( olhs_ ) + eval( orhs_ ) );
            sres_   -= -( eval( olhs_ ) + eval( orhs_ ) );
            osres_  -= -( eval( olhs_ ) + eval( orhs_ ) );
            refres_ -= -( eval( reflhs_ ) + eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the scaled dense matrix/dense matrix addition.
//
// \param scalar The scalar value.
// \return void
// \exception std::runtime_error Addition error detected.
//
// This function tests the scaled matrix addition with plain assignment, addition assignment,
// and subtraction assignment. In case any error resulting from the addition or the subsequent
// assignment is detected, a \a std::runtime_error exception is thrown.
*/
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
template< typename T >    // Type of the scalar
void OperationTest<MT1,MT2>::testScaledOperation( T scalar )
{
   BLAZE_CONSTRAINT_MUST_BE_NUMERIC_TYPE( T );

   if( scalar == T(0) )
      throw std::invalid_argument( "Invalid scalar parameter" );


#if BLAZETEST_MATHTEST_TEST_SCALED_OPERATION
   if( BLAZETEST_MATHTEST_TEST_SCALED_OPERATION > 1 )
   {
      //=====================================================================================
      // Self-scaling (M*=s)
      //=====================================================================================

      // Self-scaling (M*=s)
      {
         test_ = "Self-scaling (M*=s)";

         try {
            dres_   = lhs_ + rhs_;
            odres_  = dres_;
            sres_   = dres_;
            osres_  = dres_;
            refres_ = dres_;

            dres_   *= scalar;
            odres_  *= scalar;
            sres_   *= scalar;
            osres_  *= scalar;
            refres_ *= scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed self-scaling operation\n"
                << " Details:\n"
                << "   Scalar = " << scalar << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<MT1,MT2>();
      }


      //=====================================================================================
      // Self-scaling (M=M*s)
      //=====================================================================================

      // Self-scaling (M=M*s)
      {
         test_ = "Self-scaling (M=M*s)";

         try {
            dres_   = lhs_ + rhs_;
            odres_  = dres_;
            sres_   = dres_;
            osres_  = dres_;
            refres_ = dres_;

            dres_   = dres_   * scalar;
            odres_  = odres_  * scalar;
            sres_   = sres_   * scalar;
            osres_  = osres_  * scalar;
            refres_ = refres_ * scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed self-scaling operation\n"
                << " Details:\n"
                << "   Scalar = " << scalar << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<MT1,MT2>();
      }


      //=====================================================================================
      // Self-scaling (M=s*M)
      //=====================================================================================

      // Self-scaling (M=s*M)
      {
         test_ = "Self-scaling (M=s*M)";

         try {
            dres_   = lhs_ + rhs_;
            odres_  = dres_;
            sres_   = dres_;
            osres_  = dres_;
            refres_ = dres_;

            dres_   = scalar * dres_;
            odres_  = scalar * odres_;
            sres_   = scalar * sres_;
            osres_  = scalar * osres_;
            refres_ = scalar * refres_;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed self-scaling operation\n"
                << " Details:\n"
                << "   Scalar = " << scalar << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<MT1,MT2>();
      }


      //=====================================================================================
      // Self-scaling (M/=s)
      //=====================================================================================

      // Self-scaling (M/=s)
      {
         test_ = "Self-scaling (M/=s)";

         try {
            dres_   = lhs_ + rhs_;
            odres_  = dres_;
            sres_   = dres_;
            osres_  = dres_;
            refres_ = dres_;

            dres_   /= scalar;
            odres_  /= scalar;
            sres_   /= scalar;
            osres_  /= scalar;
            refres_ /= scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed self-scaling operation\n"
                << " Details:\n"
                << "   Scalar = " << scalar << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<MT1,MT2>();
      }


      //=====================================================================================
      // Self-scaling (M=M/s)
      //=====================================================================================

      // Self-scaling (M=M/s)
      {
         test_ = "Self-scaling (M=M/s)";

         try {
            dres_   = lhs_ + rhs_;
            odres_  = dres_;
            sres_   = dres_;
            osres_  = dres_;
            refres_ = dres_;

            dres_   = dres_   / scalar;
            odres_  = odres_  / scalar;
            sres_   = sres_   / scalar;
            osres_  = osres_  / scalar;
            refres_ = refres_ / scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed self-scaling operation\n"
                << " Details:\n"
                << "   Scalar = " << scalar << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<MT1,MT2>();
      }


      //=====================================================================================
      // Scaled addition (s*OP)
      //=====================================================================================

      // Scaled addition with the given matrices
      {
         test_  = "Scaled addition with the given matrices (s*OP)";
         error_ = "Failed addition operation";

         try {
            initResults();
            dres_   = scalar * ( lhs_ + rhs_ );
            odres_  = scalar * ( lhs_ + rhs_ );
            sres_   = scalar * ( lhs_ + rhs_ );
            osres_  = scalar * ( lhs_ + rhs_ );
            refres_ = scalar * ( reflhs_ + refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   = scalar * ( lhs_ + orhs_ );
            odres_  = scalar * ( lhs_ + orhs_ );
            sres_   = scalar * ( lhs_ + orhs_ );
            osres_  = scalar * ( lhs_ + orhs_ );
            refres_ = scalar * ( reflhs_ + refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   = scalar * ( olhs_ + rhs_ );
            odres_  = scalar * ( olhs_ + rhs_ );
            sres_   = scalar * ( olhs_ + rhs_ );
            osres_  = scalar * ( olhs_ + rhs_ );
            refres_ = scalar * ( reflhs_ + refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   = scalar * ( olhs_ + orhs_ );
            odres_  = scalar * ( olhs_ + orhs_ );
            sres_   = scalar * ( olhs_ + orhs_ );
            osres_  = scalar * ( olhs_ + orhs_ );
            refres_ = scalar * ( reflhs_ + refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Scaled addition with evaluated matrices
      {
         test_  = "Scaled addition with evaluated matrices (s*OP)";
         error_ = "Failed addition operation";

         try {
            initResults();
            dres_   = scalar * ( eval( lhs_ ) + eval( rhs_ ) );
            odres_  = scalar * ( eval( lhs_ ) + eval( rhs_ ) );
            sres_   = scalar * ( eval( lhs_ ) + eval( rhs_ ) );
            osres_  = scalar * ( eval( lhs_ ) + eval( rhs_ ) );
            refres_ = scalar * ( eval( reflhs_ ) + eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   = scalar * ( eval( lhs_ ) + eval( orhs_ ) );
            odres_  = scalar * ( eval( lhs_ ) + eval( orhs_ ) );
            sres_   = scalar * ( eval( lhs_ ) + eval( orhs_ ) );
            osres_  = scalar * ( eval( lhs_ ) + eval( orhs_ ) );
            refres_ = scalar * ( eval( reflhs_ ) + eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   = scalar * ( eval( olhs_ ) + eval( rhs_ ) );
            odres_  = scalar * ( eval( olhs_ ) + eval( rhs_ ) );
            sres_   = scalar * ( eval( olhs_ ) + eval( rhs_ ) );
            osres_  = scalar * ( eval( olhs_ ) + eval( rhs_ ) );
            refres_ = scalar * ( eval( reflhs_ ) + eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   = scalar * ( eval( olhs_ ) + eval( orhs_ ) );
            odres_  = scalar * ( eval( olhs_ ) + eval( orhs_ ) );
            sres_   = scalar * ( eval( olhs_ ) + eval( orhs_ ) );
            osres_  = scalar * ( eval( olhs_ ) + eval( orhs_ ) );
            refres_ = scalar * ( eval( reflhs_ ) + eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Scaled addition (OP*s)
      //=====================================================================================

      // Scaled addition with the given matrices
      {
         test_  = "Scaled addition with the given matrices (OP*s)";
         error_ = "Failed addition operation";

         try {
            initResults();
            dres_   = ( lhs_ + rhs_ ) * scalar;
            odres_  = ( lhs_ + rhs_ ) * scalar;
            sres_   = ( lhs_ + rhs_ ) * scalar;
            osres_  = ( lhs_ + rhs_ ) * scalar;
            refres_ = ( reflhs_ + refrhs_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   = ( lhs_ + orhs_ ) * scalar;
            odres_  = ( lhs_ + orhs_ ) * scalar;
            sres_   = ( lhs_ + orhs_ ) * scalar;
            osres_  = ( lhs_ + orhs_ ) * scalar;
            refres_ = ( reflhs_ + refrhs_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   = ( olhs_ + rhs_ ) * scalar;
            odres_  = ( olhs_ + rhs_ ) * scalar;
            sres_   = ( olhs_ + rhs_ ) * scalar;
            osres_  = ( olhs_ + rhs_ ) * scalar;
            refres_ = ( reflhs_ + refrhs_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   = ( olhs_ + orhs_ ) * scalar;
            odres_  = ( olhs_ + orhs_ ) * scalar;
            sres_   = ( olhs_ + orhs_ ) * scalar;
            osres_  = ( olhs_ + orhs_ ) * scalar;
            refres_ = ( reflhs_ + refrhs_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Scaled addition with evaluated matrices
      {
         test_  = "Scaled addition with evaluated matrices (OP*s)";
         error_ = "Failed addition operation";

         try {
            initResults();
            dres_   = ( eval( lhs_ ) + eval( rhs_ ) ) * scalar;
            odres_  = ( eval( lhs_ ) + eval( rhs_ ) ) * scalar;
            sres_   = ( eval( lhs_ ) + eval( rhs_ ) ) * scalar;
            osres_  = ( eval( lhs_ ) + eval( rhs_ ) ) * scalar;
            refres_ = ( eval( reflhs_ ) + eval( refrhs_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   = ( eval( lhs_ ) + eval( orhs_ ) ) * scalar;
            odres_  = ( eval( lhs_ ) + eval( orhs_ ) ) * scalar;
            sres_   = ( eval( lhs_ ) + eval( orhs_ ) ) * scalar;
            osres_  = ( eval( lhs_ ) + eval( orhs_ ) ) * scalar;
            refres_ = ( eval( reflhs_ ) + eval( refrhs_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   = ( eval( olhs_ ) + eval( rhs_ ) ) * scalar;
            odres_  = ( eval( olhs_ ) + eval( rhs_ ) ) * scalar;
            sres_   = ( eval( olhs_ ) + eval( rhs_ ) ) * scalar;
            osres_  = ( eval( olhs_ ) + eval( rhs_ ) ) * scalar;
            refres_ = ( eval( reflhs_ ) + eval( refrhs_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   = ( eval( olhs_ ) + eval( orhs_ ) ) * scalar;
            odres_  = ( eval( olhs_ ) + eval( orhs_ ) ) * scalar;
            sres_   = ( eval( olhs_ ) + eval( orhs_ ) ) * scalar;
            osres_  = ( eval( olhs_ ) + eval( orhs_ ) ) * scalar;
            refres_ = ( eval( reflhs_ ) + eval( refrhs_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Scaled addition (OP/s)
      //=====================================================================================

      // Scaled addition with the given matrices
      {
         test_  = "Scaled addition with the given matrices (OP/s)";
         error_ = "Failed addition operation";

         try {
            initResults();
            dres_   = ( lhs_ + rhs_ ) / scalar;
            odres_  = ( lhs_ + rhs_ ) / scalar;
            sres_   = ( lhs_ + rhs_ ) / scalar;
            osres_  = ( lhs_ + rhs_ ) / scalar;
            refres_ = ( reflhs_ + refrhs_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   = ( lhs_ + orhs_ ) / scalar;
            odres_  = ( lhs_ + orhs_ ) / scalar;
            sres_   = ( lhs_ + orhs_ ) / scalar;
            osres_  = ( lhs_ + orhs_ ) / scalar;
            refres_ = ( reflhs_ + refrhs_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   = ( olhs_ + rhs_ ) / scalar;
            odres_  = ( olhs_ + rhs_ ) / scalar;
            sres_   = ( olhs_ + rhs_ ) / scalar;
            osres_  = ( olhs_ + rhs_ ) / scalar;
            refres_ = ( reflhs_ + refrhs_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   = ( olhs_ + orhs_ ) / scalar;
            odres_  = ( olhs_ + orhs_ ) / scalar;
            sres_   = ( olhs_ + orhs_ ) / scalar;
            osres_  = ( olhs_ + orhs_ ) / scalar;
            refres_ = ( reflhs_ + refrhs_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Scaled addition with evaluated matrices
      {
         test_  = "Scaled addition with evaluated matrices (OP/s)";
         error_ = "Failed addition operation";

         try {
            initResults();
            dres_   = ( eval( lhs_ ) + eval( rhs_ ) ) / scalar;
            odres_  = ( eval( lhs_ ) + eval( rhs_ ) ) / scalar;
            sres_   = ( eval( lhs_ ) + eval( rhs_ ) ) / scalar;
            osres_  = ( eval( lhs_ ) + eval( rhs_ ) ) / scalar;
            refres_ = ( eval( reflhs_ ) + eval( refrhs_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   = ( eval( lhs_ ) + eval( orhs_ ) ) / scalar;
            odres_  = ( eval( lhs_ ) + eval( orhs_ ) ) / scalar;
            sres_   = ( eval( lhs_ ) + eval( orhs_ ) ) / scalar;
            osres_  = ( eval( lhs_ ) + eval( orhs_ ) ) / scalar;
            refres_ = ( eval( reflhs_ ) + eval( refrhs_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   = ( eval( olhs_ ) + eval( rhs_ ) ) / scalar;
            odres_  = ( eval( olhs_ ) + eval( rhs_ ) ) / scalar;
            sres_   = ( eval( olhs_ ) + eval( rhs_ ) ) / scalar;
            osres_  = ( eval( olhs_ ) + eval( rhs_ ) ) / scalar;
            refres_ = ( eval( reflhs_ ) + eval( refrhs_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   = ( eval( olhs_ ) + eval( orhs_ ) ) / scalar;
            odres_  = ( eval( olhs_ ) + eval( orhs_ ) ) / scalar;
            sres_   = ( eval( olhs_ ) + eval( orhs_ ) ) / scalar;
            osres_  = ( eval( olhs_ ) + eval( orhs_ ) ) / scalar;
            refres_ = ( eval( reflhs_ ) + eval( refrhs_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Scaled addition with addition assignment (s*OP)
      //=====================================================================================

      // Scaled addition with addition assignment with the given matrices
      {
         test_  = "Scaled addition with addition assignment with the given matrices (s*OP)";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += scalar * ( lhs_ + rhs_ );
            odres_  += scalar * ( lhs_ + rhs_ );
            sres_   += scalar * ( lhs_ + rhs_ );
            osres_  += scalar * ( lhs_ + rhs_ );
            refres_ += scalar * ( reflhs_ + refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   += scalar * ( lhs_ + orhs_ );
            odres_  += scalar * ( lhs_ + orhs_ );
            sres_   += scalar * ( lhs_ + orhs_ );
            osres_  += scalar * ( lhs_ + orhs_ );
            refres_ += scalar * ( reflhs_ + refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   += scalar * ( olhs_ + rhs_ );
            odres_  += scalar * ( olhs_ + rhs_ );
            sres_   += scalar * ( olhs_ + rhs_ );
            osres_  += scalar * ( olhs_ + rhs_ );
            refres_ += scalar * ( reflhs_ + refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   += scalar * ( olhs_ + orhs_ );
            odres_  += scalar * ( olhs_ + orhs_ );
            sres_   += scalar * ( olhs_ + orhs_ );
            osres_  += scalar * ( olhs_ + orhs_ );
            refres_ += scalar * ( reflhs_ + refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Scaled addition with addition assignment with evaluated matrices
      {
         test_  = "Scaled addition with addition assignment with evaluated matrices (s*OP)";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += scalar * ( eval( lhs_ ) + eval( rhs_ ) );
            odres_  += scalar * ( eval( lhs_ ) + eval( rhs_ ) );
            sres_   += scalar * ( eval( lhs_ ) + eval( rhs_ ) );
            osres_  += scalar * ( eval( lhs_ ) + eval( rhs_ ) );
            refres_ += scalar * ( eval( reflhs_ ) + eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   += scalar * ( eval( lhs_ ) + eval( orhs_ ) );
            odres_  += scalar * ( eval( lhs_ ) + eval( orhs_ ) );
            sres_   += scalar * ( eval( lhs_ ) + eval( orhs_ ) );
            osres_  += scalar * ( eval( lhs_ ) + eval( orhs_ ) );
            refres_ += scalar * ( eval( reflhs_ ) + eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   += scalar * ( eval( olhs_ ) + eval( rhs_ ) );
            odres_  += scalar * ( eval( olhs_ ) + eval( rhs_ ) );
            sres_   += scalar * ( eval( olhs_ ) + eval( rhs_ ) );
            osres_  += scalar * ( eval( olhs_ ) + eval( rhs_ ) );
            refres_ += scalar * ( eval( reflhs_ ) + eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   += scalar * ( eval( olhs_ ) + eval( orhs_ ) );
            odres_  += scalar * ( eval( olhs_ ) + eval( orhs_ ) );
            sres_   += scalar * ( eval( olhs_ ) + eval( orhs_ ) );
            osres_  += scalar * ( eval( olhs_ ) + eval( orhs_ ) );
            refres_ += scalar * ( eval( reflhs_ ) + eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Scaled addition with addition assignment (OP*s)
      //=====================================================================================

      // Scaled addition with addition assignment with the given matrices
      {
         test_  = "Scaled addition with addition assignment with the given matrices (OP*s)";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += ( lhs_ + rhs_ ) * scalar;
            odres_  += ( lhs_ + rhs_ ) * scalar;
            sres_   += ( lhs_ + rhs_ ) * scalar;
            osres_  += ( lhs_ + rhs_ ) * scalar;
            refres_ += ( reflhs_ + refrhs_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   += ( lhs_ + orhs_ ) * scalar;
            odres_  += ( lhs_ + orhs_ ) * scalar;
            sres_   += ( lhs_ + orhs_ ) * scalar;
            osres_  += ( lhs_ + orhs_ ) * scalar;
            refres_ += ( reflhs_ + refrhs_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   += ( olhs_ + rhs_ ) * scalar;
            odres_  += ( olhs_ + rhs_ ) * scalar;
            sres_   += ( olhs_ + rhs_ ) * scalar;
            osres_  += ( olhs_ + rhs_ ) * scalar;
            refres_ += ( reflhs_ + refrhs_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   += ( olhs_ + orhs_ ) * scalar;
            odres_  += ( olhs_ + orhs_ ) * scalar;
            sres_   += ( olhs_ + orhs_ ) * scalar;
            osres_  += ( olhs_ + orhs_ ) * scalar;
            refres_ += ( reflhs_ + refrhs_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Scaled addition with addition assignment with evaluated matrices
      {
         test_  = "Scaled addition with addition assignment with evaluated matrices (OP*s)";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += ( eval( lhs_ ) + eval( rhs_ ) ) * scalar;
            odres_  += ( eval( lhs_ ) + eval( rhs_ ) ) * scalar;
            sres_   += ( eval( lhs_ ) + eval( rhs_ ) ) * scalar;
            osres_  += ( eval( lhs_ ) + eval( rhs_ ) ) * scalar;
            refres_ += ( eval( reflhs_ ) + eval( refrhs_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   += ( eval( lhs_ ) + eval( orhs_ ) ) * scalar;
            odres_  += ( eval( lhs_ ) + eval( orhs_ ) ) * scalar;
            sres_   += ( eval( lhs_ ) + eval( orhs_ ) ) * scalar;
            osres_  += ( eval( lhs_ ) + eval( orhs_ ) ) * scalar;
            refres_ += ( eval( reflhs_ ) + eval( refrhs_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   += ( eval( olhs_ ) + eval( rhs_ ) ) * scalar;
            odres_  += ( eval( olhs_ ) + eval( rhs_ ) ) * scalar;
            sres_   += ( eval( olhs_ ) + eval( rhs_ ) ) * scalar;
            osres_  += ( eval( olhs_ ) + eval( rhs_ ) ) * scalar;
            refres_ += ( eval( reflhs_ ) + eval( refrhs_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   += ( eval( olhs_ ) + eval( orhs_ ) ) * scalar;
            odres_  += ( eval( olhs_ ) + eval( orhs_ ) ) * scalar;
            sres_   += ( eval( olhs_ ) + eval( orhs_ ) ) * scalar;
            osres_  += ( eval( olhs_ ) + eval( orhs_ ) ) * scalar;
            refres_ += ( eval( reflhs_ ) + eval( refrhs_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Scaled addition with addition assignment (OP/s)
      //=====================================================================================

      // Scaled addition with addition assignment with the given matrices
      {
         test_  = "Scaled addition with addition assignment with the given matrices (OP/s)";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += ( lhs_ + rhs_ ) / scalar;
            odres_  += ( lhs_ + rhs_ ) / scalar;
            sres_   += ( lhs_ + rhs_ ) / scalar;
            osres_  += ( lhs_ + rhs_ ) / scalar;
            refres_ += ( reflhs_ + refrhs_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   += ( lhs_ + orhs_ ) / scalar;
            odres_  += ( lhs_ + orhs_ ) / scalar;
            sres_   += ( lhs_ + orhs_ ) / scalar;
            osres_  += ( lhs_ + orhs_ ) / scalar;
            refres_ += ( reflhs_ + refrhs_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   += ( olhs_ + rhs_ ) / scalar;
            odres_  += ( olhs_ + rhs_ ) / scalar;
            sres_   += ( olhs_ + rhs_ ) / scalar;
            osres_  += ( olhs_ + rhs_ ) / scalar;
            refres_ += ( reflhs_ + refrhs_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   += ( olhs_ + orhs_ ) / scalar;
            odres_  += ( olhs_ + orhs_ ) / scalar;
            sres_   += ( olhs_ + orhs_ ) / scalar;
            osres_  += ( olhs_ + orhs_ ) / scalar;
            refres_ += ( reflhs_ + refrhs_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Scaled addition with addition assignment with evaluated matrices
      {
         test_  = "Scaled addition with addition assignment with evaluated matrices (OP/s)";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += ( eval( lhs_ ) + eval( rhs_ ) ) / scalar;
            odres_  += ( eval( lhs_ ) + eval( rhs_ ) ) / scalar;
            sres_   += ( eval( lhs_ ) + eval( rhs_ ) ) / scalar;
            osres_  += ( eval( lhs_ ) + eval( rhs_ ) ) / scalar;
            refres_ += ( eval( reflhs_ ) + eval( refrhs_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   += ( eval( lhs_ ) + eval( orhs_ ) ) / scalar;
            odres_  += ( eval( lhs_ ) + eval( orhs_ ) ) / scalar;
            sres_   += ( eval( lhs_ ) + eval( orhs_ ) ) / scalar;
            osres_  += ( eval( lhs_ ) + eval( orhs_ ) ) / scalar;
            refres_ += ( eval( reflhs_ ) + eval( refrhs_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   += ( eval( olhs_ ) + eval( rhs_ ) ) / scalar;
            odres_  += ( eval( olhs_ ) + eval( rhs_ ) ) / scalar;
            sres_   += ( eval( olhs_ ) + eval( rhs_ ) ) / scalar;
            osres_  += ( eval( olhs_ ) + eval( rhs_ ) ) / scalar;
            refres_ += ( eval( reflhs_ ) + eval( refrhs_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   += ( eval( olhs_ ) + eval( orhs_ ) ) / scalar;
            odres_  += ( eval( olhs_ ) + eval( orhs_ ) ) / scalar;
            sres_   += ( eval( olhs_ ) + eval( orhs_ ) ) / scalar;
            osres_  += ( eval( olhs_ ) + eval( orhs_ ) ) / scalar;
            refres_ += ( eval( reflhs_ ) + eval( refrhs_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Scaled addition with subtraction assignment (s*OP)
      //=====================================================================================

      // Scaled addition with subtraction assignment with the given matrices
      {
         test_  = "Scaled addition with subtraction assignment with the given matrices (s*OP)";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= scalar * ( lhs_ + rhs_ );
            odres_  -= scalar * ( lhs_ + rhs_ );
            sres_   -= scalar * ( lhs_ + rhs_ );
            osres_  -= scalar * ( lhs_ + rhs_ );
            refres_ -= scalar * ( reflhs_ + refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   -= scalar * ( lhs_ + orhs_ );
            odres_  -= scalar * ( lhs_ + orhs_ );
            sres_   -= scalar * ( lhs_ + orhs_ );
            osres_  -= scalar * ( lhs_ + orhs_ );
            refres_ -= scalar * ( reflhs_ + refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   -= scalar * ( olhs_ + rhs_ );
            odres_  -= scalar * ( olhs_ + rhs_ );
            sres_   -= scalar * ( olhs_ + rhs_ );
            osres_  -= scalar * ( olhs_ + rhs_ );
            refres_ -= scalar * ( reflhs_ + refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   -= scalar * ( olhs_ + orhs_ );
            odres_  -= scalar * ( olhs_ + orhs_ );
            sres_   -= scalar * ( olhs_ + orhs_ );
            osres_  -= scalar * ( olhs_ + orhs_ );
            refres_ -= scalar * ( reflhs_ + refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Scaled addition with subtraction assignment with evaluated matrices
      {
         test_  = "Scaled addition with subtraction assignment with evaluated matrices (s*OP)";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= scalar * ( eval( lhs_ ) + eval( rhs_ ) );
            odres_  -= scalar * ( eval( lhs_ ) + eval( rhs_ ) );
            sres_   -= scalar * ( eval( lhs_ ) + eval( rhs_ ) );
            osres_  -= scalar * ( eval( lhs_ ) + eval( rhs_ ) );
            refres_ -= scalar * ( eval( reflhs_ ) + eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   -= scalar * ( eval( lhs_ ) + eval( orhs_ ) );
            odres_  -= scalar * ( eval( lhs_ ) + eval( orhs_ ) );
            sres_   -= scalar * ( eval( lhs_ ) + eval( orhs_ ) );
            osres_  -= scalar * ( eval( lhs_ ) + eval( orhs_ ) );
            refres_ -= scalar * ( eval( reflhs_ ) + eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   -= scalar * ( eval( olhs_ ) + eval( rhs_ ) );
            odres_  -= scalar * ( eval( olhs_ ) + eval( rhs_ ) );
            sres_   -= scalar * ( eval( olhs_ ) + eval( rhs_ ) );
            osres_  -= scalar * ( eval( olhs_ ) + eval( rhs_ ) );
            refres_ -= scalar * ( eval( reflhs_ ) + eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   -= scalar * ( eval( olhs_ ) + eval( orhs_ ) );
            odres_  -= scalar * ( eval( olhs_ ) + eval( orhs_ ) );
            sres_   -= scalar * ( eval( olhs_ ) + eval( orhs_ ) );
            osres_  -= scalar * ( eval( olhs_ ) + eval( orhs_ ) );
            refres_ -= scalar * ( eval( reflhs_ ) + eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Scaled addition with subtraction assignment (OP*s)
      //=====================================================================================

      // Scaled addition with subtraction assignment with the given matrices
      {
         test_  = "Scaled addition with subtraction assignment with the given matrices (OP*s)";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= ( lhs_ + rhs_ ) * scalar;
            odres_  -= ( lhs_ + rhs_ ) * scalar;
            sres_   -= ( lhs_ + rhs_ ) * scalar;
            osres_  -= ( lhs_ + rhs_ ) * scalar;
            refres_ -= ( reflhs_ + refrhs_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   -= ( lhs_ + orhs_ ) * scalar;
            odres_  -= ( lhs_ + orhs_ ) * scalar;
            sres_   -= ( lhs_ + orhs_ ) * scalar;
            osres_  -= ( lhs_ + orhs_ ) * scalar;
            refres_ -= ( reflhs_ + refrhs_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   -= ( olhs_ + rhs_ ) * scalar;
            odres_  -= ( olhs_ + rhs_ ) * scalar;
            sres_   -= ( olhs_ + rhs_ ) * scalar;
            osres_  -= ( olhs_ + rhs_ ) * scalar;
            refres_ -= ( reflhs_ + refrhs_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   -= ( olhs_ + orhs_ ) * scalar;
            odres_  -= ( olhs_ + orhs_ ) * scalar;
            sres_   -= ( olhs_ + orhs_ ) * scalar;
            osres_  -= ( olhs_ + orhs_ ) * scalar;
            refres_ -= ( reflhs_ + refrhs_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Scaled addition with subtraction assignment with evaluated matrices
      {
         test_  = "Scaled addition with subtraction assignment with evaluated matrices (OP*s)";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= ( eval( lhs_ ) + eval( rhs_ ) ) * scalar;
            odres_  -= ( eval( lhs_ ) + eval( rhs_ ) ) * scalar;
            sres_   -= ( eval( lhs_ ) + eval( rhs_ ) ) * scalar;
            osres_  -= ( eval( lhs_ ) + eval( rhs_ ) ) * scalar;
            refres_ -= ( eval( reflhs_ ) + eval( refrhs_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   -= ( eval( lhs_ ) + eval( orhs_ ) ) * scalar;
            odres_  -= ( eval( lhs_ ) + eval( orhs_ ) ) * scalar;
            sres_   -= ( eval( lhs_ ) + eval( orhs_ ) ) * scalar;
            osres_  -= ( eval( lhs_ ) + eval( orhs_ ) ) * scalar;
            refres_ -= ( eval( reflhs_ ) + eval( refrhs_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   -= ( eval( olhs_ ) + eval( rhs_ ) ) * scalar;
            odres_  -= ( eval( olhs_ ) + eval( rhs_ ) ) * scalar;
            sres_   -= ( eval( olhs_ ) + eval( rhs_ ) ) * scalar;
            osres_  -= ( eval( olhs_ ) + eval( rhs_ ) ) * scalar;
            refres_ -= ( eval( reflhs_ ) + eval( refrhs_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   -= ( eval( olhs_ ) + eval( orhs_ ) ) * scalar;
            odres_  -= ( eval( olhs_ ) + eval( orhs_ ) ) * scalar;
            sres_   -= ( eval( olhs_ ) + eval( orhs_ ) ) * scalar;
            osres_  -= ( eval( olhs_ ) + eval( orhs_ ) ) * scalar;
            refres_ -= ( eval( reflhs_ ) + eval( refrhs_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Scaled addition with subtraction assignment (OP/s)
      //=====================================================================================

      // Scaled addition with subtraction assignment with the given matrices
      {
         test_  = "Scaled addition with subtraction assignment with the given matrices (OP/s)";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= ( lhs_ + rhs_ ) / scalar;
            odres_  -= ( lhs_ + rhs_ ) / scalar;
            sres_   -= ( lhs_ + rhs_ ) / scalar;
            osres_  -= ( lhs_ + rhs_ ) / scalar;
            refres_ -= ( reflhs_ + refrhs_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   -= ( lhs_ + orhs_ ) / scalar;
            odres_  -= ( lhs_ + orhs_ ) / scalar;
            sres_   -= ( lhs_ + orhs_ ) / scalar;
            osres_  -= ( lhs_ + orhs_ ) / scalar;
            refres_ -= ( reflhs_ + refrhs_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   -= ( olhs_ + rhs_ ) / scalar;
            odres_  -= ( olhs_ + rhs_ ) / scalar;
            sres_   -= ( olhs_ + rhs_ ) / scalar;
            osres_  -= ( olhs_ + rhs_ ) / scalar;
            refres_ -= ( reflhs_ + refrhs_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   -= ( olhs_ + orhs_ ) / scalar;
            odres_  -= ( olhs_ + orhs_ ) / scalar;
            sres_   -= ( olhs_ + orhs_ ) / scalar;
            osres_  -= ( olhs_ + orhs_ ) / scalar;
            refres_ -= ( reflhs_ + refrhs_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Scaled addition with subtraction assignment with evaluated matrices
      {
         test_  = "Scaled addition with subtraction assignment with evaluated matrices (OP/s)";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= ( eval( lhs_ ) + eval( rhs_ ) ) / scalar;
            odres_  -= ( eval( lhs_ ) + eval( rhs_ ) ) / scalar;
            sres_   -= ( eval( lhs_ ) + eval( rhs_ ) ) / scalar;
            osres_  -= ( eval( lhs_ ) + eval( rhs_ ) ) / scalar;
            refres_ -= ( eval( reflhs_ ) + eval( refrhs_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   -= ( eval( lhs_ ) + eval( orhs_ ) ) / scalar;
            odres_  -= ( eval( lhs_ ) + eval( orhs_ ) ) / scalar;
            sres_   -= ( eval( lhs_ ) + eval( orhs_ ) ) / scalar;
            osres_  -= ( eval( lhs_ ) + eval( orhs_ ) ) / scalar;
            refres_ -= ( eval( reflhs_ ) + eval( refrhs_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   -= ( eval( olhs_ ) + eval( rhs_ ) ) / scalar;
            odres_  -= ( eval( olhs_ ) + eval( rhs_ ) ) / scalar;
            sres_   -= ( eval( olhs_ ) + eval( rhs_ ) ) / scalar;
            osres_  -= ( eval( olhs_ ) + eval( rhs_ ) ) / scalar;
            refres_ -= ( eval( reflhs_ ) + eval( refrhs_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   -= ( eval( olhs_ ) + eval( orhs_ ) ) / scalar;
            odres_  -= ( eval( olhs_ ) + eval( orhs_ ) ) / scalar;
            sres_   -= ( eval( olhs_ ) + eval( orhs_ ) ) / scalar;
            osres_  -= ( eval( olhs_ ) + eval( orhs_ ) ) / scalar;
            refres_ -= ( eval( reflhs_ ) + eval( refrhs_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the transpose dense matrix/dense matrix addition.
//
// \return void
// \exception std::runtime_error Addition error detected.
//
// This function tests the transpose matrix addition with plain assignment. In case any error
// resulting from the addition or the subsequent assignment is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
void OperationTest<MT1,MT2>::testTransposeOperation()
{
#if BLAZETEST_MATHTEST_TEST_TRANSPOSE_OPERATION
   if( BLAZETEST_MATHTEST_TEST_TRANSPOSE_OPERATION > 1 )
   {
      //=====================================================================================
      // Transpose addition
      //=====================================================================================

      // Transpose addition with the given matrices
      {
         test_  = "Transpose addition with the given matrices";
         error_ = "Failed addition operation";

         try {
            initTransposeResults();
            tdres_  = trans( lhs_ + rhs_ );
            todres_ = trans( lhs_ + rhs_ );
            tsres_  = trans( lhs_ + rhs_ );
            tosres_ = trans( lhs_ + rhs_ );
            refres_ = trans( reflhs_ + refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkTransposeResults<MT1,MT2>();

         try {
            initTransposeResults();
            tdres_  = trans( lhs_ + orhs_ );
            todres_ = trans( lhs_ + orhs_ );
            tsres_  = trans( lhs_ + orhs_ );
            tosres_ = trans( lhs_ + orhs_ );
            refres_ = trans( reflhs_ + refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkTransposeResults<MT1,OMT2>();

         try {
            initTransposeResults();
            tdres_  = trans( olhs_ + rhs_ );
            todres_ = trans( olhs_ + rhs_ );
            tsres_  = trans( olhs_ + rhs_ );
            tosres_ = trans( olhs_ + rhs_ );
            refres_ = trans( reflhs_ + refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkTransposeResults<OMT1,MT2>();

         try {
            initTransposeResults();
            tdres_  = trans( olhs_ + orhs_ );
            todres_ = trans( olhs_ + orhs_ );
            tsres_  = trans( olhs_ + orhs_ );
            tosres_ = trans( olhs_ + orhs_ );
            refres_ = trans( reflhs_ + refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkTransposeResults<OMT1,OMT2>();
      }

      // Transpose addition with evaluated matrices
      {
         test_  = "Transpose addition with evaluated matrices";
         error_ = "Failed addition operation";

         try {
            initTransposeResults();
            tdres_  = trans( eval( lhs_ ) + eval( rhs_ ) );
            todres_ = trans( eval( lhs_ ) + eval( rhs_ ) );
            tsres_  = trans( eval( lhs_ ) + eval( rhs_ ) );
            tosres_ = trans( eval( lhs_ ) + eval( rhs_ ) );
            refres_ = trans( eval( reflhs_ ) + eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkTransposeResults<MT1,MT2>();

         try {
            initTransposeResults();
            tdres_  = trans( eval( lhs_ ) + eval( orhs_ ) );
            todres_ = trans( eval( lhs_ ) + eval( orhs_ ) );
            tsres_  = trans( eval( lhs_ ) + eval( orhs_ ) );
            tosres_ = trans( eval( lhs_ ) + eval( orhs_ ) );
            refres_ = trans( eval( reflhs_ ) + eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkTransposeResults<MT1,OMT2>();

         try {
            initTransposeResults();
            tdres_  = trans( eval( olhs_ ) + eval( rhs_ ) );
            todres_ = trans( eval( olhs_ ) + eval( rhs_ ) );
            tsres_  = trans( eval( olhs_ ) + eval( rhs_ ) );
            tosres_ = trans( eval( olhs_ ) + eval( rhs_ ) );
            refres_ = trans( eval( reflhs_ ) + eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkTransposeResults<OMT1,MT2>();

         try {
            initTransposeResults();
            tdres_  = trans( eval( olhs_ ) + eval( orhs_ ) );
            todres_ = trans( eval( olhs_ ) + eval( orhs_ ) );
            tsres_  = trans( eval( olhs_ ) + eval( orhs_ ) );
            tosres_ = trans( eval( olhs_ ) + eval( orhs_ ) );
            refres_ = trans( eval( reflhs_ ) + eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkTransposeResults<OMT1,OMT2>();
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the abs dense matrix/dense matrix addition.
//
// \return void
// \exception std::runtime_error Addition error detected.
//
// This function tests the abs matrix addition with plain assignment, addition assignment,
// and subtraction assignment. In case any error resulting from the addition or the subsequent
// assignment is detected, a \a std::runtime_error exception is thrown.
*/
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
void OperationTest<MT1,MT2>::testAbsOperation()
{
#if BLAZETEST_MATHTEST_TEST_ABS_OPERATION
   if( BLAZETEST_MATHTEST_TEST_ABS_OPERATION > 1 )
   {
      //=====================================================================================
      // Abs addition
      //=====================================================================================

      // Abs addition with the given matrices
      {
         test_  = "Abs addition with the given matrices";
         error_ = "Failed addition operation";

         try {
            initResults();
            dres_   = abs( lhs_ + rhs_ );
            odres_  = abs( lhs_ + rhs_ );
            sres_   = abs( lhs_ + rhs_ );
            osres_  = abs( lhs_ + rhs_ );
            refres_ = abs( reflhs_ + refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   = abs( lhs_ + orhs_ );
            odres_  = abs( lhs_ + orhs_ );
            sres_   = abs( lhs_ + orhs_ );
            osres_  = abs( lhs_ + orhs_ );
            refres_ = abs( reflhs_ + refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   = abs( olhs_ + rhs_ );
            odres_  = abs( olhs_ + rhs_ );
            sres_   = abs( olhs_ + rhs_ );
            osres_  = abs( olhs_ + rhs_ );
            refres_ = abs( reflhs_ + refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   = abs( olhs_ + orhs_ );
            odres_  = abs( olhs_ + orhs_ );
            sres_   = abs( olhs_ + orhs_ );
            osres_  = abs( olhs_ + orhs_ );
            refres_ = abs( reflhs_ + refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Abs addition with evaluated matrices
      {
         test_  = "Abs addition with evaluated matrices";
         error_ = "Failed addition operation";

         try {
            initResults();
            dres_   = abs( eval( lhs_ ) + eval( rhs_ ) );
            odres_  = abs( eval( lhs_ ) + eval( rhs_ ) );
            sres_   = abs( eval( lhs_ ) + eval( rhs_ ) );
            osres_  = abs( eval( lhs_ ) + eval( rhs_ ) );
            refres_ = abs( eval( reflhs_ ) + eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   = abs( eval( lhs_ ) + eval( orhs_ ) );
            odres_  = abs( eval( lhs_ ) + eval( orhs_ ) );
            sres_   = abs( eval( lhs_ ) + eval( orhs_ ) );
            osres_  = abs( eval( lhs_ ) + eval( orhs_ ) );
            refres_ = abs( eval( reflhs_ ) + eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   = abs( eval( olhs_ ) + eval( rhs_ ) );
            odres_  = abs( eval( olhs_ ) + eval( rhs_ ) );
            sres_   = abs( eval( olhs_ ) + eval( rhs_ ) );
            osres_  = abs( eval( olhs_ ) + eval( rhs_ ) );
            refres_ = abs( eval( reflhs_ ) + eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   = abs( eval( olhs_ ) + eval( orhs_ ) );
            odres_  = abs( eval( olhs_ ) + eval( orhs_ ) );
            sres_   = abs( eval( olhs_ ) + eval( orhs_ ) );
            osres_  = abs( eval( olhs_ ) + eval( orhs_ ) );
            refres_ = abs( eval( reflhs_ ) + eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Abs addition with addition assignment
      //=====================================================================================

      // Abs addition with addition assignment with the given matrices
      {
         test_  = "Abs addition with addition assignment with the given matrices";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += abs( lhs_ + rhs_ );
            odres_  += abs( lhs_ + rhs_ );
            sres_   += abs( lhs_ + rhs_ );
            osres_  += abs( lhs_ + rhs_ );
            refres_ += abs( reflhs_ + refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   += abs( lhs_ + orhs_ );
            odres_  += abs( lhs_ + orhs_ );
            sres_   += abs( lhs_ + orhs_ );
            osres_  += abs( lhs_ + orhs_ );
            refres_ += abs( reflhs_ + refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   += abs( olhs_ + rhs_ );
            odres_  += abs( olhs_ + rhs_ );
            sres_   += abs( olhs_ + rhs_ );
            osres_  += abs( olhs_ + rhs_ );
            refres_ += abs( reflhs_ + refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   += abs( olhs_ + orhs_ );
            odres_  += abs( olhs_ + orhs_ );
            sres_   += abs( olhs_ + orhs_ );
            osres_  += abs( olhs_ + orhs_ );
            refres_ += abs( reflhs_ + refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Abs addition with addition assignment with evaluated matrices
      {
         test_  = "Abs addition with addition assignment with evaluated matrices";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += abs( eval( lhs_ ) + eval( rhs_ ) );
            odres_  += abs( eval( lhs_ ) + eval( rhs_ ) );
            sres_   += abs( eval( lhs_ ) + eval( rhs_ ) );
            osres_  += abs( eval( lhs_ ) + eval( rhs_ ) );
            refres_ += abs( eval( reflhs_ ) + eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   += abs( eval( lhs_ ) + eval( orhs_ ) );
            odres_  += abs( eval( lhs_ ) + eval( orhs_ ) );
            sres_   += abs( eval( lhs_ ) + eval( orhs_ ) );
            osres_  += abs( eval( lhs_ ) + eval( orhs_ ) );
            refres_ += abs( eval( reflhs_ ) + eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   += abs( eval( olhs_ ) + eval( rhs_ ) );
            odres_  += abs( eval( olhs_ ) + eval( rhs_ ) );
            sres_   += abs( eval( olhs_ ) + eval( rhs_ ) );
            osres_  += abs( eval( olhs_ ) + eval( rhs_ ) );
            refres_ += abs( eval( reflhs_ ) + eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   += abs( eval( olhs_ ) + eval( orhs_ ) );
            odres_  += abs( eval( olhs_ ) + eval( orhs_ ) );
            sres_   += abs( eval( olhs_ ) + eval( orhs_ ) );
            osres_  += abs( eval( olhs_ ) + eval( orhs_ ) );
            refres_ += abs( eval( reflhs_ ) + eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Abs addition with subtraction assignment
      //=====================================================================================

      // Abs addition with subtraction assignment with the given matrices
      {
         test_  = "Abs addition with subtraction assignment with the given matrices";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= abs( lhs_ + rhs_ );
            odres_  -= abs( lhs_ + rhs_ );
            sres_   -= abs( lhs_ + rhs_ );
            osres_  -= abs( lhs_ + rhs_ );
            refres_ -= abs( reflhs_ + refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   -= abs( lhs_ + orhs_ );
            odres_  -= abs( lhs_ + orhs_ );
            sres_   -= abs( lhs_ + orhs_ );
            osres_  -= abs( lhs_ + orhs_ );
            refres_ -= abs( reflhs_ + refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   -= abs( olhs_ + rhs_ );
            odres_  -= abs( olhs_ + rhs_ );
            sres_   -= abs( olhs_ + rhs_ );
            osres_  -= abs( olhs_ + rhs_ );
            refres_ -= abs( reflhs_ + refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   -= abs( olhs_ + orhs_ );
            odres_  -= abs( olhs_ + orhs_ );
            sres_   -= abs( olhs_ + orhs_ );
            osres_  -= abs( olhs_ + orhs_ );
            refres_ -= abs( reflhs_ + refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Abs addition with subtraction assignment with evaluated matrices
      {
         test_  = "Abs addition with subtraction assignment with evaluated matrices";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= abs( eval( lhs_ ) + eval( rhs_ ) );
            odres_  -= abs( eval( lhs_ ) + eval( rhs_ ) );
            sres_   -= abs( eval( lhs_ ) + eval( rhs_ ) );
            osres_  -= abs( eval( lhs_ ) + eval( rhs_ ) );
            refres_ -= abs( eval( reflhs_ ) + eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   -= abs( eval( lhs_ ) + eval( orhs_ ) );
            odres_  -= abs( eval( lhs_ ) + eval( orhs_ ) );
            sres_   -= abs( eval( lhs_ ) + eval( orhs_ ) );
            osres_  -= abs( eval( lhs_ ) + eval( orhs_ ) );
            refres_ -= abs( eval( reflhs_ ) + eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   -= abs( eval( olhs_ ) + eval( rhs_ ) );
            odres_  -= abs( eval( olhs_ ) + eval( rhs_ ) );
            sres_   -= abs( eval( olhs_ ) + eval( rhs_ ) );
            osres_  -= abs( eval( olhs_ ) + eval( rhs_ ) );
            refres_ -= abs( eval( reflhs_ ) + eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   -= abs( eval( olhs_ ) + eval( orhs_ ) );
            odres_  -= abs( eval( olhs_ ) + eval( orhs_ ) );
            sres_   -= abs( eval( olhs_ ) + eval( orhs_ ) );
            osres_  -= abs( eval( olhs_ ) + eval( orhs_ ) );
            refres_ -= abs( eval( reflhs_ ) + eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the submatrix-wise dense matrix/dense matrix addition.
//
// \return void
// \exception std::runtime_error Addition error detected.
//
// This function tests the submatrix-wise matrix addition with plain assignment, addition
// assignment, and subtraction assignment. In case any error resulting from the addition
// or the subsequent assignment is detected, a \a std::runtime_error exception is thrown.
*/
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
void OperationTest<MT1,MT2>::testSubmatrixOperation()
{
#if BLAZETEST_MATHTEST_TEST_SUBMATRIX_OPERATION
   if( BLAZETEST_MATHTEST_TEST_SUBMATRIX_OPERATION > 1 )
   {
      if( lhs_.rows() == 0UL || lhs_.columns() == 0UL )
         return;


      //=====================================================================================
      // Submatrix-wise addition
      //=====================================================================================

      // Submatrix-wise addition with the given matrices
      {
         test_  = "Submatrix-wise addition with the given matrices";
         error_ = "Failed addition operation";

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.rows(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.rows() - row );
               for( size_t column=0UL, n=0UL; column<rhs_.columns(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, rhs_.columns() - column );
                  submatrix( dres_  , row, column, m, n ) = submatrix( lhs_ + rhs_      , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) = submatrix( lhs_ + rhs_      , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) = submatrix( lhs_ + rhs_      , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) = submatrix( lhs_ + rhs_      , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) = submatrix( reflhs_ + refrhs_, row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.rows(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.rows() - row );
               for( size_t column=0UL, n=0UL; column<orhs_.columns(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, orhs_.columns() - column );
                  submatrix( dres_  , row, column, m, n ) = submatrix( lhs_ + orhs_     , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) = submatrix( lhs_ + orhs_     , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) = submatrix( lhs_ + orhs_     , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) = submatrix( lhs_ + orhs_     , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) = submatrix( reflhs_ + refrhs_, row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.rows(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.rows() - row );
               for( size_t column=0UL, n=0UL; column<rhs_.columns(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, rhs_.columns() - column );
                  submatrix( dres_  , row, column, m, n ) = submatrix( olhs_ + rhs_     , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) = submatrix( olhs_ + rhs_     , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) = submatrix( olhs_ + rhs_     , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) = submatrix( olhs_ + rhs_     , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) = submatrix( reflhs_ + refrhs_, row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.rows(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.rows() - row );
               for( size_t column=0UL, n=0UL; column<orhs_.columns(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, orhs_.columns() - column );
                  submatrix( dres_  , row, column, m, n ) = submatrix( olhs_ + orhs_    , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) = submatrix( olhs_ + orhs_    , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) = submatrix( olhs_ + orhs_    , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) = submatrix( olhs_ + orhs_    , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) = submatrix( reflhs_ + refrhs_, row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Submatrix-wise addition with evaluated matrices
      {
         test_  = "Submatrix-wise addition with evaluated matrices";
         error_ = "Failed addition operation";

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.rows(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.rows() - row );
               for( size_t column=0UL, n=0UL; column<rhs_.columns(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, rhs_.columns() - column );
                  submatrix( dres_  , row, column, m, n ) = submatrix( eval( lhs_ ) + eval( rhs_ )      , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) = submatrix( eval( lhs_ ) + eval( rhs_ )      , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) = submatrix( eval( lhs_ ) + eval( rhs_ )      , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) = submatrix( eval( lhs_ ) + eval( rhs_ )      , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) = submatrix( eval( reflhs_ ) + eval( refrhs_ ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.rows(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.rows() - row );
               for( size_t column=0UL, n=0UL; column<orhs_.columns(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, orhs_.columns() - column );
                  submatrix( dres_  , row, column, m, n ) = submatrix( eval( lhs_ ) + eval( orhs_ )     , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) = submatrix( eval( lhs_ ) + eval( orhs_ )     , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) = submatrix( eval( lhs_ ) + eval( orhs_ )     , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) = submatrix( eval( lhs_ ) + eval( orhs_ )     , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) = submatrix( eval( reflhs_ ) + eval( refrhs_ ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.rows(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.rows() - row );
               for( size_t column=0UL, n=0UL; column<rhs_.columns(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, rhs_.columns() - column );
                  submatrix( dres_  , row, column, m, n ) = submatrix( eval( olhs_ ) + eval( rhs_ )     , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) = submatrix( eval( olhs_ ) + eval( rhs_ )     , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) = submatrix( eval( olhs_ ) + eval( rhs_ )     , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) = submatrix( eval( olhs_ ) + eval( rhs_ )     , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) = submatrix( eval( reflhs_ ) + eval( refrhs_ ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.rows(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.rows() - row );
               for( size_t column=0UL, n=0UL; column<orhs_.columns(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, orhs_.columns() - column );
                  submatrix( dres_  , row, column, m, n ) = submatrix( eval( olhs_ ) + eval( orhs_ )    , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) = submatrix( eval( olhs_ ) + eval( orhs_ )    , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) = submatrix( eval( olhs_ ) + eval( orhs_ )    , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) = submatrix( eval( olhs_ ) + eval( orhs_ )    , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) = submatrix( eval( reflhs_ ) + eval( refrhs_ ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Submatrix-wise addition with addition assignment
      //=====================================================================================

      // Submatrix-wise addition with addition assignment with the given matrices
      {
         test_  = "Submatrix-wise addition with addition assignment with the given matrices";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.rows(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.rows() - row );
               for( size_t column=0UL, n=0UL; column<rhs_.columns(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, rhs_.columns() - column );
                  submatrix( dres_  , row, column, m, n ) += submatrix( lhs_ + rhs_      , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) += submatrix( lhs_ + rhs_      , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) += submatrix( lhs_ + rhs_      , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) += submatrix( lhs_ + rhs_      , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) += submatrix( reflhs_ + refrhs_, row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.rows(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.rows() - row );
               for( size_t column=0UL, n=0UL; column<orhs_.columns(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, orhs_.columns() - column );
                  submatrix( dres_  , row, column, m, n ) += submatrix( lhs_ + orhs_     , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) += submatrix( lhs_ + orhs_     , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) += submatrix( lhs_ + orhs_     , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) += submatrix( lhs_ + orhs_     , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) += submatrix( reflhs_ + refrhs_, row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.rows(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.rows() - row );
               for( size_t column=0UL, n=0UL; column<rhs_.columns(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, rhs_.columns() - column );
                  submatrix( dres_  , row, column, m, n ) += submatrix( olhs_ + rhs_     , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) += submatrix( olhs_ + rhs_     , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) += submatrix( olhs_ + rhs_     , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) += submatrix( olhs_ + rhs_     , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) += submatrix( reflhs_ + refrhs_, row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.rows(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.rows() - row );
               for( size_t column=0UL, n=0UL; column<orhs_.columns(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, orhs_.columns() - column );
                  submatrix( dres_  , row, column, m, n ) += submatrix( olhs_ + orhs_    , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) += submatrix( olhs_ + orhs_    , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) += submatrix( olhs_ + orhs_    , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) += submatrix( olhs_ + orhs_    , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) += submatrix( reflhs_ + refrhs_, row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Submatrix-wise addition with addition assignment with evaluated matrices
      {
         test_  = "Submatrix-wise addition with addition assignment with evaluated matrices";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.rows(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.rows() - row );
               for( size_t column=0UL, n=0UL; column<rhs_.columns(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, rhs_.columns() - column );
                  submatrix( dres_  , row, column, m, n ) += submatrix( eval( lhs_ ) + eval( rhs_ )      , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) += submatrix( eval( lhs_ ) + eval( rhs_ )      , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) += submatrix( eval( lhs_ ) + eval( rhs_ )      , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) += submatrix( eval( lhs_ ) + eval( rhs_ )      , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) += submatrix( eval( reflhs_ ) + eval( refrhs_ ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.rows(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.rows() - row );
               for( size_t column=0UL, n=0UL; column<orhs_.columns(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, orhs_.columns() - column );
                  submatrix( dres_  , row, column, m, n ) += submatrix( eval( lhs_ ) + eval( orhs_ )     , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) += submatrix( eval( lhs_ ) + eval( orhs_ )     , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) += submatrix( eval( lhs_ ) + eval( orhs_ )     , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) += submatrix( eval( lhs_ ) + eval( orhs_ )     , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) += submatrix( eval( reflhs_ ) + eval( refrhs_ ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.rows(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.rows() - row );
               for( size_t column=0UL, n=0UL; column<rhs_.columns(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, rhs_.columns() - column );
                  submatrix( dres_  , row, column, m, n ) += submatrix( eval( olhs_ ) + eval( rhs_ )     , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) += submatrix( eval( olhs_ ) + eval( rhs_ )     , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) += submatrix( eval( olhs_ ) + eval( rhs_ )     , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) += submatrix( eval( olhs_ ) + eval( rhs_ )     , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) += submatrix( eval( reflhs_ ) + eval( refrhs_ ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.rows(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.rows() - row );
               for( size_t column=0UL, n=0UL; column<orhs_.columns(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, orhs_.columns() - column );
                  submatrix( dres_  , row, column, m, n ) += submatrix( eval( olhs_ ) + eval( orhs_ )    , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) += submatrix( eval( olhs_ ) + eval( orhs_ )    , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) += submatrix( eval( olhs_ ) + eval( orhs_ )    , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) += submatrix( eval( olhs_ ) + eval( orhs_ )    , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) += submatrix( eval( reflhs_ ) + eval( refrhs_ ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Submatrix-wise addition with subtraction assignment
      //=====================================================================================

      // Submatrix-wise addition with subtraction assignment with the given matrices
      {
         test_  = "Submatrix-wise addition with subtraction assignment with the given matrices";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.rows(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.rows() - row );
               for( size_t column=0UL, n=0UL; column<rhs_.columns(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, rhs_.columns() - column );
                  submatrix( dres_  , row, column, m, n ) -= submatrix( lhs_ + rhs_      , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) -= submatrix( lhs_ + rhs_      , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) -= submatrix( lhs_ + rhs_      , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) -= submatrix( lhs_ + rhs_      , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) -= submatrix( reflhs_ + refrhs_, row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.rows(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.rows() - row );
               for( size_t column=0UL, n=0UL; column<orhs_.columns(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, orhs_.columns() - column );
                  submatrix( dres_  , row, column, m, n ) -= submatrix( lhs_ + orhs_     , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) -= submatrix( lhs_ + orhs_     , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) -= submatrix( lhs_ + orhs_     , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) -= submatrix( lhs_ + orhs_     , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) -= submatrix( reflhs_ + refrhs_, row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.rows(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.rows() - row );
               for( size_t column=0UL, n=0UL; column<rhs_.columns(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, rhs_.columns() - column );
                  submatrix( dres_  , row, column, m, n ) -= submatrix( olhs_ + rhs_     , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) -= submatrix( olhs_ + rhs_     , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) -= submatrix( olhs_ + rhs_     , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) -= submatrix( olhs_ + rhs_     , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) -= submatrix( reflhs_ + refrhs_, row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.rows(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.rows() - row );
               for( size_t column=0UL, n=0UL; column<orhs_.columns(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, orhs_.columns() - column );
                  submatrix( dres_  , row, column, m, n ) -= submatrix( olhs_ + orhs_    , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) -= submatrix( olhs_ + orhs_    , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) -= submatrix( olhs_ + orhs_    , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) -= submatrix( olhs_ + orhs_    , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) -= submatrix( reflhs_ + refrhs_, row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Submatrix-wise addition with subtraction assignment with evaluated matrices
      {
         test_  = "Submatrix-wise addition with subtraction assignment with evaluated matrices";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.rows(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.rows() - row );
               for( size_t column=0UL, n=0UL; column<rhs_.columns(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, rhs_.columns() - column );
                  submatrix( dres_  , row, column, m, n ) -= submatrix( eval( lhs_ ) + eval( rhs_ )      , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) -= submatrix( eval( lhs_ ) + eval( rhs_ )      , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) -= submatrix( eval( lhs_ ) + eval( rhs_ )      , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) -= submatrix( eval( lhs_ ) + eval( rhs_ )      , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) -= submatrix( eval( reflhs_ ) + eval( refrhs_ ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.rows(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.rows() - row );
               for( size_t column=0UL, n=0UL; column<orhs_.columns(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, orhs_.columns() - column );
                  submatrix( dres_  , row, column, m, n ) -= submatrix( eval( lhs_ ) + eval( orhs_ )     , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) -= submatrix( eval( lhs_ ) + eval( orhs_ )     , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) -= submatrix( eval( lhs_ ) + eval( orhs_ )     , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) -= submatrix( eval( lhs_ ) + eval( orhs_ )     , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) -= submatrix( eval( reflhs_ ) + eval( refrhs_ ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.rows(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.rows() - row );
               for( size_t column=0UL, n=0UL; column<rhs_.columns(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, rhs_.columns() - column );
                  submatrix( dres_  , row, column, m, n ) -= submatrix( eval( olhs_ ) + eval( rhs_ )     , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) -= submatrix( eval( olhs_ ) + eval( rhs_ )     , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) -= submatrix( eval( olhs_ ) + eval( rhs_ )     , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) -= submatrix( eval( olhs_ ) + eval( rhs_ )     , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) -= submatrix( eval( reflhs_ ) + eval( refrhs_ ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.rows(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.rows() - row );
               for( size_t column=0UL, n=0UL; column<orhs_.columns(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, orhs_.columns() - column );
                  submatrix( dres_  , row, column, m, n ) -= submatrix( eval( olhs_ ) + eval( orhs_ )    , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) -= submatrix( eval( olhs_ ) + eval( orhs_ )    , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) -= submatrix( eval( olhs_ ) + eval( orhs_ )    , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) -= submatrix( eval( olhs_ ) + eval( orhs_ )    , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) -= submatrix( eval( reflhs_ ) + eval( refrhs_ ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the row-wise dense matrix/dense matrix addition.
//
// \return void
// \exception std::runtime_error Addition error detected.
//
// This function tests the row-wise matrix addition with plain assignment, addition assignment,
// and subtraction assignment. In case any error resulting from the addition or the subsequent
// assignment is detected, a \a std::runtime_error exception is thrown.
*/
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
void OperationTest<MT1,MT2>::testRowOperation()
{
#if BLAZETEST_MATHTEST_TEST_ROW_OPERATION
   if( BLAZETEST_MATHTEST_TEST_ROW_OPERATION > 1 )
   {
      if( lhs_.rows() == 0UL )
         return;


      //=====================================================================================
      // Row-wise addition
      //=====================================================================================

      // Row-wise addition with the given matrices
      {
         test_  = "Row-wise addition with the given matrices";
         error_ = "Failed addition operation";

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) = row( lhs_ + rhs_, i );
               row( odres_ , i ) = row( lhs_ + rhs_, i );
               row( sres_  , i ) = row( lhs_ + rhs_, i );
               row( osres_ , i ) = row( lhs_ + rhs_, i );
               row( refres_, i ) = row( reflhs_ + refrhs_, i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) = row( lhs_ + orhs_, i );
               row( odres_ , i ) = row( lhs_ + orhs_, i );
               row( sres_  , i ) = row( lhs_ + orhs_, i );
               row( osres_ , i ) = row( lhs_ + orhs_, i );
               row( refres_, i ) = row( reflhs_ + refrhs_, i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) = row( olhs_ + rhs_, i );
               row( odres_ , i ) = row( olhs_ + rhs_, i );
               row( sres_  , i ) = row( olhs_ + rhs_, i );
               row( osres_ , i ) = row( olhs_ + rhs_, i );
               row( refres_, i ) = row( reflhs_ + refrhs_, i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) = row( olhs_ + orhs_, i );
               row( odres_ , i ) = row( olhs_ + orhs_, i );
               row( sres_  , i ) = row( olhs_ + orhs_, i );
               row( osres_ , i ) = row( olhs_ + orhs_, i );
               row( refres_, i ) = row( reflhs_ + refrhs_, i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Row-wise addition with evaluated matrices
      {
         test_  = "Row-wise addition with evaluated matrices";
         error_ = "Failed addition operation";

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) = row( eval( lhs_ ) + eval( rhs_ ), i );
               row( odres_ , i ) = row( eval( lhs_ ) + eval( rhs_ ), i );
               row( sres_  , i ) = row( eval( lhs_ ) + eval( rhs_ ), i );
               row( osres_ , i ) = row( eval( lhs_ ) + eval( rhs_ ), i );
               row( refres_, i ) = row( eval( reflhs_ ) + eval( refrhs_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) = row( eval( lhs_ ) + eval( orhs_ ), i );
               row( odres_ , i ) = row( eval( lhs_ ) + eval( orhs_ ), i );
               row( sres_  , i ) = row( eval( lhs_ ) + eval( orhs_ ), i );
               row( osres_ , i ) = row( eval( lhs_ ) + eval( orhs_ ), i );
               row( refres_, i ) = row( eval( reflhs_ ) + eval( refrhs_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) = row( eval( olhs_ ) + eval( rhs_ ), i );
               row( odres_ , i ) = row( eval( olhs_ ) + eval( rhs_ ), i );
               row( sres_  , i ) = row( eval( olhs_ ) + eval( rhs_ ), i );
               row( osres_ , i ) = row( eval( olhs_ ) + eval( rhs_ ), i );
               row( refres_, i ) = row( eval( reflhs_ ) + eval( refrhs_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) = row( eval( olhs_ ) + eval( orhs_ ), i );
               row( odres_ , i ) = row( eval( olhs_ ) + eval( orhs_ ), i );
               row( sres_  , i ) = row( eval( olhs_ ) + eval( orhs_ ), i );
               row( osres_ , i ) = row( eval( olhs_ ) + eval( orhs_ ), i );
               row( refres_, i ) = row( eval( reflhs_ ) + eval( refrhs_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Row-wise addition with addition assignment
      //=====================================================================================

      // Row-wise addition with addition assignment with the given matrices
      {
         test_  = "Row-wise addition with addition assignment with the given matrices";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) += row( lhs_ + rhs_, i );
               row( odres_ , i ) += row( lhs_ + rhs_, i );
               row( sres_  , i ) += row( lhs_ + rhs_, i );
               row( osres_ , i ) += row( lhs_ + rhs_, i );
               row( refres_, i ) += row( reflhs_ + refrhs_, i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) += row( lhs_ + orhs_, i );
               row( odres_ , i ) += row( lhs_ + orhs_, i );
               row( sres_  , i ) += row( lhs_ + orhs_, i );
               row( osres_ , i ) += row( lhs_ + orhs_, i );
               row( refres_, i ) += row( reflhs_ + refrhs_, i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) += row( olhs_ + rhs_, i );
               row( odres_ , i ) += row( olhs_ + rhs_, i );
               row( sres_  , i ) += row( olhs_ + rhs_, i );
               row( osres_ , i ) += row( olhs_ + rhs_, i );
               row( refres_, i ) += row( reflhs_ + refrhs_, i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) += row( olhs_ + orhs_, i );
               row( odres_ , i ) += row( olhs_ + orhs_, i );
               row( sres_  , i ) += row( olhs_ + orhs_, i );
               row( osres_ , i ) += row( olhs_ + orhs_, i );
               row( refres_, i ) += row( reflhs_ + refrhs_, i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Row-wise addition with addition assignment with evaluated matrices
      {
         test_  = "Row-wise addition with addition assignment with evaluated matrices";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) += row( eval( lhs_ ) + eval( rhs_ ), i );
               row( odres_ , i ) += row( eval( lhs_ ) + eval( rhs_ ), i );
               row( sres_  , i ) += row( eval( lhs_ ) + eval( rhs_ ), i );
               row( osres_ , i ) += row( eval( lhs_ ) + eval( rhs_ ), i );
               row( refres_, i ) += row( eval( reflhs_ ) + eval( refrhs_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) += row( eval( lhs_ ) + eval( orhs_ ), i );
               row( odres_ , i ) += row( eval( lhs_ ) + eval( orhs_ ), i );
               row( sres_  , i ) += row( eval( lhs_ ) + eval( orhs_ ), i );
               row( osres_ , i ) += row( eval( lhs_ ) + eval( orhs_ ), i );
               row( refres_, i ) += row( eval( reflhs_ ) + eval( refrhs_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) += row( eval( olhs_ ) + eval( rhs_ ), i );
               row( odres_ , i ) += row( eval( olhs_ ) + eval( rhs_ ), i );
               row( sres_  , i ) += row( eval( olhs_ ) + eval( rhs_ ), i );
               row( osres_ , i ) += row( eval( olhs_ ) + eval( rhs_ ), i );
               row( refres_, i ) += row( eval( reflhs_ ) + eval( refrhs_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) += row( eval( olhs_ ) + eval( orhs_ ), i );
               row( odres_ , i ) += row( eval( olhs_ ) + eval( orhs_ ), i );
               row( sres_  , i ) += row( eval( olhs_ ) + eval( orhs_ ), i );
               row( osres_ , i ) += row( eval( olhs_ ) + eval( orhs_ ), i );
               row( refres_, i ) += row( eval( reflhs_ ) + eval( refrhs_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Row-wise addition with subtraction assignment
      //=====================================================================================

      // Row-wise addition with subtraction assignment with the given matrices
      {
         test_  = "Row-wise addition with subtraction assignment with the given matrices";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) -= row( lhs_ + rhs_, i );
               row( odres_ , i ) -= row( lhs_ + rhs_, i );
               row( sres_  , i ) -= row( lhs_ + rhs_, i );
               row( osres_ , i ) -= row( lhs_ + rhs_, i );
               row( refres_, i ) -= row( reflhs_ + refrhs_, i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) -= row( lhs_ + orhs_, i );
               row( odres_ , i ) -= row( lhs_ + orhs_, i );
               row( sres_  , i ) -= row( lhs_ + orhs_, i );
               row( osres_ , i ) -= row( lhs_ + orhs_, i );
               row( refres_, i ) -= row( reflhs_ + refrhs_, i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) -= row( olhs_ + rhs_, i );
               row( odres_ , i ) -= row( olhs_ + rhs_, i );
               row( sres_  , i ) -= row( olhs_ + rhs_, i );
               row( osres_ , i ) -= row( olhs_ + rhs_, i );
               row( refres_, i ) -= row( reflhs_ + refrhs_, i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) -= row( olhs_ + orhs_, i );
               row( odres_ , i ) -= row( olhs_ + orhs_, i );
               row( sres_  , i ) -= row( olhs_ + orhs_, i );
               row( osres_ , i ) -= row( olhs_ + orhs_, i );
               row( refres_, i ) -= row( reflhs_ + refrhs_, i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Row-wise addition with subtraction assignment with evaluated matrices
      {
         test_  = "Row-wise addition with subtraction assignment with evaluated matrices";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) -= row( eval( lhs_ ) + eval( rhs_ ), i );
               row( odres_ , i ) -= row( eval( lhs_ ) + eval( rhs_ ), i );
               row( sres_  , i ) -= row( eval( lhs_ ) + eval( rhs_ ), i );
               row( osres_ , i ) -= row( eval( lhs_ ) + eval( rhs_ ), i );
               row( refres_, i ) -= row( eval( reflhs_ ) + eval( refrhs_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) -= row( eval( lhs_ ) + eval( orhs_ ), i );
               row( odres_ , i ) -= row( eval( lhs_ ) + eval( orhs_ ), i );
               row( sres_  , i ) -= row( eval( lhs_ ) + eval( orhs_ ), i );
               row( osres_ , i ) -= row( eval( lhs_ ) + eval( orhs_ ), i );
               row( refres_, i ) -= row( eval( reflhs_ ) + eval( refrhs_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) -= row( eval( olhs_ ) + eval( rhs_ ), i );
               row( odres_ , i ) -= row( eval( olhs_ ) + eval( rhs_ ), i );
               row( sres_  , i ) -= row( eval( olhs_ ) + eval( rhs_ ), i );
               row( osres_ , i ) -= row( eval( olhs_ ) + eval( rhs_ ), i );
               row( refres_, i ) -= row( eval( reflhs_ ) + eval( refrhs_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) -= row( eval( olhs_ ) + eval( orhs_ ), i );
               row( odres_ , i ) -= row( eval( olhs_ ) + eval( orhs_ ), i );
               row( sres_  , i ) -= row( eval( olhs_ ) + eval( orhs_ ), i );
               row( osres_ , i ) -= row( eval( olhs_ ) + eval( orhs_ ), i );
               row( refres_, i ) -= row( eval( reflhs_ ) + eval( refrhs_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Row-wise addition with multiplication assignment
      //=====================================================================================

      // Row-wise addition with multiplication assignment with the given matrices
      {
         test_  = "Row-wise addition with multiplication assignment with the given matrices";
         error_ = "Failed multiplication assignment operation";

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) *= row( lhs_ + rhs_, i );
               row( odres_ , i ) *= row( lhs_ + rhs_, i );
               row( sres_  , i ) *= row( lhs_ + rhs_, i );
               row( osres_ , i ) *= row( lhs_ + rhs_, i );
               row( refres_, i ) *= row( reflhs_ + refrhs_, i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) *= row( lhs_ + orhs_, i );
               row( odres_ , i ) *= row( lhs_ + orhs_, i );
               row( sres_  , i ) *= row( lhs_ + orhs_, i );
               row( osres_ , i ) *= row( lhs_ + orhs_, i );
               row( refres_, i ) *= row( reflhs_ + refrhs_, i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) *= row( olhs_ + rhs_, i );
               row( odres_ , i ) *= row( olhs_ + rhs_, i );
               row( sres_  , i ) *= row( olhs_ + rhs_, i );
               row( osres_ , i ) *= row( olhs_ + rhs_, i );
               row( refres_, i ) *= row( reflhs_ + refrhs_, i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) *= row( olhs_ + orhs_, i );
               row( odres_ , i ) *= row( olhs_ + orhs_, i );
               row( sres_  , i ) *= row( olhs_ + orhs_, i );
               row( osres_ , i ) *= row( olhs_ + orhs_, i );
               row( refres_, i ) *= row( reflhs_ + refrhs_, i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Row-wise addition with multiplication assignment with evaluated matrices
      {
         test_  = "Row-wise addition with multiplication assignment with evaluated matrices";
         error_ = "Failed multiplication assignment operation";

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) *= row( eval( lhs_ ) + eval( rhs_ ), i );
               row( odres_ , i ) *= row( eval( lhs_ ) + eval( rhs_ ), i );
               row( sres_  , i ) *= row( eval( lhs_ ) + eval( rhs_ ), i );
               row( osres_ , i ) *= row( eval( lhs_ ) + eval( rhs_ ), i );
               row( refres_, i ) *= row( eval( reflhs_ ) + eval( refrhs_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) *= row( eval( lhs_ ) + eval( orhs_ ), i );
               row( odres_ , i ) *= row( eval( lhs_ ) + eval( orhs_ ), i );
               row( sres_  , i ) *= row( eval( lhs_ ) + eval( orhs_ ), i );
               row( osres_ , i ) *= row( eval( lhs_ ) + eval( orhs_ ), i );
               row( refres_, i ) *= row( eval( reflhs_ ) + eval( refrhs_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) *= row( eval( olhs_ ) + eval( rhs_ ), i );
               row( odres_ , i ) *= row( eval( olhs_ ) + eval( rhs_ ), i );
               row( sres_  , i ) *= row( eval( olhs_ ) + eval( rhs_ ), i );
               row( osres_ , i ) *= row( eval( olhs_ ) + eval( rhs_ ), i );
               row( refres_, i ) *= row( eval( reflhs_ ) + eval( refrhs_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) *= row( eval( olhs_ ) + eval( orhs_ ), i );
               row( odres_ , i ) *= row( eval( olhs_ ) + eval( orhs_ ), i );
               row( sres_  , i ) *= row( eval( olhs_ ) + eval( orhs_ ), i );
               row( osres_ , i ) *= row( eval( olhs_ ) + eval( orhs_ ), i );
               row( refres_, i ) *= row( eval( reflhs_ ) + eval( refrhs_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the column-wise dense matrix/dense matrix addition.
//
// \return void
// \exception std::runtime_error Addition error detected.
//
// This function tests the column-wise matrix addition with plain assignment, addition assignment,
// and subtraction assignment. In case any error resulting from the addition or the subsequent
// assignment is detected, a \a std::runtime_error exception is thrown.
*/
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
void OperationTest<MT1,MT2>::testColumnOperation()
{
#if BLAZETEST_MATHTEST_TEST_COLUMN_OPERATION
   if( BLAZETEST_MATHTEST_TEST_COLUMN_OPERATION > 1 )
   {
      if( lhs_.columns() == 0UL )
         return;


      //=====================================================================================
      // Column-wise addition
      //=====================================================================================

      // Column-wise addition with the given matrices
      {
         test_  = "Column-wise addition with the given matrices";
         error_ = "Failed addition operation";

         try {
            initResults();
            for( size_t j=0UL; j<lhs_.columns(); ++j ) {
               column( dres_  , j ) = column( lhs_ + rhs_, j );
               column( odres_ , j ) = column( lhs_ + rhs_, j );
               column( sres_  , j ) = column( lhs_ + rhs_, j );
               column( osres_ , j ) = column( lhs_ + rhs_, j );
               column( refres_, j ) = column( reflhs_ + refrhs_, j );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t j=0UL; j<lhs_.columns(); ++j ) {
               column( dres_  , j ) = column( lhs_ + orhs_, j );
               column( odres_ , j ) = column( lhs_ + orhs_, j );
               column( sres_  , j ) = column( lhs_ + orhs_, j );
               column( osres_ , j ) = column( lhs_ + orhs_, j );
               column( refres_, j ) = column( reflhs_ + refrhs_, j );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t j=0UL; j<lhs_.columns(); ++j ) {
               column( dres_  , j ) = column( olhs_ + rhs_, j );
               column( odres_ , j ) = column( olhs_ + rhs_, j );
               column( sres_  , j ) = column( olhs_ + rhs_, j );
               column( osres_ , j ) = column( olhs_ + rhs_, j );
               column( refres_, j ) = column( reflhs_ + refrhs_, j );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t j=0UL; j<lhs_.columns(); ++j ) {
               column( dres_  , j ) = column( olhs_ + orhs_, j );
               column( odres_ , j ) = column( olhs_ + orhs_, j );
               column( sres_  , j ) = column( olhs_ + orhs_, j );
               column( osres_ , j ) = column( olhs_ + orhs_, j );
               column( refres_, j ) = column( reflhs_ + refrhs_, j );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Column-wise addition with evaluated matrices
      {
         test_  = "Column-wise addition with evaluated matrices";
         error_ = "Failed addition operation";

         try {
            initResults();
            for( size_t j=0UL; j<lhs_.columns(); ++j ) {
               column( dres_  , j ) = column( eval( lhs_ ) + eval( rhs_ ), j );
               column( odres_ , j ) = column( eval( lhs_ ) + eval( rhs_ ), j );
               column( sres_  , j ) = column( eval( lhs_ ) + eval( rhs_ ), j );
               column( osres_ , j ) = column( eval( lhs_ ) + eval( rhs_ ), j );
               column( refres_, j ) = column( eval( reflhs_ ) + eval( refrhs_ ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t j=0UL; j<lhs_.columns(); ++j ) {
               column( dres_  , j ) = column( eval( lhs_ ) + eval( orhs_ ), j );
               column( odres_ , j ) = column( eval( lhs_ ) + eval( orhs_ ), j );
               column( sres_  , j ) = column( eval( lhs_ ) + eval( orhs_ ), j );
               column( osres_ , j ) = column( eval( lhs_ ) + eval( orhs_ ), j );
               column( refres_, j ) = column( eval( reflhs_ ) + eval( refrhs_ ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t j=0UL; j<lhs_.columns(); ++j ) {
               column( dres_  , j ) = column( eval( olhs_ ) + eval( rhs_ ), j );
               column( odres_ , j ) = column( eval( olhs_ ) + eval( rhs_ ), j );
               column( sres_  , j ) = column( eval( olhs_ ) + eval( rhs_ ), j );
               column( osres_ , j ) = column( eval( olhs_ ) + eval( rhs_ ), j );
               column( refres_, j ) = column( eval( reflhs_ ) + eval( refrhs_ ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t j=0UL; j<lhs_.columns(); ++j ) {
               column( dres_  , j ) = column( eval( olhs_ ) + eval( orhs_ ), j );
               column( odres_ , j ) = column( eval( olhs_ ) + eval( orhs_ ), j );
               column( sres_  , j ) = column( eval( olhs_ ) + eval( orhs_ ), j );
               column( osres_ , j ) = column( eval( olhs_ ) + eval( orhs_ ), j );
               column( refres_, j ) = column( eval( reflhs_ ) + eval( refrhs_ ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Column-wise addition with addition assignment
      //=====================================================================================

      // Column-wise addition with addition assignment with the given matrices
      {
         test_  = "Column-wise addition with addition assignment with the given matrices";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            for( size_t j=0UL; j<lhs_.columns(); ++j ) {
               column( dres_  , j ) += column( lhs_ + rhs_, j );
               column( odres_ , j ) += column( lhs_ + rhs_, j );
               column( sres_  , j ) += column( lhs_ + rhs_, j );
               column( osres_ , j ) += column( lhs_ + rhs_, j );
               column( refres_, j ) += column( reflhs_ + refrhs_, j );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t j=0UL; j<lhs_.columns(); ++j ) {
               column( dres_  , j ) += column( lhs_ + orhs_, j );
               column( odres_ , j ) += column( lhs_ + orhs_, j );
               column( sres_  , j ) += column( lhs_ + orhs_, j );
               column( osres_ , j ) += column( lhs_ + orhs_, j );
               column( refres_, j ) += column( reflhs_ + refrhs_, j );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t j=0UL; j<lhs_.columns(); ++j ) {
               column( dres_  , j ) += column( olhs_ + rhs_, j );
               column( odres_ , j ) += column( olhs_ + rhs_, j );
               column( sres_  , j ) += column( olhs_ + rhs_, j );
               column( osres_ , j ) += column( olhs_ + rhs_, j );
               column( refres_, j ) += column( reflhs_ + refrhs_, j );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t j=0UL; j<lhs_.columns(); ++j ) {
               column( dres_  , j ) += column( olhs_ + orhs_, j );
               column( odres_ , j ) += column( olhs_ + orhs_, j );
               column( sres_  , j ) += column( olhs_ + orhs_, j );
               column( osres_ , j ) += column( olhs_ + orhs_, j );
               column( refres_, j ) += column( reflhs_ + refrhs_, j );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Column-wise addition with addition assignment with evaluated matrices
      {
         test_  = "Column-wise addition with addition assignment with evaluated matrices";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            for( size_t j=0UL; j<lhs_.columns(); ++j ) {
               column( dres_  , j ) += column( eval( lhs_ ) + eval( rhs_ ), j );
               column( odres_ , j ) += column( eval( lhs_ ) + eval( rhs_ ), j );
               column( sres_  , j ) += column( eval( lhs_ ) + eval( rhs_ ), j );
               column( osres_ , j ) += column( eval( lhs_ ) + eval( rhs_ ), j );
               column( refres_, j ) += column( eval( reflhs_ ) + eval( refrhs_ ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t j=0UL; j<lhs_.columns(); ++j ) {
               column( dres_  , j ) += column( eval( lhs_ ) + eval( orhs_ ), j );
               column( odres_ , j ) += column( eval( lhs_ ) + eval( orhs_ ), j );
               column( sres_  , j ) += column( eval( lhs_ ) + eval( orhs_ ), j );
               column( osres_ , j ) += column( eval( lhs_ ) + eval( orhs_ ), j );
               column( refres_, j ) += column( eval( reflhs_ ) + eval( refrhs_ ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t j=0UL; j<lhs_.columns(); ++j ) {
               column( dres_  , j ) += column( eval( olhs_ ) + eval( rhs_ ), j );
               column( odres_ , j ) += column( eval( olhs_ ) + eval( rhs_ ), j );
               column( sres_  , j ) += column( eval( olhs_ ) + eval( rhs_ ), j );
               column( osres_ , j ) += column( eval( olhs_ ) + eval( rhs_ ), j );
               column( refres_, j ) += column( eval( reflhs_ ) + eval( refrhs_ ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t j=0UL; j<lhs_.columns(); ++j ) {
               column( dres_  , j ) += column( eval( olhs_ ) + eval( orhs_ ), j );
               column( odres_ , j ) += column( eval( olhs_ ) + eval( orhs_ ), j );
               column( sres_  , j ) += column( eval( olhs_ ) + eval( orhs_ ), j );
               column( osres_ , j ) += column( eval( olhs_ ) + eval( orhs_ ), j );
               column( refres_, j ) += column( eval( reflhs_ ) + eval( refrhs_ ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Column-wise addition with subtraction assignment
      //=====================================================================================

      // Column-wise addition with subtraction assignment with the given matrices
      {
         test_  = "Column-wise addition with subtraction assignment with the given matrices";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            for( size_t j=0UL; j<lhs_.columns(); ++j ) {
               column( dres_  , j ) -= column( lhs_ + rhs_, j );
               column( odres_ , j ) -= column( lhs_ + rhs_, j );
               column( sres_  , j ) -= column( lhs_ + rhs_, j );
               column( osres_ , j ) -= column( lhs_ + rhs_, j );
               column( refres_, j ) -= column( reflhs_ + refrhs_, j );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t j=0UL; j<lhs_.columns(); ++j ) {
               column( dres_  , j ) -= column( lhs_ + orhs_, j );
               column( odres_ , j ) -= column( lhs_ + orhs_, j );
               column( sres_  , j ) -= column( lhs_ + orhs_, j );
               column( osres_ , j ) -= column( lhs_ + orhs_, j );
               column( refres_, j ) -= column( reflhs_ + refrhs_, j );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t j=0UL; j<lhs_.columns(); ++j ) {
               column( dres_  , j ) -= column( olhs_ + rhs_, j );
               column( odres_ , j ) -= column( olhs_ + rhs_, j );
               column( sres_  , j ) -= column( olhs_ + rhs_, j );
               column( osres_ , j ) -= column( olhs_ + rhs_, j );
               column( refres_, j ) -= column( reflhs_ + refrhs_, j );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t j=0UL; j<lhs_.columns(); ++j ) {
               column( dres_  , j ) -= column( olhs_ + orhs_, j );
               column( odres_ , j ) -= column( olhs_ + orhs_, j );
               column( sres_  , j ) -= column( olhs_ + orhs_, j );
               column( osres_ , j ) -= column( olhs_ + orhs_, j );
               column( refres_, j ) -= column( reflhs_ + refrhs_, j );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Column-wise addition with subtraction assignment with evaluated matrices
      {
         test_  = "Column-wise addition with subtraction assignment with evaluated matrices";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            for( size_t j=0UL; j<lhs_.columns(); ++j ) {
               column( dres_  , j ) -= column( eval( lhs_ ) + eval( rhs_ ), j );
               column( odres_ , j ) -= column( eval( lhs_ ) + eval( rhs_ ), j );
               column( sres_  , j ) -= column( eval( lhs_ ) + eval( rhs_ ), j );
               column( osres_ , j ) -= column( eval( lhs_ ) + eval( rhs_ ), j );
               column( refres_, j ) -= column( eval( reflhs_ ) + eval( refrhs_ ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t j=0UL; j<lhs_.columns(); ++j ) {
               column( dres_  , j ) -= column( eval( lhs_ ) + eval( orhs_ ), j );
               column( odres_ , j ) -= column( eval( lhs_ ) + eval( orhs_ ), j );
               column( sres_  , j ) -= column( eval( lhs_ ) + eval( orhs_ ), j );
               column( osres_ , j ) -= column( eval( lhs_ ) + eval( orhs_ ), j );
               column( refres_, j ) -= column( eval( reflhs_ ) + eval( refrhs_ ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t j=0UL; j<lhs_.columns(); ++j ) {
               column( dres_  , j ) -= column( eval( olhs_ ) + eval( rhs_ ), j );
               column( odres_ , j ) -= column( eval( olhs_ ) + eval( rhs_ ), j );
               column( sres_  , j ) -= column( eval( olhs_ ) + eval( rhs_ ), j );
               column( osres_ , j ) -= column( eval( olhs_ ) + eval( rhs_ ), j );
               column( refres_, j ) -= column( eval( reflhs_ ) + eval( refrhs_ ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t j=0UL; j<lhs_.columns(); ++j ) {
               column( dres_  , j ) -= column( eval( olhs_ ) + eval( orhs_ ), j );
               column( odres_ , j ) -= column( eval( olhs_ ) + eval( orhs_ ), j );
               column( sres_  , j ) -= column( eval( olhs_ ) + eval( orhs_ ), j );
               column( osres_ , j ) -= column( eval( olhs_ ) + eval( orhs_ ), j );
               column( refres_, j ) -= column( eval( reflhs_ ) + eval( refrhs_ ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Column-wise addition with multiplication assignment
      //=====================================================================================

      // Column-wise addition with multiplication assignment with the given matrices
      {
         test_  = "Column-wise addition with multiplication assignment with the given matrices";
         error_ = "Failed multiplication assignment operation";

         try {
            initResults();
            for( size_t j=0UL; j<lhs_.columns(); ++j ) {
               column( dres_  , j ) *= column( lhs_ + rhs_, j );
               column( odres_ , j ) *= column( lhs_ + rhs_, j );
               column( sres_  , j ) *= column( lhs_ + rhs_, j );
               column( osres_ , j ) *= column( lhs_ + rhs_, j );
               column( refres_, j ) *= column( reflhs_ + refrhs_, j );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t j=0UL; j<lhs_.columns(); ++j ) {
               column( dres_  , j ) *= column( lhs_ + orhs_, j );
               column( odres_ , j ) *= column( lhs_ + orhs_, j );
               column( sres_  , j ) *= column( lhs_ + orhs_, j );
               column( osres_ , j ) *= column( lhs_ + orhs_, j );
               column( refres_, j ) *= column( reflhs_ + refrhs_, j );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t j=0UL; j<lhs_.columns(); ++j ) {
               column( dres_  , j ) *= column( olhs_ + rhs_, j );
               column( odres_ , j ) *= column( olhs_ + rhs_, j );
               column( sres_  , j ) *= column( olhs_ + rhs_, j );
               column( osres_ , j ) *= column( olhs_ + rhs_, j );
               column( refres_, j ) *= column( reflhs_ + refrhs_, j );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t j=0UL; j<lhs_.columns(); ++j ) {
               column( dres_  , j ) *= column( olhs_ + orhs_, j );
               column( odres_ , j ) *= column( olhs_ + orhs_, j );
               column( sres_  , j ) *= column( olhs_ + orhs_, j );
               column( osres_ , j ) *= column( olhs_ + orhs_, j );
               column( refres_, j ) *= column( reflhs_ + refrhs_, j );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Column-wise addition with multiplication assignment with evaluated matrices
      {
         test_  = "Column-wise addition with multiplication assignment with evaluated matrices";
         error_ = "Failed multiplication assignment operation";

         try {
            initResults();
            for( size_t j=0UL; j<lhs_.columns(); ++j ) {
               column( dres_  , j ) *= column( eval( lhs_ ) + eval( rhs_ ), j );
               column( odres_ , j ) *= column( eval( lhs_ ) + eval( rhs_ ), j );
               column( sres_  , j ) *= column( eval( lhs_ ) + eval( rhs_ ), j );
               column( osres_ , j ) *= column( eval( lhs_ ) + eval( rhs_ ), j );
               column( refres_, j ) *= column( eval( reflhs_ ) + eval( refrhs_ ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t j=0UL; j<lhs_.columns(); ++j ) {
               column( dres_  , j ) *= column( eval( lhs_ ) + eval( orhs_ ), j );
               column( odres_ , j ) *= column( eval( lhs_ ) + eval( orhs_ ), j );
               column( sres_  , j ) *= column( eval( lhs_ ) + eval( orhs_ ), j );
               column( osres_ , j ) *= column( eval( lhs_ ) + eval( orhs_ ), j );
               column( refres_, j ) *= column( eval( reflhs_ ) + eval( refrhs_ ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t j=0UL; j<lhs_.columns(); ++j ) {
               column( dres_  , j ) *= column( eval( olhs_ ) + eval( rhs_ ), j );
               column( odres_ , j ) *= column( eval( olhs_ ) + eval( rhs_ ), j );
               column( sres_  , j ) *= column( eval( olhs_ ) + eval( rhs_ ), j );
               column( osres_ , j ) *= column( eval( olhs_ ) + eval( rhs_ ), j );
               column( refres_, j ) *= column( eval( reflhs_ ) + eval( refrhs_ ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t j=0UL; j<lhs_.columns(); ++j ) {
               column( dres_  , j ) *= column( eval( olhs_ ) + eval( orhs_ ), j );
               column( odres_ , j ) *= column( eval( olhs_ ) + eval( orhs_ ), j );
               column( sres_  , j ) *= column( eval( olhs_ ) + eval( orhs_ ), j );
               column( osres_ , j ) *= column( eval( olhs_ ) + eval( orhs_ ), j );
               column( refres_, j ) *= column( eval( reflhs_ ) + eval( refrhs_ ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }
   }
#endif
}
//*************************************************************************************************




//=================================================================================================
//
//  ERROR DETECTION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Checking and comparing the computed results.
//
// \return void
// \exception std::runtime_error Incorrect dense result detected.
// \exception std::runtime_error Incorrect sparse result detected.
//
// This function is called after each test case to check and compare the computed results. The
// two template arguments \a LT and \a RT indicate the types of the left-hand side and right-hand
// side operands used for the computations.
*/
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
template< typename LT     // Type of the left-hand side operand
        , typename RT >   // Type of the right-hand side operand
void OperationTest<MT1,MT2>::checkResults()
{
   using blaze::IsRowMajorMatrix;

   if( !isEqual( dres_, refres_ ) || !isEqual( odres_, refres_ ) ) {
      std::ostringstream oss;
      oss.precision( 20 );
      oss << " Test : " << test_ << "\n"
          << " Error: Incorrect dense result detected\n"
          << " Details:\n"
          << "   Left-hand side " << ( IsRowMajorMatrix<LT>::value ? ( "row-major" ) : ( "column-major" ) ) << " dense matrix type:\n"
          << "     " << typeid( LT ).name() << "\n"
          << "   Right-hand side " << ( IsRowMajorMatrix<RT>::value ? ( "row-major" ) : ( "column-major" ) ) << " dense matrix type:\n"
          << "     " << typeid( RT ).name() << "\n"
          << "   Result:\n" << dres_ << "\n"
          << "   Result with opposite storage order:\n" << odres_ << "\n"
          << "   Expected result:\n" << refres_ << "\n";
      throw std::runtime_error( oss.str() );
   }

   if( !isEqual( sres_, refres_ ) || !isEqual( osres_, refres_ ) ) {
      std::ostringstream oss;
      oss.precision( 20 );
      oss << " Test : " << test_ << "\n"
          << " Error: Incorrect sparse result detected\n"
          << " Details:\n"
          << "   Left-hand side " << ( IsRowMajorMatrix<LT>::value ? ( "row-major" ) : ( "column-major" ) ) << " dense matrix type:\n"
          << "     " << typeid( LT ).name() << "\n"
          << "   Right-hand side " << ( IsRowMajorMatrix<RT>::value ? ( "row-major" ) : ( "column-major" ) ) << " dense matrix type:\n"
          << "     " << typeid( RT ).name() << "\n"
          << "   Result:\n" << sres_ << "\n"
          << "   Result with opposite storage order:\n" << osres_ << "\n"
          << "   Expected result:\n" << refres_ << "\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checking and comparing the computed transpose results.
//
// \return void
// \exception std::runtime_error Incorrect dense result detected.
// \exception std::runtime_error Incorrect sparse result detected.
//
// This function is called after each test case to check and compare the computed transpose
// results. The two template arguments \a LT and \a RT indicate the types of the left-hand
// side and right-hand side operands used for the computations.
*/
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
template< typename LT     // Type of the left-hand side operand
        , typename RT >   // Type of the right-hand side operand
void OperationTest<MT1,MT2>::checkTransposeResults()
{
   using blaze::IsRowMajorMatrix;

   if( !isEqual( tdres_, refres_ ) || !isEqual( todres_, refres_ ) ) {
      std::ostringstream oss;
      oss.precision( 20 );
      oss << " Test : " << test_ << "\n"
          << " Error: Incorrect dense result detected\n"
          << " Details:\n"
          << "   Left-hand side " << ( IsRowMajorMatrix<LT>::value ? ( "row-major" ) : ( "column-major" ) ) << " dense matrix type:\n"
          << "     " << typeid( LT ).name() << "\n"
          << "   Right-hand side " << ( IsRowMajorMatrix<RT>::value ? ( "row-major" ) : ( "column-major" ) ) << " dense matrix type:\n"
          << "     " << typeid( RT ).name() << "\n"
          << "   Transpose result:\n" << tdres_ << "\n"
          << "   Transpose result with opposite storage order:\n" << todres_ << "\n"
          << "   Expected result:\n" << refres_ << "\n";
      throw std::runtime_error( oss.str() );
   }

   if( !isEqual( tsres_, refres_ ) || !isEqual( tosres_, refres_ ) ) {
      std::ostringstream oss;
      oss.precision( 20 );
      oss << " Test : " << test_ << "\n"
          << " Error: Incorrect sparse result detected\n"
          << " Details:\n"
          << "   Left-hand side " << ( IsRowMajorMatrix<LT>::value ? ( "row-major" ) : ( "column-major" ) ) << " dense matrix type:\n"
          << "     " << typeid( LT ).name() << "\n"
          << "   Right-hand side " << ( IsRowMajorMatrix<RT>::value ? ( "row-major" ) : ( "column-major" ) ) << " dense matrix type:\n"
          << "     " << typeid( RT ).name() << "\n"
          << "   Transpose result:\n" << tsres_ << "\n"
          << "   Transpose result with opposite storage order:\n" << tosres_ << "\n"
          << "   Expected result:\n" << refres_ << "\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Initializing the non-transpose result matrices.
//
// \return void
//
// This function is called before each non-transpose test case to initialize the according result
// matrices to random values.
*/
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
void OperationTest<MT1,MT2>::initResults()
{
   const typename blaze::BaseElementType<RE>::Type min( randmin );
   const typename blaze::BaseElementType<RE>::Type max( randmax );

   randomize( dres_, min, max );
   odres_   = dres_;
   sres_    = dres_;
   osres_   = dres_;
   refres_  = dres_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Initializing the transpose result matrices.
//
// \return void
//
// This function is called before each transpose test case to initialize the according result
// matrices to random values.
*/
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
void OperationTest<MT1,MT2>::initTransposeResults()
{
   const typename blaze::BaseElementType<RE>::Type min( randmin );
   const typename blaze::BaseElementType<RE>::Type max( randmax );

   randomize( tdres_, min, max );
   todres_  = tdres_;
   tsres_   = tdres_;
   tosres_  = tdres_;
   refres_  = tdres_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Convert the given exception into a \a std::runtime_error exception.
//
// \param ex The \a std::exception to be extended.
// \return void
// \exception std::runtime_error The converted exception.
//
// This function converts the given exception to a \a std::runtime_error exception. Additionally,
// the function extends the given exception message by all available information for the failed
// test. The two template arguments \a LT and \a RT indicate the types of the left-hand side and
// right-hand side operands used for the computations.
*/
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
template< typename LT     // Type of the left-hand side operand
        , typename RT >   // Type of the right-hand side operand
void OperationTest<MT1,MT2>::convertException( const std::exception& ex )
{
   using blaze::IsRowMajorMatrix;

   std::ostringstream oss;
   oss << " Test : " << test_ << "\n"
       << " Error: " << error_ << "\n"
       << " Details:\n"
       << "   Left-hand side " << ( IsRowMajorMatrix<LT>::value ? ( "row-major" ) : ( "column-major" ) ) << " dense matrix type:\n"
       << "     " << typeid( LT ).name() << "\n"
       << "   Right-hand side " << ( IsRowMajorMatrix<LT>::value ? ( "row-major" ) : ( "column-major" ) ) << " dense matrix type:\n"
       << "     " << typeid( RT ).name() << "\n"
       << "   Error message: " << ex.what() << "\n";
   throw std::runtime_error( oss.str() );
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Testing the matrix addition between two specific matrix types.
//
// \param creator1 The creator for the left-hand side matrix.
// \param creator2 The creator for the right-hand side matrix.
// \return void
*/
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
void runTest( const Creator<MT1>& creator1, const Creator<MT2>& creator2 )
{
   for( size_t rep=0; rep<repetitions; ++rep ) {
      OperationTest<MT1,MT2>( creator1, creator2 );
   }
}
//*************************************************************************************************




//=================================================================================================
//
//  MACROS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Macro for the definition of a dense matrix/dense matrix addition test case.
*/
#define DEFINE_DMATDMATADD_OPERATION_TEST( MT1, MT2 ) \
   extern template class blazetest::mathtest::dmatdmatadd::OperationTest<MT1,MT2>
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Macro for the execution of a dense matrix/dense matrix addition test case.
*/
#define RUN_DMATDMATADD_OPERATION_TEST( C1, C2 ) \
   blazetest::mathtest::dmatdmatadd::runTest( C1, C2 )
/*! \endcond */
//*************************************************************************************************

} // namespace dmatdmatadd

} // namespace mathtest

} // namespace blazetest

#endif
