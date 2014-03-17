//=================================================================================================
/*!
//  \file blazetest/mathtest/svectdvecmult/OperationTest.h
//  \brief Header file for the sparse vector/dense vector outer product operation test
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

#ifndef _BLAZETEST_MATHTEST_SVECTDVECMULT_OPERATIONTEST_H_
#define _BLAZETEST_MATHTEST_SVECTDVECMULT_OPERATIONTEST_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <sstream>
#include <stdexcept>
#include <string>
#include <typeinfo>
#include <blaze/math/CompressedVector.h>
#include <blaze/math/constraints/Computation.h>
#include <blaze/math/constraints/DenseMatrix.h>
#include <blaze/math/constraints/DenseVector.h>
#include <blaze/math/constraints/SparseMatrix.h>
#include <blaze/math/constraints/SparseVector.h>
#include <blaze/math/constraints/StorageOrder.h>
#include <blaze/math/constraints/TransposeFlag.h>
#include <blaze/math/constraints/VecTVecMultExpr.h>
#include <blaze/math/DynamicVector.h>
#include <blaze/math/StaticVector.h>
#include <blaze/math/traits/MultExprTrait.h>
#include <blaze/math/traits/MultTrait.h>
#include <blaze/math/typetraits/BaseElementType.h>
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

namespace svectdvecmult {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Auxiliary class template for the sparse vector/dense vector multiplication operation test.
//
// This class template represents one particular outer product test between two vectors of a
// particular type. The two template arguments \a VT1 and \a VT2 represent the types of the
// left-hand side and right-hand side vector, respectively.
*/
template< typename VT1    // Type of the left-hand side sparse vector
        , typename VT2 >  // Type of the right-hand side dense vector
class OperationTest
{
 private:
   //**Type definitions****************************************************************************
   typedef typename VT1::TransposeType                TVT1;  //!< Transpose vector type 1
   typedef typename VT2::TransposeType                TVT2;  //!< Transpose vector type 2
   typedef typename blaze::MultTrait<VT1,TVT2>::Type  RE;    //!< Default result type
   typedef typename RE::OppositeType                  ORE;   //!< Default result type with opposite storage order
   typedef typename RE::TransposeType                 TRE;   //!< Transpose default result type
   typedef typename ORE::TransposeType                TORE;  //!< Transpose default result type with opposite storage order
   //**********************************************************************************************

   //**Type definitions****************************************************************************
   typedef typename VT1::ElementType          ET1;     //!< Element type 1
   typedef typename VT2::ElementType          ET2;     //!< Element type 2
   typedef typename RE::ElementType           RET;     //!< Resulting element type
   typedef blaze::DynamicVector<ET1,false>    RT1;     //!< Reference type 1
   typedef blaze::DynamicVector<ET2,true>     RT2;     //!< Reference type 2
   typedef blaze::DynamicMatrix<RET,true>     DRRE;    //!< Dense reference result type
   typedef blaze::CompressedMatrix<RET,true>  SRRE;    //!< Sparse reference result type
   typedef typename DRRE::OppositeType        ODRRE;   //!< Dense reference result type with opposite storage order
   typedef typename SRRE::OppositeType        OSRRE;   //!< Sparse reference result type with opposite storage order
   typedef typename DRRE::TransposeType       TDRRE;   //!< Transpose dense reference result type
   typedef typename SRRE::TransposeType       TSRRE;   //!< Transpose sparse reference result type
   typedef typename ODRRE::TransposeType      TODRRE;  //!< Transpose dense reference result type with opposite storage order
   typedef typename OSRRE::TransposeType      TOSRRE;  //!< Transpose sparse reference result type with opposite storage order
   typedef DRRE                               DRE;     //!< Dense result type
   typedef RE                                 SRE;     //!< Sparse result type
   typedef ODRRE                              ODRE;    //!< Dense result type with opposite storage order
   typedef ORE                                OSRE;    //!< Sparse result type with opposite storage order
   typedef TDRRE                              TDRE;    //!< Transpose dense result type
   typedef TRE                                TSRE;    //!< Transpose sparse result type
   typedef TODRRE                             TODRE;   //!< Transpose dense result type with opposite storage order
   typedef TORE                               TOSRE;   //!< Transpose sparse result type with opposite storage order

   //! Type of the outer product expression
   typedef typename blaze::MultExprTrait<VT1,TVT2>::Type  VecTVecMultExprType;
   //**********************************************************************************************

 public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit OperationTest( const Creator<VT1>& creator1, const Creator<VT2>& creator2 );
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
   void checkResults();
   void checkTransposeResults();
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   void initResults();
   void initTransposeResults();
   void convertException( const std::exception& ex );
   //@}
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   VT1   lhs_;     //!< The left-hand side sparse vector.
   TVT2  rhs_;     //!< The right-hand side dense vector.
   DRE   dres_;    //!< The dense result matrix.
   SRE   sres_;    //!< The sparse result matrix.
   ODRE  odres_;   //!< The dense result matrix with opposite storage order.
   OSRE  osres_;   //!< The sparse result matrix with opposite storage order.
   TDRE  tdres_;   //!< The transpose dense result matrix.
   TSRE  tsres_;   //!< The transpose sparse result matrix.
   TODRE todres_;  //!< The transpose dense result matrix with opposite storage order.
   TOSRE tosres_;  //!< The transpose sparse result matrix with opposite storage order.
   RT1   reflhs_;  //!< The reference left-hand side vector.
   RT2   refrhs_;  //!< The reference right-hand side vector.
   DRRE  refres_;  //!< The reference result.

   std::string test_;   //!< Label of the currently performed test.
   std::string error_;  //!< Description of the current error type.
   //@}
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( VT1    );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( VT2    );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( TVT1   );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( TVT2   );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( RT1    );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( RT2    );
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
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE      ( VT1    );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE      ( VT2    );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE         ( TVT1   );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE         ( TVT2   );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE      ( RT1    );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE         ( RT2    );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( DRE    );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( SRE    );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE   ( ODRE   );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE   ( ODRE   );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE   ( TDRE   );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE   ( TDRE   );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( TODRE  );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( TODRE  );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( DRRE   );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( SRRE   );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE   ( ODRRE  );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE   ( ODRRE  );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE   ( TDRRE  );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE   ( TDRRE  );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( TODRRE );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( TODRRE );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ET1, typename TVT1::ElementType   );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ET2, typename TVT2::ElementType   );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( VT1, typename TVT1::TransposeType );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( VT2, typename TVT2::TransposeType );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( RE , typename ORE::OppositeType   );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( RE , typename TRE::TransposeType  );
   BLAZE_CONSTRAINT_MUST_BE_VECTVECMULTEXPR_TYPE( VecTVecMultExprType );
   BLAZE_CONSTRAINT_MUST_BE_COMPUTATION_TYPE    ( VecTVecMultExprType );
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
/*!\brief Constructor for the sparse vector/dense vector outer product operation test.
//
// \param creator1 The creator for the left-hand side sparse vector of the vector outer product.
// \param creator2 The creator for the right-hand side dense vector of the vector outer product.
// \exception std::runtime_error Operation error detected.
*/
template< typename VT1    // Type of the left-hand side sparse vector
        , typename VT2 >  // Type of the right-hand side dense vector
OperationTest<VT1,VT2>::OperationTest( const Creator<VT1>& creator1, const Creator<VT2>& creator2 )
   : lhs_( creator1() )           // The left-hand side sparse vector
   , rhs_( trans( creator2() ) )  // The right-hand side dense vector
   , dres_()                      // The dense result matrix
   , sres_()                      // The sparse result matrix
   , odres_()                     // The dense result matrix with opposite storage order
   , osres_()                     // The sparse result matrix with opposite storage order
   , tdres_()                     // The transpose dense result matrix
   , tsres_()                     // The transpose sparse result matrix
   , todres_()                    // The transpose dense result matrix with opposite storage order
   , tosres_()                    // The transpose sparse result matrix with opposite storage order
   , reflhs_( lhs_ )              // The reference left-hand side vector
   , refrhs_( rhs_ )              // The reference right-hand side vector
   , refres_()                    // The reference result
   , test_()                      // Label of the currently performed test
   , error_()                     // Description of the current error type
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
/*!\brief Tests on the initial status of the vectors.
//
// \return void
// \exception std::runtime_error Initialization error detected.
//
// This function runs tests on the initial status of the vectors. In case any initialization
// error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename VT1    // Type of the left-hand side sparse vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::testInitialStatus()
{
   // Checking the size of the left-hand side operand
   if( lhs_.size() != reflhs_.size() ) {
      std::ostringstream oss;
      oss << " Test: Initial size comparison of left-hand side sparse operand\n"
          << " Error: Invalid vector size\n"
          << " Details:\n"
          << "   Sparse vector type:\n"
          << "     " << typeid( VT1 ).name() << "\n"
          << "   Detected size = " << lhs_.size() << "\n"
          << "   Expected size = " << reflhs_.size() << "\n";
      throw std::runtime_error( oss.str() );
   }

   // Checking the size of the right-hand side operand
   if( rhs_.size() != refrhs_.size() ) {
      std::ostringstream oss;
      oss << " Test: Initial size comparison of right-hand side dense operand\n"
          << " Error: Invalid vector size\n"
          << " Details:\n"
          << "   Transpose dense vector type:\n"
          << "     " << typeid( TVT2 ).name() << "\n"
          << "   Detected size = " << rhs_.size() << "\n"
          << "   Expected size = " << refrhs_.size() << "\n";
      throw std::runtime_error( oss.str() );
   }

   // Checking the initialization of the left-hand side operand
   if( !isEqual( lhs_, reflhs_ ) ) {
      std::ostringstream oss;
      oss << " Test: Initial test of initialization of left-hand side sparse operand\n"
          << " Error: Invalid vector initialization\n"
          << " Details:\n"
          << "   Sparse vector type:\n"
          << "     " << typeid( VT1 ).name() << "\n"
          << "   Current initialization:\n" << lhs_ << "\n"
          << "   Expected initialization:\n" << reflhs_ << "\n";
      throw std::runtime_error( oss.str() );
   }

   // Checking the initialization of the right-hand side operand
   if( !isEqual( rhs_, refrhs_ ) ) {
      std::ostringstream oss;
      oss << " Test: Initial test of initialization of right-hand side dense operand\n"
          << " Error: Invalid vector initialization\n"
          << " Details:\n"
          << "   Transpose dense vector type:\n"
          << "     " << typeid( TVT2 ).name() << "\n"
          << "   Current initialization:\n" << rhs_ << "\n"
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
template< typename VT1    // Type of the left-hand side sparse vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::testAssignment()
{
   try {
      lhs_ = reflhs_;
      rhs_ = refrhs_;
   }
   catch( std::exception& ex ) {
      std::ostringstream oss;
      oss << " Test: Assignment with the given vectors\n"
          << " Error: Failed assignment\n"
          << " Details:\n"
          << "   Left-hand side sparse vector type:\n"
          << "     " << typeid( VT1 ).name() << "\n"
          << "   Right-hand side transpose dense vector type:\n"
          << "     " << typeid( TVT2 ).name() << "\n"
          << "   Error message: " << ex.what() << "\n";
      throw std::runtime_error( oss.str() );
   }

   if( !isEqual( lhs_, reflhs_ ) ) {
      std::ostringstream oss;
      oss << " Test: Checking the assignment result of left-hand side sparse operand\n"
          << " Error: Invalid vector initialization\n"
          << " Details:\n"
          << "   Sparse vector type:\n"
          << "     " << typeid( VT1 ).name() << "\n"
          << "   Current initialization:\n" << lhs_ << "\n"
          << "   Expected initialization:\n" << reflhs_ << "\n";
      throw std::runtime_error( oss.str() );
   }

   if( !isEqual( rhs_, refrhs_ ) ) {
      std::ostringstream oss;
      oss << " Test: Checking the assignment result of right-hand side dense operand\n"
          << " Error: Invalid vector initialization\n"
          << " Details:\n"
          << "   Transpose dense vector type:\n"
          << "     " << typeid( TVT2 ).name() << "\n"
          << "   Current initialization:\n" << rhs_ << "\n"
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
template< typename VT1    // Type of the left-hand side sparse vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::testElementAccess()
{
   using blaze::equal;

   if( lhs_.size() > 0UL && rhs_.size() > 0UL )
   {
      if( !equal( ( lhs_ * rhs_ )(0UL,0UL), ( reflhs_ * refrhs_ )(0UL,0UL) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of outer product expression\n"
             << " Error: Unequal resulting elements at element (0,0) detected\n"
             << " Details:\n"
             << "   Left-hand side sparse vector type:\n"
             << "     " << typeid( VT1 ).name() << "\n"
             << "   Right-hand side transpose dense vector type:\n"
             << "     " << typeid( TVT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( ( lhs_ * eval( rhs_ ) )(0UL,0UL), ( reflhs_ * eval( refrhs_ ) )(0UL,0UL) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of right evaluated addition expression\n"
             << " Error: Unequal resulting elements at element (0,0) detected\n"
             << " Details:\n"
             << "   Left-hand side sparse vector type:\n"
             << "     " << typeid( VT1 ).name() << "\n"
             << "   Right-hand side transpose dense vector type:\n"
             << "     " << typeid( TVT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( ( eval( lhs_ ) * rhs_ )(0UL,0UL), ( eval( reflhs_ ) * refrhs_ )(0UL,0UL) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of left evaluated addition expression\n"
             << " Error: Unequal resulting elements at element (0,0) detected\n"
             << " Details:\n"
             << "   Left-hand side sparse vector type:\n"
             << "     " << typeid( VT1 ).name() << "\n"
             << "   Right-hand side transpose dense vector type:\n"
             << "     " << typeid( TVT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( ( eval( lhs_ ) * eval( rhs_ ) )(0UL,0UL), ( eval( reflhs_ ) * eval( refrhs_ ) )(0UL,0UL) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of fully evaluated addition expression\n"
             << " Error: Unequal resulting elements at element (0,0) detected\n"
             << " Details:\n"
             << "   Left-hand side sparse vector type:\n"
             << "     " << typeid( VT1 ).name() << "\n"
             << "   Right-hand side transpose dense vector type:\n"
             << "     " << typeid( TVT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the plain sparse vector/dense vector outer product.
//
// \return void
// \exception std::runtime_error Outer product error detected.
//
// This function tests the plain outer product with plain assignment, addition assignment,
// and subtraction assignment. In case any error resulting from the outer product or the
// subsequent assignment is detected, a \a std::runtime_error exception is thrown.
*/
template< typename VT1    // Type of the left-hand side sparse vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::testBasicOperation()
{
#if BLAZETEST_MATHTEST_TEST_BASIC_OPERATION
   if( BLAZETEST_MATHTEST_TEST_BASIC_OPERATION > 1 )
   {
      //=====================================================================================
      // Outer product
      //=====================================================================================

      // Outer product with the given vectors
      {
         test_  = "Outer product with the given vectors";
         error_ = "Failed outer product operation";

         try {
            initResults();
            dres_   = lhs_ * rhs_;
            odres_  = lhs_ * rhs_;
            sres_   = lhs_ * rhs_;
            osres_  = lhs_ * rhs_;
            refres_ = reflhs_ * refrhs_;
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Outer product with evaluated vectors
      {
         test_  = "Outer product with evaluated vectors";
         error_ = "Failed outer product operation";

         try {
            initResults();
            dres_   = eval( lhs_ ) * eval( rhs_ );
            odres_  = eval( lhs_ ) * eval( rhs_ );
            sres_   = eval( lhs_ ) * eval( rhs_ );
            osres_  = eval( lhs_ ) * eval( rhs_ );
            refres_ = eval( reflhs_ ) * eval( refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }


      //=====================================================================================
      // Outer product with addition assignment
      //=====================================================================================

      // Outer product with addition assignment with the given vectors
      {
         test_  = "Outer product with addition assignment with the given vectors";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += lhs_ * rhs_;
            odres_  += lhs_ * rhs_;
            sres_   += lhs_ * rhs_;
            osres_  += lhs_ * rhs_;
            refres_ += reflhs_ * refrhs_;
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Outer product with addition assignment with evaluated vectors
      {
         test_  = "Outer product with addition assignment with evaluated vectors";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += eval( lhs_ ) * eval( rhs_ );
            odres_  += eval( lhs_ ) * eval( rhs_ );
            sres_   += eval( lhs_ ) * eval( rhs_ );
            osres_  += eval( lhs_ ) * eval( rhs_ );
            refres_ += eval( reflhs_ ) * eval( refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }


      //=====================================================================================
      // Outer product with subtraction assignment
      //=====================================================================================

      // Outer product with subtraction assignment with the given vectors
      {
         test_  = "Outer product with subtraction assignment with the given vectors";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= lhs_ * rhs_;
            odres_  -= lhs_ * rhs_;
            sres_   -= lhs_ * rhs_;
            osres_  -= lhs_ * rhs_;
            refres_ -= reflhs_ * refrhs_;
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Outer product with subtraction assignment with evaluated vectors
      {
         test_  = "Outer product with subtraction assignment with evaluated vectors";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= eval( lhs_ ) * eval( rhs_ );
            odres_  -= eval( lhs_ ) * eval( rhs_ );
            sres_   -= eval( lhs_ ) * eval( rhs_ );
            osres_  -= eval( lhs_ ) * eval( rhs_ );
            refres_ -= eval( reflhs_ ) * eval( refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the negated sparse vector/dense vector outer product.
//
// \return void
// \exception std::runtime_error Outer product error detected.
//
// This function tests the negated outer product with plain assignment, addition assignment,
// and subtraction assignment. In case any error resulting from the outer product or the
// subsequent assignment is detected, a \a std::runtime_error exception is thrown.
*/
template< typename VT1    // Type of the left-hand side sparse vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::testNegatedOperation()
{
#if BLAZETEST_MATHTEST_TEST_NEGATED_OPERATION
   if( BLAZETEST_MATHTEST_TEST_NEGATED_OPERATION > 1 )
   {
      //=====================================================================================
      // Negated outer product
      //=====================================================================================

      // Negated outer product with the given vectors
      {
         test_  = "Negated outer product with the given vectors";
         error_ = "Failed outer product operation";

         try {
            initResults();
            dres_   = -( lhs_ * rhs_ );
            odres_  = -( lhs_ * rhs_ );
            sres_   = -( lhs_ * rhs_ );
            osres_  = -( lhs_ * rhs_ );
            refres_ = -( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Negated outer product with evaluated vectors
      {
         test_  = "Negated outer product with evaluated vectors";
         error_ = "Failed outer product operation";

         try {
            initResults();
            dres_   = -( eval( lhs_ ) * eval( rhs_ ) );
            odres_  = -( eval( lhs_ ) * eval( rhs_ ) );
            sres_   = -( eval( lhs_ ) * eval( rhs_ ) );
            osres_  = -( eval( lhs_ ) * eval( rhs_ ) );
            refres_ = -( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }


      //=====================================================================================
      // Negated outer product with addition assignment
      //=====================================================================================

      // Negated outer product with addition assignment with the given vectors
      {
         test_  = "Negated outer product with addition assignment with the given vectors";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += -( lhs_ * rhs_ );
            odres_  += -( lhs_ * rhs_ );
            sres_   += -( lhs_ * rhs_ );
            osres_  += -( lhs_ * rhs_ );
            refres_ += -( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Negated outer product with addition assignment with evaluated vectors
      {
         test_  = "Negated outer product with addition assignment with evaluated vectors";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += -( eval( lhs_ ) * eval( rhs_ ) );
            odres_  += -( eval( lhs_ ) * eval( rhs_ ) );
            sres_   += -( eval( lhs_ ) * eval( rhs_ ) );
            osres_  += -( eval( lhs_ ) * eval( rhs_ ) );
            refres_ += -( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }


      //=====================================================================================
      // Negated outer product with subtraction assignment
      //=====================================================================================

      // Negated outer product with subtraction assignment with the given vectors
      {
         test_  = "Negated outer product with subtraction assignment with the given vectors";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= -( lhs_ * rhs_ );
            odres_  -= -( lhs_ * rhs_ );
            sres_   -= -( lhs_ * rhs_ );
            osres_  -= -( lhs_ * rhs_ );
            refres_ -= -( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Negated outer product with subtraction assignment with evaluated vectors
      {
         test_  = "Negated outer product with subtraction assignment with evaluated vectors";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= -( eval( lhs_ ) * eval( rhs_ ) );
            odres_  -= -( eval( lhs_ ) * eval( rhs_ ) );
            sres_   -= -( eval( lhs_ ) * eval( rhs_ ) );
            osres_  -= -( eval( lhs_ ) * eval( rhs_ ) );
            refres_ -= -( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the scaled sparse vector/dense vector outer product.
//
// \param scalar The scalar value.
// \return void
// \exception std::runtime_error Outer product error detected.
//
// This function tests the scaled outer product with plain assignment, addition assignment,
// and subtraction assignment. In case any error resulting from the outer product or the
// subsequent assignment is detected, a \a std::runtime_error exception is thrown.
*/
template< typename VT1    // Type of the left-hand side sparse vector
        , typename VT2 >  // Type of the right-hand side dense vector
template< typename T >    // Type of the scalar
void OperationTest<VT1,VT2>::testScaledOperation( T scalar )
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
            dres_   = lhs_ * rhs_;
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

         checkResults();
      }


      //=====================================================================================
      // Self-scaling (M=M*s)
      //=====================================================================================

      // Self-scaling (M=M*s)
      {
         test_ = "Self-scaling (M=M*s)";

         try {
            dres_   = lhs_ * rhs_;
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

         checkResults();
      }


      //=====================================================================================
      // Self-scaling (M=s*M)
      //=====================================================================================

      // Self-scaling (M=s*M)
      {
         test_ = "Self-scaling (M=s*M)";

         try {
            dres_   = lhs_ * rhs_;
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

         checkResults();
      }


      //=====================================================================================
      // Self-scaling (M/=s)
      //=====================================================================================

      // Self-scaling (M/=s)
      {
         test_ = "Self-scaling (M/=s)";

         try {
            dres_   = lhs_ * rhs_;
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

         checkResults();
      }


      //=====================================================================================
      // Self-scaling (M=M/s)
      //=====================================================================================

      // Self-scaling (M=M/s)
      {
         test_ = "Self-scaling (M=M/s)";

         try {
            dres_   = lhs_ * rhs_;
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

         checkResults();
      }


      //=====================================================================================
      // Scaled outer product (s*OP)
      //=====================================================================================

      // Scaled outer product with the given vectors
      {
         test_  = "Scaled outer product with the given vectors (s*OP)";
         error_ = "Failed outer product operation";

         try {
            initResults();
            dres_   = scalar * ( lhs_ * rhs_ );
            odres_  = scalar * ( lhs_ * rhs_ );
            sres_   = scalar * ( lhs_ * rhs_ );
            osres_  = scalar * ( lhs_ * rhs_ );
            refres_ = scalar * ( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Scaled outer product with evaluated vectors
      {
         test_  = "Scaled outer product with evaluated vectors (s*OP)";
         error_ = "Failed outer product operation";

         try {
            initResults();
            dres_   = scalar * ( eval( lhs_ ) * eval( rhs_ ) );
            odres_  = scalar * ( eval( lhs_ ) * eval( rhs_ ) );
            sres_   = scalar * ( eval( lhs_ ) * eval( rhs_ ) );
            osres_  = scalar * ( eval( lhs_ ) * eval( rhs_ ) );
            refres_ = scalar * ( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }


      //=====================================================================================
      // Scaled outer product (OP*s)
      //=====================================================================================

      // Scaled outer product with the given vectors
      {
         test_  = "Scaled outer product with the given vectors (OP*s)";
         error_ = "Failed outer product operation";

         try {
            initResults();
            dres_   = ( lhs_ * rhs_ ) * scalar;
            odres_  = ( lhs_ * rhs_ ) * scalar;
            sres_   = ( lhs_ * rhs_ ) * scalar;
            osres_  = ( lhs_ * rhs_ ) * scalar;
            refres_ = ( reflhs_ * refrhs_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Scaled outer product with evaluated vectors
      {
         test_  = "Scaled outer product with evaluated vectors (OP*s)";
         error_ = "Failed outer product operation";

         try {
            initResults();
            dres_   = ( eval( lhs_ ) * eval( rhs_ ) ) * scalar;
            odres_  = ( eval( lhs_ ) * eval( rhs_ ) ) * scalar;
            sres_   = ( eval( lhs_ ) * eval( rhs_ ) ) * scalar;
            osres_  = ( eval( lhs_ ) * eval( rhs_ ) ) * scalar;
            refres_ = ( eval( reflhs_ ) * eval( refrhs_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }


      //=====================================================================================
      // Scaled outer product (OP/s)
      //=====================================================================================

      // Scaled outer product with the given vectors
      {
         test_  = "Scaled outer product with the given vectors (OP/s)";
         error_ = "Failed outer product operation";

         try {
            initResults();
            dres_   = ( lhs_ * rhs_ ) / scalar;
            odres_  = ( lhs_ * rhs_ ) / scalar;
            sres_   = ( lhs_ * rhs_ ) / scalar;
            osres_  = ( lhs_ * rhs_ ) / scalar;
            refres_ = ( reflhs_ * refrhs_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Scaled outer product with evaluated vectors
      {
         test_  = "Scaled outer product with evaluated vectors (OP/s)";
         error_ = "Failed outer product operation";

         try {
            initResults();
            dres_   = ( eval( lhs_ ) * eval( rhs_ ) ) / scalar;
            odres_  = ( eval( lhs_ ) * eval( rhs_ ) ) / scalar;
            sres_   = ( eval( lhs_ ) * eval( rhs_ ) ) / scalar;
            osres_  = ( eval( lhs_ ) * eval( rhs_ ) ) / scalar;
            refres_ = ( eval( reflhs_ ) * eval( refrhs_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }


      //=====================================================================================
      // Scaled outer product with addition assignment (s*OP)
      //=====================================================================================

      // Scaled outer product with addition assignment with the given vectors
      {
         test_  = "Scaled outer product with addition assignment with the given vectors (s*OP)";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += scalar * ( lhs_ * rhs_ );
            odres_  += scalar * ( lhs_ * rhs_ );
            sres_   += scalar * ( lhs_ * rhs_ );
            osres_  += scalar * ( lhs_ * rhs_ );
            refres_ += scalar * ( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Scaled outer product with addition assignment with evaluated vectors
      {
         test_  = "Scaled outer product with addition assignment with evaluated vectors (s*OP)";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += scalar * ( eval( lhs_ ) * eval( rhs_ ) );
            odres_  += scalar * ( eval( lhs_ ) * eval( rhs_ ) );
            sres_   += scalar * ( eval( lhs_ ) * eval( rhs_ ) );
            osres_  += scalar * ( eval( lhs_ ) * eval( rhs_ ) );
            refres_ += scalar * ( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }


      //=====================================================================================
      // Scaled outer product with addition assignment (OP*s)
      //=====================================================================================

      // Scaled outer product with addition assignment with the given vectors
      {
         test_  = "Scaled outer product with addition assignment with the given vectors (OP*s)";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += ( lhs_ * rhs_ ) * scalar;
            odres_  += ( lhs_ * rhs_ ) * scalar;
            sres_   += ( lhs_ * rhs_ ) * scalar;
            osres_  += ( lhs_ * rhs_ ) * scalar;
            refres_ += ( reflhs_ * refrhs_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Scaled outer product with addition assignment with evaluated vectors
      {
         test_  = "Scaled outer product with addition assignment with evaluated vectors (OP*s)";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += ( eval( lhs_ ) * eval( rhs_ ) ) * scalar;
            odres_  += ( eval( lhs_ ) * eval( rhs_ ) ) * scalar;
            sres_   += ( eval( lhs_ ) * eval( rhs_ ) ) * scalar;
            osres_  += ( eval( lhs_ ) * eval( rhs_ ) ) * scalar;
            refres_ += ( eval( reflhs_ ) * eval( refrhs_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }


      //=====================================================================================
      // Scaled outer product with addition assignment (OP/s)
      //=====================================================================================

      // Scaled outer product with addition assignment with the given vectors
      {
         test_  = "Scaled outer product with addition assignment with the given vectors (OP/s)";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += ( lhs_ * rhs_ ) / scalar;
            odres_  += ( lhs_ * rhs_ ) / scalar;
            sres_   += ( lhs_ * rhs_ ) / scalar;
            osres_  += ( lhs_ * rhs_ ) / scalar;
            refres_ += ( reflhs_ * refrhs_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Scaled outer product with addition assignment with evaluated vectors
      {
         test_  = "Scaled outer product with addition assignment with evaluated vectors (OP/s)";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += ( eval( lhs_ ) * eval( rhs_ ) ) / scalar;
            odres_  += ( eval( lhs_ ) * eval( rhs_ ) ) / scalar;
            sres_   += ( eval( lhs_ ) * eval( rhs_ ) ) / scalar;
            osres_  += ( eval( lhs_ ) * eval( rhs_ ) ) / scalar;
            refres_ += ( eval( reflhs_ ) * eval( refrhs_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }


      //=====================================================================================
      // Scaled outer product with subtraction assignment (s*OP)
      //=====================================================================================

      // Scaled outer product with subtraction assignment with the given vectors
      {
         test_  = "Scaled outer product with subtraction assignment with the given vectors (s*OP)";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= scalar * ( lhs_ * rhs_ );
            odres_  -= scalar * ( lhs_ * rhs_ );
            sres_   -= scalar * ( lhs_ * rhs_ );
            osres_  -= scalar * ( lhs_ * rhs_ );
            refres_ -= scalar * ( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Scaled outer product with subtraction assignment with evaluated vectors
      {
         test_  = "Scaled outer product with subtraction assignment with evaluated vectors (s*OP)";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= scalar * ( eval( lhs_ ) * eval( rhs_ ) );
            odres_  -= scalar * ( eval( lhs_ ) * eval( rhs_ ) );
            sres_   -= scalar * ( eval( lhs_ ) * eval( rhs_ ) );
            osres_  -= scalar * ( eval( lhs_ ) * eval( rhs_ ) );
            refres_ -= scalar * ( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }


      //=====================================================================================
      // Scaled outer product with subtraction assignment (OP*s)
      //=====================================================================================

      // Scaled outer product with subtraction assignment with the given vectors
      {
         test_  = "Scaled outer product with subtraction assignment with the given vectors (OP*s)";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= ( lhs_ * rhs_ ) * scalar;
            odres_  -= ( lhs_ * rhs_ ) * scalar;
            sres_   -= ( lhs_ * rhs_ ) * scalar;
            osres_  -= ( lhs_ * rhs_ ) * scalar;
            refres_ -= ( reflhs_ * refrhs_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Scaled outer product with subtraction assignment with evaluated vectors
      {
         test_  = "Scaled outer product with subtraction assignment with evaluated vectors (OP*s)";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= ( eval( lhs_ ) * eval( rhs_ ) ) * scalar;
            odres_  -= ( eval( lhs_ ) * eval( rhs_ ) ) * scalar;
            sres_   -= ( eval( lhs_ ) * eval( rhs_ ) ) * scalar;
            osres_  -= ( eval( lhs_ ) * eval( rhs_ ) ) * scalar;
            refres_ -= ( eval( reflhs_ ) * eval( refrhs_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }


      //=====================================================================================
      // Scaled outer product with subtraction assignment (OP/s)
      //=====================================================================================

      // Scaled outer product with subtraction assignment with the given vectors
      {
         test_  = "Scaled outer product with subtraction assignment with the given vectors (OP/s)";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= ( lhs_ * rhs_ ) / scalar;
            odres_  -= ( lhs_ * rhs_ ) / scalar;
            sres_   -= ( lhs_ * rhs_ ) / scalar;
            osres_  -= ( lhs_ * rhs_ ) / scalar;
            refres_ -= ( reflhs_ * refrhs_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Scaled outer product with subtraction assignment with evaluated vectors
      {
         test_  = "Scaled outer product with subtraction assignment with evaluated vectors (OP/s)";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= ( eval( lhs_ ) * eval( rhs_ ) ) / scalar;
            odres_  -= ( eval( lhs_ ) * eval( rhs_ ) ) / scalar;
            sres_   -= ( eval( lhs_ ) * eval( rhs_ ) ) / scalar;
            osres_  -= ( eval( lhs_ ) * eval( rhs_ ) ) / scalar;
            refres_ -= ( eval( reflhs_ ) * eval( refrhs_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the transpose sparse vector/dense vector outer product.
//
// \return void
// \exception std::runtime_error Outer product error detected.
//
// This function tests the transpose outer product with plain assignment, addition assignment,
// and subtraction assignment. In case any error resulting from the outer product or the
// subsequent assignment is detected, a \a std::runtime_error exception is thrown.
*/
template< typename VT1    // Type of the left-hand side sparse vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::testTransposeOperation()
{
#if BLAZETEST_MATHTEST_TEST_TRANSPOSE_OPERATION
   if( BLAZETEST_MATHTEST_TEST_TRANSPOSE_OPERATION > 1 )
   {
      //=====================================================================================
      // Transpose outer product
      //=====================================================================================

      // Transpose outer product with the given vectors
      {
         test_  = "Transpose outer product with the given vectors";
         error_ = "Failed outer product operation";

         try {
            initTransposeResults();
            tdres_  = trans( lhs_ * rhs_ );
            todres_ = trans( lhs_ * rhs_ );
            tsres_  = trans( lhs_ * rhs_ );
            tosres_ = trans( lhs_ * rhs_ );
            refres_ = trans( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkTransposeResults();
      }

      // Transpose outer product with evaluated vectors
      {
         test_  = "Transpose outer product with evaluated vectors";
         error_ = "Failed outer product operation";

         try {
            initTransposeResults();
            tdres_  = trans( eval( lhs_ ) * eval( rhs_ ) );
            todres_ = trans( eval( lhs_ ) * eval( rhs_ ) );
            tsres_  = trans( eval( lhs_ ) * eval( rhs_ ) );
            tosres_ = trans( eval( lhs_ ) * eval( rhs_ ) );
            refres_ = trans( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkTransposeResults();
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the abs sparse vector/dense vector outer product.
//
// \return void
// \exception std::runtime_error Outer product error detected.
//
// This function tests the abs outer product with plain assignment, addition assignment,
// and subtraction assignment. In case any error resulting from the outer product or the
// subsequent assignment is detected, a \a std::runtime_error exception is thrown.
*/
template< typename VT1    // Type of the left-hand side sparse vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::testAbsOperation()
{
#if BLAZETEST_MATHTEST_TEST_ABS_OPERATION
   if( BLAZETEST_MATHTEST_TEST_ABS_OPERATION > 1 )
   {
      //=====================================================================================
      // Abs outer product
      //=====================================================================================

      // Abs outer product with the given vectors
      {
         test_  = "Abs outer product with the given vectors";
         error_ = "Failed outer product operation";

         try {
            initResults();
            dres_   = abs( lhs_ * rhs_ );
            odres_  = abs( lhs_ * rhs_ );
            sres_   = abs( lhs_ * rhs_ );
            osres_  = abs( lhs_ * rhs_ );
            refres_ = abs( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Abs outer product with evaluated vectors
      {
         test_  = "Abs outer product with evaluated vectors";
         error_ = "Failed outer product operation";

         try {
            initResults();
            dres_   = abs( eval( lhs_ ) * eval( rhs_ ) );
            odres_  = abs( eval( lhs_ ) * eval( rhs_ ) );
            sres_   = abs( eval( lhs_ ) * eval( rhs_ ) );
            osres_  = abs( eval( lhs_ ) * eval( rhs_ ) );
            refres_ = abs( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }


      //=====================================================================================
      // Abs outer product with addition assignment
      //=====================================================================================

      // Abs outer product with addition assignment with the given vectors
      {
         test_  = "Abs outer product with addition assignment with the given vectors";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += abs( lhs_ * rhs_ );
            odres_  += abs( lhs_ * rhs_ );
            sres_   += abs( lhs_ * rhs_ );
            osres_  += abs( lhs_ * rhs_ );
            refres_ += abs( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Abs outer product with addition assignment with evaluated vectors
      {
         test_  = "Abs outer product with addition assignment with evaluated vectors";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += abs( eval( lhs_ ) * eval( rhs_ ) );
            odres_  += abs( eval( lhs_ ) * eval( rhs_ ) );
            sres_   += abs( eval( lhs_ ) * eval( rhs_ ) );
            osres_  += abs( eval( lhs_ ) * eval( rhs_ ) );
            refres_ += abs( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }


      //=====================================================================================
      // Abs outer product with subtraction assignment
      //=====================================================================================

      // Abs outer product with subtraction assignment with the given vectors
      {
         test_  = "Abs outer product with subtraction assignment with the given vectors";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= abs( lhs_ * rhs_ );
            odres_  -= abs( lhs_ * rhs_ );
            sres_   -= abs( lhs_ * rhs_ );
            osres_  -= abs( lhs_ * rhs_ );
            refres_ -= abs( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Abs outer product with subtraction assignment with evaluated vectors
      {
         test_  = "Abs outer product with subtraction assignment with evaluated vectors";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= abs( eval( lhs_ ) * eval( rhs_ ) );
            odres_  -= abs( eval( lhs_ ) * eval( rhs_ ) );
            sres_   -= abs( eval( lhs_ ) * eval( rhs_ ) );
            osres_  -= abs( eval( lhs_ ) * eval( rhs_ ) );
            refres_ -= abs( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the submatrix-wise sparse vector/dense vector outer product.
//
// \return void
// \exception std::runtime_error Outer product error detected.
//
// This function tests the submatrix-wise outer product with plain assignment, addition
// assignment, and subtraction assignment. In case any error resulting from the outer product
// or the subsequent assignment is detected, a \a std::runtime_error exception is thrown.
*/
template< typename VT1    // Type of the left-hand side sparse vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::testSubmatrixOperation()
{
#if BLAZETEST_MATHTEST_TEST_SUBMATRIX_OPERATION
   if( BLAZETEST_MATHTEST_TEST_SUBMATRIX_OPERATION > 1 )
   {
      if( lhs_.size() == 0UL || rhs_.size() == 0UL )
         return;


      //=====================================================================================
      // Submatrix-wise outer product
      //=====================================================================================

      // Submatrix-wise outer product with the given vectors
      {
         test_  = "Submatrix-wise outer product with the given vectors";
         error_ = "Failed outer product operation";

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.size(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.size() - row );
               for( size_t column=0UL, n=0UL; column<rhs_.size(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, rhs_.size() - column );
                  submatrix( dres_  , row, column, m, n ) = submatrix( lhs_ * rhs_      , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) = submatrix( lhs_ * rhs_      , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) = submatrix( lhs_ * rhs_      , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) = submatrix( lhs_ * rhs_      , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) = submatrix( reflhs_ * refrhs_, row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Submatrix-wise outer product with evaluated vectors
      {
         test_  = "Submatrix-wise outer product with evaluated vectors";
         error_ = "Failed outer product operation";

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.size(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.size() - row );
               for( size_t column=0UL, n=0UL; column<rhs_.size(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, rhs_.size() - column );
                  submatrix( dres_  , row, column, m, n ) = submatrix( eval( lhs_ ) * eval( rhs_ )      , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) = submatrix( eval( lhs_ ) * eval( rhs_ )      , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) = submatrix( eval( lhs_ ) * eval( rhs_ )      , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) = submatrix( eval( lhs_ ) * eval( rhs_ )      , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) = submatrix( eval( reflhs_ ) * eval( refrhs_ ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }


      //=====================================================================================
      // Submatrix-wise outer product with addition assignment
      //=====================================================================================

      // Submatrix-wise outer product with addition assignment with the given vectors
      {
         test_  = "Submatrix-wise outer product with addition assignment with the given vectors";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.size(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.size() - row );
               for( size_t column=0UL, n=0UL; column<rhs_.size(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, rhs_.size() - column );
                  submatrix( dres_  , row, column, m, n ) += submatrix( lhs_ * rhs_      , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) += submatrix( lhs_ * rhs_      , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) += submatrix( lhs_ * rhs_      , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) += submatrix( lhs_ * rhs_      , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) += submatrix( reflhs_ * refrhs_, row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Submatrix-wise outer product with addition assignment with evaluated vectors
      {
         test_  = "Submatrix-wise outer product with addition assignment with evaluated vectors";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.size(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.size() - row );
               for( size_t column=0UL, n=0UL; column<rhs_.size(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, rhs_.size() - column );
                  submatrix( dres_  , row, column, m, n ) += submatrix( eval( lhs_ ) * eval( rhs_ )      , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) += submatrix( eval( lhs_ ) * eval( rhs_ )      , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) += submatrix( eval( lhs_ ) * eval( rhs_ )      , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) += submatrix( eval( lhs_ ) * eval( rhs_ )      , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) += submatrix( eval( reflhs_ ) * eval( refrhs_ ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }


      //=====================================================================================
      // Submatrix-wise outer product with subtraction assignment
      //=====================================================================================

      // Submatrix-wise outer product with subtraction assignment with the given vectors
      {
         test_  = "Submatrix-wise outer product with subtraction assignment with the given vectors";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.size(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.size() - row );
               for( size_t column=0UL, n=0UL; column<rhs_.size(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, rhs_.size() - column );
                  submatrix( dres_  , row, column, m, n ) -= submatrix( lhs_ * rhs_      , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) -= submatrix( lhs_ * rhs_      , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) -= submatrix( lhs_ * rhs_      , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) -= submatrix( lhs_ * rhs_      , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) -= submatrix( reflhs_ * refrhs_, row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Submatrix-wise outer product with subtraction assignment with evaluated vectors
      {
         test_  = "Submatrix-wise outer product with subtraction assignment with evaluated vectors";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.size(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.size() - row );
               for( size_t column=0UL, n=0UL; column<rhs_.size(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, rhs_.size() - column );
                  submatrix( dres_  , row, column, m, n ) -= submatrix( eval( lhs_ ) * eval( rhs_ )      , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) -= submatrix( eval( lhs_ ) * eval( rhs_ )      , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) -= submatrix( eval( lhs_ ) * eval( rhs_ )      , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) -= submatrix( eval( lhs_ ) * eval( rhs_ )      , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) -= submatrix( eval( reflhs_ ) * eval( refrhs_ ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the row-wise sparse vector/dense vector outer product.
//
// \return void
// \exception std::runtime_error Outer product error detected.
//
// This function tests the row-wise outer product with plain assignment, addition assignment,
// and subtraction assignment. In case any error resulting from the outer product or the
// subsequent assignment is detected, a \a std::runtime_error exception is thrown.
*/
template< typename VT1    // Type of the left-hand side dense vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::testRowOperation()
{
#if BLAZETEST_MATHTEST_TEST_ROW_OPERATION
   if( BLAZETEST_MATHTEST_TEST_ROW_OPERATION > 1 )
   {
      //=====================================================================================
      // Row-wise outer product
      //=====================================================================================

      // Row-wise outer product with the given vectors
      {
         test_  = "Row-wise outer product with the given vectors";
         error_ = "Failed outer product operation";

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.size(); ++i ) {
               row( dres_  , i ) = row( lhs_ * rhs_, i );
               row( odres_ , i ) = row( lhs_ * rhs_, i );
               row( sres_  , i ) = row( lhs_ * rhs_, i );
               row( osres_ , i ) = row( lhs_ * rhs_, i );
               row( refres_, i ) = row( reflhs_ * refrhs_, i );
            }
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Row-wise outer product with evaluated vectors
      {
         test_  = "Row-wise outer product with the given vectors";
         error_ = "Failed outer product operation";

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.size(); ++i ) {
               row( dres_  , i ) = row( eval( lhs_ ) * eval( rhs_ ), i );
               row( odres_ , i ) = row( eval( lhs_ ) * eval( rhs_ ), i );
               row( sres_  , i ) = row( eval( lhs_ ) * eval( rhs_ ), i );
               row( osres_ , i ) = row( eval( lhs_ ) * eval( rhs_ ), i );
               row( refres_, i ) = row( eval( reflhs_ ) * eval( refrhs_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }


      //=====================================================================================
      // Row-wise outer product with addition assignment
      //=====================================================================================

      // Row-wise outer product with addition assignment with the given vectors
      {
         test_  = "Row-wise outer product with addition assignment with the given vectors";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.size(); ++i ) {
               row( dres_  , i ) += row( lhs_ * rhs_, i );
               row( odres_ , i ) += row( lhs_ * rhs_, i );
               row( sres_  , i ) += row( lhs_ * rhs_, i );
               row( osres_ , i ) += row( lhs_ * rhs_, i );
               row( refres_, i ) += row( reflhs_ * refrhs_, i );
            }
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Row-wise outer product with addition assignment with evaluated vectors
      {
         test_  = "Row-wise outer product with addition assignment with the given vectors";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.size(); ++i ) {
               row( dres_  , i ) += row( eval( lhs_ ) * eval( rhs_ ), i );
               row( odres_ , i ) += row( eval( lhs_ ) * eval( rhs_ ), i );
               row( sres_  , i ) += row( eval( lhs_ ) * eval( rhs_ ), i );
               row( osres_ , i ) += row( eval( lhs_ ) * eval( rhs_ ), i );
               row( refres_, i ) += row( eval( reflhs_ ) * eval( refrhs_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }


      //=====================================================================================
      // Row-wise outer product with subtraction assignment
      //=====================================================================================

      // Row-wise outer product with subtraction assignment with the given vectors
      {
         test_  = "Row-wise outer product with subtraction assignment with the given vectors";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.size(); ++i ) {
               row( dres_  , i ) -= row( lhs_ * rhs_, i );
               row( odres_ , i ) -= row( lhs_ * rhs_, i );
               row( sres_  , i ) -= row( lhs_ * rhs_, i );
               row( osres_ , i ) -= row( lhs_ * rhs_, i );
               row( refres_, i ) -= row( reflhs_ * refrhs_, i );
            }
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Row-wise outer product with subtraction assignment with evaluated vectors
      {
         test_  = "Row-wise outer product with subtraction assignment with the given vectors";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.size(); ++i ) {
               row( dres_  , i ) -= row( eval( lhs_ ) * eval( rhs_ ), i );
               row( odres_ , i ) -= row( eval( lhs_ ) * eval( rhs_ ), i );
               row( sres_  , i ) -= row( eval( lhs_ ) * eval( rhs_ ), i );
               row( osres_ , i ) -= row( eval( lhs_ ) * eval( rhs_ ), i );
               row( refres_, i ) -= row( eval( reflhs_ ) * eval( refrhs_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the column-wise sparse vector/dense vector outer product.
//
// \return void
// \exception std::runtime_error Outer product error detected.
//
// This function tests the column-wise outer product with plain assignment, addition assignment,
// and subtraction assignment. In case any error resulting from the outer product or the
// subsequent assignment is detected, a \a std::runtime_error exception is thrown.
*/
template< typename VT1    // Type of the left-hand side dense vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::testColumnOperation()
{
#if BLAZETEST_MATHTEST_TEST_COLUMN_OPERATION
   if( BLAZETEST_MATHTEST_TEST_COLUMN_OPERATION > 1 )
   {
      //=====================================================================================
      // Column-wise outer product
      //=====================================================================================

      // Column-wise outer product with the given vectors
      {
         test_  = "Column-wise outer product with the given vectors";
         error_ = "Failed outer product operation";

         try {
            initResults();
            for( size_t i=0UL; i<rhs_.size(); ++i ) {
               column( dres_  , i ) = column( lhs_ * rhs_, i );
               column( odres_ , i ) = column( lhs_ * rhs_, i );
               column( sres_  , i ) = column( lhs_ * rhs_, i );
               column( osres_ , i ) = column( lhs_ * rhs_, i );
               column( refres_, i ) = column( reflhs_ * refrhs_, i );
            }
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Column-wise outer product with evaluated vectors
      {
         test_  = "Column-wise outer product with the given vectors";
         error_ = "Failed outer product operation";

         try {
            initResults();
            for( size_t i=0UL; i<rhs_.size(); ++i ) {
               column( dres_  , i ) = column( eval( lhs_ ) * eval( rhs_ ), i );
               column( odres_ , i ) = column( eval( lhs_ ) * eval( rhs_ ), i );
               column( sres_  , i ) = column( eval( lhs_ ) * eval( rhs_ ), i );
               column( osres_ , i ) = column( eval( lhs_ ) * eval( rhs_ ), i );
               column( refres_, i ) = column( eval( reflhs_ ) * eval( refrhs_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }


      //=====================================================================================
      // Column-wise outer product with addition assignment
      //=====================================================================================

      // Column-wise outer product with addition assignment with the given vectors
      {
         test_  = "Column-wise outer product with addition assignment with the given vectors";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            for( size_t i=0UL; i<rhs_.size(); ++i ) {
               column( dres_  , i ) += column( lhs_ * rhs_, i );
               column( odres_ , i ) += column( lhs_ * rhs_, i );
               column( sres_  , i ) += column( lhs_ * rhs_, i );
               column( osres_ , i ) += column( lhs_ * rhs_, i );
               column( refres_, i ) += column( reflhs_ * refrhs_, i );
            }
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Column-wise outer product with addition assignment with evaluated vectors
      {
         test_  = "Column-wise outer product with addition assignment with the given vectors";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            for( size_t i=0UL; i<rhs_.size(); ++i ) {
               column( dres_  , i ) += column( eval( lhs_ ) * eval( rhs_ ), i );
               column( odres_ , i ) += column( eval( lhs_ ) * eval( rhs_ ), i );
               column( sres_  , i ) += column( eval( lhs_ ) * eval( rhs_ ), i );
               column( osres_ , i ) += column( eval( lhs_ ) * eval( rhs_ ), i );
               column( refres_, i ) += column( eval( reflhs_ ) * eval( refrhs_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }


      //=====================================================================================
      // Column-wise outer product with subtraction assignment
      //=====================================================================================

      // Column-wise outer product with subtraction assignment with the given vectors
      {
         test_  = "Column-wise outer product with subtraction assignment with the given vectors";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            for( size_t i=0UL; i<rhs_.size(); ++i ) {
               column( dres_  , i ) -= column( lhs_ * rhs_, i );
               column( odres_ , i ) -= column( lhs_ * rhs_, i );
               column( sres_  , i ) -= column( lhs_ * rhs_, i );
               column( osres_ , i ) -= column( lhs_ * rhs_, i );
               column( refres_, i ) -= column( reflhs_ * refrhs_, i );
            }
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Column-wise outer product with subtraction assignment with evaluated vectors
      {
         test_  = "Column-wise outer product with subtraction assignment with the given vectors";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            for( size_t i=0UL; i<rhs_.size(); ++i ) {
               column( dres_  , i ) -= column( eval( lhs_ ) * eval( rhs_ ), i );
               column( odres_ , i ) -= column( eval( lhs_ ) * eval( rhs_ ), i );
               column( sres_  , i ) -= column( eval( lhs_ ) * eval( rhs_ ), i );
               column( osres_ , i ) -= column( eval( lhs_ ) * eval( rhs_ ), i );
               column( refres_, i ) -= column( eval( reflhs_ ) * eval( refrhs_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
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
// This function is called after each test case to check and compare the computed results.
*/
template< typename VT1    // Type of the left-hand side sparse vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::checkResults()
{
   if( !isEqual( dres_, refres_ ) || !isEqual( odres_, refres_ ) ) {
      std::ostringstream oss;
      oss.precision( 20 );
      oss << " Test : " << test_ << "\n"
          << " Error: Incorrect dense result detected\n"
          << " Details:\n"
          << "   Left-hand side sparse vector type:\n"
          << "     " << typeid( VT1 ).name() << "\n"
          << "   Right-hand side transpose dense vector type:\n"
          << "     " << typeid( TVT2 ).name() << "\n"
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
          << "   Left-hand side sparse vector type:\n"
          << "     " << typeid( VT1 ).name() << "\n"
          << "   Right-hand side transpose dense vector type:\n"
          << "     " << typeid( TVT2 ).name() << "\n"
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
// results.
*/
template< typename VT1    // Type of the left-hand side sparse vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::checkTransposeResults()
{
   if( !isEqual( tdres_, refres_ ) || !isEqual( todres_, refres_ ) ) {
      std::ostringstream oss;
      oss.precision( 20 );
      oss << " Test : " << test_ << "\n"
          << " Error: Incorrect transpose dense result detected\n"
          << " Details:\n"
          << "   Left-hand side sparse vector type:\n"
          << "     " << typeid( VT1 ).name() << "\n"
          << "   Right-hand side transpose dense vector type:\n"
          << "     " << typeid( TVT2 ).name() << "\n"
          << "   Transpose result:\n" << tdres_ << "\n"
          << "   Transpose result with opposite storage order:\n" << todres_ << "\n"
          << "   Expected result:\n" << refres_ << "\n";
      throw std::runtime_error( oss.str() );
   }

   if( !isEqual( tsres_, refres_ ) || !isEqual( tosres_, refres_ ) ) {
      std::ostringstream oss;
      oss.precision( 20 );
      oss << " Test : " << test_ << "\n"
          << " Error: Incorrect transpose sparse result detected\n"
          << " Details:\n"
          << "   Left-hand side sparse vector type:\n"
          << "     " << typeid( VT1 ).name() << "\n"
          << "   Right-hand side transpose dense vector type:\n"
          << "     " << typeid( TVT2 ).name() << "\n"
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
template< typename VT1    // Type of the left-hand side sparse vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::initResults()
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
template< typename VT1    // Type of the left-hand side sparse vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::initTransposeResults()
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
// test.
*/
template< typename VT1    // Type of the left-hand side sparse vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::convertException( const std::exception& ex )
{
   std::ostringstream oss;
   oss << " Test : " << test_ << "\n"
       << " Error: " << error_ << "\n"
       << " Details:\n"
       << "   Left-hand side sparse vector type:\n"
       << "     " << typeid( VT1 ).name() << "\n"
       << "   Right-hand side transpose dense vector type:\n"
       << "     " << typeid( TVT2 ).name() << "\n"
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
/*!\brief Testing the vector outer product between two specific vector types.
//
// \param creator1 The creator for the left-hand side vector.
// \param creator2 The creator for the right-hand side vector.
// \return void
*/
template< typename VT1    // Type of the left-hand side sparse vector
        , typename VT2 >  // Type of the right-hand side dense vector
void runTest( const Creator<VT1>& creator1, const Creator<VT2>& creator2 )
{
   for( size_t rep=0; rep<repetitions; ++rep ) {
      OperationTest<VT1,VT2>( creator1, creator2 );
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
/*!\brief Macro for the definition of a sparse vector/dense vector outer product test case.
*/
#define DEFINE_SVECTDVECMULT_OPERATION_TEST( VT1, VT2 ) \
   extern template class blazetest::mathtest::svectdvecmult::OperationTest<VT1,VT2>
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Macro for the execution of a sparse vector/dense vector outer product test case.
*/
#define RUN_SVECTDVECMULT_OPERATION_TEST( C1, C2 ) \
   blazetest::mathtest::svectdvecmult::runTest( C1, C2 )
/*! \endcond */
//*************************************************************************************************

} // namespace svectdvecmult

} // namespace mathtest

} // namespace blazetest

#endif
