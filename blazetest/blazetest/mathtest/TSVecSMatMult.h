//=================================================================================================
/*!
//  \file blazetest/mathtest/TSVecSMatMult.h
//  \brief Header file for the sparse vector/sparse matrix multiplication math test
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

#ifndef _BLAZETEST_MATHTEST_TSVECSMATMULT_H_
#define _BLAZETEST_MATHTEST_TSVECSMATMULT_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <sstream>
#include <stdexcept>
#include <string>
#include <typeinfo>
#include <blaze/math/CompressedMatrix.h>
#include <blaze/math/constraints/DenseMatrix.h>
#include <blaze/math/constraints/DenseVector.h>
#include <blaze/math/constraints/SparseMatrix.h>
#include <blaze/math/constraints/SparseVector.h>
#include <blaze/math/constraints/StorageOrder.h>
#include <blaze/math/constraints/TransposeFlag.h>
#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/MathTrait.h>
#include <blaze/math/StaticMatrix.h>
#include <blaze/math/typetraits/IsRowMajorMatrix.h>
#include <blaze/util/constraints/Numeric.h>
#include <blaze/util/constraints/SameType.h>
#include <blazetest/system/MathTest.h>
#include <blazetest/util/Creator.h>
#include <blazetest/util/Utility.h>


namespace blazetest {

namespace mathtest {

namespace tsvecsmatmult {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Auxiliary class template for the sparse vector/sparse matrix multiplication math test.
//
// The TSVecSMatMult class template represents one particular vector/matrix multiplication test
// between a vector and a matrix of particular types. The two template arguments \a VT and \a MT
// represent the types of the left-hand side vector and right-hand side matrix, respectively.
*/
template< typename VT    // Type of the left-hand side sparse vector
        , typename MT >  // Type of the right-hand side sparse matrix
class TSVecSMatMult
{
 private:
   //**Type definitions****************************************************************************
   typedef typename VT::TransposeType                   TVT;  //!< Transpose vector type
   typedef typename MT::OppositeType                    OMT;  //!< Matrix type with opposite storage order
   typedef typename MT::TransposeType                   TMT;  //!< Transpose matrix type
   typedef typename blaze::MathTrait<TVT,MT>::MultType  RE;   //!< Default result type
   typedef typename RE::TransposeType                   TRE;  //!< Transpose default result type
   //**********************************************************************************************

   //**Type definitions****************************************************************************
   typedef typename VT::ElementType           VET;    //!< Element type of the vector type.
   typedef typename MT::ElementType           MET;    //!< Element type of the matrix type.
   typedef typename RE::ElementType           RET;    //!< Resulting element type
   typedef blaze::DynamicVector<VET,true>     VRT;    //!< Vector reference type
   typedef blaze::DynamicMatrix<MET,false>    MRT;    //!< Matrix reference type
   typedef blaze::DynamicVector<RET,true>     DRRE;   //!< Dense reference result type
   typedef blaze::CompressedVector<RET,true>  SRRE;   //!< Sparse reference result type
   typedef typename DRRE::TransposeType       TDRRE;  //!< Transpose dense reference result type
   typedef typename SRRE::TransposeType       TSRRE;  //!< Transpose sparse reference result type
   typedef DRRE                               DRE;    //!< Dense result type
   typedef RE                                 SRE;    //!< Sparse result type
   typedef TDRRE                              TDRE;   //!< Transpose dense result type
   typedef TRE                                TSRE;   //!< Transpose sparse result type
   //**********************************************************************************************

 public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit TSVecSMatMult( const Creator<VT>& creator1, const Creator<MT>& creator2 );
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
   //@}
   //**********************************************************************************************

   //**Error detection functions*******************************************************************
   /*!\name Error detection functions */
   //@{
   template< typename RT > void checkResults();
   template< typename RT > void checkTransposeResults();
   //@}
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   TVT   lhs_;      //!< The left-hand side sparse vector.
   MT    rhs_;      //!< The right-hand side sparse matrix.
   DRE   dres_;     //!< The dense result vector.
   SRE   sres_;     //!< The sparse result vector.
   VRT   reflhs_;   //!< The reference left-hand side vector.
   MRT   refrhs_;   //!< The reference right-hand side matrix.
   DRRE  refres_;   //!< The reference result.
   OMT   orhs_;     //!< The right-hand side sparse matrix with opposite storage order.
   TDRE  tdres_;    //!< The transpose dense result vector.
   TSRE  tsres_;    //!< The transpose sparse result vector.
   TDRRE trefres_;  //!< The transpose reference result.

   std::string test_;  //!< Label of the currently performed test.
   //@}
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( VT    );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( MT    );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( TVT   );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( TMT   );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( VRT   );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE ( MRT   );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( DRE   );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( SRE   );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( DRRE  );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( SRRE  );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( TDRE  );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( TSRE  );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( TDRRE );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( TSRRE );
   BLAZE_CONSTRAINT_MUST_BE_NONTRANSPOSE_VECTOR_TYPE( VT    );
   BLAZE_CONSTRAINT_MUST_BE_TRANSPOSE_VECTOR_TYPE   ( TVT   );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE   ( MT    );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( TMT   );
   BLAZE_CONSTRAINT_MUST_BE_TRANSPOSE_VECTOR_TYPE   ( VRT   );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE   ( MRT   );
   BLAZE_CONSTRAINT_MUST_BE_TRANSPOSE_VECTOR_TYPE   ( DRE   );
   BLAZE_CONSTRAINT_MUST_BE_TRANSPOSE_VECTOR_TYPE   ( SRE   );
   BLAZE_CONSTRAINT_MUST_BE_TRANSPOSE_VECTOR_TYPE   ( DRRE  );
   BLAZE_CONSTRAINT_MUST_BE_TRANSPOSE_VECTOR_TYPE   ( SRRE  );
   BLAZE_CONSTRAINT_MUST_BE_NONTRANSPOSE_VECTOR_TYPE( TDRE  );
   BLAZE_CONSTRAINT_MUST_BE_NONTRANSPOSE_VECTOR_TYPE( TSRE  );
   BLAZE_CONSTRAINT_MUST_BE_NONTRANSPOSE_VECTOR_TYPE( TDRRE );
   BLAZE_CONSTRAINT_MUST_BE_NONTRANSPOSE_VECTOR_TYPE( TSRRE );
   BLAZE_CONSTRAINT_MATRICES_MUST_HAVE_DIFFERENT_STORAGE_ORDER( MT, TMT );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MET, typename OMT::ElementType   );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MET, typename TMT::ElementType   );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( VET, typename TVT::ElementType   );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MT , typename OMT::OppositeType  );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MT , typename TMT::TransposeType );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( VT , typename TVT::TransposeType );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( RE , typename TRE::TransposeType );
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
/*!\brief Constructor for the TSVecSMatMult class template.
//
// \param creator1 The creator for the left-hand side sparse vector of the multiplication.
// \param creator2 The creator for the right-hand side sparse matrix of the multiplication.
// \exception std::runtime_error Operation error detected.
*/
template< typename VT    // Type of the left-hand side sparse vector
        , typename MT >  // Type of the right-hand side sparse matrix
TSVecSMatMult<VT,MT>::TSVecSMatMult( const Creator<VT>& creator1, const Creator<MT>& creator2 )
   : lhs_ ( trans( creator1() ) )  // The left-hand side sparse vector
   , rhs_ ( creator2() )           // The right-hand side sparse matrix
   , dres_()                       // The dense result vector
   , sres_()                       // The sparse result vector
   , reflhs_( lhs_ )               // The reference left-hand side vector
   , refrhs_( rhs_ )               // The reference right-hand side matrix
   , refres_()                     // The reference result
   , orhs_( rhs_ )                 // The right-hand side sparse matrix with transpose storage order.
   , tdres_()                      // The transpose dense result vector.
   , tsres_()                      // The transpose sparse result vector.
   , trefres_()                    // The transpose reference result.
   , test_()                       // Label of the currently performed test
{
   testInitialStatus();
   testAssignment();
   testElementAccess();
   testBasicOperation();
   testNegatedOperation();
   testScaledOperation( 2 );
   testScaledOperation( 2UL );
   testScaledOperation( 1.1F );
   testScaledOperation( 1.1 );
   testTransposeOperation();
   testAbsOperation();
}
//*************************************************************************************************




//=================================================================================================
//
//  TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Tests on the initial status of the operands.
//
// \return void
// \exception std::runtime_error Initialization error detected.
//
// This function runs tests on the initial status of the operands. In case any initialization
// error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename VT    // Type of the left-hand side sparse vector
        , typename MT >  // Type of the right-hand side sparse matrix
void TSVecSMatMult<VT,MT>::testInitialStatus()
{
   //=====================================================================================
   // Performing initial tests with the given types
   //=====================================================================================

   // Checking the size of the left-hand side operand
   if( lhs_.size() != reflhs_.size() ) {
      std::ostringstream oss;
      oss << " Test: Initial size comparison of left-hand side sparse operand\n"
          << " Error: Invalid vector size\n"
          << " Details:\n"
          << "   Transpose sparse vector type:\n"
          << "     " << typeid( TVT ).name() << "\n"
          << "   Detected size = " << lhs_.size() << "\n"
          << "   Expected size = " << reflhs_.size() << "\n";
      throw std::runtime_error( oss.str() );
   }

   // Checking the number of rows of the right-hand side operand
   if( rhs_.rows() != refrhs_.rows() ) {
      std::ostringstream oss;
      oss << " Test: Initial size comparison of right-hand side sparse operand\n"
          << " Error: Invalid number of rows\n"
          << " Details:\n"
          << "   Row-major sparse matrix type:\n"
          << "     " << typeid( MT ).name() << "\n"
          << "   Detected number of rows = " << rhs_.rows() << "\n"
          << "   Expected number of rows = " << refrhs_.rows() << "\n";
      throw std::runtime_error( oss.str() );
   }

   // Checking the number of columns of the right-hand side operand
   if( rhs_.columns() != refrhs_.columns() ) {
      std::ostringstream oss;
      oss << " Test: Initial size comparison of right-hand side sparse operand\n"
          << " Error: Invalid number of columns\n"
          << " Details:\n"
          << "   Row-major sparse matrix type:\n"
          << "     " << typeid( MT ).name() << "\n"
          << "   Detected number of columns = " << rhs_.columns() << "\n"
          << "   Expected number of columns = " << refrhs_.columns() << "\n";
      throw std::runtime_error( oss.str() );
   }

   // Checking the initialization of the left-hand side operand
   if( !isEqual( lhs_, reflhs_ ) ) {
      std::ostringstream oss;
      oss << " Test: Initial test of initialization of left-hand side sparse operand\n"
          << " Error: Invalid vector initialization\n"
          << " Details:\n"
          << "   Transpose sparse vector type:\n"
          << "     " << typeid( TVT ).name() << "\n"
          << "   Current initialization:\n" << lhs_ << "\n"
          << "   Expected initialization:\n" << reflhs_ << "\n";
      throw std::runtime_error( oss.str() );
   }

   // Checking the initialization of the right-hand side operand
   if( !isEqual( rhs_, refrhs_ ) ) {
      std::ostringstream oss;
      oss << " Test: Initial test of initialization of right-hand side sparse operand\n"
          << " Error: Invalid matrix initialization\n"
          << " Details:\n"
          << "   Row-major sparse matrix type:\n"
          << "     " << typeid( MT ).name() << "\n"
          << "   Current initialization:\n" << rhs_ << "\n"
          << "   Expected initialization:\n" << refrhs_ << "\n";
      throw std::runtime_error( oss.str() );
   }


   //=====================================================================================
   // Performing initial tests with the transpose types
   //=====================================================================================

   // Checking the number of rows of the transpose right-hand side operand
   if( orhs_.rows() != refrhs_.rows() ) {
      std::ostringstream oss;
      oss << " Test: Initial size comparison of transpose right-hand side sparse operand\n"
          << " Error: Invalid number of rows\n"
          << " Details:\n"
          << "   Column-major sparse matrix type:\n"
          << "     " << typeid( TMT ).name() << "\n"
          << "   Detected number of rows = " << orhs_.rows() << "\n"
          << "   Expected number of rows = " << refrhs_.rows() << "\n";
      throw std::runtime_error( oss.str() );
   }

   // Checking the number of columns of the transpose right-hand side operand
   if( orhs_.columns() != refrhs_.columns() ) {
      std::ostringstream oss;
      oss << " Test: Initial size comparison of transpose right-hand side sparse operand\n"
          << " Error: Invalid number of columns\n"
          << " Details:\n"
          << "   Column-major sparse matrix type:\n"
          << "     " << typeid( TMT ).name() << "\n"
          << "   Detected number of columns = " << orhs_.columns() << "\n"
          << "   Expected number of columns = " << refrhs_.columns() << "\n";
      throw std::runtime_error( oss.str() );
   }

   // Checking the initialization of the transpose right-hand side operand
   if( !isEqual( orhs_, refrhs_ ) ) {
      std::ostringstream oss;
      oss << " Test: Initial test of initialization of transpose right-hand side sparse operand\n"
          << " Error: Invalid matrix initialization\n"
          << " Details:\n"
          << "   Column-major sparse matrix type:\n"
          << "     " << typeid( TMT ).name() << "\n"
          << "   Current initialization:\n" << orhs_ << "\n"
          << "   Expected initialization:\n" << refrhs_ << "\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the vector assignment.
//
// \return void
// \exception std::runtime_error Assignment error detected.
//
// This function tests the vector assignment. In case any error is detected, a
// \a std::runtime_error exception is thrown.
*/
template< typename VT    // Type of the left-hand side sparse vector
        , typename MT >  // Type of the right-hand side sparse matrix
void TSVecSMatMult<VT,MT>::testAssignment()
{
   //=====================================================================================
   // Performing an assignment with the given types
   //=====================================================================================

   try {
      lhs_ = reflhs_;
      rhs_ = refrhs_;
   }
   catch( std::exception& ex ) {
      std::ostringstream oss;
      oss << " Test: Assignment with the given types\n"
          << " Error: Failed assignment\n"
          << " Details:\n"
          << "   Left-hand side transpose sparse vector type:\n"
          << "     " << typeid( TVT ).name() << "\n"
          << "   Right-hand side row-major sparse matrix type:\n"
          << "     " << typeid( MT ).name() << "\n"
          << "   Error message: " << ex.what() << "\n";
      throw std::runtime_error( oss.str() );
   }

   if( !isEqual( lhs_, reflhs_ ) ) {
      std::ostringstream oss;
      oss << " Test: Checking the assignment result of left-hand side sparse operand\n"
          << " Error: Invalid matrix initialization\n"
          << " Details:\n"
          << "   Transpose sparse vector type:\n"
          << "     " << typeid( TVT ).name() << "\n"
          << "   Current initialization:\n" << lhs_ << "\n"
          << "   Expected initialization:\n" << reflhs_ << "\n";
      throw std::runtime_error( oss.str() );
   }

   if( !isEqual( rhs_, refrhs_ ) ) {
      std::ostringstream oss;
      oss << " Test: Checking the assignment result of right-hand side sparse operand\n"
          << " Error: Invalid vector initialization\n"
          << " Details:\n"
          << "   Row-major sparse matrix type:\n"
          << "     " << typeid( MT ).name() << "\n"
          << "   Current initialization:\n" << rhs_ << "\n"
          << "   Expected initialization:\n" << refrhs_ << "\n";
      throw std::runtime_error( oss.str() );
   }


   //=====================================================================================
   // Performing an assignment with the transpose types
   //=====================================================================================

   try {
      orhs_ = refrhs_;
   }
   catch( std::exception& ex ) {
      std::ostringstream oss;
      oss << " Test: Assignment with the transpose types\n"
          << " Error: Failed assignment\n"
          << " Details:\n"
          << "   Right-hand side column-major sparse matrix type:\n"
          << "     " << typeid( TMT ).name() << "\n"
          << "   Error message: " << ex.what() << "\n";
      throw std::runtime_error( oss.str() );
   }

   if( !isEqual( orhs_, refrhs_ ) ) {
      std::ostringstream oss;
      oss << " Test: Checking the assignment result of transpose right-hand side sparse operand\n"
          << " Error: Invalid matrix initialization\n"
          << " Details:\n"
          << "   Column-major sparse matrix type:\n"
          << "     " << typeid( TMT ).name() << "\n"
          << "   Current initialization:\n" << orhs_ << "\n"
          << "   Expected initialization:\n" << refrhs_ << "\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the vector element access.
//
// \return void
// \exception std::runtime_error Element access error detected.
//
// This function tests the element access via the subscript operator. In case any
// error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename VT    // Type of the left-hand side sparse vector
        , typename MT >  // Type of the right-hand side sparse matrix
void TSVecSMatMult<VT,MT>::testElementAccess()
{
   using blaze::equal;


   //=====================================================================================
   // Testing the element access with the given types
   //=====================================================================================

   if( rhs_.columns() > 0UL )
   {
      if( !equal( ( lhs_ * rhs_ )[0UL], ( reflhs_ * refrhs_ )[0UL] ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of multiplication expression\n"
             << " Error: Unequal resulting elements at index 0 detected\n"
             << " Details:\n"
             << "   Left-hand side transpose sparse vector type:\n"
             << "     " << typeid( TVT ).name() << "\n"
             << "   Right-hand side row-major sparse matrix type:\n"
             << "     " << typeid( MT ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( ( lhs_ * eval( rhs_ ) )[0UL], ( reflhs_ * eval( refrhs_ ) )[0UL] ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of right evaluated multiplication expression\n"
             << " Error: Unequal resulting elements at index 0 detected\n"
             << " Details:\n"
             << "   Left-hand side transpose sparse vector type:\n"
             << "     " << typeid( TVT ).name() << "\n"
             << "   Right-hand side row-major sparse matrix type:\n"
             << "     " << typeid( MT ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( ( eval( lhs_ ) * rhs_ )[0UL], ( eval( reflhs_ ) * refrhs_ )[0UL] ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of left evaluated multiplication expression\n"
             << " Error: Unequal resulting elements at index 0 detected\n"
             << " Details:\n"
             << "   Left-hand side transpose sparse vector type:\n"
             << "     " << typeid( TVT ).name() << "\n"
             << "   Right-hand side row-major sparse matrix type:\n"
             << "     " << typeid( MT ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( ( eval( lhs_ ) * eval( rhs_ ) )[0UL], ( eval( reflhs_ ) * eval( refrhs_ ) )[0UL] ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of fully evaluated multiplication expression\n"
             << " Error: Unequal resulting elements at index 0 detected\n"
             << " Details:\n"
             << "   Left-hand side transpose sparse vector type:\n"
             << "     " << typeid( TVT ).name() << "\n"
             << "   Right-hand side row-major sparse matrix type:\n"
             << "     " << typeid( MT ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Testing the element access with the transpose types
   //=====================================================================================

   if( orhs_.columns() > 0UL )
   {
      if( !equal( ( lhs_ * orhs_ )[0UL], ( reflhs_ * refrhs_ )[0UL] ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of transpose multiplication expression\n"
             << " Error: Unequal resulting elements at index 0 detected\n"
             << " Details:\n"
             << "   Left-hand side transpose sparse vector type:\n"
             << "     " << typeid( TVT ).name() << "\n"
             << "   Right-hand side column-major sparse matrix type:\n"
             << "     " << typeid( TMT ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( ( lhs_ * eval( orhs_ ) )[0UL], ( reflhs_ * eval( refrhs_ ) )[0UL] ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of right evaluated transpose multiplication expression\n"
             << " Error: Unequal resulting elements at index 0 detected\n"
             << " Details:\n"
             << "   Left-hand side transpose sparse vector type:\n"
             << "     " << typeid( TVT ).name() << "\n"
             << "   Right-hand side column-major sparse matrix type:\n"
             << "     " << typeid( TMT ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( ( eval( lhs_ ) * orhs_ )[0UL], ( eval( reflhs_ ) * refrhs_ )[0UL] ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of left evaluated transpose multiplication expression\n"
             << " Error: Unequal resulting elements at index 0 detected\n"
             << " Details:\n"
             << "   Left-hand side transpose sparse vector type:\n"
             << "     " << typeid( TVT ).name() << "\n"
             << "   Right-hand side column-major sparse matrix type:\n"
             << "     " << typeid( TMT ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( ( eval( lhs_ ) * eval( orhs_ ) )[0UL], ( eval( reflhs_ ) * eval( refrhs_ ) )[0UL] ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of fully evaluated transpose multiplication expression\n"
             << " Error: Unequal resulting elements at index 0 detected\n"
             << " Details:\n"
             << "   Left-hand side transpose sparse vector type:\n"
             << "     " << typeid( TVT ).name() << "\n"
             << "   Right-hand side column-major sparse matrix type:\n"
             << "     " << typeid( TMT ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the plain sparse vector/sparse matrix multiplication.
//
// \return void
// \exception std::runtime_error Multiplication error detected.
//
// This function tests the plain vector/matrix multiplication with plain assignment,
// addition assignment, subtraction assignment, and multiplication assignment. In case
// any error resulting from the multiplication or the subsequent assignment is detected,
// a \a std::runtime_error exception is thrown.
*/
template< typename VT    // Type of the left-hand side sparse vector
        , typename MT >  // Type of the right-hand side sparse matrix
void TSVecSMatMult<VT,MT>::testBasicOperation()
{
#if BLAZETEST_MATHTEST_TEST_BASIC_OPERATION
   if( BLAZETEST_MATHTEST_TEST_BASIC_OPERATION > 1 )
   {
      //=====================================================================================
      // Multiplication
      //=====================================================================================

      // Multiplication with the given vector/matrix
      {
         test_ = "Multiplication with the given vector/matrix";

         try {
            dres_   = lhs_ * rhs_;
            sres_   = lhs_ * rhs_;
            refres_ = reflhs_ * refrhs_;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side row-major sparse matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<MT>();

         try {
            dres_ = lhs_ * orhs_;
            sres_ = lhs_ * orhs_;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side column-major sparse matrix type:\n"
                << "     " << typeid( TMT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<TMT>();
      }

      // Multiplication with evaluated vector/matrix
      {
         test_ = "Multiplication with evaluated vector/matrix";

         try {
            dres_ = eval( lhs_ ) * eval( rhs_ );
            sres_ = eval( lhs_ ) * eval( rhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side row-major sparse matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<MT>();

         try {
            dres_ = eval( lhs_ ) * eval( orhs_ );
            sres_ = eval( lhs_ ) * eval( orhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side column-major sparse matrix type:\n"
                << "     " << typeid( TMT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<TMT>();
      }


      //=====================================================================================
      // Multiplication with addition assignment
      //=====================================================================================

      // Multiplication with addition assignment with the given vector/matrix
      {
         test_ = "Multiplication with addition assignment with the given vector/matrix";

         try {
            dres_   += lhs_ * rhs_;
            sres_   += lhs_ * rhs_;
            refres_ += reflhs_ * refrhs_;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed addition assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side row-major sparse matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<MT>();

         try {
            dres_   += lhs_ * orhs_;
            sres_   += lhs_ * orhs_;
            refres_ += reflhs_ * refrhs_;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed addition assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side column-major sparse matrix type:\n"
                << "     " << typeid( TMT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<TMT>();
      }

      // Multiplication with addition assignment with evaluated vector/matrix
      {
         test_ = "Multiplication with addition assignment with evaluated vector/matrix";

         try {
            dres_   += eval( lhs_ ) * eval( rhs_ );
            sres_   += eval( lhs_ ) * eval( rhs_ );
            refres_ += eval( reflhs_ ) * eval( refrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed addition assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side row-major sparse matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<MT>();

         try {
            dres_   += eval( lhs_ ) * eval( orhs_ );
            sres_   += eval( lhs_ ) * eval( orhs_ );
            refres_ += eval( reflhs_ ) * eval( refrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed addition assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side column-major sparse matrix type:\n"
                << "     " << typeid( TMT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<TMT>();
      }


      //=====================================================================================
      // Multiplication with subtraction assignment
      //=====================================================================================

      // Multiplication with subtraction assignment with the given vector/matrix
      {
         test_ = "Multiplication with subtraction assignment with the given vector/matrix";

         try {
            dres_   -= lhs_ * rhs_;
            sres_   -= lhs_ * rhs_;
            refres_ -= reflhs_ * refrhs_;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side row-major sparse matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<MT>();

         try {
            dres_   -= lhs_ * orhs_;
            sres_   -= lhs_ * orhs_;
            refres_ -= reflhs_ * refrhs_;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side column-major sparse matrix type:\n"
                << "     " << typeid( TMT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<TMT>();
      }

      // Multiplication with subtraction assignment with evaluated vector/matrix
      {
         test_ = "Multiplication with subtraction assignment with evaluated vector/matrix";

         try {
            dres_   -= eval( lhs_ ) * eval( rhs_ );
            sres_   -= eval( lhs_ ) * eval( rhs_ );
            refres_ -= eval( reflhs_ ) * eval( refrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side row-major sparse matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<MT>();

         try {
            dres_   -= eval( lhs_ ) * eval( orhs_ );
            sres_   -= eval( lhs_ ) * eval( orhs_ );
            refres_ -= eval( reflhs_ ) * eval( refrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side column-major sparse matrix type:\n"
                << "     " << typeid( TMT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<TMT>();
      }


      //=====================================================================================
      // Multiplication with multiplication assignment
      //=====================================================================================

      // Multiplication with multiplication assignment with the given vector/matrix
      {
         test_ = "Multiplication with multiplication assignment with the given vector/matrix";

         try {
            dres_   *= lhs_ * rhs_;
            sres_   *= lhs_ * rhs_;
            refres_ *= reflhs_ * refrhs_;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side row-major sparse matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<MT>();

         try {
            dres_   *= lhs_ * orhs_;
            sres_   *= lhs_ * orhs_;
            refres_ *= reflhs_ * refrhs_;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side column-major sparse matrix type:\n"
                << "     " << typeid( TMT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<TMT>();
      }

      // Multiplication with multiplication assignment with evaluated vector/matrix
      {
         test_ = "Multiplication with multiplication assignment with evaluated vector/matrix";

         try {
            dres_   *= eval( lhs_ ) * eval( rhs_ );
            sres_   *= eval( lhs_ ) * eval( rhs_ );
            refres_ *= eval( reflhs_ ) * eval( refrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side row-major sparse matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<MT>();

         try {
            dres_   *= eval( lhs_ ) * eval( orhs_ );
            sres_   *= eval( lhs_ ) * eval( orhs_ );
            refres_ *= eval( reflhs_ ) * eval( refrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side column-major sparse matrix type:\n"
                << "     " << typeid( TMT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<TMT>();
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the negated sparse vector/sparse matrix multiplication.
//
// \return void
// \exception std::runtime_error Multiplication error detected.
//
// This function tests the negated vector/matrix multiplication with plain assignment, addition
// assignment, subtraction assignment, and multiplication assignment. In case any error resulting
// from the multiplication or the subsequent assignment is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename VT    // Type of the left-hand side sparse vector
        , typename MT >  // Type of the right-hand side sparse matrix
void TSVecSMatMult<VT,MT>::testNegatedOperation()
{
#if BLAZETEST_MATHTEST_TEST_NEGATED_OPERATION
   if( BLAZETEST_MATHTEST_TEST_NEGATED_OPERATION > 1 )
   {
      //=====================================================================================
      // Negated multiplication
      //=====================================================================================

      // Negated multiplication with the given vector/matrix
      {
         test_ = "Negated multiplication with the given vector/matrix";

         try {
            dres_   = -( lhs_ * rhs_ );
            sres_   = -( lhs_ * rhs_ );
            refres_ = -( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side row-major sparse matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<MT>();

         try {
            dres_ = -( lhs_ * orhs_ );
            sres_ = -( lhs_ * orhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side column-major sparse matrix type:\n"
                << "     " << typeid( TMT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<TMT>();
      }

      // Negated multiplication with evaluated vector/matrix
      {
         test_ = "Negated multiplication with evaluated vector/matrix";

         try {
            dres_ = -( eval( lhs_ ) * eval( rhs_ ) );
            sres_ = -( eval( lhs_ ) * eval( rhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side row-major sparse matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<MT>();

         try {
            dres_ = -( eval( lhs_ ) * eval( orhs_ ) );
            sres_ = -( eval( lhs_ ) * eval( orhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side column-major sparse matrix type:\n"
                << "     " << typeid( TMT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<TMT>();
      }


      //=====================================================================================
      // Negated multiplication with addition assignment
      //=====================================================================================

      // Negated multiplication with addition assignment with the given vector/matrix
      {
         test_ = "Negated multiplication with addition assignment with the given vector/matrix";

         try {
            dres_   += -( lhs_ * rhs_ );
            sres_   += -( lhs_ * rhs_ );
            refres_ += -( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed addition assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side row-major sparse matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<MT>();

         try {
            dres_   += -( lhs_ * orhs_ );
            sres_   += -( lhs_ * orhs_ );
            refres_ += -( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed addition assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side column-major sparse matrix type:\n"
                << "     " << typeid( TMT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<TMT>();
      }

      // Negated multiplication with addition assignment with the given vector/matrix
      {
         test_ = "Negated multiplication with addition assignment with evaluated vector/matrix";

         try {
            dres_   += -( eval( lhs_ ) * eval( rhs_ ) );
            sres_   += -( eval( lhs_ ) * eval( rhs_ ) );
            refres_ += -( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed addition assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side row-major sparse matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<MT>();

         try {
            dres_   += -( eval( lhs_ ) * eval( orhs_ ) );
            sres_   += -( eval( lhs_ ) * eval( orhs_ ) );
            refres_ += -( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed addition assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side column-major sparse matrix type:\n"
                << "     " << typeid( TMT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<TMT>();
      }


      //=====================================================================================
      // Negated multiplication with subtraction assignment
      //=====================================================================================

      // Negated multiplication with subtraction assignment with the given vector/matrix
      {
         test_ = "Negated multiplication with subtraction assignment with the given vector/matrix";

         try {
            dres_   -= -( lhs_ * rhs_ );
            sres_   -= -( lhs_ * rhs_ );
            refres_ -= -( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side row-major sparse matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<MT>();

         try {
            dres_   -= -( lhs_ * orhs_ );
            sres_   -= -( lhs_ * orhs_ );
            refres_ -= -( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side column-major sparse matrix type:\n"
                << "     " << typeid( TMT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<TMT>();
      }

      // Negated multiplication with subtraction assignment with the given vector/matrix
      {
         test_ = "Negated multiplication with subtraction assignment with evaluated vector/matrix";

         try {
            dres_   -= -( eval( lhs_ ) * eval( rhs_ ) );
            sres_   -= -( eval( lhs_ ) * eval( rhs_ ) );
            refres_ -= -( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side row-major sparse matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<MT>();

         try {
            dres_   -= -( eval( lhs_ ) * eval( orhs_ ) );
            sres_   -= -( eval( lhs_ ) * eval( orhs_ ) );
            refres_ -= -( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side column-major sparse matrix type:\n"
                << "     " << typeid( TMT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<TMT>();
      }


      //=====================================================================================
      // Negated multiplication with multiplication assignment
      //=====================================================================================

      // Negated multiplication with multiplication assignment with the given vector/matrix
      {
         test_ = "Negated multiplication with multiplication assignment with the given vector/matrix";

         try {
            dres_   *= -( lhs_ * rhs_ );
            sres_   *= -( lhs_ * rhs_ );
            refres_ *= -( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side row-major sparse matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<MT>();

         try {
            dres_   *= -( lhs_ * orhs_ );
            sres_   *= -( lhs_ * orhs_ );
            refres_ *= -( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side column-major sparse matrix type:\n"
                << "     " << typeid( TMT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<TMT>();
      }

      // Negated multiplication with multiplication assignment with the given vector/matrix
      {
         test_ = "Negated multiplication with multiplication assignment with evaluated vector/matrix";

         try {
            dres_   *= -( eval( lhs_ ) * eval( rhs_ ) );
            sres_   *= -( eval( lhs_ ) * eval( rhs_ ) );
            refres_ *= -( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side row-major sparse matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<MT>();

         try {
            dres_   *= -( eval( lhs_ ) * eval( orhs_ ) );
            sres_   *= -( eval( lhs_ ) * eval( orhs_ ) );
            refres_ *= -( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side column-major sparse matrix type:\n"
                << "     " << typeid( TMT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<TMT>();
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the scaled sparse vector/sparse matrix multiplication.
//
// \param scalar The scalar value.
// \return void
// \exception std::runtime_error Multiplication error detected.
//
// This function tests the scaled vector/matrix multiplication with plain assignment, addition
// assignment, subtraction assignment, and multiplication assignment. In case any error resulting
// from the multiplication or the subsequent assignment is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename VT    // Type of the left-hand side sparse vector
        , typename MT >  // Type of the right-hand side sparse matrix
template< typename T >   // Type of the scalar
void TSVecSMatMult<VT,MT>::testScaledOperation( T scalar )
{
   BLAZE_CONSTRAINT_MUST_BE_NUMERIC_TYPE( T );

   if( scalar == T(0) )
      throw std::invalid_argument( "Invalid scalar parameter" );


#if BLAZETEST_MATHTEST_TEST_SCALED_OPERATION
   if( BLAZETEST_MATHTEST_TEST_SCALED_OPERATION > 1 )
   {
      //=====================================================================================
      // Self-scaling (OP*=s)
      //=====================================================================================

      // Self-scaling (OP*=s)
      {
         test_ = "Self-scaling (OP*=s)";

         try {
            dres_   = lhs_ * rhs_;
            sres_   = dres_;
            refres_ = dres_;

            dres_   *= scalar;
            sres_   *= scalar;
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

         checkResults<MT>();
      }


      //=====================================================================================
      // Self-scaling (OP/=s)
      //=====================================================================================

      // Self-scaling (OP/=s)
      {
         test_ = "Self-scaling (OP/=s)";

         try {
            dres_   = lhs_ * rhs_;
            sres_   = dres_;
            refres_ = dres_;

            dres_   /= scalar;
            sres_   /= scalar;
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

         checkResults<MT>();
      }


      //=====================================================================================
      // Scaled multiplication (s*OP)
      //=====================================================================================

      // Scaled multiplication with the given vector/matrix
      {
         test_ = "Scaled multiplication with the given vector/matrix (s*OP)";

         try {
            dres_   = scalar * ( lhs_ * rhs_ );
            sres_   = scalar * ( lhs_ * rhs_ );
            refres_ = scalar * ( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side row-major sparse matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<MT>();

         try {
            dres_ = scalar * ( lhs_ * orhs_ );
            sres_ = scalar * ( lhs_ * orhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side column-major sparse matrix type:\n"
                << "     " << typeid( TMT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<TMT>();
      }

      // Scaled multiplication with evaluated vector/matrix
      {
         test_ = "Scaled multiplication with evaluated vector/matrix (s*OP)";

         try {
            dres_ = scalar * ( eval( lhs_ ) * eval( rhs_ ) );
            sres_ = scalar * ( eval( lhs_ ) * eval( rhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side row-major sparse matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<MT>();

         try {
            dres_ = scalar * ( eval( lhs_ ) * eval( orhs_ ) );
            sres_ = scalar * ( eval( lhs_ ) * eval( orhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side column-major sparse matrix type:\n"
                << "     " << typeid( TMT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<TMT>();
      }


      //=====================================================================================
      // Scaled multiplication (OP*s)
      //=====================================================================================

      // Scaled multiplication with the given vector/matrix
      {
         test_ = "Scaled multiplication with the given vector/matrix (OP*s)";

         try {
            dres_   = ( lhs_ * rhs_ ) * scalar;
            sres_   = ( lhs_ * rhs_ ) * scalar;
            refres_ = ( reflhs_ * refrhs_ ) * scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side row-major sparse matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<MT>();

         try {
            dres_ = ( lhs_ * orhs_ ) * scalar;
            sres_ = ( lhs_ * orhs_ ) * scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side column-major sparse matrix type:\n"
                << "     " << typeid( TMT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<TMT>();
      }

      // Scaled multiplication with evaluated vector/matrix
      {
         test_ = "Scaled multiplication with evaluated vector/matrix (OP*s)";

         try {
            dres_  = ( eval( lhs_ ) * eval( rhs_ ) ) * scalar;
            sres_  = ( eval( lhs_ ) * eval( rhs_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side row-major sparse matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<MT>();

         try {
            dres_ = ( eval( lhs_ ) * eval( orhs_ ) ) * scalar;
            sres_ = ( eval( lhs_ ) * eval( orhs_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side column-major sparse matrix type:\n"
                << "     " << typeid( TMT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<TMT>();
      }


      //=====================================================================================
      // Scaled multiplication (OP/s)
      //=====================================================================================

      // Scaled multiplication with the given vector/matrix
      {
         test_ = "Scaled multiplication with the given vector/matrix (OP/s)";

         try {
            dres_   = ( lhs_ * rhs_ ) / scalar;
            sres_   = ( lhs_ * rhs_ ) / scalar;
            refres_ = ( reflhs_ * refrhs_ ) / scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side row-major sparse matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<MT>();

         try {
            dres_ = ( lhs_ * orhs_ ) / scalar;
            sres_ = ( lhs_ * orhs_ ) / scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side column-major sparse matrix type:\n"
                << "     " << typeid( TMT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<TMT>();
      }

      // Scaled multiplication with evaluated vector/matrix
      {
         test_ = "Scaled multiplication with evaluated vector/matrix (OP/s)";

         try {
            dres_  = ( eval( lhs_ ) * eval( rhs_ ) ) / scalar;
            sres_  = ( eval( lhs_ ) * eval( rhs_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side row-major sparse matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<MT>();

         try {
            dres_ = ( eval( lhs_ ) * eval( orhs_ ) ) / scalar;
            sres_ = ( eval( lhs_ ) * eval( orhs_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side column-major sparse matrix type:\n"
                << "     " << typeid( TMT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<TMT>();
      }


      //=====================================================================================
      // Scaled multiplication with addition assignment (s*OP)
      //=====================================================================================

      // Scaled multiplication with addition assignment with the given vector/matrix
      {
         test_ = "Scaled multiplication with addition assignment with the given vector/matrix (s*OP)";

         try {
            dres_   += scalar * ( lhs_ * rhs_ );
            sres_   += scalar * ( lhs_ * rhs_ );
            refres_ += scalar * ( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed addition assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side row-major sparse matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<MT>();

         try {
            dres_   += scalar * ( lhs_ * orhs_ );
            sres_   += scalar * ( lhs_ * orhs_ );
            refres_ += scalar * ( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed addition assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side column-major sparse matrix type:\n"
                << "     " << typeid( TMT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<TMT>();
      }

      // Scaled multiplication with addition assignment with evaluated vector/matrix
      {
         test_ = "Scaled multiplication with addition assignment with evaluated vector/matrix (s*OP)";

         try {
            dres_   += scalar * ( eval( lhs_ ) * eval( rhs_ ) );
            sres_   += scalar * ( eval( lhs_ ) * eval( rhs_ ) );
            refres_ += scalar * ( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed addition assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side row-major sparse matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<MT>();

         try {
            dres_   += scalar * ( eval( lhs_ ) * eval( orhs_ ) );
            sres_   += scalar * ( eval( lhs_ ) * eval( orhs_ ) );
            refres_ += scalar * ( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed addition assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side column-major sparse matrix type:\n"
                << "     " << typeid( TMT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<TMT>();
      }


      //=====================================================================================
      // Scaled multiplication with addition assignment (OP*s)
      //=====================================================================================

      // Scaled multiplication with addition assignment with the given vector/matrix
      {
         test_ = "Scaled multiplication with addition assignment with the given vector/matrix (OP*s)";

         try {
            dres_   += ( lhs_ * rhs_ ) * scalar;
            sres_   += ( lhs_ * rhs_ ) * scalar;
            refres_ += ( reflhs_ * refrhs_ ) * scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed addition assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side row-major sparse matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<MT>();

         try {
            dres_   += ( lhs_ * orhs_ ) * scalar;
            sres_   += ( lhs_ * orhs_ ) * scalar;
            refres_ += ( reflhs_ * refrhs_ ) * scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed addition assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side column-major sparse matrix type:\n"
                << "     " << typeid( TMT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<TMT>();
      }

      // Scaled multiplication with addition assignment with evaluated vector/matrix
      {
         test_ = "Scaled multiplication with addition assignment with evaluated vector/matrix (OP*s)";

         try {
            dres_   += ( eval( lhs_ ) * eval( rhs_ ) ) * scalar;
            sres_   += ( eval( lhs_ ) * eval( rhs_ ) ) * scalar;
            refres_ += ( eval( reflhs_ ) * eval( refrhs_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed addition assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side row-major sparse matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<MT>();

         try {
            dres_   += ( eval( lhs_ ) * eval( orhs_ ) ) * scalar;
            sres_   += ( eval( lhs_ ) * eval( orhs_ ) ) * scalar;
            refres_ += ( eval( reflhs_ ) * eval( refrhs_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed addition assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side column-major sparse matrix type:\n"
                << "     " << typeid( TMT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<TMT>();
      }


      //=====================================================================================
      // Scaled multiplication with addition assignment (OP/s)
      //=====================================================================================

      // Scaled multiplication with addition assignment with the given vector/matrix
      {
         test_ = "Scaled multiplication with addition assignment with the given vector/matrix (OP/s)";

         try {
            dres_   += ( lhs_ * rhs_ ) / scalar;
            sres_   += ( lhs_ * rhs_ ) / scalar;
            refres_ += ( reflhs_ * refrhs_ ) / scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed addition assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side row-major sparse matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<MT>();

         try {
            dres_   += ( lhs_ * orhs_ ) / scalar;
            sres_   += ( lhs_ * orhs_ ) / scalar;
            refres_ += ( reflhs_ * refrhs_ ) / scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed addition assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side column-major sparse matrix type:\n"
                << "     " << typeid( TMT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<TMT>();
      }

      // Scaled multiplication with addition assignment with evaluated vector/matrix
      {
         test_ = "Scaled multiplication with addition assignment with evaluated vector/matrix (OP/s)";

         try {
            dres_   += ( eval( lhs_ ) * eval( rhs_ ) ) / scalar;
            sres_   += ( eval( lhs_ ) * eval( rhs_ ) ) / scalar;
            refres_ += ( eval( reflhs_ ) * eval( refrhs_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed addition assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side row-major sparse matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<MT>();

         try {
            dres_   += ( eval( lhs_ ) * eval( orhs_ ) ) / scalar;
            sres_   += ( eval( lhs_ ) * eval( orhs_ ) ) / scalar;
            refres_ += ( eval( reflhs_ ) * eval( refrhs_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed addition assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side column-major sparse matrix type:\n"
                << "     " << typeid( TMT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<TMT>();
      }


      //=====================================================================================
      // Scaled multiplication with subtraction assignment (s*OP)
      //=====================================================================================

      // Scaled multiplication with subtraction assignment with the given vector/matrix
      {
         test_ = "Scaled multiplication with subtraction assignment with the given vector/matrix (s*OP)";

         try {
            dres_   -= scalar * ( lhs_ * rhs_ );
            sres_   -= scalar * ( lhs_ * rhs_ );
            refres_ -= scalar * ( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side row-major sparse matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<MT>();

         try {
            dres_   -= scalar * ( lhs_ * orhs_ );
            sres_   -= scalar * ( lhs_ * orhs_ );
            refres_ -= scalar * ( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side column-major sparse matrix type:\n"
                << "     " << typeid( TMT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<TMT>();
      }

      // Scaled multiplication with subtraction assignment with evaluated vector/matrix
      {
         test_ = "Left-scaled multiplication with subtraction assignment with evaluated vector/matrix (s*OP)";

         try {
            dres_   -= scalar * ( eval( lhs_ ) * eval( rhs_ ) );
            sres_   -= scalar * ( eval( lhs_ ) * eval( rhs_ ) );
            refres_ -= scalar * ( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side row-major sparse matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<MT>();

         try {
            dres_   -= scalar * ( eval( lhs_ ) * eval( orhs_ ) );
            sres_   -= scalar * ( eval( lhs_ ) * eval( orhs_ ) );
            refres_ -= scalar * ( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side column-major sparse matrix type:\n"
                << "     " << typeid( TMT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<TMT>();
      }


      //=====================================================================================
      // Scaled multiplication with subtraction assignment (OP*s)
      //=====================================================================================

      // Scaled multiplication with subtraction assignment with the given vector/matrix
      {
         test_ = "Scaled multiplication with subtraction assignment with the given vector/matrix (OP*s)";

         try {
            dres_   -= ( lhs_ * rhs_ ) * scalar;
            sres_   -= ( lhs_ * rhs_ ) * scalar;
            refres_ -= ( reflhs_ * refrhs_ ) * scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side row-major sparse matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<MT>();

         try {
            dres_   -= ( lhs_ * orhs_ ) * scalar;
            sres_   -= ( lhs_ * orhs_ ) * scalar;
            refres_ -= ( reflhs_ * refrhs_ ) * scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side column-major sparse matrix type:\n"
                << "     " << typeid( TMT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<TMT>();
      }

      // Scaled multiplication with subtraction assignment with evaluated vector/matrix
      {
         test_ = "Scaled multiplication with subtraction assignment with evaluated vector/matrix (OP*s)";

         try {
            dres_   -= ( eval( lhs_ ) * eval( rhs_ ) ) * scalar;
            sres_   -= ( eval( lhs_ ) * eval( rhs_ ) ) * scalar;
            refres_ -= ( eval( reflhs_ ) * eval( refrhs_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side row-major sparse matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<MT>();

         try {
            dres_   -= ( eval( lhs_ ) * eval( orhs_ ) ) * scalar;
            sres_   -= ( eval( lhs_ ) * eval( orhs_ ) ) * scalar;
            refres_ -= ( eval( reflhs_ ) * eval( refrhs_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side column-major sparse matrix type:\n"
                << "     " << typeid( TMT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<TMT>();
      }


      //=====================================================================================
      // Scaled multiplication with subtraction assignment (OP/s)
      //=====================================================================================

      // Scaled multiplication with subtraction assignment with the given vector/matrix
      {
         test_ = "Scaled multiplication with subtraction assignment with the given vector/matrix (OP/s)";

         try {
            dres_   -= ( lhs_ * rhs_ ) / scalar;
            sres_   -= ( lhs_ * rhs_ ) / scalar;
            refres_ -= ( reflhs_ * refrhs_ ) / scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side row-major sparse matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<MT>();

         try {
            dres_   -= ( lhs_ * orhs_ ) / scalar;
            sres_   -= ( lhs_ * orhs_ ) / scalar;
            refres_ -= ( reflhs_ * refrhs_ ) / scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side column-major sparse matrix type:\n"
                << "     " << typeid( TMT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<TMT>();
      }

      // Scaled multiplication with subtraction assignment with evaluated vector/matrix
      {
         test_ = "Scaled multiplication with subtraction assignment with evaluated vector/matrix (OP/s)";

         try {
            dres_   -= ( eval( lhs_ ) * eval( rhs_ ) ) / scalar;
            sres_   -= ( eval( lhs_ ) * eval( rhs_ ) ) / scalar;
            refres_ -= ( eval( reflhs_ ) * eval( refrhs_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side row-major sparse matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<MT>();

         try {
            dres_   -= ( eval( lhs_ ) * eval( orhs_ ) ) / scalar;
            sres_   -= ( eval( lhs_ ) * eval( orhs_ ) ) / scalar;
            refres_ -= ( eval( reflhs_ ) * eval( refrhs_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side column-major sparse matrix type:\n"
                << "     " << typeid( TMT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<TMT>();
      }


      //=====================================================================================
      // Scaled multiplication with multiplication assignment (s*OP)
      //=====================================================================================

      // Scaled multiplication with multiplication assignment with the given vector/matrix
      {
         test_ = "Scaled multiplication with multiplication assignment with the given vector/matrix (s*OP)";

         try {
            dres_   *= scalar * ( lhs_ * rhs_ );
            sres_   *= scalar * ( lhs_ * rhs_ );
            refres_ *= scalar * ( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side row-major sparse matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<MT>();

         try {
            dres_   *= scalar * ( lhs_ * orhs_ );
            sres_   *= scalar * ( lhs_ * orhs_ );
            refres_ *= scalar * ( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side column-major sparse matrix type:\n"
                << "     " << typeid( TMT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<TMT>();
      }

      // Scaled multiplication with multiplication assignment with evaluated vector/matrix
      {
         test_ = "Scaled multiplication with multiplication assignment with evaluated vector/matrix (s*OP)";

         try {
            dres_   *= scalar * ( eval( lhs_ ) * eval( rhs_ ) );
            sres_   *= scalar * ( eval( lhs_ ) * eval( rhs_ ) );
            refres_ *= scalar * ( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side row-major sparse matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<MT>();

         try {
            dres_   *= scalar * ( eval( lhs_ ) * eval( orhs_ ) );
            sres_   *= scalar * ( eval( lhs_ ) * eval( orhs_ ) );
            refres_ *= scalar * ( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side column-major sparse matrix type:\n"
                << "     " << typeid( TMT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<TMT>();
      }


      //=====================================================================================
      // Scaled multiplication with multiplication assignment (OP*s)
      //=====================================================================================

      // Scaled multiplication with multiplication assignment with the given vector/matrix
      {
         test_ = "Scaled multiplication with multiplication assignment with the given vector/matrix (OP*s)";

         try {
            dres_   *= ( lhs_ * rhs_ ) * scalar;
            sres_   *= ( lhs_ * rhs_ ) * scalar;
            refres_ *= ( reflhs_ * refrhs_ ) * scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side row-major sparse matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<MT>();

         try {
            dres_   *= ( lhs_ * orhs_ ) * scalar;
            sres_   *= ( lhs_ * orhs_ ) * scalar;
            refres_ *= ( reflhs_ * refrhs_ ) * scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side column-major sparse matrix type:\n"
                << "     " << typeid( TMT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<TMT>();
      }

      // Scaled multiplication with multiplication assignment with evaluated vector/matrix
      {
         test_ = "Scaled multiplication with multiplication assignment with evaluated vector/matrix (OP*s)";

         try {
            dres_   *= ( eval( lhs_ ) * eval( rhs_ ) ) * scalar;
            sres_   *= ( eval( lhs_ ) * eval( rhs_ ) ) * scalar;
            refres_ *= ( eval( reflhs_ ) * eval( refrhs_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side row-major sparse matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<MT>();

         try {
            dres_   *= ( eval( lhs_ ) * eval( orhs_ ) ) * scalar;
            sres_   *= ( eval( lhs_ ) * eval( orhs_ ) ) * scalar;
            refres_ *= ( eval( reflhs_ ) * eval( refrhs_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side column-major sparse matrix type:\n"
                << "     " << typeid( TMT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<TMT>();
      }


      //=====================================================================================
      // Scaled multiplication with multiplication assignment (OP/s)
      //=====================================================================================

      // Scaled multiplication with multiplication assignment with the given vector/matrix
      {
         test_ = "Scaled multiplication with multiplication assignment with the given vector/matrix (OP/s)";

         try {
            dres_   *= ( lhs_ * rhs_ ) / scalar;
            sres_   *= ( lhs_ * rhs_ ) / scalar;
            refres_ *= ( reflhs_ * refrhs_ ) / scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side row-major sparse matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<MT>();

         try {
            dres_   *= ( lhs_ * orhs_ ) / scalar;
            sres_   *= ( lhs_ * orhs_ ) / scalar;
            refres_ *= ( reflhs_ * refrhs_ ) / scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side column-major sparse matrix type:\n"
                << "     " << typeid( TMT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<TMT>();
      }

      // Scaled multiplication with multiplication assignment with evaluated vector/matrix
      {
         test_ = "Scaled multiplication with multiplication assignment with evaluated vector/matrix (OP/s)";

         try {
            dres_   *= ( eval( lhs_ ) * eval( rhs_ ) ) / scalar;
            sres_   *= ( eval( lhs_ ) * eval( rhs_ ) ) / scalar;
            refres_ *= ( eval( reflhs_ ) * eval( refrhs_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side row-major sparse matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<MT>();

         try {
            dres_   *= ( eval( lhs_ ) * eval( orhs_ ) ) / scalar;
            sres_   *= ( eval( lhs_ ) * eval( orhs_ ) ) / scalar;
            refres_ *= ( eval( reflhs_ ) * eval( refrhs_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side column-major sparse matrix type:\n"
                << "     " << typeid( TMT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<TMT>();
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the transpose sparse vector/sparse matrix multiplication.
//
// \return void
// \exception std::runtime_error Multiplication error detected.
//
// This function tests the transpose vector/matrix multiplication with plain assignment,
// addition assignment, subtraction assignment, and multiplication assignment. In case any
// error resulting from the multiplication or the subsequent assignment is detected, a
// \a std::runtime_error exception is thrown.
*/
template< typename VT    // Type of the left-hand side sparse vector
        , typename MT >  // Type of the right-hand side sparse matrix
void TSVecSMatMult<VT,MT>::testTransposeOperation()
{
#if BLAZETEST_MATHTEST_TEST_TRANSPOSE_OPERATION
   if( BLAZETEST_MATHTEST_TEST_TRANSPOSE_OPERATION > 1 )
   {
      //=====================================================================================
      // Transpose multiplication
      //=====================================================================================

      // Transpose multiplication with the given vector/matrix
      {
         test_ = "Transpose multiplication with the given vector/matrix";

         try {
            tdres_   = trans( lhs_ * rhs_ );
            tsres_   = trans( lhs_ * rhs_ );
            trefres_ = trans( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side row-major sparse matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults<MT>();

         try {
            tdres_  = trans( lhs_ * orhs_ );
            tsres_  = trans( lhs_ * orhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side column-major sparse matrix type:\n"
                << "     " << typeid( TMT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults<TMT>();
      }

      // Transpose multiplication with evaluated vector/matrix
      {
         test_ = "Transpose multiplication with evaluated vector/matrix";

         try {
            tdres_ = trans( eval( lhs_ ) * eval( rhs_ ) );
            tsres_ = trans( eval( lhs_ ) * eval( rhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side row-major sparse matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults<MT>();

         try {
            tdres_ = trans( eval( lhs_ ) * eval( orhs_ ) );
            tsres_ = trans( eval( lhs_ ) * eval( orhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side column-major sparse matrix type:\n"
                << "     " << typeid( TMT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults<TMT>();
      }


      //=====================================================================================
      // Transpose multiplication with addition assignment
      //=====================================================================================

      // Transpose multiplication with addition assignment with the given vector/matrix
      {
         test_ = "Transpose multiplication with addition assignment with the given vector/matrix";

         try {
            tdres_   += trans( lhs_ * rhs_ );
            tsres_   += trans( lhs_ * rhs_ );
            trefres_ += trans( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed addition assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side row-major sparse matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults<MT>();

         try {
            tdres_   += trans( lhs_ * orhs_ );
            tsres_   += trans( lhs_ * orhs_ );
            trefres_ += trans( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed addition assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side column-major sparse matrix type:\n"
                << "     " << typeid( TMT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults<TMT>();
      }

      // Transpose multiplication with addition assignment with evaluated vector/matrix
      {
         test_ = "Transpose multiplication with addition assignment with evaluated vector/matrix";

         try {
            tdres_   += trans( eval( lhs_ ) * eval( rhs_ ) );
            tsres_   += trans( eval( lhs_ ) * eval( rhs_ ) );
            trefres_ += trans( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed addition assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side row-major sparse matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults<MT>();

         try {
            tdres_   += trans( eval( lhs_ ) * eval( orhs_ ) );
            tsres_   += trans( eval( lhs_ ) * eval( orhs_ ) );
            trefres_ += trans( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed addition assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side column-major sparse matrix type:\n"
                << "     " << typeid( TMT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults<TMT>();
      }


      //=====================================================================================
      // Transpose multiplication with subtraction assignment
      //=====================================================================================

      // Transpose multiplication with subtraction assignment with the given vector/matrix
      {
         test_ = "Transpose multiplication with subtraction assignment with the given vector/matrix";

         try {
            tdres_   -= trans( lhs_ * rhs_ );
            tsres_   -= trans( lhs_ * rhs_ );
            trefres_ -= trans( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side row-major sparse matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults<MT>();

         try {
            tdres_   -= trans( lhs_ * orhs_ );
            tsres_   -= trans( lhs_ * orhs_ );
            trefres_ -= trans( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side column-major sparse matrix type:\n"
                << "     " << typeid( TMT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults<TMT>();
      }

      // Transpose multiplication with subtraction assignment with evaluated vector/matrix
      {
         test_ = "Transpose multiplication with subtraction assignment with evaluated vector/matrix";

         try {
            tdres_   -= trans( eval( lhs_ ) * eval( rhs_ ) );
            tsres_   -= trans( eval( lhs_ ) * eval( rhs_ ) );
            trefres_ -= trans( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side row-major sparse matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults<MT>();

         try {
            tdres_   -= trans( eval( lhs_ ) * eval( orhs_ ) );
            tsres_   -= trans( eval( lhs_ ) * eval( orhs_ ) );
            trefres_ -= trans( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side column-major sparse matrix type:\n"
                << "     " << typeid( TMT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults<TMT>();
      }


      //=====================================================================================
      // Transpose multiplication with multiplication assignment
      //=====================================================================================

      // Transpose multiplication with multiplication assignment with the given vector/matrix
      {
         test_ = "Transpose multiplication with multiplication assignment with the given vector/matrix";

         try {
            tdres_   *= trans( lhs_ * rhs_ );
            tsres_   *= trans( lhs_ * rhs_ );
            trefres_ *= trans( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side row-major sparse matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults<MT>();

         try {
            tdres_   *= trans( lhs_ * orhs_ );
            tsres_   *= trans( lhs_ * orhs_ );
            trefres_ *= trans( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side column-major sparse matrix type:\n"
                << "     " << typeid( TMT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults<TMT>();
      }

      // Transpose multiplication with multiplication assignment with evaluated vector/matrix
      {
         test_ = "Transpose multiplication with multiplication assignment with evaluated vector/matrix";

         try {
            tdres_   *= trans( eval( lhs_ ) * eval( rhs_ ) );
            tsres_   *= trans( eval( lhs_ ) * eval( rhs_ ) );
            trefres_ *= trans( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side row-major sparse matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults<MT>();

         try {
            tdres_   *= trans( eval( lhs_ ) * eval( orhs_ ) );
            tsres_   *= trans( eval( lhs_ ) * eval( orhs_ ) );
            trefres_ *= trans( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side column-major sparse matrix type:\n"
                << "     " << typeid( TMT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults<TMT>();
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the abs sparse vector/sparse matrix multiplication.
//
// \return void
// \exception std::runtime_error Multiplication error detected.
//
// This function tests the abs vector/matrix multiplication with plain assignment, addition
// assignment, subtraction assignment, and multiplication assignment. In case any error
// resulting from the multiplication or the subsequent assignment is detected, a
// \a std::runtime_error exception is thrown.
*/
template< typename VT    // Type of the left-hand side sparse vector
        , typename MT >  // Type of the right-hand side sparse matrix
void TSVecSMatMult<VT,MT>::testAbsOperation()
{
#if BLAZETEST_MATHTEST_TEST_ABS_OPERATION
   if( BLAZETEST_MATHTEST_TEST_ABS_OPERATION > 1 )
   {
      //=====================================================================================
      // Abs multiplication
      //=====================================================================================

      // Abs multiplication with the given vector/matrix
      {
         test_ = "Abs multiplication with the given vector/matrix";

         try {
            dres_   = abs( lhs_ * rhs_ );
            sres_   = abs( lhs_ * rhs_ );
            refres_ = abs( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side row-major sparse matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<MT>();

         try {
            dres_ = abs( lhs_ * orhs_ );
            sres_ = abs( lhs_ * orhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side column-major sparse matrix type:\n"
                << "     " << typeid( TMT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<TMT>();
      }

      // Abs multiplication with evaluated vector/matrix
      {
         test_ = "Abs multiplication with evaluated vector/matrix";

         try {
            dres_  = abs( eval( lhs_ ) * eval( rhs_ ) );
            sres_  = abs( eval( lhs_ ) * eval( rhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side row-major sparse matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<MT>();

         try {
            dres_ = abs( eval( lhs_ ) * eval( orhs_ ) );
            sres_ = abs( eval( lhs_ ) * eval( orhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side column-major sparse matrix type:\n"
                << "     " << typeid( TMT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<TMT>();
      }


      //=====================================================================================
      // Abs multiplication with addition assignment
      //=====================================================================================

      // Abs multiplication with addition assignment with the given vector/matrix
      {
         test_ = "Abs multiplication with addition assignment with the given vector/matrix";

         try {
            dres_   += abs( lhs_ * rhs_ );
            sres_   += abs( lhs_ * rhs_ );
            refres_ += abs( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed addition assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side row-major sparse matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<MT>();

         try {
            dres_   += abs( lhs_ * orhs_ );
            sres_   += abs( lhs_ * orhs_ );
            refres_ += abs( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed addition assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side column-major sparse matrix type:\n"
                << "     " << typeid( TMT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<TMT>();
      }

      // Abs multiplication with addition assignment with evaluated vector/matrix
      {
         test_ = "Abs multiplication with addition assignment with evaluated vector/matrix";

         try {
            dres_   += abs( eval( lhs_ ) * eval( rhs_ ) );
            sres_   += abs( eval( lhs_ ) * eval( rhs_ ) );
            refres_ += abs( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed addition assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side row-major sparse matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<MT>();

         try {
            dres_   += abs( eval( lhs_ ) * eval( orhs_ ) );
            sres_   += abs( eval( lhs_ ) * eval( orhs_ ) );
            refres_ += abs( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed addition assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side column-major sparse matrix type:\n"
                << "     " << typeid( TMT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<TMT>();
      }


      //=====================================================================================
      // Abs multiplication with subtraction assignment
      //=====================================================================================

      // Abs multiplication with subtraction assignment with the given vector/matrix
      {
         test_ = "Abs multiplication with subtraction assignment with the given vector/matrix";

         try {
            dres_   -= abs( lhs_ * rhs_ );
            sres_   -= abs( lhs_ * rhs_ );
            refres_ -= abs( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side row-major sparse matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<MT>();

         try {
            dres_   -= abs( lhs_ * orhs_ );
            sres_   -= abs( lhs_ * orhs_ );
            refres_ -= abs( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side column-major sparse matrix type:\n"
                << "     " << typeid( TMT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<TMT>();
      }

      // Abs multiplication with subtraction assignment with evaluated vector/matrix
      {
         test_ = "Abs multiplication with subtraction assignment with evaluated vector/matrix";

         try {
            dres_   -= abs( eval( lhs_ ) * eval( rhs_ ) );
            sres_   -= abs( eval( lhs_ ) * eval( rhs_ ) );
            refres_ -= abs( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side row-major sparse matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<MT>();

         try {
            dres_   -= abs( eval( lhs_ ) * eval( orhs_ ) );
            sres_   -= abs( eval( lhs_ ) * eval( orhs_ ) );
            refres_ -= abs( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side column-major sparse matrix type:\n"
                << "     " << typeid( TMT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<TMT>();
      }


      //=====================================================================================
      // Abs multiplication with multiplication assignment
      //=====================================================================================

      // Abs multiplication with multiplication assignment with the given vector/matrix
      {
         test_ = "Abs multiplication with multiplication assignment with the given vector/matrix";

         try {
            dres_   *= abs( lhs_ * rhs_ );
            sres_   *= abs( lhs_ * rhs_ );
            refres_ *= abs( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side row-major sparse matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<MT>();

         try {
            dres_   *= abs( lhs_ * orhs_ );
            sres_   *= abs( lhs_ * orhs_ );
            refres_ *= abs( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side column-major sparse matrix type:\n"
                << "     " << typeid( TMT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<TMT>();
      }

      // Abs multiplication with multiplication assignment with evaluated vector/matrix
      {
         test_ = "Abs multiplication with multiplication assignment with evaluated vector/matrix";

         try {
            dres_   *= abs( eval( lhs_ ) * eval( rhs_ ) );
            sres_   *= abs( eval( lhs_ ) * eval( rhs_ ) );
            refres_ *= abs( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side row-major sparse matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<MT>();

         try {
            dres_   *= abs( eval( lhs_ ) * eval( orhs_ ) );
            sres_   *= abs( eval( lhs_ ) * eval( orhs_ ) );
            refres_ *= abs( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication assignment operation\n"
                << " Details:\n"
                << "   Left-hand side transpose sparse vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Right-hand side column-major sparse matrix type:\n"
                << "     " << typeid( TMT ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<TMT>();
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
// The template argument \a LT indicates the types of the left-hand side operand used for
// the computations.
*/
template< typename VT    // Type of the left-hand side sparse vector
        , typename MT >  // Type of the right-hand side sparse matrix
template< typename RT >  // Type of the right-hand side operand
void TSVecSMatMult<VT,MT>::checkResults()
{
   using blaze::IsRowMajorMatrix;

   if( !isEqual( dres_, refres_ ) ) {
      std::ostringstream oss;
      oss.precision( 20 );
      oss << " Test : " << test_ << "\n"
          << " Error: Incorrect dense result detected\n"
          << " Details:\n"
          << "   Left-hand side transpose sparse vector type:\n"
          << "     " << typeid( VT ).name() << "\n"
          << "   Right-hand side " << ( IsRowMajorMatrix<RT>::value ? ( "row-major" ) : ( "column-major" ) ) << " sparse matrix type:\n"
          << "     " << typeid( RT ).name() << "\n"
          << "   Result:\n" << dres_ << "\n"
          << "   Expected result:\n" << refres_ << "\n";
      throw std::runtime_error( oss.str() );
   }

   if( !isEqual( sres_, refres_ ) ) {
      std::ostringstream oss;
      oss.precision( 20 );
      oss << " Test : " << test_ << "\n"
          << " Error: Incorrect sparse result detected\n"
          << " Details:\n"
          << "   Left-hand side transpose sparse vector type:\n"
          << "     " << typeid( VT ).name() << "\n"
          << "   Right-hand side " << ( IsRowMajorMatrix<RT>::value ? ( "row-major" ) : ( "column-major" ) ) << " sparse matrix type:\n"
          << "     " << typeid( RT ).name() << "\n"
          << "   Result:\n" << sres_ << "\n"
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
// results. The template argument \a LT indicates the types of the left-hand side operand
// used for the computations.
*/
template< typename VT    // Type of the left-hand side sparse vector
        , typename MT >  // Type of the right-hand side sparse matrix
template< typename RT >  // Type of the right-hand side operand
void TSVecSMatMult<VT,MT>::checkTransposeResults()
{
   using blaze::IsRowMajorMatrix;

   if( !isEqual( tdres_, trefres_ ) ) {
      std::ostringstream oss;
      oss.precision( 20 );
      oss << " Test : " << test_ << "\n"
          << " Error: Incorrect dense result detected\n"
          << " Details:\n"
          << "   Left-hand side transpose sparse vector type:\n"
          << "     " << typeid( VT ).name() << "\n"
          << "   Right-hand side " << ( IsRowMajorMatrix<RT>::value ? ( "row-major" ) : ( "column-major" ) ) << " sparse matrix type:\n"
          << "     " << typeid( RT ).name() << "\n"
          << "   Transpose result:\n" << tdres_ << "\n"
          << "   Expected transpose result:\n" << trefres_ << "\n";
      throw std::runtime_error( oss.str() );
   }

   if( !isEqual( tsres_, trefres_ ) ) {
      std::ostringstream oss;
      oss.precision( 20 );
      oss << " Test : " << test_ << "\n"
          << " Error: Incorrect sparse result detected\n"
          << " Details:\n"
          << "   Left-hand side transpose sparse vector type:\n"
          << "     " << typeid( VT ).name() << "\n"
          << "   Right-hand side " << ( IsRowMajorMatrix<RT>::value ? ( "row-major" ) : ( "column-major" ) ) << " sparse matrix type:\n"
          << "     " << typeid( RT ).name() << "\n"
          << "   Transpose result:\n" << tsres_ << "\n"
          << "   Expected transpose result:\n" << trefres_ << "\n";
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
/*!\brief Testing the vector/matrix multiplication between two specific types.
//
// \param creator1 The creator for the left-hand side vector.
// \param creator2 The creator for the right-hand side matrix.
// \return void
*/
template< typename VT    // Type of the left-hand side sparse vector
        , typename MT >  // Type of the right-hand side sparse matrix
void runTest( const Creator<VT>& creator1, const Creator<MT>& creator2 )
{
   for( size_t rep=0; rep<repetitions; ++rep ) {
      TSVecSMatMult<VT,MT>( creator1, creator2 );
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
/*!\brief Macro for the execution of a sparse vector/sparse matrix multiplication test case.
*/
#define RUN_TSVECSMATMULT_TEST( C1, C2 ) \
   blazetest::mathtest::tsvecsmatmult::runTest( C1, C2 )
/*! \endcond */
//*************************************************************************************************

} // namespace tsvecsmatmult

} // namespace mathtest

} // namespace blazetest

#endif
