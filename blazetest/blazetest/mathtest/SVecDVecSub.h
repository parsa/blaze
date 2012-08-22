//=================================================================================================
/*!
//  \file blazetest/mathtest/SVecDVecSub.h
//  \brief Header file for the sparse vector/dense vector subtraction math test
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

#ifndef _BLAZETEST_MATHTEST_SVECDVECSUB_H_
#define _BLAZETEST_MATHTEST_SVECDVECSUB_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <sstream>
#include <stdexcept>
#include <string>
#include <typeinfo>
#include <blaze/math/CompressedVector.h>
#include <blaze/math/constraints/DenseVector.h>
#include <blaze/math/constraints/SparseVector.h>
#include <blaze/math/constraints/TransposeFlag.h>
#include <blaze/math/DynamicVector.h>
#include <blaze/math/MathTrait.h>
#include <blaze/math/shims/Equal.h>
#include <blaze/math/StaticVector.h>
#include <blaze/math/typetraits/IsTransposeVector.h>
#include <blaze/util/constraints/Numeric.h>
#include <blaze/util/constraints/SameType.h>
#include <blazetest/system/MathTest.h>
#include <blazetest/util/Creator.h>
#include <blazetest/util/Utility.h>


namespace blazetest {

namespace mathtest {

namespace svecdvecsub {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Auxiliary class template for the sparse vector/dense vector subtraction math test.
//
// The SVecDVecSub class template represents one particular vector subtraction test between two
// vectors of a particular type. The two template arguments \a VT1 and \a VT2 represent the
// types of the left-hand side and right-hand side vector, respectively.
*/
template< typename VT1    // Type of the left-hand side sparse vector
        , typename VT2 >  // Type of the right-hand side dense vector
class SVecDVecSub
{
 private:
   //**Type definitions****************************************************************************
   typedef typename VT1::TransposeType                    TVT1;  //!< Transpose vector type 1
   typedef typename VT2::TransposeType                    TVT2;  //!< Transpose vector type 2
   typedef typename blaze::MathTrait<VT1,VT2>::SubType    RE;    //!< Default result type
   typedef typename blaze::MathTrait<TVT1,TVT2>::SubType  TRE;   //!< Transpose default result type
   //**********************************************************************************************

   //**Enumerations********************************************************************************
   enum { TF = blaze::IsTransposeVector<VT1>::value };
   //**********************************************************************************************

   //**Type definitions****************************************************************************
   typedef typename VT1::ElementType        ET1;    //!< Element type 1
   typedef typename VT2::ElementType        ET2;    //!< Element type 2
   typedef typename RE::ElementType         RET;    //!< Resulting element type
   typedef blaze::DynamicVector<ET1,TF>     RT1;    //!< Reference type 1
   typedef blaze::DynamicVector<ET2,TF>     RT2;    //!< Reference type 2
   typedef typename RT1::TransposeType      TRT1;   //!< Transpose reference type 1
   typedef typename RT2::TransposeType      TRT2;   //!< Transpose reference type 2
   typedef blaze::DynamicVector<RET,TF>     DRRE;   //!< Dense reference result type
   typedef blaze::CompressedVector<RET,TF>  SRRE;   //!< Sparse reference result type
   typedef typename DRRE::TransposeType     TDRRE;  //!< Transpose dense reference result type
   typedef typename SRRE::TransposeType     TSRRE;  //!< Transpose sparse reference result type
   typedef RE                               DRE;    //!< Dense result type
   typedef SRRE                             SRE;    //!< Sparse result type
   typedef TRE                              TDRE;   //!< Transpose dense result type
   typedef TSRRE                            TSRE;   //!< Transpose sparse result type
   //**********************************************************************************************

 public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit SVecDVecSub( const Creator<VT1>& creator1, const Creator<VT2>& creator2 );
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
   template< typename LT, typename RT > void checkResults();
   template< typename LT, typename RT > void checkTransposeResults();
   //@}
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   VT1   lhs_;      //!< The left-hand side sparse vector.
   VT2   rhs_;      //!< The right-hand side dense vector.
   DRE   dres_;     //!< The dense vector for the result of the vector subtraction.
   SRE   sres_;     //!< The sparse vector for the result of the vector subtraction.
   RT1   reflhs_;   //!< The reference left-hand side vector.
   RT2   refrhs_;   //!< The reference right-hand side vector.
   DRRE  refres_;   //!< The reference result.
   TVT1  tlhs_;     //!< The transpose left-hand side vector.
   TVT2  trhs_;     //!< The transpose right-hand side vector.
   TDRE  tdres_;    //!< The dense vector for the result of the transpose vector subtraction.
   TSRE  tsres_;    //!< The sparse vector for the result of the transpose vector subtraction.
   TRT1  treflhs_;  //!< The reference left-hand side transpose vector.
   TRT2  trefrhs_;  //!< The reference right-hand side transpose vector.
   TDRRE trefres_;  //!< The transpose reference result.

   std::string test_;  //!< Label of the currently performed test.
   //@}
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( VT1   );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( VT2   );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( TVT1  );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( TVT2  );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( RT1   );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( RT2   );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( TRT1  );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( TRT2  );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( DRE   );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( SRE   );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( DRRE  );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( SRRE  );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( TDRE  );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( TSRE  );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( TDRRE );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( TSRRE );
   BLAZE_CONSTRAINT_VECTORS_MUST_HAVE_SAME_TRANSPOSE_FLAG( VT1 , VT2   );
   BLAZE_CONSTRAINT_VECTORS_MUST_HAVE_SAME_TRANSPOSE_FLAG( VT1 , RT1   );
   BLAZE_CONSTRAINT_VECTORS_MUST_HAVE_SAME_TRANSPOSE_FLAG( TVT1, TVT2  );
   BLAZE_CONSTRAINT_VECTORS_MUST_HAVE_SAME_TRANSPOSE_FLAG( TVT1, TRT1  );
   BLAZE_CONSTRAINT_VECTORS_MUST_HAVE_SAME_TRANSPOSE_FLAG( RT1 , RT2   );
   BLAZE_CONSTRAINT_VECTORS_MUST_HAVE_SAME_TRANSPOSE_FLAG( TRT1, TRT2  );
   BLAZE_CONSTRAINT_VECTORS_MUST_HAVE_SAME_TRANSPOSE_FLAG( VT1 , DRE   );
   BLAZE_CONSTRAINT_VECTORS_MUST_HAVE_SAME_TRANSPOSE_FLAG( VT1 , SRE   );
   BLAZE_CONSTRAINT_VECTORS_MUST_HAVE_SAME_TRANSPOSE_FLAG( RT1 , DRRE  );
   BLAZE_CONSTRAINT_VECTORS_MUST_HAVE_SAME_TRANSPOSE_FLAG( RT1 , SRRE  );
   BLAZE_CONSTRAINT_VECTORS_MUST_HAVE_SAME_TRANSPOSE_FLAG( TVT1, TDRE  );
   BLAZE_CONSTRAINT_VECTORS_MUST_HAVE_SAME_TRANSPOSE_FLAG( TVT1, TSRE  );
   BLAZE_CONSTRAINT_VECTORS_MUST_HAVE_SAME_TRANSPOSE_FLAG( TRT1, TDRRE );
   BLAZE_CONSTRAINT_VECTORS_MUST_HAVE_SAME_TRANSPOSE_FLAG( TRT1, TSRRE );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ET1, typename TVT1::ElementType   );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ET2, typename TVT2::ElementType   );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( VT1, typename TVT1::TransposeType );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( VT2, typename TVT2::TransposeType );
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
/*!\brief Constructor for the SVecDVecSub class template.
//
// \param creator1 The creator for the left-hand side sparse vector of the vector subtraction.
// \param creator2 The creator for the right-hand side dense vector of the vector subtraction.
// \exception std::runtime_error Operation error detected.
*/
template< typename VT1    // Type of the left-hand side sparse vector
        , typename VT2 >  // Type of the right-hand side dense vector
SVecDVecSub<VT1,VT2>::SVecDVecSub( const Creator<VT1>& creator1, const Creator<VT2>& creator2 )
   : lhs_( creator1() )    // The left-hand side sparse vector
   , rhs_( creator2() )    // The right-hand side dense vector
   , dres_()               // The dense vector for the result of the vector subtraction
   , sres_()               // The sparse vector for the result of the vector subtraction
   , reflhs_( lhs_ )       // The reference left-hand side vector
   , refrhs_( rhs_ )       // The reference right-hand side vector
   , refres_()             // The reference result
   , tlhs_( trans(lhs_) )  // The transpose left-hand side vector
   , trhs_( trans(rhs_) )  // The transpose right-hand side vector
   , tdres_()              // The dense vector for the result of the transpose vector subtraction
   , tsres_()              // The sparse vector for the result of the transpose vector subtraction
   , treflhs_( tlhs_ )     // The reference left-hand side transpose vector
   , trefrhs_( trhs_ )     // The reference right-hand side transpose vector
   , trefres_()            // The transpose reference result
   , test_()               // Label of the currently performed test.
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
void SVecDVecSub<VT1,VT2>::testInitialStatus()
{
   //=====================================================================================
   // Performing initial tests with the given vectors
   //=====================================================================================

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
          << "   Dense vector type:\n"
          << "     " << typeid( VT2 ).name() << "\n"
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
          << "   Dense vector type:\n"
          << "     " << typeid( VT2 ).name() << "\n"
          << "   Current initialization:\n" << rhs_ << "\n"
          << "   Expected initialization:\n" << refrhs_ << "\n";
      throw std::runtime_error( oss.str() );
   }


   //=====================================================================================
   // Performing initial tests with the transpose types
   //=====================================================================================

   // Checking the size of the left-hand side operand
   if( tlhs_.size() != treflhs_.size() ) {
      std::ostringstream oss;
      oss << " Test: Initial size comparison of transpose left-hand side sparse operand\n"
          << " Error: Invalid vector size\n"
          << " Details:\n"
          << "   Transpose sparse vector type:\n"
          << "     " << typeid( TVT1 ).name() << "\n"
          << "   Detected size = " << tlhs_.size() << "\n"
          << "   Expected size = " << treflhs_.size() << "\n";
      throw std::runtime_error( oss.str() );
   }

   // Checking the size of the right-hand side operand
   if( trhs_.size() != trefrhs_.size() ) {
      std::ostringstream oss;
      oss << " Test: Initial size comparison of transpose right-hand side dense operand\n"
          << " Error: Invalid vector size\n"
          << " Details:\n"
          << "   Transpose dense vector type:\n"
          << "     " << typeid( TVT2 ).name() << "\n"
          << "   Detected size = " << trhs_.size() << "\n"
          << "   Expected size = " << trefrhs_.size() << "\n";
      throw std::runtime_error( oss.str() );
   }

   // Checking the initialization of the left-hand side operand
   if( !isEqual( tlhs_, treflhs_ ) ) {
      std::ostringstream oss;
      oss << " Test: Initial test of initialization of transpose left-hand side sparse operand\n"
          << " Error: Invalid vector initialization\n"
          << " Details:\n"
          << "   Transpose sparse vector type:\n"
          << "     " << typeid( TVT1 ).name() << "\n"
          << "   Current initialization:\n" << tlhs_ << "\n"
          << "   Expected initialization:\n" << treflhs_ << "\n";
      throw std::runtime_error( oss.str() );
   }

   // Checking the initialization of the right-hand side operand
   if( !isEqual( trhs_, trefrhs_ ) ) {
      std::ostringstream oss;
      oss << " Test: Initial test of initialization of transpose right-hand side dense operand\n"
          << " Error: Invalid vector initialization\n"
          << " Details:\n"
          << "   Transpose dense vector type:\n"
          << "     " << typeid( TVT2 ).name() << "\n"
          << "   Current initialization:\n" << trhs_ << "\n"
          << "   Expected initialization:\n" << trefrhs_ << "\n";
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
template< typename VT1    // Type of the left-hand side sparse vector
        , typename VT2 >  // Type of the right-hand side dense vector
void SVecDVecSub<VT1,VT2>::testAssignment()
{
   //=====================================================================================
   // Performing an assignment with the given vectors
   //=====================================================================================

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
          << "   Right-hand side dense vector type:\n"
          << "     " << typeid( VT2 ).name() << "\n"
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
          << "   Dense vector type:\n"
          << "     " << typeid( VT2 ).name() << "\n"
          << "   Current initialization:\n" << rhs_ << "\n"
          << "   Expected initialization:\n" << refrhs_ << "\n";
      throw std::runtime_error( oss.str() );
   }


   //=====================================================================================
   // Performing an assignment with the transpose types
   //=====================================================================================

   try {
      tlhs_ = treflhs_;
      trhs_ = trefrhs_;
   }
   catch( std::exception& ex ) {
      std::ostringstream oss;
      oss << " Test: Assignment with the transpose types\n"
          << " Error: Failed assignment\n"
          << " Details:\n"
          << "   Transpose left-hand side sparse vector type:\n"
          << "     " << typeid( TVT1 ).name() << "\n"
          << "   Transpose right-hand side dense vector type:\n"
          << "     " << typeid( TVT2 ).name() << "\n"
          << "   Error message: " << ex.what() << "\n";
      throw std::runtime_error( oss.str() );
   }

   if( !isEqual( tlhs_, treflhs_ ) ) {
      std::ostringstream oss;
      oss << " Test: Checking the assignment result of transpose left-hand side sparse operand\n"
          << " Error: Invalid vector initialization\n"
          << " Details:\n"
          << "   Transpose sparse vector type:\n"
          << "     " << typeid( TVT1 ).name() << "\n"
          << "   Current initialization:\n" << tlhs_ << "\n"
          << "   Expected initialization:\n" << treflhs_ << "\n";
      throw std::runtime_error( oss.str() );
   }

   if( !isEqual( trhs_, trefrhs_ ) ) {
      std::ostringstream oss;
      oss << " Test: Checking the assignment result of transpose right-hand side dense operand\n"
          << " Error: Invalid vector initialization\n"
          << " Details:\n"
          << "   Transpose dense vector type:\n"
          << "     " << typeid( TVT2 ).name() << "\n"
          << "   Current initialization:\n" << trhs_ << "\n"
          << "   Expected initialization:\n" << trefrhs_ << "\n";
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
template< typename VT1    // Type of the left-hand side sparse vector
        , typename VT2 >  // Type of the right-hand side dense vector
void SVecDVecSub<VT1,VT2>::testElementAccess()
{
   using blaze::equal;


   //=====================================================================================
   // Testing the element access with the given vectors
   //=====================================================================================

   if( lhs_.size() > 0UL && rhs_.size() > 0UL )
   {
      if( !equal( ( lhs_ - rhs_ )[0UL], ( reflhs_ - refrhs_ )[0UL] ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of subtraction expression\n"
             << " Error: Unequal resulting elements at index 0 detected\n"
             << " Details:\n"
             << "   Left-hand side sparse vector type:\n"
             << "     " << typeid( VT1 ).name() << "\n"
             << "   Right-hand side dense vector type:\n"
             << "     " << typeid( VT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( ( lhs_ - eval( rhs_ ) )[0UL], ( reflhs_ - eval( refrhs_ ) )[0UL] ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of right evaluated subtraction expression\n"
             << " Error: Unequal resulting elements at index 0 detected\n"
             << " Details:\n"
             << "   Left-hand side sparse vector type:\n"
             << "     " << typeid( VT1 ).name() << "\n"
             << "   Right-hand side dense vector type:\n"
             << "     " << typeid( VT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( ( eval( lhs_ ) - rhs_ )[0UL], ( eval( reflhs_ ) - refrhs_ )[0UL] ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of left evaluated subtraction expression\n"
             << " Error: Unequal resulting elements at index 0 detected\n"
             << " Details:\n"
             << "   Left-hand side sparse vector type:\n"
             << "     " << typeid( VT1 ).name() << "\n"
             << "   Right-hand side dense vector type:\n"
             << "     " << typeid( VT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( ( eval( lhs_ ) - eval( rhs_ ) )[0UL], ( eval( reflhs_ ) - eval( refrhs_ ) )[0UL] ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of fully evaluated subtraction expression\n"
             << " Error: Unequal resulting elements at index 0 detected\n"
             << " Details:\n"
             << "   Left-hand side sparse vector type:\n"
             << "     " << typeid( VT1 ).name() << "\n"
             << "   Right-hand side dense vector type:\n"
             << "     " << typeid( VT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Testing the element access with the transpose types
   //=====================================================================================

   if( tlhs_.size() > 0UL && trhs_.size() > 0UL )
   {
      if( !equal( ( tlhs_ - trhs_ )[0UL], ( treflhs_ - trefrhs_ )[0UL] ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of transpose subtraction expression\n"
             << " Error: Unequal resulting elements at index 0 detected\n"
             << " Details:\n"
             << "   Transpose left-hand side sparse vector type:\n"
             << "     " << typeid( TVT1 ).name() << "\n"
             << "   Transpose right-hand side dense vector type:\n"
             << "     " << typeid( TVT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( ( tlhs_ - eval( trhs_ ) )[0UL], ( treflhs_ - eval( trefrhs_ ) )[0UL] ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of right evaluated transpose subtraction expression\n"
             << " Error: Unequal resulting elements at index 0 detected\n"
             << " Details:\n"
             << "   Transpose left-hand side sparse vector type:\n"
             << "     " << typeid( TVT1 ).name() << "\n"
             << "   Transpose right-hand side dense vector type:\n"
             << "     " << typeid( TVT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( ( eval( tlhs_ ) - trhs_ )[0UL], ( eval( treflhs_ ) - trefrhs_ )[0UL] ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of left evaluated transpose subtraction expression\n"
             << " Error: Unequal resulting elements at index 0 detected\n"
             << " Details:\n"
             << "   Transpose left-hand side sparse vector type:\n"
             << "     " << typeid( TVT1 ).name() << "\n"
             << "   Transpose right-hand side dense vector type:\n"
             << "     " << typeid( TVT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( ( eval( tlhs_ ) - eval( trhs_ ) )[0UL], ( eval( treflhs_ ) - eval( trefrhs_ ) )[0UL] ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of fully evaluated transpose subtraction expression\n"
             << " Error: Unequal resulting elements at index 0 detected\n"
             << " Details:\n"
             << "   Transpose left-hand side sparse vector type:\n"
             << "     " << typeid( TVT1 ).name() << "\n"
             << "   Transpose right-hand side dense vector type:\n"
             << "     " << typeid( TVT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the plain sparse vector/dense vector subtraction.
//
// \return void
// \exception std::runtime_error Subtraction error detected.
//
// This function tests the plain vector subtraction with plain assignment, addition assignment,
// subtraction assignment, and multiplication assignment. In case any error resulting from the
// subtraction or the subsequent assignment is detected, a \a std::runtime_error exception is
// thrown.
*/
template< typename VT1    // Type of the left-hand side sparse vector
        , typename VT2 >  // Type of the right-hand side dense vector
void SVecDVecSub<VT1,VT2>::testBasicOperation()
{
#if BLAZETEST_MATHTEST_TEST_BASIC_OPERATION
   if( BLAZETEST_MATHTEST_TEST_BASIC_OPERATION > 1 )
   {
      //=====================================================================================
      // Subtraction with the given vectors
      //=====================================================================================

      // Subtraction with the given vectors
      {
         test_ = " Subtraction with the given vectors";

         try {
            dres_   = lhs_ - rhs_;
            sres_   = lhs_ - rhs_;
            refres_ = reflhs_ - refrhs_;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side dense vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<VT1,VT2>();

         try {
            tdres_   = tlhs_ - trhs_;
            tsres_   = tlhs_ - trhs_;
            trefres_ = treflhs_ - trefrhs_;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction operation\n"
                << " Details:\n"
                << "   Transpose left-hand side sparse vector type:\n"
                << "     " << typeid( TVT1 ).name() << "\n"
                << "   Transpose right-hand side dense vector type:\n"
                << "     " << typeid( TVT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Subtraction with evaluated vectors
      {
         test_ = "Subtraction with evaluated vectors";

         try {
            dres_ = eval( lhs_ ) - eval( rhs_ );
            sres_ = eval( lhs_ ) - eval( rhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side dense vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<VT1,VT2>();

         try {
            tdres_ = eval( tlhs_ ) - eval( trhs_ );
            tsres_ = eval( tlhs_ ) - eval( trhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction operation\n"
                << " Details:\n"
                << "   Transpose left-hand side sparse vector type:\n"
                << "     " << typeid( TVT1 ).name() << "\n"
                << "   Transpose right-hand side dense vector type:\n"
                << "     " << typeid( TVT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Subtraction with addition assignment
      //=====================================================================================

      // Subtraction with addition assignment with the given vectors
      {
         test_ = "Subtraction with addition assignment with the given vectors";

         try {
            dres_   += lhs_ - rhs_;
            sres_   += lhs_ - rhs_;
            refres_ += reflhs_ - refrhs_;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed addition assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side dense vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<VT1,VT2>();

         try {
            tdres_   += tlhs_ - trhs_;
            tsres_   += tlhs_ - trhs_;
            trefres_ += treflhs_ - trefrhs_;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed addition assignment operation\n"
                << " Details:\n"
                << "   Transpose left-hand side sparse vector type:\n"
                << "     " << typeid( TVT1 ).name() << "\n"
                << "   Transpose right-hand side dense vector type:\n"
                << "     " << typeid( TVT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Subtraction with addition assignment with the given vectors
      {
         test_ = "Subtraction with addition assignment with evaluated vectors";

         try {
            dres_   += eval( lhs_ ) - eval( rhs_ );
            sres_   += eval( lhs_ ) - eval( rhs_ );
            refres_ += eval( reflhs_ ) - eval( refrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed addition assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side dense vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<VT1,VT2>();

         try {
            tdres_   += eval( tlhs_ ) - eval( trhs_ );
            tsres_   += eval( tlhs_ ) - eval( trhs_ );
            trefres_ += eval( treflhs_ ) - eval( trefrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed addition assignment operation\n"
                << " Details:\n"
                << "   Transpose left-hand side sparse vector type:\n"
                << "     " << typeid( TVT1 ).name() << "\n"
                << "   Transpose right-hand side dense vector type:\n"
                << "     " << typeid( TVT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Subtraction with subtraction assignment
      //=====================================================================================

      // Subtraction with subtraction assignment with the given vectors
      {
         test_ = "Subtraction with subtraction assignment with the given vectors";

         try {
            dres_   -= lhs_ - rhs_;
            sres_   -= lhs_ - rhs_;
            refres_ -= reflhs_ - refrhs_;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed addition assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side dense vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<VT1,VT2>();

         try {
            tdres_   -= tlhs_ - trhs_;
            tsres_   -= tlhs_ - trhs_;
            trefres_ -= treflhs_ - trefrhs_;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed addition assignment operation\n"
                << " Details:\n"
                << "   Transpose left-hand side sparse vector type:\n"
                << "     " << typeid( TVT1 ).name() << "\n"
                << "   Transpose right-hand side dense vector type:\n"
                << "     " << typeid( TVT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Subtraction with subtraction assignment with evaluated vectors
      {
         test_ = "Subtraction with subtraction assignment with evaluated vectors";

         try {
            dres_   -= eval( lhs_ ) - eval( rhs_ );
            sres_   -= eval( lhs_ ) - eval( rhs_ );
            refres_ -= eval( reflhs_ ) - eval( refrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed addition assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side dense vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<VT1,VT2>();

         try {
            tdres_   -= eval( tlhs_ ) - eval( trhs_ );
            tsres_   -= eval( tlhs_ ) - eval( trhs_ );
            trefres_ -= eval( treflhs_ ) - eval( trefrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed addition assignment operation\n"
                << " Details:\n"
                << "   Transpose left-hand side sparse vector type:\n"
                << "     " << typeid( TVT1 ).name() << "\n"
                << "   Transpose right-hand side dense vector type:\n"
                << "     " << typeid( TVT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Subtraction with multiplication assignment
      //=====================================================================================

      // Subtraction with multiplication assignment with the given vectors
      {
         test_ = "Subtraction with multiplication assignment with the given vectors";

         try {
            dres_   *= lhs_ - rhs_;
            sres_   *= lhs_ - rhs_;
            refres_ *= reflhs_ - refrhs_;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side dense vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<VT1,VT2>();

         try {
            tdres_   *= tlhs_ - trhs_;
            tsres_   *= tlhs_ - trhs_;
            trefres_ *= treflhs_ - trefrhs_;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication assignment operation\n"
                << " Details:\n"
                << "   Transpose left-hand side sparse vector type:\n"
                << "     " << typeid( TVT1 ).name() << "\n"
                << "   Transpose right-hand side dense vector type:\n"
                << "     " << typeid( TVT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Subtraction with multiplication assignment with evaluated vectors
      {
         test_ = "Subtraction with multiplication assignment with evaluated vectors";

         try {
            dres_   *= eval( lhs_ ) - eval( rhs_ );
            sres_   *= eval( lhs_ ) - eval( rhs_ );
            refres_ *= eval( reflhs_ ) - eval( refrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side dense vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<VT1,VT2>();

         try {
            tdres_   *= eval( tlhs_ ) - eval( trhs_ );
            tsres_   *= eval( tlhs_ ) - eval( trhs_ );
            trefres_ *= eval( treflhs_ ) - eval( trefrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication assignment operation\n"
                << " Details:\n"
                << "   Transpose left-hand side sparse vector type:\n"
                << "     " << typeid( TVT1 ).name() << "\n"
                << "   Transpose right-hand side dense vector type:\n"
                << "     " << typeid( TVT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults<TVT1,TVT2>();
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the negated sparse vector/dense vector subtraction.
//
// \return void
// \exception std::runtime_error Subtraction error detected.
//
// This function tests the negated vector subtraction with plain assignment, addition assignment,
// subtraction assignment, and multiplication assignment. In case any error resulting from the
// subtraction or the subsequent assignment is detected, a \a std::runtime_error exception is
// thrown.
*/
template< typename VT1    // Type of the left-hand side sparse vector
        , typename VT2 >  // Type of the right-hand side dense vector
void SVecDVecSub<VT1,VT2>::testNegatedOperation()
{
#if BLAZETEST_MATHTEST_TEST_NEGATED_OPERATION
   if( BLAZETEST_MATHTEST_TEST_NEGATED_OPERATION > 1 )
   {
      //=====================================================================================
      // Negated subtraction
      //=====================================================================================

      // Negated subtraction with the given vectors
      {
         test_ = "Negated subtraction with the givven types";

         try {
            dres_   = -( lhs_ - rhs_ );
            sres_   = -( lhs_ - rhs_ );
            refres_ = -( reflhs_ - refrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side dense vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<VT1,VT2>();

         try {
            tdres_   = -( tlhs_ - trhs_ );
            tsres_   = -( tlhs_ - trhs_ );
            trefres_ = -( treflhs_ - trefrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction operation\n"
                << " Details:\n"
                << "   Transpose left-hand side sparse vector type:\n"
                << "     " << typeid( TVT1 ).name() << "\n"
                << "   Transpose right-hand side dense vector type:\n"
                << "     " << typeid( TVT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Negated subtraction with evaluated vectors
      {
         test_ = "Negated subtraction with evaluated vectors";

         try {
            dres_ = -( eval( lhs_ ) - eval( rhs_ ) );
            sres_ = -( eval( lhs_ ) - eval( rhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side dense vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<VT1,VT2>();

         try {
            tdres_ = -( eval( tlhs_ ) - eval( trhs_ ) );
            tsres_ = -( eval( tlhs_ ) - eval( trhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction operation\n"
                << " Details:\n"
                << "   Transpose left-hand side sparse vector type:\n"
                << "     " << typeid( TVT1 ).name() << "\n"
                << "   Transpose right-hand side dense vector type:\n"
                << "     " << typeid( TVT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Negated subtraction with addition assignment
      //=====================================================================================

      // Negated subtraction with addition assignment with the given vectors
      {
         test_ = "Negated subtraction with addition assignment with the given vectors";

         try {
            dres_   += -( lhs_ - rhs_ );
            sres_   += -( lhs_ - rhs_ );
            refres_ += -( reflhs_ - refrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed addition assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side dense vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<VT1,VT2>();

         try {
            tdres_   += -( tlhs_ - trhs_ );
            tsres_   += -( tlhs_ - trhs_ );
            trefres_ += -( treflhs_ - trefrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed addition assignment operation\n"
                << " Details:\n"
                << "   Transpose left-hand side sparse vector type:\n"
                << "     " << typeid( TVT1 ).name() << "\n"
                << "   Transpose right-hand side dense vector type:\n"
                << "     " << typeid( TVT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Negated subtraction with addition assignment with evaluated vectors
      {
         test_ = "Negated subtraction with addition assignment with evaluated vectors";

         try {
            dres_   += -( eval( lhs_ ) - eval( rhs_ ) );
            sres_   += -( eval( lhs_ ) - eval( rhs_ ) );
            refres_ += -( eval( reflhs_ ) - eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed addition assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side dense vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<VT1,VT2>();

         try {
            tdres_   += -( eval( tlhs_ ) - eval( trhs_ ) );
            tsres_   += -( eval( tlhs_ ) - eval( trhs_ ) );
            trefres_ += -( eval( treflhs_ ) - eval( trefrhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed addition assignment operation\n"
                << " Details:\n"
                << "   Transpose left-hand side sparse vector type:\n"
                << "     " << typeid( TVT1 ).name() << "\n"
                << "   Transpose right-hand side dense vector type:\n"
                << "     " << typeid( TVT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Negated subtraction with subtraction assignment
      //=====================================================================================

      // Negated subtraction with subtraction assignment with the given vectors
      {
         test_ = "Negated subtraction with subtraction assignment with the given vectors";

         try {
            dres_   -= -( lhs_ - rhs_ );
            sres_   -= -( lhs_ - rhs_ );
            refres_ -= -( reflhs_ - refrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side dense vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<VT1,VT2>();
         try {
            tdres_   -= -( tlhs_ - trhs_ );
            tsres_   -= -( tlhs_ - trhs_ );
            trefres_ -= -( treflhs_ - trefrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction assignment operation\n"
                << " Details:\n"
                << "   Transpose left-hand side sparse vector type:\n"
                << "     " << typeid( TVT1 ).name() << "\n"
                << "   Transpose right-hand side dense vector type:\n"
                << "     " << typeid( TVT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Negated subtraction with subtraction assignment with evaluated vectors
      {
         test_ = "Negated subtraction with subtraction assignment with evaluated vectors";

         try {
            dres_   -= -( eval( lhs_ ) - eval( rhs_ ) );
            sres_   -= -( eval( lhs_ ) - eval( rhs_ ) );
            refres_ -= -( eval( reflhs_ ) - eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side dense vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<VT1,VT2>();

         try {
            tdres_   -= -( eval( tlhs_ ) - eval( trhs_ ) );
            tsres_   -= -( eval( tlhs_ ) - eval( trhs_ ) );
            trefres_ -= -( eval( treflhs_ ) - eval( trefrhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction assignment operation\n"
                << " Details:\n"
                << "   Transpose left-hand side sparse vector type:\n"
                << "     " << typeid( TVT1 ).name() << "\n"
                << "   Transpose right-hand side dense vector type:\n"
                << "     " << typeid( TVT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Negated subtraction with multiplication assignment
      //=====================================================================================

      // Negated subtraction with multiplication assignment with the given vectors
      {
         test_ = "Negated subtraction with multiplication assignment with the given vectors";

         try {
            dres_   *= -( lhs_ - rhs_ );
            sres_   *= -( lhs_ - rhs_ );
            refres_ *= -( reflhs_ - refrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side dense vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<VT1,VT2>();

         try {
            tdres_   *= -( tlhs_ - trhs_ );
            tsres_   *= -( tlhs_ - trhs_ );
            trefres_ *= -( treflhs_ - trefrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication assignment operation\n"
                << " Details:\n"
                << "   Transpose left-hand side sparse vector type:\n"
                << "     " << typeid( TVT1 ).name() << "\n"
                << "   Transpose right-hand side dense vector type:\n"
                << "     " << typeid( TVT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Negated subtraction with multiplication assignment with evaluated vectors
      {
         test_ = "Negated subtraction with multiplication assignment with evaluated vectors";

         try {
            dres_   *= -( eval( lhs_ ) - eval( rhs_ ) );
            sres_   *= -( eval( lhs_ ) - eval( rhs_ ) );
            refres_ *= -( eval( reflhs_ ) - eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side dense vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<VT1,VT2>();

         try {
            tdres_   *= -( eval( tlhs_ ) - eval( trhs_ ) );
            tsres_   *= -( eval( tlhs_ ) - eval( trhs_ ) );
            trefres_ *= -( eval( treflhs_ ) - eval( trefrhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication assignment operation\n"
                << " Details:\n"
                << "   Transpose left-hand side sparse vector type:\n"
                << "     " << typeid( TVT1 ).name() << "\n"
                << "   Transpose right-hand side dense vector type:\n"
                << "     " << typeid( TVT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults<TVT1,TVT2>();
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the scaled sparse vector/dense vector subtraction.
//
// \param scalar The scalar value.
// \return void
// \exception std::runtime_error Subtraction error detected.
//
// This function tests the scaled vector subtraction with plain assignment, addition assignment,
// subtraction assignment, and multiplication assignment. In case any error resulting from the
// subtraction or the subsequent assignment is detected, a \a std::runtime_error exception is
// thrown.
*/
template< typename VT1    // Type of the left-hand side sparse vector
        , typename VT2 >  // Type of the right-hand side dense vector
template< typename T >    // Type of the scalar
void SVecDVecSub<VT1,VT2>::testScaledOperation( T scalar )
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
            dres_   = lhs_ - rhs_;
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

         checkResults<VT1,VT2>();
      }


      //=====================================================================================
      // Self-scaling (OP/=s)
      //=====================================================================================

      // Self-scaling (OP/=s)
      {
         test_ = "Self-scaling (OP/=s)";

         try {
            dres_   = lhs_ - rhs_;
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

         checkResults<VT1,VT2>();
      }


      //=====================================================================================
      // Scaled subtraction (s*OP)
      //=====================================================================================

      // Scaled subtraction with the given vectors
      {
         test_ = "Scaled subtraction with the given vectors (s*OP)";

         try {
            dres_   = scalar * ( lhs_ - rhs_ );
            sres_   = scalar * ( lhs_ - rhs_ );
            refres_ = scalar * ( reflhs_ - refrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side dense vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<VT1,VT2>();

         try {
            tdres_   = scalar * ( tlhs_ - trhs_ );
            tsres_   = scalar * ( tlhs_ - trhs_ );
            trefres_ = scalar * ( treflhs_ - trefrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction operation\n"
                << " Details:\n"
                << "   Transpose left-hand side sparse vector type:\n"
                << "     " << typeid( TVT1 ).name() << "\n"
                << "   Transpose right-hand side dense vector type:\n"
                << "     " << typeid( TVT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Scaled subtraction with evaluated vectors
      {
         test_ = "Scaled subtraction with evaluated vectors (s*OP)";

         try {
            dres_ = scalar * ( eval( lhs_ ) - eval( rhs_ ) );
            sres_ = scalar * ( eval( lhs_ ) - eval( rhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side dense vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<VT1,VT2>();

         try {
            tdres_ = scalar * ( eval( tlhs_ ) - eval( trhs_ ) );
            tsres_ = scalar * ( eval( tlhs_ ) - eval( trhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction operation\n"
                << " Details:\n"
                << "   Transpose left-hand side sparse vector type:\n"
                << "     " << typeid( TVT1 ).name() << "\n"
                << "   Transpose right-hand side dense vector type:\n"
                << "     " << typeid( TVT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Scaled subtraction (OP*s)
      //=====================================================================================

      // Scaled subtraction with the given vectors
      {
         test_ = "Scaled subtraction with the given vectors (OP*s)";

         try {
            dres_   = ( lhs_ - rhs_ ) * scalar;
            sres_   = ( lhs_ - rhs_ ) * scalar;
            refres_ = ( reflhs_ - refrhs_ ) * scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side dense vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<VT1,VT2>();

         try {
            tdres_   = ( tlhs_ - trhs_ ) * scalar;
            tsres_   = ( tlhs_ - trhs_ ) * scalar;
            trefres_ = ( treflhs_ - trefrhs_ ) * scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction operation\n"
                << " Details:\n"
                << "   Transpose left-hand side sparse vector type:\n"
                << "     " << typeid( TVT1 ).name() << "\n"
                << "   Transpose right-hand side dense vector type:\n"
                << "     " << typeid( TVT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Scaled subtraction with evaluated vectors
      {
         test_ = "Scaled subtraction with evaluated vectors (OP*s)";

         try {
            dres_ = ( eval( lhs_ ) - eval( rhs_ ) ) * scalar;
            sres_ = ( eval( lhs_ ) - eval( rhs_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side dense vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<VT1,VT2>();

         try {
            tdres_ = ( eval( tlhs_ ) - eval( trhs_ ) ) * scalar;
            tsres_ = ( eval( tlhs_ ) - eval( trhs_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction operation\n"
                << " Details:\n"
                << "   Transpose left-hand side sparse vector type:\n"
                << "     " << typeid( TVT1 ).name() << "\n"
                << "   Transpose right-hand side dense vector type:\n"
                << "     " << typeid( TVT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Scaled subtraction (OP/s)
      //=====================================================================================

      // Scaled subtraction with the given vectors
      {
         test_ = "Scaled subtraction with the given vectors (OP/s)";

         try {
            dres_   = ( lhs_ - rhs_ ) / scalar;
            sres_   = ( lhs_ - rhs_ ) / scalar;
            refres_ = ( reflhs_ - refrhs_ ) / scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side dense vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<VT1,VT2>();

         try {
            tdres_   = ( tlhs_ - trhs_ ) / scalar;
            tsres_   = ( tlhs_ - trhs_ ) / scalar;
            trefres_ = ( treflhs_ - trefrhs_ ) / scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction operation\n"
                << " Details:\n"
                << "   Transpose left-hand side sparse vector type:\n"
                << "     " << typeid( TVT1 ).name() << "\n"
                << "   Transpose right-hand side dense vector type:\n"
                << "     " << typeid( TVT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Scaled subtraction with evaluated vectors
      {
         test_ = "Scaled subtraction with evaluated vectors (OP/s)";

         try {
            dres_ = ( eval( lhs_ ) - eval( rhs_ ) ) / scalar;
            sres_ = ( eval( lhs_ ) - eval( rhs_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side dense vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<VT1,VT2>();

         try {
            tdres_ = ( eval( tlhs_ ) - eval( trhs_ ) ) / scalar;
            tsres_ = ( eval( tlhs_ ) - eval( trhs_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction operation\n"
                << " Details:\n"
                << "   Transpose left-hand side sparse vector type:\n"
                << "     " << typeid( TVT1 ).name() << "\n"
                << "   Transpose right-hand side dense vector type:\n"
                << "     " << typeid( TVT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Scaled subtraction with addition assignment (s*OP)
      //=====================================================================================

      // Scaled subtraction with addition assignment with the given vectors
      {
         test_ = "Scaled subtraction with addition assignment with the given vectors (s*OP)";

         try {
            dres_   += scalar * ( lhs_ - rhs_ );
            sres_   += scalar * ( lhs_ - rhs_ );
            refres_ += scalar * ( reflhs_ - refrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed addition assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side dense vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<VT1,VT2>();

         try {
            tdres_   += scalar * ( tlhs_ - trhs_ );
            tsres_   += scalar * ( tlhs_ - trhs_ );
            trefres_ += scalar * ( treflhs_ - trefrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed addition assignment operation\n"
                << " Details:\n"
                << "   Transpose left-hand side sparse vector type:\n"
                << "     " << typeid( TVT1 ).name() << "\n"
                << "   Transpose right-hand side dense vector type:\n"
                << "     " << typeid( TVT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Scaled subtraction with addition assignment with evaluated vectors
      {
         test_ = "Scaled subtraction with addition assignment with evaluated vectors (s*OP)";

         try {
            dres_   += scalar * ( eval( lhs_ ) - eval( rhs_ ) );
            sres_   += scalar * ( eval( lhs_ ) - eval( rhs_ ) );
            refres_ += scalar * ( eval( reflhs_ ) - eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed addition assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side dense vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<VT1,VT2>();

         try {
            tdres_   += scalar * ( eval( tlhs_ ) - eval( trhs_ ) );
            tsres_   += scalar * ( eval( tlhs_ ) - eval( trhs_ ) );
            trefres_ += scalar * ( eval( treflhs_ ) - eval( trefrhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed addition assignment operation\n"
                << " Details:\n"
                << "   Transpose left-hand side sparse vector type:\n"
                << "     " << typeid( TVT1 ).name() << "\n"
                << "   Transpose right-hand side dense vector type:\n"
                << "     " << typeid( TVT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Scaled subtraction with addition assignment (OP*s)
      //=====================================================================================

      // Scaled subtraction with addition assignment with the given vectors
      {
         test_ = "Scaled subtraction with addition assignment with the given vectors (OP*s)";

         try {
            dres_   += ( lhs_ - rhs_ ) * scalar;
            sres_   += ( lhs_ - rhs_ ) * scalar;
            refres_ += ( reflhs_ - refrhs_ ) * scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed addition assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side dense vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<VT1,VT2>();

         try {
            tdres_   += ( tlhs_ - trhs_ ) * scalar;
            tsres_   += ( tlhs_ - trhs_ ) * scalar;
            trefres_ += ( treflhs_ - trefrhs_ ) * scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed addition assignment operation\n"
                << " Details:\n"
                << "   Transpose left-hand side sparse vector type:\n"
                << "     " << typeid( TVT1 ).name() << "\n"
                << "   Transpose right-hand side dense vector type:\n"
                << "     " << typeid( TVT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Scaled subtraction with addition assignment with evaluated vectors
      {
         test_ = "Scaled subtraction with addition assignment with evaluated vectors (OP*s)";

         try {
            dres_   += ( eval( lhs_ ) - eval( rhs_ ) ) * scalar;
            sres_   += ( eval( lhs_ ) - eval( rhs_ ) ) * scalar;
            refres_ += ( eval( reflhs_ ) - eval( refrhs_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed addition assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side dense vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<VT1,VT2>();

         try {
            tdres_   += ( eval( tlhs_ ) - eval( trhs_ ) ) * scalar;
            tsres_   += ( eval( tlhs_ ) - eval( trhs_ ) ) * scalar;
            trefres_ += ( eval( treflhs_ ) - eval( trefrhs_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed addition assignment operation\n"
                << " Details:\n"
                << "   Transpose left-hand side sparse vector type:\n"
                << "     " << typeid( TVT1 ).name() << "\n"
                << "   Transpose right-hand side dense vector type:\n"
                << "     " << typeid( TVT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Scaled subtraction with addition assignment (OP/s)
      //=====================================================================================

      // Scaled subtraction with addition assignment with the given vectors
      {
         test_ = "Scaled subtraction with addition assignment with the given vectors (OP/s)";

         try {
            dres_   += ( lhs_ - rhs_ ) / scalar;
            sres_   += ( lhs_ - rhs_ ) / scalar;
            refres_ += ( reflhs_ - refrhs_ ) / scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed addition assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side dense vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<VT1,VT2>();

         try {
            tdres_   += ( tlhs_ - trhs_ ) / scalar;
            tsres_   += ( tlhs_ - trhs_ ) / scalar;
            trefres_ += ( treflhs_ - trefrhs_ ) / scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed addition assignment operation\n"
                << " Details:\n"
                << "   Transpose left-hand side sparse vector type:\n"
                << "     " << typeid( TVT1 ).name() << "\n"
                << "   Transpose right-hand side dense vector type:\n"
                << "     " << typeid( TVT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Scaled subtraction with addition assignment with evaluated vectors
      {
         test_ = "Scaled subtraction with addition assignment with evaluated vectors (OP/s)";

         try {
            dres_   += ( eval( lhs_ ) - eval( rhs_ ) ) / scalar;
            sres_   += ( eval( lhs_ ) - eval( rhs_ ) ) / scalar;
            refres_ += ( eval( reflhs_ ) - eval( refrhs_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed addition assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side dense vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<VT1,VT2>();

         try {
            tdres_   += ( eval( tlhs_ ) - eval( trhs_ ) ) / scalar;
            tsres_   += ( eval( tlhs_ ) - eval( trhs_ ) ) / scalar;
            trefres_ += ( eval( treflhs_ ) - eval( trefrhs_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed addition assignment operation\n"
                << " Details:\n"
                << "   Transpose left-hand side sparse vector type:\n"
                << "     " << typeid( TVT1 ).name() << "\n"
                << "   Transpose right-hand side dense vector type:\n"
                << "     " << typeid( TVT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Scaled subtraction with subtraction assignment (s*OP)
      //=====================================================================================

      // Scaled subtraction with subtraction assignment with the given vectors
      {
         test_ = "Scaled subtraction with subtraction assignment with the given vectors (s*OP)";

         try {
            dres_   -= scalar * ( lhs_ - rhs_ );
            sres_   -= scalar * ( lhs_ - rhs_ );
            refres_ -= scalar * ( reflhs_ - refrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side dense vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<VT1,VT2>();

         try {
            tdres_   -= scalar * ( tlhs_ - trhs_ );
            tsres_   -= scalar * ( tlhs_ - trhs_ );
            trefres_ -= scalar * ( treflhs_ - trefrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction assignment operation\n"
                << " Details:\n"
                << "   Transpose left-hand side sparse vector type:\n"
                << "     " << typeid( TVT1 ).name() << "\n"
                << "   Transpose right-hand side dense vector type:\n"
                << "     " << typeid( TVT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Scaled subtraction with subtraction assignment with evaluated vectors
      {
         test_ = "Scaled subtraction with subtraction assignment with evaluated vectors (s*OP)";

         try {
            dres_   -= scalar * ( eval( lhs_ ) - eval( rhs_ ) );
            sres_   -= scalar * ( eval( lhs_ ) - eval( rhs_ ) );
            refres_ -= scalar * ( eval( reflhs_ ) - eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side dense vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<VT1,VT2>();

         try {
            tdres_   -= scalar * ( eval( tlhs_ ) - eval( trhs_ ) );
            tsres_   -= scalar * ( eval( tlhs_ ) - eval( trhs_ ) );
            trefres_ -= scalar * ( eval( treflhs_ ) - eval( trefrhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction assignment operation\n"
                << " Details:\n"
                << "   Transpose left-hand side sparse vector type:\n"
                << "     " << typeid( TVT1 ).name() << "\n"
                << "   Transpose right-hand side dense vector type:\n"
                << "     " << typeid( TVT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Scaled subtraction with subtraction assignment (OP*s)
      //=====================================================================================

      // Scaled subtraction with subtraction assignment with the given vectors
      {
         test_ = "Scaled subtraction with subtraction assignment with the given vectors (OP*s)";

         try {
            dres_   -= ( lhs_ - rhs_ ) * scalar;
            sres_   -= ( lhs_ - rhs_ ) * scalar;
            refres_ -= ( reflhs_ - refrhs_ ) * scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side dense vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<VT1,VT2>();

         try {
            tdres_   -= ( tlhs_ - trhs_ ) * scalar;
            tsres_   -= ( tlhs_ - trhs_ ) * scalar;
            trefres_ -= ( treflhs_ - trefrhs_ ) * scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction assignment operation\n"
                << " Details:\n"
                << "   Transpose left-hand side sparse vector type:\n"
                << "     " << typeid( TVT1 ).name() << "\n"
                << "   Transpose right-hand side dense vector type:\n"
                << "     " << typeid( TVT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Scaled subtraction with subtraction assignment with evaluated vectors
      {
         test_ = "Scaled subtraction with subtraction assignment with evaluated vectors (OP*s)";

         try {
            dres_   -= ( eval( lhs_ ) - eval( rhs_ ) ) * scalar;
            sres_   -= ( eval( lhs_ ) - eval( rhs_ ) ) * scalar;
            refres_ -= ( eval( reflhs_ ) - eval( refrhs_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side dense vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<VT1,VT2>();

         try {
            tdres_   -= ( eval( tlhs_ ) - eval( trhs_ ) ) * scalar;
            tsres_   -= ( eval( tlhs_ ) - eval( trhs_ ) ) * scalar;
            trefres_ -= ( eval( treflhs_ ) - eval( trefrhs_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction assignment operation\n"
                << " Details:\n"
                << "   Transpose left-hand side sparse vector type:\n"
                << "     " << typeid( TVT1 ).name() << "\n"
                << "   Transpose right-hand side dense vector type:\n"
                << "     " << typeid( TVT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Scaled subtraction with subtraction assignment (OP/s)
      //=====================================================================================

      // Scaled subtraction with subtraction assignment with the given vectors
      {
         test_ = "Scaled subtraction with subtraction assignment with the given vectors (OP/s)";

         try {
            dres_   -= ( lhs_ - rhs_ ) / scalar;
            sres_   -= ( lhs_ - rhs_ ) / scalar;
            refres_ -= ( reflhs_ - refrhs_ ) / scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side dense vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<VT1,VT2>();

         try {
            tdres_   -= ( tlhs_ - trhs_ ) / scalar;
            tsres_   -= ( tlhs_ - trhs_ ) / scalar;
            trefres_ -= ( treflhs_ - trefrhs_ ) / scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction assignment operation\n"
                << " Details:\n"
                << "   Transpose left-hand side sparse vector type:\n"
                << "     " << typeid( TVT1 ).name() << "\n"
                << "   Transpose right-hand side dense vector type:\n"
                << "     " << typeid( TVT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Scaled subtraction with subtraction assignment with evaluated vectors
      {
         test_ = "Scaled subtraction with subtraction assignment with evaluated vectors (OP/s)";

         try {
            dres_   -= ( eval( lhs_ ) - eval( rhs_ ) ) / scalar;
            sres_   -= ( eval( lhs_ ) - eval( rhs_ ) ) / scalar;
            refres_ -= ( eval( reflhs_ ) - eval( refrhs_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side dense vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<VT1,VT2>();

         try {
            tdres_   -= ( eval( tlhs_ ) - eval( trhs_ ) ) / scalar;
            tsres_   -= ( eval( tlhs_ ) - eval( trhs_ ) ) / scalar;
            trefres_ -= ( eval( treflhs_ ) - eval( trefrhs_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction assignment operation\n"
                << " Details:\n"
                << "   Transpose left-hand side sparse vector type:\n"
                << "     " << typeid( TVT1 ).name() << "\n"
                << "   Transpose right-hand side dense vector type:\n"
                << "     " << typeid( TVT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Scaled subtraction with multiplication assignment (s*OP)
      //=====================================================================================

      // Scaled subtraction with multiplication assignment with the given vectors
      {
         test_ = "Scaled subtraction with multiplication assignment with the given vectors (s*OP)";

         try {
            dres_   *= scalar * ( lhs_ - rhs_ );
            sres_   *= scalar * ( lhs_ - rhs_ );
            refres_ *= scalar * ( reflhs_ - refrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side dense vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<VT1,VT2>();

         try {
            tdres_   *= scalar * ( tlhs_ - trhs_ );
            tsres_   *= scalar * ( tlhs_ - trhs_ );
            trefres_ *= scalar * ( treflhs_ - trefrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication assignment operation\n"
                << " Details:\n"
                << "   Transpose left-hand side sparse vector type:\n"
                << "     " << typeid( TVT1 ).name() << "\n"
                << "   Transpose right-hand side dense vector type:\n"
                << "     " << typeid( TVT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Scaled subtraction with multiplication assignment with evaluated vectors
      {
         test_ = "Scaled subtraction with multiplication assignment with evaluated vectors (s*OP)";

         try {
            dres_   *= scalar * ( eval( lhs_ ) - eval( rhs_ ) );
            sres_   *= scalar * ( eval( lhs_ ) - eval( rhs_ ) );
            refres_ *= scalar * ( eval( reflhs_ ) - eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side dense vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<VT1,VT2>();

         try {
            tdres_   *= scalar * ( eval( tlhs_ ) - eval( trhs_ ) );
            tsres_   *= scalar * ( eval( tlhs_ ) - eval( trhs_ ) );
            trefres_ *= scalar * ( eval( treflhs_ ) - eval( trefrhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication assignment operation\n"
                << " Details:\n"
                << "   Transpose left-hand side sparse vector type:\n"
                << "     " << typeid( TVT1 ).name() << "\n"
                << "   Transpose right-hand side dense vector type:\n"
                << "     " << typeid( TVT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Scaled subtraction with multiplication assignment (OP*s)
      //=====================================================================================

      // Scaled subtraction with multiplication assignment with the given vectors
      {
         test_ = "Scaled subtraction with multiplication assignment with the given vectors (OP*s)";

         try {
            dres_   *= ( lhs_ - rhs_ ) * scalar;
            sres_   *= ( lhs_ - rhs_ ) * scalar;
            refres_ *= ( reflhs_ - refrhs_ ) * scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side dense vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<VT1,VT2>();

         try {
            tdres_   *= ( tlhs_ - trhs_ ) * scalar;
            tsres_   *= ( tlhs_ - trhs_ ) * scalar;
            trefres_ *= ( treflhs_ - trefrhs_ ) * scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication assignment operation\n"
                << " Details:\n"
                << "   Transpose left-hand side sparse vector type:\n"
                << "     " << typeid( TVT1 ).name() << "\n"
                << "   Transpose right-hand side dense vector type:\n"
                << "     " << typeid( TVT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Scaled subtraction with multiplication assignment with evaluated vectors
      {
         test_ = "Scaled subtraction with multiplication assignment with evaluated vectors (OP*s)";

         try {
            dres_   *= ( eval( lhs_ ) - eval( rhs_ ) ) * scalar;
            sres_   *= ( eval( lhs_ ) - eval( rhs_ ) ) * scalar;
            refres_ *= ( eval( reflhs_ ) - eval( refrhs_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side dense vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<VT1,VT2>();

         try {
            tdres_   *= ( eval( tlhs_ ) - eval( trhs_ ) ) * scalar;
            tsres_   *= ( eval( tlhs_ ) - eval( trhs_ ) ) * scalar;
            trefres_ *= ( eval( treflhs_ ) - eval( trefrhs_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication assignment operation\n"
                << " Details:\n"
                << "   Transpose left-hand side sparse vector type:\n"
                << "     " << typeid( TVT1 ).name() << "\n"
                << "   Transpose right-hand side dense vector type:\n"
                << "     " << typeid( TVT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Scaled subtraction with multiplication assignment (OP/s)
      //=====================================================================================

      // Scaled subtraction with multiplication assignment with the given vectors
      {
         test_ = "Scaled subtraction with multiplication assignment with the given vectors (OP/s)";

         try {
            dres_   *= ( lhs_ - rhs_ ) / scalar;
            sres_   *= ( lhs_ - rhs_ ) / scalar;
            refres_ *= ( reflhs_ - refrhs_ ) / scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side dense vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<VT1,VT2>();

         try {
            tdres_   *= ( tlhs_ - trhs_ ) / scalar;
            tsres_   *= ( tlhs_ - trhs_ ) / scalar;
            trefres_ *= ( treflhs_ - trefrhs_ ) / scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication assignment operation\n"
                << " Details:\n"
                << "   Transpose left-hand side sparse vector type:\n"
                << "     " << typeid( TVT1 ).name() << "\n"
                << "   Transpose right-hand side dense vector type:\n"
                << "     " << typeid( TVT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Scaled subtraction with multiplication assignment with evaluated vectors
      {
         test_ = "Scaled subtraction with multiplication assignment with evaluated vectors (OP/s)";

         try {
            dres_   *= ( eval( lhs_ ) - eval( rhs_ ) ) / scalar;
            sres_   *= ( eval( lhs_ ) - eval( rhs_ ) ) / scalar;
            refres_ *= ( eval( reflhs_ ) - eval( refrhs_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side dense vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<VT1,VT2>();

         try {
            tdres_   *= ( eval( tlhs_ ) - eval( trhs_ ) ) / scalar;
            tsres_   *= ( eval( tlhs_ ) - eval( trhs_ ) ) / scalar;
            trefres_ *= ( eval( treflhs_ ) - eval( trefrhs_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication assignment operation\n"
                << " Details:\n"
                << "   Transpose left-hand side sparse vector type:\n"
                << "     " << typeid( TVT1 ).name() << "\n"
                << "   Transpose right-hand side dense vector type:\n"
                << "     " << typeid( TVT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults<TVT1,TVT2>();
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the transpose sparse vector/dense vector subtraction.
//
// \return void
// \exception std::runtime_error Subtraction error detected.
//
// This function tests the transpose vector subtraction with plain assignment, addition assignment,
// subtraction assignment, and multiplication assignment. In case any error resulting from the
// subtraction or the subsequent assignment is detected, a \a std::runtime_error exception is
// thrown.
*/
template< typename VT1    // Type of the left-hand side sparse vector
        , typename VT2 >  // Type of the right-hand side dense vector
void SVecDVecSub<VT1,VT2>::testTransposeOperation()
{
#if BLAZETEST_MATHTEST_TEST_TRANSPOSE_OPERATION
   if( BLAZETEST_MATHTEST_TEST_TRANSPOSE_OPERATION > 1 )
   {
      //=====================================================================================
      // Transpose subtraction
      //=====================================================================================

      // Transpose subtraction with the given vectors
      {
         test_ = "Transpose subtraction with the given vectors";

         try {
            tdres_   = trans( lhs_ - rhs_ );
            tsres_   = trans( lhs_ - rhs_ );
            trefres_ = trans( reflhs_ - refrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side dense vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults<VT1,VT2>();

         try {
            dres_   = trans( tlhs_ - trhs_ );
            sres_   = trans( tlhs_ - trhs_ );
            refres_ = trans( treflhs_ - trefrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction operation\n"
                << " Details:\n"
                << "   Transpose left-hand side sparse vector type:\n"
                << "     " << typeid( TVT1 ).name() << "\n"
                << "   Transpose right-hand side dense vector type:\n"
                << "     " << typeid( TVT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<TVT1,TVT2>();
      }

      // Transpose subtraction with evaluated vectors
      {
         test_ = "Transpose subtraction with evaluated vectors";

         try {
            tdres_ = trans( eval( lhs_ ) - eval( rhs_ ) );
            tsres_ = trans( eval( lhs_ ) - eval( rhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side dense vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults<VT1,VT2>();

         try {
            dres_ = trans( eval( tlhs_ ) - eval( trhs_ ) );
            sres_ = trans( eval( tlhs_ ) - eval( trhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction operation\n"
                << " Details:\n"
                << "   Transpose left-hand side sparse vector type:\n"
                << "     " << typeid( TVT1 ).name() << "\n"
                << "   Transpose right-hand side dense vector type:\n"
                << "     " << typeid( TVT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Transpose subtraction with addition assignment
      //=====================================================================================

      // Transpose subtraction with addition assignment with the given vectors
      {
         test_ = "Transpose subtraction with addition assignment with the given vectors";

         try {
            tdres_   += trans( lhs_ - rhs_ );
            tsres_   += trans( lhs_ - rhs_ );
            trefres_ += trans( reflhs_ - refrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed addition assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side dense vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults<VT1,VT2>();

         try {
            dres_   += trans( tlhs_ - trhs_ );
            sres_   += trans( tlhs_ - trhs_ );
            refres_ += trans( treflhs_ - trefrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed addition assignment operation\n"
                << " Details:\n"
                << "   Transpose left-hand side sparse vector type:\n"
                << "     " << typeid( TVT1 ).name() << "\n"
                << "   Transpose right-hand side dense vector type:\n"
                << "     " << typeid( TVT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<TVT1,TVT2>();
      }

      // Transpose subtraction with addition assignment with evaluated vectors
      {
         test_ = "Transpose subtraction with addition assignment with evaluated vectors";

         try {
            tdres_   += trans( eval( lhs_ ) - eval( rhs_ ) );
            tsres_   += trans( eval( lhs_ ) - eval( rhs_ ) );
            trefres_ += trans( eval( reflhs_ ) - eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed addition assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side dense vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults<VT1,VT2>();

         try {
            dres_   += trans( eval( tlhs_ ) - eval( trhs_ ) );
            sres_   += trans( eval( tlhs_ ) - eval( trhs_ ) );
            refres_ += trans( eval( treflhs_ ) - eval( trefrhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed addition assignment operation\n"
                << " Details:\n"
                << "   Transpose left-hand side sparse vector type:\n"
                << "     " << typeid( TVT1 ).name() << "\n"
                << "   Transpose right-hand side dense vector type:\n"
                << "     " << typeid( TVT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Transpose subtraction with subtraction assignment
      //=====================================================================================

      // Transpose subtraction with subtraction assignment with the given vectors
      {
         test_ = "Transpose subtraction with subtraction assignment with the given vectors";

         try {
            tdres_   -= trans( lhs_ - rhs_ );
            tsres_   -= trans( lhs_ - rhs_ );
            trefres_ -= trans( reflhs_ - refrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side dense vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults<VT1,VT2>();

         try {
            dres_   -= trans( tlhs_ - trhs_ );
            sres_   -= trans( tlhs_ - trhs_ );
            refres_ -= trans( treflhs_ - trefrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction assignment operation\n"
                << " Details:\n"
                << "   Transpose left-hand side sparse vector type:\n"
                << "     " << typeid( TVT1 ).name() << "\n"
                << "   Transpose right-hand side dense vector type:\n"
                << "     " << typeid( TVT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<TVT1,TVT2>();
      }

      // Transpose subtraction with subtraction assignment with evaluated vectors
      {
         test_ = "Transpose subtraction with subtraction assignment with evaluated vectors";

         try {
            tdres_   -= trans( eval( lhs_ ) - eval( rhs_ ) );
            tsres_   -= trans( eval( lhs_ ) - eval( rhs_ ) );
            trefres_ -= trans( eval( reflhs_ ) - eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side dense vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults<VT1,VT2>();

         try {
            dres_   -= trans( eval( tlhs_ ) - eval( trhs_ ) );
            sres_   -= trans( eval( tlhs_ ) - eval( trhs_ ) );
            refres_ -= trans( eval( treflhs_ ) - eval( trefrhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction assignment operation\n"
                << " Details:\n"
                << "   Transpose left-hand side sparse vector type:\n"
                << "     " << typeid( TVT1 ).name() << "\n"
                << "   Transpose right-hand side dense vector type:\n"
                << "     " << typeid( TVT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Transpose subtraction with multiplication assignment
      //=====================================================================================

      // Transpose subtraction with multiplication assignment with the given vectors
      {
         test_ = "Transpose subtraction with multiplication assignment with the given vectors";

         try {
            tdres_   *= trans( lhs_ - rhs_ );
            tsres_   *= trans( lhs_ - rhs_ );
            trefres_ *= trans( reflhs_ - refrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side dense vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults<VT1,VT2>();

         try {
            dres_   *= trans( tlhs_ - trhs_ );
            sres_   *= trans( tlhs_ - trhs_ );
            refres_ *= trans( treflhs_ - trefrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication assignment operation\n"
                << " Details:\n"
                << "   Transpose left-hand side sparse vector type:\n"
                << "     " << typeid( TVT1 ).name() << "\n"
                << "   Transpose right-hand side dense vector type:\n"
                << "     " << typeid( TVT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<TVT1,TVT2>();
      }

      // Transpose subtraction with multiplication assignment with evaluated vectors
      {
         test_ = "Transpose subtraction with multiplication assignment with evaluated vectors";

         try {
            tdres_   *= trans( eval( lhs_ ) - eval( rhs_ ) );
            tsres_   *= trans( eval( lhs_ ) - eval( rhs_ ) );
            trefres_ *= trans( eval( reflhs_ ) - eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side dense vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults<VT1,VT2>();

         try {
            dres_   *= trans( eval( tlhs_ ) - eval( trhs_ ) );
            sres_   *= trans( eval( tlhs_ ) - eval( trhs_ ) );
            refres_ *= trans( eval( treflhs_ ) - eval( trefrhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication assignment operation\n"
                << " Details:\n"
                << "   Transpose left-hand side sparse vector type:\n"
                << "     " << typeid( TVT1 ).name() << "\n"
                << "   Transpose right-hand side dense vector type:\n"
                << "     " << typeid( TVT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<TVT1,TVT2>();
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the abs sparse vector/dense vector subtraction.
//
// \return void
// \exception std::runtime_error Subtraction error detected.
//
// This function tests the abs vector subtraction with plain assignment, addition assignment,
// subtraction assignment, and multiplication assignment. In case any error resulting from
// the subtraction or the subsequent assignment is detected, a \a std::runtime_error exception
// is thrown.
*/
template< typename VT1    // Type of the left-hand side sparse vector
        , typename VT2 >  // Type of the right-hand side dense vector
void SVecDVecSub<VT1,VT2>::testAbsOperation()
{
#if BLAZETEST_MATHTEST_TEST_ABS_OPERATION
   if( BLAZETEST_MATHTEST_TEST_ABS_OPERATION > 1 )
   {
      //=====================================================================================
      // Abs subtraction
      //=====================================================================================

      // Abs subtraction with the given vectors
      {
         test_ = "Abs subtraction with the given vectors";

         try {
            dres_   = abs( lhs_ - rhs_ );
            sres_   = abs( lhs_ - rhs_ );
            refres_ = abs( reflhs_ - refrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side dense vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<VT1,VT2>();

         try {
            tdres_   = abs( tlhs_ - trhs_ );
            tsres_   = abs( tlhs_ - trhs_ );
            trefres_ = abs( treflhs_ - trefrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction operation\n"
                << " Details:\n"
                << "   Transpose left-hand side sparse vector type:\n"
                << "     " << typeid( TVT1 ).name() << "\n"
                << "   Transpose right-hand side dense vector type:\n"
                << "     " << typeid( TVT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Abs subtraction with evaluated vectors
      {
         test_ = "Abs subtraction with evaluated vectors";

         try {
            dres_ = abs( eval( lhs_ ) - eval( rhs_ ) );
            sres_ = abs( eval( lhs_ ) - eval( rhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side dense vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<VT1,VT2>();

         try {
            tdres_ = abs( eval( tlhs_ ) - eval( trhs_ ) );
            tsres_ = abs( eval( tlhs_ ) - eval( trhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction operation\n"
                << " Details:\n"
                << "   Transpose left-hand side sparse vector type:\n"
                << "     " << typeid( TVT1 ).name() << "\n"
                << "   Transpose right-hand side dense vector type:\n"
                << "     " << typeid( TVT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Abs subtraction with addition assignment
      //=====================================================================================

      // Abs subtraction with addition assignment with the given vectors
      {
         test_ = "Abs subtraction with addition assignment with the given vectors";

         try {
            dres_   += abs( lhs_ - rhs_ );
            sres_   += abs( lhs_ - rhs_ );
            refres_ += abs( reflhs_ - refrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed addition assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side dense vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<VT1,VT2>();

         try {
            tdres_   += abs( tlhs_ - trhs_ );
            tsres_   += abs( tlhs_ - trhs_ );
            trefres_ += abs( treflhs_ - trefrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed addition assignment operation\n"
                << " Details:\n"
                << "   Transpose left-hand side sparse vector type:\n"
                << "     " << typeid( TVT1 ).name() << "\n"
                << "   Transpose right-hand side dense vector type:\n"
                << "     " << typeid( TVT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Abs subtraction with addition assignment with evaluated vectors
      {
         test_ = "Abs subtraction with addition assignment with evaluated vectors";

         try {
            dres_   += abs( eval( lhs_ ) - eval( rhs_ ) );
            sres_   += abs( eval( lhs_ ) - eval( rhs_ ) );
            refres_ += abs( eval( reflhs_ ) - eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed addition assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side dense vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<VT1,VT2>();

         try {
            tdres_   += abs( eval( tlhs_ ) - eval( trhs_ ) );
            tsres_   += abs( eval( tlhs_ ) - eval( trhs_ ) );
            trefres_ += abs( eval( treflhs_ ) - eval( trefrhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed addition assignment operation\n"
                << " Details:\n"
                << "   Transpose left-hand side sparse vector type:\n"
                << "     " << typeid( TVT1 ).name() << "\n"
                << "   Transpose right-hand side dense vector type:\n"
                << "     " << typeid( TVT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Abs subtraction with subtraction assignment
      //=====================================================================================

      // Abs subtraction with subtraction assignment with the given vectors
      {
         test_ = "Abs subtraction with subtraction assignment with the given types";

         try {
            dres_   -= abs( lhs_ - rhs_ );
            sres_   -= abs( lhs_ - rhs_ );
            refres_ -= abs( reflhs_ - refrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side dense vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<VT1,VT2>();

         try {
            tdres_   -= abs( tlhs_ - trhs_ );
            tsres_   -= abs( tlhs_ - trhs_ );
            trefres_ -= abs( treflhs_ - trefrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction assignment operation\n"
                << " Details:\n"
                << "   Transpose left-hand side sparse vector type:\n"
                << "     " << typeid( TVT1 ).name() << "\n"
                << "   Transpose right-hand side dense vector type:\n"
                << "     " << typeid( TVT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Abs subtraction with subtraction assignment with evaluated vectors
      {
         test_ = "Abs subtraction with subtraction assignment with evaluated vectors";

         try {
            dres_   -= abs( eval( lhs_ ) - eval( rhs_ ) );
            sres_   -= abs( eval( lhs_ ) - eval( rhs_ ) );
            refres_ -= abs( eval( reflhs_ ) - eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side dense vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<VT1,VT2>();

         try {
            tdres_   -= abs( eval( tlhs_ ) - eval( trhs_ ) );
            tsres_   -= abs( eval( tlhs_ ) - eval( trhs_ ) );
            trefres_ -= abs( eval( treflhs_ ) - eval( trefrhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction assignment operation\n"
                << " Details:\n"
                << "   Transpose left-hand side sparse vector type:\n"
                << "     " << typeid( TVT1 ).name() << "\n"
                << "   Transpose right-hand side dense vector type:\n"
                << "     " << typeid( TVT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Abs subtraction with multiplication assignment
      //=====================================================================================

      // Abs subtraction with multiplication assignment with the given vectors
      {
         test_ = "Abs subtraction with multiplication assignment with the given vectors";

         try {
            dres_   *= abs( lhs_ - rhs_ );
            sres_   *= abs( lhs_ - rhs_ );
            refres_ *= abs( reflhs_ - refrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side dense vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<VT1,VT2>();

         try {
            tdres_   *= abs( tlhs_ - trhs_ );
            tsres_   *= abs( tlhs_ - trhs_ );
            trefres_ *= abs( treflhs_ - trefrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication assignment operation\n"
                << " Details:\n"
                << "   Transpose left-hand side sparse vector type:\n"
                << "     " << typeid( TVT1 ).name() << "\n"
                << "   Transpose right-hand side dense vector type:\n"
                << "     " << typeid( TVT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Abs subtraction with multiplication assignment with evaluated vectors
      {
         test_ = "Abs subtraction with multiplication assignment with evaluated vectors";

         try {
            dres_   *= abs( eval( lhs_ ) - eval( rhs_ ) );
            sres_   *= abs( eval( lhs_ ) - eval( rhs_ ) );
            refres_ *= abs( eval( reflhs_ ) - eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side dense vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<VT1,VT2>();

         try {
            tdres_   *= abs( eval( tlhs_ ) - eval( trhs_ ) );
            tsres_   *= abs( eval( tlhs_ ) - eval( trhs_ ) );
            trefres_ *= abs( eval( treflhs_ ) - eval( trefrhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication assignment operation\n"
                << " Details:\n"
                << "   Transpose left-hand side sparse vector type:\n"
                << "     " << typeid( TVT1 ).name() << "\n"
                << "   Transpose right-hand side dense vector type:\n"
                << "     " << typeid( TVT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults<TVT1,TVT2>();
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
template< typename VT1    // Type of the left-hand side sparse vector
        , typename VT2 >  // Type of the right-hand side dense vector
template< typename LT     // Type of the left-hand side operand
        , typename RT >   // Type of the right-hand side operan
void SVecDVecSub<VT1,VT2>::checkResults()
{
   using blaze::IsTransposeVector;

   if( !isEqual( dres_, refres_ ) ) {
      std::ostringstream oss;
      oss.precision( 20 );
      oss << " Test : " << test_ << "\n"
          << " Error: Incorrect dense result vector detected\n"
          << " Details:\n"
          << "   " << ( IsTransposeVector<LT>::value ? ( "Transpose l" ) : ( "L" ) ) << "eft-hand side sparse vector type:\n"
          << "     " << typeid( VT1 ).name() << "\n"
          << "   " << ( IsTransposeVector<RT>::value ? ( "Transpose r" ) : ( "R" ) ) << "ight-hand side dense vector type:\n"
          << "     " << typeid( VT2 ).name() << "\n"
          << "   Result:\n" << dres_ << "\n"
          << "   Expected result:\n" << refres_ << "\n";
      throw std::runtime_error( oss.str() );
   }

   if( !isEqual( sres_, refres_ ) ) {
      std::ostringstream oss;
      oss.precision( 20 );
      oss << " Test : " << test_ << "\n"
          << " Error: Incorrect sparse result vector detected\n"
          << " Details:\n"
          << "   " << ( IsTransposeVector<LT>::value ? ( "Transpose l" ) : ( "L" ) ) << "eft-hand side sparse vector type:\n"
          << "     " << typeid( VT1 ).name() << "\n"
          << "   " << ( IsTransposeVector<RT>::value ? ( "Transpose r" ) : ( "R" ) ) << "ight-hand side dense vector type:\n"
          << "     " << typeid( VT2 ).name() << "\n"
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
// results. The two template arguments \a LT and \a RT indicate the types of the left-hand
// side and right-hand side operands used for the computations.
*/
template< typename VT1    // Type of the left-hand side sparse vector
        , typename VT2 >  // Type of the right-hand side dense vector
template< typename LT     // Type of the left-hand side operand
        , typename RT >   // Type of the right-hand side operand
void SVecDVecSub<VT1,VT2>::checkTransposeResults()
{
   using blaze::IsTransposeVector;

   if( !isEqual( tdres_, trefres_ ) ) {
      std::ostringstream oss;
      oss.precision( 20 );
      oss << " Test : " << test_ << "\n"
          << " Error: Incorrect dense result vector detected\n"
          << " Details:\n"
          << "   " << ( IsTransposeVector<LT>::value ? ( "Transpose l" ) : ( "L" ) ) << "eft-hand side sparse vector type:\n"
          << "     " << typeid( LT ).name() << "\n"
          << "   " << ( IsTransposeVector<RT>::value ? ( "Transpose r" ) : ( "R" ) ) << "ight-hand side dense vector type:\n"
          << "     " << typeid( RT ).name() << "\n"
          << "   Result:\n" << tdres_ << "\n"
          << "   Expected result:\n" << trefres_ << "\n";
      throw std::runtime_error( oss.str() );
   }

   if( !isEqual( tsres_, trefres_ ) ) {
      std::ostringstream oss;
      oss.precision( 20 );
      oss << " Test : " << test_ << "\n"
          << " Error: Incorrect sparse result vector detected\n"
          << " Details:\n"
          << "   " << ( IsTransposeVector<LT>::value ? ( "Transpose l" ) : ( "L" ) ) << "eft-hand side sparse vector type:\n"
          << "     " << typeid( LT ).name() << "\n"
          << "   " << ( IsTransposeVector<RT>::value ? ( "Transpose r" ) : ( "R" ) ) << "ight-hand side dense vector type:\n"
          << "     " << typeid( RT ).name() << "\n"
          << "   Result:\n" << tsres_ << "\n"
          << "   Expected result:\n" << trefres_ << "\n";
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
/*!\brief Testing the vector subtraction between two specific vector types.
//
// \param creator1 The creator for the left-hand side sparse vector.
// \param creator2 The creator for the right-hand side dense vector.
// \return void
*/
template< typename VT1    // Type of the left-hand side sparse vector
        , typename VT2 >  // Type of the right-hand side dense vector
void runTest( const Creator<VT1>& creator1, const Creator<VT2>& creator2 )
{
   for( size_t rep=0; rep<repetitions; ++rep ) {
      SVecDVecSub<VT1,VT2>( creator1, creator2 );
   }
}
//*************************************************************************************************




//=================================================================================================
//
//  MACRO DEFINITIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Macro for the definition of a sparse vector/dense vector subtraction test case.
*/
#define DEFINE_SVECDVECSUB_TEST( VT1, VT2 ) \
   extern template class blazetest::mathtest::svecdvecsub::SVecDVecSub<VT1,VT2>
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Macro for the execution of a sparse vector/dense vector subtraction test case.
*/
#define RUN_SVECDVECSUB_TEST( C1, C2 ) \
   blazetest::mathtest::svecdvecsub::runTest( C1, C2 )
/*! \endcond */
//*************************************************************************************************

} // namespace svecdvecsub

} // namespace mathtest

} // namespace blazetest

#endif
