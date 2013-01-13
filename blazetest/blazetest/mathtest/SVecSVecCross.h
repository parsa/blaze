//=================================================================================================
/*!
//  \file blazetest/mathtest/SVecSVecCross.h
//  \brief Header file for the sparse vector/sparse vector cross product math test
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

#ifndef _BLAZETEST_MATHTEST_SVECSVECCROSS_H_
#define _BLAZETEST_MATHTEST_SVECSVECCROSS_H_


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
#include <blaze/math/shims/Equal.h>
#include <blaze/math/StaticVector.h>
#include <blaze/math/traits/CrossTrait.h>
#include <blaze/math/typetraits/IsTransposeVector.h>
#include <blaze/util/constraints/Numeric.h>
#include <blaze/util/constraints/SameType.h>
#include <blazetest/system/MathTest.h>
#include <blazetest/util/Creator.h>
#include <blazetest/util/Utility.h>


namespace blazetest {

namespace mathtest {

namespace svecsveccross {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Auxiliary class template for the sparse vector/sparse vector cross product math test.
//
// The SVecSVecCross class template represents one particular vector cross product test between
// two vectors of a particular type. The two template arguments \a VT1 and \a VT2 represent the
// types of the left-hand side and right-hand side vector, respectively.
*/
template< typename VT1    // Type of the left-hand side sparse vector
        , typename VT2 >  // Type of the right-hand side sparse vector
class SVecSVecCross
{
 private:
   //**Type definitions****************************************************************************
   typedef typename blaze::CrossTrait<VT1,VT2>::Type  RE;   //!< Default result type
   typedef typename RE::TransposeType                 TRE;  //!< Transpose default result type
   //**********************************************************************************************

   //**Type definitions****************************************************************************
   typedef typename VT1::ElementType           ET1;    //!< Element type 1
   typedef typename VT2::ElementType           ET2;    //!< Element type 2
   typedef typename RE::ElementType            RET;    //!< Resulting element type
   typedef blaze::DynamicVector<ET1,false>     RT1;    //!< Reference type 1
   typedef blaze::DynamicVector<ET2,false>     RT2;    //!< Reference type 2
   typedef blaze::StaticVector<RET,3UL,false>  DRRE;   //!< Dense reference result type
   typedef blaze::CompressedVector<RET,false>  SRRE;   //!< Sparse reference result type
   typedef typename DRRE::TransposeType        TDRRE;  //!< Transpose dense reference result type
   typedef typename SRRE::TransposeType        TSRRE;  //!< Transpose sparse reference result type
   typedef RE                                  DRE;    //!< Dense result type
   typedef SRRE                                SRE;    //!< Sparse result type
   typedef TRE                                 TDRE;   //!< Transpose dense result type
   typedef TSRRE                               TSRE;   //!< Transpose sparse result type
   //**********************************************************************************************

 public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit SVecSVecCross( const Creator<VT1>& creator1, const Creator<VT2>& creator2 );
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
   void checkResults();
   void checkTransposeResults();
   //@}
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   VT1   lhs_;      //!< The left-hand side sparse vector.
   VT2   rhs_;      //!< The right-hand side sparse vector.
   RT1   reflhs_;   //!< The reference left-hand side vector.
   RT2   refrhs_;   //!< The reference right-hand side vector.
   DRE   dres_;     //!< The dense vector for the result of the vector cross product.
   SRE   sres_;     //!< The sparse vector for the result of the vector cross product.
   DRRE  refres_;   //!< The reference result.
   TDRE  tdres_;    //!< The dense vector for the result of the transpose vector cross product.
   TSRE  tsres_;    //!< The sparse vector for the result of the transpose vector cross product.
   TDRRE trefres_;  //!< The transpose reference result.

   std::string test_;  //!< Label of the currently performed test.
   //@}
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( VT1   );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( VT2   );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( RT1   );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( RT2   );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( DRE   );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( SRE   );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( DRRE  );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( SRRE  );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( TDRE  );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( TSRE  );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( TDRRE );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( TSRRE );
   BLAZE_CONSTRAINT_MUST_BE_NONTRANSPOSE_VECTOR_TYPE( VT1   );
   BLAZE_CONSTRAINT_MUST_BE_NONTRANSPOSE_VECTOR_TYPE( VT2   );
   BLAZE_CONSTRAINT_MUST_BE_NONTRANSPOSE_VECTOR_TYPE( RT1   );
   BLAZE_CONSTRAINT_MUST_BE_NONTRANSPOSE_VECTOR_TYPE( RT2   );
   BLAZE_CONSTRAINT_MUST_BE_NONTRANSPOSE_VECTOR_TYPE( DRE   );
   BLAZE_CONSTRAINT_MUST_BE_NONTRANSPOSE_VECTOR_TYPE( SRE   );
   BLAZE_CONSTRAINT_MUST_BE_NONTRANSPOSE_VECTOR_TYPE( DRRE  );
   BLAZE_CONSTRAINT_MUST_BE_NONTRANSPOSE_VECTOR_TYPE( SRRE  );
   BLAZE_CONSTRAINT_MUST_BE_TRANSPOSE_VECTOR_TYPE   ( TDRE  );
   BLAZE_CONSTRAINT_MUST_BE_TRANSPOSE_VECTOR_TYPE   ( TSRE  );
   BLAZE_CONSTRAINT_MUST_BE_TRANSPOSE_VECTOR_TYPE   ( TDRRE );
   BLAZE_CONSTRAINT_MUST_BE_TRANSPOSE_VECTOR_TYPE   ( TSRRE );
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
/*!\brief Constructor for the SVecSVecCross class template.
//
// \param creator1 The creator for the left-hand side sparse vector of the vector cross product.
// \param creator2 The creator for the right-hand side sparse vector of the vector cross product.
// \exception std::runtime_error Operation error detected.
*/
template< typename VT1    // Type of the left-hand side sparse vector
        , typename VT2 >  // Type of the right-hand side sparse vector
SVecSVecCross<VT1,VT2>::SVecSVecCross( const Creator<VT1>& creator1, const Creator<VT2>& creator2 )
   : lhs_( creator1() )  // The left-hand side sparse vector
   , rhs_( creator2() )  // The right-hand side sparse vector
   , reflhs_( lhs_ )     // The reference left-hand side vector
   , refrhs_( rhs_ )     // The reference right-hand side vector
   , dres_()             // The dense vector for the result of the vector cross product
   , sres_()             // The sparse vector for the result of the vector cross product
   , refres_()           // The reference result
   , tdres_()            // The dense vector for the result of the transpose vector cross product
   , tsres_()            // The sparse vector for the result of the transpose vector cross product
   , trefres_()          // The transpose reference result
   , test_()             // Label of the currently performed test.
{
   if( lhs_.size() != 3UL ) {
      throw std::runtime_error( "Invalid size of left-hand side operand" );
   }

   if( rhs_.size() != 3UL ) {
      throw std::runtime_error( "Invalid size of right-hand side operand" );
   }

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
        , typename VT2 >  // Type of the right-hand side sparse vector
void SVecSVecCross<VT1,VT2>::testInitialStatus()
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
      oss << " Test: Initial size comparison of right-hand side sparse operand\n"
          << " Error: Invalid vector size\n"
          << " Details:\n"
          << "   Sparse vector type:\n"
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
      oss << " Test: Initial test of initialization of right-hand side sparse operand\n"
          << " Error: Invalid vector initialization\n"
          << " Details:\n"
          << "   Sparse vector type:\n"
          << "     " << typeid( VT2 ).name() << "\n"
          << "   Current initialization:\n" << rhs_ << "\n"
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
template< typename VT1    // Type of the left-hand side sparse vector
        , typename VT2 >  // Type of the right-hand side sparse vector
void SVecSVecCross<VT1,VT2>::testAssignment()
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
          << "   Right-hand side sparse vector type:\n"
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
      oss << " Test: Checking the assignment result of right-hand side sparse operand\n"
          << " Error: Invalid vector initialization\n"
          << " Details:\n"
          << "   Sparse vector type:\n"
          << "     " << typeid( VT2 ).name() << "\n"
          << "   Current initialization:\n" << rhs_ << "\n"
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
template< typename VT1    // Type of the left-hand side sparse vector
        , typename VT2 >  // Type of the right-hand side sparse vector
void SVecSVecCross<VT1,VT2>::testElementAccess()
{
   using blaze::equal;


   //=====================================================================================
   // Testing the element access with the given vectors
   //=====================================================================================

   if( !equal( ( lhs_ % rhs_ )[0UL], ( reflhs_ % refrhs_ )[0UL] ) ) {
      std::ostringstream oss;
      oss << " Test : Element access of cross product expression\n"
          << " Error: Unequal resulting elements at index 0 detected\n"
          << " Details:\n"
          << "   Left-hand side sparse vector type:\n"
          << "     " << typeid( VT1 ).name() << "\n"
          << "   Right-hand side sparse vector type:\n"
          << "     " << typeid( VT2 ).name() << "\n";
      throw std::runtime_error( oss.str() );
   }

   if( !equal( ( lhs_ % eval( rhs_ ) )[0UL], ( reflhs_ % eval( refrhs_ ) )[0UL] ) ) {
      std::ostringstream oss;
      oss << " Test : Element access of right evaluated cross product expression\n"
          << " Error: Unequal resulting elements at index 0 detected\n"
          << " Details:\n"
          << "   Left-hand side sparse vector type:\n"
          << "     " << typeid( VT1 ).name() << "\n"
          << "   Right-hand side sparse vector type:\n"
          << "     " << typeid( VT2 ).name() << "\n";
      throw std::runtime_error( oss.str() );
   }

   if( !equal( ( eval( lhs_ ) % rhs_ )[0UL], ( eval( reflhs_ ) % refrhs_ )[0UL] ) ) {
      std::ostringstream oss;
      oss << " Test : Element access of left evaluated cross product expression\n"
          << " Error: Unequal resulting elements at index 0 detected\n"
          << " Details:\n"
          << "   Left-hand side sparse vector type:\n"
          << "     " << typeid( VT1 ).name() << "\n"
          << "   Right-hand side sparse vector type:\n"
          << "     " << typeid( VT2 ).name() << "\n";
      throw std::runtime_error( oss.str() );
   }

   if( !equal( ( eval( lhs_ ) % eval( rhs_ ) )[0UL], ( eval( reflhs_ ) % eval( refrhs_ ) )[0UL] ) ) {
      std::ostringstream oss;
      oss << " Test : Element access of fully evaluated cross product expression\n"
          << " Error: Unequal resulting elements at index 0 detected\n"
          << " Details:\n"
          << "   Left-hand side sparse vector type:\n"
          << "     " << typeid( VT1 ).name() << "\n"
          << "   Right-hand side sparse vector type:\n"
          << "     " << typeid( VT2 ).name() << "\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the plain sparse vector/sparse vector cross product.
//
// \return void
// \exception std::runtime_error Cross product error detected.
//
// This function tests the plain vector cross product with plain assignment, addition
// assignment, subtraction assignment, and multiplication assignment. In case any error
// resulting from the cros product or the subsequent assignment is detected, a
// \a std::runtime_error exception is thrown.
*/
template< typename VT1    // Type of the left-hand side sparse vector
        , typename VT2 >  // Type of the right-hand side sparse vector
void SVecSVecCross<VT1,VT2>::testBasicOperation()
{
#if BLAZETEST_MATHTEST_TEST_BASIC_OPERATION
   if( BLAZETEST_MATHTEST_TEST_BASIC_OPERATION > 1 )
   {
      //=====================================================================================
      // Cross product with the given vectors
      //=====================================================================================

      // Cross product with the given vectors
      {
         test_ = " Cross product with the given vectors";

         try {
            dres_   = lhs_ % rhs_;
            sres_   = lhs_ % rhs_;
            refres_ = reflhs_ % refrhs_;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed cross product operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side sparse vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults();
      }

      // Cross product with evaluated vectors
      {
         test_ = "Cross product with evaluated vectors";

         try {
            dres_ = eval( lhs_ ) % eval( rhs_ );
            sres_ = eval( lhs_ ) % eval( rhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed cross product operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side sparse vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults();
      }


      //=====================================================================================
      // Cross product with addition assignment
      //=====================================================================================

      // Cross product with addition assignment with the given vectors
      {
         test_ = "Cross product with addition assignment with the given vectors";

         try {
            dres_   += lhs_ % rhs_;
            sres_   += lhs_ % rhs_;
            refres_ += reflhs_ % refrhs_;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed addition assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side sparse vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults();
      }

      // Cross product with addition assignment with the given vectors
      {
         test_ = "Cross product with addition assignment with evaluated vectors";

         try {
            dres_   += eval( lhs_ ) % eval( rhs_ );
            sres_   += eval( lhs_ ) % eval( rhs_ );
            refres_ += eval( reflhs_ ) % eval( refrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed addition assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side sparse vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults();
      }


      //=====================================================================================
      // Cross product with subtraction assignment
      //=====================================================================================

      // Cross product with subtraction assignment with the given vectors
      {
         test_ = "Cross product with subtraction assignment with the given vectors";

         try {
            dres_   -= lhs_ % rhs_;
            sres_   -= lhs_ % rhs_;
            refres_ -= reflhs_ % refrhs_;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side sparse vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults();
      }

      // Cross product with subtraction assignment with evaluated vectors
      {
         test_ = "Cross product with subtraction assignment with evaluated vectors";

         try {
            dres_   -= eval( lhs_ ) % eval( rhs_ );
            sres_   -= eval( lhs_ ) % eval( rhs_ );
            refres_ -= eval( reflhs_ ) % eval( refrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side sparse vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults();
      }


      //=====================================================================================
      // Cross product with multiplication assignment
      //=====================================================================================

      // Cross product with multiplication assignment with the given vectors
      {
         test_ = "Cross product with multiplication assignment with the given vectors";

         try {
            dres_   *= lhs_ % rhs_;
            sres_   *= lhs_ % rhs_;
            refres_ *= reflhs_ % refrhs_;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side sparse vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults();
      }

      // Cross product with multiplication assignment with evaluated vectors
      {
         test_ = "Cross product with multiplication assignment with evaluated vectors";

         try {
            dres_   *= eval( lhs_ ) % eval( rhs_ );
            sres_   *= eval( lhs_ ) % eval( rhs_ );
            refres_ *= eval( reflhs_ ) % eval( refrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side sparse vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the negated sparse vector/sparse vector cross product.
//
// \return void
// \exception std::runtime_error Cross product error detected.
//
// This function tests the negated vector cross product with plain assignment, addition
// assignment, subtraction assignment, and multiplication assignment. In case any error
// resulting from the cross product or the subsequent assignment is detected, a
// \a std::runtime_error exception is thrown.
*/
template< typename VT1    // Type of the left-hand side sparse vector
        , typename VT2 >  // Type of the right-hand side sparse vector
void SVecSVecCross<VT1,VT2>::testNegatedOperation()
{
#if BLAZETEST_MATHTEST_TEST_NEGATED_OPERATION
   if( BLAZETEST_MATHTEST_TEST_NEGATED_OPERATION > 1 )
   {
      //=====================================================================================
      // Negated cross product
      //=====================================================================================

      // Negated cross product with the given vectors
      {
         test_ = "Negated cross product with the givven types";

         try {
            dres_   = -( lhs_ % rhs_ );
            sres_   = -( lhs_ % rhs_ );
            refres_ = -( reflhs_ % refrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed cross product operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side sparse vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults();
      }

      // Negated cross product with evaluated vectors
      {
         test_ = "Negated cross product with evaluated vectors";

         try {
            dres_ = -( eval( lhs_ ) % eval( rhs_ ) );
            sres_ = -( eval( lhs_ ) % eval( rhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed cross product operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side sparse vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults();
      }


      //=====================================================================================
      // Negated cross product with addition assignment
      //=====================================================================================

      // Negated cross product with addition assignment with the given vectors
      {
         test_ = "Negated cross product with addition assignment with the given vectors";

         try {
            dres_   += -( lhs_ % rhs_ );
            sres_   += -( lhs_ % rhs_ );
            refres_ += -( reflhs_ % refrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed addition assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side sparse vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults();
      }

      // Negated cross product with addition assignment with evaluated vectors
      {
         test_ = "Negated cross product with addition assignment with evaluated vectors";

         try {
            dres_   += -( eval( lhs_ ) % eval( rhs_ ) );
            sres_   += -( eval( lhs_ ) % eval( rhs_ ) );
            refres_ += -( eval( reflhs_ ) % eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed addition assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side sparse vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults();
      }


      //=====================================================================================
      // Negated cross product with subtraction assignment
      //=====================================================================================

      // Negated cross product with subtraction assignment with the given vectors
      {
         test_ = "Negated cross product with subtraction assignment with the given vectors";

         try {
            dres_   -= -( lhs_ % rhs_ );
            sres_   -= -( lhs_ % rhs_ );
            refres_ -= -( reflhs_ % refrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side sparse vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults();
      }

      // Negated cross product with subtraction assignment with evaluated vectors
      {
         test_ = "Negated cross product with subtraction assignment with evaluated vectors";

         try {
            dres_   -= -( eval( lhs_ ) % eval( rhs_ ) );
            sres_   -= -( eval( lhs_ ) % eval( rhs_ ) );
            refres_ -= -( eval( reflhs_ ) % eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side sparse vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults();
      }


      //=====================================================================================
      // Negated cross product with multiplication assignment
      //=====================================================================================

      // Negated cross product with multiplication assignment with the given vectors
      {
         test_ = "Negated cross product with multiplication assignment with the given vectors";

         try {
            dres_   *= -( lhs_ % rhs_ );
            sres_   *= -( lhs_ % rhs_ );
            refres_ *= -( reflhs_ % refrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side sparse vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults();
      }

      // Negated cross product with multiplication assignment with evaluated vectors
      {
         test_ = "Negated cross product with multiplication assignment with evaluated vectors";

         try {
            dres_   *= -( eval( lhs_ ) % eval( rhs_ ) );
            sres_   *= -( eval( lhs_ ) % eval( rhs_ ) );
            refres_ *= -( eval( reflhs_ ) % eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side sparse vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults();
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the scaled sparse vector/sparse vector cross product.
//
// \param scalar The scalar value.
// \return void
// \exception std::runtime_error Cross product error detected.
//
// This function tests the scaled vector cross product with plain assignment, addition
// assignment, subtraction assignment, and multiplication assignment. In case any error
// resulting from the cross product or the subsequent assignment is detected, a
// \a std::runtime_error exception is thrown.
*/
template< typename VT1    // Type of the left-hand side sparse vector
        , typename VT2 >  // Type of the right-hand side sparse vector
template< typename T >    // Type of the scalar
void SVecSVecCross<VT1,VT2>::testScaledOperation( T scalar )
{
   BLAZE_CONSTRAINT_MUST_BE_NUMERIC_TYPE( T );

   if( scalar == T(0) )
      throw std::invalid_argument( "Invalid scalar parameter" );


#if BLAZETEST_MATHTEST_TEST_SCALED_OPERATION
   if( BLAZETEST_MATHTEST_TEST_SCALED_OPERATION > 1 )
   {
      //=====================================================================================
      // Self-scaling (v*=s)
      //=====================================================================================

      // Self-scaling (v*=s)
      {
         test_ = "Self-scaling (v*=s)";

         try {
            dres_   = lhs_ % rhs_;
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

         checkResults();
      }
      
      
      //=====================================================================================
      // Self-scaling (v=v*s)
      //=====================================================================================

      // Self-scaling (v=v*s)
      {
         test_ = "Self-scaling (v=v*s)";

         try {
            dres_   = lhs_ % rhs_;
            sres_   = dres_;
            refres_ = dres_;

            dres_   = dres_   * scalar;
            sres_   = sres_   * scalar;
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
      // Self-scaling (v=s*v)
      //=====================================================================================

      // Self-scaling (v=s*v)
      {
         test_ = "Self-scaling (v=s*v)";

         try {
            dres_   = lhs_ % rhs_;
            sres_   = dres_;
            refres_ = dres_;

            dres_   = scalar * dres_;
            sres_   = scalar * sres_;
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
      // Self-scaling (v/=s)
      //=====================================================================================

      // Self-scaling (v/=s)
      {
         test_ = "Self-scaling (v/=s)";

         try {
            dres_   = lhs_ % rhs_;
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

         checkResults();
      }
      
      
      //=====================================================================================
      // Self-scaling (v=v/s)
      //=====================================================================================

      // Self-scaling (v=v/s)
      {
         test_ = "Self-scaling (v=v/s)";

         try {
            dres_   = lhs_ % rhs_;
            sres_   = dres_;
            refres_ = dres_;

            dres_   = dres_   / scalar;
            sres_   = sres_   / scalar;
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
      // Scaled cross product (s*OP)
      //=====================================================================================

      // Scaled cross product with the given vectors
      {
         test_ = "Scaled cross product with the given vectors (s*OP)";

         try {
            dres_   = scalar * ( lhs_ % rhs_ );
            sres_   = scalar * ( lhs_ % rhs_ );
            refres_ = scalar * ( reflhs_ % refrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed cross product operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side sparse vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults();
      }

      // Scaled cross product with evaluated vectors
      {
         test_ = "Scaled cross product with evaluated vectors (s*OP)";

         try {
            dres_ = scalar * ( eval( lhs_ ) % eval( rhs_ ) );
            sres_ = scalar * ( eval( lhs_ ) % eval( rhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed cross product operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side sparse vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults();
      }


      //=====================================================================================
      // Scaled cross product (OP*s)
      //=====================================================================================

      // Scaled cross product with the given vectors
      {
         test_ = "Scaled cross product with the given vectors (OP*s)";

         try {
            dres_   = ( lhs_ % rhs_ ) * scalar;
            sres_   = ( lhs_ % rhs_ ) * scalar;
            refres_ = ( reflhs_ % refrhs_ ) * scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed cross product operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side sparse vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults();
      }

      // Scaled cross product with evaluated vectors
      {
         test_ = "Scaled cross product with evaluated vectors (OP*s)";

         try {
            dres_ = ( eval( lhs_ ) % eval( rhs_ ) ) * scalar;
            sres_ = ( eval( lhs_ ) % eval( rhs_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed cross product operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side sparse vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults();
      }


      //=====================================================================================
      // Scaled cross product (OP/s)
      //=====================================================================================

      // Scaled cross product with the given vectors
      {
         test_ = "Scaled cross product with the given vectors (OP/s)";

         try {
            dres_   = ( lhs_ % rhs_ ) / scalar;
            sres_   = ( lhs_ % rhs_ ) / scalar;
            refres_ = ( reflhs_ % refrhs_ ) / scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed cross product operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side sparse vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults();
      }

      // Scaled cross product with evaluated vectors
      {
         test_ = "Scaled cross product with evaluated vectors (OP/s)";

         try {
            dres_ = ( eval( lhs_ ) % eval( rhs_ ) ) / scalar;
            sres_ = ( eval( lhs_ ) % eval( rhs_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed cross product operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side sparse vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults();
      }


      //=====================================================================================
      // Scaled cross product with addition assignment (s*OP)
      //=====================================================================================

      // Scaled cross product with addition assignment with the given vectors
      {
         test_ = "Scaled cross product with addition assignment with the given vectors (s*OP)";

         try {
            dres_   += scalar * ( lhs_ % rhs_ );
            sres_   += scalar * ( lhs_ % rhs_ );
            refres_ += scalar * ( reflhs_ % refrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed addition assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side sparse vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults();
      }

      // Scaled cross product with addition assignment with evaluated vectors
      {
         test_ = "Scaled cross product with addition assignment with evaluated vectors (s*OP)";

         try {
            dres_   += scalar * ( eval( lhs_ ) % eval( rhs_ ) );
            sres_   += scalar * ( eval( lhs_ ) % eval( rhs_ ) );
            refres_ += scalar * ( eval( reflhs_ ) % eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed addition assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side sparse vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults();
      }


      //=====================================================================================
      // Scaled cross product with addition assignment (OP*s)
      //=====================================================================================

      // Scaled cross product with addition assignment with the given vectors
      {
         test_ = "Scaled cross product with addition assignment with the given vectors (OP*s)";

         try {
            dres_   += ( lhs_ % rhs_ ) * scalar;
            sres_   += ( lhs_ % rhs_ ) * scalar;
            refres_ += ( reflhs_ % refrhs_ ) * scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed addition assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side sparse vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults();
      }

      // Scaled cross product with addition assignment with evaluated vectors
      {
         test_ = "Scaled cross product with addition assignment with evaluated vectors (OP*s)";

         try {
            dres_   += ( eval( lhs_ ) % eval( rhs_ ) ) * scalar;
            sres_   += ( eval( lhs_ ) % eval( rhs_ ) ) * scalar;
            refres_ += ( eval( reflhs_ ) % eval( refrhs_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed addition assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side sparse vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults();
      }


      //=====================================================================================
      // Scaled cross product with addition assignment (OP/s)
      //=====================================================================================

      // Scaled cross product with addition assignment with the given vectors
      {
         test_ = "Scaled cross product with addition assignment with the given vectors (OP/s)";

         try {
            dres_   += ( lhs_ % rhs_ ) / scalar;
            sres_   += ( lhs_ % rhs_ ) / scalar;
            refres_ += ( reflhs_ % refrhs_ ) / scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed addition assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side sparse vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults();
      }

      // Scaled cross product with addition assignment with evaluated vectors
      {
         test_ = "Scaled cross product with addition assignment with evaluated vectors (OP/s)";

         try {
            dres_   += ( eval( lhs_ ) % eval( rhs_ ) ) / scalar;
            sres_   += ( eval( lhs_ ) % eval( rhs_ ) ) / scalar;
            refres_ += ( eval( reflhs_ ) % eval( refrhs_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed addition assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side sparse vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults();
      }


      //=====================================================================================
      // Scaled cross product with subtraction assignment (s*OP)
      //=====================================================================================

      // Scaled cross product with subtraction assignment with the given vectors
      {
         test_ = "Scaled cross product with subtraction assignment with the given vectors (s*OP)";

         try {
            dres_   -= scalar * ( lhs_ % rhs_ );
            sres_   -= scalar * ( lhs_ % rhs_ );
            refres_ -= scalar * ( reflhs_ % refrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side sparse vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults();
      }

      // Scaled cross product with subtraction assignment with evaluated vectors
      {
         test_ = "Scaled cross product with subtraction assignment with evaluated vectors (s*OP)";

         try {
            dres_   -= scalar * ( eval( lhs_ ) % eval( rhs_ ) );
            sres_   -= scalar * ( eval( lhs_ ) % eval( rhs_ ) );
            refres_ -= scalar * ( eval( reflhs_ ) % eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side sparse vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults();
      }


      //=====================================================================================
      // Scaled cross product with subtraction assignment (OP*s)
      //=====================================================================================

      // Scaled cross product with subtraction assignment with the given vectors
      {
         test_ = "Scaled cross product with subtraction assignment with the given vectors (OP*s)";

         try {
            dres_   -= ( lhs_ % rhs_ ) * scalar;
            sres_   -= ( lhs_ % rhs_ ) * scalar;
            refres_ -= ( reflhs_ % refrhs_ ) * scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side sparse vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults();
      }

      // Scaled cross product with subtraction assignment with evaluated vectors
      {
         test_ = "Scaled cross product with subtraction assignment with evaluated vectors (OP*s)";

         try {
            dres_   -= ( eval( lhs_ ) % eval( rhs_ ) ) * scalar;
            sres_   -= ( eval( lhs_ ) % eval( rhs_ ) ) * scalar;
            refres_ -= ( eval( reflhs_ ) % eval( refrhs_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side sparse vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults();
      }


      //=====================================================================================
      // Scaled cross product with subtraction assignment (OP/s)
      //=====================================================================================

      // Scaled cross product with subtraction assignment with the given vectors
      {
         test_ = "Scaled cross product with subtraction assignment with the given vectors (OP/s)";

         try {
            dres_   -= ( lhs_ % rhs_ ) / scalar;
            sres_   -= ( lhs_ % rhs_ ) / scalar;
            refres_ -= ( reflhs_ % refrhs_ ) / scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side sparse vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults();
      }

      // Scaled cross product with subtraction assignment with evaluated vectors
      {
         test_ = "Scaled cross product with subtraction assignment with evaluated vectors (OP/s)";

         try {
            dres_   -= ( eval( lhs_ ) % eval( rhs_ ) ) / scalar;
            sres_   -= ( eval( lhs_ ) % eval( rhs_ ) ) / scalar;
            refres_ -= ( eval( reflhs_ ) % eval( refrhs_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side sparse vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults();
      }


      //=====================================================================================
      // Scaled cross product with multiplication assignment (s*OP)
      //=====================================================================================

      // Scaled cross product with multiplication assignment with the given vectors
      {
         test_ = "Scaled cross product with multiplication assignment with the given vectors (s*OP)";

         try {
            dres_   *= scalar * ( lhs_ % rhs_ );
            sres_   *= scalar * ( lhs_ % rhs_ );
            refres_ *= scalar * ( reflhs_ % refrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side sparse vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults();
      }

      // Scaled cross product with multiplication assignment with evaluated vectors
      {
         test_ = "Scaled cross product with multiplication assignment with evaluated vectors (s*OP)";

         try {
            dres_   *= scalar * ( eval( lhs_ ) % eval( rhs_ ) );
            sres_   *= scalar * ( eval( lhs_ ) % eval( rhs_ ) );
            refres_ *= scalar * ( eval( reflhs_ ) % eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side sparse vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults();
      }


      //=====================================================================================
      // Scaled cross product with multiplication assignment (OP*s)
      //=====================================================================================

      // Scaled cross product with multiplication assignment with the given vectors
      {
         test_ = "Scaled cross product with multiplication assignment with the given vectors (OP*s)";

         try {
            dres_   *= ( lhs_ % rhs_ ) * scalar;
            sres_   *= ( lhs_ % rhs_ ) * scalar;
            refres_ *= ( reflhs_ % refrhs_ ) * scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side sparse vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults();
      }

      // Scaled cross product with multiplication assignment with evaluated vectors
      {
         test_ = "Scaled cross product with multiplication assignment with evaluated vectors (OP*s)";

         try {
            dres_   *= ( eval( lhs_ ) % eval( rhs_ ) ) * scalar;
            sres_   *= ( eval( lhs_ ) % eval( rhs_ ) ) * scalar;
            refres_ *= ( eval( reflhs_ ) % eval( refrhs_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side sparse vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults();
      }


      //=====================================================================================
      // Scaled cross product with multiplication assignment (OP/s)
      //=====================================================================================

      // Scaled cross product with multiplication assignment with the given vectors
      {
         test_ = "Scaled cross product with multiplication assignment with the given vectors (OP/s)";

         try {
            dres_   *= ( lhs_ % rhs_ ) / scalar;
            sres_   *= ( lhs_ % rhs_ ) / scalar;
            refres_ *= ( reflhs_ % refrhs_ ) / scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side sparse vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults();
      }

      // Scaled cross product with multiplication assignment with evaluated vectors
      {
         test_ = "Scaled cross product with multiplication assignment with evaluated vectors (OP/s)";

         try {
            dres_   *= ( eval( lhs_ ) % eval( rhs_ ) ) / scalar;
            sres_   *= ( eval( lhs_ ) % eval( rhs_ ) ) / scalar;
            refres_ *= ( eval( reflhs_ ) % eval( refrhs_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed multiplication assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side sparse vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults();
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the transpose sparse vector/sparse vector cross product.
//
// \return void
// \exception std::runtime_error Cross product error detected.
//
// This function tests the transpose vector cross product with plain assignment, addition
// assignment, subtraction assignment, and multiplication assignment. In case any error resulting
// from the cross product or the subsequent assignment is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename VT1    // Type of the left-hand side sparse vector
        , typename VT2 >  // Type of the right-hand side sparse vector
void SVecSVecCross<VT1,VT2>::testTransposeOperation()
{
#if BLAZETEST_MATHTEST_TEST_TRANSPOSE_OPERATION
   if( BLAZETEST_MATHTEST_TEST_TRANSPOSE_OPERATION > 1 )
   {
      //=====================================================================================
      // Transpose cross product
      //=====================================================================================

      // Transpose cross product with the given vectors
      {
         test_ = "Transpose cross product with the given vectors";

         try {
            tdres_   = trans( lhs_ % rhs_ );
            tsres_   = trans( lhs_ % rhs_ );
            trefres_ = trans( reflhs_ % refrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed cross product operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side sparse vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults();
      }

      // Transpose cross product with evaluated vectors
      {
         test_ = "Transpose cross product with evaluated vectors";

         try {
            tdres_ = trans( eval( lhs_ ) % eval( rhs_ ) );
            tsres_ = trans( eval( lhs_ ) % eval( rhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed cross product operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side sparse vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults();
      }


      //=====================================================================================
      // Transpose cross product with addition assignment
      //=====================================================================================

      // Transpose cross product with addition assignment with the given vectors
      {
         test_ = "Transpose cross product with addition assignment with the given vectors";

         try {
            tdres_   += trans( lhs_ % rhs_ );
            tsres_   += trans( lhs_ % rhs_ );
            trefres_ += trans( reflhs_ % refrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed addition assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side sparse vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults();
      }

      // Transpose cross product with addition assignment with evaluated vectors
      {
         test_ = "Transpose cross product with addition assignment with evaluated vectors";

         try {
            tdres_   += trans( eval( lhs_ ) % eval( rhs_ ) );
            tsres_   += trans( eval( lhs_ ) % eval( rhs_ ) );
            trefres_ += trans( eval( reflhs_ ) % eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed addition assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side sparse vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults();
      }


      //=====================================================================================
      // Transpose cross product with subtraction assignment
      //=====================================================================================

      // Transpose cross product with subtraction assignment with the given vectors
      {
         test_ = "Transpose cross product with subtraction assignment with the given vectors";

         try {
            tdres_   -= trans( lhs_ % rhs_ );
            tsres_   -= trans( lhs_ % rhs_ );
            trefres_ -= trans( reflhs_ % refrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side sparse vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults();
      }

      // Transpose cross product with subtraction assignment with evaluated vectors
      {
         test_ = "Transpose cross product with subtraction assignment with evaluated vectors";

         try {
            tdres_   -= trans( eval( lhs_ ) % eval( rhs_ ) );
            tsres_   -= trans( eval( lhs_ ) % eval( rhs_ ) );
            trefres_ -= trans( eval( reflhs_ ) % eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side sparse vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults();
      }


      //=====================================================================================
      // Transpose cross product with multiplication assignment
      //=====================================================================================

      // Transpose cross product with multiplication assignment with the given vectors
      {
         test_ = "Transpose cross product with multiplication assignment with the given vectors";

         try {
            tdres_   *= trans( lhs_ % rhs_ );
            tsres_   *= trans( lhs_ % rhs_ );
            trefres_ *= trans( reflhs_ % refrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed cross product assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side sparse vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults();
      }

      // Transpose cross product with multiplication assignment with evaluated vectors
      {
         test_ = "Transpose cross product with multiplication assignment with evaluated vectors";

         try {
            tdres_   *= trans( eval( lhs_ ) % eval( rhs_ ) );
            tsres_   *= trans( eval( lhs_ ) % eval( rhs_ ) );
            trefres_ *= trans( eval( reflhs_ ) % eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed cross product assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side sparse vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkTransposeResults();
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the abs sparse vector/sparse vector cross product.
//
// \return void
// \exception std::runtime_error Cross product error detected.
//
// This function tests the abs vector cross product with plain assignment, addition assignment,
// subtraction assignment, and multiplication assignment. In case any error resulting from the
// cross product or the subsequent assignment is detected, a \a std::runtime_error exception
// is thrown.
*/
template< typename VT1    // Type of the left-hand side sparse vector
        , typename VT2 >  // Type of the right-hand side sparse vector
void SVecSVecCross<VT1,VT2>::testAbsOperation()
{
#if BLAZETEST_MATHTEST_TEST_ABS_OPERATION
   if( BLAZETEST_MATHTEST_TEST_ABS_OPERATION > 1 )
   {
      //=====================================================================================
      // Abs cross product
      //=====================================================================================

      // Abs cross product with the given vectors
      {
         test_ = "Abs cross product with the given vectors";

         try {
            dres_   = abs( lhs_ % rhs_ );
            sres_   = abs( lhs_ % rhs_ );
            refres_ = abs( reflhs_ % refrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed cross product operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side sparse vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults();
      }

      // Abs cross product with evaluated vectors
      {
         test_ = "Abs cross product with evaluated vectors";

         try {
            dres_ = abs( eval( lhs_ ) % eval( rhs_ ) );
            sres_ = abs( eval( lhs_ ) % eval( rhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed cross product operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side sparse vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults();
      }


      //=====================================================================================
      // Abs cross product with addition assignment
      //=====================================================================================

      // Abs cross product with addition assignment with the given vectors
      {
         test_ = "Abs cross product with addition assignment with the given vectors";

         try {
            dres_   += abs( lhs_ % rhs_ );
            sres_   += abs( lhs_ % rhs_ );
            refres_ += abs( reflhs_ % refrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed addition assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side sparse vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults();
      }

      // Abs cross product with addition assignment with evaluated vectors
      {
         test_ = "Abs cross product with addition assignment with evaluated vectors";

         try {
            dres_   += abs( eval( lhs_ ) % eval( rhs_ ) );
            sres_   += abs( eval( lhs_ ) % eval( rhs_ ) );
            refres_ += abs( eval( reflhs_ ) % eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed addition assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side sparse vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults();
      }


      //=====================================================================================
      // Abs cross product with subtraction assignment
      //=====================================================================================

      // Abs cross product with subtraction assignment with the given vectors
      {
         test_ = "Abs cross product with subtraction assignment with the given types";

         try {
            dres_   -= abs( lhs_ % rhs_ );
            sres_   -= abs( lhs_ % rhs_ );
            refres_ -= abs( reflhs_ % refrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side sparse vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults();
      }

      // Abs cross product with subtraction assignment with evaluated vectors
      {
         test_ = "Abs cross product with subtraction assignment with evaluated vectors";

         try {
            dres_   -= abs( eval( lhs_ ) % eval( rhs_ ) );
            sres_   -= abs( eval( lhs_ ) % eval( rhs_ ) );
            refres_ -= abs( eval( reflhs_ ) % eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed subtraction assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side sparse vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults();
      }


      //=====================================================================================
      // Abs cross product with multiplication assignment
      //=====================================================================================

      // Abs cross product with multiplication assignment with the given vectors
      {
         test_ = "Abs cross product with multiplication assignment with the given vectors";

         try {
            dres_   *= abs( lhs_ % rhs_ );
            sres_   *= abs( lhs_ % rhs_ );
            refres_ *= abs( reflhs_ % refrhs_ );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed cross product assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side sparse vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults();
      }

      // Abs cross product with multiplication assignment with evaluated vectors
      {
         test_ = "Abs cross product with multiplication assignment with evaluated vectors";

         try {
            dres_   *= abs( eval( lhs_ ) % eval( rhs_ ) );
            sres_   *= abs( eval( lhs_ ) % eval( rhs_ ) );
            refres_ *= abs( eval( reflhs_ ) % eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed cross product assignment operation\n"
                << " Details:\n"
                << "   Left-hand side sparse vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side sparse vector type:\n"
                << "     " << typeid( VT2 ).name() << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
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
// This function is called after each test case to check and compare the computed results. The
// two template arguments \a LT and \a RT indicate the types of the left-hand side and right-hand
// side operands used for the computations.
*/
template< typename VT1    // Type of the left-hand side sparse vector
        , typename VT2 >  // Type of the right-hand side sparse vector
void SVecSVecCross<VT1,VT2>::checkResults()
{
   if( !isEqual( dres_, refres_ ) ) {
      std::ostringstream oss;
      oss.precision( 20 );
      oss << " Test : " << test_ << "\n"
          << " Error: Incorrect dense result vector detected\n"
          << " Details:\n"
          << "   Left-hand side sparse vector type:\n"
          << "     " << typeid( VT1 ).name() << "\n"
          << "   Right-hand side sparse vector type:\n"
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
          << "   Left-hand side sparse vector type:\n"
          << "     " << typeid( VT1 ).name() << "\n"
          << "   Right-hand side sparse vector type:\n"
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
        , typename VT2 >  // Type of the right-hand side sparse vector
void SVecSVecCross<VT1,VT2>::checkTransposeResults()
{
   using blaze::IsTransposeVector;

   if( !isEqual( tdres_, trefres_ ) ) {
      std::ostringstream oss;
      oss.precision( 20 );
      oss << " Test : " << test_ << "\n"
          << " Error: Incorrect dense result vector detected\n"
          << " Details:\n"
          << "   Left-hand side sparse vector type:\n"
          << "     " << typeid( VT1 ).name() << "\n"
          << "   Right-hand side sparse vector type:\n"
          << "     " << typeid( VT2 ).name() << "\n"
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
          << "   Left-hand side sparse vector type:\n"
          << "     " << typeid( VT1 ).name() << "\n"
          << "   Right-hand side sparse vector type:\n"
          << "     " << typeid( VT2 ).name() << "\n"
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
/*!\brief Testing the vector cross product between two specific vector types.
//
// \param creator1 The creator for the left-hand side sparse vector.
// \param creator2 The creator for the right-hand side sparse vector.
// \return void
*/
template< typename VT1    // Type of the left-hand side sparse vector
        , typename VT2 >  // Type of the right-hand side sparse vector
void runTest( const Creator<VT1>& creator1, const Creator<VT2>& creator2 )
{
   for( size_t rep=0; rep<repetitions; ++rep ) {
      SVecSVecCross<VT1,VT2>( creator1, creator2 );
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
/*!\brief Macro for the definition of a sparse vector/sparse vector cross product test case.
*/
#define DEFINE_SVECSVECCROSS_TEST( VT1, VT2 ) \
   extern template class blazetest::mathtest::svecsveccross::SVecSVecCross<VT1,VT2>
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Macro for the execution of a sparse vector/sparse vector cross product test case.
*/
#define RUN_SVECSVECCROSS_TEST( C1, C2 ) \
   blazetest::mathtest::svecsveccross::runTest( C1, C2 )
/*! \endcond */
//*************************************************************************************************

} // namespace svecsveccross

} // namespace mathtest

} // namespace blazetest

#endif
