//=================================================================================================
/*!
//  \file blazetest/mathtest/intrinsics/OperationTest.h
//  \brief Header file for the intrinsics operation test
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

#ifndef _BLAZETEST_MATHTEST_INTRINSICS_OPERATIONTEST_H_
#define _BLAZETEST_MATHTEST_INTRINSICS_OPERATIONTEST_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <sstream>
#include <stdexcept>
#include <string>
#include <typeinfo>
#include <blaze/math/Intrinsics.h>
#include <blaze/math/shims/Equal.h>
#include <blaze/util/constraints/Numeric.h>
#include <blaze/util/Memory.h>
#include <blaze/util/NonCopyable.h>
#include <blaze/util/policies/Deallocate.h>
#include <blaze/util/Random.h>
#include <blaze/util/UniqueArray.h>


namespace blazetest {

namespace mathtest {

namespace intrinsics {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Auxiliary class template for the intrinsics operation test.
//
// This class template represents the tests of all available intrinsics operations for the given
// numeric data type \a T. In these tests both aligned and unaligned load/store operations are
// used.
*/
template< typename T >  // Data type of the intrinsic test
class OperationTest : private blaze::NonCopyable
{
 private:
   //**Type definitions****************************************************************************
   typedef blaze::IntrinsicTrait<T>  IT;             //!< Intrinsic trait for the given numeric type.
   typedef typename IT::Type         IntrinsicType;  //!< Intrinsic type for the given numeric type.

   typedef blaze::UniqueArray<T,blaze::Deallocate>  Ptr;  //!< Smart pointer for the aligned memory.
   //**********************************************************************************************

   //**********************************************************************************************
   static const size_t N  = 256;           //!< Number of numeric values to be worked on.
   static const size_t NN = N + IT::size;  //!< Total number of numeric values in each array.
   //**********************************************************************************************

 public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit OperationTest();
   //@}
   //**********************************************************************************************

 private:
   //**Test functions******************************************************************************
   /*!\name Test functions */
   //@{
   void testStore    ();
   void testStream   ();
   void testStoreu   ( size_t offset );
   void testReduction();
   //@}
   //**********************************************************************************************

   //**Error detection functions*******************************************************************
   /*!\name Error detection functions */
   //@{
   void compare( const T* a, const T* b ) const;
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
   T* a_;  //!< The first aligned array of size NN.
   T* b_;  //!< The second aligned array of size NN.

   std::string test_;  //!< Label of the currently performed test.
   //@}
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_NUMERIC_TYPE( T );
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
/*!\brief Constructor for the intrinsics operation test.
//
// \exception std::runtime_error Operation error detected.
*/
template< typename T >  // Data type of the intrinsic test
OperationTest<T>::OperationTest()
   : a_    ( blaze::allocate<T>( NN ) )  // The first aligned array of size NN
   , b_    ( blaze::allocate<T>( NN ) )  // The second aligned array of size NN
   , test_ ()                            // Label of the currently performed test
{
   testStore();
   testStream();

   for( size_t offset=0UL; offset<IT::size; ++offset ) {
      testStoreu( offset );
   }

   testReduction();
}
//*************************************************************************************************




//=================================================================================================
//
//  TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Testing the aligned store operation.
//
// \return void
// \exception std::runtime_error Load/store error detected.
//
// This function tests the aligned store operation by copying one array to another via aligned
// load and store. In case any error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename T >  // Data type of the intrinsic test
void OperationTest<T>::testStore()
{
   using blaze::load;
   using blaze::store;

   test_  = "store() operation";

   initialize();

   for( size_t i=0UL; i<N; i+=IT::size ) {
      store( b_+i, load( a_+i ) );
   }

   compare( a_, b_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the aligned, non-temporal store operation.
//
// \return void
// \exception std::runtime_error Load/store error detected.
//
// This function tests the aligned, non-temporal store operation by copying one array to another
// via aligned load and non-temporal store. In case any error is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename T >  // Data type of the intrinsic test
void OperationTest<T>::testStream()
{
   using blaze::load;
   using blaze::stream;

   test_  = "stream() operation";

   initialize();

   for( size_t i=0UL; i<N; i+=IT::size ) {
      stream( b_+i, load( a_+i ) );
   }

   compare( a_, b_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the unaligned store operation.
//
// \return void
// \exception std::runtime_error Load/store error detected.
//
// This function tests the unaligned store operation by copying one array to another via unaligned
// load and store. In case any error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename T >  // Data type of the intrinsic test
void OperationTest<T>::testStoreu( size_t offset )
{
   using blaze::loadu;
   using blaze::storeu;

   test_  = "storeu() operation";

   initialize();

   for( size_t i=0UL; i<N; i+=IT::size ) {
      storeu( b_+offset+i, loadu( a_+offset+i ) );
   }

   compare( a_+offset, b_+offset );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the reduction operation.
//
// \return void
// \exception std::runtime_error Reduction error detected.
//
// This function tests the reduction operation by comparing the results of a vectorized and a
// scalar reduction. In case any error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename T >  // Data type of the intrinsic test
void OperationTest<T>::testReduction()
{
   using blaze::load;
   using blaze::sum;

   test_ = "sum() operation";

   initialize();

   T ssum = T();
   for( size_t i=0UL; i<N; ++i ) {
      ssum += a_[i];
   }

   T vsum = T();
   for( size_t i=0UL; i<N; i+=IT::size ) {
      vsum += sum( load( a_+i ) );
   }

   if( !blaze::equal( ssum, vsum ) ) {
      std::ostringstream oss;
      oss.precision( 20 );
      oss << " Test : " << test_ << "\n"
          << " Error: Failed reduction operation\n"
          << " Details:\n"
          << "   ssum = " << ssum << "\n"
          << "   vsum = " << vsum << "\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************




//=================================================================================================
//
//  ERROR DETECTION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Comparison of the first 256 elements of the two given arrays.
//
// \param a The first array to be compared.
// \param b The second array to be compared.
// \return void
// \exception std::runtime_error Value mismatch detected.
//
// This function compares the first 256 elements of the two given arrays. In case any value of
// the two arrays differs, a \a std::runtime_error exception is thrown.
*/
template< typename T >  // Data type of the intrinsic test
void OperationTest<T>::compare( const T* a, const T* b ) const
{
   for( size_t i=0UL; i<N; ++i ) {
      if( a[i] != b[i] ) {
         std::ostringstream oss;
         oss.precision( 20 );
         oss << " Test : " << test_ << "\n"
             << " Error: Value mismatch detected at index " << i << "\n"
             << " Details:\n"
             << "   a[" << i << "] = " << a[i] << "\n"
             << "   b[" << i << "] = " << b[i] << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Initialization of all member arrays.
//
// \return void
//
// This function is called before each single test case to initialize all arrays with random
// values.
*/
template< typename T >  // Data type of the intrinsic test
void OperationTest<T>::initialize()
{
   using blaze::randomize;

   for( size_t i=0UL; i<NN; ++i ) {
      randomize( a_[i] );
      randomize( b_[i] );
   }
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Testing the intrinsics operations of a specific numeric data type.
//
// \return void
*/
template< typename T >  // Data type of the intrinsic test
void runTest()
{
   OperationTest<T>();
}
//*************************************************************************************************




//=================================================================================================
//
//  MACRO DEFINITIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Macro for the execution of an intrinsics operation test case.
*/
#define RUN_INTRINSICS_OPERATION_TEST( T ) \
   blazetest::mathtest::intrinsics::runTest<T>()
/*! \endcond */
//*************************************************************************************************

} // namespace intrinsics

} // namespace mathtest

} // namespace blazetest

#endif
