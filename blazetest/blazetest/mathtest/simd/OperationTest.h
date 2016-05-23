//=================================================================================================
/*!
//  \file blazetest/mathtest/simd/OperationTest.h
//  \brief Header file for the SIMD operation test
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

#ifndef _BLAZETEST_MATHTEST_SIMD_OPERATIONTEST_H_
#define _BLAZETEST_MATHTEST_SIMD_OPERATIONTEST_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <cmath>
#include <sstream>
#include <stdexcept>
#include <string>
#include <typeinfo>
#include <blaze/math/shims/Conjugate.h>
#include <blaze/math/shims/Equal.h>
#include <blaze/math/SIMD.h>
#include <blaze/math/typetraits/HasSIMDAbs.h>
#include <blaze/math/typetraits/HasSIMDAcos.h>
#include <blaze/math/typetraits/HasSIMDAcosh.h>
#include <blaze/math/typetraits/HasSIMDAdd.h>
#include <blaze/math/typetraits/HasSIMDAsin.h>
#include <blaze/math/typetraits/HasSIMDAsinh.h>
#include <blaze/math/typetraits/HasSIMDAtan.h>
#include <blaze/math/typetraits/HasSIMDAtanh.h>
#include <blaze/math/typetraits/HasSIMDCbrt.h>
#include <blaze/math/typetraits/HasSIMDCeil.h>
#include <blaze/math/typetraits/HasSIMDConj.h>
#include <blaze/math/typetraits/HasSIMDCos.h>
#include <blaze/math/typetraits/HasSIMDCosh.h>
#include <blaze/math/typetraits/HasSIMDDiv.h>
#include <blaze/math/typetraits/HasSIMDErf.h>
#include <blaze/math/typetraits/HasSIMDErfc.h>
#include <blaze/math/typetraits/HasSIMDExp.h>
#include <blaze/math/typetraits/HasSIMDFloor.h>
#include <blaze/math/typetraits/HasSIMDInvCbrt.h>
#include <blaze/math/typetraits/HasSIMDInvErf.h>
#include <blaze/math/typetraits/HasSIMDInvErfc.h>
#include <blaze/math/typetraits/HasSIMDInvSqrt.h>
#include <blaze/math/typetraits/HasSIMDLog.h>
#include <blaze/math/typetraits/HasSIMDLog10.h>
#include <blaze/math/typetraits/HasSIMDMult.h>
#include <blaze/math/typetraits/HasSIMDPow.h>
#include <blaze/math/typetraits/HasSIMDSin.h>
#include <blaze/math/typetraits/HasSIMDSinh.h>
#include <blaze/math/typetraits/HasSIMDSqrt.h>
#include <blaze/math/typetraits/HasSIMDSub.h>
#include <blaze/math/typetraits/HasSIMDTan.h>
#include <blaze/math/typetraits/HasSIMDTanh.h>
#include <blaze/util/constraints/Numeric.h>
#include <blaze/util/FalseType.h>
#include <blaze/util/Memory.h>
#include <blaze/util/NonCopyable.h>
#include <blaze/util/policies/Deallocate.h>
#include <blaze/util/Random.h>
#include <blaze/util/TrueType.h>


namespace blazetest {

namespace mathtest {

namespace simd {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Auxiliary class template for the SIMD operation test.
//
// This class template represents the tests of all available SIMD operations for the given
// numeric data type \a T. In these tests both aligned and unaligned load/store operations
// are used.
*/
template< typename T >  // Data type of the SIMD test
class OperationTest : private blaze::NonCopyable
{
 private:
   //**Type definitions****************************************************************************
   typedef blaze::SIMDTrait<T>  SIMD;      //!< SIMD trait for the given numeric type.
   typedef typename SIMD::Type  SIMDType;  //!< SIMD type for the given numeric type.
   //**********************************************************************************************

   //**********************************************************************************************
   enum : size_t { SIMDSIZE = SIMD::size };  //!< Number of elements in a single SIMD vector.
   enum : size_t { N = 256UL };              //!< Number of numeric values to be worked on.
   enum : size_t { NN = N + SIMDSIZE };      //!< Total number of numeric values in each array.
   //**********************************************************************************************

 public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit OperationTest();
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~OperationTest();
   //@}
   //**********************************************************************************************

 private:
   //**Test functions******************************************************************************
   /*!\name Test functions */
   //@{
   void testStorea        ();
   void testStream        ();
   void testStoreu        ( size_t offset );

   void testAddition      ( blaze::TrueType );
   void testAddition      ( blaze::FalseType );
   void testSubtraction   ( blaze::TrueType );
   void testSubtraction   ( blaze::FalseType );
   void testMultiplication( blaze::TrueType );
   void testMultiplication( blaze::FalseType );
   void testDivision      ( blaze::TrueType );
   void testDivision      ( blaze::FalseType );

   void testAbs           ( blaze::TrueType );
   void testAbs           ( blaze::FalseType );
   void testConj          ( blaze::TrueType );
   void testConj          ( blaze::FalseType );
   void testSqrt          ( blaze::TrueType );
   void testSqrt          ( blaze::FalseType );
   void testInvSqrt       ( blaze::TrueType );
   void testInvSqrt       ( blaze::FalseType );
   void testCbrt          ( blaze::TrueType );
   void testCbrt          ( blaze::FalseType );
   void testInvCbrt       ( blaze::TrueType );
   void testInvCbrt       ( blaze::FalseType );
   void testFloor         ( blaze::TrueType );
   void testFloor         ( blaze::FalseType );
   void testCeil          ( blaze::TrueType );
   void testCeil          ( blaze::FalseType );

   void testPow           ( blaze::TrueType );
   void testPow           ( blaze::FalseType );
   void testExp           ( blaze::TrueType );
   void testExp           ( blaze::FalseType );
   void testLog           ( blaze::TrueType );
   void testLog           ( blaze::FalseType );
   void testLog10         ( blaze::TrueType );
   void testLog10         ( blaze::FalseType );

   void testSin           ( blaze::TrueType );
   void testSin           ( blaze::FalseType );
   void testAsin          ( blaze::TrueType );
   void testAsin          ( blaze::FalseType );
   void testSinh          ( blaze::TrueType );
   void testSinh          ( blaze::FalseType );
   void testAsinh         ( blaze::TrueType );
   void testAsinh         ( blaze::FalseType );

   void testCos           ( blaze::TrueType );
   void testCos           ( blaze::FalseType );
   void testAcos          ( blaze::TrueType );
   void testAcos          ( blaze::FalseType );
   void testCosh          ( blaze::TrueType );
   void testCosh          ( blaze::FalseType );
   void testAcosh         ( blaze::TrueType );
   void testAcosh         ( blaze::FalseType );

   void testTan           ( blaze::TrueType );
   void testTan           ( blaze::FalseType );
   void testAtan          ( blaze::TrueType );
   void testAtan          ( blaze::FalseType );
   void testTanh          ( blaze::TrueType );
   void testTanh          ( blaze::FalseType );
   void testAtanh         ( blaze::TrueType );
   void testAtanh         ( blaze::FalseType );

   void testErf           ( blaze::TrueType );
   void testErf           ( blaze::FalseType );
   void testInvErf        ( blaze::TrueType );
   void testInvErf        ( blaze::FalseType );
   void testErfc          ( blaze::TrueType );
   void testErfc          ( blaze::FalseType );
   void testInvErfc       ( blaze::TrueType );
   void testInvErfc       ( blaze::FalseType );

   void testReduction     ();
   //@}
   //**********************************************************************************************

   //**Error detection functions*******************************************************************
   /*!\name Error detection functions */
   //@{
   void compare( const T* expected, const T* actual ) const;
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
   T* c_;  //!< The third aligned array of size NN.
   T* d_;  //!< The fourth aligned array of size NN.

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
/*!\brief Constructor for the SIMD operation test.
//
// \exception std::runtime_error Operation error detected.
*/
template< typename T >  // Data type of the SIMD test
OperationTest<T>::OperationTest()
   : a_    ( blaze::allocate<T>( NN ) )  // The first aligned array of size NN
   , b_    ( blaze::allocate<T>( NN ) )  // The second aligned array of size NN
   , c_    ( blaze::allocate<T>( NN ) )  // The third aligned array of size NN
   , d_    ( blaze::allocate<T>( NN ) )  // The fourth aligned array of size NN
   , test_ ()                            // Label of the currently performed test
{
   testStorea();
   testStream();

   for( size_t offset=0UL; offset<SIMDSIZE; ++offset ) {
      testStoreu( offset );
   }

   testAddition      ( typename blaze::HasSIMDAdd    <T,T>::Type() );
   testSubtraction   ( typename blaze::HasSIMDSub    <T,T>::Type() );
   testMultiplication( typename blaze::HasSIMDMult   <T,T>::Type() );
   testDivision      ( typename blaze::HasSIMDDiv    <T,T>::Type() );

   testAbs           ( typename blaze::HasSIMDAbs    < T >::Type() );
   testConj          ( typename blaze::HasSIMDConj   < T >::Type() );
   testSqrt          ( typename blaze::HasSIMDSqrt   < T >::Type() );
   testInvSqrt       ( typename blaze::HasSIMDInvSqrt< T >::Type() );
   testCbrt          ( typename blaze::HasSIMDCbrt   < T >::Type() );
   testInvCbrt       ( typename blaze::HasSIMDInvCbrt< T >::Type() );
   testFloor         ( typename blaze::HasSIMDFloor  < T >::Type() );
   testCeil          ( typename blaze::HasSIMDCeil   < T >::Type() );

   testPow           ( typename blaze::HasSIMDPow    < T >::Type() );
   testExp           ( typename blaze::HasSIMDExp    < T >::Type() );
   testLog           ( typename blaze::HasSIMDLog    < T >::Type() );
   testLog10         ( typename blaze::HasSIMDLog10  < T >::Type() );

   testSin           ( typename blaze::HasSIMDSin    < T >::Type() );
   testAsin          ( typename blaze::HasSIMDAsin   < T >::Type() );
   testSinh          ( typename blaze::HasSIMDSinh   < T >::Type() );
   testAsinh         ( typename blaze::HasSIMDAsinh  < T >::Type() );

   testCos           ( typename blaze::HasSIMDCos    < T >::Type() );
   testAcos          ( typename blaze::HasSIMDAcos   < T >::Type() );
   testCosh          ( typename blaze::HasSIMDCosh   < T >::Type() );
   testAcosh         ( typename blaze::HasSIMDAcosh  < T >::Type() );

   testTan           ( typename blaze::HasSIMDTan    < T >::Type() );
   testAtan          ( typename blaze::HasSIMDAtan   < T >::Type() );
   testTanh          ( typename blaze::HasSIMDTanh   < T >::Type() );
   testAtanh         ( typename blaze::HasSIMDAtanh  < T >::Type() );

   testErf           ( typename blaze::HasSIMDErf    < T >::Type() );
   testInvErf        ( typename blaze::HasSIMDInvErf < T >::Type() );
   testErfc          ( typename blaze::HasSIMDErfc   < T >::Type() );
   testInvErfc       ( typename blaze::HasSIMDInvErfc< T >::Type() );

   testReduction     ();
}
//*************************************************************************************************




//=================================================================================================
//
//  DESTRUCTOR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief The destructor for the SIMD operation test.
*/
template< typename T >  // Data type of the SIMD test
OperationTest<T>::~OperationTest()
{
   blaze::deallocate( a_ );
   blaze::deallocate( b_ );
   blaze::deallocate( c_ );
   blaze::deallocate( d_ );
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
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testStorea()
{
   using blaze::loada;
   using blaze::storea;

   test_  = "storea() operation";

   initialize();

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( b_+i, loada( a_+i ) );
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
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testStream()
{
   using blaze::loada;
   using blaze::stream;

   test_  = "stream() operation";

   initialize();

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      stream( b_+i, loada( a_+i ) );
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
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testStoreu( size_t offset )
{
   using blaze::loadu;
   using blaze::storeu;

   test_  = "storeu() operation";

   initialize();

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storeu( b_+offset+i, loadu( a_+offset+i ) );
   }

   compare( a_+offset, b_+offset );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the addition operation.
//
// \return void
// \exception std::runtime_error Addition error detected.
//
// This function tests the addition operation by comparing the results of a vectorized and a
// scalar addition. In case any error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testAddition( blaze::TrueType )
{
   using blaze::loada;
   using blaze::storea;

   test_ = "Addition operation";

   initialize();

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = a_[i] + b_[i];
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, loada( a_+i ) + loada( b_+i ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the addition operation.
//
// \return void
//
// This function is called in case the addition operation is not available for the given data
// type \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testAddition( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the subtraction operation.
//
// \return void
// \exception std::runtime_error Subtraction error detected.
//
// This function tests the subtraction operation by comparing the results of a vectorized and a
// scalar subtraction. In case any error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testSubtraction( blaze::TrueType )
{
   using blaze::loada;
   using blaze::storea;

   test_ = "Subtraction operation";

   initialize();

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = a_[i] - b_[i];
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, loada( a_+i ) - loada( b_+i ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the subtraction operation.
//
// \return void
//
// This function is called in case the subtraction operation is not available for the given data
// type \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testSubtraction( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the multiplication operation.
//
// \return void
// \exception std::runtime_error Multiplication error detected.
//
// This function tests the multiplication operation by comparing the results of a vectorized and
// a scalar multiplication. In case any error is detected, a \a std::runtime_error exception is
// thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testMultiplication( blaze::TrueType )
{
   using blaze::loada;
   using blaze::storea;

   test_ = "Multiplication operation";

   initialize();

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = a_[i] * b_[i];
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, loada( a_+i ) * loada( b_+i ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the multiplication operation.
//
// \return void
//
// This function is called in case the multiplication operation is not available for the given
// data type \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testMultiplication( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the division operation.
//
// \return void
// \exception std::runtime_error Division error detected.
//
// This function tests the division operation by comparing the results of a vectorized and a
// scalar division. In case any error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testDivision( blaze::TrueType )
{
   using blaze::loada;
   using blaze::storea;

   test_ = "Division operation";

   initialize();

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = a_[i] / b_[i];
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, loada( a_+i ) / loada( b_+i ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the division operation.
//
// \return void
//
// This function is called in case the division operation is not available for the given data
// type \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testDivision( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the absolute value operation.
//
// \return void
// \exception std::runtime_error Error in absolute value computation detected.
//
// This function tests the absolute value operation by comparing the results of a vectorized
// and a scalar absolute value operation. In case any error is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testAbs( blaze::TrueType )
{
   using std::abs;
   using blaze::abs;
   using blaze::loada;
   using blaze::storea;

   test_ = "Absolute value operation";

   initialize();

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = abs( a_[i] );
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, abs( loada( a_+i ) ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the absolute value operation.
//
// \return void
//
// This function is called in case the absolute value operation is not available for the given
// data type \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testAbs( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the conjugate operation.
//
// \return void
// \exception std::runtime_error Error in conjugate computation detected.
//
// This function tests the conjugate operation by comparing the results of a vectorized and a
// scalar conjugate operation. In case any error is detected, a \a std::runtime_error exception
// is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testConj( blaze::TrueType )
{
   using blaze::conj;
   using blaze::loada;
   using blaze::storea;

   test_ = "Conjugate operation";

   initialize();

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = conj( a_[i] );
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, conj( loada( a_+i ) ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the conjugate operation.
//
// \return void
//
// This function is called in case the conjugate operation is not available for the given data
// type \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testConj( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the square root operation.
//
// \return void
// \exception std::runtime_error Error in square root computation detected.
//
// This function tests the square root operation by comparing the results of a vectorized and
// a scalar square root operation. In case any error is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testSqrt( blaze::TrueType )
{
   using std::sqrt;
   using blaze::loada;
   using blaze::storea;

   test_ = "Square root operation";

   initialize();

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = sqrt( a_[i] );
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, sqrt( loada( a_+i ) ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the square root operation.
//
// \return void
//
// This function is called in case the square root operation is not available for the given data
// type \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testSqrt( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the inverse square root operation.
//
// \return void
// \exception std::runtime_error Error in inverse square root computation detected.
//
// This function tests the inverse square root operation by comparing the results of a vectorized
// and a scalar inverse square root operation. In case any error is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testInvSqrt( blaze::TrueType )
{
   using std::sqrt;
   using blaze::loada;
   using blaze::storea;

   test_ = "Square root operation";

   initialize();

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = T(1) / sqrt( a_[i] );
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, invsqrt( loada( a_+i ) ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the inverse square root operation.
//
// \return void
//
// This function is called in case the inverse square root operation is not available for the
// given data type \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testInvSqrt( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the cubic root operation.
//
// \return void
// \exception std::runtime_error Error in cubic root computation detected.
//
// This function tests the cubic root operation by comparing the results of a vectorized and
// a scalar cubic root operation. In case any error is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testCbrt( blaze::TrueType )
{
   using std::cbrt;
   using blaze::loada;
   using blaze::storea;

   test_ = "Cubic root operation";

   initialize();

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = cbrt( a_[i] );
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, cbrt( loada( a_+i ) ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the cubic root operation.
//
// \return void
//
// This function is called in case the cubic root operation is not available for the given data
// type \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testCbrt( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the inverse cubic root operation.
//
// \return void
// \exception std::runtime_error Error in inverse cubic root computation detected.
//
// This function tests the inverse cubic root operation by comparing the results of a vectorized
// and a scalar inverse cubic root operation. In case any error is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testInvCbrt( blaze::TrueType )
{
   using std::cbrt;
   using blaze::loada;
   using blaze::storea;

   test_ = "Cubic root operation";

   initialize();

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = T(1) / cbrt( a_[i] );
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, invcbrt( loada( a_+i ) ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the inverse cubic root operation.
//
// \return void
//
// This function is called in case the inverse cubic root operation is not available for the
// given data type \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testInvCbrt( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the floor operation.
//
// \return void
// \exception std::runtime_error Error in floor computation detected.
//
// This function tests the floor operation by comparing the results of a vectorized and a scalar
// floor operation. In case any error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testFloor( blaze::TrueType )
{
   using std::floor;
   using blaze::loada;
   using blaze::storea;

   test_ = "Floor operation";

   initialize();

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = floor( a_[i] );
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, floor( loada( a_+i ) ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the floor operation.
//
// \return void
//
// This function is called in case the floor operation is not available for the given data type
// \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testFloor( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the ceil operation.
//
// \return void
// \exception std::runtime_error Error in ceil computation detected.
//
// This function tests the ceil operation by comparing the results of a vectorized and a scalar
// ceil operation. In case any error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testCeil( blaze::TrueType )
{
   using std::ceil;
   using blaze::loada;
   using blaze::storea;

   test_ = "Ceil operation";

   initialize();

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = ceil( a_[i] );
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, ceil( loada( a_+i ) ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the ceil operation.
//
// \return void
//
// This function is called in case the ceil operation is not available for the given data type
// \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testCeil( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the power operation.
//
// \return void
// \exception std::runtime_error Error in power computation detected.
//
// This function tests the power operation by comparing the results of a vectorized and a
// scalar power operation. In case any error is detected, a \a std::runtime_error exception
// is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testPow( blaze::TrueType )
{
   using std::pow;
   using blaze::loada;
   using blaze::storea;

   test_ = "Power operation";

   initialize();

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = pow( a_[i], b_[i] );
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, pow( loada( a_+i ), loada( b_+i ) ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the power operation.
//
// \return void
//
// This function is called in case the power operation is not available for the given data type
// \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testPow( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the exponent operation.
//
// \return void
// \exception std::runtime_error Error in exponent computation detected.
//
// This function tests the exponent operation by comparing the results of a vectorized and a
// scalar exponent operation. In case any error is detected, a \a std::runtime_error exception
// is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testExp( blaze::TrueType )
{
   using std::exp;
   using blaze::loada;
   using blaze::storea;

   test_ = "Exponent operation";

   initialize();

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = exp( a_[i] );
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, exp( loada( a_+i ) ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the exp operation.
//
// \return void
//
// This function is called in case the exp operation is not available for the given data type
// \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testExp( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the natural logarithm operation.
//
// \return void
// \exception std::runtime_error Error in natural logarithm computation detected.
//
// This function tests the natural logarithm operation by comparing the results of a vectorized
// and a scalar exponent operation. In case any error is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testLog( blaze::TrueType )
{
   using std::log;
   using blaze::loada;
   using blaze::storea;

   test_ = "Natural logarithm operation";

   initialize();

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = log( a_[i] );
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, log( loada( a_+i ) ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the natural logarithm operation.
//
// \return void
//
// This function is called in case the natural logarithm operation is not available for the
// given data type \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testLog( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the common logarithm operation.
//
// \return void
// \exception std::runtime_error Error in common logarithm computation detected.
//
// This function tests the common logarithm operation by comparing the results of a vectorized
// and a scalar exponent operation. In case any error is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testLog10( blaze::TrueType )
{
   using std::log10;
   using blaze::loada;
   using blaze::storea;

   test_ = "Common logarithm operation";

   initialize();

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = log10( a_[i] );
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, log10( loada( a_+i ) ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the common logarithm operation.
//
// \return void
//
// This function is called in case the common logarithm operation is not available for the
// given data type \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testLog10( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the sine operation.
//
// \return void
// \exception std::runtime_error Error in sine computation detected.
//
// This function tests the sine operation by comparing the results of a vectorized and a scalar
// sine operation. In case any error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testSin( blaze::TrueType )
{
   using std::sin;
   using blaze::loada;
   using blaze::storea;

   test_ = "Sine operation";

   initialize();

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = sin( a_[i] );
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, sin( loada( a_+i ) ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the sine operation.
//
// \return void
//
// This function is called in case the sine operation is not available for the given data type
// \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testSin( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the inverse sine operation.
//
// \return void
// \exception std::runtime_error Error in inverse sine computation detected.
//
// This function tests the inverse sine operation by comparing the results of a vectorized and a
// scalar inverse sine operation. In case any error is detected, a \a std::runtime_error exception
// is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testAsin( blaze::TrueType )
{
   using std::asin;
   using blaze::loada;
   using blaze::storea;

   test_ = "Inverse sine operation";

   initialize();

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = asin( a_[i] );
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, asin( loada( a_+i ) ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the inverse sine operation.
//
// \return void
//
// This function is called in case the inverse sine operation is not available for the given data
// type \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testAsin( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the hyperbolic sine operation.
//
// \return void
// \exception std::runtime_error Error in hyperbolic sine computation detected.
//
// This function tests the hyperbolic sine operation by comparing the results of a vectorized
// and a scalar hyperbolic sine operation. In case any error is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testSinh( blaze::TrueType )
{
   using std::sinh;
   using blaze::loada;
   using blaze::storea;

   test_ = "Hyperbolic sine operation";

   initialize();

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = sinh( a_[i] );
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, sinh( loada( a_+i ) ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the hyperbolic sine operation.
//
// \return void
//
// This function is called in case the hyperbolic sine operation is not available for the given
// data type \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testSinh( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the inverse hyperbolic sine operation.
//
// \return void
// \exception std::runtime_error Error in inverse hyperbolic sine computation detected.
//
// This function tests the inverse hyperbolic sine operation by comparing the results of
// a vectorized and a scalar inverse hyperbolic sine operation. In case any error is detected, a
// \a std::runtime_error exception is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testAsinh( blaze::TrueType )
{
   using std::asinh;
   using blaze::loada;
   using blaze::storea;

   test_ = "Inverse hyperbolic sine operation";

   initialize();

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = asinh( a_[i] );
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, asinh( loada( a_+i ) ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the inverse hyperbolic sine operation.
//
// \return void
//
// This function is called in case the inverse hyperbolic sine operation is not available for
// the given data type \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testAsinh( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the cosine operation.
//
// \return void
// \exception std::runtime_error Error in cosine computation detected.
//
// This function tests the cosine operation by comparing the results of a vectorized and a scalar
// cosine operation. In case any error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testCos( blaze::TrueType )
{
   using std::cos;
   using blaze::loada;
   using blaze::storea;

   test_ = "Cosine operation";

   initialize();

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = cos( a_[i] );
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, cos( loada( a_+i ) ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the cosine operation.
//
// \return void
//
// This function is called in case the cosine operation is not available for the given data type
// \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testCos( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the inverse cosine operation.
//
// \return void
// \exception std::runtime_error Error in inverse cosine computation detected.
//
// This function tests the inverse cosine operation by comparing the results of a vectorized and
// a scalar inverse cosine operation. In case any error is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testAcos( blaze::TrueType )
{
   using std::acos;
   using blaze::loada;
   using blaze::storea;

   test_ = "Inverse cosine operation";

   initialize();

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = acos( a_[i] );
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, acos( loada( a_+i ) ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the inverse cosine operation.
//
// \return void
//
// This function is called in case the inverse cosine operation is not available for the given
// data type \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testAcos( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the hyperbolic cosine operation.
//
// \return void
// \exception std::runtime_error Error in hyperbolic cosine computation detected.
//
// This function tests the hyperbolic cosine operation by comparing the results of a vectorized
// and a scalar hyperbolic cosine operation. In case any error is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testCosh( blaze::TrueType )
{
   using std::cosh;
   using blaze::loada;
   using blaze::storea;

   test_ = "Hyperbolic cosine operation";

   initialize();

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = cosh( a_[i] );
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, cosh( loada( a_+i ) ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the hyperbolic cosine operation.
//
// \return void
//
// This function is called in case the hyperbolic cosine operation is not available for the given
// data type \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testCosh( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the inverse hyperbolic cosine operation.
//
// \return void
// \exception std::runtime_error Error in inverse hyperbolic cosine computation detected.
//
// This function tests the inverse hyperbolic cosine operation by comparing the results of
// a vectorized and a scalar inverse hyperbolic cosine operation. In case any error is detected,
// a \a std::runtime_error exception is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testAcosh( blaze::TrueType )
{
   using std::acosh;
   using blaze::loada;
   using blaze::storea;

   test_ = "Inverse hyperbolic cosine operation";

   initialize();

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = acosh( a_[i] );
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, acosh( loada( a_+i ) ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the inverse hyperbolic cosine operation.
//
// \return void
//
// This function is called in case the inverse hyperbolic cosine operation is not available for
// the given data type \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testAcosh( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the tangent operation.
//
// \return void
// \exception std::runtime_error Error in tangent computation detected.
//
// This function tests the tangent operation by comparing the results of a vectorized and a scalar
// tangent operation. In case any error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testTan( blaze::TrueType )
{
   using std::tan;
   using blaze::loada;
   using blaze::storea;

   test_ = "Tangent operation";

   initialize();

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = tan( a_[i] );
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, tan( loada( a_+i ) ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the tangent operation.
//
// \return void
//
// This function is called in case the tangent operation is not available for the given data type
// \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testTan( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the inverse tangent operation.
//
// \return void
// \exception std::runtime_error Error in inverse tangent computation detected.
//
// This function tests the inverse tangent operation by comparing the results of a vectorized and
// a scalar inverse tangent operation. In case any error is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testAtan( blaze::TrueType )
{
   using std::atan;
   using blaze::loada;
   using blaze::storea;

   test_ = "Inverse tangent operation";

   initialize();

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = atan( a_[i] );
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, atan( loada( a_+i ) ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the inverse tangent operation.
//
// \return void
//
// This function is called in case the inverse tangent operation is not available for the given
// data type \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testAtan( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the hyperbolic tangent operation.
//
// \return void
// \exception std::runtime_error Error in hyperbolic tangent computation detected.
//
// This function tests the hyperbolic tangent operation by comparing the results of a vectorized
// and a scalar hyperbolic tangent operation. In case any error is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testTanh( blaze::TrueType )
{
   using std::tanh;
   using blaze::loada;
   using blaze::storea;

   test_ = "Hyperbolic tangent operation";

   initialize();

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = tanh( a_[i] );
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, tanh( loada( a_+i ) ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the hyperbolic tangent operation.
//
// \return void
//
// This function is called in case the hyperbolic tangent operation is not available for the
// given data type \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testTanh( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the inverse hyperbolic tangent operation.
//
// \return void
// \exception std::runtime_error Error in inverse hyperbolic tangent computation detected.
//
// This function tests the inverse hyperbolic tangent operation by comparing the results of
// a vectorized and a scalar inverse hyperbolic tangent operation. In case any error is detected,
// a \a std::runtime_error exception is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testAtanh( blaze::TrueType )
{
   using std::atanh;
   using blaze::loada;
   using blaze::storea;

   test_ = "Inverse hyperbolic tangent operation";

   initialize();

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = atanh( a_[i] );
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, atanh( loada( a_+i ) ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the inverse hyperbolic tangent operation.
//
// \return void
//
// This function is called in case the inverse hyperbolic tangent operation is not available for
// the given data type \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testAtanh( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the error function (\c erf).
//
// \return void
// \exception std::runtime_error Error in error function computation detected.
//
// This function tests the error function (\c erf) by comparing the results of a vectorized and
// a scalar error function. In case any error is detected, a \a std::runtime_error exception is
// thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testErf( blaze::TrueType )
{
   using std::erf;
   using blaze::loada;
   using blaze::storea;

   test_ = "Error function";

   initialize();

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = erf( a_[i] );
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, erf( loada( a_+i ) ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the error function (\c erf).
//
// \return void
//
// This function is called in case the error function is not available for the given data type
// \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testErf( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the inverse error function (\c inverf).
//
// \return void
// \exception std::runtime_error Error in inverse error function computation detected.
//
// This function tests the inverse error function (\c inverf) by comparing the results of a
// vectorized and a scalar error function. In case any error is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testInvErf( blaze::TrueType )
{
   using std::erf;
   using blaze::loada;
   using blaze::storea;

   test_ = "Error function";

   initialize();

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = T(1) / erf( a_[i] );
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, inverf( loada( a_+i ) ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the inverse error function (\c inverf).
//
// \return void
//
// This function is called in case the inverse error function is not available for the given
// data type \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testInvErf( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the complementary error function (\c erfc).
//
// \return void
// \exception std::runtime_error Error in complementary error function computation detected.
//
// This function tests the complementary error function (\c erfc) by comparing the results of a
// vectorized and a scalar error function. In case any error is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testErfc( blaze::TrueType )
{
   using std::erfc;
   using blaze::loada;
   using blaze::storea;

   test_ = "Complementary error function";

   initialize();

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = erfc( a_[i] );
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, erfc( loada( a_+i ) ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the complementary error function (\c erf).
//
// \return void
//
// This function is called in case the complementary error function is not available for the
// given data type \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testErfc( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the inverse complementary error function (\c inverfc).
//
// \return void
// \exception std::runtime_error Error in inverse complementary error function computation detected.
//
// This function tests the inverse complementary error function (\c erfc) by comparing the
// results of a vectorized and a scalar error function. In case any error is detected, a
// \a std::runtime_error exception is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testInvErfc( blaze::TrueType )
{
   using std::erfc;
   using blaze::loada;
   using blaze::storea;

   test_ = "Inverse complementary error function";

   initialize();

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = T(1) / erfc( a_[i] );
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, inverfc( loada( a_+i ) ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the inverse complementary error function (\c inverf).
//
// \return void
//
// This function is called in case the inverse complementary error function is not available
// for the given data type \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testInvErfc( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the reduction operation.
//
// \return void
// \exception std::runtime_error Error in reduction computation detected.
//
// This function tests the reduction operation by comparing the results of a vectorized and a
// scalar reduction. In case any error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testReduction()
{
   using blaze::loada;
   using blaze::sum;

   test_ = "sum() operation";

   initialize();

   T ssum = T();
   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      T tmp = T();
      for( size_t j=0UL; j<SIMDSIZE; ++j ) {
         tmp += a_[i+j];
      }
      ssum += tmp;
   }

   T vsum = T();
   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      vsum += sum( loada( a_+i ) );
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
// \param expected The array of expected values.
// \param actual The array of actual values.
// \return void
// \exception std::runtime_error Value mismatch detected.
//
// This function compares the first 256 elements of the two given arrays. In case any value of
// the two arrays differs, a \a std::runtime_error exception is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::compare( const T* expected, const T* actual ) const
{
   for( size_t i=0UL; i<N; ++i ) {
      if( expected[i] != actual[i] ) {
         std::ostringstream oss;
         oss.precision( 20 );
         oss << " Test : " << test_ << "\n"
             << " Error: Value mismatch detected at index " << i << "\n"
             << " Details:\n"
             << "   expected[" << i << "] = " << expected[i] << "\n"
             << "   actual  [" << i << "] = " << actual[i] << "\n";
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
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::initialize()
{
   using blaze::randomize;

   for( size_t i=0UL; i<NN; ++i ) {
      randomize( a_[i] );
      randomize( b_[i] );
      randomize( c_[i] );
      randomize( d_[i] );
   }
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Testing the SIMD operations of a specific numeric data type.
//
// \return void
*/
template< typename T >  // Data type of the SIMD test
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
/*!\brief Macro for the execution of an SIMD operation test case.
*/
#define RUN_SIMD_OPERATION_TEST( T ) \
   blazetest::mathtest::simd::runTest<T>()
/*! \endcond */
//*************************************************************************************************

} // namespace simd

} // namespace mathtest

} // namespace blazetest

#endif
