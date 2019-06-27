//=================================================================================================
/*!
//  \file blazetest/mathtest/lapack/EigenvalueTest.h
//  \brief Header file for the LAPACK eigenvalue test
//
//  Copyright (C) 2012-2019 Klaus Iglberger - All Rights Reserved
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

#ifndef _BLAZETEST_MATHTEST_LAPACK_EIGENVALUETEST_H_
#define _BLAZETEST_MATHTEST_LAPACK_EIGENVALUETEST_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <algorithm>
#include <sstream>
#include <stdexcept>
#include <string>
#include <typeinfo>
#include <blaze/math/Aliases.h>
#include <blaze/math/Column.h>
#include <blaze/math/HermitianMatrix.h>
#include <blaze/math/LAPACK.h>
#include <blaze/math/Row.h>
#include <blaze/math/StaticMatrix.h>
#include <blaze/math/StaticVector.h>
#include <blaze/math/SymmetricMatrix.h>
#include <blaze/math/typetraits/UnderlyingElement.h>
#include <blazetest/system/LAPACK.h>
#include <blaze/util/Complex.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/Random.h>
#include <blaze/util/typetraits/IsComplex.h>


namespace blazetest {

namespace mathtest {

namespace lapack {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Auxiliary class for all tests of the LAPACK functionality.
//
// This class represents a test suite for LAPACK functionality wrapped by the Blaze library.
*/
class EigenvalueTest
{
 public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit EigenvalueTest();
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
   template< typename Type > void testGeev();
   template< typename Type, bool SOA, bool SOB, bool SOL, bool SOR > void testGges();
   template< typename Type, bool SOA, bool SOB, bool SOL, bool SOR > void testGgesSelect();
   template< typename Type > void testSyev();
   template< typename Type > void testSyevd();
   template< typename Type > void testSyevx();
   template< typename Type > void testHeev();
   template< typename Type > void testHeevd();
   template< typename Type > void testHeevx();
   //@}
   //**********************************************************************************************

   //**Error detection functions*******************************************************************
   /*!\name Error detection functions */
   //@{
   template< typename VT, typename MT, bool SO, typename ST >
   void checkEigenvector( const blaze::DenseVector<VT,false>& v,
                          const blaze::DenseMatrix<MT,SO>& A, ST w );

   template< typename VT, typename MT, bool SO, typename ST >
   void checkEigenvector( const blaze::DenseVector<VT,true>& u,
                          const blaze::DenseMatrix<MT,SO>& A, ST w );

   template< typename MT1, bool SO1, typename MT2, bool SO2, typename ST1, typename ST2 >
   void checkEigenvalue( const blaze::DenseMatrix<MT1, SO1>& A,
                         const blaze::DenseMatrix<MT2, SO2>& B, 
                         ST1 alpha, ST2 beta );
   //@}
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   std::string test_;  //!< Label of the currently performed test.
   //@}
   //**********************************************************************************************
};
//*************************************************************************************************


//=================================================================================================
//
//  HELPER FUNCTIONS / CLASSES
//
//=================================================================================================
namespace detail {

template < bool TF >
struct ConditionalTranspose;

template <>
struct ConditionalTranspose<true>
{
   template < typename MT, bool SO >
   decltype(auto) operator()( blaze::Matrix<MT,SO> const& A ) const
   {
      return trans( ~A );
   }   
};


template <>
struct ConditionalTranspose<false>
{
   template < typename MT, bool SO >
   MT const& operator()( blaze::Matrix<MT,SO> const& A ) const
   {
      return ~A;
   }   
};


template < typename ST >
inline ST singularEigenvalueThreshold();


template <>
inline double singularEigenvalueThreshold<double>()
{
   return 1e-12;
}


template <>
inline float singularEigenvalueThreshold<float>()
{
   return 1e-4f;
}


template < typename ST >
inline ST schurFactorizationTolerance();


template <>
inline double schurFactorizationTolerance<double>()
{
   return 1e-14;
}


template <>
inline float schurFactorizationTolerance<float>()
{
   return 1e-5f;
}

}


//=================================================================================================
//
//  TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Test of the eigenvalue functions for general matrices (geev).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the eigenvalue functions for general matrices for various
// data types. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void EigenvalueTest::testGeev()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   test_ = "General matrix eigenvalue computation (geev)";

   using CT = blaze::If_t< blaze::IsComplex_v<Type>, Type, blaze::complex<Type> >;

   const auto comparator = []( const CT& c1, const CT& c2 ) {
      return blaze::equal( c1, c2 );
   };

   {
      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> A1;
      randomize( A1 );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> A2( A1 );

      blaze::StaticVector<CT,3UL,blaze::rowVector> w1;
      blaze::StaticVector<CT,3UL,blaze::rowVector> w2;

      blaze::geev( A1, w1 );
      blaze::geev( A2, w2 );

      if( !std::is_permutation( w1.begin(), w1.end(), w2.begin(), comparator ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: General matrix eigenvalue computation failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major eigenvalues:\n" << w1 << "\n"
             << "   Column-major eigenvalues:\n" << w2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> A1( A );
      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> A2( A );

      blaze::StaticMatrix<CT,3UL,3UL,blaze::rowMajor> VL1;
      blaze::StaticMatrix<CT,3UL,3UL,blaze::columnMajor> VL2;

      blaze::StaticVector<CT,3UL,blaze::rowVector> w1;
      blaze::StaticVector<CT,3UL,blaze::rowVector> w2;

      blaze::geev( A1, VL1, w1 );
      blaze::geev( A2, VL2, w2 );

      if( !std::is_permutation( w1.begin(), w1.end(), w2.begin(), comparator ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: General matrix eigenvalue computation failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major eigenvalues:\n" << w1 << "\n"
             << "   Column-major eigenvalues:\n" << w2 << "\n";
         throw std::runtime_error( oss.str() );
      }

      for( size_t i=0UL; i<VL1.rows(); ++i ) {
         checkEigenvector( row( VL1, i ), A, w1[i] );
      }

      for( size_t i=0UL; i<VL2.columns(); ++i ) {
         checkEigenvector( ctrans( column( VL2, i ) ), A, w2[i] );
      }
   }

   {
      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> A1( A );
      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> A2( A );

      blaze::StaticVector<CT,3UL,blaze::rowVector> w1;
      blaze::StaticVector<CT,3UL,blaze::rowVector> w2;

      blaze::StaticMatrix<CT,3UL,3UL,blaze::rowMajor> VR1;
      blaze::StaticMatrix<CT,3UL,3UL,blaze::columnMajor> VR2;

      blaze::geev( A1, w1, VR1 );
      blaze::geev( A2, w2, VR2 );

      if( !std::is_permutation( w1.begin(), w1.end(), w2.begin(), comparator ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: General matrix eigenvalue computation failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major eigenvalues:\n" << w1 << "\n"
             << "   Column-major eigenvalues:\n" << w2 << "\n";
         throw std::runtime_error( oss.str() );
      }

      for( size_t i=0UL; i<VR1.rows(); ++i ) {
         checkEigenvector( ctrans( row( VR1, i ) ), A, w1[i] );
      }

      for( size_t i=0UL; i<VR2.columns(); ++i ) {
         checkEigenvector( column( VR2, i ), A, w2[i] );
      }
   }

   {
      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> A1( A );
      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> A2( A );

      blaze::StaticMatrix<CT,3UL,3UL,blaze::rowMajor> VL1;
      blaze::StaticMatrix<CT,3UL,3UL,blaze::columnMajor> VL2;

      blaze::StaticVector<CT,3UL,blaze::rowVector> w1;
      blaze::StaticVector<CT,3UL,blaze::rowVector> w2;

      blaze::StaticMatrix<CT,3UL,3UL,blaze::rowMajor> VR1;
      blaze::StaticMatrix<CT,3UL,3UL,blaze::columnMajor> VR2;

      blaze::geev( A1, VL1, w1, VR1 );
      blaze::geev( A2, VL2, w2, VR2 );

      if( !std::is_permutation( w1.begin(), w1.end(), w2.begin(), comparator ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: General matrix eigenvalue computation failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major eigenvalues:\n" << w1 << "\n"
             << "   Row-major left eigenvalues:\n" << VL1 << "\n"
             << "   Row-major right eigenvalues:\n" << VR1 << "\n"
             << "   Column-major eigenvalues:\n" << w2 << "\n"
             << "   Column-major left eigenvalues:\n" << VL2 << "\n"
             << "   Column-major right eigenvalues:\n" << VR2 << "\n";
         throw std::runtime_error( oss.str() );
      }

      for( size_t i=0UL; i<VL1.rows(); ++i ) {
         checkEigenvector( row( VL1, i ), A, w1[i] );
      }

      for( size_t i=0UL; i<VR1.rows(); ++i ) {
         checkEigenvector( ctrans( row( VR1, i ) ), A, w1[i] );
      }

      for( size_t i=0UL; i<VL2.columns(); ++i ) {
         checkEigenvector( ctrans( column( VL2, i ) ), A, w2[i] );
      }

      for( size_t i=0UL; i<VR2.columns(); ++i ) {
         checkEigenvector( column( VR2, i ), A, w2[i] );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the eigenvalue functions for symmetric matrices (syev).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the eigenvalue functions for symmetric matrices for various
// data types. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void EigenvalueTest::testSyev()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   test_ = "Symmetric matrix eigenvalue computation (syev)";

   {
      blaze::SymmetricMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > S;
      randomize( S );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor>    A( S );
      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> B( S );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::rowVector> wA;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::rowVector> wB;

      blaze::syev( A, wA, 'N', 'L' );
      blaze::syev( B, wB, 'N', 'U' );

      if( wA != wB ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Symmetric matrix eigenvalue computation failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major eigenvalues:\n" << wA << "\n"
             << "   Column-major eigenvalues:\n" << wB << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      blaze::SymmetricMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > S;
      randomize( S );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor>    A( S );
      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> B( S );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::rowVector> wA;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::rowVector> wB;

      blaze::syev( A, wA, 'V', 'L' );
      blaze::syev( B, wB, 'V', 'L' );

      if( wA != wB ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Symmetric matrix eigenvalue computation failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major eigenvalues:\n" << wA << "\n"
             << "   Column-major eigenvalues:\n" << wB << "\n";
         throw std::runtime_error( oss.str() );
      }

      for( size_t i=0UL; i<A.rows(); ++i ) {
         checkEigenvector( row( A, i ), S, wA[i] );
      }

      for( size_t i=0UL; i<B.columns(); ++i ) {
         checkEigenvector( column( B, i ), S, wB[i] );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the eigenvalue functions for symmetric matrices (syevd).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the eigenvalue functions for symmetric matrices for various
// data types. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void EigenvalueTest::testSyevd()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   test_ = "Symmetric matrix eigenvalue computation (syevd)";

   {
      blaze::SymmetricMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > S;
      randomize( S );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor>    A( S );
      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> B( S );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::rowVector> wA;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::rowVector> wB;

      blaze::syevd( A, wA, 'N', 'L' );
      blaze::syevd( B, wB, 'N', 'U' );

      if( wA != wB ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Symmetric matrix eigenvalue computation failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major eigenvalues:\n" << wA << "\n"
             << "   Column-major eigenvalues:\n" << wB << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      blaze::SymmetricMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > S;
      randomize( S );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor>    A( S );
      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> B( S );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::rowVector> wA;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::rowVector> wB;

      blaze::syevd( A, wA, 'V', 'L' );
      blaze::syevd( B, wB, 'V', 'U' );

      if( wA != wB ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Symmetric matrix eigenvalue computation failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major eigenvalues:\n" << wA << "\n"
             << "   Column-major eigenvalues:\n" << wB << "\n";
         throw std::runtime_error( oss.str() );
      }

      for( size_t i=0UL; i<A.rows(); ++i ) {
         checkEigenvector( row( A, i ), S, wA[i] );
      }

      for( size_t i=0UL; i<B.columns(); ++i ) {
         checkEigenvector( column( B, i ), S, wB[i] );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the eigenvalue functions for symmetric matrices (syevx).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the eigenvalue functions for symmetric matrices for various
// data types. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void EigenvalueTest::testSyevx()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   {
      test_ = "Symmetric matrix eigenvalue computation (syevx)";

      blaze::SymmetricMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > S;
      randomize( S );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor>    A( S );
      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> B( S );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::rowVector> wA;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::rowVector> wB;

      const size_t numA = blaze::syevx( A, wA, 'L' );
      const size_t numB = blaze::syevx( B, wB, 'U' );

      if( numA != numB || wA != wB ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Symmetric matrix eigenvalue computation failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major eigenvalues:\n" << wA << "\n"
             << "   Column-major eigenvalues:\n" << wB << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Symmetric matrix eigenvalue computation (syevx, floating point range)";

      blaze::SymmetricMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > S;
      randomize( S );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor>    A( S );
      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> B( S );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::rowVector> wA;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::rowVector> wB;

      const size_t numA = blaze::syevx( A, wA, 'L', 0.0, 5.0 );
      const size_t numB = blaze::syevx( B, wB, 'U', 0.0, 5.0 );

      if( numA != numB || wA != wB ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Symmetric matrix eigenvalue computation failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major eigenvalues:\n" << wA << "\n"
             << "   Column-major eigenvalues:\n" << wB << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Symmetric matrix eigenvalue computation (syevx, integral range)";

      blaze::SymmetricMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > S;
      randomize( S );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor>    A( S );
      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> B( S );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,2UL,blaze::rowVector> wA;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,2UL,blaze::rowVector> wB;

      const size_t numA = blaze::syevx( A, wA, 'L', 0, 1 );
      const size_t numB = blaze::syevx( B, wB, 'U', 0, 1 );

      if( numA != numB || wA != wB ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Symmetric matrix eigenvalue computation failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major eigenvalues:\n" << wA << "\n"
             << "   Column-major eigenvalues:\n" << wB << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Symmetric matrix eigenvalue computation (syevx)";

      blaze::SymmetricMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > S;
      randomize( S );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor>    A( S );
      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> B( S );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::rowVector> wA;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::rowVector> wB;

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor>    ZA;
      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> ZB;

      const size_t numA = blaze::syevx( A, wA, ZA, 'L' );
      const size_t numB = blaze::syevx( B, wB, ZB, 'U' );

      if( numA != numB || wA != wB ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Symmetric matrix eigenvalue computation failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major eigenvalues:\n" << wA << "\n"
             << "   Column-major eigenvalues:\n" << wB << "\n";
         throw std::runtime_error( oss.str() );
      }

      for( size_t i=0UL; i<ZA.rows(); ++i ) {
         checkEigenvector( row( ZA, i ), S, wA[i] );
      }

      for( size_t i=0UL; i<ZB.columns(); ++i ) {
         checkEigenvector( column( ZB, i ), S, wB[i] );
      }
   }

   {
      test_ = "Symmetric matrix eigenvalue computation (syevx, floating point range)";

      blaze::SymmetricMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > S;
      randomize( S );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor>    A( S );
      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> B( S );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::rowVector> wA;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::rowVector> wB;

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor>    ZA;
      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> ZB;

      const size_t numA = blaze::syevx( A, wA, ZA, 'L', 0.0, 0.5 );
      const size_t numB = blaze::syevx( B, wB, ZB, 'U', 0.0, 0.5 );

      if( numA != numB || wA != wB ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Symmetric matrix eigenvalue computation failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major eigenvalues:\n" << wA << "\n"
             << "   Column-major eigenvalues:\n" << wB << "\n";
         throw std::runtime_error( oss.str() );
      }

      for( size_t i=0UL; i<ZA.rows(); ++i ) {
         checkEigenvector( row( ZA, i ), S, wA[i] );
      }

      for( size_t i=0UL; i<ZB.columns(); ++i ) {
         checkEigenvector( column( ZB, i ), S, wB[i] );
      }
   }

   {
      test_ = "Symmetric matrix eigenvalue computation (syevx, integral range)";

      blaze::SymmetricMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > S;
      randomize( S );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor>    A( S );
      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> B( S );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,2UL,blaze::rowVector> wA;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,2UL,blaze::rowVector> wB;

      blaze::StaticMatrix<Type,2UL,3UL,blaze::rowMajor>    ZA;
      blaze::StaticMatrix<Type,3UL,2UL,blaze::columnMajor> ZB;

      const size_t numA = blaze::syevx( A, wA, ZA, 'L', 0, 1 );
      const size_t numB = blaze::syevx( B, wB, ZB, 'U', 0, 1 );

      if( numA != numB || wA != wB ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Symmetric matrix eigenvalue computation failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major eigenvalues:\n" << wA << "\n"
             << "   Column-major eigenvalues:\n" << wB << "\n";
         throw std::runtime_error( oss.str() );
      }

      for( size_t i=0UL; i<ZA.rows(); ++i ) {
         checkEigenvector( row( ZA, i ), S, wA[i] );
      }

      for( size_t i=0UL; i<ZB.columns(); ++i ) {
         checkEigenvector( column( ZB, i ), S, wB[i] );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the eigenvalue functions for Hermitian matrices (heev).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the eigenvalue functions for Hermitian matrices for various
// data types. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void EigenvalueTest::testHeev()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   test_ = "Hermitian matrix eigenvalue computation (heev)";

   {
      blaze::HermitianMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > H;
      randomize( H );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor>    A( H );
      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> B( H );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::rowVector> wA;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::rowVector> wB;

      blaze::heev( A, wA, 'N', 'L' );
      blaze::heev( B, wB, 'N', 'U' );

      if( wA != wB ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Hermitian matrix eigenvalue computation failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major eigenvalues:\n" << wA << "\n"
             << "   Column-major eigenvalues:\n" << wB << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      blaze::HermitianMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > H;
      randomize( H );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor>    A( H );
      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> B( H );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::rowVector> wA;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::rowVector> wB;

      blaze::heev( A, wA, 'V', 'L' );
      blaze::heev( B, wB, 'V', 'U' );

      if( wA != wB ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Hermitian matrix eigenvalue computation failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major eigenvalues:\n" << wA << "\n"
             << "   Column-major eigenvalues:\n" << wB << "\n";
         throw std::runtime_error( oss.str() );
      }

      for( size_t i=0UL; i<A.rows(); ++i ) {
         checkEigenvector( row( A, i ), H, wA[i] );
      }

      for( size_t i=0UL; i<B.columns(); ++i ) {
         checkEigenvector( column( B, i ), H, wB[i] );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the eigenvalue functions for Hermitian matrices (heevd).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the eigenvalue functions for Hermitian matrices for various
// data types. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void EigenvalueTest::testHeevd()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   test_ = "Hermitian matrix eigenvalue computation (heevd)";

   {
      blaze::HermitianMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > H;
      randomize( H );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor>    A( H );
      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> B( H );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::rowVector> wA;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::rowVector> wB;

      blaze::heevd( A, wA, 'N', 'L' );
      blaze::heevd( B, wB, 'N', 'U' );

      if( wA != wB ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Hermitian matrix eigenvalue computation failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major eigenvalues:\n" << wA << "\n"
             << "   Column-major eigenvalues:\n" << wB << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      blaze::HermitianMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > H;
      randomize( H );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor>    A( H );
      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> B( H );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::rowVector> wA;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::rowVector> wB;

      blaze::heevd( A, wA, 'V', 'L' );
      blaze::heevd( B, wB, 'V', 'U' );

      if( wA != wB ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Hermitian matrix eigenvalue computation failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major eigenvalues:\n" << wA << "\n"
             << "   Column-major eigenvalues:\n" << wB << "\n";
         throw std::runtime_error( oss.str() );
      }

      for( size_t i=0UL; i<A.rows(); ++i ) {
         checkEigenvector( row( A, i ), H, wA[i] );
      }

      for( size_t i=0UL; i<B.columns(); ++i ) {
         checkEigenvector( column( B, i ), H, wB[i] );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the eigenvalue functions for Hermitian matrices (heevx).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the eigenvalue functions for Hermitian matrices for various
// data types. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void EigenvalueTest::testHeevx()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   {
      test_ = "Hermitian matrix eigenvalue computation (heevx)";

      blaze::HermitianMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > H;
      randomize( H );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor>    A( H );
      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> B( H );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::rowVector> wA;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::rowVector> wB;

      const size_t numA = blaze::heevx( A, wA, 'L' );
      const size_t numB = blaze::heevx( B, wB, 'U' );

      if( numA != numB || wA != wB ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Hermitian matrix eigenvalue computation failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major eigenvalues:\n" << wA << "\n"
             << "   Column-major eigenvalues:\n" << wB << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Hermitian matrix eigenvalue computation (heevx, floating point range)";

      blaze::HermitianMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > H;
      randomize( H );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor>    A( H );
      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> B( H );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::rowVector> wA;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::rowVector> wB;

      const size_t numA = blaze::heevx( A, wA, 'L', 0.0, 5.0 );
      const size_t numB = blaze::heevx( B, wB, 'U', 0.0, 5.0 );

      if( numA != numB || wA != wB ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Hermitian matrix eigenvalue computation failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major eigenvalues:\n" << wA << "\n"
             << "   Column-major eigenvalues:\n" << wB << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Hermitian matrix eigenvalue computation (heevx, integral range)";

      blaze::HermitianMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > H;
      randomize( H );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor>    A( H );
      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> B( H );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,2UL,blaze::rowVector> wA;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,2UL,blaze::rowVector> wB;

      const size_t numA = blaze::heevx( A, wA, 'L', 0, 1 );
      const size_t numB = blaze::heevx( B, wB, 'U', 0, 1 );

      if( numA != numB || wA != wB ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Hermitian matrix eigenvalue computation failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major eigenvalues:\n" << wA << "\n"
             << "   Column-major eigenvalues:\n" << wB << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Hermitian matrix eigenvalue computation (heevx)";

      blaze::HermitianMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > H;
      randomize( H );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor>    A( H );
      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> B( H );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::rowVector> wA;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::rowVector> wB;

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor>    ZA;
      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> ZB;

      const size_t numA = blaze::heevx( A, wA, ZA, 'L' );
      const size_t numB = blaze::heevx( B, wB, ZB, 'U' );

      if( numA != numB || wA != wB ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Hermitian matrix eigenvalue computation failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major eigenvalues:\n" << wA << "\n"
             << "   Column-major eigenvalues:\n" << wB << "\n";
         throw std::runtime_error( oss.str() );
      }

      for( size_t i=0UL; i<ZA.rows(); ++i ) {
         checkEigenvector( row( ZA, i ), H, wA[i] );
      }

      for( size_t i=0UL; i<ZB.columns(); ++i ) {
         checkEigenvector( column( ZB, i ), H, wB[i] );
      }
   }

   {
      test_ = "Hermitian matrix eigenvalue computation (heevx, floating point range)";

      blaze::HermitianMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > H;
      randomize( H );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor>    A( H );
      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> B( H );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::rowVector> wA;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::rowVector> wB;

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor>    ZA;
      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> ZB;

      const size_t numA = blaze::heevx( A, wA, ZA, 'L', 0.0, 0.5 );
      const size_t numB = blaze::heevx( B, wB, ZB, 'U', 0.0, 0.5 );

      if( numA != numB || wA != wB ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Hermitian matrix eigenvalue computation failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major eigenvalues:\n" << wA << "\n"
             << "   Column-major eigenvalues:\n" << wB << "\n";
         throw std::runtime_error( oss.str() );
      }

      for( size_t i=0UL; i<ZA.rows(); ++i ) {
         checkEigenvector( row( ZA, i ), H, wA[i] );
      }

      for( size_t i=0UL; i<ZB.columns(); ++i ) {
         checkEigenvector( column( ZB, i ), H, wB[i] );
      }
   }

   {
      test_ = "Hermitian matrix eigenvalue computation (heevx, integral range)";

      blaze::HermitianMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > H;
      randomize( H );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor>    A( H );
      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> B( H );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,2UL,blaze::rowVector> wA;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,2UL,blaze::rowVector> wB;

      blaze::StaticMatrix<Type,2UL,3UL,blaze::rowMajor>    ZA;
      blaze::StaticMatrix<Type,3UL,2UL,blaze::columnMajor> ZB;

      const size_t numA = blaze::heevx( A, wA, ZA, 'L', 0, 1 );
      const size_t numB = blaze::heevx( B, wB, ZB, 'U', 0, 1 );

      if( numA != numB || wA != wB ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Hermitian matrix eigenvalue computation failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major eigenvalues:\n" << wA << "\n"
             << "   Column-major eigenvalues:\n" << wB << "\n";
         throw std::runtime_error( oss.str() );
      }

      for( size_t i=0UL; i<ZA.rows(); ++i ) {
         checkEigenvector( row( ZA, i ), H, wA[i] );
      }

      for( size_t i=0UL; i<ZB.columns(); ++i ) {
         checkEigenvector( column( ZB, i ), H, wB[i] );
      }
   }

#endif
}
//*************************************************************************************************

//*************************************************************************************************
/*!\brief Test of generalized Schur factorization functions for general matrices (gges).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the eigenvalue functions for general matrices for various
// data types. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type, bool SOA, bool SOB, bool SOL, bool SOR >
void EigenvalueTest::testGges()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   test_ = "General matrix eigenvalue and Schur form computation (gges)";

   using CT = blaze::If_t< blaze::IsComplex_v<Type>, Type, blaze::complex<Type> >;
   using detail::ConditionalTranspose;

   {
      blaze::StaticMatrix<Type,3UL,3UL,SOA> A;
      blaze::StaticMatrix<Type,3UL,3UL,SOB> B;
      randomize( A );
      randomize( B );

      blaze::StaticMatrix<Type,3UL,3UL,SOL> VL;
      blaze::StaticMatrix<Type,3UL,3UL,SOR> VR;
      
      blaze::StaticMatrix<Type,3UL,3UL,SOA> S(A);
      blaze::StaticMatrix<Type,3UL,3UL,SOB> T(B);
      blaze::StaticVector<CT,3UL,blaze::rowVector> alpha;
      blaze::StaticVector<Type,3UL,blaze::rowVector> beta;

      blaze::gges( S, T, alpha, beta, VL, VR );

      ConditionalTranspose<SOA == blaze::rowMajor> const CTA;
      ConditionalTranspose<SOB == blaze::rowMajor> const CTB;
      ConditionalTranspose<SOL == blaze::rowMajor> const CTL;
      ConditionalTranspose<SOR == blaze::rowMajor> const CTR;

      for (std::size_t i = 0; i < alpha.size(); ++i)
         checkEigenvalue(CTA(A), CTB(B), alpha[i], beta[i]);

      if( !( maxNorm(CTL(VL) * CTA(S) * CTR(trans(VR)) - CTA(A)) < detail::schurFactorizationTolerance<Type>() )
         || !( maxNorm(CTL(VL) * CTB(T) * CTR(trans(VR)) - CTB(B)) < detail::schurFactorizationTolerance<Type>() ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Matrix generalized Schur factorization failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Left Schur vectors:\n" << VL << "\n"
             << "   Right Schur vectors:\n" << VR << "\n"
             << "   alpha:\n" << alpha << "\n"
             << "   beta:\n" << beta << "\n"
             << "   Residual A:\n" << CTL(VL) * CTA(S) * CTR(trans(VR)) - CTA(A) << "\n"
             << "   Residual B:\n" << CTL(VL) * CTB(T) * CTR(trans(VR)) - CTB(B) << "\n"
             << "   ";
         throw std::runtime_error( oss.str() );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of generalized Schur factorization functions for general matrices (gges) 
// with eigenvalue selection.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the eigenvalue functions for general matrices for various
// data types. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type, bool SOA, bool SOB, bool SOL, bool SOR >
void EigenvalueTest::testGgesSelect()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   test_ = "General matrix eigenvalue and Schur form computation (gges)";

   using CT = blaze::If_t< blaze::IsComplex_v<Type>, Type, blaze::complex<Type> >;
   using RT = typename CT::value_type;
   using detail::ConditionalTranspose;

   blaze::StaticMatrix<Type,3UL,3UL,SOA> A;
   blaze::StaticMatrix<Type,3UL,3UL,SOB> B;
   randomize( A );
   randomize( B );

   blaze::StaticMatrix<Type,3UL,3UL,SOL> VL;
   blaze::StaticMatrix<Type,3UL,3UL,SOR> VR;
   
   blaze::StaticMatrix<Type,3UL,3UL,SOA> S(A);
   blaze::StaticMatrix<Type,3UL,3UL,SOB> T(B);
   blaze::StaticVector<CT,3UL,blaze::rowVector> alpha;
   blaze::StaticVector<Type,3UL,blaze::rowVector> beta;

   auto selctg = [] (RT const * alphar, RT const * alphai, RT const * beta) -> int
   {
      return *alphar > 0.;
   };

   blaze::gges( selctg, S, T, alpha, beta, VL, VR );

   ConditionalTranspose<SOA == blaze::rowMajor> const CTA;
   ConditionalTranspose<SOB == blaze::rowMajor> const CTB;
   ConditionalTranspose<SOL == blaze::rowMajor> const CTL;
   ConditionalTranspose<SOR == blaze::rowMajor> const CTR;

   for (std::size_t i = 0; i < alpha.size(); ++i)
      checkEigenvalue(CTA(A), CTB(B), alpha[i], beta[i]);

   if( !( maxNorm(CTL(VL) * CTA(S) * CTR(trans(VR)) - CTA(A)) < detail::schurFactorizationTolerance<Type>() )
      || !( maxNorm(CTL(VL) * CTB(T) * CTR(trans(VR)) - CTB(B)) < detail::schurFactorizationTolerance<Type>() ) ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
            << " Error: Matrix generalized Schur factorization failed\n"
            << " Details:\n"
            << "   Random seed = " << blaze::getSeed() << "\n"
            << "   Element type:\n"
            << "     " << typeid( Type ).name() << "\n"
            << "   Left Schur vectors:\n" << VL << "\n"
            << "   Right Schur vectors:\n" << VR << "\n"
            << "   alpha:\n" << alpha << "\n"
            << "   beta:\n" << beta << "\n"
            << "   Residual A:\n" << CTL(VL) * CTA(S) * CTR(trans(VR)) - CTA(A) << "\n"
            << "   Residual B:\n" << CTL(VL) * CTB(T) * CTR(trans(VR)) - CTB(B) << "\n"
            << "   ";
      throw std::runtime_error( oss.str() );
   }

   // Check eigenvalue order: the values selected by selctg should go first.
   std::vector<bool> selected(size(alpha));
   for (std::size_t i = 0; i < size(alpha); ++i)
   {
      RT const alphar = real(alpha[i]);
      RT const alphai = imag(alpha[i]);
      selected[i] = selctg(&alphar, &alphai, &beta[i]);

      if (i > 0 && selected[i] && !selected[i - 1])
      {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
            << " Error: wrong eigenvalue order\n"
            << " Details:\n"
            << "   Random seed = " << blaze::getSeed() << "\n"
            << "   Element type:\n"
            << "     " << typeid( Type ).name() << "\n"
            << "   Left Schur vectors:\n" << VL << "\n"
            << "   Right Schur vectors:\n" << VR << "\n"
            << "   alpha:\n" << alpha << "\n"
            << "   beta:\n" << beta << "\n"
            << "   Residual A:\n" << CTL(VL) * CTA(S) * CTR(trans(VR)) - CTA(A) << "\n"
            << "   Residual B:\n" << CTL(VL) * CTB(T) * CTR(trans(VR)) - CTB(B) << "\n"
            << "   ";
         throw std::runtime_error( oss.str() );
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
/*!\brief Checking the given right eigenvector.
//
// \param v The right eigenvector to be checked.
// \param A The corresponding dense matrix.
// \param w The corresponding eigenvalue.
// \return void
// \exception std::runtime_error Invalid right eigenvector detected.
//
// This function checks the given right eigenvector \f$v[j]\f$ by testing if it satisfies

                          \f[ A * v[j] = lambda[j] * v[j], \f]

// where \f$lambda[j]\f$ is the corresponding eigenvalue.
*/
template< typename VT    // Type of the eigenvector v
        , typename MT    // Type of the matrix A
        , bool SO        // Storage order of the matrix A
        , typename ST >  // Type of the eigenvalue w
void EigenvalueTest::checkEigenvector( const blaze::DenseVector<VT,false>& v,
                                       const blaze::DenseMatrix<MT,SO>& A, ST w )
{
   if( (~A) * (~v) != w * (~v) ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Invalid right eigenvector detected\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   System matrix:\n" << (~A) << "\n"
          << "   Eigenvalue = " << w << "\n"
          << "   Right eigenvector:\n" << (~v) << "\n"
          << "   A * v =\n" << ( (~A) * (~v) ) << "\n"
          << "   A * w =\n" << ( w * (~v) ) << "\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checking the given left eigenvector.
//
// \param u The left eigenvector to be checked.
// \param A The corresponding dense matrix.
// \param w The corresponding eigenvalue.
// \return void
// \exception std::runtime_error Invalid left eigenvector detected.
//
// This function checks the given left eigenvector \f$u[j]\f$ by testing if it satisfies

                       \f[ u[j]^{H} * A = lambda[j] * u[j]^{H}, \f]

// where \f$lambda[j]\f$ is the corresponding eigenvalue.
*/
template< typename VT    // Type of the eigenvector u
        , typename MT    // Type of the matrix A
        , bool SO        // Storage order of the matrix A
        , typename ST >  // Type of the eigenvalue w
void EigenvalueTest::checkEigenvector( const blaze::DenseVector<VT,true>& u,
                                       const blaze::DenseMatrix<MT,SO>& A, ST w )
{
   if( (~u) * (~A) != (~u) * w ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Invalid left eigenvector detected\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   System matrix:\n" << (~A) << "\n"
          << "   Eigenvalue = " << w << "\n"
          << "   Left eigenvector:\n" << (~u) << "\n"
          << "   v * A =\n" << ( (~u) * (~A) ) << "\n"
          << "   v * w =\n" << ( (~u) * w ) << "\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checking the given generalized eigenvalue.
//
// \param A The first corresponding dense matrix.
// \param B The second corresponding dense matrix.
// \param alpha The numerator of the generalized eigenvalue.
// \param beta The denominator of the generalized eigenvalue.
// \return void
// \exception std::runtime_error Invalid generalized eigenvalue detected.
//
// This function checks the given generalized eigenvalue \f$\lambda=\frac{\alpha}{\beta}\f$ by testing if it satisfies

                       \f[\det(\beta A - \alpha B) = 0 \f]

// where \f$\lambda=\frac{\alpha}{\beta}\f$ is the corresponding eigenvalue.
*/
template< typename MT1, bool SO1, typename MT2, bool SO2, typename ST1, typename ST2 >
void EigenvalueTest::checkEigenvalue( const blaze::DenseMatrix<MT1, SO1>& A,
                        const blaze::DenseMatrix<MT2, SO2>& B, 
                        ST1 alpha, ST2 beta )
{
   BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( ST1 );

   // M must be singular
   auto M = evaluate(beta * A - alpha * B);

   // Determine the type of M an the type of M elements
   using MType = decltype(M);
   using MElementType = blaze::ElementType_t<MType>;
   
   // M must be complex, because alpha is complex.
   BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( MElementType );
   using MElementScalarType = typename MElementType::value_type;

   // Evaluate eigenvalues of M
   // The initial content of M is destroyed by geev().
   blaze::DynamicVector<MElementType, blaze::rowVector> w;
   geev(M, w);

   // Check that at least one of the eigenvalues of M is close to 0 with specified tolerance.
   if (std::find_if(begin(w), end(w), 
      [] (MElementType w_i) { return abs(w_i) < detail::singularEigenvalueThreshold<MElementScalarType>(); }) == end(w))
   {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Invalid generalized eigenvalue detected (test matrix M=beta*A-alpha*B is not singular)\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Element type A:\n"
          << "     " << typeid( typename MT1::ElementType ).name() << "\n"
          << "   Element type B:\n"
          << "     " << typeid( typename MT2::ElementType ).name() << "\n"
          << "   Type alpha:\n"
          << "     " << typeid( ST1 ).name() << "\n"
          << "   Type beta:\n"
          << "     " << typeid( ST2 ).name() << "\n"
          << "   System matrix A:\n" << (~A) << "\n"
          << "   System matrix B:\n" << (~B) << "\n"
          << "   Eigenvalue numerator = " << alpha << "\n"
          << "   Eigenvalue denominator = " << beta << "\n"
          << "   Test matrix M = " << beta * A - alpha * B << "\n"
          << "   Eigenvalues of M = " << w << "\n";
      throw std::runtime_error( oss.str() );
   }
}



//=================================================================================================
//
//  GLOBAL TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Testing the LAPACK functionality.
//
// \return void
*/
void runTest()
{
   EigenvalueTest();
}
//*************************************************************************************************




//=================================================================================================
//
//  MACRO DEFINITIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Macro for the execution of the LAPACK eigenvalue test.
*/
#define RUN_LAPACK_EIGENVALUE_TEST \
   blazetest::mathtest::lapack::runTest()
/*! \endcond */
//*************************************************************************************************

} // namespace lapack

} // namespace mathtest

} // namespace blazetest

#endif
