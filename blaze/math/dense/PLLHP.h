#ifndef _BLAZE_MATH_DENSE_PLLHP_H_
#define _BLAZE_MATH_DENSE_PLLHP_H_

//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/Aliases.h>
#include <blaze/math/constraints/BLASCompatible.h>
#include <blaze/math/constraints/Hermitian.h>
#include <blaze/math/constraints/StrictlyTriangular.h>
#include <blaze/math/constraints/Symmetric.h>
#include <blaze/math/constraints/UniTriangular.h>
#include <blaze/math/constraints/Upper.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/DenseMatrix.h>
#include <blaze/math/lapack/pstrf.h>
#include <blaze/math/typetraits/IsResizable.h>
#include <blaze/math/typetraits/IsRowMajorMatrix.h>


namespace blaze {

template< typename MT1, bool SO1, typename MT2, bool SO2 >
int pllhp(const DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& L, blas_int_t* P,
		  typename RemoveComplex<ElementType_t<MT2>>::type tol );


template< typename MT1  // Type of matrix A
        , bool SO1      // Storage order of matrix A
        , typename MT2  // Type of matrix L
        , bool SO2      // Storage order of matrix L
        >
int pllhp( const DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& L, blas_int_t* P,
		   typename RemoveComplex<ElementType_t<MT2>>::type tol )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT1> );

   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT2 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT2 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UPPER_MATRIX_TYPE( MT2 );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT2> );

   if( !isSquare( *A ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid non-square matrix provided" );
   }

   const size_t n( (*A).rows() );

   if( ( !IsResizable_v<MT2> && ( (*L).rows() != n || (*L).columns() != n ) ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Dimensions of fixed size matrix do not match" );
   }

   decltype(auto) l( derestrict( *L ) );

   resize( *L, n, n, false );
   reset( l );

   if( IsRowMajorMatrix_v<MT2> ) {
      for( size_t i=0UL; i<n; ++i ) {
         for( size_t j=0UL; j<=i; ++j ) {
            l(i,j) = (*A)(i,j);
         }
      }
   }
   else {
      for( size_t j=0UL; j<n; ++j ) {
         for( size_t i=j; i<n; ++i ) {
            l(i,j) = (*A)(i,j);
         }
      }
   }

   return pstrf( l, 'L', P, tol );
}
//*************************************************************************************************

} // namespace blaze

#endif
