//=================================================================================================
/*!
//  \file blaze/math/adaptors/SymmetricMatrix.h
//  \brief Header file for the implementation of a symmetric matrix adaptor
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

#ifndef _BLAZE_MATH_ADAPTORS_SYMMETRICMATRIX_H_
#define _BLAZE_MATH_ADAPTORS_SYMMETRICMATRIX_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/adaptors/symmetricmatrix/BaseTemplate.h>
#include <blaze/math/adaptors/symmetricmatrix/DenseNonNumeric.h>
#include <blaze/math/adaptors/symmetricmatrix/DenseNumeric.h>
#include <blaze/math/adaptors/symmetricmatrix/SparseNonNumeric.h>
#include <blaze/math/adaptors/symmetricmatrix/SparseNumeric.h>
#include <blaze/math/Forward.h>
#include <blaze/math/traits/AddTrait.h>
#include <blaze/math/traits/ColumnTrait.h>
#include <blaze/math/traits/DivTrait.h>
#include <blaze/math/traits/MathTrait.h>
#include <blaze/math/traits/MultTrait.h>
#include <blaze/math/traits/RowTrait.h>
#include <blaze/math/traits/SubmatrixTrait.h>
#include <blaze/math/traits/SubTrait.h>
#include <blaze/math/typetraits/HasConstDataAccess.h>
#include <blaze/math/typetraits/Columns.h>
#include <blaze/math/typetraits/IsAdaptor.h>
#include <blaze/math/typetraits/IsResizable.h>
#include <blaze/math/typetraits/IsSquare.h>
#include <blaze/math/typetraits/IsSymmetric.h>
#include <blaze/math/typetraits/RemoveAdaptor.h>
#include <blaze/math/typetraits/Rows.h>
#include <blaze/util/constraints/Numeric.h>


namespace blaze {

//=================================================================================================
//
//  SYMMETRICMATRIX OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name SymmetricMatrix operators */
//@{
template< typename MT, bool SO, bool DF, bool NF >
inline void reset( SymmetricMatrix<MT,SO,DF,NF>& m );

template< typename MT, bool SO, bool DF, bool NF >
inline void reset( SymmetricMatrix<MT,SO,DF,NF>& m, size_t i );

template< typename MT, bool SO, bool DF, bool NF >
inline void clear( SymmetricMatrix<MT,SO,DF,NF>& m );

template< typename MT, bool SO, bool DF, bool NF >
inline bool isDefault( const SymmetricMatrix<MT,SO,DF,NF>& m );

template< typename MT, bool SO, bool DF, bool NF >
inline void swap( SymmetricMatrix<MT,SO,DF,NF>& a, SymmetricMatrix<MT,SO,DF,NF>& b ) /* throw() */;
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Resetting the given symmetric matrix.
// \ingroup symmetric_matrix
//
// \param m The symmetric matrix to be resetted.
// \return void
*/
template< typename MT  // Type of the adapted matrix
        , bool SO      // Storage order of the adapted matrix
        , bool DF      // Density flag
        , bool NF >    // Numeric flag
inline void reset( SymmetricMatrix<MT,SO,DF,NF>& m )
{
   m.reset();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Resetting the specified row/column of the given symmetric matrix.
// \ingroup symmetric_matrix
//
// \param m The symmetric matrix to be resetted.
// \param i The index of the row/column to be resetted.
// \return void
//
// This function resets the values in the specified row/column of the given symmetric matrix to
// their default value. In case the given matrix is a \a rowMajor matrix the function resets the
// values in row \a i, if it is a \a columnMajor matrix the function resets the values in column
// \a i. Note that the capacity of the row/column remains unchanged.
*/
template< typename MT  // Type of the adapted matrix
        , bool SO      // Storage order of the adapted matrix
        , bool DF      // Density flag
        , bool NF >    // Numeric flag
inline void reset( SymmetricMatrix<MT,SO,DF,NF>& m, size_t i )
{
   m.reset( i );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Clearing the given symmetric matrix.
// \ingroup symmetric_matrix
//
// \param m The symmetric matrix to be cleared.
// \return void
*/
template< typename MT  // Type of the adapted matrix
        , bool SO      // Storage order of the adapted matrix
        , bool DF      // Density flag
        , bool NF >    // Numeric flag
inline void clear( SymmetricMatrix<MT,SO,DF,NF>& m )
{
   m.clear();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the given symmetric matrix is in default state.
// \ingroup symmetric_matrix
//
// \param m The symmetric matrix to be tested for its default state.
// \return \a true in case the given matrix is component-wise zero, \a false otherwise.
//
// This function checks whether the matrix is in default state. For instance, in case the
// matrix is instantiated for a built-in integral or floating point data type, the function
// returns \a true in case all matrix elements are 0 and \a false in case any matrix element
// is not 0. The following example demonstrates the use of the \a isDefault function:

   \code
   blaze::SymmetricMatrix<int> A;
   // ... Resizing and initialization
   if( isDefault( A ) ) { ... }
   \endcode
*/
template< typename MT  // Type of the adapted matrix
        , bool SO      // Storage order of the adapted matrix
        , bool DF      // Density flag
        , bool NF >    // Numeric flag
inline bool isDefault( const SymmetricMatrix<MT,SO,DF,NF>& m )
{
   return isDefault( m.matrix_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Swapping the contents of two matrices.
// \ingroup symmetric_matrix
//
// \param a The first matrix to be swapped.
// \param b The second matrix to be swapped.
// \return void
// \exception no-throw guarantee.
*/
template< typename MT  // Type of the adapted matrix
        , bool SO      // Storage order of the adapted matrix
        , bool DF      // Density flag
        , bool NF >    // Numeric flag
inline void swap( SymmetricMatrix<MT,SO,DF,NF>& a, SymmetricMatrix<MT,SO,DF,NF>& b ) /* throw() */
{
   a.swap( b );
}
//*************************************************************************************************




//=================================================================================================
//
//  ROWS SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF, bool NF >
struct Rows< SymmetricMatrix<MT,SO,DF,NF> > : public Rows<MT>
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  COLUMNS SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF, bool NF >
struct Columns< SymmetricMatrix<MT,SO,DF,NF> > : public Columns<MT>
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISSQUARE SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF, bool NF >
struct IsSquare< SymmetricMatrix<MT,SO,DF,NF> > : public TrueType
{
   enum { value = 1 };
   typedef TrueType  Type;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISSYMMETRIC SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF, bool NF >
struct IsSymmetric< SymmetricMatrix<MT,SO,DF,NF> > : public TrueType
{
   enum { value = 1 };
   typedef TrueType  Type;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  HASCONSTDATAACCESS SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool NF >
struct HasConstDataAccess< SymmetricMatrix<MT,SO,true,NF> > : public TrueType
{
   enum { value = 1 };
   typedef TrueType  Type;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISRESIZABLE SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF, bool NF >
struct IsResizable< SymmetricMatrix<MT,SO,DF,NF> > : public IsResizable<MT>::Type
{
   enum { value = IsResizable<MT>::value };
   typedef typename IsResizable<MT>::Type  Type;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISADAPTOR SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF, bool NF >
struct IsAdaptor< SymmetricMatrix<MT,SO,DF,NF> > : public TrueType
{
   enum { value = 1 };
   typedef TrueType  Type;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  REMOVEADAPTOR SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF, bool NF >
struct RemoveAdaptor< SymmetricMatrix<MT,SO,DF,NF> >
{
   typedef MT  Type;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ADDTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO1, bool DF, bool NF, typename T, size_t M, size_t N, bool SO2 >
struct AddTrait< SymmetricMatrix<MT,SO1,DF,NF>, StaticMatrix<T,M,N,SO2> >
{
   typedef typename AddTrait< MT, StaticMatrix<T,M,N,SO2> >::Type  Type;
};

template< typename T, size_t M, size_t N, bool SO1, typename MT, bool SO2, bool DF, bool NF >
struct AddTrait< StaticMatrix<T,M,N,SO1>, SymmetricMatrix<MT,SO2,DF,NF> >
{
   typedef typename AddTrait< StaticMatrix<T,M,N,SO1>, MT >::Type  Type;
};

template< typename MT, bool SO1, bool DF, bool NF, typename T, size_t M, size_t N, bool SO2 >
struct AddTrait< SymmetricMatrix<MT,SO1,DF,NF>, HybridMatrix<T,M,N,SO2> >
{
   typedef typename AddTrait< MT, HybridMatrix<T,M,N,SO2> >::Type  Type;
};

template< typename T, size_t M, size_t N, bool SO1, typename MT, bool SO2, bool DF, bool NF >
struct AddTrait< HybridMatrix<T,M,N,SO1>, SymmetricMatrix<MT,SO2,DF,NF> >
{
   typedef typename AddTrait< HybridMatrix<T,M,N,SO1>, MT >::Type  Type;
};

template< typename MT, bool SO1, bool DF, bool NF, typename T, bool SO2 >
struct AddTrait< SymmetricMatrix<MT,SO1,DF,NF>, DynamicMatrix<T,SO2> >
{
   typedef typename AddTrait< MT, DynamicMatrix<T,SO2> >::Type  Type;
};

template< typename T, bool SO1, typename MT, bool SO2, bool DF, bool NF >
struct AddTrait< DynamicMatrix<T,SO1>, SymmetricMatrix<MT,SO2,DF,NF> >
{
   typedef typename AddTrait< DynamicMatrix<T,SO1>, MT >::Type  Type;
};

template< typename MT, bool SO1, bool DF, bool NF, typename T, bool SO2 >
struct AddTrait< SymmetricMatrix<MT,SO1,DF,NF>, CompressedMatrix<T,SO2> >
{
   typedef typename AddTrait< MT, CompressedMatrix<T,SO2> >::Type  Type;
};

template< typename T, bool SO1, typename MT, bool SO2, bool DF, bool NF >
struct AddTrait< CompressedMatrix<T,SO1>, SymmetricMatrix<MT,SO2,DF,NF> >
{
   typedef typename AddTrait< CompressedMatrix<T,SO1>, MT >::Type  Type;
};

template< typename MT1, bool SO1, bool DF1, bool NF1, typename MT2, bool SO2, bool DF2, bool NF2 >
struct AddTrait< SymmetricMatrix<MT1,SO1,DF1,NF1>, SymmetricMatrix<MT2,SO2,DF2,NF2> >
{
   typedef SymmetricMatrix< typename AddTrait<MT1,MT2>::Type >  Type;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  SUBTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO1, bool DF, bool NF, typename T, size_t M, size_t N, bool SO2 >
struct SubTrait< SymmetricMatrix<MT,SO1,DF,NF>, StaticMatrix<T,M,N,SO2> >
{
   typedef typename SubTrait< MT, StaticMatrix<T,M,N,SO2> >::Type  Type;
};

template< typename T, size_t M, size_t N, bool SO1, typename MT, bool SO2, bool DF, bool NF >
struct SubTrait< StaticMatrix<T,M,N,SO1>, SymmetricMatrix<MT,SO2,DF,NF> >
{
   typedef typename SubTrait< StaticMatrix<T,M,N,SO1>, MT >::Type  Type;
};

template< typename MT, bool SO1, bool DF, bool NF, typename T, size_t M, size_t N, bool SO2 >
struct SubTrait< SymmetricMatrix<MT,SO1,DF,NF>, HybridMatrix<T,M,N,SO2> >
{
   typedef typename SubTrait< MT, HybridMatrix<T,M,N,SO2> >::Type  Type;
};

template< typename T, size_t M, size_t N, bool SO1, typename MT, bool SO2, bool DF, bool NF >
struct SubTrait< HybridMatrix<T,M,N,SO1>, SymmetricMatrix<MT,SO2,DF,NF> >
{
   typedef typename SubTrait< HybridMatrix<T,M,N,SO1>, MT >::Type  Type;
};

template< typename MT, bool SO1, bool DF, bool NF, typename T, bool SO2 >
struct SubTrait< SymmetricMatrix<MT,SO1,DF,NF>, DynamicMatrix<T,SO2> >
{
   typedef typename SubTrait< MT, DynamicMatrix<T,SO2> >::Type  Type;
};

template< typename T, bool SO1, typename MT, bool SO2, bool DF, bool NF >
struct SubTrait< DynamicMatrix<T,SO1>, SymmetricMatrix<MT,SO2,DF,NF> >
{
   typedef typename SubTrait< DynamicMatrix<T,SO1>, MT >::Type  Type;
};

template< typename MT, bool SO1, bool DF, bool NF, typename T, bool SO2 >
struct SubTrait< SymmetricMatrix<MT,SO1,DF,NF>, CompressedMatrix<T,SO2> >
{
   typedef typename SubTrait< MT, CompressedMatrix<T,SO2> >::Type  Type;
};

template< typename T, bool SO1, typename MT, bool SO2, bool DF, bool NF >
struct SubTrait< CompressedMatrix<T,SO1>, SymmetricMatrix<MT,SO2,DF,NF> >
{
   typedef typename SubTrait< CompressedMatrix<T,SO1>, MT >::Type  Type;
};

template< typename MT1, bool SO1, bool DF1, bool NF1, typename MT2, bool SO2, bool DF2, bool NF2 >
struct SubTrait< SymmetricMatrix<MT1,SO1,DF1,NF1>, SymmetricMatrix<MT2,SO2,DF2,NF2> >
{
   typedef SymmetricMatrix< typename SubTrait<MT1,MT2>::Type >  Type;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  MULTTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF, bool NF, typename T >
struct MultTrait< SymmetricMatrix<MT,SO,DF,NF>, T >
{
   typedef SymmetricMatrix< typename MultTrait<MT,T>::Type >  Type;
   BLAZE_CONSTRAINT_MUST_BE_NUMERIC_TYPE( T );
};

template< typename T, typename MT, bool SO, bool DF, bool NF >
struct MultTrait< T, SymmetricMatrix<MT,SO,DF,NF> >
{
   typedef SymmetricMatrix< typename MultTrait<T,MT>::Type >  Type;
   BLAZE_CONSTRAINT_MUST_BE_NUMERIC_TYPE( T );
};

template< typename MT, bool SO, bool DF, bool NF, typename T, size_t N >
struct MultTrait< SymmetricMatrix<MT,SO,DF,NF>, StaticVector<T,N,false> >
{
   typedef typename MultTrait< MT, StaticVector<T,N,false> >::Type  Type;
};

template< typename T, size_t N, typename MT, bool SO, bool DF, bool NF >
struct MultTrait< StaticVector<T,N,true>, SymmetricMatrix<MT,SO,DF,NF> >
{
   typedef typename MultTrait< StaticVector<T,N,true>, MT >::Type  Type;
};

template< typename MT, bool SO, bool DF, bool NF, typename T, size_t N >
struct MultTrait< SymmetricMatrix<MT,SO,DF,NF>, HybridVector<T,N,false> >
{
   typedef typename MultTrait< MT, HybridVector<T,N,false> >::Type  Type;
};

template< typename T, size_t N, typename MT, bool SO, bool DF, bool NF >
struct MultTrait< HybridVector<T,N,true>, SymmetricMatrix<MT,SO,DF,NF> >
{
   typedef typename MultTrait< HybridVector<T,N,true>, MT >::Type  Type;
};

template< typename MT, bool SO, bool DF, bool NF, typename T >
struct MultTrait< SymmetricMatrix<MT,SO,DF,NF>, DynamicVector<T,false> >
{
   typedef typename MultTrait< MT, DynamicVector<T,false> >::Type  Type;
};

template< typename T, typename MT, bool SO, bool DF, bool NF >
struct MultTrait< DynamicVector<T,true>, SymmetricMatrix<MT,SO,DF,NF> >
{
   typedef typename MultTrait< DynamicVector<T,true>, MT >::Type  Type;
};

template< typename MT, bool SO, bool DF, bool NF, typename T >
struct MultTrait< SymmetricMatrix<MT,SO,DF,NF>, CompressedVector<T,false> >
{
   typedef typename MultTrait< MT, CompressedVector<T,false> >::Type  Type;
};

template< typename T, typename MT, bool SO, bool DF, bool NF >
struct MultTrait< CompressedVector<T,true>, SymmetricMatrix<MT,SO,DF,NF> >
{
   typedef typename MultTrait< CompressedVector<T,true>, MT >::Type  Type;
};

template< typename MT, bool SO1, bool DF, bool NF, typename T, size_t M, size_t N, bool SO2 >
struct MultTrait< SymmetricMatrix<MT,SO1,DF,NF>, StaticMatrix<T,M,N,SO2> >
{
   typedef typename MultTrait< MT, StaticMatrix<T,M,N,SO2> >::Type  Type;
};

template< typename T, size_t M, size_t N, bool SO1, typename MT, bool SO2, bool DF, bool NF >
struct MultTrait< StaticMatrix<T,M,N,SO1>, SymmetricMatrix<MT,SO2,DF,NF> >
{
   typedef typename MultTrait< StaticMatrix<T,M,N,SO1>, MT >::Type  Type;
};

template< typename MT, bool SO1, bool DF, bool NF, typename T, size_t M, size_t N, bool SO2 >
struct MultTrait< SymmetricMatrix<MT,SO1,DF,NF>, HybridMatrix<T,M,N,SO2> >
{
   typedef typename MultTrait< MT, HybridMatrix<T,M,N,SO2> >::Type  Type;
};

template< typename T, size_t M, size_t N, bool SO1, typename MT, bool SO2, bool DF, bool NF >
struct MultTrait< HybridMatrix<T,M,N,SO1>, SymmetricMatrix<MT,SO2,DF,NF> >
{
   typedef typename MultTrait< HybridMatrix<T,M,N,SO1>, MT >::Type  Type;
};

template< typename MT, bool SO1, bool DF, bool NF, typename T, bool SO2 >
struct MultTrait< SymmetricMatrix<MT,SO1,DF,NF>, DynamicMatrix<T,SO2> >
{
   typedef typename MultTrait< MT, DynamicMatrix<T,SO2> >::Type  Type;
};

template< typename T, bool SO1, typename MT, bool SO2, bool DF, bool NF >
struct MultTrait< DynamicMatrix<T,SO1>, SymmetricMatrix<MT,SO2,DF,NF> >
{
   typedef typename MultTrait< DynamicMatrix<T,SO1>, MT >::Type  Type;
};

template< typename MT, bool SO1, bool DF, bool NF, typename T, bool SO2 >
struct MultTrait< SymmetricMatrix<MT,SO1,DF,NF>, CompressedMatrix<T,SO2> >
{
   typedef typename MultTrait< MT, CompressedMatrix<T,SO2> >::Type  Type;
};

template< typename T, bool SO1, typename MT, bool SO2, bool DF, bool NF >
struct MultTrait< CompressedMatrix<T,SO1>, SymmetricMatrix<MT,SO2,DF,NF> >
{
   typedef typename MultTrait< CompressedMatrix<T,SO1>, MT >::Type  Type;
};

template< typename MT1, bool SO1, bool DF1, bool NF1, typename MT2, bool SO2, bool DF2, bool NF2 >
struct MultTrait< SymmetricMatrix<MT1,SO1,DF1,NF1>, SymmetricMatrix<MT2,SO2,DF2,NF2> >
{
   typedef typename MultTrait<MT1,MT2>::Type  Type;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DIVTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF, bool NF, typename T >
struct DivTrait< SymmetricMatrix<MT,SO,DF,NF>, T >
{
   typedef SymmetricMatrix< typename DivTrait<MT,T>::Type >  Type;
   BLAZE_CONSTRAINT_MUST_BE_NUMERIC_TYPE( T );
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  MATHTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, bool SO1, bool DF1, bool NF1, typename MT2, bool SO2, bool DF2, bool NF2 >
struct MathTrait< SymmetricMatrix<MT1,SO1,DF1,NF1>, SymmetricMatrix<MT2,SO2,DF2,NF2> >
{
   typedef SymmetricMatrix< typename MathTrait<MT1,MT2>::HighType >  HighType;
   typedef SymmetricMatrix< typename MathTrait<MT1,MT2>::LowType  >  LowType;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  SUBMATRIXTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF, bool NF >
struct SubmatrixTrait< SymmetricMatrix<MT,SO,DF,NF> >
{
   typedef typename SubmatrixTrait<MT>::Type  Type;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ROWTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF, bool NF >
struct RowTrait< SymmetricMatrix<MT,SO,DF,NF> >
{
   typedef typename RowTrait<MT>::Type  Type;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  COLUMNTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF, bool NF >
struct ColumnTrait< SymmetricMatrix<MT,SO,DF,NF> >
{
   typedef typename ColumnTrait<MT>::Type  Type;
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
