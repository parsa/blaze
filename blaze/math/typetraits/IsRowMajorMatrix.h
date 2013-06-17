//=================================================================================================
/*!
//  \file blaze/math/typetraits/IsRowMajorMatrix.h
//  \brief Header file for the IsRowMajorMatrix type trait
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

#ifndef _BLAZE_MATH_TYPETRAITS_ISROWMAJORMATRIX_H_
#define _BLAZE_MATH_TYPETRAITS_ISROWMAJORMATRIX_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <boost/type_traits/is_base_of.hpp>
#include <boost/type_traits/remove_cv.hpp>
#include <blaze/math/expressions/DenseMatrix.h>
#include <blaze/math/expressions/SparseMatrix.h>
#include <blaze/util/FalseType.h>
#include <blaze/util/SelectType.h>
#include <blaze/util/TrueType.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Auxiliary helper struct for the IsRowMajorMatrix type trait.
// \ingroup math_type_traits
*/
template< typename T >
struct IsRowMajorMatrixHelper
{
 private:
   //**********************************************************************************************
   typedef typename boost::remove_cv<T>::type  T2;
   //**********************************************************************************************

 public:
   //**********************************************************************************************
   enum { value = boost::is_base_of< DenseMatrix <T2,false>, T2 >::value ||
                  boost::is_base_of< SparseMatrix<T2,false>, T2 >::value };
   typedef typename SelectType<value,TrueType,FalseType>::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compile time check for row-major matrix types.
// \ingroup math_type_traits
//
// This type trait tests whether or not the given template argument is a row-major dense or
// sparse matrix type (i.e., a matrix whose storage order is set to \a true). In case the type
// is a row-major matrix type, the \a value member enumeration is set to 1, the nested type
// definition \a Type is \a TrueType, and the class derives from \a TrueType. Otherwise
// \a value is set to 0, \a Type is \a FalseType, and the class derives from \a FalseType.

   \code
   blaze::IsRowMajorMatrix< StaticMatrix<float,3U,3U,false> >::value   // Evaluates to 1
   blaze::IsRowMajorMatrix< const DynamicMatrix<double,false> >::Type  // Results in TrueType
   blaze::IsRowMajorMatrix< volatile CompressedMatrix<int,false> >     // Is derived from TrueType
   blaze::IsRowMajorMatrix< StaticMatrix<float,3U,3U,true> >::value    // Evaluates to 0
   blaze::IsRowMajorMatrix< const DynamicMatrix<double,true> >::Type   // Results in FalseType
   blaze::IsRowMajorMatrix< volatile CompressedMatrix<int,true> >      // Is derived from FalseType
   \endcode
*/
template< typename T >
struct IsRowMajorMatrix : public IsRowMajorMatrixHelper<T>::Type
{
 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   enum { value = IsRowMajorMatrixHelper<T>::value };
   typedef typename IsRowMajorMatrixHelper<T>::Type  Type;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************

} // namespace blaze

#endif
