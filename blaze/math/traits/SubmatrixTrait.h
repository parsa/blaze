//=================================================================================================
/*!
//  \file blaze/math/traits/SubmatrixTrait.h
//  \brief Header file for the submatrix trait
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

#ifndef _BLAZE_MATH_TRAITS_SUBMATRIXTRAIT_H_
#define _BLAZE_MATH_TRAITS_SUBMATRIXTRAIT_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/InvalidType.h>
#include <blaze/util/SelectType.h>
#include <blaze/util/typetraits/IsConst.h>
#include <blaze/util/typetraits/IsReference.h>
#include <blaze/util/typetraits/IsVolatile.h>
#include <blaze/util/typetraits/RemoveCV.h>
#include <blaze/util/typetraits/RemoveReference.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Base template for the SubmatrixTrait class.
// \ingroup math
//
// \section submatrixtrait_general General
//
// The SubmatrixTrait class template offers the possibility to select the resulting data type
// when creating a submatrix of a dense or sparse matrix. SubmatrixTrait defines the nested
// type \a Type, which represents the resulting data type of the submatrix operation. In case
// the given data type is not a dense or sparse matrix type, the resulting data type \a Type
// is set to \a INVALID_TYPE. Note that \a const and \a volatile qualifiers and reference
// modifiers are generally ignored.
//
// Per default, the SubmatrixTrait template only supports the following matrix types:
//
// <ul>
//    <li>blaze::StaticMatrix</li>
//    <li>blaze::DynamicMatrix</li>
//    <li>blaze::CompressedMatrix</li>
// </ul>
//
//
// \section submatrixtrait_specializations Creating custom specializations
//
// It is possible to specialize the SubmatrixTrait template for additional user-defined matrix
// types. The following example shows the according specialization for the DynamicMatrix class
// template:

   \code
   template< typename T1, bool SO >
   struct SubmatrixTrait< DynamicMatrix<T1,SO> >
   {
      typedef DynamicMatrix<T1,SO>  Type;
   };
   \endcode

// \n \section submatrixtrait_examples Examples
//
// The following example demonstrates the use of the SubmatrixTrait template, where depending
// on the given matrix type the according result type is selected:

   \code
   using blaze::rowMajor;
   using blaze::columnMajor;

   // Definition of the result type of a row-major dynamic matrix
   typedef blaze::DynamicMatrix<int,rowMajor>          MatrixType1;
   typedef typename SubmatrixTrait<MatrixType1>::Type  ResultType1;

   // Definition of the result type of a column-major static matrix
   typedef blaze::StaticMatrix<int,3UL,3UL,columnMajor>  MatrixType2;
   typedef typename SubmatrixTrait<MatrixType2>::Type    ResultType2;
   \endcode
*/
template< typename MT >  // Type of the matrix
struct SubmatrixTrait
{
 private:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   struct Failure { typedef INVALID_TYPE  Type; };
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   enum { qualified = IsConst<MT>::value || IsVolatile<MT>::value || IsReference<MT>::value };
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   typedef typename RemoveReference< typename RemoveCV<MT>::Type >::Type  Tmp;
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   typedef typename SelectType< qualified, SubmatrixTrait<Tmp>, Failure >::Type::Type  Type;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************

} // namespace blaze

#endif
