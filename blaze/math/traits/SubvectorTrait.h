//=================================================================================================
/*!
//  \file blaze/math/traits/SubvectorTrait.h
//  \brief Header file for the subvector trait
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

#ifndef _BLAZE_MATH_TRAITS_SUBVECTORTRAIT_H_
#define _BLAZE_MATH_TRAITS_SUBVECTORTRAIT_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/util/InvalidType.h>
#include <blaze/util/SelectType.h>
#include <blaze/util/typetraits/RemoveCV.h>
#include <blaze/util/typetraits/RemoveReference.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Base template for the SubvectorTrait class.
// \ingroup math
//
// \section subvectortrait_general General
//
// The SubvectorTrait class template offers the possibility to select the resulting data type
// when creating a subvector of a dense or sparse vector. SubvectorTrait defines the nested
// type \a Type, which represents the resulting data type of the subvector operation. In case
// the given data type is not a dense or sparse vector type, the resulting data type \a Type
// is set to \a INVALID_TYPE. Note that \a const and \a volatile qualifiers and reference
// modifiers are generally ignored.
//
// Per default, the SubvectorTrait template only supports the following vector types:
//
// <ul>
//    <li>blaze::StaticVector</li>
//    <li>blaze::DynamicVector</li>
//    <li>blaze::CompressedVector</li>
// </ul>
//
//
// \section subvectortrait_specializations Creating custom specializations
//
// It is possible to specialize the SubvectorTrait template for additional user-defined vector
// types. The following example shows the according specialization for the DynamicVector class
// template:

   \code
   template< typename T1, bool TF >
   struct SubvectorTrait< DynamicVector<T1,TF> >
   {
      typedef DynamicVector<T1,TF>  Type;
   };
   \endcode

// \n \section subvectortrait_examples Examples
//
// The following example demonstrates the use of the SubvectorTrait template, where depending
// on the given vector type the according result type is selected:

   \code
   using blaze::columnVector;
   using blaze::rowVector;

   // Definition of the result type of a dynamic column vector
   typedef blaze::DynamicVector<int,columnVector>      VectorType1;
   typedef typename SubvectorTrait<VectorType1>::Type  ResultType1;

   // Definition of the result type of the static row vector
   typedef blaze::StaticVector<int,3UL,rowVector>      VectorType2;
   typedef typename SubvectorTrait<VectorType2>::Type  ResultType2;
   \endcode
*/
template< typename VT >  // Type of the vector
struct SubvectorTrait
{
 private:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   struct Failure { typedef INVALID_TYPE  Type; };
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   enum { qualified = IsConst<VT>::value || IsVolatile<VT>::value || IsReference<VT>::value };
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   typedef typename RemoveReference< typename RemoveCV<VT>::Type >::Type  Tmp;
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   typedef typename SelectType< qualified, SubvectorTrait<Tmp>, Failure >::Type::Type  Type;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************

} // namespace blaze

#endif
