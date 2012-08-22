//=================================================================================================
/*!
//  \file blaze/math/typetraits/IsVecVecMultExpr.h
//  \brief Header file for the IsVecVecMultExpr type trait
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

#ifndef _BLAZE_MATH_TYPETRAITS_ISVECVECMULTEXPR_H_
#define _BLAZE_MATH_TYPETRAITS_ISVECVECMULTEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <boost/type_traits/remove_cv.hpp>
#include <blaze/math/expressions/Forward.h>
#include <blaze/util/FalseType.h>
#include <blaze/util/TrueType.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Auxiliary helper struct for the IsVecVecMultExpr type trait.
// \ingroup math_type_traits
*/
template< typename T >
struct IsVecVecMultExprHelper
{
   //**********************************************************************************************
   enum { value = 0 };
   typedef FalseType  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
//! Specialization of the IsVecVecMultExprHelper type trait for DVecDVecMultExpr.
template< typename VT1, typename VT2, bool TF >
struct IsVecVecMultExprHelper< DVecDVecMultExpr<VT1,VT2,TF> > : public TrueType
{
 public:
   //**********************************************************************************************
   enum { value = 1 };
   typedef TrueType  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
//! Specialization of the IsVecVecMultExprHelper type trait for DVecSVecMultExpr.
template< typename VT1, typename VT2, bool TF >
struct IsVecVecMultExprHelper< DVecSVecMultExpr<VT1,VT2,TF> > : public TrueType
{
 public:
   //**********************************************************************************************
   enum { value = 1 };
   typedef TrueType  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
//! Specialization of the IsVecVecMultExprHelper type trait for DVecTDVecMultExpr.
template< typename VT1, typename VT2 >
struct IsVecVecMultExprHelper< DVecTDVecMultExpr<VT1,VT2> > : public TrueType
{
 public:
   //**********************************************************************************************
   enum { value = 1 };
   typedef TrueType  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
//! Specialization of the IsVecVecMultExprHelper type trait for DVecTSVecMultExpr.
template< typename VT1, typename VT2 >
struct IsVecVecMultExprHelper< DVecTSVecMultExpr<VT1,VT2> > : public TrueType
{
 public:
   //**********************************************************************************************
   enum { value = 1 };
   typedef TrueType  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
//! Specialization of the IsVecVecMultExprHelper type trait for SVecDVecMultExpr.
template< typename VT1, typename VT2, bool TF >
struct IsVecVecMultExprHelper< SVecDVecMultExpr<VT1,VT2,TF> > : public TrueType
{
 public:
   //**********************************************************************************************
   enum { value = 1 };
   typedef TrueType  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
//! Specialization of the IsVecVecMultExprHelper type trait for SVecSVecMultExpr.
template< typename VT1, typename VT2, bool TF >
struct IsVecVecMultExprHelper< SVecSVecMultExpr<VT1,VT2,TF> > : public TrueType
{
 public:
   //**********************************************************************************************
   enum { value = 1 };
   typedef TrueType  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
//! Specialization of the IsVecVecMultExprHelper type trait for SVecTDVecMultExpr.
template< typename VT1, typename VT2 >
struct IsVecVecMultExprHelper< SVecTDVecMultExpr<VT1,VT2> > : public TrueType
{
 public:
   //**********************************************************************************************
   enum { value = 1 };
   typedef TrueType  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
//! Specialization of the IsVecVecMultExprHelper type trait for SVecTSVecMultExpr.
template< typename VT1, typename VT2 >
struct IsVecVecMultExprHelper< SVecTSVecMultExpr<VT1,VT2> > : public TrueType
{
 public:
   //**********************************************************************************************
   enum { value = 1 };
   typedef TrueType  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compile time check for expression types.
// \ingroup math_type_traits
//
// This type trait tests whether or not the given template parameter is a type representing a
// vector-vector multiplication. In case the type is a vector-vector multiplication expression,
// the \a value member enumeration is set to 1, the nested type definition \a Type is \a TrueType,
// and the class derives from \a TrueType. Otherwise \a value is set to 0, \a Type is \a FalseType,
// and the class derives from \a FalseType.
*/
template< typename T >
struct IsVecVecMultExpr : public IsVecVecMultExprHelper< typename boost::remove_cv<T>::type >::Type
{
 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   enum { value = IsVecVecMultExprHelper< typename boost::remove_cv<T>::type >::value };
   typedef typename IsVecVecMultExprHelper< typename boost::remove_cv<T>::type >::Type  Type;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************

} // namespace blaze

#endif
