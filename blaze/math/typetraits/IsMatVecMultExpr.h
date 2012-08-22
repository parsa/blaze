//=================================================================================================
/*!
//  \file blaze/math/typetraits/IsMatVecMultExpr.h
//  \brief Header file for the IsMatVecMultExpr type trait
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

#ifndef _BLAZE_MATH_TYPETRAITS_ISMATVECMULTEXPR_H_
#define _BLAZE_MATH_TYPETRAITS_ISMATVECMULTEXPR_H_


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
/*!\brief Auxiliary helper struct for the IsMatVecMultExpr type trait.
// \ingroup math_type_traits
*/
template< typename T >
struct IsMatVecMultExprHelper
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
//! Specialization of the IsMatVecMultExprHelper type trait for DMatDVecMultExpr.
template< typename MT, typename VT >
struct IsMatVecMultExprHelper< DMatDVecMultExpr<MT,VT> > : public TrueType
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
//! Specialization of the IsMatVecMultExprHelper type trait for TDMatDVecMultExpr.
template< typename MT, typename VT >
struct IsMatVecMultExprHelper< TDMatDVecMultExpr<MT,VT> > : public TrueType
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
//! Specialization of the IsMatVecMultExprHelper type trait for DMatSVecMultExpr.
template< typename MT, typename VT >
struct IsMatVecMultExprHelper< DMatSVecMultExpr<MT,VT> > : public TrueType
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
//! Specialization of the IsMatVecMultExprHelper type trait for TDMatSVecMultExpr.
template< typename MT, typename VT >
struct IsMatVecMultExprHelper< TDMatSVecMultExpr<MT,VT> > : public TrueType
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
//! Specialization of the IsMatVecMultExprHelper type trait for SMatDVecMultExpr.
template< typename MT, typename VT >
struct IsMatVecMultExprHelper< SMatDVecMultExpr<MT,VT> > : public TrueType
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
//! Specialization of the IsMatVecMultExprHelper type trait for TSMatDVecMultExpr.
template< typename MT, typename VT >
struct IsMatVecMultExprHelper< TSMatDVecMultExpr<MT,VT> > : public TrueType
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
//! Specialization of the IsMatVecMultExprHelper type trait for SMatSVecMultExpr.
template< typename MT, typename VT >
struct IsMatVecMultExprHelper< SMatSVecMultExpr<MT,VT> > : public TrueType
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
//! Specialization of the IsMatVecMultExprHelper type trait for TSMatSVecMultExpr.
template< typename MT, typename VT >
struct IsMatVecMultExprHelper< TSMatSVecMultExpr<MT,VT> > : public TrueType
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
// matrix-vector multiplication. In case the type is a matrix-vector multiplication expression,
// the \a value member enumeration is set to 1, the nested type definition \a Type is \a TrueType,
// and the class derives from \a TrueType. Otherwise \a value is set to 0, \a Type is \a FalseType,
// and the class derives from \a FalseType.
*/
template< typename T >
struct IsMatVecMultExpr : public IsMatVecMultExprHelper< typename boost::remove_cv<T>::type >::Type
{
 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   enum { value = IsMatVecMultExprHelper< typename boost::remove_cv<T>::type >::value };
   typedef typename IsMatVecMultExprHelper< typename boost::remove_cv<T>::type >::Type  Type;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************

} // namespace blaze

#endif
