//=================================================================================================
/*!
//  \file blaze/math/typetraits/IsMatMatAddExpr.h
//  \brief Header file for the IsMatMatAddExpr type trait
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

#ifndef _BLAZE_MATH_TYPETRAITS_ISMATMATADDEXPR_H_
#define _BLAZE_MATH_TYPETRAITS_ISMATMATADDEXPR_H_


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
/*!\brief Auxiliary helper struct for the IsMatMatAddExpr type trait.
// \ingroup math_type_traits
*/
template< typename T >
struct IsMatMatAddExprHelper
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
//! Specialization of the IsMatMatAddExprHelper type trait for DMatDMatAddExpr.
template< typename MT1, typename MT2, bool SO >
struct IsMatMatAddExprHelper< DMatDMatAddExpr<MT1,MT2,SO> > : public TrueType
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
//! Specialization of the IsMatMatAddExprHelper type trait for DMatSMatAddExpr.
template< typename MT1, typename MT2, bool SO >
struct IsMatMatAddExprHelper< DMatSMatAddExpr<MT1,MT2,SO> > : public TrueType
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
//! Specialization of the IsMatMatAddExprHelper type trait for DMatTDMatAddExpr.
template< typename MT1, typename MT2 >
struct IsMatMatAddExprHelper< DMatTDMatAddExpr<MT1,MT2> > : public TrueType
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
//! Specialization of the IsMatMatAddExprHelper type trait for DMatTSMatAddExpr.
template< typename MT1, typename MT2 >
struct IsMatMatAddExprHelper< DMatTSMatAddExpr<MT1,MT2> > : public TrueType
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
//! Specialization of the IsMatMatAddExprHelper type trait for SMatSMatAddExpr.
template< typename MT1, typename MT2 >
struct IsMatMatAddExprHelper< SMatSMatAddExpr<MT1,MT2> > : public TrueType
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
//! Specialization of the IsMatMatAddExprHelper type trait for SMatTSMatAddExpr.
template< typename MT1, typename MT2 >
struct IsMatMatAddExprHelper< SMatTSMatAddExpr<MT1,MT2> > : public TrueType
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
//! Specialization of the IsMatMatAddExprHelper type trait for TDMatSMatAddExpr.
template< typename MT1, typename MT2 >
struct IsMatMatAddExprHelper< TDMatSMatAddExpr<MT1,MT2> > : public TrueType
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
//! Specialization of the IsMatMatAddExprHelper type trait for TSMatTSMatAddExpr.
template< typename MT1, typename MT2 >
struct IsMatMatAddExprHelper< TSMatTSMatAddExpr<MT1,MT2> > : public TrueType
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
// matrix-matrix addition. In case the type is a matrix-matrix addition expression, the \a value
// member enumeration is set to 1, the nested type definition \a Type is \a TrueType, and the
// class derives from \a TrueType. Otherwise \a value is set to 0, \a Type is \a FalseType, and
// the class derives from \a FalseType.
*/
template< typename T >
struct IsMatMatAddExpr : public IsMatMatAddExprHelper< typename boost::remove_cv<T>::type >::Type
{
 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   enum { value = IsMatMatAddExprHelper< typename boost::remove_cv<T>::type >::value };
   typedef typename IsMatMatAddExprHelper< typename boost::remove_cv<T>::type >::Type  Type;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************

} // namespace blaze

#endif
