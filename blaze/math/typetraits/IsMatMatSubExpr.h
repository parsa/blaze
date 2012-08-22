//=================================================================================================
/*!
//  \file blaze/math/typetraits/IsMatMatSubExpr.h
//  \brief Header file for the IsMatMatSubExpr type trait
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

#ifndef _BLAZE_MATH_TYPETRAITS_ISMATMATSUBEXPR_H_
#define _BLAZE_MATH_TYPETRAITS_ISMATMATSUBEXPR_H_


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
/*!\brief Auxiliary helper struct for the IsMatMatSubExpr type trait.
// \ingroup math_type_traits
*/
template< typename T >
struct IsMatMatSubExprHelper
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
//! Specialization of the IsMatMatSubExprHelper type trait for DMatDMatSubExpr.
template< typename MT1, typename MT2, bool SO >
struct IsMatMatSubExprHelper< DMatDMatSubExpr<MT1,MT2,SO> > : public TrueType
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
//! Specialization of the IsMatMatSubExprHelper type trait for DMatSMatSubExpr.
template< typename MT1, typename MT2, bool SO >
struct IsMatMatSubExprHelper< DMatSMatSubExpr<MT1,MT2,SO> > : public TrueType
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
//! Specialization of the IsMatMatSubExprHelper type trait for DMatTDMatSubExpr.
template< typename MT1, typename MT2 >
struct IsMatMatSubExprHelper< DMatTDMatSubExpr<MT1,MT2> > : public TrueType
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
//! Specialization of the IsMatMatSubExprHelper type trait for DMatTSMatSubExpr.
template< typename MT1, typename MT2 >
struct IsMatMatSubExprHelper< DMatTSMatSubExpr<MT1,MT2> > : public TrueType
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
//! Specialization of the IsMatMatSubExprHelper type trait for SMatDMatSubExpr.
template< typename MT1, typename MT2, bool SO >
struct IsMatMatSubExprHelper< SMatDMatSubExpr<MT1,MT2,SO> > : public TrueType
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
//! Specialization of the IsMatMatSubExprHelper type trait for SMatSMatSubExpr.
template< typename MT1, typename MT2 >
struct IsMatMatSubExprHelper< SMatSMatSubExpr<MT1,MT2> > : public TrueType
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
//! Specialization of the IsMatMatSubExprHelper type trait for SMatTDMatSubExpr.
template< typename MT1, typename MT2 >
struct IsMatMatSubExprHelper< SMatTDMatSubExpr<MT1,MT2> > : public TrueType
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
//! Specialization of the IsMatMatSubExprHelper type trait for SMatTSMatSubExpr.
template< typename MT1, typename MT2 >
struct IsMatMatSubExprHelper< SMatTSMatSubExpr<MT1,MT2> > : public TrueType
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
//! Specialization of the IsMatMatSubExprHelper type trait for TDMatSMatSubExpr.
template< typename MT1, typename MT2 >
struct IsMatMatSubExprHelper< TDMatSMatSubExpr<MT1,MT2> > : public TrueType
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
//! Specialization of the IsMatMatSubExprHelper type trait for TSMatDMatSubExpr.
template< typename MT1, typename MT2 >
struct IsMatMatSubExprHelper< TSMatDMatSubExpr<MT1,MT2> > : public TrueType
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
//! Specialization of the IsMatMatSubExprHelper type trait for TSMatSMatSubExpr.
template< typename MT1, typename MT2 >
struct IsMatMatSubExprHelper< TSMatSMatSubExpr<MT1,MT2> > : public TrueType
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
//! Specialization of the IsMatMatSubExprHelper type trait for TSMatTSMatSubExpr.
template< typename MT1, typename MT2 >
struct IsMatMatSubExprHelper< TSMatTSMatSubExpr<MT1,MT2> > : public TrueType
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
// matrix-matrix subtraction. In case the type is a matrix-matrix subtraction expression, the
// \a value member enumeration is set to 1, the nested type definition \a Type is \a TrueType,
// and the class derives from \a TrueType. Otherwise \a value is set to 0, \a Type is \a FalseType,
// and the class derives from \a FalseType.
*/
template< typename T >
struct IsMatMatSubExpr : public IsMatMatSubExprHelper< typename boost::remove_cv<T>::type >::Type
{
 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   enum { value = IsMatMatSubExprHelper< typename boost::remove_cv<T>::type >::value };
   typedef typename IsMatMatSubExprHelper< typename boost::remove_cv<T>::type >::Type  Type;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************

} // namespace blaze

#endif
