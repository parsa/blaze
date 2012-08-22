//=================================================================================================
/*!
//  \file blaze/math/typetraits/IsMatMatMultExpr.h
//  \brief Header file for the IsMatMatMultExpr type trait
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

#ifndef _BLAZE_MATH_TYPETRAITS_ISMATMATMULTEXPR_H_
#define _BLAZE_MATH_TYPETRAITS_ISMATMATMULTEXPR_H_


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
/*!\brief Auxiliary helper struct for the IsMatMatMultExpr type trait.
// \ingroup math_type_traits
*/
template< typename T >
struct IsMatMatMultExprHelper
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
//! Specialization of the IsMatMatMultExprHelper type trait for DMatDMatMultExpr.
template< typename MT1, typename MT2 >
struct IsMatMatMultExprHelper< DMatDMatMultExpr<MT1,MT2> > : public TrueType
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
//! Specialization of the IsMatMatMultExprHelper type trait for DMatTDMatMultExpr.
template< typename MT1, typename MT2 >
struct IsMatMatMultExprHelper< DMatTDMatMultExpr<MT1,MT2> > : public TrueType
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
//! Specialization of the IsMatMatMultExprHelper type trait for TDMatDMatMultExpr.
template< typename MT1, typename MT2 >
struct IsMatMatMultExprHelper< TDMatDMatMultExpr<MT1,MT2> > : public TrueType
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
//! Specialization of the IsMatMatMultExprHelper type trait for TDMatTDMatMultExpr.
template< typename MT1, typename MT2 >
struct IsMatMatMultExprHelper< TDMatTDMatMultExpr<MT1,MT2> > : public TrueType
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
//! Specialization of the IsMatMatMultExprHelper type trait for DMatSMatMultExpr.
template< typename MT1, typename MT2 >
struct IsMatMatMultExprHelper< DMatSMatMultExpr<MT1,MT2> > : public TrueType
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
//! Specialization of the IsMatMatMultExprHelper type trait for DMatTSMatMultExpr.
template< typename MT1, typename MT2 >
struct IsMatMatMultExprHelper< DMatTSMatMultExpr<MT1,MT2> > : public TrueType
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
//! Specialization of the IsMatMatMultExprHelper type trait for TDMatSMatMultExpr.
template< typename MT1, typename MT2 >
struct IsMatMatMultExprHelper< TDMatSMatMultExpr<MT1,MT2> > : public TrueType
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
//! Specialization of the IsMatMatMultExprHelper type trait for TDMatTSMatMultExpr.
template< typename MT1, typename MT2 >
struct IsMatMatMultExprHelper< TDMatTSMatMultExpr<MT1,MT2> > : public TrueType
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
//! Specialization of the IsMatMatMultExprHelper type trait for SMatDMatMultExpr.
template< typename MT1, typename MT2 >
struct IsMatMatMultExprHelper< SMatDMatMultExpr<MT1,MT2> > : public TrueType
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
//! Specialization of the IsMatMatMultExprHelper type trait for SMatTDMatMultExpr.
template< typename MT1, typename MT2 >
struct IsMatMatMultExprHelper< SMatTDMatMultExpr<MT1,MT2> > : public TrueType
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
//! Specialization of the IsMatMatMultExprHelper type trait for TSMatDMatMultExpr.
template< typename MT1, typename MT2 >
struct IsMatMatMultExprHelper< TSMatDMatMultExpr<MT1,MT2> > : public TrueType
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
//! Specialization of the IsMatMatMultExprHelper type trait for TSMatTDMatMultExpr.
template< typename MT1, typename MT2 >
struct IsMatMatMultExprHelper< TSMatTDMatMultExpr<MT1,MT2> > : public TrueType
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
//! Specialization of the IsMatMatMultExprHelper type trait for SMatSMatMultExpr.
template< typename MT1, typename MT2 >
struct IsMatMatMultExprHelper< SMatSMatMultExpr<MT1,MT2> > : public TrueType
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
//! Specialization of the IsMatMatMultExprHelper type trait for SMatTSMatMultExpr.
template< typename MT1, typename MT2 >
struct IsMatMatMultExprHelper< SMatTSMatMultExpr<MT1,MT2> > : public TrueType
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
//! Specialization of the IsMatMatMultExprHelper type trait for TSMatSMatMultExpr.
template< typename MT1, typename MT2 >
struct IsMatMatMultExprHelper< TSMatSMatMultExpr<MT1,MT2> > : public TrueType
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
//! Specialization of the IsMatMatMultExprHelper type trait for TSMatTSMatMultExpr.
template< typename MT1, typename MT2 >
struct IsMatMatMultExprHelper< TSMatTSMatMultExpr<MT1,MT2> > : public TrueType
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
// matrix-matrix multiplication. In case the type is a matrix-matrix multiplication expression,
// the \a value member enumeration is set to 1, the nested type definition \a Type is \a TrueType,
// and the class derives from \a TrueType. Otherwise \a value is set to 0, \a Type is \a FalseType,
// and the class derives from \a FalseType.
*/
template< typename T >
struct IsMatMatMultExpr : public IsMatMatMultExprHelper< typename boost::remove_cv<T>::type >::Type
{
 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   enum { value = IsMatMatMultExprHelper< typename boost::remove_cv<T>::type >::value };
   typedef typename IsMatMatMultExprHelper< typename boost::remove_cv<T>::type >::Type  Type;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************

} // namespace blaze

#endif
