//=================================================================================================
/*!
//  \file blaze/math/typetraits/IsTVecMatMultExpr.h
//  \brief Header file for the IsTVecMatMultExpr type trait
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

#ifndef _BLAZE_MATH_TYPETRAITS_ISTVECMATMULTEXPR_H_
#define _BLAZE_MATH_TYPETRAITS_ISTVECMATMULTEXPR_H_


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
/*!\brief Auxiliary helper struct for the IsTVecMatMultExpr type trait.
// \ingroup math_type_traits
*/
template< typename T >
struct IsTVecMatMultExprHelper
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
//! Specialization of the IsTVecMatMultExprHelper type trait for TDVecDMatMultExpr.
template< typename VT, typename MT >
struct IsTVecMatMultExprHelper< TDVecDMatMultExpr<VT,MT> > : public TrueType
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
//! Specialization of the IsTVecMatMultExprHelper type trait for TDVecTDMatMultExpr.
template< typename VT, typename MT >
struct IsTVecMatMultExprHelper< TDVecTDMatMultExpr<VT,MT> > : public TrueType
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
//! Specialization of the IsTVecMatMultExprHelper type trait for TSVecDMatMultExpr.
template< typename VT, typename MT >
struct IsTVecMatMultExprHelper< TSVecDMatMultExpr<VT,MT> > : public TrueType
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
//! Specialization of the IsTVecMatMultExprHelper type trait for TSVecTDMatMultExpr.
template< typename VT, typename MT >
struct IsTVecMatMultExprHelper< TSVecTDMatMultExpr<VT,MT> > : public TrueType
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
//! Specialization of the IsTVecMatMultExprHelper type trait for TDVecSMatMultExpr.
template< typename VT, typename MT >
struct IsTVecMatMultExprHelper< TDVecSMatMultExpr<VT,MT> > : public TrueType
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
//! Specialization of the IsTVecMatMultExprHelper type trait for TDVecTSMatMultExpr.
template< typename VT, typename MT >
struct IsTVecMatMultExprHelper< TDVecTSMatMultExpr<VT,MT> > : public TrueType
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
//! Specialization of the IsTVecMatMultExprHelper type trait for TSVecSMatMultExpr.
template< typename VT, typename MT >
struct IsTVecMatMultExprHelper< TSVecSMatMultExpr<VT,MT> > : public TrueType
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
//! Specialization of the IsTVecMatMultExprHelper type trait for TSVecTSMatMultExpr.
template< typename VT, typename MT >
struct IsTVecMatMultExprHelper< TSVecTSMatMultExpr<VT,MT> > : public TrueType
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
// vector-matrix multiplication. In case the type is a vector-matrix multiplication expression,
// the \a value member enumeration is set to 1, the nested type definition \a Type is \a TrueType,
// and the class derives from \a TrueType. Otherwise \a value is set to 0, \a Type is \a FalseType,
// and the class derives from \a FalseType.
*/
template< typename T >
struct IsTVecMatMultExpr : public IsTVecMatMultExprHelper< typename boost::remove_cv<T>::type >::Type
{
 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   enum { value = IsTVecMatMultExprHelper< typename boost::remove_cv<T>::type >::value };
   typedef typename IsTVecMatMultExprHelper< typename boost::remove_cv<T>::type >::Type  Type;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************

} // namespace blaze

#endif
