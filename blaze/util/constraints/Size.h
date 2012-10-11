//=================================================================================================
/*!
//  \file blaze/util/constraints/Size.h
//  \brief Constraint on the size of a data type
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

#ifndef _BLAZE_UTIL_CONSTRAINTS_SIZE_H_
#define _BLAZE_UTIL_CONSTRAINTS_SIZE_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/constraints/ConstraintTest.h>
#include <blaze/util/Suffix.h>
#include <blaze/util/typetraits/HasSize.h>


namespace blaze {

//=================================================================================================
//
//  MUST_HAVE_SIZE CONSTRAINT
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Compile time constraint.
// \ingroup constraints
//
// Helper template class for the compile time constraint enforcement. Based on the compile time
// constant expression used for the template instantiation, either the undefined basic template
// or the specialization is selected. If the undefined basic template is selected, a compilation
// error is created.
*/
template< bool > struct CONSTRAINT_MUST_HAVE_SIZE_FAILED;
template<> struct CONSTRAINT_MUST_HAVE_SIZE_FAILED<true> { enum { value = 1 }; };
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constraint on the size of a data type.
// \ingroup constraints
//
// In case the type \a T doesn't have a size of \a S bytes, a compilation error is created.
*/
#define BLAZE_CONSTRAINT_MUST_HAVE_SIZE(T,S) \
   typedef \
      ::blaze::CONSTRAINT_TEST< \
         ::blaze::CONSTRAINT_MUST_HAVE_SIZE_FAILED< ::blaze::HasSize<T,S>::value >::value > \
      BLAZE_JOIN( CONSTRAINT_MUST_HAVE_SIZE_TYPEDEF, __LINE__ )
//*************************************************************************************************




//=================================================================================================
//
//  MUST_NOT_HAVE_SIZE CONSTRAINT
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Compile time constraint.
// \ingroup constraints
//
// Helper template class for the compile time constraint enforcement. Based on the compile time
// constant expression used for the template instantiation, either the undefined basic template
// or the specialization is selected. If the undefined basic template is selected, a compilation
// error is created.
*/
template< bool > struct CONSTRAINT_MUST_NOT_HAVE_SIZE_FAILED;
template<> struct CONSTRAINT_MUST_NOT_HAVE_SIZE_FAILED<true> { enum { value = 1 }; };
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constraint on the size of a data type.
// \ingroup constraints
//
// In case the type \a T has a size of \a S bytes, a compilation error is created.
*/
#define BLAZE_CONSTRAINT_MUST_NOT_HAVE_SIZE(T,S) \
   typedef \
      ::blaze::CONSTRAINT_TEST< \
         ::blaze::CONSTRAINT_MUST_NOT_HAVE_SIZE_FAILED< !::blaze::HasSize<T,S>::value >::value > \
      BLAZE_JOIN( CONSTRAINT_MUST_NOT_HAVE_SIZE_TYPEDEF, __LINE__ )
//*************************************************************************************************




//=================================================================================================
//
//  MUST_HAVE_1_BYTE CONSTRAINT
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Compile time constraint.
// \ingroup constraints
//
// Helper template class for the compile time constraint enforcement. Based on the compile time
// constant expression used for the template instantiation, either the undefined basic template
// or the specialization is selected. If the undefined basic template is selected, a compilation
// error is created.
*/
template< bool > struct CONSTRAINT_MUST_HAVE_1_BYTE_FAILED;
template<> struct CONSTRAINT_MUST_HAVE_1_BYTE_FAILED<true> { enum { value = 1 }; };
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constraint on the size of a data type.
// \ingroup constraints
//
// In case the type \a T doesn't have a size of exactly 1 byte, a compilation error is created.
*/
#define BLAZE_CONSTRAINT_MUST_HAVE_1_BYTE(T) \
   typedef \
      ::blaze::CONSTRAINT_TEST< \
         ::blaze::CONSTRAINT_MUST_HAVE_1_BYTE_FAILED< ::blaze::Has1Byte<T>::value >::value > \
      BLAZE_JOIN( CONSTRAINT_MUST_HAVE_1_BYTE_TYPEDEF, __LINE__ )
//*************************************************************************************************




//=================================================================================================
//
//  MUST_NOT_HAVE_1_BYTE CONSTRAINT
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Compile time constraint.
// \ingroup constraints
//
// Helper template class for the compile time constraint enforcement. Based on the compile time
// constant expression used for the template instantiation, either the undefined basic template
// or the specialization is selected. If the undefined basic template is selected, a compilation
// error is created.
*/
template< bool > struct CONSTRAINT_MUST_NOT_HAVE_1_BYTE_FAILED;
template<> struct CONSTRAINT_MUST_NOT_HAVE_1_BYTE_FAILED<true> { enum { value = 1 }; };
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constraint on the size of a data type.
// \ingroup constraints
//
// In case the type \a T has a size of exactly 1 byte, a compilation error is created.
*/
#define BLAZE_CONSTRAINT_MUST_NOT_HAVE_1_BYTE(T) \
   typedef \
      ::blaze::CONSTRAINT_TEST< \
         ::blaze::CONSTRAINT_MUST_NOT_HAVE_1_BYTE_FAILED< !::blaze::Has1Byte<T>::value >::value > \
      BLAZE_JOIN( CONSTRAINT_MUST_NOT_HAVE_1_BYTE_TYPEDEF, __LINE__ )
//*************************************************************************************************




//=================================================================================================
//
//  MUST_HAVE_2_BYTES CONSTRAINT
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Compile time constraint.
// \ingroup constraints
//
// Helper template class for the compile time constraint enforcement. Based on the compile time
// constant expression used for the template instantiation, either the undefined basic template
// or the specialization is selected. If the undefined basic template is selected, a compilation
// error is created.
*/
template< bool > struct CONSTRAINT_MUST_HAVE_2_BYTES_FAILED;
template<> struct CONSTRAINT_MUST_HAVE_2_BYTES_FAILED<true> { enum { value = 1 }; };
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constraint on the size of a data type.
// \ingroup constraints
//
// In case the type \a T doesn't have a size of exactly 2 bytes, a compilation error is created.
*/
#define BLAZE_CONSTRAINT_MUST_HAVE_2_BYTES(T) \
   typedef \
      ::blaze::CONSTRAINT_TEST< \
         ::blaze::CONSTRAINT_MUST_HAVE_2_BYTES_FAILED< ::blaze::Has2Bytes<T>::value >::value > \
      BLAZE_JOIN( CONSTRAINT_MUST_HAVE_2_BYTES_TYPEDEF, __LINE__ )
//*************************************************************************************************




//=================================================================================================
//
//  MUST_NOT_HAVE_2_BYTES CONSTRAINT
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Compile time constraint.
// \ingroup constraints
//
// Helper template class for the compile time constraint enforcement. Based on the compile time
// constant expression used for the template instantiation, either the undefined basic template
// or the specialization is selected. If the undefined basic template is selected, a compilation
// error is created.
*/
template< bool > struct CONSTRAINT_MUST_NOT_HAVE_2_BYTES_FAILED;
template<> struct CONSTRAINT_MUST_NOT_HAVE_2_BYTES_FAILED<true> { enum { value = 1 }; };
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constraint on the size of a data type.
// \ingroup constraints
//
// In case the type \a T has a size of exactly 2 bytes, a compilation error is created.
*/
#define BLAZE_CONSTRAINT_MUST_NOT_HAVE_2_BYTES(T) \
   typedef \
      ::blaze::CONSTRAINT_TEST< \
         ::blaze::CONSTRAINT_MUST_NOT_HAVE_2_BYTES_FAILED< !::blaze::Has2Bytes<T>::value >::value > \
      BLAZE_JOIN( CONSTRAINT_MUST_NOT_HAVE_2_BYTES_TYPEDEF, __LINE__ )
//*************************************************************************************************




//=================================================================================================
//
//  MUST_HAVE_4_BYTES CONSTRAINT
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Compile time constraint.
// \ingroup constraints
//
// Helper template class for the compile time constraint enforcement. Based on the compile time
// constant expression used for the template instantiation, either the undefined basic template
// or the specialization is selected. If the undefined basic template is selected, a compilation
// error is created.
*/
template< bool > struct CONSTRAINT_MUST_HAVE_4_BYTES_FAILED;
template<> struct CONSTRAINT_MUST_HAVE_4_BYTES_FAILED<true> { enum { value = 1 }; };
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constraint on the size of a data type.
// \ingroup constraints
//
// In case the type \a T doesn't have a size of exactly 4 bytes, a compilation error is created.
*/
#define BLAZE_CONSTRAINT_MUST_HAVE_4_BYTES(T) \
   typedef \
      ::blaze::CONSTRAINT_TEST< \
         ::blaze::CONSTRAINT_MUST_HAVE_4_BYTES_FAILED< ::blaze::Has4Bytes<T>::value >::value > \
      BLAZE_JOIN( CONSTRAINT_MUST_HAVE_4_BYTES_TYPEDEF, __LINE__ )
//*************************************************************************************************




//=================================================================================================
//
//  MUST_NOT_HAVE_4_BYTES CONSTRAINT
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Compile time constraint.
// \ingroup constraints
//
// Helper template class for the compile time constraint enforcement. Based on the compile time
// constant expression used for the template instantiation, either the undefined basic template
// or the specialization is selected. If the undefined basic template is selected, a compilation
// error is created.
*/
template< bool > struct CONSTRAINT_MUST_NOT_HAVE_4_BYTES_FAILED;
template<> struct CONSTRAINT_MUST_NOT_HAVE_4_BYTES_FAILED<true> { enum { value = 1 }; };
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constraint on the size of a data type.
// \ingroup constraints
//
// In case the type \a T has a size of exactly 4 bytes, a compilation error is created.
*/
#define BLAZE_CONSTRAINT_MUST_NOT_HAVE_4_BYTES(T) \
   typedef \
      ::blaze::CONSTRAINT_TEST< \
         ::blaze::CONSTRAINT_MUST_NOT_HAVE_4_BYTES_FAILED< !::blaze::Has4Bytes<T>::value >::value > \
      BLAZE_JOIN( CONSTRAINT_MUST_NOT_HAVE_4_BYTES_TYPEDEF, __LINE__ )
//*************************************************************************************************




//=================================================================================================
//
//  MUST_HAVE_8_BYTES CONSTRAINT
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Compile time constraint.
// \ingroup constraints
//
// Helper template class for the compile time constraint enforcement. Based on the compile time
// constant expression used for the template instantiation, either the undefined basic template
// or the specialization is selected. If the undefined basic template is selected, a compilation
// error is created.
*/
template< bool > struct CONSTRAINT_MUST_HAVE_8_BYTES_FAILED;
template<> struct CONSTRAINT_MUST_HAVE_8_BYTES_FAILED<true> { enum { value = 1 }; };
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constraint on the size of a data type.
// \ingroup constraints
//
// In case the type \a T doesn't have a size of exactly 8 bytes, a compilation error is created.
*/
#define BLAZE_CONSTRAINT_MUST_HAVE_8_BYTES(T) \
   typedef \
      ::blaze::CONSTRAINT_TEST< \
         ::blaze::CONSTRAINT_MUST_HAVE_8_BYTES_FAILED< ::blaze::Has8Bytes<T>::value >::value > \
      BLAZE_JOIN( CONSTRAINT_MUST_HAVE_8_BYTES_TYPEDEF, __LINE__ )
//*************************************************************************************************




//=================================================================================================
//
//  MUST_NOT_HAVE_8_BYTES CONSTRAINT
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Compile time constraint.
// \ingroup constraints
//
// Helper template class for the compile time constraint enforcement. Based on the compile time
// constant expression used for the template instantiation, either the undefined basic template
// or the specialization is selected. If the undefined basic template is selected, a compilation
// error is created.
*/
template< bool > struct CONSTRAINT_MUST_NOT_HAVE_8_BYTES_FAILED;
template<> struct CONSTRAINT_MUST_NOT_HAVE_8_BYTES_FAILED<true> { enum { value = 1 }; };
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constraint on the size of a data type.
// \ingroup constraints
//
// In case the type \a T has a size of exactly 8 bytes, a compilation error is created.
*/
#define BLAZE_CONSTRAINT_MUST_NOT_HAVE_8_BYTES(T) \
   typedef \
      ::blaze::CONSTRAINT_TEST< \
         ::blaze::CONSTRAINT_MUST_NOT_HAVE_8_BYTES_FAILED< !::blaze::Has8Bytes<T>::value >::value > \
      BLAZE_JOIN( CONSTRAINT_MUST_NOT_HAVE_8_BYTES_TYPEDEF, __LINE__ )
//*************************************************************************************************

} // namespace blaze

#endif
