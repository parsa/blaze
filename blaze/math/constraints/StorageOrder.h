//=================================================================================================
/*!
//  \file blaze/math/constraints/StorageOrder.h
//  \brief Constraints on the storage order of matrix types
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

#ifndef _BLAZE_MATH_CONSTRAINTS_STORAGEORDER_H_
#define _BLAZE_MATH_CONSTRAINTS_STORAGEORDER_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/typetraits/IsColumnMajorMatrix.h>
#include <blaze/math/typetraits/IsMatrix.h>
#include <blaze/math/typetraits/IsRowMajorMatrix.h>
#include <blaze/util/constraints/ConstraintTest.h>
#include <blaze/util/Suffix.h>


namespace blaze {

//=================================================================================================
//
//  MUST_BE_ROW_MAJOR_MATRIX_TYPE CONSTRAINT
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Compile time constraint.
// \ingroup math_constraints
//
// Helper template class for the compile time constraint enforcement. Based on the compile time
// constant expression used for the template instantiation, either the undefined basic template
// or the specialization is selected. If the undefined basic template is selected, a compilation
// error is created.
*/
template< bool > struct CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE_FAILED;
template<> struct CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE_FAILED<true> { enum { value = 1 }; };
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constraint on the data type.
// \ingroup math_constraints
//
// In case the given data type \a T is not a row-major dense or sparse matrix type (i.e. a matrix
// type whose storage order is set to \a false) a compilation error is created.
*/
#define BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE(T) \
   typedef \
      blaze::CONSTRAINT_TEST< \
         blaze::CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE_FAILED< \
            blaze::IsRowMajorMatrix<T>::value >::value > \
      BLAZE_JOIN( CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE_TYPEDEF, __LINE__ )
//*************************************************************************************************




//=================================================================================================
//
//  MUST_NOT_BE_ROW_MAJOR_MATRIX_TYPE CONSTRAINT
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Compile time constraint.
// \ingroup math_constraints
//
// Helper template class for the compile time constraint enforcement. Based on the compile time
// constant expression used for the template instantiation, either the undefined basic template
// or the specialization is selected. If the undefined basic template is selected, a compilation
// error is created.
*/
template< bool > struct CONSTRAINT_MUST_NOT_BE_ROW_MAJOR_MATRIX_TYPE_FAILED;
template<> struct CONSTRAINT_MUST_NOT_BE_ROW_MAJOR_MATRIX_TYPE_FAILED<true> { enum { value = 1 }; };
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constraint on the data type.
// \ingroup math_constraints
//
// In case the given data type \a T is a row-major dense or sparse matrix type (i.e. a matrix
// type whose storage order is set to \a false) a compilation error is created.
*/
#define BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_MAJOR_MATRIX_TYPE(T) \
   typedef \
      blaze::CONSTRAINT_TEST< \
         blaze::CONSTRAINT_MUST_NOT_BE_ROW_MAJOR_MATRIX_TYPE_FAILED< \
            !blaze::IsRowMajorMatrix<T>::value >::value \
      BLAZE_JOIN( CONSTRAINT_MUST_NOT_BE_ROW_MAJOR_MATRIX_TYPE_TYPEDEF, __LINE__ )
//*************************************************************************************************




//=================================================================================================
//
//  MUST_BE_COLUMN_MAJOR_MATRIX_TYPE CONSTRAINT
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Compile time constraint.
// \ingroup math_constraints
//
// Helper template class for the compile time constraint enforcement. Based on the compile time
// constant expression used for the template instantiation, either the undefined basic template
// or the specialization is selected. If the undefined basic template is selected, a compilation
// error is created.
*/
template< bool > struct CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE_FAILED;
template<> struct CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE_FAILED<true> { enum { value = 1 }; };
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constraint on the data type.
// \ingroup math_constraints
//
// In case the given data type \a T is not a column-major dense or sparse matrix type (i.e. a
// matrix type whose storage order is set to \a true) a compilation error is created.
*/
#define BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE(T) \
   typedef \
      blaze::CONSTRAINT_TEST< \
         blaze::CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE_FAILED< \
            blaze::IsColumnMajorMatrix<T>::value >::value > \
      BLAZE_JOIN( CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE_TYPEDEF, __LINE__ )
//*************************************************************************************************




//=================================================================================================
//
//  MUST_NOT_BE_COLUMN_MAJOR_MATRIX_TYPE CONSTRAINT
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Compile time constraint.
// \ingroup math_constraints
//
// Helper template class for the compile time constraint enforcement. Based on the compile time
// constant expression used for the template instantiation, either the undefined basic template
// or the specialization is selected. If the undefined basic template is selected, a compilation
// error is created.
*/
template< bool > struct CONSTRAINT_MUST_NOT_BE_COLUMN_MAJOR_MATRIX_TYPE_FAILED;
template<> struct CONSTRAINT_MUST_NOT_BE_COLUMN_MAJOR_MATRIX_TYPE_FAILED<true> { enum { value = 1 }; };
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constraint on the data type.
// \ingroup math_constraints
//
// In case the given data type \a T is a column-major dense or sparse matrix type (i.e. a matrix
// type whose storage order is set to \a true) a compilation error is created.
*/
#define BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_MAJOR_MATRIX_TYPE(T) \
   typedef \
      blaze::CONSTRAINT_TEST< \
         blaze::CONSTRAINT_MUST_NOT_BE_COLUMN_MAJOR_MATRIX_TYPE_FAILED< \
            !blaze::IsColumnMajorMatrix<T>::value >::value > \
      BLAZE_JOIN( CONSTRAINT_MUST_NOT_BE_COLUMN_MAJOR_MATRIX_TYPE_TYPEDEF, __LINE__ )
//*************************************************************************************************




//=================================================================================================
//
//  MUST_BE_MATRIX_WITH_STORAGE_ORDER CONSTRAINT
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Compile time constraint.
// \ingroup math_constraints
//
// Helper template class for the compile time constraint enforcement. Based on the compile time
// constant expression used for the template instantiation, either the undefined basic template
// or the specialization is selected. If the undefined basic template is selected, a compilation
// error is created.
*/
template< bool > struct CONSTRAINT_MUST_BE_MATRIX_WITH_STORAGE_ORDER_FAILED;
template<> struct CONSTRAINT_MUST_BE_MATRIX_WITH_STORAGE_ORDER_FAILED<true> { enum { value = 1 }; };
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constraint on the data type.
// \ingroup math_constraints
//
// In case the given data type \a T is not a dense or sparse matrix type and in case the
// storage order of the given dense or sparse vector type \a T is not set to \a SO, a
// compilation error is created.
*/
#define BLAZE_CONSTRAINT_MUST_BE_MATRIX_WITH_STORAGE_ORDER(T,SO) \
   typedef \
      blaze::CONSTRAINT_TEST< \
         blaze::CONSTRAINT_MUST_BE_MATRIX_WITH_STORAGE_ORDER_FAILED< \
            blaze::IsMatrix<T>::value && \
            blaze::IsColumnMajorMatrix<T>::value == SO >::value > \
      BLAZE_JOIN( CONSTRAINT_MUST_BE_MATRIX_WITH_STORAGE_ORDER_TYPEDEF, __LINE__ )
//*************************************************************************************************




//=================================================================================================
//
//  MATRICES_MUST_HAVE_SAME_STORAGE_ORDER CONSTRAINT
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Compile time constraint.
// \ingroup math_constraints
//
// Helper template class for the compile time constraint enforcement. Based on the compile time
// constant expression used for the template instantiation, either the undefined basic template
// or the specialization is selected. If the undefined basic template is selected, a compilation
// error is created.
*/
template< bool > struct CONSTRAINT_MATRICES_MUST_HAVE_SAME_STORAGE_ORDER_FAILED;
template<> struct CONSTRAINT_MATRICES_MUST_HAVE_SAME_STORAGE_ORDER_FAILED<true> { enum { value = 1 }; };
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constraint on the data type.
// \ingroup math_constraints
//
// In case either of the two given data types \a T1 or \a T2 is not a matrix type and in case
// the storage order of both matrix types doesn't match, a compilation error is created.
*/
#define BLAZE_CONSTRAINT_MATRICES_MUST_HAVE_SAME_STORAGE_ORDER(T1,T2) \
   typedef \
      blaze::CONSTRAINT_TEST< \
         blaze::CONSTRAINT_MATRICES_MUST_HAVE_SAME_STORAGE_ORDER_FAILED< \
            blaze::IsMatrix<T1>::value && \
            blaze::IsMatrix<T2>::value && \
            static_cast<int>( blaze::IsRowMajorMatrix<T1>::value ) == static_cast<int>( blaze::IsRowMajorMatrix<T2>::value ) >::value > \
      BLAZE_JOIN( CONSTRAINT_MATRICES_MUST_HAVE_SAME_STORAGE_ORDER_TYPEDEF, __LINE__ )
//*************************************************************************************************




//=================================================================================================
//
//  MATRICES_MUST_HAVE_DIFFERENT_STORAGE_ORDER CONSTRAINT
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Compile time constraint.
// \ingroup math_constraints
//
// Helper template class for the compile time constraint enforcement. Based on the compile time
// constant expression used for the template instantiation, either the undefined basic template
// or the specialization is selected. If the undefined basic template is selected, a compilation
// error is created.
*/
template< bool > struct CONSTRAINT_MATRICES_MUST_HAVE_DIFFERENT_STORAGE_ORDER_FAILED;
template<> struct CONSTRAINT_MATRICES_MUST_HAVE_DIFFERENT_STORAGE_ORDER_FAILED<true> { enum { value = 1 }; };
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constraint on the data type.
// \ingroup math_constraints
//
// In case either of the two given data types \a T1 or \a T2 is not a matrix type and in case
// the storage order of both matrix types does match, a compilation error is created.
*/
#define BLAZE_CONSTRAINT_MATRICES_MUST_HAVE_DIFFERENT_STORAGE_ORDER(T1,T2) \
   typedef \
      blaze::CONSTRAINT_TEST< \
         blaze::CONSTRAINT_MATRICES_MUST_HAVE_DIFFERENT_STORAGE_ORDER_FAILED< \
            blaze::IsMatrix<T1>::value && \
            blaze::IsMatrix<T2>::value && \
            static_cast<int>( blaze::IsRowMajorMatrix<T1>::value ) != static_cast<int>( blaze::IsRowMajorMatrix<T2>::value ) >::value > \
      BLAZE_JOIN( CONSTRAINT_MATRICES_MUST_HAVE_DIFFERENT_STORAGE_ORDER_TYPEDEF, __LINE__ )
//*************************************************************************************************

} // namespace blaze

#endif
