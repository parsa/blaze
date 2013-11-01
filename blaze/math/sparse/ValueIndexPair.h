//=================================================================================================
/*!
//  \file blaze/math/sparse/ValueIndexPair.h
//  \brief Header file for the ValueIndexPair class
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

#ifndef _BLAZE_MATH_SPARSE_VALUEINDEXPAIR_H_
#define _BLAZE_MATH_SPARSE_VALUEINDEXPAIR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/constraints/Const.h>
#include <blaze/util/constraints/Pointer.h>
#include <blaze/util/constraints/Reference.h>
#include <blaze/util/constraints/Volatile.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Index-value-pair for sparse vectors and matrices.
// \ingroup math
//
// The ValueIndexPair class represents a single index-value-pair of a sparse vector or sparse
// matrix.
*/
template< typename Type >  // Type of the value element
class ValueIndexPair
{
 public:
   //**Type definitions****************************************************************************
   typedef Type    ValueType;  //!< The value type of the value-index-pair.
   typedef size_t  IndexType;  //!< The index type of the value-index-pair.
   //**********************************************************************************************

   //**Constructors********************************************************************************
   inline ValueIndexPair();
   inline ValueIndexPair( const Type& v, size_t i );
   // No explicitly declared copy constructor.
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   // No explicitly declared destructor.
   //**********************************************************************************************

   //**Operators***********************************************************************************
   /*!\name Operators */
   //@{
   // No explicitly declared copy assignment operator.
   template< typename Other > ValueIndexPair& operator=( const Other& rhs );

   inline ValueIndexPair& operator= ( const Type& v );
   inline ValueIndexPair& operator+=( const Type& v );
   inline ValueIndexPair& operator-=( const Type& v );
   inline ValueIndexPair& operator*=( const Type& v );
   inline ValueIndexPair& operator/=( const Type& v );
   //@}
   //**********************************************************************************************

   //**Acess functions*****************************************************************************
   /*!\name Access functions */
   //@{
   inline Type&       value();
   inline const Type& value() const;
   inline size_t      index() const;
   //@}
   //**********************************************************************************************

 protected:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   Type   value_;  //!< Value of the value-index-pair.
   size_t index_;  //!< Index of the value-index-pair.
   //@}
   //**********************************************************************************************

 private:
   //**Friend declarations*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   template< typename Other > friend class ValueIndexPair;
   /*! \endcond */
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_NOT_BE_POINTER_TYPE  ( Type );
   BLAZE_CONSTRAINT_MUST_NOT_BE_REFERENCE_TYPE( Type );
   BLAZE_CONSTRAINT_MUST_NOT_BE_CONST         ( Type );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VOLATILE      ( Type );
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Default constructor for value-index-pairs.
*/
template< typename Type >  // Type of the value element
inline ValueIndexPair<Type>::ValueIndexPair()
   : value_()  // Value of the value-index-pair
   , index_()  // Index of the value-index-pair
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constructor for a direct initialization of value-index-pairs.
//
// \param v The value of the value-index-pair.
// \param i The index of the value-index-pair.
*/
template< typename Type >  // Type of the value element
inline ValueIndexPair<Type>::ValueIndexPair( const Type& v, size_t i )
   : value_( v )  // Value of the value-index-pair
   , index_( i )  // Index of the value-index-pair
{}
//*************************************************************************************************




//=================================================================================================
//
//  OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Assignment operator for different value-index-pair types.
//
// \param rhs Value-index-pair to be copied.
// \return Reference to the assigned value-index-pair.
//
// This assignment operator enables the assignment of other value-index-pair types. The given
// \a Other data type qualifies as value-index-pair type in case it provides a value() and an
// index() member function.
*/
template< typename Type >   // Type of the value element
template< typename Other >  // Data type of the right-hand side value-index-pair
inline ValueIndexPair<Type>& ValueIndexPair<Type>::operator=( const Other& rhs )
{
   value_ = rhs.value();
   index_ = rhs.index();
   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Assignment to the value of the value-index-pair.
//
// \param v The new value-index-pair value.
// \return Reference to the assigned value-index-pair.
*/
template< typename Type >  // Type of the value element
inline ValueIndexPair<Type>& ValueIndexPair<Type>::operator=( const Type& v )
{
   value_ = v;
   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Addition assignment to the value of the value-index-pair.
//
// \param v The right-hand side value to be added to the value-index-pair value.
// \return Reference to the assigned value-index-pair.
*/
template< typename Type >  // Type of the value element
inline ValueIndexPair<Type>& ValueIndexPair<Type>::operator+=( const Type& v )
{
   value_ += v;
   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Subtraction assignment to the value of the value-index-pair.
//
// \param v The right-hand side value to be subtracted from the value-index-pair value.
// \return Reference to the assigned value-index-pair.
*/
template< typename Type >  // Type of the value element
inline ValueIndexPair<Type>& ValueIndexPair<Type>::operator-=( const Type& v )
{
   value_ -= v;
   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication assignment to the value of the value-index-pair.
//
// \param v The right-hand side value for the multiplication.
// \return Reference to the assigned value-index-pair.
*/
template< typename Type >  // Type of the value element
inline ValueIndexPair<Type>& ValueIndexPair<Type>::operator*=( const Type& v )
{
   value_ *= v;
   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Division assignment to the value of the value-index-pair.
//
// \param v The right-hand side value for the division
// \return Reference to the assigned value-index-pair.
*/
template< typename Type >  // Type of the value element
inline ValueIndexPair<Type>& ValueIndexPair<Type>::operator/=( const Type& v )
{
   value_ /= v;
   return *this;
}
//*************************************************************************************************




//=================================================================================================
//
//  ACCESS FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Access to the current value of the value-index-pair.
//
// \return The current value of the value-index-pair.
*/
template< typename Type >  // Type of the value element
inline Type& ValueIndexPair<Type>::value()
{
   return value_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Access to the current value of the value-index-pair.
//
// \return The current value of the value-index-pair.
*/
template< typename Type >  // Type of the value element
inline const Type& ValueIndexPair<Type>::value() const
{
   return value_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Access to the current index of the value-index-pair.
//
// \return The current index of the value-index-pair.
*/
template< typename Type >  // Type of the value element
inline size_t ValueIndexPair<Type>::index() const
{
   return index_;
}
//*************************************************************************************************

} // namespace blaze

#endif
