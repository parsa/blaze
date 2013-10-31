//=================================================================================================
/*!
//  \file blaze/math/sparse/SparseElement.h
//  \brief Header file for the SparseElement class
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

#ifndef _BLAZE_MATH_SPARSE_SPARSEELEMENT_H_
#define _BLAZE_MATH_SPARSE_SPARSEELEMENT_H_


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
/*!\brief Index-value-pair for the sparse vectors and matrices.
// \ingroup math
//
// The SparseElement class represents a single index-value-pair of a sparse vector or sparse
// matrix.
*/
template< typename Type >  // Type of the data element
class SparseElement
{
 public:
   //**Type definitions****************************************************************************
   typedef Type    ValueType;  //!< The value type of the sparse element.
   typedef size_t  IndexType;  //!< The index type of the sparse element.
   //**********************************************************************************************

   //**Constructors********************************************************************************
   inline SparseElement();
   inline SparseElement( const Type& v, size_t i );
   // No explicitly declared copy constructor.
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   // No explicitly declared destructor.
   //**********************************************************************************************

   //**Operators***********************************************************************************
   /*!\name Operators */
   //@{
   // No explicitly declared copy assignment operator.
   template< typename Other > SparseElement& operator=( const Other& rhs );

   inline SparseElement& operator= ( const Type& v );
   inline SparseElement& operator+=( const Type& v );
   inline SparseElement& operator-=( const Type& v );
   inline SparseElement& operator*=( const Type& v );
   inline SparseElement& operator/=( const Type& v );
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
   Type   value_;  //!< Value of the sparse element.
   size_t index_;  //!< Index of the sparse element.
   //@}
   //**********************************************************************************************

 private:
   //**Friend declarations*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   template< typename Other > friend class SparseElement;
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
/*!\brief Default constructor for sparse elements.
*/
template< typename Type >  // Data type of the sparse element
inline SparseElement<Type>::SparseElement()
   : value_()  // Value of the sparse element
   , index_()  // Index of the sparse element
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constructor for a direct initialization of sparse elements.
//
// \param v The value of the sparse element.
// \param i The index of the sparse element.
*/
template< typename Type >  // Data type of the sparse element
inline SparseElement<Type>::SparseElement( const Type& v, size_t i )
   : value_( v )  // Value of the sparse element
   , index_( i )  // Index of the sparse element
{}
//*************************************************************************************************




//=================================================================================================
//
//  OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Assignment operator for different sparse element types.
//
// \param rhs Sparse element to be copied.
// \return Reference to the assigned sparse element.
//
// This assignment operator enables the assignment of other sparse element types. The given
// \a Other data type qualifies as sparse element type in case it provides a value() and an
// index() member function.
*/
template< typename Type >   // Data type of the sparse element
template< typename Other >  // Data type of the right-hand side sparse element
inline SparseElement<Type>& SparseElement<Type>::operator=( const Other& rhs )
{
   value_ = rhs.value();
   index_ = rhs.index();
   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Assignment to the value of the sparse element.
//
// \param v The new sparse element value.
// \return Reference to the assigned sparse element.
*/
template< typename Type >  // Data type of the sparse element
inline SparseElement<Type>& SparseElement<Type>::operator=( const Type& v )
{
   value_ = v;
   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Addition assignment to the value of the sparse element.
//
// \param v The right-hand side value to be added to the sparse element value.
// \return Reference to the assigned sparse element.
*/
template< typename Type >  // Data type of the sparse element
inline SparseElement<Type>& SparseElement<Type>::operator+=( const Type& v )
{
   value_ += v;
   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Subtraction assignment to the value of the sparse element.
//
// \param v The right-hand side value to be subtracted from the sparse element value.
// \return Reference to the assigned sparse element.
*/
template< typename Type >  // Data type of the sparse element
inline SparseElement<Type>& SparseElement<Type>::operator-=( const Type& v )
{
   value_ -= v;
   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication assignment to the value of the sparse element.
//
// \param v The right-hand side value for the multiplication.
// \return Reference to the assigned sparse element.
*/
template< typename Type >  // Data type of the sparse element
inline SparseElement<Type>& SparseElement<Type>::operator*=( const Type& v )
{
   value_ *= v;
   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Division assignment to the value of the sparse element.
//
// \param v The right-hand side value for the division
// \return Reference to the assigned sparse element.
*/
template< typename Type >  // Data type of the sparse element
inline SparseElement<Type>& SparseElement<Type>::operator/=( const Type& v )
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
/*!\brief Access to the current value of the sparse element.
//
// \return The current value of the sparse element.
*/
template< typename Type >  // Data type of the sparse element
inline Type& SparseElement<Type>::value()
{
   return value_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Access to the current value of the sparse element.
//
// \return The current value of the sparse element.
*/
template< typename Type >  // Data type of the sparse element
inline const Type& SparseElement<Type>::value() const
{
   return value_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Access to the current index of the sparse element.
//
// \return The current index of the sparse element.
*/
template< typename Type >  // Data type of the sparse element
inline size_t SparseElement<Type>::index() const
{
   return index_;
}
//*************************************************************************************************

} // namespace blaze

#endif
