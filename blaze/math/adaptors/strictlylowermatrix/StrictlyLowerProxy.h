//=================================================================================================
/*!
//  \file blaze/math/adaptors/strictlylowermatrix/StrictlyLowerProxy.h
//  \brief Header file for the StrictlyLowerProxy class
//
//  Copyright (C) 2013 Klaus Iglberger - All Rights Reserved
//
//  This file is part of the Blaze library. You can redistribute it and/or modify it under
//  the terms of the New (Revised) BSD License. Redistribution and use in source and binary
//  forms, with or without modification, are permitted provided that the following conditions
//  are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this list of
//     conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright notice, this list
//     of conditions and the following disclaimer in the documentation and/or other materials
//     provided with the distribution.
//  3. Neither the names of the Blaze development group nor the names of its contributors
//     may be used to endorse or promote products derived from this software without specific
//     prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
//  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
//  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
//  SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
//  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
//  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
//  BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
//  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
//  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
//  DAMAGE.
*/
//=================================================================================================

#ifndef _BLAZE_MATH_ADAPTORS_LOWERMATRIX_STRICTLYLOWERPROXY_H_
#define _BLAZE_MATH_ADAPTORS_LOWERMATRIX_STRICTLYLOWERPROXY_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <ostream>
#include <stdexcept>
#include <blaze/math/constraints/Expression.h>
#include <blaze/math/constraints/Lower.h>
#include <blaze/math/constraints/Matrix.h>
#include <blaze/math/constraints/Symmetric.h>
#include <blaze/math/constraints/Upper.h>
#include <blaze/math/proxy/Proxy.h>
#include <blaze/math/shims/Clear.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/shims/Reset.h>
#include <blaze/util/constraints/Const.h>
#include <blaze/util/constraints/Pointer.h>
#include <blaze/util/constraints/Reference.h>
#include <blaze/util/constraints/Volatile.h>
#include <blaze/util/typetraits/AddConst.h>
#include <blaze/util/typetraits/AddReference.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Access proxy for strictly lower triangular matrices.
// \ingroup strictly_lower_matrix
//
// The StrictlyLowerProxy provides controlled access to the elements of a non-const strictly lower
// triangular matrix. It guarantees that the strictly lower matrix invariant is not violated, i.e.
// that elements on the diagonal and in the upper part of the matrix remain 0. The following
// example illustrates this by means of a \f$ 3 \times 3 \f$ dense strictly lower triangular
// matrix:

   \code
   // Creating a 3x3 strictly lower triangular dense matrix
   blaze::StrictlyLowerMatrix< blaze::DynamicMatrix<int> > A( 3UL );

   A(0,1) = -2;  //        (  0 0 0 )
   A(2,1) =  3;  // => A = ( -2 0 0 )
   A(2,2) =  5;  //        (  3 5 0 )

   A(1,1) =  4;  // Invalid assignment to diagonal matrix element; results in an exception!
   A(0,2) =  7;  // Invalid assignment to upper matrix element; results in an exception!
   \endcode
*/
template< typename MT >  // Type of the adapted matrix
class StrictlyLowerProxy : public Proxy< StrictlyLowerProxy<MT>, typename MT::ElementType >
{
 private:
   //**Type definitions****************************************************************************
   //! Reference type of the underlying matrix type.
   typedef typename AddConst< typename MT::Reference >::Type  ReferenceType;
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   //! Type of the represented matrix element.
   typedef typename MT::ElementType  RepresentedType;

   //! Reference to the represented element.
   typedef typename AddReference<ReferenceType>::Type  RawReference;
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline StrictlyLowerProxy( MT& matrix, size_t row, size_t column );
            inline StrictlyLowerProxy( const StrictlyLowerProxy& ulp );
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   // No explicitly declared destructor.
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   /*!\name Assignment operators */
   //@{
                          inline const StrictlyLowerProxy& operator= ( const StrictlyLowerProxy& ulp ) const;
   template< typename T > inline const StrictlyLowerProxy& operator= ( const T& value ) const;
   template< typename T > inline const StrictlyLowerProxy& operator+=( const T& value ) const;
   template< typename T > inline const StrictlyLowerProxy& operator-=( const T& value ) const;
   template< typename T > inline const StrictlyLowerProxy& operator*=( const T& value ) const;
   template< typename T > inline const StrictlyLowerProxy& operator/=( const T& value ) const;
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline RawReference get()          const;
   inline bool         isRestricted() const;
   //@}
   //**********************************************************************************************

   //**Conversion operator*************************************************************************
   /*!\name Conversion operator */
   //@{
   inline operator RawReference() const;
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   ReferenceType value_;    //!< Reference to the accessed matrix element.
   const bool restricted_;  //!< Access flag for the accessed matrix element.
                            /*!< The flag indicates if access to the matrix element is
                                 restricted. It is \a true in case the proxy represents an
                                 element on the diagonal or in the upper part of the matrix. */
   //@}
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_MATRIX_TYPE              ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_REFERENCE_TYPE       ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_POINTER_TYPE         ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_CONST                ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VOLATILE             ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_EXPRESSION_TYPE      ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_LOWER_MATRIX_TYPE    ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UPPER_MATRIX_TYPE    ( MT );
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
/*!\brief Initialization constructor for an StrictlyLowerProxy.
//
// \param matrix Reference to the adapted matrix.
// \param row The row-index of the accessed matrix element.
// \param column The column-index of the accessed matrix element.
*/
template< typename MT >  // Type of the adapted matrix
inline StrictlyLowerProxy<MT>::StrictlyLowerProxy( MT& matrix, size_t row, size_t column )
   : value_     ( matrix( row, column ) )  // Reference to the accessed matrix element
   , restricted_( row <= column )          // Access flag for the accessed matrix element
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The copy constructor for LowerProxy.
//
// \param slp Proxy to be copied.
*/
template< typename MT >  // Type of the adapted matrix
inline StrictlyLowerProxy<MT>::StrictlyLowerProxy( const StrictlyLowerProxy& slp )
   : value_     ( slp.value_      )  // Reference to the accessed matrix element
   , restricted_( slp.restricted_ )  // Access flag for the accessed matrix element
{}
//*************************************************************************************************




//=================================================================================================
//
//  OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Copy assignment operator for StrictlyLowerProxy.
//
// \param slp Proxy to be copied.
// \return Reference to the assigned proxy.
*/
template< typename MT >  // Type of the adapted matrix
inline const StrictlyLowerProxy<MT>& StrictlyLowerProxy<MT>::operator=( const StrictlyLowerProxy& slp ) const
{
   value_ = slp.value_;

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Assignment to the accessed matrix element.
//
// \param value The new value of the matrix element.
// \return Reference to the assigned proxy.
// \exception std::invalid_argument Invalid assignment to diagonal or upper matrix element.
//
// In case the proxy represents an element on the diagonal or in the upper part of the matrix,
// a \a std::invalid_argument exception is thrown.
*/
template< typename MT >  // Type of the adapted matrix
template< typename T >   // Type of the right-hand side value
inline const StrictlyLowerProxy<MT>& StrictlyLowerProxy<MT>::operator=( const T& value ) const
{
   if( restricted_ )
      throw std::invalid_argument( "Invalid assignment to diagonal or upper matrix element" );

   value_ = value;

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Addition assignment to the accessed matrix element.
//
// \param value The right-hand side value to be added to the matrix element.
// \return Reference to the assigned proxy.
// \exception std::invalid_argument Invalid assignment to diagonal or upper matrix element.
//
// In case the proxy represents an element on the diagonal or in the upper part of the matrix,
// a \a std::invalid_argument exception is thrown.
*/
template< typename MT >  // Type of the adapted matrix
template< typename T >   // Type of the right-hand side value
inline const StrictlyLowerProxy<MT>& StrictlyLowerProxy<MT>::operator+=( const T& value ) const
{
   if( restricted_ )
      throw std::invalid_argument( "Invalid assignment to diagonal or upper matrix element" );

   value_ += value;

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Subtraction assignment to the accessed matrix element.
//
// \param value The right-hand side value to be subtracted from the matrix element.
// \return Reference to the assigned proxy.
// \exception std::invalid_argument Invalid assignment to diagonal or upper matrix element.
//
// In case the proxy represents an element on the diagonal or in the upper part of the matrix,
// a \a std::invalid_argument exception is thrown.
*/
template< typename MT >  // Type of the adapted matrix
template< typename T >   // Type of the right-hand side value
inline const StrictlyLowerProxy<MT>& StrictlyLowerProxy<MT>::operator-=( const T& value ) const
{
   if( restricted_ )
      throw std::invalid_argument( "Invalid assignment to diagonal or upper matrix element" );

   value_ -= value;

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication assignment to the accessed matrix element.
//
// \param value The right-hand side value for the multiplication.
// \return Reference to the assigned proxy.
// \exception std::invalid_argument Invalid assignment to diagonal or upper matrix element.
//
// In case the proxy represents an element on the diagonal or in the upper part of the matrix,
// a \a std::invalid_argument exception is thrown.
*/
template< typename MT >  // Type of the adapted matrix
template< typename T >   // Type of the right-hand side value
inline const StrictlyLowerProxy<MT>& StrictlyLowerProxy<MT>::operator*=( const T& value ) const
{
   if( restricted_ )
      throw std::invalid_argument( "Invalid assignment to diagonal or upper matrix element" );

   value_ *= value;

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Division assignment to the accessed matrix element.
//
// \param value The right-hand side value for the division.
// \return Reference to the assigned proxy.
// \exception std::invalid_argument Invalid assignment to diagonal or upper matrix element.
//
// In case the proxy represents an element on the diagonal or in the upper part of the matrix,
// a \a std::invalid_argument exception is thrown.
*/
template< typename MT >  // Type of the adapted matrix
template< typename T >   // Type of the right-hand side value
inline const StrictlyLowerProxy<MT>& StrictlyLowerProxy<MT>::operator/=( const T& value ) const
{
   if( restricted_ )
      throw std::invalid_argument( "Invalid assignment to diagonal or upper matrix element" );

   value_ /= value;

   return *this;
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returning the value of the accessed matrix element.
//
// \return Direct/raw reference to the accessed matrix element.
*/
template< typename MT >  // Type of the adapted matrix
inline typename StrictlyLowerProxy<MT>::RawReference StrictlyLowerProxy<MT>::get() const
{
   return value_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the proxy represents a restricted matrix element..
//
// \return \a true in case access to the matrix element is restricted, \a false if not.
*/
template< typename MT >  // Type of the adapted matrix
inline bool StrictlyLowerProxy<MT>::isRestricted() const
{
   return restricted_;
}
//*************************************************************************************************




//=================================================================================================
//
//  CONVERSION OPERATOR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Conversion to the accessed matrix element.
//
// \return Direct/raw reference to the accessed matrix element.
*/
template< typename MT >  // Type of the adapted matrix
inline StrictlyLowerProxy<MT>::operator RawReference() const
{
   return get();
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name StrictlyLowerProxy operators */
//@{
template< typename MT1, typename MT2 >
inline bool operator==( const StrictlyLowerProxy<MT1>& lhs, const StrictlyLowerProxy<MT2>& rhs );

template< typename MT, typename T >
inline bool operator==( const StrictlyLowerProxy<MT>& lhs, const T& rhs );

template< typename T, typename MT >
inline bool operator==( const T& lhs, const StrictlyLowerProxy<MT>& rhs );

template< typename MT1, typename MT2 >
inline bool operator!=( const StrictlyLowerProxy<MT1>& lhs, const StrictlyLowerProxy<MT2>& rhs );

template< typename MT, typename T >
inline bool operator!=( const StrictlyLowerProxy<MT>& lhs, const T& rhs );

template< typename T, typename MT >
inline bool operator!=( const T& lhs, const StrictlyLowerProxy<MT>& rhs );

template< typename MT1, typename MT2 >
inline bool operator<( const StrictlyLowerProxy<MT1>& lhs, const StrictlyLowerProxy<MT2>& rhs );

template< typename MT, typename T >
inline bool operator<( const StrictlyLowerProxy<MT>& lhs, const T& rhs );

template< typename T, typename MT >
inline bool operator<( const T& lhs, const StrictlyLowerProxy<MT>& rhs );

template< typename MT1, typename MT2 >
inline bool operator>( const StrictlyLowerProxy<MT1>& lhs, const StrictlyLowerProxy<MT2>& rhs );

template< typename MT, typename T >
inline bool operator>( const StrictlyLowerProxy<MT>& lhs, const T& rhs );

template< typename T, typename MT >
inline bool operator>( const T& lhs, const StrictlyLowerProxy<MT>& rhs );

template< typename MT1, typename MT2 >
inline bool operator<=( const StrictlyLowerProxy<MT1>& lhs, const StrictlyLowerProxy<MT2>& rhs );

template< typename MT, typename T >
inline bool operator<=( const StrictlyLowerProxy<MT>& lhs, const T& rhs );

template< typename T, typename MT >
inline bool operator<=( const T& lhs, const StrictlyLowerProxy<MT>& rhs );

template< typename MT1, typename MT2 >
inline bool operator>=( const StrictlyLowerProxy<MT1>& lhs, const StrictlyLowerProxy<MT2>& rhs );

template< typename MT, typename T >
inline bool operator>=( const StrictlyLowerProxy<MT>& lhs, const T& rhs );

template< typename T, typename MT >
inline bool operator>=( const T& lhs, const StrictlyLowerProxy<MT>& rhs );

template< typename MT >
inline std::ostream& operator<<( std::ostream& os, const StrictlyLowerProxy<MT>& proxy );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality comparison between two StrictlyLowerProxy objects.
// \ingroup strictly_lower_matrix
//
// \param lhs The left-hand side StrictlyLowerProxy object.
// \param rhs The right-hand side StrictlyLowerProxy object.
// \return \a true if both referenced values are equal, \a false if they are not.
*/
template< typename MT1, typename MT2 >
inline bool operator==( const StrictlyLowerProxy<MT1>& lhs, const StrictlyLowerProxy<MT2>& rhs )
{
   return ( lhs.get() == rhs.get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality comparison between a StrictlyLowerProxy object and an object of different type.
// \ingroup strictly_lower_matrix
//
// \param lhs The left-hand side StrictlyLowerProxy object.
// \param rhs The right-hand side object of other type.
// \return \a true if the referenced value and the other object are equal, \a false if they are not.
*/
template< typename MT, typename T >
inline bool operator==( const StrictlyLowerProxy<MT>& lhs, const T& rhs )
{
   return ( lhs.get() == rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality comparison between an object of different type and a StrictlyLowerProxy object.
// \ingroup strictly_lower_matrix
//
// \param lhs The left-hand side object of other type.
// \param rhs The right-hand side StrictlyLowerProxy object.
// \return \a true if the other object and the referenced value are equal, \a false if they are not.
*/
template< typename T, typename MT >
inline bool operator==( const T& lhs, const StrictlyLowerProxy<MT>& rhs )
{
   return ( lhs == rhs.get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inequality comparison between two StrictlyLowerProxy objects.
// \ingroup strictly_lower_matrix
//
// \param lhs The left-hand side StrictlyLowerProxy object.
// \param rhs The right-hand side StrictlyLowerProxy object.
// \return \a true if both referenced values are not equal, \a false if they are.
*/
template< typename MT1, typename MT2 >
inline bool operator!=( const StrictlyLowerProxy<MT1>& lhs, const StrictlyLowerProxy<MT2>& rhs )
{
   return ( lhs.get() != rhs.get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inequality comparison between a StrictlyLowerProxy object and an object of different type.
// \ingroup strictly_lower_matrix
//
// \param lhs The left-hand side StrictlyLowerProxy object.
// \param rhs The right-hand side object of other type.
// \return \a true if the referenced value and the other object are not equal, \a false if they are.
*/
template< typename MT, typename T >
inline bool operator!=( const StrictlyLowerProxy<MT>& lhs, const T& rhs )
{
   return ( lhs.get() != rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inquality comparison between an object of different type and a StrictlyLowerProxy object.
// \ingroup strictly_lower_matrix
//
// \param lhs The left-hand side object of other type.
// \param rhs The right-hand side StrictlyLowerProxy object.
// \return \a true if the other object and the referenced value are not equal, \a false if they are.
*/
template< typename T, typename MT >
inline bool operator!=( const T& lhs, const StrictlyLowerProxy<MT>& rhs )
{
   return ( lhs != rhs.get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Less-than comparison between two StrictlyLowerProxy objects.
// \ingroup strictly_lower_matrix
//
// \param lhs The left-hand side StrictlyLowerProxy object.
// \param rhs The right-hand side StrictlyLowerProxy object.
// \return \a true if the left-hand side referenced value is smaller, \a false if not.
*/
template< typename MT1, typename MT2 >
inline bool operator<( const StrictlyLowerProxy<MT1>& lhs, const StrictlyLowerProxy<MT2>& rhs )
{
   return ( lhs.get() < rhs.get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Less-than comparison between a StrictlyLowerProxy object and an object of different type.
// \ingroup strictly_lower_matrix
//
// \param lhs The left-hand side StrictlyLowerProxy object.
// \param rhs The right-hand side object of other type.
// \return \a true if the left-hand side referenced value is smaller, \a false if not.
*/
template< typename MT, typename T >
inline bool operator<( const StrictlyLowerProxy<MT>& lhs, const T& rhs )
{
   return ( lhs.get() < rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Less-than comparison between an object of different type and a StrictlyLowerProxy object.
// \ingroup strictly_lower_matrix
//
// \param lhs The left-hand side object of other type.
// \param rhs The right-hand side StrictlyLowerProxy object.
// \return \a true if the left-hand side other object is smaller, \a false if not.
*/
template< typename T, typename MT >
inline bool operator<( const T& lhs, const StrictlyLowerProxy<MT>& rhs )
{
   return ( lhs < rhs.get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Greater-than comparison between two StrictlyLowerProxy objects.
// \ingroup strictly_lower_matrix
//
// \param lhs The left-hand side StrictlyLowerProxy object.
// \param rhs The right-hand side StrictlyLowerProxy object.
// \return \a true if the left-hand side referenced value is greater, \a false if not.
*/
template< typename MT1, typename MT2 >
inline bool operator>( const StrictlyLowerProxy<MT1>& lhs, const StrictlyLowerProxy<MT2>& rhs )
{
   return ( lhs.get() > rhs.get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Greater-than comparison between a StrictlyLowerProxy object and an object of different type.
// \ingroup strictly_lower_matrix
//
// \param lhs The left-hand side StrictlyLowerProxy object.
// \param rhs The right-hand side object of other type.
// \return \a true if the left-hand side referenced value is greater, \a false if not.
*/
template< typename MT, typename T >
inline bool operator>( const StrictlyLowerProxy<MT>& lhs, const T& rhs )
{
   return ( lhs.get() > rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Greater-than comparison between an object of different type and a StrictlyLowerProxy object.
// \ingroup strictly_lower_matrix
//
// \param lhs The left-hand side object of other type.
// \param rhs The right-hand side StrictlyLowerProxy object.
// \return \a true if the left-hand side other object is greater, \a false if not.
*/
template< typename T, typename MT >
inline bool operator>( const T& lhs, const StrictlyLowerProxy<MT>& rhs )
{
   return ( lhs > rhs.get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Less-or-equal-than comparison between two StrictlyLowerProxy objects.
// \ingroup strictly_lower_matrix
//
// \param lhs The left-hand side StrictlyLowerProxy object.
// \param rhs The right-hand side StrictlyLowerProxy object.
// \return \a true if the left-hand side referenced value is smaller or equal, \a false if not.
*/
template< typename MT1, typename MT2 >
inline bool operator<=( const StrictlyLowerProxy<MT1>& lhs, const StrictlyLowerProxy<MT2>& rhs )
{
   return ( lhs.get() <= rhs.get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Less-or-equal-than comparison between a StrictlyLowerProxy object and an object of different type.
// \ingroup strictly_lower_matrix
//
// \param lhs The left-hand side StrictlyLowerProxy object.
// \param rhs The right-hand side object of other type.
// \return \a true if the left-hand side referenced value is smaller or equal, \a false if not.
*/
template< typename MT, typename T >
inline bool operator<=( const StrictlyLowerProxy<MT>& lhs, const T& rhs )
{
   return ( lhs.get() <= rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Less-or-equal-than comparison between an object of different type and a StrictlyLowerProxy object.
// \ingroup strictly_lower_matrix
//
// \param lhs The left-hand side object of other type.
// \param rhs The right-hand side StrictlyLowerProxy object.
// \return \a true if the left-hand side other object is smaller or equal, \a false if not.
*/
template< typename T, typename MT >
inline bool operator<=( const T& lhs, const StrictlyLowerProxy<MT>& rhs )
{
   return ( lhs <= rhs.get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Greater-or-equal-than comparison between two StrictlyLowerProxy objects.
// \ingroup strictly_lower_matrix
//
// \param lhs The left-hand side StrictlyLowerProxy object.
// \param rhs The right-hand side StrictlyLowerProxy object.
// \return \a true if the left-hand side referenced value is greater or equal, \a false if not.
*/
template< typename MT1, typename MT2 >
inline bool operator>=( const StrictlyLowerProxy<MT1>& lhs, const StrictlyLowerProxy<MT2>& rhs )
{
   return ( lhs.get() >= rhs.get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Greater-or-equal-than comparison between a StrictlyLowerProxy object and an object of
//        different type.
// \ingroup strictly_lower_matrix
//
// \param lhs The left-hand side StrictlyLowerProxy object.
// \param rhs The right-hand side object of other type.
// \return \a true if the left-hand side referenced value is greater or equal, \a false if not.
*/
template< typename MT, typename T >
inline bool operator>=( const StrictlyLowerProxy<MT>& lhs, const T& rhs )
{
   return ( lhs.get() >= rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Greater-or-equal-than comparison between an object of different type and a
//        StrictlyLowerProxy object.
// \ingroup strictly_lower_matrix
//
// \param lhs The left-hand side object of other type.
// \param rhs The right-hand side StrictlyLowerProxy object.
// \return \a true if the left-hand side other object is greater or equal, \a false if not.
*/
template< typename T, typename MT >
inline bool operator>=( const T& lhs, const StrictlyLowerProxy<MT>& rhs )
{
   return ( lhs >= rhs.get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Global output operator for proxies on strictly lower triangular matrices.
// \ingroup strictly_lower_matrix
//
// \param os Reference to the output stream.
// \param proxy Reference to a constant proxy object.
// \return Reference to the output stream.
*/
template< typename MT >
inline std::ostream& operator<<( std::ostream& os, const StrictlyLowerProxy<MT>& proxy )
{
   return os << proxy.get();
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\name StrictlyLowerProxy global functions */
//@{
template< typename MT >
inline void reset( const StrictlyLowerProxy<MT>& proxy );

template< typename MT >
inline void clear( const StrictlyLowerProxy<MT>& proxy );

template< typename MT >
inline bool isDefault( const StrictlyLowerProxy<MT>& proxy );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Resetting the represented element to the default initial values.
// \ingroup strictly_lower_matrix
//
// \param proxy The given access proxy.
// \return void
//
// This function resets the element represented by the access proxy to its default initial
// value.
*/
template< typename MT >
inline void reset( const StrictlyLowerProxy<MT>& proxy )
{
   using blaze::reset;

   reset( proxy.get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Clearing the represented element.
// \ingroup strictly_lower_matrix
//
// \param proxy The given access proxy.
// \return void
//
// This function clears the element represented by the access proxy to its default initial
// state.
*/
template< typename MT >
inline void clear( const StrictlyLowerProxy<MT>& proxy )
{
   using blaze::clear;

   clear( proxy.get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the represented element is in default state.
// \ingroup strictly_lower_matrix
//
// \param proxy The given access proxy
// \return \a true in case the represented element is in default state, \a false otherwise.
//
// This function checks whether the element represented by the access proxy is in default state.
// In case it is in default state, the function returns \a true, otherwise it returns \a false.
*/
template< typename MT >
inline bool isDefault( const StrictlyLowerProxy<MT>& proxy )
{
   using blaze::isDefault;

   return isDefault( proxy.get() );
}
//*************************************************************************************************

} // namespace blaze

#endif
