//=================================================================================================
/*!
//  \file blaze/math/adaptors/symmetricmatrix/NumericProxy.h
//  \brief Header file for the NumericProxy class
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

#ifndef _BLAZE_MATH_ADAPTORS_SYMMETRICMATRIX_NUMERICPROXY_H_
#define _BLAZE_MATH_ADAPTORS_SYMMETRICMATRIX_NUMERICPROXY_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <ostream>
#include <blaze/math/constraints/Expression.h>
#include <blaze/math/constraints/Lower.h>
#include <blaze/math/constraints/Matrix.h>
#include <blaze/math/constraints/Symmetric.h>
#include <blaze/math/constraints/Upper.h>
#include <blaze/math/shims/Clear.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/shims/Reset.h>
#include <blaze/util/constraints/Const.h>
#include <blaze/util/constraints/Numeric.h>
#include <blaze/util/constraints/Pointer.h>
#include <blaze/util/constraints/Reference.h>
#include <blaze/util/constraints/Volatile.h>
#include <blaze/util/InvalidType.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/IsComplex.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Access proxy for symmetric, square matrices with numeric element types.
// \ingroup symmetric_matrix
//
// The NumericProxy provides controlled access to the elements of a non-const symmetric matrix
// with numeric element type (e.g. integral values, floating point values, and complex values).
// It guarantees that a modification of element \f$ a_{ij} \f$ of the accessed matrix is also
// applied to element \f$ a_{ji} \f$. The following example illustrates this by means of a
// \f$ 3 \times 3 \f$ dense symmetric matrix:

   \code
   // Creating a 3x3 symmetric dense matrix
   blaze::SymmetricMatrix< blaze::DynamicMatrix<int> > A( 3UL );

   A(0,2) = -2;  //        (  0 0 -2 )
   A(1,1) =  3;  // => A = (  0 3  5 )
   A(1,2) =  5;  //        ( -2 5  0 )
   \endcode
*/
template< typename MT >  // Type of the adapted matrix
class NumericProxy
{
 private:
   //**struct BuiltinType**************************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Auxiliary struct to determine the value type of the represented complex element.
   */
   template< typename T >
   struct BuiltinType { typedef INVALID_TYPE  Type; };
   /*! \endcond */
   //**********************************************************************************************

   //**struct ComplexType**************************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Auxiliary struct to determine the value type of the represented complex element.
   */
   template< typename T >
   struct ComplexType { typedef typename T::value_type  Type; };
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   typedef typename MT::ElementType     RepresentedType;  //!< Type of the represented matrix element.
   typedef typename MT::Reference       Reference;        //!< Reference to the represented element.
   typedef typename MT::ConstReference  ConstReference;   //!< Reference-to-const to the represented element.
   typedef NumericProxy*                Pointer;          //!< Pointer to the represented element.
   typedef const NumericProxy*          ConstPointer;     //!< Pointer-to-const to the represented element.

   //! Value type of the represented complex element.
   typedef typename If< IsComplex<RepresentedType>
                      , ComplexType<RepresentedType>
                      , BuiltinType<RepresentedType> >::Type::Type  ValueType;
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline NumericProxy( MT& matrix, size_t row, size_t column );
            inline NumericProxy( const NumericProxy& np );
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   // No explicitly declared destructor.
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   /*!\name Assignment operators */
   //@{
                          inline NumericProxy& operator= ( const NumericProxy& sp );
   template< typename T > inline NumericProxy& operator= ( const T& value );
   template< typename T > inline NumericProxy& operator+=( const T& value );
   template< typename T > inline NumericProxy& operator-=( const T& value );
   template< typename T > inline NumericProxy& operator*=( const T& value );
   template< typename T > inline NumericProxy& operator/=( const T& value );
   //@}
   //**********************************************************************************************

   //**Access operators****************************************************************************
   /*!\name Access operators */
   //@{
   inline Pointer      operator->();
   inline ConstPointer operator->() const;
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline void reset() const;
   inline void clear() const;

   inline ConstReference get() const;
   //@}
   //**********************************************************************************************

   //**Conversion operator*************************************************************************
   /*!\name Conversion operator */
   //@{
   inline operator ConstReference() const;
   //@}
   //**********************************************************************************************

   //**Complex data access functions***************************************************************
   /*!\name Complex data access functions */
   //@{
   inline ValueType real() const;
   inline void      real( ValueType value ) const;
   inline ValueType imag() const;
   inline void      imag( ValueType value ) const;
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   MT&    matrix_;  //!< Reference to the adapted matrix.
   size_t row_;     //!< Row index of the accessed matrix element.
   size_t column_;  //!< Column index of the accessed matrix element.
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
   BLAZE_CONSTRAINT_MUST_BE_NUMERIC_TYPE             ( RepresentedType );
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
/*!\brief Initialization constructor for a NumericProxy.
//
// \param matrix Reference to the adapted matrix.
// \param row The row-index of the accessed matrix element.
// \param column The column-index of the accessed matrix element.
*/
template< typename MT >  // Type of the adapted matrix
inline NumericProxy<MT>::NumericProxy( MT& matrix, size_t row, size_t column )
   : matrix_( matrix )  // Reference to the adapted matrix
   , row_   ( row    )  // Row index of the accessed matrix element
   , column_( column )  // Column index of the accessed matrix element
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The copy constructor for NumericProxy.
//
// \param np Numeric proxy to be copied.
*/
template< typename MT >  // Type of the adapted matrix
inline NumericProxy<MT>::NumericProxy( const NumericProxy& np )
   : matrix_( np.matrix_ )  // Reference to the adapted matrix
   , row_   ( np.row_    )  // Row index of the accessed matrix element
   , column_( np.column_ )  // Column index of the accessed matrix element
{}
//*************************************************************************************************




//=================================================================================================
//
//  OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Copy assignment operator for NumericProxy.
//
// \param np Numeric proxy to be copied.
// \return Reference to the assigned proxy.
*/
template< typename MT >  // Type of the adapted matrix
inline NumericProxy<MT>& NumericProxy<MT>::operator=( const NumericProxy& np )
{
   matrix_(row_,column_) = np.matrix_(np.row_,np.column_);
   matrix_(column_,row_) = np.matrix_(np.row_,np.column_);

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Assignment to the accessed matrix element.
//
// \param value The new value of the matrix element.
// \return Reference to the assigned proxy.
*/
template< typename MT >  // Type of the adapted matrix
template< typename T >   // Type of the right-hand side value
inline NumericProxy<MT>& NumericProxy<MT>::operator=( const T& value )
{
   matrix_(row_,column_) = value;
   if( row_ != column_ )
      matrix_(column_,row_) = value;

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Addition assignment to the accessed matrix element.
//
// \param value The right-hand side value to be added to the matrix element.
// \return Reference to the assigned proxy.
*/
template< typename MT >  // Type of the adapted matrix
template< typename T >   // Type of the right-hand side value
inline NumericProxy<MT>& NumericProxy<MT>::operator+=( const T& value )
{
   matrix_(row_,column_) += value;
   if( row_ != column_ )
      matrix_(column_,row_) += value;

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Subtraction assignment to the accessed matrix element.
//
// \param value The right-hand side value to be subtracted from the matrix element.
// \return Reference to the assigned proxy.
*/
template< typename MT >  // Type of the adapted matrix
template< typename T >   // Type of the right-hand side value
inline NumericProxy<MT>& NumericProxy<MT>::operator-=( const T& value )
{
   matrix_(row_,column_) -= value;
   if( row_ != column_ )
      matrix_(column_,row_) -= value;

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication assignment to the accessed matrix element.
//
// \param value The right-hand side value for the multiplication.
// \return Reference to the assigned proxy.
*/
template< typename MT >  // Type of the adapted matrix
template< typename T >   // Type of the right-hand side value
inline NumericProxy<MT>& NumericProxy<MT>::operator*=( const T& value )
{
   matrix_(row_,column_) *= value;
   if( row_ != column_ )
      matrix_(column_,row_) *= value;

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Division assignment to the accessed matrix element.
//
// \param value The right-hand side value for the division.
// \return Reference to the assigned proxy.
*/
template< typename MT >  // Type of the adapted matrix
template< typename T >   // Type of the right-hand side value
inline NumericProxy<MT>& NumericProxy<MT>::operator/=( const T& value )
{
   matrix_(row_,column_) /= value;
   if( row_ != column_ )
      matrix_(column_,row_) /= value;

   return *this;
}
//*************************************************************************************************




//=================================================================================================
//
//  ACCESS OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Direct access to the represented matrix element.
//
// \return Pointer to the represented matrix element.
*/
template< typename MT >  // Type of the adapted matrix
inline typename NumericProxy<MT>::Pointer NumericProxy<MT>::operator->()
{
   return this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Direct access to the represented matrix element.
//
// \return Pointer to the represented matrix element.
*/
template< typename MT >  // Type of the adapted matrix
inline typename NumericProxy<MT>::ConstPointer NumericProxy<MT>::operator->() const
{
   return this;
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Reset the represented element to its default initial value.
//
// \return void
//
// This function resets the element represented by the proxy to its default initial value.
*/
template< typename MT >  // Type of the adapted matrix
inline void NumericProxy<MT>::reset() const
{
   using blaze::reset;

   reset( matrix_(row_,column_) );
   if( row_ != column_ )
      reset( matrix_(column_,row_) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Clearing the represented element.
//
// \return void
//
// This function clears the element represented by the proxy to its default initial state.
*/
template< typename MT >  // Type of the adapted matrix
inline void NumericProxy<MT>::clear() const
{
   using blaze::clear;

   clear( matrix_(row_,column_) );
   if( row_ != column_ )
      clear( matrix_(column_,row_) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returning the value of the accessed matrix element.
//
// \return Direct/raw reference to the accessed matrix element.
*/
template< typename MT >  // Type of the adapted matrix
inline typename NumericProxy<MT>::ConstReference NumericProxy<MT>::get() const
{
   return const_cast<const MT&>( matrix_ )(row_,column_);
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
inline NumericProxy<MT>::operator ConstReference() const
{
   return get();
}
//*************************************************************************************************




//=================================================================================================
//
//  COMPLEX DATA ACCESS FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns the real part of the represented complex number.
//
// \return The current real part of the represented complex number.
//
// In case the proxy represents a complex number, this function returns the current value of its
// real part.
*/
template< typename MT >  // Type of the adapted matrix
inline typename NumericProxy<MT>::ValueType NumericProxy<MT>::real() const
{
   return matrix_(row_,column_).real();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Setting the real part of the represented complex number.
//
// \param value The new value for the real part.
// \return void
//
// In case the proxy represents a complex number, this function sets a new value to its real part.
*/
template< typename MT >  // Type of the adapted matrix
inline void NumericProxy<MT>::real( ValueType value ) const
{
   matrix_(row_,column_).real( value );
   if( row_ != column_ )
      matrix_(column_,row_).real( value );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the imaginary part of the represented complex number.
//
// \return The current imaginary part of the represented complex number.
//
// In case the proxy represents a complex number, this function returns the current value of its
// imaginary part.
*/
template< typename MT >  // Type of the adapted matrix
inline typename NumericProxy<MT>::ValueType NumericProxy<MT>::imag() const
{
   return matrix_(row_,column_).imag();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Setting the imaginary part of the represented complex number.
//
// \param value The new value for the imaginary part.
// \return void
//
// In case the proxy represents a complex number, this function sets a new value to its imaginary
// part.
*/
template< typename MT >  // Type of the adapted matrix
inline void NumericProxy<MT>::imag( ValueType value ) const
{
   matrix_(row_,column_).imag( value );
   if( row_ != column_ )
      matrix_(column_,row_).imag( value );
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name NumericProxy operators */
//@{
template< typename MT1, typename MT2 >
inline bool operator==( const NumericProxy<MT1>& lhs, const NumericProxy<MT2>& rhs );

template< typename MT, typename T >
inline bool operator==( const NumericProxy<MT>& lhs, const T& rhs );

template< typename T, typename MT >
inline bool operator==( const T& lhs, const NumericProxy<MT>& rhs );

template< typename MT1, typename MT2 >
inline bool operator!=( const NumericProxy<MT1>& lhs, const NumericProxy<MT2>& rhs );

template< typename MT, typename T >
inline bool operator!=( const NumericProxy<MT>& lhs, const T& rhs );

template< typename T, typename MT >
inline bool operator!=( const T& lhs, const NumericProxy<MT>& rhs );

template< typename MT1, typename MT2 >
inline bool operator<( const NumericProxy<MT1>& lhs, const NumericProxy<MT2>& rhs );

template< typename MT, typename T >
inline bool operator<( const NumericProxy<MT>& lhs, const T& rhs );

template< typename T, typename MT >
inline bool operator<( const T& lhs, const NumericProxy<MT>& rhs );

template< typename MT1, typename MT2 >
inline bool operator>( const NumericProxy<MT1>& lhs, const NumericProxy<MT2>& rhs );

template< typename MT, typename T >
inline bool operator>( const NumericProxy<MT>& lhs, const T& rhs );

template< typename T, typename MT >
inline bool operator>( const T& lhs, const NumericProxy<MT>& rhs );

template< typename MT1, typename MT2 >
inline bool operator<=( const NumericProxy<MT1>& lhs, const NumericProxy<MT2>& rhs );

template< typename MT, typename T >
inline bool operator<=( const NumericProxy<MT>& lhs, const T& rhs );

template< typename T, typename MT >
inline bool operator<=( const T& lhs, const NumericProxy<MT>& rhs );

template< typename MT1, typename MT2 >
inline bool operator>=( const NumericProxy<MT1>& lhs, const NumericProxy<MT2>& rhs );

template< typename MT, typename T >
inline bool operator>=( const NumericProxy<MT>& lhs, const T& rhs );

template< typename T, typename MT >
inline bool operator>=( const T& lhs, const NumericProxy<MT>& rhs );

template< typename MT >
inline std::ostream& operator<<( std::ostream& os, const NumericProxy<MT>& proxy );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality comparison between two NumericProxy objects.
// \ingroup symmetric_matrix
//
// \param lhs The left-hand side NumericProxy object.
// \param rhs The right-hand side NumericProxy object.
// \return \a true if both referenced values are equal, \a false if they are not.
*/
template< typename MT1, typename MT2 >
inline bool operator==( const NumericProxy<MT1>& lhs, const NumericProxy<MT2>& rhs )
{
   return ( lhs.get() == rhs.get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality comparison between a NumericProxy object and an object of different type.
// \ingroup symmetric_matrix
//
// \param lhs The left-hand side NumericProxy object.
// \param rhs The right-hand side object of other type.
// \return \a true if the referenced value and the other object are equal, \a false if they are not.
*/
template< typename MT, typename T >
inline bool operator==( const NumericProxy<MT>& lhs, const T& rhs )
{
   return ( lhs.get() == rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality comparison between an object of different type and a NumericProxy object.
// \ingroup symmetric_matrix
//
// \param lhs The left-hand side object of other type.
// \param rhs The right-hand side NumericProxy object.
// \return \a true if the other object and the referenced value are equal, \a false if they are not.
*/
template< typename T, typename MT >
inline bool operator==( const T& lhs, const NumericProxy<MT>& rhs )
{
   return ( lhs == rhs.get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inequality comparison between two NumericProxy objects.
// \ingroup symmetric_matrix
//
// \param lhs The left-hand side NumericProxy object.
// \param rhs The right-hand side NumericProxy object.
// \return \a true if both referenced values are not equal, \a false if they are.
*/
template< typename MT1, typename MT2 >
inline bool operator!=( const NumericProxy<MT1>& lhs, const NumericProxy<MT2>& rhs )
{
   return ( lhs.get() != rhs.get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inequality comparison between a NumericProxy object and an object of different type.
// \ingroup symmetric_matrix
//
// \param lhs The left-hand side NumericProxy object.
// \param rhs The right-hand side object of other type.
// \return \a true if the referenced value and the other object are not equal, \a false if they are.
*/
template< typename MT, typename T >
inline bool operator!=( const NumericProxy<MT>& lhs, const T& rhs )
{
   return ( lhs.get() != rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inquality comparison between an object of different type and a NumericProxy object.
// \ingroup symmetric_matrix
//
// \param lhs The left-hand side object of other type.
// \param rhs The right-hand side NumericProxy object.
// \return \a true if the other object and the referenced value are not equal, \a false if they are.
*/
template< typename T, typename MT >
inline bool operator!=( const T& lhs, const NumericProxy<MT>& rhs )
{
   return ( lhs != rhs.get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Less-than comparison between two NumericProxy objects.
// \ingroup symmetric_matrix
//
// \param lhs The left-hand side NumericProxy object.
// \param rhs The right-hand side NumericProxy object.
// \return \a true if the left-hand side referenced value is smaller, \a false if not.
*/
template< typename MT1, typename MT2 >
inline bool operator<( const NumericProxy<MT1>& lhs, const NumericProxy<MT2>& rhs )
{
   return ( lhs.get() < rhs.get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Less-than comparison between a NumericProxy object and an object of different type.
// \ingroup symmetric_matrix
//
// \param lhs The left-hand side NumericProxy object.
// \param rhs The right-hand side object of other type.
// \return \a true if the left-hand side referenced value is smaller, \a false if not.
*/
template< typename MT, typename T >
inline bool operator<( const NumericProxy<MT>& lhs, const T& rhs )
{
   return ( lhs.get() < rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Less-than comparison between an object of different type and a NumericProxy object.
// \ingroup symmetric_matrix
//
// \param lhs The left-hand side object of other type.
// \param rhs The right-hand side NumericProxy object.
// \return \a true if the left-hand side other object is smaller, \a false if not.
*/
template< typename T, typename MT >
inline bool operator<( const T& lhs, const NumericProxy<MT>& rhs )
{
   return ( lhs < rhs.get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Greater-than comparison between two NumericProxy objects.
// \ingroup symmetric_matrix
//
// \param lhs The left-hand side NumericProxy object.
// \param rhs The right-hand side NumericProxy object.
// \return \a true if the left-hand side referenced value is greater, \a false if not.
*/
template< typename MT1, typename MT2 >
inline bool operator>( const NumericProxy<MT1>& lhs, const NumericProxy<MT2>& rhs )
{
   return ( lhs.get() > rhs.get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Greater-than comparison between a NumericProxy object and an object of different type.
// \ingroup symmetric_matrix
//
// \param lhs The left-hand side NumericProxy object.
// \param rhs The right-hand side object of other type.
// \return \a true if the left-hand side referenced value is greater, \a false if not.
*/
template< typename MT, typename T >
inline bool operator>( const NumericProxy<MT>& lhs, const T& rhs )
{
   return ( lhs.get() > rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Greater-than comparison between an object of different type and a NumericProxy object.
// \ingroup symmetric_matrix
//
// \param lhs The left-hand side object of other type.
// \param rhs The right-hand side NumericProxy object.
// \return \a true if the left-hand side other object is greater, \a false if not.
*/
template< typename T, typename MT >
inline bool operator>( const T& lhs, const NumericProxy<MT>& rhs )
{
   return ( lhs > rhs.get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Less-or-equal-than comparison between two NumericProxy objects.
// \ingroup symmetric_matrix
//
// \param lhs The left-hand side NumericProxy object.
// \param rhs The right-hand side NumericProxy object.
// \return \a true if the left-hand side referenced value is smaller or equal, \a false if not.
*/
template< typename MT1, typename MT2 >
inline bool operator<=( const NumericProxy<MT1>& lhs, const NumericProxy<MT2>& rhs )
{
   return ( lhs.get() <= rhs.get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Less-or-equal-than comparison between a NumericProxy object and an object of different type.
// \ingroup symmetric_matrix
//
// \param lhs The left-hand side NumericProxy object.
// \param rhs The right-hand side object of other type.
// \return \a true if the left-hand side referenced value is smaller or equal, \a false if not.
*/
template< typename MT, typename T >
inline bool operator<=( const NumericProxy<MT>& lhs, const T& rhs )
{
   return ( lhs.get() <= rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Less-or-equal-than comparison between an object of different type and a NumericProxy object.
// \ingroup symmetric_matrix
//
// \param lhs The left-hand side object of other type.
// \param rhs The right-hand side NumericProxy object.
// \return \a true if the left-hand side other object is smaller or equal, \a false if not.
*/
template< typename T, typename MT >
inline bool operator<=( const T& lhs, const NumericProxy<MT>& rhs )
{
   return ( lhs <= rhs.get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Greater-or-equal-than comparison between two NumericProxy objects.
// \ingroup symmetric_matrix
//
// \param lhs The left-hand side NumericProxy object.
// \param rhs The right-hand side NumericProxy object.
// \return \a true if the left-hand side referenced value is greater or equal, \a false if not.
*/
template< typename MT1, typename MT2 >
inline bool operator>=( const NumericProxy<MT1>& lhs, const NumericProxy<MT2>& rhs )
{
   return ( lhs.get() >= rhs.get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Greater-or-equal-than comparison between a NumericProxy object and an object of different type.
// \ingroup symmetric_matrix
//
// \param lhs The left-hand side NumericProxy object.
// \param rhs The right-hand side object of other type.
// \return \a true if the left-hand side referenced value is greater or equal, \a false if not.
*/
template< typename MT, typename T >
inline bool operator>=( const NumericProxy<MT>& lhs, const T& rhs )
{
   return ( lhs.get() >= rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Greater-or-equal-than comparison between an object of different type and a NumericProxy object.
// \ingroup symmetric_matrix
//
// \param lhs The left-hand side object of other type.
// \param rhs The right-hand side NumericProxy object.
// \return \a true if the left-hand side other object is greater or equal, \a false if not.
*/
template< typename T, typename MT >
inline bool operator>=( const T& lhs, const NumericProxy<MT>& rhs )
{
   return ( lhs >= rhs.get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Global output operator for numeric proxies.
// \ingroup symmetric_matrix
//
// \param os Reference to the output stream.
// \param proxy Reference to a constant numeric proxy object.
// \return Reference to the output stream.
*/
template< typename MT >
inline std::ostream& operator<<( std::ostream& os, const NumericProxy<MT>& proxy )
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
/*!\name NumericProxy global functions */
//@{
template< typename MT >
inline void reset( const NumericProxy<MT>& proxy );

template< typename MT >
inline void clear( const NumericProxy<MT>& proxy );

template< typename MT >
inline bool isDefault( const NumericProxy<MT>& proxy );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Resetting the represented element to the default initial values.
// \ingroup symmetric_matrix
//
// \param proxy The given access proxy.
// \return void
//
// This function resets the element represented by the numeric proxy to its default initial
// value.
*/
template< typename MT >
inline void reset( const NumericProxy<MT>& proxy )
{
   proxy.reset();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Clearing the represented element.
// \ingroup symmetric_matrix
//
// \param proxy The given access proxy.
// \return void
//
// This function clears the element represented by the numeric proxy to its default initial
// state.
*/
template< typename MT >
inline void clear( const NumericProxy<MT>& proxy )
{
   proxy.clear();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the represented element is in default state.
// \ingroup symmetric_matrix
//
// \param proxy The given access proxy
// \return \a true in case the represented element is in default state, \a false otherwise.
//
// This function checks whether the element represented by the access proxy is in default state.
// In case it is in default state, the function returns \a true, otherwise it returns \a false.
*/
template< typename MT >
inline bool isDefault( const NumericProxy<MT>& proxy )
{
   using blaze::isDefault;

   return isDefault( proxy.get() );
}
//*************************************************************************************************

} // namespace blaze

#endif
