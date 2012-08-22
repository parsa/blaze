//=================================================================================================
/*!
//  \file blaze/math/sparse/MatrixAccessProxy.h
//  \brief Header file for the MatrixAccessProxy class
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

#ifndef _BLAZE_MATH_SPARSE_MATRIXACCESSPROXY_H_
#define _BLAZE_MATH_SPARSE_MATRIXACCESSPROXY_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/typetraits/IsRowMajorMatrix.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Access proxy for sparse, \f$ M \times N \f$ matrices.
// \ingroup math
//
// The MatrixAccessProxy provides safe access to the elements of a non-const sparse matrices.\n
// The proxied access to the elements of a sparse matrix is necessary since it may be possible
// that several insertion operations happen in the same statement. The following code illustrates
// this with two examples by means of the CompressedMatrix class:

   \code
   CompressedMatrix<real> A( 4, 4 );

   // Standard usage of the function call operator to initialize a matrix element.
   // Only a single sparse matrix element is accessed!
   A(0,1) = 1.0;

   // Initialization of a matrix element via another matrix element.
   // Two sparse matrix accesses in one statement!

   // Multiple accesses to elements of the sparse matrix in one statement!
   const real result = A(0,2) + A(1,2) + A(2,2);
   \endcode

// The problem (especially with the last statement) is that several insertion operations might
// take place due to the access via the function call operator. If the function call operator
// would return a direct reference to one of the accessed elements, this reference might be
// invalidated during the evaluation of a subsequent function call operator, which results in
// undefined behavior. This class provides the necessary functionality to guarantee a safe access
// to the sparse matrix elements while preserving the intuitive use of the function call operator.
*/
template< typename MT >  // Type of the sparse matrix
class MatrixAccessProxy
{
 private:
   //**Enumerations********************************************************************************
   //! Compile time flag indicating whether the given matrix type is a row-major matrix.
   enum { rmm = IsRowMajorMatrix<MT>::value };
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   typedef MT                        MatrixType;      //!< Type of the accessed sparse matrix.
   typedef typename MT::ElementType  ElementType;     //!< Type of the accessed sparse matrix element.
   typedef ElementType&              Reference;       //!< Reference type of the accessed element.
   typedef const ElementType&        ConstReference;  //!< Reference type of the accessed constant element.
   typedef typename MT::Iterator     Iterator;        //!< Iterator type of the accessed sparse matrix.
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline MatrixAccessProxy( MT& sv, size_t i, size_t j );
            inline MatrixAccessProxy( const MatrixAccessProxy& vap );
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   // No explicitly declared destructor.
   //**********************************************************************************************

   //**Operators***********************************************************************************
   /*!\name Operators */
   //@{
                          inline MatrixAccessProxy& operator= ( const MatrixAccessProxy& vap );
   template< typename T > inline MatrixAccessProxy& operator= ( const T& value );
   template< typename T > inline MatrixAccessProxy& operator+=( const T& value );
   template< typename T > inline MatrixAccessProxy& operator-=( const T& value );
   template< typename T > inline MatrixAccessProxy& operator*=( const T& value );
   template< typename T > inline MatrixAccessProxy& operator/=( const T& value );
   //@}
   //**********************************************************************************************

   //**Conversion operator*************************************************************************
   /*!\name Conversion operator */
   //@{
   inline operator Reference() const;
   //@}
   //**********************************************************************************************

 private:
   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   Reference get()                       const;
   void      set( ConstReference value ) const;
   //@}
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   MT&    sm_;  //!< Reference to the accessed sparse matrix.
   size_t i_;   //!< Row-index of the accessed sparse matrix element.
   size_t j_;   //!< Column-index of the accessed sparse matrix element.
   //@}
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Initialization constructor for a MatrixAccessProxy.
//
// \param sv Reference to the accessed sparse matrix.
// \param i The row-index of the accessed sparse matrix element.
// \param j The column-index of the accessed sparse matrix element.
*/
template< typename MT >  // Type of the sparse matrix
inline MatrixAccessProxy<MT>::MatrixAccessProxy( MT& sv, size_t i, size_t j )
   : sm_( sv )  // Reference to the accessed sparse matrix
   , i_ ( i  )  // Row-index of the accessed sparse matrix element
   , j_ ( j  )  // Column-index of the accessed sparse matrix element
{
   const Iterator element( sm_.find( i_, j_ ) );
   if( element == ( rmm ? sm_.end(i_) : sm_.end(j_) ) )
      sm_.insert( i_, j_, ElementType() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The copy constructor for MatrixAccessProxy.
//
// \param map Sparse matrix access proxy to be copied.
*/
template< typename MT >  // Type of the sparse matrix
inline MatrixAccessProxy<MT>::MatrixAccessProxy( const MatrixAccessProxy& map )
   : sm_( map.sm_ )  // Reference to the accessed sparse matrix
   , i_ ( map.i_  )  // Row-index of the accessed sparse matrix element
   , j_ ( map.j_  )  // Column-index of the accessed sparse matrix element
{
   BLAZE_INTERNAL_ASSERT( sm_.find(i_,j_) != ( rmm ? sm_.end(i_) : sm_.end(j_) ), "Missing matrix element detected" );
}
//*************************************************************************************************




//=================================================================================================
//
//  OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Copy assignment operator for MatrixAccessProxy.
//
// \param map Sparse matrix access proxy to be copied.
// \return Reference to the assigned access proxy.
*/
template< typename MT >  // Type of the sparse matrix
inline MatrixAccessProxy<MT>& MatrixAccessProxy<MT>::operator=( const MatrixAccessProxy& map )
{
   set( map.get() );
   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Assignment to the accessed sparse matrix element.
//
// \param value The new value of the sparse matrix element.
// \return Reference to the assigned access proxy.
*/
template< typename MT >  // Type of the sparse matrix
template< typename T >   // Type of the right-hand side value
inline MatrixAccessProxy<MT>& MatrixAccessProxy<MT>::operator=( const T& value )
{
   set( value );
   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Addition assignment to the accessed sparse matrix element.
//
// \param value The right-hand side value to be added to the sparse matrix element.
// \return Reference to the assigned access proxy.
*/
template< typename MT >  // Type of the sparse matrix
template< typename T >   // Type of the right-hand side value
inline MatrixAccessProxy<MT>& MatrixAccessProxy<MT>::operator+=( const T& value )
{
   get() += value;
   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Subtraction assignment to the accessed sparse matrix element.
//
// \param value The right-hand side value to be subtracted from the sparse matrix element.
// \return Reference to the assigned access proxy.
*/
template< typename MT >  // Type of the sparse matrix
template< typename T >   // Type of the right-hand side value
inline MatrixAccessProxy<MT>& MatrixAccessProxy<MT>::operator-=( const T& value )
{
   get() -= value;
   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication assignment to the accessed sparse matrix element.
//
// \param value The right-hand side value for the multiplication.
// \return Reference to the assigned access proxy.
*/
template< typename MT >  // Type of the sparse matrix
template< typename T >   // Type of the right-hand side value
inline MatrixAccessProxy<MT>& MatrixAccessProxy<MT>::operator*=( const T& value )
{
   get() *= value;
   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Division assignment to the accessed sparse matrix element.
//
// \param value The right-hand side value for the division.
// \return Reference to the assigned access proxy.
*/
template< typename MT >  // Type of the sparse matrix
template< typename T >   // Type of the right-hand side value
inline MatrixAccessProxy<MT>& MatrixAccessProxy<MT>::operator/=( const T& value )
{
   get() /= value;
   return *this;
}
//*************************************************************************************************




//=================================================================================================
//
//  CONVERSION OPERATOR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Conversion to the accessed sparse matrix element.
//
// \return Reference to the accessed sparse matrix element.
*/
template< typename MT >  // Type of the sparse matrix
inline MatrixAccessProxy<MT>::operator Reference() const
{
   return get();
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returning the value of the accessed sparse matrix element.
//
// \return Reference to the accessed sparse matrix element.
*/
template< typename MT >  // Type of the sparse matrix
inline typename MatrixAccessProxy<MT>::Reference MatrixAccessProxy<MT>::get() const
{
   const Iterator element( sm_.find( i_, j_ ) );
   BLAZE_INTERNAL_ASSERT( element != ( rmm ? sm_.end(i_) : sm_.end(j_) ), "Missing matrix element detected" );
   return element->value();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Setting the value of the accessed sparse matrix element.
//
// \param value Reference to the new value of the sparse matrix element.
// \return void
*/
template< typename MT >  // Type of the sparse matrix
inline void MatrixAccessProxy<MT>::set( ConstReference value ) const
{
   const Iterator element( sm_.find( i_, j_ ) );
   BLAZE_INTERNAL_ASSERT( element != ( rmm ? sm_.end(i_) : sm_.end(j_) ), "Missing matrix element detected" );
   element->value() = value;
}
//*************************************************************************************************

} // namespace blaze

#endif
