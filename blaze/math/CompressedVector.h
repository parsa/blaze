//=================================================================================================
/*!
//  \file blaze/math/CompressedVector.h
//  \brief Implementation of an arbitrarily sized compressed vector
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

#ifndef _BLAZE_MATH_COMPRESSEDVECTOR_H_
#define _BLAZE_MATH_COMPRESSEDVECTOR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <algorithm>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <blaze/math/CMathTrait.h>
#include <blaze/math/CompressedMatrix.h>
#include <blaze/math/constraints/DenseVector.h>
#include <blaze/math/constraints/TransposeFlag.h>
#include <blaze/math/Functions.h>
#include <blaze/math/MathTrait.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/shims/IsNaN.h>
#include <blaze/math/sparse/SparseElement.h>
#include <blaze/math/sparse/VectorAccessProxy.h>
#include <blaze/math/SparseVector.h>
#include <blaze/math/Types.h>
#include <blaze/math/typetraits/CanAlias.h>
#include <blaze/math/typetraits/IsResizable.h>
#include <blaze/math/typetraits/IsSparseVector.h>
#include <blaze/system/Precision.h>
#include <blaze/system/TransposeFlag.h>
#include <blaze/util/Assert.h>
#include <blaze/util/constraints/Builtin.h>
#include <blaze/util/constraints/Const.h>
#include <blaze/util/constraints/FloatingPoint.h>
#include <blaze/util/constraints/Numeric.h>
#include <blaze/util/constraints/Pointer.h>
#include <blaze/util/constraints/Reference.h>
#include <blaze/util/constraints/SameSize.h>
#include <blaze/util/constraints/Volatile.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/Null.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/IsFloatingPoint.h>
#include <blaze/util/typetraits/IsNumeric.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\defgroup compressed_vector CompressedVector
// \ingroup sparse_vector
*/
/*!\brief Efficient implementation of an arbitrary sized sparse vector.
// \ingroup compressed_vector
//
// The CompressedVector class is the representation of an arbitrarily sized sparse vector,
// which stores only non-zero elements of arbitrary type. The type of the elements and the
// transpose flag of the vector can be specified via the two template parameters:

   \code
   template< typename Type, bool TF >
   class CompressedVector;
   \endcode

//  - Type: specifies the type of the vector elements. CompressedVector can be used with any
//          non-cv-qualified element type. The arithmetic operators for vector/vector and
//          vector/element operations with the same element type work for any element type
//          as long as the element type supports the arithmetic operation. Arithmetic operations
//          between vectors and elements of different element types are only supported for
//          all data types supported by the MathTrait class template (for details see the
//          MathTrait class description).
//  - TF  : specifies whether the vector is a row vector (\a blaze::rowVector) or a column
//          vector (\a blaze::columnVector). The default value is \a blaze::columnVector.
//
// Inserting/accessing elements in a compressed vector can be done by several alternative
// functions. The following example demonstrates all options:

   \code
   // Creating a non-transpose compressed vector of size 100
   CompressedVector<double,false> a( 100 );

   // The subscript operator provides access to all possible elements of the compressed vector,
   // including the zero elements. In case the subscript operator is used to access an element
   // that is currently not stored in the sparse vector, the element is inserted into the vector.
   a[42] = 2.0;

   // An alternative for inserting elements into the vector is the insert() function. However,
   // it inserts the element only in case the element is not already contained in the vector.
   A.insert( 50, 3.7 );

   // A very efficient way to add new elements to a sparse vector is the append() function.
   // Note that append() requires that the appended element's index is strictly larger than
   // the currently largest non-zero index of the vector and that the vector's capacity
   // is large enough to hold the new element.
   a.reserve( 10 );
   A.append( 51, -2.1 );

   // In order to traverse all non-zero elements currently stored in the vector, the begin()
   // and end() functions can be used. In the example, all non-zero elements of vector are
   // traversed.
   for( CompressedVector<double,false>::Iterator i=a.begin(); i!=a.end(); ++i ) {
      ... = i->value();  // Access to the value of the non-zero element
      ... = i->index();  // Access to the index of the non-zero element
   }
   \endcode

// The use of CompressedVector is very natural and intuitive. All operations (addition, subtraction,
// multiplication, scaling, ...) can be performed on all possible combinations of dense and sparse
// vectors with fitting element types. The following example gives an impression of the use of
// CompressedVector:

   \code
   using blaze::CompressedVector;
   using blaze::DynamicVector;
   using blaze::CompressedMatrix;

   CompressedVector<double> a( 2 );  // Default constructed, non-initialized 2D vectors
   a[0] = 1.0;                       // Initialization of the first element
   a[1] = 2.0;                       // Initialization of the second element

   CompressedVector<double> b( 2 );        // Empty sparse vector
   DynamicVector<float>     c( 2, 2.0F );  // Directly, homogeneously initialized dense vector
   CompressedVector<double> d;             // Default constructed dynamic vector
   CompressedMatrix<double> A;             // Default constructed row-major matrix

   d = a + b;  // Vector addition between vectors of equal element type
   d = a - c;  // Vector subtraction between a sparse and dense vector with different element types
   d = a * b;  // Component-wise vector multiplication

   a *= 2.0;      // In-place scaling of vector
   d  = a * 2.0;  // Scaling of vector a
   d  = 2.0 * a;  // Scaling of vector a

   d += a - b;  // Addition assignment
   d -= a + c;  // Subtraction assignment
   d *= a * b;  // Multiplication assignment

   double scalar = trans( a ) * b;  // Scalar/dot/inner product between two vectors

   A = a * trans( b );  // Outer product between two vectors
   \endcode


*/
template< typename Type                     // Data type of the vector
        , bool TF = defaultTransposeFlag >  // Transpose flag
class CompressedVector : public SparseVector< CompressedVector<Type,TF>, TF >
{
 private:
   //**Type definitions****************************************************************************
   typedef SparseElement<Type>  ElementBase;  //!< Base class for the compressed vector element.
   //**********************************************************************************************

   //**Private class Element***********************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Index-value-pair for the CompressedVector class.
   */
   struct Element : public ElementBase
   {
      using ElementBase::operator=;
      friend class CompressedVector;
   };
   /*! \endcond */
   //**********************************************************************************************

   //**Private class FindIndex*********************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Helper class for the lower_bound() function.
   */
   struct FindIndex : public std::binary_function<Element,size_t,bool>
   {
      inline bool operator()( const Element& element, size_t index ) const {
         return element.index() < index;
      }
      inline bool operator()( size_t index, const Element& element ) const {
         return index < element.index();
      }
      inline bool operator()( const Element& element1, const Element& element2 ) const {
         return element1.index() < element2.index();
      }
   };
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   typedef CompressedVector<Type,TF>   This;            //!< Type of this CompressedVector instance.
   typedef This                        ResultType;      //!< Result type for expression template evaluations.
   typedef CompressedVector<Type,!TF>  TransposeType;   //!< Transpose type for expression template evaluations.
   typedef Type                        ElementType;     //!< Type of the compressed vector elements.
   typedef const Type&                 ReturnType;      //!< Return type for expression template evaluations.
   typedef const CompressedVector&     CompositeType;   //!< Data type for composite expression templates.
   typedef VectorAccessProxy<This>     Reference;       //!< Reference to a non-constant vector value.
   typedef const Type&                 ConstReference;  //!< Reference to a constant vector value.
   typedef Element*                    Iterator;        //!< Iterator over non-constant elements.
   typedef const Element*              ConstIterator;   //!< Iterator over constant elements.

   //! Compressed vector length return type.
   /*! Return type of the CompressedVector<Type,TF>::length function. */
   typedef typename CMathTrait<Type>::Type  LengthType;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation flag for the detection of aliasing effects.
   /*! This compilation switch indicates whether this type potentially causes compuation errors
       due to aliasing effects. In case the type can cause aliasing effects, the \a canAlias
       switch is set to \a true, otherwise it is set to \a false. */
   enum { canAlias = 0 };
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
                           explicit inline CompressedVector();
                           explicit inline CompressedVector( size_t size );
                           explicit inline CompressedVector( size_t size, size_t nonzeros );
                                    inline CompressedVector( const CompressedVector& sv );
   template< typename VT >          inline CompressedVector( const DenseVector<VT,TF>&  dv );
   template< typename VT >          inline CompressedVector( const SparseVector<VT,TF>& sv );
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   inline ~CompressedVector();
   //@}
   //**********************************************************************************************

   //**Data access functions***********************************************************************
   /*!\name Data access functions */
   //@{
   inline Reference      operator[]( size_t index );
   inline ConstReference operator[]( size_t index ) const;
   inline Iterator       begin();
   inline ConstIterator  begin() const;
   inline Iterator       end();
   inline ConstIterator  end() const;
   //@}
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   /*!\name Assignment operators */
   //@{
                           inline CompressedVector& operator= ( const CompressedVector& rhs );
   template< typename VT > inline CompressedVector& operator= ( const DenseVector<VT,TF>&  rhs );
   template< typename VT > inline CompressedVector& operator= ( const SparseVector<VT,TF>& rhs );
   template< typename VT > inline CompressedVector& operator+=( const Vector<VT,TF>& rhs );
   template< typename VT > inline CompressedVector& operator-=( const Vector<VT,TF>& rhs );
   template< typename VT > inline CompressedVector& operator*=( const Vector<VT,TF>& rhs );

   template< typename Other >
   inline typename EnableIf< IsNumeric<Other>, CompressedVector >::Type&
      operator*=( Other rhs );

   template< typename Other >
   inline typename EnableIf< IsNumeric<Other>, CompressedVector >::Type&
      operator/=( Other rhs );
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
                              inline size_t                 size() const;
                              inline size_t                 capacity() const;
                              inline size_t                 nonZeros() const;
                              inline void                   reset();
                              inline void                   clear();
                                     void                   append( size_t index, const Type& value );
                                     Type&                  insert( size_t index, const Type& value );
                              inline Iterator               find  ( size_t index );
                              inline ConstIterator          find  ( size_t index ) const;
                              inline void                   resize( size_t n, bool preserve=true );
                              inline void                   reserve( size_t n );
                              inline LengthType             length() const;
                              inline const Type             sqrLength() const;
                              inline CompressedVector&      normalize();
                              inline const CompressedVector getNormalized() const;
   template< typename Other > inline CompressedVector&      scale( Other scalar );
                              inline void                   swap( CompressedVector& sv ) /* throw() */;
   //@}
   //**********************************************************************************************

   //**Expression template evaluation functions****************************************************
   /*!\name Expression template evaluation functions */
   //@{
   template< typename Other > inline bool isAliased( const Other* alias ) const;
   template< typename VT >    inline void assign   ( const DenseVector <VT,TF>& rhs );
   template< typename VT >    inline void assign   ( const SparseVector<VT,TF>& rhs );
   template< typename VT >    inline void addAssign( const DenseVector <VT,TF>& rhs );
   template< typename VT >    inline void addAssign( const SparseVector<VT,TF>& rhs );
   template< typename VT >    inline void subAssign( const DenseVector <VT,TF>& rhs );
   template< typename VT >    inline void subAssign( const SparseVector<VT,TF>& rhs );
   //@}
   //**********************************************************************************************

 private:
   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline size_t extendCapacity() const;
   //@}
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   size_t size_;             //!< The current size/dimension of the compressed vector.
   size_t capacity_;         //!< The maximum capacity of the compressed vector.
   Iterator begin_;          //!< Pointer to the first non-zero element of the compressed vector.
   Iterator end_;            //!< Pointer one past the last non-zero element of the compressed vector.

   static const Type zero_;  //!< Neutral element for accesses to zero elements.
   //@}
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_NOT_BE_POINTER_TYPE  ( Type );
   BLAZE_CONSTRAINT_MUST_NOT_BE_REFERENCE_TYPE( Type );
   BLAZE_CONSTRAINT_MUST_NOT_BE_CONST         ( Type );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VOLATILE      ( Type );
   BLAZE_CONSTRAINT_MUST_HAVE_SAME_SIZE( ElementBase, Element );
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  DEFINITION AND INITIALIZATION OF THE STATIC MEMBER VARIABLES
//
//=================================================================================================

template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
const Type CompressedVector<Type,TF>::zero_ = Type();




//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief The default constructor for CompressedVector.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline CompressedVector<Type,TF>::CompressedVector()
   : size_    ( 0UL )   // The current size/dimension of the compressed vector
   , capacity_( 0UL )   // The maximum capacity of the compressed vector
   , begin_   ( NULL )  // Pointer to the first non-zero element of the compressed vector
   , end_     ( NULL )  // Pointer to the last non-zero element of the compressed vector
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constructor for a compressed vector of size \a n.
//
// \param n The size of the vector.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline CompressedVector<Type,TF>::CompressedVector( size_t n )
   : size_    ( n   )   // The current size/dimension of the compressed vector
   , capacity_( 0UL )   // The maximum capacity of the compressed vector
   , begin_   ( NULL )  // Pointer to the first non-zero element of the compressed vector
   , end_     ( NULL )  // Pointer to the last non-zero element of the compressed vector
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constructor for a compressed vector of size \a n.
//
// \param n The size of the vector.
// \param nonzeros The number of expected non-zero elements.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline CompressedVector<Type,TF>::CompressedVector( size_t n, size_t nonzeros )
   : size_    ( n )                       // The current size/dimension of the compressed vector
   , capacity_( nonzeros )                // The maximum capacity of the compressed vector
   , begin_   ( new Element[capacity_] )  // Pointer to the first non-zero element of the compressed vector
   , end_     ( begin_ )                  // Pointer to the last non-zero element of the compressed vector
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The copy constructor for CompressedVector.
//
// \param sv Compressed vector to be copied.
//
// The copy constructor is explicitly defined due to the required dynamic memory management
// and in order to enable/facilitate NRV optimization.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline CompressedVector<Type,TF>::CompressedVector( const CompressedVector& sv )
   : size_    ( sv.size_ )                // The current size/dimension of the compressed vector
   , capacity_( sv.nonZeros() )           // The maximum capacity of the compressed vector
   , begin_   ( new Element[capacity_] )  // Pointer to the first non-zero element of the compressed vector
   , end_     ( begin_+capacity_ )        // Pointer to the last non-zero element of the compressed vector
{
   std::copy( sv.begin_, sv.end_, begin_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Conversion constructor from dense vectors.
//
// \param dv Dense vector to be copied.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
template< typename VT >  // Type of the foreign dense vector
inline CompressedVector<Type,TF>::CompressedVector( const DenseVector<VT,TF>& dv )
   : size_    ( (~dv).size() )  // The current size/dimension of the compressed vector
   , capacity_( 0UL  )          // The maximum capacity of the compressed vector
   , begin_   ( NULL )          // Pointer to the first non-zero element of the compressed vector
   , end_     ( NULL )          // Pointer to the last non-zero element of the compressed vector
{
   using blaze::assign;
   assign( *this, ~dv );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Conversion constructor from different sparse vectors.
//
// \param sv Sparse vector to be copied.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
template< typename VT >  // Type of the foreign sparse vector
inline CompressedVector<Type,TF>::CompressedVector( const SparseVector<VT,TF>& sv )
   : size_    ( (~sv).size() )            // The current size/dimension of the compressed vector
   , capacity_( (~sv).nonZeros() )        // The maximum capacity of the compressed vector
   , begin_   ( new Element[capacity_] )  // Pointer to the first non-zero element of the compressed vector
   , end_     ( begin_ )                  // Pointer to the last non-zero element of the compressed vector
{
   using blaze::assign;
   assign( *this, ~sv );
}
//*************************************************************************************************




//=================================================================================================
//
//  DESTRUCTOR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief The destructor for CompressedVector.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline CompressedVector<Type,TF>::~CompressedVector()
{
   delete [] begin_;
}
//*************************************************************************************************




//=================================================================================================
//
//  DATA ACCESS FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Subscript operator for the direct access to the compressed vector elements.
//
// \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
//
// This function returns a reference to the accessed value at position \a index. In case the
// compressed vector does not yet store an element for index \a index, a new element is inserted
// into the compressed vector. An alternative for traversing the non-zero elements of the sparse
// vector are the begin() and end() functions.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline typename CompressedVector<Type,TF>::Reference
   CompressedVector<Type,TF>::operator[]( size_t index )
{
   BLAZE_USER_ASSERT( index < size_, "Invalid compressed vector access index" );

   return Reference( *this, index );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Subscript operator for the direct access to the compressed vector elements.
//
// \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline typename CompressedVector<Type,TF>::ConstReference
   CompressedVector<Type,TF>::operator[]( size_t index ) const
{
   BLAZE_USER_ASSERT( index < size_, "Invalid compressed vector access index" );

   const Iterator pos( std::lower_bound( begin_, end_, index, FindIndex() ) );

   if( pos == end_ || pos->index_ != index )
      return zero_;
   else
      return pos->value_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first non-zero element of the compressed vector.
//
// \return Iterator to the first non-zero element of the compressed vector.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline typename CompressedVector<Type,TF>::Iterator CompressedVector<Type,TF>::begin()
{
   return Iterator( begin_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first non-zero element of the compressed vector.
//
// \return Iterator to the first non-zero element of the compressed vector.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline typename CompressedVector<Type,TF>::ConstIterator CompressedVector<Type,TF>::begin() const
{
   return ConstIterator( begin_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last non-zero element of the compressed vector.
//
// \return Iterator just past the last non-zero element of the compressed vector.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline typename CompressedVector<Type,TF>::Iterator CompressedVector<Type,TF>::end()
{
   return Iterator( end_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last non-zero element of the compressed vector.
//
// \return Iterator just past the last non-zero element of the compressed vector.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline typename CompressedVector<Type,TF>::ConstIterator CompressedVector<Type,TF>::end() const
{
   return ConstIterator( end_ );
}
//*************************************************************************************************




//=================================================================================================
//
//  ASSIGNMENT OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Copy assignment operator for CompressedVector.
//
// \param rhs Compressed vector to be copied.
// \return Reference to the assigned compressed vector.
//
// The compressed vector is resized according to the given compressed vector and initialized
// as a copy of this vector.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline CompressedVector<Type,TF>&
   CompressedVector<Type,TF>::operator=( const CompressedVector& rhs )
{
   if( &rhs == this ) return *this;

   const size_t nonzeros( rhs.nonZeros() );

   if( nonzeros > capacity_ ) {
      Iterator newBegin( new Element[nonzeros] );
      end_ = std::copy( rhs.begin_, rhs.end_, newBegin );
      std::swap( begin_, newBegin );
      delete [] newBegin;

      size_     = rhs.size_;
      capacity_ = nonzeros;
   }
   else {
      end_  = std::copy( rhs.begin_, rhs.end_, begin_ );
      size_ = rhs.size_;
   }

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Assignment operator for dense vectors.
//
// \param rhs Dense vector to be copied.
// \return Reference to the assigned compressed vector.
//
// The vector is resized according to the given dense vector and initialized as a copy of
// this vector.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
template< typename VT >  // Type of the right-hand side dense vector
inline CompressedVector<Type,TF>&
   CompressedVector<Type,TF>::operator=( const DenseVector<VT,TF>& rhs )
{
   using blaze::assign;

   if( CanAlias<VT>::value && (~rhs).isAliased( this ) ) {
      CompressedVector tmp( rhs );
      swap( tmp );
   }
   else {
      size_ = (~rhs).size();
      end_  = begin_;
      assign( *this, ~rhs );
   }

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Assignment operator for different sparse vectors.
//
// \param rhs Sparse vector to be copied.
// \return Reference to the assigned compressed vector.
//
// The vector is resized according to the given sparse vector and initialized as a copy of
// this vector.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
template< typename VT >  // Type of the right-hand side sparse vector
inline CompressedVector<Type,TF>&
   CompressedVector<Type,TF>::operator=( const SparseVector<VT,TF>& rhs )
{
   using blaze::assign;

   if( ( CanAlias<VT>::value && (~rhs).isAliased( this ) ) ||
       (~rhs).nonZeros() > capacity_ ) {
      CompressedVector tmp( rhs );
      swap( tmp );
   }
   else {
      size_ = (~rhs).size();
      end_  = begin_;
      assign( *this, ~rhs );
   }

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Addition assignment operator for the addition of a vector (\f$ \vec{a}+=\vec{b} \f$).
//
// \param rhs The right-hand side vector to be added to the compressed vector.
// \return Reference to the compressed vector.
// \exception std::invalid_argument Vector sizes do not match.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
template< typename VT >  // Type of the right-hand side vector
inline CompressedVector<Type,TF>& CompressedVector<Type,TF>::operator+=( const Vector<VT,TF>& rhs )
{
   using blaze::addAssign;

   if( (~rhs).size() != size_ )
      throw std::invalid_argument( "Vector sizes do not match" );

   addAssign( *this, ~rhs );

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Subtraction assignment operator for the subtraction of a vector (\f$ \vec{a}-=\vec{b} \f$).
//
// \param rhs The right-hand side vector to be subtracted from the compressed vector.
// \return Reference to the compressed vector.
// \exception std::invalid_argument Vector sizes do not match.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
template< typename VT >  // Type of the right-hand side vector
inline CompressedVector<Type,TF>& CompressedVector<Type,TF>::operator-=( const Vector<VT,TF>& rhs )
{
   using blaze::subAssign;

   if( (~rhs).size() != size_ )
      throw std::invalid_argument( "Vector sizes do not match" );

   subAssign( *this, ~rhs );

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication assignment operator for the multiplication of a vector
//        (\f$ \vec{a}*=\vec{b} \f$).
//
// \param rhs The right-hand side vector to be multiplied with the compressed vector.
// \return Reference to the compressed vector.
// \exception std::invalid_argument Vector sizes do not match.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
template< typename VT >  // Type of the right-hand side vector
inline CompressedVector<Type,TF>& CompressedVector<Type,TF>::operator*=( const Vector<VT,TF>& rhs )
{
   if( (~rhs).size() != size_ )
      throw std::invalid_argument( "Vector sizes do not match" );

   CompressedVector<Type,TF> tmp( *this * (~rhs) );
   swap( tmp );

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication assignment operator for the multiplication between a compressed vector
//        and a scalar value (\f$ \vec{a}*=s \f$).
//
// \param rhs The right-hand side scalar value for the multiplication.
// \return Reference to the compressed vector.
//
// This operator can only be used for built-in data types. Additionally, the elements of the
// compressed vector must support the multiplication assignment operator for the given scalar
// built-in data type.
*/
template< typename Type     // Data type of the vector
        , bool TF >         // Transpose flag
template< typename Other >  // Data type of the right-hand side scalar
inline typename EnableIf< IsNumeric<Other>, CompressedVector<Type,TF> >::Type&
   CompressedVector<Type,TF>::operator*=( Other rhs )
{
   for( Iterator element=begin_; element<end_; ++element )
      element->value_ *= rhs;
   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Division assignment operator for the division of a compressed vector by a scalar value
//        (\f$ \vec{a}/=s \f$).
//
// \param rhs The right-hand side scalar value for the division.
// \return Reference to the compressed vector.
//
// This operator can only be used for built-in data types. Additionally, the elements of the
// compressed vector must either support the multiplication assignment operator for the given
// floating point data type or the division assignment operator for the given integral data
// type.
*/
template< typename Type     // Data type of the vector
        , bool TF >         // Transpose flag
template< typename Other >  // Data type of the right-hand side scalar
inline typename EnableIf< IsNumeric<Other>, CompressedVector<Type,TF> >::Type&
   CompressedVector<Type,TF>::operator/=( Other rhs )
{
   BLAZE_USER_ASSERT( rhs != Other(0), "Division by zero detected" );

   typedef typename MathTrait<Type,Other>::DivType  DT;
   typedef typename If< IsNumeric<DT>, DT, Other >::Type  Tmp;

   // Depending on the two involved data types, an integer division is applied or a
   // floating point division is selected.
   if( IsNumeric<DT>::value && IsFloatingPoint<DT>::value ) {
      const Tmp tmp( Tmp(1)/static_cast<Tmp>( rhs ) );
      for( Iterator element=begin_; element!=end_; ++element )
         element->value_ *= tmp;
   }
   else {
      for( Iterator element=begin_; element!=end_; ++element )
         element->value_ /= rhs;
   }

   return *this;
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns the current size/dimension of the compressed vector.
//
// \return The size of the compressed vector.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline size_t CompressedVector<Type,TF>::size() const
{
   return size_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the maximum capacity of the compressed vector.
//
// \return The capacity of the compressed vector.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline size_t CompressedVector<Type,TF>::capacity() const
{
   return capacity_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the number of non-zero elements in the compressed vector.
//
// \return The number of non-zero elements in the compressed vector.
//
// Note that the number of non-zero elements is always smaller than the current size of the
// compressed vector.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline size_t CompressedVector<Type,TF>::nonZeros() const
{
   return end_ - begin_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Reset to the default initial values.
//
// \return void
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline void CompressedVector<Type,TF>::reset()
{
   end_ = begin_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Clearing the compressed vector.
//
// \return void
//
// After the clear() function, the size of the compressed vector is 0.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline void CompressedVector<Type,TF>::clear()
{
   size_ = 0UL;
   end_  = begin_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Appending an element to the compressed vector.
//
// \param index The index of the new element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be appended.
// \return void
//
// This function provides a very efficient way to fill a compressed vector with elements. It
// appends a new element to the end of the compressed vector without any additional check or
// memory allocation. Therefore it is strictly necessary to keep the following preconditions
// in mind:
//
//  - the index of the new element must be strictly larger than the largest index of non-zero
//    elements in the compressed vector
//  - the current number of non-zero elements must be smaller than the capacity of the vector
//
// Ignoring these preconditions might result in undefined behavior!
//
// \b Note: Although append() does not allocate new memory, it still invalidates all iterators
// returned by the end() functions!
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
void CompressedVector<Type,TF>::append( size_t index, const Type& value )
{
   BLAZE_USER_ASSERT( index < size_, "Invalid compressed vector access index" );
   BLAZE_USER_ASSERT( nonZeros() < capacity(), "Not enough reserved space" );
   BLAZE_USER_ASSERT( begin_ == end_ || (end_-1UL)->index_ < index, "Index is not strictly increasing" );

   end_->value_ = value;
   end_->index_ = index;
   ++end_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inserting an element into the compressed vector.
//
// \param index The index of the new element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be inserted.
// \return Reference to the inserted value.
// \exception std::invalid_argument Invalid compressed vector access index.
//
// This function inserts a new element into the compressed vector. However, duplicate elements
// are not allowed. In case the sparse vector already contains an element with index \a index,
// a \a std::invalid_argument exception is thrown.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
Type& CompressedVector<Type,TF>::insert( size_t index, const Type& value )
{
   BLAZE_USER_ASSERT( index < size_, "Invalid compressed vector access index" );

   const Iterator pos( std::lower_bound( begin_, end_, index, FindIndex() ) );

   if( pos != end_ && pos->index_ == index )
      throw std::invalid_argument( "Bad access index" );

   if( nonZeros() != capacity_ ) {
      std::copy_backward( pos, end_, end_+1 );
      pos->value_ = value;
      pos->index_ = index;
      ++end_;

      return pos->value_;
   }
   else {
      size_t newCapacity( extendCapacity() );

      Iterator newBegin = new Element[newCapacity];
      Iterator tmp      = std::copy( begin_, pos, newBegin );
      tmp->value_ = value;
      tmp->index_ = index;
      end_ = std::copy( pos, end_, tmp+1 );

      std::swap( newBegin, begin_ );
      delete [] newBegin;
      capacity_ = newCapacity;

      return tmp->value_;
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Searches for a specific vector element.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the element in case the index is found, end() iterator otherwise.
//
// This function can be used to check whether a specific element is contained in the sparse
// vector. It specifically searches for the element with index \a index. In case the element
// is found, the function returns an iterator to the element. Otherwise an iterator just past
// the last non-zero element of the compressed vector (the end() iterator) is returned. Note
// that the returned compressed vector iterator is subject to invalidation due to inserting
// operations via the subscript operator or the insert() function!
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline typename CompressedVector<Type,TF>::Iterator CompressedVector<Type,TF>::find( size_t index )
{
   return const_cast<Iterator>( const_cast<const This&>( *this ).find( index ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Searches for a specific vector element.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the element in case the index is found, end() iterator otherwise.
//
// This function can be used to check whether a specific element is contained in the sparse
// vector. It specifically searches for the element with index \a index. In case the element
// is found, the function returns an iterator to the element. Otherwise an iterator just past
// the last non-zero element of the compressed vector (the end() iterator) is returned. Note
// that the returned compressed vector iterator is subject to invalidation due to inserting
// operations via the subscript operator or the insert() function!
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline typename CompressedVector<Type,TF>::ConstIterator CompressedVector<Type,TF>::find( size_t index ) const
{
   const Iterator pos( std::lower_bound( begin_, end_, index, FindIndex() ) );
   if( pos != end_ && pos->index_ == index )
      return pos;
   else return end_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Changing the size of the compressed vector.
//
// \param n The new size of the compressed vector.
// \param preserve \a true if the old values of the vector should be preserved, \a false if not.
// \return void
//
// This function resizes the compressed vector using the given size to \a n. During this operation,
// new dynamic memory may be allocated in case the capacity of the compressed vector is too small.
// Therefore this function potentially changes all vector elements. In order to preserve the old
// vector values, the \a preserve flag can be set to \a true.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline void CompressedVector<Type,TF>::resize( size_t n, bool preserve )
{
   if( preserve ) {
      end_ = std::lower_bound( begin_, end_, n, FindIndex() );
   }
   else {
      end_ = begin_;
   }

   size_ = n;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Setting the minimum capacity of the compressed vector.
//
// \param n The new minimum capacity of the compressed vector.
// \return void
//
// This function increases the capacity of the compressed vector to at least \a n elements. The
// current values of the vector elements are preserved.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline void CompressedVector<Type,TF>::reserve( size_t n )
{
   if( n > capacity_ ) {
      const size_t newCapacity( n );

      // Allocating a new data and index array
      Iterator newBegin  = new Element[newCapacity];

      // Replacing the old data and index array
      end_ = std::copy( begin_, end_, newBegin );
      std::swap( newBegin, begin_ );
      capacity_ = newCapacity;
      delete [] newBegin;
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Calculation of the compressed vector length \f$|\vec{a}|\f$.
//
// \return The length of the compressed vector.
//
// This function calculates the actual length of the vector. The return type of the length()
// function depends on the actual type of the vector instance:
//
// <table border="0" cellspacing="0" cellpadding="1">
//    <tr>
//       <td width="250px"> \b Type </td>
//       <td width="100px"> \b LengthType </td>
//    </tr>
//    <tr>
//       <td>float</td>
//       <td>float</td>
//    </tr>
//    <tr>
//       <td>integral data types and double</td>
//       <td>double</td>
//    </tr>
//    <tr>
//       <td>long double</td>
//       <td>long double</td>
//    </tr>
// </table>
//
// \b Note: This operation is only defined for built-in data types. In case \a Type is a user
// defined data type the attempt to use the length() function results in a compile time error!
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
#ifndef WIN32
inline typename CompressedVector<Type,TF>::LengthType CompressedVector<Type,TF>::length() const
#else
inline typename CMathTrait<Type>::Type CompressedVector<Type,TF>::length() const
#endif
{
   BLAZE_CONSTRAINT_MUST_BE_BUILTIN_TYPE( Type );

   LengthType sum( 0 );
   for( Iterator element=begin_; element!=end_; ++element )
      sum += element->value_ * element->value_;
   return std::sqrt( sum );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Calculation of the vector square length \f$|\vec{a}|^2\f$.
//
// \return The square length of the vector.
//
// This function calculates the actual square length of the vector.
//
// \b Note: This operation is only defined for built-in data types. In case \a Type is a user
// defined data type the attempt to use the length() function results in a compile time error!
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline const Type CompressedVector<Type,TF>::sqrLength() const
{
   BLAZE_CONSTRAINT_MUST_BE_BUILTIN_TYPE( Type );

   Type sum( 0 );
   for( Iterator element=begin_; element!=end_; ++element )
      sum += element->value_ * element->value_;
   return sum;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Normalization of the compressed vector (\f$|\vec{a}|=1\f$).
//
// \return Reference to the compressed vector.
//
// Normalization of the compressed vector to a length of 1. This operation is only defined for
// floating point vectors. The attempt to use this function for an integral vector results
// in a compile time error.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline CompressedVector<Type,TF>& CompressedVector<Type,TF>::normalize()
{
   BLAZE_CONSTRAINT_MUST_BE_FLOATING_POINT_TYPE( Type );

   const Type len( length() );

   if( len == Type(0) )
      return *this;

   const Type ilen( Type(1) / len );

   for( Iterator element=begin_; element!=end_; ++element )
      element->value_ *= ilen;

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Calculation of the normalized compressed vector (\f$|\vec{a}|=1\f$).
//
// \return The normalized compressed vector.
//
// The function returns the normalized compressed vector. This operation is only defined for
// floating point vectors. The attempt to use this function for an integral vector results
// in a compile time error.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline const CompressedVector<Type,TF> CompressedVector<Type,TF>::getNormalized() const
{
   BLAZE_CONSTRAINT_MUST_BE_FLOATING_POINT_TYPE( Type );

   const Type len( length() );

   if( len == Type(0) )
      return *this;

   const Type ilen( Type(1) / len );
   CompressedVector tmp( *this );

   for( Iterator element=tmp.begin_; element!=tmp.end_; ++element ) {
      element->value_ *= ilen;
   }

   return tmp;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Scaling of the compressed vector by the scalar value \a scalar (\f$ \vec{a}=\vec{b}*s \f$).
//
// \param scalar The scalar value for the vector scaling.
// \return Reference to the compressed vector.
*/
template< typename Type     // Data type of the vector
        , bool TF >         // Transpose flag
template< typename Other >  // Data type of the scalar value
inline CompressedVector<Type,TF>& CompressedVector<Type,TF>::scale( Other scalar )
{
   for( Iterator element=begin_; element!=end_; ++element )
      element->value_ *= scalar;
   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Swapping the contents of two compressed vectors.
//
// \param sv The compressed vector to be swapped.
// \return void
// \exception no-throw guarantee.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline void CompressedVector<Type,TF>::swap( CompressedVector& sv ) /* throw() */
{
   std::swap( size_, sv.size_ );
   std::swap( capacity_, sv.capacity_ );
   std::swap( begin_, sv.begin_ );
   std::swap( end_, sv.end_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Calculating a new vector capacity.
//
// \return The new compressed vector capacity.
//
// This function calculates a new vector capacity based on the current capacity of the sparse
// vector. Note that the new capacity is restricted to the interval \f$[7..size]\f$.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline size_t CompressedVector<Type,TF>::extendCapacity() const
{
   using blaze::max;
   using blaze::min;

   size_t nonzeros( 2UL*capacity_+1UL );
   nonzeros = max( nonzeros, 7UL   );
   nonzeros = min( nonzeros, size_ );

   BLAZE_INTERNAL_ASSERT( nonzeros > capacity_, "Invalid capacity value" );

   return nonzeros;
}
//*************************************************************************************************




//=================================================================================================
//
//  EXPRESSION TEMPLATE EVALUATION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns whether the vector is aliased with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this vector, \a false if not.
*/
template< typename Type     // Data type of the vector
        , bool TF >         // Transpose flag
template< typename Other >  // Data type of the foreign expression
inline bool CompressedVector<Type,TF>::isAliased( const Other* alias ) const
{
   return static_cast<const void*>( this ) == static_cast<const void*>( alias );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the assignment of a dense vector.
//
// \param rhs The right-hand side dense vector to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
template< typename VT >  // Type of the right-hand side dense vector
inline void CompressedVector<Type,TF>::assign( const DenseVector<VT,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size_ == (~rhs).size(), "Invalid vector sizes" );

   size_t nonzeros( 0UL );

   for( size_t i=0UL; i<size_; ++i )
   {
      if( nonzeros == capacity_ )
         reserve( extendCapacity() );

      end_->value_ = (~rhs)[i];

      if( !isDefault( end_->value_ ) ) {
         end_->index_ = i;
         ++end_;
         ++nonzeros;
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the assignment of a sparse vector.
//
// \param rhs The right-hand side sparse vector to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
template< typename VT >  // Type of the right-hand side sparse vector
inline void CompressedVector<Type,TF>::assign( const SparseVector<VT,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size_ == (~rhs).size(), "Invalid vector sizes" );

   // Using the following formulation instead of a std::copy function call of the form
   //
   //          end_ = std::copy( (~rhs).begin(), (~rhs).end(), begin_ );
   //
   // results in much less requirements on the ConstIterator type provided from the right-hand
   // sparse vector type
   for( typename VT::ConstIterator element=(~rhs).begin(); element!=(~rhs).end(); ++element )
      append( element->index(), element->value() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the addition assignment of a dense vector.
//
// \param rhs The right-hand side dense vector to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
template< typename VT >  // Type of the right-hand side dense vector
inline void CompressedVector<Type,TF>::addAssign( const DenseVector<VT,TF>& rhs )
{
   typedef typename MathTrait<This,typename VT::ResultType>::AddType  AddType;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( AddType );
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( AddType, TF );
   BLAZE_CONSTRAINT_MUST_BE_REFERENCE_TYPE( typename AddType::CompositeType );

   BLAZE_INTERNAL_ASSERT( size_ == (~rhs).size(), "Invalid vector sizes" );

   AddType tmp( *this + (~rhs) );
   reset();
   assign( tmp );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the addition assignment of a sparse vector.
//
// \param rhs The right-hand side sparse vector to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
template< typename VT >  // Type of the right-hand side sparse vector
inline void CompressedVector<Type,TF>::addAssign( const SparseVector<VT,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size_ == (~rhs).size(), "Invalid vector sizes" );

   CompressedVector<Type,TF> tmp( *this + (~rhs) );
   swap( tmp );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the subtraction assignment of a dense vector.
//
// \param rhs The right-hand side dense vector to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
template< typename VT >  // Type of the right-hand side dense vector
inline void CompressedVector<Type,TF>::subAssign( const DenseVector<VT,TF>& rhs )
{
   typedef typename MathTrait<This,typename VT::ResultType>::SubType  SubType;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( SubType );
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( SubType, TF );
   BLAZE_CONSTRAINT_MUST_BE_REFERENCE_TYPE( typename SubType::CompositeType );

   BLAZE_INTERNAL_ASSERT( size_ == (~rhs).size(), "Invalid vector sizes" );

   SubType tmp( *this - (~rhs) );
   reset();
   assign( tmp );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the subtraction assignment of a sparse vector.
//
// \param rhs The right-hand side sparse vector to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
template< typename VT >  // Type of the right-hand side sparse vector
inline void CompressedVector<Type,TF>::subAssign( const SparseVector<VT,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size_ == (~rhs).size(), "Invalid vector sizes" );

   CompressedVector<Type,TF> tmp( *this - (~rhs) );
   swap( tmp );
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name CompressedVector operators */
//@{
template< typename Type, bool TF >
inline void reset( CompressedVector<Type,TF>& v );

template< typename Type, bool TF >
inline void clear( CompressedVector<Type,TF>& v );

template< typename Type, bool TF >
inline bool isnan( const CompressedVector<Type,TF>& v );

template< typename Type, bool TF >
inline bool isDefault( const CompressedVector<Type,TF>& v );

template< typename Type, bool TF >
inline const typename MathTrait<Type,Type>::MultType sq( const CompressedVector<Type,TF>& v );

template< typename Type, bool TF >
inline void swap( CompressedVector<Type,TF>& a, CompressedVector<Type,TF>& b ) /* throw() */;
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Resetting the given compressed vector.
// \ingroup compressed_vector
//
// \param v The compressed vector to be resetted.
// \return void
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline void reset( CompressedVector<Type,TF>& v )
{
   v.reset();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Clearing the given compressed vector.
// \ingroup compressed_vector
//
// \param v The compressed vector to be cleared.
// \return void
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline void clear( CompressedVector<Type,TF>& v )
{
   v.clear();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks the given compressed vector for not-a-number elements.
// \ingroup compressed_vector
//
// \param v The compressed vector to be checked for not-a-number elements.
// \return \a true if at least one element of the vector is not-a-number, \a false otherwise.
//
// This function checks the N-dimensional compressed vector for not-a-number (NaN) elements. If
// at least one element of the vector is not-a-number, the function returns \a true, otherwise
// it returns \a false.

   \code
   blaze::CompressedVector<double> a;
   // ... Resizing and initialization
   if( isnan( a ) ) { ... }
   \endcode
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline bool isnan( const CompressedVector<Type,TF>& v )
{
   typedef typename CompressedVector<Type,TF>::ConstIterator  ConstIterator;

   const ConstIterator end( v.end() );
   for( ConstIterator element=v.begin(); element!=end; ++element ) {
      if( isnan( element->value() ) ) return true;
   }
   return false;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the given compressed vector is in default state.
// \ingroup compressed_vector
//
// \param v The compressed vector to be tested for its default state.
// \return \a true in case the given vector is component-wise zero, \a false otherwise.
//
// This function checks whether the N-dimensional compressed vector is in default state. For
// instance, in case the vector is instantiated for a built-in integral or floating point data
// type, the function returns \a true in case all vector elements are 0 and \a false in case
// any vector element is not 0. The following example demonstrates the use of the \a isDefault
// function:

   \code
   blaze::CompressedVector<double> a;
   // ... Resizing and initialization
   if( isDefault( a ) ) { ... }
   \endcode
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline bool isDefault( const CompressedVector<Type,TF>& v )
{
   typedef typename CompressedVector<Type,TF>::ConstIterator  ConstIterator;

   for( ConstIterator element=v.begin(); element!=v.end(); ++element )
      if( !isDefault( element->value() ) ) return false;
   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Squaring the given compressed vector.
// \ingroup compressed_vector
//
// \param v The compressed vector to be squared.
// \return The result of the square operation.
//
// This function calculates the component product of the given compressed vector. It has the
// same effect as multiplying the vector with itself (\f$ v * v \f$).
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline const CompressedVector< typename MathTrait<Type,Type>::MultType, TF >
   sq( const CompressedVector<Type>& v )
{
   return v * v;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Swapping the contents of two compressed vectors.
// \ingroup compressed_vector
//
// \param a The first compressed vector to be swapped.
// \param b The second compressed vector to be swapped.
// \return void
// \exception no-throw guarantee.
*/
template< typename Type  // Data type of the vector
        , bool TF >      // Transpose flag
inline void swap( CompressedVector<Type,TF>& a, CompressedVector<Type,TF>& b ) /* throw() */
{
   a.swap( b );
}
//*************************************************************************************************




//=================================================================================================
//
//  ISRESIZABLE SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T, bool TF >
struct IsResizable< CompressedVector<T,TF> > : public TrueType
{
   enum { value = 1 };
   typedef TrueType  Type;
};

template< typename T, bool TF >
struct IsResizable< const CompressedVector<T,TF> > : public TrueType
{
   enum { value = 1 };
   typedef TrueType  Type;
};

template< typename T, bool TF >
struct IsResizable< volatile CompressedVector<T,TF> > : public TrueType
{
   enum { value = 1 };
   typedef TrueType  Type;
};

template< typename T, bool TF >
struct IsResizable< const volatile CompressedVector<T,TF> > : public TrueType
{
   enum { value = 1 };
   typedef TrueType  Type;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  MATHTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T1, bool TF, typename T2 >
struct MathTrait< CompressedVector<T1,TF>, T2 >
{
   typedef INVALID_TYPE                                                 HighType;
   typedef INVALID_TYPE                                                 LowType;
   typedef INVALID_TYPE                                                 AddType;
   typedef INVALID_TYPE                                                 SubType;
   typedef CompressedVector< typename MathTrait<T1,T2>::MultType, TF >  MultType;
   typedef INVALID_TYPE                                                 CrossType;
   typedef CompressedVector< typename MathTrait<T1,T2>::DivType , TF >  DivType;
   BLAZE_CONSTRAINT_MUST_BE_NUMERIC_TYPE( T2 );
};

template< typename T1, typename T2, bool TF >
struct MathTrait< T1, CompressedVector<T2,TF> >
{
   typedef INVALID_TYPE                                                 HighType;
   typedef INVALID_TYPE                                                 LowType;
   typedef INVALID_TYPE                                                 AddType;
   typedef INVALID_TYPE                                                 SubType;
   typedef CompressedVector< typename MathTrait<T1,T2>::MultType, TF >  MultType;
   typedef INVALID_TYPE                                                 CrossType;
   typedef INVALID_TYPE                                                 DivType;
   BLAZE_CONSTRAINT_MUST_BE_NUMERIC_TYPE( T1 );
};

template< typename T1, bool TF, typename T2, size_t N >
struct MathTrait< CompressedVector<T1,TF>, StaticVector<T2,N,TF> >
{
   typedef INVALID_TYPE                                                 HighType;
   typedef INVALID_TYPE                                                 LowType;
   typedef StaticVector< typename MathTrait<T1,T2>::AddType, N, TF >    AddType;
   typedef StaticVector< typename MathTrait<T1,T2>::SubType, N, TF >    SubType;
   typedef CompressedVector< typename MathTrait<T1,T2>::MultType, TF >  MultType;
   typedef INVALID_TYPE                                                 CrossType;
   typedef INVALID_TYPE                                                 DivType;
};

template< typename T1, typename T2 >
struct MathTrait< CompressedVector<T1,false>, StaticVector<T2,3UL,false> >
{
 private:
   typedef typename MathTrait<T1,T2>::MultType  T;

 public:
   typedef INVALID_TYPE                                                    HighType;
   typedef INVALID_TYPE                                                    LowType;
   typedef StaticVector< typename MathTrait<T1,T2>::AddType, 3UL, false >  AddType;
   typedef StaticVector< typename MathTrait<T1,T2>::SubType, 3UL, false >  SubType;
   typedef CompressedVector< typename MathTrait<T1,T2>::MultType, false >  MultType;
   typedef StaticVector< typename MathTrait<T,T>::SubType  , 3UL, false >  CrossType;
   typedef INVALID_TYPE                                                    DivType;
};

template< typename T1, typename T2, size_t N >
struct MathTrait< CompressedVector<T1,false>, StaticVector<T2,N,true> >
{
   typedef INVALID_TYPE                                                   HighType;
   typedef INVALID_TYPE                                                   LowType;
   typedef INVALID_TYPE                                                   AddType;
   typedef INVALID_TYPE                                                   SubType;
   typedef CompressedMatrix< typename MathTrait<T1,T2>::MultType, true >  MultType;
   typedef INVALID_TYPE                                                   CrossType;
   typedef INVALID_TYPE                                                   DivType;
};

template< typename T1, typename T2, size_t N >
struct MathTrait< CompressedVector<T1,true>, StaticVector<T2,N,false> >
{
   typedef INVALID_TYPE                         HighType;
   typedef INVALID_TYPE                         LowType;
   typedef INVALID_TYPE                         AddType;
   typedef INVALID_TYPE                         SubType;
   typedef typename MathTrait<T1,T2>::MultType  MultType;
   typedef INVALID_TYPE                         CrossType;
   typedef INVALID_TYPE                         DivType;
};

template< typename T1, size_t N, bool TF, typename T2 >
struct MathTrait< StaticVector<T1,N,TF>, CompressedVector<T2,TF> >
{
   typedef INVALID_TYPE                                                 HighType;
   typedef INVALID_TYPE                                                 LowType;
   typedef StaticVector< typename MathTrait<T1,T2>::AddType, N, TF >    AddType;
   typedef StaticVector< typename MathTrait<T1,T2>::SubType, N, TF >    SubType;
   typedef CompressedVector< typename MathTrait<T1,T2>::MultType, TF >  MultType;
   typedef INVALID_TYPE                                                 CrossType;
   typedef INVALID_TYPE                                                 DivType;
};

template< typename T1, typename T2 >
struct MathTrait< StaticVector<T1,3UL,false>, CompressedVector<T2,false> >
{
 private:
   typedef typename MathTrait<T1,T2>::MultType  T;

 public:
   typedef INVALID_TYPE                                                    HighType;
   typedef INVALID_TYPE                                                    LowType;
   typedef StaticVector< typename MathTrait<T1,T2>::AddType, 3UL, false >  AddType;
   typedef StaticVector< typename MathTrait<T1,T2>::SubType, 3UL, false >  SubType;
   typedef CompressedVector< typename MathTrait<T1,T2>::MultType, false >  MultType;
   typedef StaticVector< typename MathTrait<T,T>::SubType  , 3UL, false >  CrossType;
   typedef INVALID_TYPE                                                    DivType;
};

template< typename T1, size_t N, typename T2 >
struct MathTrait< StaticVector<T1,N,false>, CompressedVector<T2,true> >
{
   typedef INVALID_TYPE                                                    HighType;
   typedef INVALID_TYPE                                                    LowType;
   typedef INVALID_TYPE                                                    AddType;
   typedef INVALID_TYPE                                                    SubType;
   typedef CompressedMatrix< typename MathTrait<T1,T2>::MultType, false >  MultType;
   typedef INVALID_TYPE                                                    CrossType;
   typedef INVALID_TYPE                                                    DivType;
};

template< typename T1, size_t N, typename T2 >
struct MathTrait< StaticVector<T1,N,true>, CompressedVector<T2,false> >
{
   typedef INVALID_TYPE                         HighType;
   typedef INVALID_TYPE                         LowType;
   typedef INVALID_TYPE                         AddType;
   typedef INVALID_TYPE                         SubType;
   typedef typename MathTrait<T1,T2>::MultType  MultType;
   typedef INVALID_TYPE                         CrossType;
   typedef INVALID_TYPE                         DivType;
};

template< typename T1, typename T2 >
struct MathTrait< CompressedVector<T1,false>, DynamicVector<T2,false> >
{
 private:
   typedef typename MathTrait<T1,T2>::MultType  T;

 public:
   typedef INVALID_TYPE                                                    HighType;
   typedef INVALID_TYPE                                                    LowType;
   typedef DynamicVector< typename MathTrait<T1,T2>::AddType , false >     AddType;
   typedef DynamicVector< typename MathTrait<T1,T2>::SubType , false >     SubType;
   typedef CompressedVector< typename MathTrait<T1,T2>::MultType, false >  MultType;
   typedef StaticVector< typename MathTrait<T,T>::SubType, 3UL, false >    CrossType;
   typedef INVALID_TYPE                                                    DivType;
};

template< typename T1, typename T2 >
struct MathTrait< CompressedVector<T1,false>, DynamicVector<T2,true> >
{
   typedef INVALID_TYPE                                                   HighType;
   typedef INVALID_TYPE                                                   LowType;
   typedef INVALID_TYPE                                                   AddType;
   typedef INVALID_TYPE                                                   SubType;
   typedef CompressedMatrix< typename MathTrait<T1,T2>::MultType, true >  MultType;
   typedef INVALID_TYPE                                                   CrossType;
   typedef INVALID_TYPE                                                   DivType;
};

template< typename T1, typename T2 >
struct MathTrait< CompressedVector<T1,true>, DynamicVector<T2,false> >
{
   typedef INVALID_TYPE                         HighType;
   typedef INVALID_TYPE                         LowType;
   typedef INVALID_TYPE                         AddType;
   typedef INVALID_TYPE                         SubType;
   typedef typename MathTrait<T1,T2>::MultType  MultType;
   typedef INVALID_TYPE                         CrossType;
   typedef INVALID_TYPE                         DivType;
};

template< typename T1, typename T2 >
struct MathTrait< CompressedVector<T1,true>, DynamicVector<T2,true> >
{
   typedef INVALID_TYPE                                                   HighType;
   typedef INVALID_TYPE                                                   LowType;
   typedef DynamicVector< typename MathTrait<T1,T2>::AddType , true >     AddType;
   typedef DynamicVector< typename MathTrait<T1,T2>::SubType , true >     SubType;
   typedef CompressedVector< typename MathTrait<T1,T2>::MultType, true >  MultType;
   typedef INVALID_TYPE                                                   CrossType;
   typedef INVALID_TYPE                                                   DivType;
};

template< typename T1, typename T2 >
struct MathTrait< DynamicVector<T1,false>, CompressedVector<T2,false> >
{
 private:
   typedef typename MathTrait<T1,T2>::MultType  T;

 public:
   typedef INVALID_TYPE                                                    HighType;
   typedef INVALID_TYPE                                                    LowType;
   typedef DynamicVector< typename MathTrait<T1,T2>::AddType, false >      AddType;
   typedef DynamicVector< typename MathTrait<T1,T2>::SubType, false >      SubType;
   typedef CompressedVector< typename MathTrait<T1,T2>::MultType, false >  MultType;
   typedef StaticVector< typename MathTrait<T,T>::SubType, 3UL, false >    CrossType;
   typedef INVALID_TYPE                                                    DivType;
};

template< typename T1, typename T2 >
struct MathTrait< DynamicVector<T1,false>, CompressedVector<T2,true> >
{
   typedef INVALID_TYPE                                                    HighType;
   typedef INVALID_TYPE                                                    LowType;
   typedef INVALID_TYPE                                                    AddType;
   typedef INVALID_TYPE                                                    SubType;
   typedef CompressedMatrix< typename MathTrait<T1,T2>::MultType, false >  MultType;
   typedef INVALID_TYPE                                                    CrossType;
   typedef INVALID_TYPE                                                    DivType;
};

template< typename T1, typename T2 >
struct MathTrait< DynamicVector<T1,true>, CompressedVector<T2,false> >
{
   typedef INVALID_TYPE                         HighType;
   typedef INVALID_TYPE                         LowType;
   typedef INVALID_TYPE                         AddType;
   typedef INVALID_TYPE                         SubType;
   typedef typename MathTrait<T1,T2>::MultType  MultType;
   typedef INVALID_TYPE                         CrossType;
   typedef INVALID_TYPE                         DivType;
};

template< typename T1, typename T2 >
struct MathTrait< DynamicVector<T1,true>, CompressedVector<T2,true> >
{
   typedef INVALID_TYPE                                                   HighType;
   typedef INVALID_TYPE                                                   LowType;
   typedef DynamicVector< typename MathTrait<T1,T2>::AddType, true >      AddType;
   typedef DynamicVector< typename MathTrait<T1,T2>::SubType, true >      SubType;
   typedef CompressedVector< typename MathTrait<T1,T2>::MultType, true >  MultType;
   typedef INVALID_TYPE                                                   CrossType;
   typedef INVALID_TYPE                                                   DivType;
};

template< typename T1, typename T2 >
struct MathTrait< CompressedVector<T1,false>, CompressedVector<T2,false> >
{
 private:
   typedef typename MathTrait<T1,T2>::MultType  T;

 public:
   typedef CompressedVector< typename MathTrait<T1,T2>::HighType, false >  HighType;
   typedef CompressedVector< typename MathTrait<T1,T2>::LowType , false >  LowType;
   typedef CompressedVector< typename MathTrait<T1,T2>::AddType , false >  AddType;
   typedef CompressedVector< typename MathTrait<T1,T2>::SubType , false >  SubType;
   typedef CompressedVector< typename MathTrait<T1,T2>::MultType, false >  MultType;
   typedef StaticVector< typename MathTrait<T,T>::SubType, 3UL, false >    CrossType;
   typedef INVALID_TYPE                                                    DivType;
};

template< typename T1, typename T2 >
struct MathTrait< CompressedVector<T1,false>, CompressedVector<T2,true> >
{
   typedef INVALID_TYPE                                                    HighType;
   typedef INVALID_TYPE                                                    LowType;
   typedef INVALID_TYPE                                                    AddType;
   typedef INVALID_TYPE                                                    SubType;
   typedef CompressedMatrix< typename MathTrait<T1,T2>::MultType, false >  MultType;
   typedef INVALID_TYPE                                                    CrossType;
   typedef INVALID_TYPE                                                    DivType;
};

template< typename T1, typename T2 >
struct MathTrait< CompressedVector<T1,true>, CompressedVector<T2,false> >
{
   typedef INVALID_TYPE                         HighType;
   typedef INVALID_TYPE                         LowType;
   typedef INVALID_TYPE                         AddType;
   typedef INVALID_TYPE                         SubType;
   typedef typename MathTrait<T1,T2>::MultType  MultType;
   typedef INVALID_TYPE                         CrossType;
   typedef INVALID_TYPE                         DivType;
};

template< typename T1, typename T2 >
struct MathTrait< CompressedVector<T1,true>, CompressedVector<T2,true> >
{
   typedef CompressedVector< typename MathTrait<T1,T2>::HighType, true >  HighType;
   typedef CompressedVector< typename MathTrait<T1,T2>::LowType , true >  LowType;
   typedef CompressedVector< typename MathTrait<T1,T2>::AddType , true >  AddType;
   typedef CompressedVector< typename MathTrait<T1,T2>::SubType , true >  SubType;
   typedef CompressedVector< typename MathTrait<T1,T2>::MultType, true >  MultType;
   typedef INVALID_TYPE                                                   CrossType;
   typedef INVALID_TYPE                                                   DivType;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  TYPE DEFINITIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Compressed single precision vector.
// \ingroup compressed_vector
*/
typedef CompressedVector<float,false>  CVecNf;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compressed double precision vector.
// \ingroup compressed_vector
*/
typedef CompressedVector<double,false>  CVecNd;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compressed vector with system-specific precision.
// \ingroup compressed_vector
*/
typedef CompressedVector<real,false>  CVecN;
//*************************************************************************************************

} // namespace blaze

#endif
