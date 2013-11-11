//=================================================================================================
/*!
//  \file blaze/math/Quaternion.h
//  \brief Header file for the implementation of a quaternion
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

#ifndef _BLAZE_MATH_QUATERNION_H_
#define _BLAZE_MATH_QUATERNION_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <cmath>
#include <istream>
#include <ostream>
#include <blaze/math/Accuracy.h>
#include <blaze/math/Forward.h>
#include <blaze/math/shims/Equal.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/shims/IsNaN.h>
#include <blaze/math/traits/MathTrait.h>
#include <blaze/math/traits/MultTrait.h>
#include <blaze/system/Precision.h>
#include <blaze/util/Assert.h>
#include <blaze/util/constraints/Const.h>
#include <blaze/util/constraints/FloatingPoint.h>
#include <blaze/util/constraints/Volatile.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\defgroup quaternion Quaternion
// \ingroup math
*/
/*!\brief Efficient implementation of a quaternion.
// \ingroup quaternion
//
// Quaternions are a superior way to deal with rotations and orientations. This quaternion
// consists of 4 statically allocated elements, where the first element represents the real
// part and the three other elements represent the three imaginary parts. The naming
// convention of the elements is as following:

                             \f[\left(\begin{array}{*{4}{c}}
                             r & i & j & k \\
                             \end{array}\right)\f]

// These elements can be accessed directly with the subscript operator. The numbering of the
// quaternion elements is

                             \f[\left(\begin{array}{*{4}{c}}
                             0 & 1 & 2 & 3 \\
                             \end{array}\right)\f]

// \b Note: The Quaternion class can only be instantiated for non-cv-qualified floating point
// types! Therefore the only possible Quaternion instantiations are
//
//  - Quaternion<float>
//  - Quaternion<double>
//  - Quaternion<long double>
//
// The attempt to create a quaternion with an integral data type results in a compile time
// error.
*/
template< typename Type >  // Data type of the quaternion
class Quaternion
{
 public:
   //**Type definitions****************************************************************************
   typedef Type  ElementType;  //!< Type of the quaternion elements.
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline Quaternion();

   explicit inline Quaternion( Type r, Type i, Type j, Type k );

   template< typename Axis >
   explicit inline Quaternion( StaticVector<Axis,3UL,false> axis, Type angle );

   explicit inline Quaternion( Type xangle, Type yangle, Type zangle );

   template< typename Other >
   explicit inline Quaternion( const StaticVector<Other,3UL,false>& euler );

   inline Quaternion( const Quaternion& q );

   template< typename Other >
   inline Quaternion( const Quaternion<Other>& q );
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   // No explicitly declared destructor.
   //**********************************************************************************************

   //**Operators***********************************************************************************
   /*!\name Operators */
   //@{
                              inline Quaternion& operator= ( const Quaternion& rhs );
   template< typename Other > inline Quaternion& operator= ( const Quaternion<Other>& rhs );
                              inline Type        operator[]( size_t index ) const;
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
                              inline Quaternion&                set( Type r, Type i, Type j, Type k );
                              inline void                       reset();
                              inline Type                       length() const;
                              inline Quaternion&                normalize();
                              inline const Quaternion           getNormalized() const;
                              inline Quaternion&                invert();
                              inline const RotationMatrix<Type> toRotationMatrix() const;
                              inline void                       rotateX( Type angle );
                              inline void                       rotateY( Type angle );
                              inline void                       rotateZ( Type angle );
                              inline void                       swap( Quaternion& q ) /* throw() */;
   //@}
   //**********************************************************************************************

   //**Math functions******************************************************************************
   /*!\name Math functions */
   //@{
   template< typename Other, bool TF >
   inline const StaticVector<typename MultTrait<Type,Other>::Type,3UL,false>
      rotate( const StaticVector<Other,3UL,TF>& v ) const;

   template< typename Other >
   inline const StaticMatrix<typename MultTrait<Type,Other>::Type,3UL,3UL,false>
      rotate( const StaticMatrix<Other,3UL,3UL,false>& m ) const;

   template< typename Other >
   inline const StaticMatrix<typename MultTrait<Type,Other>::Type,3UL,3UL,false>
      diagRotate( const StaticMatrix<Other,3UL,3UL,false>& m ) const;

   template< typename Other >
   inline typename MathTrait<Type,Other>::HighType
      calcAngle( const StaticVector<Other,3UL,false>& axis ) const;
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   Type v_[4];  //!< The four statically allocated quaternion elements.
                /*!< Access to the quaternion values is gained via the subscript operator.
                     The order of the elements is
                     \f[\left(\begin{array}{*{4}{c}}
                     0 & 1 & 2 & 3 \\
                     \end{array}\right)\f] */
   //@}
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_FLOATING_POINT_TYPE( Type );
   BLAZE_CONSTRAINT_MUST_NOT_BE_CONST          ( Type );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VOLATILE       ( Type );
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
/*!\brief The default constructor for Quaternion.
//
// The real part of the quaternion is initialized with 1, whereas the imaginary parts are
// initialized with 0:

           \f[ \left(\begin{array}{c} 1 \\ 0 \\ 0 \\ 0 \end{array}\right) \f]
*/
template< typename Type >  // Data type of the quaternion
inline Quaternion<Type>::Quaternion()
{
   reset();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constructor for a direct initialization of all quaternion elements.
//
// \param r The initial value for the real part.
// \param i The initial value for the first imaginary part.
// \param j The initial value for the second imaginary part.
// \param k The initial value for the third imaginary part.
//
// The initial values for the quaterion have to be chosen such that the length of the
// quaternion is 1.
*/
template< typename Type >  // Data type of the quaternion
inline Quaternion<Type>::Quaternion( Type r, Type i, Type j, Type k )
{
   v_[0] = r; v_[1] = i; v_[2] = j; v_[3] = k;
   BLAZE_USER_ASSERT( std::fabs( r*r + i*i + j*j + k*k - Type(1) ) < Type(accuracy), "Invalid quaternion parameters" );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constructor for a quaternion depending on a rotation axis and angle.
//
// \param axis The rotation axis.
// \param angle The rotation angle (radian measure).
//
// This constructor creates a quaternion from the rotation axis \a axis and the rotation angle
// \a angle. \a axis may be an arbitrary, non-zero vector of any length. However, it is allowed
// to use the zero vector (0,0,0) in combination with an angle of 0. This combination results
// in a default quaternion

           \f[ \left(\begin{array}{c} 1 \\ 0 \\ 0 \\ 0 \end{array}\right) \f]
*/
template< typename Type >  // Data type of the quaternion
template< typename Axis >  // Data type of the rotation axis
inline Quaternion<Type>::Quaternion( StaticVector<Axis,3UL,false> axis, Type angle )
{
   BLAZE_CONSTRAINT_MUST_BE_FLOATING_POINT_TYPE( Axis );

   if( std::fabs(angle) < real(1E-15) ) {
      reset();
      return;
   }

   BLAZE_USER_ASSERT( axis.sqrLength() > Axis(0), "Invalid rotation axis" );

   const Type sina( std::sin( angle*Type(0.5) ) );
   const Type cosa( std::cos( angle*Type(0.5) ) );

   axis.normalize();

   v_[0] = cosa;
   v_[1] = sina * axis[0];
   v_[2] = sina * axis[1];
   v_[3] = sina * axis[2];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constructor for a quaternion rotated by the Euler angles \a xangle, \a yangle and \a zangle.
//
// \param xangle Rotation around the x-axis (radian measure).
// \param yangle Rotation around the y-axis (radian measure).
// \param zangle Rotation around the z-axis (radian measure).
//
// This constructor creates a quaternion rotated by the Euler angles \a xangle, \a yangle and
// \a zangle (all in radian measure). The rotations are applied in the order x, y, and z.
*/
template< typename Type >  // Data type of the quaternion
inline Quaternion<Type>::Quaternion( Type xangle, Type yangle, Type zangle )
{
   reset();
   rotateX( xangle );
   rotateY( yangle );
   rotateZ( zangle );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constructor for a quaternion rotated by the Euler angles \a euler.
//
// \param euler 3-dimensional vector of the three rotation angles (radian measure).
//
// This constructor creates a quaternion rotated by the Euler angles \a euler (all components
// in radian measure). The rotations are applied in the order x, y, and z.
*/
template< typename Type >   // Data type of the quaternion
template< typename Other >  // Data type of the Euler angle vector
inline Quaternion<Type>::Quaternion( const StaticVector<Other,3UL,false>& euler )
{
   reset();
   rotateX( euler[0] );
   rotateY( euler[1] );
   rotateZ( euler[2] );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The copy constructor for Quaternion.
//
// \param q Quaternion to be copied.
//
// The copy constructor is explicitly defined in order to enable/facilitate NRV optimization.
*/
template< typename Type >  // Data type of the quaternion
inline Quaternion<Type>::Quaternion( const Quaternion& q )
{
   v_[0] = q.v_[0];
   v_[1] = q.v_[1];
   v_[2] = q.v_[2];
   v_[3] = q.v_[3];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Conversion constructor from different Quaternion instances.
//
// \param q Quaternion to be copied.
*/
template< typename Type >   // Data type of the quaternion
template< typename Other >  // Data type of the foreign quaternion
inline Quaternion<Type>::Quaternion( const Quaternion<Other>& q )
{
   v_[0] = q[0];
   v_[1] = q[1];
   v_[2] = q[2];
   v_[3] = q[3];
}
//*************************************************************************************************




//=================================================================================================
//
//  OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Copy assignment operator for Quaternion.
//
// \param rhs Quaternion to be copied.
// \return Reference to the assigned quaternion.
//
// Explicit definition of a copy assignment operator for performance reasons.
*/
template< typename Type >  // Data type of the quaternion
inline Quaternion<Type>& Quaternion<Type>::operator=( const Quaternion<Type>& rhs )
{
   // This implementation is faster than the synthesized default copy assignment operator and
   // faster than an implementation with the C library function 'memcpy' in combination with a
   // protection against self-assignment. Additionally, this version goes without a protection
   // against self-assignment.
   v_[0] = rhs.v_[0];
   v_[1] = rhs.v_[1];
   v_[2] = rhs.v_[2];
   v_[3] = rhs.v_[3];
   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Assignment operator for different Quaternion instances.
//
// \param rhs Quaternion to be copied.
// \return Reference to the assigned quaternion.
*/
template< typename Type >   // Data type of the quaternion
template< typename Other >  // Data type of the foreign quaternion
inline Quaternion<Type>& Quaternion<Type>::operator=( const Quaternion<Other>& rhs )
{
   // This implementation is faster than the synthesized default copy assignment operator and
   // faster than an implementation with the C library function 'memcpy' in combination with a
   // protection against self-assignment. Additionally, this version goes without a protection
   // against self-assignment.
   v_[0] = rhs[0];
   v_[1] = rhs[1];
   v_[2] = rhs[2];
   v_[3] = rhs[3];
   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Subscript operator for the direct access to the quaternion elements.
//
// \param index Access index. The index has to be in the range \f$[0..3]\f$.
// \return Copy of the accessed element.
//
// In case BLAZE_USER_ASSERT() is active, this operator performs an index check.
*/
template< typename Type >  // Data type of the quaternion
inline Type Quaternion<Type>::operator[]( size_t index ) const
{
   BLAZE_USER_ASSERT( index < 4, "Invalid quaternion access index" );
   return v_[index];
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Setting the value of the quaternion elements.
//
// \param r The value for the real part.
// \param i The value for the first imaginary part.
// \param j The value for the second imaginary part.
// \param k The value for the third imaginary part.
*/
template< typename Type >  // Data type of the quaternion
inline Quaternion<Type>& Quaternion<Type>::set( Type r, Type i, Type j, Type k )
{
   BLAZE_USER_ASSERT( std::fabs( r*r + i*i + j*j + k*k - Type(1) ) < Type(1E-8), "Invalid quaternion parameters" );
   v_[0] = r;
   v_[1] = i;
   v_[2] = j;
   v_[3] = k;
   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Reset to the default initial values.
//
// \return void
//
// This function resets the quaternion to the default initial values. The real part of the
// quaternion is reset to 1, whereas the imaginary parts are reset to 0:

                             \f[\left(\begin{array}{*{4}{c}}
                             1 & 0 & 0 & 0 \\
                             \end{array}\right)\f]
*/
template< typename Type >  // Data type of the quaternion
inline void Quaternion<Type>::reset()
{
   v_[0] = Type(1);
   v_[1] = v_[2] = v_[3] = Type(0);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Calculation of the quaternion length \f$|\hat{q}|\f$.
//
// \return The length of the quaternion.
*/
template< typename Type >  // Data type of the quaternion
inline Type Quaternion<Type>::length() const
{
   // Although the length of the quaternion should always be exactly one, the function
   // calculates the actual length to enable length checks.
   return std::sqrt( v_[0]*v_[0] + v_[1]*v_[1] + v_[2]*v_[2] + v_[3]*v_[3] );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Normalization of the quaternion (\f$|\hat{q}|=1\f$).
//
// \return Reference to the quaternion.
*/
template< typename Type >  // Data type of the quaternion
inline Quaternion<Type>& Quaternion<Type>::normalize()
{
   const Type len( std::sqrt( v_[0]*v_[0] + v_[1]*v_[1] + v_[2]*v_[2] + v_[3]*v_[3] ) );

   if( len == Type(0) )
      return *this;

   const Type ilen( Type(1)/len );

   v_[0] *= ilen;
   v_[1] *= ilen;
   v_[2] *= ilen;
   v_[3] *= ilen;

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Calculation of the normalized quaternion (\f$|\hat{q}|=1\f$).
//
// \return The normalized quaternion.
*/
template< typename Type >  // Data type of the quaternion
inline const Quaternion<Type> Quaternion<Type>::getNormalized() const
{
   const Type len( std::sqrt( v_[0]*v_[0] + v_[1]*v_[1] + v_[2]*v_[2] + v_[3]*v_[3] ) );

   if( len == Type(0) )
      return *this;

   const Type ilen( Type(1)/len );

   return Quaternion( ilen*v_[0], ilen*v_[1], ilen*v_[2], ilen*v_[3] );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inversion of the quaternion (\f$ \hat{q} = \hat{q}^{-1} \f$).
//
// \return Reference to the inverted quaternion.
*/
template< typename Type >  // Data type of the quaternion
inline Quaternion<Type>& Quaternion<Type>::invert()
{
   v_[1] *= -Type(1);
   v_[2] *= -Type(1);
   v_[3] *= -Type(1);
   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Conversion to a rotation matrix.
//
// \return The resulting rotation matrix.
*/
template< typename Type >  // Data type of the quaternion
inline const RotationMatrix<Type> Quaternion<Type>::toRotationMatrix() const
{
   return RotationMatrix<Type>( Type(1) - Type(2)*v_[2]*v_[2] - Type(2)*v_[3]*v_[3],
                                Type(2)*( v_[1]*v_[2] - v_[0]*v_[3] ),
                                Type(2)*( v_[1]*v_[3] + v_[0]*v_[2] ),
                                Type(2)*( v_[1]*v_[2] + v_[0]*v_[3] ),
                                Type(1) - Type(2)*v_[1]*v_[1] - Type(2)*v_[3]*v_[3],
                                Type(2)*( v_[2]*v_[3] - v_[0]*v_[1] ),
                                Type(2)*( v_[1]*v_[3] - v_[0]*v_[2] ),
                                Type(2)*( v_[2]*v_[3] + v_[0]*v_[1] ),
                                Type(1) - Type(2)*v_[1]*v_[1] - Type(2)*v_[2]*v_[2] );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Rotating the quaternion around the global x-axis by \a angle degrees (radian measure).
//
// \param angle The rotation angle (radian measure).
// \return void
*/
template< typename Type >  // Data type of the quaternion
inline void Quaternion<Type>::rotateX( Type angle )
{
   const Type sina( std::sin( angle*Type(0.5) ) );
   const Type cosa( std::cos( angle*Type(0.5) ) );

   const Quaternion q( cosa*v_[0] - sina*v_[1],
                       cosa*v_[1] + sina*v_[0],
                       cosa*v_[2] - sina*v_[3],
                       cosa*v_[3] + sina*v_[2] );

   this->operator=( q );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Rotating the quaternion around the global y-axis by \a angle degrees (radian measure).
//
// \param angle The rotation angle (radian measure).
// \return void
*/
template< typename Type >  // Data type of the quaternion
inline void Quaternion<Type>::rotateY( Type angle )
{
   const Type sina( std::sin( angle*Type(0.5) ) );
   const Type cosa( std::cos( angle*Type(0.5) ) );

   const Quaternion q( cosa*v_[0] - sina*v_[2],
                       cosa*v_[1] + sina*v_[3],
                       cosa*v_[2] + sina*v_[0],
                       cosa*v_[3] - sina*v_[1] );

   this->operator=( q );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Rotating the quaternion around the global z-axis by \a angle degrees (radian measure).
//
// \param angle The rotation angle (radian measure).
// \return void
*/
template< typename Type >  // Data type of the quaternion
inline void Quaternion<Type>::rotateZ( Type angle )
{
   const Type sina( std::sin( angle*Type(0.5) ) );
   const Type cosa( std::cos( angle*Type(0.5) ) );

   const Quaternion q( cosa*v_[0] - sina*v_[3],
                       cosa*v_[1] - sina*v_[2],
                       cosa*v_[2] + sina*v_[1],
                       cosa*v_[3] + sina*v_[0] );

   this->operator=( q );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Swapping the contents of two quaternions.
//
// \param q The quaternion to be swapped.
// \return void
// \exception no-throw guarantee.
*/
template< typename Type >  // Data type of the quaternion
inline void Quaternion<Type>::swap( Quaternion& q ) /* throw() */
{
   std::swap( v_[0], q.v_[0] );
   std::swap( v_[1], q.v_[1] );
   std::swap( v_[2], q.v_[2] );
   std::swap( v_[3], q.v_[3] );
}
//*************************************************************************************************




//=================================================================================================
//
//  MATH FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Rotation of a vector v (\f$ \vec{rot} = \hat{q} \cdot \vec{v} \cdot \hat{q}^{-1} \f$).
//
// \param v The vector to be rotated.
// \return The rotated vector.
//
// The function is selected for vectors of different data type (in case \a Type and \a Other
// are supported by the MultTrait class). The function returns a vector of the higher-order
// data type of the two involved data types.
*/
template< typename Type >  // Data type of the quaternion
template< typename Other   // Data type of the vector
        , bool TF >        // Transpose flag
inline const StaticVector<typename MultTrait<Type,Other>::Type,3UL,false>
   Quaternion<Type>::rotate( const StaticVector<Other,3UL,TF>& v ) const
{
   typedef typename MultTrait<Type,Other>::Type  MT;

   // Multiplication in two steps
   const MT w( v_[1]*v[0] + v_[2]*v[1] + v_[3]*v[2] );
   const MT x( v_[0]*v[0] - v_[3]*v[1] + v_[2]*v[2] );
   const MT y( v_[0]*v[1] - v_[1]*v[2] + v_[3]*v[0] );
   const MT z( v_[0]*v[2] - v_[2]*v[0] + v_[1]*v[1] );

   return StaticVector<MT,3UL,TF>( v_[0]*x + v_[1]*w + v_[2]*z - v_[3]*y,
                                   v_[0]*y + v_[2]*w + v_[3]*x - v_[1]*z,
                                   v_[0]*z + v_[3]*w + v_[1]*y - v_[2]*x );

   // Multiplication in one step
   /*
   const MT v0( v_[0] * v_[0] );
   const MT v1( v_[1] * v_[1] );
   const MT v2( v_[2] * v_[2] );
   const MT v3( v_[3] * v_[3] );

   return StaticVector<MT,3UL,TF>( ( v[0]*( v0 + v1 - v2 - v3 ) + 2.0*v_[0]*( v_[2]*v[2] - v_[3]*v[1] ) + 2.0*v_[1]*( v_[2]*v[1] + v_[3]*v[2] ) ),
                                   ( v[1]*( v0 - v1 + v2 - v3 ) + 2.0*v_[0]*( v_[3]*v[0] - v_[1]*v[2] ) + 2.0*v_[2]*( v_[1]*v[0] + v_[3]*v[2] ) ),
                                   ( v[2]*( v0 - v1 - v2 + v3 ) + 2.0*v_[0]*( v_[1]*v[1] - v_[2]*v[0] ) + 2.0*v_[3]*( v_[1]*v[0] + v_[2]*v[1] ) ) );
   */
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Rotation of a matrix.
//
// \param m The matrix to be rotated.
// \return The rotated matrix.
//
// The function is selected for matrices of different data type (in case \a Type and \a Other
// are supported by the MultTrait class). The function returns a matrix of the higher-order
// data type of the two involved data types.
*/
template< typename Type >   // Data type of the quaternion
template< typename Other >  // Data type of the matrix
inline const StaticMatrix< typename MultTrait<Type,Other>::Type, 3UL, 3UL, false >
   Quaternion<Type>::rotate( const StaticMatrix<Other,3UL,3UL,false>& m ) const
{
   typedef typename MultTrait<Type,Other>::Type  MT;
   const RotationMatrix<MT> R( this->toRotationMatrix() );
   return R.rotate( m );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Rotation of a diagonal matrix.
//
// \param m The diagonal matrix to be rotated.
// \return The rotated matrix.
//
// The DiagRotate function is a special case of the rotate function. The matrix is assumed to
// be a diagonal matrix, which reduces the number of floating point operations of the rotation.
// The function is selected for matrices of different data type (in case \a Type and \a Other
// are supported by the MultTrait class). The function returns a matrix of the higher-order
// data type of the two involved data types.
*/
template< typename Type >   // Data type of the quaternion
template< typename Other >  // Data type of the diagonal matrix
inline const StaticMatrix< typename MultTrait<Type,Other>::Type, 3UL, 3UL, false >
   Quaternion<Type>::diagRotate( const StaticMatrix<Other,3UL,3UL,false>& m ) const
{
   typedef typename MultTrait<Type,Other>::Type  MT;
   const RotationMatrix<MT> R( this->toRotationMatrix() );
   return R.diagRotate( m );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the angle in radian measure between the quaterion and a given axis.
//
// \param axis The given axis.
// \return The angle in radian measure.
*/
template< typename Type >   // Data type of the quaternion
template< typename Other >  // Data type of the axis
inline typename MathTrait<Type,Other>::HighType
   Quaternion<Type>::calcAngle( const StaticVector<Other,3UL,false>& axis ) const
{
   typedef typename MathTrait<Type,Other>::HighType  High;

   const StaticVector<High,3UL,true> u( v_[1], v_[2], v_[3] );
   const High y  ( u.length() );
   const High x  ( v_[0] );
   const High dot( u * axis );

   return High(2) * std::atan2( y, ( dot < real(0) )?( -x ):( x ) );
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name Quaternion operators */
//@{
template< typename T1, typename T2 >
inline bool operator==( const Quaternion<T1>& lhs, const Quaternion<T2>& rhs );

template< typename T1, typename T2 >
inline bool operator!=( const Quaternion<T1>& lhs, const Quaternion<T2>& rhs );

template< typename Type >
std::ostream& operator<<( std::ostream& os, const Quaternion<Type>& q );

template< typename Type >
std::istream& operator>>( std::istream& is, Quaternion<Type>& q );

template< typename Type >
inline bool isnan( const Quaternion<Type>& q );

template< typename Type >
inline void reset( Quaternion<Type>& q );

template< typename Type >
inline void clear( Quaternion<Type>& q );

template< typename Type >
inline bool isDefault( const Quaternion<Type>& q );

template< typename Type >
inline const Quaternion<Type> inv( const Quaternion<Type>& m );

template< typename Type >
inline const Quaternion<Type> sq( const Quaternion<Type>& m );

template< typename Type >
inline void swap( Quaternion<Type>& a, Quaternion<Type>& b ) /* throw() */;
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality operator for the comparison of two quaternions.
// \ingroup quaternion
//
// \param lhs The left-hand side quaternion for the comparison.
// \param rhs The right-hand side quaternion for the comparison.
// \return \a true if the two quaternions are equal, \a false if not.
*/
template< typename T1    // Data type of the left-hand side quaternion
        , typename T2 >  // Data type of the right-hand side quaternion
inline bool operator==( const Quaternion<T1>& lhs, const Quaternion<T2>& rhs )
{
   // In order to compare the two quaternions, the data values of the lower-order data
   // type are converted to the higher-order data type within the equal function.
   if( !equal( lhs[0], rhs[0] ) ||
       !equal( lhs[1], rhs[1] ) ||
       !equal( lhs[2], rhs[2] ) ||
       !equal( lhs[2], rhs[2] ) )
      return false;
   else return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inequality operator for the comparison of two quaternions.
// \ingroup quaternion
//
// \param lhs The left-hand side quaternion for the comparison.
// \param rhs The right-hand side quaternion for the comparison.
// \return \a true if the two quaternions are not equal, \a false if they are equal.
*/
template< typename T1    // Data type of the left-hand side quaternion
        , typename T2 >  // Data type of the right-hand side quaternion
inline bool operator!=( const Quaternion<T1>& lhs, const Quaternion<T2>& rhs )
{
   return !( lhs == rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Global output operator for quaternions.
// \ingroup quaternion
//
// \param os Reference to the output stream.
// \param q Reference to a constant quaternion object.
// \return Reference to the output stream.
*/
template< typename Type >  // Data type of the quaternion
std::ostream& operator<<( std::ostream& os, const Quaternion<Type>& q )
{
   return os << "<" << q[0] << "," << q[1] << "," << q[2] << "," << q[3] << ">";
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Global input operator for quaternions.
// \ingroup quaternion
//
// \param is Reference to the input stream.
// \param q Reference to a quaternion object.
// \return The input stream.
*/
template< typename Type >  // Data type of the quaternion
std::istream& operator>>( std::istream& is, Quaternion<Type>& q )
{
   if( !is ) return is;

   char bracket1, bracket2, comma1, comma2, comma3;
   Type r(0), i(0), j(0), k(0);
   const std::istream::pos_type pos( is.tellg() );
   const std::istream::fmtflags oldFlags( is.flags() );

   // Setting the 'skip whitespaces' flag
   is >> std::skipws;

   // Extracting the quaternion
   if( !(is >> bracket1 >> r >> comma1 >> i >> comma2 >> j >> comma3 >> k >> bracket2) ||
       bracket1 != '<' || comma1 != ',' || comma2 != ',' || comma3 != ',' || bracket2 != '>' ) {
      is.clear();
      is.seekg( pos );
      is.setstate( std::istream::failbit );
      is.flags( oldFlags );
      return is;
   }

   // Transfering the input to the quaternion values
   q.set( r, i, j, k );

   // Resetting the flags
   is.flags( oldFlags );

   return is;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks the given quaternion for not-a-number elements.
// \ingroup quaternion
//
// \param q The quaternion to be checked for not-a-number elements.
// \return \a true if at least one element of the quaternion is not-a-number, \a false otherwise.
*/
template< typename Type >  // Data type of the quaternion
inline bool isnan( const Quaternion<Type>& q )
{
   if( isnan( q[0] ) || isnan( q[1] ) || isnan( q[2] ) || isnan( q[3] ) )
      return true;
   else return false;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Resetting the given quaternion.
// \ingroup quaternion
//
// \param q The quaternion to be resetted.
// \return void
*/
template< typename Type >  // Data type of the quaternion
inline void reset( Quaternion<Type>& q )
{
   q.reset();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Clearing the given quaternion.
// \ingroup quaternion
//
// \param q The quaternion to be cleared.
// \return void
//
// Clearing a quaternion is equivalent to resetting it via the reset() function.
*/
template< typename Type >  // Data type of the quaternion
inline void clear( Quaternion<Type>& q )
{
   q.reset();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the given quaternion is in default state.
// \ingroup quaternion
//
// \param q The quaternion to be tested for its default state.
// \return \a true in case the given quaternion is in default state, \a false otherwise.
//
// The function returns \a true in case the real part of the quaternion is 1 and the imaginary
// parts are 0, otherwise it returns \a false.

                             \f[\left(\begin{array}{*{4}{c}}
                             1 & 0 & 0 & 0 \\
                             \end{array}\right)\f]
*/
template< typename Type >  // Data type of the quaternion
inline bool isDefault( const Quaternion<Type>& q )
{
   return ( q[0] == Type(1) ) && ( q[1] == Type(0) ) && ( q[2] == Type(0) ) && ( q[3] == Type(0) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inverting the given quaternion (\f$ \hat{q} = \hat{p}^{-1} \f$).
// \ingroup quaternion
//
// \param q The quaternion to be inverted.
// \return The inverse quaternion.
//
// This function returns the inverse of the given quaternion.
*/
template< typename Type >  // Data type of the quaternion
inline const Quaternion<Type> inv( const Quaternion<Type>& q )
{
   return Quaternion<Type>( q[0], -q[1], -q[2], -q[3] );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Squaring the given quaternion.
// \ingroup quaternion
//
// \param q The quaternion to be squared.
// \return The result of the square operation.
//
// This function squares the given quaternion \a q. This function has the same effect as
// multiplying the quaternion with itself (\f$ q * q \f$).
*/
template< typename Type >  // Data type of the quaternion
inline const Quaternion<Type> sq( const Quaternion<Type>& q )
{
   return q * q;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Swapping the contents of two quaternions.
// \ingroup quaternion
//
// \param a The first quaternion to be swapped.
// \param b The second quaternion to be swapped.
// \return void
// \exception no-throw guarantee.
*/
template< typename Type >  // Data type of the quaternions
inline void swap( Quaternion<Type>& a, Quaternion<Type>& b ) /* throw() */
{
   a.swap( b );
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL ARITHMETIC OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name Quaternion arithmetic operators */
//@{
template< typename T1, typename T2 >
inline const Quaternion< typename MultTrait<T1,T2>::Type >
   operator*( const Quaternion<T1>& lhs, const Quaternion<T2>& rhs );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication operator for the multiplication of two quaternions
//        (\f$ \hat{q}=\hat{p}*\hat{r} \f$).
// \ingroup quaternion
//
// \param lhs The left-hand side quaternion for the multiplication.
// \param rhs The right-hand side quaternion for the multiplication.
// \return The resulting quaternion.
*/
template< typename T1    // Data type of the left-hand side quaternion
        , typename T2 >  // Data type of the right-hand side quaternion
inline const Quaternion< typename MultTrait<T1,T2>::Type >
   operator*( const Quaternion<T1>& lhs, const Quaternion<T2>& rhs )
{
   typedef typename MultTrait<T1,T2>::Type  MT;

   const MT r( lhs[0]*rhs[0] - lhs[1]*rhs[1] - lhs[2]*rhs[2] - lhs[3]*rhs[3] );
   const MT i( lhs[0]*rhs[1] + lhs[1]*rhs[0] + lhs[2]*rhs[3] - lhs[3]*rhs[2] );
   const MT j( lhs[0]*rhs[2] + lhs[2]*rhs[0] + lhs[3]*rhs[1] - lhs[1]*rhs[3] );
   const MT k( lhs[0]*rhs[3] + lhs[3]*rhs[0] + lhs[1]*rhs[2] - lhs[2]*rhs[1] );

   const MT len2( r*r + i*i + j*j + k*k );

   if( std::fabs( len2 - MT(1) ) < MT(1E-8) ) {
      return Quaternion<MT>( r, i, j, k );
   }
   else {
      const MT ilen( MT(1) / std::sqrt( len2 ) );
      return Quaternion<MT>( r*ilen, i*ilen, j*ilen, k*ilen );
   }
}
//*************************************************************************************************




//=================================================================================================
//
//  MULTTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T1, typename T2 >
struct MultTrait< Quaternion<T1>, Quaternion<T2> >
{
   typedef Quaternion< typename MultTrait<T1,T2>::Type >  MultType;
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
template< typename T1, typename T2 >
struct MathTrait< Quaternion<T1>, Quaternion<T2> >
{
   typedef Quaternion< typename MathTrait<T1,T2>::HighType >  HighType;
   typedef Quaternion< typename MathTrait<T1,T2>::LowType  >  LowType;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  TYPE DEFINITIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Quaternion of real type.
// \ingroup quaternion
*/
typedef Quaternion<real>  Quat;
//*************************************************************************************************

} // namespace blaze

#endif
