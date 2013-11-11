//=================================================================================================
/*!
//  \file blaze/math/RotationMatrix.h
//  \brief Implementation of a 3x3 rotation matrix
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

#ifndef _BLAZE_MATH_ROTATIONMATRIX_H_
#define _BLAZE_MATH_ROTATIONMATRIX_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <algorithm>
#include <cmath>
#include <ostream>
#include <limits>
#include <blaze/math/Accuracy.h>
#include <blaze/math/DenseMatrix.h>
#include <blaze/math/Forward.h>
#include <blaze/math/shims/Equal.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/shims/IsNaN.h>
#include <blaze/math/StaticMatrix.h>
#include <blaze/math/StaticVector.h>
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
//  EULER ROTATIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Order of the Euler rotation
// \ingroup math
//
// This codes are needed for the EulerAngles function in order to calculate the Euler angles
// for a specific combination of rotations.
*/
enum EulerRotation {
   XYZs =  0,  //!< Rotation order x, y, z in a static frame.
   ZYXr =  1,  //!< Rotation order z, y, x in a rotating frame.
   XYXs =  2,  //!< Rotation order x, y, x in a static frame.
   XYXr =  3,  //!< Rotation order x, y, z in a rotating frame.
   XZYs =  4,  //!< Rotation order x, z, y in a static frame.
   YZXr =  5,  //!< Rotation order y, z, x in a rotating frame.
   XZXs =  6,  //!< Rotation order x, z, x in a static frame.
   XZXr =  7,  //!< Rotation order x, z, x in a rotating frame.
   YZXs =  8,  //!< Rotation order y, z, x in a static frame.
   XZYr =  9,  //!< Rotation order x, z, y in a rotating frame.
   YZYs = 10,  //!< Rotation order y, z, y in a static frame.
   YZYr = 11,  //!< Rotation order y, z, y in a rotating frame.
   YXZs = 12,  //!< Rotation order y, x, z in a static frame.
   ZXYr = 13,  //!< Rotation order z, x, y in a rotating frame.
   YXYs = 14,  //!< Rotation order y, x, y in a static frame.
   YXYr = 15,  //!< Rotation order y, x, y in a rotating frame.
   ZXYs = 16,  //!< Rotation order z, x, y in a static frame.
   YXZr = 17,  //!< Rotation order y, x, z in a rotating frame.
   ZXZs = 18,  //!< Rotation order z, x, z in a static frame.
   ZXZr = 19,  //!< Rotation order z, x, z in a rotating frame.
   ZYXs = 20,  //!< Rotation order z, y, x in a static frame.
   XYZr = 21,  //!< Rotation order x, y, z in a rotating frame.
   ZYZs = 22,  //!< Rotation order z, y, z in a static frame.
   ZYZr = 23   //!< Rotation order z, y, z in a rotating frame.
};
//*************************************************************************************************




//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\defgroup dense_rotation_matrix RotationMatrix
// \ingroup dense_matrix
*/
/*!\brief Efficient, generic implementation of a 3x3 rotation matrix.
// \ingroup dense_rotation_matrix
//
// The RotationMatrix class is the representation of a 3x3 rotation matrix with a total of 9
// statically allocated elements of arbitrary type. The naming convention of the elements is
// as following:

                          \f[\left(\begin{array}{*{3}{c}}
                          xx & xy & xz \\
                          yx & yy & yz \\
                          zx & zy & zz \\
                          \end{array}\right)\f]\n

// These elements can be accessed directly with the 1D subscript operator or with the 2D function
// operator. The numbering of the matrix elements is

                          \f[\left(\begin{array}{*{3}{c}}
                          0 & 1 & 2 \\
                          3 & 4 & 5 \\
                          6 & 7 & 8 \\
                          \end{array}\right)\f]

// \b Note: The RotationMatrix class can only be instantiated for non-cv-qualified floating point
// types! Therefore the only possible RotationMatrix instantiations are
//
//  - RotationMatrix<float>
//  - RotationMatrix<double>
//  - RotationMatrix<long double>
//
// The attempt to create a rotation matrix with an integral data type results in a compile time
// error.
*/
template< typename Type >  // Data type of the rotation matrix
class RotationMatrix : public DenseMatrix< RotationMatrix<Type>, false >
{
 public:
   //**Type definitions****************************************************************************
   typedef RotationMatrix<Type>   This;           //!< Type of this RotationMatrix instance.
   typedef This                   ResultType;     //!< Result type for expression template evaluations.
   typedef Type                   ElementType;    //!< Type of the matrix elements.
   typedef const RotationMatrix&  CompositeType;  //!< Data type for composite expression templates.
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline RotationMatrix();

   template< typename Axis >
   explicit RotationMatrix( StaticVector<Axis,3UL> axis, Type angle );

   inline RotationMatrix( const RotationMatrix& m );

   template< typename Other >
   inline RotationMatrix( const RotationMatrix<Other>& m );
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   // No explicitly declared destructor.
   //**********************************************************************************************

   //**Operators***********************************************************************************
   /*!\name Operators */
   //@{
                              inline RotationMatrix& operator= ( const RotationMatrix& rhs );
   template< typename Other > inline RotationMatrix& operator= ( const RotationMatrix<Other>& rhs );
                              inline Type            operator[]( size_t index ) const;
                              inline Type            operator()( size_t i, size_t j ) const;
   template< typename Other > inline RotationMatrix& operator*=( const RotationMatrix<Other>& rhs );
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
                              inline size_t          rows() const;
                              inline size_t          columns() const;
                              inline void            reset();
                              inline Type            getDeterminant() const;
                              inline RotationMatrix& transpose();
                              inline RotationMatrix& invert();
                              inline void            swap( RotationMatrix& m ) /* throw() */;
   //@}
   //**********************************************************************************************

   //**Expression template evaluation functions****************************************************
   /*!\name Expression template evaluation functions */
   //@{
   template< typename Other > inline bool isAliased( const Other* alias ) const;
   //@}
   //**********************************************************************************************

   //**Math functions******************************************************************************
   /*!\name Math functions */
   //@{
   template< typename Other >
   inline const StaticMatrix< typename MultTrait<Type,Other>::Type, 3UL, 3UL, false >
      rotate( const StaticMatrix<Other,3UL,3UL,false>& m ) const;

   template< typename Other >
   inline const StaticMatrix< typename MultTrait<Type,Other>::Type, 3UL, 3UL, false >
      diagRotate( const StaticMatrix<Other,3UL,3UL,false>& m ) const;
   //@}
   //**********************************************************************************************

   //**Euler rotations*****************************************************************************
   /*!\name Euler rotations
   //
   // For the classification of the Euler rotation, the following characteristics are
   // defined:\n
   //  - Inner axis: the inner axis is the axis of the first rotation matrix multiplied
   //    to a vector.
   //  - Parity: the parity is even, if the inner axis X is followed by the middle axis
   //    Y, or Y is followed by Z, or Z is followed by X; otherwise parity is odd.
   //  - Repetition: repetition tells, if the first and last axes are the same or different.
   //  - Frame: the frame refers to the frame from which the Euler angles are calculated.
   //
   // Altogether, there are 24 possible Euler rotations. The possibilities are consisting
   // of the choice of the inner axis (X,Y or Z), the parity (even or odd), repetition
   // (yes or no) and the frame (static or rotating). E.g., an Euler order of XYZs stands
   // for the rotation order of x-, y- and z-axis in a static frame (inner axis: X, parity:
   // even, repetition: no, frame: static), whereas YXYr stands for the rotation order y-,
   // x- and y-axis in a rotating frame ( inner axis: Y, parity: odd, repetition: yes,
   // frame: rotating).
   */
   //@{
   inline const StaticVector<Type,3UL> getEulerAnglesXYZ()                   const;
          const StaticVector<Type,3UL> getEulerAngles( EulerRotation order ) const;
   //@}
   //**********************************************************************************************

 private:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline RotationMatrix( Type xx, Type xy, Type xz,
                                   Type yx, Type yy, Type yz,
                                   Type zx, Type zy, Type zz );
   //@}
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   Type v_[9];  //!< The nine statically allocated matrix elements.
                /*!< Access to the matrix elements is gained via the subscript or function call
                     operator. The order of the elements is
                     \f[\left(\begin{array}{*{3}{c}}
                     0 & 1 & 2 \\
                     3 & 4 & 5 \\
                     6 & 7 & 8 \\
                     \end{array}\right)\f] */
   //@}
   //**********************************************************************************************

   //**Friend declarations*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   template< typename Other > friend class Quaternion;

   template< typename Other >
   friend const RotationMatrix<Other> trans( const RotationMatrix<Other>& m );

   template< typename Other >
   friend const RotationMatrix<Other> inv( const RotationMatrix<Other>& m );

   template< typename T1, typename T2 >
   friend const RotationMatrix< typename MultTrait<T1,T2>::Type >
      operator*( const RotationMatrix<T1>& lhs, const RotationMatrix<T2>& rhs );
   /*! \endcond */
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
/*!\brief The default constructor for RotationMatrix.
//
// The diagonal matrix elements are initialized with 1, all other elements are initialized
// with 0.
*/
template< typename Type >  // Data type of the rotation matrix
inline RotationMatrix<Type>::RotationMatrix()
{
   v_[0] = v_[4] = v_[8] = Type(1);
   v_[1] = v_[2] = v_[3] = v_[5] = v_[6] = v_[7] = Type(0);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Rotation matrix constructor.
//
// \param axis The rotation axis.
// \param angle The rotation angle (radian measure).
//
// This constructor creates a rotation matrix from the rotation axis \a axis and the rotation
// angle \a angle. \a axis may be an arbitrary, non-zero vector of any length. However, it is
// allowed to use the zero vector (0,0,0) in combination with an angle of 0. This combination
// results in the default rotation matrix

                          \f[\left(\begin{array}{*{3}{c}}
                          1 & 0 & 0 \\
                          0 & 1 & 0 \\
                          0 & 0 & 1 \\
                          \end{array}\right)\f]
*/
template< typename Type >  // Data type of the rotation matrix
template< typename Axis >  // Data type of the rotation axis
RotationMatrix<Type>::RotationMatrix( StaticVector<Axis,3UL> axis, Type angle )
{
   BLAZE_CONSTRAINT_MUST_BE_FLOATING_POINT_TYPE( Axis );

   BLAZE_USER_ASSERT( ( axis.sqrLength() > Axis(0) || angle == Type(0) ), "Invalid matrix parameters" );

   const Type sina( std::sin(angle) );
   const Type cosa( std::cos(angle) );
   const Type tmp( Type(1)-cosa );

   axis.normalize();

   v_[0] = cosa + axis[0]*axis[0]*tmp;
   v_[1] = axis[0]*axis[1]*tmp - axis[2]*sina;
   v_[2] = axis[0]*axis[2]*tmp + axis[1]*sina;
   v_[3] = axis[1]*axis[0]*tmp + axis[2]*sina;
   v_[4] = cosa + axis[1]*axis[1]*tmp;
   v_[5] = axis[1]*axis[2]*tmp - axis[0]*sina;
   v_[6] = axis[2]*axis[0]*tmp - axis[1]*sina;
   v_[7] = axis[2]*axis[1]*tmp + axis[0]*sina;
   v_[8] = cosa + axis[2]*axis[2]*tmp;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The copy constructor for RotationMatrix.
//
// \param m Rotation matrix to be copied.
//
// The copy constructor is explicitly defined in order to enable/facilitate NRV optimization.
*/
template< typename Type >  // Data type of the rotation matrix
inline RotationMatrix<Type>::RotationMatrix( const RotationMatrix& m )
{
   v_[0] = m.v_[0];
   v_[1] = m.v_[1];
   v_[2] = m.v_[2];
   v_[3] = m.v_[3];
   v_[4] = m.v_[4];
   v_[5] = m.v_[5];
   v_[6] = m.v_[6];
   v_[7] = m.v_[7];
   v_[8] = m.v_[8];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Conversion constructor from different RotationMatrix instances.
//
// \param m Rotation matrix to be copied.
*/
template< typename Type >   // Data type of the rotation matrix
template< typename Other >  // Data type of the foreign rotation matrix
inline RotationMatrix<Type>::RotationMatrix( const RotationMatrix<Other>& m )
{
   v_[0] = m[0];
   v_[1] = m[1];
   v_[2] = m[2];
   v_[3] = m[3];
   v_[4] = m[4];
   v_[5] = m[5];
   v_[6] = m[6];
   v_[7] = m[7];
   v_[8] = m[8];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constructor for a direct initialization of all rotation matrix elements.
//
// \param xx The initial value for the xx-component.
// \param xy The initial value for the xy-component.
// \param xz The initial value for the xz-component.
// \param yx The initial value for the yx-component.
// \param yy The initial value for the yy-component.
// \param yz The initial value for the yz-component.
// \param zx The initial value for the zx-component.
// \param zy The initial value for the zy-component.
// \param zz The initial value for the zz-component.
*/
template< typename Type >  // Data type of the rotation matrix
inline RotationMatrix<Type>::RotationMatrix( Type xx, Type xy, Type xz,
                                             Type yx, Type yy, Type yz,
                                             Type zx, Type zy, Type zz )
{
   v_[0] = xx; v_[1] = xy; v_[2] = xz;
   v_[3] = yx; v_[4] = yy; v_[5] = yz;
   v_[6] = zx; v_[7] = zy; v_[8] = zz;
}
//*************************************************************************************************




//=================================================================================================
//
//  OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Copy assignment operator for RotationMatrix.
//
// \param rhs Rotation matrix to be copied.
// \return Reference to the assigned rotation matrix.
//
// Explicit definition of a copy assignment operator for performance reasons.
*/
template< typename Type >  // Data type of the rotation matrix
inline RotationMatrix<Type>& RotationMatrix<Type>::operator=( const RotationMatrix& rhs )
{
   // This implementation is faster than the synthesized default copy assignment operator and
   // faster than an implementation with the C library function 'memcpy' in combination with a
   // protection against self-assignment. Additionally, this version goes without a protection
   // against self-assignment.
   v_[0] = rhs.v_[0];
   v_[1] = rhs.v_[1];
   v_[2] = rhs.v_[2];
   v_[3] = rhs.v_[3];
   v_[4] = rhs.v_[4];
   v_[5] = rhs.v_[5];
   v_[6] = rhs.v_[6];
   v_[7] = rhs.v_[7];
   v_[8] = rhs.v_[8];
   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Assignment operator for different RotationMatrix instances.
//
// \param rhs Rotation matrix to be copied.
// \return Reference to the assigned rotation matrix.
*/
template< typename Type >   // Data type of the rotation matrix
template< typename Other >  // Data type of the foreign rotation matrix
inline RotationMatrix<Type>& RotationMatrix<Type>::operator=( const RotationMatrix<Other>& rhs )
{
   // This implementation is faster than the synthesized default copy assignment operator and
   // faster than an implementation with the C library function 'memcpy' in combination with a
   // protection against self-assignment. Additionally, this version goes without a protection
   // against self-assignment.
   v_[0] = rhs[0];
   v_[1] = rhs[1];
   v_[2] = rhs[2];
   v_[3] = rhs[3];
   v_[4] = rhs[4];
   v_[5] = rhs[5];
   v_[6] = rhs[6];
   v_[7] = rhs[7];
   v_[8] = rhs[8];
   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief 1D-access to the rotation matrix elements.
//
// \param index Access index. The index has to be in the range \f$[0..8]\f$.
// \return Copy of the accessed element.
//
// In case BLAZE_USER_ASSERT() is active, this operator performs an index check.
*/
template< typename Type >  // Data type of the rotation matrix
inline Type RotationMatrix<Type>::operator[]( size_t index ) const
{
   BLAZE_USER_ASSERT( index < 9, "Invalid rotation matrix access index" );
   return v_[index];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief 2D-access to the rotation matrix elements.
//
// \param i Access index for the row. The index has to be in the range [0..2].
// \param j Access index for the column. The index has to be in the range [0..2].
// \return Copy of the accessed element.
*/
template< typename Type >  // Data type of the rotation matrix
inline Type RotationMatrix<Type>::operator()( size_t i, size_t j ) const
{
   BLAZE_USER_ASSERT( i<3 && j<3, "Invalid rotation matrix access index" );
   return v_[i*3+j];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication assignment operator for the multiplication between two rotation matrices
//        (\f$ A*=B \f$).
//
// \param rhs The right-hand side rotation matrix for the multiplication.
// \return Reference to the rotation matrix.
*/
template< typename Type >   // Data type of the rotation matrix
template< typename Other >  // Data type of the right-hand side rotation matrix
inline RotationMatrix<Type>& RotationMatrix<Type>::operator*=( const RotationMatrix<Other>& rhs )
{
   // Creating a temporary due to data dependencies
   const RotationMatrix tmp( v_[0]*rhs[0] + v_[1]*rhs[3] + v_[2]*rhs[6],
                             v_[0]*rhs[1] + v_[1]*rhs[4] + v_[2]*rhs[7],
                             v_[0]*rhs[2] + v_[1]*rhs[5] + v_[2]*rhs[8],
                             v_[3]*rhs[0] + v_[4]*rhs[3] + v_[5]*rhs[6],
                             v_[3]*rhs[1] + v_[4]*rhs[4] + v_[5]*rhs[7],
                             v_[3]*rhs[2] + v_[4]*rhs[5] + v_[5]*rhs[8],
                             v_[6]*rhs[0] + v_[7]*rhs[3] + v_[8]*rhs[6],
                             v_[6]*rhs[1] + v_[7]*rhs[4] + v_[8]*rhs[7],
                             v_[6]*rhs[2] + v_[7]*rhs[5] + v_[8]*rhs[8] );

   return this->operator=( tmp );
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns the current number of rows of the rotation matrix.
//
// \return The number of rows of the rotation matrix.
*/
template< typename Type >  // Data type of the rotation matrix
inline size_t RotationMatrix<Type>::rows() const
{
   return size_t(3);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the current number of columns of the rotation matrix.
//
// \return The number of columns of the rotation matrix.
*/
template< typename Type >  // Data type of the rotation matrix
inline size_t RotationMatrix<Type>::columns() const
{
   return size_t(3);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Reset to the default initial values.
//
// \return void
//
// This function resets the rotation matrix to the default initial values:

                          \f[\left(\begin{array}{*{3}{c}}
                          1 & 0 & 0 \\
                          0 & 1 & 0 \\
                          0 & 0 & 1 \\
                          \end{array}\right)\f]
*/
template< typename Type >  // Data type of the rotation matrix
inline void RotationMatrix<Type>::reset()
{
   v_[0] = v_[4] = v_[8] = Type(1);
   v_[1] = v_[2] = v_[3] = v_[5] = v_[6] = v_[7] = Type(0);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Calculation of the determinant of the rotation matrix.
//
// \return The determinant of the rotation matrix.
*/
template< typename Type >  // Data type of the rotation matrix
inline Type RotationMatrix<Type>::getDeterminant() const
{
   // Although the determinant of a rotation matrix should always be exactly one, the
   // function calculates the actual determinant to enable checks.
   return v_[0]*v_[4]*v_[8] + v_[1]*v_[5]*v_[6] + v_[2]*v_[3]*v_[7] -
          v_[6]*v_[4]*v_[2] - v_[7]*v_[5]*v_[0] - v_[8]*v_[3]*v_[1];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Transposing the rotation matrix.
//
// \return Reference to the transposed rotation matrix.
//
// This function has the same effect as the invert() function (\f$ R^T = R^-1 \f$).
*/
template< typename Type >  // Data type of the rotation matrix
inline RotationMatrix<Type>& RotationMatrix<Type>::transpose()
{
   std::swap( v_[1], v_[3] );
   std::swap( v_[2], v_[6] );
   std::swap( v_[5], v_[7] );
   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inverting the matrix.
//
// \return Reference to the inverted matrix.
//
// This function has the same effect as the transpose() function (\f$ R^-1 = R^T \f$).
*/
template< typename Type >  // Data type of the rotation matrix
inline RotationMatrix<Type>& RotationMatrix<Type>::invert()
{
   std::swap( v_[1], v_[3] );
   std::swap( v_[2], v_[6] );
   std::swap( v_[5], v_[7] );
   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Swapping the contents of two 3x3 matrices.
//
// \param m The matrix to be swapped.
// \return void
// \exception no-throw guarantee.
*/
template< typename Type >  // Data type of the rotation matrix
inline void RotationMatrix<Type>::swap( RotationMatrix& m ) /* throw() */
{
   std::swap( v_[0], m.v_[0] );
   std::swap( v_[1], m.v_[1] );
   std::swap( v_[2], m.v_[2] );
   std::swap( v_[3], m.v_[3] );
   std::swap( v_[4], m.v_[4] );
   std::swap( v_[5], m.v_[5] );
   std::swap( v_[6], m.v_[6] );
   std::swap( v_[7], m.v_[7] );
   std::swap( v_[8], m.v_[8] );
}
//*************************************************************************************************




//=================================================================================================
//
//  EXPRESSION TEMPLATE EVALUATION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns whether the rotation matrix is aliased with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this matrix, \a false if not.
*/
template< typename Type >   // Data type of the matrix
template< typename Other >  // Data type of the foreign expression
inline bool RotationMatrix<Type>::isAliased( const Other* alias ) const
{
   return static_cast<const void*>( this ) == static_cast<const void*>( alias );
}
//*************************************************************************************************




//=================================================================================================
//
//  MATH FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Rotation of a matrix M (\f$ ROT=R*M*R^{-1} \f$).
//
// \param m The matrix to be rotated.
// \return The rotated matrix.
//
// The function is selected for matrices of different data type (in case \a Type and \a Other
// are supported by the MultTrait class). The function returns a matrix of the higher-order
// data type of the two involved data types.
//
// \b Note: This function is only defined for matrices of floating point type. The attempt to
// use this function with matrices of integral data type will result in a compile time error.
*/
template< typename Type >   // Data type of the rotation matrix
template< typename Other >  // Data type of the standard matrix
inline const StaticMatrix< typename MultTrait<Type,Other>::Type, 3UL, 3UL, false >
   RotationMatrix<Type>::rotate( const StaticMatrix<Other,3UL,3UL,false>& m ) const
{
   BLAZE_CONSTRAINT_MUST_BE_FLOATING_POINT_TYPE( Other );

   typedef StaticMatrix<typename MultTrait<Type,Other>::Type,3UL,3UL,false>  MT;

   //--Multiplication in two steps (number of FLOP = 90, 1 additional temporary matrix)------------

   // Precalculation of tmp = m * R(-1)
   const MT tmp( m.v_[0]*v_[0] + m.v_[1]*v_[1] + m.v_[2]*v_[2],
                 m.v_[0]*v_[3] + m.v_[1]*v_[4] + m.v_[2]*v_[5],
                 m.v_[0]*v_[6] + m.v_[1]*v_[7] + m.v_[2]*v_[8],
                 m.v_[3]*v_[0] + m.v_[4]*v_[1] + m.v_[5]*v_[2],
                 m.v_[3]*v_[3] + m.v_[4]*v_[4] + m.v_[5]*v_[5],
                 m.v_[3]*v_[6] + m.v_[4]*v_[7] + m.v_[5]*v_[8],
                 m.v_[6]*v_[0] + m.v_[7]*v_[1] + m.v_[8]*v_[2],
                 m.v_[6]*v_[3] + m.v_[7]*v_[4] + m.v_[8]*v_[5],
                 m.v_[6]*v_[6] + m.v_[7]*v_[7] + m.v_[8]*v_[8] );

   // Calculating ROT = R * tmp
   return MT( v_[0]*tmp.v_[0] + v_[1]*tmp.v_[3] + v_[2]*tmp.v_[6],
              v_[0]*tmp.v_[1] + v_[1]*tmp.v_[4] + v_[2]*tmp.v_[7],
              v_[0]*tmp.v_[2] + v_[1]*tmp.v_[5] + v_[2]*tmp.v_[8],
              v_[3]*tmp.v_[0] + v_[4]*tmp.v_[3] + v_[5]*tmp.v_[6],
              v_[3]*tmp.v_[1] + v_[4]*tmp.v_[4] + v_[5]*tmp.v_[7],
              v_[3]*tmp.v_[2] + v_[4]*tmp.v_[5] + v_[5]*tmp.v_[8],
              v_[6]*tmp.v_[0] + v_[7]*tmp.v_[3] + v_[8]*tmp.v_[6],
              v_[6]*tmp.v_[1] + v_[7]*tmp.v_[4] + v_[8]*tmp.v_[7],
              v_[6]*tmp.v_[2] + v_[7]*tmp.v_[5] + v_[8]*tmp.v_[8] );

   //--Multiplication in one step (number of FLOP = 180, no additional temporary matrix)-----------
   /*
   return MT( m.v_[0]*v_[0]*v_[0] + m.v_[4]*v_[1]*v_[1] + m.v_[8]*v_[2]*v_[2] + v_[0]*v_[1]*( m.v_[1]+m.v_[3] ) + v_[0]*v_[2]*( m.v_[2]+m.v_[6] ) + v_[1]*v_[2]*( m.v_[5]+m.v_[7] ),
              v_[0]*( m.v_[0]*v_[3] + m.v_[1]*v_[4] + m.v_[2]*v_[5] ) + v_[1]*( m.v_[3]*v_[3] + m.v_[4]*v_[4] + m.v_[5]*v_[5] ) + v_[2]*( m.v_[6]*v_[3] + m.v_[7]*v_[4] + m.v_[8]*v_[5] ),
              v_[0]*( m.v_[0]*v_[6] + m.v_[1]*v_[7] + m.v_[2]*v_[8] ) + v_[1]*( m.v_[3]*v_[6] + m.v_[4]*v_[7] + m.v_[5]*v_[8] ) + v_[2]*( m.v_[6]*v_[6] + m.v_[7]*v_[7] + m.v_[8]*v_[8] ),
              v_[3]*( m.v_[0]*v_[0] + m.v_[1]*v_[1] + m.v_[2]*v_[2] ) + v_[4]*( m.v_[3]*v_[0] + m.v_[4]*v_[1] + m.v_[5]*v_[2] ) + v_[5]*( m.v_[6]*v_[0] + m.v_[7]*v_[1] + m.v_[8]*v_[2] ),
              m.v_[0]*v_[3]*v_[3] + m.v_[4]*v_[4]*v_[4] + m.v_[8]*v_[5]*v_[5] + v_[3]*v_[4]*( m.v_[1]+m.v_[3] ) + v_[3]*v_[5]*( m.v_[2]+m.v_[6] ) + v_[4]*v_[5]*( m.v_[5]+m.v_[7] ),
              v_[3]*( m.v_[0]*v_[6] + m.v_[1]*v_[7] + m.v_[2]*v_[8] ) + v_[4]*( m.v_[3]*v_[6] + m.v_[4]*v_[7] + m.v_[5]*v_[8] ) + v_[5]*( m.v_[6]*v_[6] + m.v_[7]*v_[7] + m.v_[8]*v_[8] ),
              v_[6]*( m.v_[0]*v_[0] + m.v_[1]*v_[1] + m.v_[2]*v_[2] ) + v_[7]*( m.v_[3]*v_[0] + m.v_[4]*v_[1] + m.v_[5]*v_[2] ) + v_[8]*( m.v_[6]*v_[0] + m.v_[7]*v_[1] + m.v_[8]*v_[2] ),
              v_[6]*( m.v_[0]*v_[3] + m.v_[1]*v_[4] + m.v_[2]*v_[5] ) + v_[7]*( m.v_[3]*v_[3] + m.v_[4]*v_[4] + m.v_[5]*v_[5] ) + v_[8]*( m.v_[6]*v_[3] + m.v_[7]*v_[4] + m.v_[8]*v_[5] ),
              m.v_[0]*v_[6]*v_[6] + m.v_[4]*v_[7]*v_[7] + m.v_[8]*v_[8]*v_[8] + v_[6]*v_[7]*( m.v_[1]+m.v_[3] ) + v_[6]*v_[8]*( m.v_[2]+m.v_[6] ) + v_[7]*v_[8]*( m.v_[5]+m.v_[7] ) );
   */
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Rotation of a diagonal matrix M (\f$ ROT=R*M*R^{-1} \f$).
//
// \param m The diagonal matrix to be rotated.
// \return The rotated matrix.
//
// The DiagRotate function is a special case of the rotate function. The matrix is assumed to
// be a diagonal matrix, which reduces the number of floating point operations of the rotation.
// The function is selected for matrices of different data type (in case \a Type and \a Other
// are supported by the MultTrait class). The function returns a matrix of the higher-order
// data type of the two involved data types.
//
// \b Note: This function is only defined for matrices of floating point type. The attempt to
// use this function with matrices of integral data type will result in a compile time error.
*/
template< typename Type >   // Data type of the rotation matrix
template< typename Other >  // Data type of the diagonal standard matrix
inline const StaticMatrix< typename MultTrait<Type,Other>::Type, 3UL, 3UL, false >
   RotationMatrix<Type>::diagRotate( const StaticMatrix<Other,3UL,3UL,false>& m ) const
{
   BLAZE_CONSTRAINT_MUST_BE_FLOATING_POINT_TYPE( Other );

   typedef StaticMatrix<typename MultTrait<Type,Other>::Type,3UL,3UL,false>  MT;

   // Precalculating tmp = m * R(-1)
   const MT tmp( m.v_[0]*v_[0], m.v_[0]*v_[3], m.v_[0]*v_[6],
                 m.v_[4]*v_[1], m.v_[4]*v_[4], m.v_[4]*v_[7],
                 m.v_[8]*v_[2], m.v_[8]*v_[5], m.v_[8]*v_[8] );

   // Calculating ROT = R * tmp
   return MT( v_[0]*tmp.v_[0] + v_[1]*tmp.v_[3] + v_[2]*tmp.v_[6],
              v_[0]*tmp.v_[1] + v_[1]*tmp.v_[4] + v_[2]*tmp.v_[7],
              v_[0]*tmp.v_[2] + v_[1]*tmp.v_[5] + v_[2]*tmp.v_[8],
              v_[3]*tmp.v_[0] + v_[4]*tmp.v_[3] + v_[5]*tmp.v_[6],
              v_[3]*tmp.v_[1] + v_[4]*tmp.v_[4] + v_[5]*tmp.v_[7],
              v_[3]*tmp.v_[2] + v_[4]*tmp.v_[5] + v_[5]*tmp.v_[8],
              v_[6]*tmp.v_[0] + v_[7]*tmp.v_[3] + v_[8]*tmp.v_[6],
              v_[6]*tmp.v_[1] + v_[7]*tmp.v_[4] + v_[8]*tmp.v_[7],
              v_[6]*tmp.v_[2] + v_[7]*tmp.v_[5] + v_[8]*tmp.v_[8] );
}
//*************************************************************************************************




//=================================================================================================
//
//  EULER ROTATIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Calculation of the Euler angles (in radian measure).
//
// \return The Euler angles for a rotation order of x, y, z (radian measure).
//
// The Euler angles are calculated for a rotation order of x-, y- and z-axis.
*/
template< typename Type >  // Data type of the rotation matrix
inline const StaticVector<Type,3UL> RotationMatrix<Type>::getEulerAnglesXYZ() const
{
   const Type cy( std::sqrt( v_[0]*v_[0] + v_[3]*v_[3] ) );

   if( cy > accuracy ) {
      return StaticVector<Type,3UL>( std::atan2( v_[7], v_[8] ), std::atan2( -v_[6], cy ), std::atan2( v_[3], v_[0] ) );
   }
   else {
      return StaticVector<Type,3UL>( std::atan2( -v_[5], v_[4] ), std::atan2( -v_[6], cy ), Type(0) );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Calculation of the Euler angles for a specific rotation order.
//
// \param order The specific rotation order.
// \return The specific Euler angles (radian measure).
*/
template< typename Type >  // Data type of the rotation matrix
const StaticVector<Type,3UL> RotationMatrix<Type>::getEulerAngles( EulerRotation order ) const
{
   static const unsigned int eulSafe[4] = { 0, 1, 2, 0 };
   static const unsigned int eulNext[4] = { 1, 2, 0, 1 };

   StaticVector<Type,3UL> ea;

   // Unpacking the euler order
   const unsigned int frame( order&1 );
   const unsigned int repetition( (order&2)>>1 );
   const unsigned int parity( (order&4)>>2 );
   const unsigned int i( eulSafe[(order&24)>>3] );
   const unsigned int j( eulNext[i+parity] );
   const unsigned int k( eulNext[i+1-parity] );

   // Treatment of rotations with repetition
   if( repetition ) {
      const Type sy( std::sqrt( v_[i*3+j]*v_[i*3+j] + v_[i*3+k]*v_[i*3+k] ) );
      if( sy > accuracy ) {
         ea[0] = std::atan2( v_[i*3+j], v_[i*3+k] );
         ea[1] = std::atan2( sy, v_[i*3+i] );
         ea[2] = std::atan2( v_[j*3+i], -v_[k*3+i] );
      }
      else {
         ea[0] = std::atan2( -v_[j*3+k], v_[j*3+j] );
         ea[1] = std::atan2( sy, v_[i*3+i] );
         ea[2] = Type(0);
      }
   }

   // Treatment of rotations without repetition
   else {
      const Type cy( std::sqrt( v_[i*3+i]*v_[i*3+i] + v_[j*3+i]*v_[j*3+i] ) );
      if( cy > accuracy ) {
         ea[0] = std::atan2( v_[k*3+j], v_[k*3+k] );
         ea[1] = std::atan2( -v_[k*3+i], cy );
         ea[2] = std::atan2( v_[j*3+i], v_[i*3+i] );
      }
      else {
         ea[0] = std::atan2( -v_[j*3+k], v_[j*3+j] );
         ea[1] = std::atan2( -v_[k*3+i], cy );
         ea[2] = Type(0);
      }
   }

   // Treatment of an odd partity
   if( parity ) {
      ea[0] = -ea[0];
      ea[1] = -ea[1];
      ea[2] = -ea[2];
   }

   // Treatment of a rotating frame
   if( frame ) {
      Type tmp = ea[0];
      ea[0] = ea[2];
      ea[2] = tmp;
   }

   return ea;
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name RotationMatrix operators */
//@{
template< typename T1, typename T2 >
inline bool operator==( const RotationMatrix<T1>& lhs, const RotationMatrix<T2>& rhs );

template< typename T1, typename T2 >
inline bool operator!=( const RotationMatrix<T1>& lhs, const RotationMatrix<T2>& rhs );

template< typename Type >
std::ostream& operator<<( std::ostream& os, const RotationMatrix<Type>& m );

template< typename Type >
inline bool isnan( const RotationMatrix<Type>& m );

template< typename Type >
inline const StaticMatrix<Type,3UL,3UL,false> abs( const RotationMatrix<Type>& m );

template< typename Type >
inline const StaticMatrix<Type,3UL,3UL,false> fabs( const RotationMatrix<Type>& m );

template< typename Type >
inline void reset( RotationMatrix<Type>& m );

template< typename Type >
inline void clear( RotationMatrix<Type>& m );

template< typename Type >
inline bool isDefault( const RotationMatrix<Type>& m );

template< typename Type >
inline const RotationMatrix<Type> trans( const RotationMatrix<Type>& m );

template< typename Type >
inline const RotationMatrix<Type> inv( const RotationMatrix<Type>& m );

template< typename Type >
inline const RotationMatrix<Type> sq( const RotationMatrix<Type>& m );

template< typename Type >
inline void swap( RotationMatrix<Type>& a, RotationMatrix<Type>& b ) /* throw() */;
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality operator for the comparison of two rotation matrices.
// \ingroup dense_rotation_matrix
//
// \param lhs The left-hand side rotation matrix for the comparison.
// \param rhs The right-hand side rotation matrix for the comparison.
// \return \a true if the two rotation matrices are equal, \a false if not.
*/
template< typename T1    // Data type of the left-hand side rotation matrix
        , typename T2 >  // Data type of the right-hand side rotation matrix
inline bool operator==( const RotationMatrix<T1>& lhs, const RotationMatrix<T2>& rhs )
{
   // In order to compare the two matrices, the data values of the lower-order data
   // type are converted to the higher-order data type within the equal function.
   if( !equal( lhs[0], rhs[0] ) ||
       !equal( lhs[1], rhs[1] ) ||
       !equal( lhs[2], rhs[2] ) ||
       !equal( lhs[3], rhs[3] ) ||
       !equal( lhs[4], rhs[4] ) ||
       !equal( lhs[5], rhs[5] ) ||
       !equal( lhs[6], rhs[6] ) ||
       !equal( lhs[7], rhs[7] ) ||
       !equal( lhs[8], rhs[8] ) )
      return false;
   else return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inequality operator for the comparison of two rotation matrices.
// \ingroup dense_rotation_matrix
//
// \param lhs The left-hand side rotation matrix for the comparison.
// \param rhs The right-hand side rotation matrix for the comparison.
// \return \a true if the two rotation matrices are not equal, \a false if they are equal.
*/
template< typename T1    // Data type of the left-hand side rotation matrix
        , typename T2 >  // Data type of the right-hand side rotation matrix
inline bool operator!=( const RotationMatrix<T1>& lhs, const RotationMatrix<T2>& rhs )
{
   return !( lhs == rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Global output operator for 3x3 rotation matrices.
// \ingroup dense_rotation_matrix
//
// \param os Reference to the output stream.
// \param m Reference to a constant rotation matrix object.
// \return Reference to the output stream.
*/
template< typename Type >  // Data type of the rotation matrix
std::ostream& operator<<( std::ostream& os, const RotationMatrix<Type>& m )
{
   return os << " ( " << m[0] << " , " << m[1] << " , " << m[2] << " )\n"
             << " ( " << m[3] << " , " << m[4] << " , " << m[5] << " )\n"
             << " ( " << m[6] << " , " << m[7] << " , " << m[8] << " )\n";
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks the given rotation matrix for not-a-number elements.
// \ingroup dense_rotation_matrix
//
// \param m The rotation matrix to be checked for not-a-number elements.
// \return \a true if at least one element of the matrix is not-a-number, \a false otherwise.
*/
template< typename Type >  // Data type of the rotation matrix
inline bool isnan( const RotationMatrix<Type>& m )
{
   if( isnan( m[0] ) || isnan( m[1] ) || isnan( m[2] ) ||
       isnan( m[3] ) || isnan( m[4] ) || isnan( m[5] ) ||
       isnan( m[6] ) || isnan( m[7] ) || isnan( m[8] ) )
      return true;
   else return false;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns a matrix containing the absolute values of each single element of \a m.
// \ingroup dense_rotation_matrix
//
// \param m The input rotation matrix.
// \return The absolute value of each single element of \a m.
//
// The \a abs function calculates the absolute value of each element of the input rotation
// matrix \a m.
*/
template< typename Type >  // Data type of the rotation matrix
inline const StaticMatrix<Type,3UL,3UL,false> abs( const RotationMatrix<Type>& m )
{
   using std::abs;
   return StaticMatrix<Type,3UL,3UL,false>( abs(m[0]), abs(m[1]), abs(m[2]),
                                            abs(m[3]), abs(m[4]), abs(m[5]),
                                            abs(m[6]), abs(m[7]), abs(m[8]) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns a matrix containing the absolute values of each single element of \a m.
// \ingroup dense_rotation_matrix
//
// \param m The input rotation matrix.
// \return The absolute value of each single element of \a m.
//
// The \a fabs function calculates the absolute value of each element of the input rotation
// matrix \a m.
*/
template< typename Type >  // Data type of the rotation matrix
inline const StaticMatrix<Type,3UL,3UL,false> fabs( const RotationMatrix<Type>& m )
{
   using std::fabs;
   return StaticMatrix<Type,3UL,3UL,false>( fabs(m[0]), fabs(m[1]), fabs(m[2]),
                                            fabs(m[3]), fabs(m[4]), fabs(m[5]),
                                            fabs(m[6]), fabs(m[7]), fabs(m[8]) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Resetting the given rotation matrix.
// \ingroup dense_rotation_matrix
//
// \param m The rotation matrix to be resetted.
// \return void
*/
template< typename Type >  // Data type of the rotation matrix
inline void reset( RotationMatrix<Type>& m )
{
   m.reset();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Clearing the given rotation matrix.
// \ingroup dense_rotation_matrix
//
// \param m The rotation matrix to be cleared.
// \return void
//
// Clearing a rotation matrix is equivalent to resetting it via the reset() function.
*/
template< typename Type >  // Data type of the rotation matrix
inline void clear( RotationMatrix<Type>& m )
{
   m.reset();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the given rotation matrix is in default state.
// \ingroup dense_rotation_matrix
//
// \param m The rotation matrix to be tested for its default state.
// \return \a true in case the given matrix is component-wise zero, \a false otherwise.
*/
template< typename Type >  // Data type of the rotation matrix
inline bool isDefault( const RotationMatrix<Type>& m )
{
   return ( m[0] == Type(1) ) && ( m[1] == Type(0) ) && ( m[2] == Type(0) ) &&
          ( m[3] == Type(0) ) && ( m[4] == Type(1) ) && ( m[5] == Type(0) ) &&
          ( m[6] == Type(0) ) && ( m[7] == Type(0) ) && ( m[8] == Type(1) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Calculation of the transpose of the rotation matrix.
// \ingroup dense_rotation_matrix
//
// \param m The rotation matrix to be transposed.
// \return The transpose of the rotation matrix.
//
// This function returns the transpose of the given rotation matrix:

   \code
   blaze::Rot3 R1, R2;
   // ... Resizing and initialization
   R1 = trans( R2 );
   \endcode

// Note that this function has the same effect as the inv() function (\f$ R^T = R^-1 \f$).
*/
template< typename Type >  // Data type of the rotation matrix
inline const RotationMatrix<Type> trans( const RotationMatrix<Type>& m )
{
   return RotationMatrix<Type>( m[0], m[3], m[6],
                                m[1], m[4], m[7],
                                m[2], m[5], m[8] );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inverting the given rotation matrix.
// \ingroup dense_rotation_matrix
//
// \param m The rotation matrix to be inverted.
// \return The inverse rotation matrix.
//
// This function returns the inverse of the given rotation matrix:

   \code
   blaze::Rot3 R1, R2;
   // ... Resizing and initialization
   R1 = inv( R2 );
   \endcode

// Note that this function has the same effect as the trans() function (\f$ R^-1 = R^T \f$).
*/
template< typename Type >  // Data type of the rotation matrix
inline const RotationMatrix<Type> inv( const RotationMatrix<Type>& m )
{
   return RotationMatrix<Type>( m[0], m[3], m[6], m[1], m[4], m[7], m[2], m[5], m[8] );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Squaring the given rotation matrix.
// \ingroup dense_rotation_matrix
//
// \param m The rotation matrix to be squared.
// \return The result of the square operation.
//
// This function squares the given rotation matrix \a m. This function has the same effect as
// multiplying the rotation matrix with itself (\f$ m * m \f$).
*/
template< typename Type >  // Data type of the rotation matrix
inline const RotationMatrix<Type> sq( const RotationMatrix<Type>& m )
{
   return m * m;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Swapping the contents of two rotation matrices.
// \ingroup dense_rotation_matrix
//
// \param a The first rotation matrix to be swapped.
// \param b The second rotation matrix to be swapped.
// \return void
// \exception no-throw guarantee.
*/
template< typename Type >  // Data type of the rotation matrices
inline void swap( RotationMatrix<Type>& a, RotationMatrix<Type>& b ) /* throw() */
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
/*!\name RotationMatrix arithmetic operators */
//@{
template< typename T1, typename T2 >
inline const StaticVector< typename MultTrait<T1,T2>::Type, 3UL, false >
   operator*( const RotationMatrix<T1>& lhs, const StaticVector<T2,3UL,false>& rhs );

template< typename T1, typename T2 >
inline const StaticVector< typename MultTrait<T1,T2>::Type, 3UL, true >
   operator*( const StaticVector<T1,3UL,true>& lhs, const RotationMatrix<T2>& rhs );

template< typename T1, typename T2 >
inline const StaticMatrix< typename MultTrait<T1,T2>::Type, 3UL, 3UL, false >
   operator*( const RotationMatrix<T1>& lhs, const StaticMatrix<T2,3UL,3UL,false>& rhs );

template< typename T1, typename T2 >
inline const StaticMatrix< typename MultTrait<T1,T2>::Type, 3UL, 3UL, false >
   operator*( const StaticMatrix<T1,3UL,3UL,false>& lhs, const RotationMatrix<T2>& rhs );

template< typename T1, typename T2 >
inline const RotationMatrix< typename MultTrait<T1,T2>::Type >
   operator*( const RotationMatrix<T1>& lhs, const RotationMatrix<T2>& rhs );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication operator for the multiplication of a rotation matrix and a vector
//        (\f$ \vec{a}=B*\vec{c} \f$).
// \ingroup dense_rotation_matrix
//
// \param lhs The left-hand side rotation matrix for the multiplication.
// \param rhs The right-hand side vector for the multiplication.
// \return The resulting vector.
//
// This operator is selected for multiplications between rotation matrices and vectors of two
// different data types \a T1 and \a T2, which are supported by the MultTrait class. The operator
// returns a vector of the higher-order data type of the two involved data types.
*/
template< typename T1    // Data type of the left-hand side rotation matrix
        , typename T2 >  // Data type of the right-hand side vector
inline const StaticVector< typename MultTrait<T1,T2>::Type, 3UL, false >
   operator*( const RotationMatrix<T1>& lhs, const StaticVector<T2,3UL,false>& rhs )
{
   typedef typename MultTrait<T1,T2>::Type  MT;
   return StaticVector<MT,3UL,false>( lhs[0]*rhs[0] + lhs[1]*rhs[1] + lhs[2]*rhs[2],
                                      lhs[3]*rhs[0] + lhs[4]*rhs[1] + lhs[5]*rhs[2],
                                      lhs[6]*rhs[0] + lhs[7]*rhs[1] + lhs[8]*rhs[2] );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication operator for the multiplication of a vector and a rotation matrix
//        (\f$ \vec{a}=\vec{b}^T*B \f$).
// \ingroup dense_rotation_matrix
//
// \param lhs The left-hand side transpose vector for the multiplication.
// \param rhs The right-hand side rotation matrix for the multiplication.
// \return The resulting vector.
//
// This operator is selected for multiplications between rotation matrices and vectors of two
// different data types \a T1 and \a T2, which are supported by the MultTrait class. The operator
// returns a vector of the higher-order data type of the two involved data types.
*/
template< typename T1    // Data type of the left-hand side vector
        , typename T2 >  // Data type of the right-hand side rotation matrix
inline const StaticVector< typename MultTrait<T1,T2>::Type, 3UL, true >
   operator*( const StaticVector<T1,3UL,true>& lhs, const RotationMatrix<T2>& rhs )
{
   typedef typename MultTrait<T1,T2>::Type  MT;
   return StaticVector<MT,3UL,true>( lhs[0]*rhs[0] + lhs[1]*rhs[3] + lhs[2]*rhs[6],
                                     lhs[0]*rhs[1] + lhs[1]*rhs[4] + lhs[2]*rhs[7],
                                     lhs[0]*rhs[2] + lhs[1]*rhs[5] + lhs[2]*rhs[8] );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication operator for the multiplication of a rotation matrix and a standard
//        matrix (\f$ A=R*B \f$).
// \ingroup dense_rotation_matrix
//
// \param lhs The left-hand side rotation matrix for the multiplication.
// \param rhs The right-hand side standard matrix for the multiplication.
// \return The resulting matrix.
//
// This operator is selected for multiplications between matrices of two different data types
// \a T1 and \a T2, which are supported by the MultTrait class. The operator returns a matrix
// of the higher-order data type of the two involved matrix data types.
*/
template< typename T1    // Data type of the left-hand side rotation matrix
        , typename T2 >  // Data type of the right-hand side standard matrix
inline const StaticMatrix< typename MultTrait<T1,T2>::Type, 3UL, 3UL, false >
   operator*( const RotationMatrix<T1>& lhs, const StaticMatrix<T2,3UL,3UL,false>& rhs )
{
   typedef StaticMatrix<typename MultTrait<T1,T2>::Type,3UL,3UL,false>  MT;
   return MT( lhs[0]*rhs[0] + lhs[1]*rhs[3] + lhs[2]*rhs[6],
              lhs[0]*rhs[1] + lhs[1]*rhs[4] + lhs[2]*rhs[7],
              lhs[0]*rhs[2] + lhs[1]*rhs[5] + lhs[2]*rhs[8],
              lhs[3]*rhs[0] + lhs[4]*rhs[3] + lhs[5]*rhs[6],
              lhs[3]*rhs[1] + lhs[4]*rhs[4] + lhs[5]*rhs[7],
              lhs[3]*rhs[2] + lhs[4]*rhs[5] + lhs[5]*rhs[8],
              lhs[6]*rhs[0] + lhs[7]*rhs[3] + lhs[8]*rhs[6],
              lhs[6]*rhs[1] + lhs[7]*rhs[4] + lhs[8]*rhs[7],
              lhs[6]*rhs[2] + lhs[7]*rhs[5] + lhs[8]*rhs[8] );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication operator for the multiplication of a standard matrix and a rotation
//        matrix (\f$ A=B*R \f$).
// \ingroup dense_rotation_matrix
//
// \param lhs The left-hand side standard matrix for the multiplication.
// \param rhs The right-hand side rotation matrix for the multiplication.
// \return The resulting matrix.
//
// This operator is selected for multiplications between matrices of two different data types
// \a T1 and \a T2, which are supported by the MultTrait class. The operator returns a matrix
// of the higher-order data type of the two involved matrix data types.
*/
template< typename T1    // Data type of the left-hand side standard matrix
        , typename T2 >  // Data type of the right-hand side rotation matrix
inline const StaticMatrix< typename MultTrait<T1,T2>::Type, 3UL, 3UL, false >
   operator*( const StaticMatrix<T1,3UL,3UL,false>& lhs, const RotationMatrix<T2>& rhs )
{
   typedef StaticMatrix<typename MultTrait<T1,T2>::Type,3UL,3UL,false>  MT;
   return MT( lhs[0]*rhs[0] + lhs[1]*rhs[3] + lhs[2]*rhs[6],
              lhs[0]*rhs[1] + lhs[1]*rhs[4] + lhs[2]*rhs[7],
              lhs[0]*rhs[2] + lhs[1]*rhs[5] + lhs[2]*rhs[8],
              lhs[3]*rhs[0] + lhs[4]*rhs[3] + lhs[5]*rhs[6],
              lhs[3]*rhs[1] + lhs[4]*rhs[4] + lhs[5]*rhs[7],
              lhs[3]*rhs[2] + lhs[4]*rhs[5] + lhs[5]*rhs[8],
              lhs[6]*rhs[0] + lhs[7]*rhs[3] + lhs[8]*rhs[6],
              lhs[6]*rhs[1] + lhs[7]*rhs[4] + lhs[8]*rhs[7],
              lhs[6]*rhs[2] + lhs[7]*rhs[5] + lhs[8]*rhs[8] );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication operator for the multiplication of two rotation matrices (\f$ A=B*C \f$).
// \ingroup dense_rotation_matrix
//
// \param lhs The left-hand side rotation matrix for the multiplication.
// \param rhs The right-hand side rotation matrix for the multiplication.
// \return The resulting rotation matrix.
//
// This operator is selected for multiplications between rotation matrices of two different
// data types \a T1 and \a T2, which are supported by the MultTrait class. The operator
// returns a matrix of the higher-order data type of the two involved matrix data types.
*/
template< typename T1    // Data type of the left-hand side rotation matrix
        , typename T2 >  // Data type of the right-hand side rotation matrix
inline const RotationMatrix< typename MultTrait<T1,T2>::Type >
   operator*( const RotationMatrix<T1>& lhs, const RotationMatrix<T2>& rhs )
{
   typedef typename MultTrait<T1,T2>::Type  MT;
   return RotationMatrix<MT>( lhs.v_[0]*rhs.v_[0] + lhs.v_[1]*rhs.v_[3] + lhs.v_[2]*rhs.v_[6],
                              lhs.v_[0]*rhs.v_[1] + lhs.v_[1]*rhs.v_[4] + lhs.v_[2]*rhs.v_[7],
                              lhs.v_[0]*rhs.v_[2] + lhs.v_[1]*rhs.v_[5] + lhs.v_[2]*rhs.v_[8],
                              lhs.v_[3]*rhs.v_[0] + lhs.v_[4]*rhs.v_[3] + lhs.v_[5]*rhs.v_[6],
                              lhs.v_[3]*rhs.v_[1] + lhs.v_[4]*rhs.v_[4] + lhs.v_[5]*rhs.v_[7],
                              lhs.v_[3]*rhs.v_[2] + lhs.v_[4]*rhs.v_[5] + lhs.v_[5]*rhs.v_[8],
                              lhs.v_[6]*rhs.v_[0] + lhs.v_[7]*rhs.v_[3] + lhs.v_[8]*rhs.v_[6],
                              lhs.v_[6]*rhs.v_[1] + lhs.v_[7]*rhs.v_[4] + lhs.v_[8]*rhs.v_[7],
                              lhs.v_[6]*rhs.v_[2] + lhs.v_[7]*rhs.v_[5] + lhs.v_[8]*rhs.v_[8] );
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
struct MultTrait< RotationMatrix<T1>, StaticVector<T2,3UL,false> >
{
   typedef StaticVector< typename MultTrait<T1,T2>::Type, 3UL, false >  Type;
};

template< typename T1, typename T2 >
struct MultTrait< StaticVector<T1,3UL,true>, RotationMatrix<T2> >
{
   typedef StaticVector< typename MultTrait<T1,T2>::Type, 3UL, true >  Type;
};

template< typename T1, typename T2 >
struct MultTrait< RotationMatrix<T1>, DynamicVector<T2,false> >
{
   typedef StaticVector< typename MultTrait<T1,T2>::Type, 3UL, false >  Type;
};

template< typename T1, typename T2 >
struct MultTrait< DynamicVector<T1,true>, RotationMatrix<T2> >
{
   typedef StaticVector< typename MultTrait<T1,T2>::Type, 3UL, true >  Type;
};

template< typename T1, typename T2 >
struct MultTrait< RotationMatrix<T1>, CompressedVector<T2,false> >
{
   typedef StaticVector< typename MultTrait<T1,T2>::Type, 3UL, false >  Type;
};

template< typename T1, typename T2 >
struct MultTrait< CompressedVector<T1,true>, RotationMatrix<T2> >
{
   typedef StaticVector< typename MultTrait<T1,T2>::Type, 3UL, true >  Type;
};

template< typename T1, typename T2 >
struct MultTrait< RotationMatrix<T1>, StaticMatrix<T2,3UL,3UL,false> >
{
   typedef StaticMatrix< typename MultTrait<T1,T2>::Type, 3UL, 3UL, false >  Type;
};

template< typename T1, typename T2 >
struct MultTrait< StaticMatrix<T1,3UL,3UL,false>, RotationMatrix<T2> >
{
   typedef StaticMatrix< typename MultTrait<T1,T2>::Type, 3UL, 3UL, false >  Type;
};

template< typename T1, typename T2, bool SO >
struct MultTrait< RotationMatrix<T1>, DynamicMatrix<T2,SO> >
{
   typedef DynamicMatrix< typename MultTrait<T1,T2>::Type, false >  Type;
};

template< typename T1, bool SO, typename T2 >
struct MultTrait< DynamicMatrix<T1,SO>, RotationMatrix<T2> >
{
   typedef DynamicMatrix< typename MultTrait<T1,T2>::Type, false >  Type;
};

template< typename T1, typename T2, bool SO >
struct MultTrait< RotationMatrix<T1>, CompressedMatrix<T2,SO> >
{
   typedef DynamicMatrix< typename MultTrait<T1,T2>::Type, false >  Type;
};

template< typename T1, bool SO, typename T2 >
struct MultTrait< CompressedMatrix<T1,SO>, RotationMatrix<T2> >
{
   typedef DynamicMatrix< typename MultTrait<T1,T2>::Type, false >  Type;
};

template< typename T1, typename T2 >
struct MultTrait< RotationMatrix<T1>, RotationMatrix<T2> >
{
   typedef RotationMatrix< typename MultTrait<T1,T2>::Type >  Type;
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
struct MathTrait< RotationMatrix<T1>, RotationMatrix<T2> >
{
   typedef RotationMatrix< typename MathTrait<T1,T2>::HighType >  HighType;
   typedef RotationMatrix< typename MathTrait<T1,T2>::LowType  >  LowType;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  TYPE DEFINITIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Rotation matrix of real type.
// \ingroup dense_rotation_matrix
*/
typedef RotationMatrix<real>  Rot3;
//*************************************************************************************************

} // namespace blaze

#endif
