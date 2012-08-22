//=================================================================================================
/*!
//  \file blaze/math/Matrix.h
//  \brief Header file for all basic Matrix functionality
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

#ifndef _BLAZE_MATH_MATRIX_H_
#define _BLAZE_MATH_MATRIX_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/expressions/Matrix.h>
#include <blaze/util/Assert.h>


namespace blaze {

//=================================================================================================
//
//  GLOBAL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\name Matrix global functions */
//@{
template< typename MT, bool SO >
inline size_t rows( const Matrix<MT,SO>& m );

template< typename MT, bool SO >
inline size_t columns( const Matrix<MT,SO>& m );

template< typename MT1, bool SO1, typename MT2, bool SO2 >
inline void assign( Matrix<MT1,SO1>& lhs, const Matrix<MT2,SO2>& rhs );

template< typename MT1, bool SO1, typename MT2, bool SO2 >
inline void addAssign( Matrix<MT1,SO1>& lhs, const Matrix<MT2,SO2>& rhs );

template< typename MT1, bool SO1, typename MT2, bool SO2 >
inline void subAssign( Matrix<MT1,SO1>& lhs, const Matrix<MT2,SO2>& rhs );

template< typename MT1, bool SO1, typename MT2, bool SO2 >
inline void multAssign( Matrix<MT1,SO1>& lhs, const Matrix<MT2,SO2>& rhs );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the current number of rows of the matrix.
// \ingroup matrix
//
// \param m The given matrix.
// \return The number of rows of the matrix.
*/
template< typename MT  // Type of the matrix
        , bool SO >    // Storage order of the matrix
inline size_t rows( const Matrix<MT,SO>& m )
{
   return (~m).rows();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the current number of columns of the matrix.
// \ingroup matrix
//
// \param m The given matrix.
// \return The number of columns of the matrix.
*/
template< typename MT  // Type of the matrix
        , bool SO >    // Storage order of the matrix
inline size_t columns( const Matrix<MT,SO>& m )
{
   return (~m).columns();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the assignment of a matrix to a matrix.
// \ingroup matrix
//
// \param lhs The target left-hand side matrix.
// \param rhs The right-hand side matrix to be assigned.
// \return void
//
// This function implements the default assignment of a matrix to a matrix.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT1  // Type of the left-hand side matrix
        , bool SO1      // Storage order of the left-hand side matrix
        , typename MT2  // Type of the right-hand side matrix
        , bool SO2 >    // Storage order of the right-hand side matrix
inline void assign( Matrix<MT1,SO1>& lhs, const Matrix<MT2,SO2>& rhs )
{
   BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( (~lhs).columns() == (~rhs).columns(), "Invalid number of columns" );

   (~lhs).assign( ~rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the addition assignment of a matrix to a matrix.
// \ingroup matrix
//
// \param lhs The target left-hand side matrix.
// \param rhs The right-hand side matrix to be added.
// \return void
//
// This function implements the default addition assignment of a matrix to a matrix.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT1  // Type of the left-hand side matrix
        , bool SO1      // Storage order of the left-hand side matrix
        , typename MT2  // Type of the right-hand side matrix
        , bool SO2 >    // Storage order of the right-hand side matrix
inline void addAssign( Matrix<MT1,SO1>& lhs, const Matrix<MT2,SO2>& rhs )
{
   BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( (~lhs).columns() == (~rhs).columns(), "Invalid number of columns" );

   (~lhs).addAssign( ~rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the subtraction assignment of a matrix to matrix.
// \ingroup matrix
//
// \param lhs The target left-hand side matrix.
// \param rhs The right-hand side matrix to be subtracted.
// \return void
//
// This function implements the default subtraction assignment of a matrix to a matrix.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT1  // Type of the left-hand side matrix
        , bool SO1      // Storage order of the left-hand side matrix
        , typename MT2  // Type of the right-hand side matrix
        , bool SO2 >    // Storage order of the right-hand side matrix
inline void subAssign( Matrix<MT1,SO1>& lhs, const Matrix<MT2,SO2>& rhs )
{
   BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( (~lhs).columns() == (~rhs).columns(), "Invalid number of columns" );

   (~lhs).subAssign( ~rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the multiplication assignment of a matrix to a matrix.
// \ingroup matrix
//
// \param lhs The target left-hand side matrix.
// \param rhs The right-hand side matrix to be multiplied.
// \return void
//
// This function implements the default multiplication assignment of a matrix to a matrix.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT1  // Type of the left-hand side matrix
        , bool SO1      // Storage order of the left-hand side matrix
        , typename MT2  // Type of the right-hand side matrix
        , bool SO2 >    // Storage order of the right-hand side matrix
inline void multAssign( Matrix<MT1,SO1>& lhs, const Matrix<MT2,SO2>& rhs )
{
   BLAZE_INTERNAL_ASSERT( (~lhs).columns() == (~rhs).rows(), "Invalid matrix sizes" );

   (~lhs).multAssign( ~rhs );
}
//*************************************************************************************************

} // namespace blaze

#endif
