//=================================================================================================
/*!
//  \file blaze/math/views/submatrix/SubmatrixData.h
//  \brief Header file for the implementation of the SubmatrixData class template
//
//  Copyright (C) 2012-2017 Klaus Iglberger - All Rights Reserved
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

#ifndef _BLAZE_MATH_VIEWS_SUBMATRIX_SUBMATRIXDATA_H_
#define _BLAZE_MATH_VIEWS_SUBMATRIX_SUBMATRIXDATA_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/Exception.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Auxiliary class template for the data members of the Submatrix class.
// \ingroup submatrix
//
// The auxiliary SubmatrixData class template represents an abstraction of the data members of
// the Submatrix class template. The necessary set of data members is selected depending on the
// number of compile time submatrix arguments.
*/
template< typename MT       // Type of the matrix
        , size_t... CSAs >  // Compile time submatrix arguments
struct SubmatrixData
{};
//*************************************************************************************************




//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR ZERO COMPILE TIME ARGUMENTS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Specialization of the SubmatrixData class template for zero compile time submatrix
//        arguments.
// \ingroup submatrix
//
// This specialization of SubmatrixData adapts the class template to the requirements of zero
// compile time submatrix arguments.
*/
template< typename MT >  // Type of the matrix
struct SubmatrixData<MT>
{
 public:
   //**Type definitions****************************************************************************
   //! Composite data type of the matrix expression.
   using Operand = If_< IsExpression<MT>, MT, MT& >;
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline SubmatrixData( Operand matrix, size_t rindex, size_t cindex, size_t m, size_t n );
   // No explicitly declared copy constructor.
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   // No explicitly declared destructor.
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   SubmatrixData& operator=( const SubmatrixData& ) = delete;
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline Operand operand() const noexcept;
   inline size_t  row    () const noexcept;
   inline size_t  column () const noexcept;
   inline size_t  rows   () const noexcept;
   inline size_t  columns() const noexcept;
   //@}
   //**********************************************************************************************

 protected:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   Operand      matrix_;  //!< The matrix containing the submatrix.
   const size_t row_;     //!< The first row of the submatrix.
   const size_t column_;  //!< The first column of the submatrix.
   const size_t m_;       //!< The number of rows of the submatrix.
   const size_t n_;       //!< The number of columns of the submatrix.
   //@}
   //**********************************************************************************************
};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The constructor for SubmatrixData.
//
// \param matrix The matrix containing the submatrix.
// \param rindex The index of the first row of the submatrix in the given matrix.
// \param cindex The index of the first column of the submatrix in the given matrix.
// \param m The number of rows of the submatrix.
// \param n The number of columns of the submatrix.
// \exception std::invalid_argument Invalid submatrix specification.
*/
template< typename MT >  // Type of the matrix
inline SubmatrixData<MT>::SubmatrixData( Operand matrix, size_t rindex, size_t cindex, size_t m, size_t n )
   : matrix_( matrix )  // The matrix containing the submatrix
   , row_   ( rindex )  // The first row of the submatrix
   , column_( cindex )  // The first column of the submatrix
   , m_     ( m      )  // The number of rows of the submatrix
   , n_     ( n      )  // The number of columns of the submatrix
{
   if( ( row_ + m_ > matrix_.rows() ) || ( column_ + n_ > matrix_.columns() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid submatrix specification" );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the matrix containing the submatrix.
//
// \return The matrix containing the submatrix.
*/
template< typename MT >  // Type of the matrix
inline typename SubmatrixData<MT>::Operand SubmatrixData<MT>::operand() const noexcept
{
   return matrix_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the index of the first row of the submatrix in the underlying matrix.
//
// \return The index of the first row.
*/
template< typename MT >  // Type of the matrix
inline size_t SubmatrixData<MT>::row() const noexcept
{
   return row_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the index of the first column of the submatrix in the underlying matrix.
//
// \return The index of the first column.
*/
template< typename MT >  // Type of the matrix
inline size_t SubmatrixData<MT>::column() const noexcept
{
   return column_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the number of rows of the submatrix.
//
// \return The number of rows of the submatrix.
*/
template< typename MT >  // Type of the matrix
inline size_t SubmatrixData<MT>::rows() const noexcept
{
   return m_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the number of columns of the submatrix.
//
// \return The number of columns of the submatrix.
*/
template< typename MT >  // Type of the matrix
inline size_t SubmatrixData<MT>::columns() const noexcept
{
   return n_;
}
//*************************************************************************************************




//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR FOUR COMPILE TIME ARGUMENTS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Specialization of the SubmatrixData class template for four compile time submatrix
//        arguments.
// \ingroup submatrix
//
// This specialization of SubmatrixData adapts the class template to the requirements of two
// compile time arguments.
*/
template< typename MT  // Type of the matrix
        , size_t I     // Index of the first row
        , size_t J     // Index of the first column
        , size_t M     // Number of rows
        , size_t N >   // Number of columns
struct SubmatrixData<MT,I,J,M,N>
{
 public:
   //**Type definitions****************************************************************************
   //! Composite data type of the matrix expression.
   using Operand = If_< IsExpression<MT>, MT, MT& >;
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline SubmatrixData( Operand matrix );
   // No explicitly declared copy constructor.
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   // No explicitly declared destructor.
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   SubmatrixData& operator=( const SubmatrixData& ) = delete;
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline Operand operand() const noexcept;
   inline size_t  row    () const noexcept;
   inline size_t  column () const noexcept;
   inline size_t  rows   () const noexcept;
   inline size_t  columns() const noexcept;
   //@}
   //**********************************************************************************************

 protected:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   Operand matrix_;  //!< The matrix containing the submatrix.
   //@}
   //**********************************************************************************************
};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The constructor for SubmatrixData.
//
// \param matrix The matrix containing the submatrix.
// \exception std::invalid_argument Invalid submatrix specification.
*/
template< typename MT  // Type of the matrix
        , size_t I     // Index of the first row
        , size_t J     // Index of the first column
        , size_t M     // Number of rows
        , size_t N >   // Number of columns
inline SubmatrixData<MT,I,J,M,N>::SubmatrixData( Operand matrix )
   : matrix_( matrix )  // The matrix containing the submatrix
{
   if( ( I + M > matrix_.rows() ) || ( J + N > matrix_.columns() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid submatrix specification" );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the matrix containing the submatrix.
//
// \return The matrix containing the submatrix.
*/
template< typename MT  // Type of the matrix
        , size_t I     // Index of the first row
        , size_t J     // Index of the first column
        , size_t M     // Number of rows
        , size_t N >   // Number of columns
inline typename SubmatrixData<MT,I,J,M,N>::Operand
   SubmatrixData<MT,I,J,M,N>::operand() const noexcept
{
   return matrix_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the index of the first row of the submatrix in the underlying matrix.
//
// \return The index of the first row.
*/
template< typename MT  // Type of the matrix
        , size_t I     // Index of the first row
        , size_t J     // Index of the first column
        , size_t M     // Number of rows
        , size_t N >   // Number of columns
inline size_t SubmatrixData<MT,I,J,M,N>::row() const noexcept
{
   return I;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the index of the first column of the submatrix in the underlying matrix.
//
// \return The index of the first column.
*/
template< typename MT  // Type of the matrix
        , size_t I     // Index of the first row
        , size_t J     // Index of the first column
        , size_t M     // Number of rows
        , size_t N >   // Number of columns
inline size_t SubmatrixData<MT,I,J,M,N>::column() const noexcept
{
   return J;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the number of rows of the submatrix.
//
// \return The number of rows of the submatrix.
*/
template< typename MT  // Type of the matrix
        , size_t I     // Index of the first row
        , size_t J     // Index of the first column
        , size_t M     // Number of rows
        , size_t N >   // Number of columns
inline size_t SubmatrixData<MT,I,J,M,N>::rows() const noexcept
{
   return M;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the number of columns of the submatrix.
//
// \return The number of columns of the submatrix.
*/
template< typename MT  // Type of the matrix
        , size_t I     // Index of the first row
        , size_t J     // Index of the first column
        , size_t M     // Number of rows
        , size_t N >   // Number of columns
inline size_t SubmatrixData<MT,I,J,M,N>::columns() const noexcept
{
   return N;
}
//*************************************************************************************************

} // namespace blaze

#endif
