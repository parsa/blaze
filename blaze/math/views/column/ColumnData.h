//=================================================================================================
/*!
//  \file blaze/math/views/column/ColumnData.h
//  \brief Header file for the implementation of the ColumnData class template
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

#ifndef _BLAZE_MATH_VIEWS_COLUMN_COLUMNDATA_H_
#define _BLAZE_MATH_VIEWS_COLUMN_COLUMNDATA_H_


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
/*!\brief Auxiliary class template for the data members of the Column class.
// \ingroup column
//
// The auxiliary ColumnData class template represents an abstraction of the data members of
// the Column class template. The necessary set of data members is selected depending on the
// number of compile time column arguments.
*/
template< typename MT       // Type of the matrix
        , size_t... CCAs >  // Compile time column arguments
struct ColumnData
{};
//*************************************************************************************************




//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR ZERO COMPILE TIME COLUMN INDICES
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Specialization of the ColumnData class template for zero compile time column arguments.
// \ingroup column
//
// This specialization of ColumnData adapts the class template to the requirements of zero compile
// time column arguments.
*/
template< typename MT >  // Type of the matrix
struct ColumnData<MT>
{
 public:
   //**Type definitions****************************************************************************
   //! Composite data type of the dense matrix expression.
   using Operand = If_< IsExpression<MT>, MT, MT& >;
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline ColumnData( Operand matrix, size_t index );
   // No explicitly declared copy constructor.
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   // No explicitly declared destructor.
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   ColumnData& operator=( const ColumnData& ) = delete;
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline Operand operand() const noexcept;
   inline size_t  column () const noexcept;
   //@}
   //**********************************************************************************************

 protected:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   Operand      matrix_;  //!< The matrix containing the column.
   const size_t column_;  //!< The index of the column in the matrix.
   //@}
   //**********************************************************************************************
};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The constructor for ColumnData.
//
// \param matrix The matrix containing the column.
// \param index The index of the column.
// \exception std::invalid_argument Invalid column access index.
*/
template< typename MT >  // Type of the matrix
inline ColumnData<MT>::ColumnData( Operand matrix, size_t index )
   : matrix_( matrix )  // The matrix containing the column
   , column_( index  )  // The index of the column in the matrix
{
   if( matrix_.columns() <= column() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid column access index" );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the matrix containing the column.
//
// \return The matrix containing the column.
*/
template< typename MT >  // Type of the matrix
inline typename ColumnData<MT>::Operand ColumnData<MT>::operand() const noexcept
{
   return matrix_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the index of the column of the underlying dense matrix.
//
// \return The index of the column.
*/
template< typename MT >  // Type of the matrix
inline size_t ColumnData<MT>::column() const noexcept
{
   return column_;
}
//*************************************************************************************************




//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR ONE COMPILE TIME COLUMN INDEX
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Specialization of the ColumnData class template for a single compile time column argument.
// \ingroup column
//
// This specialization of ColumnData adapts the class template to the requirements of a single
// compile time column argument.
*/
template< typename MT  // Type of the matrix
        , size_t I >   // Compile time column index
struct ColumnData<MT,I>
{
 public:
   //**Type definitions****************************************************************************
   //! Composite data type of the dense matrix expression.
   using Operand = If_< IsExpression<MT>, MT, MT& >;
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline ColumnData( Operand matrix );
   // No explicitly declared copy constructor.
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   // No explicitly declared destructor.
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   ColumnData& operator=( const ColumnData& ) = delete;
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline           Operand   operand() const noexcept;
   inline constexpr size_t    column () const noexcept;
   //@}
   //**********************************************************************************************

 protected:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   Operand matrix_;  //!< The matrix containing the column.
   //@}
   //**********************************************************************************************
};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The constructor for ColumnData.
//
// \param matrix The matrix containing the column.
// \exception std::invalid_argument Invalid column access index.
*/
template< typename MT  // Type of the matrix
        , size_t I >   // Compile time column index
inline ColumnData<MT,I>::ColumnData( Operand matrix )
   : matrix_( matrix )  // The matrix containing the column
{
   if( matrix_.columns() <= column() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid column access index" );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the matrix containing the column.
//
// \return The matrix containing the column.
*/
template< typename MT  // Type of the matrix
        , size_t I >   // Compile time column index
inline typename ColumnData<MT,I>::Operand ColumnData<MT,I>::operand() const noexcept
{
   return matrix_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the index of the column of the underlying dense matrix.
//
// \return The index of the column.
*/
template< typename MT  // Type of the matrix
        , size_t I >   // Compile time column index
inline constexpr size_t ColumnData<MT,I>::column() const noexcept
{
   return I;
}
//*************************************************************************************************

} // namespace blaze

#endif
