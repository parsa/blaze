//=================================================================================================
/*!
//  \file blaze/math/views/row/RowData.h
//  \brief Header file for the implementation of the RowData class template
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

#ifndef _BLAZE_MATH_VIEWS_ROW_ROWDATA_H_
#define _BLAZE_MATH_VIEWS_ROW_ROWDATA_H_


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
/*!\brief Auxiliary class template for the data members of the Row class.
// \ingroup row
//
// The auxiliary RowData class template represents an abstraction of the data members of the
// Row class template. The necessary set of data members is selected depending on the number
// of compile time row arguments.
*/
template< typename MT       // Type of the matrix
        , size_t... CRAs >  // Compile time row arguments
struct RowData
{};
//*************************************************************************************************




//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR ZERO COMPILE TIME ROW INDICES
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Specialization of the RowData class template for zero compile time row arguments.
// \ingroup row
//
// This specialization of RowData adapts the class template to the requirements of zero compile
// time row arguments.
*/
template< typename MT >  // Type of the matrix
struct RowData<MT>
{
 public:
   //**Type definitions****************************************************************************
   //! Composite data type of the dense matrix expression.
   using Operand = If_< IsExpression<MT>, MT, MT& >;
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline RowData( Operand matrix, size_t index );
   // No explicitly declared copy constructor.
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   // No explicitly declared destructor.
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   RowData& operator=( const RowData& ) = delete;
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline Operand operand() const noexcept;
   inline size_t  row    () const noexcept;
   //@}
   //**********************************************************************************************

 protected:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   Operand      matrix_;  //!< The matrix containing the row.
   const size_t row_;     //!< The index of the row in the matrix.
   //@}
   //**********************************************************************************************
};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The constructor for RowData.
//
// \param matrix The matrix containing the row.
// \param index The index of the row.
// \exception std::invalid_argument Invalid row access index.
*/
template< typename MT >  // Type of the matrix
inline RowData<MT>::RowData( Operand matrix, size_t index )
   : matrix_( matrix )  // The matrix containing the row
   , row_   ( index  )  // The index of the row in the matrix
{
   if( matrix_.rows() <= row() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid row access index" );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the matrix containing the row.
//
// \return The matrix containing the row.
*/
template< typename MT >  // Type of the matrix
inline typename RowData<MT>::Operand RowData<MT>::operand() const noexcept
{
   return matrix_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the index of the row of the underlying dense matrix.
//
// \return The index of the row.
*/
template< typename MT >  // Type of the matrix
inline size_t RowData<MT>::row() const noexcept
{
   return row_;
}
//*************************************************************************************************




//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR ONE COMPILE TIME ROW INDEX
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Specialization of the RowData class template for a single compile time row argument.
// \ingroup row
//
// This specialization of RowData adapts the class template to the requirements of a single
// compile time row argument.
*/
template< typename MT     // Type of the matrix
        , size_t Index >  // Compile time row index
struct RowData<MT,Index>
{
 public:
   //**Type definitions****************************************************************************
   //! Composite data type of the dense matrix expression.
   using Operand = If_< IsExpression<MT>, MT, MT& >;
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline RowData( Operand matrix );
   // No explicitly declared copy constructor.
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   // No explicitly declared destructor.
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   RowData& operator=( const RowData& ) = delete;
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline           Operand   operand() const noexcept;
   inline constexpr size_t    row    () const noexcept;
   //@}
   //**********************************************************************************************

 protected:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   Operand matrix_;  //!< The matrix containing the row.
   //@}
   //**********************************************************************************************
};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The constructor for RowData.
//
// \param matrix The matrix containing the row.
// \exception std::invalid_argument Invalid row access index.
*/
template< typename MT     // Type of the matrix
        , size_t Index >  // Compile time row index
inline RowData<MT,Index>::RowData( Operand matrix )
   : matrix_( matrix )  // The matrix containing the row
{
   if( matrix_.rows() <= row() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid row access index" );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the matrix containing the row.
//
// \return The matrix containing the row.
*/
template< typename MT     // Type of the matrix
        , size_t Index >  // Compile time row index
inline typename RowData<MT,Index>::Operand RowData<MT,Index>::operand() const noexcept
{
   return matrix_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the index of the row of the underlying dense matrix.
//
// \return The index of the row.
*/
template< typename MT     // Type of the matrix
        , size_t Index >  // Compile time row index
inline constexpr size_t RowData<MT,Index>::row() const noexcept
{
   return Index;
}
//*************************************************************************************************

} // namespace blaze

#endif
