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
template< size_t... CCAs >  // Compile time column arguments
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
template<>
struct ColumnData<>
{
 public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline constexpr ColumnData( size_t index );
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
   inline constexpr size_t column() const noexcept;
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   const size_t column_;  //!< The index of the column in the matrix.
   //@}
   //**********************************************************************************************
};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The constructor for ColumnData.
//
// \param index The index of the column.
*/
inline constexpr ColumnData<>::ColumnData( size_t index )
   : column_( index )  // The index of the column in the matrix
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the index of the column of the underlying dense matrix.
//
// \return The index of the column.
*/
inline constexpr size_t ColumnData<>::column() const noexcept
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
template< size_t I >   // Compile time column index
struct ColumnData<I>
{
 public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline constexpr ColumnData();
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
   inline constexpr size_t column() const noexcept;
   //@}
   //**********************************************************************************************
};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The constructor for ColumnData.
*/
template< size_t I >  // Compile time column index
inline constexpr ColumnData<I>::ColumnData()
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the index of the column of the underlying dense matrix.
//
// \return The index of the column.
*/
template< size_t I >  // Compile time column index
inline constexpr size_t ColumnData<I>::column() const noexcept
{
   return I;
}
//*************************************************************************************************

} // namespace blaze

#endif
