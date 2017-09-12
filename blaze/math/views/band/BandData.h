//=================================================================================================
/*!
//  \file blaze/math/views/band/BandData.h
//  \brief Header file for the implementation of the BandData class template
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

#ifndef _BLAZE_MATH_VIEWS_BAND_BANDDATA_H_
#define _BLAZE_MATH_VIEWS_BAND_BANDDATA_H_


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
/*!\brief Auxiliary class template for the data members of the Band class.
// \ingroup band
//
// The auxiliary BandData class template represents an abstraction of the data members of the
// Band class template. The necessary set of data member is selected depending on the number
// of compile time band arguments.
*/
template< typename MT          // Type of the matrix
        , ptrdiff_t... CBAs >  // Compile time band arguments
struct BandData
{};
//*************************************************************************************************




//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR ZERO COMPILE TIME BAND INDICES
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Specialization of the BandData class template for zero compile time band arguments.
// \ingroup band
//
// This specialization of BandData adapts the class template to the requirements of zero compile
// time band arguments.
*/
template< typename MT >  // Type of the matrix
struct BandData<MT>
{
 public:
   //**Type definitions****************************************************************************
   //! Composite data type of the dense matrix expression.
   using Operand = If_< IsExpression<MT>, MT, MT& >;
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline BandData( Operand matrix, ptrdiff_t index );
   // No explicitly declared copy constructor.
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   // No explicitly declared destructor.
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   BandData& operator=( const BandData& ) = delete;
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline Operand   operand() const noexcept;
   inline ptrdiff_t band   () const noexcept;
   inline size_t    row    () const noexcept;
   inline size_t    column () const noexcept;
   //@}
   //**********************************************************************************************

 protected:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   Operand         matrix_;  //!< The matrix containing the band.
   const ptrdiff_t band_;    //!< The band index.
   const size_t    row_;     //!< The index of the row containing the first element of the band.
   const size_t    column_;  //!< The index of the column containing the first element of the band.
   //@}
   //**********************************************************************************************
};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The constructor for BandData.
//
// \param matrix The matrix containing the band.
// \param index The index of the band.
// \exception std::invalid_argument Invalid band access index.
*/
template< typename MT >  // Type of the matrix
inline BandData<MT>::BandData( Operand matrix, ptrdiff_t index )
   : matrix_( matrix )                        // The matrix containing the band
   , band_  ( index  )                        // The band index
   , row_   ( index >= 0L ?   0UL : -index )  // The index of the row containing the first element of the band
   , column_( index >= 0L ? index :    0UL )  // The index of the column containing the first element of the band
{
   if( ( band() > 0L && column() >= matrix.columns() ) ||
       ( band() < 0L && row() >= matrix.rows() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid band access index" );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the matrix containing the band.
//
// \return The matrix containing the band.
*/
template< typename MT >  // Type of the matrix
inline typename BandData<MT>::Operand BandData<MT>::operand() const noexcept
{
   return matrix_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the index of the band of the underlying dense matrix.
//
// \return The index of the band.
*/
template< typename MT >  // Type of the matrix
inline ptrdiff_t BandData<MT>::band() const noexcept
{
   return band_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the index of the row containing the first element of the band.
//
// \return The first row index.
*/
template< typename MT >  // Type of the matrix
inline size_t BandData<MT>::row() const noexcept
{
   return row_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the index of the column containing the first element of the band.
//
// \return The first column index.
*/
template< typename MT >  // Type of the matrix
inline size_t BandData<MT>::column() const noexcept
{
   return column_;
}
//*************************************************************************************************




//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR ONE COMPILE TIME BAND INDEX
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Specialization of the BandData class template for a single compile time band argument.
// \ingroup band
//
// This specialization of BandData adapts the class template to the requirements of a single
// compile time band argument.
*/
template< typename MT    // Type of the matrix
        , ptrdiff_t I >  // Compile time band index
struct BandData<MT,I>
{
 public:
   //**Type definitions****************************************************************************
   //! Composite data type of the dense matrix expression.
   using Operand = If_< IsExpression<MT>, MT, MT& >;
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline BandData( Operand matrix );
   // No explicitly declared copy constructor.
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   // No explicitly declared destructor.
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   BandData& operator=( const BandData& ) = delete;
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline           Operand   operand() const noexcept;
   inline constexpr ptrdiff_t band   () const noexcept;
   inline constexpr size_t    row    () const noexcept;
   inline constexpr size_t    column () const noexcept;
   //@}
   //**********************************************************************************************

 protected:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   Operand matrix_;  //!< The matrix containing the band.
   //@}
   //**********************************************************************************************
};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The constructor for BandData.
//
// \param matrix The matrix containing the band.
// \exception std::invalid_argument Invalid band access index.
*/
template< typename MT    // Type of the matrix
        , ptrdiff_t I >  // Compile time band index
inline BandData<MT,I>::BandData( Operand matrix )
   : matrix_( matrix )  // The matrix containing the band
{
   if( ( band() > 0L && column() >= matrix.columns() ) ||
       ( band() < 0L && row() >= matrix.rows() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid band access index" );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the matrix containing the band.
//
// \return The matrix containing the band.
*/
template< typename MT    // Type of the matrix
        , ptrdiff_t I >  // Compile time band index
inline typename BandData<MT,I>::Operand BandData<MT,I>::operand() const noexcept
{
   return matrix_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the index of the band of the underlying dense matrix.
//
// \return The index of the band.
*/
template< typename MT    // Type of the matrix
        , ptrdiff_t I >  // Compile time band index
inline constexpr ptrdiff_t BandData<MT,I>::band() const noexcept
{
   return I;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the index of the row containing the first element of the band.
//
// \return The first row index.
*/
template< typename MT    // Type of the matrix
        , ptrdiff_t I >  // Compile time band index
inline constexpr size_t BandData<MT,I>::row() const noexcept
{
   return ( I >= 0L ? 0UL : -I );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the index of the column containing the first element of the band.
//
// \return The first column index.
*/
template< typename MT    // Type of the matrix
        , ptrdiff_t I >  // Compile time band index
inline constexpr size_t BandData<MT,I>::column() const noexcept
{
   return ( I >= 0L ? I : 0UL );
}
//*************************************************************************************************

} // namespace blaze

#endif
