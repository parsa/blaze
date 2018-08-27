//=================================================================================================
/*!
//  \file blaze/math/views/columns/ColumnsData.h
//  \brief Header file for the implementation of the ColumnsData class template
//
//  Copyright (C) 2012-2018 Klaus Iglberger - All Rights Reserved
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

#ifndef _BLAZE_MATH_VIEWS_COLUMNS_COLUMNSDATA_H_
#define _BLAZE_MATH_VIEWS_COLUMNS_COLUMNSDATA_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/Assert.h>
#include <blaze/util/SmallArray.h>
#include <blaze/util/Types.h>
#include <blaze/util/Unused.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Auxiliary class template for the data members of the Columns class.
// \ingroup columns
//
// The auxiliary ColumnsData class template represents an abstraction of the data members of the
// Columns class template. The necessary set of data members is selected depending on the number
// of compile time column arguments. The basic implementation of ColumnsData adapts the class
// template to the requirements of multiple compile time column arguments.
*/
template< size_t... CCAs >  // Compile time column arguments
struct ColumnsData
{
 public:
   //**Type definitions****************************************************************************
   using Indices = std::array<size_t,sizeof...(CCAs)>;  //!< Type of the container for column indices.
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   template< typename... RCAs >
   explicit inline ColumnsData( RCAs... args ) noexcept;

   ColumnsData( const ColumnsData& ) = default;
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~ColumnsData() = default;
   //@}
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   /*!\name Assignment operators */
   //@{
   ColumnsData& operator=( const ColumnsData& ) = delete;
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   static inline constexpr const Indices& idces  () noexcept;
   static inline constexpr size_t         idx    ( size_t i ) noexcept;
   static inline constexpr size_t         columns() noexcept;
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   static constexpr Indices indices_{ { CCAs... } };  //!< The indices of the columns in the matrix.
   //@}
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
// Definition and initialization of the static member variables
template< size_t... CCAs >  // Compile time column arguments
constexpr typename ColumnsData<CCAs...>::Indices ColumnsData<CCAs...>::indices_;
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief The constructor for ColumnsData.
//
// \param args The optional column arguments.
*/
template< size_t... CCAs >    // Compile time column arguments
template< typename... RCAs >  // Optional column arguments
inline ColumnsData<CCAs...>::ColumnsData( RCAs... args ) noexcept
{
   UNUSED_PARAMETER( args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the indices of the specified columns in the underlying matrix.
//
// \return The indices of the specified columns.
*/
template< size_t... CCAs >  // Compile time column arguments
inline constexpr const typename ColumnsData<CCAs...>::Indices& ColumnsData<CCAs...>::idces() noexcept
{
   return indices_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the index of the specified column in the underlying matrix.
//
// \param i Access index for the column.
// \return The index of the specified column.
*/
template< size_t... CCAs >  // Compile time column arguments
inline constexpr size_t ColumnsData<CCAs...>::idx( size_t i ) noexcept
{
   BLAZE_USER_ASSERT( i < columns(), "Invalid column access index" );
   return indices_[i];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of columns.
//
// \return The number of columns.
*/
template< size_t... CCAs >  // Compile time column arguments
inline constexpr size_t ColumnsData<CCAs...>::columns() noexcept
{
   return sizeof...( CCAs );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR ZERO COMPILE TIME COLUMN ARGUMENTS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the ColumnsData class template for zero compile time column arguments.
// \ingroup columns
//
// This specialization of ColumnsData adapts the class template to the requirements of zero compile
// time column arguments.
*/
template<>
struct ColumnsData<>
{
 public:
   //**Type definitions****************************************************************************
   using Indices = SmallArray<size_t,8UL>;  //!< Type of the container for column indices.
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   template< typename T, typename... RCAs >
   explicit inline ColumnsData( const T* indices, size_t n, RCAs... args );

   ColumnsData( const ColumnsData& ) = default;
   ColumnsData( ColumnsData&& ) = default;
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~ColumnsData() = default;
   //@}
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   /*!\name Assignment operators */
   //@{
   ColumnsData& operator=( const ColumnsData& ) = delete;
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline const Indices& idces  () const noexcept;
   inline size_t         idx    ( size_t i ) const noexcept;
   inline size_t         columns() const noexcept;
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   Indices indices_;  //!< The indices of the columns in the matrix.
   //@}
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief The constructor for ColumnsData.
//
// \param indices Pointer to the first index of the selected columns.
// \param n The total number of indices.
// \param args The optional column arguments.
*/
template< typename T          // Type of the column indices
        , typename... RCAs >  // Optional column arguments
inline ColumnsData<>::ColumnsData( const T* indices, size_t n, RCAs... args )
   : indices_( indices, indices+n )  // The indices of the columns in the matrix
{
   UNUSED_PARAMETER( args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the indices of the specified columns in the underlying matrix.
//
// \return The indices of the specified columns.
*/
inline const ColumnsData<>::Indices& ColumnsData<>::idces() const noexcept
{
   return indices_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the index of the specified column in the underlying matrix.
//
// \param i Access index for the column.
// \return The index of the specified column.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active.
*/
inline size_t ColumnsData<>::idx( size_t i ) const noexcept
{
   BLAZE_USER_ASSERT( i < columns(), "Invalid column access index" );
   return indices_[i];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of columns.
//
// \return The number of columns.
*/
inline size_t ColumnsData<>::columns() const noexcept
{
   return indices_.size();
}
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
