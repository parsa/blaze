//=================================================================================================
/*!
//  \file blaze/math/views/elements/ElementsData.h
//  \brief Header file for the implementation of the ElementsData class template
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

#ifndef _BLAZE_MATH_VIEWS_ELEMENTS_ELEMENTSDATA_H_
#define _BLAZE_MATH_VIEWS_ELEMENTS_ELEMENTSDATA_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/Assert.h>
#include <blaze/util/SmallVector.h>
#include <blaze/util/Types.h>
#include <blaze/util/Unused.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR TWO COMPILE TIME ARGUMENTS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Auxiliary class template for the data members of the Elements class.
// \ingroup elements
//
// The auxiliary ElementsData class template represents an abstraction of the data members of the
// Elements class template. The necessary set of data members is selected depending on the number
// of compile time element arguments. The basic implementation of ElementsData adapts the class
// template to the requirements of multiple compile time element arguments.
*/
template< size_t... CEAs >  // Compile time element arguments
struct ElementsData
{
 public:
   //**Type definitions****************************************************************************
   using Indices = std::array<size_t,sizeof...(CEAs)>;  //!< Type of the container for element indices.
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   template< typename... REAs >
   explicit inline ElementsData( REAs... args ) noexcept;
   // No explicitly declared copy constructor.
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   // No explicitly declared destructor.
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   ElementsData& operator=( const ElementsData& ) = delete;
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   static inline constexpr const Indices& idces() noexcept;
   static inline constexpr size_t         idx  ( size_t i ) noexcept;
   static inline constexpr size_t         size () noexcept;
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   static constexpr Indices indices_{ { CEAs... } };  //!< The indices of the elements in the vector.
   //@}
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
// Definition and initialization of the static member variables
template< size_t... CEAs >  // Compile time element arguments
constexpr typename ElementsData<CEAs...>::Indices ElementsData<CEAs...>::indices_;
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief The constructor for ElementsData.
//
// \param args The optional element arguments.
*/
template< size_t... CEAs >    // Compile time element arguments
template< typename... REAs >  // Optional element arguments
inline ElementsData<CEAs...>::ElementsData( REAs... args ) noexcept
{
   UNUSED_PARAMETER( args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the indices of the specified elements in the underlying vector.
//
// \return The indices of the specified elements.
*/
template< size_t... CEAs >  // Compile time element arguments
inline constexpr const typename ElementsData<CEAs...>::Indices& ElementsData<CEAs...>::idces() noexcept
{
   return indices_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the index of the specified element in the underlying vector.
//
// \param i Access index for the element.
// \return The index of the specified element.
*/
template< size_t... CEAs >  // Compile time element arguments
inline constexpr size_t ElementsData<CEAs...>::idx( size_t i ) noexcept
{
   BLAZE_USER_ASSERT( i < size(), "Invalid element access index" );
   return indices_[i];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of elements.
//
// \return The number of elements.
*/
template< size_t... CEAs >  // Compile time element arguments
inline constexpr size_t ElementsData<CEAs...>::size() noexcept
{
   return sizeof...( CEAs );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR ZERO COMPILE TIME ELEMENT ARGUMENTS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the ElementsData class template for zero compile time element arguments.
// \ingroup elements
//
// This specialization of ElementsData adapts the class template to the requirements of zero
// compile time element arguments.
*/
template<>
struct ElementsData<>
{
 public:
   //**Type definitions****************************************************************************
   using Indices = SmallVector<size_t,8UL>;  //!< Type of the container for element indices.
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   template< typename T, typename... REAs >
   explicit inline ElementsData( const T* indices, size_t n, REAs... args );

   inline ElementsData( const ElementsData& ) = default;
   inline ElementsData( ElementsData&& ) = default;
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   // No explicitly declared destructor.
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   ElementsData& operator=( const ElementsData& ) = delete;
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline const Indices& idces() const noexcept;
   inline size_t         idx  ( size_t i ) const noexcept;
   inline size_t         size () const noexcept;
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   Indices indices_;  //!< The indices of the elements in the vector.
   //@}
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief The constructor for ElementsData.
//
// \param indices Pointer to the first index of the selected elements.
// \param n The total number of indices.
// \param args The optional element arguments.
*/
template< typename T          // Type of the element indices
        , typename... REAs >  // Optional element arguments
inline ElementsData<>::ElementsData( const T* indices, size_t n, REAs... args )
   : indices_( indices, indices+n )  // The indices of the elements in the vector
{
   UNUSED_PARAMETER( args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the indices of the specified elements in the underlying vector.
//
// \return The indices of the specified elements.
*/
inline const ElementsData<>::Indices& ElementsData<>::idces() const noexcept
{
   return indices_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the index of the specified element in the underlying vector.
//
// \param i Access index for the element.
// \return The index of the specified element.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active.
*/
inline size_t ElementsData<>::idx( size_t i ) const noexcept
{
   BLAZE_USER_ASSERT( i < size(), "Invalid element access index" );
   return indices_[i];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of elements.
//
// \return The number of elements.
*/
inline size_t ElementsData<>::size() const noexcept
{
   return indices_.size();
}
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
