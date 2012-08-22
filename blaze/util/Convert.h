//=================================================================================================
/*!
//  \file blaze/util/Convert.h
//  \brief Conversion functionality
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

#ifndef _BLAZE_UTIL_CONVERT_H_
#define _BLAZE_UTIL_CONVERT_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <cstdlib>
#include <sstream>
#include <stdexcept>
#include <string>
#include <typeinfo>
#include <blaze/util/NonCreatable.h>
#include <blaze/util/typetraits/IsBaseOf.h>


namespace blaze {

//=================================================================================================
//
//  CLASS TEMPLATE CASTCONVERTER
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Base class for the CastConverter specializations.
// \ingroup util
*/
template< typename To, typename From, int IsBase >
struct CastConverter
{};
/*! \endcond */
//*************************************************************************************************





//=================================================================================================
//
//  PARTIAL TEMPLATE SPECIALIZATION OF CASTCONVERTER
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Dynamic cast converter between pointers of type \a From and \a To.
// \ingroup util
*/
template< typename To, typename From >
struct CastConverter<To*,From*,0> : private NonCreatable
{
 public:
   //**Conversion functions************************************************************************
   /*!\name Conversion functions */
   //@{
   static inline To* convert( From* from );
   //@}
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Dynamic cast between pointers of type \a From and type \a To.
//
// \param from The pointer of type \a From to be converted.
// \return The converted pointer of type \a To.
*/
template< typename To, typename From >
inline To* CastConverter<To*,From*,0>::convert( From* from )
{
   return dynamic_cast<To*>( from );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  PARTIAL TEMPLATE SPECIALIZATION OF CASTCONVERTER
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Static cast converter between pointers of type \a From and \a To.
// \ingroup util
*/
template< typename To, typename From >
struct CastConverter<To*,From*,1> : private NonCreatable
{
 public:
   //**Conversion functions************************************************************************
   /*!\name Conversion functions */
   //@{
   static inline To* convert( From* from );
   //@}
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Static cast between pointers of type \a From and type \a To.
//
// \param from The pointer of type \a From to be converted.
// \return The converted pointer of type \a To.
*/
template< typename To, typename From >
inline To* CastConverter<To*,From*,1>::convert( From* from )
{
   return static_cast<To*>( from );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  CLASS TEMPLATE CONVERTER
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Converter class between type \a From and type \a To.
// \ingroup util
*/
template< typename To, typename From >
struct Converter : private NonCreatable
{
 public:
   //**Conversion functions************************************************************************
   /*!\name Conversion functions */
   //@{
   static inline To convert( const From& from );
   //@}
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Conversion from type \a From to type \a To.
//
// \param from The data value to be converted.
// \return The converted data value.
*/
template< typename To, typename From >
inline To Converter<To,From>::convert( const From& from )
{
   return static_cast<To>( from );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  TOTAL TEMPLATE SPECIALIZATION OF CONVERTER
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Converter class between pointers of type \a From and type \a To.
// \ingroup util
*/
template< typename To, typename From >
struct Converter<To*,From*> : private NonCreatable
{
 public:
   //**Conversion functions************************************************************************
   /*!\name Conversion functions */
   //@{
   static inline To* convert( From* from );
   //@}
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Conversion between pointers of type \a From and type \a To.
//
// \param from The pointer of type \a From to be converted.
// \return The converted pointer of type \a To.
*/
template< typename To, typename From >
inline To* Converter<To*,From*>::convert( From* from )
{
   return CastConverter<To*,From*,blaze::IsBaseOf<To,From>::yes>::convert( from );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  PARTIAL TEMPLATE SPECIALIZATION OF CONVERTER
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Converter class between type \a std::string and type \a To.
// \ingroup util
*/
template< typename To >
struct Converter<To,std::string> : private NonCreatable
{
 public:
   //**Conversion functions************************************************************************
   /*!\name Conversion functions */
   //@{
   static inline To convert( const std::string& from );
   //@}
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Conversion from type \a std::string to type \a To.
//
// \param from The \a std::string value to be converted.
// \return The converted data value.
*/
template< typename To >
inline To Converter<To,std::string>::convert( const std::string& from )
{
   To to;
   std::istringstream iss( from );
   if( !(iss >> to) ) {
      std::ostringstream error;
      error << "Invalid cast from std::string to " << typeid(to).name() << "\n";
      throw std::runtime_error( error.str() );
   }
   return to;
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  PARTIAL TEMPLATE SPECIALIZATION OF CONVERTER
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Converter class between type \a From and type \a std::string.
// \ingroup util
*/
template< typename From >
struct Converter<std::string,From> : private NonCreatable
{
 public:
   //**Conversion functions************************************************************************
   /*!\name Conversion functions */
   //@{
   static inline std::string convert( const From& from );
   //@}
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Conversion from type \a From to type \a std::string.
//
// \param from The data value to be converted.
// \return The converted data value.
*/
template< typename From >
inline std::string Converter<std::string,From>::convert( const From& from )
{
   std::ostringstream oss;
   if( !(oss << from) ) {
      std::ostringstream error;
      error << "Invalid cast from " << typeid(from).name() << " to std::string\n";
      throw std::runtime_error( error.str() );
   }
   return oss.str();
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  TOTAL TEMPLATE SPECIALIZATION OF CONVERTER
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Resolution of the ambiguity for the conversion \a std::string to \a std::string
// \ingroup util
//
// This converter resolves the ambiguity for the instantiation of Converter<std::string,std::string>.
*/
template<>
struct Converter<std::string,std::string> : private NonCreatable
{
 public:
   //**Conversion functions************************************************************************
   /*!\name Conversion functions */
   //@{
   static inline std::string convert( const std::string& from );
   //@}
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Conversion from \a std::string to \a std::string.
//
// \param from The string to be converted.
// \return The converted string.
*/
inline std::string Converter<std::string,std::string>::convert( const std::string& from )
{
   return from;
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  CONVERSION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Conversion from type \a From to type \a To.
//
// \param from The data value to be converted.
// \return The converted data value.
//
// The \a convert function transforms the data value \a from of type \a From to the data type
// \a To. The syntax for this operation is similar to the C++ cast operators. For example, in
// order to convert a built-in integer value \p integer to a \a std::string, use
//
//                    \code convert<std::string>( integer ) \endcode
//
// The \a convert function supports any possible type conversion in the most efficient way.
*/
template< typename To, typename From >
inline To convert( const From& from )
{
   return Converter<To,From>::convert( from );
}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Conversion from a \a std::string to an integer value.
//
// \param from The string to be converted.
// \return The converted integer value.
*/
template<>
inline int convert<int,std::string>( const std::string& from )
{
   return std::atoi( from.c_str() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Conversion from a \a std::string to an unsigned integer value.
//
// \param from The string to be converted.
// \return The converted unsigned integer value.
*/
template<>
inline unsigned int convert<unsigned int,std::string>( const std::string& from )
{
   return static_cast<unsigned int>( std::atoi( from.c_str() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Conversion from a \a std::string to a float value.
//
// \param from The string to be converted.
// \return The converted float value.
*/
template<>
inline float convert<float,std::string>( const std::string& from )
{
   return static_cast<float>( std::atof( from.c_str() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Conversion from a \a std::string to a double value.
//
// \param from The string to be converted.
// \return The converted double value.
*/
template<>
inline double convert<double,std::string>( const std::string& from )
{
   return std::atof( from.c_str() );
}
/*! \endcond */
//*************************************************************************************************




//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Conversion from a constant character array to type \a To.
//
// \param from The constant character array to be converted.
// \return The converted data value.
*/
template< typename To >
inline To convert( const char* const from )
{
   return Converter<To,std::string>::convert( from );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Conversion from a constant character array to an integer value.
//
// \param from The constant character array to be converted.
// \return The converted integer value.
*/
template<>
inline int convert<int>( const char* const from )
{
   return std::atoi( from );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Conversion from a constant character array to an unsigned integer value.
//
// \param from The constant character array to be converted.
// \return The converted unsigned integer value.
*/
template<>
inline unsigned int convert<unsigned int>( const char* const from )
{
   return static_cast<unsigned int>( std::atof( from ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Conversion from a constant character array to a \a float value.
//
// \param from The constant character array to be converted.
// \return The converted float value.
*/
template<>
inline float convert<float>( const char* const from )
{
   return static_cast<float>( std::atof( from ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Conversion from a constant character array to a \a double value.
//
// \param from The constant character array to be converted.
// \return The converted double value.
*/
template<>
inline double convert<double>( const char* const from )
{
   return std::atof( from );
}
/*! \endcond */
//*************************************************************************************************




//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Conversion from a non-constant character array to type \a To.
//
// \param from The non-constant character array to be converted.
// \return The converted data value.
*/
template< typename To >
inline To convert( char* const from )
{
   return Converter<To,std::string>::convert( from );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Conversion from a non-constant character array to an integer value.
//
// \param from The non-constant character array to be converted.
// \return The converted integer value.
*/
template<>
inline int convert<int>( char* const from )
{
   return std::atoi( from );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Conversion from a non-constant character array to an unsigned integer value.
//
// \param from The non-constant character array to be converted.
// \return The converted unsigned integer value.
*/
template<>
inline unsigned int convert<unsigned int>( char* const from )
{
   return static_cast<unsigned int>( std::atoi( from ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Conversion from a non-constant character array to a \a float value.
//
// \param from The non-constant character array to be converted.
// \return The converted float value.
*/
template<>
inline float convert<float>( char* const from )
{
   return static_cast<float>( std::atof( from ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Conversion from a non-constant character array to a \a double value.
//
// \param from The non-constant character array to be converted.
// \return The converted double value.
*/
template<>
inline double convert<double>( char* const from )
{
   return std::atof( from );
}
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
