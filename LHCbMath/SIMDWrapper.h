/*****************************************************************************\
* (c) Copyright 2000-2018 CERN for the benefit of the LHCb Collaboration      *
*                                                                             *
* This software is distributed under the terms of the GNU General Public      *
* Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   *
*                                                                             *
* In applying this licence, CERN does not waive the privileges and immunities *
* granted to it by virtue of its status as an Intergovernmental Organization  *
* or submit itself to any jurisdiction.                                       *
\*****************************************************************************/

#pragma once

#include "LHCbMath/TypeMapping.h"
#include "LHCbMath/bit_cast.h"

#include <immintrin.h>
#include <limits>

#include <algorithm>
#include <array>
#include <bitset>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>

template <typename T>
inline T approx_log2( const T& x ) {
  auto vx = castToInt( x );
  auto mx = castToFloat( ( vx & 0x007FFFFF ) | 0x3F000000 );
  auto y  = T( vx ) * 1.1920928955078125e-7f;
  return y - 124.22551499f - 1.498030302f * mx - 1.72587999f / ( 0.3520887068f + mx );
}

template <typename T>
inline T approx_log( const T& x ) {
  return approx_log2( x ) * 0.69314718f;
}

namespace SIMDWrapper {
  template <typename T>
  inline T log( T const& x ) {
    return approx_log( x );
  }

  enum InstructionSet { Scalar = 0, AVX2, AVX256, AVX512, EndOfList };
  inline std::string instructionSetName( InstructionSet set ) {
    switch ( set ) {
    case Scalar:
      return "Scalar";
    case AVX2:
      return "AVX2";
    case AVX256:
      return "AVX256";
    case AVX512:
      return "AVX512";
    default:
      return "Unknown";
    }
  }

  template <typename bare, std::size_t size, typename T>
  inline std::ostream& print_vector( std::ostream& os, T const& x, const char* name ) {
    std::array<bare, size> tmp;
    x.store( tmp.data() );
    os << name << "{";
    for ( std::size_t i = 0; i < size - 1; ++i ) { os << tmp[i] << ", "; }
    return os << tmp[size - 1] << "}";
  }

  template <std::size_t size, typename T>
  inline std::ostream& print_bitset( std::ostream& os, T const& x, const char* name ) {
    auto s = std::bitset<size>{x}.to_string();
    std::reverse( s.begin(), s.end() );
    return os << name << "{" << s << "}";
  }

  // Define some helper functions that will always return which instruction set
  // was used for the types in a given namespace when the stack was built.
  namespace scalar {
    using SIMDWrapper::log;
  } // namespace scalar
  namespace avx2 {
    using SIMDWrapper::log;
    InstructionSet stackInstructionSet();
  } // namespace avx2
  namespace avx256 {
    using SIMDWrapper::log;
    InstructionSet stackInstructionSet();
  } // namespace avx256
  namespace avx512 {
    using SIMDWrapper::log;
    InstructionSet stackInstructionSet();
  } // namespace avx512
  namespace best {
    using SIMDWrapper::log;
    InstructionSet stackInstructionSet();
  } // namespace best

  namespace scalar {
    constexpr InstructionSet instructionSet() { return Scalar; }
    class mask_v {
    public:
      mask_v() {} // Constructor must be empty
      mask_v( int m ) : data( m ) {}

      mask_v& operator=( const int& m ) {
        data = m;
        return *this;
      }

      operator int() const { return data; }

      friend mask_v operator&&( const mask_v& lhs, const mask_v& rhs ) { return lhs.data & rhs.data; }
      friend mask_v operator||( const mask_v& lhs, const mask_v& rhs ) { return lhs.data | rhs.data; }
      friend mask_v operator!( const mask_v& x ) { return !x.data; }

      friend bool all( const mask_v& mask ) { return mask == 1; }
      friend bool none( const mask_v& mask ) { return mask == 0; }
      friend bool any( const mask_v& mask ) { return mask != 0; }

    private:
      int data;
    };

    inline int popcount( const mask_v& mask ) { return (int)mask; }

    class int_v;

    class float_v {
    public:
      float_v() {} // Constructor must be empty
      float_v( float f ) : data( f ) {}
      float_v( const float* f ) : data( *f ) {}

      float_v& operator=( const float& f ) {
        data = f;
        return *this;
      }

      inline operator int_v() const;

      float cast() const { return data; }

      void store( float* ptr ) const { *ptr = data; }

      void compressstore( int mask, float* ptr ) const {
        if ( mask ) *ptr = data;
      }

      float hmax() const { return data; }
      float hmin() const { return data; }
      float hadd() const { return data; }
      float hmax( float mask ) const { return ( mask ) ? data : std::numeric_limits<float>::min(); }
      float hmin( float mask ) const { return ( mask ) ? data : std::numeric_limits<float>::max(); }
      float hadd( float mask ) const { return ( mask ) ? data : 0; }

      friend float_v operator+( const float_v& lhs, const float_v& rhs ) { return ( lhs.data + rhs.data ); }
      friend float_v operator-( const float_v& lhs, const float_v& rhs ) { return ( lhs.data - rhs.data ); }
      friend float_v operator*( const float_v& lhs, const float_v& rhs ) { return ( lhs.data * rhs.data ); }
      friend float_v operator/( const float_v& lhs, const float_v& rhs ) { return ( lhs.data / rhs.data ); }
      friend float_v operator-( const float_v& x ) { return -1.f * x; }

      friend float_v abs( const float_v& v ) { return std::abs( v.data ); }
      friend float_v copysign( const float_v& x, const float_v& y ) { return std::copysign( x.data, y.data ); }

      friend float_v sqrt( const float_v& v ) { return std::sqrt( v.data ); }

      friend float_v min( const float_v& lhs, const float_v& rhs ) { return std::min( lhs, rhs ); }
      friend float_v max( const float_v& lhs, const float_v& rhs ) { return std::max( lhs, rhs ); }

      friend float_v signselect( const float_v& s, const float_v& a, const float_v& b ) {
        return ( s.data > 0 ) ? a : b;
      }
      friend float_v select( int mask, const float_v& a, const float_v& b ) { return ( mask ) ? a : b; }

      friend mask_v        operator<( const float_v& lhs, const float_v& rhs ) { return lhs.data < rhs.data; }
      friend mask_v        operator>( const float_v& lhs, const float_v& rhs ) { return lhs.data > rhs.data; }
      friend mask_v        operator==( const float_v& lhs, const float_v& rhs ) { return lhs.data == rhs.data; }
      friend std::ostream& operator<<( std::ostream& os, float_v const& x ) { return os << x.cast(); }

    private:
      float data;
    };

    class int_v {
    public:
      int_v() {} // Constructor must be empty
      int_v( int f ) : data( f ) {}
      int_v( const int* f ) : data( *f ) {}

      int_v& operator=( const int& f ) {
        data = f;
        return *this;
      }

          operator float_v() const { return float_v( float( data ) ); }
      int cast() const { return data; } // don't cast directly to avoid conflict

      void store( int* ptr ) const { *ptr = data; }

      void compressstore( int mask, int* ptr ) const {
        if ( mask ) *ptr = data;
      }

      int hmax() const { return data; }
      int hmin() const { return data; }
      int hadd() const { return data; }
      int hmax( int mask ) const { return ( mask ) ? data : std::numeric_limits<int>::min(); }
      int hmin( int mask ) const { return ( mask ) ? data : std::numeric_limits<int>::max(); }
      int hadd( int mask ) const { return ( mask ) ? data : 0; }

      friend int_v operator+( const int_v& lhs, const int_v& rhs ) { return lhs.data + rhs.data; }
      friend int_v operator-( const int_v& lhs, const int_v& rhs ) { return lhs.data - rhs.data; }
      friend int_v operator*( const int_v& lhs, const int_v& rhs ) { return lhs.data * rhs.data; }

      friend int_v operator&( const int_v& lhs, const int_v& rhs ) { return lhs.data & rhs.data; }
      friend int_v operator|( const int_v& lhs, const int_v& rhs ) { return lhs.data | rhs.data; }

      friend int_v operator<<( const int_v& lhs, const int_v& rhs ) { return lhs.data << rhs.data; }
      friend int_v operator>>( const int_v& lhs, const int_v& rhs ) { return lhs.data >> rhs.data; }

      friend int_v signselect( const float_v& s, const int_v& a, const int_v& b ) {
        return ( s > float_v( 0.f ) ) ? a : b;
      }
      friend int_v select( int mask, const int_v& a, const int_v& b ) { return ( mask ) ? a : b; }

      friend mask_v operator<( const int_v& lhs, const int_v& rhs ) { return lhs.data < rhs.data; }
      friend mask_v operator>( const int_v& lhs, const int_v& rhs ) { return lhs.data > rhs.data; }
      friend mask_v operator==( const int_v& lhs, const int_v& rhs ) { return lhs.data == rhs.data; }

      friend std::ostream& operator<<( std::ostream& os, int_v const& x ) { return os << x.cast(); }

    private:
      int data;
    };

    inline float_v::operator int_v() const { return int_v( int( data ) ); }

    inline int_v castToInt( const float_v& x ) { return bit_cast<int>( x.cast() ); }

    inline float_v castToFloat( const int_v& x ) { return bit_cast<float>( x.cast() ); }

    inline int_v gather( const int* base, const int_v& idx ) { return base[idx.cast()]; }

    inline float_v gather( const float* base, const int_v& idx ) { return base[idx.cast()]; }

    inline int_v maskgather( const int* base, const int_v& idx, const mask_v& mask, const int_v& source ) {
      return select( mask, base[idx.cast()], source );
    }

    inline float_v maskgather( const float* base, const int_v& idx, const mask_v& mask, const float_v& source ) {
      return select( mask, base[idx.cast()], source );
    }

    struct types {
      static const size_t size = 1;
      using int_v              = scalar::int_v;
      using float_v            = scalar::float_v;
      using mask_v             = scalar::mask_v;
      static mask_v mask_true() { return true; }
      static mask_v mask_false() { return false; }
      static int_v  indices() { return 0; }
      static int_v  indices( int start ) { return start; }
      static int    popcount( mask_v const& mask ) { return scalar::popcount( mask ); }
      static mask_v loop_mask( int, int ) { return true; }
    };
  } // namespace scalar
} // namespace SIMDWrapper

template <>
struct LHCb::type_map<SIMDWrapper::scalar::float_v> {
  using int_t = SIMDWrapper::scalar::int_v;
};

namespace SIMDWrapper {
#ifndef __AVX2__
  namespace avx2 {
    constexpr InstructionSet instructionSet() { return scalar::instructionSet(); }
    using float_v = scalar::float_v;
    using int_v   = scalar::int_v;
    using types   = scalar::types;
  } // namespace avx2
#else
  namespace avx2 {
    constexpr InstructionSet instructionSet() { return AVX2; }
    // Permutation for avx2 compress, from https://github.com/lemire/simdprune
    alignas( 32 ) const uint32_t compress_mask256_epi32[] = {
        0, 1, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 7, 0, 2, 3, 4, 5, 6, 7, 7, 2, 3, 4, 5, 6, 7, 7, 7, 0, 1, 3, 4, 5,
        6, 7, 7, 1, 3, 4, 5, 6, 7, 7, 7, 0, 3, 4, 5, 6, 7, 7, 7, 3, 4, 5, 6, 7, 7, 7, 7, 0, 1, 2, 4, 5, 6, 7, 7, 1, 2,
        4, 5, 6, 7, 7, 7, 0, 2, 4, 5, 6, 7, 7, 7, 2, 4, 5, 6, 7, 7, 7, 7, 0, 1, 4, 5, 6, 7, 7, 7, 1, 4, 5, 6, 7, 7, 7,
        7, 0, 4, 5, 6, 7, 7, 7, 7, 4, 5, 6, 7, 7, 7, 7, 7, 0, 1, 2, 3, 5, 6, 7, 7, 1, 2, 3, 5, 6, 7, 7, 7, 0, 2, 3, 5,
        6, 7, 7, 7, 2, 3, 5, 6, 7, 7, 7, 7, 0, 1, 3, 5, 6, 7, 7, 7, 1, 3, 5, 6, 7, 7, 7, 7, 0, 3, 5, 6, 7, 7, 7, 7, 3,
        5, 6, 7, 7, 7, 7, 7, 0, 1, 2, 5, 6, 7, 7, 7, 1, 2, 5, 6, 7, 7, 7, 7, 0, 2, 5, 6, 7, 7, 7, 7, 2, 5, 6, 7, 7, 7,
        7, 7, 0, 1, 5, 6, 7, 7, 7, 7, 1, 5, 6, 7, 7, 7, 7, 7, 0, 5, 6, 7, 7, 7, 7, 7, 5, 6, 7, 7, 7, 7, 7, 7, 0, 1, 2,
        3, 4, 6, 7, 7, 1, 2, 3, 4, 6, 7, 7, 7, 0, 2, 3, 4, 6, 7, 7, 7, 2, 3, 4, 6, 7, 7, 7, 7, 0, 1, 3, 4, 6, 7, 7, 7,
        1, 3, 4, 6, 7, 7, 7, 7, 0, 3, 4, 6, 7, 7, 7, 7, 3, 4, 6, 7, 7, 7, 7, 7, 0, 1, 2, 4, 6, 7, 7, 7, 1, 2, 4, 6, 7,
        7, 7, 7, 0, 2, 4, 6, 7, 7, 7, 7, 2, 4, 6, 7, 7, 7, 7, 7, 0, 1, 4, 6, 7, 7, 7, 7, 1, 4, 6, 7, 7, 7, 7, 7, 0, 4,
        6, 7, 7, 7, 7, 7, 4, 6, 7, 7, 7, 7, 7, 7, 0, 1, 2, 3, 6, 7, 7, 7, 1, 2, 3, 6, 7, 7, 7, 7, 0, 2, 3, 6, 7, 7, 7,
        7, 2, 3, 6, 7, 7, 7, 7, 7, 0, 1, 3, 6, 7, 7, 7, 7, 1, 3, 6, 7, 7, 7, 7, 7, 0, 3, 6, 7, 7, 7, 7, 7, 3, 6, 7, 7,
        7, 7, 7, 7, 0, 1, 2, 6, 7, 7, 7, 7, 1, 2, 6, 7, 7, 7, 7, 7, 0, 2, 6, 7, 7, 7, 7, 7, 2, 6, 7, 7, 7, 7, 7, 7, 0,
        1, 6, 7, 7, 7, 7, 7, 1, 6, 7, 7, 7, 7, 7, 7, 0, 6, 7, 7, 7, 7, 7, 7, 6, 7, 7, 7, 7, 7, 7, 7, 0, 1, 2, 3, 4, 5,
        7, 7, 1, 2, 3, 4, 5, 7, 7, 7, 0, 2, 3, 4, 5, 7, 7, 7, 2, 3, 4, 5, 7, 7, 7, 7, 0, 1, 3, 4, 5, 7, 7, 7, 1, 3, 4,
        5, 7, 7, 7, 7, 0, 3, 4, 5, 7, 7, 7, 7, 3, 4, 5, 7, 7, 7, 7, 7, 0, 1, 2, 4, 5, 7, 7, 7, 1, 2, 4, 5, 7, 7, 7, 7,
        0, 2, 4, 5, 7, 7, 7, 7, 2, 4, 5, 7, 7, 7, 7, 7, 0, 1, 4, 5, 7, 7, 7, 7, 1, 4, 5, 7, 7, 7, 7, 7, 0, 4, 5, 7, 7,
        7, 7, 7, 4, 5, 7, 7, 7, 7, 7, 7, 0, 1, 2, 3, 5, 7, 7, 7, 1, 2, 3, 5, 7, 7, 7, 7, 0, 2, 3, 5, 7, 7, 7, 7, 2, 3,
        5, 7, 7, 7, 7, 7, 0, 1, 3, 5, 7, 7, 7, 7, 1, 3, 5, 7, 7, 7, 7, 7, 0, 3, 5, 7, 7, 7, 7, 7, 3, 5, 7, 7, 7, 7, 7,
        7, 0, 1, 2, 5, 7, 7, 7, 7, 1, 2, 5, 7, 7, 7, 7, 7, 0, 2, 5, 7, 7, 7, 7, 7, 2, 5, 7, 7, 7, 7, 7, 7, 0, 1, 5, 7,
        7, 7, 7, 7, 1, 5, 7, 7, 7, 7, 7, 7, 0, 5, 7, 7, 7, 7, 7, 7, 5, 7, 7, 7, 7, 7, 7, 7, 0, 1, 2, 3, 4, 7, 7, 7, 1,
        2, 3, 4, 7, 7, 7, 7, 0, 2, 3, 4, 7, 7, 7, 7, 2, 3, 4, 7, 7, 7, 7, 7, 0, 1, 3, 4, 7, 7, 7, 7, 1, 3, 4, 7, 7, 7,
        7, 7, 0, 3, 4, 7, 7, 7, 7, 7, 3, 4, 7, 7, 7, 7, 7, 7, 0, 1, 2, 4, 7, 7, 7, 7, 1, 2, 4, 7, 7, 7, 7, 7, 0, 2, 4,
        7, 7, 7, 7, 7, 2, 4, 7, 7, 7, 7, 7, 7, 0, 1, 4, 7, 7, 7, 7, 7, 1, 4, 7, 7, 7, 7, 7, 7, 0, 4, 7, 7, 7, 7, 7, 7,
        4, 7, 7, 7, 7, 7, 7, 7, 0, 1, 2, 3, 7, 7, 7, 7, 1, 2, 3, 7, 7, 7, 7, 7, 0, 2, 3, 7, 7, 7, 7, 7, 2, 3, 7, 7, 7,
        7, 7, 7, 0, 1, 3, 7, 7, 7, 7, 7, 1, 3, 7, 7, 7, 7, 7, 7, 0, 3, 7, 7, 7, 7, 7, 7, 3, 7, 7, 7, 7, 7, 7, 7, 0, 1,
        2, 7, 7, 7, 7, 7, 1, 2, 7, 7, 7, 7, 7, 7, 0, 2, 7, 7, 7, 7, 7, 7, 2, 7, 7, 7, 7, 7, 7, 7, 0, 1, 7, 7, 7, 7, 7,
        7, 1, 7, 7, 7, 7, 7, 7, 7, 0, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 0, 1, 2, 3, 4, 5, 6, 6, 1, 2, 3, 4,
        5, 6, 6, 6, 0, 2, 3, 4, 5, 6, 6, 6, 2, 3, 4, 5, 6, 6, 6, 6, 0, 1, 3, 4, 5, 6, 6, 6, 1, 3, 4, 5, 6, 6, 6, 6, 0,
        3, 4, 5, 6, 6, 6, 6, 3, 4, 5, 6, 6, 6, 6, 6, 0, 1, 2, 4, 5, 6, 6, 6, 1, 2, 4, 5, 6, 6, 6, 6, 0, 2, 4, 5, 6, 6,
        6, 6, 2, 4, 5, 6, 6, 6, 6, 6, 0, 1, 4, 5, 6, 6, 6, 6, 1, 4, 5, 6, 6, 6, 6, 6, 0, 4, 5, 6, 6, 6, 6, 6, 4, 5, 6,
        6, 6, 6, 6, 6, 0, 1, 2, 3, 5, 6, 6, 6, 1, 2, 3, 5, 6, 6, 6, 6, 0, 2, 3, 5, 6, 6, 6, 6, 2, 3, 5, 6, 6, 6, 6, 6,
        0, 1, 3, 5, 6, 6, 6, 6, 1, 3, 5, 6, 6, 6, 6, 6, 0, 3, 5, 6, 6, 6, 6, 6, 3, 5, 6, 6, 6, 6, 6, 6, 0, 1, 2, 5, 6,
        6, 6, 6, 1, 2, 5, 6, 6, 6, 6, 6, 0, 2, 5, 6, 6, 6, 6, 6, 2, 5, 6, 6, 6, 6, 6, 6, 0, 1, 5, 6, 6, 6, 6, 6, 1, 5,
        6, 6, 6, 6, 6, 6, 0, 5, 6, 6, 6, 6, 6, 6, 5, 6, 6, 6, 6, 6, 6, 6, 0, 1, 2, 3, 4, 6, 6, 6, 1, 2, 3, 4, 6, 6, 6,
        6, 0, 2, 3, 4, 6, 6, 6, 6, 2, 3, 4, 6, 6, 6, 6, 6, 0, 1, 3, 4, 6, 6, 6, 6, 1, 3, 4, 6, 6, 6, 6, 6, 0, 3, 4, 6,
        6, 6, 6, 6, 3, 4, 6, 6, 6, 6, 6, 6, 0, 1, 2, 4, 6, 6, 6, 6, 1, 2, 4, 6, 6, 6, 6, 6, 0, 2, 4, 6, 6, 6, 6, 6, 2,
        4, 6, 6, 6, 6, 6, 6, 0, 1, 4, 6, 6, 6, 6, 6, 1, 4, 6, 6, 6, 6, 6, 6, 0, 4, 6, 6, 6, 6, 6, 6, 4, 6, 6, 6, 6, 6,
        6, 6, 0, 1, 2, 3, 6, 6, 6, 6, 1, 2, 3, 6, 6, 6, 6, 6, 0, 2, 3, 6, 6, 6, 6, 6, 2, 3, 6, 6, 6, 6, 6, 6, 0, 1, 3,
        6, 6, 6, 6, 6, 1, 3, 6, 6, 6, 6, 6, 6, 0, 3, 6, 6, 6, 6, 6, 6, 3, 6, 6, 6, 6, 6, 6, 6, 0, 1, 2, 6, 6, 6, 6, 6,
        1, 2, 6, 6, 6, 6, 6, 6, 0, 2, 6, 6, 6, 6, 6, 6, 2, 6, 6, 6, 6, 6, 6, 6, 0, 1, 6, 6, 6, 6, 6, 6, 1, 6, 6, 6, 6,
        6, 6, 6, 0, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 0, 1, 2, 3, 4, 5, 5, 5, 1, 2, 3, 4, 5, 5, 5, 5, 0, 2,
        3, 4, 5, 5, 5, 5, 2, 3, 4, 5, 5, 5, 5, 5, 0, 1, 3, 4, 5, 5, 5, 5, 1, 3, 4, 5, 5, 5, 5, 5, 0, 3, 4, 5, 5, 5, 5,
        5, 3, 4, 5, 5, 5, 5, 5, 5, 0, 1, 2, 4, 5, 5, 5, 5, 1, 2, 4, 5, 5, 5, 5, 5, 0, 2, 4, 5, 5, 5, 5, 5, 2, 4, 5, 5,
        5, 5, 5, 5, 0, 1, 4, 5, 5, 5, 5, 5, 1, 4, 5, 5, 5, 5, 5, 5, 0, 4, 5, 5, 5, 5, 5, 5, 4, 5, 5, 5, 5, 5, 5, 5, 0,
        1, 2, 3, 5, 5, 5, 5, 1, 2, 3, 5, 5, 5, 5, 5, 0, 2, 3, 5, 5, 5, 5, 5, 2, 3, 5, 5, 5, 5, 5, 5, 0, 1, 3, 5, 5, 5,
        5, 5, 1, 3, 5, 5, 5, 5, 5, 5, 0, 3, 5, 5, 5, 5, 5, 5, 3, 5, 5, 5, 5, 5, 5, 5, 0, 1, 2, 5, 5, 5, 5, 5, 1, 2, 5,
        5, 5, 5, 5, 5, 0, 2, 5, 5, 5, 5, 5, 5, 2, 5, 5, 5, 5, 5, 5, 5, 0, 1, 5, 5, 5, 5, 5, 5, 1, 5, 5, 5, 5, 5, 5, 5,
        0, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 0, 1, 2, 3, 4, 4, 4, 4, 1, 2, 3, 4, 4, 4, 4, 4, 0, 2, 3, 4, 4,
        4, 4, 4, 2, 3, 4, 4, 4, 4, 4, 4, 0, 1, 3, 4, 4, 4, 4, 4, 1, 3, 4, 4, 4, 4, 4, 4, 0, 3, 4, 4, 4, 4, 4, 4, 3, 4,
        4, 4, 4, 4, 4, 4, 0, 1, 2, 4, 4, 4, 4, 4, 1, 2, 4, 4, 4, 4, 4, 4, 0, 2, 4, 4, 4, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4,
        4, 0, 1, 4, 4, 4, 4, 4, 4, 1, 4, 4, 4, 4, 4, 4, 4, 0, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0, 1, 2, 3,
        3, 3, 3, 3, 1, 2, 3, 3, 3, 3, 3, 3, 0, 2, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 3, 0, 1, 3, 3, 3, 3, 3, 3, 1,
        3, 3, 3, 3, 3, 3, 3, 0, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 0, 1, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 2, 2,
        2, 2, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    class int_v;

    class float_v {
    public:
      float_v() {} // Constructor must be empty
      float_v( scalar::float_v& f ) : data( _mm256_set1_ps( f.cast() ) ) {}
      float_v( float f ) : data( _mm256_set1_ps( f ) ) {}
      float_v( const float* f ) : data( _mm256_loadu_ps( f ) ) {}
      float_v( __m256 f ) : data( f ) {}

      float_v& operator=( const __m256& f ) {
        data = f;
        return *this;
      }

      float_v& operator=( const float_v& f ) {
        data = f.data;
        return *this;
      }

      operator __m256() const { return data; }

      inline operator int_v() const;

      void store( float* ptr ) const { _mm256_storeu_ps( ptr, data ); }

      void compressstore( __m256 mask, float* ptr ) const {
        __m256i perm =
            _mm256_load_si256( (const __m256i*)compress_mask256_epi32 + ( _mm256_movemask_ps( mask ) ^ 0xFF ) );
        _mm256_storeu_ps( ptr, _mm256_permutevar8x32_ps( data, perm ) );
      }

      friend float_v operator+( const float_v& lhs, const float_v& rhs ) { return _mm256_add_ps( lhs, rhs ); }
      friend float_v operator-( const float_v& lhs, const float_v& rhs ) { return _mm256_sub_ps( lhs, rhs ); }
      friend float_v operator*( const float_v& lhs, const float_v& rhs ) { return _mm256_mul_ps( lhs, rhs ); }
      friend float_v operator/( const float_v& lhs, const float_v& rhs ) { return _mm256_div_ps( lhs, rhs ); }
      friend float_v operator-( const float_v& x ) { return -1.f * x; }

      friend float_v operator&( const float_v& lhs, const float_v& rhs ) { return _mm256_and_ps( lhs, rhs ); }
      friend float_v operator|( const float_v& lhs, const float_v& rhs ) { return _mm256_or_ps( lhs, rhs ); }
      friend float_v operator^( const float_v& lhs, const float_v& rhs ) { return _mm256_xor_ps( lhs, rhs ); }

      friend float_v operator&&( const float_v& lhs, const float_v& rhs ) { return _mm256_and_ps( lhs, rhs ); }
      friend float_v operator||( const float_v& lhs, const float_v& rhs ) { return _mm256_or_ps( lhs, rhs ); }
      friend float_v operator!( const float_v& x ) { return x ^ _mm256_castsi256_ps( _mm256_set1_epi32( -1 ) ); }

      friend std::ostream& operator<<( std::ostream& os, float_v const& x ) {
        return print_vector<float, 8>( os, x, "avx2" );
      }

      friend float_v min( const float_v& lhs, const float_v& rhs ) { return _mm256_min_ps( lhs, rhs ); }
      friend float_v max( const float_v& lhs, const float_v& rhs ) { return _mm256_max_ps( lhs, rhs ); }
      friend float_v abs( const float_v& v ) { return v & _mm256_castsi256_ps( _mm256_set1_epi32( 0x7FFFFFFF ) ); }
      friend float_v copysign( const float_v& x, const float_v& y ) {
        return x ^ ( y & _mm256_castsi256_ps( _mm256_set1_epi32( 0x80000000 ) ) );
      }

      friend float_v signselect( const float_v& s, const float_v& a, const float_v& b ) {
        return _mm256_blendv_ps( a, b, s );
      }
      friend float_v select( __m256 mask, const float_v& a, const float_v& b ) {
        return _mm256_blendv_ps( b, a, mask );
      }

      friend float_v sqrt( const float_v& v ) { return _mm256_sqrt_ps( v ); }

      friend float_v operator<( const float_v& lhs, const float_v& rhs ) {
        return _mm256_cmp_ps( lhs, rhs, _CMP_LT_OS );
      }
      friend float_v operator>( const float_v& lhs, const float_v& rhs ) {
        return _mm256_cmp_ps( lhs, rhs, _CMP_GT_OS );
      }
      friend float_v operator==( const float_v& lhs, const float_v& rhs ) {
        return _mm256_cmp_ps( lhs, rhs, _CMP_EQ_OS );
      }

      friend bool all( const float_v& mask ) { return _mm256_movemask_ps( mask ) == 0xFF; }
      friend bool none( const float_v& mask ) { return _mm256_movemask_ps( mask ) == 0x00; }
      friend bool any( const float_v& mask ) { return _mm256_movemask_ps( mask ) != 0x00; }

      float hmax() const {
        __m128 r = _mm_max_ps( _mm256_extractf128_ps( data, 0 ), _mm256_extractf128_ps( data, 1 ) );
        r        = _mm_max_ps( r, _mm_shuffle_ps( r, r, _MM_SHUFFLE( 2, 3, 0, 1 ) ) );
        r        = _mm_max_ps( r, _mm_shuffle_ps( r, r, _MM_SHUFFLE( 1, 0, 3, 2 ) ) );
        return _mm_cvtss_f32( r );
      }
      float hmin() const {
        __m128 r = _mm_min_ps( _mm256_extractf128_ps( data, 0 ), _mm256_extractf128_ps( data, 1 ) );
        r        = _mm_min_ps( r, _mm_shuffle_ps( r, r, _MM_SHUFFLE( 2, 3, 0, 1 ) ) );
        r        = _mm_min_ps( r, _mm_shuffle_ps( r, r, _MM_SHUFFLE( 1, 0, 3, 2 ) ) );
        return _mm_cvtss_f32( r );
      }
      float hadd() const {
        __m128 r = _mm_add_ps( _mm256_extractf128_ps( data, 0 ), _mm256_extractf128_ps( data, 1 ) );
        r        = _mm_add_ps( r, _mm_shuffle_ps( r, r, _MM_SHUFFLE( 2, 3, 0, 1 ) ) );
        r        = _mm_add_ps( r, _mm_shuffle_ps( r, r, _MM_SHUFFLE( 1, 0, 3, 2 ) ) );
        return _mm_cvtss_f32( r );
      }
      float hmax( const __m256 mask ) const { return select( mask, *this, std::numeric_limits<float>::min() ).hmax(); }
      float hmin( const __m256 mask ) const { return select( mask, *this, std::numeric_limits<float>::max() ).hmin(); }
      float hadd( const __m256 mask ) const { return select( mask, *this, 0.f ).hadd(); }

    private:
      __m256 data;
    };

    inline int popcount( const float_v& mask ) { return _mm_popcnt_u32( _mm256_movemask_ps( mask ) ); }

    class int_v {
    public:
      int_v() {} // Constructor must be empty
      int_v( int f ) : data( _mm256_set1_epi32( f ) ) {}
      int_v( const int* f ) : data( _mm256_loadu_si256( (__m256i*)f ) ) {}
      int_v( __m256i f ) : data( f ) {}

      int_v& operator=( const __m256i& f ) {
        data = f;
        return *this;
      }

      operator __m256i() const { return data; }
      operator float_v() const { return float_v( _mm256_cvtepi32_ps( data ) ); }

      void store( int* ptr ) const { _mm256_storeu_si256( (__m256i*)ptr, data ); }

      void compressstore( __m256 mask, int* ptr ) const {
        __m256i perm =
            _mm256_load_si256( (const __m256i*)compress_mask256_epi32 + ( _mm256_movemask_ps( mask ) ^ 0xFF ) );
        _mm256_storeu_si256( (__m256i*)ptr, _mm256_permutevar8x32_epi32( data, perm ) );
      }

      int hmax() const {
        __m128i r = _mm_max_epi32( _mm256_extractf128_si256( data, 0 ), _mm256_extractf128_si256( data, 1 ) );
        r         = _mm_max_epi32( r, _mm_shuffle_epi32( r, _MM_SHUFFLE( 2, 3, 0, 1 ) ) );
        r         = _mm_max_epi32( r, _mm_shuffle_epi32( r, _MM_SHUFFLE( 1, 0, 3, 2 ) ) );
        return _mm_extract_epi32( r, 0 );
      }
      int hmin() const {
        __m128i r = _mm_min_epi32( _mm256_extractf128_si256( data, 0 ), _mm256_extractf128_si256( data, 1 ) );
        r         = _mm_min_epi32( r, _mm_shuffle_epi32( r, _MM_SHUFFLE( 2, 3, 0, 1 ) ) );
        r         = _mm_min_epi32( r, _mm_shuffle_epi32( r, _MM_SHUFFLE( 1, 0, 3, 2 ) ) );
        return _mm_extract_epi32( r, 0 );
      }
      int hadd() const {
        __m128i r = _mm_add_epi32( _mm256_extractf128_si256( data, 0 ), _mm256_extractf128_si256( data, 1 ) );
        r         = _mm_add_epi32( r, _mm_shuffle_epi32( r, _MM_SHUFFLE( 2, 3, 0, 1 ) ) );
        r         = _mm_add_epi32( r, _mm_shuffle_epi32( r, _MM_SHUFFLE( 1, 0, 3, 2 ) ) );
        return _mm_extract_epi32( r, 0 );
      }
      int hmax( const __m256 mask ) const { return select( mask, *this, std::numeric_limits<int>::min() ).hmax(); }
      int hmin( const __m256 mask ) const { return select( mask, *this, std::numeric_limits<int>::max() ).hmin(); }
      int hadd( const __m256 mask ) const { return select( mask, *this, 0 ).hadd(); }

      friend int_v operator+( const int_v& lhs, const int_v& rhs ) { return _mm256_add_epi32( lhs, rhs ); }
      friend int_v operator-( const int_v& lhs, const int_v& rhs ) { return _mm256_sub_epi32( lhs, rhs ); }
      friend int_v operator*( const int_v& lhs, const int_v& rhs ) { return _mm256_mul_epi32( lhs, rhs ); }

      friend int_v operator&( const int_v& lhs, const int_v& rhs ) { return _mm256_and_si256( lhs, rhs ); }
      friend int_v operator|( const int_v& lhs, const int_v& rhs ) { return _mm256_or_si256( lhs, rhs ); }

      friend int_v operator<<( const int_v& lhs, const int_v& rhs ) { return _mm256_sllv_epi32( lhs, rhs ); }
      friend int_v operator>>( const int_v& lhs, const int_v& rhs ) { return _mm256_srlv_epi32( lhs, rhs ); }

      friend std::ostream& operator<<( std::ostream& os, int_v const& x ) {
        return print_vector<int, 8>( os, x, "avx2" );
      }

      friend int_v signselect( const float_v& s, const int_v& a, const int_v& b ) {
        return _mm256_castps_si256( _mm256_blendv_ps( _mm256_castsi256_ps( a ), _mm256_castsi256_ps( b ), s ) );
      }
      friend int_v select( const __m256 mask, const int_v& a, const int_v& b ) {
        return _mm256_castps_si256( _mm256_blendv_ps( _mm256_castsi256_ps( b ), _mm256_castsi256_ps( a ), mask ) );
      }

      friend float_v operator<( const int_v& lhs, const int_v& rhs ) {
        return _mm256_castsi256_ps( _mm256_cmpgt_epi32( rhs, lhs ) );
      }
      friend float_v operator>( const int_v& lhs, const int_v& rhs ) {
        return _mm256_castsi256_ps( _mm256_cmpgt_epi32( lhs, rhs ) );
      }
      friend float_v operator==( const int_v& lhs, const int_v& rhs ) {
        return _mm256_castsi256_ps( _mm256_cmpeq_epi32( lhs, rhs ) );
      }

    private:
      __m256i data;
    };

    inline float_v::operator int_v() const { return int_v( _mm256_cvttps_epi32( data ) ); }

    inline int_v castToInt( const float_v& x ) { return int_v( _mm256_castps_si256( x ) ); }

    inline float_v castToFloat( const int_v& x ) { return float_v( _mm256_castsi256_ps( x ) ); }

    inline int_v gather( const int* base, const int_v& idx ) {
      return _mm256_i32gather_epi32( base, idx, sizeof( int ) );
    }

    inline float_v gather( const float* base, const int_v& idx ) {
      return _mm256_i32gather_ps( base, idx, sizeof( float ) );
    }

    inline int_v maskgather( const int* base, const int_v& idx, const float_v& mask, const int_v& source ) {
      return _mm256_mask_i32gather_epi32( source, base, idx, castToInt( mask ), sizeof( int ) );
    }

    inline float_v maskgather( const float* base, const int_v& idx, const float_v& mask, const float_v& source ) {
      return _mm256_mask_i32gather_ps( source, base, idx, mask, sizeof( float ) );
    }

    struct types {
      static const size_t size = 8;
      using int_v              = avx2::int_v;
      using float_v            = avx2::float_v;
      using mask_v             = avx2::float_v;
      static mask_v  mask_true() { return mask_v( -1.f ); }
      static mask_v  mask_false() { return mask_v( 0.f ); }
      static int_v   indices() { return _mm256_setr_epi32( 0, 1, 2, 3, 4, 5, 6, 7 ); }
      static int_v   indices( int start ) { return indices() + start; }
      static int     popcount( mask_v const& mask ) { return avx2::popcount( mask ); }
      static float_v loop_mask( int i, int n ) {
        return _mm256_cmp_ps( _mm256_setr_ps( 0, 1, 2, 3, 4, 5, 6, 7 ), _mm256_set1_ps( n - i ), _CMP_LT_OS );
      }
    };
  } // namespace avx2
} // namespace SIMDWrapper

template <>
struct LHCb::type_map<SIMDWrapper::avx2::float_v> {
  using int_t = SIMDWrapper::avx2::int_v;
};

namespace SIMDWrapper {
#endif

#ifndef __AVX512F__
  namespace avx512 {
    constexpr InstructionSet instructionSet() { return avx2::instructionSet(); }
    using float_v = avx2::float_v;
    using int_v   = avx2::int_v;
    using types   = avx2::types;
  } // namespace avx512

  namespace avx256 {
    constexpr InstructionSet instructionSet() { return avx2::instructionSet(); }
    using float_v = avx2::float_v;
    using int_v   = avx2::int_v;
    using types   = avx2::types;
  } // namespace avx256
#else
  namespace avx256 {
    constexpr InstructionSet instructionSet() { return AVX256; }
    class mask_v {
    public:
      mask_v() {} // Constructor must be empty
      mask_v( __mmask8 m ) : data( m ) {}

      mask_v& operator=( const __mmask8& m ) {
        data = m;
        return *this;
      }

      operator __mmask8() const { return data; }

      friend mask_v operator&&( const mask_v& lhs, const mask_v& rhs ) { return lhs & rhs; }
      friend mask_v operator||( const mask_v& lhs, const mask_v& rhs ) { return lhs | rhs; }
      friend mask_v operator!( const mask_v& x ) { return ~x; }

      friend bool all( const mask_v& mask ) { return mask == 0xFF; }
      friend bool none( const mask_v& mask ) { return mask == 0x00; }
      friend bool any( const mask_v& mask ) { return mask != 0x00; }

      friend std::ostream& operator<<( std::ostream& os, mask_v const& x ) {
        return print_bitset<8>( os, x, "avx256" );
      }

    private:
      __mmask8 data;
    };

    inline int popcount( mask_v const& mask ) { return _mm_popcnt_u32( mask ); }

    class int_v;

    class float_v {
    public:
      float_v() {} // Constructor must be empty
      float_v( scalar::float_v& f ) : data( _mm256_set1_ps( f.cast() ) ) {}
      float_v( float f ) : data( _mm256_set1_ps( f ) ) {}
      float_v( const float* f ) : data( _mm256_loadu_ps( f ) ) {}
      float_v( __m256 f ) : data( f ) {}

      float_v& operator=( const __m256& f ) {
        data = f;
        return *this;
      }

      operator __m256() const { return data; }
      inline operator int_v() const;

      void store( float* ptr ) const { _mm256_storeu_ps( ptr, data ); }

      void compressstore( const mask_v& mask, float* ptr ) const { _mm256_mask_compressstoreu_ps( ptr, mask, data ); }

      float hmax() const {
        __m128 r = _mm_max_ps( _mm256_extractf128_ps( data, 0 ), _mm256_extractf128_ps( data, 1 ) );
        r = _mm_max_ps( r, _mm_shuffle_ps( r, r, _MM_SHUFFLE( 2, 3, 0, 1 ) ) );
        r = _mm_max_ps( r, _mm_shuffle_ps( r, r, _MM_SHUFFLE( 1, 0, 3, 2 ) ) );
        return _mm_cvtss_f32( r );
      }
      float hmin() const {
        __m128 r = _mm_min_ps( _mm256_extractf128_ps( data, 0 ), _mm256_extractf128_ps( data, 1 ) );
        r = _mm_min_ps( r, _mm_shuffle_ps( r, r, _MM_SHUFFLE( 2, 3, 0, 1 ) ) );
        r = _mm_min_ps( r, _mm_shuffle_ps( r, r, _MM_SHUFFLE( 1, 0, 3, 2 ) ) );
        return _mm_cvtss_f32( r );
      }
      float hadd() const {
        __m128 r = _mm_add_ps( _mm256_extractf128_ps( data, 0 ), _mm256_extractf128_ps( data, 1 ) );
        r = _mm_add_ps( r, _mm_shuffle_ps( r, r, _MM_SHUFFLE( 2, 3, 0, 1 ) ) );
        r = _mm_add_ps( r, _mm_shuffle_ps( r, r, _MM_SHUFFLE( 1, 0, 3, 2 ) ) );
        return _mm_cvtss_f32( r );
      }
      float hmax( const mask_v& mask ) const { return select( mask, *this, std::numeric_limits<float>::min() ).hmax(); }
      float hmin( const mask_v& mask ) const { return select( mask, *this, std::numeric_limits<float>::max() ).hmin(); }
      float hadd( const mask_v& mask ) const { return select( mask, *this, 0.f ).hadd(); }

      friend float_v operator+( const float_v& lhs, const float_v& rhs ) { return _mm256_add_ps( lhs, rhs ); }
      friend float_v operator-( const float_v& lhs, const float_v& rhs ) { return _mm256_sub_ps( lhs, rhs ); }
      friend float_v operator*( const float_v& lhs, const float_v& rhs ) { return _mm256_mul_ps( lhs, rhs ); }
      friend float_v operator/( const float_v& lhs, const float_v& rhs ) { return _mm256_div_ps( lhs, rhs ); }
      friend float_v operator-( const float_v& x ) { return -1.f * x; }

      friend float_v operator&( const float_v& lhs, const float_v& rhs ) { return _mm256_and_ps( lhs, rhs ); }
      friend float_v operator|( const float_v& lhs, const float_v& rhs ) { return _mm256_or_ps( lhs, rhs ); }
      friend float_v operator^( const float_v& lhs, const float_v& rhs ) { return _mm256_xor_ps( lhs, rhs ); }

      friend float_v operator&&( const float_v& lhs, const float_v& rhs ) { return _mm256_and_ps( lhs, rhs ); }
      friend float_v operator||( const float_v& lhs, const float_v& rhs ) { return _mm256_or_ps( lhs, rhs ); }
      friend float_v operator!( const float_v& x ) { return x ^ _mm256_castsi256_ps( _mm256_set1_epi32( -1 ) ); }

      friend float_v min( const float_v& lhs, const float_v& rhs ) { return _mm256_min_ps( lhs, rhs ); }
      friend float_v max( const float_v& lhs, const float_v& rhs ) { return _mm256_max_ps( lhs, rhs ); }
      friend float_v abs( const float_v& v ) { return v & _mm256_castsi256_ps( _mm256_set1_epi32( 0x7FFFFFFF ) ); }
      friend float_v copysign( const float_v& x, const float_v& y ) {
        return x ^ ( y & _mm256_castsi256_ps( _mm256_set1_epi32( 0x80000000 ) ) );
      }

      friend float_v signselect( const float_v& s, const float_v& a, const float_v& b ) {
        return _mm256_mask_mov_ps( a, s < float_v( 0.f ), b );
      }
      friend float_v select( const mask_v& mask, const float_v& a, const float_v& b ) {
        return _mm256_mask_mov_ps( b, mask, a );
      }

      friend float_v sqrt( const float_v& v ) { return _mm256_sqrt_ps( v ); }

      friend mask_v operator<( const float_v& lhs, const float_v& rhs ) {
        return _mm256_cmp_ps_mask( lhs, rhs, _CMP_LT_OS );
      }
      friend mask_v operator>( const float_v& lhs, const float_v& rhs ) {
        return _mm256_cmp_ps_mask( lhs, rhs, _CMP_GT_OS );
      }
      friend mask_v operator==( const float_v& lhs, const float_v& rhs ) {
        return _mm256_cmp_ps_mask( lhs, rhs, _CMP_EQ_OS );
      }
      friend std::ostream& operator<<( std::ostream& os, float_v const& x ) {
        return print_vector<float, 8>( os, x, "avx256" );
      }

    private:
      __m256 data;
    };

    class int_v {
    public:
      int_v() {} // Constructor must be empty
      int_v( int f ) : data( _mm256_set1_epi32( f ) ) {}
      int_v( const int* f ) : data( _mm256_loadu_si256( (__m256i*)f ) ) {}
      constexpr int_v( __m256i f ) : data( f ) {}

      int_v& operator=( const __m256i& f ) {
        data = f;
        return *this;
      }

      operator __m256i() const { return data; }
      operator float_v() const { return float_v( _mm256_cvtepi32_ps( data ) ); }

      void store( int* ptr ) const { _mm256_storeu_si256( (__m256i*)ptr, data ); }

      void compressstore( const mask_v& mask, int* ptr ) const { _mm256_mask_compressstoreu_epi32( ptr, mask, data ); }

      int hmax() const {
        __m128i r = _mm_max_epi32( _mm256_extractf128_si256( data, 0 ), _mm256_extractf128_si256( data, 1 ) );
        r = _mm_max_epi32( r, _mm_shuffle_epi32( r, _MM_SHUFFLE( 2, 3, 0, 1 ) ) );
        r = _mm_max_epi32( r, _mm_shuffle_epi32( r, _MM_SHUFFLE( 1, 0, 3, 2 ) ) );
        return _mm_extract_epi32( r, 0 );
      }
      int hmin() const {
        __m128i r = _mm_min_epi32( _mm256_extractf128_si256( data, 0 ), _mm256_extractf128_si256( data, 1 ) );
        r = _mm_min_epi32( r, _mm_shuffle_epi32( r, _MM_SHUFFLE( 2, 3, 0, 1 ) ) );
        r = _mm_min_epi32( r, _mm_shuffle_epi32( r, _MM_SHUFFLE( 1, 0, 3, 2 ) ) );
        return _mm_extract_epi32( r, 0 );
      }
      int hadd() const {
        __m128i r = _mm_add_epi32( _mm256_extractf128_si256( data, 0 ), _mm256_extractf128_si256( data, 1 ) );
        r = _mm_add_epi32( r, _mm_shuffle_epi32( r, _MM_SHUFFLE( 2, 3, 0, 1 ) ) );
        r = _mm_add_epi32( r, _mm_shuffle_epi32( r, _MM_SHUFFLE( 1, 0, 3, 2 ) ) );
        return _mm_extract_epi32( r, 0 );
      }
      int hmax( const mask_v& mask ) const { return select( mask, *this, std::numeric_limits<int>::min() ).hmax(); }
      int hmin( const mask_v& mask ) const { return select( mask, *this, std::numeric_limits<int>::max() ).hmin(); }
      int hadd( const mask_v& mask ) const { return select( mask, *this, 0 ).hadd(); }

      friend int_v operator+( const int_v& lhs, const int_v& rhs ) { return _mm256_add_epi32( lhs, rhs ); }
      friend int_v operator-( const int_v& lhs, const int_v& rhs ) { return _mm256_sub_epi32( lhs, rhs ); }
      friend int_v operator*( const int_v& lhs, const int_v& rhs ) { return _mm256_mul_epi32( lhs, rhs ); }

      friend int_v operator&( const int_v& lhs, const int_v& rhs ) { return _mm256_and_si256( lhs, rhs ); }
      friend int_v operator|( const int_v& lhs, const int_v& rhs ) { return _mm256_or_si256( lhs, rhs ); }

      friend int_v signselect( const float_v& s, const int_v& a, const int_v& b ) {
        return _mm256_mask_mov_epi32( a, s < float_v( 0.f ), b );
      }
      friend int_v select( const mask_v& mask, const int_v& a, const int_v& b ) {
        return _mm256_mask_mov_epi32( b, mask, a );
      }

      friend mask_v operator<( const int_v& lhs, const int_v& rhs ) { return _mm256_cmplt_epi32_mask( lhs, rhs ); }
      friend mask_v operator>( const int_v& lhs, const int_v& rhs ) { return _mm256_cmpgt_epi32_mask( lhs, rhs ); }
      friend mask_v operator==( const int_v& lhs, const int_v& rhs ) { return _mm256_cmpeq_epi32_mask( lhs, rhs ); }
      friend std::ostream& operator<<( std::ostream& os, int_v const& x ) {
        return print_vector<int, 8>( os, x, "avx256" );
      }

    private:
      __m256i data;
    };

    inline float_v::operator int_v() const { return int_v( _mm256_cvttps_epi32( data ) ); }

    inline int_v castToInt( const float_v& x ) { return int_v( _mm256_castps_si256( x ) ); }

    inline float_v castToFloat( const int_v& x ) { return float_v( _mm256_castsi256_ps( x ) ); }

    inline int_v gather( const int* base, const int_v& idx ) {
      return _mm256_i32gather_epi32( base, idx, sizeof( int ) );
    }

    inline float_v gather( const float* base, const int_v& idx ) {
      return _mm256_i32gather_ps( base, idx, sizeof( float ) );
    }

    inline int_v maskgather( const int* base, const int_v& idx, const mask_v& mask, const int_v& source ) {
      return _mm256_mmask_i32gather_epi32( source, mask, idx, base, sizeof( int ) );
    }

    inline float_v maskgather( const float* base, const int_v& idx, const mask_v& mask, const float_v& source ) {
      return _mm256_mmask_i32gather_ps( source, mask, idx, base, sizeof( float ) );
    }

    struct types {
      static const size_t size = 8;
      using int_v = avx256::int_v;
      using float_v = avx256::float_v;
      using mask_v = avx256::mask_v;
      static mask_v mask_true() { return 0xFF; }
      static mask_v mask_false() { return 0x00; }
      static int_v indices() { return _mm256_setr_epi32( 0, 1, 2, 3, 4, 5, 6, 7 ); }
      static int_v indices( int start ) { return indices() + start; }
      static int popcount( mask_v const& mask ) { return avx256::popcount( mask ); }
      static mask_v loop_mask( int i, int n ) { return ( ( i + 8 ) > n ) ? ~( 0xFF << ( n & 7 ) ) : 0xFF; }
    };
  } // namespace avx256

  namespace avx512 {
    constexpr InstructionSet instructionSet() { return AVX512; }
    class mask_v {
    public:
      mask_v() {} // Constructor must be empty
      mask_v( __mmask16 m ) : data( m ) {}

      mask_v& operator=( const __mmask16& m ) {
        data = m;
        return *this;
      }

      operator __mmask16() const { return data; }

      friend mask_v operator&&( const mask_v& lhs, const mask_v& rhs ) { return lhs & rhs; }
      friend mask_v operator||( const mask_v& lhs, const mask_v& rhs ) { return lhs | rhs; }
      friend mask_v operator!( const mask_v& x ) { return ~x; }

      friend bool all( const mask_v& mask ) { return mask == 0xFFFF; }
      friend bool none( const mask_v& mask ) { return mask == 0x0000; }
      friend bool any( const mask_v& mask ) { return mask != 0x0000; }
      friend std::ostream& operator<<( std::ostream& os, mask_v const& x ) {
        return print_bitset<16>( os, x, "avx512" );
      }

    private:
      __mmask16 data;
    };

    inline int popcount( mask_v const& mask ) { return _mm_popcnt_u32( mask ); }

    class int_v;

    class float_v {
    public:
      float_v() {} // Constructor must be empty
      float_v( scalar::float_v& f ) : data( _mm512_set1_ps( f.cast() ) ) {}
      float_v( float f ) : data( _mm512_set1_ps( f ) ) {}
      float_v( const float* f ) : data( _mm512_loadu_ps( f ) ) {}
      float_v( __m512 f ) : data( f ) {}

      float_v& operator=( const __m512& f ) {
        data = f;
        return *this;
      }

      operator __m512() const { return data; }
      inline operator int_v() const;

      void store( float* ptr ) const { _mm512_storeu_ps( ptr, data ); }

      void compressstore( const mask_v& mask, float* ptr ) const { _mm512_mask_compressstoreu_ps( ptr, mask, data ); }

      float hmax() const { return _mm512_reduce_max_ps( data ); }
      float hmin() const { return _mm512_reduce_min_ps( data ); }
      float hadd() const { return _mm512_reduce_add_ps( data ); }
      float hmax( const mask_v& mask ) const { return _mm512_mask_reduce_max_ps( mask, data ); }
      float hmin( const mask_v& mask ) const { return _mm512_mask_reduce_min_ps( mask, data ); }
      float hadd( const mask_v& mask ) const { return _mm512_mask_reduce_add_ps( mask, data ); }

      friend float_v operator+( const float_v& lhs, const float_v& rhs ) { return _mm512_add_ps( lhs, rhs ); }
      friend float_v operator-( const float_v& lhs, const float_v& rhs ) { return _mm512_sub_ps( lhs, rhs ); }
      friend float_v operator*( const float_v& lhs, const float_v& rhs ) { return _mm512_mul_ps( lhs, rhs ); }
      friend float_v operator/( const float_v& lhs, const float_v& rhs ) { return _mm512_div_ps( lhs, rhs ); }
      friend float_v operator-( const float_v& x ) { return -1.f * x; }

      friend float_v operator&( const float_v& lhs, const float_v& rhs ) { return _mm512_and_ps( lhs, rhs ); }
      friend float_v operator|( const float_v& lhs, const float_v& rhs ) { return _mm512_or_ps( lhs, rhs ); }
      friend float_v operator^( const float_v& lhs, const float_v& rhs ) { return _mm512_xor_ps( lhs, rhs ); }

      friend float_v operator&&( const float_v& lhs, const float_v& rhs ) { return _mm512_and_ps( lhs, rhs ); }
      friend float_v operator||( const float_v& lhs, const float_v& rhs ) { return _mm512_or_ps( lhs, rhs ); }
      friend float_v operator!( const float_v& x ) { return x ^ _mm512_castsi512_ps( _mm512_set1_epi32( -1 ) ); }

      friend float_v min( const float_v& lhs, const float_v& rhs ) { return _mm512_min_ps( lhs, rhs ); }
      friend float_v max( const float_v& lhs, const float_v& rhs ) { return _mm512_max_ps( lhs, rhs ); }
      friend float_v abs( const float_v& v ) { return _mm512_abs_ps( v ); }
      friend float_v copysign( const float_v& x, const float_v& y ) {
        return x ^ ( y & _mm512_castsi512_ps( _mm512_set1_epi32( 0x80000000 ) ) );
      }

      friend float_v signselect( const float_v& s, const float_v& a, const float_v& b ) {
        return _mm512_mask_mov_ps( a, s < float_v( 0.f ), b );
      }
      friend float_v select( const mask_v& mask, const float_v& a, const float_v& b ) {
        return _mm512_mask_mov_ps( b, mask, a );
      }

      friend float_v sqrt( const float_v& v ) { return _mm512_sqrt_ps( v ); }

      friend mask_v operator<( const float_v& lhs, const float_v& rhs ) {
        return _mm512_cmp_ps_mask( lhs, rhs, _CMP_LT_OS );
      }
      friend mask_v operator>( const float_v& lhs, const float_v& rhs ) {
        return _mm512_cmp_ps_mask( lhs, rhs, _CMP_GT_OS );
      }
      friend mask_v operator==( const float_v& lhs, const float_v& rhs ) {
        return _mm512_cmp_ps_mask( lhs, rhs, _CMP_EQ_OS );
      }
      friend std::ostream& operator<<( std::ostream& os, float_v const& x ) {
        return print_vector<float, 16>( os, x, "avx512" );
      }

    private:
      __m512 data;
    };

    class int_v {
    public:
      int_v() {} // Constructor must be empty
      int_v( int f ) : data( _mm512_set1_epi32( f ) ) {}
      int_v( const int* f ) : data( _mm512_loadu_si512( f ) ) {}
      constexpr int_v( __m512i f ) : data( f ) {}

      int_v& operator=( const __m512i& f ) {
        data = f;
        return *this;
      }

      operator __m512i() const { return data; }
      operator float_v() const { return float_v( _mm512_cvtepi32_ps( data ) ); }

      void store( int* ptr ) const { _mm512_storeu_si512( ptr, data ); }

      void compressstore( const mask_v& mask, int* ptr ) const { _mm512_mask_compressstoreu_epi32( ptr, mask, data ); }

      int hmax() const { return _mm512_reduce_max_epi32( data ); }
      int hmin() const { return _mm512_reduce_min_epi32( data ); }
      int hadd() const { return _mm512_reduce_add_epi32( data ); }
      int hmax( const mask_v& mask ) const { return _mm512_mask_reduce_max_epi32( mask, data ); }
      int hmin( const mask_v& mask ) const { return _mm512_mask_reduce_min_epi32( mask, data ); }
      int hadd( const mask_v& mask ) const { return _mm512_mask_reduce_add_epi32( mask, data ); }

      friend int_v operator+( const int_v& lhs, const int_v& rhs ) { return _mm512_add_epi32( lhs, rhs ); }
      friend int_v operator-( const int_v& lhs, const int_v& rhs ) { return _mm512_sub_epi32( lhs, rhs ); }
      friend int_v operator*( const int_v& lhs, const int_v& rhs ) { return _mm512_mul_epi32( lhs, rhs ); }

      friend int_v operator&( const int_v& lhs, const int_v& rhs ) { return _mm512_and_si512( lhs, rhs ); }
      friend int_v operator|( const int_v& lhs, const int_v& rhs ) { return _mm512_or_si512( lhs, rhs ); }

      friend int_v signselect( const float_v& s, const int_v& a, const int_v& b ) {
        return _mm512_mask_mov_epi32( a, s < float_v( 0.f ), b );
      }
      friend int_v select( const mask_v& mask, const int_v& a, const int_v& b ) {
        return _mm512_mask_mov_epi32( b, mask, a );
      }

      friend mask_v operator<( const int_v& lhs, const int_v& rhs ) { return _mm512_cmplt_epi32_mask( lhs, rhs ); }
      friend mask_v operator>( const int_v& lhs, const int_v& rhs ) { return _mm512_cmpgt_epi32_mask( lhs, rhs ); }
      friend mask_v operator==( const int_v& lhs, const int_v& rhs ) { return _mm512_cmpeq_epi32_mask( lhs, rhs ); }
      friend std::ostream& operator<<( std::ostream& os, int_v const& x ) {
        return print_vector<int, 16>( os, x, "avx512" );
      }

    private:
      __m512i data;
    };

    inline float_v::operator int_v() const { return int_v( _mm512_cvttps_epi32( data ) ); }

    inline int_v castToInt( const float_v& x ) { return int_v( _mm512_castps_si512( x ) ); }

    inline float_v castToFloat( const int_v& x ) { return float_v( _mm512_castsi512_ps( x ) ); }

    inline int_v gather( const int* base, const int_v& idx ) {
      return _mm512_i32gather_epi32( idx, base, sizeof( int ) );
    }

    inline float_v gather( const float* base, const int_v& idx ) {
      return _mm512_i32gather_ps( idx, base, sizeof( float ) );
    }

    inline int_v maskgather( const int* base, const int_v& idx, const mask_v& mask, const int_v& source ) {
      return _mm512_mask_i32gather_epi32( source, mask, idx, base, sizeof( int ) );
    }

    inline float_v maskgather( const float* base, const int_v& idx, const mask_v& mask, const float_v& source ) {
      return _mm512_mask_i32gather_ps( source, mask, idx, base, sizeof( float ) );
    }

    struct types {
      static const size_t size = 16;
      using int_v = avx512::int_v;
      using float_v = avx512::float_v;
      using mask_v = avx512::mask_v;
      static mask_v mask_true() { return 0xFFFF; }
      static mask_v mask_false() { return 0x0000; }
      static int_v indices() { return _mm512_setr_epi32( 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 ); }
      static int_v indices( int start ) { return indices() + start; }
      static int popcount( mask_v const& mask ) { return avx512::popcount( mask ); }
      static mask_v loop_mask( int i, int n ) { return ( ( i + 16 ) > n ) ? ~( 0xFFFF << ( n & 15 ) ) : 0xFFFF; }
    };
  } // namespace avx512
} // namespace SIMDWrapper

template <>
struct LHCb::type_map<SIMDWrapper::avx256::float_v> {
  using int_t = SIMDWrapper::avx256::int_v;
};

template <>
struct LHCb::type_map<SIMDWrapper::avx512::float_v> {
  using int_t = SIMDWrapper::avx512::int_v;
};

namespace SIMDWrapper {
#endif

  namespace best {
    using float_v = avx256::float_v;
    using int_v   = avx256::int_v;
    using types   = avx256::types;
    constexpr InstructionSet instructionSet() { return avx256::instructionSet(); }
  } // namespace best
} // namespace SIMDWrapper

// Helper to build SIMDWrapper accessors:
#define SOA_ACCESSOR( name, location )                                                                                 \
  template <typename T>                                                                                                \
  inline auto name( const int key ) const {                                                                            \
    return T( location + key );                                                                                        \
  }                                                                                                                    \
  template <typename T>                                                                                                \
  inline auto gather_##name( const T& key ) const {                                                                    \
    return gather( location, key );                                                                                    \
  }                                                                                                                    \
  template <typename T, typename KeyT, typename MaskT>                                                                 \
  inline auto maskgather_##name( const KeyT& key, const MaskT& mask, const T& src ) const {                            \
    return maskgather( location, key, mask, src );                                                                     \
  }                                                                                                                    \
  template <typename T>                                                                                                \
  inline auto store_##name( const int key, const T& v ) {                                                              \
    v.store( location + key );                                                                                         \
  }                                                                                                                    \
  template <typename T, typename MaskT>                                                                                \
  inline auto compressstore_##name( const int key, const MaskT& mask, const T& v ) {                                   \
    v.compressstore( mask, location + key );                                                                           \
  }

// With variadic parameters:
#define SOA_ACCESSOR_VAR( name, location, ... )                                                                        \
  template <typename T>                                                                                                \
  inline auto name( const int key, __VA_ARGS__ ) const {                                                               \
    return T( location + key );                                                                                        \
  }                                                                                                                    \
  template <typename T>                                                                                                \
  inline auto gather_##name( const T& key, __VA_ARGS__ ) const {                                                       \
    return gather( location, key );                                                                                    \
  }                                                                                                                    \
  template <typename T, typename KeyT, typename MaskT>                                                                 \
  inline auto maskgather_##name( const T& key, const MaskT& mask, const T& src, __VA_ARGS__ ) const {                  \
    return maskgather( location, key, mask, src );                                                                     \
  }                                                                                                                    \
  template <typename T>                                                                                                \
  inline auto store_##name( const int key, __VA_ARGS__, const T& v ) {                                                 \
    v.store( location + key );                                                                                         \
  }                                                                                                                    \
  template <typename T, typename MaskT>                                                                                \
  inline auto compressstore_##name( const int key, __VA_ARGS__, MaskT mask, const T& v ) {                             \
    v.compressstore( mask, location + key );                                                                           \
  }

// Align size to AVX512 boundary
constexpr int align_size( int n ) { return ( n / 16 + 1 ) * 16; }

// Maths :

template <typename T>
inline T faster_atan2( const T& y, const T& x ) { // error < 0.07 rad, no 0/0 security
  const T c1    = M_PI / 4.0;
  const T c2    = 3.0 * M_PI / 4.0;
  T       abs_y = abs( y );

  T x_plus_y = x + abs_y;
  T x_sub_y  = x - abs_y;
  T y_sub_x  = abs_y - x;

  T nom = signselect( x, x_sub_y, x_plus_y );
  T den = signselect( x, x_plus_y, y_sub_x );

  T r     = nom / den;
  T angle = signselect( x, c1, c2 ) - c1 * r;

  return copysign( angle, y );
}

template <typename T>
inline T vapprox_log( const T& x ) {
  return castToInt( x ) * 8.2629582881927490e-8f - 87.989971088f;
}

template <typename M, typename T>
inline void swap( const M& mask, T& a, T& b ) {
  T tmp = select( mask, b, a );
  b     = select( mask, a, b );
  a     = tmp;
}
