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

#include <cstring>
#include <type_traits>

// type punning in way which does not rely on C99 compatiblity...
// (i.e. do _not_ use a union, but use std::memcpy, which does
// a reinterpret_cast<unsigned char*> behind the scenes, and thus
// avoids undefined behaviour -- see eg. the 'notes' section of
// http://en.cppreference.com/w/cpp/string/byte/memcpy)
// [ note: gcc explicitly supports the use of a union for type-punning
// as an extension -- and generates the same code as the std::memcpy version ]
//
// TODO: C++20: replace with std::bit_cast

template <typename Out, typename In>
constexpr std::enable_if_t<
    ( sizeof( Out ) == sizeof( In ) ) && std::is_trivially_copyable_v<In> && std::is_trivial_v<Out>, Out>
bit_cast( In src ) {
  static_assert( !std::is_same_v<std::decay_t<Out>, std::decay_t<In>>,
                 "Identity punning not allowed -- just remove the punning!" );
  Out dest{0}; // constexpr forbids uninitialized variables....
  std::memcpy( &dest, &src, sizeof( Out ) );
  return dest;
}
