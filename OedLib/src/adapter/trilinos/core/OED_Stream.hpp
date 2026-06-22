#pragma once
#ifndef OED_STREAM_HPP
#define OED_STREAM_HPP

#include <ostream>
#include <string>

namespace OED
{


  namespace Trilinos_Adapter
  {

  template <typename _CharT, typename _Traits>
  class basic_nullstream : virtual public std::basic_ostream<_CharT, _Traits>
  {
  public:
    explicit basic_nullstream() : std::basic_ostream<_CharT, _Traits>(NULL) {}
  };

  using nullstream = basic_nullstream<char, std::char_traits<char>>;

  }

}

#endif
