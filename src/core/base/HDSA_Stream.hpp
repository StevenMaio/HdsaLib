#pragma once
#ifndef HDSA_STREAM_HPP
#define HDSA_STREAM_HPP

#include <ostream>
#include <string>

namespace HDSA
{

  template <typename _CharT, typename _Traits>
  class basic_nullstream : virtual public std::basic_ostream<_CharT, _Traits>
  {
  public:
    explicit basic_nullstream() : std::basic_ostream<_CharT, _Traits>(NULL) {}
  };

  using nullstream = basic_nullstream<char, std::char_traits<char>>;

}

#endif
