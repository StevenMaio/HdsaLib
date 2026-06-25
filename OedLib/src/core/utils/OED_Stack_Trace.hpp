#pragma once
#ifndef OED_STACKTRACE_HPP
#define OED_STACKTRACE_HPP

#ifdef MrHyDE_ENABLE_OED
#include <iostream>
#include <stdexcept>
#include "Teuchos_TestForException.hpp"

#define OED_TEST_FOR_EXCEPTION(throw_exception_test, Exception, msg) \
  TEUCHOS_TEST_FOR_EXCEPTION(throw_exception_test, Exception, msg)
#else
// TODO: fix this or w/e
#define OED_TEST_FOR_EXCEPTION(throw_exception_test, Exception, msg)
// { \
//   const bool throw_exception = (throw_exception_test); \
//   if (throw_exception) \
//   { \
//     throw Exception(msg); \
//   } \
// }
#endif

#endif
