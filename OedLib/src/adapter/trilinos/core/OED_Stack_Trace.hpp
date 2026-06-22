#pragma once
#ifndef OED_STACKTRACE_HPP
#define OED_STACKTRACE_HPP

#include <iostream>
#include <stdexcept>
#include "Teuchos_TestForException.hpp"

#define OED_TEST_FOR_EXCEPTION(throw_exception_test, Exception, msg) \
  TEUCHOS_TEST_FOR_EXCEPTION(throw_exception_test, Exception, msg)

#endif
