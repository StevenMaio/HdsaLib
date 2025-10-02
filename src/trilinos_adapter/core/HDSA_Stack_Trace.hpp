#pragma once
#ifndef HDSA_STACKTRACE_HPP
#define HDSA_STACKTRACE_HPP

#include <iostream>
#include <stdexcept>
#include "Teuchos_TestForException.hpp"

#define HDSA_TEST_FOR_EXCEPTION(throw_exception_test, Exception, msg) \
  TEUCHOS_TEST_FOR_EXCEPTION(throw_exception_test, Exception, msg)

#endif
