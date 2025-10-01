#pragma once

#include <memory>
#include <type_traits>

#include <cstddef>
#include <utility>

#include "Teuchos_RCP.hpp"

namespace HDSA
{

  template <class T>
  using Ptr = Teuchos::RCP<T>;

  static const Teuchos::ENull nullPtr = Teuchos::null;

}

namespace HDSA
{

  template <class T, class... Args>
  inline Ptr<T> makePtr(Args &&...args)
  {
    return Teuchos::rcp(new T(std::forward<Args>(args)...));
  }

  template <class T>
  inline Ptr<T> makePtrFromRef(T &obj)
  {
    return Teuchos::rcpFromRef(obj);
  }

  template <class T, class U>
  inline Ptr<T> dynamicPtrCast(const Ptr<U> &r) noexcept
  {
    return Teuchos::rcp_dynamic_cast<T>(r);
  }

}
