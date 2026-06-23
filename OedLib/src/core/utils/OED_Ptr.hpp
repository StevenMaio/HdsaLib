#pragma once

#include <memory>

#ifdef MrHyDE_ENABLE_OED

#include "Teuchos_RCP.hpp"

namespace OED
{

  template <class T>
  using Ptr = Teuchos::RCP<T>;

  static const Teuchos::ENull nullPtr = Teuchos::null;

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
#else // use std::shared_ptr
namespace OED
{
  template <class T>
  using Ptr = std::shared_ptr<T>;

  static const nullptr_t nullPtr = nullptr;

  template <class T, class... Args>
  inline Ptr<T> makePtr(Args &&...args)
  {
    return std::make_shared<T>(std::forward<Args>(args)...);
  }

  template <class T>
  inline Ptr<T> makePtrFromRef(T &obj)
  {
    // this function isn't used in OedLib
    return OED::nullPtr;
  }

  template <class T, class U>
  inline Ptr<T> dynamicPtrCast(const Ptr<U> &r) noexcept
  {
    return dynamic_cast<T>(r);
  }

}
#endif