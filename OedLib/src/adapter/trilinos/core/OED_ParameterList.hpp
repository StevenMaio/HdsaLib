#pragma once
#include "OED_Ptr.hpp"
#include "Teuchos_XMLParameterListCoreHelpers.hpp"

namespace OED
{

  namespace Trilinos_Adapter
  {
    using ParameterList = Teuchos::ParameterList;

    inline void updateParametersFromXmlFile(const std::string &filename, ParameterList &parlist)
    {
      Teuchos::Ptr<ParameterList> p{&parlist};
      Teuchos::updateParametersFromXmlFile(filename, p);
    }

  }

} // namespace OED
