/***********************************************************************
 HdsaLib - A library for Hyper-differential Sensitivity Analysis
 
 Questions? Contact Joseph Hart (joshart@sandia.gov)
************************************************************************/

#pragma once
#include "HDSA_Ptr.hpp"
#include "Teuchos_XMLParameterListCoreHelpers.hpp"

namespace HDSA
{
  using ParameterList = Teuchos::ParameterList;

  inline void updateParametersFromXmlFile(const std::string &filename, ParameterList &parlist)
  {
    Teuchos::Ptr<ParameterList> p{&parlist};
    Teuchos::updateParametersFromXmlFile(filename, p);
  }

} // namespace HDSA
