#ifndef HDSA_ROL_PDEOPT_STDVECTOR_HPP
#define HDSA_ROL_PDEOPT_STDVECTOR_HPP

#include <algorithm>
#include <cstdlib>

#include "ROL_Vector.hpp"
#include "ROL_StdVector.hpp"
#include "../../../PDE-OPT/TOOLS/pdevector.hpp"

template <class RealT>
class ROL_PDEOPT_Std_Vector : public ROL_Vector<RealT> {
  
private:
  int dim_;
  HDSA::Ptr<std::vector<RealT> > vec_ptr_;
  
public:
  
  ROL_PDEOPT_Std_Vector(int dim): dim_(dim)
  {
    vec_ptr_ = HDSA::makePtr<std::vector<RealT> >(dim,0.0);
    ROL_Vector<RealT>::rol_vec_ = HDSA::makePtr<PDE_OptVector<RealT> >(HDSA::makePtr<ROL::StdVector<RealT> >(vec_ptr_));
  }
  
  ~ROL_PDEOPT_Std_Vector()
  { }
  
  // Access the (i,j) element
  RealT operator () (int k) const 
  {
    return (*vec_ptr_)[k];
  }

  // Replace the kth element of the vector by val
  void Replace_Element(int k, RealT val)
  {
    (*vec_ptr_)[k] = val;
  } 

  // Get the data on this processor
  std::vector<RealT> Get_Data_on_Processor(void) const
  {
    return *vec_ptr_;
  }

  // Get the indices on this processor
  std::vector<int> Get_Indices_on_Processor(void) const 
  {
    std::vector<int> indices = std::vector<int>(vec_ptr_->size());
    for(unsigned int k = 0; k < indices.size(); k++)
      {
	indices[k] = k;
      }
    return indices;
  }

  HDSA::Ptr<HDSA::Vector<RealT> > Define_Clone() const
  {
    HDSA::Ptr<HDSA::Vector<RealT> > vec = HDSA::makePtr<ROL_PDEOPT_Std_Vector<RealT> >(dim_);
    return vec;
  }

};


#endif
