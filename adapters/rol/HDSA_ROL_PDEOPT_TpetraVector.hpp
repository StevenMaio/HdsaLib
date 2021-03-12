#ifndef HDSA_ROL_PDEOPT_TPETRAVECTOR_HPP
#define HDSA_ROL_PDEOPT_TPETRAVECTOR_HPP

#include "ROL_Vector.hpp"
#include "../../../PDE-OPT/TOOLS/pdevector.hpp"
#include "../../../PDE-OPT/TOOLS/pde.hpp"
#include "../../../PDE-OPT/TOOLS/dynpde.hpp"
#include "../../../PDE-OPT/TOOLS/assembler_def.hpp"

template <class RealT>
class ROL_PDEOPT_Tpetra_Vector : public ROL_Vector<RealT> {
  
private:
  HDSA::Ptr<Tpetra::MultiVector<> > vec_ptr_;
  const HDSA::Ptr<PDE<RealT> > pde_;
  const HDSA::Ptr<DynamicPDE<RealT> > dyn_pde_;
  const HDSA::Ptr<Assembler<RealT> > assembler_;
  HDSA::Ptr<const Teuchos::Comm<int> > comm_;

public:
  
  ROL_PDEOPT_Tpetra_Vector(const HDSA::Ptr<PDE<RealT> > & pde, const HDSA::Ptr<Assembler<RealT> > & assembler): pde_(pde), assembler_(assembler)
  {
    vec_ptr_  = assembler->createControlVector();  
    vec_ptr_->putScalar(0.0);
    ROL_Vector<RealT>::rol_vec_ = HDSA::makePtr<PDE_PrimalOptVector<RealT> >(vec_ptr_,pde_,assembler_);
    comm_ = vec_ptr_->getMap()->getComm();
  }
 
  ROL_PDEOPT_Tpetra_Vector(const HDSA::Ptr<DynamicPDE<RealT> > & dyn_pde, const HDSA::Ptr<Assembler<RealT> > & assembler): dyn_pde_(dyn_pde), assembler_(assembler)
  {
    vec_ptr_  = assembler->createControlVector();  
    vec_ptr_->putScalar(0.0);
    ROL_Vector<RealT>::rol_vec_ = HDSA::makePtr<PDE_PrimalOptVector<RealT> >(vec_ptr_,dyn_pde_,*assembler_);
    comm_ = vec_ptr_->getMap()->getComm();
  }
 
  ~ROL_PDEOPT_Tpetra_Vector()
  { }
  
  // Access the kth element of the vector
  RealT operator () (int k) const 
  {
    RealT val;
    (*vec_ptr_).sync<Kokkos::HostSpace> ();
    auto e_vec_2d = (*vec_ptr_).getLocalView<Kokkos::HostSpace> ();
    auto e_vec_1d = Kokkos::subview (e_vec_2d, Kokkos::ALL (), 0);
    if(vec_ptr_->getMap()->isNodeGlobalElement(k))
      {
  	val = e_vec_1d(vec_ptr_->getMap()->getLocalElement(k));
      }    
    char *buff = (char*)&val;
    long long int k_ll = (long long int)k;
    Teuchos::ArrayView<long long int> node = Teuchos::ArrayView<long long int>(&k_ll,1);
    int proc = 0;
    Teuchos::ArrayView<int> node_id = Teuchos::ArrayView<int>(&proc,1);
    auto lookup = vec_ptr_->getMap()->getRemoteIndexList(node,node_id);
    if(lookup==Tpetra::IDNotPresent)
      {
  	std::cout << "There is an error in vector elementwise access" << std::endl;
      }
    comm_->barrier();
    comm_->broadcast(proc,8,buff);
    
    return val;
  }
  

  // Replace the kth element of the vector by val
  void Replace_Element(int k, RealT val)
  {
    vec_ptr_->replaceGlobalValue(k,0,val);
  } 

  // Get the data on this processor
  std::vector<RealT> Get_Data_on_Processor(void) const
  {
    // Get data from Tpetra vector into data
    Teuchos::ArrayRCP<RealT> array = vec_ptr_->getDataNonConst(0);
    std::vector<RealT> data = std::vector<RealT>(array.size());
    for(unsigned int k = 0; k < data.size(); k++)
      {
	data[k] = array[k];
      }
    return data;
  }

  // Get the indices on this processor
  std::vector<int> Get_Indices_on_Processor(void) const 
  {
    auto array = vec_ptr_->getMap()->getMyGlobalIndices();
    std::vector<int> indices = std::vector<int>(array.size());
    for(unsigned int k = 0; k < indices.size(); k++)
      {
	indices[k] = array[k];
      }
    return indices;
  }

  HDSA::Ptr<HDSA::Vector<RealT> > Define_Clone() const
  {
    HDSA::Ptr<HDSA::Vector<RealT> > vec;
    if(pde_ != HDSA::nullPtr)
      {
	vec = HDSA::makePtr<ROL_PDEOPT_Tpetra_Vector<RealT> >(pde_,assembler_);
      }
    else
      {
	vec = HDSA::makePtr<ROL_PDEOPT_Tpetra_Vector<RealT> >(dyn_pde_,assembler_);
      }
    return vec;
  }
 
  HDSA::Ptr<Tpetra::MultiVector<> > get_tpetra_vec(void)
  {
    return vec_ptr_;
  }

};


#endif


