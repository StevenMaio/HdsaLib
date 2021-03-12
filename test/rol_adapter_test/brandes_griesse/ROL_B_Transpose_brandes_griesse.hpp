#ifndef ROL_B_TRANSPOSE_BRANDES_GRIESSE_HPP
#define ROL_B_TRANSPOSE_BRANDES_GRIESSE_HPP

#include "../../../../PDE-OPT/TOOLS/pdeconstraint.hpp"

template <class RealT>
class ROL_B_Transpose_brandes_griesse: public ROL_B_Transpose_Steady_State<RealT> {

private:
  HDSA::Ptr<PDE_Brandes_Griesse_Parameter<RealT> > pde_param_;
  HDSA::Ptr<PDE_Constraint<RealT> > con_param_;
  HDSA::Ptr<Assembler<RealT> > assembler_;
  HDSA::Ptr<Intrepid::FieldContainer<RealT> > z_coeff_;

public:

  ROL_B_Transpose_brandes_griesse(const HDSA::Ptr<MeshManager<RealT> > & meshMgr, const HDSA::Ptr<Assembler<RealT> > & assembler,
				  const HDSA::Ptr<Teuchos::ParameterList > & parlist, const HDSA::Ptr<std::ostream> & outStream, const HDSA::Ptr<const Teuchos::Comm<int> > & comm,
				  const HDSA::Ptr<ROL::Objective_SimOpt<RealT> > & obj, const HDSA::Ptr<ROL::Constraint_SimOpt<RealT> > & con, const HDSA::Ptr<ROL::Vector<RealT> > & u,
				  const HDSA::Ptr<ROL::Vector<RealT> > & z, const HDSA::Ptr<HDSA::Vector<RealT> > & theta): ROL_B_Transpose_Steady_State<RealT>(obj,con, u, z, theta)
  {
    assembler_ = assembler;
    pde_param_ = HDSA::makePtr<PDE_Brandes_Griesse_Parameter<RealT> >(*parlist);
    con_param_ = HDSA::makePtr<PDE_Constraint<RealT> >(pde_param_,meshMgr,comm,*parlist,*outStream);
  }

  virtual ~ROL_B_Transpose_brandes_griesse()
  { }

  void applyThetaJacobianTranspose(HDSA::Vector<RealT> & grad, const ROL::Vector<RealT> & lambda, const ROL::Vector<RealT> & u, const ROL::Vector<RealT> & z, const HDSA::Vector<RealT> & theta)
  {
    if(z_coeff_ == HDSA::nullPtr)
      {
	z_coeff_ = get_z_Field_Container(z);
	pde_param_->Update_Z_input(z_coeff_);
      }

    const Std_Vector<RealT> &etheta = dynamic_cast<const Std_Vector<RealT>&>(theta);
    HDSA::Ptr<std::vector<RealT> > theta_std = etheta.get_std_vec();
    Std_Vector<RealT> &egrad = dynamic_cast<Std_Vector<RealT>&>(grad);
    HDSA::Ptr<std::vector<RealT> > grad_std = egrad.get_std_vec();

    HDSA::Ptr<ROL::Vector<RealT> > theta_rol = HDSA::makePtr<PDE_OptVector<RealT> >(HDSA::makePtr<ROL::StdVector<RealT> >(theta_std));
    HDSA::Ptr<ROL::Vector<RealT> > grad_rol = HDSA::makePtr<PDE_OptVector<RealT> >(HDSA::makePtr<ROL::StdVector<RealT> >(grad_std));

    RealT tol = 1.e-8;
    con_param_->update(u,*theta_rol);
    con_param_->applyAdjointJacobian_2(*grad_rol,lambda,u,*theta_rol,tol);
  } 

  void Update_Constraint(const HDSA::Ptr<ROL::Constraint_SimOpt<RealT> > & con, const ROL::Vector<RealT> & u,
			 const ROL::Vector<RealT> & z, const HDSA::Vector<RealT> & theta) const
  {
    HDSA::Ptr<Linear_PDE_Constraint<RealT> > pdecon = HDSA::dynamicPtrCast<Linear_PDE_Constraint<RealT> >(con);
    const Std_Vector<RealT> &etheta = dynamic_cast<const Std_Vector<RealT>&>(theta);
    HDSA::Ptr<std::vector<RealT> > theta_std = etheta.get_std_vec();
    pdecon->setParameter(*theta_std);
    con->update(u,z);
  }
  
private:

  HDSA::Ptr<Intrepid::FieldContainer<RealT> > get_z_Field_Container(const ROL::Vector<RealT> & z)
  {
    HDSA::Ptr<Intrepid::FieldContainer<RealT> > z_coeff = assembler_->get_z_field_container(dynamic_cast<const ROL::TpetraMultiVector<RealT>&>(z).getVector());
    return z_coeff;
  }

};


#endif
