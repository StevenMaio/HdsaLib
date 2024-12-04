/***********************************************************************
 Hyper-Differential Sensitivity Analysis Library (HDSAlib)
 
 Copyright 2018 National Technology & Engineering Solutions of Sandia,
 LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the
 U.S. Government retains certain rights in this software.”
 
 Questions? Contact Joseph Hart (joshart@sandia.gov) and/or
 Bart van Bloemen Waanders (bartv@sandia.gov)
 ************************************************************************/

#ifndef OBJ_MRHYDE_HPP
#define OBJ_MRHYDE_HPP

#include "ROL_StdVector.hpp"
#include "ROL_RiskVector.hpp"
#include "ROL_Objective.hpp"
#include "ROL_BoundConstraint.hpp"
#include "Teuchos_ParameterList.hpp"

#include "solverManager.hpp"
#include "postprocessManager.hpp"
#include "parameterManager.hpp"

//#include "ROL_RiskVector.hpp"

#include <iostream>
#include <fstream>
#include <string>

//#include <random> //for normal noise...not sure if this is necessary...

namespace HDSA {

  //  using namespace MrHyDE;  
  
  template<class Real>
  class Objective_MrHyDE {
    
  private:
    
    Real noise_;                                            //standard deviation of normal additive noise to add to data (0 for now)
    Teuchos::RCP<MrHyDE::SolverManager<SolverNode> > solver_;                                     // Solver object for MILO (solves FWD, ADJ, computes gradient, etc.)
    Teuchos::RCP<MrHyDE::PostprocessManager<SolverNode> > postproc_;                              // Postprocessing object for MILO (write solution, computes response, etc.)
    Teuchos::RCP<MrHyDE::ParameterManager<SolverNode> > params_;

  public:
    
    /*!
     \brief A constructor generating data
     */
    Objective_MrHyDE(Teuchos::RCP<MrHyDE::SolverManager<SolverNode> > solver,
                   Teuchos::RCP<MrHyDE::PostprocessManager<SolverNode> > postproc,
                   Teuchos::RCP<MrHyDE::ParameterManager<SolverNode> > & params) :
    solver_(solver), postproc_(postproc), params_(params) {
      
    } //end constructor
    
    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    
    
    //! Compute gradient of objective function with respect to parameters
    void gradient(HDSA::Vector<Real> &g, const HDSA::Vector<Real> &Params, Real &tol){

      bool newparams = this->checkNewParams(Params);

      if (newparams) {
	MrHyDE_OptVector curr_params = params_->getCurrentVector();
	ROL::Ptr<ROL::Vector<Real> > rparams = curr_params.clone();

	// correct cast
        MrHyDE_OptVector Paramsp = 
          Teuchos::dyn_cast<MrHyDE_OptVector >(dynamic_cast<ROL::Vector<Real> &>(*rparams));

	const HDSA::Vector_MrHyDE<Real> &eParams = dynamic_cast<const HDSA::Vector_MrHyDE<Real>&>(Params);
	rparams->set(*eParams.mrhyde_vec);
      
        params_->updateParams(Paramsp);
        DFAD val = 0.0;
        solver_->forwardModel(val);
	Real value = val.val();
	std::cout << "value = " << value << std::endl;
      }
       HDSA::Vector_MrHyDE<Real> sens = 
	 Teuchos::dyn_cast<HDSA::Vector_MrHyDE<Real> >(dynamic_cast<HDSA::Vector<Real> &>(g));

      sens.zeros();
      MrHyDE_OptVector esens = 
          Teuchos::dyn_cast<MrHyDE_OptVector >(dynamic_cast<ROL::Vector<Real> &>(*sens.mrhyde_vec));
      
      solver_->adjointModel(esens);
      //      g.scale(0.5); //bvbw test
    }

        //! Compute gradient of objective function with respect to parameters
    void Misfit_Gradient(HDSA::Vector<Real> &g, const HDSA::Vector<Real> &state, const HDSA::Vector<Real> &Params, Real &tol){

      bool newparams = this->checkNewParams(Params);

      if (newparams) {
        MrHyDE_OptVector Paramsp = 
        Teuchos::dyn_cast<MrHyDE_OptVector >(const_cast<HDSA::Vector<Real> &>(Params));
      
        params_->updateParams(Paramsp);
        DFAD val = 0.0;
        solver_->forwardModel(val);

      }
      // MrHyDE_OptVector sens = 
      // Teuchos::dyn_cast<MrHyDE_OptVector >(const_cast<HDSA::Vector<Real> &>(g));
      // sens.zero();
      
      // solver->adjointModel(sens);

    }
    
    bool checkNewParams(const HDSA::Vector<Real> &Params) {
      MrHyDE_OptVector curr_params = params_->getCurrentVector();
      //      std::cout << "curr_params" << curr_params << std::endl;
      // curr_params.print(std::cout);
      //      Kokkos::View<ScalarT***,HostDevice> cparams("curr_params", curr_params);
      // KokkosTools::print(cparams,"curr_params");

      ROL::Ptr<ROL::Vector<Real> > diff = curr_params.clone();
      MrHyDE_OptVector ediff = 
      Teuchos::dyn_cast<MrHyDE_OptVector >(const_cast<ROL::Vector<Real> &>(*diff));
      HDSA::Vector_MrHyDE<Real> eParams = 
	Teuchos::dyn_cast<HDSA::Vector_MrHyDE<Real> >(const_cast<HDSA::Vector<Real> &>(Params));
      ediff.zero();
      ediff.set(curr_params);
      ediff.axpy(-1.0,*eParams.mrhyde_vec);
      ScalarT dnorm = ediff.norm();
      ScalarT refnorm = curr_params.norm();
      dnorm = dnorm/refnorm;
      ScalarT reltol = 1.0e-12;
      bool newparams = false;
      if (dnorm > reltol) {
        newparams = true;
      }
      return newparams;
    }

    //! Compute the Hessian-vector product of the objective function
    void hessVec(HDSA::Vector<Real> &hv, const HDSA::Vector<Real> &v, const HDSA::Vector<Real> &Params, Real &tol ){
      //this->ROL::Objective<Real>::hessVec(hv,v,Params,tol);  need FD routine
      HDSA::Ptr<HDSA::Vector<Real> > g = Params.clone();
      HDSA::Ptr<HDSA::Vector<Real> > Params_pert = Params.clone();
      this->gradient(*g,Params,tol);
      Params_pert->set(Params);
      Real h = 1.e-4;
      Params_pert->axpy(h,v);
      this->gradient(hv,*Params_pert,tol);
      hv.axpy(-1.0,*g);
      hv.scale(1/h);
    }
    
    //print out Hessian (estimated via component-wise FD; to get inverse covariance in linear-Gaussian Bayesian inverse problem)
    // void printHess(const string & filename, const Vector<Real> & xin, const int & commrank){
    //   StdVector<Real> x = Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(xin));
    //   int paramDim = x.getVector()->size();
      
    //   Teuchos::RCP<StdVector<Real> > g = Teuchos::rcp(new StdVector<Real>(Teuchos::rcp(new vector<Real>(paramDim,0.0))));
    //   ScalarT gtol = sqrt(ROL_EPSILON<Real>());
    //   this->gradient(*g,x,gtol);
    //   Teuchos::RCP<StdVector<Real> > gnew = Teuchos::rcp(new StdVector<Real>(Teuchos::rcp(new vector<Real>(paramDim,0.0))));
      
    //   vector<vector<ScalarT> > hessStash(paramDim);
      
    //   //Real h = 1.e-3*x.norm(); //step length
    //   Real h = std::max(static_cast<Real>(1.0),
    //                     x.norm())*sqrt(ROL_EPSILON<Real>()); ///step length...more like what ROL has...
      
    //   //perturb each component
    //   for(int i=0; i<x.dimension(); i++){
    //     //compute new step
    //     Teuchos::RCP<StdVector<Real> > xnew = Teuchos::rcp(new StdVector<Real>(Teuchos::rcp(new vector<Real>(paramDim,0.0))));
    //     xnew->set(x);
    //     xnew->axpy(h,*x.basis(i));
        
    //     //gradient at new step
    //     gnew->zero();
    //     this->gradient(*gnew,*xnew,gtol);
        
    //     //i-th column (or row...) of Hessian
    //     gnew->axpy(static_cast<Real>(-1.0),*g);
    //     gnew->scale(static_cast<Real>(1.0)/h);
    //     Teuchos::RCP<vector<ScalarT> > gnewv = gnew->getVector();
        
    //     vector<ScalarT> row(paramDim);
    //     for(int j=0; j<paramDim; j++)
    //       row[j] = (*gnewv)[j];
    //     hessStash[i] = row;
    //   }
      
    //   if(commrank == 0){
    //     std::ofstream respOUT(filename);
    //     respOUT.precision(16);
    //     for(int i=0; i<paramDim; i++){
    //       for(int j=0; j<paramDim; j++)
    //         respOUT << hessStash[i][j] << " ";
    //       respOUT << endl;
    //     }
    //     respOUT.close();
    //   }
    // }
    
    /*!
     \brief Generate data to plot objective function
     */
    void generate_plot(Real difflo, Real diffup, Real diffstep){
      
      Teuchos::RCP<std::vector<ScalarT> > Params_rcp = Teuchos::rcp(new std::vector<ScalarT>(1,0.0) );
      ROL::StdVector<ScalarT> Params(Params_rcp);
      std::ofstream output ("Objective.dat");
      
      Real diff = 0.0;
      Real val = 0.0;
      Real tol = 1.e-16;
      int n = (diffup-difflo)/diffstep + 1;
      for(int i=0;i<n;i++){
        diff = difflo + i*diffstep;
        (*Params_rcp)[0] = diff;
        val = this->value(Params,tol);
        if(output.is_open()){
          output << std::scientific << diff << " " << val << "\n";
        }
      }
      output.close();
    }
    
    /*!
     \brief Generate data to plot objective function
     */
    void generate_plot(Real alo, Real aup, Real astep, Real blo, Real bup, Real bstep){
      
      Teuchos::RCP<std::vector<ScalarT> > Params_rcp = Teuchos::rcp(new std::vector<ScalarT>(2,0.0) );
      ROL::StdVector<ScalarT> Params(Params_rcp);
      std::ofstream output ("Objective.dat");
      
      Real a = 0.0;
      Real b = 0.0;
      Real val = 0.0;
      Real tol = 1.e-16;
      int n = (aup-alo)/astep + 1;
      int m = (bup-blo)/bstep + 1;
      for(int i=0;i<n;i++){
        a = alo + i*astep;
        for(int j=0;j<m;j++){
          b = blo + j*bstep;
          (*Params_rcp)[0] = a;
          (*Params_rcp)[1] = b;
          val = this->value(Params,tol);
          if(output.is_open()){
            output << std::scientific << a << " " << b << " " << val << "\n";
          }
        }
      }
      output.close();
    }

    void do_solop(bool solop_flag)
    {
      postproc_->hdsa_solop = solop_flag;
    }
    
  }; //end description of Objective class


  /*!
   \brief Inequality constraints on optimization parameters
   
   
   ---
   */
  
}//end namespace ROL

#endif

