#ifndef PDE_SHALLOW_ICE_PARAM_HPP
#define PDE_SHALLOW_ICE_PARAM_HPP

#include "../../../PDE-OPT/TOOLS/pde.hpp"
#include "../../../PDE-OPT/TOOLS/fe.hpp"

#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_CellTools.hpp"

#include "ROL_Ptr.hpp"

template <class Real>
class PDE_shallow_ice_param : public PDE<Real> {
private:
  // Finite element basis information
  ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > basisPtr_;
  std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > > basisPtrs_;
  // Cell cubature information
  ROL::Ptr<Intrepid::Cubature<Real> > cellCub_;
  ROL::Ptr<Intrepid::Cubature<Real> > bdryCub_;
  // Cell node information
  ROL::Ptr<Intrepid::FieldContainer<Real> > volCellNodes_;
  std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > bdryCellNodes_;
  std::vector<std::vector<std::vector<int> > > bdryCellLocIds_;
  // Finite element definition
  ROL::Ptr<FE<Real> > fe_;
  std::vector<std::vector<ROL::Ptr<FE<Real> > > > feBdry_;
  // Local degrees of freedom on boundary, for each side of the reference cell (first index).
  std::vector<std::vector<int> > fidx_;

  // Coordinates of degrees freedom on boundary cells.
  // Indexing:  [sideset number][local side id](cell number, value at dof)
  std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > bdryCellDofValues_;

  // Field pattern, offsets, etc.
  std::vector<std::vector<int> > fieldPattern_;  // local Field/DOF pattern; set from DOF manager 
  int numFields_;                                // number of fields (equations in the PDE)
  int numDofs_;                                  // total number of degrees of freedom for all (local) fields
  std::vector<int> offset_;                      // for each field, a counting offset
  std::vector<int> numFieldDofs_;                // for each field, number of degrees of freedom

  ROL::Ptr<FieldHelper<Real> > fieldHelper_;

  int L_;
  Real a_;
  std::vector<Real> uncertain_basis_grid_;

  Real T_;
  Real rho_;
  Real g_;
  Real A_;
  Real width_;
  Real height_;

  int num_coeff_load_;
  std::vector<Real> basal_coeff_;
  std::vector<Real> forcing_coeff_;
  int N_;
  Real h_;

  ROL::Ptr<Intrepid::FieldContainer<Real> > Z_input;

public:
  PDE_shallow_ice_param(Teuchos::ParameterList &parlist) {
    // Finite element fields -- NOT DIMENSION INDEPENDENT!
    basisPtr_ = ROL::makePtr<Intrepid::Basis_HGRAD_QUAD_C1_FEM<Real, Intrepid::FieldContainer<Real> >>();
    // Volume quadrature rules.
    shards::CellTopology cellType = basisPtr_->getBaseCellTopology();         // get the cell type from any basis
    Intrepid::DefaultCubatureFactory<Real> cubFactory;                           // create cubature factory
    int cubDegree = parlist.sublist("Problem").get("Cubature Degree", 2);        // set cubature degree, e.g., 2
    cellCub_ = cubFactory.create(cellType, cubDegree);                           // create default cubature
    // Boundary quadrature rules.
    int d = cellType.getDimension();
    shards::CellTopology bdryCellType = cellType.getCellTopologyData(d-1, 0);
    int bdryCubDegree = parlist.sublist("Problem").get("Boundary Cubature Degree",2); // set cubature degree, e.g., 2
    bdryCub_ = cubFactory.create(bdryCellType, bdryCubDegree);
    basisPtrs_.clear(); basisPtrs_.push_back(basisPtr_); basisPtrs_.push_back(basisPtr_); basisPtrs_.push_back(basisPtr_);

    numDofs_ = 0;
    numFields_ = basisPtrs_.size();
    offset_.resize(numFields_);
    numFieldDofs_.resize(numFields_);
    for (int i=0; i<numFields_; ++i) {
      if (i==0) {
        offset_[i] = 0;
      }
      else {
        offset_[i] = offset_[i-1] + basisPtrs_[i-1]->getCardinality();
      }
      numFieldDofs_[i] = basisPtrs_[i]->getCardinality();
      numDofs_ += numFieldDofs_[i];
    }

    T_  = parlist.sublist("Time Discretization").get("End Time",1.0);

    L_ = parlist.sublist("Problem").get("Number of Uncertain Basis Functions", 10);    
    a_ = parlist.sublist("Problem").get("Noise Level", .2);
    uncertain_basis_grid_.resize(L_+1);
    for(int i = 0; i < L_+1; i++)
      {
	uncertain_basis_grid_[i] = static_cast<Real>(i)/static_cast<Real>(L_);
      }

    width_ = parlist.sublist("Geometry").get("Width",1.0);
    height_ = parlist.sublist("Geometry").get("Height",1.0);

    Real A_scale = parlist.sublist("Problem").get("Thermal Scaling", 1.0);

    rho_ = 910.0;
    g_ = 9.81;
    A_ = std::pow(10.0,-16.0)*A_scale;

    num_coeff_load_ = parlist.sublist("Problem").get("Number of Coefficients in Loaded Fields", 10); 
    N_ = std::sqrt(num_coeff_load_)-1;
    h_ = 1.0/static_cast<Real>(N_);

    Real beta_scale = parlist.sublist("Problem").get("Sliding Scale", 1.0);
    basal_coeff_.resize(num_coeff_load_);
    // read in data
    std::ifstream in("Log_Basal_Sliding.txt");          
    // read the elements in the file into a vector  
    // test file open   
    if (in) 
      {   
	for(int j = 0; j < num_coeff_load_; j++)
	  {
	    in >> basal_coeff_[j];
	  }
      }
    else
      {
	std::cout << "Error loading the data from Log_Basal_Sliding.txt" << std::endl;
      }   
    for(int j = 0; j < num_coeff_load_; j++)
      {
	basal_coeff_[j] += std::log(beta_scale);
      }

    forcing_coeff_.resize(num_coeff_load_);
    // read in data
    std::ifstream in_2("Forcing.txt");       
    // read the elements in the file into a vector  
    // test file open   
    if (in_2) 
      {   
	for(int j = 0; j < num_coeff_load_; j++)
	  {
	    in_2 >> forcing_coeff_[j];
	  }
      }
    else
      {
	std::cout << "Error loading the data from Forcing.txt" << std::endl;
      }     
    
  }

   void residual(ROL::Ptr<Intrepid::FieldContainer<Real> > & res,
                const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    
     // Retrieve dimensions.
    int c  = u_coeff->dimension(0);
    int p  = cellCub_->getNumPoints();
    int f = basisPtr_->getCardinality();
    int d = fe_->gradN()->dimension(3);
 
    // Initialize residuals.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > R(3);
    R[0]   = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, f);
    R[1] = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, f);
    R[2] = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, f);

    // Split u_coeff into components.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > U;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > Z;
    fieldHelper_->splitFieldCoeff(Z, Z_input);
    //fieldHelper_->splitFieldCoeff(Z, z_coeff);

    PDE<Real>::setParameter(*z_param);

    ROL::Ptr<Intrepid::FieldContainer<Real> > valBeta_eval
      = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    ROL::Ptr<Intrepid::FieldContainer<Real> > valHflux_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    compute_parameters(valBeta_eval,valHflux_eval);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        (*valBeta_eval)(i,j) = std::exp((*valBeta_eval)(i,j));
      }
    }

    ROL::Ptr<Intrepid::FieldContainer<Real> > valS_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    fe_->evaluateValue(valS_eval, U[0]);

    ROL::Ptr<Intrepid::FieldContainer<Real> > valB_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    fe_->evaluateValue(valB_eval, Z[0]);

    ROL::Ptr<Intrepid::FieldContainer<Real> > gradS_eval
      = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p, d);
    fe_->evaluateGradient(gradS_eval, U[0]);

    // COMPUTE BASAL VELOCITY
    ROL::Ptr<Intrepid::FieldContainer<Real> > kappa
      = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    computeBasalVelocity(kappa,valBeta_eval,valS_eval,valB_eval,0);

    // MULTIPLY kappa * grad(S)
    Intrepid::FieldContainer<Real> kappa_gradS(c, p, d);
    Intrepid::FunctionSpaceTools::tensorMultiplyDataData<Real>(kappa_gradS,
                                                               *kappa,
                                                               *gradS_eval);
    // INTEGRATE (kappa * grad(S)) . grad(N)
    Intrepid::FunctionSpaceTools::integrate<Real>(*R[0],
                                                  kappa_gradS,
                                                  *(fe_->gradNdetJ()),
                                                  Intrepid::COMP_CPP, false);

    // COMPUTE THERMAL VELOCITY
    Real A = A_;
    ROL::Ptr<Intrepid::FieldContainer<Real> > h
      = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    computeThermalVelocity(h,valS_eval,valB_eval,A,0);
    
    // MULTIPLY grad(S) . grad(S)
    Intrepid::FieldContainer<Real> gradS_gradS(c, p);
    Intrepid::FunctionSpaceTools::dotMultiplyDataData<Real>(gradS_gradS,
							    *gradS_eval,
							    *gradS_eval);
    // MULTIPLY h*grad(S) . grad(S)
    Intrepid::FieldContainer<Real> h_gradS_gradS(c, p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(h_gradS_gradS,
                                                               *h,
                                                               gradS_gradS);

   // MULTIPLY h * grad(S) . grad(S) * grad(S)
    Intrepid::FieldContainer<Real> h_gradS_gradS_gradS(c, p, d);
    Intrepid::FunctionSpaceTools::tensorMultiplyDataData<Real>(h_gradS_gradS_gradS,
                                                               h_gradS_gradS,
                                                               *gradS_eval);
    // INTEGRATE h * grad(S) . grad(S) * grad(S) . grad(N)
    Intrepid::FunctionSpaceTools::integrate<Real>(*R[0],
                                                  h_gradS_gradS_gradS,
                                                  *(fe_->gradNdetJ()),
                                                  Intrepid::COMP_CPP, true);


    Intrepid::RealSpaceTools<Real>::scale(*valHflux_eval,-1.0);
    Intrepid::FunctionSpaceTools::integrate<Real>(*R[0], *valHflux_eval, *(fe_->NdetJ()), Intrepid::COMP_CPP, true); 

    // Surface velocity terms

    ROL::Ptr<Intrepid::FieldContainer<Real> > v = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    computeSurfaceVelocity(v,valS_eval,valB_eval,A_,0);
    ROL::Ptr<Intrepid::FieldContainer<Real> > Sx = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    ROL::Ptr<Intrepid::FieldContainer<Real> > Sy = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    for(int i = 0; i < c; i++)
      {
        for(int j = 0; j < p; j++)
          {
            (*Sx)(i,j) = (*gradS_eval)(i,j,0);
            (*Sy)(i,j) = (*gradS_eval)(i,j,1);
          }
      }

    // MULTIPLY v * (grad(S) . grad(S))
    Intrepid::FieldContainer<Real> v_gradS_gradS(c, p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(v_gradS_gradS,
                                                               *v,
                                                               gradS_gradS);    
   // MULTIPLY v * (grad(S) . grad(S)) * Sx
    Intrepid::FieldContainer<Real> v_gradS_gradS_Sx(c, p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(v_gradS_gradS_Sx,
                                                               *Sx,
                                                               v_gradS_gradS); 

    Intrepid::FunctionSpaceTools::integrate<Real>(*R[1], v_gradS_gradS_Sx, *(fe_->NdetJ()), Intrepid::COMP_CPP, false);

    // MULTIPLY v * (grad(S) . grad(S)) * Sy
    Intrepid::FieldContainer<Real> v_gradS_gradS_Sy(c, p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(v_gradS_gradS_Sy,
                                                               *Sy,
                                                               v_gradS_gradS); 

    Intrepid::FunctionSpaceTools::integrate<Real>(*R[2], v_gradS_gradS_Sy, *(fe_->NdetJ()), Intrepid::COMP_CPP, false);

    // Combine the residuals.
    fieldHelper_->combineFieldCoeff(res, R);
   }

  void Jacobian_1(ROL::Ptr<Intrepid::FieldContainer<Real> > & jac,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // Retrieve dimensions.
    int c  = u_coeff->dimension(0);
    int p  = cellCub_->getNumPoints();
    int f = basisPtr_->getCardinality();
    int d = fe_->gradN()->dimension(3);
    
    // INITILAIZE JACOBIAN
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > J(3);
    for (int i = 0; i < 3; ++i) {
      J[i].resize(3);
      for (int j = 0; j < 3; j++)
	{
	  J[i][j] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
	}
    }

    // Split u_coeff into components.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > U;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > Z;
    fieldHelper_->splitFieldCoeff(Z, Z_input);
    //fieldHelper_->splitFieldCoeff(Z, z_coeff);

    PDE<Real>::setParameter(*z_param);

    ROL::Ptr<Intrepid::FieldContainer<Real> > valBeta_eval
      = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    ROL::Ptr<Intrepid::FieldContainer<Real> > valHflux_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    compute_parameters(valBeta_eval,valHflux_eval);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        (*valBeta_eval)(i,j) = std::exp((*valBeta_eval)(i,j));
      }
    }

    ROL::Ptr<Intrepid::FieldContainer<Real> > valS_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    fe_->evaluateValue(valS_eval, U[0]);

    ROL::Ptr<Intrepid::FieldContainer<Real> > valB_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    fe_->evaluateValue(valB_eval, Z[0]);

    ROL::Ptr<Intrepid::FieldContainer<Real> > gradS_eval
      = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p, d);
    fe_->evaluateGradient(gradS_eval, U[0]);

    // COMPUTE BASAL VELOCITY
    ROL::Ptr<Intrepid::FieldContainer<Real> > kappa
      = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    computeBasalVelocity(kappa,valBeta_eval,valS_eval,valB_eval,0);
    ROL::Ptr<Intrepid::FieldContainer<Real> > d_kappa
      = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    computeBasalVelocity(d_kappa,valBeta_eval,valS_eval,valB_eval,1);
    // MULTIPLY kappa * grad(N)
    Intrepid::FieldContainer<Real> kappa_gradN(c, f, p, d);
    Intrepid::FunctionSpaceTools::tensorMultiplyDataField<Real>(kappa_gradN,
                                                                *kappa,
                                                                *(fe_->gradN()));
    // INTEGRATE (kappa * grad(N)) . grad(N)
    Intrepid::FunctionSpaceTools::integrate<Real>(*J[0][0],
                                                  kappa_gradN,
                                                  *(fe_->gradNdetJ()),
                                                  Intrepid::COMP_CPP, false);
    // MULTIPLY d_kappa * grad(S)
    Intrepid::FieldContainer<Real> d_kappa_gradS(c, p, d);
    Intrepid::FunctionSpaceTools::tensorMultiplyDataData<Real>(d_kappa_gradS,
                                                               *d_kappa,
                                                               *gradS_eval);
    // MULTIPLY (d_kappa * grad(S)) . grad(N)
    Intrepid::FieldContainer<Real> d_kappa_gradS_gradN(c, f, p);
    Intrepid::FunctionSpaceTools::dotMultiplyDataField<Real>(d_kappa_gradS_gradN,
                                                             d_kappa_gradS,
                                                             *(fe_->gradNdetJ()));
    // INTEGRATE (d_kappa * grad(S)) . grad(N) * N
    Intrepid::FunctionSpaceTools::integrate<Real>(*J[0][0],
                                                  d_kappa_gradS_gradN,
                                                  *(fe_->N()),
                                                  Intrepid::COMP_CPP, true);

    // COMPUTE THERMAL VELOCITY
    Real A = A_;
    ROL::Ptr<Intrepid::FieldContainer<Real> > h
      = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    computeThermalVelocity(h,valS_eval,valB_eval,A,0);
    ROL::Ptr<Intrepid::FieldContainer<Real> > d_h
      = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    computeThermalVelocity(d_h,valS_eval,valB_eval,A,1);

    // MULTIPLY grad(S) . grad(S)
    Intrepid::FieldContainer<Real> gradS_gradS(c, p);
    Intrepid::FunctionSpaceTools::dotMultiplyDataData<Real>(gradS_gradS,
							    *gradS_eval,
							    *gradS_eval);
    // MULTIPLY d_h*grad(S) . grad(S)
    Intrepid::FieldContainer<Real> dh_gradS_gradS(c, p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(dh_gradS_gradS,
                                                               *d_h,
                                                               gradS_gradS);

    // MULTIPLY dh * grad(S) . grad(S) * grad(S)
    Intrepid::FieldContainer<Real> dh_gradS_gradS_gradS(c, p, d);
    Intrepid::FunctionSpaceTools::tensorMultiplyDataData<Real>(dh_gradS_gradS_gradS,
                                                               dh_gradS_gradS,
                                                               *gradS_eval);
    // MULTIPLY dh * grad(S) . grad(S) * grad(S) . grad(N)
    Intrepid::FieldContainer<Real> dh_gradS_gradS_gradS_gradN(c, f, p);
    Intrepid::FunctionSpaceTools::dotMultiplyDataField<Real>(dh_gradS_gradS_gradS_gradN,
							     dh_gradS_gradS_gradS,
							     *fe_->gradNdetJ());
    // INTEGRATE dh * grad(S) . grad(S) * grad(S) . grad(N) * N
    Intrepid::FunctionSpaceTools::integrate<Real>(*J[0][0],
                                                  dh_gradS_gradS_gradS_gradN,
                                                  *(fe_->N()),
                                                  Intrepid::COMP_CPP, true);

    // MULTIPLY grad(S) . grad(N)
    Intrepid::FieldContainer<Real> gradS_gradN(c, f, p);
    Intrepid::FunctionSpaceTools::dotMultiplyDataField<Real>(gradS_gradN,
    							    *gradS_eval,
    							     *(fe_->gradN()));
    // MULTIPLY h*grad(S) . grad(N)
    Intrepid::FieldContainer<Real> h_gradS_gradN(c, f, p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(h_gradS_gradN,
    								*h,
    								gradS_gradN);

    // MULTIPLY h * grad(S) . grad(N) * grad(S)
    Intrepid::FieldContainer<Real> h_gradS_gradN_gradS(c, f, p, d);
    for(int i = 0; i < c; i++)
      {
    	for(int j = 0; j < p; j++)
    	  {
    	    for(int k = 0; k < f; k++)
    	      {
    		h_gradS_gradN_gradS(i,k,j,0) = h_gradS_gradN(i,k,j)*(*gradS_eval)(i,j,0);
    		h_gradS_gradN_gradS(i,k,j,1) = h_gradS_gradN(i,k,j)*(*gradS_eval)(i,j,1);
    	      }
    	  }
      }
        						        
    // MULTIPLY 2.0 * h * grad(S) . grad(N) * grad(S)
    Intrepid::RealSpaceTools<Real>::scale(h_gradS_gradN_gradS,2.0);

    // INTEGRATE 2.0 * h * grad(S) . grad(N) * grad(S) . grad(N)
    Intrepid::FunctionSpaceTools::integrate<Real>(*J[0][0],
                                                  h_gradS_gradN_gradS,
                                                  *(fe_->gradNdetJ()),
                                                  Intrepid::COMP_CPP, true);

    // MULTIPLY h*grad(S) . grad(S)
    Intrepid::FieldContainer<Real> h_gradS_gradS(c, p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(h_gradS_gradS,
                                                               *h,
                                                               gradS_gradS);

   // MULTIPLY h * grad(S) . grad(S) * grad(N)
    Intrepid::FieldContainer<Real> h_gradS_gradS_gradN(c, f, p, d);
    Intrepid::FunctionSpaceTools::tensorMultiplyDataField<Real>(h_gradS_gradS_gradN,
								h_gradS_gradS,
								*(fe_->gradN()));
    // INTEGRATE h * grad(S) . grad(S) * grad(N) . grad(N)
    Intrepid::FunctionSpaceTools::integrate<Real>(*J[0][0],
                                                  h_gradS_gradS_gradN,
                                                  *(fe_->gradNdetJ()),
                                                  Intrepid::COMP_CPP, true);

    // Surface velocity

    ROL::Ptr<Intrepid::FieldContainer<Real> > Sx = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    ROL::Ptr<Intrepid::FieldContainer<Real> > Sy = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    for(int i = 0; i < c; i++)
      {
        for(int j = 0; j < p; j++)
          {
            (*Sx)(i,j) = (*gradS_eval)(i,j,0);
            (*Sy)(i,j) = (*gradS_eval)(i,j,1);
          }
      }

    ROL::Ptr<Intrepid::FieldContainer<Real> > v = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    computeSurfaceVelocity(v,valS_eval,valB_eval,A,0);
    ROL::Ptr<Intrepid::FieldContainer<Real> > d1_v = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    computeSurfaceVelocity(d1_v,valS_eval,valB_eval,A,1);
    
    // J[1][0]

    // d1_v * (grad(S) . grad(S)) * Sx * N * N                                                                                                                                                                   

    // MULTIPLY d1_v * (grad(S) . grad(S))                                                                                                                                                                     

    Intrepid::FieldContainer<Real> dv1_gradS_gradS(c, p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(dv1_gradS_gradS,
                                                               *d1_v,
                                                               gradS_gradS);

    // MULTIPLY d1_v * (grad(S) . grad(S)) * Sx                                                                              
    Intrepid::FieldContainer<Real> dv1_gradS_gradS_Sx(c, p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(dv1_gradS_gradS_Sx,
                                                               dv1_gradS_gradS,
                                                               *Sx);

    // MULTIPLY d1_v * (grad(S) . grad(S)) * Sx * N                                                                                                                                 
    Intrepid::FieldContainer<Real> dv1_gradS_gradS_Sx_N(c, f, p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(dv1_gradS_gradS_Sx_N,
                                                                dv1_gradS_gradS_Sx,
                                                                *(fe_->N()));

    // INTEGRATE d1_v * (grad(S) . grad(S)) * Sx * N * N                                                                                                                               
    Intrepid::FunctionSpaceTools::integrate<Real>(*J[1][0],
                                                  dv1_gradS_gradS_Sx_N,
                                                  *(fe_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);

    ////////////////////////////////////////////////////////////////////////////////////////////////

    // v * 2.0 * (grad(S) . grad(N)) * Sx * N

    // MULTIPLY v * Sx                                                                                                                                                                                           
    Intrepid::FieldContainer<Real> v_Sx(c, p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(v_Sx,
                                                               *v,
                                                               *Sx);

    // MULTIPLY 2.0 * v * Sx
    Intrepid::RealSpaceTools<Real>::scale(v_Sx,2.0);

    // MULTIPLY 2.0 * v * Sx * grad(S)                                                                                                                                                                      
    Intrepid::FieldContainer<Real> v_Sx_gradS(c, p, d);
    Intrepid::FunctionSpaceTools::tensorMultiplyDataData<Real>(v_Sx_gradS,
                                                               v_Sx,
                                                               *gradS_eval);

    // MULTIPLY 2.0 * v * Sx * (grad(S) . grad(N))                                                                                                                                                          
    Intrepid::FieldContainer<Real> v_Sx_gradS_gradN(c, f, p);
    Intrepid::FunctionSpaceTools::dotMultiplyDataField<Real>(v_Sx_gradS_gradN,
    							     v_Sx_gradS,
    							     *(fe_->gradNdetJ()));

    // INTEGRATE 2.0 * v * Sx * (grad(S) . grad(N)) * N                                                                                                                                                   
    Intrepid::FunctionSpaceTools::integrate<Real>(*J[1][0],
    						  *(fe_->N()),
                                                  v_Sx_gradS_gradN,
                                                  Intrepid::COMP_CPP, true);      
                                                                                                                         
    //////////////////////////////////////////////////////////////////////////////////////////////////

    // v * (grad(S) . grad(S)) * Nx * N

    // MULTIPLY v * (grad(S) . grad(S))                                                                                                                                                                   
    Intrepid::FieldContainer<Real> v_gradS_gradS(c, p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(v_gradS_gradS,
                                                               *v,
                                                               gradS_gradS);

    // MULTIPLY v * (grad(S) . grad(S)) * Nx                                                                                                                                                               
    Intrepid::FieldContainer<Real> v_gradS_gradS_Nx(c, f, p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(v_gradS_gradS_Nx,
    								v_gradS_gradS,
    								*(fe_->DNDdetJ(0)));

    // INTEGRATE v * (grad(S) . grad(S)) * Nx * N                                                                                                                                                        
    Intrepid::FunctionSpaceTools::integrate<Real>(*J[1][0],
    						  *(fe_->N()),
                                                  v_gradS_gradS_Nx,
    						  Intrepid::COMP_CPP, true);   
    
    ////////////////////////////////////////////////////////////////////////////

    // J[2][0]

    // d1_v * (grad(S) . grad(S)) * Sy * N * N                                                                                                                                                                   

    // MULTIPLY d1_v * (grad(S) . grad(S)) * Sy                                                                              
    Intrepid::FieldContainer<Real> dv1_gradS_gradS_Sy(c, p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(dv1_gradS_gradS_Sy,
                                                               dv1_gradS_gradS,
                                                               *Sy);

    // MULTIPLY d1_v * (grad(S) . grad(S)) * Sy * N                                                                                                                                 
    Intrepid::FieldContainer<Real> dv1_gradS_gradS_Sy_N(c, f, p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(dv1_gradS_gradS_Sy_N,
                                                                dv1_gradS_gradS_Sy,
                                                                *(fe_->N()));

    // INTEGRATE d1_v * (grad(S) . grad(S)) * Sy * N * N                                                                                                                               
    Intrepid::FunctionSpaceTools::integrate<Real>(*J[2][0],
                                                  dv1_gradS_gradS_Sy_N,
                                                  *(fe_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);

    ////////////////////////////////////////////////////////////////////////////////////////////////

    // v * 2.0 * (grad(S) . grad(N)) * Sy * N

    // MULTIPLY v * Sy                                                                                                                                                                                           
    Intrepid::FieldContainer<Real> v_Sy(c, p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(v_Sy,
                                                               *v,
                                                               *Sy);

    // MULTIPLY 2.0 * v * Sy
    Intrepid::RealSpaceTools<Real>::scale(v_Sy,2.0);

    // MULTIPLY 2.0 * v * Sy * grad(S)                                                                                                                                                                      
    Intrepid::FieldContainer<Real> v_Sy_gradS(c, p, d);
    Intrepid::FunctionSpaceTools::tensorMultiplyDataData<Real>(v_Sy_gradS,
                                                               v_Sy,
                                                               *gradS_eval);

    // MULTIPLY 2.0 * v * Sy * (grad(S) . grad(N))                                                                                                                                                          
    Intrepid::FieldContainer<Real> v_Sy_gradS_gradN(c, f, p);
    Intrepid::FunctionSpaceTools::dotMultiplyDataField<Real>(v_Sy_gradS_gradN,
    							     v_Sy_gradS,
    							     *(fe_->gradNdetJ()));

    // INTEGRATE 2.0 * v * Sy * (grad(S) . grad(N)) * N                                                                                                                                                   
    Intrepid::FunctionSpaceTools::integrate<Real>(*J[2][0],
    						  *(fe_->N()),
                                                  v_Sy_gradS_gradN,
                                                  Intrepid::COMP_CPP, true);   
    
                                                                                                                         
    //////////////////////////////////////////////////////////////////////////////////////////////////

    // v * (grad(S) . grad(S)) * Ny * N

    // MULTIPLY v * (grad(S) . grad(S)) * Ny                                                                                                                                                               
    Intrepid::FieldContainer<Real> v_gradS_gradS_Ny(c, f, p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(v_gradS_gradS_Ny,
    								v_gradS_gradS,
    								*(fe_->DNDdetJ(1)));

    // INTEGRATE v * (grad(S) . grad(S)) * Ny * N                                                                                                                                                        
    Intrepid::FunctionSpaceTools::integrate<Real>(*J[2][0],
    						  *(fe_->N()),
                                                  v_gradS_gradS_Ny,
    						  Intrepid::COMP_CPP, true);   
    
    ////////////////////////////////////////////////////////////////////////////


    // Combine the jacobians.
    fieldHelper_->combineFieldCoeff(jac, J);
  }

  void Jacobian_2(ROL::Ptr<Intrepid::FieldContainer<Real> > & jac,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_shallow_ice::Jacobian_2): Jacobian is zero.");
  }
  
  void Jacobian_3(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > & jac,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // GET DIMENSIONS
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
    int p = fe_->gradN()->dimension(2);
    int d = fe_->gradN()->dimension(3);
    // ADD CONTROL TERM TO RESIDUAL
    int size = z_param->size();
    int size_field = size/2;

    PDE<Real>::setParameter(*z_param);

    // Split u_coeff into components.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > U;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > Z;
    fieldHelper_->splitFieldCoeff(Z, Z_input);

    ROL::Ptr<Intrepid::FieldContainer<Real> > val_expBeta_eval
      = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    ROL::Ptr<Intrepid::FieldContainer<Real> > valHflux_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    compute_parameters(val_expBeta_eval,valHflux_eval);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        (*val_expBeta_eval)(i,j) = std::exp((*val_expBeta_eval)(i,j));
      }
    }

    ROL::Ptr<Intrepid::FieldContainer<Real> > valBeta_eval_nom
      = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    ROL::Ptr<Intrepid::FieldContainer<Real> > valHflux_eval_nom =
      ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    compute_parameters(valBeta_eval_nom,valHflux_eval_nom,true);

    Intrepid::RealSpaceTools<Real>::scale(*valHflux_eval_nom,-1.0);

    ROL::Ptr<Intrepid::FieldContainer<Real> > valS_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    fe_->evaluateValue(valS_eval, U[0]);

    ROL::Ptr<Intrepid::FieldContainer<Real> > valB_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    fe_->evaluateValue(valB_eval, Z[0]);

    ROL::Ptr<Intrepid::FieldContainer<Real> > gradS_eval
      = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p, d);
    fe_->evaluateGradient(gradS_eval, U[0]);

    // COMPUTE BASAL VELOCITY
    ROL::Ptr<Intrepid::FieldContainer<Real> > kappa
      = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    computeBasalVelocity(kappa,val_expBeta_eval,valS_eval,valB_eval,0);

    ROL::Ptr<Intrepid::FieldContainer<Real> > beta_kappa_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*beta_kappa_eval,
							       *valBeta_eval_nom,
							       *kappa);

    // Basal sliding parameters
    for (int i = 0; i < size_field; ++i) 
      {
	std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > J(3);
	J[0] = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, f);
	J[1] = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, f);
	J[2] = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, f);

	ROL::Ptr<Intrepid::FieldContainer<Real> > phi_i_eval =
	  ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
	compute_phi_k(phi_i_eval, i);
	
	ROL::Ptr<Intrepid::FieldContainer<Real> > beta_kappa_phi_i_eval =
	  ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
	Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*beta_kappa_phi_i_eval,
								   *beta_kappa_eval,
								   *phi_i_eval);

	ROL::Ptr<Intrepid::FieldContainer<Real> > beta_kappa_phi_i_eval_gradS =
	  ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p, d);
	Intrepid::FunctionSpaceTools::tensorMultiplyDataData<Real>(*beta_kappa_phi_i_eval_gradS,
								   *beta_kappa_phi_i_eval,
								   *gradS_eval);
	Intrepid::FunctionSpaceTools::integrate<Real>(*J[0],
						      *beta_kappa_phi_i_eval_gradS,
						      *(fe_->gradNdetJ()),
						      Intrepid::COMP_CPP, false);
	
	fieldHelper_->combineFieldCoeff(jac[i], J);
      }

    // Forcing term parameters
    for (int i = 0; i < size_field; ++i) 
      {
	std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > J(3);
	J[0] = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, f);
	J[1] = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, f);
	J[2] = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, f);

	ROL::Ptr<Intrepid::FieldContainer<Real> > phi_i_eval =
	  ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
	compute_phi_k(phi_i_eval, i);
	
	ROL::Ptr<Intrepid::FieldContainer<Real> > hflux_phi_i_eval =
	  ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
	Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*hflux_phi_i_eval,
								   *valHflux_eval_nom,
								   *phi_i_eval);
	Intrepid::FunctionSpaceTools::integrate<Real>(*J[0],
						      *hflux_phi_i_eval,
						      *(fe_->NdetJ()),
						      Intrepid::COMP_CPP, false);
	
	fieldHelper_->combineFieldCoeff(jac[size_field+i], J);
      }
  }

  void Hessian_11(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_shallow_ice::Hessian_11): Hessian is zero.");
  }

  void Hessian_12(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_shallow_ice::Hessian_12): Hessian is zero.");
  }

  void Hessian_13(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_shallow_ice::Hessian_13): Hessian is zero.");
  }

  void Hessian_21(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_shallow_ice::Hessian_21): Hessian is zero.");
  }

  void Hessian_22(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_shallow_ice::Hessian_22): Hessian is zero.");
  }

  void Hessian_23(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_shallow_ice::Hessian_23): Hessian is zero.");
  }

  void Hessian_31(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_shallow_ice::Hessian_31): Hessian is zero.");
  }

  void Hessian_32(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_shallow_ice::Hessian_32): Hessian is zero.");
  }

  void Hessian_33(std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_shallow_ice::Hessian_33): Hessian is zero.");
  }

   void RieszMap_1(ROL::Ptr<Intrepid::FieldContainer<Real> > & riesz) {
     throw Exception::NotImplemented(">>> (PDE_shallow_ice::RieszMap_1): Not implemented.");
  }

  void RieszMap_2(ROL::Ptr<Intrepid::FieldContainer<Real> > & riesz) {
    throw Exception::NotImplemented(">>> (PDE_shallow_ice::RieszMap_2): Not implemented.");
  }

  void setFieldPattern(const std::vector<std::vector<int> > & fieldPattern) {
    fieldPattern_ = fieldPattern;
    fieldHelper_ = ROL::makePtr<FieldHelper<Real>>(numFields_, numDofs_, numFieldDofs_, fieldPattern_);
  }

 void setCellNodes(const ROL::Ptr<Intrepid::FieldContainer<Real> > &volCellNodes,
                    const std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > &bdryCellNodes,
                    const std::vector<std::vector<std::vector<int> > > &bdryCellLocIds) {
    volCellNodes_ = volCellNodes;
    bdryCellNodes_ = bdryCellNodes;
    bdryCellLocIds_ = bdryCellLocIds;
    // Finite element definition.
    fe_ = ROL::makePtr<FE<Real>>(volCellNodes_,basisPtr_,cellCub_);
    // Set local boundary DOFs.
    fidx_ = fe_->getBoundaryDofs();
    // Compute Dirichlet values at DOFs.
    int d = basisPtr_->getBaseCellTopology().getDimension();
    int numSidesets = bdryCellLocIds_.size();
    bdryCellDofValues_.resize(numSidesets);
    feBdry_.resize(numSidesets); 
    for (int i=0; i<numSidesets; ++i) {
      int numLocSides = bdryCellLocIds_[i].size();
      bdryCellDofValues_[i].resize(numLocSides);
      feBdry_[i].resize(numLocSides);
      for (int j=0; j<numLocSides; ++j) {
        int c = bdryCellLocIds_[i][j].size();
        int f = basisPtr_->getCardinality();
	if(c > 0)
	  {
	    feBdry_[i][j] = ROL::makePtr<FE<Real> >(bdryCellNodes_[i][j],basisPtr_,bdryCub_,j);
	  }
	bdryCellDofValues_[i][j] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
        ROL::Ptr<Intrepid::FieldContainer<Real> > coords =
          ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, d);
        if (c > 0) {
          fe_->computeDofCoords(coords, bdryCellNodes_[i][j]);
        }
        for (int k=0; k<c; ++k) {
          for (int l=0; l<f; ++l) {
            std::vector<Real> dofpoint(d);
            for (int m=0; m<d; ++m) {
              dofpoint[m] = (*coords)(k, l, m);
            }
            (*bdryCellDofValues_[i][j])(k, l) = 0.0;
          }
        }
      }
    }

  }

  ROL::Ptr<Intrepid::FieldContainer<Real> > getBoundaryCoeff(const Intrepid::FieldContainer<Real> & cell_coeff, int sideSet, int cell) const {
    std::vector<int> bdryCellLocId = bdryCellLocIds_[sideSet][cell];
    const int numCellsSide = bdryCellLocId.size();
    const int f = basisPtr_->getCardinality();
    
    ROL::Ptr<Intrepid::FieldContainer<Real > > bdry_coeff = 
      ROL::makePtr<Intrepid::FieldContainer<Real > >(numCellsSide, f);
    for (int i = 0; i < numCellsSide; ++i) {
      for (int j = 0; j < f; ++j) {
	(*bdry_coeff)(i, j) = cell_coeff(bdryCellLocId[i], j);
      }
    }
    return bdry_coeff;
  }

  const ROL::Ptr<FE<Real> > getFE(void) const {
    return fe_;
  }

  const ROL::Ptr<FieldHelper<Real> > getFieldHelper(void) const {
    return fieldHelper_;
  }

  std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>>> getFields() {
    return basisPtrs_;
  }

  const std::vector<std::vector<std::vector<int> > > getBdryCellLocIds(void) const {
    return bdryCellLocIds_;
  }

private:
  
  void computeBasalVelocity(ROL::Ptr<Intrepid::FieldContainer<Real> > &kappa,
			    const ROL::Ptr<Intrepid::FieldContainer<Real> > &beta,
			    const ROL::Ptr<Intrepid::FieldContainer<Real> > &u,
			    const ROL::Ptr<Intrepid::FieldContainer<Real> > &z,
			    const int deriv = 0 ) const {
    // GET DIMENSIONS
    int c = fe_->gradN()->dimension(0);
    int p = fe_->gradN()->dimension(2);
    int d = fe_->gradN()->dimension(3);
    std::vector<Real> pt(d);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
	pt[0] = (*fe_->cubPts())(i,j,0)/width_;
	pt[1] = (*fe_->cubPts())(i,j,1)/height_;
        // Compute basal velocity
        (*kappa)(i,j) = evaluateBasalVelocity((*beta)(i,j),(*u)(i,j),(*z)(i,j),pt,deriv);
      }
    }
  }

  Real evaluateBasalVelocity(const Real beta, const Real u, const Real z, const std::vector<Real> & x, const int deriv = 0) const {
    Real val = 0.0;
    if(deriv == 0)
      {
	val = beta*rho_*g_*std::pow(u-z,2.0);
      }
    if(deriv == 1)
      {
	val = beta*rho_*g_*2.0*(u-z);
      }
    if(deriv == 2)
      {
	val = beta*rho_*g_*2.0;
      }
    return val;
  }

  void computeThermalVelocity(ROL::Ptr<Intrepid::FieldContainer<Real> > &h,
			      const ROL::Ptr<Intrepid::FieldContainer<Real> > &u,
			      const ROL::Ptr<Intrepid::FieldContainer<Real> > &z,
			      const Real A,
			      const int deriv = 0 ) const {
    // GET DIMENSIONS
    int c = fe_->gradN()->dimension(0);
    int p = fe_->gradN()->dimension(2);
    int d = fe_->gradN()->dimension(3);
    std::vector<Real> pt(d);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
	pt[0] = (*fe_->cubPts())(i,j,0)/width_;
	pt[1] = (*fe_->cubPts())(i,j,1)/height_;
        // Compute thermal velocity
        (*h)(i,j) = evaluateThermalVelocity((*u)(i,j),(*z)(i,j),pt,A,deriv);
      }
    }
  }

  Real evaluateThermalVelocity(const Real u, const Real z, const std::vector<Real> & x, const Real A, const int deriv = 0) const {
    Real val = 0.0;
    if(deriv == 0)
      {
	val = 0.4*A*std::pow(rho_*g_,3.0)*std::pow(u-z,5.0);
      }
    if(deriv == 1)
      {
	val = 0.4*A*std::pow(rho_*g_,3.0)*5.0*std::pow(u-z,4.0);
      }
    if(deriv == 2)
      {
        val = 0.4*A*std::pow(rho_*g_,3.0)*20.0*std::pow(u-z,3.0);
      }
    return val;
  }

  void computeSurfaceVelocity(ROL::Ptr<Intrepid::FieldContainer<Real> > &v,
			      const ROL::Ptr<Intrepid::FieldContainer<Real> > &u,
			      const ROL::Ptr<Intrepid::FieldContainer<Real> > &z,
			      const Real A,
			      const int deriv = 0 ) const {
    // GET DIMENSIONS
    int c = fe_->gradN()->dimension(0);
    int p = fe_->gradN()->dimension(2);
    int d = fe_->gradN()->dimension(3);
    std::vector<Real> pt(d);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
 	pt[0] = (*fe_->cubPts())(i,j,0)/width_;
	pt[1] = (*fe_->cubPts())(i,j,1)/height_;
        // Compute thermal velocity
        (*v)(i,j) = evaluateSurfaceVelocity((*u)(i,j),(*z)(i,j),pt,A,deriv);
      }
    }
  }

  Real evaluateSurfaceVelocity(const Real u, const Real z, const std::vector<Real> & x, const Real A, const int deriv = 0) const {
    Real val = 0.0;
    if(deriv == 0)
      {
	val = 0.5*A*std::pow(rho_*g_,3.0)*std::pow(u-z,4.0);
      }
    if(deriv == 1)
      {
	val = 0.5*A*std::pow(rho_*g_,3.0)*4.0*std::pow(u-z,3.0);
      }
    if(deriv == 2)
      {
        val = 0.5*A*std::pow(rho_*g_,3.0)*12.0*std::pow(u-z,2.0);
      }
    return 1000.0*val;
  }

  void compute_parameters(ROL::Ptr<Intrepid::FieldContainer<Real> > & beta, ROL::Ptr<Intrepid::FieldContainer<Real> > & h_flux,  bool nominal_param = false) const {
    int c = fe_->gradN()->dimension(0);
    int p = fe_->gradN()->dimension(2);
    int d = fe_->gradN()->dimension(3);
    std::vector<Real> pt(d);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
	pt[0] = (*fe_->cubPts())(i,j,0)/width_;
	pt[1] = (*fe_->cubPts())(i,j,1)/height_;
        (*beta)(i,j) = basal_sliding(pt,nominal_param);
        (*h_flux)(i,j) = rhs(pt,nominal_param);
      }
    }

  }

  Real basal_sliding(const std::vector<Real> & pt, const bool nominal_param = false) const {
    Real val = 0.0;
    Real x = pt[0];
    Real y = pt[1];

    int i = std::floor(x/h_);
    int j = std::floor(y/h_);
    for(int ii = -1; ii < 2; ii++)
      {
	for(int jj = -1; jj < 2; jj++)
	  {
	    int i_index = i+ii;
	    int j_index = j+jj;
	    if(i_index>=0 && j_index>=0 && i_index<N_+1 && j_index<N_+1)
	      {
		int k = j_index*(N_+1)+i_index;
		val += basal_coeff_[k]*Linear_FE_Basis_Fun_Eval(x,i_index)*Linear_FE_Basis_Fun_Eval(y,j_index);
	      }
	  }
      }
    if(!nominal_param)
      {
	const std::vector<Real> param = PDE<Real>::getParameter();
	if((int)param.size()>0)
	  {
	    Real val_uncertain = 0.0;
	    int count = 0;
	    Real basis_fun_x = 0.0;
	    Real basis_fun_y = 0.0;
	    for(int j = 0; j < L_+1; j++)
	      {
	    	basis_fun_y = Parameter_Basis_Fun_Eval(y,j);
	    	for(int i = 0; i < L_+1; i++)
	    	  {
	    	    basis_fun_x = Parameter_Basis_Fun_Eval(x,i);
	    	    val_uncertain += param[count]*basis_fun_x*basis_fun_y;
	    	    count = count+1;
	    	  }
	      }
	    val_uncertain = a_*val_uncertain;
	    val_uncertain = 1+val_uncertain;
	    val = val*val_uncertain;
	  }
      }

    return val;
  }

  Real rhs(const std::vector<Real> & pt, bool & nominal_param = false) const {
    Real val = 0.0;
    Real x = pt[0];
    Real y = pt[1];
    Real t = PDE<Real>::getTime();
    Real pi(M_PI);
    Real time_coeff = 1.0;

    int i = std::floor(x/h_);
    int j = std::floor(y/h_);
    for(int ii = -1; ii < 2; ii++)
      {
	for(int jj = -1; jj < 2; jj++)
	  {
	    int i_index = i+ii;
	    int j_index = j+jj;
	    if(i_index>=0 && j_index>=0 && i_index<N_+1 && j_index<N_+1)
	      {
		int k = j_index*(N_+1)+i_index;
		val += time_coeff*forcing_coeff_[k]*Linear_FE_Basis_Fun_Eval(x,i_index)*Linear_FE_Basis_Fun_Eval(y,j_index);
	      }
	  }
      }
    if(!nominal_param)
      {
	const std::vector<Real> param = PDE<Real>::getParameter();
	if((int)param.size()>0)
	  {
	    Real val_uncertain = 0.0;
	    int count = (L_+1)*(L_+1);
	    Real basis_fun_x = 0.0;
	    Real basis_fun_y = 0.0;
	    for(int j = 0; j < L_+1; j++)
	      {
		basis_fun_y = Parameter_Basis_Fun_Eval(y,j);
		for(int i = 0; i < L_+1; i++)
		  {
		    basis_fun_x = Parameter_Basis_Fun_Eval(x,i);
		    val_uncertain += param[count]*basis_fun_x*basis_fun_y;
		    count = count+1;
		  }
	      }
	    val_uncertain = a_*val_uncertain;
	    val = val*val_uncertain;
	  }
      }

    return val;
  }

  void compute_phi_k(ROL::Ptr<Intrepid::FieldContainer<Real> > & phi, int k) const {
    int k_x = k%(L_+1);
    int k_y = std::floor(k/(L_+1));
    int c = fe_->gradN()->dimension(0);
    int p = fe_->gradN()->dimension(2);
    int d = fe_->gradN()->dimension(3);
    std::vector<Real> pt(d);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
	pt[0] = (*fe_->cubPts())(i,j,0)/width_;
	pt[1] = (*fe_->cubPts())(i,j,1)/height_;
        (*phi)(i,j) = a_*Parameter_Basis_Fun_Eval(pt[0],k_x)*Parameter_Basis_Fun_Eval(pt[1],k_y);
      }
    }

  }

  Real Linear_FE_Basis_Fun_Eval(const Real & x, int i) const
  {
    // Evaluates ith 1D FE basis function at point x
    Real val = 0.0;
    Real xi = static_cast<Real>(i)*h_;
    Real dist = std::abs(x-xi);
    if( dist < h_ )
      {
	val = 1.0 - dist/h_;
      }
    return val;
  }

  Real Parameter_Basis_Fun_Eval(const Real & x, int i) const
  {
    // Evaluates ith 1D basis function on a grid of L+1 nodes at point x
    Real val = 0.0;
    Real dist = std::abs(x-uncertain_basis_grid_[i]);
    if( dist < 1.0/static_cast<Real>(L_) )
      {
	val = 1.0 - dist*static_cast<Real>(L_);
      }
    return val;
  }

public:
  void Update_Z_input(ROL::Ptr<Intrepid::FieldContainer<Real> > & Z)
  {
    Z_input = Z;
  }

}; // PDE_shallow_ice

#endif



// #ifndef PDE_SHALLOW_ICE_PARAM_HPP
// #define PDE_SHALLOW_ICE_PARAM_HPP

// #include "../../../PDE-OPT/TOOLS/pde.hpp"
// #include "../../../PDE-OPT/TOOLS/fe.hpp"

// #include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"
// #include "Intrepid_HGRAD_QUAD_C2_FEM.hpp"
// #include "Intrepid_DefaultCubatureFactory.hpp"
// #include "Intrepid_FunctionSpaceTools.hpp"
// #include "Intrepid_CellTools.hpp"

// #include "ROL_Ptr.hpp"

// template <class Real>
// class PDE_shallow_ice_param : public PDE<Real> {
// private:
//   // Finite element basis information
//   ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > basisPtr_;
//   std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > > basisPtrs_;
//   // Cell cubature information
//   ROL::Ptr<Intrepid::Cubature<Real> > cellCub_;
//   ROL::Ptr<Intrepid::Cubature<Real> > bdryCub_;
//   // Cell node information
//   ROL::Ptr<Intrepid::FieldContainer<Real> > volCellNodes_;
//   std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > bdryCellNodes_;
//   std::vector<std::vector<std::vector<int> > > bdryCellLocIds_;
//   // Finite element definition
//   ROL::Ptr<FE<Real> > fe_;
//   std::vector<std::vector<ROL::Ptr<FE<Real> > > > feBdry_;
//   // Local degrees of freedom on boundary, for each side of the reference cell (first index).
//   std::vector<std::vector<int> > fidx_;

//   // Coordinates of degrees freedom on boundary cells.
//   // Indexing:  [sideset number][local side id](cell number, value at dof)
//   std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > bdryCellDofValues_;

//   // Field pattern, offsets, etc.
//   std::vector<std::vector<int> > fieldPattern_;  // local Field/DOF pattern; set from DOF manager 
//   int numFields_;                                // number of fields (equations in the PDE)
//   int numDofs_;                                  // total number of degrees of freedom for all (local) fields
//   std::vector<int> offset_;                      // for each field, a counting offset
//   std::vector<int> numFieldDofs_;                // for each field, number of degrees of freedom

//   ROL::Ptr<FieldHelper<Real> > fieldHelper_;

//   int L_;
//   Real a_;
//   std::vector<Real> uncertain_basis_grid_;

//   Real T_;
//   Real rho_;
//   Real g_;
//   Real A_;
//   Real width_;
//   Real height_;

//   int num_coeff_load_;
//   std::vector<Real> basal_coeff_;
//   std::vector<Real> forcing_coeff_;
//   int N_;
//   Real h_;

//   ROL::Ptr<Intrepid::FieldContainer<Real> > Z_input;

// public:
//   PDE_shallow_ice_param(Teuchos::ParameterList &parlist) {
//     // Finite element fields -- NOT DIMENSION INDEPENDENT!
//     basisPtr_ = ROL::makePtr<Intrepid::Basis_HGRAD_QUAD_C1_FEM<Real, Intrepid::FieldContainer<Real> >>();
//     // Volume quadrature rules.
//     shards::CellTopology cellType = basisPtr_->getBaseCellTopology();         // get the cell type from any basis
//     Intrepid::DefaultCubatureFactory<Real> cubFactory;                           // create cubature factory
//     int cubDegree = parlist.sublist("Problem").get("Cubature Degree", 2);        // set cubature degree, e.g., 2
//     cellCub_ = cubFactory.create(cellType, cubDegree);                           // create default cubature
//     // Boundary quadrature rules.
//     int d = cellType.getDimension();
//     shards::CellTopology bdryCellType = cellType.getCellTopologyData(d-1, 0);
//     int bdryCubDegree = parlist.sublist("Problem").get("Boundary Cubature Degree",2); // set cubature degree, e.g., 2
//     bdryCub_ = cubFactory.create(bdryCellType, bdryCubDegree);
//     basisPtrs_.clear(); basisPtrs_.push_back(basisPtr_); basisPtrs_.push_back(basisPtr_); basisPtrs_.push_back(basisPtr_);

//     numDofs_ = 0;
//     numFields_ = basisPtrs_.size();
//     offset_.resize(numFields_);
//     numFieldDofs_.resize(numFields_);
//     for (int i=0; i<numFields_; ++i) {
//       if (i==0) {
//         offset_[i] = 0;
//       }
//       else {
//         offset_[i] = offset_[i-1] + basisPtrs_[i-1]->getCardinality();
//       }
//       numFieldDofs_[i] = basisPtrs_[i]->getCardinality();
//       numDofs_ += numFieldDofs_[i];
//     }

//     T_  = parlist.sublist("Time Discretization").get("End Time",1.0);

//     L_ = parlist.sublist("Problem").get("Number of Uncertain Basis Functions", 10);    
//     a_ = parlist.sublist("Problem").get("Noise Level", .2);
//     uncertain_basis_grid_.resize(L_+1);
//     for(int i = 0; i < L_+1; i++)
//       {
// 	uncertain_basis_grid_[i] = static_cast<Real>(i)/static_cast<Real>(L_);
//       }

//     width_ = parlist.sublist("Geometry").get("Width",1.0);
//     height_ = parlist.sublist("Geometry").get("Height",1.0);

//     Real A_scale = parlist.sublist("Problem").get("Thermal Scaling", 1.0);

//     rho_ = 910.0;
//     g_ = 9.81;
//     A_ = std::pow(10.0,-16.0)*A_scale;

//     num_coeff_load_ = parlist.sublist("Problem").get("Number of Coefficients in Loaded Fields", 10); 
//     N_ = std::sqrt(num_coeff_load_)-1;
//     h_ = 1.0/static_cast<Real>(N_);

//     Real beta_scale = parlist.sublist("Problem").get("Sliding Scale", 1.0);
//     basal_coeff_.resize(num_coeff_load_);
//     // read in data
//     std::ifstream in("Log_Basal_Sliding.txt");          
//     // read the elements in the file into a vector  
//     // test file open   
//     if (in) 
//       {   
// 	for(int j = 0; j < num_coeff_load_; j++)
// 	  {
// 	    in >> basal_coeff_[j];
// 	  }
//       }
//     else
//       {
// 	std::cout << "Error loading the data from Log_Basal_Sliding.txt" << std::endl;
//       }   
//     for(int j = 0; j < num_coeff_load_; j++)
//       {
// 	basal_coeff_[j] += std::log(beta_scale);
//       }

//     forcing_coeff_.resize(num_coeff_load_);
//     // read in data
//     std::ifstream in_2("Forcing.txt");       
//     // read the elements in the file into a vector  
//     // test file open   
//     if (in_2) 
//       {   
// 	for(int j = 0; j < num_coeff_load_; j++)
// 	  {
// 	    in_2 >> forcing_coeff_[j];
// 	  }
//       }
//     else
//       {
// 	std::cout << "Error loading the data from Forcing.txt" << std::endl;
//       }     
    
//   }

//    void residual(ROL::Ptr<Intrepid::FieldContainer<Real> > & res,
//                 const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
//                 const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
//                 const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    
//      // Retrieve dimensions.
//     int c  = u_coeff->dimension(0);
//     int p  = cellCub_->getNumPoints();
//     int f = basisPtr_->getCardinality();
//     int d = fe_->gradN()->dimension(3);
 
//     // Initialize residuals.
//     std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > R(3);
//     R[0]   = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, f);
//     R[1] = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, f);
//     R[2] = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, f);

//     // Split u_coeff into components.
//     std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > U;
//     fieldHelper_->splitFieldCoeff(U, u_coeff);
//     std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > Z;
//     fieldHelper_->splitFieldCoeff(Z, Z_input);
//     //fieldHelper_->splitFieldCoeff(Z, z_coeff);

//     PDE<Real>::setParameter(*z_param);

//     ROL::Ptr<Intrepid::FieldContainer<Real> > valBeta_eval
//       = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
//     ROL::Ptr<Intrepid::FieldContainer<Real> > valHflux_eval =
//       ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
//     compute_parameters(valBeta_eval,valHflux_eval);
//     for (int i = 0; i < c; ++i) {
//       for (int j = 0; j < p; ++j) {
//         (*valBeta_eval)(i,j) = std::exp((*valBeta_eval)(i,j));
//       }
//     }

//     ROL::Ptr<Intrepid::FieldContainer<Real> > valS_eval =
//       ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
//     fe_->evaluateValue(valS_eval, U[0]);

//     ROL::Ptr<Intrepid::FieldContainer<Real> > valB_eval =
//       ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
//     fe_->evaluateValue(valB_eval, Z[0]);

//     ROL::Ptr<Intrepid::FieldContainer<Real> > gradS_eval
//       = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p, d);
//     fe_->evaluateGradient(gradS_eval, U[0]);

//     // COMPUTE BASAL VELOCITY
//     ROL::Ptr<Intrepid::FieldContainer<Real> > kappa
//       = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
//     computeBasalVelocity(kappa,valBeta_eval,valS_eval,valB_eval,0);

//     // MULTIPLY kappa * grad(S)
//     Intrepid::FieldContainer<Real> kappa_gradS(c, p, d);
//     Intrepid::FunctionSpaceTools::tensorMultiplyDataData<Real>(kappa_gradS,
//                                                                *kappa,
//                                                                *gradS_eval);
//     // INTEGRATE (kappa * grad(S)) . grad(N)
//     Intrepid::FunctionSpaceTools::integrate<Real>(*R[0],
//                                                   kappa_gradS,
//                                                   *(fe_->gradNdetJ()),
//                                                   Intrepid::COMP_CPP, false);

//     // COMPUTE THERMAL VELOCITY
//     Real A = A_;
//     const std::vector<Real> param = PDE<Real>::getParameter();
//     if((int)param.size()>0)
//       {
//     	A = A*(1.0+a_*param[2*(L_+1)*(L_+1)]);
//       }
//     ROL::Ptr<Intrepid::FieldContainer<Real> > h
//       = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
//     computeThermalVelocity(h,valS_eval,valB_eval,A,0);
    
//     // MULTIPLY grad(S) . grad(S)
//     Intrepid::FieldContainer<Real> gradS_gradS(c, p);
//     Intrepid::FunctionSpaceTools::dotMultiplyDataData<Real>(gradS_gradS,
// 							    *gradS_eval,
// 							    *gradS_eval);
//     // MULTIPLY h*grad(S) . grad(S)
//     Intrepid::FieldContainer<Real> h_gradS_gradS(c, p);
//     Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(h_gradS_gradS,
//                                                                *h,
//                                                                gradS_gradS);

//    // MULTIPLY h * grad(S) . grad(S) * grad(S)
//     Intrepid::FieldContainer<Real> h_gradS_gradS_gradS(c, p, d);
//     Intrepid::FunctionSpaceTools::tensorMultiplyDataData<Real>(h_gradS_gradS_gradS,
//                                                                h_gradS_gradS,
//                                                                *gradS_eval);
//     // INTEGRATE h * grad(S) . grad(S) * grad(S) . grad(N)
//     Intrepid::FunctionSpaceTools::integrate<Real>(*R[0],
//                                                   h_gradS_gradS_gradS,
//                                                   *(fe_->gradNdetJ()),
//                                                   Intrepid::COMP_CPP, true);


//     Intrepid::RealSpaceTools<Real>::scale(*valHflux_eval,-1.0);
//     Intrepid::FunctionSpaceTools::integrate<Real>(*R[0], *valHflux_eval, *(fe_->NdetJ()), Intrepid::COMP_CPP, true); 

//     // Surface velocity terms

//     ROL::Ptr<Intrepid::FieldContainer<Real> > v = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
//     computeSurfaceVelocity(v,valS_eval,valB_eval,A_,0);
//     ROL::Ptr<Intrepid::FieldContainer<Real> > Sx = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
//     ROL::Ptr<Intrepid::FieldContainer<Real> > Sy = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
//     for(int i = 0; i < c; i++)
//       {
//         for(int j = 0; j < p; j++)
//           {
//             (*Sx)(i,j) = (*gradS_eval)(i,j,0);
//             (*Sy)(i,j) = (*gradS_eval)(i,j,1);
//           }
//       }

//     // MULTIPLY v * (grad(S) . grad(S))
//     Intrepid::FieldContainer<Real> v_gradS_gradS(c, p);
//     Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(v_gradS_gradS,
//                                                                *v,
//                                                                gradS_gradS);    
//    // MULTIPLY v * (grad(S) . grad(S)) * Sx
//     Intrepid::FieldContainer<Real> v_gradS_gradS_Sx(c, p);
//     Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(v_gradS_gradS_Sx,
//                                                                *Sx,
//                                                                v_gradS_gradS); 

//     Intrepid::FunctionSpaceTools::integrate<Real>(*R[1], v_gradS_gradS_Sx, *(fe_->NdetJ()), Intrepid::COMP_CPP, false);

//     // MULTIPLY v * (grad(S) . grad(S)) * Sy
//     Intrepid::FieldContainer<Real> v_gradS_gradS_Sy(c, p);
//     Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(v_gradS_gradS_Sy,
//                                                                *Sy,
//                                                                v_gradS_gradS); 

//     Intrepid::FunctionSpaceTools::integrate<Real>(*R[2], v_gradS_gradS_Sy, *(fe_->NdetJ()), Intrepid::COMP_CPP, false);

//     // Combine the residuals.
//     fieldHelper_->combineFieldCoeff(res, R);
//    }

//   void Jacobian_1(ROL::Ptr<Intrepid::FieldContainer<Real> > & jac,
//                   const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
//                   const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
//                   const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
//     // Retrieve dimensions.
//     int c  = u_coeff->dimension(0);
//     int p  = cellCub_->getNumPoints();
//     int f = basisPtr_->getCardinality();
//     int d = fe_->gradN()->dimension(3);
    
//     // INITILAIZE JACOBIAN
//     std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > J(3);
//     for (int i = 0; i < 3; ++i) {
//       J[i].resize(3);
//       for (int j = 0; j < 3; j++)
// 	{
// 	  J[i][j] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
// 	}
//     }

//     // Split u_coeff into components.
//     std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > U;
//     fieldHelper_->splitFieldCoeff(U, u_coeff);
//     std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > Z;
//     fieldHelper_->splitFieldCoeff(Z, Z_input);
//     //fieldHelper_->splitFieldCoeff(Z, z_coeff);

//     PDE<Real>::setParameter(*z_param);

//     ROL::Ptr<Intrepid::FieldContainer<Real> > valBeta_eval
//       = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
//     ROL::Ptr<Intrepid::FieldContainer<Real> > valHflux_eval =
//       ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
//     compute_parameters(valBeta_eval,valHflux_eval);
//     for (int i = 0; i < c; ++i) {
//       for (int j = 0; j < p; ++j) {
//         (*valBeta_eval)(i,j) = std::exp((*valBeta_eval)(i,j));
//       }
//     }

//     ROL::Ptr<Intrepid::FieldContainer<Real> > valS_eval =
//       ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
//     fe_->evaluateValue(valS_eval, U[0]);

//     ROL::Ptr<Intrepid::FieldContainer<Real> > valB_eval =
//       ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
//     fe_->evaluateValue(valB_eval, Z[0]);

//     ROL::Ptr<Intrepid::FieldContainer<Real> > gradS_eval
//       = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p, d);
//     fe_->evaluateGradient(gradS_eval, U[0]);

//     // COMPUTE BASAL VELOCITY
//     ROL::Ptr<Intrepid::FieldContainer<Real> > kappa
//       = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
//     computeBasalVelocity(kappa,valBeta_eval,valS_eval,valB_eval,0);
//     ROL::Ptr<Intrepid::FieldContainer<Real> > d_kappa
//       = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
//     computeBasalVelocity(d_kappa,valBeta_eval,valS_eval,valB_eval,1);
//     // MULTIPLY kappa * grad(N)
//     Intrepid::FieldContainer<Real> kappa_gradN(c, f, p, d);
//     Intrepid::FunctionSpaceTools::tensorMultiplyDataField<Real>(kappa_gradN,
//                                                                 *kappa,
//                                                                 *(fe_->gradN()));
//     // INTEGRATE (kappa * grad(N)) . grad(N)
//     Intrepid::FunctionSpaceTools::integrate<Real>(*J[0][0],
//                                                   kappa_gradN,
//                                                   *(fe_->gradNdetJ()),
//                                                   Intrepid::COMP_CPP, false);
//     // MULTIPLY d_kappa * grad(S)
//     Intrepid::FieldContainer<Real> d_kappa_gradS(c, p, d);
//     Intrepid::FunctionSpaceTools::tensorMultiplyDataData<Real>(d_kappa_gradS,
//                                                                *d_kappa,
//                                                                *gradS_eval);
//     // MULTIPLY (d_kappa * grad(S)) . grad(N)
//     Intrepid::FieldContainer<Real> d_kappa_gradS_gradN(c, f, p);
//     Intrepid::FunctionSpaceTools::dotMultiplyDataField<Real>(d_kappa_gradS_gradN,
//                                                              d_kappa_gradS,
//                                                              *(fe_->gradNdetJ()));
//     // INTEGRATE (d_kappa * grad(S)) . grad(N) * N
//     Intrepid::FunctionSpaceTools::integrate<Real>(*J[0][0],
//                                                   d_kappa_gradS_gradN,
//                                                   *(fe_->N()),
//                                                   Intrepid::COMP_CPP, true);

//     // COMPUTE THERMAL VELOCITY
//     Real A = A_;
//     const std::vector<Real> param = PDE<Real>::getParameter();
//     if((int)param.size()>0)
//       {
// 	A = A*(1.0+a_*param[2*(L_+1)*(L_+1)]);
//       }
//     ROL::Ptr<Intrepid::FieldContainer<Real> > h
//       = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
//     computeThermalVelocity(h,valS_eval,valB_eval,A,0);
//     ROL::Ptr<Intrepid::FieldContainer<Real> > d_h
//       = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
//     computeThermalVelocity(d_h,valS_eval,valB_eval,A,1);

//     // MULTIPLY grad(S) . grad(S)
//     Intrepid::FieldContainer<Real> gradS_gradS(c, p);
//     Intrepid::FunctionSpaceTools::dotMultiplyDataData<Real>(gradS_gradS,
// 							    *gradS_eval,
// 							    *gradS_eval);
//     // MULTIPLY d_h*grad(S) . grad(S)
//     Intrepid::FieldContainer<Real> dh_gradS_gradS(c, p);
//     Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(dh_gradS_gradS,
//                                                                *d_h,
//                                                                gradS_gradS);

//     // MULTIPLY dh * grad(S) . grad(S) * grad(S)
//     Intrepid::FieldContainer<Real> dh_gradS_gradS_gradS(c, p, d);
//     Intrepid::FunctionSpaceTools::tensorMultiplyDataData<Real>(dh_gradS_gradS_gradS,
//                                                                dh_gradS_gradS,
//                                                                *gradS_eval);
//     // MULTIPLY dh * grad(S) . grad(S) * grad(S) . grad(N)
//     Intrepid::FieldContainer<Real> dh_gradS_gradS_gradS_gradN(c, f, p);
//     Intrepid::FunctionSpaceTools::dotMultiplyDataField<Real>(dh_gradS_gradS_gradS_gradN,
// 							     dh_gradS_gradS_gradS,
// 							     *fe_->gradNdetJ());
//     // INTEGRATE dh * grad(S) . grad(S) * grad(S) . grad(N) * N
//     Intrepid::FunctionSpaceTools::integrate<Real>(*J[0][0],
//                                                   dh_gradS_gradS_gradS_gradN,
//                                                   *(fe_->N()),
//                                                   Intrepid::COMP_CPP, true);

//     // MULTIPLY grad(S) . grad(N)
//     Intrepid::FieldContainer<Real> gradS_gradN(c, f, p);
//     Intrepid::FunctionSpaceTools::dotMultiplyDataField<Real>(gradS_gradN,
//     							    *gradS_eval,
//     							     *(fe_->gradN()));
//     // MULTIPLY h*grad(S) . grad(N)
//     Intrepid::FieldContainer<Real> h_gradS_gradN(c, f, p);
//     Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(h_gradS_gradN,
//     								*h,
//     								gradS_gradN);

//     // MULTIPLY h * grad(S) . grad(N) * grad(S)
//     Intrepid::FieldContainer<Real> h_gradS_gradN_gradS(c, f, p, d);
//     for(int i = 0; i < c; i++)
//       {
//     	for(int j = 0; j < p; j++)
//     	  {
//     	    for(int k = 0; k < f; k++)
//     	      {
//     		h_gradS_gradN_gradS(i,k,j,0) = h_gradS_gradN(i,k,j)*(*gradS_eval)(i,j,0);
//     		h_gradS_gradN_gradS(i,k,j,1) = h_gradS_gradN(i,k,j)*(*gradS_eval)(i,j,1);
//     	      }
//     	  }
//       }
        						        
//     // MULTIPLY 2.0 * h * grad(S) . grad(N) * grad(S)
//     Intrepid::RealSpaceTools<Real>::scale(h_gradS_gradN_gradS,2.0);

//     // INTEGRATE 2.0 * h * grad(S) . grad(N) * grad(S) . grad(N)
//     Intrepid::FunctionSpaceTools::integrate<Real>(*J[0][0],
//                                                   h_gradS_gradN_gradS,
//                                                   *(fe_->gradNdetJ()),
//                                                   Intrepid::COMP_CPP, true);

//     // MULTIPLY h*grad(S) . grad(S)
//     Intrepid::FieldContainer<Real> h_gradS_gradS(c, p);
//     Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(h_gradS_gradS,
//                                                                *h,
//                                                                gradS_gradS);

//    // MULTIPLY h * grad(S) . grad(S) * grad(N)
//     Intrepid::FieldContainer<Real> h_gradS_gradS_gradN(c, f, p, d);
//     Intrepid::FunctionSpaceTools::tensorMultiplyDataField<Real>(h_gradS_gradS_gradN,
// 								h_gradS_gradS,
// 								*(fe_->gradN()));
//     // INTEGRATE h * grad(S) . grad(S) * grad(N) . grad(N)
//     Intrepid::FunctionSpaceTools::integrate<Real>(*J[0][0],
//                                                   h_gradS_gradS_gradN,
//                                                   *(fe_->gradNdetJ()),
//                                                   Intrepid::COMP_CPP, true);

//     // Surface velocity

//     ROL::Ptr<Intrepid::FieldContainer<Real> > Sx = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
//     ROL::Ptr<Intrepid::FieldContainer<Real> > Sy = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
//     for(int i = 0; i < c; i++)
//       {
//         for(int j = 0; j < p; j++)
//           {
//             (*Sx)(i,j) = (*gradS_eval)(i,j,0);
//             (*Sy)(i,j) = (*gradS_eval)(i,j,1);
//           }
//       }

//     ROL::Ptr<Intrepid::FieldContainer<Real> > v = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
//     computeSurfaceVelocity(v,valS_eval,valB_eval,A,0);
//     ROL::Ptr<Intrepid::FieldContainer<Real> > d1_v = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
//     computeSurfaceVelocity(d1_v,valS_eval,valB_eval,A,1);
    
//     // J[1][0]

//     // d1_v * (grad(S) . grad(S)) * Sx * N * N                                                                                                                                                                   

//     // MULTIPLY d1_v * (grad(S) . grad(S))                                                                                                                                                                     

//     Intrepid::FieldContainer<Real> dv1_gradS_gradS(c, p);
//     Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(dv1_gradS_gradS,
//                                                                *d1_v,
//                                                                gradS_gradS);

//     // MULTIPLY d1_v * (grad(S) . grad(S)) * Sx                                                                              
//     Intrepid::FieldContainer<Real> dv1_gradS_gradS_Sx(c, p);
//     Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(dv1_gradS_gradS_Sx,
//                                                                dv1_gradS_gradS,
//                                                                *Sx);

//     // MULTIPLY d1_v * (grad(S) . grad(S)) * Sx * N                                                                                                                                 
//     Intrepid::FieldContainer<Real> dv1_gradS_gradS_Sx_N(c, f, p);
//     Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(dv1_gradS_gradS_Sx_N,
//                                                                 dv1_gradS_gradS_Sx,
//                                                                 *(fe_->N()));

//     // INTEGRATE d1_v * (grad(S) . grad(S)) * Sx * N * N                                                                                                                               
//     Intrepid::FunctionSpaceTools::integrate<Real>(*J[1][0],
//                                                   dv1_gradS_gradS_Sx_N,
//                                                   *(fe_->NdetJ()),
//                                                   Intrepid::COMP_CPP, false);

//     ////////////////////////////////////////////////////////////////////////////////////////////////

//     // v * 2.0 * (grad(S) . grad(N)) * Sx * N

//     // MULTIPLY v * Sx                                                                                                                                                                                           
//     Intrepid::FieldContainer<Real> v_Sx(c, p);
//     Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(v_Sx,
//                                                                *v,
//                                                                *Sx);

//     // MULTIPLY 2.0 * v * Sx
//     Intrepid::RealSpaceTools<Real>::scale(v_Sx,2.0);

//     // MULTIPLY 2.0 * v * Sx * grad(S)                                                                                                                                                                      
//     Intrepid::FieldContainer<Real> v_Sx_gradS(c, p, d);
//     Intrepid::FunctionSpaceTools::tensorMultiplyDataData<Real>(v_Sx_gradS,
//                                                                v_Sx,
//                                                                *gradS_eval);

//     // MULTIPLY 2.0 * v * Sx * (grad(S) . grad(N))                                                                                                                                                          
//     Intrepid::FieldContainer<Real> v_Sx_gradS_gradN(c, f, p);
//     Intrepid::FunctionSpaceTools::dotMultiplyDataField<Real>(v_Sx_gradS_gradN,
//     							     v_Sx_gradS,
//     							     *(fe_->gradNdetJ()));

//     // INTEGRATE 2.0 * v * Sx * (grad(S) . grad(N)) * N                                                                                                                                                   
//     Intrepid::FunctionSpaceTools::integrate<Real>(*J[1][0],
//     						  *(fe_->N()),
//                                                   v_Sx_gradS_gradN,
//                                                   Intrepid::COMP_CPP, true);      
                                                                                                                         
//     //////////////////////////////////////////////////////////////////////////////////////////////////

//     // v * (grad(S) . grad(S)) * Nx * N

//     // MULTIPLY v * (grad(S) . grad(S))                                                                                                                                                                   
//     Intrepid::FieldContainer<Real> v_gradS_gradS(c, p);
//     Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(v_gradS_gradS,
//                                                                *v,
//                                                                gradS_gradS);

//     // MULTIPLY v * (grad(S) . grad(S)) * Nx                                                                                                                                                               
//     Intrepid::FieldContainer<Real> v_gradS_gradS_Nx(c, f, p);
//     Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(v_gradS_gradS_Nx,
//     								v_gradS_gradS,
//     								*(fe_->DNDdetJ(0)));

//     // INTEGRATE v * (grad(S) . grad(S)) * Nx * N                                                                                                                                                        
//     Intrepid::FunctionSpaceTools::integrate<Real>(*J[1][0],
//     						  *(fe_->N()),
//                                                   v_gradS_gradS_Nx,
//     						  Intrepid::COMP_CPP, true);   
    
//     ////////////////////////////////////////////////////////////////////////////

//     // J[2][0]

//     // d1_v * (grad(S) . grad(S)) * Sy * N * N                                                                                                                                                                   

//     // MULTIPLY d1_v * (grad(S) . grad(S)) * Sy                                                                              
//     Intrepid::FieldContainer<Real> dv1_gradS_gradS_Sy(c, p);
//     Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(dv1_gradS_gradS_Sy,
//                                                                dv1_gradS_gradS,
//                                                                *Sy);

//     // MULTIPLY d1_v * (grad(S) . grad(S)) * Sy * N                                                                                                                                 
//     Intrepid::FieldContainer<Real> dv1_gradS_gradS_Sy_N(c, f, p);
//     Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(dv1_gradS_gradS_Sy_N,
//                                                                 dv1_gradS_gradS_Sy,
//                                                                 *(fe_->N()));

//     // INTEGRATE d1_v * (grad(S) . grad(S)) * Sy * N * N                                                                                                                               
//     Intrepid::FunctionSpaceTools::integrate<Real>(*J[2][0],
//                                                   dv1_gradS_gradS_Sy_N,
//                                                   *(fe_->NdetJ()),
//                                                   Intrepid::COMP_CPP, false);

//     ////////////////////////////////////////////////////////////////////////////////////////////////

//     // v * 2.0 * (grad(S) . grad(N)) * Sy * N

//     // MULTIPLY v * Sy                                                                                                                                                                                           
//     Intrepid::FieldContainer<Real> v_Sy(c, p);
//     Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(v_Sy,
//                                                                *v,
//                                                                *Sy);

//     // MULTIPLY 2.0 * v * Sy
//     Intrepid::RealSpaceTools<Real>::scale(v_Sy,2.0);

//     // MULTIPLY 2.0 * v * Sy * grad(S)                                                                                                                                                                      
//     Intrepid::FieldContainer<Real> v_Sy_gradS(c, p, d);
//     Intrepid::FunctionSpaceTools::tensorMultiplyDataData<Real>(v_Sy_gradS,
//                                                                v_Sy,
//                                                                *gradS_eval);

//     // MULTIPLY 2.0 * v * Sy * (grad(S) . grad(N))                                                                                                                                                          
//     Intrepid::FieldContainer<Real> v_Sy_gradS_gradN(c, f, p);
//     Intrepid::FunctionSpaceTools::dotMultiplyDataField<Real>(v_Sy_gradS_gradN,
//     							     v_Sy_gradS,
//     							     *(fe_->gradNdetJ()));

//     // INTEGRATE 2.0 * v * Sy * (grad(S) . grad(N)) * N                                                                                                                                                   
//     Intrepid::FunctionSpaceTools::integrate<Real>(*J[2][0],
//     						  *(fe_->N()),
//                                                   v_Sy_gradS_gradN,
//                                                   Intrepid::COMP_CPP, true);   
    
                                                                                                                         
//     //////////////////////////////////////////////////////////////////////////////////////////////////

//     // v * (grad(S) . grad(S)) * Ny * N

//     // MULTIPLY v * (grad(S) . grad(S)) * Ny                                                                                                                                                               
//     Intrepid::FieldContainer<Real> v_gradS_gradS_Ny(c, f, p);
//     Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(v_gradS_gradS_Ny,
//     								v_gradS_gradS,
//     								*(fe_->DNDdetJ(1)));

//     // INTEGRATE v * (grad(S) . grad(S)) * Ny * N                                                                                                                                                        
//     Intrepid::FunctionSpaceTools::integrate<Real>(*J[2][0],
//     						  *(fe_->N()),
//                                                   v_gradS_gradS_Ny,
//     						  Intrepid::COMP_CPP, true);   
    
//     ////////////////////////////////////////////////////////////////////////////


//     // Combine the jacobians.
//     fieldHelper_->combineFieldCoeff(jac, J);
//   }

//   void Jacobian_2(ROL::Ptr<Intrepid::FieldContainer<Real> > & jac,
//                   const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
//                   const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
//                   const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
//     throw Exception::Zero(">>> (PDE_shallow_ice::Jacobian_2): Jacobian is zero.");
//   }
  
//   void Jacobian_3(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > & jac,
//                   const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
//                   const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
//                   const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
//     // GET DIMENSIONS
//     int c = fe_->gradN()->dimension(0);
//     int f = fe_->gradN()->dimension(1);
//     int p = fe_->gradN()->dimension(2);
//     int d = fe_->gradN()->dimension(3);
//     // ADD CONTROL TERM TO RESIDUAL
//     int size = z_param->size();
//     int size_field = (size-1)/2;

//     PDE<Real>::setParameter(*z_param);

//     // Split u_coeff into components.
//     std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > U;
//     fieldHelper_->splitFieldCoeff(U, u_coeff);
//     std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > Z;
//     fieldHelper_->splitFieldCoeff(Z, Z_input);

//     ROL::Ptr<Intrepid::FieldContainer<Real> > val_expBeta_eval
//       = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
//     ROL::Ptr<Intrepid::FieldContainer<Real> > valHflux_eval =
//       ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
//     compute_parameters(val_expBeta_eval,valHflux_eval);
//     for (int i = 0; i < c; ++i) {
//       for (int j = 0; j < p; ++j) {
//         (*val_expBeta_eval)(i,j) = std::exp((*val_expBeta_eval)(i,j));
//       }
//     }

//     ROL::Ptr<Intrepid::FieldContainer<Real> > valBeta_eval_nom
//       = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
//     ROL::Ptr<Intrepid::FieldContainer<Real> > valHflux_eval_nom =
//       ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
//     compute_parameters(valBeta_eval_nom,valHflux_eval_nom,true);

//     Intrepid::RealSpaceTools<Real>::scale(*valHflux_eval_nom,-1.0);

//     ROL::Ptr<Intrepid::FieldContainer<Real> > valS_eval =
//       ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
//     fe_->evaluateValue(valS_eval, U[0]);

//     ROL::Ptr<Intrepid::FieldContainer<Real> > valB_eval =
//       ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
//     fe_->evaluateValue(valB_eval, Z[0]);

//     ROL::Ptr<Intrepid::FieldContainer<Real> > gradS_eval
//       = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p, d);
//     fe_->evaluateGradient(gradS_eval, U[0]);

//     // COMPUTE BASAL VELOCITY
//     ROL::Ptr<Intrepid::FieldContainer<Real> > kappa
//       = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
//     computeBasalVelocity(kappa,val_expBeta_eval,valS_eval,valB_eval,0);

//     ROL::Ptr<Intrepid::FieldContainer<Real> > beta_kappa_eval =
//       ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
//     Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*beta_kappa_eval,
// 							       *valBeta_eval_nom,
// 							       *kappa);

//     // Basal sliding parameters
//     for (int i = 0; i < size_field; ++i) 
//       {
// 	std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > J(3);
// 	J[0] = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, f);
// 	J[1] = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, f);
// 	J[2] = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, f);

// 	ROL::Ptr<Intrepid::FieldContainer<Real> > phi_i_eval =
// 	  ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
// 	compute_phi_k(phi_i_eval, i);
	
// 	ROL::Ptr<Intrepid::FieldContainer<Real> > beta_kappa_phi_i_eval =
// 	  ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
// 	Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*beta_kappa_phi_i_eval,
// 								   *beta_kappa_eval,
// 								   *phi_i_eval);

// 	ROL::Ptr<Intrepid::FieldContainer<Real> > beta_kappa_phi_i_eval_gradS =
// 	  ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p, d);
// 	Intrepid::FunctionSpaceTools::tensorMultiplyDataData<Real>(*beta_kappa_phi_i_eval_gradS,
// 								   *beta_kappa_phi_i_eval,
// 								   *gradS_eval);
// 	Intrepid::FunctionSpaceTools::integrate<Real>(*J[0],
// 						      *beta_kappa_phi_i_eval_gradS,
// 						      *(fe_->gradNdetJ()),
// 						      Intrepid::COMP_CPP, false);
	
// 	fieldHelper_->combineFieldCoeff(jac[i], J);
//       }

//     // Forcing term parameters
//     for (int i = 0; i < size_field; ++i) 
//       {
// 	std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > J(3);
// 	J[0] = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, f);
// 	J[1] = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, f);
// 	J[2] = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, f);

// 	ROL::Ptr<Intrepid::FieldContainer<Real> > phi_i_eval =
// 	  ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
// 	compute_phi_k(phi_i_eval, i);
	
// 	ROL::Ptr<Intrepid::FieldContainer<Real> > hflux_phi_i_eval =
// 	  ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
// 	Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*hflux_phi_i_eval,
// 								   *valHflux_eval_nom,
// 								   *phi_i_eval);
// 	Intrepid::FunctionSpaceTools::integrate<Real>(*J[0],
// 						      *hflux_phi_i_eval,
// 						      *(fe_->NdetJ()),
// 						      Intrepid::COMP_CPP, false);
	
// 	fieldHelper_->combineFieldCoeff(jac[size_field+i], J);
//       }

//     // Flow rate factor
//     std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > J_A(3);
//     J_A[0] = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, f);
//     J_A[1] = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, f);
//     J_A[2] = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, f);
//     ROL::Ptr<Intrepid::FieldContainer<Real> > h
//       = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
//     computeThermalVelocity(h,valS_eval,valB_eval,A_*a_,0);
//     // MULTIPLY grad(S) . grad(S)
//     Intrepid::FieldContainer<Real> gradS_gradS(c, p);
//     Intrepid::FunctionSpaceTools::dotMultiplyDataData<Real>(gradS_gradS,
// 							    *gradS_eval,
// 							    *gradS_eval);
//     // MULTIPLY h*grad(S) . grad(S)
//     Intrepid::FieldContainer<Real> h_gradS_gradS(c, p);
//     Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(h_gradS_gradS,
//                                                                *h,
//                                                                gradS_gradS);
//     ROL::Ptr<Intrepid::FieldContainer<Real> > h_gradS_gradS_gradS =
//       ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p, d);
//     Intrepid::FunctionSpaceTools::tensorMultiplyDataData<Real>(*h_gradS_gradS_gradS,
// 							       h_gradS_gradS,
// 							       *gradS_eval);
//     Intrepid::FunctionSpaceTools::integrate<Real>(*J_A[0],
// 						  *h_gradS_gradS_gradS,
// 						  *(fe_->gradNdetJ()),
// 						  Intrepid::COMP_CPP, false);
//     // Surface velocity terms

//     ROL::Ptr<Intrepid::FieldContainer<Real> > v = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
//     computeSurfaceVelocity(v,valS_eval,valB_eval,A_*a_,0);
//     ROL::Ptr<Intrepid::FieldContainer<Real> > Sx = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
//     ROL::Ptr<Intrepid::FieldContainer<Real> > Sy = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
//     for(int i = 0; i < c; i++)
//       {
//         for(int j = 0; j < p; j++)
//           {
//             (*Sx)(i,j) = (*gradS_eval)(i,j,0);
//             (*Sy)(i,j) = (*gradS_eval)(i,j,1);
//           }
//       }

//     // MULTIPLY v * (grad(S) . grad(S))
//     Intrepid::FieldContainer<Real> v_gradS_gradS(c, p);
//     Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(v_gradS_gradS,
//                                                                *v,
//                                                                gradS_gradS);    
//    // MULTIPLY v * (grad(S) . grad(S)) * Sx
//     Intrepid::FieldContainer<Real> v_gradS_gradS_Sx(c, p);
//     Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(v_gradS_gradS_Sx,
//                                                                *Sx,
//                                                                v_gradS_gradS); 

//     Intrepid::FunctionSpaceTools::integrate<Real>(*J_A[1], v_gradS_gradS_Sx, *(fe_->NdetJ()), Intrepid::COMP_CPP, false);

//     // MULTIPLY v * (grad(S) . grad(S)) * Sy
//     Intrepid::FieldContainer<Real> v_gradS_gradS_Sy(c, p);
//     Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(v_gradS_gradS_Sy,
//                                                                *Sy,
//                                                                v_gradS_gradS); 

//     Intrepid::FunctionSpaceTools::integrate<Real>(*J_A[2], v_gradS_gradS_Sy, *(fe_->NdetJ()), Intrepid::COMP_CPP, false);

//     fieldHelper_->combineFieldCoeff(jac[2*size_field], J_A);
//   }

//   void Hessian_11(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
//                   const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
//                   const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
//                   const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
//                   const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
//     throw Exception::Zero(">>> (PDE_shallow_ice::Hessian_11): Hessian is zero.");
//   }

//   void Hessian_12(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
//                   const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
//                   const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
//                   const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
//                   const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
//     throw Exception::Zero(">>> (PDE_shallow_ice::Hessian_12): Hessian is zero.");
//   }

//   void Hessian_13(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > & hess,
//                   const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
//                   const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
//                   const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
//                   const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
//     throw Exception::Zero(">>> (PDE_shallow_ice::Hessian_13): Hessian is zero.");
//   }

//   void Hessian_21(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
//                   const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
//                   const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
//                   const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
//                   const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
//     throw Exception::Zero(">>> (PDE_shallow_ice::Hessian_21): Hessian is zero.");
//   }

//   void Hessian_22(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
//                   const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
//                   const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
//                   const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
//                   const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
//     throw Exception::Zero(">>> (PDE_shallow_ice::Hessian_22): Hessian is zero.");
//   }

//   void Hessian_23(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > & hess,
//                   const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
//                   const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
//                   const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
//                   const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
//     throw Exception::Zero(">>> (PDE_shallow_ice::Hessian_23): Hessian is zero.");
//   }

//   void Hessian_31(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > & hess,
//                   const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
//                   const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
//                   const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
//                   const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
//     throw Exception::Zero(">>> (PDE_shallow_ice::Hessian_31): Hessian is zero.");
//   }

//   void Hessian_32(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > & hess,
//                   const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
//                   const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
//                   const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
//                   const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
//     throw Exception::Zero(">>> (PDE_shallow_ice::Hessian_32): Hessian is zero.");
//   }

//   void Hessian_33(std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > & hess,
//                   const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
//                   const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
//                   const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
//                   const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
//     throw Exception::Zero(">>> (PDE_shallow_ice::Hessian_33): Hessian is zero.");
//   }

//    void RieszMap_1(ROL::Ptr<Intrepid::FieldContainer<Real> > & riesz) {
//      throw Exception::NotImplemented(">>> (PDE_shallow_ice::RieszMap_1): Not implemented.");
//   }

//   void RieszMap_2(ROL::Ptr<Intrepid::FieldContainer<Real> > & riesz) {
//     throw Exception::NotImplemented(">>> (PDE_shallow_ice::RieszMap_2): Not implemented.");
//   }

//   void setFieldPattern(const std::vector<std::vector<int> > & fieldPattern) {
//     fieldPattern_ = fieldPattern;
//     fieldHelper_ = ROL::makePtr<FieldHelper<Real>>(numFields_, numDofs_, numFieldDofs_, fieldPattern_);
//   }

//  void setCellNodes(const ROL::Ptr<Intrepid::FieldContainer<Real> > &volCellNodes,
//                     const std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > &bdryCellNodes,
//                     const std::vector<std::vector<std::vector<int> > > &bdryCellLocIds) {
//     volCellNodes_ = volCellNodes;
//     bdryCellNodes_ = bdryCellNodes;
//     bdryCellLocIds_ = bdryCellLocIds;
//     // Finite element definition.
//     fe_ = ROL::makePtr<FE<Real>>(volCellNodes_,basisPtr_,cellCub_);
//     // Set local boundary DOFs.
//     fidx_ = fe_->getBoundaryDofs();
//     // Compute Dirichlet values at DOFs.
//     int d = basisPtr_->getBaseCellTopology().getDimension();
//     int numSidesets = bdryCellLocIds_.size();
//     bdryCellDofValues_.resize(numSidesets);
//     feBdry_.resize(numSidesets); 
//     for (int i=0; i<numSidesets; ++i) {
//       int numLocSides = bdryCellLocIds_[i].size();
//       bdryCellDofValues_[i].resize(numLocSides);
//       feBdry_[i].resize(numLocSides);
//       for (int j=0; j<numLocSides; ++j) {
//         int c = bdryCellLocIds_[i][j].size();
//         int f = basisPtr_->getCardinality();
// 	if(c > 0)
// 	  {
// 	    feBdry_[i][j] = ROL::makePtr<FE<Real> >(bdryCellNodes_[i][j],basisPtr_,bdryCub_,j);
// 	  }
// 	bdryCellDofValues_[i][j] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
//         ROL::Ptr<Intrepid::FieldContainer<Real> > coords =
//           ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, d);
//         if (c > 0) {
//           fe_->computeDofCoords(coords, bdryCellNodes_[i][j]);
//         }
//         for (int k=0; k<c; ++k) {
//           for (int l=0; l<f; ++l) {
//             std::vector<Real> dofpoint(d);
//             for (int m=0; m<d; ++m) {
//               dofpoint[m] = (*coords)(k, l, m);
//             }
//             (*bdryCellDofValues_[i][j])(k, l) = 0.0;
//           }
//         }
//       }
//     }

//   }

//   ROL::Ptr<Intrepid::FieldContainer<Real> > getBoundaryCoeff(const Intrepid::FieldContainer<Real> & cell_coeff, int sideSet, int cell) const {
//     std::vector<int> bdryCellLocId = bdryCellLocIds_[sideSet][cell];
//     const int numCellsSide = bdryCellLocId.size();
//     const int f = basisPtr_->getCardinality();
    
//     ROL::Ptr<Intrepid::FieldContainer<Real > > bdry_coeff = 
//       ROL::makePtr<Intrepid::FieldContainer<Real > >(numCellsSide, f);
//     for (int i = 0; i < numCellsSide; ++i) {
//       for (int j = 0; j < f; ++j) {
// 	(*bdry_coeff)(i, j) = cell_coeff(bdryCellLocId[i], j);
//       }
//     }
//     return bdry_coeff;
//   }

//   const ROL::Ptr<FE<Real> > getFE(void) const {
//     return fe_;
//   }

//   const ROL::Ptr<FieldHelper<Real> > getFieldHelper(void) const {
//     return fieldHelper_;
//   }

//   std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>>> getFields() {
//     return basisPtrs_;
//   }

//   const std::vector<std::vector<std::vector<int> > > getBdryCellLocIds(void) const {
//     return bdryCellLocIds_;
//   }

// private:
  
//   void computeBasalVelocity(ROL::Ptr<Intrepid::FieldContainer<Real> > &kappa,
// 			    const ROL::Ptr<Intrepid::FieldContainer<Real> > &beta,
// 			    const ROL::Ptr<Intrepid::FieldContainer<Real> > &u,
// 			    const ROL::Ptr<Intrepid::FieldContainer<Real> > &z,
// 			    const int deriv = 0 ) const {
//     // GET DIMENSIONS
//     int c = fe_->gradN()->dimension(0);
//     int p = fe_->gradN()->dimension(2);
//     int d = fe_->gradN()->dimension(3);
//     std::vector<Real> pt(d);
//     for (int i = 0; i < c; ++i) {
//       for (int j = 0; j < p; ++j) {
// 	pt[0] = (*fe_->cubPts())(i,j,0)/width_;
// 	pt[1] = (*fe_->cubPts())(i,j,1)/height_;
//         // Compute basal velocity
//         (*kappa)(i,j) = evaluateBasalVelocity((*beta)(i,j),(*u)(i,j),(*z)(i,j),pt,deriv);
//       }
//     }
//   }

//   Real evaluateBasalVelocity(const Real beta, const Real u, const Real z, const std::vector<Real> & x, const int deriv = 0) const {
//     Real val = 0.0;
//     if(deriv == 0)
//       {
// 	val = beta*rho_*g_*std::pow(u-z,2.0);
//       }
//     if(deriv == 1)
//       {
// 	val = beta*rho_*g_*2.0*(u-z);
//       }
//     if(deriv == 2)
//       {
// 	val = beta*rho_*g_*2.0;
//       }
//     return val;
//   }

//   void computeThermalVelocity(ROL::Ptr<Intrepid::FieldContainer<Real> > &h,
// 			      const ROL::Ptr<Intrepid::FieldContainer<Real> > &u,
// 			      const ROL::Ptr<Intrepid::FieldContainer<Real> > &z,
// 			      const Real A,
// 			      const int deriv = 0 ) const {
//     // GET DIMENSIONS
//     int c = fe_->gradN()->dimension(0);
//     int p = fe_->gradN()->dimension(2);
//     int d = fe_->gradN()->dimension(3);
//     std::vector<Real> pt(d);
//     for (int i = 0; i < c; ++i) {
//       for (int j = 0; j < p; ++j) {
// 	pt[0] = (*fe_->cubPts())(i,j,0)/width_;
// 	pt[1] = (*fe_->cubPts())(i,j,1)/height_;
//         // Compute thermal velocity
//         (*h)(i,j) = evaluateThermalVelocity((*u)(i,j),(*z)(i,j),pt,A,deriv);
//       }
//     }
//   }

//   Real evaluateThermalVelocity(const Real u, const Real z, const std::vector<Real> & x, const Real A, const int deriv = 0) const {
//     Real val = 0.0;
//     if(deriv == 0)
//       {
// 	val = 0.4*A*std::pow(rho_*g_,3.0)*std::pow(u-z,5.0);
//       }
//     if(deriv == 1)
//       {
// 	val = 0.4*A*std::pow(rho_*g_,3.0)*5.0*std::pow(u-z,4.0);
//       }
//     if(deriv == 2)
//       {
//         val = 0.4*A*std::pow(rho_*g_,3.0)*20.0*std::pow(u-z,3.0);
//       }
//     return val;
//   }

//   void computeSurfaceVelocity(ROL::Ptr<Intrepid::FieldContainer<Real> > &v,
// 			      const ROL::Ptr<Intrepid::FieldContainer<Real> > &u,
// 			      const ROL::Ptr<Intrepid::FieldContainer<Real> > &z,
// 			      const Real A,
// 			      const int deriv = 0 ) const {
//     // GET DIMENSIONS
//     int c = fe_->gradN()->dimension(0);
//     int p = fe_->gradN()->dimension(2);
//     int d = fe_->gradN()->dimension(3);
//     std::vector<Real> pt(d);
//     for (int i = 0; i < c; ++i) {
//       for (int j = 0; j < p; ++j) {
//  	pt[0] = (*fe_->cubPts())(i,j,0)/width_;
// 	pt[1] = (*fe_->cubPts())(i,j,1)/height_;
//         // Compute thermal velocity
//         (*v)(i,j) = evaluateSurfaceVelocity((*u)(i,j),(*z)(i,j),pt,A,deriv);
//       }
//     }
//   }

//   Real evaluateSurfaceVelocity(const Real u, const Real z, const std::vector<Real> & x, const Real A, const int deriv = 0) const {
//     Real val = 0.0;
//     if(deriv == 0)
//       {
// 	val = 0.5*A*std::pow(rho_*g_,3.0)*std::pow(u-z,4.0);
//       }
//     if(deriv == 1)
//       {
// 	val = 0.5*A*std::pow(rho_*g_,3.0)*4.0*std::pow(u-z,3.0);
//       }
//     if(deriv == 2)
//       {
//         val = 0.5*A*std::pow(rho_*g_,3.0)*12.0*std::pow(u-z,2.0);
//       }
//     return 1000.0*val;
//   }

//   void compute_parameters(ROL::Ptr<Intrepid::FieldContainer<Real> > & beta, ROL::Ptr<Intrepid::FieldContainer<Real> > & h_flux,  bool nominal_param = false) const {
//     int c = fe_->gradN()->dimension(0);
//     int p = fe_->gradN()->dimension(2);
//     int d = fe_->gradN()->dimension(3);
//     std::vector<Real> pt(d);
//     for (int i = 0; i < c; ++i) {
//       for (int j = 0; j < p; ++j) {
// 	pt[0] = (*fe_->cubPts())(i,j,0)/width_;
// 	pt[1] = (*fe_->cubPts())(i,j,1)/height_;
//         (*beta)(i,j) = basal_sliding(pt,nominal_param);
//         (*h_flux)(i,j) = rhs(pt,nominal_param);
//       }
//     }

//   }

//   Real basal_sliding(const std::vector<Real> & pt, const bool nominal_param = false) const {
//     Real val = 0.0;
//     Real x = pt[0];
//     Real y = pt[1];

//     int i = std::floor(x/h_);
//     int j = std::floor(y/h_);
//     for(int ii = -1; ii < 2; ii++)
//       {
// 	for(int jj = -1; jj < 2; jj++)
// 	  {
// 	    int i_index = i+ii;
// 	    int j_index = j+jj;
// 	    if(i_index>=0 && j_index>=0 && i_index<N_+1 && j_index<N_+1)
// 	      {
// 		int k = j_index*(N_+1)+i_index;
// 		val += basal_coeff_[k]*Linear_FE_Basis_Fun_Eval(x,i_index)*Linear_FE_Basis_Fun_Eval(y,j_index);
// 	      }
// 	  }
//       }
//     if(!nominal_param)
//       {
// 	const std::vector<Real> param = PDE<Real>::getParameter();
// 	if((int)param.size()>0)
// 	  {
// 	    Real val_uncertain = 0.0;
// 	    int count = 0;
// 	    Real basis_fun_x = 0.0;
// 	    Real basis_fun_y = 0.0;
// 	    for(int j = 0; j < L_+1; j++)
// 	      {
// 	    	basis_fun_y = Parameter_Basis_Fun_Eval(y,j);
// 	    	for(int i = 0; i < L_+1; i++)
// 	    	  {
// 	    	    basis_fun_x = Parameter_Basis_Fun_Eval(x,i);
// 	    	    val_uncertain += param[count]*basis_fun_x*basis_fun_y;
// 	    	    count = count+1;
// 	    	  }
// 	      }
// 	    val_uncertain = a_*val_uncertain;
// 	    val_uncertain = 1+val_uncertain;
// 	    val = val*val_uncertain;
// 	  }
//       }

//     return val;
//   }

//   Real rhs(const std::vector<Real> & pt, bool & nominal_param = false) const {
//     Real val = 0.0;
//     Real x = pt[0];
//     Real y = pt[1];
//     Real t = PDE<Real>::getTime();
//     Real pi(M_PI);
//     Real time_coeff = 0.0;

//     int i = std::floor(x/h_);
//     int j = std::floor(y/h_);
//     for(int ii = -1; ii < 2; ii++)
//       {
// 	for(int jj = -1; jj < 2; jj++)
// 	  {
// 	    int i_index = i+ii;
// 	    int j_index = j+jj;
// 	    if(i_index>=0 && j_index>=0 && i_index<N_+1 && j_index<N_+1)
// 	      {
// 		int k = j_index*(N_+1)+i_index;
// 		val += time_coeff*forcing_coeff_[k]*Linear_FE_Basis_Fun_Eval(x,i_index)*Linear_FE_Basis_Fun_Eval(y,j_index);
// 	      }
// 	  }
//       }
//     if(!nominal_param)
//       {
// 	const std::vector<Real> param = PDE<Real>::getParameter();
// 	if((int)param.size()>0)
// 	  {
// 	    Real val_uncertain = 0.0;
// 	    int count = (L_+1)*(L_+1);
// 	    Real basis_fun_x = 0.0;
// 	    Real basis_fun_y = 0.0;
// 	    for(int j = 0; j < L_+1; j++)
// 	      {
// 		basis_fun_y = Parameter_Basis_Fun_Eval(y,j);
// 		for(int i = 0; i < L_+1; i++)
// 		  {
// 		    basis_fun_x = Parameter_Basis_Fun_Eval(x,i);
// 		    val_uncertain += param[count]*basis_fun_x*basis_fun_y;
// 		    count = count+1;
// 		  }
// 	      }
// 	    val_uncertain = a_*val_uncertain;
// 	    val = val*val_uncertain;
// 	  }
//       }

//     return val;
//   }

//   void compute_phi_k(ROL::Ptr<Intrepid::FieldContainer<Real> > & phi, int k) const {
//     int k_x = k%(L_+1);
//     int k_y = std::floor(k/(L_+1));
//     int c = fe_->gradN()->dimension(0);
//     int p = fe_->gradN()->dimension(2);
//     int d = fe_->gradN()->dimension(3);
//     std::vector<Real> pt(d);
//     for (int i = 0; i < c; ++i) {
//       for (int j = 0; j < p; ++j) {
// 	pt[0] = (*fe_->cubPts())(i,j,0)/width_;
// 	pt[1] = (*fe_->cubPts())(i,j,1)/height_;
//         (*phi)(i,j) = a_*Parameter_Basis_Fun_Eval(pt[0],k_x)*Parameter_Basis_Fun_Eval(pt[1],k_y);
//       }
//     }

//   }

//   Real Linear_FE_Basis_Fun_Eval(const Real & x, int i) const
//   {
//     // Evaluates ith 1D FE basis function at point x
//     Real val = 0.0;
//     Real xi = static_cast<Real>(i)*h_;
//     Real dist = std::abs(x-xi);
//     if( dist < h_ )
//       {
// 	val = 1.0 - dist/h_;
//       }
//     return val;
//   }

//   Real Parameter_Basis_Fun_Eval(const Real & x, int i) const
//   {
//     // Evaluates ith 1D basis function on a grid of L+1 nodes at point x
//     Real val = 0.0;
//     Real dist = std::abs(x-uncertain_basis_grid_[i]);
//     if( dist < 1.0/static_cast<Real>(L_) )
//       {
// 	val = 1.0 - dist*static_cast<Real>(L_);
//       }
//     return val;
//   }

// public:
//   void Update_Z_input(ROL::Ptr<Intrepid::FieldContainer<Real> > & Z)
//   {
//     Z_input = Z;
//   }

// }; // PDE_shallow_ice

// #endif
