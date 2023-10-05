#ifndef PDE_SHALLOW_ICE_HPP
#define PDE_SHALLOW_ICE_HPP

#include "../../../PDE-OPT/TOOLS/pde.hpp"
#include "../../../PDE-OPT/TOOLS/fe.hpp"

#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_CellTools.hpp"

#include "ROL_Ptr.hpp"

template <class Real>
class PDE_shallow_ice : public PDE<Real> {
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
  std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > bdryCellDofValues_xvel;
  std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > bdryCellDofValues_yvel;

  // Field pattern, offsets, etc.
  std::vector<std::vector<int> > fieldPattern_;  // local Field/DOF pattern; set from DOF manager 
  int numFields_;                                // number of fields (equations in the PDE)
  int numDofs_;                                  // total number of degrees of freedom for all (local) fields
  std::vector<int> offset_;                      // for each field, a counting offset
  std::vector<int> numFieldDofs_;                // for each field, number of degrees of freedom
  ROL::Ptr<FieldHelper<Real> > fieldHelper_;

  Real width_;
  Real height_;

  Real B_;
  Real rho_;
  Real g_;
  Real nob1_;
  Real nob2_;
  Real nob3_;

  int num_coeff_load_;
  std::vector<Real> bedrock_topo_coeff_;
  std::vector<Real> surface_height_coeff_;
  std::vector<Real> xvel_coeff_;
  std::vector<Real> yvel_coeff_;
  int N_;
  Real h_;

public:
  PDE_shallow_ice(Teuchos::ParameterList &parlist) {
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
    basisPtrs_.clear(); basisPtrs_.push_back(basisPtr_); basisPtrs_.push_back(basisPtr_);

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

    width_ = parlist.sublist("Geometry").get("Width",1.0);
    height_ = parlist.sublist("Geometry").get("Height",1.0);

    rho_ = 910.0;
    g_ = 9.81;
    B_ = parlist.sublist("Problem").get("B",1.0);
    nob1_ = parlist.sublist("Problem").get("nob1",1.0);
    nob2_ = parlist.sublist("Problem").get("nob2",1.0);
    nob3_ = parlist.sublist("Problem").get("nob3",1.0);

    num_coeff_load_ = parlist.sublist("Problem").get("Number of Coefficients in Loaded Fields", 10); 
    N_ = std::sqrt(num_coeff_load_)-1;
    h_ = 1.0/static_cast<Real>(N_);

    bedrock_topo_coeff_.resize(num_coeff_load_);
    // read in data
    std::ifstream in_bedrock("Bedrock_Topography.txt");          
    // read the elements in the file into a vector  
    // test file open   
    if (in_bedrock) 
      {   
	for(int j = 0; j < num_coeff_load_; j++)
	  {
	    in_bedrock >> bedrock_topo_coeff_[j];
	  }
      }
    else
      {
	std::cout << "Error loading the data from Bedrock_Topography.txt" << std::endl;
      }   

    surface_height_coeff_.resize(num_coeff_load_);
    // read in data
    std::ifstream in_surface("Surface_Height.txt");          
    // read the elements in the file into a vector  
    // test file open   
    if (in_surface) 
      {   
	for(int j = 0; j < num_coeff_load_; j++)
	  {
	    in_surface >> surface_height_coeff_[j];
	  }
      }
    else
      {
	std::cout << "Error loading the data from Surface_Height.txt" << std::endl;
      }  

    xvel_coeff_.resize(num_coeff_load_);
    // read in data
    std::ifstream in_xvel("Horizontal_Velocity.txt");          
    // read the elements in the file into a vector  
    // test file open   
    if (in_xvel) 
      {   
	for(int j = 0; j < num_coeff_load_; j++)
	  {
	    in_xvel >> xvel_coeff_[j];
	  }
      }
    else
      {
	std::cout << "Error loading the data from Horzontal_Velocity.txt" << std::endl;
      }  

    yvel_coeff_.resize(num_coeff_load_);
    // read in data
    std::ifstream in_yvel("Vertical_Velocity.txt");          
    // read the elements in the file into a vector  
    // test file open   
    if (in_yvel) 
      {   
	for(int j = 0; j < num_coeff_load_; j++)
	  {
	    in_yvel >> yvel_coeff_[j];
	  }
      }
    else
      {
	std::cout << "Error loading the data from Vertical_Velocity.txt" << std::endl;
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
     std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > R(2);
     R[0]   = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, f);
     R[1] = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, f);

     // Split u_coeff into components.
     std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > U;
     fieldHelper_->splitFieldCoeff(U, u_coeff);
     std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > Z;
     fieldHelper_->splitFieldCoeff(Z, z_coeff);
     
     ROL::Ptr<Intrepid::FieldContainer<Real> > valB_eval
       = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
     ROL::Ptr<Intrepid::FieldContainer<Real> > valS_eval
       = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
     ROL::Ptr<Intrepid::FieldContainer<Real> > valSx_eval
       = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
     ROL::Ptr<Intrepid::FieldContainer<Real> > valSy_eval
       = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
     compute_parameters(valB_eval,valS_eval,valSx_eval,valSy_eval);

     ROL::Ptr<Intrepid::FieldContainer<Real> > valH_eval
       = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
     Intrepid::RealSpaceTools<Real>::add(*valH_eval, *valS_eval);
     Intrepid::RealSpaceTools<Real>::subtract(*valH_eval, *valB_eval);

     ROL::Ptr<Intrepid::FieldContainer<Real> > valBeta_eval =
       ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
     fe_->evaluateValue(valBeta_eval, Z[0]);
     for (int i = 0; i < c; ++i) {
       for (int j = 0; j < p; ++j) {
	 (*valBeta_eval)(i,j) = nob1_*std::exp((*valBeta_eval)(i,j));
       }
     }
     
     ROL::Ptr<Intrepid::FieldContainer<Real> > val_xvel_eval =
       ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
     fe_->evaluateValue(val_xvel_eval, U[0]);
     ROL::Ptr<Intrepid::FieldContainer<Real> > grad_xvel_eval
       = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p, d);
     fe_->evaluateGradient(grad_xvel_eval, U[0]);
     ROL::Ptr<Intrepid::FieldContainer<Real> > val_yvel_eval =
       ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
     fe_->evaluateValue(val_yvel_eval, U[1]);
     ROL::Ptr<Intrepid::FieldContainer<Real> > grad_yvel_eval
       = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p, d);
     fe_->evaluateGradient(grad_yvel_eval, U[1]);
     
    // MULTIPLY beta * xvel
    Intrepid::FieldContainer<Real> beta_xvel(c, p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(beta_xvel,
                                                               *val_xvel_eval,
                                                               *valBeta_eval);
    Intrepid::FunctionSpaceTools::integrate<Real>(*R[0], beta_xvel, *(fe_->NdetJ()), Intrepid::COMP_CPP, false);

    // MULTIPLY rho*g*H*Sx
    Intrepid::FieldContainer<Real> H_Sx(c, p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(H_Sx,
                                                               *valH_eval,
                                                               *valSx_eval);
    Intrepid::RealSpaceTools<Real>::scale(H_Sx, rho_*g_);
    Intrepid::FunctionSpaceTools::integrate<Real>(*R[0], H_Sx, *(fe_->NdetJ()), Intrepid::COMP_CPP, true);

    // MULTIPLY beta * yvel
    Intrepid::FieldContainer<Real> beta_yvel(c, p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(beta_yvel,
                                                               *val_yvel_eval,
                                                               *valBeta_eval);
    Intrepid::FunctionSpaceTools::integrate<Real>(*R[1], beta_yvel, *(fe_->NdetJ()), Intrepid::COMP_CPP, false);

    // MULTIPLY rho*g*H*Sy
    Intrepid::FieldContainer<Real> H_Sy(c, p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(H_Sy,
                                                               *valH_eval,
                                                               *valSy_eval);
    Intrepid::RealSpaceTools<Real>::scale(H_Sy, rho_*g_);
    Intrepid::FunctionSpaceTools::integrate<Real>(*R[1], H_Sy, *(fe_->NdetJ()), Intrepid::COMP_CPP, true);

    // Nonlinear terms
    Intrepid::FieldContainer<Real> epsilon(c, p);
    for(int i = 0; i < c; i++)
      {
	for(int j = 0; j < p; j++)
	  {
	    //epsilon(i,j) = std::pow((*grad_xvel_eval)(i,j,0),2.0) + std::pow((*grad_yvel_eval)(i,j,1),2.0) + 0.25*( (*grad_xvel_eval)(i,j,1) + (*grad_yvel_eval)(i,j,0) ) + (*grad_xvel_eval)(i,j,0)*(*grad_yvel_eval)(i,j,1);
	    epsilon(i,j) = 1.0 + nob2_*std::pow((*grad_xvel_eval)(i,j,0),2.0) + nob2_*std::pow((*grad_yvel_eval)(i,j,1),2.0) + nob3_*0.25*( (*grad_xvel_eval)(i,j,1) );
	  }
      }

    Intrepid::FieldContainer<Real> nu(c, p);
    for(int i = 0; i < c; i++)
      {
	for(int j = 0; j < p; j++)
	  {
	    nu(i,j) = 0.5*B_*(1.0/std::cbrt(epsilon(i,j)));
	    //nu(i,j) = 0.5*B_*std::pow(epsilon(i,j),2.0);
	  }
      }

    Intrepid::FieldContainer<Real> zx(c, p, 2);
    Intrepid::FieldContainer<Real> zy(c, p, 2);
    for(int i = 0; i < c; i++)
      {
	for(int j = 0; j < p; j++)
	  {
	    zx(i,j,0) = 4.0*(*grad_xvel_eval)(i,j,0) + 2.0*(*grad_yvel_eval)(i,j,1);
	    zx(i,j,1) = (*grad_xvel_eval)(i,j,1) + (*grad_yvel_eval)(i,j,0);
	    zy(i,j,0) = 4.0*(*grad_yvel_eval)(i,j,1) + 2.0*(*grad_xvel_eval)(i,j,0);
	    zy(i,j,1) = (*grad_xvel_eval)(i,j,1) + (*grad_yvel_eval)(i,j,0);
	    // zx(i,j,0) = 1.0;
	    // zx(i,j,1) = 1.0;
	    // zy(i,j,0) = 1.0;
	    // zy(i,j,1) = 1.0;
	  }
      }

    Intrepid::FieldContainer<Real> wx(c, p, 2);
    Intrepid::FieldContainer<Real> wy(c, p, 2);
    for(int i = 0; i < c; i++)
      {
	for(int j = 0; j < p; j++)
	  {
	    wx(i,j,0) = nu(i,j)*(*valH_eval)(i,j)*zx(i,j,0);
	    wx(i,j,1) = nu(i,j)*(*valH_eval)(i,j)*zx(i,j,1);
	    wy(i,j,0) = nu(i,j)*(*valH_eval)(i,j)*zy(i,j,0);
	    wy(i,j,1) = nu(i,j)*(*valH_eval)(i,j)*zy(i,j,1);
	  }
      }

    Intrepid::FunctionSpaceTools::integrate<Real>(*R[0], wx, *(fe_->gradNdetJ()), Intrepid::COMP_CPP, true);
    Intrepid::FunctionSpaceTools::integrate<Real>(*R[1], wy, *(fe_->gradNdetJ()), Intrepid::COMP_CPP, true);

    // APPLY DIRICHLET CONDITIONS
    int numLocalSideIds = bdryCellLocIds_[0].size();
    for (int j = 0; j < numLocalSideIds; ++j) {
      int numCellsSide = bdryCellLocIds_[0][j].size();
      int numBdryDofs = fidx_[j].size();
      for (int k = 0; k < numCellsSide; ++k) {
	int cidx = bdryCellLocIds_[0][j][k];
	for (int l = 0; l < numBdryDofs; ++l) {
	  (*R[0])(cidx,fidx_[j][l]) = (*U[0])(cidx,fidx_[j][l]) - (*bdryCellDofValues_xvel[0][j])(k,fidx_[j][l]);
	  (*R[1])(cidx,fidx_[j][l]) = (*U[1])(cidx,fidx_[j][l]) - (*bdryCellDofValues_yvel[0][j])(k,fidx_[j][l]);
	}
      }
    }

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
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > J(2);
    for (int i = 0; i < 2; ++i) {
      J[i].resize(2);
      for (int j = 0; j < 2; j++)
	{
	  J[i][j] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
	}
    }
    
    // Split u_coeff into components.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > U;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > Z;
    fieldHelper_->splitFieldCoeff(Z, z_coeff);
    
    ROL::Ptr<Intrepid::FieldContainer<Real> > valB_eval
      = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    ROL::Ptr<Intrepid::FieldContainer<Real> > valS_eval
      = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    ROL::Ptr<Intrepid::FieldContainer<Real> > valSx_eval
      = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    ROL::Ptr<Intrepid::FieldContainer<Real> > valSy_eval
      = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    compute_parameters(valB_eval,valS_eval,valSx_eval,valSy_eval);

    ROL::Ptr<Intrepid::FieldContainer<Real> > valH_eval
      = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    Intrepid::RealSpaceTools<Real>::add(*valH_eval, *valS_eval);
    Intrepid::RealSpaceTools<Real>::subtract(*valH_eval, *valB_eval);
    
    ROL::Ptr<Intrepid::FieldContainer<Real> > valBeta_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    fe_->evaluateValue(valBeta_eval, Z[0]);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
	(*valBeta_eval)(i,j) = nob1_*std::exp((*valBeta_eval)(i,j));
      }
    }
    
    ROL::Ptr<Intrepid::FieldContainer<Real> > val_xvel_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    fe_->evaluateValue(val_xvel_eval, U[0]);
    ROL::Ptr<Intrepid::FieldContainer<Real> > grad_xvel_eval
      = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p, d);
    fe_->evaluateGradient(grad_xvel_eval, U[0]);
    ROL::Ptr<Intrepid::FieldContainer<Real> > val_yvel_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    fe_->evaluateValue(val_yvel_eval, U[1]);
    ROL::Ptr<Intrepid::FieldContainer<Real> > grad_yvel_eval
      = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p, d);
    fe_->evaluateGradient(grad_yvel_eval, U[1]);
    
    // MULTIPLY beta * xvel
    Intrepid::FieldContainer<Real> beta_phi(c, f, p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(beta_phi,
								*valBeta_eval,
								*(fe_->N()));
    Intrepid::FunctionSpaceTools::integrate<Real>(*J[0][0], beta_phi, *(fe_->NdetJ()), Intrepid::COMP_CPP, false);

    Intrepid::FunctionSpaceTools::integrate<Real>(*J[1][1], beta_phi, *(fe_->NdetJ()), Intrepid::COMP_CPP, false);

    // Nonlinear terms
    Intrepid::FieldContainer<Real> epsilon(c, p);
    for(int i = 0; i < c; i++)
      {
	for(int j = 0; j < p; j++)
	  {
	    //epsilon(i,j) = std::pow((*grad_xvel_eval)(i,j,0),2.0) + std::pow((*grad_yvel_eval)(i,j,1),2.0) + 0.25*( (*grad_xvel_eval)(i,j,1) + (*grad_yvel_eval)(i,j,0) ) + (*grad_xvel_eval)(i,j,0)*(*grad_yvel_eval)(i,j,1);
	    epsilon(i,j) = 1.0 + nob2_*std::pow((*grad_xvel_eval)(i,j,0),2.0) + nob2_*std::pow((*grad_yvel_eval)(i,j,1),2.0) + nob3_*0.25*( (*grad_xvel_eval)(i,j,1) );
	  }
      }
    
    Intrepid::FieldContainer<Real> nu(c, p);
    for(int i = 0; i < c; i++)
      {
	for(int j = 0; j < p; j++)
	  {
	    nu(i,j) = 0.5*B_*(1.0/std::cbrt(epsilon(i,j))); 
	    //nu(i,j) = 0.5*B_*std::pow(epsilon(i,j),2.0);
	  }
      }

    Intrepid::FieldContainer<Real> zx(c, p, 2);
    Intrepid::FieldContainer<Real> zy(c, p, 2);
    for(int i = 0; i < c; i++)
      {
	for(int j = 0; j < p; j++)
	  {
	    zx(i,j,0) = 4.0*(*grad_xvel_eval)(i,j,0) + 2.0*(*grad_yvel_eval)(i,j,1);
	    zx(i,j,1) = (*grad_xvel_eval)(i,j,1) + (*grad_yvel_eval)(i,j,0);
	    zy(i,j,0) = 4.0*(*grad_yvel_eval)(i,j,1) + 2.0*(*grad_xvel_eval)(i,j,0);
	    zy(i,j,1) = (*grad_xvel_eval)(i,j,1) + (*grad_yvel_eval)(i,j,0);
	    // zx(i,j,0) = 1.0;
	    // zx(i,j,1) = 1.0;
	    // zy(i,j,0) = 1.0;
	    // zy(i,j,1) = 1.0;
	  }
      }

    Intrepid::FieldContainer<Real> d_phi(*fe_->gradN()); // (c,f,p,d) 

    Intrepid::FieldContainer<Real> nu_diff_x(c,f,p);
    Intrepid::FieldContainer<Real> nu_diff_y(c,f,p);
    Intrepid::FieldContainer<Real> zx_diff_x(c,f,p,d); 
    Intrepid::FieldContainer<Real> zx_diff_y(c,f,p,d); 
    Intrepid::FieldContainer<Real> zy_diff_x(c,f,p,d); 
    Intrepid::FieldContainer<Real> zy_diff_y(c,f,p,d);    

    for(int i = 0; i < c; i++)
      {
	for(int j = 0; j < p; j++)
	  {
	    Real val1 = (-1.0/3.0)*nu(i,j)*(1.0/epsilon(i,j)); 
	    //Real val1 = (2.0)*nu(i,j)*(1.0/epsilon(i,j));
	    for(int k = 0; k < f; k++)
	      {
		//nu_diff_x(i,k,j) = val1*(2.0*(*grad_xvel_eval)(i,j,0)+(*grad_yvel_eval)(i,j,1))*d_phi(i,k,j,0) + val1*0.25*d_phi(i,k,j,1);
		//nu_diff_y(i,k,j) = val1*0.25*d_phi(i,k,j,0) + val1*(2.0*(*grad_yvel_eval)(i,j,1)+(*grad_xvel_eval)(i,j,0))*d_phi(i,k,j,1);
		nu_diff_x(i,k,j) = val1*( nob2_*2.0*(*grad_xvel_eval)(i,j,0)*d_phi(i,k,j,0) + nob3_*0.25*d_phi(i,k,j,1) );
		nu_diff_y(i,k,j) = val1*( nob2_*2.0*(*grad_yvel_eval)(i,j,1)*d_phi(i,k,j,1) );

		zx_diff_x(i,k,j,0) = 4.0*d_phi(i,k,j,0);
		zx_diff_x(i,k,j,1) = d_phi(i,k,j,1);
		zx_diff_y(i,k,j,0) = 2.0*d_phi(i,k,j,1);
		zx_diff_y(i,k,j,1) = d_phi(i,k,j,0);

		zy_diff_x(i,k,j,0) = 2.0*d_phi(i,k,j,0);
		zy_diff_x(i,k,j,1) = d_phi(i,k,j,1);
		zy_diff_y(i,k,j,0) = 4.0*d_phi(i,k,j,1);
		zy_diff_y(i,k,j,1) = d_phi(i,k,j,0);

		// zx_diff_x(i,k,j,0) = 0.0;
		// zx_diff_x(i,k,j,1) = 0.0;
		// zx_diff_y(i,k,j,0) = 0.0;
		// zx_diff_y(i,k,j,1) = 0.0;

		// zy_diff_x(i,k,j,0) = 0.0;
		// zy_diff_x(i,k,j,1) = 0.0;
		// zy_diff_y(i,k,j,0) = 0.0;
		// zy_diff_y(i,k,j,1) = 0.0;
	      }
	  }
      }

    Intrepid::FieldContainer<Real> wx1_diff_x_1(c,f,p); 
    Intrepid::FieldContainer<Real> wx1_diff_x_2(c,f,p); 
    Intrepid::FieldContainer<Real> wx2_diff_x_1(c,f,p); 
    Intrepid::FieldContainer<Real> wx2_diff_x_2(c,f,p); 
    Intrepid::FieldContainer<Real> wx1_diff_y_1(c,f,p); 
    Intrepid::FieldContainer<Real> wx1_diff_y_2(c,f,p); 
    Intrepid::FieldContainer<Real> wx2_diff_y_1(c,f,p); 
    Intrepid::FieldContainer<Real> wx2_diff_y_2(c,f,p); 
    Intrepid::FieldContainer<Real> wy1_diff_x_1(c,f,p); 
    Intrepid::FieldContainer<Real> wy1_diff_x_2(c,f,p); 
    Intrepid::FieldContainer<Real> wy2_diff_x_1(c,f,p); 
    Intrepid::FieldContainer<Real> wy2_diff_x_2(c,f,p); 
    Intrepid::FieldContainer<Real> wy1_diff_y_1(c,f,p); 
    Intrepid::FieldContainer<Real> wy1_diff_y_2(c,f,p); 
    Intrepid::FieldContainer<Real> wy2_diff_y_1(c,f,p); 
    Intrepid::FieldContainer<Real> wy2_diff_y_2(c,f,p); 
 

    for(int i = 0; i < c; i++)
      {
	for(int j = 0; j < p; j++)
	  {
	    Real h = (*valH_eval)(i,j);
	    for(int k = 0; k < f; k++)
	      {
		wx1_diff_x_1(i,k,j) = nu_diff_x(i,k,j)*h*zx(i,j,0);
		wx1_diff_x_2(i,k,j) = nu(i,j)*h*zx_diff_x(i,k,j,0);

		wx2_diff_x_1(i,k,j) = nu_diff_x(i,k,j)*h*zx(i,j,1);
		wx2_diff_x_2(i,k,j) = nu(i,j)*h*zx_diff_x(i,k,j,1);

		wx1_diff_y_1(i,k,j) = nu_diff_y(i,k,j)*h*zx(i,j,0);
		wx1_diff_y_2(i,k,j) = nu(i,j)*h*zx_diff_y(i,k,j,0);

		wx2_diff_y_1(i,k,j) = nu_diff_y(i,k,j)*h*zx(i,j,1);
		wx2_diff_y_2(i,k,j) = nu(i,j)*h*zx_diff_y(i,k,j,1);

		wy1_diff_x_1(i,k,j) = nu_diff_x(i,k,j)*h*zy(i,j,0);
		wy1_diff_x_2(i,k,j) = nu(i,j)*h*zy_diff_x(i,k,j,0);

		wy2_diff_x_1(i,k,j) = nu_diff_x(i,k,j)*h*zy(i,j,1);
		wy2_diff_x_2(i,k,j) = nu(i,j)*h*zy_diff_x(i,k,j,1);

		wy1_diff_y_1(i,k,j) = nu_diff_y(i,k,j)*h*zy(i,j,0);
		wy1_diff_y_2(i,k,j) = nu(i,j)*h*zy_diff_y(i,k,j,0);

		wy2_diff_y_1(i,k,j) = nu_diff_y(i,k,j)*h*zy(i,j,1);
		wy2_diff_y_2(i,k,j) = nu(i,j)*h*zy_diff_y(i,k,j,1);
	      }
	  }
      }

    Intrepid::FunctionSpaceTools::integrate<Real>(*J[0][0], wx1_diff_x_1, *(fe_->DNDdetJ(0)), Intrepid::COMP_CPP, true);
    Intrepid::FunctionSpaceTools::integrate<Real>(*J[0][0], wx1_diff_x_2, *(fe_->DNDdetJ(0)), Intrepid::COMP_CPP, true); // Do not change this
    Intrepid::FunctionSpaceTools::integrate<Real>(*J[0][0], *(fe_->DNDdetJ(1)), wx2_diff_x_1, Intrepid::COMP_CPP, true);
    Intrepid::FunctionSpaceTools::integrate<Real>(*J[0][0], wx2_diff_x_2, *(fe_->DNDdetJ(1)), Intrepid::COMP_CPP, true); // Do not change this

    Intrepid::FunctionSpaceTools::integrate<Real>(*J[0][1], *(fe_->DNDdetJ(0)), wx1_diff_y_1, Intrepid::COMP_CPP, true);
    Intrepid::FunctionSpaceTools::integrate<Real>(*J[0][1], *(fe_->DNDdetJ(0)), wx1_diff_y_2, Intrepid::COMP_CPP, true); // Do not change this
    Intrepid::FunctionSpaceTools::integrate<Real>(*J[0][1], *(fe_->DNDdetJ(1)), wx2_diff_y_1, Intrepid::COMP_CPP, true);
    Intrepid::FunctionSpaceTools::integrate<Real>(*J[0][1], *(fe_->DNDdetJ(1)), wx2_diff_y_2, Intrepid::COMP_CPP, true); // Do not change this

    Intrepid::FunctionSpaceTools::integrate<Real>(*J[1][0], wy1_diff_x_1, *(fe_->DNDdetJ(0)), Intrepid::COMP_CPP, true);
    Intrepid::FunctionSpaceTools::integrate<Real>(*J[1][0], wy1_diff_x_2, *(fe_->DNDdetJ(0)), Intrepid::COMP_CPP, true); // Do not change this
    Intrepid::FunctionSpaceTools::integrate<Real>(*J[1][0], *(fe_->DNDdetJ(1)), wy2_diff_x_1, Intrepid::COMP_CPP, true);
    Intrepid::FunctionSpaceTools::integrate<Real>(*J[1][0], wy2_diff_x_2, *(fe_->DNDdetJ(1)), Intrepid::COMP_CPP, true); // Do not change this

    Intrepid::FunctionSpaceTools::integrate<Real>(*J[1][1], *(fe_->DNDdetJ(0)), wy1_diff_y_1, Intrepid::COMP_CPP, true);
    Intrepid::FunctionSpaceTools::integrate<Real>(*J[1][1], *(fe_->DNDdetJ(0)), wy1_diff_y_2, Intrepid::COMP_CPP, true); // Do not change this
    Intrepid::FunctionSpaceTools::integrate<Real>(*J[1][1], *(fe_->DNDdetJ(1)), wy2_diff_y_1, Intrepid::COMP_CPP, true);
    Intrepid::FunctionSpaceTools::integrate<Real>(*J[1][1], *(fe_->DNDdetJ(1)), wy2_diff_y_2, Intrepid::COMP_CPP, true); // Do not change this

    // APPLY DIRICHLET CONDITIONS
    int numLocalSideIds = bdryCellLocIds_[0].size();
    for (int j = 0; j < numLocalSideIds; ++j) {
      int numCellsSide = bdryCellLocIds_[0][j].size();
      int numBdryDofs = fidx_[j].size();
      for (int k = 0; k < numCellsSide; ++k) {
	int cidx = bdryCellLocIds_[0][j][k];
	for (int l = 0; l < numBdryDofs; ++l) {
	  for (int m = 0; m < f; ++m) {
	    (*J[0][0])(cidx,fidx_[j][l],m) = static_cast<Real>(0);
	    (*J[1][1])(cidx,fidx_[j][l],m) = static_cast<Real>(0);
	    (*J[1][0])(cidx,fidx_[j][l],m) = static_cast<Real>(0);
	    (*J[0][1])(cidx,fidx_[j][l],m) = static_cast<Real>(0);
	  }
	  (*J[0][0])(cidx,fidx_[j][l],fidx_[j][l]) = static_cast<Real>(1);
	  (*J[1][1])(cidx,fidx_[j][l],fidx_[j][l]) = static_cast<Real>(1);
	}
      }
    }

    // Combine the jacobians
    fieldHelper_->combineFieldCoeff(jac, J);
  }

  void Jacobian_2(ROL::Ptr<Intrepid::FieldContainer<Real> > & jac,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // Retrieve dimensions.
    int c  = u_coeff->dimension(0);
    int p  = cellCub_->getNumPoints();
    int f = basisPtr_->getCardinality();
    int d = fe_->gradN()->dimension(3);
    
    // INITILAIZE JACOBIAN
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > J(2);
    for (int i = 0; i < 2; ++i) {
      J[i].resize(2);
      for (int j = 0; j < 2; j++)
	{
	  J[i][j] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
	}
    }

    // Split u_coeff into components.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > U;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > Z;
    fieldHelper_->splitFieldCoeff(Z, z_coeff);
    
    ROL::Ptr<Intrepid::FieldContainer<Real> > valB_eval
      = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    ROL::Ptr<Intrepid::FieldContainer<Real> > valS_eval
      = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    ROL::Ptr<Intrepid::FieldContainer<Real> > valSx_eval
      = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    ROL::Ptr<Intrepid::FieldContainer<Real> > valSy_eval
      = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    compute_parameters(valB_eval,valS_eval,valSx_eval,valSy_eval);
    
    ROL::Ptr<Intrepid::FieldContainer<Real> > valH_eval
      = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    Intrepid::RealSpaceTools<Real>::add(*valH_eval, *valS_eval);
    Intrepid::RealSpaceTools<Real>::subtract(*valH_eval, *valB_eval);
    
    ROL::Ptr<Intrepid::FieldContainer<Real> > valBeta_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    fe_->evaluateValue(valBeta_eval, Z[0]);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
	(*valBeta_eval)(i,j) = nob1_*std::exp((*valBeta_eval)(i,j));
      }
    }
    
    ROL::Ptr<Intrepid::FieldContainer<Real> > val_xvel_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    fe_->evaluateValue(val_xvel_eval, U[0]);
    ROL::Ptr<Intrepid::FieldContainer<Real> > grad_xvel_eval
      = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p, d);
    fe_->evaluateGradient(grad_xvel_eval, U[0]);
    ROL::Ptr<Intrepid::FieldContainer<Real> > val_yvel_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    fe_->evaluateValue(val_yvel_eval, U[1]);
    ROL::Ptr<Intrepid::FieldContainer<Real> > grad_yvel_eval
      = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p, d);
    fe_->evaluateGradient(grad_yvel_eval, U[1]);
    
    // MULTIPLY beta * xvel
    Intrepid::FieldContainer<Real> beta_xvel(c, p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(beta_xvel,
                                                               *val_xvel_eval,
                                                               *valBeta_eval);
    // MULTIPLY beta * xvel * phi
    Intrepid::FieldContainer<Real> beta_xvel_phi(c, f, p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(beta_xvel_phi,
								beta_xvel,
								*(fe_->N()));
    Intrepid::FunctionSpaceTools::integrate<Real>(*J[0][0], beta_xvel_phi, *(fe_->NdetJ()), Intrepid::COMP_CPP, false);

    // MULTIPLY beta * yvel
    Intrepid::FieldContainer<Real> beta_yvel(c, p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(beta_yvel,
                                                               *val_yvel_eval,
                                                               *valBeta_eval);
    // MULTIPLY beta * yvel * phi
    Intrepid::FieldContainer<Real> beta_yvel_phi(c, f, p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(beta_yvel_phi,
								beta_yvel,
								*(fe_->N()));
    Intrepid::FunctionSpaceTools::integrate<Real>(*J[1][0], beta_yvel_phi, *(fe_->NdetJ()), Intrepid::COMP_CPP, false);

    // APPLY DIRICHLET CONDITIONS
    int numLocalSideIds = bdryCellLocIds_[0].size();
    for (int j = 0; j < numLocalSideIds; ++j) {
      int numCellsSide = bdryCellLocIds_[0][j].size();
      int numBdryDofs = fidx_[j].size();
      for (int k = 0; k < numCellsSide; ++k) {
	int cidx = bdryCellLocIds_[0][j][k];
	for (int l = 0; l < numBdryDofs; ++l) {
	  for (int m = 0; m < f; ++m) {
	    (*J[0][0])(cidx,fidx_[j][l],m) = static_cast<Real>(0);
	    (*J[1][1])(cidx,fidx_[j][l],m) = static_cast<Real>(0);
	    (*J[1][0])(cidx,fidx_[j][l],m) = static_cast<Real>(0);
	    (*J[0][1])(cidx,fidx_[j][l],m) = static_cast<Real>(0);
	  }
	}
      }
    }
 
    // Combine the residuals.
    fieldHelper_->combineFieldCoeff(jac, J);
  }
  
  void Jacobian_3(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > & jac,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_shallow_ice::Jacobian_3): Jacobian is zero.");
  }

  void Hessian_11(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> (PDE_shallow_ice::Hessian_11): Not implemented.");
  }

  void Hessian_12(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> (PDE_shallow_ice::Hessian_12): Not implemented.");
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
    throw Exception::NotImplemented(">>> (PDE_shallow_ice::Hessian_21): Not implemented.");
  }

  void Hessian_22(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> (PDE_shallow_ice::Hessian_22): Not implemented.");
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
    bdryCellDofValues_xvel.resize(numSidesets);
    bdryCellDofValues_yvel.resize(numSidesets);
    feBdry_.resize(numSidesets); 
    for (int i=0; i<numSidesets; ++i) {
      int numLocSides = bdryCellLocIds_[i].size();
      bdryCellDofValues_xvel[i].resize(numLocSides);
      bdryCellDofValues_yvel[i].resize(numLocSides);
      feBdry_[i].resize(numLocSides);
      for (int j=0; j<numLocSides; ++j) {
        int c = bdryCellLocIds_[i][j].size();
        int f = basisPtr_->getCardinality();
	if(c > 0)
	  {
	    feBdry_[i][j] = ROL::makePtr<FE<Real> >(bdryCellNodes_[i][j],basisPtr_,bdryCub_,j);
	  }
	bdryCellDofValues_xvel[i][j] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
	bdryCellDofValues_yvel[i][j] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
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
            (*bdryCellDofValues_xvel[i][j])(k, l) = DirichletFunc_xvel(dofpoint, i, j);
	    (*bdryCellDofValues_yvel[i][j])(k, l) = DirichletFunc_yvel(dofpoint, i, j);
          }
        }
      }
    }

  }
  
  void setFieldPattern(const std::vector<std::vector<int> > & fieldPattern) {
    fieldPattern_ = fieldPattern;
    fieldHelper_ = ROL::makePtr<FieldHelper<Real>>(numFields_, numDofs_, numFieldDofs_, fieldPattern_);
  }

  const ROL::Ptr<FE<Real> > getFE(void) const {
    return fe_;
  }

  std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>>> getFields() {
    return basisPtrs_;
  }

  const std::vector<std::vector<std::vector<int> > > getBdryCellLocIds(void) const {
    return bdryCellLocIds_;
  }

  const ROL::Ptr<FieldHelper<Real> > getFieldHelper(void) const {
    return fieldHelper_;
  }

private:
  
  void compute_parameters(ROL::Ptr<Intrepid::FieldContainer<Real> > & b, ROL::Ptr<Intrepid::FieldContainer<Real> > & S,
			  ROL::Ptr<Intrepid::FieldContainer<Real> > & Sx, ROL::Ptr<Intrepid::FieldContainer<Real> > & Sy) const {
    int c = fe_->gradN()->dimension(0);
    int p = fe_->gradN()->dimension(2);
    int d = fe_->gradN()->dimension(3);
    std::vector<Real> pt(d);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for ( int k = 0; k < d; ++k) {
          pt[k] = (*fe_->cubPts())(i,j,k);
        }
        (*b)(i,j) = bedrock(pt);
	Surface_Height((*S)(i,j), (*Sx)(i,j), (*Sy)(i,j), pt);
      }
    }

  }

  Real bedrock(const std::vector<Real> & pt) const {
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
		val += bedrock_topo_coeff_[k]*Linear_FE_Basis_Fun_Eval(x,i_index,0)*Linear_FE_Basis_Fun_Eval(y,j_index,0);
	      }
	  }
      }

    return val;
  }

  void Surface_Height(Real & S, Real & Sx, Real & Sy, const std::vector<Real> & pt) const {
    S = 0.0;
    Sx = 0.0;
    Sy = 0.0;
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
		S += surface_height_coeff_[k]*Linear_FE_Basis_Fun_Eval(x,i_index,0)*Linear_FE_Basis_Fun_Eval(y,j_index,0);
		Sx += surface_height_coeff_[k]*Linear_FE_Basis_Fun_Eval(x,i_index,1)*Linear_FE_Basis_Fun_Eval(y,j_index,0);
		Sy += surface_height_coeff_[k]*Linear_FE_Basis_Fun_Eval(x,i_index,0)*Linear_FE_Basis_Fun_Eval(y,j_index,1);
	      }
	  }
      }
  }

  Real Linear_FE_Basis_Fun_Eval(const Real & x, int i, int deriv) const
  {
    // Evaluates ith 1D FE basis function at point x
    Real val = 0.0;
    Real xi = static_cast<Real>(i)*h_;
    Real dist = std::abs(x-xi);
    if( dist < h_ )
      {
	if(deriv == 0)
	  {
	    val = 1.0 - dist/h_;
	  }
	else if(deriv == 1)
	  {
	    if(x < xi)
	      {	  
		val = 1.0/h_;
	      }
	    if(x > xi)
	      {
		val = -1.0/h_;
	      }
	  }
      }
    return val;
  }

  Real DirichletFunc_xvel(const std::vector<Real> & coords, int sideset, int locSideId) const {
    std::vector<int> I = std::vector<int>(N_+1);
    Real v = 0;
    if(locSideId == 0) // y=0 boundary
      {
	for(int i = 0; i < N_+1; i++)
	  {
	    I[i] = i;
	  }
	v = coords[0]/width_;
      }
    else if(locSideId == 1) // x=1 boundary
      {
	for(int i = 0; i < N_+1; i++)
	  {
	    I[i] = (N_+1)*(i+1)-1;
	  }
	v = coords[1]/height_;
      }
    else if(locSideId == 2) // y=1 boundary
      {
      	for(int i = 0; i < N_+1; i++)
	  {
	    I[i] = (N_+1)*N_ + i;
	  }
	v = coords[0]/width_;
      }
    else if(locSideId == 3) // x=0 boundary
      {
	for(int i = 0; i < N_+1; i++)
	  {
	    I[i] = (N_+1)*i;
	  }
	v = coords[1]/height_;
      }
    
    Real val = 0.0;
    for(int i = 0; i < N_+1; i++)
      {
	val += xvel_coeff_[I[i]]*Linear_FE_Basis_Fun_Eval(v,i,0);
      }

    return val;
  }

  Real DirichletFunc_yvel(const std::vector<Real> & coords, int sideset, int locSideId) const {
     std::vector<int> I = std::vector<int>(N_+1);
    Real v = 0;
    if(locSideId == 0) // y=0 boundary
      {
	for(int i = 0; i < N_+1; i++)
	  {
	    I[i] = i;
	  }
	v = coords[0]/width_;
      }
    else if(locSideId == 1) // x=1 boundary
      {
	for(int i = 0; i < N_+1; i++)
	  {
	    I[i] = (N_+1)*(i+1)-1;
	  }
	v = coords[1]/height_;
      }
    else if(locSideId == 2) // y=1 boundary
      {
      	for(int i = 0; i < N_+1; i++)
	  {
	    I[i] = (N_+1)*N_ + i;
	  }
	v = coords[0]/width_;
      }
    else if(locSideId == 3) // x=0 boundary
      {
	for(int i = 0; i < N_+1; i++)
	  {
	    I[i] = (N_+1)*i;
	  }
	v = coords[1]/height_;
      }
    
    Real val = 0.0;
    for(int i = 0; i < N_+1; i++)
      {
	val += yvel_coeff_[I[i]]*Linear_FE_Basis_Fun_Eval(v,i,0);
      }

    return val;
  }

}; // PDE_shallow_ice

#endif
