#include "Teuchos_GlobalMPISession.hpp"

#include "ROL_Stream.hpp"
#include "ROL_ParameterList.hpp"
#include "ROL_Objective_SimOpt.hpp"
#include "ROL_OptimizationSolver.hpp"
#include "ROL_ReducedDynamicObjective.hpp"
#include "ROL_Bounds.hpp"
#include "ROL_DynamicConstraintCheck.hpp"
#include "ROL_DynamicObjectiveCheck.hpp"

#include <iostream>
#include <random>
#include <unordered_set>
#include <set>

#include "../../../src/source_file.hpp"

#include "../../../../PDE-OPT/TOOLS/solver.cpp"
#include "../../../../PDE-OPT/TOOLS/solver_def.hpp"
#include "../../../../PDE-OPT/TOOLS/assembler.cpp"
#include "../../../../PDE-OPT/TOOLS/assembler_def.hpp"
#include "../../../../PDE-OPT/TOOLS/meshmanager.hpp"
#include "../../../../PDE-OPT/TOOLS/dynconstraint.hpp"
#include "../../../../PDE-OPT/TOOLS/ltiobjective.hpp"
#include "../../../../PDE-OPT/TOOLS/pdevector.hpp"
#include "../../../../PDE-OPT/TOOLS/pdeobjective.hpp"
#include "../mesh_darcy_adv_diff.hpp"
#include "../pde_darcy_adv_diff.hpp"
#include "../dynpde_darcy_adv_diff.hpp"

template<class RealT>
RealT Permeability_Eval(std::vector<RealT> & coords)
{ 
  int T = 5;
  std::vector<RealT> amp = std::vector<RealT>(T);
  std::vector<RealT> sdx = std::vector<RealT>(T);
  std::vector<RealT> sdy = std::vector<RealT>(T);
  std::vector<RealT> c = std::vector<RealT>(T);
  std::vector<RealT> h = std::vector<RealT>(T);
  std::vector<RealT> k= std::vector<RealT>(T);
  
  amp[0] = -2.0; amp[1] = -2.2; amp[2] = 2.1; amp[3] = 2.0; amp[4] = 1.9;
  sdx[0] = 1.0; sdx[1] = 2.0; sdx[2] = 17.0; sdx[3] = 7.0; sdx[4] = 29.0;
  sdy[0] = 50.0; sdy[1] = 43.0; sdy[2] = 34.0; sdy[3] = 46.0; sdy[4] = 58.0;
  c[0] = 10.0; c[1] = 6.0; c[2] = 8.0; c[3] = 6.0; c[4] = 11.0;
  h[0] = .6; h[1] = .45; h[2] = .27; h[3] = .53; h[4] = .78;
  k[0] = .75; k[1] = .15; k[2] = .43; k[3] = .55; k[4] = .62; 

  RealT val = 0.0;
  for(int i = 0; i < T; i++)
    {
      val += amp[i]*std::exp(-1.0*( sdx[i]*(coords[0]-h[i])*(coords[0]-h[i]) + c[i]*(coords[0]-h[i])*(coords[1]-k[i]) + sdy[i]*(coords[1]-k[i])*(coords[1]-k[i]) ) );
    }

  ////////////////////////
  val = val/1.35;
  //////////////////////
  
  return val;
}

template<class RealT>
void Map_Index_to_Coords(int k, std::vector<RealT> & coords, int d, int nx, int ny, RealT nx_float, RealT ny_float)
{
  coords[0] = static_cast<RealT>( (k%(nx+1)) )*(1.0/nx_float);
  coords[1] = static_cast<RealT>( std::floor( static_cast<RealT>(k)/(nx_float+1.0) ) )*(1.0/ny_float);
}

template<class RealT>
void Set_Permeability(HDSA::Ptr<Tpetra::MultiVector<> > & z_ptr, HDSA::Ptr<HDSA::ParameterList> & parlist)
{
  int nx = parlist->sublist("Geometry").get("NX",0);
  int ny = parlist->sublist("Geometry").get("NY",0);
  int d = (nx+1)*(ny+1);
  RealT nx_float = static_cast<RealT>(nx);
  RealT ny_float = static_cast<RealT>(ny);
  std::vector<RealT> coords = std::vector<RealT>(2,0.0);
  for(int k = 0; k < d; k++)
    {
      Map_Index_to_Coords<RealT>(k,coords,d,nx,ny,nx_float,ny_float);
      z_ptr->replaceGlobalValue(2*k,0,Permeability_Eval<RealT>(coords));
      z_ptr->replaceGlobalValue(2*k+1,0,0.0);
    }
}

int main(int argc, char *argv[]) {
  using RealT = double;
  
  unsigned seed = 13420958;
  std::default_random_engine generator;
  generator.seed(seed);

  /*** Initialize communicator. ***/
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  HDSA::Ptr<const HDSA::Comm<int> > comm = HDSA::makePtr<HDSA::Comm<int> >();

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  const int myRank = comm->getRank();
  HDSA::nullstream bhs;
  HDSA::Ptr<std::ostream> outStream;
  if(myRank == 0)
    {
      outStream = HDSA::makePtrFromRef(std::cout);
    }
  else
    {	
      outStream =  HDSA::makePtrFromRef(bhs);
    }
  
  /*** Read in XML input ***/
  std::string filename = "input_data_generation.xml";
  HDSA::Ptr<HDSA::ParameterList> parlist = HDSA::makePtr<HDSA::ParameterList>();
  HDSA::updateParametersFromXmlFile( filename, *parlist );

  int nt         = parlist->sublist("Time Discretization").get("Number of Time Steps", 100);
  RealT T        = parlist->sublist("Time Discretization").get("End Time",             1.0);
  RealT dt       = T/static_cast<RealT>(nt);
  
  /*************************************************************************/
  /***************** BUILD GOVERNING PDE ***********************************/
  /*************************************************************************/
  /*** Initialize main data structure. ***/
  HDSA::Ptr<MeshManager<RealT> > meshMgr
    = HDSA::makePtr<MeshManager_darcy_adv_diff<RealT> >(*parlist);
  // Initialize PDE describing advection-diffusion equation
  HDSA::Ptr<DynamicPDE_darcy_adv_diff<RealT> > pde
    = HDSA::makePtr<DynamicPDE_darcy_adv_diff<RealT> >(*parlist);
  
  /*************************************************************************/
  /***************** BUILD CONSTRAINT **************************************/
  /*************************************************************************/
  HDSA::Ptr<DynConstraint<RealT> > dyn_con
    = HDSA::makePtr<DynConstraint<RealT> >(pde,meshMgr,comm->Get_Teuchos_Communicator(),*parlist,*outStream);
  dyn_con->getAssembler()->printMeshData(*outStream);
  
  /*************************************************************************/
  /***************** BUILD STATE VECTORS ***********************************/
  /*************************************************************************/
  HDSA::Ptr<Tpetra::MultiVector<> > u0_ptr, uo_ptr, un_ptr,z_ptr, ck_ptr;
  u0_ptr = dyn_con->getAssembler()->createStateVector();
  uo_ptr = dyn_con->getAssembler()->createStateVector();
  un_ptr = dyn_con->getAssembler()->createStateVector();
  z_ptr = dyn_con->getAssembler()->createControlVector();
  ck_ptr = dyn_con->getAssembler()->createResidualVector();
  HDSA::Ptr<ROL::Vector<RealT> > u0, uo, un, ck, z;
  u0 = HDSA::makePtr<PDE_PrimalSimVector<RealT> >(u0_ptr,pde,*dyn_con->getAssembler());
  uo = HDSA::makePtr<PDE_PrimalSimVector<RealT> >(uo_ptr,pde,*dyn_con->getAssembler());
  un = HDSA::makePtr<PDE_PrimalSimVector<RealT> >(un_ptr,pde,*dyn_con->getAssembler());
  ck = HDSA::makePtr<PDE_DualSimVector<RealT> >(ck_ptr,pde,*dyn_con->getAssembler());
  z = HDSA::makePtr<PDE_PrimalOptVector<RealT> >(z_ptr,pde,*dyn_con->getAssembler());
  u0->zero();

  dyn_con->setSolveParameters(*parlist);
  
  /*************************************************************************/
  /***************** BUILD REDUCED COST FUNCTIONAL *************************/
  /*************************************************************************/
  std::vector<ROL::TimeStamp<RealT> > timeStamp(nt);
  for( int k=0; k<nt; ++k ) {
    timeStamp.at(k).t.resize(2);
    timeStamp.at(k).t.at(0) = k*dt;
    timeStamp.at(k).t.at(1) = (k+1)*dt;
  }

  // Run derivative checks                                                                                                                                                                                
  bool checkDeriv = parlist->sublist("Problem").get("Check Derivatives",true);
  if ( checkDeriv )
    {
      uo->randomize(); un->randomize(); z->randomize();
      ROL::ValidateFunction<RealT> validate(1,13,20,11,true,*outStream);
      ROL::DynamicConstraintCheck<RealT>::check(*dyn_con,validate,*uo,*un,*z);
    }

  /*************************************************************************/
  /***************** GENERATE DATA *****************************************/
  /*************************************************************************/
  
  // 1. Generate data with a FE mesh twice the resolution desired for the inverse problem
  // 2. Run Sparcify_Data.m matlab code to remove unwanted data and write to data.txt
  RealT sigma = parlist->sublist("Problem").get("Noise sigma",0.0);  
  std::normal_distribution<RealT> distribution = std::normal_distribution<RealT>(0.0,sigma); 
  std::vector<std::vector<RealT> > data;
  std::vector<std::vector<RealT> > clean_data;
  data.resize(nt);
  clean_data.resize(nt);
  Set_Permeability<RealT>(z_ptr,parlist);
  int dim = u0->dimension();
  uo->set(*u0); un->zero();
  data[0].resize(dim);
  clean_data[0].resize(dim);
  for(int i = 0; i < dim; i++)
    {
      clean_data[0][i] = uo->dot(*un->basis(i));
      data[0][i] = clean_data[0][i]*(1.0 + distribution(generator));
    }
  for(int k = 1; k < nt; k++)
    {
      // Advance time stepper
      dyn_con->solve(*ck, *uo, *un, *z, timeStamp[k-1]);
      uo->set(*un);
      data[k].resize(dim);
      clean_data[k].resize(dim);
      for(int i = 0; i < dim; i++)
  	{
	  clean_data[k][i] = uo->dot(*un->basis(i));
  	  data[k][i] = clean_data[k][i]*(1.0 + distribution(generator));
  	}
    }

  std::string name;
  std::ofstream fout;
  name="noisy_true_state.txt";
  fout.open(name);
  for(int i=0; i < nt; i++)
    {
      for(int j=0; j < dim; j++)
  	{
  	  fout << data[i][j] << std::setw(20);
  	}
      fout<<""<<std::endl;
    }
  fout.close();

  name="clean_true_state.txt";
  fout.open(name);
  for(int i=0; i < nt; i++)
    {
      for(int j=0; j < dim; j++)
  	{
  	  fout << clean_data[i][j] << std::setw(20);
  	}
      fout<<""<<std::endl;
    }
  fout.close();

  dyn_con->outputTpetraVector(z_ptr,"true_perm.txt");
  dyn_con->outputTpetraData();
  dyn_con->getAssembler()->printMeshData(*outStream);

  return 0;
}
