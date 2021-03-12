#include "Teuchos_Comm.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_Version.hpp"

#include "ROL_Stream.hpp"
#include "ROL_ParameterList.hpp"
#include "ROL_Objective_SimOpt.hpp"
#include "ROL_OptimizationSolver.hpp"

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
#include "../../../../PDE-OPT/TOOLS/linearpdeconstraint.hpp"
#include "../../../../PDE-OPT/TOOLS/pdevector.hpp"
#include "../../../../PDE-OPT/TOOLS/pdeobjective.hpp"
#include "../mesh_cdr.hpp"
#include "../pde_cdr.hpp"

template<class RealT>
void Map_Index_to_Coords(int k, std::vector<RealT> & coords, int d, int nx, int ny, RealT nx_float, RealT ny_float)
{
  coords[0] = static_cast<RealT>( (k%(nx+1)) )*(1.0/nx_float);
  coords[1] = static_cast<RealT>( std::floor( static_cast<RealT>(k)/(nx_float+1.0) ) )*(1.0/ny_float);
}

template<class RealT>
RealT Source_Eval(RealT x, RealT y)
{
  RealT val = std::exp(-100 * ( (x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) ) );
  return val;
}

template<class RealT>
void Set_source(HDSA::Ptr<Tpetra::MultiVector<> > & z_ptr, Teuchos::RCP<Teuchos::ParameterList> & parlist)
{
  int nx = parlist->sublist("Geometry").get("NX",0);
  int ny = parlist->sublist("Geometry").get("NY",0);
  int dim = (nx+1)*(ny+1);
  RealT nx_float = static_cast<RealT>(nx);
  RealT ny_float = static_cast<RealT>(ny);
  std::vector<RealT> coords = std::vector<RealT>(2,0.0);
  for(int k = 0; k < dim; k++)
    {
      Map_Index_to_Coords<RealT>(k,coords,dim,nx,ny,nx_float,ny_float);
      z_ptr->replaceGlobalValue(k,0,Source_Eval<RealT>(coords[0],coords[1]));
    }
}

int main(int argc, char *argv[]) {
  using RealT = double;
  
  /*** Initialize communicator. ***/
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  HDSA::Ptr<const Teuchos::Comm<int> > comm
    = Tpetra::getDefaultComm();

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
  HDSA::Ptr<Teuchos::ParameterList> parlist = Teuchos::getParametersFromXmlFile("input_data_generation.xml");

  // Initialize PDE
  HDSA::Ptr<MeshManager<RealT> > meshMgr = HDSA::makePtr<MeshManager_CDR<RealT> >(*parlist);
  HDSA::Ptr<PDE<RealT> > pde = HDSA::makePtr<PDE_CDR<RealT> >(*parlist);
  HDSA::Ptr<PDE_CDR<RealT> > pde_cdr = HDSA::dynamicPtrCast<PDE_CDR<RealT> >(pde);
  HDSA::Ptr<ROL::Constraint_SimOpt<RealT> > con = HDSA::makePtr<Linear_PDE_Constraint<RealT> >(pde,meshMgr,comm,*parlist,*outStream,true);
  // Cast the constraint and get the assembler.
  HDSA::Ptr<Linear_PDE_Constraint<RealT> > pdecon = HDSA::dynamicPtrCast<Linear_PDE_Constraint<RealT> >(con);
  HDSA::Ptr<Assembler<RealT> > assembler = pdecon->getAssembler();
  // Create state vector and set to zeroes
  HDSA::Ptr<Tpetra::MultiVector<> > u_ptr, z_ptr, p_ptr, r_ptr;
  HDSA::Ptr<ROL::Vector<RealT> > up, zp, pp, rp;
  u_ptr  = assembler->createStateVector();   u_ptr->putScalar(0.0);
  z_ptr  = assembler->createControlVector(); z_ptr->putScalar(0.0);
  p_ptr  = assembler->createStateVector();   p_ptr->putScalar(0.0);
  r_ptr  = assembler->createStateVector();   r_ptr->putScalar(0.0);
  up  = HDSA::makePtr<PDE_PrimalSimVector<RealT> >(u_ptr,pde,assembler);
  zp  = HDSA::makePtr<PDE_PrimalOptVector<RealT> >(z_ptr,pde,assembler);
  pp  = HDSA::makePtr<PDE_PrimalSimVector<RealT> >(p_ptr,pde,assembler);
  rp  = HDSA::makePtr<PDE_DualSimVector<RealT> >(r_ptr,pde,assembler);

  /*************************************************************************/
  /***************** GENERATE DATA *****************************************/
  /*************************************************************************/
  
  // 1. Generate data with a FE mesh twice the resolution desired for the inverse problem
  // 2. Run Sparcify_Data.m matlab code to remove unwanted data and write to data.txt
  Set_source<RealT>(z_ptr,parlist);
  int dim = up->dimension();
  RealT tol = 1.e-8;
  con->solve(*rp, *up, *zp,tol);

  std::string name;
  std::ofstream fout;
  name="true_state.txt";
  fout.open(name);
  for(int j=0; j < dim; j++)
    {
       fout << std::setprecision(16) << up->dot(*rp->basis(j)) << "  " << std::endl;
    }
  fout.close();

  name="true_source.txt";
  fout.open(name);
  for(int j=0; j < dim; j++)
    {
       fout << std::setprecision(16) << zp->dot(*rp->basis(j)) << "  " << std::endl;
    }
  fout.close();

  pdecon->outputTpetraData();
  pdecon->getAssembler()->printMeshData(*outStream);

  return 0;
}
