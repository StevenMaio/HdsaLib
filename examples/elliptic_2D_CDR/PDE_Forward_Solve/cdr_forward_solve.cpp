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
void Set_source(HDSA::Ptr<Tpetra::MultiVector<> > & z_ptr, Teuchos::RCP<Teuchos::ParameterList> & parlist)
{
    std::ifstream inputFile("source_input.txt");          
    RealT value;
    int count = 0;
    // read the elements in the file into a vector  
    // test file open   
    if (inputFile) {   
      while ( inputFile >> value ) {
	z_ptr->replaceGlobalValue(count,0,value);
	count += 1;
      }
    }
    else
      {
	std::cout << "Error loading the source input" << std::endl;
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
  HDSA::Ptr<Teuchos::ParameterList> parlist = Teuchos::getParametersFromXmlFile("input_forward_solve.xml");

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
  Set_source<RealT>(z_ptr,parlist);
  int dim = up->dimension();
  RealT tol = 1.e-8;
  con->solve(*rp, *up, *zp,tol);

  std::string name;
  std::ofstream fout;
  name="state_solution.txt";
  fout.open(name);
  for(int j=0; j < dim; j++)
    {
       fout << std::setprecision(16) << up->dot(*rp->basis(j)) << "  " << std::endl;
    }
  fout.close();

  return 0;
}
