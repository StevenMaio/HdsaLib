#include "Teuchos_GlobalMPISession.hpp"

#include "ROL_Reduced_Objective_SimOpt.hpp"
#include "ROL_OptimizationSolver.hpp"

#include "../../../PDE-OPT/TOOLS/solver.cpp"
#include "../../../PDE-OPT/TOOLS/solver_def.hpp"
#include "../../../PDE-OPT/TOOLS/assembler.cpp"
#include "../../../PDE-OPT/TOOLS/assembler_def.hpp"
#include "../../../PDE-OPT/TOOLS/meshmanager.hpp"
#include "../../../PDE-OPT/TOOLS/linearpdeconstraint.hpp"
#include "../../../PDE-OPT/TOOLS/pdeconstraint.hpp"
#include "../../../PDE-OPT/TOOLS/pdevector.hpp"
#include "../../../PDE-OPT/TOOLS/pdeobjective.hpp"
#include "../../src/source_file.hpp"
#include "mesh.hpp"
#include "operators.hpp"

typedef double RealT;

int main(int argc, char *argv[]) {
  HDSA::nullstream bhs;
  HDSA::Ptr<std::ostream> outStream;
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &bhs);
  HDSA::Ptr<const HDSA::Comm<int> > comm = HDSA::makePtr<HDSA::Comm<int> >();
  std::string filename = "input.xml";
  HDSA::Ptr<HDSA::ParameterList> parlist = HDSA::makePtr<HDSA::ParameterList>();
  HDSA::updateParametersFromXmlFile( filename, *parlist );
  const int myRank = comm->getRank();
  if (myRank == 0)
    {
      outStream = HDSA::makePtrFromRef(std::cout);
    }
  else 
    {
      outStream = HDSA::makePtrFromRef(bhs);
    }

  HDSA::Ptr<MeshManager_Construction<RealT> > meshMgr = HDSA::makePtr<MeshManager_Construction<RealT> >(*parlist);

  HDSA::Ptr<PDE<RealT> > stiff_pde = HDSA::makePtr<PDE_Stiffness_Op<RealT> >(*parlist);
  HDSA::Ptr<ROL::Constraint_SimOpt<RealT> > stiff_con = HDSA::makePtr<Linear_PDE_Constraint<RealT> >(stiff_pde,meshMgr,comm->Get_Teuchos_Communicator(),*parlist,*outStream,true);

  HDSA::Ptr<PDE<RealT> > mass_pde = HDSA::makePtr<PDE_Mass_Op<RealT> >(*parlist);
  HDSA::Ptr<ROL::Constraint_SimOpt<RealT> > mass_con = HDSA::makePtr<Linear_PDE_Constraint<RealT> >(mass_pde,meshMgr,comm->Get_Teuchos_Communicator(),*parlist,*outStream,true);

  // Cast the constraint and get the assembler.
  HDSA::Ptr<Linear_PDE_Constraint<RealT> > pdecon = HDSA::dynamicPtrCast<Linear_PDE_Constraint<RealT> >(stiff_con);
  HDSA::Ptr<Assembler<RealT> > assembler = pdecon->getAssembler();
  HDSA::Ptr<Tpetra::MultiVector<> > zp_ptr  = assembler->createControlVector();   zp_ptr->putScalar(0.0);
  HDSA::Ptr<ROL::Vector<RealT> > zp = HDSA::makePtr<PDE_PrimalOptVector<RealT> >(zp_ptr,stiff_pde,assembler);

  HDSA::Ptr<ROL::Vector<RealT> > Kz = zp->clone();
  HDSA::Ptr<ROL::Vector<RealT> > Mz = zp->clone();
  HDSA::Ptr<ROL::Vector<RealT> > tmp1 = zp->clone();
  int n = zp->dimension();
  
  std::vector<std::vector<RealT> > K;
  std::vector<std::vector<RealT> > M;
  K.resize(n);
  M.resize(n);
  for(int i = 0; i < n; i++)
    {
      K[i].resize(n);
      M[i].resize(n);
    }
  
  RealT tol = 1.e-8;
  for(int j = 0; j < n; j++)
    {
      zp->set(*Kz->basis(j));
      stiff_con->applyJacobian_1(*Kz,*zp,*tmp1,*tmp1,tol);
      mass_con->applyJacobian_1(*Mz,*zp,*tmp1,*tmp1,tol);
      for(int i = 0; i < n; i++)
	{
	  K[i][j] = Kz->dot(*zp->basis(i));
	  M[i][j] = Mz->dot(*zp->basis(i));
	}
    }
  std::string name;
  std::ofstream fout;
  name = "Stiffness_Matrix.txt";
  fout.open(name);
  for(int i = 0; i < n; i++)
    {
      for(int j = 0; j < n; j++)
	{
	  fout << std::setprecision(16) << K[i][j] << "  ";
	}
      fout << "  " << std::endl;
    }
  fout.close();
  name = "Mass_Matrix.txt";
  fout.open(name);
  for(int i = 0; i < n; i++)
    {
      for(int j = 0; j < n; j++)
	{
	  fout << std::setprecision(16) << M[i][j] << "  ";
	}
      fout << "  " << std::endl;
    }
  fout.close();
  
  return 0;
}
