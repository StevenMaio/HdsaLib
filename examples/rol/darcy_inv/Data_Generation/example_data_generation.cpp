#include "Teuchos_GlobalMPISession.hpp"

#include "ROL_TpetraMultiVector.hpp"
#include "ROL_Algorithm.hpp"
#include "ROL_Reduced_Objective_SimOpt.hpp"
#include "ROL_Bounds.hpp"

#include "../../../../../Trilinos/packages/rol/example/PDE-OPT/TOOLS/meshmanager.hpp"
#include "../../../../../Trilinos/packages/rol/example/PDE-OPT/TOOLS/pdeconstraint.hpp"
#include "../../../../../Trilinos/packages/rol/example/PDE-OPT/TOOLS/pdeobjective.hpp"
#include "../../../../../Trilinos/packages/rol/example/PDE-OPT/TOOLS/pdevector.hpp"
#include "../pde_darcy_flow.hpp"
#include "../mesh_darcy_flow.hpp"
#include "HDSA_Stream.hpp"
#include "HDSA_Comm.hpp"
#include "HDSA_Ptr.hpp"
#include "HDSA_ParameterList.hpp"
#include "../../../../../Trilinos/packages/rol/example/PDE-OPT/TOOLS/solver.cpp"
#include "../../../../../Trilinos/packages/rol/example/PDE-OPT/TOOLS/solver_def.hpp"
#include "../../../../../Trilinos/packages/rol/example/PDE-OPT/TOOLS/assembler.cpp"
#include "../../../../../Trilinos/packages/rol/example/PDE-OPT/TOOLS/assembler_def.hpp"

typedef double RealT;

template <class RealT>
void set_perm(HDSA::Ptr<Tpetra::MultiVector<>> &z_ptr)
{
  int num_coeff_load = z_ptr->getGlobalLength();
  // read in data
  std::ifstream in("Log_Permability.txt");
  // read the elements in the file into a vector
  // test file open
  RealT val = 0.0;
  if (in)
  {
    for (int j = 0; j < num_coeff_load; j++)
    {
      in >> val;
      z_ptr->replaceGlobalValue(j, 0, val);
    }
  }
  else
  {
    std::cout << "Error loading the data from Log_Permabilty.txt" << std::endl;
  }
}

template <class RealT>
void set_Parameters(std::vector<RealT> &param)
{
  int theta_dim = param.size();

  // read in data
  std::ifstream in("theta_true.txt");
  // read the elements in the file into a vector
  // test file open
  RealT val = 0.0;
  if (in)
  {
    for (int j = 0; j < theta_dim; j++)
    {
      in >> val;
      param[j] = val;
    }
  }
  else
  {
    std::cout << "Error loading the data from theta_star.txt" << std::endl;
  }
}

int main(int argc, char *argv[])
{

  HDSA::nullstream bhs;
  Teuchos::GlobalMPISession mpiSession(&argc, &argv, &bhs);
  HDSA::Ptr<const HDSA::Comm<int>> comm = HDSA::makePtr<HDSA::Comm<int>>();
  HDSA::Ptr<std::ostream> outStream;
  int myRank = comm->getRank();
  if (myRank == 0)
  {
    outStream = HDSA::makePtrFromRef(std::cout);
  }
  else
  {
    outStream = HDSA::makePtrFromRef(bhs);
  }
  std::string filename = "input.xml";
  HDSA::Ptr<HDSA::ParameterList> parlist = HDSA::makePtr<HDSA::ParameterList>();
  HDSA::updateParametersFromXmlFile(filename, *parlist);

  /*************************************************************************/
  /***************** BUILD GOVERNING PDE ***********************************/
  /*************************************************************************/
  /*** Initialize main data structure. ***/
  HDSA::Ptr<MeshManager<RealT>> meshMgr = HDSA::makePtr<MeshManager_darcy_flow<RealT>>(*parlist);
  HDSA::Ptr<PDE_darcy_flow<RealT>> pde = HDSA::makePtr<PDE_darcy_flow<RealT>>(*parlist);
  HDSA::Ptr<ROL::Constraint_SimOpt<RealT>> con = HDSA::makePtr<PDE_Constraint<RealT>>(pde, meshMgr, comm->Get_Teuchos_Communicator(), *parlist, *outStream);
  HDSA::Ptr<PDE_Constraint<RealT>> pdecon = HDSA::dynamicPtrCast<PDE_Constraint<RealT>>(con);
  HDSA::Ptr<Assembler<RealT>> assembler = pdecon->getAssembler();
  assembler->printMeshData(*outStream);

  int theta_modes = parlist->sublist("Problem").get("Uncertain Modes per dimension", 1);
  std::vector<RealT> param = std::vector<RealT>(std::pow(theta_modes, 2), 0.0);
  set_Parameters<RealT>(param);
  pdecon->setParameter(param);
  con->setSolveParameters(*parlist);

  /*************************************************************************/
  /***************** BUILD VECTORS *****************************************/
  /*************************************************************************/
  // Create state vectors
  HDSA::Ptr<Tpetra::MultiVector<>> u_ptr, p_ptr, r_ptr, z_ptr;
  HDSA::Ptr<ROL::Vector<RealT>> up, zp, pp, rp;
  u_ptr = assembler->createStateVector();
  up = HDSA::makePtr<PDE_PrimalSimVector<RealT>>(u_ptr, pde, assembler, *parlist);
  p_ptr = assembler->createStateVector();
  pp = HDSA::makePtr<PDE_PrimalSimVector<RealT>>(p_ptr, pde, assembler, *parlist);
  // Create residual vector
  r_ptr = assembler->createStateVector();
  rp = HDSA::makePtr<PDE_PrimalSimVector<RealT>>(r_ptr, pde, assembler, *parlist);
  // Create control vectors
  z_ptr = assembler->createControlVector();
  zp = HDSA::makePtr<PDE_PrimalOptVector<RealT>>(z_ptr, pde, assembler, *parlist);

  set_perm<RealT>(z_ptr);
  RealT tol = 1.e-8;
  con->solve(*pp, *up, *zp, tol);

  pdecon->outputTpetraVector(u_ptr, "true_u.txt");
  pdecon->outputTpetraVector(z_ptr, "true_z.txt");
  pdecon->outputTpetraData();

  return 0;
}
