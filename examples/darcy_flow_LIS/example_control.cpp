#include "Teuchos_GlobalMPISession.hpp"

#include "ROL_TpetraMultiVector.hpp"
#include "ROL_Algorithm.hpp"
#include "ROL_Reduced_Objective_SimOpt.hpp"
#include "ROL_Bounds.hpp"

#include "../../../PDE-OPT/TOOLS/solver.cpp"
#include "../../../PDE-OPT/TOOLS/solver_def.hpp"
#include "../../../PDE-OPT/TOOLS/assembler.cpp"
#include "../../../PDE-OPT/TOOLS/assembler_def.hpp"
#include "../../../PDE-OPT/TOOLS/meshmanager.hpp"
#include "../../../PDE-OPT/TOOLS/linearpdeconstraint.hpp"
#include "../../../PDE-OPT/TOOLS/pdeobjective.hpp"
#include "../../../PDE-OPT/TOOLS/pdevector.hpp"
#include "pde_darcy_flow.hpp"
#include "obj_darcy_flow.hpp"
#include "mesh_darcy_flow.hpp"

#include "../../src/source_file.hpp"
#include "Opt_Problem_Objects_darcy_flow_LIS.hpp"

typedef double RealT;

int main(int argc, char *argv[]) {

  HDSA::nullstream bhs;
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &bhs);
  HDSA::Ptr<const HDSA::Comm<int> > comm = HDSA::makePtr<HDSA::Comm<int> >();
  std::string filename = "input.xml";
  HDSA::Ptr<HDSA::ParameterList> parlist = HDSA::makePtr<HDSA::ParameterList>();
  HDSA::updateParametersFromXmlFile( filename, *parlist );

  int L = parlist->sublist("Problem").get("Number of Uncertainty Basis Function", 10);
  int dim = (L+1)*(L+1);
  HDSA::Ptr<HDSA::Vector<RealT> > theta = HDSA::makePtr<Std_Vector<RealT> >(dim);

  HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > OP_Objects = HDSA::makePtr<Opt_Problem_Objects_darcy_flow_LIS<RealT> >(parlist,theta,comm);
  OP_Objects->Solve_Optimization_Problem();
  OP_Objects->Write_Optimal_Solution();
  
  return 0;
}
