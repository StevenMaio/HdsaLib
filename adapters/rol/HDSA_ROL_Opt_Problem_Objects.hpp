#ifndef OPT_PROBLEM_OBJECTS_ROL_HPP
#define OPT_PROBLEM_OBJECTS_ROL_HPP

// ROL instantiation of Opt_Problem_Objects

template <class RealT>
class Opt_Problem_Objects_ROL : public HDSA::Opt_Problem_Objects<RealT> {

private:
  const HDSA::Ptr<HDSA::ParameterList > parlist_;
  HDSA::nullstream bhs_;
  HDSA::Ptr<std::ostream> outStream_;

public:

  Opt_Problem_Objects_ROL(const HDSA::Ptr<HDSA::ParameterList > & parlist, const HDSA::Ptr<const HDSA::Comm<int> > & comm): parlist_(parlist)
  { 
    int myRank = comm->getRank();
    if(myRank == 0)
      {
	outStream_ = HDSA::makePtrFromRef(std::cout);
      }
    else
      {	
	outStream_ = HDSA::makePtrFromRef(bhs_);
      }
  }

  virtual ~Opt_Problem_Objects_ROL()
  { }

  virtual HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > Construct_Opt_Problem_Objects(const HDSA::Ptr<HDSA::Vector<RealT> > & theta, const HDSA::Ptr<const HDSA::Comm<int> > & comm) = 0;

  void Solve_Optimization_Problem()
  {
    bool  use_Full_Space = parlist_->sublist("Problem").get("Full Space",false);
    
    if(use_Full_Space)
      {
	HDSA::Ptr<ROL::Objective_SimOpt<RealT> > rol_obj;
	try {
	  rol_obj =  dynamic_cast<ROL_FS_Objective<RealT>&>(*HDSA::Opt_Problem_Objects<RealT>::fs_obj).get_objective_function();
	} catch(...) {
	  rol_obj =  dynamic_cast<ROL_FS_Objective_Model_Error<RealT>&>(*HDSA::Opt_Problem_Objects<RealT>::fs_obj).get_objective_function();
	}

	const HDSA::Ptr<ROL::Constraint_SimOpt<RealT> > rol_con =  dynamic_cast<const ROL_Constraint<RealT>&>(*HDSA::Opt_Problem_Objects<RealT>::con).get_constraint();
	HDSA::Ptr<ROL::Vector<RealT> > u_rol = dynamic_cast<const ROL_Vector<RealT>&>(*HDSA::Opt_Problem_Objects<RealT>::u).get_rol_vec();
	HDSA::Ptr<ROL::Vector<RealT> > z_rol = dynamic_cast<const ROL_Vector<RealT>&>(*HDSA::Opt_Problem_Objects<RealT>::z).get_rol_vec();
	HDSA::Ptr<ROL::Vector<RealT> > lambda_rol = dynamic_cast<const ROL_Vector<RealT>&>(*HDSA::Opt_Problem_Objects<RealT>::lambda).get_rol_vec();
	ROL::Vector_SimOpt<RealT> x_rol(u_rol,z_rol);

	HDSA::Ptr<ROL::Vector<RealT> > rp = u_rol->dual().clone();
	ROL::OptimizationProblem<RealT> optProb(rol_obj, HDSA::makePtrFromRef(x_rol), rol_con, rp);
	ROL::OptimizationSolver<RealT> optSolver(optProb, *parlist_);
	std::clock_t timer = std::clock();
	optSolver.solve(*outStream_);
	*outStream_ << "Optimization time: "
		    << static_cast<RealT>(std::clock()-timer)/static_cast<RealT>(CLOCKS_PER_SEC)
		    << " seconds." << std::endl << std::endl;
      }
    else
      {
	//const HDSA::Ptr<ROL::Objective<RealT> > rol_obj =  dynamic_cast<const ROL_RS_Objective<RealT>&>(*HDSA::Opt_Problem_Objects<RealT>::rs_obj).get_objective_function();
	HDSA::Ptr<ROL::Objective<RealT> > rol_obj;
	try {
	  rol_obj =  dynamic_cast<ROL_RS_Objective<RealT>&>(*HDSA::Opt_Problem_Objects<RealT>::rs_obj).get_objective_function();
	} catch(...) {
	  rol_obj =  dynamic_cast<ROL_RS_Objective_Model_Error<RealT>&>(*HDSA::Opt_Problem_Objects<RealT>::rs_obj).get_objective_function();
	}
	HDSA::Ptr<ROL::Vector<RealT> > z_rol = dynamic_cast<const ROL_Vector<RealT>&>(*HDSA::Opt_Problem_Objects<RealT>::z).get_rol_vec();
	
	// Build optimization problem and check derivatives
	ROL::OptimizationProblem<RealT> optProb(rol_obj,z_rol);
	// Build optimization solver and solve
	ROL::OptimizationSolver<RealT> optSolver(optProb,*parlist_);
	std::clock_t timer = std::clock();
	optSolver.solve(*outStream_);
	*outStream_ << "Trust Region Time: "
		    << static_cast<RealT>(std::clock()-timer)/static_cast<RealT>(CLOCKS_PER_SEC)
		    << " seconds." << std::endl << std::endl;
      }
  }

  virtual void Load_Optimal_Solution() = 0;

  virtual void Write_Optimal_Solution() = 0;
  
};


#endif
