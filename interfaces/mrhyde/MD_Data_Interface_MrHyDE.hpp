#ifndef HDSA_MD_DATA_INTERFACE_MRHYDE_HPP
#define HDSA_MD_DATA_INTERFACE_MRHYDE_HPP

  template <class RealT>
  class MD_Data_Interface_MrHyDE : public HDSA::MD_Data_Interface<RealT> {

  private:
    Teuchos::RCP<MrHyDE::SolverManager<SolverNode> > solve_;
  public:
    MD_Data_Interface_MrHyDE(Teuchos::RCP<MrHyDE::SolverManager<SolverNode> > &solve)
    { 
 	solve_=solve;
    }

    virtual ~MD_Data_Interface_MrHyDE()
    { }

    HDSA::Ptr<HDSA::Vector<RealT> > Load_Optimal_u(void) const{
    int num_coeff_load = 200;
    HDSA::Ptr<HDSA::Vector<RealT> > u_opt = HDSA::makePtr<HDSA::Vector_MrHyDE_State<RealT> >(solve_);

    std::vector<RealT> opt_u_coeff = std::vector<RealT>(num_coeff_load);
    std::ifstream in("u_opt.txt");
    if (in)
      {
        for(int j = 0; j < num_coeff_load; j++)
        {
          in >> opt_u_coeff[j];
        }
      }
    else
     {
        std::cout << "Error loading the data from u_opt.txt" << std::endl;
     }

    HDSA::Vector_MrHyDE_State<RealT> &eu_opt = dynamic_cast<HDSA::Vector_MrHyDE_State<RealT>&>(*u_opt);
    for(int k = 0; k < num_coeff_load; k++)
        {
          eu_opt.mrhyde_state_vec[0][0]->replaceGlobalValue(k,0,opt_u_coeff[k]);
        }
      return u_opt;
    } 

    HDSA::Ptr<HDSA::Vector<RealT> > Load_Optimal_z(void) const{
    int num_coeff_load = 200;
    ROL::Ptr<ROL::Vector<RealT> > z_opt_rol = solve_->params->getCurrentVector().clone();
    MrHyDE_OptVector &z_opt = dynamic_cast<MrHyDE_OptVector&>(*z_opt_rol);
    
    std::vector<RealT> opt_z_coeff = std::vector<RealT>(num_coeff_load);
    std::ifstream in("z_opt.txt");
    if (in)
      {
        for(int j = 0; j < num_coeff_load; j++)
        {
          in >> opt_z_coeff[j];
        }
      }
    else
     {
        std::cout << "Error loading the data from z_opt.txt" << std::endl;
     }

    for(int k = 0; k < num_coeff_load; k++)
        {
        z_opt.getField()[0]->getVector()->replaceGlobalValue(k,0,opt_z_coeff[k]);
        }
        HDSA::Ptr<HDSA::Vector<RealT> > z_opt_hdsa = HDSA::makePtr<HDSA::Vector_MrHyDE<RealT> >(z_opt);
	HDSA::Vector_MrHyDE<RealT> &ez_opt_hdsa = dynamic_cast<HDSA::Vector_MrHyDE<RealT>&>(*z_opt_hdsa);
	ez_opt_hdsa.mrhyde_vec->set(z_opt);	
        return z_opt_hdsa; 
    }

    HDSA::Ptr<HDSA::MultiVector<RealT> > Load_Z_Data(void) const{
    int num_coeff_load = 200;
    //HDSA::Ptr<HDSA::Vector<RealT> > z1 = HDSA::makePtr<HDSA::Vector_MrHyDE_Steady_State<RealT> >(solve_);
    //HDSA::Ptr<HDSA::Vector<RealT> > z2 = HDSA::makePtr<HDSA::Vector_MrHyDE_Steady_State<RealT> >(solve_);

    ROL::Ptr<ROL::Vector<RealT> > z1_rol = solve_->params->getCurrentVector().clone();
    HDSA::Ptr<HDSA::Vector<RealT> > z1 = HDSA::makePtr<HDSA::ROL_Vector<RealT> >(z1_rol);

    ROL::Ptr<ROL::Vector<RealT> > z2_rol = solve_->params->getCurrentVector().clone();
    HDSA::Ptr<HDSA::Vector<RealT> > z2 = HDSA::makePtr<HDSA::ROL_Vector<RealT> >(z2_rol);

    MrHyDE_OptVector &ez1 = dynamic_cast<MrHyDE_OptVector&>(*z1_rol);
    MrHyDE_OptVector &ez2 = dynamic_cast<MrHyDE_OptVector&>(*z2_rol);
    
    //    HDSA::Vector_MrHyDE_Steady_State<RealT> &ez1 = dynamic_cast<HDSA::Vector_MrHyDE_Steady_State<RealT>&>(*z1);
    // HDSA::Vector_MrHyDE_Steady_State<RealT> &ez2 = dynamic_cast<HDSA::Vector_MrHyDE_Steady_State<RealT>&>(*z2);
    std::ifstream in("Z.txt");
    RealT val = 0.0;
    if (in)
      {
        for(int j = 0; j < num_coeff_load; j++)
        {
          in >> val;
	  //          ez1.mrhyde_steady_state_vec[0]->replaceGlobalValue(j,0,val);
	  ez1.getField()[0]->getVector()->replaceGlobalValue(j,0,val);
          in >> val;
          // ez2.mrhyde_steady_state_vec[0]->replaceGlobalValue(j,0,val);
          ez2.getField()[0]->getVector()->replaceGlobalValue(j,0,val);
        }
      }
    else
     {
        std::cout << "Error loading the data from Z.txt" << std::endl;
     }

    HDSA::Ptr<HDSA::MultiVector<RealT> > Z = HDSA::makePtr<HDSA::MultiVector<RealT> >(2,*z1);
    (*Z)[0]->set(*z1);
    (*Z)[1]->set(*z2);
	return Z;
    }

    HDSA::Ptr<HDSA::MultiVector<RealT> > Load_D_Data(void) const{
    int num_coeff_load = 200;
    HDSA::Ptr<HDSA::Vector<RealT> > d1 = HDSA::makePtr<HDSA::Vector_MrHyDE_Steady_State<RealT> >(solve_);
    HDSA::Ptr<HDSA::Vector<RealT> > d2 = HDSA::makePtr<HDSA::Vector_MrHyDE_Steady_State<RealT> >(solve_);
    HDSA::Vector_MrHyDE_Steady_State<RealT> &ed1 = dynamic_cast<HDSA::Vector_MrHyDE_Steady_State<RealT>&>(*d1);
    HDSA::Vector_MrHyDE_Steady_State<RealT> &ed2 = dynamic_cast<HDSA::Vector_MrHyDE_Steady_State<RealT>&>(*d2);
    std::ifstream in("D.txt");
    RealT val = 0.0;
    if (in)
      {   
        for(int j = 0; j < num_coeff_load; j++)
        {   
          in >> val;
          ed1.mrhyde_steady_state_vec[0]->replaceGlobalValue(j,0,val);
          in >> val;
          ed2.mrhyde_steady_state_vec[0]->replaceGlobalValue(j,0,val);
        }   
      }   
    else
     {   
        std::cout << "Error loading the data from D.txt" << std::endl;
     }   

    HDSA::Ptr<HDSA::MultiVector<RealT> > D = HDSA::makePtr<HDSA::MultiVector<RealT> >(2,*d1);
    (*D)[0]->set(*d1);
    (*D)[1]->set(*d2);
	return D;
    }

  };

#endif
