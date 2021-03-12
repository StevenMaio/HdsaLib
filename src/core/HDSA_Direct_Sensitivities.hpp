#ifndef HDSA_DIRECT_SENSITIVITIES_HPP
#define HDSA_DIRECT_SENSITIVITIES_HPP

namespace HDSA
{
  
  template <class RealT>
  class Direct_Sensitivity_Computation
  {
  private:
    HDSA::Ptr<HDSA::Vector<RealT> > theta_;
    HDSA::Ptr<HDSA::ParameterList> parlist_sensitivity_;
    const HDSA::Ptr<const HDSA::Comm<int> >  comm_;
    HDSA::Ptr<Opt_Problem_Objects<RealT> > OP_Objects_Factory_;
    HDSA::Ptr<Weight_Matrices<RealT> > weight_matrices_factory_;
    int sample_index_;

  public:

    Direct_Sensitivity_Computation(const HDSA::Ptr<HDSA::Vector<RealT> > & theta, const HDSA::Ptr<HDSA::ParameterList> & parlist_sensitivity,
				   const HDSA::Ptr<const HDSA::Comm<int> > & comm, const HDSA::Ptr<Opt_Problem_Objects<RealT> > & OP_Objects_Factory,
				   const HDSA::Ptr<Weight_Matrices<RealT> > & weight_matrices_factory, const int & sample_index):
      theta_(theta), parlist_sensitivity_(parlist_sensitivity), comm_(comm), OP_Objects_Factory_(OP_Objects_Factory), weight_matrices_factory_(weight_matrices_factory), sample_index_(sample_index)
    { }

    void Compute(void)
    {
      // This function executes direct sensitivity computation
      std::clock_t timer_direct = std::clock();
      
      int theta_dim = theta_->dimension();
      
      HDSA::Ptr<HDSA::Processor_Distribution<RealT> > proc_dist = HDSA::makePtr<HDSA::Processor_Distribution<RealT> >(comm_,theta_dim);
      HDSA::Ptr<HDSA::Comm<int> > subcomm = comm_->createSubcommunicator(proc_dist->Get_Comm_Split_Ranks());
      
      HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > OP_Objects_subcomm = OP_Objects_Factory_->Construct_Opt_Problem_Objects(theta_,subcomm);
      OP_Objects_subcomm->Load_Optimal_Solution();
      HDSA::Ptr<HDSA::Weight_Matrices<RealT> > weight_matrices_subcomm = weight_matrices_factory_->Construct_Weight_Matrices(theta_,subcomm); 
      HDSA::Ptr<HDSA::Nominal_Data<RealT> > Nom_subcomm = HDSA::makePtr<HDSA::Nominal_Data<RealT> >(parlist_sensitivity_,OP_Objects_subcomm);
      
      bool reduced_space_sen = parlist_sensitivity_->sublist("Formulation").get("Reduced Space Sensitivities",true);
      HDSA::Ptr<HDSA::Sensitivity_Operators<RealT> > Sen_Op_subcomm;
      if(reduced_space_sen)
	{
	  Sen_Op_subcomm = HDSA::makePtr<HDSA::Sensitivity_Operators_RS<RealT> >(OP_Objects_subcomm,Nom_subcomm,subcomm, proc_dist->Get_Comm_Split_Ranks());
	}
      else
	{
	  Sen_Op_subcomm = HDSA::makePtr<HDSA::Sensitivity_Operators_FS<RealT> >(OP_Objects_subcomm,Nom_subcomm,subcomm, proc_dist->Get_Comm_Split_Ranks());
	}
      
      HDSA::Ptr<HDSA::Vector<RealT> > theta_1 = OP_Objects_subcomm->theta->Clone();
      HDSA::Ptr<HDSA::Vector<RealT> > theta_2 = OP_Objects_subcomm->theta->Clone();
      HDSA::Ptr<HDSA::Vector<RealT> > z_1 = OP_Objects_subcomm->z->Clone();
      HDSA::Ptr<HDSA::Vector<RealT> > z_2 = OP_Objects_subcomm->z->Clone();
      
      bool output_vectors = parlist_sensitivity_->sublist("Formulation").get("Output Sensitivity Operator Columns", false);
      std::vector<RealT> sensitivity_indices = std::vector<RealT>(theta_dim,0.0);
      
      for(int k = 0; k < theta_dim; k++)
	{
	  if(proc_dist->Does_Processor_Own_Vector(k))
	    {
	      theta_1->basis(k);
	      
	      Sen_Op_subcomm->Apply_Sensitivity_Operator(z_1,theta_1);
	      
	      if(output_vectors)
		{
		  // Write solutions to text files
		  std::string name;
		  std::ofstream fout;
		  name = "Sensitivity_Operator_Column_" + std::to_string(k+1) + ".txt";
		  fout.open(name);
		  for(int k = 0; k < z_1->dimension(); k++)
		    {
		      fout << (*z_1)(k) << std::setw(20);
		    }
		  fout.close();
		}
	      
	      weight_matrices_subcomm->Apply_z_Weight_Mat(z_2,z_1);
	      weight_matrices_subcomm->Apply_theta_Weight_Mat(theta_2,theta_1);
	      sensitivity_indices[k] = std::sqrt(z_2->dot(*z_1))/std::sqrt(theta_1->dot(*theta_2));
	      
	      theta_2->zero();
	      z_1->zero();
	      z_2->zero();
	      if(subcomm->getRank() == 0)
		{
		  std::cout << "Computed sensitivity index " << k+1 << " out of " << theta_dim << " with local sensitivity sample number " << sample_index_ << std::endl;
		}
	    }
	}
      
      comm_->barrier();
      for(int k = 0; k < theta_dim; k++)
	{
	  char *buff = (char*)&sensitivity_indices[k];
	  comm_->broadcast(proc_dist->Get_Procs_Loop_Distribution()[k][0],8,buff);
	}
      comm_->barrier();
      
      // Write solutions to text files
      std::string name;
      std::ofstream fout;
      name = "Sensitivity_Indices.txt";
      fout.open(name);
      for(int k = 0; k < theta_dim; k++)
	{
	  fout << sensitivity_indices[k] << std::setw(20);
	}
      fout.close();
      
      RealT Time = static_cast<RealT>(std::clock()-timer_direct)/static_cast<RealT>(CLOCKS_PER_SEC);
      if(comm_->getRank() == 0)
	{
	  std::cout << " " << std::endl;
	  std::cout << "Total time for sensitivity computation: " << Time << " seconds with local sensitivity sample number " << sample_index_ << std::endl;
	}
    }
    
  };

}

#endif
