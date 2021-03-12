#ifndef HDSA_SAMPLE_LOCAL_SENSITIVITIES_HPP
#define HDSA_SAMPLE_LOCAL_SENSITIVITIES_HPP

namespace HDSA
{

  template<typename RealT>
  void Sample_Local_Sensitivities(const HDSA::Ptr<const HDSA::Comm<int> > & comm, const HDSA::Ptr<HDSA::ParameterList> & parlist_sensitivity, 
				  const HDSA::Ptr<Opt_Problem_Objects<RealT> > & OP_Objects_Factory, const HDSA::Ptr<HDSA::Weight_Matrices<RealT> > & weight_matrices,  
				  const HDSA::Ptr<HDSA::Parameter_Sampler<RealT> > & sampler)
  {
    // Divide processors to compute local sensitivities in parallel
    int num_samp = parlist_sensitivity->sublist("Formulation").get("Number of Parameter Samples",1);
    int num_sub_comm = 0;
    if(num_samp <= comm->getSize())
      {
	num_sub_comm = num_samp;
      }
    else
      {
	num_sub_comm = comm->getSize();
	num_samp = num_sub_comm;
	std::cout << "Number of processors fewer than the number of requested samples, computing fewer samples" << std::endl;
      }
    
    int procs_per_sub_comm = floor(comm->getSize()/num_sub_comm);
    int remainder = comm->getSize() - procs_per_sub_comm*num_sub_comm;
    std::vector<std::vector<int> > proc_distribution;
    proc_distribution.resize(num_sub_comm);
    int count = 0;
    for(int k = 0; k < remainder; k++)
      {
	proc_distribution[k].resize(procs_per_sub_comm+1);
	for(int i = 0; i < procs_per_sub_comm+1; i++)
	  {
	    proc_distribution[k][i] = count;
	    count = count+1;
	  }
      }
    for(int k = remainder; k < num_sub_comm; k++)
      {
	proc_distribution[k].resize(procs_per_sub_comm);
	for(int i = 0; i < procs_per_sub_comm; i++)
	  {
	    proc_distribution[k][i] = count;
	    count = count+1;
	  }
      }
    
    std::vector<int> Comm_Split_Ranks;
    int myRank = comm->getRank();
    int myRank_position = 0;
    for(int k = 0; k < num_sub_comm; k++)
      {
	if( myRank >= proc_distribution[k].front() && myRank <= proc_distribution[k].back() )
	  {
	    myRank_position = k;
	    for(unsigned i = 0; i < proc_distribution[k].size(); i++)
	      {
		Comm_Split_Ranks.push_back(proc_distribution[k][i]);
	      }
	  }
      }
    
    if(myRank == 0)
      {
	std::cout << "Computing " << num_samp << " local sensitivity(ies)" << std::endl;
      }

    // Draw theta samples where the local sensitivities will be computed
    std::vector<HDSA::Ptr<HDSA::Vector<RealT> > > theta_samps = sampler->Draw_Samples(num_samp); 
    // Split into subcommunicators and compute local sensitivities in parallel
    HDSA::Ptr<HDSA::Comm<int> > Comm_split = comm->createSubcommunicator(Comm_Split_Ranks);
    HDSA::Ptr<HDSA::Sensitivity_Computation<RealT> > Sen_Comp = HDSA::makePtr<HDSA::Sensitivity_Computation<RealT> >(theta_samps[myRank_position], parlist_sensitivity, Comm_split, 
														     OP_Objects_Factory, weight_matrices, Comm_Split_Ranks, myRank_position); 
  }

}

#endif
