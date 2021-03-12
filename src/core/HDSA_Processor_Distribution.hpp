#ifndef HDSA_PROCESSOR_DISTRIBUTION_HPP
#define HDSA_PROCESSOR_DISTRIBUTION_HPP

namespace HDSA
{

  template <class RealT>
  class Processor_Distribution
  {
  private:
    const HDSA::Ptr<const HDSA::Comm<int> > comm_;
    int num_vecs_;
    int numSubComm_; 
    std::vector<int> Comm_Split_Ranks_; // Which processors will work with processor myRank
    std::vector<std::vector<int> > Procs_Loop_Distribution_; // Vector which identifies which processors will execute which matrix vector products in the parallel loop
    
  public:
    
    Processor_Distribution(const HDSA::Ptr<const HDSA::Comm<int> > & comm, int & num_vecs): comm_(comm), num_vecs_(num_vecs)
    {
      // Distribute matrix vector products to subcommunicators 
      int myRank = comm_->getRank();
      int numProcs = comm_->getSize();
      int myRank_Position = 0;
      std::vector<int> Subcomm_Loop_Distribution; // Vector which identifies which subcommunicators will execute which matrix vector products in the parallel loop
      Procs_Loop_Distribution_.resize(num_vecs_);
      Subcomm_Loop_Distribution.resize(num_vecs_);
      
      // Split into subcommunicators here

      if(numProcs <= num_vecs_) // Case when the number of processors is less than the number of iterations in the parallel loop 
	{
	  int num_vecs_per_Proc = floor(num_vecs_/numProcs); // Number of loop iterations per processor
	  int remainder = num_vecs_-num_vecs_per_Proc*numProcs; // How many processors get num_vecs_per_Proc+1 loop iterations
	  if(remainder == 0)
	    {
	      for(int k = 0; k < num_vecs_; k++)
		{
		  Procs_Loop_Distribution_[k].push_back(floor(k/num_vecs_per_Proc));
		  Subcomm_Loop_Distribution[k] = Procs_Loop_Distribution_[k][0];
		}
	    }
	  else
	    {
	      for(int k = 0; k < remainder*(num_vecs_per_Proc+1); k++)
		{
		  Procs_Loop_Distribution_[k].push_back(floor(k/(num_vecs_per_Proc+1)));
		  Subcomm_Loop_Distribution[k] = Procs_Loop_Distribution_[k][0];
		}
	      for(int k = remainder*(num_vecs_per_Proc+1); k < num_vecs_; k++)
		{
		  Procs_Loop_Distribution_[k].push_back(floor((k-remainder)/num_vecs_per_Proc));
		  Subcomm_Loop_Distribution[k] = Procs_Loop_Distribution_[k][0];
		}
	    }
	  Comm_Split_Ranks_.push_back(myRank);
	  numSubComm_ = numProcs;
	}
      else // Case when we have more processors than the number of iterations in the parallel loop
	{
	  int Proc_per_num_vecs = floor(numProcs/num_vecs_); // Number of processors per parallel loop iteration
	  int remainder = numProcs-Proc_per_num_vecs*num_vecs_; // How many loop iterations get Proc_per_num_vecs+a processors
	  if(remainder == 0)
	    {
	      for(int k = 0; k < num_vecs_; k++)
		{
		  Subcomm_Loop_Distribution[k] = k;
		  for(int i = 0; i < Proc_per_num_vecs; i++)
		    {
		      Procs_Loop_Distribution_[k].push_back(k*Proc_per_num_vecs+i);
		      if(k*Proc_per_num_vecs+i == myRank)
			{
			  myRank_Position = k;
			}
		    }
		}
	    }
	  else
	    {
	      for(int k = 0; k < remainder; k++)
		{
		  Subcomm_Loop_Distribution[k] = k;
		  for(int i = 0; i < Proc_per_num_vecs+1; i++)
		    {
		      Procs_Loop_Distribution_[k].push_back(k*(Proc_per_num_vecs+1)+i);
		      if(k*(Proc_per_num_vecs+1)+i == myRank)
			{
			  myRank_Position = k;
			}
		    }
		}
	      for(int k = remainder; k < num_vecs_; k++)
		{
		  Subcomm_Loop_Distribution[k] = k;
		  for(int i = 0; i < Proc_per_num_vecs; i++)
		    {
		      Procs_Loop_Distribution_[k].push_back(remainder*(Proc_per_num_vecs+1)+(k-remainder)*Proc_per_num_vecs+i);
		      if(remainder*(Proc_per_num_vecs+1)+(k-remainder)*Proc_per_num_vecs+i == myRank)
			{
			  myRank_Position = k;
			}
		    }
		}
	    }
	  
	  numSubComm_ = num_vecs_;
	  for(unsigned int i = 0; i < Procs_Loop_Distribution_[myRank_Position].size(); i++)
	    {
	      Comm_Split_Ranks_.push_back(Procs_Loop_Distribution_[myRank_Position][i]);
	    }
	}
    }

    std::vector<int> Get_Comm_Split_Ranks(void)
    {
      return Comm_Split_Ranks_;
    }

    std::vector<std::vector<int> > Get_Procs_Loop_Distribution(void)
    {
      return Procs_Loop_Distribution_;
    }

    int Get_numSubComm(void)
    {
      return numSubComm_;
    }

    bool Does_Processor_Own_Vector(int & k)
    {
      bool own = false;
      if(Procs_Loop_Distribution_[k][0] == Comm_Split_Ranks_[0])
	{
	  own = true;
	}
      return own;
    }

    void Broadcast_Matrix(HDSA::Ptr<HDSA::Dense_Matrix<RealT> > & A)
    {
      int dim = A->numRows();
      for(int k = 0; k < num_vecs_; k++)
	{
	  char *buff = (char*)A->Get_Element_Ptr(0,k);
	  comm_->broadcast(Get_Procs_Loop_Distribution()[k][0],8*dim,buff);
	}
      comm_->barrier();
    }

  };

}

#endif
