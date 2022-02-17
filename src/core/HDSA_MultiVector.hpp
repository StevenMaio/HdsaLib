#ifndef HDSA_MULTIVECTOR_HPP
#define HDSA_MULTIVECTOR_HPP

#include <fstream>

namespace HDSA
{

template <class RealT>
class MultiVector {

private:
  int num_vecs_;
  std::vector<HDSA::Ptr<HDSA::Vector<RealT> > > vecs_;

public:

  MultiVector(const int num_vecs, const HDSA::Ptr<HDSA::Vector<RealT> > & vec): num_vecs_(num_vecs)
  {
    vecs_.resize(num_vecs);
    for(int k = 0; k < num_vecs; k++)
      {
	vecs_[k] = vec->Clone();
      }
  }

  virtual ~MultiVector()
  { }

  // Access the kth vector
  HDSA::Ptr<HDSA::Vector<RealT> > operator [] (int k) const
  {
    return vecs_[k];
  }

  std::vector<RealT> dot(HDSA::Vector<RealT> & x) const
  {
    std::vector<RealT> ips = std::vector<RealT>(num_vecs_,0.0);
    for(int k = 0; k < num_vecs_; k++)
      {
	ips[k] = vecs_[k]->dot(x);
      }
    return ips;
  }


  int Number_of_Vectors(void) const
  {
    return num_vecs_;
  }

  void Write_to_File(std::string & name) const
  {
      std::ofstream fout;
      fout.open(name);
      for(int i = 0; i < vecs_[0]->dimension(); i++)
	{
	  for(int j = 0; j < num_vecs_; j++)
	    {
	      fout << std::setprecision(16) << (*vecs_[j])(i) << "  ";
	    }
	  fout << "  " << std::endl;
	}
      fout.close();
  }

};

}

#endif
