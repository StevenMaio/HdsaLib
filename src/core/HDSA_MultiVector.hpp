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

  HDSA::Ptr<HDSA::Vector<RealT> > MatVec(HDSA::Vector<RealT> & x) const
  {
    HDSA::Ptr<HDSA::Vector<RealT> > Ax = HDSA::makePtr<Std_Vector<RealT> >(num_vecs_);
    for(int k = 0; k < num_vecs_; k++)
      {
        Ax->Replace_Element(k,vecs_[k]->dot(x));
      }
    return Ax;
  }

  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > MatMat(HDSA::MultiVector<RealT> & x) const
  {
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > C = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(x.Number_of_Vectors(),num_vecs_);
    for(int i = 0; i < x.Number_of_Vectors(); i++)
      {
	for(int j = 0; j < num_vecs_; j++)
	{
	  C->Replace_Element(i,j,x[i]->dot(*vecs_[j]));
	}
      }
    return C;
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
