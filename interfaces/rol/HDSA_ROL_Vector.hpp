#ifndef HDSA_ROL_VECTOR_HPP
#define HDSA_ROL_VECTOR_HPP

#include "ROL_Vector.hpp"
#include "ROL_StdVector.hpp"

namespace HDSA
{

template <class RealT>
class ROL_Vector : public Vector<RealT> {

public:
  ROL::Ptr<ROL::Vector<RealT> > rol_vec;


  ROL_Vector(ROL::Vector<RealT> & rol_vec_in) 
  { 
    rol_vec = rol_vec_in.clone();
  }

  virtual ~ROL_Vector()
  { }

  // Clone the vector
  HDSA::Ptr<HDSA::Vector<RealT> > clone() const
  {
    HDSA::Ptr<HDSA::Vector<RealT> > hdsa_vector = HDSA::makePtr<HDSA::ROL_Vector<RealT> >(*rol_vec);
    return hdsa_vector;
  }

  // compute the dot product of this and x
  RealT dot( const HDSA::Vector<RealT> &x ) const
  {
    const HDSA::ROL_Vector<RealT> &ex = dynamic_cast<const HDSA::ROL_Vector<RealT>&>(x);
    RealT val = ex.rol_vec->dot(*rol_vec);
    return val;
  }

  // add alpha*x to this
  void axpy( const RealT alpha, const HDSA::Vector<RealT> &x )
  {
    const HDSA::ROL_Vector<RealT> &ex = dynamic_cast<const HDSA::ROL_Vector<RealT>&>(x);
    rol_vec->axpy(alpha,*ex.rol_vec);
  }

  // return vector dimension
  int dimension() const
  {
    return rol_vec->dimension();
  }
  
  // set this=val elementwise
  void setScalar( const RealT val )
  {
    rol_vec->setScalar(val);
  }
  
  void randomize_standard_normal( )
  {
    ROL::Elementwise::NormalRandom<RealT> nr;
    rol_vec->applyUnary(nr);
  }

  void Write_to_File(std::string & name) const
  {
    try{
      ROL::Ptr<std::vector<RealT> > vec = dynamic_cast<ROL::StdVector<RealT>&>(*rol_vec).getVector();
      std::ofstream fout;
      fout.open(name);
      for(int i = 0; i < rol_vec->dimension(); i++)
	{
	  fout << std::setprecision(16) << (*vec)[i] << "  ";
	}
      fout.close();
    }
    catch (...)
      {
	std::cout << "Write_to_File is currently not supported for this vector type" << std::endl;
      }
  }
 
};

}

#endif
