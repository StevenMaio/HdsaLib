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
  static ROL::Elementwise::NormalRandom<RealT> nr;


  ROL_Vector(ROL::Vector<RealT> & rol_vec_in) 
  { 
    rol_vec = rol_vec_in.clone();
    rol_vec->set(rol_vec_in);
  }

  ROL_Vector(const ROL::Ptr<ROL::Vector<RealT>>& rol_vec_in) : rol_vec(rol_vec_in) {};

  virtual ~ROL_Vector()
  { }

  // Clone the vector
  HDSA::Ptr<HDSA::Vector<RealT> > clone() const override
  {
    ROL::Ptr<ROL::Vector<RealT> > rol_vec_clone = rol_vec->clone();
    rol_vec_clone->zero(); //ROL clone() vector is not initialized
    return Teuchos::rcp(new HDSA::ROL_Vector<RealT>(rol_vec_clone)); 
  }

  // compute the dot product of this and x
  RealT dot( const HDSA::Vector<RealT> &x ) const override
  {
    const HDSA::ROL_Vector<RealT> &ex = dynamic_cast<const HDSA::ROL_Vector<RealT>&>(x);
    RealT val = ex.rol_vec->dot(*rol_vec);
    return val;
  }

  // add alpha*x to this
  void axpy( const RealT alpha, const HDSA::Vector<RealT> &x ) override
  {
    const HDSA::ROL_Vector<RealT> &ex = dynamic_cast<const HDSA::ROL_Vector<RealT>&>(x);
    rol_vec->axpy(alpha,*ex.rol_vec);
  }

  // return vector dimension
  int dimension() const override
  {
    return rol_vec->dimension();
  }
  
  // set this=val elementwise
  void setScalar( const RealT val ) override
  {
    rol_vec->setScalar(val);
  }

  // set this=val elementwise
  void set(ROL::Vector<RealT> & rol_vec_in)
  {
    rol_vec->set(rol_vec_in);
  }
  
  void randomize_standard_normal( ) override
  {
    rol_vec->applyUnary(nr);
  }

  void Write_to_File(const std::string & name) const override
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

template<>
ROL::Elementwise::NormalRandom<double> ROL_Vector<double>::nr = ROL::Elementwise::NormalRandom<double>(0.0,1.0);

}

#endif
