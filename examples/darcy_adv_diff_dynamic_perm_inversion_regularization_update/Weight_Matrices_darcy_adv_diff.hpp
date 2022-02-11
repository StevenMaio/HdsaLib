#ifndef WEIGHT_MATRICES_DARCY_ADV_DIFF_HPP
#define WEIGHT_MATRICES_DARCY_ADV_DIFF_HPP

template <class RealT>
class Weight_Matrices_darcy_adv_diff : public HDSA::Weight_Matrices<RealT> {

private:
  HDSA::Ptr<HDSA::ParameterList> parlist_;
  HDSA::Ptr<HDSA::ParameterList> parlist_sensitivity_;

public:

  Weight_Matrices_darcy_adv_diff(const HDSA::Ptr<HDSA::ParameterList> & parlist, const HDSA::Ptr<HDSA::ParameterList> & parlist_sensitivity):
    HDSA::Weight_Matrices<RealT>(parlist_sensitivity), parlist_(parlist), parlist_sensitivity_(parlist_sensitivity)
  { }

  HDSA::Ptr<HDSA::Weight_Matrices<RealT> > Construct_Weight_Matrices(const HDSA::Ptr<HDSA::Vector<RealT> > & theta, const HDSA::Ptr<const HDSA::Comm<int> > & comm) 
  {
    HDSA::Ptr<HDSA::Weight_Matrices<RealT> > weight_matrices = HDSA::makePtr<Weight_Matrices_darcy_adv_diff<RealT> >(parlist_,parlist_sensitivity_);
    return weight_matrices;
  }

  void Apply_theta_Weight_Mat(HDSA::Ptr<HDSA::Vector<RealT> > & theta_out, const HDSA::Ptr<HDSA::Vector<RealT> > & theta_in) const 
  {
    int L = parlist_->sublist("Problem").get("Number of Uncertainty Basis Function", 10);

    // Entries [0] to [9*16-1] is the source
    for(int i = 0; i < 9*16; i++)
      {
	RealT val = 0.0;
	int i_box_num = i/9;
	int i_node_index = i%9;

	for(int j = 0; j < 9*16; j++)
	  {
	    int j_box_num = j/9;
	    int j_node_index = j%9;
	    
	    RealT x_ip = 0.0;
	    RealT y_ip = 0.0;
	    if(i_box_num == j_box_num)
	      {
		int x1 = i_node_index%3;
		int x2 = j_node_index%3;
		int y1 = i_node_index/3;
		int y2 = j_node_index/3;
		x_ip = Evaluate_Box_Basis_Function_Inner_Products(x1,x2);
		y_ip = Evaluate_Box_Basis_Function_Inner_Products(y1,y2);
	      }
	    
	    val += (*theta_in)(j)*x_ip*y_ip;
	  }
	
	theta_out->Replace_Element(i,val);
      }

    // Entry [9*16] is the diffusion coefficient
    theta_out->Replace_Element(9*16,(*theta_in)(9*16));
    
    // Entry [9*16+1] to [9*16+(L+1)] is the left boundary condition
    int offset = 9*16+1;
    for(int i = 0; i < (L+1); i++)
      {
	RealT val = 0.0;
	int xi = Evaluate_x_coordinate(i,L);

	for(int j = 0; j < (L+1); j++)
	  {
	    int xj = Evaluate_x_coordinate(j,L);
	    RealT x_ip = Evaluate_Basis_Function_Inner_Products(xi,xj,L);
	    val += (*theta_in)(offset+j)*x_ip;
	  }
	theta_out->Replace_Element(offset+i,val);
      }

    // Entry [9*16+(L+1)+1] to [9*16+(L+1)+(L+1)] is the right boundary condition
    offset = 9*16+1+(L+1);
    for(int i = 0; i < (L+1); i++)
      {
	RealT val = 0.0;
	int xi = Evaluate_x_coordinate(i,L);

	for(int j = 0; j < (L+1); j++)
	  {
	    int xj = Evaluate_x_coordinate(j,L);
	    RealT x_ip = Evaluate_Basis_Function_Inner_Products(xi,xj,L);
	    val += (*theta_in)(offset+j)*x_ip;
	  }
	theta_out->Replace_Element(offset+i,val);
      }
  }
    
  void Apply_z_Weight_Mat(HDSA::Ptr<HDSA::Vector<RealT> > & z_out, const HDSA::Ptr<HDSA::Vector<RealT> > & z_in) const
  {
    int dim = z_in->dimension();
    int nx = parlist_->sublist("Geometry").get("NX", 10);
    int ny = parlist_->sublist("Geometry").get("NY", 10);
    RealT nx_float = static_cast<RealT>(nx);
    RealT ny_float = static_cast<RealT>(ny);
    RealT delta_x = 1.0/nx_float;
    RealT delta_y = 1.0/ny_float;

    HDSA::Ptr<Tpetra::MultiVector<> > z_tpetra = HDSA::dynamicPtrCast<ROL_PDEOPT_Tpetra_Vector<RealT> >(z_in)->get_tpetra_vec(); 
    Teuchos::Array<int> Map_myGlobalIds;
    for(int i = 0; i < dim; i++)
      {
	Map_myGlobalIds.push_back(i);
      }  
    HDSA::Ptr<const Tpetra::Map<> > Map = HDSA::makePtr<Tpetra::Map<> >(dim,Map_myGlobalIds,0,z_tpetra->getMap()->getComm());
    HDSA::Ptr<Tpetra::MultiVector<> > vec_global = HDSA::makePtr<Tpetra::MultiVector<> >(Map, 1, true);
    Tpetra::Import<> Importer = Tpetra::Import<>(z_tpetra->getMap(),vec_global->getMap());
    vec_global->doImport(*z_tpetra,Importer,Tpetra::REPLACE);
    (*vec_global).sync<Kokkos::HostSpace> ();
    auto e_vec_global_2d = (*vec_global).getLocalView<Kokkos::HostSpace> ();
    auto e_vec_global_1d = Kokkos::subview (e_vec_global_2d, Kokkos::ALL (), 0);

    // Store values of inner products which populate the mass matrix
    std::vector<RealT> ip_corner;
    ip_corner.resize(3); 
    // First entry is the inner product of basis function with itself, second is with its neighbor, third is with its adjacent
    ip_corner[0] = (1.0/9.0)*delta_x*delta_y;
    ip_corner[1] = (1.0/18.0)*delta_x*delta_y;
    ip_corner[2] = (1.0/36.0)*delta_x*delta_y;
    std::vector<RealT> ip_bdry_not_corner;
    ip_bdry_not_corner.resize(4);
    // 1st entry is the inner product of basis function with itself, 2nd is with its interior neighbor, 3rd is with its side neighbor, 4th is with its adjacent
    ip_bdry_not_corner[0] = (2.0/9.0)*delta_x*delta_y;
    ip_bdry_not_corner[1] = (1.0/9.0)*delta_x*delta_y;
    ip_bdry_not_corner[2] = (1.0/18.0)*delta_x*delta_y;
    ip_bdry_not_corner[3] = (1.0/36.0)*delta_x*delta_y;
    std::vector<RealT> ip_interior;
    ip_interior.resize(3); 
    // First entry is the inner product of basis function with itself, second is with its neighbor, third is with its adjacent
    ip_interior[0] = (4.0/9.0)*delta_x*delta_y;
    ip_interior[1] = (1.0/9.0)*delta_x*delta_y;
    ip_interior[2] = (1.0/36.0)*delta_x*delta_y;
    
    RealT val = 0.0;
    z_out->zero();
    auto mat_vec_map = HDSA::dynamicPtrCast<ROL_PDEOPT_Tpetra_Vector<RealT> >(z_out)->get_tpetra_vec()->getMap();

    for(int i = 0; i < dim; i+=2)
      {
	if(mat_vec_map->isNodeGlobalElement(i))
	  {
	    val = 0.0;
	    int j = i/2;
	    
	    if(j == 0)
	      {
		val += e_vec_global_1d(2*0)*ip_corner[0];
		val += e_vec_global_1d(2*1)*ip_corner[1];
		val += e_vec_global_1d(2*(nx+1))*ip_corner[1];
		val += e_vec_global_1d(2*(nx+2))*ip_corner[2];
	      }
	    else if(j == dim/2-1)
	      {
		val += e_vec_global_1d(2*j)*ip_corner[0];
		val += e_vec_global_1d(2*(j-1))*ip_corner[1];
		val += e_vec_global_1d(2*(j-nx-1))*ip_corner[1];
		val += e_vec_global_1d(2*(j-nx-2))*ip_corner[2];
	      }
	    else if(j == nx)
	      {
		val += e_vec_global_1d(2*j)*ip_corner[0];
		val += e_vec_global_1d(2*(j-1))*ip_corner[1];
		val += e_vec_global_1d(2*(j+nx+1))*ip_corner[1];
		val += e_vec_global_1d(2*(j+nx))*ip_corner[2];
	      }
	    else if(j == ny*(nx+1))
	      {
		val += e_vec_global_1d(2*j)*ip_corner[0];
		val += e_vec_global_1d(2*(j+1))*ip_corner[1];
		val += e_vec_global_1d(2*(j-nx-1))*ip_corner[1];
		val += e_vec_global_1d(2*(j-nx))*ip_corner[2];	
	      }
	    else if(j < nx)
	      {
		val += e_vec_global_1d(2*j)*ip_bdry_not_corner[0];
		val += e_vec_global_1d(2*(j+1))*ip_bdry_not_corner[2];
		val += e_vec_global_1d(2*(j-1))*ip_bdry_not_corner[2];
		val += e_vec_global_1d(2*(j+nx+1))*ip_bdry_not_corner[1];
		val += e_vec_global_1d(2*(j+nx))*ip_bdry_not_corner[3];
		val += e_vec_global_1d(2*(j+nx+2))*ip_bdry_not_corner[3];
	      }
	    else if(j > ny*(nx+1))
	      {
		val += e_vec_global_1d(2*j)*ip_bdry_not_corner[0];
		val += e_vec_global_1d(2*(j+1))*ip_bdry_not_corner[2];
		val += e_vec_global_1d(2*(j-1))*ip_bdry_not_corner[2];
		val += e_vec_global_1d(2*(j-nx-1))*ip_bdry_not_corner[1];
		val += e_vec_global_1d(2*(j-nx))*ip_bdry_not_corner[3];
		val += e_vec_global_1d(2*(j-nx-2))*ip_bdry_not_corner[3];
	      }
	    else if(j%(nx+1) == 0)
	      {
		val += e_vec_global_1d(2*j)*ip_bdry_not_corner[0];
		val += e_vec_global_1d(2*(j+1))*ip_bdry_not_corner[1];
		val += e_vec_global_1d(2*(j-nx-1))*ip_bdry_not_corner[2];
		val += e_vec_global_1d(2*(j+nx+1))*ip_bdry_not_corner[2];
		val += e_vec_global_1d(2*(j-nx))*ip_bdry_not_corner[3];
		val += e_vec_global_1d(2*(j+nx+2))*ip_bdry_not_corner[3];
	      }
	    else if((j+1)%(nx+1) == 0)
	      {
		val += e_vec_global_1d(2*j)*ip_bdry_not_corner[0];
		val += e_vec_global_1d(2*(j-1))*ip_bdry_not_corner[1];
		val += e_vec_global_1d(2*(j-nx-1))*ip_bdry_not_corner[2];
		val += e_vec_global_1d(2*(j+nx+1))*ip_bdry_not_corner[2];
		val += e_vec_global_1d(2*(j-nx-2))*ip_bdry_not_corner[3];
		val += e_vec_global_1d(2*(j+nx))*ip_bdry_not_corner[3];
	      }
	    else
	      {
		val += e_vec_global_1d(2*j)*ip_interior[0];
		val += e_vec_global_1d(2*(j+1))*ip_interior[1];
		val += e_vec_global_1d(2*(j-1))*ip_interior[1];
		val += e_vec_global_1d(2*(j+nx+1))*ip_interior[1];
		val += e_vec_global_1d(2*(j-nx-1))*ip_interior[1];
		val += e_vec_global_1d(2*(j+nx+2))*ip_interior[2];
		val += e_vec_global_1d(2*(j+nx))*ip_interior[2];
		val += e_vec_global_1d(2*(j-nx))*ip_interior[2];
		val += e_vec_global_1d(2*(j-nx-2))*ip_interior[2];
	      }
	  
	  }
	z_out->Replace_Element(i,val);   
      }

  }
   
private:

  RealT Evaluate_Box_Basis_Function_Inner_Products(int i, int j) const
  {
    RealT val = 0.0;
    
    if(i == j)
      {
	val = 0.1/3.0;
      }
    else if(std::abs(i-j) == 1)
      {
	val = 0.05/6.0;
      }
    
    return val;
  }

  RealT Evaluate_Basis_Function_Inner_Products(int i, int j, int L) const
  {
    RealT val = 0.0;
    RealT Ld = static_cast<RealT>(L);
    
    if(i == j)
      {
	if(i == 0)
	  {
	    val = 1.0/(3.0*Ld);
	  }
	else if(i == L)
	  {
	    val = 1.0/(3.0*Ld);
	  }
	else
	  {
	    val = 2.0/(3.0*Ld);
	  }
      }
    else if(i+1 == j || i-1 == j)
      {
	val = 1.0/(6.0*Ld);
      }
    
    return val;
  }
  
  int Evaluate_x_coordinate(int i, int L) const
  {
    int xi = i%(L+1);
    return xi;
  }
  
  int Evaluate_y_coordinate(int i, int L) const
  {
    int yi = std::ceil(static_cast<RealT>(i+1)/static_cast<RealT>(L+1)) - 1;
    return yi;
  }
  
};


#endif
