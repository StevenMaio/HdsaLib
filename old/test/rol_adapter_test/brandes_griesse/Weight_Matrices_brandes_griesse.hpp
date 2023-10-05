#ifndef WEIGHT_MATRICES_BRANDES_GRIESSE_HPP
#define WEIGHT_MATRICES_BRANDES_GRIESSE_HPP

template <class RealT>
class Weight_Matrices_brandes_griesse : public HDSA::Weight_Matrices<RealT> {

private:
  HDSA::Ptr<HDSA::ParameterList> parlist_sensitivity_;
  HDSA::Ptr<HDSA::ParameterList> parlist_;
  std::vector<int> z_vec_zeros_;

public:

  Weight_Matrices_brandes_griesse(const HDSA::Ptr<HDSA::ParameterList> & parlist, const HDSA::Ptr<HDSA::ParameterList> & parlist_sensitivity):
    HDSA::Weight_Matrices<RealT>(parlist_sensitivity), parlist_sensitivity_(parlist_sensitivity), parlist_(parlist)
  { 
    int nx = parlist_->sublist("Geometry").get("NX", 10);
    int ny = parlist_->sublist("Geometry").get("NY", 10);
    int z_dim = (nx+1)*(ny+1);
    std::vector<RealT> z = std::vector<RealT>(z_dim,0.0);
    // find controller zeros
    std::ifstream inputFile("control_read.txt");          
    RealT value;
    int count = 0;
    // read the elements in the file into a vector  
    // test file open   
    if (inputFile) {   
      while ( inputFile >> value ) {
        z[count] = value;
	count += 1;
      }
    }
    else
      {
	std::cout << "Error loading the optimal control solution" << std::endl;
      }    

    z_vec_zeros_ = std::vector<int>(z_dim,0);
    for(int k = 0; k < z_dim; k++)
      {
	if(z[k] != 0.0)
	  {
	    z_vec_zeros_[k] = 1;
	  }
      }

  }

  HDSA::Ptr<HDSA::Weight_Matrices<RealT> > Construct_Weight_Matrices(const HDSA::Ptr<HDSA::Vector<RealT> > & theta, const HDSA::Ptr<const HDSA::Comm<int> > & comm) 
  {
    HDSA::Ptr<HDSA::Weight_Matrices<RealT> > weight_matrices = HDSA::makePtr<Weight_Matrices_brandes_griesse<RealT> >(parlist_,parlist_sensitivity_);
    return weight_matrices;
  }

  void Apply_theta_Weight_Mat(HDSA::Ptr<HDSA::Vector<RealT> > & theta_out, const HDSA::Ptr<HDSA::Vector<RealT> > & theta_in) const 
  {
    theta_out->Replace_Element(0,0.75*(*theta_in)(0));
    theta_out->Replace_Element(1,1.50*(*theta_in)(1));
    for(int k = 2; k < theta_in->dimension(); k++)
    {
      theta_out->Replace_Element(k,(*theta_in)(k));
    }
  }
    
  void Apply_z_Weight_Mat(HDSA::Ptr<HDSA::Vector<RealT> > & z_out, const HDSA::Ptr<HDSA::Vector<RealT> > & z_in) const
  {
    int dim = z_in->dimension();
    int nx = parlist_->sublist("Geometry").get("NX", 10);
    int ny = parlist_->sublist("Geometry").get("NY", 10);
    RealT nx_float = static_cast<RealT>(nx);
    RealT ny_float = static_cast<RealT>(ny);
    RealT delta_x = 2.0/nx_float;
    RealT delta_y = 2.0/ny_float;

    for(int i = 0; i < dim; i++)
      {
	if(z_vec_zeros_[i] == 0.0)
	  {
	    z_in->Replace_Element(i,0.0);
	  }
      }
     
    HDSA::Ptr<Tpetra::MultiVector<> > z_tpetra = HDSA::dynamicPtrCast<ROL_PDEOPT_Tpetra_Vector<RealT> >(z_in)->get_tpetra_vec(); 
    //Teuchos::Array<long long int> Map_myGlobalIds;
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

    for(int i = 0; i < dim; i++)
      {
	if(mat_vec_map->isNodeGlobalElement(i))
	  {
	    val = 0.0;
	   
	    if(z_vec_zeros_[i] != 0.0)
	      {	
		if(i == 0)
		  {
		    val += e_vec_global_1d(0)*ip_corner[0];
		    val += e_vec_global_1d(1)*ip_corner[1];
		    val += e_vec_global_1d(nx+1)*ip_corner[1];
		    val += e_vec_global_1d(nx+2)*ip_corner[2];
		  }
		else if(i == dim-1)
		  {
		    val += e_vec_global_1d(i)*ip_corner[0];
		    val += e_vec_global_1d(i-1)*ip_corner[1];
		    val += e_vec_global_1d(i-nx-1)*ip_corner[1];
		    val += e_vec_global_1d(i-nx-2)*ip_corner[2];
		  }
		else if(i == nx)
		  {
		    val += e_vec_global_1d(i)*ip_corner[0];
		    val += e_vec_global_1d(i-1)*ip_corner[1];
		    val += e_vec_global_1d(i+nx+1)*ip_corner[1];
		    val += e_vec_global_1d(i+nx)*ip_corner[2];
		  }
		else if(i == ny*(nx+1))
		  {
		    val += e_vec_global_1d(i)*ip_corner[0];
		    val += e_vec_global_1d(i+1)*ip_corner[1];
		    val += e_vec_global_1d(i-nx-1)*ip_corner[1];
		    val += e_vec_global_1d(i-nx)*ip_corner[2];	
		  }
		else if(i < nx)
		  {
		    val += e_vec_global_1d(i)*ip_bdry_not_corner[0];
		    val += e_vec_global_1d(i+1)*ip_bdry_not_corner[2];
		    val += e_vec_global_1d(i-1)*ip_bdry_not_corner[2];
		    val += e_vec_global_1d(i+nx+1)*ip_bdry_not_corner[1];
		    val += e_vec_global_1d(i+nx)*ip_bdry_not_corner[3];
		    val += e_vec_global_1d(i+nx+2)*ip_bdry_not_corner[3];
		  }
		else if(i > ny*(nx+1))
		  {
		    val += e_vec_global_1d(i)*ip_bdry_not_corner[0];
		    val += e_vec_global_1d(i+1)*ip_bdry_not_corner[2];
		    val += e_vec_global_1d(i-1)*ip_bdry_not_corner[2];
		    val += e_vec_global_1d(i-nx-1)*ip_bdry_not_corner[1];
		    val += e_vec_global_1d(i-nx)*ip_bdry_not_corner[3];
		    val += e_vec_global_1d(i-nx-2)*ip_bdry_not_corner[3];
		  }
		else if(i%(nx+1) == 0)
		  {
		    val += e_vec_global_1d(i)*ip_bdry_not_corner[0];
		    val += e_vec_global_1d(i+1)*ip_bdry_not_corner[1];
		    val += e_vec_global_1d(i-nx-1)*ip_bdry_not_corner[2];
		    val += e_vec_global_1d(i+nx+1)*ip_bdry_not_corner[2];
		    val += e_vec_global_1d(i-nx)*ip_bdry_not_corner[3];
		    val += e_vec_global_1d(i+nx+2)*ip_bdry_not_corner[3];
		  }
		else if((i+1)%(nx+1) == 0)
		  {
		    val += e_vec_global_1d(i)*ip_bdry_not_corner[0];
		    val += e_vec_global_1d(i-1)*ip_bdry_not_corner[1];
		    val += e_vec_global_1d(i-nx-1)*ip_bdry_not_corner[2];
		    val += e_vec_global_1d(i+nx+1)*ip_bdry_not_corner[2];
		    val += e_vec_global_1d(i-nx-2)*ip_bdry_not_corner[3];
		    val += e_vec_global_1d(i+nx)*ip_bdry_not_corner[3];
		  }
		else
		  {
		    val += e_vec_global_1d(i)*ip_interior[0];
		    val += e_vec_global_1d(i+1)*ip_interior[1];
		    val += e_vec_global_1d(i-1)*ip_interior[1];
		    val += e_vec_global_1d(i+nx+1)*ip_interior[1];
		    val += e_vec_global_1d(i-nx-1)*ip_interior[1];
		    val += e_vec_global_1d(i+nx+2)*ip_interior[2];
		    val += e_vec_global_1d(i+nx)*ip_interior[2];
		    val += e_vec_global_1d(i-nx)*ip_interior[2];
		    val += e_vec_global_1d(i-nx-2)*ip_interior[2];
		  }	
	      
		z_out->Replace_Element(i,val);   
	      }
	    
	  }
      }
  }
    
};


#endif
