#include "../../../../Trilinos/packages/rol/example/PDE-OPT/TOOLS/meshmanager.hpp"

template <class Real>
class MeshManager_shallow_ice : public MeshManager_Rectangle<Real>
{

private:
  int nx_;
  int ny_;
  ROL::Ptr<std::vector<std::vector<std::vector<int>>>> meshSidesets_;

public:
  MeshManager_shallow_ice(Teuchos::ParameterList &parlist) : MeshManager_Rectangle<Real>(parlist)
  {
    nx_ = parlist.sublist("Geometry").get("NX", 3);
    ny_ = parlist.sublist("Geometry").get("NY", 3);
    computeSidesets();
  }

  void computeSidesets()
  {

    int numSidesets = 2;
    meshSidesets_ = ROL::makePtr<std::vector<std::vector<std::vector<int>>>>(numSidesets);

    // Dirichlet
    (*meshSidesets_)[0].resize(4);
    (*meshSidesets_)[0][0].resize(0);
    (*meshSidesets_)[0][1].resize(0);
    (*meshSidesets_)[0][2].resize(0);
    (*meshSidesets_)[0][3].resize(ny_);
    // Neumann
    (*meshSidesets_)[1].resize(4);
    (*meshSidesets_)[1][0].resize(nx_);
    (*meshSidesets_)[1][1].resize(ny_);
    (*meshSidesets_)[1][2].resize(nx_);
    (*meshSidesets_)[1][3].resize(0);

    for (int i = 0; i < nx_; ++i)
    {
      (*meshSidesets_)[1][0][i] = i;
    }
    for (int i = 0; i < ny_; ++i)
    {
      (*meshSidesets_)[1][1][i] = (i + 1) * nx_ - 1;
    }
    for (int i = 0; i < nx_; ++i)
    {
      (*meshSidesets_)[1][2][i] = i + nx_ * (ny_ - 1);
    }
    for (int i = 0; i < ny_; ++i)
    {
      (*meshSidesets_)[0][3][i] = i * nx_;
    }

  } // computeSidesets

  ROL::Ptr<std::vector<std::vector<std::vector<int>>>> getSidesets(
      const bool verbose = false,
      std::ostream &outStream = std::cout) const
  {
    if (verbose)
    {
      outStream << "Mesh_shallow_ice: getSidesets called" << std::endl;
      outStream << "Mesh_shallow_ice: numSidesets = " << meshSidesets_->size() << std::endl;
    }
    return meshSidesets_;
  }

}; // MeshManager_shallow_ice
