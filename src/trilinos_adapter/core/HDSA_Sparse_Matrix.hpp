#ifndef HDSA_SPARSE_MATRIX_HPP
#define HDSA_SPARSE_MATRIX_HPP

#include "Tpetra_CrsMatrix_decl.hpp"

namespace HDSA
{

    template <class RealT,
              class LO = Tpetra::Map<>::local_ordinal_type,
              class GO = Tpetra::Map<>::global_ordinal_type,
              class Node = Tpetra::Map<>::node_type>
    class Sparse_Matrix
    {

    private:
        HDSA::Ptr<Tpetra::CrsMatrix<RealT, LO, GO, Node> > A_;

    public:
        // Null constructor
        Sparse_Matrix(void)
        {
        }

        Sparse_Matrix(HDSA::Ptr<Tpetra::CrsMatrix<RealT, LO, GO, Node> > & A): A_(A)
        {
        }

        ~Sparse_Matrix()
        {
        }

        HDSA::Ptr<HDSA::Sparse_Matrix<RealT> > clone(void) const
        {
            HDSA::Ptr<Tpetra::CrsMatrix<RealT, LO, GO, Node> > B = HDSA::makePtr<Tpetra::CrsMatrix<RealT, LO, GO, Node> >(*A_,Teuchos::Copy);
            HDSA::Ptr<HDSA::Sparse_Matrix<RealT> > B_sm = HDSA::makePtr<HDSA::Sparse_Matrix<RealT> > (B);
            return B_sm;
        }
    };

}

#endif
