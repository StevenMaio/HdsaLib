#ifndef HDSA_SPARSE_MATRIX_HPP
#define HDSA_SPARSE_MATRIX_HPP

#include "Tpetra_CrsMatrix_decl.hpp"
#include "HDSA_Tpetra_Vector.hpp"

namespace HDSA
{

    template <class RealT,
              class LO = Tpetra::Map<>::local_ordinal_type,
              class GO = Tpetra::Map<>::global_ordinal_type,
              class Node = Tpetra::Map<>::node_type>
    class Sparse_Matrix
    {

    private:
        HDSA::Ptr<Tpetra::CrsMatrix<RealT, LO, GO, Node>> A_;

    public:
        // Null constructor
        Sparse_Matrix(void)
        {
        }

        Sparse_Matrix(HDSA::Ptr<Tpetra::CrsMatrix<RealT, LO, GO, Node>> &A) : A_(A)
        {
        }

        ~Sparse_Matrix()
        {
        }

        HDSA::Ptr<HDSA::Sparse_Matrix<RealT>> Clone(void) const
        {
            HDSA::Ptr<Tpetra::CrsMatrix<RealT, LO, GO, Node>> B = HDSA::makePtr<Tpetra::CrsMatrix<RealT, LO, GO, Node>>(*A_, Teuchos::Copy);
            HDSA::Ptr<HDSA::Sparse_Matrix<RealT>> B_sm = HDSA::makePtr<HDSA::Sparse_Matrix<RealT>>(B);
            return B_sm;
        }

        RealT Frobenius_Norm(void) const
        {
            return A_->getFrobeniusNorm();
        }

        HDSA::Ptr<Tpetra::CrsMatrix<RealT, LO, GO, Node>> Get_Tpetra_Matrix(void) const
        {
            return A_;
        }

        void Apply(HDSA::Vector<RealT> &x_out, const HDSA::Vector<RealT> &x_in) const
        {
            const HDSA::Tpetra_Vector<RealT> &ex_in = dynamic_cast<const HDSA::Tpetra_Vector<RealT> &>(x_in);
            HDSA::Tpetra_Vector<RealT> &ex_out = dynamic_cast<HDSA::Tpetra_Vector<RealT> &>(x_out);
            A_->apply(*ex_in.getVector(), *ex_out.getVector());
        }

        void Set(HDSA::Sparse_Matrix<RealT> &B)
        {
            // Prepare A for updates
            A_->resumeFill();

            // Loop over each row of B and copy its entries to A
            for (Tpetra::global_size_t row = 0; row < B.Get_Tpetra_Matrix()->getGlobalNumRows(); ++row)
            {
                size_t numEntries = B.Get_Tpetra_Matrix()->getNumEntriesInGlobalRow(row);
                typename Tpetra::CrsMatrix<>::nonconst_global_inds_host_view_type indices;
                typename Tpetra::CrsMatrix<>::nonconst_values_host_view_type values;
                Kokkos::resize(indices, numEntries);
                Kokkos::resize(values, numEntries);

                // Get the global row copy from B
                B.Get_Tpetra_Matrix()->getGlobalRowCopy(row, indices, values, numEntries);
                // Replace the entries in A with those from B
                A_->replaceGlobalValues(row, indices, values);
            }

            // Complete the fill process
            A_->fillComplete();
        }

        void Scaled_Plus(RealT &alpha, HDSA::Sparse_Matrix<RealT> &B)
        {
            A_->resumeFill(); // Prepare A for updates
            // Loop over each row of B and copy its entries to A
            for (Tpetra::global_size_t row = 0; row < B.Get_Tpetra_Matrix()->getGlobalNumRows(); ++row)
            {
                size_t numEntries = B.Get_Tpetra_Matrix()->getNumEntriesInGlobalRow(row);
                typename Tpetra::CrsMatrix<>::nonconst_global_inds_host_view_type indices;
                typename Tpetra::CrsMatrix<>::nonconst_values_host_view_type values;
                Kokkos::resize(indices, numEntries);
                Kokkos::resize(values, numEntries);

                // Get the global row copy from B
                B.Get_Tpetra_Matrix()->getGlobalRowCopy(row, indices, values, numEntries);
                for (int k = 0; k < numEntries; k++)
                {
                    values[k] = alpha * values[k];
                }
                // Replace the entries in A with those from B
                Teuchos::ArrayView<const GO> indicesView(indices.data(), numEntries);
                Teuchos::ArrayView<const RealT> valuesView(values.data(), numEntries);
                A_->sumIntoGlobalValues(row, indicesView, valuesView);
            }
            A_->fillComplete(); // Complete the fill process
        }
    };

}

#endif
