/***********************************************************************
 HdsaLib - A library for Hyper-differential Sensitivity Analysis
 
 Questions? Contact Joseph Hart (joshart@sandia.gov)
************************************************************************/

#ifndef HDSA_SPARSE_MATRIX_HPP
#define HDSA_SPARSE_MATRIX_HPP

#include "Tpetra_CrsMatrix_decl.hpp"
#include "TpetraExt_MatrixMatrix.hpp"
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
        bool is_symmetric_;

    public:
        // Null constructor
        Sparse_Matrix(void)
        {
        }

        Sparse_Matrix(HDSA::Ptr<Tpetra::CrsMatrix<RealT, LO, GO, Node>> &A, bool is_symmetric = false) : A_(A), is_symmetric_(is_symmetric)
        {
        }

        Sparse_Matrix(HDSA::Vector<RealT> &vec, bool reciprocate_diag)
        {
            HDSA::Tpetra_Vector<RealT> tpetra_vec = dynamic_cast<HDSA::Tpetra_Vector<RealT> &>(vec);
            HDSA::Ptr<Tpetra::MultiVector<RealT, LO, GO, Node>> d = tpetra_vec.getVector();

            HDSA::Ptr<const Tpetra::Map<LO, GO, Node>> map = d->getMap();
            A_ = HDSA::makePtr<Tpetra::CrsMatrix<RealT, LO, GO, Node>>(map, 1);

            auto sLocal = d->getLocalViewHost(Tpetra::Access::ReadOnly);
            const size_t localNum = map->getLocalNumElements();
            for (size_t i = 0; i < localNum; ++i)
            {
                GO gid = map->getGlobalElement(i);
                RealT val = sLocal(i, 0);
                if (reciprocate_diag)
                {
                    val = Teuchos::ScalarTraits<RealT>::one() / val;
                }
                A_->insertGlobalValues(gid, Teuchos::tuple(gid), Teuchos::tuple(val));
            }

            A_->fillComplete(map, map);
        }

        ~Sparse_Matrix()
        {
        }

        HDSA::Ptr<const Teuchos::Comm<int>> Get_Comm(void) const
        {
            return A_->getComm();
        }

        void Set_Symmetric(void)
        {
            is_symmetric_ = true;
        }

        bool Is_Symmetric(void) const
        {
            return is_symmetric_;
        }

        HDSA::Ptr<HDSA::Sparse_Matrix<RealT>> Clone(int max_entries_per_row = 0) const
        {
            if (max_entries_per_row == 0)
            {
                max_entries_per_row = A_->getGlobalMaxNumRowEntries();
            }
            HDSA::Ptr<Tpetra::CrsMatrix<RealT, LO, GO, Node>> B = HDSA::makePtr<Tpetra::CrsMatrix<RealT, LO, GO, Node>>(A_->getRowMap(), max_entries_per_row);
            HDSA::Ptr<HDSA::Sparse_Matrix<RealT>> B_sm = HDSA::makePtr<HDSA::Sparse_Matrix<RealT>>(B);
            B_sm->Scale(0.0);
            return B_sm;
        }

        // Compute C = this * B, with options for transposes
        void Matrix_Matrix_Multiply(HDSA::Sparse_Matrix<RealT> &C, const HDSA::Sparse_Matrix<RealT> &B, bool A_trans = false, bool B_trans = false) const
        {
            HDSA::Ptr<Tpetra::CrsMatrix<RealT, LO, GO, Node>> B_tpetra = B.Get_Tpetra_Matrix();
            HDSA::Ptr<Tpetra::CrsMatrix<RealT, LO, GO, Node>> C_tpetra = C.Get_Tpetra_Matrix();
            Tpetra::MatrixMatrix::Multiply(*A_, A_trans, *B_tpetra, B_trans, *C_tpetra);
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

        void Scale(RealT alpha)
        {
            A_->scale(alpha);
        }

        void Set(HDSA::Sparse_Matrix<RealT> &B)
        {
            // Prepare A for updates
            A_->resumeFill();

            // Loop over each row of B and copy its entries to A
            for (Tpetra::global_size_t row = 0; row < B.Get_Tpetra_Matrix()->getGlobalNumRows(); ++row)
            {
                if (A_->getRowMap()->isNodeGlobalElement(row))
                {
                    size_t numEntries = B.Get_Tpetra_Matrix()->getNumEntriesInGlobalRow(row);
                    typename Tpetra::CrsMatrix<>::nonconst_global_inds_host_view_type indices;
                    typename Tpetra::CrsMatrix<>::nonconst_values_host_view_type values;
                    Kokkos::resize(indices, numEntries);
                    Kokkos::resize(values, numEntries);

                    // Get the global row copy from B
                    B.Get_Tpetra_Matrix()->getGlobalRowCopy(row, indices, values, numEntries);
                    // Replace the entries in A with those from B
                    A_->insertGlobalValues(row, numEntries, &values[0], &indices[0]);
                }
            }

            // Complete the fill process
            A_->fillComplete();
        }

        void Scaled_Plus(const RealT &alpha, const HDSA::Sparse_Matrix<RealT> &B)
        {
            A_->resumeFill(); // Prepare A for updates
            // Loop over each row of B and copy its entries to A
            for (Tpetra::global_size_t row = 0; row < B.Get_Tpetra_Matrix()->getGlobalNumRows(); ++row)
            {
                if (A_->getRowMap()->isNodeGlobalElement(row))
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
            }
            A_->fillComplete(); // Complete the fill process
        }

        void Begin_Fill(void)
        {
            A_->resumeFill();
        }

        void End_Fill(void)
        {
            A_->fillComplete(); // Complete the fill process
        }

        void Set_Entry(const int &i, const int &j, const RealT &val)
        {
            const long long col = static_cast<long long>(j);
            A_->insertGlobalValues(i, 1, &val, &col);
        }

        void Get_Global_Row(const int &row, std::vector<int> &col_indices, std::vector<RealT> &vals)
        {
            size_t numEntries = A_->getNumEntriesInGlobalRow(row);
            typename Tpetra::CrsMatrix<>::nonconst_global_inds_host_view_type indices;
            typename Tpetra::CrsMatrix<>::nonconst_values_host_view_type values;
            Kokkos::resize(indices, numEntries);
            Kokkos::resize(values, numEntries);
            A_->getGlobalRowCopy(row, indices, values, numEntries);

            col_indices.resize(numEntries);
            vals.resize(numEntries);
            for (int k = 0; k < numEntries; k++)
            {
                col_indices[k] = indices[k];
                vals[k] = values[k];
            }
        }

        int Get_Number_of_Rows(void) const
        {
            return A_->getGlobalNumRows();
        }

        int Get_Number_of_Columns(void) const
        {
            return A_->getGlobalNumCols();
        }

        int Get_Max_Nonzeros_Per_Row(void) const
        {
            return A_->getGlobalMaxNumRowEntries();
        }

        bool Is_Row_Owned(int row) const
        {
            return A_->getRowMap()->isNodeGlobalElement(row);
        }
    };

}

#endif
