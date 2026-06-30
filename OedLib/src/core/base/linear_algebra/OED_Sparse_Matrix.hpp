#pragma once

#include <unordered_map>
#include <utility>
#include <cmath>

#include "OED_Vector.hpp"
#include "OED_Linear_Operator.hpp"

namespace OED
{
    template <class RealT>
    class Sparse_Matrix : public Linear_Operator<RealT>
    {
    private:
        std::unordered_map<std::pair<int, int>, RealT> entries_;

    public:
        Sparse_Matrix() : Linear_Operator<RealT>() {}

        void Add_Entry(int row, int col, RealT val)
        {
            if (std::abs(val) > 1e-12)
            {
                std::pair<int, int> coord{row, col};
                this->entries_[coord] = val;
                std::cout << "Sparse_Matrix::Add_Entry added (" << row << ", " << col << ")=" << val << std::endl;
            }
        }

        void Apply(OED::Vector<RealT> &y, const OED::Vector<RealT> &x) const override
        {
            y.Zeros();
            for (auto &entry: this->entries_)
            {
                std::pair<int, int> coord = entry.first;
                RealT val = entry.second;

                RealT y_i = y.Get_Entry(coord.first);
                y_i += val * x.Get_Entry(coord.second);

                y.Set_Entry(coord.first, y_i);
            }
        }

        void Transpose_Apply(OED::Vector<RealT> &x, const OED::Vector<RealT> &y) const 
        {
            x.Zeros();
            for (auto &entry: this->entries_)
            {
                std::pair<int, int> coord = entry.first;
                RealT val = entry.second;

                RealT x_j = x.Get_Entry(coord.second);
                x_j += val * y.Get_Entry(coord.first);

                x.Set_Entry(coord.second, x_j);
            }
        }
    };
}
