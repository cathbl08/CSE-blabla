#ifndef FILE_TILE_HPP
#define FILE_TILE_HPP

#include <array>
#include "simd.hpp"
#include "ordering.hpp" // For ORDERING enum

namespace ASC_bla
{
    using namespace ASC_HPC;

    template <typename T, size_t H, size_t W, ORDERING ORD = RowMajor>
    class Tile;

    // Addition for Tiles
    template <typename T, size_t H, size_t W, ORDERING ORD>
    Tile<T, H, W, ORD> operator+(const Tile<T, H, W, ORD>& a, const Tile<T, H, W, ORD>& b)
    {
        Tile<T, H, W, ORD> result;
        if constexpr (ORD == RowMajor)
        {
            for (size_t i = 0; i < H; ++i)
            {
                result.regs[i] = a.regs[i] + b.regs[i];
            }
        }
        else // ColMajor
        {
            for (size_t j = 0; j < W; ++j)
            {
                result.regs[j] = a.regs[j] + b.regs[j];
            }
        }
        return result;
    }

    // Scalar multiplication for Tiles
    template <typename S, typename T, size_t H, size_t W, ORDERING ORD>
    auto operator*(S scalar, const Tile<T, H, W, ORD>& tile)
    {
        Tile<T, H, W, ORD> result;
        if constexpr (ORD == RowMajor)
        {
            for (size_t i = 0; i < H; ++i)
            {
                result.regs[i] = scalar * tile.regs[i];
            }
        }
        else // ColMajor
        {
            for (size_t j = 0; j < W; ++j)
            {
                result.regs[j] = scalar * tile.regs[j];
            }
        }
        return result;
    }

    template <typename T, size_t H, size_t W, ORDERING ORD>
    class Tile
    {
    public:
        static constexpr size_t RegSize = (ORD == RowMajor) ? W : H;
        static constexpr size_t NumRegs = (ORD == RowMajor) ? H : W;

        std::array<SIMD<T, RegSize>, NumRegs> regs;

        Tile() = default;

        // Load from memory
        template <typename MAT>
        Tile(const MAT& mat, size_t row_start, size_t col_start)
        {
            if constexpr (ORD == RowMajor)
            {
                for (size_t i = 0; i < H; ++i)
                {
                    // uses that mat has a data() pointer and is row-major by default
                    regs[i] = SIMD<T, W>(&mat(row_start + i, col_start));
                }
            }
            else // ColMajor
            {
                for (size_t j = 0; j < W; ++j)
                {
                    // This would be more complex for a generic matrix expression.
                    // For a contiguous ColMajor MatrixView, it would be efficient.
                    // For simplicity, we do a gather operation here.
                    std::array<T, H> col_data;
                    for(size_t i=0; i<H; ++i)
                        col_data[i] = mat(row_start+i, col_start+j);
                    regs[j] = SIMD<T, H>(col_data.data());
                }
            }
        }

        // Store to memory
        template <typename MAT>
        void Store(MAT& mat, size_t row_start, size_t col_start) const
        {
            if constexpr (ORD == RowMajor)
            {
                for (size_t i = 0; i < H; ++i)
                {
                    regs[i].store(&mat(row_start + i, col_start));
                }
            }
            // ColMajor store would be similar but column-wise.
        }
    };
} // namespace ASC_bla

#endif // FILE_TILE_HPP
