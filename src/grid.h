#ifndef GRID_H_INCLUDED
#define GRID_H_INCLUDED

#include "vector.h"

namespace eufield {

struct molecule
{
        int i, j;
        int p, xm, xp, ym, yp;

        molecule()
        {
            p = xm = xp = ym = yp = i = j = 0;
        }

        molecule(int ai, int aj, int ap, int axm, int axp, int aym, int ayp)
        {
            i = ai;
            j = aj;
            p = ap;
            xm = axm;
            xp = axp;
            ym = aym;
            yp = ayp;
        }
};

class grid
{
    public:
        int** idxs;
        int** bcs;
        int nx, ny, n;
        float h;
        float startx, starty;

        grid(int nx, int ny, int n)
        {
            this->nx = nx;
            this->ny = ny;
            this->n = n;

            idxs = new int*[ny];
            bcs = new int*[ny];
            for (int i = 0; i < ny; ++i)
            {
                idxs[i] = new int[ny];
                std::fill_n(idxs[i], ny, 0);
                bcs[i] = new int[ny];
                std::fill_n(bcs[i], ny, 0);
            }
        }

        ~grid()
        {
            for (int i = 0; i < ny; ++i)
            {
                delete[] idxs[i];
                delete[] bcs[i];
            }
            delete[] idxs;
            delete[] bcs;
        }

        int& operator()(int i, int j)
        {
            return idxs[i][j];
        }

        molecule get_molecule(int i, int j)
        {
            molecule result;
            result.i = i;
            result.j = j;
            result.p = (*this)(i, j);
            if (result.p < 0) return result;
            result.xp = (*this)(i + 1, j);
            result.xm = (*this)(i - 1, j);
            result.yp = (*this)(i, j + 1);
            result.ym = (*this)(i, j - 1);
            return result;
        }

        geometry::vector getpos(int i, int j)
        {
            return geometry::vector(startx + i * h, starty + j * h, 0);
        }

        bool getidx(geometry::vector v, int& i, int& j)
        {
            i = floor((v[0] - startx) / h) + 1.e-5;
            j = floor((v[1] - starty) / h) + 1.e-5;

            if ((i < 0) || (i >= nx) || (j < 0) || (j >= ny)) return false;

            return bcs[i][j] == 0;
        }

        bool isinside(geometry::vector v)
        {
            int i, j;
            return getidx(v, i, j);
        }
};

}
#endif /*GRID_H_INCLUDED*/
