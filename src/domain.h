#ifndef DOMAIN_H_INCLUDED
#define DOMAIN_H_INCLUDED

#include "grid.h"
#include "data.h"
#include "vector.h"
#include <list>
#include <fstream>

namespace eufield {

class domain
{
        typedef std::list<data*> datalist;
        datalist dtlst;

        void get_bilinear_params(geometry::vector v, geometry::vector& p, geometry::vector& q, int& ll, int& lr, int& tl, int& tr)
        {
            int i, j;
            if (!grid->getidx(v, i, j)) throw -1;
            q = (v - grid->getpos(i, j)) / grid->h;
            p = geometry::vector(1) - q;

            ll = (*grid)(i, j);
            lr = (*grid)(i + 1, j);
            tl = (*grid)(i, j + 1);
            tr = (*grid)(i + 1, j + 1);
        }
    public:
        eufield::grid *grid;
        float *phi;
        float *grad2phi;
        float **u;
        float **gradphi;
        float **omega;
        float ***gradu;
        float* dummy;

        domain(eufield::grid *grid)
        {
            this->grid = grid;

            dtlst.push_back(new data1("phi", grid, &phi, true));
            dtlst.push_back(new data1("dummy", grid, &dummy, false));
            dtlst.push_back(new data1("grad2phi", grid, &grad2phi, false));
            dtlst.push_back(new data3("u", grid, &u, false));
            dtlst.push_back(new data3("gradphi", grid, &gradphi, false));
            dtlst.push_back(new data3("omega", grid, &omega, false));

            gradu = new float**[3];
            dtlst.push_back(new data3("gradu", grid, &gradu[0], false));
            dtlst.push_back(new data3("gradv", grid, &gradu[1], false));
            dtlst.push_back(new data3("gradw", grid, &gradu[2], false));
        }
        ~domain()
        {
            for (datalist::iterator iter = dtlst.begin(); iter != dtlst.end(); ++iter)
                delete *iter;

            delete[] gradu;
        }

        void write(std::string fname)
        {
            std::ofstream stream(fname.c_str());

            stream << "# vtk DataFile Version 2.0" << std::endl;
            stream << "Data Output" << std::endl;
            stream << "ASCII" << std::endl;
            stream << "DATASET UNSTRUCTURED_GRID" << std::endl;

            int npoints = 0;
            int* visidxs = new int[grid->n];
            std::fill_n(visidxs, grid->n, -1);
            for (int j = 0; j < grid->ny; ++j)
                for (int i = 0; i < grid->nx; ++i)
                    if (grid->bcs[i][j] == 0) visidxs[(*grid)(i, j)] = npoints++;

            stream << "POINTS " << npoints << " float" << std::endl;
            for (int j = 0; j < grid->ny; ++j)
                for (int i = 0; i < grid->nx; ++i)
                    if (grid->bcs[i][j] == 0) stream << grid->getpos(i, j) << std::endl;

            int ncells = 0;
            for (int i = 0; i < grid->nx; ++i)
                for (int j = 0; j < grid->ny; ++j)
                {
                    if (grid->bcs[i][j] != 0) continue;
                    if (grid->bcs[i + 1][j] != 0 || grid->bcs[i][j + 1] != 0 || grid->bcs[i + 1][j + 1] != 0) continue;

                    ++ncells;
                }

            int idx0, idx1, idx2, idx3;
            stream << "CELLS " << ncells << " " << ncells * 5 << std::endl;
            for (int i = 0; i < grid->nx; ++i)
                for (int j = 0; j < grid->ny; ++j)
                {
                    if (grid->bcs[i][j] != 0) continue;
                    if (grid->bcs[i + 1][j] != 0 || grid->bcs[i][j + 1] != 0 || grid->bcs[i + 1][j + 1] != 0) continue;

                    idx0 = visidxs[(*grid)(i, j)];
                    idx1 = visidxs[(*grid)(i + 1, j)];
                    idx2 = visidxs[(*grid)(i, j + 1)];
                    idx3 = visidxs[(*grid)(i + 1, j + 1)];

                    stream << "4 " << idx0 << " " << idx1 << " " << idx2 << " " << idx3 << std::endl;
                }

            stream << "CELL_TYPES " << ncells << std::endl;
            for (int i = 0; i < ncells; ++i)
                stream << "8" << std::endl;

            stream << "POINT_DATA " << npoints << std::endl;

            for (datalist::iterator iter = dtlst.begin(); iter != dtlst.end(); ++iter)
                if ((*iter)->writable) (*iter)->write(stream);

            stream.close();
            delete[] visidxs;
        }

        float getval(float* d, geometry::vector v)
        {
            geometry::vector p, q;
            int ll, lr, tl, tr;

            get_bilinear_params(v, p, q, ll, lr, tl, tr);

            return p[0] * p[1] * d[ll] + q[0] * p[1] * d[lr] + p[0] * q[1] * d[tl] + q[0] * q[1] * d[tr];
        }
        geometry::vector getvec(float** d, int idx)
        {
            return geometry::vector(d[0][idx], d[1][idx], d[2][idx]);
        }
        geometry::vector getvec(float** d, geometry::vector v)
        {
            geometry::vector p, q;
            int ll, lr, tl, tr;

            get_bilinear_params(v, p, q, ll, lr, tl, tr);

            float x = p[0] * p[1] * d[0][ll] + q[0] * p[1] * d[0][lr] + p[0] * q[1] * d[0][tl] + q[0] * q[1] * d[0][tr];
            float y = p[0] * p[1] * d[1][ll] + q[0] * p[1] * d[1][lr] + p[0] * q[1] * d[1][tl] + q[0] * q[1] * d[1][tr];
            float z = p[0] * p[1] * d[2][ll] + q[0] * p[1] * d[2][lr] + p[0] * q[1] * d[2][tl] + q[0] * q[1] * d[2][tr];

            return geometry::vector(x, y, z);
        }

        void calcgrad(float* p, float** gradp)
        {
            float hh = 2. * grid->h;

            for (int i = 0; i < grid->nx; ++i)
                for (int j = 0; j < grid->ny; ++j)
                {
                    if (grid->bcs[i][j] != 0) continue;

                    molecule m = grid->get_molecule(i, j);
                    gradp[0][m.p] = (p[m.xp] - p[m.xm]) / hh;
                    gradp[1][m.p] = (p[m.yp] - p[m.ym]) / hh;
                    gradp[2][m.p] = 0.;
                }
        }
        void calcgrad2(float* p, float* grad2p)
        {
            float h2 = grid->h * grid->h;

            for (int i = 0; i < grid->nx; ++i)
                for (int j = 0; j < grid->ny; ++j)
                {
                    if (grid->bcs[i][j] != 0) continue;

                    molecule m = grid->get_molecule(i, j);
                    grad2p[m.p] = (p[m.xp] + p[m.xm] + p[m.yp] + p[m.ym] - 4. * p[m.p]) / h2;
                }
        }
        void calcgradphi()
        {
            calcgrad(phi, gradphi);
        }
        void calcgradu()
        {
            calcgrad(u[0], gradu[0]);
            calcgrad(u[1], gradu[1]);
            calcgrad(u[2], gradu[2]);
        }
        void calcomega()
        {
            for (int i = 0; i < grid->nx; ++i)
                for (int j = 0; j < grid->ny; ++j)
                {
                    if (grid->bcs[i][j] != 0) continue;

                    int idx = (*grid)(i, j);
                    omega[0][idx] = 0.5 * (gradu[2][1][idx] - gradu[1][2][idx]);
                    omega[1][idx] = 0.5 * (gradu[0][2][idx] - gradu[2][0][idx]);
                    omega[2][idx] = 0.5 * (gradu[1][0][idx] - gradu[0][1][idx]);

                    /*
                     int idx = (*grid)(i, j);
                     float u, v, a, b, p, q;
                     u = this->u[0][idx];
                     v = this->u[1][idx];
                     a = gradu[0][0][idx];
                     b = gradu[0][1][idx];
                     p = gradu[1][0][idx];
                     q = gradu[1][1][idx];
                     omega[0][idx] = 0.;
                     omega[1][idx] = 0.;
                     omega[2][idx] = (u * (p * u + q * v) - v  (a * u + b * v)) / (u * u + v * v);
                     */
                }
        }
        void calcgrad2phi()
        {
            calcgrad2(phi, grad2phi);
        }
};
}
#endif /* DOMAIN_H_INCLUDED */
