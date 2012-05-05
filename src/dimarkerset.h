#ifndef DIMARKERSET_H_INCLUDED
#define DIMARKERSET_H_INCLUDED

#include "vector.h"
#include "domain.h"
#include "dimarker.h"
#include <list>
#include <map>
#include <algorithm>

namespace lfield {
class dimarkerset: public std::list<dimarker>
{
        eufield::domain* domain;
        float lmax;

        void init_I(float *I)
        {
            for (int i = 0; i < domain->grid->nx; ++i)
                for (int j = 0; j < domain->grid->ny; ++j)
                {
                    if (domain->grid->bcs[i][j] != 0) continue;

                    int idx = domain->grid->operator ()(i, j);
                    I[idx] = domain->phi[idx] > 0 ? domain->phi[idx] : 0.;
                }

        }

        void init_deltas(float **gradI, float *grad2I, float *delta)
        {
            for (int i = 0; i < domain->grid->nx; ++i)
                for (int j = 0; j < domain->grid->ny; ++j)
                {
                    if (domain->grid->bcs[i][j] != 0) continue;

                    int idx = domain->grid->operator ()(i, j);
                    geometry::vector gI = domain->getvec(gradI, idx);
                    geometry::vector gPhi = domain->getvec(domain->gradphi, idx);
                    float g2I = grad2I[idx], g2phi = domain->grad2phi[idx], l2phi = gPhi.length2();
                    delta[idx] = g2I / l2phi - geometry::dot(gI, gPhi) * g2phi / l2phi / l2phi;
                }

        }

        bool check_intersect(geometry::vector r1, geometry::vector r2, float phi1, float phi2, geometry::vector& intersect)
        {
            if (phi1 * phi2 > 0 || (fabs(phi1) < 1e-6 && fabs(phi2) < 1e-6)) return false;

            float delta = fabs(phi1 - phi2);
            intersect = r1 * (fabs(phi2) / delta) + r2 * (fabs(phi1) / delta);

            return true;
        }

        void insert_markers(geometry::vector sec1, geometry::vector sec2, int current_count)
        {
            geometry::vector s = sec1 - sec2;
            /*
             if (s.length2() < domain->grid->h * domain->grid->h / 4.)
             return;
             */

            dimarker result;

            result.pos = 0.5 * (sec1 + sec2);
            s.normalize();
            result.normal[0] = +s[1];
            result.normal[1] = -s[0];
            result.normal[2] = s[2];
            geometry::vector local_norm = domain->getvec(domain->gradphi, result.pos);
            /*
             if (geometry::dot(local_norm, result.normal) < 0)
             result.normal *= -1.;
             */

            int i, j;
            domain->grid->getidx(result.pos, i, j);
            int idx = domain->grid->operator ()(i, j);
            float phi = domain->phi[idx];
            if (result.getphi(domain->grid->getpos(i, j)) * phi < 0.) result.normal *= -1;

            //result.normal = domain->getvec(domain->gradphi, result.pos);
            /*
             if (false && result.normal.length2() < 0.95)
             {
             return;

             int i, j;
             domain->grid->getidx(result.pos, i, j);
             int idx = domain->grid->operator ()(i, j);
             int idxi = domain->grid->operator ()(i + 1, j);
             int idxj = domain->grid->operator ()(i, j + 1);
             int idxij = domain->grid->operator ()(i + 1, j + 1);

             float phi = domain->phi[idx];
             float phii = domain->phi[idxi];
             float phij = domain->phi[idxj];
             float phiij = domain->phi[idxij];

             result.normal[0] = (phii + phiij - phi - phij) / 2. / domain->grid->h;
             result.normal[0] = (phiij + phij - phi - phii) / 2. / domain->grid->h;
             result.normal[0] = 0.;
             }
             result.normal.normalize();
             */
            /*
             if (fabs(geometry::dot(result.normal, s)) > 0.9)
             return;
             */

            //push_back(result);
            float l = (sec2 - sec1).length();
            static const float nmax = 4;

            int n = (nmax * l) / lmax;
            //n = n > 0 ? n : 1;
            n -= current_count;
            if (n <= 0) return;
            //n = n % 2 == 1 ? n + 1 : n;
            float dalpha = 1. / ((float) n);
            float alpha = dalpha / 2.;
            for (int i = 0; i < n; ++i)
            {
                dimarker d;
                d.pos = alpha * sec1 + (1. - alpha) * sec2;
                d.normal = result.normal;
                push_back(d);
                //push_back(dimarker::create(domain, alpha * sec1 + (1. - alpha) * sec2));
                alpha += dalpha;
            }

        }

        bool check_normals(int idx, int idxi, int idxj, int idxij)
        {
            geometry::vector n = domain->getvec(domain->gradphi, idx);
            geometry::vector ni = domain->getvec(domain->gradphi, idxi);
            geometry::vector nj = domain->getvec(domain->gradphi, idxj);
            geometry::vector nij = domain->getvec(domain->gradphi, idxij);

            return false;
        }

    public:

        dimarkerset(eufield::domain* adomain)
        {
            domain = adomain;
            lmax = domain->grid->h * sqrt(2.);
        }

        void init_markers_delta()
        {
            float *I = new float[domain->grid->n];
            float *grad2I = new float[domain->grid->n];
            float *delta = new float[domain->grid->n];
            float **gradI = new float*[3];
            gradI[0] = new float[domain->grid->n];
            gradI[1] = new float[domain->grid->n];
            gradI[2] = new float[domain->grid->n];

            init_I(I);
            domain->calcgrad(I, gradI);
            domain->calcgrad2(I, grad2I);
            init_deltas(gradI, grad2I, delta);

            for (int i = 0; i < domain->grid->nx; ++i)
                for (int j = 0; j < domain->grid->ny; ++j)
                {
                    if (domain->grid->bcs[i][j] != 0) continue;

                    int idx = domain->grid->operator ()(i, j);
                    int idxi = domain->grid->operator ()(i + 1, j);
                    int idxj = domain->grid->operator ()(i, j + 1);
                    int idxij = domain->grid->operator ()(i + 1, j + 1);

                    if (!(domain->phi[idx] * domain->phi[idxi] <= 0 || domain->phi[idxi] * domain->phi[idxij] <= 0 || domain->phi[idxij] * domain->phi[idxj] <= 0 || domain->phi[idxj] * domain->phi[idx] <= 0)) continue;

                    geometry::vector r = domain->grid->getpos(i, j);
                    geometry::vector ri = domain->grid->getpos(i + 1, j);
                    geometry::vector rj = domain->grid->getpos(i, j + 1);
                    geometry::vector rij = domain->grid->getpos(i + 1, j + 1);

                    float dd = delta[idx];
                    float ddi = delta[idxi];
                    float ddj = delta[idxj];
                    float ddij = delta[idxij];

                    float lg = domain->getvec(domain->gradphi, idx).length();
                    float lgi = domain->getvec(domain->gradphi, idxi).length();
                    float lgj = domain->getvec(domain->gradphi, idxj).length();
                    float lgij = domain->getvec(domain->gradphi, idxij).length();

                    geometry::vector x = (dd * r * lg + ddi * ri * lgi + ddj * rj * lgj + ddij * rij * lgij) / (dd * lg + ddi * lgi + ddj * lgj + ddij * lgij);

                    push_back(dimarker(x, geometry::vector()));
                }

            delete[] I;
            delete[] grad2I;
            delete[] gradI[0];
            delete[] gradI[1];
            delete[] gradI[2];
            delete gradI;
            delete delta;
        }

        void init_markers_intersects()
        {
            for (int i = 0; i < domain->grid->nx; ++i)
                for (int j = 0; j < domain->grid->ny; ++j)
                {

                    if (domain->grid->bcs[i][j] != 0) continue;

                    int idx = domain->grid->operator ()(i, j);
                    int idxi = domain->grid->operator ()(i + 1, j);
                    int idxj = domain->grid->operator ()(i, j + 1);
                    int idxij = domain->grid->operator ()(i + 1, j + 1);

                    geometry::vector r = domain->grid->getpos(i, j);
                    geometry::vector ri = domain->grid->getpos(i + 1, j);
                    geometry::vector rj = domain->grid->getpos(i, j + 1);
                    geometry::vector rij = domain->grid->getpos(i + 1, j + 1);

                    float dd = domain->phi[idx];
                    float ddi = domain->phi[idxi];
                    float ddj = domain->phi[idxj];
                    float ddij = domain->phi[idxij];

                    geometry::vector sectB, sectR, sectL, sectT;

                    bool sectedB = check_intersect(r, ri, dd, ddi, sectB);
                    bool sectedR = check_intersect(ri, rij, ddi, ddij, sectR);
                    bool sectedT = check_intersect(rij, rj, ddij, ddj, sectT);
                    bool sectedL = check_intersect(rj, r, ddj, dd, sectL);

                    if (!(sectedB || sectedR || sectedT || sectedL) || (sectedB && sectedR && sectedT && sectedL)) continue;

                    if (sectedB && sectedR) insert_markers(sectB, sectR, 0);
                    if (sectedR && sectedT) insert_markers(sectR, sectT, 0);
                    if (sectedT && sectedL) insert_markers(sectT, sectL, 0);
                    if (sectedL && sectedB) insert_markers(sectL, sectB, 0);

                    if (sectedL && sectedR && !(sectedT || sectedB)) insert_markers(sectL, sectR, 0);
                    if (sectedT && sectedB && !(sectedL || sectedR)) insert_markers(sectT, sectB, 0);

                }
        }

        void advect(float dt)
        {
            for (iterator i = begin(); i != end(); ++i)
                i->advect_RK2(domain, dt);
        }

        void correct()
        {
            bool* correcteds = new bool[domain->grid->n];
            dimarker** correctors = new dimarker*[domain->grid->n];
            float* oldls = new float[domain->grid->n];

            std::copy(domain->phi, domain->phi + domain->grid->n, oldls);
            std::fill_n(correcteds, domain->grid->n, false); // * domain->grid->h);

            for (iterator i = begin(); i != end(); ++i)
                i->correct(domain, correcteds, correctors);

            /*
             float* phi = domain->phi;
             for (int i = 0; i < domain->grid->nx && false; ++i)
             for (int j = 0; j < domain->grid->ny && false; ++j)
             {
             if (domain->grid->bcs[i][j] != 0)
             continue;

             int idx = (*domain->grid)(i, j);
             if (phi[idx] * oldls[idx] > 0.)
             continue;

             eufield::molecule m = domain->grid->get_molecule(i, j);
             if (phi[idx] * oldls[m.xm] < 0. && phi[idx] * oldls[m.xp] < 0. && phi[idx] * oldls[m.ym] < 0. && phi[idx] * oldls[m.yp] < 0.)
             phi[idx] = oldls[idx];//(phi[m.xm] + phi[m.xp] + phi[m.ym] + phi[m.yp]) / 4.;
             if (phi[idx] * phi[m.xm] < 0. && phi[idx] * phi[m.xp] < 0. && phi[idx] * phi[m.ym] < 0. && phi[idx] * phi[m.yp] < 0.)
             phi[idx] = (phi[m.xm] + phi[m.xp] + phi[m.ym] + phi[m.yp]) / 4.;
             }
             */

            delete[] correcteds;
            delete[] correctors;
            delete[] oldls;
        }

        void add_markers()
        {
            int* count = new int[domain->grid->n];
            std::fill_n(count, domain->grid->n, 0);

            int i, j;

            for (iterator iter = begin(); iter != end(); ++iter)
            {
                domain->grid->getidx(iter->pos, i, j);
                ++count[domain->grid->operator ()(i, j)];
            }

            for (int i = 0; i < domain->grid->nx; ++i)
                for (int j = 0; j < domain->grid->ny; ++j)
                {
                    if (domain->grid->bcs[i][j] != 0) continue;

                    int idx = domain->grid->operator ()(i, j);
                    int idxi = domain->grid->operator ()(i + 1, j);
                    int idxj = domain->grid->operator ()(i, j + 1);
                    int idxij = domain->grid->operator ()(i + 1, j + 1);

                    int c = count[idx];
                    if (c > 0) continue;
                    //				if (1.0 / fabs(domain->grad2phi[idx]) < domain->grid->h)
                    //					continue;

                    float dd = domain->phi[idx];
                    float ddi = domain->phi[idxi];
                    float ddj = domain->phi[idxj];
                    float ddij = domain->phi[idxij];

                    //				if (fabs(dd) > 1. || fabs(ddi) > 1. || fabs(ddj) > 1. || fabs(ddij) > 1.)
                    //					continue;

                    geometry::vector r = domain->grid->getpos(i, j);
                    geometry::vector ri = domain->grid->getpos(i + 1, j);
                    geometry::vector rj = domain->grid->getpos(i, j + 1);
                    geometry::vector rij = domain->grid->getpos(i + 1, j + 1);

                    geometry::vector sectB, sectR, sectL, sectT;

                    bool sectedB = check_intersect(r, ri, dd, ddi, sectB);
                    bool sectedR = check_intersect(ri, rij, ddi, ddij, sectR);
                    bool sectedT = check_intersect(rij, rj, ddij, ddj, sectT);
                    bool sectedL = check_intersect(rj, r, ddj, dd, sectL);

                    if (!(sectedB || sectedR || sectedT || sectedL) || (sectedB && sectedR && sectedT && sectedL)) continue;

                    if (sectedB && sectedR) insert_markers(sectB, sectR, c);
                    if (sectedR && sectedT) insert_markers(sectR, sectT, c);
                    if (sectedT && sectedL) insert_markers(sectT, sectL, c);
                    if (sectedL && sectedB) insert_markers(sectL, sectB, c);

                    if (sectedL && sectedR && !(sectedT || sectedB)) insert_markers(sectL, sectR, c);
                    if (sectedT && sectedB && !(sectedL || sectedR)) insert_markers(sectT, sectB, c);

                }
        }

        void delete_bad_vecs(std::list<iterator>* lst)
        {
            int ndel = 0;
            std::list<iterator>::iterator* del = new std::list<iterator>::iterator[lst->size()];
            bool* ok = new bool[lst->size()];
            std::fill_n(ok, lst->size(), false);

            int ii = 0, jj;
            for (std::list<iterator>::iterator i = lst->begin(); i != lst->end(); ++i)
            {
                jj = ii - 1;
                for (std::list<iterator>::iterator j = i; j != lst->end(); ++j)
                {
                    ++jj;
                    if (j == i) continue;

                    if (j == lst->end()) continue;

                    if (geometry::dot((*i)->pos, (*j)->pos) > 0.9) ok[ii] = ok[jj] = true;
                }
                if (!ok[ii]) del[ndel++] = i;
            }

            for (int i = 0; i < ndel; ++i)
            {
                erase(*del[i]);
                lst->erase(del[i]);
            }
        }

        void delete_last_n(std::list<iterator>* lst, int n)
        {
            while (lst->size() > n)
                lst->erase(lst->begin());

            for (std::list<iterator>::iterator i = lst->begin(); i != lst->end(); ++i)
                erase(*i);
        }

        void delete_markers()
        {
            static const uint NMAX = 5;
            std::list<iterator>** cellmarkers = new std::list<iterator>*[domain->grid->n];
            for (int i = 0; i < domain->grid->n; ++i)
                cellmarkers[i] = new std::list<lfield::dimarkerset::iterator>();

            for (iterator iter = begin(); iter != end(); ++iter)
            {
                int i, j;
                domain->grid->getidx(iter->pos, i, j);
                cellmarkers[domain->grid->operator ()(i, j)]->push_back(iter);
            }

            for (int i = 0; i < domain->grid->n; ++i)
            {
                /*
                 int idx = domain->grid->operator ()(i, j);
                 int idxi = domain->grid->operator ()(i + 1, j);
                 int idxj = domain->grid->operator ()(i, j + 1);
                 int idxij = domain->grid->operator ()(i + 1, j + 1);

                 float dd = domain->phi[idx];
                 float ddi = domain->phi[idxi];
                 float ddj = domain->phi[idxj];
                 float ddij = domain->phi[idxij];
                 */

                /*
                 for (uint jj = 0; jj < NMAX; ++jj)
                 cellmarkers[ii]->erase(cellmarkers[ii]->begin());
                 */

                if (cellmarkers[i]->size() > NMAX) delete_bad_vecs(cellmarkers[i]);

                if (false & cellmarkers[i]->size() > NMAX) delete_last_n(cellmarkers[i], NMAX - cellmarkers[i]->size());

            }

            for (int i = 0; i < domain->grid->n; ++i)
                delete cellmarkers[i];
            delete[] cellmarkers;
        }

        void write(std::string fname)
        {
            std::ofstream stream(fname.c_str());

            stream << "# vtk DataFile Version 2.0" << std::endl;
            stream << "Data Output" << std::endl;
            stream << "ASCII" << std::endl;
            stream << "DATASET UNSTRUCTURED_GRID" << std::endl;

            stream << "POINTS " << size() << " float" << std::endl;
            for (iterator i = begin(); i != end(); ++i)
                stream << i->pos << std::endl;

            stream << "CELLS " << size() << " " << size() * 2 << std::endl;
            int n = 0;
            for (iterator i = begin(); i != end(); ++i)
                stream << "1 " << n++ << std::endl;

            stream << "CELL_TYPES " << size() << std::endl;
            for (uint i = 0; i < size(); ++i)
                stream << "1" << std::endl;

            stream << "POINT_DATA " << size() << std::endl;
            stream << "VECTORS normals float" << std::endl;
            for (iterator i = begin(); i != end(); ++i)
                stream << i->normal << std::endl;

            stream.close();
        }
};
}

#endif /* DIMARKERSET_H_INCLUDED */

