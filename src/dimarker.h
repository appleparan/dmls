#ifndef DIMARKER_H_INCLUDED
#define DIMARKER_H_INCLUDED

#include "vector.h"
#include "domain.h"

namespace lfield {

class dimarker
{
        void apply_values(int idx, float phi, eufield::domain* domain, bool* correcteds, dimarker** correctors)
        {
            domain->phi[idx] = phi;

            domain->gradphi[0][idx] = normal[0];
            domain->gradphi[1][idx] = normal[1];
            domain->gradphi[2][idx] = normal[2];

            correcteds[idx] = true;
            correctors[idx] = this;
        }

        bool can_correct(int i, int j, float phinode, eufield::domain* domain, int r)
        {

            float phi2;
            for (int p = 1 - r; p <= r; ++p)
                for (int q = 1 - r; q <= r; ++q)
                {
                    phi2 = domain->phi[domain->grid->operator()(i + p, j + q)];
                    if (phinode * phi2 < 0.) // && (fabs(phinode) > 1e-6 || fabs(phi2) > 1e-6))
                    return true;
                }
            return false;
        }

        void correct_bisector(int i, int j, eufield::domain* domain, bool* correcteds, dimarker** correctors)
        {
            if (i <= 0 || j <= 0 || i >= domain->grid->nx || j >= domain->grid->ny) return;

            /*
             float h = domain->grid->h;
             if (fabs(domain->getval(domain->phi, pos)) > h)
             return;
             */

            int idx = domain->grid->operator ()(i, j);
            float phinode = domain->phi[idx];

            geometry::vector xnode = domain->grid->getpos(i, j);
            float phi = getphi(xnode);

            geometry::vector marker_node = xnode - pos;

            if (!correcteds[idx])
            {
                apply_values(idx, phi, domain, correcteds, correctors);
                return;
            }

            float phi12 = correctors[idx]->getphi(pos);
            float phi21 = getphi(correctors[idx]->pos);

            float stronger = fabs(phi12) > fabs(phi21) ? phi12 : phi21;
            if ((stronger < 0. && phinode < phi) || (stronger > 0. && phinode > phi)) apply_values(idx, phi, domain, correcteds, correctors);
        }

    public:
        geometry::vector pos, normal;

        dimarker()
        {

        }

        dimarker(geometry::vector apos, geometry::vector anormal)
        {
            pos = apos;
            normal = anormal;
        }

        static dimarker create(eufield::domain* domain, geometry::vector pos)
        {
            geometry::vector normal = domain->getvec(domain->gradphi, pos);
            normal.normalize();
            /*
             int i, j;
             if (!domain->grid->getidx(pos, i, j))
             throw -1;

             int idx = domain->grid->operator ()(i, j);
             int idxi = domain->grid->operator ()(i + 1, j);
             int idxj = domain->grid->operator ()(i, j + 1);
             int idxij = domain->grid->operator ()(i + 1, j + 1);

             float phi = domain->phi[idx];
             float phii = domain->phi[idxi];
             float phij = domain->phi[idxj];
             float phiij = domain->phi[idxij];

             geometry::vector normal(phii + phiij - phi - phij, phij + phiij - phi - phii, 0);
             normal /= domain->grid->h;
             normal.normalize();
             */

            return dimarker(pos, normal);
        }

        static geometry::vector get_omega(eufield::domain* domain, geometry::vector pos, geometry::vector n)
        {
            geometry::vector omega;
            float phix, phiy, a, b, p, q;
            a = domain->getval(domain->gradu[0][0], pos);
            b = domain->getval(domain->gradu[0][1], pos);
            p = domain->getval(domain->gradu[1][0], pos);
            q = domain->getval(domain->gradu[1][1], pos);
            phix = n[0];
            phiy = n[1];
            omega[0] = 0.;
            omega[1] = 0.;
            omega[2] = (p * phiy * phiy - q * phix * phiy + a * phix * phiy - b * phix * phix) / (phix * phix + phiy * phiy);

            return omega;
        }

        void advect_RK4(eufield::domain* domain, float dt)
        {
            geometry::vector k1, k2, k3, k4, a, b, c, d, nn, pp;
            a = get_omega(domain, pos, normal);
            nn = normal + geometry::cross(a, normal) * dt / 2.;
            nn.normalize();
            k1 = domain->getvec(domain->u, pos);
            pp = pos + k1 * dt / 2;

            b = get_omega(domain, pp, nn);
            nn = normal + geometry::cross(b, normal) * dt / 2.;
            nn.normalize();
            k2 = domain->getvec(domain->u, pp);
            pp = pos + k2 * dt / 2;

            c = get_omega(domain, pp, nn);
            nn = normal + geometry::cross(c, normal) * dt;
            nn.normalize();
            k3 = domain->getvec(domain->u, pp);
            pp = pos + k3 * dt;

            d = get_omega(domain, pp, nn);
            k4 = domain->getvec(domain->u, pp);

            pos += (k1 + 2. * k2 + 2. * k3 + k4) * (dt / 6.);
            normal += geometry::cross((a + 2. * b + 2. * c + d) / 6., normal) * dt;
            normal.normalize();
        }

        void advect_RK2(eufield::domain* domain, float dt)
        {
            geometry::vector k1, k2, a, b, nn, pp;
            a = get_omega(domain, pos, normal);
            nn = normal + geometry::cross(a, normal) * dt / 2.;
            nn.normalize();
            k1 = domain->getvec(domain->u, pos);
            pp = pos + k1 * dt / 2;

            b = get_omega(domain, pp, nn);
            k2 = domain->getvec(domain->u, pp);

            pos += k2 * dt;
            normal += geometry::cross(b, normal) * dt;
            normal.normalize();
        }

        float getphi(geometry::vector v)
        {
            return geometry::dot(normal, v - pos);
        }

        void correct(eufield::domain* domain, bool* correcteds, dimarker** correctors)
        {
            int i, j;
            if (!domain->grid->getidx(pos, i, j)) throw -1;

            /*
             correct_point(i, j, domain, dists);
             correct_point(i + 1, j, domain, dists);
             correct_point(i + 1, j + 1, domain, dists);
             correct_point(i, j + 1, domain, dists);
             */

            int r = 1;

            if (!can_correct(i, j, domain->phi[domain->grid->operator ()(i, j)], domain, r)) return;

            for (int p = 1 - r; p <= r; ++p)
                for (int q = 1 - r; q <= r; ++q)
                    correct_bisector(i + p, j + q, domain, correcteds, correctors);

        }
};
}

#endif /* DIMARKER_H_INCLUDED */
