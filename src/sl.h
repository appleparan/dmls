#ifndef SL_H_INCLUDED
#define SL_H_INCLUDED

#include <algorithm>

namespace ls {
class sl
{
    public:
        static void advect(eufield::domain *d, eufield::grid *g, float* phi, float dt)
        {
            float* result = new float[g->n];
            //std::fill_n(result,g->n,100000);
            std::copy(phi, phi + g->n, result);

            for (int i = 0; i < g->nx; ++i)
                for (int j = 0; j < g->ny; ++j)
                {
                    if (g->bcs[i][j] != 0) continue;

                    int idx = (*g)(i, j);

                    geometry::vector vel = d->getvec(d->u, idx);
                    geometry::vector pos = g->getpos(i, j);
                    geometry::vector newpos = pos - vel * dt;

                    if (g->isinside(newpos)) result[idx] = d->getval(d->phi, newpos);
                }

            std::copy(result, result + g->n, phi);
            delete[] result;
        }
};
}

#endif /* SL_H_INCLUDED */
