#ifndef FMM_H_INCLUDED
#define FMM_H_INCLUDED

#include "grid.h"
#include <algorithm>
#include <math.h>
#include <map>

namespace ls
{
class fmm
{
	eufield::grid* grid;
	enum node_conds
	{
		none, alive, accepted
	};
	node_conds* conds;
	float *result, *data;
	float h2;

	typedef std::multimap<float, eufield::molecule> phi_idx_map;
	phi_idx_map *alive_nodes;
	phi_idx_map::iterator *alive_iters;

	fmm(eufield::grid* grid, float *data, float *result)
	{
		this->grid = grid;
		this->data = data;
		this->result = result;
		h2 = grid->h * grid->h;
		conds = new node_conds[grid->n];
		std::fill_n(conds, grid->n, none);
		alive_nodes = new phi_idx_map();
		alive_iters = new phi_idx_map::iterator[grid->n];
	}

	~fmm()
	{
		delete[] conds;
		delete alive_nodes;
		delete[] alive_iters;
	}

	static float select(float x, float y)
	{
		if (x < 0)
			return y;
		if (y < 0)
			return x;
		return std::min(x, y);
	}

	void init_interface_adjacents(bool anchored)
	{
		for (int i = 0; i < grid->nx; ++i)
			for (int j = 0; j < grid->ny; ++j)
			{
				if (grid->bcs[i][j] != 0)
					continue;

				eufield::molecule m = grid->get_molecule(i, j);

				float val = data[m.p];
				if (val < 0)
					continue;

				float xm_val = -1, xp_val = -1, yp_val = -1, ym_val = -1;

				if (data[m.xm] < 0)
					xm_val = val / (val - data[m.xm]) * grid->h;
				if (data[m.xp] < 0)
					xp_val = val / (val - data[m.xp]) * grid->h;
				if (data[m.yp] < 0)
					yp_val = val / (val - data[m.yp]) * grid->h;
				if (data[m.ym] < 0)
					ym_val = val / (val - data[m.ym]) * grid->h;

				if (xm_val >= 0 || xp_val >= 0 || yp_val >= 0 || ym_val >= 0)
				{
					if (anchored)
						result[m.p] = data[m.p];
					else
						result[m.p] = select(select(xm_val, xp_val), select(ym_val, yp_val));
					accept(m);
				}
			}
	}

	void change_sign()
	{
		int n = grid->n;
		for (int i = 0; i < n; ++i)
		{
			data[i] *= -1.;
			result[i] *= -1.;
		}
	}

	void reinit_side()
	{
		while (!alive_nodes->empty())
		{
			phi_idx_map::iterator top = alive_nodes->begin();
			accept(top->second);
		}
	}

	void accept(const eufield::molecule m)
	{
		if (conds[m.p] == alive)
			alive_nodes->erase(alive_iters[m.p]);
		conds[m.p] = accepted;

		if (m.xm >= 0 && data[m.xm] > 0. && conds[m.xm] != accepted)
			update(m.i - 1, m.j, m.xm);
		if (m.xp >= 0 && data[m.xp] > 0. && conds[m.xp] != accepted)
			update(m.i + 1, m.j, m.xp);
		if (m.ym >= 0 && data[m.ym] > 0. && conds[m.ym] != accepted)
			update(m.i, m.j - 1, m.ym);
		if (m.yp >= 0 && data[m.yp] > 0. && conds[m.yp] != accepted)
			update(m.i, m.j + 1, m.yp);
	}

	void update(int i, int j, int idx)
	{
		if (grid->bcs[i][j] != 0)
			return;

		eufield::molecule m;
		if (conds[idx] == alive)
		{
			m = alive_iters[idx]->second;
			alive_nodes->erase(alive_iters[idx]);
		}
		else
		{
			m = grid->get_molecule(i, j);
			conds[m.p] = alive;
		}

		float xp_val = -1, xm_val = -1, yp_val = -1, ym_val = -1;

		if (conds[m.xm] == accepted)
			xm_val = data[m.xm];
		if (conds[m.xp] == accepted)
			xp_val = data[m.xp];
		if (conds[m.ym] == accepted)
			ym_val = data[m.ym];
		if (conds[m.yp] == accepted)
			yp_val = data[m.yp];

		result[m.p] = calc_phi(select(xm_val, xp_val), select(ym_val, yp_val));

		alive_iters[m.p] = alive_nodes->insert(std::make_pair(result[m.p], m));

		if (m.p != idx)
			throw -1;
	}

	float calc_phi(float phi1, float phi2)
	{
		if (phi1 < 0)
			return calc_phi(phi2);
		if (phi2 < 0)
			return calc_phi(phi1);

		float phi_max = std::max(phi1, phi2);

		if (((phi_max - phi1) * (phi_max - phi1) + (phi_max - phi2) * (phi_max - phi2)) / h2 > 1)
		{
			return calc_phi(std::min(phi1, phi2));
		}
		else
			return solve_parabol(2., -2. * (phi1 + phi2), phi1 * phi1 + phi2 * phi2 - h2);

	}

	float calc_phi(float phi)
	{
		return solve_parabol(1., -2. * phi, phi * phi - h2);
	}

	float solve_parabol(float a, float b, float c)
	{
		return (-b + sqrt(b * b - 4. * a * c)) / 2. / a;
	}

public:
	static void reinit(eufield::grid *g, float* phi)
	{
		float *result = new float[g->n];
		std::copy(phi, phi + g->n, result);
		fmm m(g, phi, result);

		m.init_interface_adjacents(true);
		m.reinit_side();
		m.change_sign();

		m.init_interface_adjacents(true);
		m.reinit_side();
		m.change_sign();

		std::copy(result, result + g->n, phi);

		delete[] result;
	}
};
}
#endif /*FMM_H_INCLUDED*/
