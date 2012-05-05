#ifndef ZALESAK_H_INCLUDED
#define ZALESAK_H_INCLUDED

#include "vector.h"
#include "grid.h"
#include "domain.h"
#include "fmm.h"
#include "sl.h"
#include "dimarkerset.h"

#include <iostream>
#include <math.h>
#include <sstream>

class zalesak
{
	static void make_sphere(eufield::grid *grid, float* phi, geometry::vector pos, float radius)
	{
		for (int i = 0; i < grid->nx; ++i)
			for (int j = 0; j < grid->ny; ++j)
			{
				int idx = (*grid)(i, j);

				float left_wall = pos[0] - radius / 6;
				float right_wall = pos[0] + radius / 6;
				float top_wall = pos[1] + radius / 1.5;
				float bottom_wall = pos[1] - radius;

				float x = i * grid->h + grid->startx;
				float y = j * grid->h + grid->starty;

				float val1 = (pos - geometry::vector(x, y, 0)).length() - radius;

				float val2 = INFINITY;

				if (bottom_wall <= y && y <= top_wall && left_wall <= x && x <= right_wall)//inside void
				{
					float top = top_wall - y;
					float bottom = y - bottom_wall;
					float left = x - left_wall;
					float right = right_wall - x;
					if (fabs(val2) > fabs(top))
						val2 = top;
					if (fabs(val2) > fabs(bottom))
						val2 = bottom;
					if (fabs(val2) > fabs(left))
						val2 = left;
					if (fabs(val2) > fabs(right))
						val2 = right;
				}
				else if (left_wall <= x && x <= right_wall)
				{
					if (fabs(val2) > fabs(y - bottom_wall))
						val2 = y - bottom_wall;
					if (fabs(val2) > fabs(y - top_wall))
						val2 = top_wall - y;
				}
				else if (bottom_wall <= y && y <= top_wall)
				{
					if (fabs(val2) > fabs(x - left_wall))
						val2 = x - left_wall;
					if (fabs(val2) > fabs(x - right_wall))
						val2 = right_wall - x;
				}
				else
				{
					float ul = sqrt((x - left_wall) * (x - left_wall) + (y - top_wall) * (y - top_wall));
					float ur = sqrt((x - right_wall) * (x - right_wall) + (y - top_wall) * (y - top_wall));
					float bl = sqrt((x - left_wall) * (x - left_wall) + (y - bottom_wall) * (y - bottom_wall));
					float br = sqrt((x - right_wall) * (x - right_wall) + (y - bottom_wall) * (y - bottom_wall));

					val2 = -std::min(std::min(std::min(ul, ur), bl), br);
				}
				//val2 *= grid->h;

				phi[idx] = std::max(val1, val2);
			}
	}

	static void init_domain(eufield::grid *&grid, eufield::domain *&dom, lfield::dimarkerset *&dms, float omega)
	{
		int n = 50;
		int nx = n + 2, ny = n + 2;

		grid = new eufield::grid(nx, ny, nx * ny);
		grid->startx = grid->starty = 0;
		grid->h = 100. / n;
		for (int i = 0, k = 0; i < nx; ++i)
			for (int j = 0; j < ny; ++j)
			{
				grid->idxs[i][j] = k++;
				if (i == 0 || i == nx - 1 || j == 0 || j == ny - 1)
					grid->bcs[i][j] = 1;
			}
		geometry::vector center(grid->h * (nx - 1.) / 2., grid->h * (ny - 1.) / 2., 0);

		dom = new eufield::domain(grid);
		dms = new lfield::dimarkerset(dom);

		make_sphere(grid, dom->phi, center + geometry::vector(0, 25, 0), 15);
		ls::fmm::reinit(dom->grid, dom->phi);

		for (int i = 0; i < grid->nx; ++i)
			for (int j = 0; j < grid->ny; ++j)
			{
				int idx = (*grid)(i, j);
				if (idx < 0)
					continue;

				float x = i * grid->h + grid->startx;
				float y = j * grid->h + grid->starty;

				dom->u[0][idx] = omega * (center.y() - y);
				dom->u[1][idx] = omega * (x - center.x());
			}

		dom->calcgradphi();
		dom->calcgradu();
		dom->calcomega();
		dom->calcgrad2phi();

		dms->init_markers_intersects();
	}

public:
	static void run()
	{
		float PI = 3.14159265;
		float omega = PI / 314.;
		float dt = (2. * PI / omega) / 628.;

		/////////////////////////////// initializing
		eufield::grid *grid;
		eufield::domain *dom;
		lfield::dimarkerset *dms;
		init_domain(grid, dom, dms, omega);

		////////////////////////////// ls
		std::ostringstream filename;

		filename.str("");
		filename << "out/ls" << 0 << ".vtk";
		dom->write(filename.str());
		filename.str("");
		filename << "out/dm" << 0 << ".vtk";
		dms->write(filename.str());

		for (int i = 1; i <= 2 * 628; ++i)
		{
			ls::sl::advect(dom, dom->grid, dom->phi, dt);

			dms->advect(dt);
			dms->correct();

			ls::fmm::reinit(dom->grid, dom->phi);

			//dms->add_markers();

			if (true)//i % 157 == 0)
			{
				filename.str("");
				filename << "out/ls" << i << ".vtk";
				dom->write(filename.str());

				filename.str("");
				filename << "out/dm" << i << ".vtk";
				dms->write(filename.str());
			}

			std::cout << "time step " << i << " done." << std::endl << std::flush;
		}
	}
};

#endif /* ZALESAK_H_INCLUDED */
