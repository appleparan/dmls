#ifndef VORTEX_H_INCLUDED
#define VORTEX_H_INCLUDED

#include "vector.h"
#include "grid.h"
#include "domain.h"
#include "fmm.h"
#include "sl.h"
#include "dimarkerset.h"

#include <iostream>
#include <math.h>
#include <sstream>

class vortex
{
	static void make_sphere(eufield::grid *grid, float* phi, geometry::vector pos, float radius)
	{
		for (int i = 0; i < grid->nx; ++i)
			for (int j = 0; j < grid->ny; ++j)
			{
				int idx = (*grid)(i, j);

				float x = i * grid->h + grid->startx;
				float y = j * grid->h + grid->starty;

				phi[idx] = (pos - geometry::vector(x, y, 0)).length() - radius;
			}
	}

	static void init_domain(eufield::grid *&grid, eufield::domain *&dom, lfield::dimarkerset *&dms, float& dt)
	{
		float PI = 3.14159265;

		int n = 128;
		int nx = n, ny = n;

		grid = new eufield::grid(nx, ny, nx * ny);
		grid->startx = grid->starty = 0;
		grid->h = 1. / (n - 1.);
		for (int i = 0, k = 0; i < nx; ++i)
			for (int j = 0; j < ny; ++j)
			{
				grid->idxs[i][j] = k++;
				if (i == 0 || i == nx - 1 || j == 0 || j == ny - 1)
					grid->bcs[i][j] = 1;
			}
		geometry::vector center(grid->h * nx / 2., grid->h * ny / 2., 0);

		dom = new eufield::domain(grid);
		dms = new lfield::dimarkerset(dom);

		make_sphere(grid, dom->phi, center + geometry::vector(0, 0.25, 0), 0.15);
		ls::fmm::reinit(dom->grid, dom->phi);

		float maxu = -1;
		for (int i = 0; i < grid->nx; ++i)
			for (int j = 0; j < grid->ny; ++j)
			{
				int idx = (*grid)(i, j);
				if (idx < 0)
					continue;

				float x = i * grid->h + grid->startx;
				float y = j * grid->h + grid->starty;

				float sx = sin(PI * x), cx = cos(PI * x), sy = sin(PI * y), cy = cos(PI * y);

				dom->u[0][idx] = -2. * sx * sx * sy * cy;
				dom->u[1][idx] = +2. * sx * cx * sy * sy;

				dom->dummy[idx] = sx * sx * sy * sy / PI;

				float u = sqrt(dom->u[0][idx] * dom->u[0][idx] + dom->u[1][idx] * dom->u[1][idx]);
				if (u > maxu)
					maxu = u;
			}

		dom->calcgradphi();
		dom->calcgradu();
		dom->calcomega();
		dom->calcgrad2phi();

		dms->init_markers_intersects();

		dt = grid->h / maxu;
	}

	static void reverse_domain(eufield::grid *grid, eufield::domain *dom, float tT)
	{
		float PI = 3.14159265;

		for (int i = 0; i < grid->nx; ++i)
			for (int j = 0; j < grid->ny; ++j)
			{
				int idx = (*grid)(i, j);

				geometry::vector v = grid->getpos(i, j);
				float x = v[0];
				float y = v[1];

				float sx = sin(PI * x), cx = cos(PI * x), sy = sin(PI * y), cy = cos(PI * y);
				float ct = cos(PI * tT);

				dom->u[0][idx] = -2. * sx * sx * sy * cy * ct;
				dom->u[1][idx] = +2. * sx * cx * sy * sy * ct;
			}
	}

public:

	static void run()
	{
		float dt;
		float T = 8;

		/////////////////////////////// initializing
		eufield::grid *grid;
		eufield::domain *dom;
		lfield::dimarkerset *dms;
		init_domain(grid, dom, dms, dt);
		//dt /= 10;
		int nsteps = ceil(T / dt / 2.);
		float t = -dt;

		std::cout << "h: " << dom->grid->h << std::endl << std::flush;
		std::cout << "dt: " << dt << std::endl << std::flush;
		std::cout << "T: " << T << std::endl << std::flush;
		std::cout << "nsteps: " << nsteps << std::endl << std::flush;

		////////////////////////////// ls
		std::ostringstream filename;

		filename.str("");
		filename << "out/ls" << 0 << ".vtk";
		dom->write(filename.str());
		filename.str("");
		filename << "out/dm" << 0 << ".vtk";
		dms->write(filename.str());

		int i = 0;
		while (t / T <= 0.5)
		{
			t += dt;

			++i;
			reverse_domain(grid, dom, t / T);
			dom->calcgradu();

			ls::sl::advect(dom, dom->grid, dom->phi, dt);

			dms->advect(dt);
			dms->correct();

			//ls::fmm::reinit(dom->grid, dom->phi);

			//if (t / T < 0.5)
			dms->add_markers();
			//if (t / T > 0.5)
			//dms->delete_markers();

			if (i % 1 == 0)// == 0 and i>479)
			//if (fabs(t - floor(t)) < dt)
			{
				std::cout << "*";

				filename.str("");
				filename << "out/ls" << i << ".vtk";
				dom->write(filename.str());

				filename.str("");
				filename << "out/dm" << i << ".vtk";
				dms->write(filename.str());
			}

			std::cout << "time step " << i << " done.  t/T= " <<t/T<< std::endl << std::flush;
		}

		std::cout << t << std::endl;
	}
};

#endif /* VORTEX_H_INCLUDED */
