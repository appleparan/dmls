#ifndef DATA_H_INCLUDED
#define DATA_H_INCLUDED

#include "grid.h"
#include <string>
#include <ostream>

namespace eufield
{

class data
{
protected:
	std::string name;
	eufield::grid* grid;
public:
	bool writable;
	data(std::string aname, eufield::grid* agrid, bool awritable = true)
	{
		name = aname;
		grid = agrid;
		writable = awritable;
	}

	virtual ~data()
	{
	}

	virtual void write(std::ostream& s)=0;
};

class data1: public data
{
	float **d;
public:
	data1(std::string aname, eufield::grid* agrid, float** ad, bool awritable) :
		data(aname, agrid, awritable), d(ad)
	{
		*d = new float[agrid->n];
		std::fill_n(*d, agrid->n, 0);
	}
	virtual ~data1()
	{
		delete[] *d;
	}
	void write(std::ostream& s)
	{
		s << "SCALARS " << name << " float 1" << std::endl << "LOOKUP_TABLE default" << std::endl;
		for (int j = 0; j < grid->ny; ++j)
			for (int i = 0; i < grid->nx; ++i)
				if (grid->bcs[i][j] == 0)
					s << (*d)[(*grid)(i, j)] << std::endl;

	}
};

class data3: public data
{
	float ***d;
public:
	data3(std::string aname, eufield::grid* agrid, float*** ad, bool awritable) :
		data(aname, agrid, awritable), d(ad)
	{
		*d = new float*[3];
		(*d)[0] = new float[agrid->n];
		(*d)[1] = new float[agrid->n];
		(*d)[2] = new float[agrid->n];

		std::fill_n((*d)[0], agrid->n, 0);
		std::fill_n((*d)[1], agrid->n, 0);
		std::fill_n((*d)[2], agrid->n, 0);
	}
	virtual ~data3()
	{
		delete[] *d[0];
		delete[] *d[1];
		delete[] *d[2];
		delete[] *d;
	}
	void write(std::ostream& s)
	{
		s << "VECTORS " << name << " float" << std::endl;
		for (int j = 0; j < grid->ny; ++j)
			for (int i = 0; i < grid->nx; ++i)
				if (grid->bcs[i][j] == 0)
				{
					int idx = (*grid)(i, j);
					s << (*d)[0][idx] << " " << (*d)[1][idx] << " " << (*d)[2][idx] << std::endl;
				}
	}
};
}

#endif /* DATA_H_INCLUDED */
