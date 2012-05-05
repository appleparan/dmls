#ifndef VECTOR2D_H_INCLUDED
#define VECTOR2D_H_INCLUDED

#include <iostream>
#include <math.h>

namespace geometry
{
class vector
{
public:
	vector()
	{
		init(0., 0., 0.);
	}
	vector(float c)
	{
		init(c, c, c);
	}
	vector(float x, float y, float z)
	{
		init(x, y, z);
	}
	vector(const vector& vi)
	{
		init(vi[0], vi[1], vi[2]);
	}
	vector(float vec[])
	{
		init(vec[0], vec[1], vec[2]);
	}

	void init(float ax, float ay, float az)
	{
		v[0] = ax;
		v[1] = ay;
		v[2] = az;
	}

	float& operator[](int index)
	{
		return v[index];
	}
	const float& operator[](int index) const
	{
		return v[index];
	}

	float& x()
	{
		return v[0];
	}
	float& y()
	{
		return v[1];
	}
	float& z()
	{
		return v[2];
	}

	vector& operator=(const vector &vi)
	{
		v[0] = vi[0];
		v[1] = vi[1];
		v[2] = vi[2];
		return *this;
	}
	vector& operator+=(const vector &vi)
	{
		v[0] += vi[0];
		v[1] += vi[1];
		v[2] += vi[2];
		return *this;
	}
	vector& operator-=(const vector &vi)
	{
		v[0] -= vi[0];
		v[1] -= vi[1];
		v[2] -= vi[2];
		return *this;
	}
	vector& operator*=(const vector &vi)
	{
		v[0] *= vi[0];
		v[1] *= vi[1];
		v[2] *= vi[2];
		return *this;
	}
	vector& operator+=(float c)
	{
		v[0] += c;
		v[1] += c;
		v[2] += c;
		return *this;
	}
	vector& operator-=(float c)
	{
		v[0] -= c;
		v[1] -= c;
		v[2] -= c;
		return *this;
	}
	vector& operator*=(float c)
	{
		v[0] *= c;
		v[1] *= c;
		v[2] *= c;
		return *this;
	}
	vector& operator/=(float c)
	{
		float cInv = 1. / c;
		v[0] *= cInv;
		v[1] *= cInv;
		v[2] *= cInv;
		return *this;
	}

	vector operator+(const vector &vi) const
	{
		return vector(*this) += vi;
	}
	vector operator-(const vector &vi) const
	{
		return vector(*this) -= vi;
	}
	vector operator*(const vector &vi) const
	{
		return vector(*this) *= vi;
	}
	vector operator+(float c) const
	{
		return vector(*this) += c;
	}
	vector operator-(float c) const
	{
		return vector(*this) -= c;
	}
	vector operator*(float c) const
	{
		return vector(*this) *= c;
	}
	vector operator/(float c) const
	{
		return vector(*this) /= c;
	}
	vector operator-() const
	{
		return vector(-v[0], -v[1], -v[2]);
	}

	float length() const
	{
		return sqrt(length2());
	}
	float length2() const
	{
		return v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
	}
	void normalize()
	{
		(*this) /= length();
	}
	vector hat() const
	{
		return (*this) / length();
	}
	//vector cross(const vector& vi) const;
	//float dot(const vector& vi) const;

private:
	float v[3];
};

inline vector operator*(float c, const vector &vi)
{
	return vector(c * vi[0], c * vi[1], c * vi[2]);
}

inline vector operator/(float c, const vector &vi)
{
	float cInv = 1.0 / c;
	return vector(cInv * vi[0], cInv * vi[1], cInv * vi[2]);
}

/*
 inline float vector::dot(const vector& vi) const
 {
 return v[0] * vi[0] + v[1] * vi[1] + v[2] * vi[2];
 }
 */

inline float dot(const vector &v1, const vector &v2)
{
	return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

inline vector cross(const vector &v1, const vector &v2)
{
	return vector(v1[1] * v2[2] - v1[2] * v2[1], -(v1[0] * v2[2] - v1[2] * v2[0]), v1[0] * v2[1] - v1[1] * v2[0]);
}

inline std::ostream& operator<<(std::ostream& out, const vector& v)
{
	return out << v[0] << " " << v[1] << " " << v[2];
}

inline std::istream& operator>>(std::istream& in, vector& v)
{
	return in >> v[0] >> v[1] >> v[2];
}
}
#endif /* VECTOR2D_H_INCLUDED */
