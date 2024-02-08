///////////////////////////////////////////////////////////////////////////////////////////////
//                                VOROCRUST-MESHING 1.0                                      //
// Copyright 2022 National Technology & Engineering Solutions of Sandia, LLC (NTESS).        //
// Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain  //
// rights in this software.                                                                  //
//                                                                                           //
// Redistribution and use in source and binary forms, with or without modification, are      //
// permitted provided that the following conditions are met:                                 //  
//                                                                                           //
// 1. Redistributions of source code must retain the above copyright notice, this list of    //
// conditions and the following disclaimer.                                                  //
//                                                                                           //
// 2. Redistributions in binary form must reproduce the above copyright notice, this list    //
// of conditions and the following disclaimer in the // documentation and/or other materials //
// provided with the distribution.                                                           //
//                                                                                           //
// 3. Neither the name of the copyright holder nor the names of its contributors may be      //
// used to endorse or promote products derived from this software without specific prior     //
// written permission.                                                                       //
//-------------------------------------------------------------------------------------------//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY       //
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF   //
// MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE//
// COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, //
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        //
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)    //
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR  //
// TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS        //
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                              //
///////////////////////////////////////////////////////////////////////////////////////////////
//                                     Author                                                //
//                                Mohamed S. Ebeida                                          //
//                                msebeid@sandia.gov                                         //
///////////////////////////////////////////////////////////////////////////////////////////////
//  MeshingGeometricalMethods.cpp                                 Last modified (07/11/2022) //
///////////////////////////////////////////////////////////////////////////////////////////////

#include "MeshingGeometricalMethods.h"

MeshingGeometricalMethods::MeshingGeometricalMethods()
{

}

MeshingGeometricalMethods::~MeshingGeometricalMethods()
{

}

double MeshingGeometricalMethods::get_point_angle(double x, double y)
{
	#pragma region Get Point Angle:
	if (fabs(x) < 1E-10 && fabs(y) < 1E-10) return DBL_MAX; // this is a singularity

	double angle(0.0);
	if (fabs(x) > fabs(y)) angle = atan(fabs(y) / fabs(x));
	else                   angle = 0.5 * PI - atan(fabs(x) / fabs(y));

	if (x > -1E-10 && y > -1E-10) return angle;
	else if (x > -1E-10)          return 2 * PI - angle;
	else if (y > -1E-10)          return PI - angle;
	else                          return PI + angle;
	return 0.0;
	#pragma endregion
}

double MeshingGeometricalMethods::distance(size_t num_dim, double* x, double* y)
{
	#pragma region Distance between two points:
	double dst = 0.0;
	for (size_t idim = 0; idim < num_dim; idim++)
	{
		double dx = y[idim] - x[idim];
		dst += dx * dx;
	}
	dst = sqrt(dst);
	return dst;
	#pragma endregion
}

double MeshingGeometricalMethods::cos_angle(size_t num_dim, double* x, double* y, double* z)
{
	#pragma region Distance between two points:
	double hxy = distance(num_dim, x, y);
	double hyz = distance(num_dim, y, z);
	if (hxy < 1E-10 || hyz < 1E-10) return 1.0;
	double dot(0.0);
	for (size_t idim = 0; idim < num_dim; idim++)
	{
		double dxy = x[idim] - y[idim];
		double dzy = z[idim] - y[idim];
		dot += dxy * dzy;
	}
	return dot /= (hxy * hyz);
	#pragma endregion
}

void MeshingGeometricalMethods::get_power_vertex(size_t num_dim, double* co, double ro, double* c1, double r1, double* pv)
{
	#pragma region Get Power Vertex:
	double h(0.0);
	for (size_t idim = 0; idim < num_dim; idim++)
	{
		double dx = c1[idim] - co[idim];
		h += dx * dx;
	}
	h = sqrt(h);

	double ho = (h * h + ro * ro - r1 * r1) / (2 * h);

	for (size_t idim = 0; idim < num_dim; idim++) pv[idim] = co[idim] + ho * (c1[idim] - co[idim]) / h;

	#pragma endregion
}

bool MeshingGeometricalMethods::get_power_vertex(size_t num_dim, size_t num_points, double** centers, double* radii, double* pv)
{
	#pragma region Power Vertex of n spheres:

	double** c = new double*[num_points]; double* r = new double[num_points];
	for (size_t i = 0; i < num_points; i++)
	{
		c[i] = centers[i]; r[i] = radii[i];
	}

	double** basis = new double*[num_dim];
	for (size_t ibasis = 0; ibasis < num_dim; ibasis++)
	{
		basis[ibasis] = new double[num_dim];
		for (size_t idim = 0; idim < num_dim; idim++) basis[ibasis][idim] = 0.0;
		basis[ibasis][ibasis] = 1.0;
	}

	double* xst = new double[num_dim];  double* dir = new double[num_dim];
	double* dart = new double[num_dim]; double* vect = new double[num_dim];

	for (size_t idim = 0; idim < num_dim; idim++) xst[idim] = c[0][idim];

	double closest_distance = 0.0; bool trimming_failed(false);
	for (size_t ispoke = 0; ispoke < num_points - 1; ispoke++)
	{
		// throw a random spoke
		_rsampler.sample_uniformly_from_unit_sphere(dart, num_dim - ispoke);

		for (size_t idim = 0; idim < num_dim; idim++) dir[idim] = 0.0;

		for (size_t jbasis = ispoke; jbasis < num_dim; jbasis++)
		{
			for (size_t idim = 0; idim < num_dim; idim++) dir[idim] += dart[jbasis - ispoke] * basis[jbasis][idim];
		}

		// trim using all remaining points
		double alpha_min(DBL_MAX); size_t neighbor_index(0);
		for (size_t ipoint = ispoke + 1; ipoint < num_points; ipoint++)
		{
			double* qv = new double[num_dim]; double* mv = new double[num_dim];
			get_power_vertex(num_dim, c[0], r[0], c[ipoint], r[ipoint], qv);
			for (size_t idim = 0; idim < num_dim; idim++) mv[idim] = c[ipoint][idim] - c[0][idim];

			double dot_1(0.0), dot_2(0.0);
			for (size_t idim = 0; idim < num_dim; idim++)
			{
				dot_1 += (qv[idim] - xst[idim]) * mv[idim];
				dot_2 += dir[idim] * mv[idim];
			}
			delete[] qv; delete[] mv;
			if (fabs(dot_2) < 1E-10) continue; // spoke is parallel to Power Hyperplane

			double alpha = dot_1 / dot_2;
			if (fabs(alpha) < fabs(alpha_min))
			{
				alpha_min = alpha; neighbor_index = ipoint;
			}
		}
		if (neighbor_index == 0)
		{
			for (size_t ibasis = 0; ibasis < num_dim; ibasis++) delete[] basis[ibasis];
			delete[] basis; delete[] xst; delete[] dir; delete[] dart; delete[] vect;
			delete[] c; delete[] r;
			return false;
		}

		// update spoke start
		for (size_t idim = 0; idim < num_dim; idim++) xst[idim] += alpha_min * dir[idim];

		double* tmp = c[neighbor_index]; c[neighbor_index] = c[ispoke + 1]; c[ispoke + 1] = tmp;
		double rtmp = r[neighbor_index]; r[neighbor_index] = r[ispoke + 1]; r[ispoke + 1] = rtmp;

		// update basis
		for (size_t idim = 0; idim < num_dim; idim++) vect[idim] = c[ispoke + 1][idim] - c[0][idim];

		double norm;
		get_normal_component(num_dim, ispoke, basis, vect, norm);
		for (size_t idim = 0; idim < num_dim; idim++) basis[ispoke][idim] = vect[idim];

		// update remaining basis
		for (size_t jbasis = ispoke + 1; jbasis < num_dim; jbasis++)
		{
			for (size_t idim = 0; idim < num_dim; idim++)
			{
				for (size_t jdim = 0; jdim < num_dim; jdim++) vect[jdim] = 0.0;
				vect[idim] = 1.0;
				get_normal_component(num_dim, jbasis, basis, vect, norm);
				if (norm > 0.1) break;
			}
			for (size_t idim = 0; idim < num_dim; idim++) basis[jbasis][idim] = vect[idim];
		}
	}

	// project the Voronoi vertex to the space of the points
	for (size_t idim = 0; idim < num_dim; idim++) vect[idim] = xst[idim] - c[0][idim];

	for (size_t jbasis = num_points - 1; jbasis < num_dim; jbasis++)
	{
		double dot(0.0);
		for (size_t idim = 0; idim < num_dim; idim++) dot += vect[idim] * basis[jbasis][idim];
		for (size_t idim = 0; idim < num_dim; idim++) vect[idim] -= dot * basis[jbasis][idim];
	}

	for (size_t idim = 0; idim < num_dim; idim++) pv[idim] = c[0][idim] + vect[idim];

	for (size_t ibasis = 0; ibasis < num_dim; ibasis++) delete[] basis[ibasis];
	delete[] basis; delete[] xst; delete[] dir; delete[] dart; delete[] vect;
	delete[] c; delete[] r;
	return true;
	#pragma endregion
}

int MeshingGeometricalMethods::get_normal_component(size_t num_dim, size_t num_basis, double** basis, double* vect, double &norm)
{
	#pragma region Get Normal component to some basis:

	double* comp = new double[num_basis];

	// project point to current basis
	for (size_t ibasis = 0; ibasis < num_basis; ibasis++)
	{
		comp[ibasis] = 0.0;
		for (size_t idim = 0; idim < num_dim; idim++) comp[ibasis] += vect[idim] * basis[ibasis][idim];
	}

	// get vector component orthogonal to current basis
	for (size_t ibasis = 0; ibasis < num_basis; ibasis++)
	{
		for (size_t idim = 0; idim < num_dim; idim++) vect[idim] -= comp[ibasis] * basis[ibasis][idim];
	}

	delete[] comp;

	norm = 0.0;
	for (size_t idim = 0; idim < num_dim; idim++) norm += vect[idim] * vect[idim];

	if (fabs(norm) < 1E-10) return 1;

	norm = 1.0 / sqrt(norm);
	for (size_t idim = 0; idim < num_dim; idim++) vect[idim] *= norm;

	return 0;
	#pragma endregion
}

double MeshingGeometricalMethods::dot_product(size_t num_dim, double* v1, double* v2)
{
	#pragma region dot product:
	double dot = 0.0;
	for (size_t idim = 0; idim < num_dim; idim++)
	{
		dot += v1[idim] * v2[idim];
	}
	return dot;
	#pragma endregion
}

bool MeshingGeometricalMethods::normalize_vector(size_t num_dim, double* vec)
{
	#pragma region Normalize Vecotr:
	double norm(0.0);
	for (size_t idim = 0; idim < num_dim; idim++)
	{
		norm += vec[idim] * vec[idim];
	}
	norm = sqrt(norm);
	if (norm < 1E-10) return false;
	for (size_t idim = 0; idim < num_dim; idim++) vec[idim] /= norm;
	return true;
	#pragma endregion
}

bool MeshingGeometricalMethods::overlapping_spheres(size_t num_dim, double* sphere_i, double* sphere_j)
{
	#pragma region Overlapping spheres check:
	double h = distance(num_dim, sphere_i, sphere_j);
	double ri = sphere_i[num_dim];
	double rj = sphere_j[num_dim];
	if (h < ri + rj + 1E-10) return true;
	return false;
	#pragma endregion
}

double MeshingGeometricalMethods::get_3d_triangle_area(double** corners)
{
	#pragma region Get Area of a triangle:
	double* a = new double[3]; double* b = new double[3];
	for (size_t idim = 0; idim < 3; idim++)
	{
		a[idim] = corners[1][idim] - corners[0][idim];
		b[idim] = corners[2][idim] - corners[0][idim];
	}
	double nx = a[1] * b[2] - a[2] * b[1];
	double ny = a[2] * b[0] - a[0] * b[2];
	double nz = a[0] * b[1] - a[1] * b[0];
	double norm = sqrt(nx * nx + ny * ny + nz * nz);
	delete[] a; delete[] b;
	return 0.5 * norm;
	#pragma endregion
}

bool MeshingGeometricalMethods::get_3d_triangle_normal(double** corners, double* normal)
{
	#pragma region Get Face Area and Normal:
	double* a = new double[3]; double* b = new double[3];
	for (size_t idim = 0; idim < 3; idim++)
	{
		a[idim] = corners[1][idim] - corners[0][idim];
		b[idim] = corners[2][idim] - corners[0][idim];
	}
	if (!normalize_vector(3, a) || !normalize_vector(3, b))
	{
		for (size_t idim = 0; idim < 3; idim++) normal[idim] = 0.0;
		delete[] a; delete[] b;
		return false; // a degenerate face due to edge length
	}
	normal[0] = a[1] * b[2] - a[2] * b[1];
	normal[1] = a[2] * b[0] - a[0] * b[2];
	normal[2] = a[0] * b[1] - a[1] * b[0];
	double norm = normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2];
	norm = sqrt(norm);
	if (fabs(norm) < 1E-10)
	{
		delete[] a; delete[] b;
		return false; // a degenerate face due to angle
	}
	for (size_t idim = 0; idim < 3; idim++) normal[idim] /= norm;
	delete[] a; delete[] b;
	return true;
	#pragma endregion
}

int MeshingGeometricalMethods::project_to_3d_triangle(double* p, double** corners, double* q, double &proj_dist)
{
	#pragma region Project to 3d Triangle:
	double* n = new double[3];
	get_3d_triangle_normal(corners, n);
	double* a = new double[3]; double* b = new double[3];
	for (size_t idim = 0; idim < 3; idim++) a[idim] = corners[0][idim] - p[idim];
	double dot = dot_product(3, a, n);
	for (size_t idim = 0; idim < 3; idim++) q[idim] = p[idim] + dot * n[idim];

	bool in_triangle(true);
	double* m = new double[3];
	for (size_t i = 0; i < 3; i++)
	{
		size_t ip = i + 1; if (ip == 3) ip = 0;
		size_t im = 2; if (i > 0) im = i - 1;
		double* tmp = corners[im];
		corners[im] = q;
		if (!get_3d_triangle_normal(corners, m))
		{
			corners[im] = tmp;
			continue;
		}
		dot = dot_product(3, m, n);
		corners[im] = tmp;
		if (dot < -1E-10)
		{
			in_triangle = false; break;
		}
	}
	if (in_triangle)
	{
		proj_dist = distance(3, p, q);
		delete[] n; delete[] m; delete[] a; delete[] b;
		return 0;
	}

	// project to edges/corners
	double hclosest(DBL_MAX);
	for (size_t i = 0; i < 3; i++)
	{
		size_t ip = i + 1; if (ip == 3) ip = 0;
		size_t im = 2; if (i > 0) im = i - 1;
		for (size_t idim = 0; idim < 3; idim++)
		{
			a[idim] = p[idim] - corners[i][idim];
			b[idim] = corners[ip][idim] - corners[i][idim];
		}
		dot = dot_product(3, a, b);

		double b_norm(0.0);
		for (size_t idim = 0; idim < 3; idim++) b_norm += b[idim] * b[idim];

		double u = dot / b_norm;
		if (u < 0.0) u = 0.0;
		if (u > 1.0) u = 1.0;
		for (size_t idim = 0; idim < 3; idim++) a[idim] = corners[i][idim] + u * (corners[ip][idim] - corners[i][idim]);

		double h = distance(3, a, p);
		if (h < hclosest)
		{
			hclosest =  h;
			for (size_t idim = 0; idim < 3; idim++) q[idim] = a[idim];
		}
	}
	proj_dist = distance(3, p, q);
	delete[] n; delete[] m; delete[] a; delete[] b;
	return 0;
	#pragma endregion
}

int MeshingGeometricalMethods::project_to_3d_line_segment(double* p, double* xo, double* xn, double* q, double& proj_dist)
{
	#pragma region Project to 3d line segement:
	double* a = new double[3];
	double* b = new double[3];
	for (size_t idim = 0; idim < 3; idim++)
	{
		a[idim] = p[idim] - xo[idim];
		b[idim] = xn[idim] - xo[idim];
	}
	double dot = dot_product(3, a, b);

	double b_norm(0.0);
	for (size_t idim = 0; idim < 3; idim++) b_norm += b[idim] * b[idim];

	double u = dot / b_norm;
	if (u < 0.0) u = 0.0;
	if (u > 1.0) u = 1.0;
	for (size_t idim = 0; idim < 3; idim++) q[idim] = xo[idim] + u * (xn[idim] - xo[idim]);
	proj_dist = distance(3, p, q);

	delete[] a; delete[] b;
	return 0;
	#pragma endregion
}

int MeshingGeometricalMethods::project_to_3d_line(double* p, double* xo, double* edir, double* q, double &proj_dist)
{
	#pragma region Project to 3d line segement:
	double* a = new double[3];
	for (size_t idim = 0; idim < 3; idim++) a[idim] = p[idim] - xo[idim];

	double u = dot_product(3, a, edir);
	for (size_t idim = 0; idim < 3; idim++) q[idim] = xo[idim] + u * edir[idim];
	proj_dist = distance(3, p, q);

	delete[] a;
	return 0;
	#pragma endregion
}
