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
// MeshingSpheresMethods.cpp                                      Last modified (07/12/2018) //
///////////////////////////////////////////////////////////////////////////////////////////////

#include "MeshingSpheresMethods.h"

MeshingSpheresMethods::MeshingSpheresMethods()
{	

}

MeshingSpheresMethods::~MeshingSpheresMethods()
{
   
}

bool MeshingSpheresMethods::point_covered(double* point, double* sphere, double alpha_coverage)
{
	double beta(1.0);
	if (alpha_coverage > 1E-10) beta = sqrt(1.0 - alpha_coverage * alpha_coverage);
	double h = _geom.distance(3, point, sphere);
	if (h < (beta + 1E-8) * sphere[3] + 1E-8) return true;
	return false;
}

bool MeshingSpheresMethods::point_covered(double* point, MeshingSmartTree* spheres, double Lip, double alpha_coverage)
{
	#pragma region Point Cover Check:
	if (spheres->get_num_tree_points() == 0) return false;

	double beta(1.0);
	if (alpha_coverage > 1E-10) beta = sqrt(1.0 - alpha_coverage * alpha_coverage);

	size_t closest_sphere; double closest_dst(DBL_MAX);
	spheres->get_closest_tree_point(point, closest_sphere, closest_dst);

	double closest_sphere_radius = spheres->get_tree_point_attrib(closest_sphere, 0);

	double estimated_face_center_radius = closest_sphere_radius + Lip * closest_dst;

	size_t num(0), cap(100);
	size_t* neighbor_spheres = new size_t[cap];
	double R = estimated_face_center_radius / (1.0 - Lip);
	spheres->get_tree_points_in_sphere(point, R, num, neighbor_spheres, cap);

	double* sphere = new double[4];
	for (size_t i = 0; i < num; i++)
	{
		size_t isphere = neighbor_spheres[i];
		spheres->get_tree_point(isphere, 4, sphere);

		double h = _geom.distance(3, point, sphere);
		if (h < (beta + 1E-8) * sphere[3] + 1E-8)
		{
			delete[] sphere; delete[] neighbor_spheres;
			return true;
		}
	}

	delete[] sphere; delete[] neighbor_spheres;
	return false;
	#pragma endregion
}


bool MeshingSpheresMethods::edge_covered(double** edge_corners, MeshingSmartTree* spheres, double Lip, double alpha_coverage)
{
	#pragma region Edge Cover Check:
	if (spheres->get_num_tree_points() == 0) return false;

	double beta = sqrt(1.0 - alpha_coverage * alpha_coverage);
	double* edge_center = new double[3];
	for (size_t idim = 0; idim < 3; idim++)
	{
		edge_center[idim] = 0.5 * (edge_corners[0][idim] + edge_corners[1][idim]);
	}

	size_t closest_sphere; double closest_dst(DBL_MAX);
	spheres->get_closest_tree_point(edge_center, closest_sphere, closest_dst);
	double closest_sphere_radius = spheres->get_tree_point_attrib(closest_sphere, 0);
	double estimated_edge_center_radius = closest_sphere_radius + Lip * closest_dst;

	size_t num(0), cap(100);
	size_t* neighbor_spheres = new size_t[cap];
	double R = estimated_edge_center_radius * (1.0 + Lip) / (1.0 - Lip);
	R += estimated_edge_center_radius;
	spheres->get_tree_points_in_sphere(edge_center, R, num, neighbor_spheres, cap);

	double* sphere_center = new double[3];
	for (size_t i = 0; i < num; i++)
	{
		size_t isphere = neighbor_spheres[i];
		double isphere_radius = spheres->get_tree_point_attrib(isphere, 0);
		spheres->get_tree_point(isphere, sphere_center);

		double ho = _geom.distance(3, edge_corners[0], sphere_center);
		double h1 = _geom.distance(3, edge_corners[1], sphere_center);
		if (ho < (beta + 1E-8)* isphere_radius + 1E-8 && h1 <  (beta + 1E-8)* isphere_radius + 1E-8)
		{
			delete[] edge_center; delete[] sphere_center; delete[] neighbor_spheres;
			return true;
		}
	}
	delete[] edge_center; delete[] sphere_center; delete[] neighbor_spheres;
	return false;
	#pragma endregion
}

bool MeshingSpheresMethods::face_covered(double** face_corners, MeshingSmartTree* spheres, double Lip, double alpha_coverage)
{
	#pragma region Face is Covered:

	if (spheres->get_num_tree_points() == 0) return false;

	double beta = sqrt(1.0 - alpha_coverage * alpha_coverage);

	double* face_center = new double[3];
	for (size_t idim = 0; idim < 3; idim++)
	{
		face_center[idim] = (face_corners[0][idim] + face_corners[1][idim] + face_corners[2][idim]) / 3.0;
	}

	size_t closest_sphere; double closest_dst(DBL_MAX);
	spheres->get_closest_tree_point(face_center, closest_sphere, closest_dst);

	double closest_sphere_radius = spheres->get_tree_point_attrib(closest_sphere, 0);

	double estimated_face_center_radius = closest_sphere_radius + Lip * closest_dst;

	size_t num(0), cap(100);
	size_t* neighbor_spheres = new size_t[cap];
	double R = estimated_face_center_radius * (1 + Lip) / (1.0 - Lip);
	R += estimated_face_center_radius; 
	spheres->get_tree_points_in_sphere(face_center, R, num, neighbor_spheres, cap);

	double* sphere = new double[4];
	for (size_t i = 0; i < num; i++)
	{
		size_t isphere = neighbor_spheres[i];
		spheres->get_tree_point(isphere, 4, sphere);

		double ho = _geom.distance(3, face_corners[0], sphere);
		double h1 = _geom.distance(3, face_corners[1], sphere);
		double h2 = _geom.distance(3, face_corners[2], sphere);
		if (ho < (beta + 1E-8) * sphere[3] + 1E-8 && h1 < (beta + 1E-8) * sphere[3] + 1E-8 && h2 < (beta + 1E-8) * sphere[3] + 1E-8)
		{
			delete[] face_center; delete[] sphere; delete[] neighbor_spheres;
			return true;
		}
	}

	if (false)
	{
		// validate Sphere collection
		size_t num_spheres = spheres->get_num_tree_points();
		for (size_t isphere = 0; isphere < num_spheres; isphere++)
		{
			spheres->get_tree_point(isphere, 4, sphere);

			double ho = _geom.distance(3, face_corners[0], sphere);
			double h1 = _geom.distance(3, face_corners[1], sphere);
			double h2 = _geom.distance(3, face_corners[2], sphere);
			if (ho < (beta + 1E-8) * sphere[3] + 1E-8 && h1 < (beta + 1E-8) * sphere[3] + 1E-8 && h2 < (beta + 1E-8) * sphere[3] + 1E-8)
			{
				delete[] face_center; delete[] sphere; delete[] neighbor_spheres;
				return true;
			}
		}
	}
	delete[] face_center; delete[] sphere; delete[] neighbor_spheres;
	return false;
	#pragma endregion
}


int MeshingSpheresMethods::impose_lipschitz_continuity_brute(MeshingSmartTree* spheres, double Lip, bool &shrunk)
{
	#pragma region Impose Lipschitz Continuity Brute Force:
	size_t num_spheres = spheres->get_num_tree_points();

	double* sphere_i = new double[4];
	double* sphere_j = new double[4];

	shrunk = false;
	while (true)
	{
		bool done = true;
		for (size_t i = 0; i < num_spheres; i++)
		{
			for (size_t j = i + 1; j < num_spheres; j++)
			{
				spheres->get_tree_point(i, 4, sphere_i);
				spheres->get_tree_point(j, 4, sphere_j);
				double hij = _geom.distance(3, sphere_i, sphere_j);
				double rLip = fmin(sphere_i[3], sphere_j[3]) + Lip * hij;
				if (rLip < sphere_i[3])
				{
					done = false;
					shrunk = true;
					spheres->set_tree_point_attrib(i, 0, rLip);
				}
				if (rLip < sphere_j[3])
				{
					done = false;
					shrunk = true;
					spheres->set_tree_point_attrib(j, 0, rLip);
				}
			}
		}
		if (done) break;
	}
	delete[] sphere_i; delete[] sphere_j;
	return 0;
	#pragma endregion
}

int MeshingSpheresMethods::impose_lipschitz_continuity(MeshingSmartTree* spheres, double Lip, bool &sphere_shrunk)
{
	#pragma region Impose Lipschitz Continuity:

	size_t num_spheres = spheres->get_num_tree_points();

	size_t* I = new size_t[num_spheres];
	double* R = new double[num_spheres];
	double* new_radii = new double[num_spheres];
	for (size_t i = 0; i < num_spheres; i++)
	{
		I[i] = i;
		new_radii[i] = spheres->get_tree_point_attrib(i, 0);
		R[i] = -new_radii[i];
	}
	_memo.quicksort(R, I, 0, num_spheres - 1);

	size_t num_iter(0);
	while (true)
	{
		bool done = true;
		for (size_t i = 0; i < num_spheres; i++)
		{
			size_t isphere = I[i];

			double* sphere_i = new double[4];
			spheres->get_tree_point(isphere, 4, sphere_i);

			size_t num(0), cap(100);
			size_t* neighbor_spheres = new size_t[cap];
			spheres->get_tree_points_in_sphere(sphere_i, 2.0 * sphere_i[3], num, neighbor_spheres, cap);

			// adjust radius of current sphere
			for (size_t ii = 0; ii < num; ii++)
			{
				#pragma region Shrink large Spheres:
				size_t sphere_index = neighbor_spheres[ii];
				if (sphere_index == isphere) continue;

				double* sphere_j = new double[4];
				spheres->get_tree_point(sphere_index, 4, sphere_j);

				double h = _geom.distance(3, sphere_i, sphere_j);

				bool overlapped = _geom.overlapping_spheres(3, sphere_i, sphere_j);

				double Lip_radius = sphere_j[3] + Lip * h;
				if (new_radii[isphere] > Lip_radius)
				{
					// Shrink larger sphere to satisfy Global lipschitz constant
					new_radii[isphere] = Lip_radius;
					done = false; sphere_shrunk = true;
				}
				delete[] sphere_j;
				#pragma endregion
			}

			delete[] sphere_i; delete[] neighbor_spheres;
		}

		for (size_t i = 0; i < num_spheres; i++)
		{
			spheres->set_tree_point_attrib(i, 0, new_radii[i]);
		}

		num_iter++;
		if (done) break;
	}
	delete[] I; delete[] R; delete[] new_radii;
	return 0;
	#pragma endregion
}

int MeshingSpheresMethods::get_overlapping_spheres(double* sphere, MeshingSmartTree* spheres, double Lip,
	                                        size_t &num_overlapping_spheres, size_t* &overlapping_spheres)
{
	#pragma region Get Overlapping Spheres:
	if (spheres->get_num_tree_points() == 0)
	{
		num_overlapping_spheres = 0;
		overlapping_spheres = 0;
		return 0;
	}

	size_t iclosest; double hclosest(DBL_MAX);
	spheres->get_closest_tree_point(sphere, iclosest, hclosest);
	double rclosest = spheres->get_tree_point_attrib(iclosest, 0);

	double r = rclosest + Lip * hclosest;
	if (r < sphere[3]) r = sphere[3];

	double R = r * (1 + Lip) / (1.0 - Lip);
	R += r;

	size_t num(0), cap(100);
	size_t* neighbor_spheres = new size_t[cap];
	spheres->get_tree_points_in_sphere(sphere, R, num, neighbor_spheres, cap);

	num_overlapping_spheres = 0;
	for (size_t i = 0; i < num; i++)
	{
		double* neighbor_sphere = new double[4];
		spheres->get_tree_point(neighbor_spheres[i], 4, neighbor_sphere);
		if (_geom.overlapping_spheres(3, sphere, neighbor_sphere)) num_overlapping_spheres++;
		delete[] neighbor_sphere;
	}

	overlapping_spheres = new size_t[num_overlapping_spheres];
	num_overlapping_spheres = 0;
	for (size_t i = 0; i < num; i++)
	{
		double* neighbor_sphere = new double[4];
		spheres->get_tree_point(neighbor_spheres[i], 4, neighbor_sphere);
		if (_geom.overlapping_spheres(3, sphere, neighbor_sphere))
		{
			overlapping_spheres[num_overlapping_spheres] = neighbor_spheres[i];
			num_overlapping_spheres++;
		}
		delete[] neighbor_sphere;
	}
	delete[] neighbor_spheres;
	return 0;
	#pragma endregion
}

int MeshingSpheresMethods::resolve_sliver(size_t num_spheres, double** spheres, bool* fixed, size_t &sliver_sphere_index, double &sliver_sphere_radius)
{
	#pragma region Resolve a sliver by minimal shrinkage:
	double* pv = new double[4];
	double** triplet = new double*[3];
	double* normal = new double[3];
	double* upper_int = new double[3];
	double* lower_int = new double[3];

	sliver_sphere_index = num_spheres; double dr_min(DBL_MAX);
	for (size_t i = 0; i < num_spheres; i++)
	{
		triplet[0] = spheres[i];
		for (size_t j = i + 1; j < num_spheres; j++)
		{
			triplet[1] = spheres[j];
			for (size_t k = j + 1; k < num_spheres; k++)
			{
				triplet[2] = spheres[k];
				get_power_vertex(3, 3, triplet, pv);

				_geom.get_3d_triangle_normal(triplet, normal);

				double h = _geom.distance(3, pv, triplet[0]);
				double v = sqrt(triplet[0][3] * triplet[0][3] - h * h);
				for (size_t idim = 0; idim < 3; idim++)
				{
					upper_int[idim] = pv[idim] + v * normal[idim];
					lower_int[idim] = pv[idim] - v * normal[idim];
				}

				for (size_t l = 0; l < num_spheres; l++)
				{
					if (fixed[l]) continue; // Corner and edge sphere has fixed radius
					if (l == i || l == j || l == k) continue;
					double hup = _geom.distance(3, spheres[l], upper_int);
					double hlo = _geom.distance(3, spheres[l], lower_int);
					double hmin = fmin(hup, hlo);
					if (hmin > spheres[l][3] - 1E-10) continue;

					hmin *= 0.95;

					double dr = spheres[l][3] - hmin;
					if (dr < dr_min)
					{
						dr_min = dr;
						sliver_sphere_index = l;
						sliver_sphere_radius = hmin;
					}
				}
			}
		}
	}

	delete[] pv; delete[] triplet; delete[] normal;
	delete[] upper_int; delete[] lower_int;
	return 0;
	#pragma endregion
}


int MeshingSpheresMethods::get_power_vertex(size_t num_dim, double* sphere_i, double* sphere_j, double* pv)
{
	#pragma region Get Power Vertex:
	double h(0.0);
	for (size_t idim = 0; idim < num_dim; idim++)
	{
		double dx = sphere_j[idim] - sphere_i[idim];
		h += dx * dx;
	}
	h = sqrt(h);

	if (h < 1E-10) return 1;

	double ri = sphere_i[num_dim];
	double rj = sphere_j[num_dim];
	double ho = (h * h + ri * ri - rj * rj) / (2 * h);

	for (size_t idim = 0; idim < num_dim; idim++) pv[idim] = sphere_i[idim] + ho * (sphere_j[idim] - sphere_i[idim]) / h;
	
	return 0;
	#pragma endregion
}


bool MeshingSpheresMethods::get_power_vertex(size_t num_dim, size_t num_spheres, double** spheres, double* pv)
{
	#pragma region Power Vertex of n spheres:

	double** s = new double*[num_spheres];
	for (size_t i = 0; i < num_spheres; i++) s[i] = spheres[i];

	double** basis = new double*[num_dim];
	for (size_t ibasis = 0; ibasis < num_dim; ibasis++)
	{
		basis[ibasis] = new double[num_dim];
		for (size_t idim = 0; idim < num_dim; idim++) basis[ibasis][idim] = 0.0;
		basis[ibasis][ibasis] = 1.0;
	}

	double* xst = new double[num_dim];  double* dir = new double[num_dim];
	double* dart = new double[num_dim]; double* vect = new double[num_dim];

	for (size_t idim = 0; idim < num_dim; idim++) xst[idim] = s[0][idim];

	double closest_distance = 0.0; bool trimming_failed(false);
	for (size_t ispoke = 0; ispoke < num_spheres - 1; ispoke++)
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
		for (size_t ipoint = ispoke + 1; ipoint < num_spheres; ipoint++)
		{
			double* qv = new double[num_dim]; double* mv = new double[num_dim];
			get_power_vertex(num_dim, s[0], s[ipoint], qv);
			for (size_t idim = 0; idim < num_dim; idim++) mv[idim] = s[ipoint][idim] - s[0][idim];

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
			delete[] s; 
			return false;
		}

		// update spoke start
		for (size_t idim = 0; idim < num_dim; idim++) xst[idim] += alpha_min * dir[idim];

		double* tmp = s[neighbor_index]; s[neighbor_index] = s[ispoke + 1]; s[ispoke + 1] = tmp;

		// update basis
		for (size_t idim = 0; idim < num_dim; idim++) vect[idim] = s[ispoke + 1][idim] - s[0][idim];

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
	for (size_t idim = 0; idim < num_dim; idim++) vect[idim] = xst[idim] - s[0][idim];

	for (size_t jbasis = num_spheres - 1; jbasis < num_dim; jbasis++)
	{
		double dot(0.0);
		for (size_t idim = 0; idim < num_dim; idim++) dot += vect[idim] * basis[jbasis][idim];
		for (size_t idim = 0; idim < num_dim; idim++) vect[idim] -= dot * basis[jbasis][idim];
	}

	for (size_t idim = 0; idim < num_dim; idim++) pv[idim] = s[0][idim] + vect[idim];

	for (size_t ibasis = 0; ibasis < num_dim; ibasis++) delete[] basis[ibasis];
	delete[] basis; delete[] xst; delete[] dir; delete[] dart; delete[] vect;
	delete[] s;
	return true;
	#pragma endregion
}

int MeshingSpheresMethods::get_normal_component(size_t num_dim, size_t num_basis, double** basis, double* vect, double &norm)
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
