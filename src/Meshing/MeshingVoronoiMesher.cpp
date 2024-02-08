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
//  MeshingVoronoiMesher.cpp                                      Last modified (11/02/2020) //
///////////////////////////////////////////////////////////////////////////////////////////////

#include "MeshingVoronoiMesher.h"
#include "MeshingRandomSampler.h"

//RMMTMP
#include <chrono>
#include <ctime>
#include <numeric>

#ifdef USE_EXODUS
#include "MeshingExodusIO.h"
#endif

MeshingVoronoiMesher::MeshingVoronoiMesher()
{

}

int MeshingVoronoiMesher::ensure_sharp_edge_spheres_are_Delaunay(size_t& num_seeds, double*& seeds, double*& seeds_sizing, size_t*& seeds_region_id,
	                                                             size_t num_surface_spheres, double* surface_spheres, double* surface_spheres_sizing,
	                                                             size_t num_edge_spheres, double* edge_spheres, double* edge_spheres_sizing,
	                                                             size_t num_corner_spheres, double* corner_spheres, double* corner_spheres_sizing,
	                                                             size_t num_sizing_points, double* sizing_points, double* sizing_value,
	                                                             int num_threads, double Lip, double rmax)
{
	#pragma region Add seeds to Sharp Spheres to ensure they are Delaunay:
	// Seeds tree
	
	size_t  seeds_tree_origin(0);
	size_t* seeds_tree_right = new size_t[num_seeds];
	size_t* seeds_tree_left = new size_t[num_seeds];
	build_balanced_kd_tree(num_seeds, 3, seeds, seeds_tree_origin, seeds_tree_right, seeds_tree_left);

	double* seeds_sorted = new double[num_seeds * 3];
	size_t* seeds_old_index = new size_t[num_seeds];
	size_t* seeds_new_index = new size_t[num_seeds];
	re_enumerate_points_for_better_memory_access(num_seeds, 3, seeds,
		                                         seeds_tree_origin, seeds_tree_right, seeds_tree_left,
		                                         seeds_sorted, seeds_old_index, seeds_new_index);

	size_t num_lip(0);
	if (true)
	{
		double* seed_sizing_sorted = new double[num_seeds];
		for (size_t i = 0; i < num_seeds; i++) seed_sizing_sorted[i] = seeds_sizing[seeds_old_index[i]];
		impose_lipschitz_continuity(num_seeds, 3, seeds_sorted, seed_sizing_sorted, num_threads, Lip, seeds_tree_origin, seeds_tree_right, seeds_tree_left);
		for (size_t i = 0; i < num_seeds; i++) seeds_sizing[seeds_old_index[i]] = seed_sizing_sorted[i];
		delete[] seed_sizing_sorted;
		num_lip = num_seeds;
	}

	// Surface sphere tree
	size_t  surface_spheres_tree_origin(0);
	size_t* surface_spheres_tree_right(0);
	size_t* surface_spheres_tree_left(0);
	double* surface_spheres_sorted(0);
	size_t* surface_spheres_old_index(0);
	size_t* surface_spheres_new_index(0);

	if (num_surface_spheres > 0)
	{
		surface_spheres_tree_right = new size_t[num_surface_spheres];
		surface_spheres_tree_left = new size_t[num_surface_spheres];
		build_balanced_kd_tree(num_surface_spheres, 3, surface_spheres, surface_spheres_tree_origin, surface_spheres_tree_right, surface_spheres_tree_left);

		surface_spheres_sorted = new double[num_surface_spheres * 3];
		surface_spheres_old_index = new size_t[num_surface_spheres];
		surface_spheres_new_index = new size_t[num_surface_spheres];
		re_enumerate_points_for_better_memory_access(num_surface_spheres, 3, surface_spheres,
			surface_spheres_tree_origin, surface_spheres_tree_right, surface_spheres_tree_left,
			surface_spheres_sorted, surface_spheres_old_index, surface_spheres_new_index);
	}

	// Edge sphere tree
	size_t  edge_spheres_tree_origin(0);
	size_t* edge_spheres_tree_right(0);
	size_t* edge_spheres_tree_left(0);
	double* edge_spheres_sorted(0);
	size_t* edge_spheres_old_index(0);
	size_t* edge_spheres_new_index(0);

	if (num_edge_spheres > 0)
	{
		edge_spheres_tree_right = new size_t[num_edge_spheres];
		edge_spheres_tree_left = new size_t[num_edge_spheres];
		build_balanced_kd_tree(num_edge_spheres, 3, edge_spheres, edge_spheres_tree_origin, edge_spheres_tree_right, edge_spheres_tree_left);

		edge_spheres_sorted = new double[num_edge_spheres * 3];
		edge_spheres_old_index = new size_t[num_edge_spheres];
		edge_spheres_new_index = new size_t[num_edge_spheres];
		re_enumerate_points_for_better_memory_access(num_edge_spheres, 3, edge_spheres,
			edge_spheres_tree_origin, edge_spheres_tree_right, edge_spheres_tree_left,
			edge_spheres_sorted, edge_spheres_old_index, edge_spheres_new_index);
	}

	// Corner sphere tree
	size_t  corner_spheres_tree_origin(0);
	size_t* corner_spheres_tree_right(0);
	size_t* corner_spheres_tree_left(0);
	double* corner_spheres_sorted(0);
	size_t* corner_spheres_old_index(0);
	size_t* corner_spheres_new_index(0);

	if (num_corner_spheres > 0)
	{
		corner_spheres_tree_right = new size_t[num_surface_spheres];
		corner_spheres_tree_left = new size_t[num_surface_spheres];
		build_balanced_kd_tree(num_corner_spheres, 3, corner_spheres, corner_spheres_tree_origin, corner_spheres_tree_right, corner_spheres_tree_left);

		corner_spheres_sorted = new double[num_corner_spheres * 3];
		corner_spheres_old_index = new size_t[num_corner_spheres];
		corner_spheres_new_index = new size_t[num_corner_spheres];
		re_enumerate_points_for_better_memory_access(num_corner_spheres, 3, corner_spheres,
			corner_spheres_tree_origin, corner_spheres_tree_right, corner_spheres_tree_left,
			corner_spheres_sorted, corner_spheres_old_index, corner_spheres_new_index);
	}

	// Sizing spheres tree
	size_t  sizing_points_tree_origin(0);
	size_t* sizing_points_tree_right(0);
	size_t* sizing_points_tree_left(0);
	double* sizing_points_sorted(0);
	size_t* sizing_points_old_index(0);
	size_t* sizing_points_new_index(0);

	if (num_sizing_points > 0)
	{
		sizing_points_tree_right = new size_t[num_surface_spheres];
		sizing_points_tree_left = new size_t[num_surface_spheres];
		build_balanced_kd_tree(num_sizing_points, 3, sizing_points, sizing_points_tree_origin, sizing_points_tree_right, sizing_points_tree_left);

		sizing_points_sorted = new double[num_sizing_points * 3];
		sizing_points_old_index = new size_t[num_sizing_points];
		sizing_points_new_index = new size_t[num_sizing_points];
		re_enumerate_points_for_better_memory_access(num_sizing_points, 3, sizing_points,
			sizing_points_tree_origin, sizing_points_tree_right, sizing_points_tree_left,
			sizing_points_sorted, sizing_points_old_index, sizing_points_new_index);
	}

	MeshingRandomSampler rsampler;

	double* xo = new double[4];
	double* xn = new double[4];
	double* pv = new double[4];
	double* y = new double[2];
	double* x = new double[4];
	double* x_best = new double[4];

	double* eo = new double[3];
	double* e1 = new double[3];
	double* e2 = new double[3];

	size_t added_seeds_cap(100), num_added_seeds(0);
	double* added_seeds = new double[added_seeds_cap * 3];
	double* added_seeds_sizing = new double[added_seeds_cap];
	size_t* added_seeds_id = new size_t[added_seeds_cap];

	for (int isphere = 0; isphere < num_edge_spheres; isphere++)
	{
		xo[0] = edge_spheres[isphere * 3];
		xo[1] = edge_spheres[isphere * 3 + 1];
		xo[2] = edge_spheres[isphere * 3 + 2];
		xo[3] = edge_spheres_sizing[isphere];
		
		for (size_t ii = 0; ii < 2; ii++)
		{
			double* e_dir(0);
			if (ii == 1) e_dir = eo;

			size_t iclosest(SIZE_MAX); double hclosest(DBL_MAX);
			get_closest_tree_point(num_edge_spheres, 3, edge_spheres_sorted,
				                   edge_spheres_tree_origin, edge_spheres_tree_right, edge_spheres_tree_left,
				                   edge_spheres_new_index[isphere], e_dir, iclosest, hclosest);

			if (iclosest != SIZE_MAX && hclosest < xo[3] + edge_spheres_sizing[edge_spheres_old_index[iclosest]])
			{
				size_t jsphere = edge_spheres_old_index[iclosest];
				xn[0] = edge_spheres[jsphere * 3];
				xn[1] = edge_spheres[jsphere * 3 + 1];
				xn[2] = edge_spheres[jsphere * 3 + 2];
				xn[3] = edge_spheres_sizing[jsphere];
			}
			else
			{
				iclosest = SIZE_MAX; hclosest = DBL_MAX;
				get_closest_tree_point(num_corner_spheres, 3, corner_spheres_sorted,
					                   corner_spheres_tree_origin, corner_spheres_tree_right, corner_spheres_tree_left,
					                   xo, e_dir, iclosest, hclosest);

				size_t jsphere = corner_spheres_old_index[iclosest];
				xn[0] = corner_spheres[jsphere * 3];
				xn[1] = corner_spheres[jsphere * 3 + 1];
				xn[2] = corner_spheres[jsphere * 3 + 2];
				xn[3] = corner_spheres_sizing[jsphere];
			}
			
			for (size_t idim = 0; idim < 3; idim++) eo[idim] = xo[idim] - xn[idim];

			double norm(0.0);
			for (size_t idim = 0; idim < 3; idim++) norm += eo[idim] * eo[idim];
			norm = sqrt(norm);

			for (size_t idim = 0; idim < 3; idim++) eo[idim] /= norm;
			
			while (true)
			{
				rsampler.sample_uniformly_from_unit_sphere(e1, 3);
				double dot(0.0);
				for (size_t idim = 0; idim < 3; idim++) dot += e1[idim] * eo[idim];
				for (size_t idim = 0; idim < 3; idim++) e1[idim] -= dot * eo[idim];

				norm = 0.0;
				for (size_t idim = 0; idim < 3; idim++) norm += e1[idim] * e1[idim];
				norm = sqrt(norm);

				if (norm < 0.001) continue;

				for (size_t idim = 0; idim < 3; idim++) e1[idim] /= norm;

				break;
			}
			
			e2[0] = eo[1] * e1[2] - eo[2] * e1[1];
			e2[1] = eo[2] * e1[0] - eo[0] * e1[2];
			e2[2] = eo[0] * e1[1] - eo[1] * e1[0];

			double dst(hclosest);
			double h = (xo[3] * xo[3] - xn[3] * xn[3] + dst * dst) / (2.0 * dst);

			for (size_t idim = 0; idim < 3; idim++) pv[idim] = xo[idim] + (h / dst) * (xn[idim] - xo[idim]);

			pv[3] = sqrt(xo[3] * xo[3] - h * h);

			size_t num_neighbor_seeds(0);
			size_t* neighbor_seeds(0);
			get_tree_points_in_sphere(num_seeds, 3, seeds_sorted, seeds_tree_origin, seeds_tree_right, seeds_tree_left, pv, pv[3], num_neighbor_seeds, neighbor_seeds);
			delete[] neighbor_seeds;

			if (num_neighbor_seeds > 2) continue;			

			double h_best(0.0); size_t num_mc(10);
			size_t id_best(SIZE_MAX);
			for (size_t imc = 0; imc < num_mc; imc++)
			{
				rsampler.sample_uniformly_from_unit_sphere(y, 2);
				for (size_t idim = 0; idim < 3; idim++) x[idim] = pv[idim] + pv[3] * (y[0] * e1[idim] + y[1] * e2[idim]);
				
				if (num_surface_spheres > 0)
				{
					#pragma region Make sure seed is outside surface spheres:
					size_t num_neighbor_spheres(0);
					size_t* neighbor_spheres(0);

					size_t iclosest(num_surface_spheres); double hclosest(DBL_MAX);
					get_closest_tree_point(num_surface_spheres, 3, surface_spheres_sorted,
						                   surface_spheres_tree_origin, surface_spheres_tree_right, surface_spheres_tree_left,
						                   x, iclosest, hclosest);

					size_t sphere_index = surface_spheres_old_index[iclosest];
					double rs = surface_spheres_sizing[sphere_index] + Lip * hclosest;

					double R = rs / (1 - Lip);

					get_tree_points_in_sphere(num_surface_spheres, 3, surface_spheres_sorted,
						                      surface_spheres_tree_origin, surface_spheres_tree_right, surface_spheres_tree_left,
						                      x, R, num_neighbor_spheres, neighbor_spheres);

					bool bad_seed(false);
					for (size_t i = 0; i < num_neighbor_spheres; i++)
					{
						size_t isphere = neighbor_spheres[i];
						double h(0.0);
						for (size_t idim = 0; idim < 3; idim++)
						{
							double dx = x[idim] - surface_spheres_sorted[isphere * 3 + idim];
							h += dx * dx;
						}
						h = sqrt(h);
						isphere = surface_spheres_old_index[isphere];
						if (h < surface_spheres_sizing[isphere] - 1E-10)
						{
							bad_seed = true;
							break;
						}
					}

					if (neighbor_spheres != 0) delete[] neighbor_spheres;
					if (bad_seed) continue;
					#pragma endregion
				}
				
				if (true)
				{
					#pragma region Make sure seed sphere is not too close to an existing seed:
					size_t iclosest(num_seeds); double hclosest(DBL_MAX);
					get_closest_tree_point(num_seeds, 3, seeds_sorted,
						                   seeds_tree_origin, seeds_tree_right, seeds_tree_left,
						                   x, iclosest, hclosest);

					size_t sphere_index = seeds_old_index[iclosest];
					if (seeds_region_id[sphere_index] == 0) continue; // An exterior seed

					size_t new_seed_id = seeds_region_id[sphere_index];

					if (hclosest > h_best)
					{
						h_best = hclosest;
						for (size_t idim = 0; idim < 3; idim++) x_best[idim] = x[idim];
						x_best[3] = pv[3];
						id_best = new_seed_id;
					}					
					#pragma endregion
				}				
			}

			if (h_best < 1E-10)
			{
				vcm_cout << "VoroCrust Warning!, could not find a proper seed location to ensure Edge Sphere is Delaunay!!" << std::endl;
				vcm_cout << "Sphere i: " << xo[0] << " " << xo[1] << " " << xo[2] << " " << xo[3] << std::endl;
				vcm_cout << "Sphere j: " << xn[0] << " " << xn[1] << " " << xn[2] << " " << xn[3] << std::endl;				

				if (e_dir != 0)
				{
					vcm_cout << "e_dir: " << e_dir[0] << " " << e_dir[1] << " " << e_dir[2] << std::endl;
				}				
				continue;
			}

			for (size_t idim = 0; idim < 3; idim++) added_seeds[num_added_seeds * 3 + idim] = x_best[idim];
			added_seeds_sizing[num_added_seeds] = x_best[3];
			added_seeds_id[num_added_seeds] = id_best;
			num_added_seeds++;

			if (num_added_seeds == added_seeds_cap)
			{
				added_seeds_cap *= 2;
				double* tmp_added_seeds = new double[added_seeds_cap * 3];
				double* tmp_added_seeds_sizing = new double[added_seeds_cap];
				size_t* tmp_added_seeds_id = new size_t[added_seeds_cap];

				for (size_t ii = 0; ii < num_added_seeds * 3; ii++) tmp_added_seeds[ii] = added_seeds[ii];
				for (size_t ii = 0; ii < num_added_seeds; ii++) tmp_added_seeds_sizing[ii] = added_seeds_sizing[ii];
				for (size_t ii = 0; ii < num_added_seeds; ii++) tmp_added_seeds_id[ii] = added_seeds_id[ii];

				delete[] added_seeds; added_seeds = tmp_added_seeds;
				delete[] added_seeds_sizing; added_seeds_sizing = tmp_added_seeds_sizing;
				delete[] added_seeds_id; added_seeds_id = tmp_added_seeds_id;
			}
		}		
	}
	vcm_cout << "  * Number of added seeds " << num_added_seeds << std::endl;

	if (num_added_seeds > 0)
	{
		#pragma region Adding New Seeds
		double* new_seeds = new double[(num_seeds + num_added_seeds) * 3];
		double* new_seeds_sizing = new double[num_seeds + num_added_seeds];
		size_t* new_seeds_region_id = new size_t[num_seeds + num_added_seeds];
		for (size_t i = 0; i < num_seeds * 3; i++) new_seeds[i] = seeds[i];
		for (size_t i = 0; i < num_seeds; i++) new_seeds_sizing[i] = seeds_sizing[i];
		for (size_t i = 0; i < num_seeds; i++) new_seeds_region_id[i] = seeds_region_id[i];

		
		for (size_t ii = 0; ii < num_added_seeds * 3; ii++) new_seeds[num_seeds * 3 + ii] = added_seeds[ii];
		for (size_t ii = 0; ii < num_added_seeds; ii++) new_seeds_sizing[num_seeds + ii] = added_seeds_sizing[ii];
		for (size_t ii = 0; ii < num_added_seeds; ii++) new_seeds_region_id[num_seeds + ii] = added_seeds_id[ii];

		delete[] seeds; delete[] seeds_sizing; delete[] seeds_region_id;
		num_seeds += num_added_seeds;
		seeds = new_seeds; seeds_sizing = new_seeds_sizing;
		seeds_region_id = new_seeds_region_id;
		#pragma endregion
	}
	

	delete[] xo; delete[] xn; delete[] pv;
	delete[] y; delete[] x; delete[] x_best;

	delete[] eo; delete[] e1; delete[] e2;

	delete[] added_seeds; delete[] added_seeds_sizing; delete[] added_seeds_id;

	delete[] seeds_tree_right; delete[] seeds_tree_left;
	delete[] seeds_sorted; delete[] seeds_old_index; delete[] seeds_new_index;

	if (num_surface_spheres > 0)
	{
		delete[] surface_spheres_tree_right; delete[] surface_spheres_tree_left;
		delete[] surface_spheres_sorted; delete[] surface_spheres_old_index; delete[] surface_spheres_new_index;
	}

	if (num_edge_spheres > 0)
	{
		delete[] edge_spheres_tree_right; delete[] edge_spheres_tree_left;
		delete[] edge_spheres_sorted; delete[] edge_spheres_old_index; delete[] edge_spheres_new_index;
	}

	if (num_corner_spheres > 0)
	{
		delete[] corner_spheres_tree_right; delete[] corner_spheres_tree_left;
		delete[] corner_spheres_sorted; delete[] corner_spheres_old_index; delete[] corner_spheres_new_index;
	}

	if (num_sizing_points > 0)
	{
		delete[] sizing_points_tree_right; delete[] sizing_points_tree_left;
		delete[] sizing_points_sorted; delete[] sizing_points_old_index; delete[] sizing_points_new_index;
	}
	return 0;
	#pragma endregion
}

int MeshingVoronoiMesher::impose_interior_seeds(size_t& num_seeds, double*& seeds, double*& seeds_sizing, size_t*& seeds_region_id,
	                                     size_t num_surface_spheres, double* surface_spheres, double* surface_spheres_sizing,
	                                     size_t num_edge_spheres, double* edge_spheres, double* edge_spheres_sizing,
	                                     size_t num_corner_spheres, double* corner_spheres, double* corner_spheres_sizing,
	                                     size_t num_sizing_points, double* sizing_points, double* sizing_value,
	                                     size_t num_imposed_seeds, double* imposed_seeds,
	                                     int num_threads, double Lip, double rmax)
{
	#pragma region Impose Interior Seeds:
	// Seeds tree
	size_t num_seeds_reported(num_seeds);
	vcm_cout << "  * Number of seeds " << num_seeds << std::endl;
	vcm_cout << "  * Number of seeds to impose " << num_imposed_seeds << std::endl;

	size_t  seeds_tree_origin(0);
	size_t* seeds_tree_right = new size_t[num_seeds];
	size_t* seeds_tree_left = new size_t[num_seeds];
	build_balanced_kd_tree(num_seeds, 3, seeds, seeds_tree_origin, seeds_tree_right, seeds_tree_left);

	double* seeds_sorted = new double[num_seeds * 3];
	size_t* seeds_old_index = new size_t[num_seeds];
	size_t* seeds_new_index = new size_t[num_seeds];
	re_enumerate_points_for_better_memory_access(num_seeds, 3, seeds,
		                                         seeds_tree_origin, seeds_tree_right, seeds_tree_left,
		                                         seeds_sorted, seeds_old_index, seeds_new_index);

	size_t num_lip(0);
	if (true)
	{
		double* seed_sizing_sorted = new double[num_seeds];
		for (size_t i = 0; i < num_seeds; i++) seed_sizing_sorted[i] = seeds_sizing[seeds_old_index[i]];
		impose_lipschitz_continuity(num_seeds, 3, seeds_sorted, seed_sizing_sorted, num_threads, Lip, seeds_tree_origin, seeds_tree_right, seeds_tree_left);
		for (size_t i = 0; i < num_seeds; i++) seeds_sizing[seeds_old_index[i]] = seed_sizing_sorted[i];
		delete[] seed_sizing_sorted;
		num_lip = num_seeds;
	}

	// Surface sphere tree
	size_t  surface_spheres_tree_origin(0);
	size_t* surface_spheres_tree_right(0);
	size_t* surface_spheres_tree_left(0);
	double* surface_spheres_sorted(0);
	size_t* surface_spheres_old_index(0);
	size_t* surface_spheres_new_index(0);

	if (num_surface_spheres > 0)
	{
		surface_spheres_tree_right = new size_t[num_surface_spheres];
		surface_spheres_tree_left = new size_t[num_surface_spheres];
		build_balanced_kd_tree(num_surface_spheres, 3, surface_spheres, surface_spheres_tree_origin, surface_spheres_tree_right, surface_spheres_tree_left);

		surface_spheres_sorted = new double[num_surface_spheres * 3];
		surface_spheres_old_index = new size_t[num_surface_spheres];
		surface_spheres_new_index = new size_t[num_surface_spheres];
		re_enumerate_points_for_better_memory_access(num_surface_spheres, 3, surface_spheres,
			surface_spheres_tree_origin, surface_spheres_tree_right, surface_spheres_tree_left,
			surface_spheres_sorted, surface_spheres_old_index, surface_spheres_new_index);
	}

	// Edge sphere tree
	size_t  edge_spheres_tree_origin(0);
	size_t* edge_spheres_tree_right(0);
	size_t* edge_spheres_tree_left(0);
	double* edge_spheres_sorted(0);
	size_t* edge_spheres_old_index(0);
	size_t* edge_spheres_new_index(0);

	if (num_edge_spheres > 0)
	{
		edge_spheres_tree_right = new size_t[num_edge_spheres];
		edge_spheres_tree_left = new size_t[num_edge_spheres];
		build_balanced_kd_tree(num_edge_spheres, 3, edge_spheres, edge_spheres_tree_origin, edge_spheres_tree_right, edge_spheres_tree_left);

		edge_spheres_sorted = new double[num_edge_spheres * 3];
		edge_spheres_old_index = new size_t[num_edge_spheres];
		edge_spheres_new_index = new size_t[num_edge_spheres];
		re_enumerate_points_for_better_memory_access(num_edge_spheres, 3, edge_spheres,
			edge_spheres_tree_origin, edge_spheres_tree_right, edge_spheres_tree_left,
			edge_spheres_sorted, edge_spheres_old_index, edge_spheres_new_index);
	}

	// Corner sphere tree
	size_t  corner_spheres_tree_origin(0);
	size_t* corner_spheres_tree_right(0);
	size_t* corner_spheres_tree_left(0);
	double* corner_spheres_sorted(0);
	size_t* corner_spheres_old_index(0);
	size_t* corner_spheres_new_index(0);

	if (num_corner_spheres > 0)
	{
		corner_spheres_tree_right = new size_t[num_surface_spheres];
		corner_spheres_tree_left = new size_t[num_surface_spheres];
		build_balanced_kd_tree(num_corner_spheres, 3, corner_spheres, corner_spheres_tree_origin, corner_spheres_tree_right, corner_spheres_tree_left);

		corner_spheres_sorted = new double[num_corner_spheres * 3];
		corner_spheres_old_index = new size_t[num_corner_spheres];
		corner_spheres_new_index = new size_t[num_corner_spheres];
		re_enumerate_points_for_better_memory_access(num_corner_spheres, 3, corner_spheres,
			corner_spheres_tree_origin, corner_spheres_tree_right, corner_spheres_tree_left,
			corner_spheres_sorted, corner_spheres_old_index, corner_spheres_new_index);
	}

	// Sizing spheres tree
	size_t  sizing_points_tree_origin(0);
	size_t* sizing_points_tree_right(0);
	size_t* sizing_points_tree_left(0);
	double* sizing_points_sorted(0);
	size_t* sizing_points_old_index(0);
	size_t* sizing_points_new_index(0);

	if (num_sizing_points > 0)
	{
		sizing_points_tree_right = new size_t[num_surface_spheres];
		sizing_points_tree_left = new size_t[num_surface_spheres];
		build_balanced_kd_tree(num_sizing_points, 3, sizing_points, sizing_points_tree_origin, sizing_points_tree_right, sizing_points_tree_left);

		sizing_points_sorted = new double[num_sizing_points * 3];
		sizing_points_old_index = new size_t[num_sizing_points];
		sizing_points_new_index = new size_t[num_sizing_points];
		re_enumerate_points_for_better_memory_access(num_sizing_points, 3, sizing_points,
			sizing_points_tree_origin, sizing_points_tree_right, sizing_points_tree_left,
			sizing_points_sorted, sizing_points_old_index, sizing_points_new_index);
	}

	size_t num_valid_seeds(0);
	bool* valid_seed = new bool[num_imposed_seeds];
	size_t* imposed_seed_region_id = new size_t[num_imposed_seeds];

#if defined USE_OPEN_MP
	omp_set_num_threads(num_threads);
#pragma omp parallel for
#endif
	for (int iseed = 0; iseed < num_imposed_seeds; iseed++)
	{
		int thread_id(0);
#if defined USE_OPEN_MP
		thread_id = omp_get_thread_num();
#endif
		double* x = new double[4];

		x[0] = imposed_seeds[iseed * 4];
		x[1] = imposed_seeds[iseed * 4 + 1];
		x[2] = imposed_seeds[iseed * 4 + 2];
		x[3] = imposed_seeds[iseed * 4 + 3];

		valid_seed[iseed] = false;

		if (true)
		{
			#pragma region Make sure seed is inside domain:
			size_t num_neighbor_spheres(0);
			size_t* neighbor_spheres(0);

			size_t iclosest(num_seeds); double hclosest(DBL_MAX);
			get_closest_tree_point(num_seeds, 3, seeds_sorted,
				seeds_tree_origin, seeds_tree_right, seeds_tree_left,
				x, iclosest, hclosest);

			size_t sphere_index = seeds_old_index[iclosest];
			if (seeds_region_id[sphere_index] == 0)
			{
				delete[] x;
				continue;
			}

			imposed_seed_region_id[iseed] = seeds_region_id[sphere_index];

			#pragma endregion
		}

		if (num_surface_spheres > 0)
		{
			#pragma region Make sure seed is outside surface spheres:
			size_t num_neighbor_spheres(0);
			size_t* neighbor_spheres(0);

			size_t iclosest(num_surface_spheres); double hclosest(DBL_MAX);
			get_closest_tree_point(num_surface_spheres, 3, surface_spheres_sorted,
				surface_spheres_tree_origin, surface_spheres_tree_right, surface_spheres_tree_left,
				x, iclosest, hclosest);

			size_t sphere_index = surface_spheres_old_index[iclosest];
			double rs = surface_spheres_sizing[sphere_index] + Lip * hclosest;

			double R = rs / (1 - Lip);

			get_tree_points_in_sphere(num_surface_spheres, 3, surface_spheres_sorted,
				surface_spheres_tree_origin, surface_spheres_tree_right, surface_spheres_tree_left,
				x, R, num_neighbor_spheres, neighbor_spheres);

			bool bad_seed(false);
			for (size_t i = 0; i < num_neighbor_spheres; i++)
			{
				size_t isphere = neighbor_spheres[i];
				double h(0.0);
				for (size_t idim = 0; idim < 3; idim++)
				{
					double dx = x[idim] - surface_spheres_sorted[isphere * 3 + idim];
					h += dx * dx;
				}
				h = sqrt(h);
				isphere = surface_spheres_old_index[isphere];
				if (h < surface_spheres_sizing[isphere] - 1E-10)
				{
					bad_seed = true;
					break;
				}
			}
			
			delete[] neighbor_spheres;
			if (bad_seed)
			{
				delete[] x;
				continue;
			}
			#pragma endregion
		}

		if (num_edge_spheres)
		{
			#pragma region Make sure seed is outside edge spheres:
			size_t num_neighbor_spheres(0);
			size_t* neighbor_spheres(0);

			size_t iclosest(num_edge_spheres); double hclosest(DBL_MAX);
			get_closest_tree_point(num_edge_spheres, 3, edge_spheres_sorted,
				edge_spheres_tree_origin, edge_spheres_tree_right, edge_spheres_tree_left,
				x, iclosest, hclosest);

			size_t sphere_index = edge_spheres_old_index[iclosest];
			double rs = edge_spheres_sizing[sphere_index] + Lip * hclosest;

			double R = rs / (1 - Lip);

			get_tree_points_in_sphere(num_edge_spheres, 3, edge_spheres_sorted,
				edge_spheres_tree_origin, edge_spheres_tree_right, edge_spheres_tree_left,
				x, R, num_neighbor_spheres, neighbor_spheres);

			bool bad_seed(false);
			for (size_t i = 0; i < num_neighbor_spheres; i++)
			{
				size_t isphere = neighbor_spheres[i];
				double h(0.0);
				for (size_t idim = 0; idim < 3; idim++)
				{
					double dx = x[idim] - edge_spheres_sorted[isphere * 3 + idim];
					h += dx * dx;
				}
				h = sqrt(h);
				isphere = edge_spheres_old_index[isphere];
				if (h < edge_spheres_sizing[isphere] - 1E-10)
				{
					bad_seed = true;
					break;
				}
			}
			
			delete[] neighbor_spheres;
			if (bad_seed)
			{
				delete[] x;
				continue;
			}
			#pragma endregion
		}

		if (num_corner_spheres)
		{
			#pragma region Make sure seed is outside corner spheres:
			size_t num_neighbor_spheres(0);
			size_t* neighbor_spheres(0);

			size_t iclosest(num_corner_spheres); double hclosest(DBL_MAX);
			get_closest_tree_point(num_corner_spheres, 3, corner_spheres_sorted,
				corner_spheres_tree_origin, corner_spheres_tree_right, corner_spheres_tree_left,
				x, iclosest, hclosest);

			size_t sphere_index = corner_spheres_old_index[iclosest];
			double rs = corner_spheres_sizing[sphere_index] + Lip * hclosest;

			double R = rs / (1 - Lip);

			get_tree_points_in_sphere(num_corner_spheres, 3, corner_spheres_sorted,
				corner_spheres_tree_origin, corner_spheres_tree_right, corner_spheres_tree_left,
				x, R, num_neighbor_spheres, neighbor_spheres);

			bool bad_seed(false);
			for (size_t i = 0; i < num_neighbor_spheres; i++)
			{
				size_t isphere = neighbor_spheres[i];
				double h(0.0);
				for (size_t idim = 0; idim < 3; idim++)
				{
					double dx = x[idim] - corner_spheres_sorted[isphere * 3 + idim];
					h += dx * dx;
				}
				h = sqrt(h);
				isphere = corner_spheres_old_index[isphere];
				if (h < corner_spheres_sizing[isphere] - 1E-10)
				{
					bad_seed = true;
					break;
				}
			}
			
			delete[] neighbor_spheres;
			if (bad_seed)
			{
				delete[] x;
				continue;
			}
			#pragma endregion
		}
		delete[] x;
		valid_seed[iseed] = true;
		num_valid_seeds++;
	}

	vcm_cout << "  * Number of imposed seeds " << num_valid_seeds << std::endl;


	if (num_valid_seeds > 0)
	{
		#pragma region Adding New Seeds
		double* new_seeds = new double[(num_seeds + num_valid_seeds) * 3];
		double* new_seeds_sizing = new double[num_seeds + num_valid_seeds];
		size_t* new_seeds_region_id = new size_t[num_seeds + num_valid_seeds];
		for (size_t i = 0; i < num_seeds * 3; i++) new_seeds[i] = seeds[i];
		for (size_t i = 0; i < num_seeds; i++) new_seeds_sizing[i] = seeds_sizing[i];
		for (size_t i = 0; i < num_seeds; i++) new_seeds_region_id[i] = seeds_region_id[i];

		num_valid_seeds = 0;
		for (size_t iseed = 0; iseed < num_imposed_seeds; iseed++)
		{
			if (!valid_seed[iseed]) continue;

			double* x = new double[4];
			x[0] = imposed_seeds[iseed * 4];
			x[1] = imposed_seeds[iseed * 4 + 1];
			x[2] = imposed_seeds[iseed * 4 + 2];
			x[3] = imposed_seeds[iseed * 4 + 3];

			new_seeds[(num_seeds + num_valid_seeds) * 3] = x[0];
			new_seeds[(num_seeds + num_valid_seeds) * 3 + 1] = x[1];
			new_seeds[(num_seeds + num_valid_seeds) * 3 + 2] = x[2];
			new_seeds_sizing[num_seeds + num_valid_seeds] = x[3];
			new_seeds_region_id[num_seeds + num_valid_seeds] = imposed_seed_region_id[iseed];
			delete[] x;
			num_valid_seeds++;
		}
		delete[] seeds; delete[] seeds_sizing; delete[] seeds_region_id;
		num_seeds += num_valid_seeds; 
		seeds = new_seeds; seeds_sizing = new_seeds_sizing;
		seeds_region_id = new_seeds_region_id;
		#pragma endregion
	}
	delete[] valid_seed; delete[] imposed_seed_region_id;
	
	delete[] seeds_tree_right; delete[] seeds_tree_left;
	delete[] seeds_sorted; delete[] seeds_old_index; delete[] seeds_new_index;

	if (num_surface_spheres > 0)
	{
		delete[] surface_spheres_tree_right; delete[] surface_spheres_tree_left;
		delete[] surface_spheres_sorted; delete[] surface_spheres_old_index; delete[] surface_spheres_new_index;
	}

	if (num_edge_spheres > 0)
	{
		delete[] edge_spheres_tree_right; delete[] edge_spheres_tree_left;
		delete[] edge_spheres_sorted; delete[] edge_spheres_old_index; delete[] edge_spheres_new_index;
	}

	if (num_corner_spheres > 0)
	{
		delete[] corner_spheres_tree_right; delete[] corner_spheres_tree_left;
		delete[] corner_spheres_sorted; delete[] corner_spheres_old_index; delete[] corner_spheres_new_index;
	}

	if (num_sizing_points > 0)
	{
		delete[] sizing_points_tree_right; delete[] sizing_points_tree_left;
		delete[] sizing_points_sorted; delete[] sizing_points_old_index; delete[] sizing_points_new_index;
	}
	return 0;
	#pragma endregion
}


int MeshingVoronoiMesher::generate_interior_seeds(size_t& num_seeds, double*& seeds, double*& seeds_sizing, size_t*& seeds_region_id,
	                                       size_t num_surface_spheres, double* surface_spheres, double* surface_spheres_sizing,
	                                       size_t num_edge_spheres, double* edge_spheres, double* edge_spheres_sizing,
	                                       size_t num_corner_spheres, double* corner_spheres, double* corner_spheres_sizing,
	                                       size_t num_sizing_points, double* sizing_points, double* sizing_value,
	                                       int num_threads, double Lip, double rmax)
{
	#pragma region Generate Interior Seeds:
	// Seeds tree
	size_t num_seeds_reported(num_seeds);
	vcm_cout << "  * Number of seeds " << num_seeds << std::endl;

	size_t  seeds_tree_origin(0);
	size_t* seeds_tree_right = new size_t[num_seeds];
	size_t* seeds_tree_left = new size_t[num_seeds];
	build_balanced_kd_tree(num_seeds, 3, seeds, seeds_tree_origin, seeds_tree_right, seeds_tree_left);

	double* seeds_sorted = new double[num_seeds * 3];
	size_t* seeds_old_index = new size_t[num_seeds];
	size_t* seeds_new_index = new size_t[num_seeds];
	re_enumerate_points_for_better_memory_access(num_seeds, 3, seeds,
		                                         seeds_tree_origin, seeds_tree_right, seeds_tree_left,
		                                         seeds_sorted, seeds_old_index, seeds_new_index);

	size_t num_lip(0);
	if (true)
	{
		double* seed_sizing_sorted = new double[num_seeds];
		for (size_t i = 0; i < num_seeds; i++) seed_sizing_sorted[i] = seeds_sizing[seeds_old_index[i]];
		impose_lipschitz_continuity(num_seeds, 3, seeds_sorted, seed_sizing_sorted, num_threads, Lip, seeds_tree_origin, seeds_tree_right, seeds_tree_left);
		for (size_t i = 0; i < num_seeds; i++) seeds_sizing[seeds_old_index[i]] = seed_sizing_sorted[i];
		delete[] seed_sizing_sorted;
		num_lip = num_seeds;
	}

	// Surface sphere tree
	size_t  surface_spheres_tree_origin(0);
	size_t* surface_spheres_tree_right(0);
	size_t* surface_spheres_tree_left(0);
	double* surface_spheres_sorted(0);
	size_t* surface_spheres_old_index(0);
	size_t* surface_spheres_new_index(0);

	if (num_surface_spheres > 0)
	{
		surface_spheres_tree_right = new size_t[num_surface_spheres];
		surface_spheres_tree_left = new size_t[num_surface_spheres];
		build_balanced_kd_tree(num_surface_spheres, 3, surface_spheres, surface_spheres_tree_origin, surface_spheres_tree_right, surface_spheres_tree_left);

		surface_spheres_sorted = new double[num_surface_spheres * 3];
		surface_spheres_old_index = new size_t[num_surface_spheres];
		surface_spheres_new_index = new size_t[num_surface_spheres];
		re_enumerate_points_for_better_memory_access(num_surface_spheres, 3, surface_spheres,
			                                         surface_spheres_tree_origin, surface_spheres_tree_right, surface_spheres_tree_left,
			                                         surface_spheres_sorted, surface_spheres_old_index, surface_spheres_new_index);
	}
	
	// Edge sphere tree
	size_t  edge_spheres_tree_origin(0);
	size_t* edge_spheres_tree_right(0);
	size_t* edge_spheres_tree_left(0);
	double* edge_spheres_sorted(0);
	size_t* edge_spheres_old_index(0);
	size_t* edge_spheres_new_index(0);

	if (num_edge_spheres > 0)
	{
		edge_spheres_tree_right = new size_t[num_edge_spheres];
		edge_spheres_tree_left = new size_t[num_edge_spheres];
		build_balanced_kd_tree(num_edge_spheres, 3, edge_spheres, edge_spheres_tree_origin, edge_spheres_tree_right, edge_spheres_tree_left);

		edge_spheres_sorted = new double[num_edge_spheres * 3];
		edge_spheres_old_index = new size_t[num_edge_spheres];
		edge_spheres_new_index = new size_t[num_edge_spheres];
		re_enumerate_points_for_better_memory_access(num_edge_spheres, 3, edge_spheres,
			                                         edge_spheres_tree_origin, edge_spheres_tree_right, edge_spheres_tree_left,
			                                         edge_spheres_sorted, edge_spheres_old_index, edge_spheres_new_index);
	}
	

	// Corner sphere tree
	size_t  corner_spheres_tree_origin(0);
	size_t* corner_spheres_tree_right(0);
	size_t* corner_spheres_tree_left(0);
	double* corner_spheres_sorted(0);
	size_t* corner_spheres_old_index(0);
	size_t* corner_spheres_new_index(0);

	if (num_corner_spheres > 0)
	{
		corner_spheres_tree_right = new size_t[num_surface_spheres];
		corner_spheres_tree_left = new size_t[num_surface_spheres];
		build_balanced_kd_tree(num_corner_spheres, 3, corner_spheres, corner_spheres_tree_origin, corner_spheres_tree_right, corner_spheres_tree_left);

		corner_spheres_sorted = new double[num_corner_spheres * 3];
		corner_spheres_old_index = new size_t[num_corner_spheres];
		corner_spheres_new_index = new size_t[num_corner_spheres];
		re_enumerate_points_for_better_memory_access(num_corner_spheres, 3, corner_spheres,
			                                         corner_spheres_tree_origin, corner_spheres_tree_right, corner_spheres_tree_left,
			                                         corner_spheres_sorted, corner_spheres_old_index, corner_spheres_new_index);
	}
	
	// Sizing spheres tree
	size_t  sizing_points_tree_origin(0);
	size_t* sizing_points_tree_right(0);
	size_t* sizing_points_tree_left(0);
	double* sizing_points_sorted(0);
	size_t* sizing_points_old_index(0);
	size_t* sizing_points_new_index(0);

	if (num_sizing_points > 0)
	{
		sizing_points_tree_right = new size_t[num_surface_spheres];
		sizing_points_tree_left = new size_t[num_surface_spheres];
		build_balanced_kd_tree(num_sizing_points, 3, sizing_points, sizing_points_tree_origin, sizing_points_tree_right, sizing_points_tree_left);

		sizing_points_sorted = new double[num_sizing_points * 3];
		sizing_points_old_index = new size_t[num_sizing_points];
		sizing_points_new_index = new size_t[num_sizing_points];
		re_enumerate_points_for_better_memory_access(num_sizing_points, 3, sizing_points,
			                                         sizing_points_tree_origin, sizing_points_tree_right, sizing_points_tree_left,
			                                         sizing_points_sorted, sizing_points_old_index, sizing_points_new_index);
	}

	std::vector<MeshingRandomSampler*> rsamplers(num_threads);
	for (size_t ithread = 0; ithread < num_threads; ithread++)
	{
		rsamplers[ithread] = new MeshingRandomSampler(ithread);
	}

	size_t num_front_seeds(0);
	for (size_t i = 0; i < num_seeds; i++)
	{
		if (seeds_region_id[i] == 0) continue;
		num_front_seeds++;
	}
	size_t* front_indices = new size_t[num_front_seeds];
	num_front_seeds = 0;
	for (size_t i = 0; i < num_seeds; i++)
	{
		if (seeds_region_id[i] == 0) continue;
		front_indices[num_front_seeds] = i;
		num_front_seeds++;
	}

	double deactivating_time(0.0), spokes_time(0.0), lip_time(0.0), update_seeds_time(0.0);
	
	size_t interior_spheres_cap(1000 * num_threads);
	size_t num_interior_spheres(0);
	double* interior_spheres = new double[interior_spheres_cap * 4];
	size_t* interior_spheres_region_id = new size_t[interior_spheres_cap];

	MeshingTimer timer;

	size_t nit(0);
	while (true)
	{
		nit++;
		// randomize front seeds
		for (size_t i = 0; i < num_front_seeds; i++)
		{
			double u = rsamplers[0]->generate_uniform_random_number();
			size_t j = size_t(u * num_front_seeds);
			if (j == num_front_seeds) j--;

			size_t tmp = front_indices[i];
			front_indices[i] = front_indices[j];
			front_indices[j] = tmp;
		}

		bool* is_active_seed = new bool[num_front_seeds];

		timer.start();
		
		double* front_spheres = new double[num_front_seeds * 4];
		for (size_t iseed = 0; iseed < num_front_seeds; iseed++)
		{
			size_t seed_index(front_indices[iseed]);
			if (seed_index < num_seeds)
			{
				front_spheres[iseed * 4] = seeds[seed_index * 3];
				front_spheres[iseed * 4 + 1] = seeds[seed_index * 3 + 1];
				front_spheres[iseed * 4 + 2] = seeds[seed_index * 3 + 2];
				front_spheres[iseed * 4 + 3] = seeds_sizing[seed_index];
			}
			else
			{
				size_t sphere_index = seed_index - num_seeds;
				front_spheres[iseed * 4] = interior_spheres[sphere_index * 4];
				front_spheres[iseed * 4 + 1] = interior_spheres[sphere_index * 4 + 1];
				front_spheres[iseed * 4 + 2] = interior_spheres[sphere_index * 4 + 2];
				front_spheres[iseed * 4 + 3] = interior_spheres[sphere_index * 4 + 3];
			}
		}

		for (size_t iactive = 0; iactive < num_front_seeds; iactive++) is_active_seed[iactive] = false;

		int num_good_seeds(1); is_active_seed[0] = true;
		for (size_t i = 1; i < num_front_seeds; i++)
		{
			size_t si_index = front_indices[i];

			bool good_seed(true);
			for (size_t j = 0; j < i - 1; j++)
			{
				if (is_active_seed[j])
				{
					size_t sj_index = front_indices[j];

					double h(0.0);
					for (size_t idim = 0; idim < 3; idim++)
					{
						double dx = front_spheres[i * 4 + idim] - front_spheres[j * 4 + idim];
						h += dx * dx;
					}
					h = sqrt(h);
					if (h < (3.0 + 2.0 * Lip) * (front_spheres[i * 4 + 3] + front_spheres[j * 4 + 3]) - 1E-10)
					{
						good_seed = false;
						break;
					}
				}
			}

			if (good_seed)
			{
				is_active_seed[i] = true;
				num_good_seeds++;
				if (num_good_seeds == 100 * num_threads) break;
			}
		}
		delete[] front_spheres;
		
		deactivating_time += timer.report_timing();
		
		size_t num_active_seeds(0);
		for (size_t iactive = 0; iactive < num_front_seeds; iactive++)
		{
			if (is_active_seed[iactive]) num_active_seeds++;
		}

		size_t* active_front_seeds = new size_t[num_active_seeds];
		size_t* inactive_front_seeds = new size_t[num_front_seeds - num_active_seeds];
		
		num_active_seeds = 0; size_t num_inactive_seeds(0);
		for (size_t iactive = 0; iactive < num_front_seeds; iactive++)
		{
			if (is_active_seed[iactive])
			{
				active_front_seeds[num_active_seeds] = front_indices[iactive];
				num_active_seeds++;
			}
			else
			{
				inactive_front_seeds[num_inactive_seeds] = front_indices[iactive];
				num_inactive_seeds++;
			}
		}
		delete[] is_active_seed; delete[] front_indices; front_indices = 0;

		timer.reset_timer();

		size_t* num_new_spheres = new size_t[num_active_seeds];
		double** new_spheres = new double*[num_active_seeds];

#if defined USE_OPEN_MP
		omp_set_num_threads(num_threads);
		#pragma omp parallel for
#endif
		for (int iactive = 0; iactive < num_active_seeds; iactive++)
		{
			// sample new front using spoke darts
			new_spheres[iactive] = 0;
			size_t seed_index = active_front_seeds[iactive];

			int thread_id(0);
#if defined USE_OPEN_MP
			thread_id = omp_get_thread_num();
#endif
			double* x = new double[3]; double r(0.0);

			if (seed_index < num_seeds)
			{
				x[0] = seeds[seed_index * 3];
				x[1] = seeds[seed_index * 3 + 1];
				x[2] = seeds[seed_index * 3 + 2];
				r = seeds_sizing[seed_index];
			}
			else
			{
				size_t sphere_index = seed_index - num_seeds;
				x[0] = interior_spheres[sphere_index * 4];
				x[1] = interior_spheres[sphere_index * 4 + 1];
				x[2] = interior_spheres[sphere_index * 4 + 2];
				r = interior_spheres[sphere_index * 4 + 3];
			}
			
			double alpha = 3.0 / (1 - Lip);
			size_t num_neighbor_seeds(0);
			size_t* neighbor_seeds(0);

			get_tree_points_in_sphere(num_seeds, 3, seeds_sorted, seeds_tree_origin, seeds_tree_right, seeds_tree_left,
				                      x, alpha * r, num_neighbor_seeds, neighbor_seeds);

			size_t num_neighbor_interior_spheres(0);
			size_t* neighbor_interior_spheres(0);
			if (num_interior_spheres > 0)
			{
				size_t cap(10);
				neighbor_interior_spheres = new size_t[cap];
				for (size_t i = 0; i < num_interior_spheres; i++)
				{
					double h(0.0);
					for (size_t idim = 0; idim < 3; idim++)
					{
						double dx = interior_spheres[i * 4 + idim] - x[idim];
						h += dx * dx;
					}
					h = sqrt(h);
					if (h > alpha * r) continue;

					neighbor_interior_spheres[num_neighbor_interior_spheres] = i;
					num_neighbor_interior_spheres++;
					if (num_neighbor_interior_spheres == cap)
					{
						cap *= 2;
						size_t* tmp = new size_t[cap];
						for (size_t j = 0; j < num_neighbor_interior_spheres; j++) tmp[j] = neighbor_interior_spheres[j];
						delete[] neighbor_interior_spheres;
						neighbor_interior_spheres = tmp;
					}
				}
				if (num_neighbor_interior_spheres == 0)
				{
					delete[] neighbor_interior_spheres;
					neighbor_interior_spheres = 0;
				}
			}

			size_t num_neighbor_surface_spheres(0);
			size_t* neighbor_surface_spheres(0);
			if (num_surface_spheres > 0)
			{
				size_t iclosest(num_surface_spheres); double hclosest(DBL_MAX);
				get_closest_tree_point(num_surface_spheres, 3, surface_spheres_sorted,
					                   surface_spheres_tree_origin, surface_spheres_tree_right, surface_spheres_tree_left,
					                   x, iclosest, hclosest);

				size_t sphere_index = surface_spheres_old_index[iclosest];
				double rs = surface_spheres_sizing[sphere_index] + Lip * hclosest;

				double R = 2 * r + (rs + 2 * r * Lip) / (1 - Lip);

				get_tree_points_in_sphere(num_surface_spheres, 3, surface_spheres_sorted, 
					                      surface_spheres_tree_origin, surface_spheres_tree_right, surface_spheres_tree_left,
					                      x, R, num_neighbor_surface_spheres, neighbor_surface_spheres);
			}
	
			size_t num_neighbor_edge_spheres(0);
			size_t* neighbor_edge_spheres(0);
			if (num_edge_spheres > 0)
			{
				size_t iclosest(num_edge_spheres); double hclosest(DBL_MAX);
				get_closest_tree_point(num_edge_spheres, 3, edge_spheres_sorted,
					                   edge_spheres_tree_origin, edge_spheres_tree_right, edge_spheres_tree_left,
					                   x, iclosest, hclosest);

				size_t sphere_index = edge_spheres_old_index[iclosest];
				double rs = edge_spheres_sizing[sphere_index] + Lip * hclosest;

				double R = 2 * r + (rs + 2 * r * Lip) / (1 - Lip);

				get_tree_points_in_sphere(num_edge_spheres, 3, edge_spheres_sorted,
					                      edge_spheres_tree_origin, edge_spheres_tree_right, edge_spheres_tree_left,
					                      x, R, num_neighbor_edge_spheres, neighbor_edge_spheres);
			}

			size_t num_neighbor_corner_spheres(0);
			size_t* neighbor_corner_spheres(0);
			if (num_corner_spheres > 0)
			{
				size_t iclosest(num_corner_spheres); double hclosest(DBL_MAX);
				get_closest_tree_point(num_corner_spheres, 3, corner_spheres_sorted,
					                   corner_spheres_tree_origin, corner_spheres_tree_right, corner_spheres_tree_left,
					                   x, iclosest, hclosest);

				size_t sphere_index = corner_spheres_old_index[iclosest];
				double rs = corner_spheres_sizing[sphere_index] + Lip * hclosest;

				double R = 2 * r + (rs + 2 * r * Lip) / (1 - Lip);

				get_tree_points_in_sphere(num_corner_spheres, 3, corner_spheres_sorted,
					                      corner_spheres_tree_origin, corner_spheres_tree_right, corner_spheres_tree_left,
					                      x, R, num_neighbor_corner_spheres, neighbor_corner_spheres);
			}

			size_t num_spheres(num_neighbor_seeds + num_neighbor_interior_spheres + num_neighbor_surface_spheres + num_neighbor_edge_spheres + num_neighbor_corner_spheres);
			double* spheres = new double[num_spheres * 4];
			for (size_t i = 0; i < num_neighbor_seeds; i++)
			{
				size_t seed_index = seeds_old_index[neighbor_seeds[i]];
				spheres[i * 4] = seeds[seed_index * 3];
				spheres[i * 4 + 1] = seeds[seed_index * 3 + 1];
				spheres[i * 4 + 2] = seeds[seed_index * 3 + 2];
				spheres[i * 4 + 3] = seeds_sizing[seed_index];
			}
			num_spheres = num_neighbor_seeds;

			for (size_t i = 0; i < num_neighbor_interior_spheres; i++)
			{
				size_t sphere_index = neighbor_interior_spheres[i];
				spheres[num_spheres * 4 + i * 4] = interior_spheres[sphere_index * 4];
				spheres[num_spheres * 4 + i * 4 + 1] = interior_spheres[sphere_index * 4 + 1];
				spheres[num_spheres * 4 + i * 4 + 2] = interior_spheres[sphere_index * 4 + 2];
				spheres[num_spheres * 4 + i * 4 + 3] = interior_spheres[sphere_index * 4 + 3];
			}
			num_spheres += num_neighbor_interior_spheres;

			for (size_t i = 0; i < num_neighbor_surface_spheres; i++)
			{
				size_t sphere_index = surface_spheres_old_index[neighbor_surface_spheres[i]];
				spheres[num_spheres * 4 + i * 4] = surface_spheres[sphere_index * 3];
				spheres[num_spheres * 4 + i * 4 + 1] = surface_spheres[sphere_index * 3 + 1];
				spheres[num_spheres * 4 + i * 4 + 2] = surface_spheres[sphere_index * 3 + 2];
				spheres[num_spheres * 4 + i * 4 + 3] = surface_spheres_sizing[sphere_index];
			}
			if (neighbor_surface_spheres != 0) delete[] neighbor_surface_spheres;
			num_spheres += num_neighbor_surface_spheres;

			for (size_t i = 0; i < num_neighbor_edge_spheres; i++)
			{
				size_t sphere_index = edge_spheres_old_index[neighbor_edge_spheres[i]];
				spheres[num_spheres * 4 + i * 4] = edge_spheres[sphere_index * 3];
				spheres[num_spheres * 4 + i * 4 + 1] = edge_spheres[sphere_index * 3 + 1];
				spheres[num_spheres * 4 + i * 4 + 2] = edge_spheres[sphere_index * 3 + 2];
				spheres[num_spheres * 4 + i * 4 + 3] = edge_spheres_sizing[sphere_index];
			}
			if (neighbor_edge_spheres != 0) delete[] neighbor_edge_spheres;
			num_spheres += num_neighbor_edge_spheres;

			for (size_t i = 0; i < num_neighbor_corner_spheres; i++)
			{
				size_t sphere_index = corner_spheres_old_index[neighbor_corner_spheres[i]];
				spheres[num_spheres * 4 + i * 4] = corner_spheres[sphere_index * 3];
				spheres[num_spheres * 4 + i * 4 + 1] = corner_spheres[sphere_index * 3 + 1];
				spheres[num_spheres * 4 + i * 4 + 2] = corner_spheres[sphere_index * 3 + 2];
				spheres[num_spheres * 4 + i * 4 + 3] = corner_spheres_sizing[sphere_index];
			}
			if (neighbor_corner_spheres != 0) delete[] neighbor_corner_spheres;
			num_spheres += num_neighbor_corner_spheres;

			size_t num_darts(20);
			double* dart = new double[3];
			num_new_spheres[iactive] = 0;
			new_spheres[iactive] = new double[num_darts * 4];
			for (size_t idart = 0; idart < num_darts; idart++)
			{
				rsamplers[thread_id]->sample_uniformly_from_unit_sphere(dart, 3);
				double u = 1.0 + rsamplers[thread_id]->generate_uniform_random_number();
				for (size_t idim = 0; idim < 3; idim++)
				{
					double dx = u * r * dart[idim];
					dart[idim] = x[idim] + dx;
				}

				bool good_dart(true); 
				double closest_dst_sq(0.0); size_t closest_seed(seed_index);
				for (size_t idim = 0; idim < 3; idim++)
				{
					double dx = dart[idim] - x[idim];
					closest_dst_sq += dx * dx;
				}

				for (size_t isphere = 0; isphere < num_spheres; isphere++)
				{
					double hsq(0.0);
					for (size_t idim = 0; idim < 3; idim++)
					{
						double dx = dart[idim] - spheres[isphere * 4 + idim];
						hsq += dx * dx;
					}
					if (isphere < num_neighbor_seeds + num_neighbor_interior_spheres && hsq < closest_dst_sq)
					{
						closest_dst_sq = hsq;
						if (isphere < num_neighbor_seeds) closest_seed = seeds_old_index[neighbor_seeds[isphere]];
						else                              closest_seed = num_seeds + neighbor_interior_spheres[isphere - num_neighbor_seeds];
					}
					if (hsq < spheres[isphere * 4 + 3] * spheres[isphere * 4 + 3])
					{
						good_dart = false;
						break;
					}
				}

				for (size_t isphere = 0; isphere < num_new_spheres[iactive]; isphere++)
				{
					double hsq(0.0);
					for (size_t idim = 0; idim < 3; idim++)
					{
						double dx = dart[idim] - new_spheres[iactive][isphere * 4 + idim];
						hsq += dx * dx;
					}
					if (hsq < new_spheres[iactive][isphere * 4 + 3] * new_spheres[iactive][isphere * 4 + 3])
					{
						good_dart = false; // dart is inside a new sphere around the current seed
					}
				}

				if (seed_index < num_seeds && closest_seed < num_seeds && 
					seeds_region_id[seed_index] != seeds_region_id[closest_seed]) good_dart = false;
				
				if (seed_index >= num_seeds && closest_seed < num_seeds && 
					interior_spheres_region_id[seed_index - num_seeds] != seeds_region_id[closest_seed]) good_dart = false;

				if (seed_index >= num_seeds && closest_seed >= num_seeds &&
					interior_spheres_region_id[seed_index - num_seeds] != interior_spheres_region_id[closest_seed - num_seeds]) good_dart = false;

				if (seed_index < num_seeds && closest_seed >= num_seeds &&
					seeds_region_id[seed_index] != interior_spheres_region_id[closest_seed - num_seeds]) good_dart = false;

				if (!good_dart) continue;

				new_spheres[iactive][num_new_spheres[iactive] * 4] = dart[0];
				new_spheres[iactive][num_new_spheres[iactive] * 4 + 1] = dart[1];
				new_spheres[iactive][num_new_spheres[iactive] * 4 + 2] = dart[2];

				double sz(0.0);
				if (closest_seed < num_seeds) sz = seeds_sizing[closest_seed] + Lip * sqrt(closest_dst_sq);
				else                          sz = interior_spheres[(closest_seed - num_seeds) * 4 + 3] + Lip * sqrt(closest_dst_sq);
				sz = fmin(sz, rmax);

				if (num_sizing_points > 0)
				{
					size_t iclosest(0); double hclosest(DBL_MAX);
					get_closest_tree_point(num_sizing_points, 3, sizing_points_sorted,
							               sizing_points_tree_origin, sizing_points_tree_right, sizing_points_tree_left,
							               dart, iclosest, hclosest);
					sz = fmin(sz, sizing_value[sizing_points_old_index[iclosest]]);
				}

				new_spheres[iactive][num_new_spheres[iactive] * 4 + 3] = sz;

				num_new_spheres[iactive]++;
			} // end darts

			if (neighbor_interior_spheres != 0) delete[] neighbor_interior_spheres;
			delete[] x; delete[] dart; delete[] spheres; delete[] neighbor_seeds;
		}
		
		spokes_time += timer.report_timing();

		timer.reset_timer();

		// add new seeds, update front_indeces and recreate seeds_tree
		size_t num_new_seeds(0);
		for (size_t iactive = 0; iactive < num_active_seeds; iactive++) num_new_seeds += num_new_spheres[iactive];

		if (num_new_seeds == 0 && num_inactive_seeds == 0)
		{
			if (num_interior_spheres > 0)
			{
				size_t num_all_seeds = num_seeds + num_interior_spheres;
				double* new_seeds = new double[num_all_seeds * 3];
				size_t* new_seeds_region_id = new size_t[num_all_seeds];
				double* new_seeds_sizing = new double[num_all_seeds];

				for (size_t i = 0; i < num_seeds * 3; i++) new_seeds[i] = seeds[i];
				for (size_t i = 0; i < num_seeds; i++) new_seeds_region_id[i] = seeds_region_id[i];
				for (size_t i = 0; i < num_seeds; i++) new_seeds_sizing[i] = seeds_sizing[i];

				for (size_t i = 0; i < num_interior_spheres; i++)
				{
					for (size_t idim = 0; idim < 3; idim++)  new_seeds[(num_seeds + i) * 3 + idim] = interior_spheres[i * 4 + idim];
					new_seeds_sizing[num_seeds + i] = interior_spheres[i * 4 + 3];
					new_seeds_region_id[num_seeds + i] = interior_spheres_region_id[i];
				}
				delete[] seeds; delete[] seeds_region_id; delete[] seeds_sizing;
				num_seeds = num_all_seeds;
				seeds = new_seeds; seeds_region_id = new_seeds_region_id; seeds_sizing = new_seeds_sizing;

				// rebuild seeds tree
				delete[] seeds_tree_right; delete[] seeds_tree_left;
				seeds_tree_right = new size_t[num_seeds];
				seeds_tree_left = new size_t[num_seeds];
				build_balanced_kd_tree(num_seeds, 3, seeds, seeds_tree_origin, seeds_tree_right, seeds_tree_left);

				delete[] seeds_sorted; delete[] seeds_old_index; delete[] seeds_new_index;
				seeds_sorted = new double[num_seeds * 3];
				seeds_old_index = new size_t[num_seeds];
				seeds_new_index = new size_t[num_seeds];
				re_enumerate_points_for_better_memory_access(num_seeds, 3, seeds,
					                                         seeds_tree_origin, seeds_tree_right, seeds_tree_left,
					                                         seeds_sorted, seeds_old_index, seeds_new_index);
			}
			delete[] active_front_seeds;
			for (size_t iactive = 0; iactive < num_active_seeds; iactive++) delete[] new_spheres[iactive];
			delete[] new_spheres; delete[] num_new_spheres;
			break; //Spoke darts terminates here
		}
		
		if (num_new_seeds == 0)
		{
			delete[] active_front_seeds;
			for (size_t iactive = 0; iactive < num_active_seeds; iactive++) delete[] new_spheres[iactive];
			delete[] new_spheres; delete[] num_new_spheres;

			num_front_seeds = num_inactive_seeds;
			front_indices = inactive_front_seeds;
			continue;
		}

		if (num_seeds + num_interior_spheres + num_new_seeds > 1.2 * num_seeds_reported)
		{
			num_seeds_reported = num_seeds + num_interior_spheres + num_new_seeds;
			vcm_cout << "  * Number of seeds " << num_seeds_reported << std::endl;
		}

		// update front indices
		num_front_seeds = num_inactive_seeds + num_new_seeds;
		front_indices = new size_t[num_front_seeds];
		for (size_t i = 0; i < num_inactive_seeds; i++) front_indices[i] = inactive_front_seeds[i];
		for (size_t i = 0; i < num_new_seeds; i++) front_indices[num_inactive_seeds + i] = num_seeds + num_interior_spheres + i;
		delete[] inactive_front_seeds;

		if (num_interior_spheres + num_new_seeds < interior_spheres_cap)
		{
			for (size_t iactive = 0; iactive < num_active_seeds; iactive++)
			{
				size_t seed_index = active_front_seeds[iactive];
				for (size_t ii = 0; ii < num_new_spheres[iactive]; ii++)
				{
					for (size_t idim = 0; idim <= 3; idim++) interior_spheres[num_interior_spheres * 4 + idim] = new_spheres[iactive][ii * 4 + idim];
					if (seed_index < num_seeds) interior_spheres_region_id[num_interior_spheres] = seeds_region_id[seed_index];
					else                        interior_spheres_region_id[num_interior_spheres] = interior_spheres_region_id[seed_index - num_seeds];

					num_interior_spheres++;
				}
				delete[] new_spheres[iactive];
			}
			delete[] active_front_seeds;
			delete[] new_spheres; delete[] num_new_spheres;
			continue;
		}

		size_t num_all_seeds = num_seeds + num_interior_spheres + num_new_seeds;
		double* new_seeds = new double[num_all_seeds * 3];
		size_t* new_seeds_region_id = new size_t[num_all_seeds];
		double* new_seeds_sizing = new double[num_all_seeds];

		for (size_t i = 0; i < num_seeds * 3; i++) new_seeds[i] = seeds[i];
		for (size_t i = 0; i < num_seeds; i++) new_seeds_region_id[i] = seeds_region_id[i];
		for (size_t i = 0; i < num_seeds; i++) new_seeds_sizing[i] = seeds_sizing[i];

		for (size_t i = 0; i < num_interior_spheres; i++)
		{
			for (size_t idim = 0; idim < 3; idim++)  new_seeds[(num_seeds + i) * 3 + idim] = interior_spheres[i * 4 + idim];
			new_seeds_sizing[num_seeds + i] = interior_spheres[i * 4 + 3];
			new_seeds_region_id[num_seeds + i] = interior_spheres_region_id[i];
		}

		// save interior spheres here 
		save_spheres_csv("interior_spheres.csv", 3, num_interior_spheres, interior_spheres, 0, 0);
		
		num_new_seeds = num_interior_spheres;
		for (size_t iactive = 0; iactive < num_active_seeds; iactive++)
		{
			size_t seed_index = active_front_seeds[iactive];
			for (size_t ii = 0; ii < num_new_spheres[iactive]; ii++)
			{
				for (size_t idim = 0; idim < 3; idim++) new_seeds[num_seeds * 3 + num_new_seeds * 3 + idim] = new_spheres[iactive][ii * 4 + idim];
				new_seeds_sizing[num_seeds + num_new_seeds] = new_spheres[iactive][ii * 4 + 3];
				
				if (seed_index < num_seeds) new_seeds_region_id[num_seeds + num_new_seeds] = seeds_region_id[seed_index];
				else                        new_seeds_region_id[num_seeds + num_new_seeds] = interior_spheres_region_id[seed_index - num_seeds];
			
				num_new_seeds++;
			}
			delete[] new_spheres[iactive];
		}

		num_interior_spheres = 0;

		delete[] active_front_seeds;
		delete[] new_spheres; delete[] num_new_spheres;

		delete[] seeds; delete[] seeds_region_id; delete[] seeds_sizing;

		num_seeds += num_new_seeds;
		seeds = new_seeds; seeds_region_id = new_seeds_region_id; seeds_sizing = new_seeds_sizing;
		
		// rebuild seeds tree
		delete[] seeds_tree_right; delete[] seeds_tree_left;
		seeds_tree_right = new size_t[num_seeds];
		seeds_tree_left = new size_t[num_seeds];
		build_balanced_kd_tree(num_seeds, 3, seeds, seeds_tree_origin, seeds_tree_right, seeds_tree_left);

		delete[] seeds_sorted; delete[] seeds_old_index; delete[] seeds_new_index;
		seeds_sorted = new double[num_seeds * 3];
		seeds_old_index = new size_t[num_seeds];
		seeds_new_index = new size_t[num_seeds];
		re_enumerate_points_for_better_memory_access(num_seeds, 3, seeds,
			                                         seeds_tree_origin, seeds_tree_right, seeds_tree_left,
			                                         seeds_sorted, seeds_old_index, seeds_new_index);

		update_seeds_time += timer.report_timing();

		if (num_seeds > 2 * num_lip)
		{
			timer.reset_timer();

			double* seed_sizing_sorted = new double[num_seeds];
			for (size_t i = 0; i < num_seeds; i++) seed_sizing_sorted[i] = seeds_sizing[seeds_old_index[i]];
			impose_lipschitz_continuity(num_seeds, 3, seeds_sorted, seed_sizing_sorted, num_threads, Lip, seeds_tree_origin, seeds_tree_right, seeds_tree_left);
			for (size_t i = 0; i < num_seeds; i++) seeds_sizing[seeds_old_index[i]] = seed_sizing_sorted[i];
			delete[] seed_sizing_sorted;
			num_lip = num_seeds;
	
			lip_time += timer.report_timing();
		}
	}

	if (true)
	{
		timer.reset_timer();

		double* seed_sizing_sorted = new double[num_seeds];
		for (size_t i = 0; i < num_seeds; i++) seed_sizing_sorted[i] = seeds_sizing[seeds_old_index[i]];
		impose_lipschitz_continuity(num_seeds, 3, seeds_sorted, seed_sizing_sorted, num_threads, Lip, seeds_tree_origin, seeds_tree_right, seeds_tree_left);
		for (size_t i = 0; i < num_seeds; i++) seeds_sizing[seeds_old_index[i]] = seed_sizing_sorted[i];
		delete[] seed_sizing_sorted;
		
		lip_time += timer.report_timing();
	}

	vcm_cout << "  * Number of seeds " << num_seeds << std::endl;
	vcm_cout << "  * Number of iterations = " << nit << std::endl;
	vcm_cout << "  * Deactivation time = " << deactivating_time << std::endl;
	vcm_cout << "  * Spoke throwing time = " << spokes_time << std::endl;
	vcm_cout << "  * Lipschitz time = " << lip_time << std::endl;
	vcm_cout << "  * Seed update time = " << update_seeds_time << std::endl;

	delete[] seeds_tree_right; delete[] seeds_tree_left;
	delete[] seeds_sorted; delete[] seeds_old_index; delete[] seeds_new_index;

	delete[] interior_spheres;
	delete[] interior_spheres_region_id;


	if (num_surface_spheres > 0)
	{
		delete[] surface_spheres_tree_right; delete[] surface_spheres_tree_left;
		delete[] surface_spheres_sorted; delete[] surface_spheres_old_index; delete[] surface_spheres_new_index;
	}

	if (num_edge_spheres > 0)
	{
		delete[] edge_spheres_tree_right; delete[] edge_spheres_tree_left;
		delete[] edge_spheres_sorted; delete[] edge_spheres_old_index; delete[] edge_spheres_new_index;
	}

	if (num_corner_spheres > 0)
	{
		delete[] corner_spheres_tree_right; delete[] corner_spheres_tree_left;
		delete[] corner_spheres_sorted; delete[] corner_spheres_old_index; delete[] corner_spheres_new_index;
	}

	if (num_sizing_points > 0)
	{
		delete[] sizing_points_tree_right; delete[] sizing_points_tree_left;
		delete[] sizing_points_sorted; delete[] sizing_points_old_index; delete[] sizing_points_new_index;
	}

	for (size_t ithread = 0; ithread < num_threads; ithread++) delete rsamplers[ithread];

	timer.reset_timer();

	
	return 0;
	#pragma endregion
}

int MeshingVoronoiMesher::generate_3d_voronoi_mesh(int num_threads, size_t num_seeds, double* seeds, size_t* seed_region_id, double* seed_sizing,
	                                        size_t &num_vertices, double* &vertices, size_t &num_faces, size_t** &faces)
{
    #pragma region Generate a general 3d Voronoi Mesh:

	size_t tree_origin(SIZE_MAX);
	size_t* tree_right = new size_t[num_seeds]; size_t* tree_left = new size_t[num_seeds];
	build_balanced_kd_tree(num_seeds, 3, seeds, tree_origin, tree_right, tree_left);

	// enumerate 
	double* seeds_sorted = new double[num_seeds * 3]; size_t* seed_old_index = new size_t[num_seeds]; size_t* seed_new_index = new size_t[num_seeds];
	re_enumerate_points_for_better_memory_access(num_seeds, 3, seeds, tree_origin, tree_right, tree_left, seeds_sorted, seed_old_index, seed_new_index);

	bool* ghost_seed = new bool[num_seeds];
	for (int iseed = 0; iseed < num_seeds; iseed++)
	{
		size_t old_seed_index = seed_old_index[iseed];
		if (seed_region_id[old_seed_index] == 0) ghost_seed[iseed] = true;
		else                                     ghost_seed[iseed] = false;
	}

	size_t* num_cell_corners = new size_t[num_seeds];
	double** cell_corners = new double*[num_seeds];
	size_t** cell_corner_neighbors = new size_t*[num_seeds];
	size_t** cell_corner_seeds = new size_t*[num_seeds];
	size_t** cell_corner_indices = new size_t*[num_seeds]; // seed and local indices

#if defined USE_OPEN_MP
	omp_set_num_threads(num_threads);
#pragma omp parallel for
#endif
	for (int iseed = 0; iseed < num_seeds; iseed++)
	{
		size_t old_seed_index = seed_old_index[iseed];
		if (seed_region_id[old_seed_index] == 0) continue;

		get_3d_Voronoi_cell(num_seeds, seeds_sorted, tree_origin, tree_right, tree_left, 
			               iseed, 2.0 * seed_sizing[old_seed_index], ghost_seed,
			               num_cell_corners[old_seed_index], cell_corners[old_seed_index], 
			               cell_corner_neighbors[old_seed_index], cell_corner_seeds[old_seed_index], cell_corner_indices[old_seed_index]);

		for (size_t icorner = 0; icorner < num_cell_corners[old_seed_index]; icorner++)
		{
			for (size_t i = 0; i < 3; i++) cell_corner_seeds[old_seed_index][icorner * 3 + i] = seed_old_index[cell_corner_seeds[old_seed_index][icorner * 3 + i]];
			cell_corner_indices[old_seed_index][icorner * 3] = seed_old_index[cell_corner_indices[old_seed_index][icorner * 3]];
		}
	}

	delete[] ghost_seed; delete[] seeds_sorted; delete[] tree_right; delete[] tree_left;
	delete[] seed_old_index; delete[] seed_new_index;
	
	for (int iseed = 0; iseed < num_seeds; iseed++)
	{
		#pragma region Tag local index of every vertex:
		if (seed_region_id[iseed] == 0) continue;

		for (size_t icorner = 0; icorner < num_cell_corners[iseed]; icorner++)
		{
			double vtx_x = cell_corners[iseed][icorner * 3 + 0];
			double vtx_y = cell_corners[iseed][icorner * 3 + 1];
			double vtx_z = cell_corners[iseed][icorner * 3 + 2];

			size_t vtx_seed_index = cell_corner_indices[iseed][icorner * 3];

			cell_corner_indices[iseed][icorner * 3 + 1] = SIZE_MAX;
			double hmin(DBL_MAX);
			for (size_t jcorner = 0; jcorner < num_cell_corners[vtx_seed_index]; jcorner++)
			{
				double dx = cell_corners[vtx_seed_index][jcorner * 3 + 0] - vtx_x;
				double dy = cell_corners[vtx_seed_index][jcorner * 3 + 1] - vtx_y;
				double dz = cell_corners[vtx_seed_index][jcorner * 3 + 2] - vtx_z;
				double dst(sqrt(dx * dx + dy * dy + dz * dz));
				if (dst < hmin)
				{
					hmin = dst;
					cell_corner_indices[iseed][icorner * 3 + 1] = jcorner;
				}
			}
			if (hmin > 1E-6)
			{
				// unweld this corner // NEED TO VERIFY MESH CONFORMITY AFTER GENERATING FACES
				cell_corner_indices[iseed][icorner * 3] = iseed;
				cell_corner_indices[iseed][icorner * 3 + 1] = icorner;
				//vcm_cout << "WARNING: (MeshingVoronoiMesher.cpp): Distance between two welded vertices = " << hmin << ", seed index = " << iseed << ", icorner = " << icorner << std::endl;
			}
		}
		#pragma endregion
	}

	// making sure that each corner points to the chosen one
	for (int iseed = 0; iseed < num_seeds; iseed++)
	{
		if (seed_region_id[iseed] == 0) continue;
		for (size_t icorner = 0; icorner < num_cell_corners[iseed]; icorner++)
		{
			size_t vtx_seed = cell_corner_indices[iseed][icorner * 3];
			size_t vtx_corner = cell_corner_indices[iseed][icorner * 3 + 1];
			while (!(cell_corner_indices[vtx_seed][vtx_corner * 3] == vtx_seed &&
				     cell_corner_indices[vtx_seed][vtx_corner * 3 + 1] == vtx_corner))
			{
				size_t new_vtx_seed = cell_corner_indices[vtx_seed][vtx_corner * 3];
				size_t new_vtx_corner = cell_corner_indices[vtx_seed][vtx_corner * 3 + 1];
				vtx_seed = new_vtx_seed; vtx_corner = new_vtx_corner;
			}
			cell_corner_indices[iseed][icorner * 3] = vtx_seed;
			cell_corner_indices[iseed][icorner * 3 + 1] = vtx_corner;
		}
	}

	// count meshed cells and vertices
	size_t num_meshed_cells(0); num_vertices = 0;
	for (int iseed = 0; iseed < num_seeds; iseed++)
	{
		if (seed_region_id[iseed] == 0) continue;
		num_meshed_cells++;
		for (size_t icorner = 0; icorner < num_cell_corners[iseed]; icorner++)
		{
			if (cell_corner_indices[iseed][icorner * 3] != iseed) continue; // vertex is assigned to another cell
			if (cell_corner_indices[iseed][icorner * 3 + 1] != icorner) continue; // vertex is assigned to another corner
			cell_corner_indices[iseed][icorner * 3 + 2] = num_vertices;
			num_vertices++;
		}
	}

	vertices = new double[num_vertices * 3]; size_t vtx_index(0);
	for (int iseed = 0; iseed < num_seeds; iseed++)
	{
		#pragma region Extract vertices:
		if (seed_region_id[iseed] == 0) continue;
		num_meshed_cells++;
		for (size_t icorner = 0; icorner < num_cell_corners[iseed]; icorner++)
		{
			if (cell_corner_indices[iseed][icorner * 3] != iseed) continue; // vertex is assigned to another cell
			if (cell_corner_indices[iseed][icorner * 3 + 1] != icorner) continue; // vertex is assigned to another corner
			vertices[vtx_index * 3] = cell_corners[iseed][icorner * 3];
			vertices[vtx_index * 3 + 1] = cell_corners[iseed][icorner * 3 + 1];
			vertices[vtx_index * 3 + 2] = cell_corners[iseed][icorner * 3 + 2];
			vtx_index++;
		}
		delete[] cell_corners[iseed];
		#pragma endregion
	}
	delete[] cell_corners;

	// count non_degenerate faces
	num_faces= 0;
	std::vector<size_t> face_corners;
	for (int iseed = 0; iseed < num_seeds; iseed++)
	{
		#pragma region Count Faces
		if (seed_region_id[iseed] == 0) continue;
		std::set<size_t> neighbors;
		for (size_t icorner = 0; icorner < num_cell_corners[iseed]; icorner++)
		{
			for (size_t i = 0; i < 3; i++) neighbors.insert(cell_corner_seeds[iseed][icorner * 3 + i]);
		}

		size_t num_corners = num_cell_corners[iseed];
		size_t num_neigbors(neighbors.size());
		for (std::set<size_t>::iterator it = neighbors.begin(); it != neighbors.end(); ++it)
		{
			size_t neighbor_seed = *it;
			if (neighbor_seed < iseed && seed_region_id[neighbor_seed] != 0) continue; // this face is assigned with another seed

			size_t first_corner(num_corners);
			for (size_t icorner = 0; icorner < num_corners; icorner++)
			{
				for (size_t i = 0; i < 3; i++)
				{
					if (cell_corner_seeds[iseed][icorner * 3 + i] == neighbor_seed)
					{
						first_corner = icorner;
						break;
					}
				}
				if (first_corner != num_corners) break;
			}

			size_t num_face_corners(0);
			face_corners.clear();
			size_t current_corner(first_corner);

			size_t num_iter(0);
			while (true)
			{
				if (num_iter > 100 * num_face_corners)
				{
					vcm_cout << "WARNING::Singular face between seeds " << iseed << " and " << neighbor_seed << std::endl;
					break;
				}

				num_iter++;

				size_t vtx_seed = cell_corner_indices[iseed][current_corner * 3];
				size_t vtx_corner = cell_corner_indices[iseed][current_corner * 3 + 1];
				size_t vtx_index = cell_corner_indices[vtx_seed][vtx_corner * 3 + 2];

				bool redundant(false);
				for (size_t k = 0; k < num_face_corners; k++)
				{
					if (vtx_index == face_corners[k]) redundant = true;
				}

				if (num_face_corners > 2 && vtx_index == face_corners[0]) break;

				if (num_face_corners == 0 || !redundant)
				{
					face_corners.push_back(vtx_index);
					num_face_corners++;
				}
				 
				for (size_t i = 0; i < 3; i++)
				{
					if (cell_corner_seeds[iseed][current_corner * 3 + i] == neighbor_seed)
					{
						current_corner = cell_corner_neighbors[iseed][current_corner * 3 + i];
						break;
					}
				}
				if (current_corner == first_corner) break;
			}
			if (num_face_corners >= 3) num_faces++;
		}
		#pragma endregion
	}

	faces = new size_t*[num_faces]; size_t face_index(0);
	for (int iseed = 0; iseed < num_seeds; iseed++)
	{
		#pragma region Extract Faces
		if (seed_region_id[iseed] == 0) continue;
		std::set<size_t> neighbors;
		for (size_t icorner = 0; icorner < num_cell_corners[iseed]; icorner++)
		{
			for (size_t i = 0; i < 3; i++) neighbors.insert(cell_corner_seeds[iseed][icorner * 3 + i]);
		}

		size_t num_corners = num_cell_corners[iseed];
		size_t num_neigbors(neighbors.size());
		for (std::set<size_t>::iterator it = neighbors.begin(); it != neighbors.end(); ++it)
		{
			size_t neighbor_seed = *it;
			if (neighbor_seed < iseed && seed_region_id[neighbor_seed] != 0) continue; // this face is assigned with another seed

			size_t first_corner(num_corners);
			for (size_t icorner = 0; icorner < num_corners; icorner++)
			{
				for (size_t i = 0; i < 3; i++)
				{
					if (cell_corner_seeds[iseed][icorner * 3 + i] == neighbor_seed)
					{
						first_corner = icorner;
						break;
					}
				}
				if (first_corner != num_corners) break;
			}

			size_t num_face_corners(0);
			face_corners.clear();
			size_t current_corner(first_corner);
			size_t num_iter(0);
			while (true)
			{
				if (num_iter > 100 * num_face_corners)
				{
					//vcm_cout << "WARNING::Singular face between seeds " << iseed << " and " << neighbor_seed << std::endl;
					break;
				}

				num_iter++;

				size_t vtx_seed = cell_corner_indices[iseed][current_corner * 3];
				size_t vtx_corner = cell_corner_indices[iseed][current_corner * 3 + 1];
				size_t vtx_index = cell_corner_indices[vtx_seed][vtx_corner * 3 + 2];

				bool redundant(false);
				for (size_t k = 0; k < num_face_corners; k++)
				{
					if (vtx_index == face_corners[k]) redundant = true;
				}

				if (num_face_corners > 2 && vtx_index == face_corners[0]) break;

				if (num_face_corners == 0 || !redundant)
				{
					face_corners.push_back(vtx_index);
					num_face_corners++;
				}

				for (size_t i = 0; i < 3; i++)
				{
					if (cell_corner_seeds[iseed][current_corner * 3 + i] == neighbor_seed)
					{
						current_corner = cell_corner_neighbors[iseed][current_corner * 3 + i];
						break;
					}
				}
				if (current_corner == first_corner) break;
			}

			if (num_face_corners >= 3)
			{
				faces[face_index] = new size_t[num_face_corners + 3];
				faces[face_index][0] = num_face_corners;
				for (size_t i = 0; i < num_face_corners; i++) faces[face_index][i + 1] = face_corners[i];
				faces[face_index][num_face_corners + 1] = iseed;
				faces[face_index][num_face_corners + 2] = neighbor_seed;
				face_index++;
			}
		}
		delete[] cell_corner_seeds[iseed];
		delete[] cell_corner_neighbors[iseed];
		#pragma endregion
	}

	delete[] num_cell_corners; delete[] cell_corner_seeds; delete[] cell_corner_neighbors;

	for (int iseed = 0; iseed < num_seeds; iseed++)
	{
		if (seed_region_id[iseed] == 0) continue;
		delete[] cell_corner_indices[iseed];
	}
	delete[] cell_corner_indices;

	// adjust regions id
	
 	return 0;
	#pragma endregion
}

int MeshingVoronoiMesher::get_3d_Voronoi_cell(size_t num_seeds, double* seeds, size_t tree_origin, size_t* tree_right, size_t* tree_left,
	                                   size_t seed_index, double r, bool* ghost_seed,
	                                   size_t& num_cell_corners, double* &cell_corners, size_t* &cell_corner_neighbors, 
	                                   size_t*& cell_corner_seeds, size_t*& cell_corner_indices)
{
	#pragma region Get 3d Voronoi Cell:
	double st_time = 0.0; // omp_get_wtime();
	double* x = new double[3];
	for (size_t idim = 0; idim < 3; idim++) x[idim] = seeds[seed_index * 3 + idim];

	double* xo = new double[3]; double* n = new double[3];

	while (true)
	{
		r *= 2.0;

		size_t num_local_neighbors(0); size_t* local_neighbors(0);
		get_tree_points_in_sphere(num_seeds, 3, seeds, tree_origin, tree_right, tree_left, 
			                      x, r, num_local_neighbors, local_neighbors);

		if (num_local_neighbors == 0) continue;

		double* local_seeds = new double[num_local_neighbors * 3];
		size_t iseed(0);
		for (size_t ineighbor = 0; ineighbor < num_local_neighbors; ineighbor++)
		{
			size_t neighbor_seed_index = local_neighbors[ineighbor];
			for (size_t idim = 0; idim < 3; idim++) local_seeds[ineighbor * 3 + idim] = seeds[neighbor_seed_index * 3 + idim];
		}
		
		size_t corners_cap;
		construct_initial_polytope(x, 10 * r, num_cell_corners, corners_cap, cell_corners, cell_corner_neighbors, cell_corner_seeds);

		for (size_t icorner = 0; icorner < num_cell_corners; icorner++)
		{
			// get closest point from tree
			if (cell_corners[icorner * 3] == DBL_MAX) continue;

			double lambda_min(DBL_MAX); size_t closest_seed(num_local_neighbors);
			for (size_t ipoint = 0; ipoint < num_local_neighbors; ipoint++)
			{
				double xy_dot_xy(0.0), xv_dot_xy(0.0);
				for (size_t idim = 0; idim < 3; idim++)
				{
					double dxy = local_seeds[ipoint * 3 + idim] - x[idim];
					double dxv = cell_corners[icorner * 3 + idim] - x[idim];
					xy_dot_xy += dxy * dxy;
					xv_dot_xy += dxy * dxv;
				}
				if (xv_dot_xy < 1E-10) continue;
				double lambda = xy_dot_xy / xv_dot_xy;
				if (lambda < lambda_min)
				{
					lambda_min = lambda;
					closest_seed = ipoint;
				}
			}
			if (lambda_min > 2.0) continue; // a valid voronoi vertex

			double norm(0.0);
			for (size_t idim = 0; idim < 3; idim++)
			{
				xo[idim] = 0.5 * (local_seeds[closest_seed * 3 + idim] + x[idim]);
				n[idim] = local_seeds[closest_seed * 3 + idim] - x[idim];
				norm += n[idim] * n[idim];
			}
			norm = sqrt(norm);
			for (size_t idim = 0; idim < 3; idim++) n[idim] /= norm;

			double xi_dot_n(0.0);
			for (size_t idim = 0; idim < 3; idim++) xi_dot_n += n[idim] * (cell_corners[icorner * 3 + idim] - xo[idim]);
			if (xi_dot_n <= 1E-10) continue; // a valid voronoi vertex, accounting for the same round-off error in the trimming method

			trim_vertex(closest_seed, xo, n, icorner, icorner, num_cell_corners, corners_cap, cell_corners, cell_corner_neighbors, cell_corner_seeds);
		}
		
		double r_v(0.0);
		for (size_t icorner = 0; icorner < num_cell_corners; icorner++)
		{
			if (cell_corners[icorner * 3 + 0] == DBL_MAX) continue;

			double dst(0.0);
			for (size_t idim = 0; idim < 3; idim++)
			{
				double dx = cell_corners[icorner * 3 + idim] - x[idim];
				dst += dx * dx;
			}
			if (dst > r_v) r_v = dst;
		}
		r_v = sqrt(r_v);
		
		if (2 * r_v < r + 1E-10)
		{
			size_t num_good_corners(0);
			size_t* corner_new_index = new size_t[num_cell_corners];
			for (size_t icorner = 0; icorner < num_cell_corners; icorner++)
			{
				if (cell_corners[icorner * 3 + 0] == DBL_MAX) continue;
				for (size_t i = 0; i < 3; i++) cell_corner_seeds[icorner * 3 + i] = local_neighbors[cell_corner_seeds[icorner * 3 + i]];
				corner_new_index[icorner] = num_good_corners;
				num_good_corners++;
			}

			double* good_corners = new double[num_good_corners * 3];
			size_t* good_corner_neighbors = new size_t[num_good_corners * 3];
			size_t* good_corner_seeds = new size_t[num_good_corners * 3];

			cell_corner_indices = new size_t[num_good_corners * 3];

			size_t corner_index(0);
			for (size_t icorner = 0; icorner < num_cell_corners; icorner++)
			{
				if (cell_corners[icorner * 3 + 0] == DBL_MAX) continue;
				for (size_t i = 0; i < 3; i++) good_corners[corner_index * 3 + i] = cell_corners[icorner * 3 + i];
				for (size_t i = 0; i < 3; i++) good_corner_neighbors[corner_index * 3 + i] = corner_new_index[cell_corner_neighbors[icorner * 3 + i]];
				for (size_t i = 0; i < 3; i++) good_corner_seeds[corner_index * 3 + i] = cell_corner_seeds[icorner * 3 + i];

				// pick the neighbor with the lowest index to associate this corner to
				double dx(cell_corners[icorner * 3] - x[0]);
				double dy(cell_corners[icorner * 3 + 1] - x[1]);
				double dz(cell_corners[icorner * 3 + 2] - x[2]);
				double rv(sqrt(dx* dx + dy * dy + dz * dz));

				cell_corner_indices[corner_index * 3] = SIZE_MAX;
				for (size_t ineighbor = 0; ineighbor < num_local_neighbors; ineighbor++)
				{
					if (ghost_seed[local_neighbors[ineighbor]]) continue;

					dx = cell_corners[icorner * 3] - local_seeds[ineighbor * 3];
					dy = cell_corners[icorner * 3 + 1] - local_seeds[ineighbor * 3 + 1];
					dz = cell_corners[icorner * 3 + 2] - local_seeds[ineighbor * 3 + 2];
					double sv(sqrt(dx * dx + dy * dy + dz * dz));

					if (fabs(rv - sv) < 1E-10 && local_neighbors[ineighbor] < cell_corner_indices[corner_index * 3])
						cell_corner_indices[corner_index * 3] = local_neighbors[ineighbor];
				}
				corner_index++;
			}
			delete[] cell_corners; delete[] cell_corner_neighbors; delete[] cell_corner_seeds;
			delete[] corner_new_index;

			num_cell_corners = num_good_corners;
			cell_corners = good_corners;
			cell_corner_neighbors = good_corner_neighbors;
			cell_corner_seeds = good_corner_seeds;

			delete[] local_neighbors; delete[] local_seeds;
			break;
		}

		delete[] local_neighbors; delete[] local_seeds;
		delete[] cell_corners; delete[] cell_corner_neighbors; delete[] cell_corner_seeds;
	}
	delete[] x;  delete[] xo; delete[] n;
	return 0;
	#pragma endregion
}

int MeshingVoronoiMesher::construct_initial_polytope(double* x, double r, size_t& num_corners, size_t& corners_cap, double*& corners, size_t*& corners_neighbors, size_t*& corners_seeds)
{
	#pragma region Construct Initial Polytope:
	num_corners = 8; corners_cap = 100;

	double dst = r;

	corners = new double[corners_cap * 3];
	corners_neighbors = new size_t[corners_cap * 3];
	corners_seeds = new size_t[corners_cap * 3];
	for (size_t icorner = 0; icorner < num_corners; icorner++)
	{
		for (size_t idim = 0; idim < 3; idim++) corners[icorner * 3 + idim] = x[idim];

		for (size_t idim = 0; idim < 3; idim++) corners_neighbors[icorner * 3 + idim] = icorner;

		for (size_t idim = 0; idim < 3; idim++) corners_seeds[icorner * 3 + idim] = SIZE_MAX;
	}

	corners[0 * 3 + 0] += dst; corners[1 * 3 + 0] += dst; corners[2 * 3 + 0] += dst; corners[3 * 3 + 0] += dst;
	corners[4 * 3 + 0] -= dst; corners[5 * 3 + 0] -= dst; corners[6 * 3 + 0] -= dst; corners[7 * 3 + 0] -= dst;

	corners[2 * 3 + 1] += dst; corners[3 * 3 + 1] += dst; corners[6 * 3 + 1] += dst; corners[7 * 3 + 1] += dst;
	corners[0 * 3 + 1] -= dst; corners[1 * 3 + 1] -= dst; corners[4 * 3 + 1] -= dst; corners[5 * 3 + 1] -= dst;

	corners[0 * 3 + 2] += dst; corners[2 * 3 + 2] += dst; corners[4 * 3 + 2] += dst; corners[6 * 3 + 2] += dst;
	corners[1 * 3 + 2] -= dst; corners[3 * 3 + 2] -= dst; corners[5 * 3 + 2] -= dst; corners[7 * 3 + 2] -= dst;

	corners_neighbors[0 * 3 + 0] = 1; corners_neighbors[0 * 3 + 1] = 2; corners_neighbors[0 * 3 + 2] = 4;
	corners_neighbors[1 * 3 + 0] = 5; corners_neighbors[1 * 3 + 1] = 3; corners_neighbors[1 * 3 + 2] = 0;
	corners_neighbors[2 * 3 + 0] = 6; corners_neighbors[2 * 3 + 1] = 0; corners_neighbors[2 * 3 + 2] = 3;
	corners_neighbors[3 * 3 + 0] = 2; corners_neighbors[3 * 3 + 1] = 1; corners_neighbors[3 * 3 + 2] = 7;
	corners_neighbors[4 * 3 + 0] = 0; corners_neighbors[4 * 3 + 1] = 6; corners_neighbors[4 * 3 + 2] = 5;
	corners_neighbors[5 * 3 + 0] = 4; corners_neighbors[5 * 3 + 1] = 7; corners_neighbors[5 * 3 + 2] = 1;
	corners_neighbors[6 * 3 + 0] = 7; corners_neighbors[6 * 3 + 1] = 4; corners_neighbors[6 * 3 + 2] = 2;
	corners_neighbors[7 * 3 + 0] = 3; corners_neighbors[7 * 3 + 1] = 5; corners_neighbors[7 * 3 + 2] = 6;
	
	corners_seeds[0 * 3 + 0] -= 0; corners_seeds[0 * 3 + 1] -= 4; corners_seeds[0 * 3 + 2] -= 3;
	corners_seeds[1 * 3 + 0] -= 5; corners_seeds[1 * 3 + 1] -= 0; corners_seeds[1 * 3 + 2] -= 3;
	corners_seeds[2 * 3 + 0] -= 4; corners_seeds[2 * 3 + 1] -= 0; corners_seeds[2 * 3 + 2] -= 1;
	corners_seeds[3 * 3 + 0] -= 0; corners_seeds[3 * 3 + 1] -= 5; corners_seeds[3 * 3 + 2] -= 1;
	corners_seeds[4 * 3 + 0] -= 4; corners_seeds[4 * 3 + 1] -= 2; corners_seeds[4 * 3 + 2] -= 3;
	corners_seeds[5 * 3 + 0] -= 2; corners_seeds[5 * 3 + 1] -= 5; corners_seeds[5 * 3 + 2] -= 3;
	corners_seeds[6 * 3 + 0] -= 2; corners_seeds[6 * 3 + 1] -= 4; corners_seeds[6 * 3 + 2] -= 1;
	corners_seeds[7 * 3 + 0] -= 5; corners_seeds[7 * 3 + 1] -= 2; corners_seeds[7 * 3 + 2] -= 1;

	return 0;
	#pragma endregion
}




int MeshingVoronoiMesher::trim_vertex(size_t trimming_seed_index, double* xo, double* n, size_t corner_index, size_t old_corner,
	                           size_t& num_corners, size_t& corners_cap, double*& corners, size_t*& corners_neighbors, size_t*& corners_seeds)
{
	#pragma region Trim Polytope (Recursively):

	size_t num_old_corners(num_corners);
	size_t io(0);
	if (old_corner != corner_index)
	{
		for (size_t i = 0; i < 3; i++)
		{
			if (corners_neighbors[corner_index * 3 + i] == old_corner)
			{
				io = i + 1;
				if (io == 3) io = 0;
				break;
			}
		}
	}

	for (size_t ii = 0; ii < 3; ii++)
	{
		size_t i = io + ii;
		if (i >= 3) i -= 3;
		size_t icorner = corners_neighbors[corner_index * 3 + i];
		if (icorner == SIZE_MAX)
			continue;

		bool visited_before(false);
		for (size_t j = 0; j < 3; j++)
		{
			if (corners_neighbors[icorner * 3 + j] != SIZE_MAX) continue;
			visited_before = true;
			break;
		}
		if (visited_before) continue;


		double xj_dot_n(0.0);
		for (size_t idim = 0; idim < 3; idim++) xj_dot_n += n[idim] * (corners[icorner * 3 + idim] - xo[idim]);

		// break that edge
		corners_neighbors[corner_index * 3 + i] = SIZE_MAX;
		if (xj_dot_n > 1E-10)
		{
			trim_vertex(trimming_seed_index, xo, n, icorner, corner_index, num_corners, corners_cap, corners, corners_neighbors, corners_seeds);
		}
		else
		{
			// create a new corner
			double xi_dot_n(0.0);
			for (size_t idim = 0; idim < 3; idim++) xi_dot_n += n[idim] * (corners[corner_index * 3 + idim] - xo[idim]);
			double lambda = xi_dot_n / (xi_dot_n - xj_dot_n);

			// adjusting round-off error to preserve polytope convexity
			if (lambda < 0.0) lambda = 0.0;
			if (lambda > 1.0) lambda = 1.0;

			double* xk = new double[3];
			for (size_t idim = 0; idim < 3; idim++) xk[idim] = corners[corner_index * 3 + idim] + lambda * (corners[icorner * 3 + idim] - corners[corner_index * 3 + idim]);

			size_t* seeds = new size_t[3]; seeds[i] = corners_seeds[corner_index * 3 + i];
			size_t im = 2; if (i > 0) im = i - 1;
			seeds[im] = corners_seeds[corner_index * 3 + im];
			size_t ip = i + 1; if (ip == 3) ip = 0;
			seeds[ip] = trimming_seed_index;

			size_t* neighbors = new size_t[3]; neighbors[i] = icorner;
			neighbors[im] = SIZE_MAX; neighbors[ip] = SIZE_MAX;

			for (size_t idim = 0; idim < 3; idim++) corners[num_corners * 3 + idim] = xk[idim];
			for (size_t idim = 0; idim < 3; idim++) corners_neighbors[num_corners * 3 + idim] = neighbors[idim];
			for (size_t idim = 0; idim < 3; idim++) corners_seeds[num_corners * 3 + idim] = seeds[idim];

			delete[] xk; delete[] seeds; delete[] neighbors;

			for (size_t j = 0; j < 3; j++)
			{
				if (corners_neighbors[icorner * 3 + j] == corner_index)
				{
					corners_neighbors[icorner * 3 + j] = num_corners;
					break;
				}
			}

			num_corners++;

			if (num_corners == corners_cap)
			{
				corners_cap *= 2;
				double* new_corners = new double[corners_cap * 3];
				size_t* new_corners_neighbors = new size_t[corners_cap * 3];
				size_t* new_corners_seeds = new size_t[corners_cap * 3];
				for (size_t j = 0; j < num_corners * 3; j++) new_corners[j] = corners[j];
				for (size_t j = 0; j < num_corners * 3; j++) new_corners_neighbors[j] = corners_neighbors[j];
				for (size_t j = 0; j < num_corners * 3; j++) new_corners_seeds[j] = corners_seeds[j];
				delete[] corners; delete[] corners_neighbors; delete[] corners_seeds;
				corners = new_corners; corners_neighbors = new_corners_neighbors; corners_seeds = new_corners_seeds;
			}
		}
	}

	// deleting a vertex
	for (size_t idim = 0; idim < 3; idim++)
	{
		corners[corner_index * 3 + idim] = DBL_MAX;
		corners_neighbors[corner_index * 3 + idim] = SIZE_MAX;
		corners_seeds[corner_index * 3 + idim] = SIZE_MAX;
	}

	if (corner_index == old_corner)
	{
		for (size_t i = num_old_corners; i < num_corners; i++)
		{
			size_t im(num_corners - 1); if (i > num_old_corners) im = i - 1;
			size_t ip(i + 1); if (ip == num_corners) ip = num_old_corners;
			for (size_t j = 0; j < 3; j++)
			{
				if (corners_neighbors[i * 3 + j] == SIZE_MAX) continue;
				size_t jp = j + 1; if (jp == 3) jp = 0;
				size_t jm = 2; if (j > 0) jm = j - 1;
				corners_neighbors[i * 3 + jp] = ip;
				corners_neighbors[i * 3 + jm] = im;
				break;
			}
		}
	}
	return 0;
	#pragma endregion
}

int MeshingVoronoiMesher::save_Voronoi_tessellation(std::string file_name,
	                                                size_t num_vertices, double* vertices, size_t num_faces, size_t** faces)
{
	#pragma region Save PLY mesh:

	std::fstream file(file_name.c_str(), std::ios::out);

	file << "ply" << std::endl; // header

	file << "format ascii 1.0" << std::endl;

	file << "comment VoroCrust generated" << std::endl;

	file << "element vertex " << num_vertices << std::endl;

	file << "property float x" << std::endl;

	file << "property float y" << std::endl;

	file << "property float z" << std::endl;

	file << "element face " << num_faces << std::endl;

	file << "property list uchar int vertex_indices" << std::endl;

	file << "end_header" << std::endl;

	// Points
	for (size_t ivtx = 0; ivtx < num_vertices; ivtx++)
	{
		file << vertices[ivtx * 3] << " " << vertices[ivtx * 3 + 1] << " " << vertices[ivtx * 3 + 2] << std::endl;
	}

	// face
	size_t corner_index(0);
	for (size_t iface = 0; iface < num_faces; iface++)
	{
		size_t num_face_corners(faces[iface][0]);

		file << num_face_corners;
		for (size_t i = 1; i <= num_face_corners; i++)
		{
			size_t vtx_index = faces[iface][i];
			file << " " << vtx_index;
		}
		file << std::endl;
	}
	return 0;
	#pragma endregion
}

int MeshingVoronoiMesher::save_Voronoi_tessellation(std::string file_name,
	                                                size_t num_vertices, double* vertices, size_t num_faces, size_t** faces,
	                                                double* xo, double* edir)
{
	#pragma region Save PLY mesh:

	// count faces to be plotted

	size_t num_valid_faces(0);
	bool* valid_faces = new bool[num_faces];
	for (size_t iface = 0; iface < num_faces; iface++)
	{
		size_t num_face_corners(faces[iface][0]);
		valid_faces[iface] = true;
		for (size_t i = 1; i <= num_face_corners; i++)
		{
			size_t vtx_index = faces[iface][i];
			double dot(0.0);
			for (size_t idim = 0; idim < 3; idim++)
			{
				double dx = vertices[vtx_index * 3 + idim] - xo[idim];
				dot += dx * edir[idim];
			}
			if (dot < -1E-10)
			{
				valid_faces[iface] = false;
				break;
			}
		}
		if (valid_faces[iface]) num_valid_faces++;
	}


	std::fstream file(file_name.c_str(), std::ios::out);

	file << "ply" << std::endl; // header

	file << "format ascii 1.0" << std::endl;

	file << "comment VoroCrust generated" << std::endl;

	file << "element vertex " << num_vertices << std::endl;

	file << "property float x" << std::endl;

	file << "property float y" << std::endl;

	file << "property float z" << std::endl;

	file << "element face " << num_valid_faces << std::endl;

	file << "property list uchar int vertex_indices" << std::endl;

	file << "end_header" << std::endl;

	// Points
	for (size_t ivtx = 0; ivtx < num_vertices; ivtx++)
	{
		file << vertices[ivtx * 3] << " " << vertices[ivtx * 3 + 1] << " " << vertices[ivtx * 3 + 2] << std::endl;
	}

	// face
	size_t corner_index(0);
	for (size_t iface = 0; iface < num_faces; iface++)
	{
		if (!valid_faces[iface]) continue;

		size_t num_face_corners(faces[iface][0]);

		file << num_face_corners;
		for (size_t i = 1; i <= num_face_corners; i++)
		{
			size_t vtx_index = faces[iface][i];
			file << " " << vtx_index;
		}
		file << std::endl;
	}
	delete[] valid_faces;
	return 0;
	#pragma endregion
}

int MeshingVoronoiMesher::save_vcg_file_pflotran(std::string outFile, size_t num_seeds, double* seeds, size_t* seed_region_id,
	                                      size_t num_vertices, double* vertices, size_t num_faces, size_t** faces,
	                                      size_t* cell_num_faces, size_t** cell_faces)
{
	#pragma region Save VCG FILE for PFLOTRAN:

	vcm_cout << "\nVoroCrust::Saving VoroCrust Geomerty (.vcg) file:" << std::endl;

	// delete ghost seeds
	size_t num_meshed_cells(0);
	size_t* seed_new_index = new size_t[num_seeds];
	for (size_t iseed = 0; iseed < num_seeds; iseed++)
	{
		seed_new_index[iseed] = num_seeds;
		if (seed_region_id[iseed] == 0) continue; // a ghost seed

		seed_new_index[iseed] = num_meshed_cells;
		num_meshed_cells++;
	}

	double* face_area_height = new double[num_faces * 2];
	for (size_t iface = 0; iface < num_faces; iface++)
	{
		#pragma region Estimating Face Area and Height:
		face_area_height[iface * 2] = 0.0;
		size_t ko = faces[iface][1];
		size_t num_face_corners(faces[iface][0]);
		for (size_t i = 2; i < num_face_corners; i++)
		{
			size_t ip = i + 1;

			size_t kk = faces[iface][i];
			size_t kp = faces[iface][ip];

			double ax = vertices[kk * 3] - vertices[ko * 3];
			double ay = vertices[kk * 3 + 1] - vertices[ko * 3 + 1];
			double az = vertices[kk * 3 + 2] - vertices[ko * 3 + 2];

			double bx = vertices[kp * 3] - vertices[ko * 3];
			double by = vertices[kp * 3 + 1] - vertices[ko * 3 + 1];
			double bz = vertices[kp * 3 + 2] - vertices[ko * 3 + 2];

			double nx = ay * bz - az * by;
			double ny = az * bx - ax * bz;
			double nz = ax * by - ay * bx;

			double norm(sqrt(nx * nx + ny * ny + nz * nz));
			face_area_height[iface * 2] += 0.5 * norm;
		}

		size_t seed_i = faces[iface][num_face_corners + 1];
		size_t seed_j = faces[iface][num_face_corners + 2];
		double dx = seeds[seed_i * 3] - seeds[seed_j * 3];
		double dy = seeds[seed_i * 3 + 1] - seeds[seed_j * 3 + 1];
		double dz = seeds[seed_i * 3 + 2] - seeds[seed_j * 3 + 2];
		face_area_height[iface * 2 + 1] = 0.5 * sqrt(dx * dx + dy * dy + dz * dz);
		#pragma endregion
	}

	std::fstream file(outFile.c_str(), std::ios::out);
	file << "CELLS " << num_meshed_cells << std::endl;
	double vol(0.0), min_vol(DBL_MAX), max_vol(0.0);
	for (size_t iseed = 0; iseed < num_seeds; iseed++)
	{
		if (seed_region_id[iseed] == 0) continue;

		double cell_volume(0.0);
		for (size_t iface = 0; iface < cell_num_faces[iseed]; iface++)
		{
			size_t face_index = cell_faces[iseed][iface];
			cell_volume += face_area_height[face_index * 2] * face_area_height[face_index * 2 + 1] / 3;
		}

		file << seed_new_index[iseed] + 1 << " " << std::setprecision(16)
			 << seeds[iseed * 3] << " " << seeds[iseed * 3 + 1] << " " << seeds[iseed * 3 + 2] 
			 << " " << cell_volume << " " << seed_region_id[iseed] << std::endl;

		vol += cell_volume;
		if (cell_volume < min_vol) min_vol = cell_volume;
		if (cell_volume > max_vol) max_vol = cell_volume;
	}
	vcm_cout << "  * Total Volume = " << vol << std::endl;
	vcm_cout << "    Minimum cell volume = " << min_vol << std::endl;
	vcm_cout << "    Maximum cell volume = " << max_vol << std::endl;
	vcm_cout << "    Average cell volume = " << vol / num_meshed_cells << std::endl;

	file << "CONNECTIONS " << num_faces << std::endl;
	for (size_t iface = 0; iface < num_faces; iface++)
	{
		size_t num_face_corners = faces[iface][0];
		size_t i1 = faces[iface][num_face_corners + 1];
		size_t i2 = faces[iface][num_face_corners + 2];

		double nx = seeds[i2 * 3 + 0] - seeds[i1 * 3 + 0];
		double ny = seeds[i2 * 3 + 1] - seeds[i1 * 3 + 1];
		double nz = seeds[i2 * 3 + 2] - seeds[i1 * 3 + 2];
	
		double mx = 0.5 * (seeds[i2 * 3 + 0] + seeds[i1 * 3 + 0]);
		double my = 0.5 * (seeds[i2 * 3 + 1] + seeds[i1 * 3 + 1]);
		double mz = 0.5 * (seeds[i2 * 3 + 2] + seeds[i1 * 3 + 2]);

		double norm(sqrt(nx * nx + ny * ny + nz * nz));

		if (norm < 1E-6)
		{
			vcm_cout << "WARNING::SAVE VCG FILE:: Two Face seeds are too close" << std::endl;
		}
		nx /= norm; ny /= norm; nz /= norm;

		file << seed_new_index[i1] + 1 << " ";
		if (seed_region_id[i2] != 0)	file << seed_new_index[i2] + 1 << " ";
		else                            file << seed_new_index[i1] + 1 << " "; // duplicate seeds indicate a boundary connection
		file << mx << " " << my << " " << mz << " " << face_area_height[iface * 2] << " ";
		file << nx << " " << ny << " " << nz << std::endl;
	}

	size_t min_num_connections(SIZE_MAX), max_num_connections(0), num_connections(0);
	for (size_t iseed = 0; iseed < num_seeds; iseed++)
	{
		if (seed_region_id[iseed] == 0) continue;
		num_connections += cell_num_faces[iseed];
		if (cell_num_faces[iseed] < min_num_connections)  min_num_connections = cell_num_faces[iseed];
		if (cell_num_faces[iseed] > max_num_connections)  max_num_connections = cell_num_faces[iseed];
	}

	vcm_cout << "  * Total Number of connections = " << num_connections << std::endl;
	vcm_cout << "    Minimum cell connections = " << min_num_connections << std::endl;
	vcm_cout << "    Maximum cell connections = " << max_num_connections << std::endl;
	vcm_cout << "    Average cell connections = " << num_connections / num_meshed_cells << std::endl;
	delete[] seed_new_index; delete[] face_area_height;
	return 0;
	#pragma endregion
}

int MeshingVoronoiMesher::save_spheres_csv(std::string file_name, size_t num_dim, size_t num_spheres, double* spheres, double* spheres_sizing, size_t* spheres_region_id)
{
	#pragma region Save tree to CSV file:
	std::fstream file(file_name.c_str(), std::ios::out);
	// Spheres
	file << "x1coord";
	for (size_t idim = 1; idim < num_dim; idim++) file << ", x" << idim + 1 << "coord";
	file << ", radius" << std::endl;

	for (size_t i = 0; i < num_spheres; i++)
	{
		if (spheres_sizing != 0)
		{
			//if (fabs(spheres[i * num_dim] - 0.5) > spheres_sizing[i]) continue;
			file << std::setprecision(16) << spheres[i * num_dim];
			for (size_t idim = 1; idim < num_dim; idim++) file << ", " << spheres[i * num_dim + idim];
			file << ", " << spheres_sizing[i];
			if (spheres_region_id != 0) file << ", " << spheres_region_id[i];
		}
		else
		{
			file << std::setprecision(16) << spheres[i * (num_dim + 1)];
			for (size_t idim = 1; idim <= num_dim; idim++) file << ", " << spheres[i * (num_dim + 1) + idim];
			if (spheres_region_id != 0) file << ", " << spheres_region_id[i];
		}
		file << std::endl;
	}
	return 0;
	#pragma endregion
}

int MeshingVoronoiMesher::load_spheres_csv(std::string file_name, size_t num_dim, size_t &num_spheres, double* &spheres)
{
	#pragma region Load tree from CSV file:
	
	size_t cap(100); num_spheres = 0;
	spheres = new double[cap * (num_dim + 1)];

	//open file
	std::ifstream tmpfile(file_name.c_str());

	if (tmpfile.is_open() && tmpfile.good())
	{
		std::string line = "";
		size_t iline(0);
		while (getline(tmpfile, line))
		{
			std::istringstream iss(line);
			std::vector<std::string> tokens;
			get_tokens(line, ',', tokens);

			if (iline == 0)
			{
				iline++;
				continue;
			}
			
			if (tokens.size() <= num_dim) continue;

			for (size_t i = 0; i <= num_dim; i++) spheres[num_spheres * (num_dim + 1) + i] = _memo.string_to_double(tokens[i]);
			num_spheres++;

			if (num_spheres == cap)
			{
				cap *= 2;
				double* tmp_spheres = new double[cap * (num_dim + 1)];
				for (size_t i = 0; i < num_spheres * (num_dim + 1); i++) tmp_spheres[i] = spheres[i];
				delete[] spheres;
				spheres = tmp_spheres;
			}
			iline++;
		}
	}
	return 0;
	#pragma endregion
}


int MeshingVoronoiMesher::build_balanced_kd_tree(size_t num_points, size_t num_dim, double* points, size_t& tree_origin, size_t* tree_right, size_t* tree_left)
{
	#pragma region Build Balanced kd-tree:
	
	size_t* tree_nodes_sorted = new size_t[num_points];
	for (size_t i = 0; i < num_points; i++) tree_nodes_sorted[i] = i;

	for (size_t iseed = 0; iseed < num_points; iseed++)
	{
		tree_left[iseed] = iseed; tree_right[iseed] = iseed;
	}
	tree_origin = SIZE_MAX;

	size_t tree_height(0);

	size_t target_pos = num_points / 2; 
	kd_tree_balance_quicksort(num_points, num_dim, points, tree_origin, tree_right, tree_left, tree_height, 
		                      target_pos, 0, num_points - 1, 0, tree_nodes_sorted);

	delete[] tree_nodes_sorted;

	//vcm_cout << "      * Executed in " << cpu_time << " seconds." << std::endl;
	//vcm_cout << "      * Number of tree levels = " << tree_max_height << std::endl;
	//vcm_cout << "      * Number of perfectly balanced tree levels = " << ceil(log2(1 + num_points)) << std::endl;

	return 0;
	#pragma endregion
}

int MeshingVoronoiMesher::re_enumerate_points_for_better_memory_access(size_t num_points, size_t num_dim, double* points,
	                                                            size_t& tree_origin, size_t* tree_right, size_t* tree_left,
	                                                            double* points_sorted, size_t* point_old_index, size_t* point_new_index)
{
	#pragma region Re enumerate Points for better memory access:
	size_t num_traversed(0);
	kd_tree_get_nodes_order(tree_right,tree_left, 0, tree_origin, num_traversed, point_old_index);

	for (size_t i = 0; i < num_points; i++)
	{
		size_t current_point = point_old_index[i];
		point_new_index[current_point] = i;
	}

	// new tree containers with the current order

	size_t* tree_right_sorted = new size_t[num_points];
	size_t* tree_left_sorted = new size_t[num_points];
	
	for (size_t i = 0; i < num_points; i++)
	{
		size_t current_point = point_old_index[i];
		for (size_t idim = 0; idim < 3; idim++) points_sorted[i * num_dim + idim] = points[current_point * num_dim + idim];
		tree_right_sorted[i] = point_new_index[tree_right[current_point]];
		tree_left_sorted[i] = point_new_index[tree_left[current_point]];
	}
	tree_origin = point_new_index[tree_origin];

	for (size_t i = 0; i < num_points; i++)
	{
		tree_right[i] = tree_right_sorted[i];
		tree_left[i] = tree_left_sorted[i];
	}
	delete[] tree_right_sorted; delete[] tree_left_sorted;
	return 0;
	#pragma endregion
}


int MeshingVoronoiMesher::get_closest_tree_point(size_t num_points, size_t num_dim, double* points, size_t tree_origin, size_t* tree_right, size_t* tree_left, 
	                                             double* x, size_t& closest_tree_point, double& closest_distance)
{
	#pragma region Closest Neighbor Search using kd tree:
	closest_tree_point = num_points;
	size_t num_nodes_visited = 0;
	if (num_points == 0) return 1;
	kd_tree_get_closest_seed(num_points, num_dim, points, tree_origin, tree_right, tree_left, 
		                     x, 0, tree_origin, closest_tree_point, closest_distance, num_nodes_visited);
	return 0;
	#pragma endregion
}

int MeshingVoronoiMesher::get_closest_tree_point(size_t num_points, size_t num_dim, double* points, size_t tree_origin, size_t* tree_right, size_t* tree_left,
	                                             double* x, double* e_dir, size_t& closest_tree_point, double& closest_distance)
{
	#pragma region Closest Neighbor Search using kd tree:
	closest_tree_point = num_points;
	size_t num_nodes_visited = 0;
	if (num_points == 0) return 1;
	kd_tree_get_closest_seed(num_points, num_dim, points, tree_origin, tree_right, tree_left,
		x, e_dir, 0, tree_origin, closest_tree_point, closest_distance, num_nodes_visited);
	return 0;
	#pragma endregion
}

int MeshingVoronoiMesher::get_closest_tree_point(size_t num_points, size_t num_dim, double* points, size_t tree_origin, size_t* tree_right, size_t* tree_left,
	                                             size_t tree_point_index, double* e_dir, size_t& closest_tree_point, double& closest_distance)
{
	#pragma region Closest Neighbor Search using kd tree:
	closest_tree_point = num_points;
	size_t num_nodes_visited = 0;
	if (num_points == 0) return 1;
	kd_tree_get_closest_seed(num_points, num_dim, points, tree_origin, tree_right, tree_left,
		                     tree_point_index, e_dir, 0, tree_origin, closest_tree_point, closest_distance, num_nodes_visited);
	return 0;
	#pragma endregion
}


int MeshingVoronoiMesher::get_tree_points_in_sphere(size_t num_points, size_t num_dim, double* points, size_t tree_origin, size_t* tree_right, size_t* tree_left,
	                                         double* x, double r, size_t& num_points_in_sphere, size_t*& points_in_sphere)
{
	#pragma region tree sphere neighbor search:
	num_points_in_sphere = 0;
	size_t capacity = 10;
	points_in_sphere = new size_t[capacity];
	kd_tree_get_seeds_in_sphere(num_points, num_dim, points, tree_origin, tree_right, tree_left, 
		                        x, r, 0, tree_origin, num_points_in_sphere, points_in_sphere, capacity);
	if (num_points_in_sphere == 0)
	{
		delete[] points_in_sphere;
		points_in_sphere = 0;
	}
	return 0;
	#pragma endregion
}


int MeshingVoronoiMesher::kd_tree_balance_quicksort(size_t num_points, size_t num_dim, double* points,
	                                         size_t &tree_origin, size_t* tree_right, size_t* tree_left, size_t& tree_height,
	                                         size_t target_pos, size_t left, size_t right, size_t active_dim, size_t* tree_nodes_sorted)
{
	#pragma region kd tree balance:
	kd_tree_quicksort_adjust_target_position(num_points, num_dim, points, target_pos, left, right, active_dim, tree_nodes_sorted);

	// target position is correct .. add to tree
	if (tree_origin == SIZE_MAX)    tree_origin = tree_nodes_sorted[target_pos];
	else                            kd_tree_add_point(num_points, num_dim, points, tree_origin, tree_right, tree_left, tree_height, tree_nodes_sorted[target_pos]);

	/* recursion */
	active_dim++;
	if (active_dim == num_dim) active_dim = 0;

	if (target_pos + 1 < right)
	{
		kd_tree_balance_quicksort(num_points, num_dim, points, tree_origin, tree_right, tree_left, tree_height,
			(target_pos + 1 + right) / 2, target_pos + 1, right, active_dim, tree_nodes_sorted);
	}
	else if (right > target_pos)
	{
		kd_tree_add_point(num_points, num_dim, points, tree_origin, tree_right, tree_left, tree_height, tree_nodes_sorted[right]);
	}

	if (target_pos > left + 1)
	{
		kd_tree_balance_quicksort(num_points, num_dim, points, tree_origin, tree_right, tree_left, tree_height,
			(left + target_pos - 1) / 2, left, target_pos - 1, active_dim, tree_nodes_sorted);
	}
	else if (left < target_pos)
	{
		kd_tree_add_point(num_points, num_dim, points, tree_origin, tree_right, tree_left, tree_height, tree_nodes_sorted[left]);
	}
	return 0;
	#pragma endregion
}

int MeshingVoronoiMesher::kd_tree_quicksort_adjust_target_position(size_t num_points, size_t num_dim, double* points, 
	                                                        size_t target_pos, size_t left, size_t right, size_t active_dim, size_t* tree_nodes_sorted)
{
	#pragma region kd tree Quick sort pivot:
	size_t i = left, j = right;

	size_t pivot_seed = tree_nodes_sorted[(left + right) / 2];
	double pivot = points[pivot_seed * num_dim + active_dim];

	/* partition */
	while (i <= j)
	{
		while (points[tree_nodes_sorted[i] * num_dim + active_dim] < pivot)
			i++;
		while (points[tree_nodes_sorted[j] * num_dim + active_dim] > pivot)
			j--;

		if (i <= j)
		{
			size_t tmp_index = tree_nodes_sorted[i];
			tree_nodes_sorted[i] = tree_nodes_sorted[j];
			tree_nodes_sorted[j] = tmp_index;

			i++;
			if (j > 0) j--;
		}
	};

	/* recursion */
	if (i < right && i <= target_pos && right >= target_pos)
		kd_tree_quicksort_adjust_target_position(num_points, num_dim, points, target_pos, i, right, active_dim, tree_nodes_sorted);
	if (j > 0 && left < j && left <= target_pos && j >= target_pos)
		kd_tree_quicksort_adjust_target_position(num_points, num_dim, points, target_pos, left, j, active_dim, tree_nodes_sorted);
	return 0;
	#pragma endregion
}

int MeshingVoronoiMesher::kd_tree_add_point(size_t num_points, size_t num_dim, double* points,
	                                 size_t tree_origin, size_t* tree_right, size_t* tree_left, size_t &tree_height,
	                                 size_t seed_index)
{
	#pragma region kd tree add point:
	// insert sphere into tree
	size_t parent_index(tree_origin); size_t d_index(0);
	size_t branch_height(1);
	while (true)
	{
		if (points[seed_index * num_dim + d_index] > points[parent_index * num_dim + d_index])
		{
			if (tree_right[parent_index] == parent_index)
			{
				tree_right[parent_index] = seed_index;
				branch_height++;
				break;
			}
			else
			{
				parent_index = tree_right[parent_index];
				branch_height++;
			}
		}
		else
		{
			if (tree_left[parent_index] == parent_index)
			{
				tree_left[parent_index] = seed_index;
				branch_height++;
				break;
			}
			else
			{
				parent_index = tree_left[parent_index];
				branch_height++;
			}
		}
		d_index++;
		if (d_index == num_dim) d_index = 0;
	}
	if (branch_height > tree_height) tree_height = branch_height;
	return 0;
	#pragma endregion
}

int MeshingVoronoiMesher::kd_tree_get_closest_seed(size_t num_points, size_t num_dim, double* points,
	                                        size_t tree_origin, size_t* tree_right, size_t* tree_left,
	                                        double* x, size_t d_index, size_t node_index,
	                                        size_t& closest_seed, double& closest_distance,
	                                        size_t& num_nodes_visited)
{
	#pragma region kd tree closest neighbor search:
	if (d_index == num_dim) d_index = 0;

	double dst(0.0);
	for (size_t idim = 0; idim < num_dim; idim++)
	{
		double dx = x[idim] - points[node_index * num_dim + idim];
		dst += dx * dx;
	}
	dst = sqrt(dst);
	num_nodes_visited++;
	if (dst < closest_distance)
	{
		// add to neighbors
		closest_seed = node_index;
		closest_distance = dst;
	}

	double neighbor_max = x[d_index] + closest_distance;
	if (tree_right[node_index] != node_index && neighbor_max > points[node_index * num_dim + d_index])
	{
		kd_tree_get_closest_seed(num_points, num_dim, points, tree_origin, tree_right, tree_left,
			                     x, d_index + 1, tree_right[node_index], closest_seed, closest_distance, num_nodes_visited);
	}

	double neighbor_min = x[d_index] - closest_distance;
	if (tree_left[node_index] != node_index && neighbor_min < points[node_index * num_dim + d_index])
	{
		kd_tree_get_closest_seed(num_points, num_dim, points, tree_origin, tree_right, tree_left, 
			                     x, d_index + 1, tree_left[node_index], closest_seed, closest_distance, num_nodes_visited);
	}
	return 0;
	#pragma endregion
}

int MeshingVoronoiMesher::kd_tree_get_closest_seed(size_t num_points, size_t num_dim, double* points,
	                                               size_t tree_origin, size_t* tree_right, size_t* tree_left,
	                                               double* x, double* e_dir, size_t d_index, size_t node_index,
	                                               size_t& closest_seed, double& closest_distance,
	                                               size_t& num_nodes_visited)
{
	#pragma region kd tree closest neighbor search:
	if (d_index == num_dim) d_index = 0;

	double dst(0.0), dot(0.0);
	for (size_t idim = 0; idim < num_dim; idim++)
	{
		double dx = points[node_index * num_dim + idim] - x[idim];
		dst += dx * dx;

		if (e_dir != 0) dot += dx * e_dir[idim];
	}
	dst = sqrt(dst);
	num_nodes_visited++;
	if (dst < closest_distance && dot >= 0.0)
	{
		// update closest seed
		closest_seed = node_index;
		closest_distance = dst;
	}

	double neighbor_max = x[d_index] + closest_distance;
	if (tree_right[node_index] != node_index && neighbor_max > points[node_index * num_dim + d_index])
	{
		kd_tree_get_closest_seed(num_points, num_dim, points, tree_origin, tree_right, tree_left,
			x, e_dir, d_index + 1, tree_right[node_index], closest_seed, closest_distance, num_nodes_visited);
	}

	double neighbor_min = x[d_index] - closest_distance;
	if (tree_left[node_index] != node_index && neighbor_min < points[node_index * num_dim + d_index])
	{
		kd_tree_get_closest_seed(num_points, num_dim, points, tree_origin, tree_right, tree_left,
			x, e_dir, d_index + 1, tree_left[node_index], closest_seed, closest_distance, num_nodes_visited);
	}
	return 0;
	#pragma endregion
}

int MeshingVoronoiMesher::kd_tree_get_closest_seed(size_t num_points, size_t num_dim, double* points,
	                                               size_t tree_origin, size_t* tree_right, size_t* tree_left,
	                                               size_t tree_point_index, double* e_dir,
	                                               size_t d_index, size_t node_index,
	                                               size_t& closest_seed, double& closest_distance,
	                                               size_t& num_nodes_visited)
{
	#pragma region kd tree closest neighbor search:
	if (d_index == num_dim) d_index = 0;

	double dst(0.0), dot(0.0);
	for (size_t idim = 0; idim < num_dim; idim++)
	{
		double dx = points[node_index * num_dim + idim] - points[tree_point_index * num_dim + idim];
		dst += dx * dx;

		if (e_dir != 0) dot += dx * e_dir[idim];
	}
	dst = sqrt(dst);
	num_nodes_visited++;
	if (dst < closest_distance && dot >= 0.0 && tree_point_index != node_index)
	{
		// update closest seed
		closest_seed = node_index;
		closest_distance = dst;
	}

	double neighbor_max = points[tree_point_index * num_dim + d_index] + closest_distance;
	if (tree_right[node_index] != node_index && neighbor_max > points[node_index * num_dim + d_index])
	{
		kd_tree_get_closest_seed(num_points, num_dim, points, tree_origin, tree_right, tree_left,
			tree_point_index, e_dir, d_index + 1, tree_right[node_index], closest_seed, closest_distance, num_nodes_visited);
	}

	double neighbor_min = points[tree_point_index * num_dim + d_index] - closest_distance;
	if (tree_left[node_index] != node_index && neighbor_min < points[node_index * num_dim + d_index])
	{
		kd_tree_get_closest_seed(num_points, num_dim, points, tree_origin, tree_right, tree_left,
			tree_point_index, e_dir, d_index + 1, tree_left[node_index], closest_seed, closest_distance, num_nodes_visited);
	}
	return 0;
	#pragma endregion
}

int MeshingVoronoiMesher::kd_tree_get_seeds_in_sphere(size_t num_points, size_t num_dim, double* points,
	                                           size_t tree_origin, size_t* tree_right, size_t* tree_left,
	                                           double* x, double r, size_t d_index, size_t node_index,
	                                           size_t& num_points_in_sphere, size_t*& points_in_sphere, size_t& capacity)
{
	#pragma region kd tree recursive sphere neighbor search:
	if (d_index == num_dim) d_index = 0;

	double dst_sq(0.0);
	for (size_t idim = 0; idim < num_dim; idim++)
	{
		double dx = points[node_index * num_dim + idim] - x[idim];
		dst_sq += dx * dx;
	}
		
	if (dst_sq < (1 + 2.0E-6) * r * r)
	{
		points_in_sphere[num_points_in_sphere] = node_index;
		num_points_in_sphere++;

		if (num_points_in_sphere == capacity)
		{
			capacity *= 2;
			size_t* new_points = new size_t[capacity];
			for (size_t ipoint = 0; ipoint < num_points_in_sphere; ipoint++) new_points[ipoint] = points_in_sphere[ipoint];
			delete[] points_in_sphere;
			points_in_sphere = new_points;
		}
	}

	bool check_right(false), check_left(false);
	double neighbor_min(x[d_index] - r), neighbor_max(x[d_index] + r);

	if (tree_right[node_index] != node_index && neighbor_max > points[node_index * num_dim + d_index])
	{
		kd_tree_get_seeds_in_sphere(num_points, num_dim, points, tree_origin, tree_right, tree_left, 
			                        x, r, d_index + 1, tree_right[node_index], num_points_in_sphere, points_in_sphere, capacity);
	}

	if (tree_left[node_index] != node_index && neighbor_min < points[node_index * num_dim + d_index])
	{
		kd_tree_get_seeds_in_sphere(num_points, num_dim, points, tree_origin, tree_right, tree_left, 
			                        x, r, d_index + 1, tree_left[node_index], num_points_in_sphere, points_in_sphere, capacity);
	}
	return 0;
	#pragma endregion
}

int MeshingVoronoiMesher::kd_tree_get_nodes_order(size_t* tree_right, size_t* tree_left, 
	                                       size_t d_index, size_t node_index,
	                                       size_t& num_traversed, size_t* ordered_indices)
{
	ordered_indices[num_traversed] = node_index; num_traversed++;
	if (tree_right[node_index] != node_index) kd_tree_get_nodes_order(tree_right, tree_left, d_index + 1, tree_right[node_index], num_traversed, ordered_indices);
	if (tree_left[node_index] != node_index) kd_tree_get_nodes_order(tree_right, tree_left, d_index + 1, tree_left[node_index], num_traversed, ordered_indices);
	return 0;
}


int MeshingVoronoiMesher::impose_lipschitz_continuity(size_t num_spheres, size_t num_dim, double* spheres, double* spheres_sizing, 
	                                           int num_threads, double Lip,
	                                           size_t tree_origin, size_t* tree_right, size_t* tree_left)
{
	#pragma region Impose Lipschitz Continuity:
	
	double* x = new double[num_dim];
	size_t* num_neighbors = new size_t[num_spheres];
	size_t** neighbors = new size_t*[num_spheres];
	
#if defined USE_OPEN_MP
	omp_set_num_threads(num_threads);
#pragma omp parallel for
#endif
	for (int i = 0; i < num_spheres; i++)
	{ 
		for (size_t idim = 0; idim < num_dim; idim++) x[idim] = spheres[i * num_dim + idim];
		get_tree_points_in_sphere(num_spheres, num_dim, spheres, tree_origin, tree_right, tree_left, x, 2.0 * spheres_sizing[i], num_neighbors[i], neighbors[i]);
	}

	size_t num_iter(0);
	while (true)
	{
		bool done = true;
		for (size_t isphere = 0; isphere < num_spheres; isphere++)
		{
			for (size_t j = 0; j < num_neighbors[isphere]; j++)
			{
				size_t jsphere = neighbors[isphere][j];
				if (spheres_sizing[jsphere] > spheres_sizing[isphere]) continue;

				double h(0.0);
				for (size_t idim = 0; idim < num_dim; idim++)
				{
					double dx = spheres[jsphere * num_dim + idim] - spheres[isphere * num_dim + idim];
					h += dx * dx;
				}
				h = sqrt(h);
				double rLip = spheres_sizing[jsphere] + h * Lip;
				if (spheres_sizing[isphere] > rLip + 1E-10)
				{
					spheres_sizing[isphere] = rLip;
					done = false;
				}
			}
		}
		if (done) break;
	}

	for (size_t i = 0; i < num_spheres; i++)
	{
		if (num_neighbors[i] > 0) delete[] neighbors[i];
	}
	delete[] neighbors; delete[] num_neighbors; delete[] x;
	return 0;
	#pragma endregion

}

int MeshingVoronoiMesher::get_tokens(std::string line, char separator, std::vector<std::string>& tokens)
{
	#pragma region convert line to tokens:
	tokens.clear();
	size_t start(0);
	while (true)
	{
		size_t end = line.find_first_of(separator, start);
		if (end == std::string::npos) end = line.size();

		tokens.push_back(line.substr(start, end - start));
		start = end + 1;

		if (end == line.size())
		{
			// end of string, exit loop
			break;
		}
	}
	return 0;
	#pragma endregion
}

int MeshingVoronoiMesher::save_exodus_file(std::string file_name, size_t num_seeds, double* seeds, size_t* seed_region_id,
                                    size_t num_vertices, double* vertices, size_t num_faces, size_t** faces,
                                    size_t* cell_num_faces, size_t** cell_faces)
{
	#pragma region Save Exodus Files:
#ifdef USE_EXODUS
    int num_nodes = int(num_vertices);
    float* x = new float[num_nodes];
    float* y = new float[num_nodes];
    float* z = new float[num_nodes];
    for (size_t ivtx = 0; ivtx < num_vertices; ivtx++)
    {
        x[ivtx] = float(vertices[ivtx * 3]);
        y[ivtx] = float(vertices[ivtx * 3 + 1]);
        z[ivtx] = float(vertices[ivtx * 3 + 2]);
    }
    
    size_t faces_num_entries(0);
    for (size_t iface = 0; iface< num_faces; iface++) faces_num_entries += faces[iface][0];
    
    int* num_face_corners = new int[num_faces];
    int* faces_entries = new int[faces_num_entries];
    faces_num_entries = 0;
    for (size_t iface = 0; iface < num_faces; iface++)
    {
        num_face_corners[iface] = int(faces[iface][0]);
        for (size_t icorner = 1; icorner <= faces[iface][0]; icorner++)
        {
            faces_entries[faces_num_entries++] = int(faces[iface][icorner]) + 1;
        }
    }
    
    size_t num_valid_cells(0);
    size_t cells_num_entries(0);
    for (size_t icell = 0; icell < num_seeds; icell++)
    {
        if (seed_region_id[icell] == 0) continue; /* A ghost node */
        num_valid_cells++;
        cells_num_entries+= cell_num_faces[icell];
    }
    int* num_cell_faces = new int[num_valid_cells];
    int* cells_entries = new int[cells_num_entries];
    num_valid_cells = 0; cells_num_entries = 0;
    for (size_t icell = 0; icell < num_seeds; icell++)
    {
        if (seed_region_id[icell] == 0) continue; /* A ghost node */
        
        for (size_t iface = 0; iface < cell_num_faces[icell]; iface++)
        {
            cells_entries[cells_num_entries++] = int(cell_faces[icell][iface]) + 1;
        }
        num_cell_faces[num_valid_cells] = cell_num_faces[icell];
        num_valid_cells++;
    }
    
    
    int exoid, num_dim, num_elem, num_elem_blk;
    int num_elem_in_block[10], num_total_nodes_per_blk[10];
    int num_face_in_block[10], num_total_faces_per_blk[10];
    int num_node_sets, error;
    int i, j;
    int bids;
    int num_qa_rec, num_info;
    int CPU_word_size, IO_word_size;
    
    num_dim       = 3;
    num_elem      = num_valid_cells;
    num_elem_blk  = 1;
    num_node_sets = 0;
    
    char *coord_names[3], *qa_record[2][4], *info[3];
    char *block_names[10];
    char *title = "This is a test";
    ex_opts(EX_VERBOSE | EX_ABORT);

    /* Specify compute and i/o word size */

    CPU_word_size = 0; /* sizeof(float) */
    IO_word_size  = 4; /* (4 bytes) */
    
    /* create EXODUS II file */

    exoid = ex_create(file_name.c_str(), /* filename path */
                      EX_CLOBBER,        /* create mode */
                      &CPU_word_size,    /* CPU float word size in bytes */
                      &IO_word_size);    /* I/O float word size in bytes */
    printf("after ex_create for test.exo, exoid = %d\n", exoid);
    printf(" cpu word size: %d io word size: %d\n", CPU_word_size, IO_word_size);

    /* initialize file with parameters */
    {
      ex_init_params par;

      ex_copy_string(par.title, title, MAX_LINE_LENGTH + 1);
      par.num_dim       = num_dim;
      par.num_nodes     = num_nodes;
      par.num_edge      = 0;
      par.num_edge_blk  = 0;
      par.num_face      = int(num_faces);
      par.num_face_blk  = 1;
      par.num_elem      = num_elem;
      par.num_elem_blk  = num_elem_blk;
      par.num_node_sets = num_node_sets;
      par.num_edge_sets = 0;
      par.num_face_sets = 0;
      par.num_side_sets = 0;
      par.num_elem_sets = 0;
      par.num_node_maps = 0;
      par.num_edge_maps = 0;
      par.num_face_maps = 0;
      par.num_elem_maps = 0;

      error = ex_put_init_ext(exoid, &par);

      printf("after ex_put_init_ext, error = %d\n", error);

      if (error) {
        ex_close(exoid);
        exit(-1);
      }
    }

    /* write nodal coordinates values and names to database */
    
    error = ex_put_coord(exoid, x, y, z);
    printf("after ex_put_coord, error = %d\n", error);

    if (error) {
      ex_close(exoid);
      exit(-1);
    }

    coord_names[0] = "x";
    coord_names[1] = "y";
    coord_names[2] = "z";

    error = ex_put_coord_names(exoid, coord_names);
    printf("after ex_put_coord_names, error = %d\n", error);

    if (error) {
      ex_close(exoid);
      exit(-1);
    }

    /* Write the face block parameters */
    block_names[0]             = "face_block_1";
    num_face_in_block[0]       = int(num_faces);    // all faces in one block
    num_total_nodes_per_blk[0] = faces_num_entries; // number of all faces corners (size of connect array)
    bids                       = 10;                // I am not sure what is this

    error = ex_put_block(exoid, EX_FACE_BLOCK, bids, "nsided", num_face_in_block[0],
                         num_total_nodes_per_blk[0], 0, 0, 0);
    printf("after ex_put_block, error = %d\n", error);

    if (error) {
      ex_close(exoid);
      exit(-1);
    }

    /* Write face block names */
    error = ex_put_names(exoid, EX_FACE_BLOCK, block_names);
    printf("after ex_put_names, error = %d\n", error);

    if (error) {
      ex_close(exoid);
      exit(-1);
    }

    /* write face connectivity */
    error = ex_put_conn(exoid, EX_FACE_BLOCK, bids, faces_entries, NULL, NULL);
    printf("after ex_put_conn, error = %d\n", error);

    if (error) {
      ex_close(exoid);
      exit(-1);
    }

    error = ex_put_entity_count_per_polyhedra(exoid, EX_FACE_BLOCK, bids, num_face_corners);
    printf("after ex_put_entity_count_per_polyhedra, error = %d\n", error);

    if (error) {
      ex_close(exoid);
      exit(-1);
    }


    /* write element block parameters */
    block_names[0] = "nfaced_1";

    num_elem_in_block[0]       = num_valid_cells;
    num_total_faces_per_blk[0] = cells_num_entries;

    bids = 10;

    error = ex_put_block(exoid, EX_ELEM_BLOCK, bids, "nfaced", num_elem_in_block[0], 0, 0,
                             num_total_faces_per_blk[0], 0);
    printf("after ex_put_block, error = %d\n", error);

    if (error) {
      ex_close(exoid);
      exit(-1);
    }

    /* Write element block names */
    error = ex_put_names(exoid, EX_ELEM_BLOCK, block_names);
    printf("after ex_put_names, error = %d\n", error);

    if (error) {
      ex_close(exoid);
      exit(-1);
    }

    /* write element-face connectivity */
    error = ex_put_conn(exoid, EX_ELEM_BLOCK, bids, NULL, NULL, cells_entries);
    printf("after ex_put_conn, error = %d\n", error);

    if (error) {
      ex_close(exoid);
      exit(-1);
    }

    error = ex_put_entity_count_per_polyhedra(exoid, EX_ELEM_BLOCK, bids, num_cell_faces);
    printf("after ex_put_entity_count_per_polyhedra, error = %d\n", error);

    if (error) {
      ex_close(exoid);
      exit(-1);
    }

    /* write QA records; test empty and just blank-filled records */
    num_qa_rec = 2;

    qa_record[0][0] = "TESTWT-NFACED";
    qa_record[0][1] = "testwt-nfaced";
    qa_record[0][2] = "2010/02/15";
    qa_record[0][3] = "06:35:15";
    qa_record[1][0] = "";
    qa_record[1][1] = "                            ";
    qa_record[1][2] = "";
    qa_record[1][3] = "                        ";

    error = ex_put_qa(exoid, num_qa_rec, qa_record);
    printf("after ex_put_qa, error = %d\n", error);

    if (error) {
      ex_close(exoid);
      exit(-1);
    }

    /* write information records; test empty and just blank-filled records */
    num_info = 3;

    info[0] = "This is the first information record.";
    info[1] = "";
    info[2] = "                                     ";

    error = ex_put_info(exoid, num_info, info);
    printf("after ex_put_info, error = %d\n", error);

    if (error) {
      ex_close(exoid);
      exit(-1);
    }

    /* close the EXODUS files
     */
    error = ex_close(exoid);
    printf("after ex_close, error = %d\n", error);
    if (error) {
      ex_close(exoid);
      exit(-1);
    }
    return 0;
#endif
    return 1;
	#pragma endregion
}
