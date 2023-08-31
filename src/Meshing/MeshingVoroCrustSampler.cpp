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
// MeshingVoroCrustSampler.cpp                                    Last modified (07/12/2022) //
///////////////////////////////////////////////////////////////////////////////////////////////

#include "MeshingVoroCrustSampler.h"
#include "MeshingVoronoiMesher.h"

#define alpha_sz 0.4

MeshingVoroCrustSampler::MeshingVoroCrustSampler(size_t num_threads)
{
	vcm_cout << "VoroCrust::Initiating VoroCrust Sampler .. " << std::endl;

	_input_sizing_function = 0;

	_r_min_edges = 0.0; _r_min_surface = 0.0;
	
	_ymin = DBL_MAX;

	// _r_min_edges = 0.1; _r_min_surface = 0.02;

	_face_normal = 0;

	_active_pool_size = 0; _active_pool_parent_face = 0; _active_pool_children = 0; _active_pool_cdf = 0;

	_pool_sampling_time = 0.0;  _pool_refinement_time = 0.0; _sizing_time = 0.0;

	size_t num_available_threads = 0;
#if defined USE_OPEN_MP
	#pragma omp parallel
	{
		int thread_id = omp_get_thread_num();
		#pragma omp critical(num_threads)
		if (thread_id > num_available_threads) num_available_threads = thread_id;
	}
#else
	vcm_cout << "  * OPNE_MP is disabled, VoroCrust is using a single thread!!" << std::endl;
#endif
	num_available_threads++;

	if (num_threads < num_available_threads) _num_threads = num_threads;
	else                                     _num_threads = num_available_threads;

	_rsamplers.resize(_num_threads);

#if defined USE_OPEN_MP
	vcm_cout << "  * VoroCrust is using " << _num_threads << " / "<< num_available_threads << " available threads!" << std::endl;
	omp_set_num_threads(int(_num_threads));
	#pragma omp parallel
	{
		int thread_id = omp_get_thread_num();
#else
	for (size_t thread_id = 0; thread_id < _num_threads; thread_id++)
	{
#endif
		_rsamplers[thread_id] = new MeshingRandomSampler(thread_id);
	}
}

MeshingVoroCrustSampler::~MeshingVoroCrustSampler()
{
	clear_memory();
}

int MeshingVoroCrustSampler::dont_detect_features_generate_point_clouds(size_t num_points, double** points, size_t num_faces, size_t** faces,
	                                                             MeshingSmartTree* surface_point_cloud, MeshingSmartTree* edge_point_cloud,
	                                                             size_t &num_sharp_edges, size_t** &sharp_edges, size_t &num_sharp_corners, size_t*& sharp_corners)
{
	#pragma region Generate point clouds for smooth models:
	vcm_cout << "VoroCrust::Point Clouds Generation for Smooth models:" << std::endl;

	clock_t start_time, end_time; double cpu_time;
	start_time = clock();

	num_sharp_corners = 0; num_sharp_edges = 0;

	surface_point_cloud->disable_auto_balance();
	double* dart = new double[4]; double* proj = new double[3];
	double** corners = new double* [3];
	_face_normal = new double* [num_faces];

	bool repeat_index(false);
	for (size_t iface = 0; iface < num_faces; iface++)
	{
		#pragma region Generate Surface Point Cloud uniform on each face:
		if (repeat_index)
		{
			iface--;
			repeat_index = false;
		}

		size_t* attrib = new size_t[2];
		attrib[0] = 2; attrib[1] = iface;

		for (size_t i = 1; i <= 3; i++)
		{
			size_t corner_index = faces[iface][i];
			corners[i - 1] = points[corner_index];
		}

		_face_normal[iface] = new double[3];
		if (!_geom.get_3d_triangle_normal(corners, _face_normal[iface]))
		{
			vcm_cout << "Warning::A degenerate face has been detected and removed!" << std::endl;
			delete[] _face_normal[iface]; delete[] faces[iface];
			faces[iface] = faces[num_faces - 1]; num_faces--;
			repeat_index = true; continue;
		}

		for (size_t i = 1; i <= 3; i++)
		{
			size_t corner_index = faces[iface][i];
			surface_point_cloud->add_tree_point(3, points[corner_index], _face_normal[iface], attrib);
		}

		for (size_t i = 1; i <= 3; i++)
		{
			size_t ip = i + 1; if (ip > 3) ip = 1;
			size_t i_index = faces[iface][i];
			size_t ip_index = faces[iface][ip];

			for (size_t idim = 0; idim < 3; idim++) dart[idim] = 0.5 * (points[i_index][idim] + points[ip_index][idim]);
			surface_point_cloud->add_tree_point(3, dart, _face_normal[iface], attrib);
		}

		/*
		size_t num_mc_points(1000000 / num_faces);
		if (num_mc_points > 3)
		{
			num_mc_points -= 3;
			for (size_t ipoint = 0; ipoint < num_mc_points; ipoint++)
			{
				_rsampler.sample_uniformly_from_simplex(dart, 3, 3, corners);
				surface_point_cloud->add_tree_point(3, dart, face_normal[iface], attrib);
			}
		}
		*/
		
		delete[] attrib;
		#pragma endregion
	}

	if (true)
	{
		double* cdf = new double[num_faces]; double total_area(0.0);
		for (size_t iface = 0; iface < num_faces; iface++)
		{
			size_t i1(faces[iface][1]), i2(faces[iface][2]), i3(faces[iface][3]);
			corners[0] = points[i1]; corners[1] = points[i2]; corners[2] = points[i3];
			cdf[iface] = _geom.get_3d_triangle_area(corners);
			total_area += cdf[iface];
		}
		for (size_t iface = 1; iface < num_faces; iface++) cdf[iface] += cdf[iface - 1];
		for (size_t iface = 0; iface < num_faces; iface++) cdf[iface] /= cdf[num_faces - 1]; 

		_r_min_surface = sqrt(total_area / 1000000);

		size_t num_mc_points(1000000);
		double min_feature_size_surface = 10 * sqrt(total_area / num_mc_points);
		for (size_t ipoint = 0; ipoint < num_mc_points; ipoint++)
		{
			#pragma region Generate Surface Point Cloud uniform all faces:
			double u = _rsampler.generate_uniform_random_number();
			size_t ist(0), iend(num_faces - 1);
			while (iend > ist)
			{
				size_t imid = (ist + iend) / 2;
				if (cdf[imid] > u)
				{
					iend = imid;
				}
				else if (cdf[imid] < u)
				{
					ist = imid + 1;
				}
				else
				{
					ist = imid; iend = ist;
				}
			}

			size_t face_index = ist;

			size_t i1(faces[face_index][1]), i2(faces[face_index][2]), i3(faces[face_index][3]);
			corners[0] = points[i1]; corners[1] = points[i2]; corners[2] = points[i3];

			_rsampler.sample_uniformly_from_simplex(dart, 3, 3, corners);

			size_t* attrib = new size_t[2];
			attrib[0] = 2; attrib[1] = face_index;
			surface_point_cloud->add_tree_point(3, dart, _face_normal[face_index], attrib);
			delete[] attrib;
			#pragma endregion
		}
		delete[] cdf;
	}

	surface_point_cloud->kd_tree_build_balanced();
	surface_point_cloud->save_tree_csv("surface_point_cloud.csv", 3);
	vcm_cout << "  * Surface point cloud has " << surface_point_cloud->get_num_tree_points() << " points" << std::endl;

	_num_faces = num_faces;
	
	delete[] corners;
	delete[] dart; delete[] proj;

	end_time = clock();
	cpu_time = ((double)(end_time - start_time)) / CLOCKS_PER_SEC;
	vcm_cout << "  * executed in " << cpu_time << " seconds!" << std::endl;
	return 0;
	#pragma endregion
}

int MeshingVoroCrustSampler::detect_features_generate_point_clouds_clean(size_t num_points, double** points, size_t num_faces, size_t** faces,
	                                                              double smooth_angle_threshold,
	                                                              MeshingSmartTree* surface_point_cloud, MeshingSmartTree* edge_point_cloud,
	                                                              size_t &num_sharp_edges, size_t** &sharp_edges, size_t &num_sharp_corners, size_t*& sharp_corners)
{
	#pragma region Generate Surface seeds for a clean model:
	vcm_cout << "VoroCrust::Feature Detection and Point Clouds Generation (Clean Input):" << std::endl;

	clock_t start_time, end_time; double cpu_time;
	start_time = clock();

	size_t** point_faces;
	_mesh.build_point_faces(num_points, num_faces, faces, point_faces);

	_mesh.build_face_normals(points, num_faces, faces, _face_normal);
	_num_faces = num_faces;

	// extract sharp edges
	size_t cap_sharp_edges(100);
	num_sharp_edges = 0;
	sharp_edges = new size_t*[cap_sharp_edges];
	for (size_t iface = 0; iface < num_faces; iface++)
	{
		#pragma region Sharp Edges:
		for (size_t ic = 1; ic <= 3; ic++)
		{
			size_t icp = ic + 1; if (icp == 4) icp = 1;
			size_t ic_index = faces[iface][ic];
			size_t icp_index = faces[iface][icp];

			size_t num_edge_faces = _mesh.get_num_edge_faces(ic_index, icp_index, faces, point_faces);

			if (num_edge_faces != 2)
			{
				size_t* sharp_edge = new size_t[2];
				if (ic_index < icp_index) { sharp_edge[0] = ic_index; sharp_edge[1] = icp_index; }
				else { sharp_edge[0] = icp_index; sharp_edge[1] = ic_index; }
				_memo.add_entry(sharp_edge, num_sharp_edges, sharp_edges, cap_sharp_edges);
			}
			else if (ic_index < icp_index)
			{
				size_t* edge_faces;
				_mesh.get_edge_faces(ic_index, icp_index, faces, point_faces, edge_faces);
				size_t fi = edge_faces[1]; size_t fj = edge_faces[2];
				delete[] edge_faces;
				double cos_ang = _geom.dot_product(3, _face_normal[fi], _face_normal[fj]);
				if (cos_ang < smooth_angle_threshold)
				{
					size_t* sharp_edge = new size_t[2];
					if (ic_index < icp_index) { sharp_edge[0] = ic_index; sharp_edge[1] = icp_index; }
					else { sharp_edge[0] = icp_index; sharp_edge[1] = ic_index; }
					_memo.add_entry(sharp_edge, num_sharp_edges, sharp_edges, cap_sharp_edges);
				}
			}
		}
		#pragma endregion
	}

	size_t cap_sharp_corners(100);
	num_sharp_corners = 0;
	sharp_corners = new size_t[cap_sharp_corners];
	size_t* num_point_sharp_edges = new size_t[num_points];
	size_t** point_sharp_edges = new size_t*[num_points];
	for (size_t ipoint = 0; ipoint < num_points; ipoint++) num_point_sharp_edges[ipoint] = 0;

	for (size_t iedge = 0; iedge < num_sharp_edges; iedge++)
	{
		size_t i1 = sharp_edges[iedge][0];
		size_t i2 = sharp_edges[iedge][1];

		num_point_sharp_edges[i1]++;
		num_point_sharp_edges[i2]++;
	}

	for (size_t ipoint = 0; ipoint < num_points; ipoint++)
	{
		size_t num = num_point_sharp_edges[ipoint];
		point_sharp_edges[ipoint] = new size_t[num];
		num_point_sharp_edges[ipoint] = 0;
	}

	for (size_t iedge = 0; iedge < num_sharp_edges; iedge++)
	{
		size_t i1 = sharp_edges[iedge][0];
		size_t i2 = sharp_edges[iedge][1];

		size_t num_1 = num_point_sharp_edges[i1];
		point_sharp_edges[i1][num_1] = iedge;

		size_t num_2 = num_point_sharp_edges[i2];
		point_sharp_edges[i2][num_2] = iedge;

		num_point_sharp_edges[i1]++;
		num_point_sharp_edges[i2]++;
	}

	for (size_t ipoint = 0; ipoint < num_points; ipoint++)
	{
		#pragma region recount number of sharp edge without redundant edges:
		size_t num_sharp_edges(0);
		if (num_point_sharp_edges[ipoint] > 0)
		{
			num_sharp_edges = 1;
			for (size_t ied = 1; ied < num_point_sharp_edges[ipoint]; ied++)
			{
				size_t iedge_index = point_sharp_edges[ipoint][ied];
				size_t i1 = sharp_edges[iedge_index][0]; size_t i2 = sharp_edges[iedge_index][1];
				bool redundant(false);
				for (size_t jed = 0; jed < ied; jed++)
				{
					size_t jedge_index = point_sharp_edges[ipoint][jed];
					size_t j1 = sharp_edges[jedge_index][0]; size_t j2 = sharp_edges[jedge_index][1];
					if (i1 == j1 && i2 == j2) { redundant = true; break; }
				}
				if (redundant) continue;
				num_sharp_edges++;
			}
		}

		if (num_sharp_edges > 0 && num_sharp_edges != 2)
		{
			// A sharp corner is connected to one sharp edge or more than two sharp edges
			_memo.add_entry(ipoint, num_sharp_corners, sharp_corners, cap_sharp_corners);
		}
		else if (num_sharp_edges == 2)
		{
			size_t iedge = point_sharp_edges[ipoint][0];
			size_t jedge(0);
			for (size_t ied = 1; ied < num_point_sharp_edges[ipoint]; ied++)
			{
				size_t iedge_index = point_sharp_edges[ipoint][ied];
				size_t i1 = sharp_edges[iedge_index][0]; size_t i2 = sharp_edges[iedge_index][1];
				bool redundant(false);
				for (size_t jed = 0; jed < ied; jed++)
				{
					size_t jedge_index = point_sharp_edges[ipoint][jed];
					size_t j1 = sharp_edges[jedge_index][0]; size_t j2 = sharp_edges[jedge_index][1];
					if (i1 == j1 && i2 == j2) { redundant = true; break; }
				}
				if (redundant) continue;
				jedge = iedge_index;
				break;
			}
			size_t i1 = sharp_edges[iedge][0]; size_t i2 = sharp_edges[iedge][1];
			if (i1 != ipoint)
			{
				size_t tmp = i1; i1 = i2; i2 = tmp;
			}
			size_t j1 = sharp_edges[jedge][0]; size_t j2 = sharp_edges[jedge][1];
			if (j1 != ipoint)
			{
				size_t tmp = j1; j1 = j2; j2 = tmp;
			}
			double dot = -_geom.cos_angle(3, points[i2], points[ipoint], points[j2]);
			if (dot < smooth_angle_threshold)
			{
				// a sharp corner connected to sharp edges enclosing a sharp angle
				_memo.add_entry(ipoint, num_sharp_corners, sharp_corners, cap_sharp_corners);
			}
		}
		#pragma endregion
	}
	

	vcm_cout << "  * Number of Sharp Corners = " << num_sharp_corners << std::endl;
	vcm_cout << "  * Number of Sharp Edges = " << num_sharp_edges << std::endl;

	delete[] num_point_sharp_edges;

	if (num_sharp_edges > 0)
	{
		#pragma region Sampling edges:
		double* cdf = new double[num_sharp_edges];
		double* dart = new double[4];

		for (size_t iedge = 0; iedge < num_sharp_edges; iedge++)
		{
			size_t i1(sharp_edges[iedge][0]), i2(sharp_edges[iedge][1]);
			cdf[iedge] = _geom.distance(3, points[i1], points[i2]);
		}
		for (size_t iedge = 1; iedge < num_sharp_edges; iedge++) cdf[iedge] += cdf[iedge - 1];
		for (size_t iedge = 0; iedge < num_sharp_edges; iedge++) cdf[iedge] /= cdf[num_sharp_edges - 1];

		size_t num_mc_points(100000);
		for (size_t ipoint = 0; ipoint < num_mc_points; ipoint++)
		{
			#pragma region Uniform Sampling of Sharp Edges:
			double u = _rsampler.generate_uniform_random_number();
			size_t ist(0), iend(num_sharp_edges - 1);
			while (iend > ist)
			{
				size_t imid = (ist + iend) / 2;
				if (cdf[imid] > u)
				{
					iend = imid;
				}
				else if (cdf[imid] < u)
				{
					ist = imid + 1;
				}
				else
				{
					ist = imid; iend = ist;
				}
			}

			size_t edge_index = ist;
			size_t i1(sharp_edges[edge_index][0]), i2(sharp_edges[edge_index][1]);
			u = _rsampler.generate_uniform_random_number();
			for (size_t idim = 0; idim < 3; idim++) dart[idim] = points[i1][idim] + u * (points[i2][idim] - points[i1][idim]);

			double* tangent = new double[3];
			for (size_t idim = 0; idim < 3; idim++) tangent[idim] = points[i2][idim] - points[i1][idim];
			_geom.normalize_vector(3, tangent);

			size_t* attrib = new size_t[2];
			attrib[0] = 2; attrib[1] = edge_index;

			dart[3] = 1.0;

			edge_point_cloud->add_tree_point(4, dart, tangent, attrib);

			delete[] attrib; delete[] tangent;
			#pragma endregion
		}
		delete[] cdf; delete[] dart;
		#pragma endregion
	}

	size_t num_mc_face_points(100000);
	if (true)
	{
		#pragma region Sampling Surface:
		double* cdf = new double[num_faces];
		double** corners = new double*[3];
		double* dart = new double[4];

		for (size_t iface = 0; iface < num_faces; iface++)
		{
			size_t i1(faces[iface][1]), i2(faces[iface][2]), i3(faces[iface][3]);
			corners[0] = points[i1]; corners[1] = points[i2]; corners[2] = points[i3];
			cdf[iface] = _geom.get_3d_triangle_area(corners);
		}
		for (size_t iface = 1; iface < num_faces; iface++) cdf[iface] += cdf[iface - 1];
		for (size_t iface = 0; iface < num_faces; iface++) cdf[iface] /= cdf[num_faces - 1];
		
		for (size_t ipoint = 0; ipoint < num_mc_face_points; ipoint++)
		{
			#pragma region Uniform Sampling of Faces:
			double u = _rsampler.generate_uniform_random_number();
			size_t ist(0), iend(num_faces - 1);
			while (iend > ist)
			{
				size_t imid = (ist + iend) / 2;
				if (cdf[imid] > u)
				{
					iend = imid;
				}
				else if (cdf[imid] < u)
				{
					ist = imid + 1;
				}
				else
				{
					ist = imid; iend = ist;
				}
			}

			size_t face_index = ist;

			size_t i1(faces[face_index][1]), i2(faces[face_index][2]), i3(faces[face_index][3]);
			corners[0] = points[i1]; corners[1] = points[i2]; corners[2] = points[i3];

			_rsampler.sample_uniformly_from_simplex(dart, 3, 3, corners);

			dart[3] = 1.0;

			size_t* attrib = new size_t[2];
			attrib[0] = 2; attrib[1] = face_index;
			surface_point_cloud->add_tree_point(4, dart, _face_normal[face_index], attrib);
			delete[] attrib;
			#pragma endregion
		}
		
		delete[] cdf; delete[] corners; delete[] dart;
		#pragma endregion
	}
	
	for (size_t ipoint = 0; ipoint < num_points; ipoint++) delete[] point_faces[ipoint];
	delete[] point_faces;

	for (size_t ipoint = 0; ipoint < num_points; ipoint++) delete[] point_sharp_edges[ipoint];
	delete[] point_sharp_edges;

	end_time = clock();
	cpu_time = ((double)(end_time - start_time)) / CLOCKS_PER_SEC;
	vcm_cout << "  * executed in " << cpu_time << " seconds!" << std::endl;

	return 0;
	#pragma endregion
}

int MeshingVoroCrustSampler::detect_features_generate_point_clouds(size_t num_points, double** points, size_t num_faces, size_t** faces,
	                                                        double smooth_angle_threshold, double rmin,
	                                                        MeshingSmartTree* surface_point_cloud, MeshingSmartTree* edge_point_cloud,
	                                                        size_t &num_sharp_edges, size_t** &sharp_edges, size_t &num_sharp_corners, size_t* &sharp_corners)
{
	#pragma region Detect Sharp features and generate point clouds:
	vcm_cout << "VoroCrust::Feature Detection and Point Clouds Generation:" << std::endl;

	clock_t start_time, end_time; double cpu_time;
	start_time = clock();

	num_sharp_corners = 0; num_sharp_edges = 0;

	surface_point_cloud->disable_auto_balance();
	double* dart = new double[4]; double* proj = new double[3];
	double** corners = new double*[3];
	double** face_normal = new double*[num_faces];

	size_t* attrib = new size_t[2];

	bool repeat_index(false);
	for (size_t iface = 0; iface < num_faces; iface++)
	{
		#pragma region Generate Surface Point Cloud uniform on each face:
		if (repeat_index)
		{
			iface--;
			repeat_index = false;
		}

		attrib[0] = 2; attrib[1] = iface;

		for (size_t i = 1; i <= 3; i++)
		{
			size_t corner_index = faces[iface][i];
			corners[i - 1] = points[corner_index];
		}

		face_normal[iface] = new double[3];
		if (!_geom.get_3d_triangle_normal(corners, face_normal[iface]))
		{
			vcm_cout << "Warning::A degenerate face has been detected and removed!" << std::endl;
			delete[] face_normal[iface]; delete[] faces[iface];
			faces[iface] = faces[num_faces - 1]; num_faces--;
			repeat_index = true; continue;
		}

		for (size_t i = 1; i <= 3; i++)
		{
			size_t corner_index = faces[iface][i];
			surface_point_cloud->add_tree_point(3, points[corner_index], face_normal[iface], attrib);
		}

		for (size_t i = 1; i <= 3; i++)
		{
			size_t ip = i + 1; if (ip > 3) ip = 1;
			size_t i_index = faces[iface][i];
			size_t ip_index = faces[iface][ip];

			for (size_t idim = 0; idim < 3; idim++) dart[idim] = 0.5 * (points[i_index][idim] + points[ip_index][idim]);
			surface_point_cloud->add_tree_point(3, dart, face_normal[iface], attrib);
		}

		size_t num_mc_points(1000000 / num_faces);
		if (num_mc_points > 3)
		{
			num_mc_points -= 3;
			for (size_t ipoint = 0; ipoint < num_mc_points; ipoint++)
			{
				_rsampler.sample_uniformly_from_simplex(dart, 3, 3, corners);
				surface_point_cloud->add_tree_point(3, dart, face_normal[iface], attrib);
			}
		}
		#pragma endregion
	}

	double* cdf = new double[num_faces]; double total_area(0.0);
	for (size_t iface = 0; iface < num_faces; iface++)
	{
		size_t i1(faces[iface][1]), i2(faces[iface][2]), i3(faces[iface][3]);
		corners[0] = points[i1]; corners[1] = points[i2]; corners[2] = points[i3];
		cdf[iface] = _geom.get_3d_triangle_area(corners);
		total_area += cdf[iface];
	}
	for (size_t iface = 1; iface < num_faces; iface++) cdf[iface] += cdf[iface - 1];
	for (size_t iface = 0; iface < num_faces; iface++) cdf[iface] /= cdf[num_faces - 1];

	size_t num_mc_points(1000000);
	double min_feature_size_surface = 10 * sqrt(total_area / num_mc_points);
	for (size_t ipoint = 0; ipoint < num_mc_points; ipoint++)
	{
		#pragma region Generate Surface Point Cloud uniform all faces:
		double u = _rsampler.generate_uniform_random_number();
		size_t ist(0), iend(num_faces - 1);
		while (iend > ist)
		{
			size_t imid = (ist + iend) / 2;
			if (cdf[imid] > u)
			{
				iend = imid;
			}
			else if (cdf[imid] < u)
			{
				ist = imid + 1;
			}
			else
			{
				ist = imid; iend = ist;
			}
		}

		size_t face_index = ist;

		size_t i1(faces[face_index][1]), i2(faces[face_index][2]), i3(faces[face_index][3]);
		corners[0] = points[i1]; corners[1] = points[i2]; corners[2] = points[i3];

		_rsampler.sample_uniformly_from_simplex(dart, 3, 3, corners);


		attrib[0] = 2; attrib[1] = face_index;
		surface_point_cloud->add_tree_point(3, dart, face_normal[face_index], attrib);
		#pragma endregion
	}
	delete[] cdf;
	delete[] attrib;

	surface_point_cloud->kd_tree_build_balanced();
	//surface_point_cloud.write_points_csv("surface_point_cloud.csv", 3);
	vcm_cout << "  * Surface point cloud has " << surface_point_cloud->get_num_tree_points() << " points" << std::endl;

	size_t* cap_normals = new size_t[num_points];
	size_t* num_normals = new size_t[num_points];
	double*** point_normals = new double**[num_points];
	for (size_t ipoint = 0; ipoint < num_points; ipoint++)
	{
		size_t cap(10), num(0);
		cap_normals[ipoint] = cap; num_normals[ipoint] = num;
		point_normals[ipoint] = new double*[cap];
	}

	bool* point_visited = new bool[num_points];
	for (size_t ipoint = 0; ipoint < num_points; ipoint++) point_visited[ipoint] = false;

	for (size_t iface = 0; iface < num_faces; iface++)
	{
		#pragma region Collecting point normals:

		for (size_t i = 1; i <= 3; i++)
		{
			size_t point_index = faces[iface][i];

			if (point_visited[point_index]) continue;

			point_visited[point_index] = true;

			bool redundant(false);
			for (size_t inormal = 0; inormal < num_normals[point_index]; inormal++)
			{
				double dot = _geom.dot_product(3, point_normals[point_index][inormal], face_normal[iface]);
				if (dot > smooth_angle_threshold)
				{
					redundant = true;
					break;
				}
			}
			if (!redundant)
			{
				for (size_t idim = 0; idim < 3; idim++) dart[idim] = face_normal[iface][idim];
				_memo.add_entry(dart, num_normals[point_index], point_normals[point_index], cap_normals[point_index]);
				dart = new double[3];
			}

			size_t num_mc = 10;
			for (size_t imc = 0; imc < num_mc; imc++)
			{
				#pragma region Capturing Sharp Neighbors Faces with narrow angles:
				for (size_t j = 1; j <= 3; j++)
				{
					size_t corner_index = faces[iface][j];
					corners[j - 1] = points[corner_index];
				}
				_rsampler.sample_uniformly_from_simplex(dart, 3, 3, corners);
				size_t iclosest(0); double hclosest(DBL_MAX);
				surface_point_cloud->get_closest_non_smooth_tree_point(dart, rmin, smooth_angle_threshold, iclosest, hclosest);
				size_t* attrib = surface_point_cloud->get_tree_point_attrib(iclosest);
				size_t face_index = attrib[1];
				double dot = _geom.dot_product(3, face_normal[iface], face_normal[face_index]);
				if (dot <= smooth_angle_threshold) continue; // faces do not form a narrow angle

				for (size_t j = 1; j <= 3; j++)
				{
					size_t corner_index = faces[face_index][j];
					corners[j - 1] = points[corner_index];
				}
				double proj_dst(0.0);
				_geom.project_to_3d_triangle(points[point_index], corners, proj, proj_dst);
				if (proj_dst < rmin)
				{
					for (size_t idim = 0; idim < 3; idim++) dart[idim] = -face_normal[face_index][idim];
					bool redundant(false);
					for (size_t inormal = 0; inormal < num_normals[point_index]; inormal++)
					{
						double dot = _geom.dot_product(3, point_normals[point_index][inormal], dart);
						if (dot > smooth_angle_threshold)
						{
							redundant = true;
							break;
						}
					}
					if (!redundant)
					{
						_memo.add_entry(dart, num_normals[point_index], point_normals[point_index], cap_normals[point_index]);
						dart = new double[3];
					}
				}
				#pragma endregion
			}

			while (true)
			{
				#pragma region Collecting normals around point_index:
				bool done = true;
				size_t iclosest(0); double hclosest(DBL_MAX);
				surface_point_cloud->get_closest_non_smooth_tree_point(points[point_index], num_normals[point_index], point_normals[point_index], smooth_angle_threshold, iclosest, hclosest);

				if (hclosest < min_feature_size_surface)
				{
					size_t* attrib = surface_point_cloud->get_tree_point_attrib(iclosest);
					size_t face_index = attrib[1];
					size_t j1(faces[face_index][1]), j2(faces[face_index][2]), j3(faces[face_index][3]);
					corners[0] = points[j1]; corners[1] = points[j2]; corners[2] = points[j3];

					double proj_dst(0.0);
					_geom.project_to_3d_triangle(points[point_index], corners, proj, proj_dst);

					if (proj_dst < rmin)
					{
						// Another normal
						for (size_t idim = 0; idim < 3; idim++) dart[idim] = face_normal[face_index][idim];

						bool redundant(false);
						for (size_t inormal = 0; inormal < num_normals[point_index]; inormal++)
						{
							double dot = _geom.dot_product(3, point_normals[point_index][inormal], dart);
							if (dot > smooth_angle_threshold)
							{
								redundant = true;
								break;
							}
						}
						if (!redundant)
						{
							_memo.add_entry(dart, num_normals[point_index], point_normals[point_index], cap_normals[point_index]);
							dart = new double[3];
							done = false;
						}
					}
				}
				if (done) break;
				#pragma endregion
			}
		}
		#pragma endregion
	}

	// detecting sharp edges
	size_t cap_sharp_edges(100); num_sharp_edges = 0;
	sharp_edges = new size_t*[cap_sharp_edges];
	double* mid_edge = new double[3];
	for (size_t iface = 0; iface < num_faces; iface++)
	{
		#pragma region Detect Sharp Edges:

		if (true)
		{
			size_t i1(faces[iface][1]), i2(faces[iface][2]), i3(faces[iface][3]);
			for (size_t idim = 0; idim < 3; idim++) dart[idim] = (points[i1][idim] + points[i2][idim] + points[i3][idim]) / 3;
			double d1 = _geom.distance(3, dart, points[i1]);
			double d2 = _geom.distance(3, dart, points[i2]);
			double d3 = _geom.distance(3, dart, points[i3]);
			double dmax = fmax(d1, d2);
			dmax = fmax(dmax, d3);

			size_t iclosest(0); double hclosest(DBL_MAX);
			surface_point_cloud->get_closest_non_smooth_tree_point(dart, rmin, smooth_angle_threshold, iclosest, hclosest);
			if (hclosest > dmax) continue;
		}


		for (size_t i = 1; i <= 3; i++)
		{
			size_t ip = i + 1; if (ip == 4) ip = 1;
			size_t i_index = faces[iface][i];
			size_t ip_index = faces[iface][ip];

			for (size_t idim = 0; idim < 3; idim++) mid_edge[idim] = 0.5 * (points[i_index][idim] + points[ip_index][idim]);

			size_t num_mc(100);
			for (size_t imc = 0; imc < num_mc; imc++)
			{
				size_t i1(faces[iface][1]), i2(faces[iface][2]), i3(faces[iface][3]);
				corners[0] = points[i1]; corners[1] = points[i2]; corners[2] = points[i3];
				_rsampler.sample_uniformly_from_simplex(dart, 3, 3, corners);

				size_t iclosest(0); double hclosest(DBL_MAX);
				surface_point_cloud->get_closest_non_smooth_tree_point(dart, rmin, smooth_angle_threshold, iclosest, hclosest);
				size_t* attrib = surface_point_cloud->get_tree_point_attrib(iclosest);
				size_t face_index = attrib[1];

				size_t j1(faces[face_index][1]), j2(faces[face_index][2]), j3(faces[face_index][3]);
				corners[0] = points[j1]; corners[1] = points[j2]; corners[2] = points[j3];
				double proj_dst(0.0), d_mid(0.0);
				_geom.project_to_3d_triangle(dart, corners, proj, proj_dst);
				_geom.project_to_3d_triangle(mid_edge, corners, proj, d_mid);
				if (proj_dst > rmin && d_mid < 0.1 * rmin)
				{
					// A sharp edge: Two non-cosmooth faces that have a significant gap inbetween yet the distance betwen the edge mid point and jface is not that large
					size_t* sharp_edge = new size_t[2];
					if (i_index < ip_index) { sharp_edge[0] = i_index; sharp_edge[1] = ip_index; }
					else { sharp_edge[0] = ip_index; sharp_edge[1] = i_index; }
					_memo.add_entry(sharp_edge, num_sharp_edges, sharp_edges, cap_sharp_edges);
					break;
				}
			}
		}
		#pragma endregion
	}
	delete[] mid_edge;

	vcm_cout << "  * " << num_sharp_edges << " sharp edges has been detected" << std::endl;

	// Edge Point Cloud
	double total_length(0.0);
	for (size_t iedge = 0; iedge < num_sharp_edges; iedge++)
	{
		size_t i1 = sharp_edges[iedge][0]; size_t i2 = sharp_edges[iedge][1];
		total_length += _geom.distance(3, points[i1], points[i2]);
	}

	size_t num_edge_points = 100000;
	double edge_spacing = total_length / num_edge_points;
	double min_feature_size_edges = edge_spacing;

	edge_point_cloud->disable_auto_balance();
	double* tangent = new double[3];
	for (size_t iedge = 0; iedge < num_sharp_edges; iedge++)
	{
		#pragma region Generate Edge Point Cloud:
		size_t i1 = sharp_edges[iedge][0]; size_t i2 = sharp_edges[iedge][1];
		double edge_length = _geom.distance(3, points[i1], points[i2]);
		for (size_t idim = 0; idim < 3; idim++) tangent[idim] = points[i2][idim] - points[i1][idim];
		_geom.normalize_vector(3, tangent);
		size_t* attrib = new size_t[2]; attrib[0] = 2; attrib[1] = iedge;
		edge_point_cloud->add_tree_point(3, points[i1], tangent, attrib);
		edge_point_cloud->add_tree_point(3, points[i2], tangent, attrib);

		size_t num(size_t(edge_length / edge_spacing));
		if (num < 5) num = 5;
		double dx = edge_length / (num + 1);
		for (size_t ipoint = 1; ipoint <= num; ipoint++)
		{
			for (size_t idim = 0; idim < 3; idim++) dart[idim] = points[i1][idim] + ipoint * dx * tangent[idim];
			edge_point_cloud->add_tree_point(3, dart, tangent, attrib);
		}
		delete[] attrib;
		#pragma endregion
	}
	edge_point_cloud->kd_tree_build_balanced();
	//edge_point_cloud->write_points_csv("edge_point_cloud.csv", 3);
	vcm_cout << "  * Edge point cloud has " << edge_point_cloud->get_num_tree_points() << " points" << std::endl;

	bool* point_on_sharp_corner = new bool[num_points];
	for (size_t ipoint = 0; ipoint < num_points; ipoint++) point_on_sharp_corner[ipoint] = false;

	bool* point_on_sharp_edge = new bool[num_points];
	for (size_t ipoint = 0; ipoint < num_points; ipoint++) point_on_sharp_edge[ipoint] = false;

	// Sharp corners
	size_t cap_sharp_corners(100); num_sharp_corners = 0;
	sharp_corners = new size_t[cap_sharp_corners];

	for (size_t iedge = 0; iedge < num_sharp_edges; iedge++)
	{
		#pragma region Identify Corner points due to sharp edges:
		size_t i1 = sharp_edges[iedge][0]; size_t i2 = sharp_edges[iedge][1];
		point_on_sharp_edge[i1] = true; point_on_sharp_edge[i2] = true;

		double edge_length = _geom.distance(3, points[i1], points[i2]);

		if (true)
		{
			for (size_t idim = 0; idim < 3; idim++) dart[idim] = 0.5 * (points[i1][idim] + points[i2][idim]);
			size_t iclosest(0); double hclosest(DBL_MAX);
			edge_point_cloud->get_closest_non_smooth_tree_edge_point(dart, 1E-10, smooth_angle_threshold, iclosest, hclosest);
			if (hclosest > 0.5 * edge_length) continue;
		}

		if (2.0 * edge_spacing < edge_length) edge_length = 2.0 * edge_spacing;



		for (size_t idim = 0; idim < 3; idim++) tangent[idim] = points[i2][idim] - points[i1][idim];
		_geom.normalize_vector(3, tangent);

		size_t num_mc = 100;
		for (size_t imc = 0; imc < num_mc; imc++)
		{
			double u = _rsampler.generate_uniform_random_number();
			if (imc < num_mc / 2)
			{
				if (imc == 0) u = 0;
				for (size_t idim = 0; idim < 3; idim++) dart[idim] = points[i1][idim] + u * edge_length * tangent[idim];
			}
			else
			{
				if (imc == num_mc / 2) u = 0;
				for (size_t idim = 0; idim < 3; idim++) dart[idim] = points[i2][idim] - u * edge_length * tangent[idim];
			}

			if (true)
			{
				#pragma region Check if a non-cosmooth edge passes too close to the edge endpoints:
				size_t iclosest(0); double hclosest(DBL_MAX);
				edge_point_cloud->get_closest_non_smooth_tree_edge_point(dart, 1E-10, smooth_angle_threshold, iclosest, hclosest);

				size_t* attrib = edge_point_cloud->get_tree_point_attrib(iclosest);
				size_t edge_index = attrib[1];
				size_t j1(sharp_edges[edge_index][0]), j2(sharp_edges[edge_index][1]);

				double proj_dst(0.0), d_1(0.0), d_2(0.0);
				_geom.project_to_3d_line_segment(dart, points[j1], points[j2], proj, proj_dst);
				_geom.project_to_3d_line_segment(points[i1], points[j1], points[j2], proj, d_1);
				_geom.project_to_3d_line_segment(points[i2], points[j1], points[j2], proj, d_2);

				if (proj_dst > rmin && d_1 < 0.1 * rmin && !point_on_sharp_corner[i1])
				{
					// Two non-cosmooth sharp edges meeting at a sharp Corner i1
					_memo.add_entry(i1, num_sharp_corners, sharp_corners, cap_sharp_corners);
					point_on_sharp_corner[i1] = true;
				}
				if (proj_dst > rmin && d_2 < 0.1 * rmin && !point_on_sharp_corner[i2])
				{
					// Two non-cosmooth sharp edges meeting at a sharp Corner i2
					_memo.add_entry(i2, num_sharp_corners, sharp_corners, cap_sharp_corners);
					point_on_sharp_corner[i2] = true;
				}
				#pragma endregion
			}

			if (true)
			{
				#pragma region Check if a non-cosmooth face passes too close to the edge endpointss:
				size_t iclosest(0); double hclosest(DBL_MAX);
				surface_point_cloud->get_closest_non_smooth_tree_point(dart, rmin, smooth_angle_threshold, iclosest, hclosest);

				size_t* attrib = surface_point_cloud->get_tree_point_attrib(iclosest);
				size_t face_index = attrib[1];

				size_t j1(faces[face_index][1]), j2(faces[face_index][2]), j3(faces[face_index][3]);
				corners[0] = points[j1]; corners[1] = points[j2]; corners[2] = points[j3];

				double proj_dst(0.0), d_1(0.0), d_2(0.0);
				_geom.project_to_3d_triangle(dart, corners, proj, proj_dst);
				_geom.project_to_3d_triangle(points[i1], corners, proj, d_1);
				_geom.project_to_3d_triangle(points[i2], corners, proj, d_2);

				if (proj_dst > rmin && d_1 < 0.1 * rmin && !point_on_sharp_corner[i1])
				{
					// A sharp Corner
					_memo.add_entry(i1, num_sharp_corners, sharp_corners, cap_sharp_corners);
					point_on_sharp_corner[i1] = true;
				}
				if (proj_dst > rmin && d_2 < 0.1 * rmin && !point_on_sharp_corner[i2])
				{
					// A sharp Corner
					_memo.add_entry(i2, num_sharp_corners, sharp_corners, cap_sharp_corners);
					point_on_sharp_corner[i2] = true;
				}
				#pragma endregion
			}
		}
		#pragma endregion
	}

	for (size_t ipoint = 0; ipoint < num_points; ipoint++)
	{
		#pragma region Identify Corner Points due to Surface Normals:
		if (point_on_sharp_corner[ipoint]) continue;

		if (num_normals[ipoint] > 2)
		{
			// A sharp Corner
			_memo.add_entry(ipoint, num_sharp_corners, sharp_corners, cap_sharp_corners);
			point_on_sharp_corner[ipoint] = true;
		}
		else if (!point_on_sharp_edge[ipoint] && num_normals[ipoint] > 1)
		{
			// A sharp Corner that has no sharp edges
			_memo.add_entry(ipoint, num_sharp_corners, sharp_corners, cap_sharp_corners);
			point_on_sharp_corner[ipoint] = true;
		}
		#pragma endregion
	}

	vcm_cout << "  * " << num_sharp_corners << " sharp corners has been detected" << std::endl;

	for (size_t ipoint = 0; ipoint < num_points; ipoint++)
	{
		size_t num = num_normals[ipoint];
		for (size_t i = 0; i < num; i++) delete[] point_normals[ipoint][i];
		delete[] point_normals[ipoint];
	}
	delete[] point_normals; delete[] cap_normals; delete[] num_normals;  delete[] point_visited;
	delete[] point_on_sharp_corner; delete[] point_on_sharp_edge;
	for (size_t iface = 0; iface < num_faces; iface++) delete[] face_normal[iface];
	delete[] face_normal; delete[] corners;
	delete[] dart; delete[] proj; delete[] tangent;

	end_time = clock();
	cpu_time = ((double)(end_time - start_time)) / CLOCKS_PER_SEC;
	vcm_cout << "  * executed in " << cpu_time << " seconds!" << std::endl;
	return 0;
	#pragma endregion
}



int MeshingVoroCrustSampler::generate_corner_spheres(size_t num_points, double** points, size_t num_sharp_corners, size_t* sharp_corners,
	                                          MeshingSmartTree* surface_point_cloud, MeshingSmartTree* edge_point_cloud, double feature_TOL,
	                                          MeshingSmartTree* corner_spheres,
	                                          double smooth_angle_threshold, double rmax, double Lip)
{
	#pragma region Generating Corner Spheres:
	vcm_cout << "VoroCrust::Generating Corner Spheres:" << std::endl;

	clock_t start_time, end_time; double cpu_time;
	start_time = clock();

	size_t num_samples(0);
	double* new_sphere = new double[4];
	size_t* attrib = new size_t[2]; attrib[0] = 2;
	double rmin(DBL_MAX); size_t imin(0);
	for (size_t icorner = 0; icorner < num_sharp_corners; icorner++)
	{
		size_t corner_index = sharp_corners[icorner];

		for (size_t idim = 0; idim < 3; idim++) new_sphere[idim] = points[corner_index][idim];

		if (icorner > 0)
		{
			size_t iclosest(0); double hclosest(DBL_MAX);
			corner_spheres->get_closest_tree_point(new_sphere, iclosest, hclosest);
			if (hclosest < feature_TOL) continue; // redundant point
		}

		size_t iclosest(0); double hclosest(DBL_MAX);
		surface_point_cloud->get_closest_non_smooth_tree_point(new_sphere, feature_TOL, smooth_angle_threshold, iclosest, hclosest);
		double surf_sz = alpha_sz * hclosest;

		hclosest = DBL_MAX;
		edge_point_cloud->get_closest_non_smooth_tree_edge_point(new_sphere, feature_TOL, smooth_angle_threshold, iclosest, hclosest);
		double edge_sz = alpha_sz * hclosest;

		new_sphere[3] = fmin(surf_sz, edge_sz);

		if (new_sphere[3] > rmax) new_sphere[3] = rmax;

		if (_input_sizing_function != 0)
		{
			size_t iclosest; double hclosest(DBL_MAX);
			_input_sizing_function->get_closest_tree_point(new_sphere, iclosest, hclosest);
			double* x_sz = _input_sizing_function->get_tree_point(iclosest);
			double r_sz = x_sz[3];
			if (new_sphere[3] > r_sz) new_sphere[3] = r_sz;
		}

		if (new_sphere[3] < rmin)
		{
			rmin = new_sphere[3];
			imin = num_samples;
		}

		attrib[1] = corner_index;
		corner_spheres->add_tree_point(4, new_sphere, 0, attrib);
		num_samples++;
	}

	for (size_t icorner = 0; icorner < num_samples; icorner++)
	{
		double* corner_sphere = corner_spheres->get_tree_point(icorner);
		size_t iclosest(0); double hclosest(DBL_MAX);
		corner_spheres->get_closest_tree_point(icorner, iclosest, hclosest);
		if (corner_sphere[3] > alpha_sz * hclosest) corner_sphere[3] = alpha_sz * hclosest;
	}
	delete[] attrib; delete[] new_sphere;

	bool shrunk;
	_smethods.impose_lipschitz_continuity_brute(corner_spheres, Lip, shrunk);


	double min_radius(DBL_MAX), max_radius(0.0);
	for (size_t i = 0; i < num_samples; i++)
	{
		double* sphere = corner_spheres->get_tree_point(i);
		if (sphere[3] < min_radius)
			min_radius = sphere[3];
		if (sphere[3] > max_radius)
			max_radius = sphere[3];
	}

	end_time = clock();
	cpu_time = ((double)(end_time - start_time)) / CLOCKS_PER_SEC;

	corner_spheres->save_tree_csv("corner_spheres.csv", 4);

	vcm_cout << "  * Number of Spheres = " << corner_spheres->get_num_tree_points() << std::endl;
	vcm_cout << "  * Minimum Corner Sphere Radius = " << min_radius << std::endl;
	vcm_cout << "  * Maximum Corner Sphere Radius = " << max_radius << std::endl;
	vcm_cout << "  * Corner spheres are saved corner_spheres.csv" << std::endl;
	vcm_cout << "  * executed in " << cpu_time << " seconds!" << std::endl;
	return 0;
	#pragma endregion
}


int MeshingVoroCrustSampler::generate_edge_spheres(size_t num_points, double** points, size_t num_edges, size_t** edges,
	                                        MeshingSmartTree* surface_point_cloud, MeshingSmartTree* edge_point_cloud, double feature_TOL,
	                                        MeshingSmartTree* edge_spheres, MeshingSmartTree* corner_spheres,
	                                        double smooth_angle_threshold, double rmax, double Lip, double coverage_radius_ratio)
{
	#pragma region Generate Edge Spheres:

	vcm_cout << "VoroCrust::Generating Edge Spheres:" << std::endl;

	clock_t start_time, end_time; double cpu_time;
	start_time = clock();

	double min_distance_ratio = sqrt(1.0 - coverage_radius_ratio * coverage_radius_ratio);

	size_t num_attrib = 2;
	size_t* attrib = new size_t[num_attrib];
	double* new_sphere = new double[4];
	double* proj = new double[3];
	double* tangent = new double[3];
	double* vec = new double[3];

	size_t num_samples(0);
	double rmin(DBL_MAX); size_t imin(0);

	size_t nit(0);
	while (true)
	{
		nit++;
		initiate_active_pool_edges(num_points, points, num_edges, edges);

		size_t num_successive_misses(0), max_num_successive_misses(100);
		while (true)
		{
			#pragma region Simple MPS on Edges:
			size_t sphere_edge;
			sample_active_edges_uniformly(num_points, points, num_edges, edges, new_sphere, sphere_edge);

			bool covered(false);
			if (!covered) covered = _smethods.point_covered(new_sphere, corner_spheres, Lip, 0.0);
			if (!covered) covered = _smethods.point_covered(new_sphere, edge_spheres, Lip, coverage_radius_ratio);

			bool valid_dart(true);
			if (covered) valid_dart = false;

			if (valid_dart)
			{
				#pragma region Create an edge sphere:

				size_t iclosest(0); double hclosest(DBL_MAX);
				surface_point_cloud->get_closest_non_smooth_tree_point(new_sphere, feature_TOL, smooth_angle_threshold, iclosest, hclosest);
				double surf_sz = alpha_sz * hclosest;

				hclosest = DBL_MAX;
				edge_point_cloud->get_closest_non_smooth_tree_edge_point(new_sphere, feature_TOL, smooth_angle_threshold, iclosest, hclosest);
				double edge_sz = alpha_sz * hclosest;

				hclosest = DBL_MAX;
				corner_spheres->get_closest_non_smooth_tree_edge_point(new_sphere, feature_TOL, smooth_angle_threshold, iclosest, hclosest);
				edge_sz = fmin(edge_sz, alpha_sz * hclosest);

				if (num_samples > 1)
				{
					size_t num_nearby_spheres_s(0); size_t* nearby_spheres_s(0);
					_smethods.get_overlapping_spheres(new_sphere, edge_spheres, Lip, num_nearby_spheres_s, nearby_spheres_s);
					for (size_t i = 0; i < num_nearby_spheres_s; i++)
					{
						size_t sphere_index = nearby_spheres_s[i];
						double* near_by_sphere = edge_spheres->get_tree_point(sphere_index);
						double dst = _geom.distance(3, new_sphere, near_by_sphere);
						double coverage_sz = dst / min_distance_ratio;
						edge_sz = fmin(edge_sz, coverage_sz);            // new sphere should not alpha cover an old sphere center
						double Lip_sz = near_by_sphere[3] + Lip * dst;
						edge_sz = fmin(edge_sz, Lip_sz);                 // respect Lipschitz constant locally
					}
					delete[] nearby_spheres_s;
				}

				new_sphere[3] = fmin(surf_sz, edge_sz);
				if (new_sphere[3] > rmax) new_sphere[3] = rmax;

				if (new_sphere[3] < _r_min_edges) new_sphere[3] = _r_min_edges;

				if (_input_sizing_function != 0)
				{
					size_t iclosest; double hclosest(DBL_MAX);
					_input_sizing_function->get_closest_tree_point(new_sphere, iclosest, hclosest);
					double* x_sz = _input_sizing_function->get_tree_point(iclosest);
					double r_sz = x_sz[3];
					if (new_sphere[3] > r_sz) new_sphere[3] = r_sz;
				}

				attrib[0] = 2; attrib[1] = sphere_edge;
				edge_spheres->add_tree_point(4, new_sphere, 0, attrib);

				if (new_sphere[3] < rmin)
				{
					rmin = new_sphere[3];
					imin = num_samples;
				}

				num_samples++;

				if (num_samples / 1000 * 1000 == num_samples)
				{
					vcm_cout << "  * Number of Spheres = (" << corner_spheres->get_num_tree_points();
					vcm_cout << ", " << edge_spheres->get_num_tree_points() << ")" << std::endl;
				}
				#pragma endregion
			}

			if (!valid_dart) num_successive_misses++;
			else             num_successive_misses = 0;

			if (num_successive_misses == max_num_successive_misses)
			{
				refine_active_pool_edges(num_points, points, num_edges, edges, edge_spheres, corner_spheres, Lip, coverage_radius_ratio);
				if (_active_pool_size == 0) break;
				num_successive_misses = 0;
			}
			#pragma endregion
		}

		bool sphere_shrunk_Lip(false);
		_smethods.impose_lipschitz_continuity(edge_spheres, Lip, sphere_shrunk_Lip);
		if (sphere_shrunk_Lip) continue;


		bool spheres_shrunk(false);
		for (size_t sphere_i_index = 0; sphere_i_index < num_samples; sphere_i_index++)
		{
			#pragma region Validate proper sizing in extremely narrow regions:
			double* sphere_i = edge_spheres->get_tree_point(sphere_i_index);
			size_t* sphere_i_attrib = edge_spheres->get_tree_point_attrib(sphere_i_index);
			size_t edge_i_index = sphere_i_attrib[1];

			size_t i1(edges[edge_i_index][0]), i2(edges[edge_i_index][1]);

			 // Collect spheres overlapping with isphere:
			size_t num_nearby_spheres(0); size_t* nearby_spheres(0);
			_smethods.get_overlapping_spheres(sphere_i, edge_spheres, Lip, num_nearby_spheres, nearby_spheres);

			for (size_t ii = 0; ii < num_nearby_spheres; ii++)
			{
				size_t sphere_j_index = nearby_spheres[ii];
				size_t* sphere_attrib = edge_spheres->get_tree_point_attrib(sphere_j_index);
				size_t edge_j_index = sphere_attrib[1];
				if (edge_i_index == edge_j_index) continue;

				size_t j1(edges[edge_j_index][0]), j2(edges[edge_j_index][1]);
				double proj_dst(0.0);

				_geom.project_to_3d_line_segment(sphere_i, points[j1], points[j2], proj, proj_dst);

				for (size_t idim = 0; idim < 3; idim++) tangent[idim] = points[j2][idim] - points[j1][idim];
				for (size_t idim = 0; idim < 3; idim++) vec[idim] = sphere_i[idim] - proj[idim];

				if (_geom.normalize_vector(3, tangent) && _geom.normalize_vector(3, vec))
				{
					double dot = _geom.dot_product(3, tangent, vec);
					if (fabs(dot) > smooth_angle_threshold) continue;

					proj_dst *= alpha_sz;
					if (proj_dst < feature_TOL) continue;
					if (proj_dst < _r_min_edges) proj_dst = _r_min_edges;
					if (sphere_i[3] < (1 + 1E-8) * proj_dst) continue;

					sphere_i[3] = proj_dst; spheres_shrunk = true;
				}
			}

			delete[] nearby_spheres;
			#pragma endregion
		}
		if (spheres_shrunk) continue;
		break;
	}

	double min_radius(DBL_MAX), max_radius(0.0);
	for (size_t i = 0; i < num_samples; i++)
	{
		double* sphere = edge_spheres->get_tree_point(i);
		if (sphere[3] < min_radius)
			min_radius = sphere[3];
		if (sphere[3] > max_radius)
			max_radius = sphere[3];
	}

	end_time = clock();
	cpu_time = ((double)(end_time - start_time)) / CLOCKS_PER_SEC;

	edge_spheres->save_tree_csv("edges_spheres.csv", 4);

	vcm_cout << "  * Number of Spheres = (" << corner_spheres->get_num_tree_points();
	vcm_cout << ", " << edge_spheres->get_num_tree_points() << ")" << std::endl;
	vcm_cout << "  * Minimum Edge Sphere Radius = " << min_radius << std::endl;
	vcm_cout << "  * Maximum Edge Sphere Radius = " << max_radius << std::endl;
	vcm_cout << "  * Edge spheres are saved edges_spheres.csv" << std::endl;
	vcm_cout << "  * executed in " << cpu_time << " seconds!" << std::endl;

	delete[] new_sphere; delete[] attrib; delete[] proj; delete[] vec; delete[] tangent;
	return 0;
	#pragma endregion
}



int MeshingVoroCrustSampler::generate_surface_spheres(size_t num_points, double** points, size_t num_faces, size_t** faces,
	                                           MeshingSmartTree* surface_point_cloud, MeshingSmartTree* edge_point_cloud, double R_MIN,
	                                           MeshingSmartTree* &surface_spheres, MeshingSmartTree* edge_spheres, MeshingSmartTree* corner_spheres,
	                                           MeshingSmartTree* &ext_surf_seeds, MeshingSmartTree* &int_surf_seeds,
	                                           double smooth_angle_threshold, double rmax, double Lip, double coverage_radius_ratio)
{
	#pragma region Solve MPS over Triangular Surfaces:

	vcm_cout << "VoroCrust::Generating Surface Spheres:" << std::endl;

	clock_t start_time, end_time; double cpu_time;
	start_time = clock();

	clock_t mps_start_time, mps_end_time; double mps_cpu_time(0.0);
	clock_t sliver_start_time, sliver_end_time; double sliver_cpu_time(0.0);
	clock_t Lip_start_time, Lip_end_time; double Lip_cpu_time(0.0);

	size_t num_patches(100);

	size_t* new_sphere_face = new size_t[_num_threads];
	bool* covered_sample = new bool[_num_threads];
	double** new_sphere = new double*[_num_threads];
	
#if defined USE_OPEN_MP
	omp_set_num_threads(int(_num_threads));
	#pragma omp parallel
	{
		int thread_id = omp_get_thread_num();
#else
	for (size_t thread_id = 0; thread_id < _num_threads; thread_id++)
	{
#endif

		new_sphere[thread_id] = new double[4];
	}

	size_t num_samples(0);

	double** corners = new double*[3]; double* proj = new double[3];
	double* vec = new double[3]; double* normal = new double[3];

	double min_radius(DBL_MAX); size_t imin(0);
	while (true)
	{
		bool surface_sphere_shrunk(false), edge_sphere_shrunk(false), corner_sphere_shrunk(false);

		mps_start_time = clock();
		for (size_t ipatch = 0; ipatch < num_patches; ipatch++)
		{
			initiate_active_pool(num_points, points, num_faces, faces, ipatch, num_patches);

			size_t num_misses(0), max_num_misses(100);
			while (true)
			{
				#pragma region Simple MPS:
				
				sample_active_pool(num_samples, num_points, points, num_faces, faces,
					               surface_point_cloud, edge_point_cloud, R_MIN,
					               surface_spheres, edge_spheres, corner_spheres,
					               smooth_angle_threshold, rmax, Lip, coverage_radius_ratio,
					               new_sphere, new_sphere_face, covered_sample);

				// adding samples
				for (size_t thread_id = 0; thread_id < _num_threads; thread_id++)
				{
					if (covered_sample[thread_id])
					{
						num_misses++;
						if (num_misses >= max_num_misses) break;
					}
					else if (!covered_sample[thread_id])
					{
						size_t num_attrib = 2;
						size_t* attrib = new size_t[num_attrib];
						attrib[0] = num_attrib;

						attrib[1] = new_sphere_face[thread_id];
						surface_spheres->add_tree_point(4, new_sphere[thread_id], 0, attrib);

						if (new_sphere[thread_id][1] < _ymin)
						{
							_ymin = new_sphere[thread_id][1];
						}


						delete[] attrib;

						num_samples++; 
						num_misses = 0; // any new sample rest number of misses

						if (new_sphere[thread_id][3] < min_radius)
						{
							min_radius = new_sphere[thread_id][3];
							imin = num_samples;
						}

						if (num_samples / 5000 * 5000 == num_samples)
						{
							vcm_cout << "  * Number of Spheres = (" << corner_spheres->get_num_tree_points();
							vcm_cout << ", " << edge_spheres->get_num_tree_points();
							vcm_cout << ", " << surface_spheres->get_num_tree_points() << ")" << std::endl;
						}
					}
				}

				if (num_misses >= max_num_misses)
				{
					refine_active_pool(num_points, points, num_faces, faces, surface_spheres, edge_spheres, corner_spheres, Lip, coverage_radius_ratio);
					if (_active_pool_size == 0) break;
					num_misses = 0;
				}
				#pragma endregion
			}
		}
		mps_end_time = clock();
		mps_cpu_time += ((double)(mps_end_time - mps_start_time)) / CLOCKS_PER_SEC;

		if (_r_min_surface > 0.0 && false)
		{
			#pragma region process singular spheres:
			// DEBUG VC WILD
			surface_spheres->save_tree_csv("surface_spheres.csv", 4);
			for (size_t i = 0; i < surface_spheres->get_num_tree_points(); i++)
			{
				double* ss = surface_spheres->get_tree_point(i);
				if (fabs(ss[3] - _r_min_surface) < 1E-10) ss[3] = 0.0;
			}
			surface_spheres->save_tree_csv("good_surface_spheres.csv", 4);
			for (size_t i = 0; i < surface_spheres->get_num_tree_points(); i++)
			{
				double* ss = surface_spheres->get_tree_point(i);
				if (fabs(ss[3]) < 1E-10) ss[3] = _r_min_surface;
				else ss[3] = 0.0;
			}
			surface_spheres->save_tree_csv("singular_surface_spheres.csv", 4);
			#pragma endregion
		}
		
		if (R_MIN > 0.0)
		{
			bool spheres_shrunk_sizing(false);
			for (size_t sphere_i_index = 0; sphere_i_index < num_samples; sphere_i_index++)
			{
				#pragma region Validate proper sizing in extremely narrow regions:
				double* sphere_i = surface_spheres->get_tree_point(sphere_i_index);
				size_t* sphere_i_attrib = surface_spheres->get_tree_point_attrib(sphere_i_index);
				size_t face_i_index = sphere_i_attrib[1];

				size_t i1(faces[face_i_index][1]), i2(faces[face_i_index][2]), i3(faces[face_i_index][2]);

				// Collect spheres overlapping with isphere:
				size_t num_nearby_spheres(0); size_t* nearby_spheres(0);
				_smethods.get_overlapping_spheres(sphere_i, surface_spheres, Lip, num_nearby_spheres, nearby_spheres);

				for (size_t ii = 0; ii < num_nearby_spheres; ii++)
				{
					size_t sphere_j_index = nearby_spheres[ii];
					size_t* sphere_attrib = surface_spheres->get_tree_point_attrib(sphere_j_index);
					size_t face_j_index = sphere_attrib[1];
					if (face_i_index == face_j_index) continue;

					size_t j1(faces[face_j_index][1]), j2(faces[face_j_index][2]), j3(faces[face_j_index][2]);
					corners[0] = points[j1]; corners[1] = points[j2]; corners[2] = points[j3];

					double proj_dst(0.0);
					_geom.project_to_3d_triangle(sphere_i, corners, proj, proj_dst);

					_geom.get_3d_triangle_normal(corners, normal);

					for (size_t idim = 0; idim < 3; idim++) vec[idim] = sphere_i[idim] - proj[idim];

					if (_geom.normalize_vector(3, vec))
					{
						double dot = _geom.dot_product(3, normal, vec);
						dot = sqrt(1 - dot * dot);
						if (fabs(dot) > smooth_angle_threshold) continue;

						proj_dst *= alpha_sz;
						if (proj_dst < R_MIN) continue;

						if (proj_dst < _r_min_surface) proj_dst = _r_min_surface;

						if (sphere_i[3] < (1 + 1E-8) * proj_dst) continue;
						sphere_i[3] = proj_dst; spheres_shrunk_sizing = true;
					}
				}
				delete[] nearby_spheres;
				#pragma endregion
			}
			if (spheres_shrunk_sizing) continue;
		}

		bool sphere_shrunk_Lip(false);
		Lip_start_time = clock();
		_smethods.impose_lipschitz_continuity(surface_spheres, Lip, sphere_shrunk_Lip);
		Lip_end_time = clock();
		Lip_cpu_time += ((double)(Lip_end_time - Lip_start_time)) / CLOCKS_PER_SEC;
		if (sphere_shrunk_Lip) continue;

		// We have maximal coverage and Lipschitz continuity at this point
		vcm_cout << "  * Number of Spheres = (" << corner_spheres->get_num_tree_points();
		vcm_cout << ", " << edge_spheres->get_num_tree_points();
		vcm_cout << ", " << surface_spheres->get_num_tree_points() << ")" << std::endl;

		bool sphere_shrunk_slivers;
		sliver_start_time = clock();
		generate_surface_seeds(num_points, points, num_faces, faces, surface_spheres, edge_spheres, corner_spheres, Lip, sphere_shrunk_slivers, ext_surf_seeds, int_surf_seeds);
		sliver_end_time = clock();
		sliver_cpu_time += ((double)(sliver_end_time - sliver_start_time)) / CLOCKS_PER_SEC;
		if (sphere_shrunk_slivers) continue;


		size_t num_surf_faces = int_surf_seeds->get_num_tree_points();
		size_t** surf_faces = new size_t*[num_surf_faces];
		for (size_t iface = 0; iface < num_surf_faces; iface++)
		{
			size_t* attrib = int_surf_seeds->get_tree_point_attrib(iface);
			surf_faces[iface] = new size_t[4];
			surf_faces[iface][0] = 3;
			surf_faces[iface][1] = attrib[2];
			surf_faces[iface][2] = attrib[3];
			surf_faces[iface][3] = attrib[4];
		}
		bool water_tight = validate_surface_mesh(surface_spheres, edge_spheres, corner_spheres, num_surf_faces, surf_faces);
		for (size_t iface = 0; iface < num_surf_faces; iface++) delete[] surf_faces[iface];
		delete[] surf_faces;
		if (water_tight) break;

	}

#if defined USE_OPEN_MP
	omp_set_num_threads(int(_num_threads));
	#pragma omp parallel
	{
	int thread_id = omp_get_thread_num();
#else
	for (size_t thread_id = 0; thread_id < _num_threads; thread_id++)
	{
#endif
		delete[] new_sphere[thread_id];
	}
	delete[] new_sphere_face; delete[] covered_sample; delete[] new_sphere;

	min_radius = DBL_MAX; double max_radius(0.0);
	for (size_t i = 0; i < num_samples; i++)
	{
		double* sphere = surface_spheres->get_tree_point(i);
		if (sphere[3] < min_radius) min_radius = sphere[3];
		if (sphere[3] > max_radius) max_radius = sphere[3];
	}
	end_time = clock();
	cpu_time = ((double)(end_time - start_time)) / CLOCKS_PER_SEC;

	surface_spheres->save_tree_csv("surface_spheres.csv", 4);

	vcm_cout << "  * Minimum Surface Sphere Radius = " << min_radius << std::endl;
	vcm_cout << "  * Maximum Surface Sphere Radius = " << max_radius << std::endl;
	vcm_cout << "  * Surface spheres are saved surface_spheres.csv" << std::endl;
	vcm_cout << "  * executed in " << cpu_time << " seconds!" << std::endl;
	vcm_cout << "    -> MPS Time = " << mps_cpu_time << " seconds!" << std::endl;
	vcm_cout << "       -> Pool Sampling Time = " << _pool_sampling_time << " seconds!" << std::endl;
	vcm_cout << "          -> Sizing Time = " << _sizing_time << " seconds!" << std::endl;
	vcm_cout << "       -> Pool Refinement Time = " << _pool_refinement_time << " seconds!" << std::endl;
	vcm_cout << "    -> Lip Time = " << Lip_cpu_time << " seconds!" << std::endl;
	vcm_cout << "    -> Sliver Time = " << sliver_cpu_time << " seconds!" << std::endl;

	delete[] corners; delete[] proj; delete[] vec; delete[] normal;
	return 0;
	#pragma endregion
}


int MeshingVoroCrustSampler::color_surface_seeds(MeshingSmartTree* surface_spheres, MeshingSmartTree* edge_spheres, MeshingSmartTree* corner_spheres, MeshingSmartTree* upper_seeds, MeshingSmartTree* lower_seeds,
	                                      MeshingSmartTree* spheres, MeshingSmartTree* seeds, size_t &num_subregions)
{
	#pragma region Color Surface Seeds:
	vcm_cout << "VoroCrust::Coloring Surface Seeds:" << std::endl;

	clock_t start_time, end_time; double cpu_time;
	start_time = clock();

	size_t num_surface_spheres(surface_spheres->get_num_tree_points());
	size_t num_edge_spheres(edge_spheres->get_num_tree_points());
	size_t num_corner_spheres(corner_spheres->get_num_tree_points());

	// merging spheres into one Tree
	for (size_t isphere = 0; isphere < num_surface_spheres; isphere++)
	{
		double* x = surface_spheres->get_tree_point(isphere);
		double* normal = surface_spheres->get_tree_point_normal(isphere);
		size_t* attrib = surface_spheres->get_tree_point_attrib(isphere);
		spheres->add_tree_point(4, x, normal, attrib);
	}
	//surface_spheres->clear_memory();

	for (size_t isphere = 0; isphere < num_edge_spheres; isphere++)
	{
		double* x = edge_spheres->get_tree_point(isphere);
		double* normal = edge_spheres->get_tree_point_normal(isphere);
		size_t* attrib = edge_spheres->get_tree_point_attrib(isphere);
		spheres->add_tree_point(4, x, normal, attrib);
	}
	//edge_spheres->clear_memory();

	for (size_t isphere = 0; isphere < num_corner_spheres; isphere++)
	{
		double* x = corner_spheres->get_tree_point(isphere);
		double* normal = corner_spheres->get_tree_point_normal(isphere);
		size_t* attrib = corner_spheres->get_tree_point_attrib(isphere);
		spheres->add_tree_point(4, x, normal, attrib);
	}
	//corner_spheres->clear_memory();


	// merging seeds into one Tree

	size_t num_upper_seeds(upper_seeds->get_num_tree_points());
	size_t num_lower_seeds(lower_seeds->get_num_tree_points());

	// adjust id of seed pair for upper seeds
	for (size_t iseed = 0; iseed < num_upper_seeds; iseed++)
	{
		// copy point location, normal and attributes
		double* x = upper_seeds->get_tree_point(iseed);
		double* normal = upper_seeds->get_tree_point_normal(iseed);
		size_t* attrib = upper_seeds->get_tree_point_attrib(iseed);
		attrib[1] += num_upper_seeds;
		seeds->add_tree_point(4, x, normal, attrib);
	}
	upper_seeds->clear_memory();

	for (size_t iseed = 0; iseed < num_lower_seeds; iseed++)
	{
		// copy point location, normal and attributes
		double* x = lower_seeds->get_tree_point(iseed);
		double* normal = lower_seeds->get_tree_point_normal(iseed);
		size_t* attrib = lower_seeds->get_tree_point_attrib(iseed);
		seeds->add_tree_point(4, x, normal, attrib);
	}
	lower_seeds->clear_memory();

	size_t num_seeds(seeds->get_num_tree_points());
	for (size_t iseed = 0; iseed < num_seeds; iseed++)
	{
		// Establish Sphere -> Seeds directional graph
		size_t* attrib = seeds->get_tree_point_attrib(iseed);
		size_t si(attrib[2]), sj(attrib[3]), sk(attrib[4]);
		spheres->graph_connect(si, iseed); spheres->graph_connect(sj, iseed); spheres->graph_connect(sk, iseed);
	}

	double* face_normal = new double[3];
	double** face_corners = new double*[3];
	size_t num_spheres(spheres->get_num_tree_points());
	for (size_t isphere = 0; isphere < num_surface_spheres; isphere++)
	{
		#pragma region Connect Seeds of Surface Spheres:
		size_t * sphere_seeds(0);
		spheres->graph_get_neighbors(isphere, sphere_seeds);
		if (sphere_seeds == 0) continue;

		size_t num_sphere_seeds = sphere_seeds[1];
		for (size_t i = 0; i < num_sphere_seeds; i++)
		{
			size_t seed_i(sphere_seeds[2 + i]);
			size_t* attrib = seeds->get_tree_point_attrib(seed_i);
			attrib[5] = 0;
		}

		size_t first_neighbor(isphere), jsphere(isphere);
		if (true)
		{
			#pragma region first face dictates orientation:
			size_t seed_1(sphere_seeds[2]);
			size_t* attrib = seeds->get_tree_point_attrib(seed_1);

			size_t seed_2(attrib[1]);
			size_t si(attrib[2]), sj(attrib[3]), sk(attrib[4]);

			while (si != isphere)
			{
				size_t tmp = si; si = sj; sj = sk; sk = tmp; // rotate indices
			}

			double* x_1 = seeds->get_tree_point(seed_1);
			double* x_2 = seeds->get_tree_point(seed_2);

			face_corners[0] = spheres->get_tree_point(si);
			face_corners[1] = spheres->get_tree_point(sj);
			face_corners[2] = spheres->get_tree_point(sk);

			_geom.get_3d_triangle_normal(face_corners, face_normal);

			double dot(0.0);
			for (size_t idim = 0; idim < 3; idim++) dot += (x_2[idim] - x_1[idim]) * face_normal[idim];

			if (dot > 0.0)
			{
				attrib[5] = 1; // inside
				attrib = seeds->get_tree_point_attrib(seed_2);
				attrib[5] = 2; // outside
			}
			else
			{
				attrib[5] = 2; // outside
				attrib = seeds->get_tree_point_attrib(seed_2);
				attrib[5] = 1; // inside
			}
			first_neighbor = sj;
			jsphere = sk;
			#pragma endregion
		}

		while (true)
		{
			#pragma region Mark loop seeds:
			bool done(true);

			for (size_t i = 0; i < num_sphere_seeds; i++)
			{
				size_t seed_1(sphere_seeds[2 + i]);
				size_t* attrib = seeds->get_tree_point_attrib(seed_1);

				size_t seed_2(attrib[1]);
				size_t si(attrib[2]), sj(attrib[3]), sk(attrib[4]);

				while (si != isphere)
				{
					size_t tmp = si; si = sj; sj = sk; sk = tmp; // rotate indices
				}

				if (attrib[5] != 0) continue; // pair is already marked

				if (sj != jsphere && sk != jsphere) continue; // wrong triangle

				done = false;
				break;
			}
			if (done) 
				break;

			for (size_t i = 0; i < num_sphere_seeds; i++)
			{
				size_t seed_1(sphere_seeds[2 + i]);
				size_t* attrib = seeds->get_tree_point_attrib(seed_1);

				size_t seed_2(attrib[1]);
				size_t si(attrib[2]), sj(attrib[3]), sk(attrib[4]);

				while (si != isphere)
				{
					size_t tmp = si; si = sj; sj = sk; sk = tmp; // rotate indices
				}

				if (attrib[5] != 0) continue; // pair is already marked

				if (sj != jsphere && sk != jsphere) continue; // wrong triangle

				if (sj != jsphere)
				{
					size_t tmp = sj; sj = sk; sk = tmp; // flip triangle
				}

				double* x_1 = seeds->get_tree_point(seed_1);
				double* x_2 = seeds->get_tree_point(seed_2);

				face_corners[0] = spheres->get_tree_point(si);
				face_corners[1] = spheres->get_tree_point(sj);
				face_corners[2] = spheres->get_tree_point(sk);

				_geom.get_3d_triangle_normal(face_corners, face_normal);

				double dot(0.0);
				for (size_t idim = 0; idim < 3; idim++) dot += (x_2[idim] - x_1[idim]) * face_normal[idim];

				if (dot > 0.0)
				{
					attrib[5] = 1; // inside
					attrib = seeds->get_tree_point_attrib(seed_2);
					attrib[5] = 2; // outside
				}
				else
				{
					attrib[5] = 2; // outside
					attrib = seeds->get_tree_point_attrib(seed_2);
					attrib[5] = 1; // inside
				}
				jsphere = sk;
				if (jsphere == first_neighbor) done = true;
				break;
			}
			if (done) break;
			#pragma endregion
		}

		for (size_t i = 0; i < num_sphere_seeds; i++)
		{
			#pragma region connect Sphere seeds:
			size_t seed_i(sphere_seeds[2 + i]);
			size_t* attrib_i = seeds->get_tree_point_attrib(seed_i);
			for (size_t j = i + 1; j < num_sphere_seeds; j++)
			{
				size_t seed_j(sphere_seeds[2 + j]);
				size_t* attrib_j = seeds->get_tree_point_attrib(seed_j);

				if (attrib_i[5] == attrib_j[5])
				{
					seeds->graph_connect_nodes(seed_i, seed_j);
				}
			}
			#pragma endregion
		}

		for (size_t i = 0; i < num_sphere_seeds; i++)
		{
			size_t seed_i(sphere_seeds[2 + i]);
			size_t* attrib = seeds->get_tree_point_attrib(seed_i);
			attrib[5] = 0;
		}
		#pragma endregion
	}

	double* e1 = new double[3];
	double* e2 = new double[3];
	double* e3 = new double[3];

	for (size_t s_index = 0; s_index < num_edge_spheres; s_index++)
	{
		#pragma region Connect Seeds of Edge Spheres:
		size_t isphere(s_index + num_surface_spheres);

		size_t* sphere_seeds(0);
		spheres->graph_get_neighbors(isphere, sphere_seeds);
		if (sphere_seeds == 0) continue;

		size_t num_sphere_seeds = sphere_seeds[1];
		for (size_t i = 0; i < num_sphere_seeds; i++)
		{
			size_t seed_i(sphere_seeds[2 + i]);
			size_t* attrib = seeds->get_tree_point_attrib(seed_i);
			attrib[5] = 0;
		}

		size_t isphere_neighbor(isphere);
		for (size_t i = 0; i < num_sphere_seeds; i++)
		{
			size_t seed_1(sphere_seeds[2 + i]);
			size_t* attrib = seeds->get_tree_point_attrib(seed_1);

			size_t seed_2(attrib[1]);
			size_t si(attrib[2]), sj(attrib[3]), sk(attrib[4]);

			while (si != isphere)
			{
				size_t tmp = si; si = sj; sj = sk; sk = tmp; // rotate indices
			}

			if (sj >= num_surface_spheres)
			{
				isphere_neighbor = sj; break;
			}

			if (sk >= num_surface_spheres)
			{
				isphere_neighbor = sk; break;
			}
		}

		size_t num_loops(0);
		for (size_t i = 0; i < num_sphere_seeds; i++)
		{
			size_t seed_1(sphere_seeds[2 + i]);
			size_t* attrib = seeds->get_tree_point_attrib(seed_1);

			size_t seed_2(attrib[1]);
			size_t si(attrib[2]), sj(attrib[3]), sk(attrib[4]);

			while (si != isphere)
			{
				size_t tmp = si; si = sj; sj = sk; sk = tmp; // rotate indices
			}

			if (sj == isphere_neighbor || sk == isphere_neighbor) num_loops++;
		}
		num_loops /= 2; // because every triangle was counted twice

		size_t* jsphere = new size_t[num_loops]; num_loops = 0;
		for (size_t i = 0; i < num_sphere_seeds; i++)
		{
			#pragma region First loop Layer:
			size_t seed_1(sphere_seeds[2 + i]);
			size_t* attrib = seeds->get_tree_point_attrib(seed_1);

			size_t seed_2(attrib[1]);
			size_t si(attrib[2]), sj(attrib[3]), sk(attrib[4]);

			while (si != isphere)
			{
				size_t tmp = si; si = sj; sj = sk; sk = tmp; // rotate indices
			}

			if (sj != isphere_neighbor && sk != isphere_neighbor) continue;

			if (sj != isphere_neighbor)
			{
				size_t tmp = sj; sj = sk; sk = tmp; // flip indices
			}

			bool found(false);
			for (size_t j = 0; j < num_loops; j++)
			{
				if (jsphere[j] == sk)
				{
					found = true; break;
				}
			}
			if (found) continue;

			jsphere[num_loops] = sk; num_loops++;
			#pragma endregion
		}

		// form orthogonal basis
		double norm(0.0);
		for (size_t idim = 0; idim < 3; idim++)
		{
			double* xi = spheres->get_tree_point(isphere);
			double* xj = spheres->get_tree_point(isphere_neighbor);
			e1[idim] = xj[idim] - xi[idim];
			norm += e1[idim] * e1[idim];
		}
		norm = sqrt(norm);
		for (size_t idim = 0; idim < 3; idim++) e1[idim] /= norm;

		for (size_t idim = 0; idim < 3; idim++)
		{
			double* xi = spheres->get_tree_point(isphere);
			double* xj = spheres->get_tree_point(jsphere[0]);
			e2[idim] = xj[idim] - xi[idim];
		}
		double dot(0.0);
		for (size_t idim = 0; idim < 3; idim++) dot += e2[idim] * e1[idim];
		for (size_t idim = 0; idim < 3; idim++) e2[idim] -= dot * e1[idim];

		norm = 0.0;
		for (size_t idim = 0; idim < 3; idim++) norm += e2[idim] * e2[idim];
		norm = sqrt(norm);
		for (size_t idim = 0; idim < 3; idim++) e2[idim] /= norm;

		e3[0] = e1[1] * e2[2] - e1[2] * e2[1];
		e3[1] = e1[2] * e2[0] - e1[0] * e2[2];
		e3[2] = e1[0] * e2[1] - e1[1] * e2[0];

		// sort loops
		double* theta = new double[num_loops];
		for (size_t i = 0; i < num_loops; i++)
		{
			double* xi = spheres->get_tree_point(isphere);
			double* xj = spheres->get_tree_point(jsphere[i]);

			double dot_2(0.0), dot_3(0.0);
			for (size_t idim = 0; idim < 3; idim++)
			{
				dot_2 += (xj[idim] - xi[idim]) * e2[idim];
				dot_3 += (xj[idim] - xi[idim]) * e3[idim];
			}
			theta[i] = _geom.get_point_angle(dot_2, dot_3);
		}
		_memo.quicksort(theta, jsphere, 0, num_loops - 1);

		for (size_t iloop = 0; iloop < num_loops; iloop++)
		{
			// Mark first triangle in the loop

			for (size_t i = 0; i < num_sphere_seeds; i++)
			{
				#pragma region Mark first triangle:
				size_t seed_1(sphere_seeds[2 + i]);
				size_t* attrib = seeds->get_tree_point_attrib(seed_1);

				if (attrib[5] != 0) continue; // pair is already marked

				size_t seed_2(attrib[1]);
				size_t si(attrib[2]), sj(attrib[3]), sk(attrib[4]);

				while (si != isphere)
				{
					size_t tmp = si; si = sj; sj = sk; sk = tmp; // rotate indices
				}

				if (sj != isphere_neighbor && sj != jsphere[iloop]) continue; // wrong triangle
				if (sk != isphere_neighbor && sk != jsphere[iloop]) continue; // wrong triangle


				if (sj != isphere_neighbor)
				{
					size_t tmp = sj; sj = sk; sk = tmp; // flip triangle
				}

				double* x_1 = seeds->get_tree_point(seed_1);
				double* x_2 = seeds->get_tree_point(seed_2);

				face_corners[0] = spheres->get_tree_point(si);
				face_corners[1] = spheres->get_tree_point(sj);
				face_corners[2] = spheres->get_tree_point(sk);

				_geom.get_3d_triangle_normal(face_corners, face_normal);

				double dot(0.0);
				for (size_t idim = 0; idim < 3; idim++) dot += (x_2[idim] - x_1[idim]) * face_normal[idim];

				size_t im(iloop + 1);
				size_t ip(im + 1); if (ip > num_loops) ip = 1;

				if (dot > 0.0)
				{
					attrib[5] = im; // inside
					attrib = seeds->get_tree_point_attrib(seed_2);
					attrib[5] = ip; // outside
				}
				else
				{
					attrib[5] = ip; // outside
					attrib = seeds->get_tree_point_attrib(seed_2);
					attrib[5] = im; // inside
				}
				jsphere[iloop] = sk;
				break;
				#pragma endregion
			}


			while (true)
			{
				#pragma region Mark loop seeds:
				bool done(true);
				for (size_t i = 0; i < num_sphere_seeds; i++)
				{
					size_t seed_1(sphere_seeds[2 + i]);
					size_t* attrib = seeds->get_tree_point_attrib(seed_1);

					if (attrib[5] != 0) continue; // pair is already marked

					size_t seed_2(attrib[1]);
					size_t si(attrib[2]), sj(attrib[3]), sk(attrib[4]);

					while (si != isphere)
					{
						size_t tmp = si; si = sj; sj = sk; sk = tmp; // rotate indices
					}

					if (sj != jsphere[iloop] && sk != jsphere[iloop]) continue; // wrong triangle

					done = false;
				}

				if (done) break; // all pairs are marked

				for (size_t i = 0; i < num_sphere_seeds; i++)
				{
					size_t seed_1(sphere_seeds[2 + i]);
					size_t* attrib = seeds->get_tree_point_attrib(seed_1);

					if (attrib[5] != 0) continue; // pair is already marked

					size_t seed_2(attrib[1]);
					size_t si(attrib[2]), sj(attrib[3]), sk(attrib[4]);

					while (si != isphere)
					{
						size_t tmp = si; si = sj; sj = sk; sk = tmp; // rotate indices
					}

					if (sj != jsphere[iloop] && sk != jsphere[iloop]) continue; // wrong triangle

					if (sj != jsphere[iloop])
					{
						size_t tmp = sj; sj = sk; sk = tmp; // flip triangle
					}

					double* x_1 = seeds->get_tree_point(seed_1);
					double* x_2 = seeds->get_tree_point(seed_2);

					face_corners[0] = spheres->get_tree_point(si);
					face_corners[1] = spheres->get_tree_point(sj);
					face_corners[2] = spheres->get_tree_point(sk);

					_geom.get_3d_triangle_normal(face_corners, face_normal);

					double dot(0.0);
					for (size_t idim = 0; idim < 3; idim++) dot += (x_2[idim] - x_1[idim]) * face_normal[idim];

					size_t im(iloop + 1);
					size_t ip(im + 1); if (ip > num_loops) ip = 1;

					if (dot > 0.0)
					{
						attrib[5] = im; // inside
						attrib = seeds->get_tree_point_attrib(seed_2);
						attrib[5] = ip; // outside
					}
					else
					{
						attrib[5] = ip; // outside
						attrib = seeds->get_tree_point_attrib(seed_2);
						attrib[5] = im; // inside
					}
					jsphere[iloop] = sk;
					if (jsphere[iloop] >= num_surface_spheres) done = true;
					break;
				}
				if (done) break;
				#pragma endregion
			}
		}

		for (size_t i = 0; i < num_sphere_seeds; i++)
		{
			#pragma region connect Sphere seeds:
			size_t seed_i(sphere_seeds[2 + i]);
			size_t* attrib_i = seeds->get_tree_point_attrib(seed_i);
			for (size_t j = i + 1; j < num_sphere_seeds; j++)
			{
				size_t seed_j(sphere_seeds[2 + j]);
				size_t* attrib_j = seeds->get_tree_point_attrib(seed_j);

				if (attrib_i[5] == attrib_j[5]) seeds->graph_connect_nodes(seed_i, seed_j);
			}
			#pragma endregion
		}

		for (size_t i = 0; i < num_sphere_seeds; i++)
		{
			size_t seed_i(sphere_seeds[2 + i]);
			size_t* attrib = seeds->get_tree_point_attrib(seed_i);
			attrib[5] = 0;
		}

		delete[] theta; delete[] jsphere;

		#pragma endregion
	}

	size_t region_id(0);
	while (true)
	{
		#pragma region Coloring Disjoint Subgraphs:
		bool done = true;
		for (size_t iseed = 0; iseed < num_seeds; iseed++)
		{
			size_t* attrib = seeds->get_tree_point_attrib(iseed);
			if (attrib[5] != 0) continue;
			region_id++; attrib[5] = region_id;
			done = false;
			break;
		}
		if (done) break;

		while (true)
		{
			#pragma region Flooding Subregion:
			bool subregion_done(true);
			for (size_t iseed = 0; iseed < num_seeds; iseed++)
			{
				// Establish Sphere -> Seeds directional graph
				size_t* attrib = seeds->get_tree_point_attrib(iseed);
				if (attrib[5] != region_id) continue;

				size_t* neighbor_seeds(0);
				seeds->graph_get_neighbors(iseed, neighbor_seeds);
				if (neighbor_seeds == 0) continue;

				size_t num_neighbors(neighbor_seeds[1]);
				for (size_t i = 0; i < num_neighbors; i++)
				{
					size_t neighbor_seed = neighbor_seeds[2 + i];

					size_t* neighbor_attrib = seeds->get_tree_point_attrib(neighbor_seed);

					if (neighbor_attrib[5] != 0)
					{
						if (neighbor_attrib[5] != region_id)
						{
							vcm_cout << "Error::Invalid coloring!!!" << std::endl;
						}
						continue;
					}
					neighbor_attrib[5] = region_id;
					subregion_done = false;
				}
			}
			if (subregion_done) break;
			#pragma endregion
		}
		#pragma endregion
	}

	num_subregions = region_id - 1;

	if (true)
	{
		#pragma region marking ghost seeds:
		double* xmin = new double[3];
		double* xmax = new double[3];
		for (size_t idim = 0; idim < 3; idim++) xmin[idim] = DBL_MAX;
		for (size_t idim = 0; idim < 3; idim++) xmax[idim] = -DBL_MAX;

		for (size_t iseed = 0; iseed < num_seeds; iseed++)
		{
			double* xseed = seeds->get_tree_point(iseed);
			for (size_t idim = 0; idim < 3; idim++)
			{
				if (xseed[idim] < xmin[idim]) xmin[idim] = xseed[idim];
				if (xseed[idim] > xmax[idim]) xmax[idim] = xseed[idim];
			}
		}
		for (size_t idim = 0; idim < 3; idim++) xmax[idim] += 2.0* (xmax[idim] - xmin[idim]);

		size_t iclosest(0); double hclosest(DBL_MAX);
		seeds->get_closest_tree_point(xmax, iclosest, hclosest);

		size_t* attrib = seeds->get_tree_point_attrib(iclosest);
		size_t ghost_region_id = attrib[5];

		for (size_t iseed = 0; iseed < num_seeds; iseed++)
		{
			size_t* attrib = seeds->get_tree_point_attrib(iseed);
			if (attrib[5] == ghost_region_id) attrib[5] = 0;
		}
		for (size_t iseed = 0; iseed < num_seeds; iseed++)
		{
			size_t* attrib = seeds->get_tree_point_attrib(iseed);
			if (attrib[5] > ghost_region_id) attrib[5]--;
		}
		delete[] xmin; delete[] xmax;
		#pragma endregion
	}

	delete[] e1; delete[] e2; delete[] e3;
	delete[] face_corners; delete[] face_normal;

	end_time = clock();
	end_time = clock();
	cpu_time = ((double)(end_time - start_time)) / CLOCKS_PER_SEC;

	vcm_cout << "  * Number of subregions is " << num_subregions << std::endl;
	vcm_cout << "  * executed in " << cpu_time << " seconds!" << std::endl;
	return 0;
	#pragma endregion
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Private Methods:
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int MeshingVoroCrustSampler::clear_memory()
{
	#pragma region Clear Memory:

	for (size_t iface = 0; iface < _num_faces; iface++) delete[] _face_normal[iface];
	delete[] _face_normal;

#if defined USE_OPEN_MP
	omp_set_num_threads(int(_num_threads));
#pragma omp parallel
	{
		int thread_id = omp_get_thread_num();
#else
	for (size_t thread_id = 0; thread_id < _num_threads; thread_id++)
	{
#endif
		delete _rsamplers[thread_id];
	}
	_rsamplers.clear();
	return 0;
	#pragma endregion
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int MeshingVoroCrustSampler::initiate_active_pool_edges(size_t num_points, double** points, size_t num_edges, size_t** edges)
{
	#pragma region Reset Edges Active Pool:
	_active_pool_size = num_edges;
	_active_pool_parent_face = new size_t[_active_pool_size];
	_active_pool_children = new bool*[_active_pool_size];
	_active_pool_cdf = new double[_active_pool_size];
	for (size_t ied = 0; ied < _active_pool_size; ied++)
	{
		_active_pool_parent_face[ied] = ied;
		_active_pool_children[ied] = 0;
		size_t ip = edges[ied][0];
		size_t iq = edges[ied][1];
		double* p = points[ip];
		double* q = points[iq];
		_active_pool_cdf[ied] = _geom.distance(3, p, q);
	}
	for (size_t ied = 1; ied < _active_pool_size; ied++) _active_pool_cdf[ied] += _active_pool_cdf[ied - 1];
	for (size_t ied = 0; ied < _active_pool_size; ied++) _active_pool_cdf[ied] /= _active_pool_cdf[_active_pool_size - 1];
	_active_pool_ref_level = 0;
	_active_pool_ref_level = 0;
	return 0;
	#pragma endregion
}

int MeshingVoroCrustSampler::refine_active_pool_edges(size_t num_points, double** points, size_t num_edges, size_t** edges,
	                                           MeshingSmartTree* edge_spheres, MeshingSmartTree* corner_spheres, double Lip, double coverage_radius_ratio)
{
	#pragma region Refine Active Pool Of Edges:

	size_t     active_pool_capacity(100); size_t     active_pool_size(0);
	size_t*    active_pool_parent_face = new size_t[active_pool_capacity];
	bool**     active_pool_children = new bool*[active_pool_capacity];
	double*    active_pool_cdf = new double[active_pool_capacity];
	double** edge_corners = new double*[3];
	for (size_t icorner = 0; icorner < 3; icorner++) edge_corners[icorner] = new double[3];

	for (size_t active_edge_index = 0; active_edge_index < _active_pool_size; active_edge_index++)
	{
		get_active_edge_corners(num_points, points, num_edges, edges, active_edge_index, edge_corners);
		bool covered = _smethods.edge_covered(edge_corners, corner_spheres, Lip, 0.0);
		if (!covered) covered = _smethods.edge_covered(edge_corners, edge_spheres, Lip, coverage_radius_ratio);
		if (covered) continue;
		refine_active_edge(active_edge_index, edge_corners, active_pool_size, active_pool_capacity, active_pool_parent_face, active_pool_children, active_pool_cdf);
	}

	for (size_t icorner = 0; icorner < 3; icorner++) delete[] edge_corners[icorner];
	delete[] edge_corners;

	// adjusting cdf
	for (size_t iface = 1; iface < active_pool_size; iface++) active_pool_cdf[iface] += active_pool_cdf[iface - 1];
	for (size_t iface = 0; iface < active_pool_size; iface++) active_pool_cdf[iface] /= active_pool_cdf[active_pool_size - 1];

	// Update active pool
	delete_active_pool();
	if (active_pool_size > 0)
	{
		_active_pool_size = active_pool_size;
		_active_pool_cdf = active_pool_cdf;
		_active_pool_parent_face = active_pool_parent_face;
		_active_pool_children = active_pool_children;
		_active_pool_ref_level++;
	}
	else
	{
		delete[] active_pool_cdf;
		delete[] active_pool_parent_face;
		delete[] active_pool_children;
	}
	return 0;
	#pragma endregion
}

int MeshingVoroCrustSampler::sample_active_edges_uniformly(size_t num_points, double** points, size_t num_edges, size_t** edges, double* dart, size_t &dart_edge)
{
	#pragma region Sample Sharp Edges OF CAD Model Uniformly:
	double u = _rsampler.generate_uniform_random_number();
	size_t ist(0), iend(_active_pool_size - 1);
	while (iend > ist)
	{
		size_t imid = (ist + iend) / 2;
		if (_active_pool_cdf[imid] > u)
		{
			iend = imid;
		}
		else if (_active_pool_cdf[imid] < u)
		{
			ist = imid + 1;
		}
		else
		{
			ist = imid; iend = ist;
		}
	}

	dart_edge = _active_pool_parent_face[ist];

	double** edge_corners = new double*[3];
	for (size_t i = 0; i < 3; i++) edge_corners[i] = new double[3];
	get_active_edge_corners(num_points, points, num_edges, edges, ist, edge_corners);

	u = _rsampler.generate_uniform_random_number();
	for (size_t idim = 0; idim < 3; idim++) dart[idim] = edge_corners[0][idim] + u * (edge_corners[1][idim] - edge_corners[0][idim]);

	for (size_t i = 0; i < 3; i++) delete[] edge_corners[i];
	delete[] edge_corners;
	return 0;
	#pragma endregion
}

int MeshingVoroCrustSampler::get_active_edge_corners(size_t num_points, double** points, size_t num_edges, size_t** edges, size_t active_edge_index, double** edge_corners)
{
	#pragma region Retrieve Active Face Corners:
	size_t iedge = _active_pool_parent_face[active_edge_index];
	for (size_t i = 0; i < 2; i++)
	{
		size_t corner_index = edges[iedge][i];
		for (size_t idim = 0; idim < 3; idim++) edge_corners[i][idim] = points[corner_index][idim];
	}

	for (size_t idim = 0; idim < 3; idim++)
	{
		edge_corners[2][idim] = 0.5 * (edge_corners[0][idim] + edge_corners[1][idim]);
	}

	for (size_t iref = 0; iref < _active_pool_ref_level; iref++)
	{
		bool bo = _active_pool_children[active_edge_index][iref];
		if (!bo)
		{
			// element 0 2
			for (size_t idim = 0; idim < 3; idim++)
			{
				edge_corners[1][idim] = edge_corners[2][idim];
			}
		}
		else
		{
			// element 2 1
			for (size_t idim = 0; idim < 3; idim++)
			{
				edge_corners[0][idim] = edge_corners[2][idim];
			}
		}

		for (size_t idim = 0; idim < 3; idim++)
		{
			edge_corners[2][idim] = 0.5 * (edge_corners[0][idim] + edge_corners[1][idim]);
		}
	}

	return 0;
	#pragma endregion
}

int MeshingVoroCrustSampler::refine_active_edge(size_t active_edge_index, double** edge_corners,
	                                     size_t &active_pool_size, size_t &active_pool_capacity,
	                                     size_t* &active_pool_parent_face, bool** &active_pool_children, double* &active_pool_cdf)
{
	#pragma region Refine Active Edge:
	// Refine and Add uncovered children to pool
	double** child_corners = new double*[2];
	for (size_t icorner = 0; icorner < 2; icorner++) child_corners[icorner] = new double[3];

	for (size_t ichild = 0; ichild < 2; ichild++)
	{
		#pragma region Get Child Corners:
		for (size_t idim = 0; idim < 3; idim++)
		{
			if (ichild == 0)
			{
				// Child 0 2
				child_corners[0][idim] = edge_corners[0][idim];
				child_corners[1][idim] = edge_corners[2][idim];
			}
			else if (ichild == 1)
			{
				// Child 2 1
				child_corners[0][idim] = edge_corners[2][idim];
				child_corners[1][idim] = edge_corners[1][idim];
			}
		}
		#pragma endregion

		#pragma region Add child to active pool:
		size_t child_level = _active_pool_ref_level + 1;
		bool* child_index = new bool[child_level * 2];
		for (size_t ilev = 0; ilev < _active_pool_ref_level; ilev++)
		{
			child_index[ilev] = _active_pool_children[active_edge_index][ilev];
		}
		if (ichild == 0)
		{
			child_index[_active_pool_ref_level] = false;
		}
		else if (ichild == 1)
		{
			child_index[_active_pool_ref_level] = true;
		}

		size_t num(active_pool_size), cap(active_pool_capacity);
		_memo.add_entry(_active_pool_parent_face[active_edge_index], num, active_pool_parent_face, cap);

		num = active_pool_size;  cap = active_pool_capacity;
		_memo.add_entry(child_index, num, active_pool_children, cap);

		double child_edge_length = _geom.distance(3, child_corners[0], child_corners[1]);
		num = active_pool_size;  cap = active_pool_capacity;
		_memo.add_entry(child_edge_length, num, active_pool_cdf, cap);

		active_pool_size = num; active_pool_capacity = cap;
		#pragma endregion
	}

	for (size_t icorner = 0; icorner < 2; icorner++) delete[] child_corners[icorner];
	delete[] child_corners;

	return 0;
	#pragma endregion
}

int MeshingVoroCrustSampler::initiate_active_pool(size_t num_points, double** points, size_t num_faces, size_t** faces, size_t active_patch, size_t &num_patches)
{
	#pragma region Initiate Active Pool:
	size_t bin_size = num_faces / num_patches;
	if (bin_size == 0)
	{
		num_patches = 1;
		bin_size = num_faces;
	}
	size_t min_index = active_patch * bin_size;
	size_t max_index = min_index + bin_size;
	if (active_patch == num_patches - 1) max_index = num_faces;

	size_t num(max_index - min_index);

	_active_pool_size = num; num = 0;
	_active_pool_parent_face = new size_t[_active_pool_size];
	_active_pool_children = new bool*[_active_pool_size];
	_active_pool_cdf = new double[_active_pool_size];
	for (size_t iface = 0; iface < num_faces; iface++)
	{
		if (iface >= min_index && iface < max_index)
		{
			_active_pool_parent_face[num] = iface;
			_active_pool_children[num] = 0;

			double** corners = new double*[3];
			for (size_t icorner = 0; icorner < 3; icorner++)
			{
				size_t corner_index = faces[iface][1 + icorner];
				corners[icorner] = points[corner_index];
			}
			_active_pool_cdf[num] = _geom.get_3d_triangle_area(corners);
			num++;
			delete[] corners;
		}
	}
	for (size_t iface = 1; iface < _active_pool_size; iface++) _active_pool_cdf[iface] += _active_pool_cdf[iface - 1];
	for (size_t iface = 0; iface < _active_pool_size; iface++) _active_pool_cdf[iface] /= _active_pool_cdf[_active_pool_size - 1];
	_active_pool_ref_level = 0;
	return 0;
	#pragma endregion
}

int MeshingVoroCrustSampler::sample_active_pool(size_t num_samples, size_t num_points, double** points, size_t num_faces, size_t** faces,
	                                     MeshingSmartTree* surface_point_cloud, MeshingSmartTree* edge_point_cloud, double rmin,
	                                     MeshingSmartTree* surface_spheres, MeshingSmartTree* edge_spheres, MeshingSmartTree* corner_spheres,
	                                     double smooth_angle_threshold, double rmax, double Lip, double coverage_radius_ratio,
	                                     double** new_sphere, size_t* new_sphere_face, bool* covered_sample)
{
	#pragma region Samping Active Pool:

#if defined USE_OPEN_MP
	double start_omp(0.0), end_omp(0.0);
	start_omp = omp_get_wtime();
	omp_set_num_threads(int(_num_threads));
	#pragma omp parallel
	{
		int thread_id = omp_get_thread_num();
#else
	clock_t start_time, end_time;
	start_time = clock();
	for (size_t thread_id = 0; thread_id < _num_threads; thread_id++)
	{
#endif
		sample_active_faces_uniformly(thread_id, num_points, points, num_faces, faces, new_sphere[thread_id], new_sphere_face[thread_id]);

		covered_sample[thread_id] = false;
		if (!covered_sample[thread_id]) covered_sample[thread_id] = _smethods.point_covered(new_sphere[thread_id], corner_spheres, Lip, 0.0);
		if (!covered_sample[thread_id]) covered_sample[thread_id] = _smethods.point_covered(new_sphere[thread_id], edge_spheres, Lip, 0.0);
		if (!covered_sample[thread_id]) covered_sample[thread_id] = _smethods.point_covered(new_sphere[thread_id], surface_spheres, Lip, coverage_radius_ratio);
	}

	double min_distance_ratio = sqrt(1.0 - coverage_radius_ratio * coverage_radius_ratio);

	

#if defined USE_OPEN_MP
	double start_omp_sz(0.0), end_omp_sz(0.0);
	start_omp_sz = omp_get_wtime();
	omp_set_num_threads(int(_num_threads));
	#pragma omp parallel
	{
		int thread_id = omp_get_thread_num();
#else
	clock_t start_time_sz, end_time_sz;
	start_time_sz = clock();
	for (size_t thread_id = 0; thread_id < _num_threads; thread_id++)
	{
#endif
		if (!covered_sample[thread_id])
		{
			for (size_t jthread = 0; jthread < thread_id; jthread++)
			{
				if (covered_sample[jthread]) continue;
				covered_sample[thread_id] = _smethods.point_covered(new_sphere[thread_id], new_sphere[jthread], coverage_radius_ratio);
				if (covered_sample[thread_id]) break;
			}

			if (!covered_sample[thread_id])
			{
				size_t sphere_face_index = new_sphere_face[thread_id];

				size_t iclosest(0); double hclosest(DBL_MAX);
				surface_point_cloud->get_closest_non_smooth_tree_point(new_sphere[thread_id], _face_normal[sphere_face_index], smooth_angle_threshold, iclosest, hclosest);
				double surf_sz = fmin(rmax, alpha_sz * hclosest);

				if (_input_sizing_function != 0)
				{
					size_t iclosest; double hclosest(DBL_MAX);
					_input_sizing_function->get_closest_tree_point(new_sphere[thread_id], iclosest, hclosest);
					double* x_sz = _input_sizing_function->get_tree_point(iclosest);
					double r_sz = x_sz[3];
					if (surf_sz > r_sz) surf_sz = r_sz;
				}

				if (edge_point_cloud->get_num_tree_points() > 0)
				{
					hclosest = DBL_MAX;
					edge_point_cloud->get_closest_non_smooth_tree_edge_point(new_sphere[thread_id], rmin, smooth_angle_threshold, iclosest, hclosest);
					surf_sz = fmin(surf_sz, alpha_sz * hclosest);
				}

				if (corner_spheres->get_num_tree_points() > 0)
				{
					hclosest = DBL_MAX;
					corner_spheres->get_closest_non_smooth_tree_edge_point(new_sphere[thread_id], rmin, smooth_angle_threshold, iclosest, hclosest);
					surf_sz = fmin(surf_sz, alpha_sz * hclosest);
				}

				if (num_samples > 0)
				{
					size_t num_nearby_spheres_s(0); size_t* nearby_spheres_s(0);
					_smethods.get_overlapping_spheres(new_sphere[thread_id], surface_spheres, Lip, num_nearby_spheres_s, nearby_spheres_s);
					for (size_t i = 0; i < num_nearby_spheres_s; i++)
					{
						size_t sphere_index = nearby_spheres_s[i];
						double* near_by_sphere = surface_spheres->get_tree_point(sphere_index);
						double dst = _geom.distance(3, new_sphere[thread_id], near_by_sphere);
						double coverage_sz = dst / min_distance_ratio;
						surf_sz = fmin(surf_sz, coverage_sz);            // new sphere should not alpha cover an old sphere center
						double Lip_sz = near_by_sphere[3] + Lip * dst;
						surf_sz = fmin(surf_sz, Lip_sz);                 // respect Lipschitz constant locally
					}
					delete[] nearby_spheres_s;
				}

				surf_sz = fmax(surf_sz, _r_min_surface);

				new_sphere[thread_id][3] = surf_sz;
			}
		}
	}

#if defined USE_OPEN_MP
	end_omp_sz = omp_get_wtime();
	_sizing_time += end_omp_sz - start_omp_sz;
#else
	end_time_sz = clock();
	double cpu_time = ((double)(end_time_sz - start_time_sz)) / CLOCKS_PER_SEC;
	_sizing_time += cpu_time;
#endif	

#if defined USE_OPEN_MP
	end_omp = omp_get_wtime();
	_pool_sampling_time += end_omp - start_omp;
#else
	end_time = clock();
	cpu_time = ((double)(end_time - start_time)) / CLOCKS_PER_SEC;
	_pool_sampling_time += cpu_time;
#endif
	
	return 0;
	#pragma endregion
}

int MeshingVoroCrustSampler::refine_active_pool(size_t num_points, double** points, size_t num_faces, size_t** faces,
	                                     MeshingSmartTree* surface_spheres, MeshingSmartTree* edge_spheres, MeshingSmartTree* corner_spheres, double Lip, double coverage_radius_ratio)
{
	#pragma region Generate Active Pool:

#if defined USE_OPEN_MP
	double start_omp(0.0), end_omp(0.0);
	start_omp = omp_get_wtime();
#else
	clock_t start_time, end_time;
	start_time = clock();
#endif
	
	size_t num_thread_faces = _active_pool_size / _num_threads;
	num_thread_faces++;

	size_t* active_pool_capacity = new size_t[_num_threads];
	size_t* active_pool_size = new size_t[_num_threads];
	size_t** active_pool_parent_face = new size_t*[_num_threads];
	bool*** active_pool_children = new bool**[_num_threads];
	double** active_pool_cdf = new double*[_num_threads];

#if defined USE_OPEN_MP
	omp_set_num_threads(int(_num_threads));
	#pragma omp parallel
	{
		int thread_id = omp_get_thread_num();
#else
	for (size_t thread_id = 0; thread_id < _num_threads; thread_id++)
	{
#endif
		active_pool_capacity[thread_id] = 100;
		active_pool_size[thread_id] = 0;
		active_pool_parent_face[thread_id] = new size_t[active_pool_capacity[thread_id]];
		active_pool_children[thread_id] = new bool*[active_pool_capacity[thread_id]];
		active_pool_cdf[thread_id] = new double[active_pool_capacity[thread_id]];
	
		double** face_corners = new double*[6];
		for (size_t icorner = 0; icorner < 6; icorner++) face_corners[icorner] = new double[3];

		for (size_t iface = 0; iface < num_thread_faces; iface++)
		{
			size_t active_face_index = thread_id * num_thread_faces + iface;
			if (active_face_index >= _active_pool_size) continue;

			get_active_face_corners(num_points, points, num_faces, faces, active_face_index, face_corners);
			bool covered = _smethods.face_covered(face_corners, surface_spheres, Lip, coverage_radius_ratio);
			if (!covered) covered = _smethods.face_covered(face_corners, edge_spheres, Lip, 0.0);
			if (!covered) covered = _smethods.face_covered(face_corners, corner_spheres, Lip, 0.0);
			if (covered) continue;
			refine_active_face(active_face_index, face_corners, active_pool_size[thread_id], active_pool_capacity[thread_id], active_pool_parent_face[thread_id], active_pool_children[thread_id], active_pool_cdf[thread_id]);
		}

		for (size_t icorner = 0; icorner < 6; icorner++) delete[] face_corners[icorner];
		delete[] face_corners;
	}

	// Update active pool
	delete_active_pool();

	for (size_t ithread = 1; ithread < _num_threads; ithread++)  active_pool_size[ithread] += active_pool_size[ithread - 1];
	

	if (active_pool_size[_num_threads - 1] > 0)
	{
		_active_pool_size = active_pool_size[_num_threads - 1];
		_active_pool_cdf = new double[_active_pool_size];
		_active_pool_parent_face = new size_t[_active_pool_size];
		_active_pool_children = new bool*[_active_pool_size];

#if defined USE_OPEN_MP
		omp_set_num_threads(int(_num_threads));
		#pragma omp parallel
		{
			int thread_id = omp_get_thread_num();
#else
		for (size_t thread_id = 0; thread_id < _num_threads; thread_id++)
		{
#endif
			size_t num_thread_faces = active_pool_size[thread_id];
			size_t num_prior_faces = 0;
			if (thread_id > 0) num_prior_faces = active_pool_size[thread_id - 1];
			num_thread_faces -= num_prior_faces;

			for (size_t iface = 0; iface < num_thread_faces; iface++)
			{
				size_t active_face_index = num_prior_faces + iface;
				_active_pool_cdf[active_face_index] = active_pool_cdf[thread_id][iface];
				_active_pool_parent_face[active_face_index] = active_pool_parent_face[thread_id][iface];
				_active_pool_children[active_face_index] = active_pool_children[thread_id][iface];
			}
			delete[] active_pool_cdf[thread_id];
			delete[] active_pool_parent_face[thread_id];
			delete[] active_pool_children[thread_id];
		}
		// adjusting cdf
		for (size_t iface = 1; iface < _active_pool_size; iface++) _active_pool_cdf[iface] += _active_pool_cdf[iface - 1];
		for (size_t iface = 0; iface < _active_pool_size; iface++) _active_pool_cdf[iface] /= _active_pool_cdf[_active_pool_size - 1];

		_active_pool_ref_level++;
	}
	else
	{
#if defined USE_OPEN_MP
		omp_set_num_threads(int(_num_threads));
		#pragma omp parallel
		{
			int thread_id = omp_get_thread_num();
#else
		for (size_t thread_id = 0; thread_id < _num_threads; thread_id++)
		{
#endif
			delete[] active_pool_cdf[thread_id];
			delete[] active_pool_parent_face[thread_id];
			delete[] active_pool_children[thread_id];
		}
	}

	delete[] active_pool_size; delete[] active_pool_capacity;
	delete[] active_pool_cdf; delete[] active_pool_parent_face; delete[] active_pool_children;

#if defined USE_OPEN_MP
	end_omp = omp_get_wtime();
	_pool_refinement_time += end_omp - start_omp;
#else
	end_time = clock();
	double cpu_time = ((double)(end_time - start_time)) / CLOCKS_PER_SEC;
	_pool_refinement_time += cpu_time;
#endif

	return 0;
	#pragma endregion
}


int MeshingVoroCrustSampler::sample_active_faces_uniformly(size_t thread_id, size_t num_points, double** points, size_t num_faces, size_t** faces, double* dart, size_t&dart_face)
{
	#pragma region Sample CAD Model Uniformly:
	double u = _rsamplers[thread_id]->generate_uniform_random_number();
	size_t ist(0), iend(_active_pool_size - 1);
	while (iend > ist)
	{
		size_t imid = (ist + iend) / 2;
		if (_active_pool_cdf[imid] > u)
		{
			iend = imid;
		}
		else if (_active_pool_cdf[imid] < u)
		{
			ist = imid + 1;
		}
		else
		{
			ist = imid; iend = ist;
		}
	}

	dart_face = _active_pool_parent_face[ist];

	double** face_corners = new double* [6];
	for (size_t i = 0; i < 6; i++) face_corners[i] = new double[3];
	get_active_face_corners(num_points, points, num_faces, faces, ist, face_corners);

	_rsamplers[thread_id]->sample_uniformly_from_simplex(dart, 3, 3, face_corners);

	for (size_t i = 0; i < 6; i++) delete[] face_corners[i];
	delete[] face_corners;
	return 0;
	#pragma endregion
}

int MeshingVoroCrustSampler::get_active_face_corners(size_t num_points, double** points, size_t num_faces, size_t** faces, size_t active_face_index, double** face_corners)
{
	#pragma region Retrieve Active Face Corners:
	size_t iface = _active_pool_parent_face[active_face_index];
	for (size_t i = 0; i < 3; i++)
	{
		size_t corner_index = faces[iface][i + 1];
		for (size_t idim = 0; idim < 3; idim++) face_corners[i][idim] = points[corner_index][idim];
	}

	for (size_t idim = 0; idim < 3; idim++)
	{
		face_corners[3][idim] = 0.5 * (face_corners[0][idim] + face_corners[1][idim]);
		face_corners[4][idim] = 0.5 * (face_corners[1][idim] + face_corners[2][idim]);
		face_corners[5][idim] = 0.5 * (face_corners[2][idim] + face_corners[0][idim]);
	}

	for (size_t iref = 0; iref < _active_pool_ref_level; iref++)
	{
		bool bo = _active_pool_children[active_face_index][2 * iref];
		bool b1 = _active_pool_children[active_face_index][2 * iref + 1];
		if (!bo && !b1)
		{
			// element 0 3 5
			for (size_t idim = 0; idim < 3; idim++)
			{
				face_corners[1][idim] = face_corners[3][idim];
				face_corners[2][idim] = face_corners[5][idim];
			}
		}
		else if (!bo && b1)
		{
			// element 3 1 4
			for (size_t idim = 0; idim < 3; idim++)
			{
				face_corners[0][idim] = face_corners[3][idim];
				face_corners[2][idim] = face_corners[4][idim];
			}
		}
		else if (bo && !b1)
		{
			// element 5 4 2
			for (size_t idim = 0; idim < 3; idim++)
			{
				face_corners[0][idim] = face_corners[5][idim];
				face_corners[1][idim] = face_corners[4][idim];
			}
		}
		else
		{
			// element 3 4 5
			for (size_t idim = 0; idim < 3; idim++)
			{
				face_corners[0][idim] = face_corners[3][idim];
				face_corners[1][idim] = face_corners[4][idim];
				face_corners[2][idim] = face_corners[5][idim];
			}
		}
		for (size_t idim = 0; idim < 3; idim++)
		{
			face_corners[3][idim] = 0.5 * (face_corners[0][idim] + face_corners[1][idim]);
			face_corners[4][idim] = 0.5 * (face_corners[1][idim] + face_corners[2][idim]);
			face_corners[5][idim] = 0.5 * (face_corners[2][idim] + face_corners[0][idim]);
		}
	}
	return 0;
	#pragma endregion
}

int MeshingVoroCrustSampler::refine_active_face(size_t active_face_index, double** face_corners,
	                                     size_t &active_pool_size, size_t &active_pool_capacity,
	                                     size_t* &active_pool_parent_face, bool** &active_pool_children, double* &active_pool_cdf)
{
	#pragma region Refine Active Face:
	double** child_corners = new double*[3];
	for (size_t icorner = 0; icorner < 3; icorner++) child_corners[icorner] = new double[3];

	for (size_t ichild = 0; ichild < 4; ichild++)
	{
		#pragma region Get Child Corners:
		for (size_t idim = 0; idim < 3; idim++)
		{
			if (ichild == 0)
			{
				// Child 0 3 5
				child_corners[0][idim] = face_corners[0][idim];
				child_corners[1][idim] = face_corners[3][idim];
				child_corners[2][idim] = face_corners[5][idim];
			}
			else if (ichild == 1)
			{
				// Child 3 1 4
				child_corners[0][idim] = face_corners[3][idim];
				child_corners[1][idim] = face_corners[1][idim];
				child_corners[2][idim] = face_corners[4][idim];
			}
			else if (ichild == 2)
			{
				// Child 5 4 2
				child_corners[0][idim] = face_corners[5][idim];
				child_corners[1][idim] = face_corners[4][idim];
				child_corners[2][idim] = face_corners[2][idim];
			}
			else if (ichild == 3)
			{
				// Child 3 4 5
				child_corners[0][idim] = face_corners[3][idim];
				child_corners[1][idim] = face_corners[4][idim];
				child_corners[2][idim] = face_corners[5][idim];
			}
		}
		#pragma endregion

		#pragma region Add child to active pool:
		size_t child_level = _active_pool_ref_level + 1;
		bool* child_index = new bool[child_level * 2];
		for (size_t ilev = 0; ilev < _active_pool_ref_level; ilev++)
		{
			child_index[ilev * 2] = _active_pool_children[active_face_index][ilev * 2];
			child_index[ilev * 2 + 1] = _active_pool_children[active_face_index][ilev * 2 + 1];
		}
		if (ichild == 0)
		{
			child_index[_active_pool_ref_level * 2] = false;
			child_index[_active_pool_ref_level * 2 + 1] = false;
		}
		else if (ichild == 1)
		{
			child_index[_active_pool_ref_level * 2] = false;
			child_index[_active_pool_ref_level * 2 + 1] = true;
		}
		else if (ichild == 2)
		{
			child_index[_active_pool_ref_level * 2] = true;
			child_index[_active_pool_ref_level * 2 + 1] = false;
		}
		else if (ichild == 3)
		{
			child_index[_active_pool_ref_level * 2] = true;
			child_index[_active_pool_ref_level * 2 + 1] = true;
		}

		size_t num(active_pool_size), cap(active_pool_capacity);
		_memo.add_entry(_active_pool_parent_face[active_face_index], num, active_pool_parent_face, cap);

		num = active_pool_size;  cap = active_pool_capacity;
		_memo.add_entry(child_index, num, active_pool_children, cap);

		double child_face_area = _geom.get_3d_triangle_area(child_corners);
		num = active_pool_size;  cap = active_pool_capacity;
		_memo.add_entry(child_face_area, num, active_pool_cdf, cap);

		active_pool_size = num; active_pool_capacity = cap;

		#pragma endregion
	}

	for (size_t icorner = 0; icorner < 3; icorner++) delete[] child_corners[icorner];
	delete[] child_corners;

	return 0;
	#pragma endregion
}

int MeshingVoroCrustSampler::delete_active_pool()
{
	#pragma region Delete Active Pool:
	if (_active_pool_size > 0)
	{
		for (size_t i = 0; i < _active_pool_size; i++)
		{
			if (_active_pool_children[i] == 0) continue;
			delete[] _active_pool_children[i];
		}
		delete[] _active_pool_children; delete[] _active_pool_parent_face; delete[] _active_pool_cdf;
		_active_pool_size = 0; _active_pool_children = 0; _active_pool_parent_face = 0; _active_pool_cdf = 0;
	}
	return 0;
	#pragma endregion
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Retrive overlapping spheres with an overlapp witness from the surface
int MeshingVoroCrustSampler::generate_surface_seeds(size_t num_points, double** points, size_t num_faces, size_t** faces,
	                                         MeshingSmartTree* surface_spheres, MeshingSmartTree* edge_spheres, MeshingSmartTree* corner_spheres,
	                                         double Lip, bool &shrunk_s,
	                                         MeshingSmartTree* upper_seeds, MeshingSmartTree* lower_seeds)
{
	#pragma region Generate Surface Seeds:
	upper_seeds->clear_memory(); lower_seeds->clear_memory();
	shrunk_s = false;
	size_t num_spheres_s = surface_spheres->get_num_tree_points();
	size_t num_spheres_e = edge_spheres->get_num_tree_points();
	size_t num_spheres_c = corner_spheres->get_num_tree_points();

	size_t num_sliver_spheres(0), cap_sliver_spheres(100);
	size_t* sliver_spheres = new size_t[cap_sliver_spheres];

	size_t num_sliver_spheres_radii(0), cap_sliver_spheres_radii(100);
	double* sliver_spheres_radii = new double[cap_sliver_spheres_radii];

	double* sphere = new double[4];

	double* sphere_i = new double[4];
	double* sphere_j = new double[4];
	double* sphere_k = new double[4];

	double* c_ijk = new double[3];
	double** centers = new double*[3];
	double* radii = new double[3];
	double* triplet_normal = new double[3];
	double* upper_seed = new double[4];
	double* lower_seed = new double[4];

	size_t* attrib = new size_t[6];
	attrib[0] = 6; attrib[5] = 0; // Seed Marker

	for (size_t isphere = 0; isphere < num_spheres_s; isphere++)
	{
		#pragma region Form intersection Pairs:
		size_t sphere_index_i(isphere);
		surface_spheres->get_tree_point(isphere, 4, sphere_i);

		size_t sphere_i_face;
		surface_spheres->get_tree_point_attrib(isphere, 0, sphere_i_face);

		size_t num_near_by_spheres(0); size_t* near_by_spheres(0);
		if (true)
		{
			#pragma region Collect spheres overlapping with isphere:

			size_t num_nearby_spheres_s(0); size_t* nearby_spheres_s(0);
			_smethods.get_overlapping_spheres(sphere_i, surface_spheres, Lip, num_nearby_spheres_s, nearby_spheres_s);

			size_t num_nearby_spheres_e(0); size_t* nearby_spheres_e(0);
			_smethods.get_overlapping_spheres(sphere_i, edge_spheres, Lip, num_nearby_spheres_e, nearby_spheres_e);

			size_t num_nearby_spheres_c(0); size_t* nearby_spheres_c(0);
			_smethods.get_overlapping_spheres(sphere_i, corner_spheres, Lip, num_nearby_spheres_c, nearby_spheres_c);

			num_near_by_spheres = num_nearby_spheres_s + num_nearby_spheres_e + num_nearby_spheres_c;
			near_by_spheres = new size_t[num_near_by_spheres];
			for (size_t i = 0; i < num_nearby_spheres_s; i++)
				near_by_spheres[i] = nearby_spheres_s[i];
			for (size_t i = 0; i < num_nearby_spheres_e; i++)
				near_by_spheres[num_nearby_spheres_s + i] = nearby_spheres_e[i] + num_spheres_s;
			for (size_t i = 0; i < num_nearby_spheres_c; i++)
				near_by_spheres[num_nearby_spheres_s + num_nearby_spheres_e + i] = nearby_spheres_c[i] + num_spheres_s + num_spheres_e;

			delete[] nearby_spheres_s; delete[] nearby_spheres_e; delete[] nearby_spheres_c;
			#pragma endregion
		}

		if (false)
		{
			#pragma region Plot Overlapping Spheres:
			double** spheres = new double*[num_near_by_spheres];
			for (size_t jsphere = 0; jsphere < num_near_by_spheres; jsphere++)
			{
				spheres[jsphere] = new double[4];
				size_t sphere_index_j = near_by_spheres[jsphere];
				if (sphere_index_j < num_spheres_s)
					surface_spheres->get_tree_point(sphere_index_j, 4, spheres[jsphere]);
				else if (sphere_index_j < num_spheres_s + num_spheres_e)
					edge_spheres->get_tree_point(sphere_index_j - num_spheres_s, 4, spheres[jsphere]);
				else
					corner_spheres->get_tree_point(sphere_index_j - num_spheres_s - num_spheres_e, 4, spheres[jsphere]);
			}

			for (size_t jsphere = 0; jsphere < num_near_by_spheres; jsphere++) delete[] spheres[jsphere];
			delete[] spheres;
			#pragma endregion
		}

		for (size_t jsphere = 0; jsphere < num_near_by_spheres; jsphere++)
		{
			size_t sphere_index_j = near_by_spheres[jsphere];

			if (sphere_index_j <= sphere_index_i) continue;

			if (sphere_index_j < num_spheres_s)
				surface_spheres->get_tree_point(sphere_index_j, 4, sphere_j);
			else if (sphere_index_j < num_spheres_s + num_spheres_e)
				edge_spheres->get_tree_point(sphere_index_j - num_spheres_s, 4, sphere_j);
			else
				corner_spheres->get_tree_point(sphere_index_j - num_spheres_s - num_spheres_e, 4, sphere_j);

			double dst_ij = _geom.distance(3, sphere_i, sphere_j);

			if (dst_ij > sphere_i[3] + sphere_j[3] - 1E-10) continue;

			for (size_t ksphere = 0; ksphere < num_near_by_spheres; ksphere++)
			{
				size_t sphere_index_k = near_by_spheres[ksphere];

				if (sphere_index_k <= sphere_index_i) continue;

				if (sphere_index_k <= sphere_index_j) continue;

				if (sphere_index_k < num_spheres_s)
					surface_spheres->get_tree_point(sphere_index_k, 4, sphere_k);
				else if (sphere_index_k < num_spheres_s + num_spheres_e)
					edge_spheres->get_tree_point(sphere_index_k - num_spheres_s, 4, sphere_k);
				else
					corner_spheres->get_tree_point(sphere_index_k - num_spheres_s - num_spheres_e, 4, sphere_k);

				double dst_ik = _geom.distance(3, sphere_i, sphere_k);
				if (dst_ik > sphere_i[3] + sphere_k[3] - 1E-10) continue;

				double dst_jk = _geom.distance(3, sphere_j, sphere_k);
				if (dst_jk > sphere_j[3] + sphere_k[3] - 1E-10) continue;

				// Three overlapping spheres
				centers[0] = sphere_i; centers[1] = sphere_j; centers[2] = sphere_k;
				radii[0] = sphere_i[3]; radii[1] = sphere_j[3]; radii[2] = sphere_k[3];
				_geom.get_power_vertex(3, 3, centers, radii, c_ijk);

				double hi = _geom.distance(3, c_ijk, sphere_i);
				if (hi > sphere_i[3] - 1E-10) continue;

				double vi = sqrt(sphere_i[3] * sphere_i[3] - hi * hi);

				// Three overlapping spheres with an intersection pair
				double area = _geom.get_3d_triangle_area(centers);

				if (area < 1E-10)
					continue; // spheres are colinear

				_geom.get_3d_triangle_normal(centers, triplet_normal);

				// adjust triplet normal

				size_t fi_i1 = faces[sphere_i_face][1];
				size_t fi_i2 = faces[sphere_i_face][2];
				size_t fi_i3 = faces[sphere_i_face][3];
				double** fi_corners = new double*[3];
				fi_corners[0] = points[fi_i1]; fi_corners[1] = points[fi_i2], fi_corners[2] = points[fi_i3];
				double* fi_normal = new double[3];
				_geom.get_3d_triangle_normal(fi_corners, fi_normal);

				double dot = _geom.dot_product(3, fi_normal, triplet_normal);
				delete[] fi_normal; delete[] fi_corners;
				if (dot < 0.0)
				{
					attrib[2] = sphere_index_i;
					attrib[3] = sphere_index_k;
					attrib[4] = sphere_index_j;
					for (size_t idim = 0; idim < 3; idim++) triplet_normal[idim] = -triplet_normal[idim];
				}
				else
				{
					attrib[2] = sphere_index_i;
					attrib[3] = sphere_index_j;
					attrib[4] = sphere_index_k;
				}

				for (size_t idim = 0; idim < 3; idim++)
				{
					upper_seed[idim] = c_ijk[idim] + vi * triplet_normal[idim];
					lower_seed[idim] = c_ijk[idim] - vi * triplet_normal[idim];
				}

				size_t num(0), cap(10);
				size_t* covering_spheres = new size_t[cap];

				size_t si(sphere_index_i), sj(sphere_index_j), sk(sphere_index_k);
				if (si >= num_spheres_s) si = num_spheres_s;
				if (sj >= num_spheres_s) sj = num_spheres_s;
				if (sk >= num_spheres_s) sk = num_spheres_s;
				bool upper_covered = _mesh.point_covered(upper_seed, surface_spheres, Lip, 0.0, si, sj, sk, num, cap, covering_spheres);
				bool lower_covered = _mesh.point_covered(lower_seed, surface_spheres, Lip, 0.0, si, sj, sk, num, cap, covering_spheres);

				si = sphere_index_i; sj = sphere_index_j; sk = sphere_index_k;
				if (si < num_spheres_s || si >= num_spheres_s + num_spheres_e) si = num_spheres_e;
				else si -= num_spheres_s;
				if (sj < num_spheres_s || sj >= num_spheres_s + num_spheres_e) sj = num_spheres_e;
				else sj -= num_spheres_s;
				if (sk < num_spheres_s || sk >= num_spheres_s + num_spheres_e) sk = num_spheres_e;
				else sk -= num_spheres_s;

				size_t old_num(num);
				if (!upper_covered) upper_covered = _mesh.point_covered(upper_seed, edge_spheres, Lip, 0.0, si, sj, sk, num, cap, covering_spheres);
				if (!lower_covered) lower_covered = _mesh.point_covered(lower_seed, edge_spheres, Lip, 0.0, si, sj, sk, num, cap, covering_spheres);
				for (size_t ii = old_num; ii < num; ii++) covering_spheres[ii] += num_spheres_s;

				si = sphere_index_i; sj = sphere_index_j; sk = sphere_index_k;
				if (si < num_spheres_s + num_spheres_e) si = num_spheres_c;
				else si -= (num_spheres_s + num_spheres_e);
				if (sj < num_spheres_s + num_spheres_e) sj = num_spheres_c;
				else sj -= (num_spheres_s + num_spheres_e);
				if (sk < num_spheres_s + num_spheres_e) sk = num_spheres_c;
				else sk -= (num_spheres_s + num_spheres_e);

				old_num = num;
				if (!upper_covered) upper_covered = _mesh.point_covered(upper_seed, corner_spheres, Lip, 0.0, si, sj, sk, num, cap, covering_spheres);
				if (!lower_covered) lower_covered = _mesh.point_covered(lower_seed, corner_spheres, Lip, 0.0, si, sj, sk, num, cap, covering_spheres);
				for (size_t ii = old_num; ii < num; ii++) covering_spheres[ii] += num_spheres_s + num_spheres_e;

				if (upper_covered != lower_covered)
				{
					#pragma region A sliver:

					_memo.add_entry(sphere_index_i, num, covering_spheres, cap);
					_memo.add_entry(sphere_index_j, num, covering_spheres, cap);
					_memo.add_entry(sphere_index_k, num, covering_spheres, cap);

					double** spheres = new double*[num]; bool* fixed = new bool[num];
					for (size_t ii = 0; ii < num; ii++)
					{
						spheres[ii] = new double[4];
						if (covering_spheres[ii] < num_spheres_s)
						{
							surface_spheres->get_tree_point(covering_spheres[ii], 4, spheres[ii]);
							fixed[ii] = false;
						}
						else if (covering_spheres[ii] < num_spheres_s + num_spheres_e)
						{
							edge_spheres->get_tree_point(covering_spheres[ii] - num_spheres_s, 4, spheres[ii]);
							fixed[ii] = true;
						}
						else
						{
							corner_spheres->get_tree_point(covering_spheres[ii] - num_spheres_s - num_spheres_e, 4, spheres[ii]);
							fixed[ii] = true;
						}
					}

					size_t sphere_index; double sphere_radius;
					_smethods.resolve_sliver(num, spheres, fixed, sphere_index, sphere_radius);

					if (sphere_index < num)
					{
						sphere_index = covering_spheres[sphere_index];
						bool found(false);
						for (size_t jj = 0; jj < num_sliver_spheres; jj++)
						{
							if (sliver_spheres[jj] != sphere_index) continue;

							if (sliver_spheres_radii[jj] > sphere_radius) sliver_spheres_radii[jj] = sphere_radius;
							found = true;
							break;
						}
						if (!found)
						{
							_memo.add_entry(sphere_index, num_sliver_spheres, sliver_spheres, cap_sliver_spheres);
							_memo.add_entry(sphere_radius, num_sliver_spheres_radii, sliver_spheres_radii, cap_sliver_spheres_radii);
						}
					}
					for (size_t ii = 0; ii < num; ii++) delete[] spheres[ii];
					delete[] spheres; delete[] fixed;
					#pragma endregion
				}
				else if (num_sliver_spheres == 0 && !upper_covered && !lower_covered)
				{
					#pragma region A valid Intersection Pairs:
					size_t iclosest; double hclosest(DBL_MAX);
					size_t upper_seed_index(upper_seeds->get_num_tree_points());
					size_t lower_seed_index(lower_seeds->get_num_tree_points());

					lower_seeds->get_closest_tree_point(lower_seed, iclosest, hclosest);
					if (hclosest < 1E-10) lower_seed_index = iclosest;

					hclosest = DBL_MAX;
					upper_seeds->get_closest_tree_point(upper_seed, iclosest, hclosest);
					if (hclosest < 1E-10) upper_seed_index = iclosest;

					if (lower_seed_index == lower_seeds->get_num_tree_points())
					{
						lower_seed[3] = fmin(sphere_i[3], sphere_j[3]);
						lower_seed[3] = fmin(lower_seed[3], sphere_k[3]);

						attrib[1] = upper_seed_index;
						lower_seeds->add_tree_point(4, lower_seed, triplet_normal, attrib);
					}

					if (upper_seed_index == upper_seeds->get_num_tree_points())
					{
						upper_seed[3] = fmin(sphere_i[3], sphere_j[3]);
						upper_seed[3] = fmin(upper_seed[3], sphere_k[3]);

						attrib[1] = lower_seed_index;
						upper_seeds->add_tree_point(4, upper_seed, triplet_normal, attrib);
					}
					#pragma endregion
				}

				delete[] covering_spheres;
			}
		}
		delete[] near_by_spheres;
		#pragma endregion
	}

	delete[] sphere;
	delete[] sphere_i; delete[] sphere_j; delete[] sphere_k;
	delete[] c_ijk; delete[] centers; delete[] radii;
	delete[] triplet_normal; delete[] upper_seed; delete[] lower_seed;
	delete[] attrib;

	if (num_sliver_spheres > 0)
	{
		#pragma region Shrink Sliver spheres:
		vcm_cout << "  * " << num_sliver_spheres << " sliver spheres were detected and shrunk!" << std::endl;

		double* sphere = new double[4];
		for (size_t isphere = 0; isphere < num_sliver_spheres; isphere++)
		{
			size_t sphere_index = sliver_spheres[isphere];
			if (sphere_index < num_spheres_s)
			{
				surface_spheres->get_tree_point(sphere_index, 4, sphere);
				surface_spheres->set_tree_point_attrib(sphere_index, 0, sliver_spheres_radii[isphere]);
				shrunk_s = true;
			}
			else
			{
				// only surface sphere may shrink here
				vcm_cout << "Warning:: A feature sphere is trying to shrink to eliminate slivers!!" << std::endl;
			}
		}
		upper_seeds->clear_memory();
		lower_seeds->clear_memory();
		delete[] sphere;
		#pragma endregion
	}
	else
	{
		vcm_cout << "  * No sliver spheres were detected!" << std::endl;
	
		// report min distance between upper seeds, lower seeds, lower and upper seeds

		bool* shrink_sphere = new bool[num_spheres_s];
		for (size_t i = 0; i < num_spheres_s; i++) shrink_sphere[i] = false;

		double min_dst_uu(DBL_MAX), min_dst_ll(DBL_MAX), min_dst_ul(DBL_MAX);
		size_t num_upper_seeds = upper_seeds->get_num_tree_points();
		size_t num_lower_seeds = lower_seeds->get_num_tree_points();
		for (size_t iseed = 0; iseed < num_upper_seeds; iseed++)
		{
			size_t* attrib = upper_seeds->get_tree_point_attrib(iseed);
			size_t si(attrib[2]), sj(attrib[3]), sk(attrib[3]);
			double min_radius(DBL_MAX);
			if (si < num_spheres_s)
			{
				double* sphere = surface_spheres->get_tree_point(si);
				if (sphere[3] < min_radius) min_radius = sphere[3];
			}
			if (sj < num_spheres_s)
			{
				double* sphere = surface_spheres->get_tree_point(sj);
				if (sphere[3] < min_radius) min_radius = sphere[3];
			}
			if (sk < num_spheres_s)
			{
				double* sphere = surface_spheres->get_tree_point(sk);
				if (sphere[3] < min_radius) min_radius = sphere[3];
			}

			size_t iclosest(iseed); double hclosest(DBL_MAX);
			upper_seeds->get_closest_tree_point(iseed, iclosest, hclosest);
			hclosest /= min_radius;

			if (hclosest < min_dst_uu) min_dst_uu = hclosest;

			if (hclosest < 1E-2)
			{
				// shrink surface spheres associated with this pair
				shrunk_s = true;
				if (si < num_spheres_s) shrink_sphere[si] = true;
				if (sj < num_spheres_s) shrink_sphere[sj] = true;
				if (sk < num_spheres_s) shrink_sphere[sk] = true;
			}

			double* seed = upper_seeds->get_tree_point(iseed);
			size_t jclosest = num_lower_seeds; double vclosest = DBL_MAX;
			lower_seeds->get_closest_tree_point(seed, jclosest, vclosest);
			vclosest /= min_radius;
			if (vclosest < min_dst_ul) min_dst_ul = vclosest;
		}

		for (size_t iseed = 0; iseed < num_lower_seeds; iseed++)
		{
			size_t* attrib = lower_seeds->get_tree_point_attrib(iseed);
			size_t si(attrib[2]), sj(attrib[3]), sk(attrib[3]);
			double min_radius(DBL_MAX);
			if (si < num_spheres_s)
			{
				double* sphere = surface_spheres->get_tree_point(si);
				if (sphere[3] < min_radius) min_radius = sphere[3];
			}
			if (sj < num_spheres_s)
			{
				double* sphere = surface_spheres->get_tree_point(sj);
				if (sphere[3] < min_radius) min_radius = sphere[3];
			}
			if (sk < num_spheres_s)
			{
				double* sphere = surface_spheres->get_tree_point(sk);
				if (sphere[3] < min_radius) min_radius = sphere[3];
			}

			size_t iclosest(iseed); double hclosest(DBL_MAX);
			lower_seeds->get_closest_tree_point(iseed, iclosest, hclosest);
			hclosest /= min_radius;
			if (hclosest < min_dst_ll) min_dst_ll = hclosest;

			if (hclosest < 1E-2)
			{
				// shrink surface spheres associated with this pair
				shrunk_s = true;
				if (si < num_spheres_s) shrink_sphere[si] = true;
				if (sj < num_spheres_s) shrink_sphere[sj] = true;
				if (sk < num_spheres_s) shrink_sphere[sk] = true;
			}

			double* seed = lower_seeds->get_tree_point(iseed);
			size_t jclosest = num_lower_seeds; double vclosest = DBL_MAX;
			upper_seeds->get_closest_tree_point(seed, jclosest, vclosest);
			vclosest /= min_radius;
			if (vclosest < min_dst_ul) min_dst_ul = vclosest;
		}
		
		vcm_cout << "  * Min. relative distance between upper seeds = " << min_dst_uu << std::endl;
		vcm_cout << "  * Min. relative distance between lower seeds = " << min_dst_ll << std::endl;
		vcm_cout << "  * Min. relative distance between lower and upper seeds = " << min_dst_ul << std::endl;

		if (shrunk_s)
		{
			vcm_cout << "  * shrinking spheres to eliminate too close seeds!" << std::endl;
			for (size_t isphere = 0; isphere < num_spheres_s; isphere++)
			{
				if (!shrink_sphere[isphere]) continue;
				
				double* sphere = surface_spheres->get_tree_point(isphere);
				sphere[3] *= 0.95;
			}
			upper_seeds->clear_memory();
			lower_seeds->clear_memory();
		}

		delete[] shrink_sphere;

		
	}
	delete[] sliver_spheres;
	delete[] sliver_spheres_radii;
	return 0;
	#pragma endregion
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool MeshingVoroCrustSampler::validate_surface_mesh(MeshingSmartTree* surface_spheres, MeshingSmartTree* edge_spheres, MeshingSmartTree* corner_spheres, size_t num_faces, size_t** faces)
{
	#pragma region Validate VC Surface Mesh:
	size_t ns = surface_spheres->get_num_tree_points();
	size_t ne = edge_spheres->get_num_tree_points();
	size_t nc = corner_spheres->get_num_tree_points();

	size_t num_points = ns + ne + nc;
	double** points = new double*[num_points];
	for (size_t i = 0; i < num_points; i++) points[i] = new double[3];

	bool* shrunk = new bool[num_points];
	for (size_t i = 0; i < num_points; i++) shrunk[i] = false;

	for (size_t i = 0; i < ns; i++) surface_spheres->get_tree_point(i, points[i]);
	for (size_t i = 0; i < ne; i++) edge_spheres->get_tree_point(i, points[ns + i]);
	for (size_t i = 0; i < nc; i++) corner_spheres->get_tree_point(i, points[ns + ne + i]);

	MeshingPolyMesh mesh;

	size_t num_single(0);
	bool water_tight(true);
	if (true)
	{
		#pragma region Check Watertightness:
		size_t * num_point_faces = new size_t[num_points];
		for (size_t i = 0; i < num_points; i++) num_point_faces[i] = 0;

		for (size_t iface = 0; iface < num_faces; iface++)
		{
			for (size_t i = 1; i <= 3; i++)
			{
				size_t ipoint = faces[iface][i];
				num_point_faces[ipoint]++;
			}
		}
		size_t** point_faces = new size_t*[num_points];
		for (size_t ipoint = 0; ipoint < num_points; ipoint++)
		{
			size_t num_faces = num_point_faces[ipoint];
			point_faces[ipoint] = new size_t[num_faces + 1];
			point_faces[ipoint][0] = 0;
		}

		delete[] num_point_faces;

		for (size_t iface = 0; iface < num_faces; iface++)
		{
			for (size_t i = 1; i <= 3; i++)
			{
				size_t ipoint = faces[iface][i];
				size_t num_faces = point_faces[ipoint][0];
				point_faces[ipoint][1 + num_faces] = iface;
				point_faces[ipoint][0]++;
			}
		}

		for (size_t iface = 0; iface < num_faces; iface++)
		{
			for (size_t i = 1; i <= 3; i++)
			{
				size_t ip = i + 1; if (ip == 4) ip = 1;

				size_t i1 = faces[iface][i]; size_t i2 = faces[iface][ip];

				if (mesh.get_num_edge_faces(i1, i2, faces, point_faces) == 1)
				{
					if (i1 < ns)
					{
						double* sphere = surface_spheres->get_tree_point(i1);
						sphere[3] *= 0.95;
						water_tight = false;
						num_single++;
					}
					if (i2 < ns)
					{
						double* sphere = surface_spheres->get_tree_point(i2);
						sphere[3] *= 0.95;
						water_tight = false;
						num_single++;
					}					
				}
				
				if (mesh.get_num_edge_faces(i1, i2, faces, point_faces) > 2 && (i1 < ns || i2 < ns))
				{
					// a surface sphere cannot be part of a non-manifold edge
					if (i1 < ns)
					{
						double* sphere = surface_spheres->get_tree_point(i1);
						sphere[3] *= 0.95;
					}
					if (i2 < ns)
					{
						double* sphere = surface_spheres->get_tree_point(i2);
						sphere[3] *= 0.95;
					}
					water_tight = false; num_single++;
				}
			}
		}
		for (size_t ipoint = 0; ipoint < num_points; ipoint++) delete[] point_faces[ipoint];
		delete[] point_faces;
		#pragma endregion
	}

	if (water_tight)
	{
		vcm_cout << "  * Correctness of VoroCrust surface has been verified!" << std::endl;
		mesh.save_mesh_obj("surface_mesh.obj", num_points, points, num_faces, faces);
		vcm_cout << "  * Surface Mesh is saved in surface_mesh.obj!" << std::endl;
	}
	else
	{
		//mesh.save_mesh_obj("surface_mesh_nonwatertight.obj", num_points, points, num_faces, faces);
		vcm_cout << " *** VoroCrust surface is not topologically equivalent to VoroCrust spheres, Nummber of invalid edges = " << num_single <<  std::endl;
	}

	for (size_t i = 0; i < num_points; i++) delete[] points[i];
	delete[] points; delete[] shrunk;
	return water_tight;
	#pragma endregion
}
