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
// MeshingPolyMesh.cpp                                            Last modified (07/12/2022) //
///////////////////////////////////////////////////////////////////////////////////////////////

#include "MeshingPolyMesh.h"

MeshingPolyMesh::MeshingPolyMesh()
{
    reset_global_variables();
}

MeshingPolyMesh::MeshingPolyMesh(size_t num_points, double** points, size_t num_faces, size_t** faces, double smooth_angle_threshold)
{
	#pragma region Constructor:
	vcm_cout << "VoroCrust::New input mesh has been created:" << std::endl;

    reset_global_variables();

	_smooth_angle_threshold = cos(smooth_angle_threshold * PI / 180.0);

	_num_points = num_points;
	_points = new double*[_num_points];
	for (size_t ipoint = 0; ipoint < _num_points; ipoint++)
	{
		_points[ipoint] = new double[3];
		for (size_t idim = 0; idim < 3; idim++) _points[ipoint][idim] = points[ipoint][idim];
	}

	_num_faces = num_faces;
	_faces = new size_t*[_num_faces];
	for (size_t iface = 0; iface < _num_faces; iface++)
	{
		size_t num = faces[iface][0];
		_faces[iface] = new size_t[num + 1];
		_faces[iface][0] = num;
		for (size_t i = 1; i <= num; i++) _faces[iface][i] = faces[iface][i];
	}

	process_model();

	vcm_cout << "  * Number of Input mesh points = " << _num_points << std::endl;
	vcm_cout << "  * Number of Input mesh faces = " << _num_faces << std::endl;
	vcm_cout << "  * Number of Sharp Corners = " << _num_sharp_corners << std::endl;
	vcm_cout << "  * Number of Sharp Edges = " << _num_sharp_edges << std::endl;
	#pragma endregion
}

MeshingPolyMesh::~MeshingPolyMesh()
{
	clear_memory();
}

int MeshingPolyMesh::reset_global_variables()
{
    #pragma region Reset Global Variables:
    _diag = 0.0; _xmin = 0; _xmax = 0;

    _num_points = 0; _num_faces = 0;

    _num_sharp_corners = 0; _num_sharp_edges = 0;

    _smooth_angle_threshold = 0.0;

    _points = 0;  _faces = 0; _point_sharp_edges = 0; _point_faces = 0;

    _face_normal = 0; _face_neighbors = 0;

    _point_basis = 0;

    _sharp_corners = 0; _sharp_edges = 0; _sharp_edge_neighbors = 0;
    return 0;
    #pragma endregion
}

int MeshingPolyMesh::clear_memory()
{
	#pragma region Clear Memory:
	if (_points != 0)
	{
		for (size_t ipoint = 0; ipoint < _num_points; ipoint++) delete[] _points[ipoint];
		delete[] _points; _points = 0;
	}

	if (_faces != 0)
	{
		for (size_t iface = 0; iface < _num_faces; iface++) delete[] _faces[iface];
		delete[] _faces; _faces = 0;
	}

	if (_xmin != 0)
	{
		delete[] _xmin; _xmin = 0;
	}

	if (_xmax != 0)
	{
		delete[] _xmax; _xmax = 0;
	}

	if (_point_sharp_edges != 0)
	{
		for (size_t ipoint = 0; ipoint < _num_points; ipoint++) delete[] _point_sharp_edges[ipoint];
		delete[] _point_sharp_edges; _point_sharp_edges = 0;
	}

	if (_point_faces != 0)
	{
		for (size_t ipoint = 0; ipoint < _num_points; ipoint++) delete[] _point_faces[ipoint];
		delete[] _point_faces; _point_faces = 0;
	}

	if (_face_normal != 0)
	{
		for (size_t iface = 0; iface < _num_faces; iface++) delete[] _face_normal[iface];
		delete[] _face_normal; _face_normal = 0;
	}

	if (_face_neighbors != 0)
	{
		for (size_t iface = 0; iface < _num_faces; iface++) delete[] _face_neighbors[iface];
		delete[] _face_neighbors; _face_neighbors = 0;
	}

	if (_sharp_corners != 0)
	{
		delete[] _sharp_corners; _sharp_corners = 0;
	}

	if (_sharp_edges != 0)
	{
		for (size_t iedge = 0; iedge < _num_sharp_edges; iedge++) delete[] _sharp_edges[iedge];
		delete[] _sharp_edges; _sharp_edges = 0;
	}

	if (_sharp_edge_neighbors != 0)
	{
		for (size_t iedge = 0; iedge < _num_sharp_edges; iedge++) delete[] _sharp_edge_neighbors[iedge];
		delete[] _sharp_edge_neighbors; _sharp_edge_neighbors = 0;
	}

	if (_point_basis != 0)
	{
		for (size_t ipoint = 0; ipoint < _num_points; ipoint++) delete[] _point_basis[ipoint];
		delete[] _point_basis;
	}

	_plc_corners_point_cloud.clear_memory();
	_plc_edge_point_cloud.clear_memory();
	_plc_surface_point_cloud.clear_memory();

    reset_global_variables();

	return 0;
	#pragma endregion
}

int MeshingPolyMesh::read_input_obj_file(std::string filename, size_t &num_points, double** &points, size_t &num_faces, size_t** &faces)
{
	#pragma region Reading Obj File:
	vcm_cout << "VoroCrust::Reading obj file:" << std::endl;

	//open file
	std::ifstream tmpfile(filename.c_str());

	// count vertices and faces
	num_points = 0; num_faces = 0;
	if (tmpfile.is_open() && tmpfile.good())
	{
		std::string line = "";
		while (getline(tmpfile, line))
		{
			std::istringstream iss(line);
			std::vector<std::string> tokens;
			copy(std::istream_iterator<std::string>(iss), std::istream_iterator<std::string>(),
				std::back_inserter<std::vector<std::string> >(tokens));

			if (tokens.size() == 0) continue;

			if (tokens[0] == "v") num_points++;
			if (tokens[0] == "f") num_faces++;
		}
	}

	points = new double*[num_points];
	faces = new size_t*[num_faces];

	std::ifstream myfile(filename.c_str());

	num_points = 0; num_faces = 0;
	if (myfile.is_open() && myfile.good())
	{
		std::string line = "";
		while (getline(myfile, line))
		{
			std::istringstream iss(line);
			std::vector<std::string> tokens;
			copy(std::istream_iterator<std::string>(iss), std::istream_iterator<std::string>(),
				std::back_inserter<std::vector<std::string> >(tokens));

			if (tokens.size() == 0) continue;

			if (tokens[0] == "v")
			{
				points[num_points] = new double[3];
				for (size_t idim = 0; idim < 3; idim++) points[num_points][idim] = _memo.string_to_double(tokens[1 + idim]);
				num_points++;
			}
			if (tokens[0] == "f")
			{

				for (size_t i = 1; i <= 3; i++)
				{
					std::size_t pos = tokens[i].find("/");
					tokens[i] = tokens[i].substr(0, pos);
				}

				faces[num_faces] = new size_t[4];
				faces[num_faces][0] = 3;
				for (size_t i = 1; i <= 3; i++) faces[num_faces][i] = size_t(_memo.string_to_double(tokens[i]) - 1);

				bool redundant(false);
				if (false)
				{
					// check if redundant facet
					for (size_t jface = 0; jface < num_faces; jface++)
					{
						size_t num_common(0);
						for (size_t ii = 1; ii <= 3; ii++)
						{
							for (size_t jj = 1; jj <= 3; jj++)
							{
								if (faces[num_faces][ii] == faces[jface][jj]) num_common++;
							}
						}
						if (num_common == 3)
						{
							redundant = true;
							break;
						}
					}
				}

				if (redundant)
				{
					delete faces[num_faces];
					continue;
				}
				num_faces++;
			}
		}
	}

	if (false)
	{
		#pragma region Pleace the origin at the lower left corner of the bounding box of the model:
		double* xmin = new double[3];
		double* xmax = new double[3];

		for (size_t idim = 0; idim < 3; idim++)
		{
			xmin[idim] = DBL_MAX;
			xmax[idim] = -DBL_MAX;
		}

		for (size_t i = 0; i < num_points; i++)
		{
			for (size_t idim = 0; idim < 3; idim++)
			{
				if (points[i][idim] < xmin[idim]) xmin[idim] = points[i][idim];
				if (points[i][idim] > xmax[idim]) xmax[idim] = points[i][idim];
			}
		}

		double* dx = new double[3];
		for (size_t idim = 0; idim < 3; idim++)
		{
			dx[idim] = xmax[idim] - xmin[idim];
		}

		// place the origin in the lower left corner of the bounding box of the model
		for (size_t i = 0; i < num_points; i++)
		{
			for (size_t idim = 0; idim < 3; idim++) points[i][idim] -= xmin[idim];
		}
		delete[] xmin;
		delete[] xmax;
		delete[] dx;
		#pragma endregion
	}

	weld_nearby_points(num_points, points, num_faces, faces, 1E-8);
	save_mesh_obj("clean.obj", num_points, points, num_faces, faces);

	vcm_cout << "  * Number of Input mesh points = " << num_points << std::endl;
	vcm_cout << "  * Number of Input mesh faces = " << num_faces << std::endl;

	

	return 0;
	#pragma endregion
}

int MeshingPolyMesh::save_mesh_obj(std::string file_name, size_t num_points, double** points, size_t num_faces, size_t** faces)
{
	#pragma region Save Triangular mesh:

	std::fstream file(file_name.c_str(), std::ios::out);

	// Points
	for (size_t i = 0; i < num_points; i++)
	{
		file << std::setprecision(16) << "v " << points[i][0] << " " << points[i][1] << " " << points[i][2] << std::endl;
	}

	// faces
	for (size_t i = 0; i < num_faces; i++)
	{
		file << "f ";

		for (size_t j = 1; j <= 3; j++) file << " " << faces[i][j] + 1;

		file << std::endl;

	}
	return 0;
	#pragma endregion
}

int MeshingPolyMesh::weld_nearby_points(size_t &num_points, double** &points, size_t &num_faces, size_t** &faces, double hmin)
{
	#pragma region Weld Nearby Points:
	MeshingSmartTree tmptree; size_t num(0);
	size_t* point_map = new size_t[num_points];
	for (size_t i = 0; i < num_points; i++)
	{
		if (i == 0)
		{
			tmptree.add_tree_point(3, points[i], 0, 0);
			point_map[i] = num;
			num++;
		}
		else
		{
			size_t iclosest; double hclosest(DBL_MAX);
			tmptree.get_closest_tree_point(points[i], iclosest, hclosest);
			if (hclosest < hmin)
			{
				point_map[i] = iclosest;
				continue;
			}
			else
			{
				tmptree.add_tree_point(3, points[i], 0, 0);
				point_map[i] = num;
				num++;
			}
		}
	}

	if (num == num_points)
	{
		delete[] point_map;
		return 0;
	}

	double** new_points = new double*[num];
	for (size_t i = 0; i < num; i++)
	{
		double* x = tmptree.get_tree_point(i);
		new_points[i] = new double[3];
		for (size_t idim = 0; idim < 3; idim++) new_points[i][idim] = x[idim];
	}

	for (size_t iface = 0; iface < num_faces; iface++)
	{
		size_t num_corners = faces[iface][0];
		for (size_t i = 1; i <= num_corners; i++)
		{
			size_t corner_index = faces[iface][i];
			faces[iface][i] = point_map[corner_index];
		}
	}

	if (true)
	{
		#pragma region Remove degenerate faces:
		size_t iface(0);
		while (iface < num_faces)
		{
			size_t num_corners = faces[iface][0];
			bool degenerate(false);
			for (size_t i = 1; i <= num_corners; i++)
			{
				size_t ip = i + 1;
				if (ip > num_corners) ip = 1;

				if (faces[iface][i] == faces[iface][ip])
				{
					// A degenerate face;
					degenerate = true;
					delete[] faces[iface];
					num_faces--;
					faces[iface] = faces[num_faces];
					break;
				}
			}
			if (degenerate) continue;
			iface++;
		}
		#pragma endregion
	}

	if (false)
	{
		#pragma region Remove redundant faces:
		size_t iface(0);
		while (iface < num_faces)
		{
			bool redundant(false);
			for (size_t jface = iface + 1; jface < num_faces; jface++)
			{
				if (faces[iface][0] != faces[jface][0]) continue;

				bool same(true);
				for (size_t i = 1; i <= faces[iface][0]; i++)
				{
					bool found(false);
					for (size_t j = 1; j <= faces[jface][0]; j++)
					{
						if (faces[iface][i] == faces[jface][j])
						{
							found = true; break;
						}
					}
					if (!found)
					{
						same = false;
						break;
					}
				}
				if (same)
				{
					redundant = true; break;
				}
			}
			if (redundant)
			{
				// A redundant face;
				delete[] faces[iface];
				num_faces--;
				faces[iface] = faces[num_faces];
			}
			else
			{
				iface++;
			}
		}
		#pragma endregion
	}

	for (size_t i = 0; i < num_points; i++) delete[] points[i];
	delete[] points; delete[] point_map;
	tmptree.clear_memory();
	num_points = num; points = new_points;

	return 0;
	#pragma endregion
}

// Retrive overlapping spheres with an overlapp witness from the surface
bool MeshingPolyMesh::generate_surface_mesh(MeshingSmartTree* surface_spheres, MeshingSmartTree* edge_spheres, MeshingSmartTree* corner_spheres,
	                                 double Lip, bool &shrunk_s, bool &shrunk_e, bool &shrunk_c)
{
	#pragma region Generate Surface Seeds:

	MeshingSmartTree upper_seeds; MeshingSmartTree lower_seeds;
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
	size_t* attrib = new size_t[5]; attrib[0] = 5;

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

			write_spheres_csv("nearby_spheres.csv", num_near_by_spheres, spheres);

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

				double* si_normal = surface_spheres->get_tree_point_normal(sphere_index_i);

				// adjust triplet normal
				double dot = _geom.dot_product(3, si_normal, triplet_normal);
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
				bool upper_covered = point_covered(upper_seed, surface_spheres, Lip, 0.0, si, sj, sk, num, cap, covering_spheres);
				bool lower_covered = point_covered(lower_seed, surface_spheres, Lip, 0.0, si, sj, sk, num, cap, covering_spheres);

				si = sphere_index_i; sj = sphere_index_j; sk = sphere_index_k;
				if (si < num_spheres_s || si >= num_spheres_s + num_spheres_e) si = num_spheres_e;
				else si -= num_spheres_s;
				if (sj < num_spheres_s || sj >= num_spheres_s + num_spheres_e) sj = num_spheres_e;
				else sj -= num_spheres_s;
				if (sk < num_spheres_s || sk >= num_spheres_s + num_spheres_e) sk = num_spheres_e;
				else sk -= num_spheres_s;

				size_t old_num(num);
				if (!upper_covered) upper_covered = point_covered(upper_seed, edge_spheres, Lip, 0.0, si, sj, sk, num, cap, covering_spheres);
				if (!lower_covered) lower_covered = point_covered(lower_seed, edge_spheres, Lip, 0.0, si, sj, sk, num, cap, covering_spheres);
				for (size_t ii = old_num; ii < num; ii++) covering_spheres[ii] += num_spheres_s;

				si = sphere_index_i; sj = sphere_index_j; sk = sphere_index_k;
				if (si < num_spheres_s + num_spheres_e) si = num_spheres_c;
				else si -= (num_spheres_s + num_spheres_e);
				if (sj < num_spheres_s + num_spheres_e) sj = num_spheres_c;
				else sj -= (num_spheres_s + num_spheres_e);
				if (sk < num_spheres_s + num_spheres_e) sk = num_spheres_c;
				else sk -= (num_spheres_s + num_spheres_e);

				old_num = num;
				if (!upper_covered) upper_covered = point_covered(upper_seed, corner_spheres, Lip, 0.0, si, sj, sk, num, cap, covering_spheres);
				if (!lower_covered) lower_covered = point_covered(lower_seed, corner_spheres, Lip, 0.0, si, sj, sk, num, cap, covering_spheres);
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
						spheres[ii] = new double[4]; fixed[ii] = false;
						if (covering_spheres[ii] < num_spheres_s)
							surface_spheres->get_tree_point(covering_spheres[ii], 4, spheres[ii]);
						else if (covering_spheres[ii] < num_spheres_s + num_spheres_e)
							edge_spheres->get_tree_point(covering_spheres[ii] - num_spheres_s, 4, spheres[ii]);
						else
							corner_spheres->get_tree_point(covering_spheres[ii] - num_spheres_s - num_spheres_e, 4, spheres[ii]);
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
					size_t upper_seed_index(upper_seeds.get_num_tree_points());
					size_t lower_seed_index(lower_seeds.get_num_tree_points());

					lower_seeds.get_closest_tree_point(lower_seed, iclosest, hclosest);
					if (hclosest < 1E-10) lower_seed_index = iclosest;

					hclosest = DBL_MAX;
					upper_seeds.get_closest_tree_point(upper_seed, iclosest, hclosest);
					if (hclosest < 1E-10) upper_seed_index = iclosest;

					if (lower_seed_index == lower_seeds.get_num_tree_points())
					{
						lower_seed[3] = fmin(sphere_i[3], sphere_j[3]);
						lower_seed[3] = fmin(lower_seed[3], sphere_k[3]);

						attrib[1] = upper_seed_index;
						lower_seeds.add_tree_point(4, lower_seed, triplet_normal, attrib);
					}

					if (upper_seed_index == upper_seeds.get_num_tree_points())
					{
						upper_seed[3] = fmin(sphere_i[3], sphere_j[3]);
						upper_seed[3] = fmin(upper_seed[3], sphere_k[3]);

						attrib[1] = lower_seed_index;
						upper_seeds.add_tree_point(4, upper_seed, triplet_normal, attrib);
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

	double rmin(DBL_MAX);
	if (num_sliver_spheres > 0)
	{
		#pragma region Shrink Sliver spheres:
		vcm_cout << "  * " << num_sliver_spheres << " sliver spheres were detected and shrunk!" << std::endl;

		double* sphere = new double[4];
		for (size_t isphere = 0; isphere < num_sliver_spheres; isphere++)
		{
			size_t sphere_index = sliver_spheres[isphere];
			if (sliver_spheres_radii[isphere] < rmin) rmin = sliver_spheres_radii[isphere];
			if (sphere_index < num_spheres_s)
			{
				surface_spheres->get_tree_point(sphere_index, 4, sphere);
				surface_spheres->set_tree_point_attrib(sphere_index, 0, sliver_spheres_radii[isphere]);
				shrunk_s = true;
			}
			else if (sphere_index < num_spheres_s + num_spheres_e)
			{
				edge_spheres->get_tree_point(sphere_index - num_spheres_s, 4, sphere);
				edge_spheres->set_tree_point_attrib(sphere_index - num_spheres_s, 0, sliver_spheres_radii[isphere]);
				shrunk_e = true;
			}
			else
			{
				corner_spheres->get_tree_point(sphere_index - num_spheres_s - num_spheres_e, 4, sphere);
				corner_spheres->set_tree_point_attrib(sphere_index - num_spheres_s - num_spheres_e, 0, sliver_spheres_radii[isphere]);
				shrunk_c = true;
			}
		}
		vcm_cout << "  * Min radius = " << rmin << std::endl;

		delete[] sphere;
		#pragma endregion
	}
	else
	{
		vcm_cout << "  * No sliver spheres were detected!" << std::endl;
	}
	delete[] sliver_spheres;
	delete[] sliver_spheres_radii;

	if (num_sliver_spheres > 0)
	{
		upper_seeds.clear_memory(); lower_seeds.clear_memory();
		generate_surface_mesh(surface_spheres, edge_spheres, corner_spheres, 0.25, shrunk_s, shrunk_e, shrunk_c);
	}

	size_t num_surf_faces = lower_seeds.get_num_tree_points();
	size_t** surf_faces = new size_t*[num_surf_faces];
	for (size_t iface = 0; iface < num_surf_faces; iface++)
	{
		size_t* attrib = lower_seeds.get_tree_point_attrib(iface);
		surf_faces[iface] = new size_t[4];
		surf_faces[iface][0] = 3;
		surf_faces[iface][1] = attrib[2];
		surf_faces[iface][2] = attrib[3];
		surf_faces[iface][3] = attrib[4];
	}

	valid_surface_mesh(surface_spheres, edge_spheres, corner_spheres, num_surf_faces, surf_faces, shrunk_s, shrunk_e, shrunk_c);

	return 0;
	#pragma endregion
}

bool MeshingPolyMesh::valid_surface_mesh(MeshingSmartTree* surface_spheres, MeshingSmartTree* edge_spheres, MeshingSmartTree* corner_spheres, size_t num_faces, size_t** faces, bool &surface_sphere_shrunk, bool &edge_sphere_shrunk, bool &corner_sphere_shrunk)
{
	#pragma region Plot Active Pool:
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

	size_t num_shrunk(0);
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
					#pragma region Shrink spheres of problematic edge:
					water_tight = false;
					double h = _geom.distance(3, points[i1], points[i2]);
					double rmax = 0.49 * h;

					for (size_t ii = 0; ii < 2; ii++)
					{
						size_t sphere_index(i1);
						if (ii == 1) sphere_index = i2;
						if (sphere_index < ns)
						{
							double r = surface_spheres->get_tree_point_attrib(sphere_index, 0);
							if (r > rmax)
							{
								surface_spheres->set_tree_point_attrib(sphere_index, 0, rmax);
								surface_sphere_shrunk = true;
								if (!shrunk[sphere_index])
								{
									shrunk[sphere_index] = true;
									num_shrunk++;
								}
							}
						}
						else if (sphere_index < ns + ne)
						{
							double r = edge_spheres->get_tree_point_attrib(sphere_index - ns, 0);
							if (r > rmax)
							{
								edge_spheres->set_tree_point_attrib(sphere_index - ns, 0, rmax);
								edge_sphere_shrunk = true;
								if (!shrunk[sphere_index])
								{
									shrunk[sphere_index] = true;
									num_shrunk++;
								}
							}
						}
						else if (sphere_index < ns + ne + nc)
						{
							double r = corner_spheres->get_tree_point_attrib(sphere_index - ns - ne, 0);
							if (r > rmax)
							{
								corner_spheres->set_tree_point_attrib(sphere_index - ns - ne, 0, rmax);
								edge_sphere_shrunk = true;
								if (!shrunk[sphere_index])
								{
									shrunk[sphere_index] = true;
									num_shrunk++;
								}
							}
						}
					}
					#pragma endregion
				}
			}
		}
		for (size_t ipoint = 0; ipoint < num_points; ipoint++) delete[] point_faces[ipoint];
		delete[] point_faces;
		#pragma endregion
	}

	if (true)
	{
		mesh.save_mesh_obj("clean_mesh.obj", num_points, points, num_faces, faces);
		vcm_cout << "  * Surface Mesh is saved in clean_mesh.obj!" << std::endl;
	}
	else
	{
		vcm_cout << "  * VoroCrust surface is not watertight, " << num_shrunk << " spheres have shrunk!" << std::endl;
	}

	for (size_t i = 0; i < num_points; i++) delete[] points[i];
	delete[] points; delete[] shrunk;
	return 0;
	#pragma endregion
}

bool MeshingPolyMesh::point_covered(double* point, MeshingSmartTree* spheres, double Lip, double alpha_coverage, size_t si, size_t sj, size_t sk,
	                         size_t &num_covering_spheres, size_t& cap_covering_spheres, size_t* &covering_spheres)
{
	#pragma region Point Cover Check:
	if (spheres->get_num_tree_points() == 0) return false;

	double beta = sqrt(1.0 - alpha_coverage * alpha_coverage);

	size_t closest_sphere; double closest_dst(DBL_MAX);
	spheres->get_closest_tree_point(point, closest_sphere, closest_dst);

	double closest_sphere_radius = spheres->get_tree_point_attrib(closest_sphere, 0);

	double estimated_face_center_radius = closest_sphere_radius + Lip * closest_dst;

	size_t num(0), cap(100);
	size_t* neighbor_spheres = new size_t[cap];
	double R = estimated_face_center_radius / (1.0 - Lip);
	spheres->get_tree_points_in_sphere(point, R, num, neighbor_spheres, cap);

	bool covered(false);
	double* sphere = new double[4];
	for (size_t i = 0; i < num; i++)
	{
		size_t isphere = neighbor_spheres[i];
		if (isphere == si || isphere == sj || isphere == sk) continue;
		spheres->get_tree_point(isphere, 4, sphere);

		double h = _geom.distance(3, point, sphere);
		if (h < beta * sphere[3] + 1E-10)
		{
			covered = true;
			_memo.add_entry(isphere, num_covering_spheres, covering_spheres, cap_covering_spheres);
		}
	}

	delete[] sphere; delete[] neighbor_spheres;

	return covered;
	#pragma endregion
}


int MeshingPolyMesh::process_model()
{
	#pragma region Process Model:

	build_bounding_box();

	build_point_faces();

	remove_isolated_points();

	if (!water_tight_model())
	{
		//vcm_cout << "  * ... aborting, Good Bye!" << std::endl;
		//clear_memory();
		//abort();
	}

	build_face_normals();

	impose_consistent_orientation();

	extract_model_sharp_features();

	// report min dihedral angle
	double dot_min(1.0);
	for (size_t iface = 0; iface < _num_faces; iface++)
	{
		for (size_t i = 1; i <= 3; i++)
		{
			size_t ip = i + 1; if (ip == 4) ip = 1;

			size_t i1 = _faces[iface][i];
			size_t i2 = _faces[iface][ip];

			if (sharp_edge(i1, i2)) continue;

			size_t* edge_faces;
			get_edge_faces(i1, i2, edge_faces);

			size_t fi = edge_faces[1]; size_t fj = edge_faces[2];
			if (fi != iface)
			{
				fi = edge_faces[2]; fj = edge_faces[1];
			}
			delete[] edge_faces;

			double dot = _geom.dot_product(3, _face_normal[fi], _face_normal[fj]);
			if (dot < dot_min)
			{
				dot_min = dot;
				double ang_min = 180.0 - acos(dot_min) * 180.0 / PI;
			}

		}
	}
	double ang = 180.0 - acos(dot_min) * 180.0 / PI;

	vcm_cout << "  * Min. dihedral angle between smooth neighbors = " << size_t(ang) << " degrees" << std::endl;

	return 0;
	#pragma endregion
}


int MeshingPolyMesh::build_point_sharp_edges()
{
	#pragma region build point sharp edges:
	_point_sharp_edges = new size_t*[_num_points];

	for (size_t ipoint = 0; ipoint < _num_points; ipoint++)
	{
		_point_sharp_edges[ipoint] = new size_t[1];
		_point_sharp_edges[ipoint][0] = 0;
	}

	for (size_t iedge = 0; iedge < _num_sharp_edges; iedge++)
	{
		size_t i1 = _sharp_edges[iedge][0];
		size_t i2 = _sharp_edges[iedge][1];
		_point_sharp_edges[i1][0]++;
		_point_sharp_edges[i2][0]++;
	}

	for (size_t ipoint = 0; ipoint < _num_points; ipoint++)
	{
		size_t num = _point_sharp_edges[ipoint][0];
		delete[] _point_sharp_edges[ipoint];
		_point_sharp_edges[ipoint] = new size_t[1 + num];
		_point_sharp_edges[ipoint][0] = 0;
	}

	for (size_t iedge = 0; iedge < _num_sharp_edges; iedge++)
	{
		size_t i1 = _sharp_edges[iedge][0];
		size_t i2 = _sharp_edges[iedge][1];
		size_t num_1 = _point_sharp_edges[i1][0];
		size_t num_2 = _point_sharp_edges[i2][0];
		_point_sharp_edges[i1][1 + num_1] = iedge;
		_point_sharp_edges[i2][1 + num_2] = iedge;
		_point_sharp_edges[i1][0]++;
		_point_sharp_edges[i2][0]++;
	}
	return 0;
	#pragma endregion
}

int MeshingPolyMesh::build_point_faces()
{
	#pragma region Buiding Mesh data structure:
	// building transpose of connectivity
	size_t* num_point_faces = new size_t[_num_points];
	for (size_t i = 0; i < _num_points; i++) num_point_faces[i] = 0;

	for (size_t iface = 0; iface < _num_faces; iface++)
	{
		for (size_t i = 1; i <= 3; i++)
		{
			size_t ipoint = _faces[iface][i];
			num_point_faces[ipoint]++;
		}
	}
	_point_faces = new size_t*[_num_points];
	for (size_t ipoint = 0; ipoint < _num_points; ipoint++)
	{
		size_t num_faces = num_point_faces[ipoint];
		_point_faces[ipoint] = new size_t[num_faces + 1];
		_point_faces[ipoint][0] = 0;
	}

	delete[] num_point_faces;

	for (size_t iface = 0; iface < _num_faces; iface++)
	{
		for (size_t i = 1; i <= 3; i++)
		{
			size_t ipoint = _faces[iface][i];
			size_t num_faces = _point_faces[ipoint][0];
			_point_faces[ipoint][1 + num_faces] = iface;
			_point_faces[ipoint][0]++;
		}
	}
	return 0;
	#pragma endregion
}

int MeshingPolyMesh::build_point_faces(size_t num_points, size_t num_faces, size_t** faces, size_t** &point_faces)
{
	#pragma region Buiding Mesh data structure:
	// building transpose of connectivity
	size_t* num_point_faces = new size_t[num_points];
	for (size_t i = 0; i < num_points; i++) num_point_faces[i] = 0;

	for (size_t iface = 0; iface < num_faces; iface++)
	{
		for (size_t i = 1; i <= 3; i++)
		{
			size_t ipoint = faces[iface][i];
			num_point_faces[ipoint]++;
		}
	}
	point_faces = new size_t*[num_points];
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
	return 0;
	#pragma endregion
}

int MeshingPolyMesh::build_face_normals()
{
	#pragma region Build Face Normals
	_face_normal = new double*[_num_faces];
	for (size_t iface = 0; iface < _num_faces; iface++)
	{
		double** corners = new double*[3];
		for (size_t i = 0; i < 3; i++)
		{
			size_t ipoint = _faces[iface][1 + i];
			corners[i] = _points[ipoint];
		}
		_face_normal[iface] = new double[3];
		_geom.get_3d_triangle_normal(corners, _face_normal[iface]);
		delete[] corners;
	}
	return 0;
	#pragma endregion
}

int MeshingPolyMesh::build_face_normals(double** points, size_t num_faces, size_t** faces, double** &face_normal)
{
	#pragma region Build Face Normals
	face_normal = new double*[num_faces];
	for (size_t iface = 0; iface < num_faces; iface++)
	{
		double** corners = new double*[3];
		for (size_t i = 0; i < 3; i++)
		{
			size_t ipoint = faces[iface][1 + i];
			corners[i] = points[ipoint];
		}
		face_normal[iface] = new double[3];
		_geom.get_3d_triangle_normal(corners, face_normal[iface]);
		delete[] corners;
	}
	return 0;
	#pragma endregion
}

int MeshingPolyMesh::impose_consistent_orientation()
{
	#pragma region Consistent Orientation:

	bool* fixed = new bool[_num_faces];
	for (size_t iface = 1; iface < _num_faces; iface++) fixed[iface] = false;

	while (true)
	{
		bool done = true;
		for (size_t iface = 0; iface < _num_faces; iface++)
		{
			if (fixed[iface]) continue;
			fixed[iface] = true; done = false;
			break;
		}

		if (done) break; // all manifold patches are consistently oriented now

		done = true;

		for (size_t iface = 0; iface < _num_faces; iface++)
		{
			if (!fixed[iface]) continue;

			for (size_t i = 1; i <= 3; i++)
			{
				size_t ip = i + 1; if (ip == 4) ip = 1;
				size_t pi = _faces[iface][i]; size_t pj = _faces[iface][ip];

				//if (sharp_edge(pi, pj)) continue;

				if (get_num_edge_faces(pi, pj) != 2) continue; // a non manifold edge

				size_t* edge_faces;
				get_edge_faces(pi, pj, edge_faces);

				size_t jface = edge_faces[1];
				if (jface == iface) jface = edge_faces[2];
				delete[] edge_faces;

				if (fixed[jface]) continue;

				for (size_t j = 1; j <= 3; j++)
				{
					size_t jp = j + 1; if (jp == 4) jp = 1;
					size_t qi = _faces[jface][j]; size_t qj = _faces[jface][jp];
					if (qi == pi && qj == pj)
					{
						// adjust orientation of jface
						_faces[jface][j] = qj;
						_faces[jface][jp] = qi;
						done = false;
						break;
					}
				}
				fixed[jface] = true;
			}
		}
		if (done) break;
	}
	delete[] fixed;
	adjust_cad_orientation();
	return 0;
	#pragma endregion
}

bool MeshingPolyMesh::water_tight_model()
{
	#pragma region Check Watertightness:

	for (size_t iface = 0; iface < _num_faces; iface++)
	{
		for (size_t i = 1; i <= 3; i++)
		{
			size_t ip = i + 1; if (ip == 4) ip = 1;

			size_t num_edge_faces = get_num_edge_faces(_faces[iface][i], _faces[iface][ip]);
			if (num_edge_faces == 1)
			{
				vcm_cout << "VoroCrust Error: Input mesh is NOT watertight!" << std::endl;
				vcm_cout << "  * Edge ("<< _faces[iface][i] << "," << _faces[iface][ip] << ") has a single face!" << std::endl;
				return false;
			}
		}
	}

	bool no_skinny(true);
	double* edir = new double[3]; double* q = new double[3];
	while (false)
	{
		bool done = true;
		for (size_t iface = 0; iface < _num_faces; iface++)
		{
			for (size_t i = 1; i <= 3; i++)
			{
				size_t im = i - 1; if (im == 0) im = 3;
				size_t ip = i + 1; if (ip == 4) ip = 1;
				size_t j = _faces[iface][i];
				size_t jm = _faces[iface][im];
				size_t jp = _faces[iface][ip];
				for (size_t idim = 0; idim < 3; idim++) edir[idim] = _points[jp][idim] - _points[jm][idim];
				_geom.normalize_vector(3, edir);
				double dst(DBL_MAX);
				_geom.project_to_3d_line(_points[j], _points[jm], edir, q, dst);
				double h = _geom.distance(3, _points[jm], _points[jp]);
				if (dst < 1E-8 * h)
				{
					no_skinny = false;
					vcm_cout << "VoroCrust Error: Input mesh has invalid faces!" << std::endl;
					vcm_cout << "  * Face (" << j << "," << jm << "," << jp << ") is too skinny!" << std::endl;
				}
				if (dst < 1E-3 * h)
				{
					done = false;
					// Move j away from its projection
					for (size_t idim = 0; idim < 3; idim++)
					{
						double dx = _points[j][idim] - q[idim];
						_points[j][idim] += 1E-3 * h;
					}
				}
			}
		}
		if (done) break;
	}
	delete[] edir; delete[] q;
	return no_skinny;
	#pragma endregion
}

int MeshingPolyMesh::extract_model_sharp_features()
{
	#pragma region Extract Model Sharp Features:
	size_t cap(100), num(0);
	size_t** sharp_edges = new size_t*[cap];
	for (size_t iface = 0; iface < _num_faces; iface++)
	{
		#pragma region Sharp Edges:
		for (size_t ic = 1; ic <= 3; ic++)
		{
			size_t icp = ic + 1; if (icp == 4) icp = 1;
			size_t ic_index = _faces[iface][ic];
			size_t icp_index = _faces[iface][icp];

			size_t num_edge_faces = get_num_edge_faces(ic_index, icp_index);

			size_t* edge_faces;
			get_edge_faces(ic_index, icp_index, edge_faces);
			size_t min_face_index(_num_faces);
			for (size_t ii = 1; ii <= num_edge_faces; ii++)
			{
				if (edge_faces[ii] < min_face_index) min_face_index = edge_faces[ii];
			}

			size_t fi(_num_faces); size_t fj(_num_faces);
			if (num_edge_faces == 2)
			{
				fi = edge_faces[1]; fj = edge_faces[2];
			}

			delete[] edge_faces;
			if (min_face_index != iface) continue; // the face with the lowest index only decides on a sharp edge

			if (num_edge_faces != 2)
			{
				// a non-manifold edge
				size_t* sharp_edge = new size_t[2];
				sharp_edge[0] = ic_index; sharp_edge[1] = icp_index;
				_memo.add_entry(sharp_edge, num, sharp_edges, cap);
			}
			else
			{
				double cos_ang = _geom.dot_product(3, _face_normal[fi], _face_normal[fj]);

				bool reversed(false);

				for (size_t jc = 1; jc <= 3; jc++)
				{
					size_t jcp = jc + 1; if (jcp == 4) jcp = 1;
					size_t jc_index = _faces[fj][ic];
					size_t jcp_index = _faces[fj][icp];
					if (jc_index == ic_index && jcp_index == icp_index)
					{
						reversed = true;
					}
				}
				if (reversed) cos_ang = - cos_ang;

				if (cos_ang < _smooth_angle_threshold)
				{
					size_t* sharp_edge = new size_t[2];
					sharp_edge[0] = ic_index; sharp_edge[1] = icp_index;
					_memo.add_entry(sharp_edge, num, sharp_edges, cap);
				}
			}
		}
		#pragma endregion
	}

	if (num > 0)
	{
		_num_sharp_edges = num;
		_sharp_edges = new size_t*[num];
		for (size_t i = 0; i < num; i++) _sharp_edges[i] = sharp_edges[i];
	}
	delete[] sharp_edges;

	size_t** point_edges = new size_t*[_num_points];
	for (size_t ipoint = 0; ipoint < _num_points; ipoint++)
	{
		point_edges[ipoint] = new size_t[1];
		point_edges[ipoint][0] = 0;
	}

	for (size_t iedge = 0; iedge < _num_sharp_edges; iedge++)
	{
		size_t i1 = _sharp_edges[iedge][0];
		size_t i2 = _sharp_edges[iedge][1];
		point_edges[i1][0]++;
		point_edges[i2][0]++;
	}

	for (size_t ipoint = 0; ipoint < _num_points; ipoint++)
	{
		size_t num = point_edges[ipoint][0];
		delete[] point_edges[ipoint];
		point_edges[ipoint] = new size_t[1 + num];
		point_edges[ipoint][0] = 0;
	}

	for (size_t iedge = 0; iedge < _num_sharp_edges; iedge++)
	{
		size_t i1 = _sharp_edges[iedge][0];
		size_t i2 = _sharp_edges[iedge][1];
		size_t num_1 = point_edges[i1][0];
		size_t num_2 = point_edges[i2][0];
		point_edges[i1][1 + num_1] = iedge;
		point_edges[i2][1 + num_2] = iedge;
		point_edges[i1][0]++;
		point_edges[i2][0]++;
	}

	cap = 100; num = 0;
	size_t* sharp_corners = new size_t[cap];
	for (size_t ipoint = 0; ipoint < _num_points; ipoint++)
	{
		#pragma region Sharp Corners:

		size_t num_point_faces(_point_faces[ipoint][0]);

		if (point_edges[ipoint][0] == 1 || point_edges[ipoint][0] > 2)
		{
			_memo.add_entry(ipoint, num, sharp_corners, cap);
			continue;
		}

		if (point_edges[ipoint][0] == 2)
		{
			size_t ied = point_edges[ipoint][1];
			size_t jed = point_edges[ipoint][2];
			size_t i1 = _sharp_edges[ied][0];
			if (i1 == ipoint) i1 = _sharp_edges[ied][1];
			size_t j1 = _sharp_edges[jed][0];
			if (j1 == ipoint) j1 = _sharp_edges[jed][1];
			double* vi = new double[3];
			double* vj = new double[3];
			for (size_t idim = 0; idim < 3; idim++)
			{
				vi[idim] = _points[i1][idim] - _points[ipoint][idim];
				vj[idim] = _points[ipoint][idim] - _points[j1][idim];
			}
			_geom.normalize_vector(3, vi);
			_geom.normalize_vector(3, vj);
			double dot = _geom.dot_product(3, vi, vj);
			delete[] vi; delete[] vj;
			if (dot < _smooth_angle_threshold)
			{
				_memo.add_entry(ipoint, num, sharp_corners, cap);
				continue;
			}
			continue;
		}

		bool sharp_corner(false);
		for (size_t iface = 0; iface < num_point_faces; iface++)
		{
			size_t fi = _point_faces[ipoint][1 + iface];
			for (size_t jface = iface + 1; jface < num_point_faces; jface++)
			{
				if (iface == jface) continue;
				size_t fj = _point_faces[ipoint][1 + jface];
				double cos_ang = _geom.dot_product(3, _face_normal[fi], _face_normal[fj]);
				if (cos_ang >= _smooth_angle_threshold) continue;

				size_t num_common_vertices(0);
				for (size_t i = 1; i <= _faces[fi][0]; i++)
				{
					for (size_t j = 1; j <= _faces[fj][0]; j++)
					{
						if (_faces[fi][i] == _faces[fj][j]) num_common_vertices++;
					}
				}
				if (num_common_vertices == 1)
				{
					sharp_corner = true; break;
				}
			}
			if (sharp_corner) break;
		}
		if (sharp_corner)
		{
			_memo.add_entry(ipoint, num, sharp_corners, cap);
		}
		#pragma endregion
	}

	for (size_t ipoint = 0; ipoint < _num_points; ipoint++) delete[] point_edges[ipoint];
	delete[] point_edges;

	if (num > 0)
	{
		_sharp_corners = new size_t[num];
		for (size_t i = 0; i < num; i++) _sharp_corners[i] = sharp_corners[i];
		_num_sharp_corners = num;
	}
	delete[] sharp_corners;

	if (_num_sharp_corners > 0) _memo.quicksort(_sharp_corners, 0, _num_sharp_corners - 1);

	build_point_sharp_edges();

	return 0;
	#pragma endregion
}

size_t MeshingPolyMesh::get_num_edge_faces(size_t point_i, size_t point_j)
{
	#pragma region get_edge_faces:
	size_t num_edge_faces = 0;
	size_t num = _point_faces[point_i][0];
	for (size_t i = 1; i <= num; i++)
	{
		size_t face_index = _point_faces[point_i][i];

		for (size_t j = 1; j <= _faces[face_index][0]; j++)
		{
			size_t jp = j + 1; if (jp > _faces[face_index][0]) jp = 1;

			size_t point_index = _faces[face_index][j];
			size_t next_index = _faces[face_index][jp];

			if (point_index == point_i && next_index == point_j) num_edge_faces++;

			if (point_index == point_j && next_index == point_i) num_edge_faces++;
		}
	}
	return num_edge_faces;
	#pragma endregion
}

size_t MeshingPolyMesh::get_num_edge_faces(size_t point_i, size_t point_j, size_t** faces, size_t** point_faces)
{
	#pragma region get_edge_faces:
	size_t num_edge_faces = 0;
	size_t num = point_faces[point_i][0];
	for (size_t i = 1; i <= num; i++)
	{
		size_t face_index = point_faces[point_i][i];

		for (size_t j = 1; j <= faces[face_index][0]; j++)
		{
			size_t jp = j + 1; if (jp > faces[face_index][0]) jp = 1;

			size_t point_index = faces[face_index][j];
			size_t next_index = faces[face_index][jp];

			if (point_index == point_i && next_index == point_j) num_edge_faces++;

			if (point_index == point_j && next_index == point_i) num_edge_faces++;
		}
	}
	return num_edge_faces;
	#pragma endregion
}

int MeshingPolyMesh::get_edge_faces(size_t point_i, size_t point_j, size_t* &edge_faces)
{
	#pragma region get_edge_faces:
	size_t num_edge_faces = get_num_edge_faces(point_i, point_j);
	edge_faces = new size_t[num_edge_faces + 1];
	edge_faces[0] = num_edge_faces; num_edge_faces = 0;

	size_t num = _point_faces[point_i][0];
	for (size_t i = 1; i <= num; i++)
	{
		size_t face_index = _point_faces[point_i][i];

		for (size_t j = 1; j <= _faces[face_index][0]; j++)
		{
			size_t jp = j + 1; if (jp > _faces[face_index][0]) jp = 1;

			size_t point_index = _faces[face_index][j];
			size_t next_index = _faces[face_index][jp];

			if ((point_index == point_i && next_index == point_j) || (point_index == point_j && next_index == point_i))
			{
				edge_faces[1 + num_edge_faces] = face_index;
				num_edge_faces++;
				break;
			}
		}
	}
	return 0;
	#pragma endregion
}

int MeshingPolyMesh::get_edge_faces(size_t point_i, size_t point_j, size_t** faces, size_t** point_faces, size_t* &edge_faces)
{
	#pragma region get_edge_faces:
	size_t num_edge_faces = get_num_edge_faces(point_i, point_j, faces, point_faces);
	edge_faces = new size_t[num_edge_faces + 1];
	edge_faces[0] = num_edge_faces; num_edge_faces = 0;

	size_t num = point_faces[point_i][0];
	for (size_t i = 1; i <= num; i++)
	{
		size_t face_index = point_faces[point_i][i];

		for (size_t j = 1; j <= faces[face_index][0]; j++)
		{
			size_t jp = j + 1; if (jp > faces[face_index][0]) jp = 1;

			size_t point_index = faces[face_index][j];
			size_t next_index = faces[face_index][jp];

			if ((point_index == point_i && next_index == point_j) || (point_index == point_j && next_index == point_i))
			{
				edge_faces[1 + num_edge_faces] = face_index;
				num_edge_faces++;
				break;
			}
		}
	}
	return 0;
	#pragma endregion
}

bool MeshingPolyMesh::sharp_edge(size_t pi, size_t pj)
{
	#pragma region Check If Sharp Edge exists between two mesh points:
	if (_num_sharp_edges == 0) return false;
	if (_point_sharp_edges[pi][0] == 0) return false;
	for (size_t ied = 1; ied <= _point_sharp_edges[pi][0]; ied++)
	{
		size_t edge_index = _point_sharp_edges[pi][ied];
		size_t i1 = _sharp_edges[edge_index][0];
		size_t i2 = _sharp_edges[edge_index][1];
		if (i1 == pj || i2 == pj) return true;
	}
	return false;
	#pragma endregion
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////

int MeshingPolyMesh::remove_isolated_points()
{
	#pragma region Remove isolated points:
	size_t num_isolated(0);
	for (size_t ipoint = 0; ipoint < _num_points; ipoint++)
	{
		if (_point_faces[ipoint][0] == 0) num_isolated++;
	}
	if (num_isolated == 0) return 0;

	size_t num_new_points(_num_points - num_isolated);
	double** new_points = new double*[num_new_points];
	size_t* point_new_index = new size_t[_num_points];

	num_new_points = 0;
	for (size_t ipoint = 0; ipoint < _num_points; ipoint++)
	{
		if (_point_faces[ipoint][0] == 0)
		{
			delete[] _points[ipoint];
			continue;
		}
		point_new_index[ipoint] = num_new_points;
		new_points[num_new_points] = _points[ipoint];
		num_new_points++;
	}

	for (size_t iface = 0; iface < _num_faces; iface++)
	{
		size_t num(_faces[iface][0]);
		for (size_t i = 1; i <= num; i++)
		{
			size_t old_index = _faces[iface][i];
			_faces[iface][i] = point_new_index[old_index];
		}
	}

	for (size_t ipoint = 0; ipoint < _num_points; ipoint++) delete[] _point_faces[ipoint];
	delete[] _point_faces;

	delete[] _points; _points = new_points; _num_points = num_new_points;
	delete[] point_new_index;



	build_point_faces();
	return 0;
	#pragma endregion
}

int MeshingPolyMesh::build_bounding_box()
{
	#pragma region Build Bounding Box:
	_xmin = new double[3];
	_xmax = new double[3];
	for (size_t idim = 0; idim < 3; idim++)
	{
		_xmin[idim] = DBL_MAX; _xmax[idim] = -DBL_MAX;
	}
	for (size_t ipoint = 0; ipoint < _num_points; ipoint++)
	{
		for (size_t idim = 0; idim < 3; idim++)
		{
			if (_points[ipoint][idim] < _xmin[idim]) _xmin[idim] = _points[ipoint][idim];
			if (_points[ipoint][idim] > _xmax[idim]) _xmax[idim] = _points[ipoint][idim];
		}
	}
	_diag = 0.0;
	for (size_t idim = 0; idim < 3; idim++)
	{
		double DX = _xmax[idim] - _xmin[idim];
		_diag += DX * DX;
	}
	_diag = sqrt(_diag);
	return 0;
	#pragma endregion
}

int MeshingPolyMesh::adjust_cad_orientation()
{
	#pragma region Adjust Orientation:

	double* xcorner = new double[3];
	for (size_t idim = 0; idim < 3; idim++)
	{
		xcorner[idim] = _xmin[idim] - (_xmax[idim] - _xmin[idim]);
	}

	double DMIN(DBL_MAX); size_t min_cad_point;
	for (size_t ipoint = 0; ipoint < _num_points; ipoint++)
	{
		double h = _geom.distance(3, _points[ipoint], xcorner);
		if (h < DMIN)
		{
			DMIN = h; min_cad_point = ipoint;
		}
	}
	size_t min_cad_face = _point_faces[min_cad_point][1];

	double dot(0.0);
	for (size_t idim = 0; idim < 3; idim++)
	{
		double dx = xcorner[idim] - _points[min_cad_point][idim];
		dot += _face_normal[min_cad_face][idim] * dx;
	}

	if (dot < 0)
	{
		// reverse orientation
		for (size_t iface = 0; iface < _num_faces; iface++)
		{
			size_t tmp = _faces[iface][2];
			_faces[iface][2] = _faces[iface][3];
			_faces[iface][3] = tmp;
			for (size_t idim = 0; idim < 3; idim++) _face_normal[iface][idim] = -_face_normal[iface][idim];
		}
	}
	delete[] xcorner;
	return 0;
	#pragma endregion
}

int MeshingPolyMesh::write_spheres_csv(std::string file_name, size_t num_spheres, double** spheres)
{
	#pragma region Plot Spheres:
	std::fstream file(file_name.c_str(), std::ios::out);
	// Spheres
	file << "x coord,y coord,z coord,radius" << std::endl;
	for (size_t i = 0; i < num_spheres; i++)
	{
		file << std::setprecision(12) << spheres[i][0] << ", " << spheres[i][1] << ", " << spheres[i][2] << ", " << spheres[i][3] << std::endl;
	}
	return 0;
	#pragma endregion
}
