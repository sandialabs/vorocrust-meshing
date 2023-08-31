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
//  MeshingVoroCrust.cpp                                          Last modified (07/12/2022) //
///////////////////////////////////////////////////////////////////////////////////////////////

#include "MeshingVoroCrust.h"


MeshingVoroCrust::MeshingVoroCrust(std::string input_filename)
: input_filename(input_filename)
{
	// generate_mesh_from_seeds_file(1, false, false);
}


MeshingVoroCrust::~MeshingVoroCrust()
{

}


int MeshingVoroCrust::execute()
{
	#pragma region Start Here:
	// open the file
	std::ifstream infile(this->input_filename);
	//std::ifstream infile( options->get_filename() );

	std::vector<std::string> filenames;
	double Lip_const(0.1);
	double vc_ang_tol(20.0);
	double r_min(0.0), r_max(DBL_MAX);
	double feature_TOL(0.0);
	int num_threads(1);

	size_t num_loop_ref(0);
	double ref_ang_tol(20.0);

	bool geneate_vcg_file(false);
    bool generate_exodus_file(false);
	bool generate_monitoring_points(false), impose_monitoring_points(false);

	bool generate_mesh_from_seeds(false);

	std::string sizing_file;

	if (infile.is_open() && infile.good())
	{
		std::string line = "";
		while (getline(infile, line))
		{
			std::istringstream iss(line);
			std::vector<std::string> tokens;
			copy(std::istream_iterator<std::string>(iss), std::istream_iterator<std::string>(),
				std::back_inserter<std::vector<std::string> >(tokens));

			if (tokens.size() == 0) continue;

			if (tokens[0] == "INPUT_MESH_FILE") filenames.push_back(tokens[2]);
			if (tokens[0] == "SIZING_FILE") sizing_file = tokens[2];
			if (tokens[0] == "R_MIN") r_min = _memo.string_to_double(tokens[2]);
			if (tokens[0] == "R_MAX") r_max = _memo.string_to_double(tokens[2]);
			if (tokens[0] == "LIP_CONST") Lip_const = _memo.string_to_double(tokens[2]);
			if (tokens[0] == "VC_ANGLE") vc_ang_tol = _memo.string_to_double(tokens[2]);
			if (tokens[0] == "NUM_THREADS") num_threads = int(_memo.string_to_double(tokens[2]));
			if (tokens[0] == "GENERATE_VCG_FILE") geneate_vcg_file = true;
            if (tokens[0] == "GENERATE_EXODUS_FILE") generate_exodus_file = true;
			if (tokens[0] == "GENERATE_MONITORING_POINTS") generate_monitoring_points = true;
			if (tokens[0] == "IMPOSE_MONITORING_POINTS")   impose_monitoring_points = true;
			if (tokens[0] == "GENERATE_MESH_FROM_SEEDS")   generate_mesh_from_seeds = true;
			if (tokens[0] == "SEEDS_FILE") filenames.push_back(tokens[2]);
		}
	}
	else
	{
		vcm_cout << "\n*** Failed to open input file, good bye!";
		char ii;
		std::cin >> ii;
		return 1;
	}

	if (generate_mesh_from_seeds)
	{
		vcm_cout << "MeshingVoroCrust::Input Data:" << std::endl;
		vcm_cout << "  * Seeds file = " << filenames[0] << std::endl;
		if (geneate_vcg_file)
			vcm_cout << "    Generate VCG file!" << std::endl;
		if (generate_exodus_file)
			vcm_cout << "    Generate Exodus file!" << std::endl;
		vcm_cout << "    Number of OpenMP threads = " << num_threads << std::endl;
		generate_mesh_from_seeds_file(filenames[0], num_threads, geneate_vcg_file, generate_exodus_file);
		vcm_cout << "\n*** MeshingVoroCrust::Mission Accomplished!" << std::endl;
		return 0;
	}

	size_t num_input_files(filenames.size());

	// Echo input parameters
	//vcm_cout << "MeshingVoroCrust::Input Data:" << std::endl;
	vcm_cout << "MeshingVoroCrust::Input Data:" << std::endl;
	for (size_t i = 0; i < num_input_files; i++)
		vcm_cout << "  * Input Mesh = " << filenames[i] << std::endl;
	if (sizing_file != "")
		vcm_cout << "    Sizing File = " << sizing_file << std::endl;

	vcm_cout << "    Minimum Sphere Radius = " << r_min << std::endl;
	vcm_cout << "    Maximum Sphere Radius = " << r_max << std::endl;
	vcm_cout << "    Lipschitz Const = " << Lip_const << std::endl;
	vcm_cout << "    VC Smooth Angle Threshold = " << vc_ang_tol << std::endl;
	if (feature_TOL == 0.0)
		vcm_cout << "    A clean input model " << std::endl;
	else
		vcm_cout << "    Feature Tolernce = " << feature_TOL << std::endl;
	vcm_cout << "    Number of OpenMP threads = " << num_threads << std::endl;
	if (geneate_vcg_file)
		vcm_cout << "    Generate VCG file!" << std::endl;
	if (generate_exodus_file)
		vcm_cout << "    Generate Exodus file!" << std::endl;
	if (generate_monitoring_points)
		vcm_cout << "    Generate Monitoring Points!" << std::endl;
	if (impose_monitoring_points)
		vcm_cout << "    Impose Monitoring Points!" << std::endl;

    //If Exodus is not enabled but the user requested to generate an Exodus file, we generate a warning.
    #ifdef NO_EXODUS
    if(generate_exodus_file){
            vcm_cout << "WARNING: The input file requests that MeshingVoroCrust generates an Exodus file, but Exodus support is disabled." <<
            " No Exodus file will be generated." << std::endl;
    }
    #endif //NO_EXODUS



	MeshingSmartTree* input_sizing_function(0);

	if (sizing_file != "")
	{
		input_sizing_function = new MeshingSmartTree(3);
		input_sizing_function->load_tree_csv(sizing_file, 4);
	}

	size_t num_points(0); double** points(0); size_t num_faces(0); size_t** faces(0);
	if (num_input_files > 1)
	{
		#pragma region Merging input models:
		MeshingPolyMesh tmp_plc;
		for (size_t i = 0; i < num_input_files; i++)
		{
			size_t tmp_num_points; double** tmp_points; size_t tmp_num_faces; size_t** tmp_faces;
			tmp_plc.read_input_obj_file(filenames[i].c_str(), tmp_num_points, tmp_points, tmp_num_faces, tmp_faces);

			for (size_t iface = 0; iface < tmp_num_faces; iface++)
			{
				for (size_t ic = 1; ic <= 3; ic++) tmp_faces[iface][ic] += num_points;
			}

			size_t combined_num_points(num_points + tmp_num_points);
			double** combined_points = new double* [combined_num_points];
			for (size_t ipoint = 0; ipoint < num_points; ipoint++) combined_points[ipoint] = points[ipoint];
			for (size_t ipoint = 0; ipoint < tmp_num_points; ipoint++) combined_points[num_points + ipoint] = tmp_points[ipoint];

			size_t combined_num_faces(num_faces + tmp_num_faces);
			size_t** combined_faces = new size_t * [combined_num_faces];
			for (size_t iface = 0; iface < num_faces; iface++) combined_faces[iface] = faces[iface];
			for (size_t iface = 0; iface < tmp_num_faces; iface++) combined_faces[num_faces + iface] = tmp_faces[iface];

			delete[] tmp_points; delete[] tmp_faces;
			if (i > 0)
			{
				delete[] points; delete[] faces;
			}

			num_points = combined_num_points; num_faces = combined_num_faces;
			points = combined_points; faces = combined_faces;
		}
		tmp_plc.weld_nearby_points(num_points, points, num_faces, faces, 1E-8);
		tmp_plc.save_mesh_obj("model_combined.obj", num_points, points, num_faces, faces);
		#pragma endregion
	}
	else
	{
		MeshingPolyMesh tmp_plc;
		tmp_plc.read_input_obj_file(filenames[0].c_str(), num_points, points, num_faces, faces);
	}

	//clock_t start_time, end_time; double cpu_time;
	MeshingTimer timer;

	execute_vc(num_threads, num_points, points, num_faces, faces, input_sizing_function, r_min, r_max, Lip_const, vc_ang_tol, geneate_vcg_file,
		       generate_exodus_file,generate_monitoring_points, impose_monitoring_points);

	vcm_cout << "\n*** MeshingVoroCrust::Mission Accomplished in " << timer.report_timing() << " seconds! ***" << std::endl;

	for (size_t ipoint = 0; ipoint < num_points; ipoint++) delete[] points[ipoint];
	delete[] points;

	for (size_t iface = 0; iface < num_faces; iface++) delete[] faces[iface];
	delete[] faces;


	return 0;
	#pragma endregion
}


int MeshingVoroCrust::execute_vc(int num_threads, size_t num_points, double** points, size_t num_faces, size_t** faces, MeshingSmartTree* input_sizing_function,
	                             double rmin, double rmax, double Lip, double smooth_angle_threshold, bool generate_vcg_file,bool generate_exodus_file,
	                             bool generate_monitoring_points, bool impose_monitoring_points)
{
	#pragma region Execute MeshingVoroCrust Wild:

	bool smooth_model(false);

	smooth_angle_threshold = cos(smooth_angle_threshold * PI / 180.0);

	MeshingSmartTree* surface_point_cloud = new MeshingSmartTree();
	MeshingSmartTree* edge_point_cloud = new MeshingSmartTree();
	size_t num_sharp_corners(0); size_t* sharp_corners(0);
	size_t num_sharp_edges(0); size_t** sharp_edges(0);

	MeshingSmartTree* surface_spheres = new MeshingSmartTree();
	MeshingSmartTree* edge_spheres = new MeshingSmartTree();
	MeshingSmartTree* corner_spheres = new MeshingSmartTree();
	MeshingSmartTree* surf_ext_seeds = new MeshingSmartTree();
	MeshingSmartTree* surf_int_seeds = new MeshingSmartTree();

	size_t num_subregions(0);
	MeshingSmartTree* interface_spheres = new MeshingSmartTree();
	MeshingSmartTree* seeds = new MeshingSmartTree();

	MeshingVoroCrustSampler* vc_sampler = new MeshingVoroCrustSampler(num_threads);

	vc_sampler->set_input_sizing_function(input_sizing_function);

	if (smooth_model)
		vc_sampler->dont_detect_features_generate_point_clouds(num_points, points, num_faces, faces, surface_point_cloud, edge_point_cloud, num_sharp_edges, sharp_edges, num_sharp_corners, sharp_corners);
	else if (rmin == 0.0) // A clean input
		vc_sampler->detect_features_generate_point_clouds_clean(num_points, points, num_faces, faces, smooth_angle_threshold, surface_point_cloud, edge_point_cloud, num_sharp_edges, sharp_edges, num_sharp_corners, sharp_corners);
	else
		vc_sampler->detect_features_generate_point_clouds(num_points, points, num_faces, faces, smooth_angle_threshold, rmin, surface_point_cloud, edge_point_cloud, num_sharp_edges, sharp_edges, num_sharp_corners, sharp_corners);

	if (num_sharp_corners > 0)
		vc_sampler->generate_corner_spheres(num_points, points, num_sharp_corners, sharp_corners, surface_point_cloud, edge_point_cloud, rmin, corner_spheres, smooth_angle_threshold, rmax, Lip);
	if (num_sharp_edges > 0)
		vc_sampler->generate_edge_spheres(num_points, points, num_sharp_edges, sharp_edges, surface_point_cloud, edge_point_cloud, rmin, edge_spheres, corner_spheres, smooth_angle_threshold, rmax, Lip, 0.5);

	vc_sampler->generate_surface_spheres(num_points, points, num_faces, faces, surface_point_cloud, edge_point_cloud, rmin, surface_spheres, edge_spheres, corner_spheres, surf_ext_seeds, surf_int_seeds, smooth_angle_threshold, rmax, Lip, 0.4);

	vc_sampler->color_surface_seeds(surface_spheres, edge_spheres, corner_spheres, surf_ext_seeds, surf_int_seeds, interface_spheres, seeds, num_subregions);
	delete surf_ext_seeds; delete surf_int_seeds; delete interface_spheres;

	delete surface_point_cloud; delete edge_point_cloud;
	for (size_t iedge = 0; iedge < num_sharp_edges; iedge++) delete[] sharp_edges[iedge];
	if (num_sharp_edges > 0) delete[] sharp_edges;
	if (num_sharp_corners > 0) delete[] sharp_corners;

	delete vc_sampler;

	size_t num_surface_seeds(seeds->get_num_tree_points());

	size_t num_seed_points(0);
	double* seed_points(0); size_t* seed_points_region_id(0); double* seed_points_sizing(0);
	generate_interior_seeds(seeds, surface_spheres, edge_spheres, corner_spheres, input_sizing_function, num_threads, Lip, rmax,
		                    impose_monitoring_points, generate_monitoring_points,
		                    num_seed_points, seed_points, seed_points_region_id, seed_points_sizing);

	if (generate_vcg_file || generate_exodus_file)
	{
		generate_explicit_mesh(num_threads, num_seed_points, seed_points, seed_points_region_id, seed_points_sizing, num_surface_seeds, generate_vcg_file, generate_exodus_file);
	}
	return 0;
	#pragma endregion
}

int MeshingVoroCrust::generate_interior_seeds(MeshingSmartTree* seeds_tree, MeshingSmartTree* surface_spheres_tree, MeshingSmartTree* edge_spheres_tree, MeshingSmartTree* corner_spheres_tree,
	                                          MeshingSmartTree* sz_function_tree, int num_threads, double Lip, double rmax,
	                                          bool impose_monitoring_points, bool generate_monitoring_points,
	                                          size_t &num_seeds, double* &seeds, size_t* &seeds_region_id, double* &seeds_sizing)
{
	#pragma region Generate Interior Seeds:

	num_seeds = seeds_tree->get_num_tree_points();
	size_t num_surface_seeds(num_seeds);
	seeds = new double[num_seeds * 3];
	seeds_region_id = new size_t[num_seeds];
	seeds_sizing = new double[num_seeds];

	for (size_t iseed = 0; iseed < num_seeds; iseed++)
	{
		double* x = seeds_tree->get_tree_point(iseed);
		size_t* attrib = seeds_tree->get_tree_point_attrib(iseed);
		for (size_t idim = 0; idim < 3; idim++) seeds[iseed * 3 + idim] = x[idim];
		seeds_region_id[iseed] = attrib[5];
		seeds_sizing[iseed] = x[3];
	}
	delete seeds_tree;

	size_t num_surface_spheres = surface_spheres_tree->get_num_tree_points();
	double* surface_spheres(0);
	double* surface_spheres_sizing(0);
	if (num_surface_spheres > 0)
	{
		surface_spheres = new double[num_surface_spheres * 3];
		surface_spheres_sizing = new double[num_surface_spheres];
		for (size_t isphere = 0; isphere < num_surface_spheres; isphere++)
		{
			double* x = surface_spheres_tree->get_tree_point(isphere);
			for (size_t idim = 0; idim < 3; idim++) surface_spheres[isphere * 3 + idim] = x[idim];
			surface_spheres_sizing[isphere] = x[3];
		}
	}
	delete surface_spheres_tree;

	size_t num_edge_spheres = edge_spheres_tree->get_num_tree_points();
	double* edge_spheres(0);
	double* edge_spheres_sizing(0);
	if (num_edge_spheres > 0)
	{
		edge_spheres = new double[num_edge_spheres * 3];
		edge_spheres_sizing = new double[num_edge_spheres];
		for (size_t isphere = 0; isphere < num_edge_spheres; isphere++)
		{
			double* x = edge_spheres_tree->get_tree_point(isphere);
			for (size_t idim = 0; idim < 3; idim++) edge_spheres[isphere * 3 + idim] = x[idim];
			edge_spheres_sizing[isphere] = x[3];
		}
	}
	delete edge_spheres_tree;

	size_t num_corner_spheres = corner_spheres_tree->get_num_tree_points();
	double* corner_spheres(0);
	double* corner_spheres_sizing(0);
	if (num_corner_spheres > 0)
	{
		corner_spheres = new double[num_corner_spheres * 3];
		corner_spheres_sizing = new double[num_corner_spheres];
		for (size_t isphere = 0; isphere < num_corner_spheres; isphere++)
		{
			double* x = corner_spheres_tree->get_tree_point(isphere);
			for (size_t idim = 0; idim < 3; idim++) corner_spheres[isphere * 3 + idim] = x[idim];
			corner_spheres_sizing[isphere] = x[3];
		}
	}
	delete corner_spheres_tree;

	size_t num_sz_points(0);
	if (sz_function_tree != 0) num_sz_points = sz_function_tree->get_num_tree_points();
	double* sz_points(0);
	double* sz_value(0);
	if (num_sz_points > 0)
	{
		sz_points = new double[num_sz_points * 3];
		sz_value = new double[num_sz_points];
		for (size_t ipoint = 0; ipoint < num_sz_points; ipoint++)
		{
			double* x = sz_function_tree->get_tree_point(ipoint);
			for (size_t idim = 0; idim < 3; idim++) sz_points[ipoint * 3 + idim] = x[idim];
			sz_value[ipoint] = x[3];
		}
	}
	delete sz_function_tree;

	MeshingVoronoiMesher mesher;

	MeshingTimer timer;
	size_t num_imposed_seeds(0); double* imposed_seeds(0);

	if (true)
	{
		vcm_cout << "MeshingVoroCrust::Ensuring Sharp Feature Spheres are Delaunay:" << std::endl;

		timer.start();

		mesher.ensure_sharp_edge_spheres_are_Delaunay(num_seeds, seeds, seeds_sizing, seeds_region_id,
			                                             num_surface_spheres, surface_spheres, surface_spheres_sizing,
			                                             num_edge_spheres, edge_spheres, edge_spheres_sizing,
			                                             num_corner_spheres, corner_spheres, corner_spheres_sizing,
			                                             num_sz_points, sz_points, sz_value,
			                                             num_threads, Lip, rmax);

		vcm_cout << "  * executed in " << timer.report_timing() << " seconds!" << std::endl;
	}

	size_t num_critical_seeds(num_seeds);

	size_t num_imposed_seeds_added(0);
	if (impose_monitoring_points)
	{
		vcm_cout << "MeshingVoroCrust::Imposing Monitoring Points:" << std::endl;

		timer.start();
		
		mesher.load_spheres_csv("monitoring_points.csv", 3, num_imposed_seeds, imposed_seeds);

		mesher.impose_interior_seeds(num_seeds, seeds, seeds_sizing, seeds_region_id,
			                         num_surface_spheres, surface_spheres, surface_spheres_sizing,
			                         num_edge_spheres, edge_spheres, edge_spheres_sizing,
			                         num_corner_spheres, corner_spheres, corner_spheres_sizing,
			                         num_sz_points, sz_points, sz_value,
			                         num_imposed_seeds, imposed_seeds, num_threads, Lip, rmax);

		generate_monitoring_points = false;

		num_imposed_seeds_added = num_seeds - num_surface_seeds;
		vcm_cout << "  * added " << num_imposed_seeds_added << " / " << num_imposed_seeds << " desired imposed points!" << std::endl;
		vcm_cout << "  * executed in " << timer.report_timing() << " seconds!" << std::endl;
	}


	vcm_cout << "MeshingVoroCrust::Generating Interior Seeds:" << std::endl;

	timer.start();

	mesher.generate_interior_seeds(num_seeds, seeds, seeds_sizing, seeds_region_id,
		                           num_surface_spheres, surface_spheres, surface_spheres_sizing,
		                           num_edge_spheres, edge_spheres, edge_spheres_sizing,
		                           num_corner_spheres, corner_spheres, corner_spheres_sizing,
		                           num_sz_points, sz_points, sz_value, num_threads, Lip, rmax);
	

	if (impose_monitoring_points)
	{
		vcm_cout << "MeshingVoroCrust::Retrieving closest seeds to Monitoring Points:" << std::endl;

		double** seeds_ = new double* [num_seeds];
		for (size_t iseed = 0; iseed < num_seeds; iseed++)
		{
			seeds_[iseed] = new double[3];
			for (size_t idim = 0; idim < 3; idim++) seeds_[iseed][idim] = seeds[iseed * 3 + idim];
		}

		double* x = new double[4];

		MeshingSmartTree seeds_tree(3, num_seeds, seeds_);

		size_t num_invalidated(0);
		bool* valid_seed = new bool[num_seeds];
		for (size_t iseed = 0; iseed < num_seeds; iseed++) valid_seed[iseed] = true;

		for (size_t iseed = 0; iseed < num_imposed_seeds; iseed++)
		{
			for (size_t idim = 0; idim < 4; idim++) x[idim] = imposed_seeds[iseed * 4 + idim];

			// find closest point in seeds .. if distance in zero, make sure this sphere has no interior seeds
			size_t iclosest(SIZE_MAX); double hclosest(DBL_MAX);
			seeds_tree.get_closest_tree_point(x, iclosest, hclosest);

			if (hclosest < 1E-10)
			{
				size_t num(0), cap(100);
				size_t* points_in_sphere = new size_t[cap];
				seeds_tree.get_tree_points_in_sphere(x, x[3], num, points_in_sphere, cap);
				for (size_t i = 0; i < num; i++)
				{
					size_t seed_index = points_in_sphere[i];

					if (valid_seed[seed_index] && seed_index >= num_critical_seeds + num_imposed_seeds_added)
					{
						valid_seed[seed_index] = false;
						num_invalidated++;
					}
				}
				seeds_sizing[iclosest] = x[3];
			}
		}

		for (size_t iseed = 0; iseed < num_imposed_seeds; iseed++)
		{
			for (size_t idim = 0; idim < 4; idim++) x[idim] = imposed_seeds[iseed * 4 + idim];

			// find closest point in seeds .. if distance in zero, make sure this sphere has no interior seeds
			size_t iclosest(SIZE_MAX); double hclosest(DBL_MAX);
			seeds_tree.get_closest_tree_point(x, iclosest, hclosest);

			if (hclosest < 1E-10)
				continue;
			else
			{
				for (size_t idim = 0; idim < 3; idim++)
				{
					imposed_seeds[iseed * 4 + idim] = seeds[iclosest * 3 + idim];
				}
				imposed_seeds[iseed * 4 + 3] = seeds_sizing[iclosest];
			}
		}

		if (num_invalidated > 0)
		{
			size_t num_valid_seeds = num_seeds - num_invalidated;
			double* tmp = new double[num_valid_seeds * 3];
			double* tmps = new double[num_valid_seeds];
			size_t* tmp_id = new size_t[num_valid_seeds];
			size_t ivalid = 0;
			for (size_t iseed = 0; iseed < num_seeds; iseed++)
			{
				if (!valid_seed[iseed]) continue;
				for (size_t idim = 0; idim < 3; idim++) tmp[ivalid * 3 + idim] = seeds[iseed * 3 + idim];
				tmps[ivalid] = seeds_sizing[iseed];
				tmp_id[ivalid] = seeds_region_id[iseed];
				ivalid++;
			}
			num_seeds = num_valid_seeds;
			delete[] seeds; seeds = tmp;
			delete[] seeds_sizing; seeds_sizing = tmps;
			delete[] seeds_region_id; seeds_region_id = tmp_id;
		}
		

		mesher.save_spheres_csv("monitoring_points_closest_seeds.csv", 3, num_imposed_seeds, imposed_seeds, 0, 0);
		vcm_cout << "  * Updated monitoring points are saved in monitoring_points_closest_seeds.csv" << std::endl;

		delete[] imposed_seeds;

		delete[] x;

		delete[] valid_seed;

		vcm_cout << "  * executed in " << timer.report_timing() << " seconds!" << std::endl;
	}

	mesher.save_spheres_csv("seeds.csv", 3, num_seeds, seeds, seeds_sizing, seeds_region_id);
	vcm_cout << "  * All seeds are stored in seeds.csv, saved in " << timer.report_timing() << " seconds!" << std::endl;

	vcm_cout << "  * Number of added seeds " << num_seeds - num_surface_seeds << std::endl;
	vcm_cout << "  * executed in " << timer.report_timing() << " seconds!" << std::endl;

	if (num_surface_spheres > 0)
	{
		delete[] surface_spheres; delete[] surface_spheres_sizing;
	}
	if (num_edge_spheres > 0)
	{
		delete[] edge_spheres; delete[] edge_spheres_sizing;
	}
	if (num_corner_spheres > 0)
	{
		delete[] corner_spheres; delete[] corner_spheres_sizing;
	}
	if (num_sz_points > 0)
	{
		delete[] sz_points; delete[] sz_value;
	}

	if (generate_monitoring_points)
	{
		#pragma region Generate Monitoring Points:
		vcm_cout << "MeshingVoroCrust::Generating Monitoring Points:" << std::endl;
		size_t num_monitoring_points(num_seeds - num_surface_seeds);
		if (num_monitoring_points > 100) num_monitoring_points = 100;
		vcm_cout << "  * Number of Monitoring Points = " << num_monitoring_points << std::endl;

		size_t* seed_new_index = new size_t[num_seeds];
		for (size_t iseed = 0; iseed < num_seeds; iseed++) seed_new_index[iseed] = iseed;

		// shuffle interior indices
		for (size_t iseed = num_surface_seeds; iseed < num_seeds; iseed++)
		{
			double u = _rsampler.generate_uniform_random_number();
			size_t jseed = num_surface_seeds + size_t(u * (num_seeds - num_surface_seeds));
			if (jseed == num_seeds) jseed--;
			size_t tmp = seed_new_index[iseed];
			seed_new_index[iseed] = seed_new_index[jseed];
			seed_new_index[jseed] = tmp;
		}

		double* mon_pts = new double[num_monitoring_points * 4];
		for (size_t i = 0; i < num_monitoring_points; i++)
		{
			size_t seed_index = seed_new_index[num_surface_seeds + i];
			mon_pts[i * 4] = seeds[seed_index * 3];
			mon_pts[i * 4 + 1] = seeds[seed_index * 3 + 1];
			mon_pts[i * 4 + 2] = seeds[seed_index * 3 + 2];
			mon_pts[i * 4 + 3] = seeds_sizing[seed_index];
		}
		mesher.save_spheres_csv("monitoring_points.csv", 3, num_monitoring_points, mon_pts, 0, 0);
		vcm_cout << "  * Monitoring points are saved in monitoring_points.csv" << std::endl;

		delete[] seed_new_index; delete[] mon_pts;
		#pragma endregion
	}
	return 0;
	#pragma endregion
}


int MeshingVoroCrust::generate_explicit_mesh(int num_threads, size_t num_seeds, double* seeds, size_t* seeds_region_id, double* seeds_sizing,
                                      size_t num_surface_seeds, bool generate_vcg_file, bool generate_exodus_file)
{
	#pragma region Generate Explicit Mesh:
	vcm_cout << "MeshingVoroCrust::Generating Explicit Voronoi Mesh:" << std::endl;

	MeshingTimer timer;

	size_t num_regions(0);
	for (size_t iseed = 0; iseed < num_seeds; iseed++)
	{
		if (seeds_region_id[iseed] > num_regions) num_regions = seeds_region_id[iseed];
	}
	num_regions++;

	vcm_cout << "  * Number of surface seeds patches =  " << num_regions << std::endl;

	MeshingVoronoiMesher mesher;
	size_t num_vertices(0); double* vertices(0);
	size_t num_faces(0); size_t** faces(0);

	mesher.generate_3d_voronoi_mesh(num_threads, num_seeds, seeds, seeds_region_id, seeds_sizing,
		                            num_vertices, vertices, num_faces, faces);

	if (false)
	{
		double* xo = new double[3]; double* edir = new double[3];
		xo[0] = 0.4; xo[1] = 0.4; xo[2] = 0.4;
		edir[0] = 1.0; edir[1] = 1.0; edir[2] = 1.0;
		mesher.save_Voronoi_tessellation("mesh.ply", num_vertices, vertices, num_faces, faces, xo, edir);
		delete[] xo; delete[] edir;
	}
	

	// adjusting region id
	size_t* region_new_index = new size_t[num_regions];
	for (size_t iregion = 0; iregion < num_regions; iregion++) region_new_index[iregion] = iregion;

	for (size_t iface = 0; iface < num_faces; iface++)
	{
		size_t num_corners = faces[iface][0];
		size_t seed_i = faces[iface][num_corners + 1];
		size_t seed_j = faces[iface][num_corners + 2];
		if (seed_i < num_surface_seeds || seed_j < num_surface_seeds) continue;
		if (region_new_index[seeds_region_id[seed_i]] == region_new_index[seeds_region_id[seed_j]]) continue;

		// two internal neighbor seeds from different region
		if (region_new_index[seeds_region_id[seed_i]] < region_new_index[seeds_region_id[seed_j]])
			region_new_index[seeds_region_id[seed_j]] = region_new_index[seeds_region_id[seed_i]];
		else
			region_new_index[seeds_region_id[seed_i]] = region_new_index[seeds_region_id[seed_j]];
	}

	// adjust new region ids
	for (size_t iregion = 0; iregion < num_regions; iregion++)
	{
		if (region_new_index[iregion] == iregion) continue;
		size_t id = region_new_index[iregion];
		while (region_new_index[id] != id) id = region_new_index[id]; // CHECK THAT LINE
		region_new_index[iregion] = id;
	}

	size_t num_active_regions(0);
	for (size_t iregion = 0; iregion < num_regions; iregion++)
	{
		if (region_new_index[iregion] == iregion) num_active_regions++;
	}
	num_active_regions--;

	vcm_cout << "  * Number of volumetric regions =  " << num_active_regions << std::endl;

	for (size_t iseed = 0; iseed < num_seeds; iseed++)
	{
		seeds_region_id[iseed] = region_new_index[seeds_region_id[iseed]];
	}

	delete[] region_new_index;

	size_t* cell_num_faces = new size_t[num_seeds];
	for (size_t iseed = 0; iseed < num_seeds; iseed++) cell_num_faces[iseed] = 0;

	for (size_t iface = 0; iface < num_faces; iface++)
	{
		size_t num_corners = faces[iface][0];
		size_t seed_i = faces[iface][num_corners + 1];
		size_t seed_j = faces[iface][num_corners + 2];
		if (seeds_region_id[seed_i] != 0) cell_num_faces[seed_i]++;
		if (seeds_region_id[seed_j] != 0) cell_num_faces[seed_j]++;
	}

	size_t num_meshed_cells(0);
	size_t** cell_faces = new size_t*[num_seeds];
	for (size_t iseed = 0; iseed < num_seeds; iseed++)
	{
		if (seeds_region_id[iseed] == 0) continue;
		cell_faces[iseed] = new size_t[cell_num_faces[iseed]];
		cell_num_faces[iseed] = 0;
		num_meshed_cells++;
	}

	for (size_t iface = 0; iface < num_faces; iface++)
	{
		size_t num_corners = faces[iface][0];
		size_t seed_i = faces[iface][num_corners + 1];
		size_t seed_j = faces[iface][num_corners + 2];
		if (seeds_region_id[seed_i] != 0)
		{
			cell_faces[seed_i][cell_num_faces[seed_i]] = iface; cell_num_faces[seed_i]++;
		}
		if (seeds_region_id[seed_j] != 0)
		{
			cell_faces[seed_j][cell_num_faces[seed_j]] = iface; cell_num_faces[seed_j]++;
		}
	}

	vcm_cout << "  * Number of seeds =  " << num_seeds << std::endl;
	vcm_cout << "  * Number of meshed cells =  " << num_meshed_cells << std::endl;
	vcm_cout << "  * executed in " << timer.report_timing() << " seconds!" << std::endl;

	if (generate_vcg_file)
		mesher.save_vcg_file_pflotran("mesh.vcg", num_seeds, seeds, seeds_region_id, num_vertices, vertices, num_faces, faces, cell_num_faces, cell_faces);

	if (generate_exodus_file)
		mesher.save_exodus_file("mesh.exo", num_seeds, seeds, seeds_region_id, num_vertices, vertices, num_faces, faces, cell_num_faces, cell_faces);

	for (size_t iface = 0; iface < num_faces; iface++) delete[] faces[iface];
	delete[] faces; delete[] vertices; delete[] cell_num_faces;

	for (size_t iseed = 0; iseed < num_seeds; iseed++)
	{
		if (seeds_region_id[iseed] == 0) continue;
		delete[] cell_faces[iseed];
	}
	delete[] cell_faces;

	delete[] seeds; delete[] seeds_region_id; delete[] seeds_sizing;
	return 0;
	#pragma endregion
}


int MeshingVoroCrust::generate_mesh_from_seeds_file(std::string file_name, int num_threads, bool generate_vcg_file, bool generate_exodus_file)
{
	size_t num_seeds;
	double* seeds; size_t* seeds_region_id; double* seeds_sizing;
	read_seeds_csv(file_name, num_seeds, seeds, seeds_region_id, seeds_sizing);
	generate_explicit_mesh(num_threads, num_seeds, seeds, seeds_region_id, seeds_sizing, num_seeds, generate_vcg_file, generate_exodus_file);
	return 0;
}


int MeshingVoroCrust::read_seeds_csv(std::string file_name, size_t& num_seeds, double* &seeds, size_t* &seeds_region_id, double* &seeds_sizing)
{
	#pragma region Read Seeds from seeds.csv:
	vcm_cout << "MeshingVoroCrust::Reeding seeds from " << file_name.c_str() << std::endl;

	MeshingTimer timer;

	std::ifstream tmpfile(file_name.c_str());
	num_seeds = 0;
	if (tmpfile.is_open() && tmpfile.good())
	{
		std::string line = "";
		getline(tmpfile, line);
		while (getline(tmpfile, line))
		{
			std::istringstream iss(line);
			std::vector<std::string> tokens;
			get_tokens(line, ',', tokens);

			if (tokens.size() == 0) continue;

			num_seeds++;
		}
	}

	std::ifstream tmpfile_2(file_name.c_str());

	seeds = new double[num_seeds * 3];
	seeds_sizing = new double[num_seeds];
	seeds_region_id = new size_t[num_seeds];

	num_seeds = 0;
	if (tmpfile_2.is_open() && tmpfile_2.good())
	{
		std::string line = "";
		getline(tmpfile_2, line);
		while (getline(tmpfile_2, line))
		{
			std::istringstream iss(line);
			std::vector<std::string> tokens;
			get_tokens(line, ',', tokens);

			if (tokens.size() == 0) continue;

			for (size_t idim = 0; idim < 3; idim++) seeds[num_seeds * 3 + idim] = _memo.string_to_double(tokens[idim]);
			seeds_sizing[num_seeds] = _memo.string_to_double(tokens[3]);
			seeds_region_id[num_seeds] = _memo.string_to_size_t(tokens[4]);
			num_seeds++;
		}
	}

	vcm_cout << "  * executed in = " << timer.report_timing() << " seconds! ***" << std::endl;
	return 0;
	#pragma endregion
}

int MeshingVoroCrust::get_tokens(std::string line, char separator, std::vector<std::string>& tokens)
{
	#pragma region Get Tokens from a line:
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



