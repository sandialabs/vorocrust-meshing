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

#include <MeshingOptionParser.h>
#include <MeshingVoroCrust.h>


#ifdef USE_KOKKOS
	#include <Kokkos_Core.hpp>
#endif



class OptionParserVCMesh : public MeshingOptionParser
{
public:


	OptionParserVCMesh(const int argc, char* argv[])
	: MeshingOptionParser(argc, argv)
	{
		this->parse();
	}


	virtual bool parse_options(void)
	{
		if( this->option_exists("-vc") )
		{
			this->method   = this->method_type_vc();
			this->filename = this->get_option_value("-vc", true);
			return true;
		}
		return false;
	}


	virtual void PrettyPrint()
	{
		std::cout << "method: ";

		if(this->method == method_type_vc())
		{
			std::cout << "VC";
		}
		else
		{
			std::cout << "UNKNOWN";
		}
		std::cout << std::endl;

		std::cout << "Filename: " << this->filename << std::endl;
	}


	virtual void usage() const
	{
		std::cout << "Usage:" << std::endl
			      << "  " << this->appname << " [-h] <method> <input_filename>" << std::endl
			      << std::endl
			      << "  " << "method:" << std::endl
			      << "  " << "  -vc    : Use the vc method" << std::endl
				  << "  " << "  -h     : Display help message and exit." << std::endl
				  << std::endl;
	}


	// Method Type Helpers
	int method_type_vc(void) { return 1; }

};    // class OptionParserVCMesh



int main(int argc, char* argv[])
{
	#ifdef USE_KOKKOS
	Kokkos::initialize(argc, argv);
	{
	#endif // USE_KOKKOS

		vcm_cout << vcm_cout.generate_banner();

		OptionParserVCMesh options(argc, argv);

		// Select and run the appropriate method
		if( options.get_method() == options.method_type_vc())
		{
			MeshingVoroCrust vc(options.get_filename());
			vc.execute();
		}

	#ifdef USE_KOKKOS
	}
	Kokkos::finalize();
	#endif // USE_KOKKOS

	return 0;
}

