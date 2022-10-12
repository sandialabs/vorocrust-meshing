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

#include<iostream>
#include<cmath>
#include<cstdio>
#include<cstdlib>
#include<cstring>

#include <Kokkos_Core.hpp>
#include <KokkosBlas.hpp>

#define N 10

typedef Kokkos::View<int*> ViewIntVector_T;

void print_view(const std::string& label, const ViewIntVector_T& view)
{
    std::cout << label << ": ";
    for(int i=0; i<N; i++)
    {
        std::cout << view(i) << " ";
    }
    std::cout << std::endl;
}

int main(int argc, char** argv)
{
    Kokkos::initialize(argc, argv);
    {
        std::cout << "-- Kokkos Configuration --" << std::endl;
        Kokkos::print_configuration(std::cout);

        std::cout << "-- Initialize Data --" << std::endl;

        ViewIntVector_T input("input",   N);
        ViewIntVector_T output("output", N);

        for(int i=0; i<N; i++)
        {
            input(i)  = -1;
            output(i) =  99;
        }

        print_view("input",  input);
        print_view("output", output);

        std::cout << "-- compute abs(output, input) --" << std::endl;

        // See: https://github.com/kokkos/kokkos-kernels/blob/master/src/blas/KokkosBlas1_abs.hpp
        KokkosBlas::abs(output, input);

        print_view("input",  input);
        print_view("output", output);

        std::cout << "-- validate results --" << std::endl;

        for(int i=0; i<N; i++)
        {
            const int actual = output(i);
            const int expect = std::abs(input(i));
            assert( actual == expect );
        }
        std::cout << "OK" << std::endl;
    }

    Kokkos::finalize();
    return 0;
}


