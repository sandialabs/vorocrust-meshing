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


