
#include <iostream>

#include "MeshingOptionParser.h"

int main(const int argc, char* argv[])
{
    MeshingOptionParser op(argc, argv);

    op.PrettyPrint();

    std::cout << "Done." << std::endl;
    return 0;
}
