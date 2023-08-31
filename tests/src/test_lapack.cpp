#include <iostream>

extern "C"
{
    void dgesv_(int* n, int* nrhs, double* a, int* lda, int* ipivot, double* b, int* ldb, int* info);
}

#define MAX 10

int
main()
{
    // Values needed for dgesv
    int    n;
    int    nrhs = 1;
    double a[ MAX ][ MAX ];
    double b[ 1 ][ MAX ];
    int    lda = MAX;
    int    ldb = MAX;
    int    ipiv[ MAX ];
    int    info;
    // Other values
    int i, j;


    // Let's make a 4x4 matrix
    n = 4;

    // Initialize it to 0's
    for(i=0; i<n; i++)
    {
        for(j=0; j<n; j++)
        {
            a[j][i] = 0;
        }
        b[0][i] = 0;
    }

    // Set identity matrix
    for(i=0; i<n; i++)
        a[i][i] = 1;

    // set up b
    for(i=0; i<n; i++)
        b[0][i] = 1;

    // print out the matrix
    for(i = 0; i < n; i++)
    {
        std::cout << "row " << i << ":  [";
        for(j = 0; j < n; j++)
        {
            std::cout << a[j][i] << " ";
        }
        std::cout << "]   [ " << b[0][i] << " ]" << std::endl;
    }

    // Solve the linear system
    dgesv_(&n, &nrhs, &a[ 0 ][ 0 ], &lda, ipiv, &b[ 0 ][ 0 ], &ldb, &info);

    // Check for success
    if(info == 0)
    {
        // Write the answer
        std::cout << "The answer is\n";
        for(i = 0; i < n; i++)
        {
            std::cout << "b[" << i << "]\t" << b[ 0 ][ i ] << "\n";
        }
    }
    else
    {
        // Write an error message
        std::cerr << "dgesv returned error " << info << "\n";
    }
    return info;
}