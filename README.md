VoroCrust-Meshing
=========
VoroCrust-Meshing is the first provably correct algorithm for conforming Voronoi meshing
of non-convex and non-manifold domains with guarantees on the quality of both
surface and volume elements.

SRC Number
----------------------
VoroCrust-Meshing is SCR 2820.0. It was approved for open source release 8/30/2022.


How to Build VoroCrust-Meshing
----------------------
### Prerequisites for Building

- **CMake** : Version 3.12.x or later preferred.
  - `$ cmake --version` can be run to get the version of the currently installed
    version of CMake.
  - This can be downloaded at https://www.cmake.org
  - On RHEL / CENTOS systems: `sudo yum install cmake`
  - On OSX if using homebrew: `brew install cmake`. 
  - On Sandia systems with SEMS modules installed:
    `module load sems-env && module load sems-cmake/3.12.2`
- **C++ Compiler** : This can vary from platform to platform. Windows systems
  will likely be using Visual Studio. *nix systems and OSX (terminal) will
  have their own compilers.
- **OpenMP** : The [OpenMP library][3] is required by VoroCrust. Most modern
  compilers will have openmp installed already.
  If it is installed then CMake should detect it automatically.
  - _Warning for OSX users_: The default compiler that comes with XCode (Clang)
    does not come from Apple with OpenMP installed. The best way to typically
    deal with this is to install GCC via [homebrew][1] or [macports][2]. If you
    really wish to use the Apple compiler, you will need to build and install
    OpenMP somehow.  Since Apple uses Clang (LLVM), the best reference to how
    to do this can be found at the [LLVM OpenMP Support page][4].
    
    NOTE: If using homebrew on computer with Monterey and an M1 chip: Go into gcc 
    bin directory (e.g /opt/homebrew/Cellar/gcc/12.1.0/bin)  and set up aliases 
    for c++, cpp, g++, and gcc using the commands: 
    
    ln -s c++-12 c++

    ln -s cpp-12 cpp 

    etc., 

    or vc-meshing will not find the correct compilers. 

Once you have all the required prerequisites installed for VoroCrust, you can
build it.

### Preferred directory structure
We will use the following directory structure for our example of how to
properly build VoroCrust-meshing. You are not _required_ to use this folder 
structure but some automated testing may not work if you choose a different 
layout.

Generally speaking, we would create a directory called `vorocrust-meshing` into 
which we will clone the VoroCrust git repositories, and we create a `build` 
directory in the vorocrust source folder.

The directory hierarchy should look something like this:

```
vorocrust-meshing/
     |
     +--- vorocrust-meshing/
     |    |
     |    + build/
     |
     +--- vorocrust-data/
     |
     +--- vorocrust-workshop/
```

_Note_: the vorocrust-data and vorocrust-workshop directories aren't required
for building at this time. We are currently working on populating the
vorocrust-data repository with sample meshes, etc.
The vorocrust-workshop repository contains some powerpoint presentations from
previous workshops, etc. We split this into a separate repository because the
files can be quite large.


### Building on *nix / OSX
1. Clone Vorocrust-Meshing
    ```bash
    $ mkdir vorocrust-meshing
    $ cd vorocrust-meshing
    $ git clone git@github.com:sandialabs/vorocrust-meshing.git
    ```
2. Create your build directory
    ```bash
    $ cd vorocrust-meshing
    $ mkdir build
    $ cd build
    ```
3. Configure
    ```bash
    $ cmake ../.
    ```
    On Macs with M1 chips it may be necessary to force cmake to find the correct
    compilers using
    ```bash
    $ cmake -D CMAKE_CXX_COMPILER=g++ -D CMAKE_C_COMPILER=gcc ../.
    ```
    
4. Build
    ```bash
    $ make
    ```
    Optionally you can add the arguments `-j N` to use N cores in a parallel
    build on your computer. This can speed up builds.

    If you would like to see the exact commands from the build, you can also
    run `make VERBOSE=1` to increase the verbosity of the makefile.

5. Test
    ```
    $ make test
    ```
    Currently, cmake expects vorocrust-data to be checked out in the same
    directory as the vorocrust repository. See the preferred directory structure
    above.

### Building on Windows
We are working on a website that provides more detailed steps on building
VoroCrust-Meshing for Windows. When that is available we will add a link here.


CMake Options for VoroCrust-Meshing
---------------------------
There are quite a number of CMake options for VoroCrust-Meshing that can be used to configure
the build. Here is a list of them:

| CMake Option                         | Default Value | Advanced | Description                                                                  |
|--------------------------------------|:-------------:|:--------:|------------------------------------------------------------------------------|                                          |
| `VOROCRUST_BUILD_MESH`               |      ON       |    NO    | Include `src/Mesh` codes                                                     |
| `VOROCRUST_ENABLE_EXODUS`            |      OFF      |    NO    | Build EXODUS mesh output capability (for Mesh)                               |
| `VOROCRUST_ENABLE_LONG_TESTS`        |      OFF      |    NO    | Enables some longer-running tests.                                           |
| `VOROCRUST_ENABLE_MPI`               |      OFF      |    NO    | Enable MPI support and build with MPI compilers.                             |
| `VOROCRUST_ENABLE_OPENMP`            |      ON       |    NO    | Enable OpenMP support.                                                       |
| `VOROCRUST_TPL_BUILD_OPENBLAS`       |      OFF      |    NO    | Build OpenBLAS (ExternalProject)                                             |
| `VOROCRUST_TPL_ENABLE_KOKKOS`        |      OFF      |    NO    | Build with Kokkos (FetchContent)                                             |
| `VOROCRUST_TPL_ENABLE_KOKKOSKERNELS` |      OFF      |    NO    | Build with KokkosKernels (FetchContent)                                      |
| `VOROCRUST_TPL_USE_LAPACK`           |      OFF      |    NO    | Build with LAPACK support.                                                   |
| `VOROCRUST_VERBOSE_CMAKE`            |      OFF      |    NO    | Enable extra verbosity in CMake.                                             |
| `VOROCRUST_ENABLE_TPL_TESTS`         |      OFF      |   YES    | Run test suites for TPLs built via ExternalProject                           |
| `VOROCRUST_ENABLE_TPL_LOG_BUILD`     |      OFF      |   YES    | Enable ExternalProject TPL *build* logging.                                  |
| `VOROCRUST_ENABLE_TPL_LOG_CONF`      |      OFF      |   YES    | Enable ExternalProject TPL *configure* logging.                              |
| `VOROCRUST_ENABLE_TPL_LOG_DOWNLOAD`  |      OFF      |   YES    | Enable ExternalProject TPL *download* logging.                               |
| `VOROCRUST_ENABLE_TPL_LOG_INST`      |      OFF      |   YES    | Enable ExternalProject TPL *installation* logging.                           |
| `VOROCRUST_ENABLE_TPL_LOG_TEST`      |      OFF      |   YES    | Enable ExternalProject TPL *test* logging.                                   |
| `VOROCRUST_ENABLE_WARNINGS`          |      OFF      |   YES    | Enable extra C++ warnings (i.e., -Wall)                                      |
| `VOROCRUST_ENABLE_TUCKERMPI`         |      OFF      |   YES    | Enable TuckerMPI (EXPERIMENTAL!!!)                                           |
| `VOROCRUST_TPL_BUILD_HDF5`           |      OFF      |   YES    | Builds the HDF5 TPL. Requires ZLib. (ExternalProject)                        |
| `VOROCRUST_TPL_BUILD_NETCDF`         |      OFF      |   YES    | Builds the NetCDF TPL. Requires HDF5. (ExternalProject)                      |
| `VOROCRUST_TPL_BUILD_SEACAS`         |      OFF      |   YES    | Builds the SEACAS TPL. Requires NetCDF. (ExternalProject)                    |
| `VOROCRUST_TPL_BUILD_OPENBLAS`       |      OFF      |   YES    | Builds the OpenBLAS TPL (ExternalProject)                                    |
| `VOROCRUST_TPL_BUILD_ZLIB`           |      OFF      |   YES    | Builds the ZLib TPL (ExternalProject)                                        |
| `VOROCRUST_TPL_KOKKOS_VERSION`       |    3.5.00     |   YES    | The 'release' version of Kokkos and Kokkos-Kernels to include.               |
| `VOROCRUST_VERSION_USES_GIT_SHA1`    |      OFF      |   YES    | Attempt to find and use the current Git SHA1, commit date for version info.  |
| `VOROCRUST_TPL_NO_DOWNLOAD`          |      ON       |   YES    | If a TPL is available in the `tpls/` dir use it rather than downloading.     |
| `VOROCRUST_ENABLE_DEVTEST`           |      OFF      |   YES    | Enable EXPERIMENTAL features that are guarded. Sometimes platform dependent. |

*Notes*
- Enabling `LAPACK` will make VoroCrust try and _find_ installed BLAS and LAPACK libraries but if they can't
  be found then it will attempt to build OpenBLAS.
- Enabling `OPENBLAS` will also enable `LAPACK`.
- If `LAPACK` and `OPENBLAS` are enabled, then VoroCrust will skip detection of an installed LAPACK and will
  build its own copy of OpenBLAS.
- Enabling `EXODUS` will cause the following TPLs to be built and installed in the Build directory:
  - zlib
  - hdf5 - requires zlib
  - netcdf - requires hdf5
  - SEACS - just the Exodus library part. This depends on HDF5 and NetCDF.
    This will generate a libexodus file which will be linked into VoroCrust.
- Some features are enabled on some platforms and disabled on others.
  We use the `ENABLE_DEVTEST` flag to guard this in many cases so Devops can
  enable a feature while testing it out while retaining cmake errors that will
  prevent an experimental feature from being used outside of Devops testing.

Citing VoroCrust
----------------
Please use the following BibTeX entry to cite VoroCrust:

```bibtex
@article{10.1145/3337680,
    author = {Abdelkader, Ahmed and
              Bajaj, Chandrajit L. and
              Ebeida, Mohamed S. and
              Mahmoud, Ahmed H. and
              Mitchell, Scott A. and
              Owens, John D. and
              Rushdi, Ahmad A.},
    title = {VoroCrust: Voronoi Meshing Without Clipping},
    year = {2020},
    issue_date = {May 2020},
    publisher = {Association for Computing Machinery},
    address = {New York, NY, USA},
    volume = {39},
    number = {3},
    issn = {0730-0301},
    url = {https://doi.org/10.1145/3337680},
    doi = {10.1145/3337680},
    journal = {ACM Trans. Graph.},
    month = may,
    articleno = {23},
    numpages = {16},
    keywords = {meshing, Voronoi,
                sharp features,
                refinement,
                union of balls,
                slivers,
                Poisson-disk sampling}
}
```


[1]: https://brew.sh
[2]: https://www.macports.org
[3]: https://www.openmp.org
[4]: https://openmp.llvm.org

<!-- TEST COMMENT -->


