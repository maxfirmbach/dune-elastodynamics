/* config.h.  Generated from config_collected.h.cmake by CMake.
   It was generated from config_collected.h.cmake which in turn is generated automatically
   from the config.h.cmake files of modules this module depends on. */

/* Define to 1 if you have module dune-elastodynamics available */
#define HAVE_DUNE_ELASTODYNAMICS 1


/* Define to 1 if you have module dune-common available */
#define HAVE_DUNE_COMMON 1


/* Define to 1 if you have module dune-uggrid available */
#define HAVE_DUNE_UGGRID 1


/* Define to 1 if you have module dune-geometry available */
#define HAVE_DUNE_GEOMETRY 1


/* Define to 1 if you have module dune-python available */
#define HAVE_DUNE_PYTHON 0


/* Define to 1 if you have module dune-typetree available */
#define HAVE_DUNE_TYPETREE 1


/* Define to 1 if you have module dune-istl available */
#define HAVE_DUNE_ISTL 1


/* Define to 1 if you have module dune-grid available */
#define HAVE_DUNE_GRID 1


/* Define to 1 if you have module dune-localfunctions available */
#define HAVE_DUNE_LOCALFUNCTIONS 1


/* Define to 1 if you have module dune-foamgrid available */
#define HAVE_DUNE_FOAMGRID 1


/* Define to 1 if you have module dune-functions available */
#define HAVE_DUNE_FUNCTIONS 1


/* begin private */
/* Define to the version of dune-common */
#define DUNE_COMMON_VERSION "2.8-git"

/* Define to the major version of dune-common */
#define DUNE_COMMON_VERSION_MAJOR 2

/* Define to the minor version of dune-common */
#define DUNE_COMMON_VERSION_MINOR 8

/* Define to the revision of dune-common */
#define DUNE_COMMON_VERSION_REVISION 0

/* Standard debug streams with a level below will collapse to doing nothing */
#define DUNE_MINIMAL_DEBUG_LEVEL 4

/* does the compiler support __attribute__((deprecated))? */
#define HAS_ATTRIBUTE_DEPRECATED 1

/* does the compiler support __attribute__((deprecated("message"))? */
#define HAS_ATTRIBUTE_DEPRECATED_MSG 1

/* does the compiler support __attribute__((unused))? */
#define HAS_ATTRIBUTE_UNUSED 1

/* does the standard library provide experimental::make_array() ? */
#define DUNE_HAVE_CXX_EXPERIMENTAL_MAKE_ARRAY 1

/* does the standard library provide experimental::is_detected ? */
#define DUNE_HAVE_CXX_EXPERIMENTAL_IS_DETECTED 1

/* does the standard library provide identity ? */
/* #undef DUNE_HAVE_CXX_STD_IDENTITY */

/* Define if you have a BLAS library. */
#define HAVE_BLAS 1

/* Define if you have LAPACK library. */
#define HAVE_LAPACK 1

/* Define if you have the MPI library.  */
#define HAVE_MPI ENABLE_MPI

/* Deactivate cxx bindings for MPI */
#if HAVE_MPI
#define MPICH_SKIP_MPICXX 1
#define OMPI_SKIP_MPICXX 1
#define MPI_NO_CPPBIND 1
#define MPIPP_H
#define _MPICC_H
#endif

/* Define if you have the GNU GMP library. The value should be ENABLE_GMP
   to facilitate activating and deactivating GMP using compile flags. */
#define HAVE_GMP ENABLE_GMP

/* Define if you have the GCC Quad-Precision library. The value should be ENABLE_QUADMATH
   to facilitate activating and deactivating QuadMath using compile flags. */
#define HAVE_QUADMATH ENABLE_QUADMATH

/* Define if you have the Vc library. The value should be ENABLE_VC
   to facilitate activating and deactivating Vc using compile flags. */
#define HAVE_VC ENABLE_VC

/* Define to 1 if you have the Threading Building Blocks (TBB) library */
#define HAVE_TBB 1




/* old feature support macros which were tested until 2.7, kept around for one more release */
/* As these are now always supported due to the new compiler requirements, they are directly */
/* defined without an explicit test. */
#define DUNE_HAVE_CXX_CLASS_TEMPLATE_ARGUMENT_DEDUCTION 1
#define DUNE_HAVE_CXX_OPTIONAL 1
#define DUNE_HAVE_CXX_VARIANT 1
#define DUNE_SUPPORTS_CXX_THROW_IN_CONSTEXPR 1
#define DUNE_HAVE_C_ALIGNED_ALLOC 1
#define DUNE_HAVE_CXX_BOOL_CONSTANT 1
#define DUNE_HAVE_CXX_EXPERIMENTAL_BOOL_CONSTANT 0
#define DUNE_HAVE_HEADER_EXPERIMENTAL_TYPE_TRAITS 0
#define DUNE_HAVE_CXX_APPLY 1
#define DUNE_HAVE_CXX_EXPERIMENTAL_APPLY 0
#define HAVE_IS_INDEXABLE_SUPPORT 1

/* Define to ENABLE_UMFPACK if the UMFPack library is available */
#define HAVE_UMFPACK ENABLE_SUITESPARSE

/* Define to ENABLE_SUITESPARSE if the SuiteSparse library is available */
#define HAVE_SUITESPARSE ENABLE_SUITESPARSE

/* Define to ENABLE_SUITESPARSE if the SuiteSparse's AMD library is available */
/* #undef HAVE_SUITESPARSE_AMD */

/* Define to ENABLE_SUITESPARSE if the SuiteSparse's BTF library is available */
/* #undef HAVE_SUITESPARSE_BTF */

/* Define to ENABLE_SUITESPARSE if the SuiteSparse's CAMD library is available */
/* #undef HAVE_SUITESPARSE_CAMD */

/* Define to ENABLE_SUITESPARSE if the SuiteSparse's CCOLAMD library is available */
/* #undef HAVE_SUITESPARSE_CCOLAMD */

/* Define to ENABLE_SUITESPARSE if the SuiteSparse's CHOLMOD library is available */
#define HAVE_SUITESPARSE_CHOLMOD ENABLE_SUITESPARSE

/* Define to ENABLE_SUITESPARSE if the SuiteSparse's COLAMD library is available */
/* #undef HAVE_SUITESPARSE_COLAMD */

/* Define to ENABLE_SUITESPARSE if the SuiteSparse's CXSPARSE library is available */
/* #undef HAVE_SUITESPARSE_CXSPARSE */

/* Define to ENABLE_SUITESPARSE if the SuiteSparse's KLU library is available */
/* #undef HAVE_SUITESPARSE_KLU */

/* Define to ENABLE_SUITESPARSE if the SuiteSparse's LDL library is available */
#define HAVE_SUITESPARSE_LDL ENABLE_SUITESPARSE

/* Define to ENABLE_SUITESPARSE if the SuiteSparse's RBIO library is available */
/* #undef HAVE_SUITESPARSE_RBIO */

/* Define to ENABLE_SUITESPARSE if the SuiteSparse's SPQR library is available
   and if it's version is at least 4.3 */
#define HAVE_SUITESPARSE_SPQR ENABLE_SUITESPARSE

/* Define to ENABLE_SUITESPARSE if the SuiteSparse's UMFPACK library is available */
#define HAVE_SUITESPARSE_UMFPACK ENABLE_SUITESPARSE

/* Define to 1 if METIS is available */
#define HAVE_METIS 1

/* Define to 1 if the Scotch replacement for METIS is used. */
/* #undef HAVE_SCOTCH_METIS */

/* Define to 1 if you have the ParMETIS library. */
#define HAVE_PARMETIS 1

/* Define to 1 if the PTScotch replacement for ParMETIS is used. */
/* #undef HAVE_PTSCOTCH_PARMETIS */

/* Define to 1 if PT-Scotch is available */
#define HAVE_PTSCOTCH 1

/* Used to call lapack functions */
#define LAPACK_NEEDS_UNDERLINE



/* Define to the version of dune-common */
#define DUNE_UGGRID_VERSION "2.8-git"

/* Define to the major version of dune-common */
#define DUNE_UGGRID_VERSION_MAJOR 2

/* Define to the minor version of dune-common */
#define DUNE_UGGRID_VERSION_MINOR 8

/* Define to the revision of dune-common */
#define DUNE_UGGRID_VERSION_REVISION 0

/* begin private section */

/* see parallel/ddd/dddi.h */
/* #undef DDD_MAX_PROCBITS_IN_GID */

/* Define to 1 if you can safely include both <sys/time.h> and <time.h>. */
/* #undef TIME_WITH_SYS_TIME */

/* Define to 1 if UGGrid should use the complete set of green refinement rules for tetrahedra */
/* #undef DUNE_UGGRID_TET_RULESET */

/* end private section */





/* Define to the version of dune-geometry */
#define DUNE_GEOMETRY_VERSION "2.8-git"

/* Define to the major version of dune-geometry */
#define DUNE_GEOMETRY_VERSION_MAJOR 2

/* Define to the minor version of dune-geometry */
#define DUNE_GEOMETRY_VERSION_MINOR 8

/* Define to the revision of dune-geometry */
#define DUNE_GEOMETRY_VERSION_REVISION 0





/* Define to the version of dune-typetree */
#define DUNE_TYPETREE_VERSION "2.8-git"

/* Define to the major version of dune-typetree */
#define DUNE_TYPETREE_VERSION_MAJOR 2

/* Define to the minor version of dune-typetree */
#define DUNE_TYPETREE_VERSION_MINOR 8

/* Define to the revision of dune-typetree */
#define DUNE_TYPETREE_VERSION_REVISION 0






/* Define to ENABLE_SUPERLU if the SuperLU library is available */
#define HAVE_SUPERLU ENABLE_SUPERLU

/* Define to the integer type that SuperLU was compiled for
   See e.g. what int_t is defined to in slu_sdefs.h */
#define SUPERLU_INT_TYPE int

/* Define to ENABLE_ARPACKPP if the ARPACK++ library is available */
#define HAVE_ARPACKPP ENABLE_ARPACKPP

/* Define to the version of dune-istl */
#define DUNE_ISTL_VERSION "2.8-git"

/* Define to the major version of dune-istl */
#define DUNE_ISTL_VERSION_MAJOR 2

/* Define to the minor version of dune-istl */
#define DUNE_ISTL_VERSION_MINOR 8

/* Define to the revision of dune-istl */
#define DUNE_ISTL_VERSION_REVISION 0

/* Enable/Disable the backwards compatibility of the category enum/method in dune-istl solvers, preconditioner, etc. */
#define DUNE_ISTL_SUPPORT_OLD_CATEGORY_INTERFACE 1





/* Define to the version of dune-grid */
#define DUNE_GRID_VERSION "2.8-git"

/* Define to the major version of dune-grid */
#define DUNE_GRID_VERSION_MAJOR 2

/* Define to the minor version of dune-grid */
#define DUNE_GRID_VERSION_MINOR 8

/* Define to the revision of dune-grid */
#define DUNE_GRID_VERSION_REVISION 0

/* Define to 1 if psurface library is found */
#define HAVE_PSURFACE 1

/* Define to 1 if AmiraMesh library is found */
/* #undef HAVE_AMIRAMESH */

/* Define to 1 if you have at least psurface version 2.0 */
/* #undef HAVE_PSURFACE_2_0 */

/* Alberta version found by configure, either 0x200 for 2.0 or 0x300 for 3.0 */
/* #undef DUNE_ALBERTA_VERSION */

/* This is only true if alberta-library was found by configure _and_ if the
   application uses the ALBERTA_CPPFLAGS */
#define HAVE_ALBERTA ENABLE_ALBERTA

/* This is only true if UG was found by configure _and_ if the application
   uses the UG_CPPFLAGS */
#define HAVE_UG ENABLE_UG

/* Define to 1 if you have mkstemp function */
#define HAVE_MKSTEMP 1







/* Define to the version of dune-localfunctions */
#define DUNE_LOCALFUNCTIONS_VERSION "2.8-git"

/* Define to the major version of dune-localfunctions */
#define DUNE_LOCALFUNCTIONS_VERSION_MAJOR 2

/* Define to the minor version of dune-localfunctions */
#define DUNE_LOCALFUNCTIONS_VERSION_MINOR 8

/* Define to the revision of dune-localfunctions */
#define DUNE_LOCALFUNCTIONS_VERSION_REVISION 0






/* Define to the version of dune-functions */
#define DUNE_FUNCTIONS_VERSION "2.8-git"

/* Define to the major version of dune-functions */
#define DUNE_FUNCTIONS_VERSION_MAJOR 2

/* Define to the minor version of dune-functions */
#define DUNE_FUNCTIONS_VERSION_MINOR 8

/* Define to the revision of dune-functions */
#define DUNE_FUNCTIONS_VERSION_REVISION 0



/* begin dune-elastodynamics
   put the definitions for config.h specific to
   your project here. Everything above will be
   overwritten
*/

/* begin private */
/* Name of package */
#define PACKAGE "dune-elastodynamics"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "max.firmbach@unibw.de"

/* Define to the full name of this package. */
#define PACKAGE_NAME "dune-elastodynamics"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "dune-elastodynamics 1.0"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "dune-elastodynamics"

/* Define to the home page for this package. */
#define PACKAGE_URL ""

/* Define to the version of this package. */
#define PACKAGE_VERSION "1.0"

/* end private */

/* Define to the version of dune-elastodynamics */
#define DUNE_ELASTODYNAMICS_VERSION "1.0"

/* Define to the major version of dune-elastodynamics */
#define DUNE_ELASTODYNAMICS_VERSION_MAJOR 1

/* Define to the minor version of dune-elastodynamics */
#define DUNE_ELASTODYNAMICS_VERSION_MINOR 0

/* Define to the revision of dune-elastodynamics */
#define DUNE_ELASTODYNAMICS_VERSION_REVISION 0

/* end dune-elastodynamics
   Everything below here will be overwritten
*/ 

/* Grid type magic for DGF parser */

/* ALBERTAGRID not available, enable with cmake variable DUNE_GRID_GRIDTYPE_SELECTOR=ON */
/* UGGRID not available, enable with cmake variable DUNE_GRID_GRIDTYPE_SELECTOR=ON */
/* ONEDGRID not available, enable with cmake variable DUNE_GRID_GRIDTYPE_SELECTOR=ON */
/* YASPGRID not available, enable with cmake variable DUNE_GRID_GRIDTYPE_SELECTOR=ON */

