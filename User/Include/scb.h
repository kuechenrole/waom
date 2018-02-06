/*
** svn $Id$
*******************************************************************************
<<<<<<< HEAD
** Copyright (c) 2002-2016 The ROMS/TOMS Group                               **
=======
** Copyright (c) 2002-2013 The ROMS/TOMS Group                               **
>>>>>>> master
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for Southern California Bight.
**
** Application flag:   SCB
** Input script:       ocean_scb.in
*/


#undef  AFT_EIGENMODES          /* Adjoint Finite Time Eigenmodes */
#undef  CORRELATION             /* Background-error Correlation Check */
#undef  FORCING_SV              /* Forcing Singular Vectors */
#undef  FT_EIGENMODES           /* Finite Time Eigenmodes */
#undef  IS4DVAR                 /* Incremental, strong constraint 4DVAR */
#define NLM_DRIVER              /* Nonlinear Basic State trajectory */
#undef  OPT_PERTURBATION        /* Optimal perturbations */
#undef  PICARD_TEST             /* Picard Iterations Test */
#undef  R_SYMMETRY              /* Representer Matrix Symmetry Test */
#undef  S4DVAR                  /* Strong constraint 4DVAR */
#undef  SANITY_CHECK            /* Sanity Check */
#undef  SO_SEMI                 /* Stochastic Optimals: Semi-norm */
#undef  TLM_CHECK               /* Tangent Linear Model Check */
#undef  W4DPSAS                 /* Weak constraint 4D-PSAS */
#undef  W4DVAR                  /* Weak constraint 4DVAR */

/*
**-----------------------------------------------------------------------------
**  Nonlinear basic state tracjectory.
**-----------------------------------------------------------------------------
*/

#if defined NLM_DRIVER
# define UV_ADV
# define UV_COR
# define UV_LDRAG
# define UV_VIS2
# define MIX_S_UV
<<<<<<< HEAD
# define SPLINES_VDIFF
# define SPLINES_VVISC
=======
>>>>>>> master
# define DJ_GRADPS
# define TS_U3HADVECTION
# define TS_C4VADVECTION
# define TS_DIF2
# define MIX_S_TS
# undef  MIX_GEO_TS
# define SALINITY
# define NONLIN_EOS
# define CURVGRID
# define MASKING
# define SOLVE3D
<<<<<<< HEAD
=======
# define SPLINES
>>>>>>> master
# define AVERAGES
# define AVERAGES_QUADRATIC
# define DIAGNOSTICS_UV
# define DIAGNOSTICS_TS

<<<<<<< HEAD
=======
# undef  MY25_MIXING
>>>>>>> master
# define LMD_MIXING
# ifdef LMD_MIXING
#  define SOLAR_SOURCE
#  define LMD_RIMIX
#  define LMD_CONVEC
#  define LMD_NONLOCAL
#  define LMD_SKPP
<<<<<<< HEAD
#  define RI_SPLINES
=======
>>>>>>> master
# endif

# define ANA_BSFLUX
# define ANA_BTFLUX

<<<<<<< HEAD
=======
# undef  M2CLIMATOLOGY
# undef  M3CLIMATOLOGY
# undef  TCLIMATOLOGY
# undef  M2CLM_NUDGING
# undef  M3CLM_NUDGING
# undef  TCLM_NUDGING

>>>>>>> master
# define FORWARD_MIXING
# define FORWARD_READ
# define FORWARD_WRITE
# define OUT_DOUBLE

#else

/*
**-----------------------------------------------------------------------------
**  Adjoint-based drivers.
**-----------------------------------------------------------------------------
*/

# define UV_ADV
# define UV_COR
# define UV_LDRAG
# define UV_VIS2
# define MIX_S_UV
<<<<<<< HEAD
# define SPLINES_VDIFF
# define SPLINES_VVISC
=======
>>>>>>> master
# define DJ_GRADPS
# define TS_U3HADVECTION
# define TS_C4VADVECTION
# define TS_DIF2
# define MIX_S_TS
# undef  MIX_GEO_TS
# define SALINITY
# define NONLIN_EOS
# define CURVGRID
# define MASKING
# define SOLVE3D
<<<<<<< HEAD

=======
# define SPLINES

# undef  MY25_MIXING
>>>>>>> master
# define LMD_MIXING
# ifdef LMD_MIXING
#  define SOLAR_SOURCE
#  define LMD_RIMIX
#  define LMD_CONVEC
#  define LMD_NONLOCAL
#  define LMD_SKPP
<<<<<<< HEAD
#  define RI_SPLINES
# endif

# define ANA_BSFLUX
# define ANA_BTFLUX

=======
# endif

# undef CLOSED_OBC
# ifdef CLOSED_OBC
#  define NORTHERN_WALL
#  define SOUTHERN_WALL
#  define EASTERN_WALL
#  define WESTERN_WALL
# else
#  define EASTERN_WALL
#  define WEST_FSCHAPMAN
#  define WEST_M2FLATHER
#  define WEST_M3CLAMPED
#  define WEST_TCLAMPED
#  define SOUTH_FSCHAPMAN
#  define SOUTH_M2FLATHER
#  define SOUTH_M3CLAMPED
#  define SOUTH_TCLAMPED
#  define NORTH_FSCHAPMAN
#  define NORTH_M2FLATHER
#  define NORTH_M3CLAMPED
#  define NORTH_TCLAMPED
# endif
# define ANA_BSFLUX
# define ANA_BTFLUX

# undef  M2CLIMATOLOGY
# undef  M3CLIMATOLOGY
# undef  TCLIMATOLOGY
# undef  M2CLM_NUDGING
# undef  M3CLM_NUDGING
# undef  TCLM_NUDGING

>>>>>>> master
# if defined W4DPSAS || defined W4DVAR
#  define CONVOLVE
# endif

# define VCONVOLUTION
# define IMPLICIT_VCONV
# undef  FULL_GRID

# define FORWARD_MIXING
# define FORWARD_READ
# define FORWARD_WRITE
# define OUT_DOUBLE
#endif
