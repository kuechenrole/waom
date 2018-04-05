/*
** svn $Id: icetest.h 1307 2008-01-10 00:22:36Z bgalton $
*******************************************************************************
** Copyright (c) 2002-2008 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for Simplified Ice Shelf Ocean Cavity Model Test.
**
** Application flag:   WAOM10
** Input script:       ocean_waom10.in
*/

/*********** momentum *************/

#define UV_COR
#define UV_VIS2
#define UV_QDRAG
#define UV_ADV

/*********** tracer ***************/

#define TS_A4HADVECTION
#define TS_A4VADVECTION
#define TS_DIF2
#define SALINITY
#define NONLIN_EOS

/*********** pressure *************/

#define DJ_GRADPS

/*********** model config *********/

#define SOLVE3D
#define CURVGRID
#define SPHERICAL

/*********** splines reconstruction ***********/

#define SPLINES_VDIFF
#define SPLINES_VVISC

/*********** ana fields ***********/

#define ANA_BSFLUX
#define ANA_BTFLUX
#define ANA_SRFLUX
/*
#define ANA_INITIAL
#define ANA_SMFLUX
#define ANA_STFLUX
#define ANA_SSFLUX
*/

/*********** mixing hor mom *******/

#define MIX_S_UV  #MIX_GEO_UV better said kait, but just works with visc4 for kaitlin


/*********** mixing hor tracer ****/

#define MIX_ISO_TS  #MIX_S_TS you had before


/*********** mixing vertical turbulent momentum and tracer ************/

#define LMD_MIXING

#ifdef LMD_MIXING
# define LMD_CONVEC
# define RI_SPLINES
# define LMD_DDMIX
# define LMD_RIMIX
# define LMD_SKPP
# define LMD_BKPP
# define LMD_NONLOCAL #;..try not big diff expected
# define LMD_SHAPIRO  #..optional ben instabilty probs
#endif



/*********** ice shelf options **********/

#define ICESHELF
#ifdef ICESHELF
# define LIMIT_ICESTRESS
# undef ANA_SEAICE
# define ICESHELF_3EQN_VBC
#endif

/*********** tides **********************/

/*SSH RED POT (Nodal included in tds file)*/
  
#define SSH_TIDES
#define ADD_FSOBC
#define UV_TIDES
#define ADD_M2OBC
#define RAMP_TIDES
/*
#define FSOBC_REDUCED
#define POT_TIDES
#define TIMELESS_DATA
*/

/*SSH RED POT NoNodal*/
/*
#define RAMP_TIDES
#define SSH_TIDES
#define FSOBC_REDUCED
#define ADD_M2OBC
#define ADD_FSOBC
*/

/*SSH UV POT NoNodal*/
/*
#define RAMP_TIDES
#define SSH_TIDES
#define UV_TIDES
#define ADD_M2OBC
#define ADD_FSOBC
#define POT_TIDES
#define TIMELESS_DATA
*/


/************ NETCDF Input output *************/

#define PERFECT_RESTART
#define AVERAGES

/************ other ************************/

#define MASKING
#undef WET_DRY
#define LIMIT_BSTRESS
#undef LIMIT_STFLX

#define SURFACE_OVERFLUX_FIX
#define QCORRECTION
#define SCORRECTION

