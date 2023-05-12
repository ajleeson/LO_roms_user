/*
** svn $Id: bio_toy.h 1001 2020-01-10 22:41:16Z arango $
*******************************************************************************
** Copyright (c) 2002-2020 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Based on bio_toy.h
**
** Options 1-dimensional vertical mixing example for course
** Coastal Ocean Dynamics (16:712:503) John Wilkin and Bob Chant
**
** Application flag: IDEAL_2PER
** Input script: roms_V-less.in
*/

#undef ANA_PSOURCE
#define ANA_GRID
#define ANA_INITIAL

#define UV_VIS2
#define MIX_S_UV
#define TS_DIF2
#define MIX_S_TS

#define UV_ADV
#define UV_COR
#define UV_QDRAG
#define DJ_GRADPS
#define SPLINES_VDIFF
#define SPLINES_VVISC
#define SOLAR_SOURCE
#undef  NONLIN_EOS
#define SALINITY
#define AVERAGES
#define AVERAGES_AKV        /* for COD */
#define AVERAGES_AKT        /* for COD */
#define SOLVE3D

/*                      output options */
#define DIAGNOSTICS_TS
#define DIAGNOSTICS_UV

/*  hadvection options are set roms.in */

/* vertical turbulence closure options */
#undef  LMD_MIXING
#define GLS_MIXING          /* for COD */

#ifdef LMD_MIXING
# define LMD_RIMIX
# define LMD_CONVEC
# define LMD_SKPP
# define LMD_BKPP
# define LMD_NONLOCAL
# define RI_SPLINES
#endif

#ifdef GLS_MIXING           /* for COD */
# define N2S2_HORAVG
# define RI_SPLINES
# undef  CANUTO_A
# undef  CANUTO_B
# define KANTHA_CLAYSON
# undef  CHARNOK
# undef  ZOS_HSIG
# undef  CRAIG_BANNER
# undef  TKE_WAVEDISS
#endif

#undef BULK_FLUXES   
#ifdef BULK_FLUXES
# define ANA_PAIR
# define ANA_HUMIDITY
# define ANA_TAIR
# define ANA_SRFLUX
# define LONGWAVE
# define ANA_RAIN
# define ANA_CLOUD
# define ANA_WINDS
#else
/* for COD press grad periodic forcing */
# undef  ATM_PRESS  /* grad(p) forcing for COD_1DMIX */
# undef  ANA_PAIR  /* periodic forcing for COD_1DMIX */
# undef  RAMP_TIDES
# define ANA_SRFLUX
# define ANA_SMFLUX
# define ANA_STFLUX
#endif

#define ANA_SSFLUX
#define ANA_BSFLUX
#define ANA_BTFLUX
#undef  ANA_WWAVE           /* for COD */

/* to use volume-less vertical sources */
#define LWSRC_MASS_ONLY


