#ifndef _AMBPRO_HEDEAER_20150326_
#define _AMBPRO_HEDEAER_20150326_
#include "../Common/BaseFunction.h"
#include "../Common/Common.h"
#include "../Common/GNSSERROR.h"

#ifdef AMBPRO_EXPORTS
#define AMBPRO_API _declspec(dllexport)
#else
#define AMBPRO_API _declspec(dllimport)
#endif

AMBPRO_API   double probability(double ave, double sig);
AMBPRO_API  void Fix_EWL_Amb(const GNSSDATA &Obs, SSat_t *Ssat, double Azel[],char *Station);
#endif