#ifndef _RTCLKPROCESS_HEADER_20150603_
#define _RTCLKPROCESS_HEADER_20150603_

#include "../GNSSProcess/SPPProcess.h"


GNSSPROCESS_API    void RTCLKInitial(const GNSSDATA Obs,const double SatXYZ[], const double SatClk[], double SatAzEl[], double ephRMS[]);


#endif