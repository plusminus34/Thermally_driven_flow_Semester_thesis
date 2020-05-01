#pragma once

#include <cmath>

const double TO_RAD = 3.141592653589793 / 180.0;//shouldn't that exist somewhere?
const double POLE_PHI = 90.0;//mNorthPole.Lat;
const double POLE_LAMBDA = 0;//mNorthPole.Lon;
const double ZRPI18 = 57.2957795;
const double ZPIR18 = 0.0174532925;
const double SINPOL = sin(ZPIR18*POLE_PHI * TO_RAD);
const double COSPOL = cos(ZPIR18*POLE_PHI * TO_RAD);

void RlatRlonToLatLon(double rlat, double rlon, double& lat, double& lon) {
	double zphis = ZPIR18 * rlat;
	double zlams = rlon > 180 ? rlon - 360 : rlon;
	zlams *= ZPIR18;
	double arg = COSPOL * cos(zphis*TO_RAD)*cos(zlams*TO_RAD) + SINPOL * sin(zphis*TO_RAD);

	// lat
	lat = ZRPI18 * asin(arg) / TO_RAD;

	// lon
	double ZLAMPOL = ZPIR18 * POLE_LAMBDA;
	const double tmp1 = sin(ZLAMPOL*TO_RAD) * (-SINPOL * cos(zlams*TO_RAD) * cos(zphis*TO_RAD) + COSPOL * sin(zphis*TO_RAD));
	const double tmp2 = cos(ZLAMPOL*TO_RAD) * sin(zlams*TO_RAD) * cos(zphis*TO_RAD);
	double ZARG1 = tmp1 - tmp2;
	double ZARG2 = tmp1 + tmp2;
	lon = 0;
	if (std::abs(ZARG2) < 1E-30)
	{
		if (std::abs(ZARG1) < 1E-30)
			lon = 0.0;
		else if (ZARG1 < 0.)
			lon = 90.0;
		else
			lon = -90.0;
	}
	else
		lon = ZRPI18 * atan2(ZARG1, ZARG2) / TO_RAD;
}

void LatLonToRlanRlon(double lat, double lon, double& rlat, double& rlon) {
	double ZLAMPOL = ZPIR18 * POLE_LAMBDA;
	double ZPHI = ZPIR18 * lat;
	double ZLAM = lon;
	if (ZLAM > 180.0)
		ZLAM = ZLAM - 360.0;
	ZLAM = ZPIR18 * ZLAM;

	// lat
	double ZARG = COSPOL * cos(ZPHI*TO_RAD) * cos((ZLAM - ZLAMPOL) *TO_RAD) + SINPOL * sin(ZPHI*TO_RAD);
	rlat = ZRPI18 * asin(ZARG) / TO_RAD;


	// lon
	double ZARG1 = -sin((ZLAM - ZLAMPOL)*TO_RAD) * cos(ZPHI*TO_RAD);
	double ZARG2 = -SINPOL * cos(ZPHI*TO_RAD) * cos((ZLAM - ZLAMPOL)*TO_RAD) + COSPOL * sin(ZPHI*TO_RAD);
	if (std::abs(ZARG2) < 1E-30)
	{
		if (std::abs(ZARG1) < 1E-30)
		{
			rlon = 0.0;
		}
		else
		{
			if (ZARG1 > 0)
				rlon = 90.0;
			else
				rlon = -90.0;
		}
	}
	else
		rlon = ZRPI18 * atan2(ZARG1, ZARG2) / TO_RAD;
}

void degreeLengths(double lat, double& lat_length, double& lon_length) {
	//From https://en.wikipedia.org/wiki/Geographic_coordinate_system
	//Length of a degree in meters
	lat_length = 111132.92 - 559.82 * cos(2 * lat*TO_RAD) + 1.175 * cos(4 * lat*TO_RAD) - 0.0023*cos(6 * lat*TO_RAD);
	lon_length = 111412.84 * cos(lat*TO_RAD) - 93.5 * cos(3 * lat*TO_RAD) + 0.118*cos(5 * lat*TO_RAD);
}

