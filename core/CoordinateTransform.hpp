#pragma once

#include <cmath>
#include <iostream>

// TODO obtain rotated north pole from .nc file?
const double POLE_PHI = 43;// latitude of rotated north pole
const double POLE_LAMBDA = -170;// longitude of rotated noth pole
const double ZRPI18 = 180.0 / 3.141592653589793;//   57.2957795;// 180/PI
const double ZPIR18 = 1 / ZRPI18;//0.0174532925;// PI/180
const double SINPOL = sin(ZPIR18*POLE_PHI);
const double COSPOL = cos(ZPIR18*POLE_PHI);

class CoordinateTransform {
public:

	static void RlatRlonToLatLon(double rlat, double rlon, double& lat, double& lon) {
		double zphis = ZPIR18 * rlat;
		double zlams = rlon > 180 ? rlon - 360 : rlon;
		zlams *= ZPIR18;
		double arg = COSPOL * cos(zphis) * cos(zlams) + SINPOL * sin(zphis);

		// lat
		lat = ZRPI18 * asin(arg);

		// lon
		double ZLAMPOL = ZPIR18 * POLE_LAMBDA;
		const double tmpa = -SINPOL * cos(zlams) + COSPOL * sin(zphis);
		const double tmpb = sin(zlams)*cos(zphis);
		double ZARG1 = sin(ZLAMPOL) * tmpa - cos(ZLAMPOL)*tmpb;
		double ZARG2 = cos(ZLAMPOL) * tmpa - sin(ZLAMPOL)*tmpb;
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
			lon = ZRPI18 * atan2(ZARG1, ZARG2);
	}

	static void LatLonToRlanRlon(double lat, double lon, double& rlat, double& rlon) {
		double ZLAMPOL = ZPIR18 * POLE_LAMBDA;
		double ZPHI = ZPIR18 * lat;
		double ZLAM = lon;
		if (ZLAM > 180.0)
			ZLAM = ZLAM - 360.0;
		ZLAM = ZPIR18 * ZLAM;

		// lat
		double ZARG = COSPOL * cos(ZPHI) * cos((ZLAM - ZLAMPOL)) + SINPOL * sin(ZPHI);
		rlat = ZRPI18 * asin(ZARG);


		// lon
		double ZARG1 = -sin((ZLAM - ZLAMPOL)) * cos(ZPHI);
		double ZARG2 = -SINPOL * cos(ZPHI) * cos((ZLAM - ZLAMPOL)) + COSPOL * sin(ZPHI);
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
			rlon = ZRPI18 * atan2(ZARG1, ZARG2);
	}

	static void degreeLengths(double lat, double& lat_length, double& lon_length) {
		//From https://en.wikipedia.org/wiki/Geographic_coordinate_system
		//Length of a degree in meters
		lat_length = 111132.92 - 559.82 * cos(2 * lat) + 1.175 * cos(4 * lat) - 0.0023*cos(6 * lat);
		lon_length = 111412.84 * cos(lat) - 93.5 * cos(3 * lat) + 0.118*cos(5 * lat);
	}

	static void degreeLengthsSimple(double lat, double& lat_length, double& lon_length) {
		lat_length = 6367449;//*cos(0)
		lon_length = lat_length * cos(lat);
	}
};