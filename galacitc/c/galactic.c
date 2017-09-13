#include <math.h>
#include <stdio.h>

#define PI 3.14159265358979323846
#define deg2rad (PI / 180.0)
#define rad2deg (180.0 / PI)
#define gp_ra 192.859508
#define gp_dec 27.128336
#define lp 122.932
//#define sin_gp_dec (sin(gp_dec * deg2rad))
//#define cos_gp_dec (sqrt(1.0_DP - sin_gp_dec * sin_gp_dec))
#define sin_gp_dec 0.45598511203168229
#define cos_gp_dec 0.88998740305998381

void galactic1(double ra, double dec, double *x, double *y)
{
	double l, b;
	double d_ra, l0;
	double sin_b, cos_b;
	double sin_dec, cos_dec;
	double sin_l0, cos_l0;

	d_ra = (ra - gp_ra) * deg2rad;
	sin_dec = sin(dec * deg2rad);
	cos_dec = sqrt(1.0 - sin_dec * sin_dec);
	sin_b = sin_gp_dec * sin_dec + cos_gp_dec * cos_dec * cos(d_ra);
	b = asin(sin_b) * rad2deg;
	cos_b = sqrt(1.0 - sin_b * sin_b);

	sin_l0 = cos_dec * sin(d_ra) / cos_b;
	cos_l0 = (sin_dec - sin_gp_dec * sin_b) / (cos_gp_dec * cos_b);
	l0 = atan2(sin_l0, cos_l0) * rad2deg;
	l  = lp - l0;
	if (l < 0) l = l + 360;
	*x = l;
	*y = b;
	return;
}

void galactic2(double l, double b, double *x, double *y)
{
	double ra, dec;
	double d_l, ra0;
	double sin_b, cos_b;
	double sin_dec, cos_dec;
	double sin_ra0, cos_ra0;

	sin_b = sin(b * deg2rad);
	cos_b = sqrt(1.0 - sin_b * sin_b);
	d_l = (lp - l) * deg2rad;
	sin_dec = sin_b * sin_gp_dec + cos_b * cos_gp_dec * cos(d_l);
	dec = asin(sin_dec) * rad2deg;
	cos_dec = sqrt(1.0 - sin_dec * sin_dec);
	sin_ra0 = cos_b * sin(d_l) / cos_dec;
	cos_ra0 = (sin_b - sin_gp_dec * sin_dec) / (cos_gp_dec * cos_dec);
	ra0 = atan2(sin_ra0, cos_ra0) * rad2deg;
	ra  = ra0 + gp_ra;
	*x = ra;
	*y = dec;
	return;
}

void galactic0(int j, double lon1, double lat1, double *lon2, double *lat2)
{
	if (j == 1) {
		galactic1(lon1, lat1, lon2, lat2);
	} else if (j == 2) {
		galactic2(lon1, lat1, lon2, lat2);
	} else {
		fprintf(stderr,"warning:\n");
		fprintf(stderr,"j == 1, convert (ra, dec) to (l, b)\n");
		fprintf(stderr,"j == 2, convert (l, b) to (ra, dec)\n");
	}
	return;
}
