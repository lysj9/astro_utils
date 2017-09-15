Module galactic
	implicit none
	integer,parameter::DP=kind(1.0d0)
	real(kind=DP),parameter::PI=3.14159265358979323846_DP
	real(kind=DP),parameter::deg2rad=PI/180.0_DP
	real(kind=DP),parameter::rad2deg=180.0_DP/PI
!	gp_ra, gp_dec, lp: given by Reid & Brunthaler 2004
!	http://adsabs.harvard.edu/abs/2004ApJ...616..872R
!	real(kind=DP),parameter::gc_ra=266.405100_DP,gc_dec=-28.936175_DP
	real(kind=DP),parameter::gp_ra=192.859508_DP,gp_dec=27.128336_DP
	real(kind=DP),parameter::lp=122.932_DP
!	real(kind=DP),parameter::sin_gp_dec = sin(gp_dec * deg2rad)
!	real(kind=DP),parameter::cos_gp_dec = sqrt(1.0_DP - sin_gp_dec * sin_gp_dec)
	real(kind=DP),parameter::sin_gp_dec=0.45598511203168229_DP
	real(kind=DP),parameter::cos_gp_dec=0.88998740305998381_DP

	contains
	subroutine galactic1(ra,dec,l,b)
	implicit none
	real(kind=DP),intent(in)::ra,dec
	real(kind=DP),intent(out)::l,b
	real(kind=DP)::d_ra,l0
	real(kind=DP)::sin_b,cos_b
	real(kind=DP)::sin_dec,cos_dec
	real(kind=DP)::sin_l0,cos_l0
!	real(kind=DP)::sin_gp_dec,cos_gp_dec
	
	d_ra = (ra - gp_ra) * deg2rad
	sin_dec = sin(dec * deg2rad)
	cos_dec = sqrt(1.0_DP - sin_dec * sin_dec) ! cos(dec/gc_dec) in range [0,1]
!	sin_gp_dec = sin(gp_dec * deg2rad)
!	cos_gp_dec = sqrt(1.0_DP - sin_gp_dec * sin_gp_dec)
	sin_b = sin_gp_dec * sin_dec + cos_gp_dec * cos_dec * cos(d_ra)
	b = asin(sin_b) * rad2deg
	cos_b = sqrt(1.0_DP - sin_b * sin_b)
	sin_l0 = cos_dec * sin(d_ra) / cos_b
	cos_l0 = (sin_dec - sin_gp_dec * sin_b) / (cos_gp_dec * cos_b)
	l0 = atan2(sin_l0, cos_l0) * rad2deg
	l  = lp - l0
	if (l < 0) l = l + 360
	return
	end subroutine galactic1

	subroutine galactic2(l,b,ra,dec)
	implicit none
	real(kind=DP),intent(in)::l,b
	real(kind=DP),intent(out)::ra,dec
	real(kind=DP)::d_l,ra0
	real(kind=DP)::sin_b,cos_b
	real(kind=DP)::sin_dec,cos_dec
	real(kind=DP)::sin_ra0,cos_ra0
!	real(kind=DP)::sin_gp_dec,cos_gp_dec
	
	sin_b = sin(b * deg2rad)
	cos_b = sqrt(1.0_DP - sin_b * sin_b)
!	sin_gp_dec = sin(gp_dec * deg2rad)
!	cos_gp_dec = sqrt(1.0_DP - sin_gp_dec * sin_gp_dec)
	d_l = (lp - l) * deg2rad
	sin_dec = sin_b * sin_gp_dec + cos_b * cos_gp_dec * cos(d_l)
	dec = asin(sin_dec) * rad2deg
	cos_dec = sqrt(1.0_DP - sin_dec * sin_dec)
	sin_ra0 = cos_b * sin(d_l) / cos_dec
	cos_ra0 = (sin_b - sin_gp_dec * sin_dec) / (cos_gp_dec * cos_dec)
	ra0 = atan2(sin_ra0, cos_ra0) * rad2deg
	ra  = ra0 + gp_ra
	return
	end subroutine galactic2

	subroutine galactic0(j, lon1, lat1, lon2, lat2)
	use ISO_FORTRAN_ENV
	integer,intent(in)::j
	real(kind=DP),intent(in)::lon1,lat1
	real(kind=DP),intent(out)::lon2, lat2
	if (j == 1) then
		call galactic1(lon1, lat1, lon2, lat2)
	else if (j == 2) then
		call galactic2(lon1, lat1, lon2, lat2)
	else
		write (ERROR_UNIT, *) "warning:"
		write (ERROR_UNIT, *) "j == 1, convert (ra, dec) to (l, b)"
		write (ERROR_UNIT, *) "j == 2, convert (l, b) to (ra, dec)"
	end if
	return
	end subroutine galactic0
end Module galactic
