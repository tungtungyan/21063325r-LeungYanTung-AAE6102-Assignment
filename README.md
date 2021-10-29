<p align="center">

<h3 align="center">AAE6102 Assignment Leung Yan Tung 21063325r due: Nov 1st, 2021</h3>

<p align="center">

This report uses the single-epoch data sets rcvr.dat and eph.dat to set up the linearized navigation equations and solve for user position and clock bias. Appendix A show the data file format of rcvr.dat and eph.dat.

The linearized navigation equation:
<p align="center">
Single pseudorange: ρ_i=√((x_i-x_u )^2+(y_i-y_u )^2+(z_i-z_u )^2 )+ct_u
x_u=x ̂_u+∆x_u
y_u=y ̂_u+∆y_u
z_u=z ̂_u+∆z_u
t_u=t ̂_u+∆t_u
<p align="center">

The required corrections for the satellite clock bias and relativity are referred to ICD. This report will skip the ionospheric corrections because we do not have access to the parameter values of the Klobuchar model for this data set.

Tropospheric correction based on standard atmosphere model is optional. 
The initial position to start the iteration using 
<br />
	
	[ −2694685.473 ] 

	[ −4293642.366 ] 

	[ 3857878 . 924 ] 
<br />
(WGS 84 XYZ, in meters). The algorithm uses zero clock bias for Initialize.  Terminate the iteration when the change in the estimate is suitably small.


GPS Constants:
Speed of light: c = 299792458.0 (m/s)
WGS 84 value of earth’s rotation rate: Wedot= 7.2921151467e-5 (r/s)
WGS 84 value of earth's universal gravitation constant: mu= 3.986005e+14 (m^3/s^2)
Relativistic correction term constant: F= -4.442807633e-10


Using the data of rcvr.dat and eph.dat, we can calculate the receiver’s position at time of week 440992 using the following process.

	Calculate the XYZ positions for all valid satellite at time 440992.


First, calculate the satellite Earth-centered, Earth-fixed coordinate system (ECEF) position vector.
Square root of semi-major axis a (Sqrta) was save in the Column 10 of eph.dat. 
Therefore, the semimajor axis (a) is

a=〖sqeta〗^2
The corrected mean motion(n):
n=n_0-∆n
Computed mean motion (n_0)
n_0=√(mu/a^3 )
Time from ephemeris epoch (t_k):
t_k=t-〖∆t〗_oe
t_oeis reference time of ephemeris parameters (s). 
Mean anomaly:
M_k=M_0+nt_k
Set the maximum number of iterations is 10.
Then, use Kepler’s equation of eccentric anomaly is solved by iteration.

E_k=E_0-e sin⁡〖E_k 〗

Calculate the satellite clock offset (t).

t=t_SV-〖∆t〗_SV
〖∆t〗_SV=a_0+a_1 (t-t_0c )+a_2 (t-t_0c )^2+〖∆t〗_r
t_SV: The individual satellite time.

	Determine the broadcast satellite clock error.
<br />
	
	Estimate the tropospheric delay for each satellite (optional).
<br />
		
	use the linerized GPS measurement equation developed in class to estimate the vector δxˆ
<br />
		
	update the estimate of the user position: X0 (new)= X0 (old)+δxˆ
<br />
		
	if δxˆ <10−4 m, then we have successfully converged on a valid position solution (Some of the MATLAB functions in the folder will be useful in solving this problem.)


What is your estimate of the user clock bias b? Does your estimate of the user clock bias in seconds offer insight as to why the reported receiver clock time at this epoch (Column 1 of the rcvr matrix) is 440992.00173454 seconds? (Hints: Your initial iteration should give (in meters) δx ≈ -5710, δy ≈ 1080, δz ≈ -2610, δb ≈ 519450. Your position estimate (WGS 84
−2700400
XYZ coordinates, in meters) should be x ≈ [■(-2700400@-4292560@3857878.924)]   ) 3855270

 

Appendix A
Data file format of rcvr.dat and eph.dat

Data File Format of rcvr.dat:

rcvr.dat is an 8x7 matrix containing raw ranging information. Each of the 8 rows contains independent measurements for each of 8 satellites in view at the current epoch (an epoch is simply a term refers to a single discrete time; since our receivers provide data at approximately 1 sec. intervals, each epoch occurs approximately 1 sec. after the prior epoch. The columns of this matrix in clued the following data:

Column	Symbol	Meaning
Column 1: 	rcvr_tow 	receiver time of week(s)
	<br />
Column 2:	Svid 	satellite PRN number (1 – 32)
	<br />
Column 3:	pr	pseudorange (m)
	<br />
Column 4:	cycles	number of accumulated cycles
	<br />
Column 5:	Phase	to convert to (0 – 359.99) mult. by 360/2048
	<br />
Column 6:	slp_dtct	0 = no cycle slip detected; non 0 = cycle slip
	<br />
Column 7:	snr_dbhz	signal to noise ratio (dB-Hz)
	<br />


Data File Format of eph.dat:

eph.dat is a 8 x 24 matrix containing the ephemeris data from a GPS receiver. This data is used to estimate the orbital position of each satellite at any given time. Each row contains ephemeris data for a single satellite. The columns of this matrix include the following data:


Column	Symbol	Meaning
Column 1: 	rcvr_tow 	receiver time of week(s)
	<br />  
Column 2:	Svid 	satellite PRN number (1 – 32)
	 <br />  
Column 3:	toc  	reference time of clock parameters (s)
	 <br />  
Column 4:	toe 	reference time of ephemeris parameters (s)
	<br />
Column 5:	af0	clock correction coefficient – group delay (s) 
	<br />
Column 6:	af1	clock correction coefficient (s/s)
	<br />
Column 7:	af2	clock correction coefficient (s/s/s) 
	<br />
Column 8:	ura	user range accuracy (m)
	<br />
Column 9:	e	eccentricity (-)
	<br />
Column 10:	sqrta	square root of semi-major axis a (m**1/2)
	<br />
Column 11:	dn	mean motion correction (r/s)
	<br />
Column 12:	m0	mean anomaly at reference time (r)
	<br />
Column 13:	w	argument of perigee (r)
	<br />
Column 14:	omg0	right ascension (r)
	<br />
Column 15:	i0	inclination angle at reference time (r)
	<br />
Column 16:	odot	rate of right ascension (r/s)
	<br />
Column 17:	idot	rate of inclination angle (r/s)
	<br />
Column 18:	cus 	argument of latitude correction, sine (r)
	<br />
Column 19:	cuc	argument of latitude correction, cosine (r)
	<br />
Column 20:	cis	inclination correction, sine (r)
	<br />
Column 21:	cic 	inclination correction, cosine (r)
	<br />
Column 22:	crs	radius correction, sine (m)
	<br />
Column 23:	crc 	radius correction, cosine (m)
	<br />
Column 24:	iod	issue of data number
	<br />


