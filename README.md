<p align="center">

<h3 align="center">AAE6102 Assignment <br />Leung Yan Tung 21063325r <br /> due: Nov 1st, 2021</h3>

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
<br />Speed of light: c = 299792458.0 (m/s)
<br />WGS 84 value of earth’s rotation rate: Wedot= 7.2921151467e-5 (r/s)
<br />WGS 84 value of earth's universal gravitation constant: mu= 3.986005e+14 (m^3/s^2)
<br />Relativistic correction term constant: F= -4.442807633e-10

Using the data of rcvr.dat and eph.dat, we can calculate the receiver’s position at time of week 440992 using the following process.

	Calculate the XYZ positions for all valid satellite at time 440992.


First, calculate the satellite Earth-centered, Earth-fixed coordinate system (ECEF) position vector.
Square root of semi-major axis a (Sqrta) was save in the Column 10 of eph.dat. 
Therefore, the semimajor axis (a) is a=〖sqeta〗^2<img width="128" alt="Screenshot 2021-10-29 at 4 50 39 PM" src="https://user-images.githubusercontent.com/71690213/139407961-a6c07131-4723-4909-a11d-1709382fe2d7.png">

<br />The corrected mean motion(n):
<br />n=n_0-∆n<img width="127" alt="Screenshot 2021-10-29 at 4 50 45 PM" src="https://user-images.githubusercontent.com/71690213/139408058-e0d4f313-77b9-4fb8-8db9-efe4133368d0.png">

<br />Computed mean motion (n_0)
<br />n_0=√(mu/a^3 )<img width="146" alt="Screenshot 2021-10-29 at 4 50 55 PM" src="https://user-images.githubusercontent.com/71690213/139408090-f0b05a1b-6199-4a8f-b3de-87be78e7378d.png">

<br />n_0=√(mu/a^3 )<img width="146" alt="Screenshot 2021-10-29 at 4 50 55 PM" src="https://user-images.githubusercontent.com/71690213/139408228-396bd7d3-b34c-41e5-b73f-be06245e6d24.png">

<br />Time from ephemeris epoch (t_k):
<br />t_k=t-〖∆t〗_oe<img width="135" alt="Screenshot 2021-10-29 at 4 51 02 PM" src="https://user-images.githubusercontent.com/71690213/139408243-021110ba-8c6a-48b9-b6f0-ae97f42bc160.png">

<br />t_oeis reference time of ephemeris parameters (s). 
<br />Mean anomaly:
<br />M_k=M_0+nt_k<img width="153" alt="Screenshot 2021-10-29 at 4 51 07 PM" src="https://user-images.githubusercontent.com/71690213/139408255-0fd371b8-f0d8-427d-9f44-58128baf3569.png">

<br />Set the maximum number of iterations is 10.
<br />Then, use Kepler’s equation of eccentric anomaly is solved by iteration.

E_k=E_0-e sin⁡〖E_k 〗<img width="179" alt="Screenshot 2021-10-29 at 4 51 17 PM" src="https://user-images.githubusercontent.com/71690213/139408278-2c561019-8690-4299-8998-1e376c14742b.png">


Calculate the satellite clock offset (t).

<br />t=t_SV-〖∆t〗_SV
<br />〖∆t〗_SV=a_0+a_1 (t-t_0c )+a_2 (t-t_0c )^2+〖∆t〗_r
	<img width="289" alt="Screenshot 2021-10-29 at 4 51 27 PM" src="https://user-images.githubusercontent.com/71690213/139408321-b6fe2752-a780-4a30-a99e-a87148305b09.png">

<br />t_SV: The individual satellite time.

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
		     <br />
<strong>[eph.dat](https://github.com/tungtungyan/21063325r-LeungYanTung-AAE6102-Assignment/blob/main/eph.dat)</strong>
<strong>[rcvr.dat](https://github.com/tungtungyan/21063325r-LeungYanTung-AAE6102-Assignment/blob/main/rcvr.dat)</strong>
