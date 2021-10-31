<p align="center">

<h3 align="center">AAE6102 Assignment <br />Leung Yan Tung 21063325r <br /> </h3>

<p align="center">


<!-- Introdution -->
# Introdution
<details open="open">
This report uses the single-epoch data sets rcvr.dat and eph.dat to set up the linearized navigation equations and solve for user position and clock bias. Appendix A show the data file format of rcvr.dat and eph.dat.

	
The linearized navigation equation:
<p align="center">
<img width="473" alt="Screenshot 2021-10-31 at 10 42 10 AM" src="https://user-images.githubusercontent.com/71690213/139564718-86bf7ebf-0037-4fae-8ddd-c58975619690.png">
<p align="center">

The required corrections for the satellite clock bias and relativity are referred to ICD. This report will skip the ionospheric corrections because we do not have access to the parameter values of the Klobuchar model for this data set.

The required corrections for the satellite clock bias and relativity are referred to ICD. This report will skip the ionospheric corrections because we do not have access to the parameter values of the Klobuchar model for this data set. This report will consider the tropospheric correction based on standard atmosphere model. 


	[ −2694685.473 ] 
	[ −4293642.366 ]   (WGS 84 XYZ, in meters)
	[ 3857878 .924 ] 

	
The algorithm uses zero clock bias for Initialize. Terminate the iteration when the change in the estimate is suitably small.
</details>

# Constants
<details open="open">
	
## GPS Constants:
Speed of light:                                         c = 299792458.0 (m/s)
<br />WGS 84 value of earth’s rotation rate:                  Wedot= 7.2921151467e-5 (r/s)
<br />WGS 84 value of earth's universal gravitation constant: mu= 3.986005e+14 (m^3/s^2)
<br />Relativistic correction term constant:                  F= -4.442807633e-10

## Tropospheric constant:
Pressure: 		  Pr = 1013.25(mbar)
<br />Temperature: 		  Tr = 291.15(K)
<br />Temperature at sea-level: temp_0 = 15
</details>
	
# Calculation of processing
<details open="open">
Using the data of rcvr.dat and eph.dat, it can calculate the receiver’s position at time of week 440992 using the following process.

## 1. Calculate the XYZ positions for all valid satellite at time 440992.

<details open="open">
Calculate the satellite Earth-centered, Earth-fixed coordinate system (ECEF) position vector.
<br />
Square root of semi-major axis a (Sqrta) was saved in the Column 10 of eph.dat.
<br />
Therefore, the semimajor axis (a) is
<br />
<p align="center">
<img width="128" alt="Screenshot 2021-10-29 at 4 50 39 PM" src="https://user-images.githubusercontent.com/71690213/139407961-a6c07131-4723-4909-a11d-1709382fe2d7.png">
<p align="center">
	
<br />The corrected mean motion(n):
<p align="center">
<img width="127" alt="Screenshot 2021-10-29 at 4 50 45 PM" src="https://user-images.githubusercontent.com/71690213/139408058-e0d4f313-77b9-4fb8-8db9-efe4133368d0.png">
<p align="center">
	
<br />Computed mean motion (n_0)
<p align="center">
<img width="146" alt="Screenshot 2021-10-29 at 4 50 55 PM" src="https://user-images.githubusercontent.com/71690213/139408090-f0b05a1b-6199-4a8f-b3de-87be78e7378d.png">
<br /><img width="146" alt="Screenshot 2021-10-29 at 4 50 55 PM" src="https://user-images.githubusercontent.com/71690213/139408228-396bd7d3-b34c-41e5-b73f-be06245e6d24.png">
<p align="center">

<br />Time from ephemeris epoch (t_k):
<p align="center">
<img width="135" alt="Screenshot 2021-10-29 at 4 51 02 PM" src="https://user-images.githubusercontent.com/71690213/139408243-021110ba-8c6a-48b9-b6f0-ae97f42bc160.png">
<p align="center">	
Where <img width="28" alt="Screenshot 2021-10-31 at 11 20 26 AM" src="https://user-images.githubusercontent.com/71690213/139565825-60844995-0bd7-43c7-a8ef-854db1cf38d5.png">
 reference time of ephemeris parameters (s). 

<br />Mean anomaly:
<p align="center">
<img width="153" alt="Screenshot 2021-10-29 at 4 51 07 PM" src="https://user-images.githubusercontent.com/71690213/139408255-0fd371b8-f0d8-427d-9f44-58128baf3569.png">
<p align="center">

<br />Assume that the maximum number of iterations is 10.
<br />Then, using Kepler’s equation of eccentric anomaly is solved by iteration.
<p align="center">
<img width="179" alt="Screenshot 2021-10-29 at 4 51 17 PM" src="https://user-images.githubusercontent.com/71690213/139408278-2c561019-8690-4299-8998-1e376c14742b.png">
<p align="center">
	
<br />Calculate the true anomaly:
<p align="center">
<img width="289" alt="Screenshot 2021-10-29 at 4 51 27 PM" src="https://user-images.githubusercontent.com/71690213/139410511-9387819f-9c11-454c-a1c0-7506a31485d2.png">
<p align="center">
	
<br />Argument of latitude:
<p align="center">
<img width="132" alt="Screenshot 2021-10-29 at 4 51 37 PM" src="https://user-images.githubusercontent.com/71690213/139410530-44e28dc8-b853-4c0f-bbea-2e2b3d708fc5.png">
<p align="center">

<br />Use the Second Harmonic Perturb to find to Correction value:
<br />
Argument of latitude Correction:
<p align="center">
<img width="306" alt="Screenshot 2021-10-29 at 4 51 47 PM" src="https://user-images.githubusercontent.com/71690213/139410551-83e449c2-f57c-42c2-851d-bdb240e07a27.png">
<p align="center">


<br />Radius Correction:
<p align="center">
<img width="295" alt="Screenshot 2021-10-29 at 4 51 58 PM" src="https://user-images.githubusercontent.com/71690213/139410597-09dc2425-67c4-4b83-afd7-d752548ca7b1.png">
<p align="center">


<br />Inclination correction:
<p align="center">
<img width="290" alt="Screenshot 2021-10-29 at 4 52 06 PM" src="https://user-images.githubusercontent.com/71690213/139410625-c2c9f946-fae7-4943-8b0c-02d8c4fc6571.png">
<p align="center">
	
<br />Use the Correction value for Corrected:
<br />
Corrected argument of latitude:
<p align="center"><img width="124" alt="Screenshot 2021-10-31 at 11 27 57 AM" src="https://user-images.githubusercontent.com/71690213/139566032-fba07064-13a9-435c-8a21-fb8fee5a09f8.png">
<p align="center">


<br />Corrected radius:
<br />
<p align="center"><img width="251" alt="Screenshot 2021-10-29 at 4 52 17 PM" src="https://user-images.githubusercontent.com/71690213/139415014-db907346-4e22-4e3c-b205-fb6b818397c7.png">
<p align="center">

<br />
Corrected Inclination:
<p align="center">
<img width="243" alt="Screenshot 2021-10-29 at 4 52 24 PM" src="https://user-images.githubusercontent.com/71690213/139415027-27acc67c-78bc-40d2-8b9b-ee019ef1c91d.png">
<p align="center">
<br /> 
Where IDOT in the Column 17 of eph.dat.

<br />Corrected longitude of ascending node:
<p align="center">
<img width="287" alt="Screenshot 2021-10-29 at 4 52 31 PM" src="https://user-images.githubusercontent.com/71690213/139415057-0bcf22e6-d19a-4354-be26-5ee00f066a82.png">
<p align="center">
	
<br />Position in orbital plane:
<br />x-coordinate:
<p align="center"><img width="155" alt="Screenshot 2021-10-29 at 4 52 37 PM" src="https://user-images.githubusercontent.com/71690213/139415344-f630aea4-a16f-4805-97e4-d1ec60b12e40.png">
<p align="center">
<br />y-coordinate:
<img width="140" alt="Screenshot 2021-10-29 at 4 52 45 PM" src="https://user-images.githubusercontent.com/71690213/139415355-617f71d2-a0df-49fb-b9a6-d93ff29dd427.png">


Therefor,
The Earth-fixed coordinates calculate by 

<img width="345" alt="Screenshot 2021-10-29 at 4 54 00 PM" src="https://user-images.githubusercontent.com/71690213/139415371-a008871b-474e-4702-b3d0-4c7bf778eb7c.png">


However, the coordinates have some rotation error when transformed the coordinate from ECEF coordinate to Earth-Centered, Inertial (ECI) coordinate system.

The rotational matrix for coordinate transform correction:
</details>

	Determine the broadcast satellite clock error.
<br />Calculate the satellite clock offset (t).

<img width="452" alt="Screenshot 2021-10-29 at 4 55 01 PM" src="https://user-images.githubusercontent.com/71690213/139415999-7a3c8b68-b4ef-43bc-bac2-5f1e7881669c.png">

t_SV: The individual satellite time.

	
###### Estimate the tropospheric delay for each satellite (optional).
<br />First, resolve the tropospheric error of each iteration. 
Using iterative algorithm to transform the initial position of ECEF coordinates to Geodetic coordinates (ϕ,λ,h).
ECEF coordinates=
<img width="169" alt="Screenshot 2021-10-29 at 4 57 13 PM" src="https://user-images.githubusercontent.com/71690213/139416027-c2b1ea05-8dd0-4024-b263-4050cb7e06a0.png">

Then, Transform the Geodetic coordinates to east, north, up (ENU) coordinates.


Using Saastamoinen model to enale tropospheric error estimation.
<img width="261" alt="Screenshot 2021-10-29 at 4 57 35 PM" src="https://user-images.githubusercontent.com/71690213/139416037-0ec706c5-c046-456f-9a2e-62625d2e020a.png">

Tropospheric delay include two part which are the hydrostatic delay (trph) and wet delay (trpw).

<img width="688" alt="Screenshot 2021-10-29 at 4 58 12 PM" src="https://user-images.githubusercontent.com/71690213/139416048-80ff8d26-9d48-433f-8e26-07de28f74ae2.png">

<img width="459" alt="Screenshot 2021-10-29 at 4 58 34 PM" src="https://user-images.githubusercontent.com/71690213/139416067-569f1a56-79e4-474d-abf0-21d063c15d13.png">

where p, T and H are pressure, temperature and humidity at height h. 
where e is the partial water vapor pressure. 

		
	use the linerized GPS measurement equation developed in class to estimate the vector δxˆ
<br />
		
	update the estimate of the user position: X0 (new)= X0 (old)+δxˆ
<br />
		
	if δxˆ <10−4 m, then we have successfully converged on a valid position solution (Some of the MATLAB functions in the folder will be useful in solving this problem.)


What is your estimate of the user clock bias b? Does your estimate of the user clock bias in seconds offer insight as to why the reported receiver clock time at this epoch (Column 1 of the rcvr matrix) is 440992.00173454 seconds? (Hints: Your initial iteration should give (in meters) δx ≈ -5710, δy ≈ 1080, δz ≈ -2610, δb ≈ 519450. Your position estimate (WGS 84
−2700400
XYZ coordinates, in meters) should be x ≈ [■(-2700400@-4292560@3857878.924)]   ) 385527
</details> 
 

Appendix A
		     <br />
<strong>[eph.dat](https://github.com/tungtungyan/21063325r-LeungYanTung-AAE6102-Assignment/blob/main/eph.dat)</strong>
<strong>[rcvr.dat](https://github.com/tungtungyan/21063325r-LeungYanTung-AAE6102-Assignment/blob/main/rcvr.dat)</strong>
