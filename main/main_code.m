%AAE6102 Assignment 1 
%Leung Yan Tung 21063325r

clc;
clear;
close all;
warning('off', 'all');

%% Settings
ENA_TROPO_ERR_CORR = true; % identify to enable tropospheric error estimation using Saastamoinen model

%% Read given data
% ephemeris: rcvr_tow, svid, toc, toe, af0, af1, af2, ura, e, sqrta, dn, m0, w, omg0, i0, odot, idot, cus, cuc, cis, cic, crs, crc, iod
eph = importdata(fullfile('Data','eph.dat')); %import ephemeris data
eph = sortrows(eph,[1, 2]);% sorts the rows of eph in ascending order based on the rcvr_tow and svid

% receiver: rcvr_tow, svid, pr, cycle, phase,cycle slp, snr
rcvr = importdata(fullfile('Data','rcvr.dat')); %import receiver data
rcvr = sortrows(rcvr,[1, 2]);% sorts the rows of rcvr in ascending order based on the rcvr_tow and svid

% extract the all the value of the seconde column
eph_c2=eph(:,2); 
rcvr_c2=rcvr(:,2);

% select common satellite
[~, eph_ord] = ismember(eph_c2, rcvr_c2);
[~, rcvr_ord] = ismember(rcvr_c2, eph_c2);

% extract the value of common satellite
eph = eph(eph_ord, :);
rcvr = rcvr(rcvr_ord, :);

%% Set given information
ini_pos = [-2694685.473; -4293642.366; 3857878.924]; % initial position
c = 299792458.0; % speed of light, m/s
wedot = 7.2921151467e-5; % WGS 84 value of earth’s rotation rate, r/s
mu = 3.986005e+14; % WGS 84 value of earth's universal gravitation constant, m^3/s^2
F = -4.442807633e-10; % Relativistic correction term constant
tar_XR = [-2700400; -4292560; 3855270];  % target receiver position

[wlat, wlon, walt] = wgsxyz2lla(tar_XR);
tar_XR_lla = [wlat, wlon, walt]'; % target receiver position

%% Tropospheric constant
% Standard atmosphere - Berg, 1948 (Bernese)
% pressure [mbar]
Pr = 1013.25;
% temperature [K]
Tr = 291.15;
% numerical constants for the algorithm [-] [m] [mbar]
Hr = 50.0;
% temperature at sea-level
temp_0 = 15;

%% Main program
% 1. Calculate satellite Earth-centered, Earth-fixed Coordinate system(ECEF) position vector
XS = nan(size(eph,1), 3); % matrix storing satellite ECEF position
dtS = nan(size(eph,1), 1); % array storing satellite clock offset
for i = 1:size(eph,1)
    % assign requried variables 
    pr      =   rcvr(i, 3);
    t       =   eph(i, 1) - pr/c; %initial time
    svid    =   eph(i, 2);
    toc     =   eph(i, 3);
    toe     =   eph(i, 4);
    af0     =   eph(i, 5);
    af1     =   eph(i, 6);
    af2     =   eph(i, 7);
    e       =   eph(i, 9);
    sqrta   =   eph(i, 10);
    dn      =   eph(i, 11);
    m0      =   eph(i, 12);
    w       =   eph(i, 13);
    omg0 	=   eph(i, 14);
    i0      =   eph(i, 15);
    odot    =   eph(i, 16);
    idot    =   eph(i, 17);
    cus     =   eph(i, 18);
    cuc     =   eph(i, 19);
    cis     =   eph(i, 20);
    cic     =   eph(i, 21);
    crs     =   eph(i, 22);
    crc     =   eph(i, 23);    
    
    
    % satellite ECEF position calculation
    a = sqrta^2;   % 1. semimajor axis
    
    n = sqrt(mu / a^3) + dn;   % 2. corrected mean motion (rad/sec)
    
    t_k = t - toe;   % 3. time from ephemeris epoch
    if t_k > 302400  % half week compensation
        t_k = t_k - 604800;
    end
    if t_k < -302400  % half week compensation
        t_k = t_k + 604800;
    end
    
    M_k = m0 + n*(t_k);   % 4. mean anomaly
    
    % 5. eccentric anomaly, Kepler's equation (solve by iterations)
    max_iter = 10; % set maximum number of iterations  
    iter = 0;
    E_k = M_k;
    E_k0 = E_k;
    while iter < max_iter || abs(E_k-E_k0)>1e-12
        E_k0 = E_k;
        E_k = E_k - (E_k - e * sin(E_k) - M_k) / (1 - e * cos(E_k));
        % fprintf('dE_k: %e\n', abs(E_k0-E_k));
        iter = iter + 1;
    end    
  

    % Relativistic clock correction
    dt_r =F*e*sqrta*sin(E_k);

    % satellite clock offset calculation
    dt_sv = af0 + af1 * (t - toc) + af2 * (t - toc).^2+ dt_r;
    dtS(i,1) = dt_sv;

    f_k = atan2(sqrt(1 - e^2) * sin(E_k), cos(E_k) - e);   % 6. true anomaly
    
    phi_k = f_k + w;   % 7. argument of latitude
    
    du_k = cuc * cos(2 * phi_k) + cus * sin(2 * phi_k);   % 8. argument of latitude correction
    dr_k = crc * cos(2 * phi_k) + crs * sin(2 * phi_k);   % 9. radius correction
    di_k = cic * cos(2 * phi_k) + cis * sin(2 * phi_k);   % 10. inclination correction
    
    u_k = phi_k + du_k;   % 11. corrected argument of latitude
    r_k = a * (1 - e * cos(E_k)) + dr_k;   % 12. corrected radius
    i_k = i0 + di_k + idot * t_k;   % 13. corrected inclination 
    
    omg_k = omg0 + (odot - wedot) * t_k - wedot * toe; % 14. corrected longitude of node
    
    x_p = r_k * cos(u_k);   % 15. in-plane x position
    y_p = r_k * sin(u_k);   % 16. in-plane y position
    
    x_s = x_p * cos(omg_k) - y_p * cos(i_k) * sin(omg_k);   % ECEF x-coordinate
    y_s = x_p * sin(omg_k) + y_p * cos(i_k) * cos(omg_k);   % ECEF y-coordinate
    z_s = y_p * sin(i_k);   % ECEF z-coordinate
    
    %coordinates transform correction
    Theta= wedot*pr/c;
    S_c = [cos(Theta) -sin(Theta) 0; sin(Theta) cos(Theta) 0; 0 0 1] * [x_s; y_s; z_s];
    XS(i, :) = [ S_c];  % store corresponding satellite position
end

% least squares to estimate receiver location 
XR = ini_pos;  % initial LS solution: ECEF position, meters
dtR = 0;  % initial LS solution: receiver clock offset, second
[wlat, wlon, walt] = wgsxyz2lla(XR); % initial LS solution: LLA position, deg,deg,meter
x = [XR; dtR]; % solution; 1-3: position in ECEF, meters; 4: receiver clock offset, second;
iter = 0;  % iteration counter

posErr = norm(XR-tar_XR); % positioning error to given target position

dryRunTbl = [iter, nan, nan, nan, nan, nan, XR', dtR, dtR*c, wlat, wlon, walt, nan]; % dryrun table to trace variable change

fprintf('Time of Iteration: %d ----(initial) \n', iter);
fprintf('Initial position: ECEF(m): \nX \t\t  Y\t\t    Z\n%.4f(m), %.4f(m), %.4f(m)\n\n', XR);
fprintf('WGS84 LLA: \nLatitude\t\tLongitude\t\tAltitude \n%.9f(degree), %.9f(degree), %.4f(m))\n', wlat, wlon, walt)
fprintf('\nTotal position error:%.4f(m)\n', posErr);

while 1
    % tropospheric error (saastamoinen model), resolve the tropospheric error each iteration
    if ENA_TROPO_ERR_CORR
        humi = 1.0;
        [wlat, wlon, walt] = wgsxyz2lla(XR);
        tropo_err = zeros(size(XS, 1), 1);
        for i = 1:size(XS, 1)
            enu = wgsxyz2enu(XS(i,:)', wlat, wlon, walt);
            el = asin(enu(3)/sqrt(enu(1)^2 + enu(2)^2));
            if walt < 0
                hgt = 0;
            else
                hgt = walt;
            end
            pres = Pr * (1 - 2.2557e-5 * hgt)^5.2568;
            temp = temp_0 - 6.5e-3 * hgt + 273.16;
            e = 6.108 * humi * exp((17.15 * temp - 4684.0) / (temp - 38.45));
            
            % saastamoninen model
            z = pi / 2.0 - el;
            trph = 0.0022768 * pres / (1.0-0.00266*cos(2.0 * wlat) - 0.00028 * hgt / 1e3) / cos(z);
            trpw = 0.002277 * (1255.0 / temp + 0.05) * e / cos(z);
            tropo_err(i, 1) = trph + trpw;
        end
    end
   
   % tropospheric error (saastamoinen model)
    
    pr = rcvr(:, 3); % extract measured pseudorange

    dist = sqrt(sum((XS - XR').^2, 2)); % geometric distance
    b = dist + c .* (dtR - dtS) + tropo_err; % approx pseudorange

    H = [[XR' - XS]./dist, ones(size(XS, 1),1).*c];  % design matrix
    
    dx = inv(H' * H) * H' * (pr - b); % LS solution
    
    residual = pr - (b + H*dx); % calculate residual
    residualSE = residual' * residual; % residual, squared error
    
    % update solution
    x = x + dx; % update solution (All)
    XR = XR + dx(1:3); % update ECEF position
    [wlat, wlon, walt] = wgsxyz2lla(XR); % convert updated position to wgs84 LLA
    dtR = dtR + dx(4); % update receiver clock offset, second
    posErr = norm(XR-tar_XR); % positioning error to given target position
    
    iter = iter + 1; % update iteration counter
    
    dryRunTbl = [dryRunTbl; iter, dx', residualSE, XR', dtR, dtR*c, wlat, wlon, walt, posErr];
    fprintf('\n----------------\n');
    fprintf('\nTime of Iteration: %d\n', iter);
    fprintf('∆x:\nX: %.4f(m), Y: %.4f(m), Z: %.4f(m), T: %.7f(s), Distance(T*c): %.4f(m)\n', dx, dx(4)*c);
    fprintf('\nUpdated position:\nECEF(m): X:%.4fm, Y:%.4fm, Z:%.4fm\nWGS84 LLA: Latitude:%.9f (degree), Longitude:%.9f (degree), Altitude:%.4f(m)\n', XR, wlat, wlon, walt);
    fprintf('\nUpdated receiver clock offset: R_t: %.7f(s), Distance(R_t*c): %.4f(m)\n', dtR, dtR*c);
    fprintf('\nTotal position error: %.4f(m)\n', posErr);
    fprintf('\nLS residual, squared error: %.4f(m^2)\n', residualSE);
    
    if norm(dx(1:3)) < 1e-4 || ...  % check delta position solution smaller than threshold
            iter > 10  % check iteration number
       break; 
    end
end

% LS summary
fprintf('\n----- Summary ----\n');
fprintf('Init. position ECEF: X: %.4f(m), Y:%.4f(m), Z:%.4f(m) \n -> Final. position ECEF: X: %.4f(m),Y: %.4f(m),Z: %.4f(m)\n', ini_pos, XR);
fprintf('\nTarget position ECEF: X: %.4f(m), Y:%.4f(m), Z:%.4f(m)\n\nPositioning error: %.4f(m)\n', tar_XR, norm([tar_XR-XR]));
fprintf('\nReceiver clock offset: %.7f(s), %.4f(m)\n', dtR, dtR*c);