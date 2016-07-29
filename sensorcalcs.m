%Point B-field plotting
%George Roberts 27/07/2016

clear;

%Add function paths
addpath('vfunc');
addpath('sfunc');

load('data/spherepts.mat');  %Load sphere points from EEGMesh - called 
%EEGPts1, EEGPts2, EEGPts3

Q = [1 1 0];  %Point-like current dipole
R0 = [0 0.09 0]; %Position of current dipole

xpts = EEGPts1(:,1);  
ypts = EEGPts1(:,2);
zpts = EEGPts1(:,3);

xp = 0.106 * (xpts/20); %Set EEG locations to radius of 10.6cm from origin
yp = 0.106 * (ypts/20);
zp = 0.106 * (zpts/20);

R = sqrt(xp.^2 + yp.^2 + zp.^2);  %calculate |r| at each point

%getting normal vectors
erx = xp./R;
ery = yp./R;
erz = zp./R;

%Get theta and phi for each point in grid
theta = acos(zp./R);
phi = atan2(yp,xp);

%Get theta unit vector at every point
thx = cos(theta).*cos(phi);
thy = cos(theta).*sin(phi);
thz = -sin(theta);

%Get phi unit vector at every point
phx = -sin(phi);
phy = cos(phi);
phz = zeros(size(xp));

%Calculate components of B-field at every point in EEG montage
[Bx_tot, By_tot, Bz_tot] = pointsBfield(Q,R0,[xp, yp, zp]);

%Get radial, theta, phi components
Br = Bx_tot.*erx + By_tot.*ery + Bz_tot.*erz;
Bt = Bx_tot.*thx + By_tot.*thy + Bz_tot.*thz;
Bp = Bx_tot.*phx + By_tot.*phy + Bz_tot.*phz;

%Plot these as scatter3
scattersize = 50*ones(size(Br));
scatter3(xp, yp, zp, scattersize, Bp, 'filled');
colorbar;

