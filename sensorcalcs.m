%Point B-field plotting
%George Roberts 27/07/2016

clear;

%Add function paths
addpath('vfunc');
addpath('sfunc');

load('data/spherepts.mat');  %Load sphere points from EEGMesh - called 
%EEGPts1, EEGPts2, EEGPts3
ntheta = 180;
theta_t = linspace(0,pi,ntheta);
[rx, ry, rz] = meshgrid(linspace(-0.1,0.1,8));
[vtx, vty, vtz] = dipolefangrid(rx,ry,rz,theta_t);

% p = 5; q = 4; s = 6;
% 
% rvect = [rx(p,q,s), ry(p,q,s), rz(p,q,s)];
% rmult = repmat(rvect, [ntheta 1]);

% R0 = [rmult(:,1), rmult(:,2), rmult(:,3)];
% Q = [squeeze(vtx(p,q,s,:)), squeeze(vty(p,q,s,:)), squeeze(vtz(p,q,s,:))];


% Q = [0.7 1 1];  %Point-like current dipole
% Q = normrows(Q);
% R0 = [0 0.045 0.07]; %Position of current dipole

Q = [1 0 0; -1 0 0; 1 1 1];  %Point-like current dipole
Q = normrows(Q);
R0 = [0.02 0.045 0.07; -0.02 -0.045 0.07; 0.04 -0.040 0.073]; %Position of current dipole

% proj_angle = getprojangle(Q,R0) + 90;
% disp(['Projection angle is ', num2str(proj_angle)]);

xpts = EEGPts1(:,1);  
ypts = EEGPts1(:,2);
zpts = EEGPts1(:,3);

xp = 0.106 * (xpts/20); %Set sensor locations to radius of 10.6cm from origin
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

%Plot origin and dipoles
plot3(0, 0, 0,'go','LineWidth',3);  %origin
hold on;

for ndip = 1:size(Q,1)
    quiver3(R0(ndip,1), R0(ndip,2), R0(ndip,3), 0.05*Q(ndip,1), 0.05*Q(ndip,2), 0.05*Q(ndip,3), 'b-', 'LineWidth', 1);
    hold on;
    plot3(R0(ndip,1), R0(ndip,2), R0(ndip,3),'ro','LineWidth',2);
hold on;
    
end

%Plot sensors as scatter3
scattersize = 50*ones(size(Br));
scatter3(xp, yp, zp, scattersize, Br, 'filled');
colorbar;

