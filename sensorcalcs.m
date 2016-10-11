%Point B-field plotting
%George Roberts 27/07/2016

clear;

%Add function paths
addpath('vfunc');
addpath('sfunc');

user_save_data = 1;  %Change this to 1 if you want to save the Bdata - will
%save dipole location, orientation, and Bdata (with xp, yp, zp)

load('data/spherepts.mat');  %Load sphere points from EEGMesh - called 
%EEGPts1, EEGPts2, EEGPts3
ntheta = 40;
theta_t = linspace(0,pi,ntheta);
[rx, ry, rz] = meshgrid(linspace(-0.085,0.085,30));
[vtx, vty, vtz] = dipolefangrid(rx,ry,rz,theta_t);

%debugging the dipole fan

p = 15; q = 23; s = 27;

rvect = [rx(p,q,s), ry(p,q,s), rz(p,q,s)];
rmult = repmat(rvect, [ntheta 1]);

R0 = [rmult(:,1), rmult(:,2), rmult(:,3)];
Q = [squeeze(vtx(p,q,s,:)), squeeze(vty(p,q,s,:)), squeeze(vtz(p,q,s,:))];


% Q = [0.9397, 0, -0.3420];  %Point-like current dipole
% Q = rand(20, 3) - 0.5;
% Q = normrows(Q);

% Q = [0 0 1];

% R0 =  [0.0274,  0, 0.0752];
% Rrand = 0.12 * (rand(500, 3) - 0.5);
% Rrand = Rrand(Rrand(:,3) > 0, :);
% R0 = Rrand(1:20, :);
% R0 = [0 0 0.09];

anglestr = [' '];

%Get projection angles and display
% proj_angle = getprojangle(Q,R0) + 90;
% for n = 1:size(proj_angle,1)
% anglestr = [anglestr, '  ', sprintf('%g',proj_angle(n,1))];
% end

% disp(anglestr);

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
axis equal;
xlabel('x axis');
ylabel('y axis');
zlabel('z axis');

% for npt = 1:size(xp,1)
%    quiver3(xp(npt), yp(npt), zp(npt), 0.02*thx(npt), 0.02*thy(npt), 0.02*thz(npt), 'r-');
%    hold on;
%     
%     
% end

% quiver3(xp, yp, zp, 0.02*thx, 0.02*thy, 0.02*thz, 'r-');
quiver3(0, 0, 0, 1.1*R0(1,1), 1.1*R0(1,2), 1.1*R0(1,3), 'r-');

Br_fT = Br * 1e15;
Bp_fT = Bp * 1e15;

Br_rms = sqrt(mean(Br_fT.^2, 1))
Bp_rms = sqrt(mean(Bp_fT.^2, 1))

%SAVE THAT DATA
if (user_save_data == 1)
    ntime = now;
    datastr = ['sequencedata/', 'B_', datestr(ntime, 1), '-', datestr(ntime,13), '.mat'];
   save(datastr, 'Br', 'Bt', 'xp', 'yp', 'zp', 'R0', 'Q');
end
