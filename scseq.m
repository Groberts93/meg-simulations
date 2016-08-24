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

xpts = EEGPts1(:,1);  
ypts = EEGPts1(:,2);
zpts = EEGPts1(:,3);
xp = 0.106 * (xpts/20); %Set sensor locations to radius of 10.6cm from origin
yp = 0.106 * (ypts/20);
zp = 0.106 * (zpts/20);
[rx, ry, rz] = meshgrid(linspace(-0.085,0.085,30));

% Q = [1 0 0; 1 0 0; 1 0 0];  %Point-like current dipole
ndips = 2;
% Q = repmat([1 0 0], ndips, 1);
r_const = 0.08;
t_const = 85;
Q = [sind(t_const), 0, -cosd(t_const); sind(t_const), 0, cosd(t_const)];
Q = normrows(Q);
Br_seq = zeros(size(xpts,1), ndips);
Bt_seq = zeros(size(xpts,1), ndips);


randR =  0.2 * (rand(10000,3) - 0.5);
randR = randR(randR(:,3) > 0, :);
loc_r = randR(:,1).^2 + randR(:,2).^2 + randR(:,3).^2;
loc_r = sqrt(loc_r);
randR = randR(loc_r < 0.1, :);
randR = randR(1:ndips,:);
% R0 = [-0.03 0.04 0.07; 0.05 0.00 0.07];
R0 = [r_const*cosd(t_const), 0.00, r_const*sind(t_const); -r_const*cosd(t_const), 0.00, r_const*sind(t_const)]; 
% R0 = randR;

% anglestr = [' '];
% 
% %Get projection angles and display
% proj_angle = getprojangle(Q,R0) + 90;
% for n = 1:size(proj_angle,1)
% anglestr = [anglestr, '  ', sprintf('%g',proj_angle(n,1))];
% end

% disp(anglestr);





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


for ndip = 1:ndips
%Calculate components of B-field at every point in EEG montage
[Bx_tot, By_tot, Bz_tot] = pointsBfield(Q(ndip,:),R0(ndip,:),[xp, yp, zp]);

%Get radial, theta, phi components
Br = Bx_tot.*erx + By_tot.*ery + Bz_tot.*erz;
Bt = Bx_tot.*thx + By_tot.*thy + Bz_tot.*thz;
Bp = Bx_tot.*phx + By_tot.*phy + Bz_tot.*phz;

%Plot origin and dipoles

figure;
plot3(0, 0, 0,'go','LineWidth',3);  %origin
hold on;
quiver3(R0(ndip,1), R0(ndip,2), R0(ndip,3), 0.05*Q(ndip,1), 0.05*Q(ndip,2), 0.05*Q(ndip,3), 'b-', 'LineWidth', 1);
hold on;
plot3(R0(ndip,1), R0(ndip,2), R0(ndip,3),'ro','LineWidth',2);
hold on;

%Plot sensors as scatter3

scattersize = 50*ones(size(Br));
scatter3(xp, yp, zp, scattersize, Br, 'filled');
colorbar;
axis equal;
xlabel('x axis');
ylabel('y axis');
zlabel('z axis');

Br_seq(:, ndip) = Br;
Bt_seq(:, ndip) = Bt;

end

%SAVE THAT DATA
if (user_save_data == 1)
    ntime = now;
    datastr = ['sequencedata/', 'B_', datestr(ntime, 1), '-', datestr(ntime,13), '.mat'];
   save(datastr, 'Br_seq', 'Bt_seq', 'xp', 'yp', 'zp', 'R0', 'Q');
end
