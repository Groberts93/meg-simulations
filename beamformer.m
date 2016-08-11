%George Roberts - 28/07/16
%Beamforming reconstruction of source power

clear;

addpath('bfunc');
addpath('sfunc');
addpath('vfunc');

% Q = [1 0 0];  %Point-like current dipole
% R0 = [0 0.045 0.07]; %Position of current dipole

%load data, initialise grid
% load('data/bdata0_0p045_0p07.mat');
% load('data/data_2dips.mat');
load('data/b_test_102.mat');
Br1 = Br;
load('data/b_test_100_556.mat');
Br2 = Br;
load('data/b_test_110_m546.mat');
Br3 = Br;

Br = [Br1, Br2, Br3];


[rx, ry, rz] = meshgrid(linspace(-0.085,0.085,30));
sizr = size(rx);

%get dipole orientations at every point on grid
n_theta = 5;
theta_t = linspace(0,pi,n_theta);
[vtx, vty, vtz] = dipolefangrid(rx, ry, rz, theta_t);  %ALSO SWAPPED HERE

Br = Br.*1e15;

nch = size(Br,1);
f = 600;
nt = 300*f;
t1 = randn(1, nt);  %time - assume normal distributed around zero
t2 = randn(1, nt);
t3 = randn(1, nt);

%init points and get normal vectors
R = sqrt(xp.^2 + yp.^2 + zp.^2);  %calculate |r| at each point

%getting normal vectors
erx = xp./R;
ery = yp./R;
erz = zp./R;


% 
B = Br(:,1)*t1 + Br(:,2)*t2 + Br(:,3)*t3 + 100*randn(nch,nt);
% B = Br(:,1)*t1 + Br(:,2)*t2 + 100*randn(nch,nt);
C = cov(B');
Cinv = inv(C);
Z = zeros(sizr(1),sizr(1),sizr(1),n_theta);
maxz = zeros(sizr(1),sizr(1),sizr(1));

dlocx = 0.045;
dlocy = 0.00;
dlocz = 0.07;

% dlocx = 0.02;
% dlocy = 0.03;
% dlocz = 0.06;


[~,Ix] = min(abs(rx(:) - dlocx));
[~,Iy] = min(abs(ry(:) - dlocy));
[~,Iz] = min(abs(rz(:) - dlocz));

[xi, xj, xk] = ind2sub(sizr,Ix);
[yi, yj, yk] = ind2sub(sizr,Iy);
[zi, zj, zk] = ind2sub(sizr,Iz);

xpt = max([xi, xj, xk]);
ypt = max([yi, yj, yk]);
zpt = max([zi, zj, zk]);

for xprb = 1:30
    for yprb = 1:30
        for zprb = 1:30
            
            for tprb = 1:n_theta
               
                Q = [vtx(xprb, yprb, zprb, tprb), ...
                    vty(xprb, yprb, zprb, tprb), ...
                    vtz(xprb, yprb, zprb, tprb)];
                
%                 if (tprb == 90)
%                     Q
%                 end
%                 
                R0 = [rx(xprb, yprb, zprb), ry(xprb, yprb, zprb), rz(xprb, yprb, zprb)];

%                 Q = [vtx(xi, xprb, xk, tprb), ...
%                     vty(yprb, yj, yk, tprb), ...
%                     vtz(zi, zj, zprb, tprb)];
                
%                 R0 = [rx(xi, xprb, xk), ry(yprb, yj, yk), rz(zi, zj, zprb)];

%                 Q = [vtx(xprb, yprb, zprb, tprb), ...
%                     vty(xprb, yprb, zprb, tprb), ...
%                     vtz(xprb, yprb, zprb, tprb)];
%                 
%                 R0 = [rx(xprb, yprb, zprb), ry(xprb, yprb, zprb), rz(xprb, yprb, zprb)];

                %need to get lead field here
                %Calculate components of B-field at every point in EEG montage
                [Bx_tot, By_tot, Bz_tot] = pointsBfield(Q,R0,[xp, yp, zp]);
                
                %Get radial, theta, phi components
                Lr = Bx_tot.*erx + By_tot.*ery + Bz_tot.*erz;
                
                Wtr = (Lr'*Cinv)/(Lr'*Cinv*Lr);
                Z(xprb,yprb,zprb,tprb) = (Wtr*C*(Wtr'))/(Wtr*Wtr');
              
                
            end
            maxz(xprb, yprb, zprb) = max(abs(Z(xprb, yprb, zprb,:)));
        end
    end
end

theta_data_max_z = squeeze(Z(xpt,ypt,zpt,:));
plot(theta_data_max_z);
[~, theta_max] = max(theta_data_max_z);
disp(['maximum power at theta = ', num2str(theta_max), ' degrees']);

figure;
h = imagesc(maxz(:,:,zpt));
minz = min(min(min(min(Z))));
maxmaxz = max(max(max(maxz)));
v = [minz maxmaxz];
colorbar;

figure;
for i = 1:30
subplot(6,5,i);
imagesc(maxz(:,:,i),v);
end

% figure;
% plot(squeeze(maxz(15,:,27)));
% [~, ind] = max(maxz)