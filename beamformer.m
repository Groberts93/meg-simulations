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
load('data/data_2dips.mat');
% load('data/bdata6.mat');
[rx, ry, rz] = meshgrid(linspace(-0.085,0.085,30));
sizr = size(rx);

%get dipole orientations at every point on grid
n_theta = 4;
theta_t = linspace(0,pi,n_theta);
[vtx, vty, vtz] = dipolefangrid(rx, ry, rz, theta_t);

Br = Br.*1e15;

nch = size(Br,1);
f = 600;
nt = 300*f;
t = randn(1, nt);  %time - assume normal distributed around zero

%init points and get normal vectors
R = sqrt(xp.^2 + yp.^2 + zp.^2);  %calculate |r| at each point

%getting normal vectors
erx = xp./R;
ery = yp./R;
erz = zp./R;



B = Br*t + 100*randn(nch,nt);
C = cov(B');
Cinv = inv(C);
Z = zeros(sizr(1),sizr(1),sizr(1),n_theta);
maxz = zeros(sizr(1),sizr(1));

dlocx = 0;
dlocy = 0.045;
dlocz = 0.07;

[~,Ix] = min(abs(rx(:) - dlocx));
[~,Iy] = min(abs(ry(:) - dlocy));
[~,Iz] = min(abs(rz(:) - dlocz));

[xi, xj, xk] = ind2sub(sizr,Ix);
[yi, yj, yk] = ind2sub(sizr,Iy);
[zi, zj, zk] = ind2sub(sizr,Iz);



for xprb = 1:30
    for yprb = 1:30
        for zprb = 1:30
            
            for tprb = 1:n_theta
               
%                 Q = [vtx(xprb, yprb, zprb, tprb), ...
%                     vty(xprb, yprb, zprb, tprb), ...
%                     vtz(xprb, yprb, zprb, tprb)];
%                 
%                 R0 = [rx(xprb, yprb, zprb), ry(xprb, yprb, zprb), rz(xprb, yprb, zprb)];

%                 Q = [vtx(xi, xj, xk, tprb), ...
%                     vty(yi, yj, yk, tprb), ...
%                     vtz(zi, zj, zk, tprb)];
%                 
%                 R0 = [rx(xi, xj, xk), ry(yi, yj, yk), rz(zi, zj, zk)];

                Q = [vtx(xprb, yprb, zprb, tprb), ...
                    vty(xprb, yprb, zprb, tprb), ...
                    vtz(xprb, yprb, zprb, tprb)];
                
                R0 = [rx(xprb, yprb, zprb), ry(xprb, yprb, zprb), rz(xprb, yprb, zprb)];

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

% theta_data_max_z = squeeze(Z(15,23,27,:));
% plot(theta_data_max_z);
% [~, theta_max] = max(theta_data_max_z);
% disp(['maximum power at theta = ', num2str(theta_max), ' degrees']);

figure;
h = imagesc(maxz(:,:,27));
v = caxis;

figure;
for i = 1:30
subplot(5,6,i);
imagesc(maxz(:,:,i),v);
end

% figure;
% plot(squeeze(maxz(15,:,27)));
% [~, ind] = max(maxz)