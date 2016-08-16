%George Roberts - 28/07/16
%Beamforming reconstruction of source power

clear;

addpath('bfunc');
addpath('sfunc');
addpath('vfunc');

%load data, initialise grid
load('newdata/B_15-Aug-2016-15:33:08.mat');
Br1 = Br;
Bt1 = Bt;
load('newdata/B_15-Aug-2016-15:33:37.mat');
Br2 = Br;
Bt2 = Bt;
load('newdata/B_15-Aug-2016-15:33:47.mat');
Br3 = Br;
Bt3 = Bt;

Br = [Br1, Br2, Br3];
Bt = [Bt1, Bt2, Bt3];

Brt1 = [Br1; Bt1];
Brt2 = [Br2; Bt2];
Brt3 = [Br3; Bt3];


[rx, ry, rz] = meshgrid(linspace(-0.085,0.085,31));
sizr = size(rx);

%get dipole orientations at every point on grid
n_theta = 180;
theta_t = linspace(0,pi,n_theta);
[vtx, vty, vtz] = dipolefangrid(rx, ry, rz, theta_t);  %ALSO SWAPPED HERE
% 
% Br = Br.*1e15;
% Bt = Bt.*1e15;
Brt1 = Brt1.*1e15;
Brt2 = Brt2.*1e15;
Brt3 = Brt3.*1e15;

nch = size(Brt1,1);
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

%Get theta and phi for each point
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


B = Brt1*t1 + Brt2*t2 + Brt3*t3 + 200*randn(nch,nt);
C = cov(B');
Cinv = inv(C);
Z = zeros(sizr(1),sizr(1),sizr(1),n_theta);
maxz = zeros(sizr(1),sizr(1),sizr(1));

dlocx = 0.045;
dlocy = 0.00;
dlocz = 0.07;

%x and y coordinates get switched here.  watch out

[~,Ix] = min(abs(rx(:) - dlocx));
[~,Iy] = min(abs(ry(:) - dlocy));
[~,Iz] = min(abs(rz(:) - dlocz));

[xi, xj, xk] = ind2sub(sizr,Ix);
[yi, yj, yk] = ind2sub(sizr,Iy);
[zi, zj, zk] = ind2sub(sizr,Iz);

xpt = max([xi, xj, xk]);
ypt = max([yi, yj, yk]);
zpt = max([zi, zj, zk]);

for xprb = 1:31
    for yprb = 1:31
        for zprb = zpt:zpt
            
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
                Lt = Bx_tot.*thx + By_tot.*thy + Bz_tot.*thz;
                
                L = [Lr; Lt];
                
                Wtr = (L'*Cinv)/(L'*Cinv*L);
                Z(xprb,yprb,zprb,tprb) = (Wtr*C*(Wtr'))/(Wtr*Wtr');
              
                
            end
            maxz(xprb, yprb, zprb) = max(abs(Z(xprb, yprb, zprb,:)));
        end
    end
end

figure;
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

% figure;
% for i = 1:30
% subplot(6,5,i);
% imagesc(maxz(:,:,i),v);
% end

% figure;
% plot(squeeze(maxz(15,:,27)));
% [~, ind] = max(maxz)