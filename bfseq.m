%George Roberts - 28/07/16
%Beamforming reconstruction of source power

clear;

addpath('bfunc');
addpath('sfunc');
addpath('vfunc');

%load data, initialise grid
% load('newdata/B_15-Aug-2016-15:33:08.mat');
% Br1 = Br;
% Bt1 = Bt;
% load('newdata/B_15-Aug-2016-15:33:37.mat');
% Br2 = Br;
% Bt2 = Bt;
% load('newdata/B_15-Aug-2016-15:33:47.mat');
% Br3 = Br;
% Bt3 = Bt;

% load('sequencedata/B_22-Aug-2016-13:15:45.mat');
% load('sequencedata/B_22-Aug-2016-13:22:44.mat');
% load('sequencedata/B_22-Aug-2016-12:46:49.mat');
load('sequencedata/B_22-Aug-2016-13:33:20.mat');

ndips = size(R0,1);

% Br = [Br1, Br2, Br3];
% Bt = [Bt1, Bt2, Bt3];
% Brt = [Br; Bt];

Br = Br_seq;
Bt = Bt_seq;
Brt = [Br_seq; Bt_seq];

dlocx = -0.03;
dlocy = 0.04;
dlocz = 0.07;

n_theta = 180;


[rx, ry, rz] = meshgrid(linspace(-0.085,0.085,31));
sizr = size(rx);

rtmode = 0;

while (rtmode > 3 || rtmode < 1)
    rtmode = input('Enter 1 to use radial components, 2 for tangential components, or 3 for both: ');
end

bfmode = 0;
while (bfmode > 3 || bfmode < 1)
    bfmode = input('Choose mode.  Enter 1 for single voxel, 2 for xy slice at z plane, 3 for full xyz sweep: ');
end

if (bfmode == 1)
prompt = {'x coordinate:','y coordinate:' 'z coordinate'};
dlg_title = 'Single voxel mode.  Enter coordinates';
num_lines = 1;
defaultans = {num2str(dlocx),num2str(dlocy),num2str(dlocz)};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
dlocx = str2num(answer{1});
dlocy = str2num(answer{2});
dlocz = str2num(answer{3});

elseif (bfmode == 2)
prompt = {'z layer', 'ntheta'};
dlg_title = ['xy slice mode.  Enter z-layer from 1-', num2str(sizr(1)), ' and ntheta'];
num_lines = [1 50];
defaultans = {num2str(28), num2str(n_theta)};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
zslice = str2num(answer{1});
n_theta = str2num(answer{2});

elseif (bfmode == 3)
    prompt = {'ntheta'};
dlg_title = ['Full xyz sweep mode.  Enter ntheta '];
num_lines = [1 50];
defaultans = {num2str(n_theta)};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
n_theta = str2num(answer{1});
end

dlocy = -dlocy;  %Negate y coordinate.  Now you can enter dlocx, dlocy and dlocz as
%if they were the xyz coordinates of the dipole in sensorcalcs.m and
%everything should square up with the original coordinates.  I think.




%get dipole orientations at every point on grid

theta_t = linspace(0,pi,n_theta);
[vtx, vty, vtz] = dipolefangrid(rx, ry, rz, theta_t);  %ALSO SWAPPED HERE
% 
Br = Br.*1e15;
Bt = Bt.*1e15;
Brt = Brt.*1e15;

if (rtmode == 1)
    nch = size(Br,1);
elseif (rtmode == 2)
    nch = size(Bt,1);
else
    nch = size(Brt,1);
end

% nch = size(Brt,1);
f = 600;
nt = 300*f;
t1 = randn(1, nt);  %time - assume normal distributed around zero
t2 = randn(1, nt);
t3 = randn(1, nt);

tn = randn(ndips, nt);

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

B = zeros(nch, nt);

for ndip = 1:ndips
    if (rtmode == 1)
        B = B + Br(:,ndip)*tn(ndip,:);
    elseif (rtmode == 2)
        B = B + Bt(:,ndip)*tn(ndip,:);
    else
        B = B + Brt(:,ndip)*tn(ndip,:);
    end
end

B = B + 200*randn(nch, nt);

C = cov(B');
Cinv = inv(C);
Z = zeros(sizr(1),sizr(1),sizr(1),n_theta);
maxz = zeros(sizr(1),sizr(1),sizr(1));



%x and y coordinates get switched here.  watch out

[~,Iy] = min(abs(rx(:) - dlocx));
[~,Ix] = min(abs(ry(:) - dlocy));
[~,Iz] = min(abs(rz(:) - dlocz));

[xi, xj, xk] = ind2sub(sizr,Ix);
[yi, yj, yk] = ind2sub(sizr,Iy);
[zi, zj, zk] = ind2sub(sizr,Iz);

xpt = max([xi, xj, xk]);
ypt = max([yi, yj, yk]);
zpt = max([zi, zj, zk]);

if (bfmode == 1)
    xrange = [xpt, xpt];
    yrange = [ypt, ypt];
    zrange = [zpt, zpt];
    zslice = zpt;
elseif (bfmode == 2)
    xrange = [1, sizr(1)];
    yrange = [1, sizr(1)];
    zrange = [zslice, zslice];
elseif (bfmode == 3)
    xrange = [1, sizr(1)];
    yrange = [1, sizr(1)];
    zrange = [1, sizr(1)];
    zslice = zpt;
end

tic;
for xprb = xrange(1):xrange(2)
    for yprb = yrange(1):yrange(2)
        for zprb = zrange(1):zrange(2)
            
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
                if (rtmode == 1)
                Lr = Bx_tot.*erx + By_tot.*ery + Bz_tot.*erz;
                L = Lr;
                elseif (rtmode == 2)
                Lt = Bx_tot.*thx + By_tot.*thy + Bz_tot.*thz;
                L = Lt;
                else
                    Lr = Bx_tot.*erx + By_tot.*ery + Bz_tot.*erz;
                    Lt = Bx_tot.*thx + By_tot.*thy + Bz_tot.*thz;
                    L = [Lr; Lt];
                end
                
                
                Wtr = (L'*Cinv)/(L'*Cinv*L);
                Z(xprb,yprb,zprb,tprb) = (Wtr*C*(Wtr'))/(Wtr*Wtr');
              
                
            end
            maxz(xprb, yprb, zprb) = max(abs(Z(xprb, yprb, zprb,:)));
        end
    end
    disp(['x slice = ', num2str(xprb)]);
end
toc;

% for zn = 1:31
%     for tn = 1:n_theta
%    maxz(:,:,zn) = maxz(:,:,zn)';
%     end
% end

figure;
theta_data_max_z = squeeze(Z(xpt,ypt,zslice,:));
plot(theta_data_max_z);
[~, theta_max] = max(theta_data_max_z);
disp(['maximum power at theta = ', num2str(theta_max), ' degrees']);

figure;
minz = min(min(min(min(Z))));
maxmaxz = max(max(max(maxz)));
v = [minz maxmaxz];
h = imagesc(maxz(:,:,zslice), v);
colorbar;

figure;
for i = 1:31
    subplot(6,6,i);
    imagesc(maxz(:,:,i),v);
end

% figure;
% plot(squeeze(maxz(15,:,27)));
% [~, ind] = max(maxz)