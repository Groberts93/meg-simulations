%George Roberts 27/07/2016
%Calculating dipole strengths at surface

clear;
addpath('vfunc');
addpath('sfunc');
Q = [1 0 0];  %Point-like current dipole
tic
n_d = 20;
% Q = randn(n_d,3);
Q = normrows(Q);

R0 = [0 0.09 0]; %Position of current dipole
% R0 = rand(n_d,3) - 0.5*ones(n_d,3);

x = linspace(-2, 2, 100);  %Grid size - should be symmetrical

% GRIDS
% Create r and start grid centred on 0
[sx, sy, sz] = sphere(150);
sx = 0.106 * sx;
sy = 0.106 * sy;
sz = 0.106 * sz;

R = sqrt(sx.^2 + sy.^2 + sz.^2);  %calculate |r| at each point

%getting normal vectors
erx = sx./R;
ery = sy./R;
erz = sz./R;

%Get theta and phi for each point in grid
theta = acos(sz./R);
phi = atan2(sy,sx);

%Get theta unit vector at every point
thx = cos(theta).*cos(phi);
thy = cos(theta).*sin(phi);
thz = -sin(theta);

%Get phi unit vector at every point
phx = -sin(phi);
phy = cos(phi);
phz = zeros(size(sx));

%Calculate B-field from all dipoles
[Bx, By, Bz] = surfaceBfield(Q,R0,sx,sy,sz);

%Calculate components of B-field in spherical coordinates
Br = surfdot2(Bx,By,Bz,erx,ery,erz);
Bt = surfdot2(Bx,By,Bz,thx,thy,thz);
Bp = surfdot2(Bx,By,Bz,phx,phy,phz);
toc

for b = 1:3
    figure;
for n = 1:size(R0,1)
% quiver3(R0(n,1),R0(n,2),R0(n,3),Q(n,1),Q(n,2),Q(n,3),'r-','LineWidth',2);
% hold on;
hold on;
plot3(R0(n,1),R0(n,2),R0(n,3),'bx','LineWidth',8);

end

% [xsp,ysp,zsp] = sphere(150); %Get sphere with 150^2 points
% xsp = 0.106*xsp;
% ysp = 0.106*ysp;
% zsp = 0.106*zsp;
if b == 1
    B = Br; bstring = 'Br';
elseif b == 2
    B = Bt; bstring = 'Bt';
else
    B = Bp; bstring = 'Bp';
end
       
% h = slice(sx,sy,sz,B,xsp,ysp,zsp); %Get slice of B-field at sphere points
h = surf(sx,sy,sz);
set(h,'CData', B);
set(h,'EdgeAlpha',0);
set(h,'FaceAlpha',0.6);
shading interp;
% caxis([-5.6e-16, 5.6e-16]);
% caxis([-8 8]);
axis tight;

colormap(jet);
c = colorbar;
title(['Magnitude of ', bstring, ' component on surface']);
ylabel(c,'Magnitude of component of B-field /T');
end