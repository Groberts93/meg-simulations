clear;
addpath('vfunc');
%George Roberts 18/07/2016
%Spherically symmetric B-field calculation from Sarvas

%Constants
mu_0 = 4*pi*1e-7;  %vacuum permeability
Qmag = 1e-9;  %|Q| = 1nAm

%Vectors
Q = [1 0 0];  %Point-like current dipole
Q = Q/norm(Q);
R0 = [0 -0.75 0]; %Position of current dipole
x = linspace(-2, 2, 30);  %Grid size - should be symmetrical
y = linspace(-2, 2, 8);  %secondary grid for plotting streamlines

%GRIDS
%Create r and start grid centred on 0
[rx, ry, rz] = meshgrid(x,x,x);
[sx, sy, sz] = meshgrid(y,y,y);

%Get grid of a = r - R0
[ax,ay,az] = meshtrans(-R0,rx,ry,rz);

%Constructing the B-field
A = sqrt(ax.^2 + ay.^2 + az.^2);  %calculate |a| = |r - Rq| at each point
R = sqrt(rx.^2 + ry.^2 + rz.^2);  %calculate |r| at each point
[R0x, R0y, R0z] = meshrep(R0,rx); %make R0 field
R0r = meshdot(R0,rx,ry,rz);  %get r0 . r at each point
R0r2 = meshdot2(R0x,R0y,R0z,rx,ry,rz);
AR = meshdot2(rx,ry,rz,ax,ay,az);
F = A.*(R.*A + R.^2 - R0r);
rcoeff = (R.^-1).*(A.^2) + (A.^-1).*(AR) + 2*A + 2*R;
r0coeff = A + 2*R + (A.^-1).*AR;
[rpx, rpy, rpz] = meshscale(rcoeff,rx,ry,rz);
[r0px, r0py, r0pz] = meshscale(r0coeff,R0x,R0y,R0z);
dFx = rpx - r0px;  dFy = rpy - r0py;  dFz = rpz - r0pz;

[cqr0x, cqr0y, cqr0z] = meshrep(cross(Q,R0),rx);
[Fcqr0x, Fcqr0y, Fcqr0z] = meshscale(F,cqr0x,cqr0y,cqr0z);

qr0r = meshdot(cross(Q,R0),rx,ry,rz);

%getting normal vectors
erx = rx./R;
ery = ry./R;
erz = rz./R;

%Full B-field calculation
Bx = (Qmag*mu_0/(4*pi))*(F.^-2).*(Fcqr0x - qr0r.*dFx);
By = (Qmag*mu_0/(4*pi))*(F.^-2).*(Fcqr0y - qr0r.*dFy);
Bz = (Qmag*mu_0/(4*pi))*(F.^-2).*(Fcqr0z - qr0r.*dFz);

%Get theta and phi for each point in grid
theta = acos(rz./R);
phi = atan2(ry,rx);

%Get theta unit vector at every point
thx = cos(theta).*cos(phi);
thy = cos(theta).*sin(phi);
thz = -sin(theta);

%Get phi unit vector at every point
phx = -sin(phi);
phy = cos(phi);
phz = zeros(size(rx));

%Calculate components of B-field in spherical coordinates
Br = meshdot2(Bx,By,Bz,erx,ery,erz);
Bt = meshdot2(Bx,By,Bz,thx,thy,thz);
Bp = meshdot2(Bx,By,Bz,phx,phy,phz);


% quiver3(rx,ry,rz,thx,thy,thz);
% hold on;

quiver3(R0(1),R0(2),R0(3),Q(1),Q(2),Q(3),'r-','LineWidth',2);
hold on;
% streamline(rx,ry,rz,phihatx,phihaty,phihatz,sx,sy,sz);
% hold on;


[xsp,ysp,zsp] = sphere(150); %Get sphere with 150^2 points
h = slice(rx,ry,rz,Br,xsp,ysp,zsp); %Get slice of B-field at sphere points
set(h,'EdgeAlpha',0);
set(h,'FaceAlpha',0.6);
shading interp;
caxis([-5.6e-16, 5.6e-16]);
axis tight;
colormap(jet);
c = colorbar;
title('Vector plot of B-field and magnitude of component on surface');
ylabel(c,'Magnitude of component of B-field /T');