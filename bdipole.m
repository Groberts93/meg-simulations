clear;
addpath('vfunc');
%George Roberts 18/07/2016
%Radial component of B-field

%Constants
mu_0 = 4*pi*1e-7;  %vacuum permeability
Qmag = 1e-9;  %|Q| = 1nAm
const = mu_0 * Qmag / (4*pi);

%Vectors
Q = [1; 0; 0];  %Point-like current dipole
Q = Q/norm(Q); %normalise dipole vector, now it's a unit vector
Rq = [0; -0.75; 0]; %Position of current dipole
x = linspace(-4, 4, 100);  %High density grid for slicing

%Create XYZ grid centred on 0
[Rx, Ry, Rz] = meshgrid(x,x,x);

%Get grid of a = r - Rq
[Ax,Ay,Az] = meshtrans(-Rq,Rx,Ry,Rz);

A = sqrt(Ax.^2 + Ay.^2 + Az.^2);  %and calculate A = |r - Rq| at each point
R = sqrt(Rx.^2 + Ry.^2 + Rz.^2);  %calculate 
erx = Rx./R;
ery = Ry./R;
erz = Rz./R;

%Calculate the scalar triple product Q x R . er
[qrqx, qrqy, qrqz] = meshrep(cross(Q,Rq),Rx);
% dprod = meshdot2(qrqx,qrqy,qrqz,erx,ery,erz);
% Br = -const*(dprod./(A.^3));  %Finally, calculate scalar field Br

Bx = -const*qrqx.*(A.^-3);
By = -const*qrqy.*(A.^-3);
Bz = -const*qrqz.*(A.^-3);



%Vector B-field calculation
x = linspace(-0.9, 0.9, 4);  %Grid size - should be symmetrical
[X2, Y2, Z2] = meshgrid(x,x,x);
[Xr2,Yr2,Zr2] = meshtrans(-Rq,X2,Y2,Z2);
R2 = sqrt(Xr2.^2 + Yr2.^2 + Zr2.^2);
[U2,V2,W2] = meshcross(Q,Xr2,Yr2,Zr2);  %Calculate Q x (r - Rq) at each point

%Normalise vectors according to |r - Rq|
U2 = const*U2./(R2.^3);
V2 = const*V2./(R2.^3);
W2 = const*W2./(R2.^3);

%Get theta and phi for each point in grid
theta = acos(Rz./R);
phi = atan2(Ry,Rx);

%Get theta unit vector at every point
thx = cos(theta).*cos(phi);
thy = cos(theta).*sin(phi);
thz = -sin(theta);

%Get phi unit vector at every point
phx = -sin(phi);
phy = cos(phi);
phz = zeros(size(Rx));

%Get theta and phi components of B field
Br = meshdot2(Bx,By,Bz,erx,ery,erz);
Bt = meshdot2(Bx,By,Bz,thx,thy,thz);
Bp = meshdot2(Bx,By,Bz,phx,phy,phz);

quiver3(X2,Y2,Z2,U2,V2,W2,'LineWidth',2);  %Plot field in blue
% quiver3(Rx,Ry,Rz,thx,thy,thz,'LineWidth',2);  %Plot field in blue
hold on;
quiver3(Rq(1),Rq(2),Rq(3),Q(1),Q(2),Q(3),1.0,'r-','LineWidth',2);
%Plot dipole in red
hold on;
[xsp,ysp,zsp] = sphere(150); %Get sphere with 150^2 points
h = slice(Rx,Ry,Rz,Br,xsp,ysp,zsp); %Get slice of B-field at sphere points
set(h,'EdgeAlpha',0);
set(h,'FaceAlpha',0.6);
shading interp
caxis([-5.6e-16, 5.6e-16]);
axis equal;
colormap(jet);
c = colorbar;
title('Vector plot of B-field and magnitude of component on surface');
ylabel(c,'Magnitude of component of B-field /T');