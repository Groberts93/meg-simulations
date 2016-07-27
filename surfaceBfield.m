function [Bx_tot, By_tot, Bz_tot] = surfaceBfield(Q,R0,xsp,ysp,zsp)

%Sanitization of input
sizQ = size(Q);
sizR0 = size(R0);
sizrx = size(xsp);
sizry = size(ysp);
% sizrz = size(rz);

assert(sum(sizQ == sizR0) == length(sizQ), 'Dipole and position arrays do not match in size');
assert(sum(sizrx == sizry) == length(sizrx), 'Supplied grid arrays do not match in size');
assert(length(sizQ) == 1 || length(sizQ) == 2);

%Constants
mu_0 = 4*pi*1e-7;  %vacuum permeability
Qmag = 1e-9;  %|Q| = 1nAm



%GRIDS
%Create r and start grid centred on 0
R = sqrt(xsp.^2 + ysp.^2 + zsp.^2);  %calculate |r| at each point
Bx_tot = zeros(size(xsp));
By_tot = zeros(size(xsp));
Bz_tot = zeros(size(xsp));

for n = 1:sizQ(1)
    
%Get grid of a = r - R0
[aspx,aspy,aspz] = surftrans(-R0(n,:),xsp,ysp,zsp);
A = sqrt(aspx.^2 + aspy.^2 + aspz.^2);  %calculate |a| = |r - Rq| at each point
[R0x, R0y, R0z] = surfrep(R0(n,:),xsp); %make R0 field


R0r = surfdot(R0(n,:),xsp,ysp,zsp);  %get r0 . r at each point
AR = surfdot2(xsp,ysp,zsp,aspx,aspy,aspz);
F = A.*(R.*A + R.^2 - R0r);
rcoeff = (R.^-1).*(A.^2) + (A.^-1).*(AR) + 2*A + 2*R;
r0coeff = A + 2*R + (A.^-1).*AR;
[rpx, rpy, rpz] = surfscale(rcoeff,xsp,ysp,zsp);
[r0px, r0py, r0pz] = surfscale(r0coeff,R0x,R0y,R0z);
dFx = rpx - r0px;  dFy = rpy - r0py;  dFz = rpz - r0pz;
[cqr0x, cqr0y, cqr0z] = surfrep(cross(Q(n,:),R0(n,:)),xsp);
[Fcqr0x, Fcqr0y, Fcqr0z] = surfscale(F,cqr0x,cqr0y,cqr0z);

qr0r = surfdot(cross(Q(n,:),R0(n,:)),xsp,ysp,zsp);

%Full B-field calculation
Bx = (F.^-2).*(Fcqr0x - qr0r.*dFx);
By = (F.^-2).*(Fcqr0y - qr0r.*dFy);
Bz = (F.^-2).*(Fcqr0z - qr0r.*dFz);

Bx_tot = Bx_tot + Bx;
By_tot = By_tot + By;
Bz_tot = Bz_tot + Bz;

    
end

Bx_tot = (Qmag*mu_0/(4*pi))*Bx_tot;
By_tot = (Qmag*mu_0/(4*pi))*By_tot;
Bz_tot = (Qmag*mu_0/(4*pi))*Bz_tot;



end