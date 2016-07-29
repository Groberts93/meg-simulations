%George Roberts 22/07/2016
%Inverse problem in MEG

clear;

load('data/Br_sparse.mat');  %load Br, rx, ry, rz
sizr = size(rx);
sqch = 16;
Nch = sqch^2;
[xsp, ysp, zsp] = sphere(sqch-1);
h = slice(rx, ry, rz, Br, xsp, ysp, zsp);

for n = 1:sqch
hold on;
plot3(xsp(n,:), ysp(n,:), zsp(n,:),'rx');

end

L = get(h,'CData');
sizl = size(L);
shading interp;
% caxis([-5.6e-16, 5.6e-16]);
caxis([-8 8]);
set(h,'EdgeAlpha',0);
set(h,'FaceAlpha',0.6);

L = reshape(L,1,sizl(1)*sizl(2));
t = linspace(0, 0.5, 150);
f = 18;
q = 1e-9 * sin(2*pi*f*t);


Bt = L'*q;

% GRIDS
% Get back R and unit vectors from passed grid
R = sqrt(rx.^2 + ry.^2 + rz.^2);  %calculate |r| at each point

%getting normal vectors
erx = rx./R;
ery = ry./R;
erz = rz./R;

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

p = 4; q = 4; s = 4;

ntheta = 180;
theta_t = linspace(0, pi, ntheta);

[vtx, vty, vtz] = dipolefangrid(rx,ry,rz,theta_t);

% sizt = [sizr ntheta];
% vtx = zeros(sizt);
% vty = zeros(sizt);
% vtz = zeros(sizt);
% 
% for xi = 1:sizr(1)
%     for yi = 1:sizr(2)
%         for zi = 1:sizr(3)
% thx1 = thx(xi,yi,zi);
% thy1 = thy(xi,yi,zi);
% thz1 = thz(xi,yi,zi);
% phx1 = phx(xi,yi,zi);
% phy1 = phy(xi,yi,zi);
% phz1 = phz(xi,yi,zi);
% 
% 
% vtx(xi,yi,zi,:) = cos(theta_t).*thx1 + sin(theta_t).*phx1;
% vty(xi,yi,zi,:) = cos(theta_t).*thy1 + sin(theta_t).*phy1;
% vtz(xi,yi,zi,:) = cos(theta_t).*thz1 + sin(theta_t).*phz1;
%         
%         end
%     end
% end

% nt = 3;
% quiver3(rx, ry, rz, vtx(:,:,:,nt), vty(:,:,:,nt), vtz(:,:,:,nt));

rvect = [rx(p,q,s), ry(p,q,s), rz(p,q,s)];
rmult = repmat(rvect, [ntheta 1]);

hold on;
quiver3( rmult(:,1), rmult(:,2), rmult(:,3), squeeze(vtx(p,q,s,:)), ...
    squeeze(vty(p,q,s,:)), squeeze(vtz(p,q,s,:)));

hold on;
plot3( rvect(1), rvect(2), rvect(3),'r*','LineWidth', 4);

hold on;
plot3( 0, 0, 0,'c*','LineWidth', 4);

